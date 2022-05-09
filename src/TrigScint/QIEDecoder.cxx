#include <bitset>
#include <iostream>
#include <fstream>

#include <TMath.h>
#include <TTimeStamp.h>

#include "Framework/Configure/Parameters.h"  // Needed to import parameters from configuration file
#include "Framework/Event.h"
#include "Framework/EventProcessor.h"  //Needed to declare processor
#include "Packing/Utility/Reader.h"
#include "Packing/Utility/Hex.h"
#include "Packing/Utility/Mask.h"
#include "TrigScint/Event/TrigScintQIEDigis.h"
#include "TrigScint/Event/QIEStream.h"

namespace trigscint {

/**
 * Smallest "chunk" of data coming from TS frontend
 *
 * Represents a single sample of data collected from the
 * scintillator bars.
 *
 * Similar to qie_frame from legacy python script
 */
struct TimeSample {
  std::vector<uint8_t> adcs;
  std::vector<uint8_t> tdcs;
  int capid, ce, bc0;

  /**
   * Parse the input stream of bytes into this time sample
   *
   * @note The input buffer needs to be _exactly_ 11 bytes.
   *
   * ## Format
   * Each row is a BYTE (two hex digits, 8 bits)
   *
   *  | ???? | capid (2 bits) | ce (1 bit) | bc0 (1 bit) |
   *  | adc0 |
   *  | adc1 |
   *  | .... |
   *  | adc8 |
   *  | letdc0 | letdc1 | letdc2 | letdc3 |
   *  | letdc4 | letdc5 | letdc6 | letdc7 |
   *
   * this is LETDC; only the two most significant bits included
   *    they are shipped as least significant bits --> shift them 
   *    want LE TDC = 3 to correspond to 64 > 49 (which is maxTDC in sim)
   */
  TimeSample(const std::vector<uint8_t>& buffer) {
    capid = (buffer.at(0) >> 2) & 0x3;
    ce    = (buffer.at(0) >> 1) & 0x1;
    bc0   = (buffer.at(0) >> 0) & 0x1;
    for (std::size_t i_channel{0}; i_channel < 8; i_channel++) {
      adcs.push_back(buffer.at(1+i_channel));
      uint8_t tdc_pack = buffer.at(9 + 1*(i_channel>3));
      uint8_t letdc = ((tdc_pack >> (2*(i_channel%4))) & 0xff);
      tdcs.push_back(16*(letdc+1));
    }
  } // decoding constructor

};  // TimeSample

class QIEDecoder : public framework::Producer {
 public:
  QIEDecoder(const std::string &name, framework::Process &process)
      : Producer(name, process) {}
  ~QIEDecoder() = default;
  virtual void configure(framework::config::Parameters &ps) final override;
  virtual void produce(framework::Event &event) final override;
 private:
  std::vector<TimeSample> extract_timesamples(const std::vector<uint32_t>& stream);
 private:
  /// the channel mapping
  std::string channelMapFileName_;
  std::ifstream channelMapFile_;
  std::map< int, int> channelMap_;

  // output collection name
  std::string outputCollection_;
  
  // verbosity for very specific printouts that don't play well with logger format
  bool verbose_{false};

  // number of channels in the pad 
  int nChannels_{50};
  // configurable flag, to set the isRealData bit in the event header
  bool isRealData_{false};
  /// binary file reader
  packing::utility::Reader reader_;
}; // QIEDecoder

void QIEDecoder::configure(framework::config::Parameters &ps) {
  // Configure this instance of the encoder
  outputCollection_ = ps.getParameter<std::string>("output_collection");
  channelMapFileName_ = ps.getParameter<std::string>("channel_map_file");
  nChannels_ = ps.getParameter<int>("number_channels");
  isRealData_ = ps.getParameter<bool>("is_real_data");
  verbose_ = ps.getParameter<bool>("verbose");
  std::string input_file{ps.getParameter<std::string>("input_file")};

  ldmx_log(debug) << "In configure, got parameters:" <<
    "\noutput_collection = " << outputCollection_ <<
    "\ninput_file = " << input_file <<
    "\nchannel_map_file = " <<  channelMapFileName_ <<
    "\nnumber_channels  = " <<  nChannels_ <<
    "\nis_real_data  = " <<  isRealData_ <<
    "\nverbose          = " << verbose_ ;

  channelMapFile_.open(channelMapFileName_, std::ios::in);
  if (! channelMapFile_.is_open() ) {
    // appears that the exception is not being thrown
    EXCEPTION_RAISE("BadMapFile", "The channel mapping file cannot be opened.");
    ldmx_log(fatal) <<  "The channel mapping file cannot be opened.";
    return;
  }
  int chID, elID;
  while (! channelMapFile_.eof() ){
    channelMapFile_ >> elID >> chID ;
    // make the map based on electronics ID, to look up channel ID.
    // the reason is, we will only know the elecID from the position
    // of the word in the stream. so these need to be strictly ordered.
    // writing the right bar in the right position is easier if we can just
    // read this map from beginning to end.
    // barID can always be set, or looked up, as a property of the digi.
  
    // here make the elecID the key (other way around when encoding)
    channelMap_.insert(std::pair<int,int>(elID, chID));
    ldmx_log(debug) << elID << "  chID " << chID ;
  }
  channelMapFile_.close();
  if (elID != nChannels_ - 1)
    ldmx_log(fatal) << "The set number of channels " << nChannels_ 
      << " seems not to match the number from the map (+1) :" << elID; 

  reader_.open(input_file);
  return;
}  // configure

/**
 * clean and split the raw fiber data into time samples
 *
 * The fiber data is "contaminated" with different signal phrases.
 * - 0xfbf7 : needs to be removed
 * - 0x7cbc : needs to be removed
 * - 0xbc   : separates time samples within a fiber
 *
 * After this cleaning and splitting, we provide the encoded time sample
 * bytes to the TimeSample constructor to decode into data.
 */
std::vector<TimeSample> QIEDecoder::extract_timesamples(const std::vector<uint32_t>& stream) {
  // clean
  std::vector<uint16_t> cleaned;
  for (const uint32_t& word : stream) {
    for (std::size_t i{0}; i < 2; i++) {
      uint16_t subword = ((word >> 16*i) & 0xffff);
      if (subword != 0xf7fb and subword != 0xbc7c) 
        cleaned.push_back(subword);
    }
  }

  // split
  std::vector<std::vector<uint8_t>> timesamples;
  timesamples.emplace_back(); // need first one
  for (const uint16_t& word : cleaned) {
    uint8_t b0 = ((word >> 8) & 0xff);
    uint8_t b1 = ((word >> 0) & 0xff);
    for (uint8_t byte : {b0,b1}) {
      if (byte == 0xbc) {
        // new timesample
        ldmx_log(debug) << "New Timesample, previous one has " 
          << timesamples.back().size() << " bytes";
        timesamples.emplace_back();
      } else {
        timesamples.back().push_back(byte);
      }
    }
  }
  ldmx_log(debug) << "Last time sample has " << timesamples.back().size() << " bytes";

  // decode
  std::vector<TimeSample> decoded_ts;
  for (const auto& raw_ts : timesamples) {
    if (raw_ts.size() == 11) {
      decoded_ts.emplace_back(raw_ts);
    } else {
      ldmx_log(debug) << "Received " << raw_ts.size() << " != 11 bytes for a time sample";
    }
  }

  return decoded_ts;
}

void QIEDecoder::produce(framework::Event &event) {
  static uint64_t w; // dummy word for decoding

  ldmx_log(debug) << "QIEDecoder: produce() starts! Event number: "
    << event.getEventHeader().getEventNumber();

  if (event.getEventNumber() < 2) {
    ldmx_log(debug) << "First event, throwing away words until first event signal word.";
    w = 0;
    while (reader_ and w != 0xffffffffffffffff) reader_ >> w;
  }

  /**
   * The next 64-bit word has the timestamp burried within it.
   *
   * The second half of the next 64-bits contains the timestamp
   * but with endian-ness flopped so we need to read in the bytes
   * in reverse order.
   */
  uint32_t dummy;
  reader_ >> dummy;
  uint8_t ts_bytes[4];
  reader_ >> ts_bytes[3] >> ts_bytes[2] >> ts_bytes[1] >> ts_bytes[0];
  uint32_t *ts = reinterpret_cast<uint32_t*>(ts_bytes);
  ldmx_log(debug) << "Header timestamp word: " << packing::utility::hex<uint32_t>(*ts);
  uint32_t timeSpill = (*ts & packing::utility::mask<26>);
  /// no longer part of event info, set to 0 to keep compatibility
  uint32_t timeEpoch = 0; // UTC seconds
  uint32_t timeClock = 0; // ns clock ticks since last whole second
  ldmx_log(debug) << "time stamp words are : " << timeEpoch 
    << " (" << std::bitset<32>(timeEpoch) << ") and " 
    << timeClock << " (" << std::bitset<32>(timeClock) << ") clock ticks, and " 
    << timeSpill << " (" << std::bitset<32>(timeSpill) << ", or, " 
    << std::hex << timeSpill << std::dec << ") counts since start of spill";
  event.getEventHeader().setIntParameter("timeSinceSpill", timeSpill);
  //  event.getEventHeader().setTimeSinceSpill(timeSpill); //not working 

  TTimeStamp *timeStamp = new TTimeStamp(timeEpoch);
  //  timeStamp->SetNanoSec(timeClock); //not sure about this one...
  event.getEventHeader().setTimestamp(*timeStamp);

  //  ldmx::EventHeader *eh = (ldmx::EventHeader*)event.getEventHeaderPtr();
  //eh->setIsRealData(isRealData_); //doesn't help
  
  /*   // comment for now, requires pushing change ot EventHeader
  event.getEventHeader().setIsRealData(isRealData_);
  ldmx_log(debug) << "Decoder bool isRealData = " << isRealData_ 
    << " and after attempt of setting it, event header returns " 
    <<  event.getEventHeader().isRealData();
  */

  ldmx_log(debug) << "Done with timestamp header, moving onto raw fiber data.";

  /**
   * We read 64-bit words until the boundary word (all F's) is encountered
   *
   * We assume the last boundary word was read by the last event.
   * @note This means the first event needs to handle initializing the reader.
   *
   * ## Format Notes
   * Each row is a 64-bit word.
   *
   *  |    ?nonsense?     | Time Since Spill (32 bits) |
   *  | Fiber 1 (32 bits) | Fiber 2 (32 bits) |
   *  |  ...continue...   |  ...continue...   |
   *  | End of Event Signal Word (All 1s)     |
   *
   * Look at TimeSample for how to decode the
   * two "columns" of raw fiber data.
   * The fiber data needs to be byte flopped like the time
   * stamp word above.
   */
  std::vector<uint32_t> fiber_1_raw, fiber_2_raw;
  w = 0;
  while (reader_ and w != 0xffffffffffffffff) {
    reader_ >> w;
    uint8_t *bytes = reinterpret_cast<uint8_t*>(&w);
    fiber_1_raw.push_back(
        (bytes[4]<<24)+(bytes[5]<<16)+(bytes[6]<<8)+bytes[7]);
    fiber_2_raw.push_back(
        (bytes[0]<<24)+(bytes[1]<<16)+(bytes[2]<<8)+bytes[3]);
  }

  ldmx_log(debug) << "Done reading raw fiber data "
    << fiber_1_raw.size() << " fiber 1 words and "
    << fiber_2_raw.size() << " fiber 2 words";
  ldmx_log(debug) << "extracting time samples from them";

  // clean_kchar, clean_BC7C, split across 'BC' for both fibers
  std::vector<TimeSample> fiber_1{extract_timesamples(fiber_1_raw)},
                          fiber_2{extract_timesamples(fiber_2_raw)};

  if (fiber_1.size() != fiber_2.size()) {
    ldmx_log(fatal) << "Non-matching number of time samples between fibers!";
    return;
  }

  ldmx_log(debug) << fiber_1.size() << " time samples in this event.";

  /// error codes, calculated during decoding
  bool crc1_error{false}, crc2_error{false}, cid_unsync{false}, cid_skip{false};
  /**
   * adc and tdc readout map in event
   *
   * Resort the samples into by-channel rather than
   * by time sample. While sorting, we also check for
   * various error codes.
   *
   * key: int electronics ID of channel
   * val: vector of adc/tdc indexed by time sample
   */
  std::map<int, std::vector<int>> adc_map, tdc_map;
  for (std::size_t i_ts{0}; i_ts < fiber_1.size(); i_ts++) {
    if (fiber_1.at(i_ts).capid != fiber_2.at(i_ts).capid) cid_unsync = true;
    if (fiber_1.at(i_ts).ce != 0) crc1_error = true;
    if (fiber_2.at(i_ts).ce != 0) crc2_error = true;
    if (i_ts > 0) {
      // logic needs some work, what happens if corrupt time sample word in middle?
      if ((fiber_1.at(i_ts-1).capid+1)%4 != fiber_1.at(i_ts).capid%4) {
        ldmx_log(debug) << "Found CIDSkip in fiber 1";
        cid_skip = true;
      }
      if ((fiber_2.at(i_ts-1).capid+1)%4 != fiber_2.at(i_ts).capid%4) {
        ldmx_log(debug) << "Found CIDSkip in fiber 2";
        cid_skip = true;
      }
    } // check if corrupted word
    for (std::size_t i_c{0}; i_c < 8; i_c++) {
      adc_map[i_c  ].push_back(fiber_1.at(i_ts).adcs.at(i_c));
      adc_map[i_c+8].push_back(fiber_2.at(i_ts).adcs.at(i_c));
      tdc_map[i_c  ].push_back(fiber_1.at(i_ts).tdcs.at(i_c));
      tdc_map[i_c+8].push_back(fiber_2.at(i_ts).tdcs.at(i_c));
    } // loop over channels in a fiber
  } // loop over time samples

  if (crc1_error) ldmx_log(debug) << "CRC Error on fiber 1"; 
  if (crc2_error) ldmx_log(debug) << "CRC Error on fiber 2";
  if (cid_unsync) ldmx_log(debug) << "CID unsync between fibers";
  if (cid_skip  ) ldmx_log(warn) << "At least one CIDSkip found in event";

  /**
   * Now that the digis are sorted correctly,
   * we can translate the electronics IDs into bar numbers
   * put the information into our C++ class
   * and then ship them out
   */
  std::vector<trigscint::TrigScintQIEDigis> outDigis;
  for (const auto& [ elec_id, adcs ] : adc_map) {
    TrigScintQIEDigis  digi;
    digi.setADC(adcs);
    if (channelMap_.find(elec_id) == channelMap_.end()) {
      ldmx_log(fatal) << "Couldn't find the bar ID corresponding to electronics ID " 
        << elec_id << "!! Skipping." ;
      continue;
    }
    int bar = channelMap_[elec_id];
    digi.setElecID(elec_id);
    digi.setChanID( bar );
    digi.setTDC(tdc_map.at(elec_id));
    digi.setTimeSinceSpill(timeSpill);
    if (bar == 0)
      ldmx_log(debug) << "for bar 0, got time since spill "<< digi.getTimeSinceSpill();
    outDigis.push_back( digi );
    ldmx_log(debug) << "Made digi with elecID = " << digi.getElecID()
            << ", barID = " << digi.getChanID()
            << ", third adc value " << digi.getADC().at(2)
            << " and third tdc " << digi.getTDC().at(2) ;
  }
 
  event.add(outputCollection_, outDigis);
} // produce

}  // namespace trigscint                                                                                                                           
DECLARE_PRODUCER_NS(trigscint, QIEDecoder);

