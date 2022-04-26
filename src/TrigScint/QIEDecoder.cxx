#include <bitset>
#include <iostream>
#include <fstream>

#include <TMath.h>
#include <TTimeStamp.h>

#include "Framework/Configure/Parameters.h"  // Needed to import parameters from configuration file
#include "Framework/Event.h"
#include "Framework/EventProcessor.h"  //Needed to declare processor
#include "Packing/Utility/Reader.h"
#include "TrigScint/Event/TrigScintQIEDigis.h"
#include "TrigScint/Event/QIEStream.h"

namespace trigscint {

class QIEDecoder : public framework::Producer {
 public:
  QIEDecoder(const std::string &name, framework::Process &process)
      : Producer(name, process) {}

  /**
   * Default destructor, closes up boost archive and input stream
   */
  ~QIEDecoder() = default;

  /**
   * Configure our converter based off the configuration parameters
   * decoded from the passed python script
   */
  virtual void configure(framework::config::Parameters &ps) final override;

  virtual void produce(framework::Event &event) final override;

  virtual void onProcessStart() final override;

  virtual void onProcessEnd() final override;

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
  // number of time samples making up the event
  int nSamples_{5};
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
  nSamples_ = ps.getParameter<int>("number_time_samples");
  isRealData_ = ps.getParameter<bool>("is_real_data");
  verbose_ = ps.getParameter<bool>("verbose");
  std::string input_file{ps.getParameter<std::string>("input_file")};

  ldmx_log(debug) << "In configure, got parameters:" <<
    "\noutput_collection = " << outputCollection_ <<
    "\ninput_file = " << input_file <<
    "\nchannel_map_file = " <<  channelMapFileName_ <<
    "\nnumber_channels  = " <<  nChannels_ <<
    "\nnumber_time_samples  = " <<  nSamples_ <<
    "\nis_real_data  = " <<  isRealData_ <<
    "\nverbose          = " << verbose_ ;

  channelMapFile_.open(channelMapFileName_, std::ios::in);
  if (! channelMapFile_.is_open() ) {
    EXCEPTION_RAISE("BadMapFile",   "The channel mapping file cannot be opened."); // <-- appears this needs implementing first 
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
  if (elID != nChannels_ - 1 )
    ldmx_log(fatal) << "The set number of channels " << nChannels_ << " seems not to match the number from the map (+1) :" << elID; 

  reader_.open(input_file);
  return;
}  // configure

/**
 * Similar to qie_frame from legacy python script
 */
struct TimeSample {
  std::vector<uint8_t> adcs;
  std::vector<uint8_t> tdcs;
  int capid, ce, bc0;

  /**
   * Parse the input stream of bytes into this time sample
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
    if (buffer.size() != 11) {
      std::cerr << "TimeSample : was expecting 11 bytes, got " << buffer.size() << std::endl;
      return;
    }
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

class EventPacket {
  uint32_t time_since_spill_;
  /// no longer part of event info, set to 0 to keep compatibility
  uint32_t time_stamp_{0}, time_stamp_tick_{0};
  /// error codes, calculated during decoding
  bool crc0_error_{false}, crc1_error_{false}, cid_unsync_{false}, cid_skip_{false};
  /**
   * adc readout map in event
   *
   * key: int electronics ID of channel
   * val: vector of adc indexed by time sample
   */
  std::map<int, std::vector<int>> adc_map_;

  /**
   * tdc readout map in event
   *
   * key: int electronics ID of channel
   * val: vector of tdc indexed by time sample
   */
  std::map<int, std::vector<int>> tdc_map_;

  /**
   * clean and split the raw fiber data into time samples
   *
   * The fiber data is "contaminated" with different signal phrases.
   * - 0xfbf7 : needs to be removed
   * - 0x7cbc : needs to be removed
   * - 0xbc   : separates time samples within a fiber
   *
   */
  std::vector<TimeSample> extract_timesamples(const std::vector<uint32_t>& stream) {
    // clean
    std::vector<uint16_t> cleaned;
    cleaned.reserve(stream.size()*2);
    for (const uint32_t& word : stream) {
      for (std::size_t i{0}; i < 2; i++) {
        uint16_t subword = ((word >> 16*i) & 0xffff);
        if (subword != 0xfbf7 and subword != 0x7cbc) 
          cleaned.push_back(subword);
      }
    }

    // split
    std::vector<std::vector<uint8_t>> timesamples;
    for (const uint16_t& word : cleaned) {
      for (std::size_t i{0}; i < 2; i++) {
        uint8_t byte = ((word >> 8*i) & 0xff);
        if (byte == 0xbc) {
          // new timesample
          timesamples.emplace_back();
        } else {
          timesamples.back().push_back(byte);
        }
      }
    }

    // decode
    std::vector<TimeSample> decoded_ts;
    for (const auto& raw_ts : timesamples) {
      decoded_ts.emplace_back(raw_ts);
    }

    return decoded_ts;
  }

 public:
  uint32_t getTimeSinceSpill() const { return time_since_spill_; }
  uint32_t getTimeStamp() const { return time_stamp_; }
  uint32_t getTimeStampTick() const { return time_stamp_tick_; }
  bool crc0Error() const { return crc0_error_; }
  bool crc1Error() const { return crc1_error_; }
  bool cidUnsync() const { return cid_unsync_; }
  bool cidSkip() const { return cid_skip_; }
  const std::map<int,std::vector<int>>& adcs() const {
    return adc_map_;
  }
  const std::vector<int>& adcs(int elec_id) const {
    return adcs().at(elec_id);
  }
  const std::map<int,std::vector<int>>& tdcs() const {
    return tdc_map_;
  }
  const std::vector<int>& tdcs(int elec_id) const {
    return tdcs().at(elec_id);
  }
  /**
   * Function called by Reader when streaming
   *
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
   * Look at extract_timesamples and TimeSample for how to decode the
   * two "columns" of raw fiber data.
   */
  packing::utility::Reader& read(packing::utility::Reader& r) {
    static uint64_t w; // dummy word for reading

    // first word is timestamp word
    r >> w;
    time_since_spill_ = w & 0xffffffff;

    // the rest of the words are split across fibers
    std::vector<uint32_t> fiber_1_raw, fiber_2_raw;
    while (r and w != 0xffffffffffffffff) {
      r >> w;
      fiber_1_raw.push_back(w & 0xffffffff);
      fiber_2_raw.push_back((w >> 32) & 0xffffffff);
    }

    // clean_kchar, clean_BC7C, split across 'BC' for both fibers
    std::vector<TimeSample> fiber_1{extract_timesamples(fiber_1_raw)},
                            fiber_2{extract_timesamples(fiber_2_raw)};

    if (fiber_1.size() != fiber_2.size()) {
      std::cout << "Non-matching number of time samples :( " << std::endl;
      return r;
    }

    adc_map_.clear();
    tdc_map_.clear();
    for (std::size_t i_ts{0}; i_ts < fiber_1.size(); i_ts++) {
      if (fiber_1.at(i_ts).capid != fiber_2.at(i_ts).capid) cid_unsync_ = true;
      if (fiber_1.at(i_ts).ce != 0) crc0_error_ = true;
      if (fiber_2.at(i_ts).ce != 0) crc1_error_ = true;
      if (i_ts > 0) {
        // logic needs some work, what happens if corrupt time sample word in middle?
        if ((fiber_1.at(i_ts-1).capid+1)%4 != fiber_1.at(i_ts).capid%4) {
          std::cout << "Found CIDSkip in fiber 1" << std::endl;
          cid_skip_ = true;
        }
        if ((fiber_2.at(i_ts-1).capid+1)%4 != fiber_2.at(i_ts).capid%4) {
          std::cout << "Found CIDSkip in fiber 1" << std::endl;
          cid_skip_ = true;
        }
      } // check if corrupted word
      for (std::size_t i_c{0}; i_c < 8; i_c++) {
        adc_map_[i_c  ].push_back(fiber_1.at(i_ts).adcs.at(i_c));
        adc_map_[i_c+8].push_back(fiber_2.at(i_ts).adcs.at(i_c));
        tdc_map_[i_c  ].push_back(fiber_1.at(i_ts).tdcs.at(i_c));
        tdc_map_[i_c+8].push_back(fiber_2.at(i_ts).tdcs.at(i_c));
      } // loop over channels in a fiber
    } // loop over time samples

    return r;
  } // read
}; // EventPacket

void QIEDecoder::produce(framework::Event &event) {
  // dummy packet for transient decoding
  static EventPacket event_packet;

  ldmx_log(debug) << "QIEDecoder: produce() starts! Event number: "
    << event.getEventHeader().getEventNumber();

  ldmx_log(debug) << "Reading next event packet";
  if (!(reader_ >> event_packet)) {
    ldmx_log(debug) << "no more events";
    return; //abortEvent();
  }

  //  ldmx::EventHeader *eh = (ldmx::EventHeader*)event.getEventHeaderPtr();
  //eh->setIsRealData(isRealData_); //doesn't help
  
  /*   // comment for now, requires pushing change ot EventHeader
  event.getEventHeader().setIsRealData(isRealData_);
  ldmx_log(debug) << "Decoder bool isRealData = " << isRealData_ 
    << " and after attempt of setting it, event header returns " <<  event.getEventHeader().isRealData();
  */
    //turns out this need to be configurable for now, to read real data
  int nSamp = nSamples_; //QIEStream::NUM_SAMPLES ;  
  ldmx_log(debug) << "num samples = " << nSamp;
  
  uint32_t timeEpoch = event_packet.getTimeStamp();
  uint32_t timeClock = event_packet.getTimeStampTick();
  uint32_t timeSpill = event_packet.getTimeSinceSpill();
  ldmx_log(debug) << "time stamp words are : " << timeEpoch 
    << " (" << std::bitset<64>(timeEpoch) << ") and " 
    << timeClock << " (" << std::bitset<64>(timeClock) << ") clock ticks, and " 
    << timeSpill << " (" << std::bitset<64>(timeSpill) << ", or, " 
    << std::hex << timeSpill << std::dec << ") counts since start of spill";

  int sigBitsSkip=6; //the first 6 bits are part of something else. 
  int divisor=TMath::Power(2, 32-sigBitsSkip);
  timeSpill=timeSpill%divisor; //remove them by taking remainder in division by the values of the last skipped bit
  //timeSpill=timeSpill/spillTimeConv_;
  ldmx_log(debug) << "After taking it mod 2^" << 32-sigBitsSkip << " (which is " << divisor << ", spill time is " << timeSpill;
  event.getEventHeader().setIntParameter("timeSinceSpill", timeSpill);
  //  event.getEventHeader().setTimeSinceSpill(timeSpill); //not working 

  TTimeStamp *timeStamp = new TTimeStamp(timeEpoch);
  //  timeStamp->SetNanoSec(timeClock); //not sure about this one...
  event.getEventHeader().setTimestamp(*timeStamp);

  /* -- TS event header done; read the channel contents -- */

  ldmx_log(debug) << "Done reading in header, ADC and TDC for event " << event.getEventNumber();
  std::vector<trigscint::TrigScintQIEDigis> outDigis;
  for (const auto& [ elec_id, adcs ] : event_packet.adcs()) {
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
    digi.setTDC(event_packet.tdcs(elec_id));
    digi.setTimeSinceSpill(timeSpill);
    if (bar == 0)
      ldmx_log(debug) << "for bar 0, got time since spill "<< digi.getTimeSinceSpill();
    outDigis.push_back( digi );
    ldmx_log(debug) << "Iterator points to key " << elec_id
            << " and mapped channel supposedly  is " << bar ;
    ldmx_log(debug) << "Made digi with elecID = " << digi.getElecID()
            << ", barID = " << digi.getChanID()
            << ", third adc value " << digi.getADC().at(2)
            << " and third tdc " << digi.getTDC().at(2) ;
  }
 
  event.add(outputCollection_, outDigis);
} // produce

void QIEDecoder::onProcessStart() {
  ldmx_log(debug) << "Process starts!";
  return;
}

void QIEDecoder::onProcessEnd() {
  ldmx_log(debug) << "Process ends!";
  return;
}

}  // namespace trigscint                                                                                                                           

DECLARE_PRODUCER_NS(trigscint, QIEDecoder);

