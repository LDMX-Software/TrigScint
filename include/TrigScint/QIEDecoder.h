#ifndef TRIGSCINT_QIEDECODER_H
#define TRIGSCINT_QIEDECODER_H

#include <iostream>
#include <fstream>
                                             
#include "Framework/Configure/Parameters.h"  // Needed to import parameters from configuration file
#include "Framework/Event.h"
#include "Framework/EventProcessor.h"  //Needed to declare processor
#include "TrigScint/Event/TrigScintQIEDigis.h"
#include "TrigScint/Event/QIEStream.h"

#include <TTimeStamp.h>
namespace trigscint {


  class QIEDecoder : public framework::Producer{
 public:
  QIEDecoder(const std::string &name, framework::Process &process)
      : Producer(name, process) {}

  /**
   * Default destructor, closes up boost archive and input stream
   */
  virtual ~QIEDecoder() = default;

	/**
	 * Configure our converter based off the configuration parameters
	 * decoded from the passed python script
	 */
	virtual void configure(framework::config::Parameters &ps);

  virtual void produce(framework::Event &event);



  virtual void onProcessStart();

  virtual void onProcessEnd();

  private:

	/// the channel mapping
	std::string channelMapFileName_;
	std::ifstream channelMapFile_;
	std::map< int, int> channelMap_;

	// input/output collection and pass name 
	std::string inputCollection_;
	std::string outputCollection_;
	std::string inputPassName_;
	
	// verbosity for very specific printouts that don't play well with logger format
	bool verbose_{false};

	// number of channels in the pad 
	int nChannels_{50};
	// number of time samples making up the event
	int nSamples_{5};
	// configurable flag, to set the isRealData bit in the event header
	bool isRealData_{false};

  }; //encoder

}//namespace trigscint 
                                                                                                 

#endif // TRIGSCINT_QIEDECODER_H
