#include "TrigScint/TrigScintRecHitProducer.h"
#include "Framework/Exception/Exception.h"
#include "Framework/RandomNumberSeedService.h"

#include <iostream>

namespace trigscint {

TrigScintRecHitProducer::TrigScintRecHitProducer(const std::string &name,
                                                 framework::Process &process)
    : Producer(name, process) {}

TrigScintRecHitProducer::~TrigScintRecHitProducer() {}

void TrigScintRecHitProducer::configure(
    framework::config::Parameters &parameters) {
  // Configure this instance of the producer
  pedestal_ = parameters.getParameter<double>("pedestal");
  gain_ = parameters.getParameter<double>("gain");
  mevPerMip_ = parameters.getParameter<double>("mev_per_mip");
  pePerMip_ = parameters.getParameter<double>("pe_per_mip");
  inputCollection_ = parameters.getParameter<std::string>("input_collection");
  inputPassName_ = parameters.getParameter<std::string>("input_pass_name");
  outputCollection_ = parameters.getParameter<std::string>("output_collection");
  verbose_ = parameters.getParameter<bool>("verbose");
  sample_of_interest_ = parameters.getParameter<int>("sample_of_interest");
}

void TrigScintRecHitProducer::produce(framework::Event &event) {
  // initialize QIE object for linearizing ADCs
  SimQIE qie;

  // Ensure the sample of interest <4
  if(sample_of_interest_>3) {
    ldmx_log(error)<<"sample_of_interest_ should be one of 0,1,2,3\n"
		   <<"Currently, sample_of_interest = "<<sample_of_interest_
		   <<"\n";
    return;
  }

  // looper over sim hits and aggregate energy depositions
  // for each detID
  const auto digis{event.getCollection<trigscint::TrigScintQIEDigis>(
      inputCollection_, inputPassName_)};

  std::vector<ldmx::TrigScintHit> trigScintHits;
  for (const auto &digi : digis) {
    ldmx::TrigScintHit hit;
    auto adc{digi.getADC()};
    auto tdc{digi.getTDC()};

    hit.setModuleID(0);
    hit.setBarID(digi.getChanID());
    hit.setBeamEfrac(-1.);

    hit.setAmplitude(qie.ADC2Q(adc[sample_of_interest_]) +
                     qie.ADC2Q(adc[sample_of_interest_+1]));  // femptocoulombs

    if (tdc[sample_of_interest_] > 49)
      hit.setTime(-999.);
    else
      hit.setTime(tdc[sample_of_interest_] * 0.5);

    hit.setEnergy((qie.ADC2Q(adc[sample_of_interest_])+
		   qie.ADC2Q(adc[sample_of_interest_+1])-pedestal_)
		  *6250. / gain_ * mevPerMip_ / pePerMip_);  // MeV
    hit.setPE((qie.ADC2Q(adc[sample_of_interest_])+
	       qie.ADC2Q(adc[sample_of_interest_+1]) - pedestal_)
	      *6250. / gain_);
    trigScintHits.push_back(hit);
  }
  // Create the container to hold the
  // digitized trigger scintillator hits.

  event.add(outputCollection_, trigScintHits);
}
}  // namespace trigscint

DECLARE_PRODUCER_NS(trigscint, TrigScintRecHitProducer);
