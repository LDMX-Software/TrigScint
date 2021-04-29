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
  En_Reco_Option_ = parameters.getParameter<int>("En_Reco_Option");

  input_pulse_shape_ =
      parameters.getParameter<std::string>("input_pulse_shape");
  if (input_pulse_shape_ == "Expo") {
    pulse_params_.clear();
    pulse_params_.push_back(parameters.getParameter<double>("expo_k"));
    pulse_params_.push_back(parameters.getParameter<double>("expo_tmax"));

    ldmx_log(debug) << "expo_k =" << pulse_params_[0];
    ldmx_log(debug) << "expo_tmax =" << pulse_params_[1];
  }

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


    if (tdc[sample_of_interest_] > 49)
      hit.setTime(-999.);
    else
      hit.setTime(tdc[sample_of_interest_] * 0.5);

    if(En_Reco_Option_==1) {
      auto Charge = ChargeReconstruction(adc,tdc);
      hit.setAmplitude(Charge);
      hit.setEnergy((Charge - pedestal_) * 6250. /
		    gain_ * mevPerMip_ / pePerMip_);  // MeV
      hit.setPE((Charge - pedestal_) * 6250. /gain_);
    }
    else {
      // femptocoulombs
      hit.setAmplitude(qie.ADC2Q(adc[sample_of_interest_]) +
		       qie.ADC2Q(adc[sample_of_interest_+1]));
      // MeV
      hit.setEnergy((qie.ADC2Q(adc[sample_of_interest_])+
		     qie.ADC2Q(adc[sample_of_interest_+1])-pedestal_)
		    *6250. / gain_ * mevPerMip_ / pePerMip_);
      hit.setPE((qie.ADC2Q(adc[sample_of_interest_])+
		 qie.ADC2Q(adc[sample_of_interest_+1]) - pedestal_)
		*6250. / gain_);
      trigScintHits.push_back(hit);
    }
  }
  // Create the container to hold the
  // digitized trigger scintillator hits.

  event.add(outputCollection_, trigScintHits);
}

Double_t ChargeReconstruction(std::vector<int>adc
                              ,std::vector<int>tdc
                              ,int sample_of_interest=2) {

  std::vector<Double_t> Qdep;   // charge deposited by each pulse
  auto Qdata = new Double_t[5]; // Linearized charge
  auto QErr = new Double_t[5];  // Quantization error
  int npulses = 0;              // No. of true pulses
  TVectorD params;              // stores the reco amplitudes
  TVectorD errors;              // Associated errors
  Double_t chisquare;           // Chi-square error
  int poi=0;                    // The pulse of interest
  SimQIE qie;

  for(int i=0;i<tdc.size();i++) {
    Qdata[i] = qie.ADC2Q(adc[i]);
    QErr[i] = qie.QErr(Qdata[i]);
    if(tdc[i]<50) {
      if(i==sample_of_interest)
        poi =i;
      npulses++;
      // auto pulse = new Expo(pulse_params_[0], pulse_params_[1]);
      auto pulse = new Expo(0.1,5.0);
      pulse->AddPulse(25*i+2*tdc[i],1);
      for(float ii=0;ii<5;ii++)
        Qdep.push_back(pulse->Integrate(25*ii,25*ii+25));
    }
  }

  TLinearFitter *lf=new TLinearFitter(npulses);
  TString formula;
  formula.Form("hyp%d",npulses);
  lf->SetFormula(formula);
  lf->FixParameter(0,0);        // fix the constant (par0) to 0

  lf->AssignData(5,npulses,&Qdep[0],Qdata,QErr);
  lf->Eval();

  lf->GetParameters(params);
  lf->GetErrors(errors);
  chisquare=lf->GetChisquare();
  return params(1+poi);
}
}  // namespace trigscint

DECLARE_PRODUCER_NS(trigscint, TrigScintRecHitProducer);
