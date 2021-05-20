#include "TrigScint/NumericalRecHitProducer.h"
#include "Framework/Exception/Exception.h"
#include "Framework/RandomNumberSeedService.h"
#include "TMath.h"

#include <iostream>
#include "TLinearFitter.h"
#include<iomanip>

namespace trigscint {

NumericalRecHitProducer::NumericalRecHitProducer(const std::string &name,
                                                 framework::Process &process)
    : Producer(name, process) {}

NumericalRecHitProducer::~NumericalRecHitProducer() {}

void NumericalRecHitProducer::configure(
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
  tdc_thr_ = parameters.getParameter<double>("tdc_thr");
  qie_sf_ = parameters.getParameter<double>("qie_sf");

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

void NumericalRecHitProducer::produce(framework::Event &event) {
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

    auto Charge =
      ChargeReconstruction(adc,tdc,sample_of_interest_);

    hit.setAmplitude(Charge);
    hit.setEnergy((Charge - pedestal_) * 6250. /
		  gain_ * mevPerMip_ / pePerMip_);  // MeV
    hit.setPE((Charge - pedestal_) * 6250. /gain_);
    trigScintHits.push_back(hit);
  }
  // Create the container to hold the
  // digitized trigger scintillator hits.

  event.add(outputCollection_, trigScintHits);
}
  Double_t NumericalRecHitProducer::ChargeReconstruction
  (std::vector<int>adc,std::vector<int>tdc,int sample) {
    int npulses = 0;              // No. of true pulses
    int poi=0;                    // The pulse of interest
    std::vector<float> Charge_;	  // stores pulse amplitudes
    auto Qdata = new Double_t[5]; // Linearized charge
    float tend = 1000/qie_sf_;	  // 1 time sample (in ns)
    float k_ = pulse_params_[0];
    float tmax_ = pulse_params_[1];
    float par0 = (exp(k_*tmax_)-1)/(k_*tmax_)*exp(-k_*tend);
    SimQIE qie;
    auto pulse = new Expo(k_,tmax_);
    
    for(int i=0;i<tdc.size();i++) {
      Qdata[i] = qie.ADC2Q(adc[i]);
      if(tdc[i]<50) {
	if(i==sample_of_interest_)
	  poi=npulses;
	double CostFunction(double* params);
	npulses++;
      }
    }
  /////////////// For Debigging purposes
  if(verbose_) {
    std::cout<<"TS \t|\t0\t|\t1\t|\t2\t|\t3\t|\t4\t|\n"
	     <<"---------------------------------------------"
	     <<"--------------------------------------------\n"
	     <<"tdc \t|";
    for(int i=0;i<5;i++)
      std::cout<<std::setw(10)<<tdc[i]<<"\t|";
    
    std::cout<<"\nadc \t|";
    for(int i=0;i<5;i++)
      std::cout<<std::setw(10)<<adc[i]<<"\t|";
    
    std::cout<<"\nQdata\t|";
    for(int i=0;i<5;i++)
      std::cout<<std::setw(10)<<Qdata[i]<<"\t|";

    // std::cout<<"\nQErr\t|";
    // for(int i=0;i<5;i++)
    //   std::cout<<std::setw(10)<<QErr[i]<<"\t|";

    std::cout<<"\n---------------------------------------------"
	     <<"--------------------------------------------";
    // for(int n = 0;n<npulses;n++){
    //   std::cout<<std::setw(10)<<"\nPulse"<<n<<"\t|";
    // for(int i=0;i<5;i++)
    // 	std::cout<<std::setw(10)<<params(1+n)*Qdep[n*5+i]<<"\t|";
    // }
    
    std::cout<<"\n"
	     <<"\nnpulses = "<<npulses<<std::endl
	     <<"poi = "<<poi<<std::endl;
	     // <<"chisquare = "<<chisquare<<std::endl
	     // <<"params:\n";
    
    // params.Print();
    // std::cout<<"errors:\n";
    // errors.Print();

    std::cout<<"Charge_[]: ";
    for(int i=0;i<Charge_.size();i++)
      std::cout<<std::setw(10)<<Charge_[i]<<"\t|";
  }
  std::cout<<"\n";

  return Charge_[poi];
  }

  double CostFunction(double* params) {
    auto pulse = new Expo(pulse_params_[0],pulse_params_[1]);
    pulse->AddPulse(params[0],params[1]);

    double Qpred = pulse->Integrate(0,1000/qie_sf);
    double Tpred = qie.TDC(pulse,T0)/2;
    
    double cost1 = pow(,2);
    double cost2 = pow(,2);
    return 0;
  }
}  // namespace trigscint

DECLARE_PRODUCER_NS(trigscint, NumericalRecHitProducer);
