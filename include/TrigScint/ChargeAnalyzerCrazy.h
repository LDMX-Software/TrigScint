/**
 * @file ChargeAnalyzerCrazy.h
 * @brief
 * @author
 */

#ifndef TRIGSCINT_CHARGEANALYZERCRAZY_H
#define TRIGSCINT_CHARGEANALYZERCRAZY_H

//LDMX Framework                               
#include "Framework/Configure/Parameters.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "TrigScint/Event/EventReadout.h"
#include "TrigScint/Event/TrigScintQIEDigis.h"
#include "TrigScint/Event/TestBeamHit.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"

namespace trigscint {

  /**
   * @class ChargeAnayzerCrazy
   * @brief
   */
  class ChargeAnalyzerCrazy : public framework::Analyzer {
  public:

	ChargeAnalyzerCrazy(const std::string& name, framework::Process& process); // : framework::Analyzer(name, process) {}
	virtual ~ChargeAnalyzerCrazy();
    virtual void configure(framework::config::Parameters &parameters);

    virtual void analyze( const framework::Event &event) final override;

	//
	virtual void onFileOpen();

	//
	virtual void onFileClose();

    virtual void onProcessStart() final override;

    virtual void onProcessEnd() final override;

    float convertToID( float yVal ) { return (yVal+yOffset_)*yToIDfactor_; }

  private:

	std::vector <std::vector <TH1F*> > vChargeVsTime;
	

	//configurable parameters
    std::string inputCol_;
    std::string inputPassName_{""};
	std::vector<double> peds_;
	std::vector<double> gain_;
	int startSample_{0};
	int EventDisplayNumber_{0};
	float PE;

	//plotting stuff 
	int evNb;
    int nEv{200};
    int nChannels{16};
    int nTrkMax{100};
    	int fillNb{0};
    	
    	int Qlow_thr = 2.e2;
    	int Qmed_thr = 5.e3;
    	int Qmed_thr2 = 30.e3;
    	int Qhigh_thr = 10.e3;
    	int Qhigh_thr2 = 64.e3;
    	
    	int PE_low = 1;
    	int PE_med1 = 10;
    	int PE_med2 = 150;
    	int PE_high = 400;
    	
    	int Si=0;
    	int Sbar0=0;
    	int Sbar1=0;
    	
    	//std::vector<std::vector<float>> QChan0;
  	//std::vector<std::vector<float>> QChan1;
  	//std::vector<std::vector<float>> QChan2;
	//std::vector<std::vector<float>> QChan3;
	//std::vector<std::vector<float>> QChan4;
	//std::vector<std::vector<float>> QChan5;
	//std::vector<std::vector<float>> QChan6;
	//std::vector<std::vector<float>> QChan7;
	//std::vector<std::vector<float>> QChan8;
	//std::vector<std::vector<float>> QChan9;
	//std::vector<std::vector<float>> QChan10;
	//std::vector<std::vector<float>> QChan11;
	//std::vector<std::vector<float>> QChan12;
	//std::vector<std::vector<float>> QChan13;
	//std::vector<std::vector<float>> QChan14;
	//std::vector<std::vector<float>> QChan15;
	std::vector<float> QChan0;
	std::vector<float> QChan1;
	std::vector<float> QChan2;
	std::vector<float> QChan3;
	std::vector<float> QChan4;
	std::vector<float> QChan5;
	std::vector<float> QChan6;
	std::vector<float> QChan7;
	std::vector<float> QChan8;
	std::vector<float> QChan9;
	std::vector<float> QChan10;
	std::vector<float> QChan11;
	std::vector<float> QChan12;
	std::vector<float> QChan13;
	std::vector<float> QChan14;
	std::vector<float> QChan15;
	std::vector<std::vector<float>> QAll;
	std::vector<std::vector<float>> TSAll;
	std::vector<std::vector<float>> FlagList;
    	
    	float QL=0.0;
    	float QH=0.0;
	//TH2F* hTDCfireChanvsEvent;
    double yOffset_{35.};
    double yToIDfactor_{50./80.};

  };
}

#endif /* TRIGSCINT_CHARGEANALYZERCRAZY_H */
