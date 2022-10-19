/**
 * @file HitAnalyzer.h
 * @brief
 * @author
 */

#ifndef TRIGSCINT_HITANALYZER_H
#define TRIGSCINT_HITANALYZER_H

//LDMX Framework                               
#include "Framework/Configure/Parameters.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "TrigScint/Event/TestBeamHit.h"
//#include "TrigScint/Event/TrigScintQIEDigis.h"
#include "TH1.h"
#include "TH2.h"

namespace trigscint {

  /**
   * @class HitAnalyzer
   * @brief
   */
  class HitAnalyzer : public framework::Analyzer {
  public:

	HitAnalyzer(const std::string& name, framework::Process& process); // : framework::Analyzer(name, process) {}
	virtual ~HitAnalyzer();
    virtual void configure(framework::config::Parameters &parameters);

    virtual void analyze(const framework::Event &event) final override;

	//
	virtual void onFileOpen();

	//
	virtual void onFileClose();

    virtual void onProcessStart() final override;

    virtual void onProcessEnd() final override;


  private:

	std::vector <std::vector <TH1F*> > vChargeVsTime;
	

	//configurable parameters
    	std::string inputCol_;
    	std::string inputPassName_{""};
	std::vector<double> peds_;
	int startSample_{0};

	//plotting stuff 
	int evNb;
    	int nEv{200};
    	int nChannels{16};
    	int S=0;
    	
    	int PE_low = 1;
    	int PE_med1 = 10;
    	int PE_med2 = 150;
    	
    	//int nTrkMax{100};
	//TH2F* hEvDisp;
	//TH2F* hEvDispPE;

	int fillNb{0};
	
	//match nev, nchan above
	TH1F* hBar;
	TH2F* hPEvsbar;
	TH1F* hPE_all;
	TH1F* hPE_low[16]; //Pedestal
	TH1F* hPE_med1[16]; //SinglePE
	TH1F* hPE_med2[16]; //MIP
	TH1F* hPE_high[16]; //Region X
	TH2F* hPEvsbar_low;
	TH2F* hPEvsbar_med1;
	TH2F* hPEvsbar_med2;
	TH2F* hPEvsbar_high;
  };
}

#endif /* TRIGSCINT_HITANALYZER_H */
