/**
 * @file ChargeAnalyzer.h
 * @brief
 * @author
 */

#ifndef TRIGSCINT_CHARGEANALYZER_H
#define TRIGSCINT_CHARGEANALYZER_H

//LDMX Framework                               
#include "Framework/Configure/Parameters.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "TrigScint/Event/EventReadout.h"
#include "TrigScint/Event/TrigScintQIEDigis.h"
#include "TH1.h"
#include "TH2.h"

namespace trigscint {

  /**
   * @class QIEAnalyzer
   * @brief
   */
  class ChargeAnalyzer : public framework::Analyzer {
  public:

	ChargeAnalyzer(const std::string& name, framework::Process& process); // : framework::Analyzer(name, process) {}
	virtual ~ChargeAnalyzer();
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

	//plotting stuff 
	int evNb;
    int nEv{200};
    int nChannels{16};
    int nTrkMax{100};
    	int fillNb{0};
    	
    	int Qlow_thr = 2.e2;
    	int Qmed_thr = 1.e3;
    	int Qhigh_thr = 10.e3;
    	
	TH1F* hQ_low[16];
	TH1F* hQ_med[16];
	TH1F* hQ_high[16];
	TH2F* hQvsTS_low[16];
	TH2F* hQvsTS_med[16];
	TH2F* hQvsTS_high[16];
	TH1F* hQTot_low[16];
	TH1F* hQTot_med[16];
	TH1F* hQTot_high[16];
	TH1F* hQAvg_low[16];
	TH1F* hQAvg_med[16];
	TH1F* hQAvg_high[16];
	TH1F* hQMin_low[16];
	TH1F* hQMin_med[16];
	TH1F* hQMin_high[16];
	TH1F* hQMax_low[16];
	TH1F* hQMax_med[16];
	TH1F* hQMax_high[16];
	TH1F* hQMed_low[16];
	TH1F* hQMed_med[16];
	TH1F* hQMed_high[16];
	TH2F* hQTotvschan_low;
	TH2F* hQTotvschan_med;
	TH2F* hQTotvschan_high;
	TH2F* hQAvgvschan_low;
	TH2F* hQAvgvschan_med;
	TH2F* hQAvgvschan_high;
	TH2F* hQMinvschan_low;
	TH2F* hQMinvschan_med;
	TH2F* hQMinvschan_high;
	TH2F* hQMaxvschan_low;
	TH2F* hQMaxvschan_med;
	TH2F* hQMaxvschan_high;
	TH2F* hQMedvschan_low;
	TH2F* hQMedvschan_med;
	TH2F* hQMedvschan_high;
	TH1F* hQTot;
	
	//TH2F* hTDCfireChanvsEvent;
    double yOffset_{35.};
    double yToIDfactor_{50./80.};

  };
}

#endif /* TRIGSCINT_QIEANALYZER_H */