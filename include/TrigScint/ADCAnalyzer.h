/**
 * @file ADCAnalyzer.h
 * @brief
 * @author
 */

#ifndef TRIGSCINT_ADCANALYZER_H
#define TRIGSCINT_ADCANALYZER_H

//LDMX Framework                               
#include "Framework/Configure/Parameters.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "TrigScint/Event/EventReadout.h"
#include "TrigScint/Event/TrigScintQIEDigis.h"
#include "TH1.h"
#include "TH2.h"

namespace trigscint {

  /**
   * @class ADCAnalyzer
   * @brief
   */
  class ADCAnalyzer : public framework::Analyzer {
  public:

	ADCAnalyzer(const std::string& name, framework::Process& process); // : framework::Analyzer(name, process) {}
	virtual ~ADCAnalyzer();
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
    	int S=0;
    	int S1=0;

	//match nev, nchan above
	//TH1F* hOut[200][16];;
	//TH1F* hPE[16];
	//TH2F* hPEvsT[16];
	//TH2F* hPedSubtractedAvgQvsT[16];
	//TH2F* hPedSubtractedTotQvsPed[16];
	//TH2F* hPedSubtractedTotQvsN[16];
	//TH2F* hTotQvsPed[16];
	//TH2F* hPedSubtractedPEvsN[16];
	//TH2F* hPedSubtractedPEvsT[16];
	//TH2F* hAvgQvsT[16];
	TH1F* hADC[16];
	TH1F* hTDC[16];
	TH1F* hADC_Total;
	TH1F* hTDC_Total;
	TH2F* hADCvsTS[16];
	TH2F* hTDCvsTS[16];
	//TH2F* hTDCfireChanvsEvent;
    double yOffset_{35.};
    double yToIDfactor_{50./80.};

  };
}

#endif /* TRIGSCINT_ADCANALYZER_H */
