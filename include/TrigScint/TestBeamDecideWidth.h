/**
 * @file TestBeamDecideWidth.h
 * @brief
 * @author
 */

#ifndef TRIGSCINT_TESTBEAMDECIDEWIDTH_H
#define TRIGSCINT_TESTBEAMDECIDEWIDTH_H

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
  class TestBeamDecideWidth : public framework::Analyzer {
  public:

	TestBeamDecideWidth(const std::string& name, framework::Process& process); // : framework::Analyzer(name, process) {}
	virtual ~TestBeamDecideWidth();
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
	float PE;

	//plotting stuff 
	int evNb;
    	int nEv{200};
    	int nChannels{16};
    	int nTrkMax{100};
    	int fillNb{0};
    	int widthlimit = 10;
    	int maximumpulsewidth = 7;
	float totChargelimit = 5000.0;    	
    	std::vector <float> QChan0;
	std::vector <float> QChan1;
	std::vector <float> QChan2;
	std::vector <float> QChan3;
	std::vector <float> QChan4;
	std::vector <float> QChan5;
	std::vector <float> QChan6;
	std::vector <float> QChan7;
	std::vector <float> QChan8;
	std::vector <float> QChan9;
	std::vector <float> QChan10;
	std::vector <float> QChan11;
    	
    	std::vector<std::vector<float>> TotQChan0;
	std::vector<std::vector<float>> TotQChan1;
	std::vector<std::vector<float>> TotQChan2;
	std::vector<std::vector<float>> TotQChan3;
	std::vector<std::vector<float>> TotQChan4;
	std::vector<std::vector<float>> TotQChan5;
	std::vector<std::vector<float>> TotQChan6;
	std::vector<std::vector<float>> TotQChan7;
	std::vector<std::vector<float>> TotQChan8;
	std::vector<std::vector<float>> TotQChan9;
	std::vector<std::vector<float>> TotQChan10;
	std::vector<std::vector<float>> TotQChan11;

	std::vector<int> PulseWidth;
	
	TH1F* ChannelWidthData[12]; 

	//TH2F* hTDCfireChanvsEvent;
    double yOffset_{35.};
    double yToIDfactor_{50./80.};

  };
}

#endif /* TRIGSCINT_TESTBEAMDECIDEWIDTH_H */
