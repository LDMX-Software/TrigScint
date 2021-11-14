/**
 * @file QIEAnalyzer.h
 * @brief
 * @author
 */

#ifndef TRIGSCINT_QIEANALYZER_H
#define TRIGSCINT_QIEANALYZER_H

//LDMX Framework                               
#include "Framework/Configure/Parameters.h"
#include "Framework/EventProcessor.h" //Needed to declare processor
#include "TrigScint/Event/EventReadout.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"

namespace trigscint {

  /**
   * @class QIEAnalyzer
   * @brief
   */
  class QIEAnalyzer : public framework::Analyzer {
  public:

    QIEAnalyzer(const std::string& name, framework::Process& process); // : framework::Analyzer(name, process) {}
    virtual ~QIEAnalyzer();
    virtual void configure(framework::config::Parameters &parameters);

    virtual void analyze(const framework::Event &event) final override;

    //
    virtual void onFileOpen();

    //
    virtual void onFileClose();

    virtual void onProcessStart() final override;

    virtual void onProcessEnd() final override;

    float convertToID( float yVal ) { return (yVal+yOffset_)*yToIDfactor_; }

  private:

    std::vector <std::vector <TH1F*> > vChargeVsTime;
	
    // 
    // TH1* hId;
    // TH1* hNtracksEvnb;
    // TH1* hNtracks;
    // TH1* hFindableTracks;
    // TH2* hTrackMatrix;
    // TH1* hNhits;
    // TH2 *hSimYEvnb;
    // TH2 *hSimIDEvnb;
    // TH2 *   hIdEvnb;
	
    std::string inputCol_;
    std::string inputPassName_{""};
    std::vector<double> peds_;

    int evNb;
    // const
    static const int maxTS{30};	// Maximum time samples recorded per channel
    static const int nEv{1000};
    static const int nChannels{16}; // No. of channels available
    int nTrkMax{100};

    //match nev, nchan above
    TH1F* hOut[200][16];;
    TH1F* hPE[16];
    TH2F* hPEvsT[16];

    TH2F* hAvgQiQj[nChannels*(nChannels-1)/2]; // Integrated charge correlation plots
    TH2F* hQiQj[nChannels*(nChannels-1)/2]; // charge correlation plots
    TH2F* hQiQj_ts[nChannels*(nChannels-1)/2][maxTS]; // charge correlation plots per time sample
    TH2F* PulseShape;				      // Pulse shape in the high charge event
    TH1F* Rel_Dev;				      // Relative deviation in charge deposition
    TH2F* Abs_Dev;				      // Absoulte deviation w.r.t. median
    TH2F* hGoodPulses;	// Pulse shapes with Qmed<300
    TH1F* AllNoise;	// Charge in 1st 15 ts
    TH1F* hQAvg15;	// Average Q in 1dt 15 ts
    
    TH1F* hPhot1Pulse[100];	// per event pulse shape
    int photcount{0};		// No. of single photon events passed
    TH1F* hQ15Phot1;		// charge distribution in single photon event
    
    double yOffset_{35.};
    double yToIDfactor_{50./80.};

  };
}

#endif /* TRIGSCINT_QIEANALYZER_H */
