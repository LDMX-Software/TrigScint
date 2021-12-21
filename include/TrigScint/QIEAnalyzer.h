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
#include "TrigScint/QIEInputPulse.h"
#include "TrigScint/SimQIE.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "QIEInputPulse.h"

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

    // /*
    //  * Integrated charge given by configurable pulse shape
    //  * @par ts time sample to be evaluated
    //  * @par par array of parameters of the pulse.
    //  * par[0] = Total integral of the pulse (Q0)
    //  * par[1] = start time of the pulse (t0)
    //  * par[2] = 1/RC time constant (k)
    //  * par[3] = time when pulse attains maximum (tmax)
    //  * par[4] = pedestal
    //  */
    // Double_t SinglePulseShape(Double_t ts, Double_t* par);

    // /// array of integers 0,1,2,..,maxTS-1
    // float* time;
    /// array of zeroes
    float* zeroes;

    // const
    static const int maxTS{30};	// Maximum time samples recorded per channel
    static const int nEv{200};
    static const int nChannels{16}; // No. of channels available
    
    // // Temporarily store charge per time sample
    // float charge[maxTS];

  private:

    std::vector <std::vector <TH1F*> > vChargeVsTime;
    SimQIE* smq;
	
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
    int nTrkMax{100};

    //match nev, nchan above
    TH1F* hOut[nEv][nChannels];;
    TH1F* hPE[nChannels];
    TH2F* hPEvsT[nChannels];
    TH2F* hAllPulses[nChannels];

    TH2F* hAvgQiQj[nChannels*(nChannels-1)/2]; // Integrated charge correlation plots
    TH2F* hQiQj[nChannels*(nChannels-1)/2]; // charge correlation plots
    TH2F* hQiQj_ts[nChannels*(nChannels-1)/2][maxTS]; // charge correlation plots per time sample
    TH2F* PulseShape;				      // Pulse shape in the high charge event
    TH1F* Rel_Dev;				      // Relative deviation in charge deposition
    TH2F* Abs_Dev;				      // Absoulte deviation w.r.t. median
    TH2F* hGoodPulses;	// Pulse shapes with Qmed<300
    
    TH1F* AllNoise;	// Charge in 1st 15 ts
    TH1F* Q15_weighted;	// Charge in 1st 15 ts, weighted by sensitivity
    TH1F* hQAvg15[nChannels];	// Average Q in 1st 15 ts
    TH1F* hQAvg1530[nChannels];	// Average Q in last 15 ts
    
    TH1F* hPhot1Pulse[100];	// per event pulse shape
    int photcount{0};		// No. of single photon events passed
    TH1F* hQ15Phot1;		// charge distribution in single photon event
    TH2F* AllPulses;

    // Single-PE fitting
    TH1F* hPhot1Q0;
    TH1F* hPhot1T0;
    TH1F* hPhot1K;
    TH1F* hPhot1TMax;
    TH1F* hPhot1Ped;
    TH1F* hPhot1Chi;

    TH1F* hPhot1Trial3K[3];	// k's calculated from trial 3
    TH1F* hPhot1Trial3KAvg;	// Average k's calculated from trial 3

    TH1F* hQ18Avg[2];		// 2 hists for LYSO and Plastic respectively.

    int Q400Count{0};		// No. of single photon events passed
    TH1F* hQ400Freq[2];		// Frequency of getting 400fC charge per time sample
    TH1F* hQ400Pulses[100];	// Pulses with avg. charge around 400fC.
    TH2F* hQ400Pulse_Distribution[2]; // Collection of all the pulses. LYSO and Plastic
    TH1F* hQ400QAvg[2];		      // Average charge distribution for 420<Q15Avg<430

    double yOffset_{35.};
    double yToIDfactor_{50./80.};

    void FitSinglePulse(Float_t* Qi);	// Fit a single pulse using 6 time samples
  };
}

#endif /* TRIGSCINT_QIEANALYZER_H */
