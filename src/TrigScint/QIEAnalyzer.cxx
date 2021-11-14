/**
 * @file QIEAnalyzer.cxx
 * @brief An analyzer drawing the most basic quanities of Trigger Scintillator bars
 * @author Lene Kristian Bryngemark, Stanford University
 */

#include "TrigScint/QIEAnalyzer.h"

namespace trigscint {

  QIEAnalyzer::QIEAnalyzer(const std::string& name,
			   framework::Process& process)
    : Analyzer(name, process) {}
  QIEAnalyzer::~QIEAnalyzer() {}
  
  void QIEAnalyzer::configure(framework::config::Parameters &parameters){

    inputCol_  = parameters.getParameter< std::string >("inputCollection");
    inputPassName_  = parameters.getParameter< std::string >("inputPassName");
    peds_  = parameters.getParameter< std::vector<double> >("pedestals");

    std::cout << " [ QIEAnalyzer ] In configure(), got parameters " 
	      << "\n\t inputCollection = " << inputCol_
	      << "\n\t inputPassName = " << inputPassName_
	      << "\n\t pedestals[0] = " << peds_[0]
	      << "\t." << std::endl;

    return;
  }

  void QIEAnalyzer::analyze(const framework::Event &event) {

    const auto channels{event.getCollection<trigscint::EventReadout >(inputCol_, inputPassName_)};

    int evNb = event.getEventNumber();
    int nChan = channels.size();

    std::vector<float> AllQs[nChannels];
    float AvgQs[nChannels];
    for (auto chan : channels) {
      std::vector<float> q = chan.getQ();
      int nTimeSamp = q.size();
      int bar = chan.getChanID();
      AllQs[bar] = chan.getQ();
      AvgQs[bar] = chan.getAvgQ();
      float qTot = 0;
      int firstT = -1;

      if(chan.getAvgQ()>100)
	for(int ts=0;ts<AllQs[bar].size();ts++)
	  PulseShape->Fill(ts,AllQs[bar][ts]);

      float QAvg15=0;
      for(int ts=0;ts<15;ts++){
	AllNoise->Fill(AllQs[bar][ts]);
	QAvg15+=AllQs[bar][ts];
      }
      hQAvg15->Fill(QAvg15/15.0);
      
      // Single photon event
      bool GoodPhot1Event=false;
      if(25<QAvg15/15.0 && QAvg15/15.0<45){
	for(int ts=0;ts<15;ts++)
	  GoodPhot1Event |= (AllQs[bar][ts]>190);
	if(GoodPhot1Event)
	  for(int ts=0;ts<15;ts++)
	    hQ15Phot1->Fill(AllQs[bar][ts]);

	if(photcount<100){
	  for(int ts=0;ts<15;ts++)
	    hPhot1Pulse[photcount]->Fill(ts,AllQs[bar][ts]);
	  photcount++;
	}
      }
      
      float Tailfall = AllQs[bar][25]+AllQs[bar][26];
      Tailfall += AllQs[bar][27];
      Tailfall += AllQs[bar][28];
      Tailfall += AllQs[bar][29];
      
      // bool Tailfall = AllQs[bar][25]>AllQs[bar][26];
      // Tailfall &= AllQs[bar][26]>AllQs[bar][27];
      // Tailfall &= AllQs[bar][27]>AllQs[bar][28];
      // Tailfall &= AllQs[bar][28]>AllQs[bar][29];
      
      // if((chan.getMedQ()<300)&&(Tailfall<1400))
      if(Tailfall<1400)
	for(int ts=0;ts<AllQs[bar].size();ts++)
	  hGoodPulses->Fill(ts,AllQs[bar][ts]);
      Abs_Dev->Fill(chan.getMedQ(),chan.getMaxQ()-chan.getMinQ());
    }//over channels

    // Charge Correlation
    for(int i=1;i<nChannels;i++){
      for(int j=0;j<i;j++){
	for(int ts=0;(ts<AllQs[i].size())&&(ts<AllQs[j].size());ts++){
	  hQiQj[(i-1)*i/2+j]->Fill(AllQs[i][ts],AllQs[j][ts]);
	  hQiQj_ts[(i-1)*i/2+j][ts]->Fill(AllQs[i][ts],AllQs[j][ts]);
	}
	hAvgQiQj[(i-1)*i/2+j]->Fill(AvgQs[i],AvgQs[j]);
      }
    }

    return;
  }

  // /*
  void QIEAnalyzer::onFileOpen() {
    
    return;
  }

  void QIEAnalyzer::onFileClose() {

    return;
  }
  //  */
  
  void QIEAnalyzer::onProcessStart() {

    getHistoDirectory();

    int nTimeSamp=40;
    int PEmax=500;
    int Qmax=3.5e5;

    // hQ15Phot1 = new TH1F("hQ15Phot1","Single photon events charge statistics;Q [fC];",520,-20,500);
    hQ15Phot1 = new TH1F("hQ15Phot1","Single photon events charge statistics;Q [fC];",43,-16,500);

    for (int i=0;i<100;i++)
      hPhot1Pulse[i] = new TH1F(Form("hPhot1Pulse_%i",i),"Event with 25< QAvg15 <45;time sample;Q [fC]",15,0,15);
    
    AllNoise = new TH1F("AllNoise","Charge distribution in 1st 15 time samples;Charge [fC];",620,-20,600);
    hQAvg15 = new TH1F("hQAvg15","Average Charge distribution in 1st 15 time samples;Charge [fC];",620,-20,600);
    PulseShape = new TH2F("hPulse","Events with Qavg>100;time sample;Q [fC]",30,0,30,100,-20,1000);
    hGoodPulses = new TH2F("hGoodPulses","Events with Qmed<300;time sample;Q [fC]",30,0,30,100,-20,5000);
    Abs_Dev = new TH2F("hAbs_Dev",";Qmed[fC];Qmax-Qmin [fC]",100,-500,500,100,0,5000);
    
    int counter=0;
    for(int c1=1;c1<nChannels;c1++){
      for(int c2=0;c2<c1;c2++){
    	hQiQj[counter] = new TH2F(Form("hQ%iQ%i",c1,c2),Form(";Charge [fC], channel %i;Charge [fC], channel %i;",c1,c2),200,-100,1000,200,-100,1000);
	hAvgQiQj[counter] = new TH2F(Form("hAvgQ%iQ%i",c1,c2),Form(";Average Charge [fC], channel %i;Average Charge [fC], channel %i;",c1,c2),300,-100,500,300,-100,500);
	for(int ts=0;ts<maxTS;ts++)
	  hQiQj_ts[counter][ts] = new TH2F(Form("hQ%iQ%i_ts%i",c1,c2,ts),Form("Time sample %i;Charge [fC], channel %i;Charge [fC], channel %i;",ts,c1,c2),200,-100,1000,200,-100,1000);
	
    	counter++;
    	if(counter>200) return;
      }
    }
    evNb = 0;
	
    return;
  }
  

  void QIEAnalyzer::onProcessEnd() {
    return;
  }


}

DECLARE_ANALYZER_NS(trigscint, QIEAnalyzer)
