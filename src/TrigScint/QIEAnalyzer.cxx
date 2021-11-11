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


      // for(int chan = 0;chan<nChannels;chan++){
      // 	for(int ts=0;ts<maxTS;ts++){
      // 	  Qi[][].push_back();
      // 	}
      // }
      
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

    int counter=0;
    for(int c1=1;c1<nChannels;c1++){
      for(int c2=0;c2<c1;c2++){
    	// hQiQj[counter] = new TH2F(Form("hQ%iQ%i",c1,c2),"",100,0,Qmax/10,100,0,Qmax/10);
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
