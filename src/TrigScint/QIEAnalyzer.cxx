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
      for (int iT = 0; iT < q.size() ; iT++) {
	ldmx_log(debug) << "in event " << evNb << "; channel " << bar << ", got charge[" << iT << "] = " << q.at(iT);
	if ( evNb < nEv && bar < nChannels ) //stick within the predefined histogram array
	  hOut[evNb][bar]->Fill(iT, q.at(iT));
	if ( q.at(iT) > 2*fabs(peds_[ bar ]) ) 	{ //integrate all charge well above ped to convert to a PE count
	  qTot+=q.at(iT);
	  if (firstT =-1) //keep track of first time sample above threshold
	    firstT=iT;
	}//if above threshold
      }//over time samples
      float PE = qTot*6250./4.e6;
      hPE[ bar ]->Fill( PE );
      hPEvsT[ bar ]->Fill( firstT, PE );
	  
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

    /*
      int yMax = 50;
      int yMin = -yMax;
      int nBinsY = (yMax-yMin)/1; //1 mm resolution
      //    
      hSimYEvnb = new TH2F("hSimYEvnb","Beam electron y for each event;Event number;y",nEv,0,nEv, nBinsY,yMin,yMax);
      hSimIDEvnb = new TH2F("hSimIDEvnb","Beam electron y converted to channel ID, for each event;Event number;ID",nEv,0,nEv, nBinsY,convertToID(yMin),convertToID(yMax));
      hIdEvnb = new TH2F("hIdEvnb","Id track vs event, filled with N_{PE};Event;Channel ID;N_{PE}", nEv,0,nEv, nChannels+2,0,nChannels+2);
      hId=new TH1F("hId","channel id;Channel ID", nChannels,0,nChannels);
      hNtracksEvnb=new TH1F("hNtracksEvnb","N_tracks for each event;Event;N_{tracks}", nEv,0,nEv);
      hNtracks=new TH1F("hNtracks","N_tracks per event;Event;N_{tracks}", nTrkMax,0,nTrkMax);
      hFindableTracks=new TH1F("hFindableTracks","N_findableTracks per event;Event;N_{tracks}^{findable}", nTrkMax,0,nTrkMax);
      hTrackMatrix=new TH2F("hTrackMatrix","Found vs findable tracks per event;N_{tracks}^{findable};N_{tracks}", nTrkMax,-0.5,nTrkMax-0.5, nTrkMax,-0.5,nTrkMax-0.5);
      hNhits=new TH1F("hNhits","Nhits in a track;Track width [channel Nb];N_{tracks}", 60, 0, 6);
    */
    int nTimeSamp=40;
    for (int iE = 0; iE<nEv; iE++) {
      for (int iB = 0; iB<nChannels; iB++) { 
	//TH1F * hOut = new TH1F(Form("hCharge_chan%i_ev%i", iB, iE), Form(";time sample; Q, channel %i, event %i [fC]", iB, iE), nTimeSamp, -0.5, nTimeSamp-0.5);
	hOut[iE][iB] = new TH1F(Form("hCharge_chan%i_ev%i", iB, iE), Form(";time sample; Q, channel %i, event %i [fC]", iB, iE), nTimeSamp,-0.5,nTimeSamp-0.5);
      }
    }
    int PEmax=500;
    int Qmax=3.5e5;
    for (int iB = 0; iB<nChannels; iB++) {
      hPE[iB]=new TH1F(Form("hPE_chan%i", iB), Form(";PE, chan%i", iB),5*PEmax,0,PEmax);
      hPEvsT[iB]=new TH2F(Form("hPEvsT_chan%i", iB), Form(";First time sample above summing threshold;PE, chan%i", iB),nTimeSamp+1,-1.5,nTimeSamp-0.5, 5*PEmax,0,PEmax);
    }
    
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
