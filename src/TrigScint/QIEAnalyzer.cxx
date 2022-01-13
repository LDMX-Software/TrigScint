/**
 * @file QIEAnalyzer.cxx
 * @brief An analyzer drawing the most basic quanities of Trigger Scintillator bars
 * @author Lene Kristian Bryngemark, Stanford University
 */

#include "TrigScint/QIEAnalyzer.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"

#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"


namespace trigscint {

  /*                                                    
   * Integrated charge given by configurable pulse shape
   * @par ts time sample to be evaluated                
   * @par par array of parameters of the pulse.         
   * par[0] = Total integral of the pulse (Q0)          
   * par[1] = start time of the pulse (t0)              
   * par[2] = 1/RC time constant (k)                    
   * par[3] = time when pulse attains maximum (tmax)    
   * par[4] = pedestal                                  
   */
  Double_t SinglePulseShape(Double_t* ts, Double_t* par);
  Double_t SinglePulseShapeNoPed(Double_t* ts_, Double_t* par);

  // Temporarily store charge per time sample
  float charge[15];
  /// array of integers 0,1,2,..,maxTS-1
  float time[15];
  // to manually calculate the pedestal
  Float_t PEDESTAL{0};
  
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
    std::vector<float> Digi_weight[nChannels];
    std::vector<float> QErr[nChannels];
    float AvgQs[nChannels];
    
    for (auto chan : channels) {
      std::vector<float> q = chan.getQ();
      int nTimeSamp = q.size();
      int bar = chan.getChanID();
      AllQs[bar] = chan.getQ();
      Digi_weight[bar] = chan.getQ();
      QErr[bar] = chan.getQ();
      AvgQs[bar] = chan.getAvgQ();
      float qTot = 0;
      int firstT = -1;
      int nQ400{0};
      for (int iT = 0; iT < q.size() ; iT++) {
	Double_t temp = smq->QErr(AllQs[bar][iT])*AllQs[bar][iT];
	if(420<AllQs[bar][iT] && AllQs[bar][iT]<430)
	  nQ400++;
	QErr[bar][iT] = temp;
	if(temp>0) Digi_weight[bar][iT] = 1/temp;
	else Digi_weight[bar][iT] = 0;
	ldmx_log(debug) << "in event " << evNb << "; channel " << bar << ", got charge[" << iT << "] = " << q.at(iT);
	// if ( evNb < nEv && bar < nChannels ) //stick within the predefined histogram array
	//   hOut[evNb][bar]->Fill(iT, q.at(iT));
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
      float QAvg18=0;
      float QAvg1530=0;
      
      for(int ts=0;ts<30;ts++){
	AllPulses->Fill(ts,AllQs[bar][ts]);
	hAllPulses[bar]->Fill(ts,AllQs[bar][ts],Digi_weight[bar][ts]);
      }
      for(int ts=0;ts<15;ts++){
	AllNoise->Fill(AllQs[bar][ts]);
	
	// Double_t err = smq->QErr(AllQs[bar][ts])*AllQs[bar][ts];
	// if(err>0)
	//   Q15_weighted->Fill(AllQs[bar][ts],1.0/err);
	Q15_weighted->Fill(AllQs[bar][ts],Digi_weight[bar][ts]);
	QAvg15+=AllQs[bar][ts];
	QAvg18+=AllQs[bar][ts];
	QAvg1530+=AllQs[bar][ts+15];
      }

      for(int ts=15;ts<18;ts++)
	QAvg18+=AllQs[bar][ts];

      hQAvg15[bar]->Fill(QAvg15/15.0);
      hQAvg1530[bar]->Fill(QAvg1530/15.0);
      hQ18Avg[1-bar%2]->Fill(QAvg18/18.0);

      if(420<(QAvg15/15.0) && (QAvg15/15.0)<430) {
	hQ400Freq[1-bar%2]->Fill(nQ400);
	hQ400QAvg[1-bar%2]->Fill(AvgQs[bar]);
	
	if(Q400Count<100){
	  for(int ts=0;ts<nTimeSamp;ts++)
	    hQ400Pulses[Q400Count]->Fill(ts,AllQs[bar][ts]);
	  Q400Count++;
	}
	
	for(int ts=0;ts<nTimeSamp;ts++)
	  hQ400Pulse_Distribution[1-bar%2]->Fill(ts,AllQs[bar][ts],Digi_weight[bar][ts]);
      }
      
      // Single photon event
      bool GoodPhot1Event=false;
      if(25<QAvg15/15.0 && QAvg15/15.0<45){
	// Define a good photon event
	for(int ts=0;ts<15;ts++)
	  GoodPhot1Event |= (AllQs[bar][ts]>190);
	
	if(GoodPhot1Event){
	  if(photcount == 3000) return;
	  
	  for(int ts=0;ts<15;ts++)
	    hQ15Phot1->Fill(AllQs[bar][ts]);
	  
	    // Fit pulse to single-PE events
	    for(int i=0;i<15;i++)
	      charge[i] = AllQs[bar][i];
	    
	    // Find maximum of pulse
	    int t1 = 0;		// First time sample to consider
	    int t2 = 5;		// Last tim sample to consider
	    float QMax15=AllQs[bar][0]; // Pulse maximum
	    
	    for(int i=1;i<15;i++) {
	      if(AllQs[bar][i]>QMax15) {
	    	QMax15 = AllQs[bar][i];
	    	t1 = i-1;
	    	t2 = i+4;
	      }
	    }
	    if(t1>10) continue;
	    // Find pedestal by averaging charge collected -1,+4 time samples
	    // away from the pulse maximum
	    Float_t nts=0;
	    PEDESTAL=0;
	    for(int i=0;i<t1;i++){
	      PEDESTAL+=AllQs[bar][i];
	      nts++;
	    }
	    for(int i=14;i>t2;i--){
	      PEDESTAL+=AllQs[bar][i];
	      nts++;
	    }
	    PEDESTAL = PEDESTAL/nts;
	  
	    // printf("photcount = %d\n",photcount);
	    
	    TF1* func = new TF1("myfunc",SinglePulseShape,0,380,5);
	    TF1* func2 = new TF1("Pedestal-assumed",SinglePulseShapeNoPed,0,380,4); 	// Fit 2
	    TGraphErrors* tg = new TGraphErrors(maxTS,time,&AllQs[bar][0],0,&QErr[bar][0]);

	    GuessPar(func, &AllQs[bar][0],t1+1,true);
	    GuessPar(func2, &AllQs[bar][0],t1+1,false);
	    
	    // Print event info for debugging
	    ldmx_log(debug)<<"\n*** This prints before fitting";
	    ldmx_log(debug)<<"Good PE event no.: "<<photcount;

	    ldmx_log(debug)<<"\nid\t";
	    for(int i=0;i<15;i++)
	      ldmx_log(debug)<<i<<"\t";

	    ldmx_log(debug)<<"\nQ[]\t";
	    for(int i=0;i<15;i++)
	      ldmx_log(debug)<<AllQs[bar][i]<<"\t";

	    ldmx_log(debug)<<"\nt1=  "<<t1<<"\tt2= "<<t2<<"\n";
	    
	    TCanvas* tc = new TCanvas(Form("c_evt_%i",photcount),"bb",800,600);
	    tg->Draw("APE");
	    tg->Fit(func,"NQ");
	    tg->Fit(func2,"NQ");
	    
	    hPhot1Q0->Fill(func->GetParameter(0));
	    hPhot1T0->Fill(func->GetParameter(1));
	    hPhot1K->Fill(func->GetParameter(2));
	    hPhot1TMax->Fill(func->GetParameter(3));
	    hPhot1Ped->Fill(func->GetParameter(4));
	    hPhot1Chi->Fill(func->GetChisquare());

	    hPhot1Q0NoPed->Fill(func2->GetParameter(0));
	    hPhot1T0NoPed->Fill(func2->GetParameter(1));
	    hPhot1KNoPed->Fill(func2->GetParameter(2));
	    hPhot1TMaxNoPed->Fill(func2->GetParameter(3));
	    hPhot1PedNoPed->Fill(PEDESTAL);
	    hPhot1ChiNoPed->Fill(func2->GetChisquare());
	    
	    delete func;
	    delete func2;
	    delete tg;
	    delete tc;
	}
	photcount++;
      }	// End of good photon loop
    
      Abs_Dev->Fill(chan.getMedQ(),chan.getMaxQ()-chan.getMinQ());
    }
  
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
  // }//over channels
  
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

    smq = new SimQIE();
    getHistoDirectory();

    AllNoise = new TH1F("AllNoise","Charge distribution in 1st 15 time samples;Charge [fC];",620,-20,600);
    AllPulses = new TH2F("AllPulses","All pulses seen;Time sample; Charge[fC]",30,0,30,100,-20,5000);
    Q15_weighted = new TH1F("Q15_weighted","Charge distribution in 1st 15 time samples;Charge [fC];",620,-20,600);
    PulseShape = new TH2F("hPulse","Events with Qavg>100;time sample;Q [fC]",30,0,30,100,-20,1000);
    hGoodPulses = new TH2F("hGoodPulses","Events with Qmed<300;time sample;Q [fC]",30,0,30,100,-20,5000);
    Abs_Dev = new TH2F("hAbs_Dev",";Qmed[fC];Qmax-Qmin [fC]",100,-500,500,100,0,5000);
    
    hQ18Avg[0] = new TH1F("hQ18Avg_LYSO","Average charge in 1st 18ts (LYSO);Charge [fC];",620,-20,600);
    hQ18Avg[1] = new TH1F("hQ18Avg_Plast","Average charge in 1st 18ts (Plastic);Charge [fC];",620,-20,600);

    hQ400Freq[0] = new TH1F("hQ400Freq_LYSO","Occurance of charge 424fC in LYSO;No. of occurances;",31,0,31);
    hQ400Freq[1] = new TH1F("hQ400Freq_Plast","Occurance of charge 424fC in Plastic;No. of occurances;",31,0,31);
    hQ400Pulse_Distribution[0] = new TH2F("hQ400Pulse_Distribution_LYSO","Overall pulse shape for 420<Q15Avg<430 [LYSO];time sample;Charge [fC]",31,0,31,102,-20,1000);
    hQ400Pulse_Distribution[1] = new TH2F("hQ400Pulse_Distribution_Plast","Overall pulse shape for 420<Q15Avg<430 [Plastic];time sample;Charge [fC]",31,0,31,102,-20,1000);
    hQ400QAvg[0] = new TH1F("hQ400QAvg_LYSO","Average Charge for events with 420<Q15Avg<430 [LYSO];Average Charge [fC];",1000,300,500);
    hQ400QAvg[1] = new TH1F("hQ400QAvg_Plast","Average Charge for events with 420<Q15Avg<430 [Plastic];Average Charge [fC];",1000,300,500);

    // Trial 3
    hPhot1Trial3K[0] = new TH1F("hPhot1Trial3K_0","",100,0,0.3);
    hPhot1Trial3K[1] = new TH1F("hPhot1Trial3K_1","",100,0,0.3);
    hPhot1Trial3K[2] = new TH1F("hPhot1Trial3K_2","",100,0,0.3);
    hPhot1Trial3KAvg = new TH1F("hPhot1Trial3KAvg","",100,0,0.3);

    
    // Single-PE fitting related histograms
    hPhot1Q0 = new TH1F("hPhot1Q0","",100,-1000,1000);
    hPhot1T0 = new TH1F("hPhot1T0","",100,-1000,1000);
    hPhot1K = new TH1F("hPhot1K","",100,-20,30);
    hPhot1TMax = new TH1F("hPhot1TMax","",100,-1000,1000);
    hPhot1Ped = new TH1F("hPhot1Ped","",500,-500,500);
    hPhot1Chi = new TH1F("hPhot1Chi","",100,0,1e5);

    // Single-PE fitting related histograms
    hPhot1Q0NoPed = new TH1F("hPhot1Q0NoPed","",100,-1000,1000);
    hPhot1T0NoPed = new TH1F("hPhot1T0NoPed","",100,-1000,1000);
    hPhot1KNoPed = new TH1F("hPhot1KNoPed","",100,-20,30);
    hPhot1TMaxNoPed = new TH1F("hPhot1TMaxNoPed","",100,-1000,1000);
    hPhot1PedNoPed = new TH1F("hPhot1PedNoPed","",500,-500,500);
    hPhot1ChiNoPed = new TH1F("hPhot1ChiNoPed","",100,0,1e5);
    
    for(int i=0;i<nChannels;i++){
      hQAvg15[i] = new TH1F(Form("hQAvg15_%i",i),Form("Average Charge distribution in 1st 15 time samples [Channel %i];Charge [fC];",i),620,-20,600);
      hQAvg1530[i] = new TH1F(Form("hQAvg1530_%i",i),Form("Average Charge distribution in last 15 time samples [Channel %i];Charge [fC];",i),620,-20,600);
    }
    
    time[0]=0.;
    zeroes = new float[15];
    for (int ii=1;ii<15;ii++){
      time[ii] = time[ii-1]+25.0;
      zeroes[ii] = 0.0;
    }
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
    // for (int iE = 0; iE<nEv; iE++) {
    //   for (int iB = 0; iB<nChannels; iB++) {
    // 	//TH1F * hOut = new TH1F(Form("hCharge_chan%i_ev%i", iB, iE), Form(";time sample; Q, channel %i, event %i [fC]", iB, iE), nTimeSamp, -0.5, nTimeSamp-0.5);
    // 	hOut[iE][iB] = new TH1F(Form("hCharge_chan%i_ev%i", iB, iE), Form(";time sample; Q, channel %i, event %i [fC]", iB, iE), nTimeSamp,-0.5,nTimeSamp-0.5);
    //   }
    // }
    int PEmax=500;
    int Qmax=3.5e5;
    for (int iB = 0; iB<nChannels; iB++) {
      hPE[iB]=new TH1F(Form("hPE_chan%i", iB), Form(";PE, chan%i", iB),5*PEmax,0,PEmax);
      hPEvsT[iB]=new TH2F(Form("hPEvsT_chan%i", iB), Form(";First time sample above summing threshold;PE, chan%i", iB),nTimeSamp+1,-1.5,nTimeSamp-0.5, 5*PEmax,0,PEmax);
    }
    
    hQ15Phot1 = new TH1F("hQ15Phot1","Single photon events charge statistics;Q [fC];",43,-16,500);

    for (int i=0;i<100;i++)
      hPhot1Pulse[i] = new TH1F(Form("hPhot1Pulse_%i",i),"Event with 25< QAvg15 <45;time sample;Q [fC]",15,0,15);    
    
    for (int i=0;i<100;i++)
      hQ400Pulses[i] = new TH1F(Form("hQ400Pulses_%i",i),Form("Event %i with 420< QAvg15 <430;time sample;Q [fC]",i),31,0,31);
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

    for(int i=0;i<nChannels;i++)
      hAllPulses[i] = new TH2F(Form("hAllPulses_%i",i),Form("Channel %i;Time sample; Charge[fC]",i),30,0,30,100,-20,4000);
    return;
  }
  

  void QIEAnalyzer::onProcessEnd() {
    return;
  }


  Double_t SinglePulseShape(Double_t* ts, Double_t* par) {
    auto pls = new Expo(par[2],par[3]);
    pls->AddPulse(par[1],par[0]);
    Double_t integral = par[4]+pls->Integrate(ts[0],ts[0]+25);
    delete pls;
    return integral;
  }
  
  Double_t SinglePulseShapeNoPed(Double_t* ts_, Double_t* par) {
    auto pls = new Expo(par[2],par[3]); // Initiate the pulse object with k and tmax
    pls->AddPulse(par[1],par[0]);	      // add a pulse with t0 and integral Q0

    Double_t integral = PEDESTAL+pls->Integrate(ts_[0],ts_[0]+25);
    delete pls;
    return integral;
  }
  
  void QIEAnalyzer::GuessPar(TF1* ff, float* charge,int max_id,bool CalculatePed)
  {
    float T=25.0;                   // time period
    float k_ = log((charge[max_id+1]-charge[max_id+2])/(charge[max_id+2]-charge[max_id+3]))/T;

    float const1 = 1/(1-exp(-k_*T));

    float ped_{0};
    if(CalculatePed)
      ped_ = (charge[max_id+2]-charge[max_id+1]*exp(-k_*T))*const1;
    else
      ped_ = PEDESTAL;

    float Q0_ = charge[max_id-1]+charge[max_id]-2*ped_+(charge[max_id+1]-charge[max_id+2])*const1*const1;
    float t0_ = time[max_id-1];
    float tmax_ = 60.0;

    if(CalculatePed){
      ff->SetParameters(Q0_,t0_,k_,tmax_,ped_);
      ff->SetParNames("Q0","T0","K","TMax","Ped");
    }
    else{
      ff->SetParameters(Q0_,t0_,k_,tmax_);
      ff->SetParNames("Q0","T0","K","TMax");
    }

    ldmx_log(debug)<<"The guessed values:\n";
    ldmx_log(debug)<<"Q0\tt0\tK\ttmax\tpedestal\n";
    ldmx_log(debug)<<Q0_<<"\t"<<t0_<<"\t"<<k_<<"\t"<<tmax_<<"\t"<<ped_<<"\n";
    return;
  }
  
}

DECLARE_ANALYZER_NS(trigscint, QIEAnalyzer)
