#include "TrigScint/ADCAnalyzer.h"

namespace trigscint {

  ADCAnalyzer::ADCAnalyzer(const std::string& name,
						   framework::Process& process)
    : Analyzer(name, process) {}
  ADCAnalyzer::~ADCAnalyzer() {}
  
  void ADCAnalyzer::configure(framework::config::Parameters &parameters){

    inputCol_  = parameters.getParameter< std::string >("inputCollection");
    inputPassName_  = parameters.getParameter< std::string >("inputPassName");
    peds_  = parameters.getParameter< std::vector<double> >("pedestals");
    gain_  = parameters.getParameter< std::vector<double> >("gain");
    startSample_  = parameters.getParameter< int >("startSample");

    std::cout << " [ ADCAnalyzer ] In configure(), got parameters " 
	      << "\n\t inputCollection = " << inputCol_
	      << "\n\t inputPassName = " << inputPassName_
	      << "\n\t startSample = " << startSample_
	      << "\n\t pedestals[0] = " << peds_[0]
	      << "\n\t gain[0] = " << gain_[0]
	      << "\t." << std::endl;

    return;
  }

  void ADCAnalyzer::analyze(const framework::Event &event) {

    //std::cout << "getting collection" << std::endl;
    const auto channels{event.getCollection<trigscint::TrigScintQIEDigis>(inputCol_, inputPassName_)};
    //std::cout << "got..." << std::endl;
  int evNb = event.getEventNumber();
  int nChannels = channels.size();
  //int S=0;
  
  for (auto chan : channels) {
    //int bar = chan.getBarID();
    int bar = chan.getChanID();
    std::vector<int> adc = chan.getADC();
    std::vector<int> tdc = chan.getTDC();
    
    //for(unsigned int i=0; i< adc.size(); i++) {
      //std::cout << "Element number " << "=" << i << "~" << adc[i] << "corresponding to channel number" << "==" << bar << std::endl;
    //}
    for(unsigned int i=0; i<adc.size(); i++) {
      //std::cout << "Element number " << "=" << i << "~" << adc[i] << "corresponding to channel number" << "==" << bar << std::endl;
      hADC[bar]->Fill(adc[i]);
      hADC_Total->Fill(adc[i]);
      hADCvsTS[bar]->SetMarkerStyle(45);
      hADCvsTS[bar]->SetMarkerColor(kBlue);
      hADCvsTS[bar]->SetMarkerSize(0.3);
      hADCvsTS[bar]->Fill(i,adc[i]);
      //if (adc[i] > 0) hADC_GT0[bar]->Fill(adc[i]);
      S=S+1;
    }

    for(unsigned int i=0; i<tdc.size(); i++) {
      //std::cout << "Element number " << "=" << i << "~" << adc[i] << "corresponding to channel number" << "==" << bar << std::endl;
      hTDC[bar]->Fill(tdc[i]);
      hTDC_Total->Fill(tdc[i]);
      hTDCvsTS[bar]->SetMarkerStyle(45);
      hTDCvsTS[bar]->SetMarkerColor(kBlue);
      hTDCvsTS[bar]->SetMarkerSize(0.3);
      hTDCvsTS[bar]->Fill(i,tdc[i]);
      S1=S1+1;
    }
    
    //std::cout << adc.size() << ":" << tdc.size() << std::endl;
    //std::cout << bar << ":" << adc[bar] << std::endl;
    //std::cout << "Random ADC counts for channel number = " << bar << ":" << adc[bar] << std::endl;
    //std::cout << adc[0] << std::endl; 
    //int adc = chan.getADC();
  }
  //hPEmaxVsDelta->Fill(leadBar-subleadBar, peLead);
  //std::cout << "Total number of entries in file based on ADC : " << S << std::endl;
  //std::cout << "Total number of entries in file based on TDC : " << S1 << std::endl;
  return;
  }

  void ADCAnalyzer::onFileOpen() {
    std::cout << "\n\n File is opening! My analyzer should do something -- like print this \n\n" << std::endl;
    //std::cout << "Total number of entries in file : " << S << std::endl;
    return;
  }

  void ADCAnalyzer::onFileClose() {
    //std::cout << "Total number of entries in file : " << S << std::endl;
    return;
  }
  
  void ADCAnalyzer::onProcessStart() {
    std::cout << "\n\n Process starts! My analyzer should do something -- like print this \n\n" << std::endl;

    getHistoDirectory();

  
  //int nTimeSamp=40;
  //int PEmax=400;
  //int nPEbins=2*PEmax;
  //float Qmax=PEmax/(6250./4.e6);
  //float Qmin=-100;
  //int nQbins=(Qmax-Qmin)/4;
  
  for (int iB=0; iB<nChannels; iB++) {
    //hADC[iB] = new TH1F(Form("hADC_chan%i",iB), Form("; ADC, chan%i", iB),20,0,240);
    hADC[iB] = new TH1F(Form("hADC_chan%i",iB), Form("ADC counts for chan%i; ADC", iB),80,0,240);
    hADCvsTS[iB] = new TH2F(Form("hADCvsTS_chan%i",iB),Form("ADC counts vs timesample for chan%i; Timesample ; ADC", iB),30,0,29,120,0,240);
    hTDCvsTS[iB] = new TH2F(Form("hTDCvsTS_chan%i",iB),Form("TDC counts vs timesample for chan%i; Timesample ; TDC", iB),30,0,29,40,60,70);
    hTDC[iB] = new TH1F(Form("hTDC_chan%i",iB), Form("; TDC, chan%i", iB),20,60,70);
    hADC_Total = new TH1F(Form("hADC_Total"),Form("; Total ADC,"),80,0,240);
    hTDC_Total = new TH1F(Form("hTDC_Total"),Form("; Total TDC,"),20,60,70);
  }

  fillNb=0;
  evNb=0;

  //for (int iB = 0; iB<nChannels; iB++) {
    //hPE[iB]=new TH1F(Form("hPE_chan%i", iB), Form(";PE, chan%i", iB),nPEbins,0,PEmax);
    //hPEVsDelta[iB]=new TH2F(Form("hPEVsDelta_chan%i", iB), Form(";#Delta_{barID};PE, chan%i has max PE", iB),nChannels+1,-nChannels/2-0.5,nChannels/2+0.5, nPEbins,0,PEmax);
    //hDeltaPEVsDelta[iB]=new TH2F(Form("hDeltaPEVsDelta_chan%i", iB), Form(";#Delta_{barID};#Delta_PE, chan%i has max PE", iB),nChannels+1,-nChannels/2-0.5,nChannels/2+0.5, nPEbins,0,PEmax);
  //}

  //make event displays for events where there are channels with weird event hit PE counts 
  //for (int iE = 0; iE<nEv; iE++) {
    //for (int iB = 0; iB<nChannels; iB++) { 
    //hOut[iE][iB] = new TH1F(Form("hCharge_chan%i_nv%i", iB, iE), Form(";time sample; Q, chan %i, ev %i [fC]", iB, iE), nTimeSamp,-0.5,nTimeSamp-0.5);
    //}
  //}
  
  //hPEmaxVsDelta=new TH2F("hPEmaxVsDelta",";#Delta_{barID};PE, max hit",nChannels,-nChannels/2,nChannels/2, nPEbins,0,PEmax);
  //hEvDisp = new TH2F(Form("hEvDisp_ev%i", nEv), ";Event number; Bar ID; PE", nEv,0.5,nEv+0.5, nChannels,-0.5,nChannels-0.5);
  //hEvDispPE = new TH2F("hEvDispPEcut", ";Event number; Bar ID; PE", nEv,0.5,nEv+0.5, nChannels,-0.5,nChannels-0.5);
  

  //fillNb = 0;
  //evNb = 0;
  
    return;
  }
  

  void ADCAnalyzer::onProcessEnd() {

    return;
  }


}

DECLARE_ANALYZER_NS(trigscint, ADCAnalyzer)
