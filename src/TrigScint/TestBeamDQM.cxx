/**
 * @file TestBeamDQM.cxx
 * @brief An analyzer drawing QIE digi quantities
 * @author Andrew Whitbeck, Texas Tech University
 */

#include "TrigScint/TestBeamDQM.h"

namespace trigscint {

  TestBeamDQM::TestBeamDQM(const std::string& name,
						   framework::Process& process)
    : Analyzer(name, process) {}
  TestBeamDQM::~TestBeamDQM() {}
  
  void TestBeamDQM::configure(framework::config::Parameters &parameters){

    inputCol_  = parameters.getParameter< std::string >("inputCollection");
    inputPassName_  = parameters.getParameter< std::string >("inputPassName");
    peds_  = parameters.getParameter< std::vector<double> >("pedestals");
    startSample_  = parameters.getParameter< int >("startSample");

    std::cout << " [ TestBeamDQM ] In configure(), got parameters " 
	      << "\n\t inputCollection = " << inputCol_
	      << "\n\t inputPassName = " << inputPassName_
	      << "\n\t startSample = " << startSample_
	      << "\n\t pedestals[0] = " << peds_[0]
	      << "\t." << std::endl;

    return;
  }

  void TestBeamDQM::analyze(const framework::Event &event) {

    const auto channels{event.getCollection<trigscint::TrigScintQIEDigis >(inputCol_, inputPassName_)};

	int evNb = event.getEventNumber();
	int nChan = channels.size();
	
	for (auto chan : channels) {
	  int bar = chan.getBarID();
	  float PE = chan.getPE();
	  hPEVsDelta[leadBar]->Fill(leadBar-subleadBar, peLead);
	  hDeltaPEVsDelta[leadBar]->Fill(leadBar-subleadBar, peLead-peSublead);
	}
	hPEmaxVsDelta->Fill(leadBar-subleadBar, peLead);

	return;
  }

  void TestBeamDQM::onFileOpen() {
    std::cout << "\n\n File is opening! My analyzer should do something -- like print this \n\n" << std::endl;

    return;
  }

  void TestBeamDQM::onFileClose() {

    return;
  }
  
  void TestBeamDQM::onProcessStart() {
    std::cout << "\n\n Process starts! My analyzer should do something -- like print this \n\n" << std::endl;

    getHistoDirectory();

	
	int nTimeSamp=40;
	int PEmax=400;
	int nPEbins=2*PEmax;
	float Qmax=PEmax/(6250./4.e6);
	float Qmin=-10;
	int nQbins=(Qmax-Qmin)/4;
	
	for (int iB = 0; iB<nChannels; iB++) {
	  hPE[iB]=new TH1F(Form("hPE_chan%i", iB), Form(";PE, chan%i", iB),nPEbins,0,PEmax);
	  hPEVsDelta[iB]=new TH2F(Form("hPEVsDelta_chan%i", iB), Form(";#Delta_{barID};PE, chan%i has max PE", iB),nChannels+1,-nChannels/2-0.5,nChannels/2+0.5, nPEbins,0,PEmax);
	  hDeltaPEVsDelta[iB]=new TH2F(Form("hDeltaPEVsDelta_chan%i", iB), Form(";#Delta_{barID};#Delta_PE, chan%i has max PE", iB),nChannels+1,-nChannels/2-0.5,nChannels/2+0.5, nPEbins,0,PEmax);
	}

	//make event displays for events where there are channels with weird event hit PE counts 
	for (int iE = 0; iE<nEv; iE++) {
	  for (int iB = 0; iB<nChannels; iB++) { 
		hOut[iE][iB] = new TH1F(Form("hCharge_chan%i_nv%i", iB, iE), Form(";time sample; Q, chan %i, ev %i [fC]", iB, iE), nTimeSamp,-0.5,nTimeSamp-0.5);
	  }
	}
	
	hPEmaxVsDelta=new TH2F("hPEmaxVsDelta",";#Delta_{barID};PE, max hit",nChannels,-nChannels/2,nChannels/2, nPEbins,0,PEmax);
	hEvDisp = new TH2F(Form("hEvDisp_ev%i", nEv), ";Event number; Bar ID; PE", nEv,0.5,nEv+0.5, nChannels,-0.5,nChannels-0.5);
	hEvDispPE = new TH2F("hEvDispPEcut", ";Event number; Bar ID; PE", nEv,0.5,nEv+0.5, nChannels,-0.5,nChannels-0.5);
  

	fillNb = 0;
	evNb = 0;
	
    return;
  }
  

  void TestBeamDQM::onProcessEnd() {

    return;
  }


}

DECLARE_ANALYZER_NS(trigscint, TestBeamDQM)
