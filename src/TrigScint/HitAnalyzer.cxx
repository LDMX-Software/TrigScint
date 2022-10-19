/**
 * @file HitAnalyzer.cxx
 * @brief An analyzer to analyse the TestbeamHit collection
 * @author Dhruvanshu Parmar, Texas Tech University
*/

#include "TrigScint/HitAnalyzer.h"

namespace trigscint {

  HitAnalyzer::HitAnalyzer(const std::string& name,
						   framework::Process& process)
    : Analyzer(name, process) {}
  HitAnalyzer::~HitAnalyzer() {}
  
  void HitAnalyzer::configure(framework::config::Parameters &parameters){

    inputCol_  = parameters.getParameter< std::string >("inputCollection");
    inputPassName_  = parameters.getParameter< std::string >("inputPassName");
    peds_  = parameters.getParameter< std::vector<double> >("pedestals");
    startSample_  = parameters.getParameter< int >("startSample");

    std::cout << " [ HitAnalyzer ] In configure(), got parameters " 
	      << "\n\t inputCollection = " << inputCol_
	      << "\n\t inputPassName = " << inputPassName_
	      << "\n\t startSample = " << startSample_
	      << "\n\t pedestals[0] = " << peds_[0]
	      << "\t." << std::endl;

    return;
  }

void HitAnalyzer::analyze(const framework::Event &event) {
	//gStyle->SetPalette(kBird);

	const auto channels{event.getCollection<trigscint::TestBeamHit>(inputCol_, inputPassName_)};

	int evNb = event.getEventNumber();
	int nChan = channels.size();

	for (auto chan : channels) {

		//int bar = chan.getChanID();
		int bar = chan.getBarID();
		float PE = chan.getPE();

		hBar->Fill(bar);
		hPE_all->Fill(PE);

		if (PE < PE_low) {
			hPE_low[bar]->Fill(PE);
			hPEvsbar_low->SetMarkerStyle(45);
			hPEvsbar_low->SetMarkerColor(kBlue);
			hPEvsbar_low->SetMarkerSize(0.3);
			hPEvsbar_low->Fill(bar,PE);
		}

		if (PE < PE_med1 && PE > PE_low) {
			hPE_med1[bar]->Fill(PE);
			hPEvsbar_med1->SetMarkerStyle(45);
			hPEvsbar_med1->SetMarkerColor(kBlue);
			hPEvsbar_med1->SetMarkerSize(0.3);
			hPEvsbar_med1->Fill(bar,PE);
		}

		if (PE < PE_med2 && PE > PE_med1) {
			hPE_med2[bar]->Fill(PE);
			hPEvsbar_med2->SetMarkerStyle(45);
			hPEvsbar_med2->SetMarkerColor(kBlue);
			hPEvsbar_med2->SetMarkerSize(0.3);
			hPEvsbar_med2->Fill(bar,PE);
		}

		if (PE > PE_med2) {
			hPE_high[bar]->Fill(PE);
			hPEvsbar_high->SetMarkerStyle(45);
			hPEvsbar_high->SetMarkerColor(kBlue);
			hPEvsbar_high->SetMarkerSize(0.3);
			hPEvsbar_high->Fill(bar,PE);
		}

		hPEvsbar->SetMarkerStyle(45);
		hPEvsbar->SetMarkerColor(kBlue);
		hPEvsbar->SetMarkerSize(0.3);
		hPEvsbar->Fill(bar,PE);

	}
	S=S+1;
	return;	
}

void HitAnalyzer::onFileOpen() {
	std::cout << "\n\n File is opening! My Hit Analyzer will start analysing ASAP \n\n" << std::endl;
	return;
}

void HitAnalyzer::onFileClose() {
	std::cout << "\n\n File analysis complete \n\n" << std::endl;
	std::cout << "\n\n Total numbe of events analysed : " << S << std::endl;
	return;
}

void HitAnalyzer::onProcessStart() {
	std::cout << "\n\n Process starts! My analyser should start plotting data for you \n\n" << std::endl;

	getHistoDirectory();

	int nPE_bins = 200;
	int PEmax = 1200;
	for (int iB=0; iB < nChannels; iB++) {
		hPE_low[iB] = new TH1F(Form("hPE_low_bar%i",iB),Form("PE for bar %i (Pedestal); PE ", iB),nPE_bins,0,PE_low);
		hPE_med1[iB] = new TH1F(Form("hPE_med1_bar%i",iB),Form("PE for bar %i (SinglePE); PE ", iB),nPE_bins,PE_low,PE_med1);
		hPE_med2[iB] = new TH1F(Form("hPE_med2_bar%i",iB),Form("PE for bar %i (MIP); PE ", iB),nPE_bins,PE_med1,PE_med2);
		hPE_high[iB] = new TH1F(Form("hPE_high_bar%i",iB),Form("PE for bar %i (RegionX); PE ", iB),nPE_bins,PE_med2,PEmax);
		hPE_all = new TH1F(Form("hPE"),Form("PE for all bars; PE"),nPE_bins,0,PEmax);
		hPEvsbar = new TH2F(Form("hPE2D"),Form("PE vs bars; bar; PE"),15,0,15,nPE_bins,0,PEmax);
		hPEvsbar_low = new TH2F(Form("hPE2D_low"),Form("PE vs bars (Pedestal); bar; PE"),15,0,15,nPE_bins,0,PE_low);
		hPEvsbar_med1 = new TH2F(Form("hPE2D_med1"),Form("PE vs bars (SinglePE); bar; PE"),15,0,15,nPE_bins,PE_low,PE_med1);
		hPEvsbar_med2 = new TH2F(Form("hPE2D_med2"),Form("PE vs bars (MIP); bar; PE"),15,0,15,nPE_bins,PE_med1,PE_med2);
		hPEvsbar_high = new TH2F(Form("hPE2D_high"),Form("PE vs bars (RegionX); bar; PE"),15,0,15,nPE_bins,PE_med2,PEmax);
		hBar = new TH1F(Form("Bar_ID"),Form("Bar ID ; bar"),15,0,15);
	}
	fillNb=0;
	evNb=0;

	return;
	}

void HitAnalyzer::onProcessEnd() {
	std::cout << "\n\n Process analysis completed. Finishing up \n\n" << std::endl;
	return;
}

}

DECLARE_ANALYZER_NS(trigscint, HitAnalyzer)