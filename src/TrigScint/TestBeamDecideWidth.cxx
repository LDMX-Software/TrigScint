//Analyzer File to find threshold pulsewidth for integration

#include "TrigScint/TestBeamDecideWidth.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLatex.h"
#include <numeric>
#include <fstream>
using namespace std;

namespace trigscint {

	TestBeamDecideWidth::TestBeamDecideWidth(const std::string& name,
						   framework::Process& process)
    : Analyzer(name, process) {}
  TestBeamDecideWidth::~TestBeamDecideWidth() {}
  
  void TestBeamDecideWidth::configure(framework::config::Parameters &parameters){

    inputCol_  = parameters.getParameter< std::string >("inputCollection");
    inputPassName_  = parameters.getParameter< std::string >("inputPassName");
    peds_  = parameters.getParameter< std::vector<double> >("pedestals");
    gain_  = parameters.getParameter< std::vector<double> >("gain");  //to do: vector
    //startSample_  = parameters.getParameter< int >("startSample");
    //nInstrumentedChannels_  = parameters.getParameter< int >("nInstrumentedChannels");

    std::cout << " [ TestBeamDecideWidth ] In configure(), got parameters " 
	      << "\n\t inputCollection = " << inputCol_
	      << "\n\t inputPassName = " << inputPassName_
	      //<< "\n\t outputCollection = " << outputCollection_
	      //<< "\n\t startSample = " << startSample_
	      //<< "\n\t pulseWidth = " << pulseWidth_
	     // << "\n\t pulseWidthLYSO = " << pulseWidthLYSO_
			  << "\n\t gain[0] = " << gain_[0] 
	      //<< "\n\t nInstrumentedChannels = " << nInstrumentedChannels_
	      //<< "\n\t doCleanHits = " << doCleanHits_
	      << "\n\t pedestals[0] = " << peds_[0]
	      //<< "\n\t MIPresponse[0] = " << MIPresponse_[0]
	      << "\t." << std::endl;

    return;
  }

std::pair<std::vector<float>,std::vector<float>> AnalyseWidthData(std::vector<std::vector<float>> X)   //A function defined which follows my initial approach of finding maximum of a charge vector
{
	std::vector<float> Widths;
	std::vector<float> TotQ;
	for(int iT=0;iT<X.size();iT++) {
		Widths.push_back(iT+12);
		std::vector<float> TQ = X[iT];
		float QMAX = *std::max_element(std::begin(TQ),std::end(TQ));
		TotQ.push_back(QMAX);
	}
	return make_pair(Widths,TotQ);
}

void DrawGraph(std::pair<std::vector<float>,std::vector<float>> Y, int title) //Takes the timesample and total charge calculated as input for different integration windows and plots a graph
{
	std::vector<float> PW = Y.first;
	std::vector<float> TotQ = Y.second;
	TCanvas *c1 = new TCanvas();
	auto gr1 = new TGraph();
	for(unsigned int k=0;k<PW.size();k++) {
		gr1->SetMarkerColor(4);
		gr1->SetMarkerStyle(20);
		gr1->SetTitle("Total Charge vs sample starting points;StartSample;Total Charge");
		gr1->SetPoint(gr1->GetN(),PW[k],TotQ[k]);
	}
	gr1->GetXaxis()->SetRangeUser(1,18);
	gr1->Draw("");
	TString outs = "/home/dhruvanshu/ldmx-sw/Channel" + std::to_string(title) + "PulseVariation.png";
	c1->Print(outs);

}

float MakeHistograms(std::vector<std::vector<float>> Z, TH1F* hist, int histtitle) //Plots histograms to depict number of events favoring a given integration window, and lastly returning the average window as output.
{
	TCanvas *c2 = new TCanvas();
	gStyle->SetLabelSize(0.048,"x");
	gStyle->SetLabelSize(0.048,"y");
	gStyle->SetTitleSize(0.045,"x");
	gStyle->SetTitleSize(0.045,"y");
	gStyle->SetOptStat(1111111);
	std::vector<float> AvgIndex;
	std::vector<float> WQ = Z[0];
	for(int qi=0; qi < WQ.size()-1; qi++){ 	//Looping over events corresponding to a channel
		std::vector<float> QIVector;
		for(int qj=0; qj < Z.size(); qj++) {  //Looping over pulse width data given an event selected by qi loop
			//std::cout << "CHecking for event number " << qi << " for width value " << qj << std::endl; 
			std::vector<float> QJdata = Z[qj];
			float QJValue = QJdata[qi+1];  
			//std::cout << "Charge deposited is " << QJValue << std::endl;
			QIVector.push_back(QJValue); //Storing the QTot data corresponding to a given event for different widths.
		}
		float QImax = *std::max_element(std::begin(QIVector),std::end(QIVector));
		int index = std::distance(std::begin(QIVector),std::max_element(std::begin(QIVector),std::end(QIVector))); //Returns index of the charge considered for a given event.
		//std::cout << "Channel " << histtitle << " data for event number " << qi << " with maximum total charge of " << QImax << " fC at startsample " << index << std::endl;
		hist->Fill(index+12);
		AvgIndex.push_back(index+12);
		QIVector.clear();
	}
	Float_t ymax = hist->GetMaximum();
	TLatex *text = new TLatex(19,ymax/1.5,"PulseWidth = 7");
	//text.SetTextSize(0.1);
	hist->SetLineWidth(2);
	hist->SetLineColor(kBlue);	//Attempts to beautify the histograms
	hist->Draw();
	hist->SetStats(111111);
	text->SetTextSize(0.04);
  text->Draw();
	TString histout = "/home/dhruvanshu/ldmx-sw/HistogramDataPulseWidth_Channel_" + std::to_string(histtitle) + ".png";
	c2->Print(histout);
	float Total = std::accumulate(std::begin(AvgIndex),std::end(AvgIndex),0.0);
	return Total/AvgIndex.size();

}
void TestBeamDecideWidth::analyze(const framework::Event &event) {

for(int width=0;width<widthlimit;width++) { //Starting with unit pulsewidth. Increasing the width and looping over it. The looping is weird, I assumed it will loop through
	//all events for a given integration window. However, it is actually picking one event and looping through all integration windows first :}
	const auto channels{event.getCollection<trigscint::EventReadout>(inputCol_,inputPassName_)};
	int evNb = event.getEventNumber();
	//std::cout << "Pulse width set to" << "=" << width << std::endl;
	//std::cout << "Processing event number:" << evNb << std::endl;
	int nChannels = channels.size();
	//std::cout << "Number of channels to analyse in this event = " << nChannels << std::endl;
	for (auto chan : channels) {
		int bar = chan.getChanID();
		if (bar >= 12) { continue; } 
		float ped = peds_.at(bar);
		int startT = 12 + width;  
		int elec = chan.getElecID();
		float maxQ = -999.;
		float totSubtrQ = 0.;
		std::vector<float> q = chan.getQ();
		for(int iT = startT; iT < q.size(); iT++) {
			float subQ=q.at(iT)-ped;
		  	                                                              //Loop for timesample of charge
			if(iT - startT < maximumpulsewidth) {  //Looping condition is different. We fix the pulsewidth, and vary the starting sample.
				if (subQ > maxQ)
				maxQ = subQ; //q.ait(iT)
				if(subQ > 0) 
				totSubtrQ+=subQ;
			}
			else if (subQ < 0 || q.at(iT) < 0 ) {
				break;
			}
		}
		//Because the looping criteria is fix event and vary width, I am storing the data as vectors corresponding to different integration window for a given channel. These input vectors
		//corresponding to a given window has respective total charge for a given event as inputs.
		if (bar==0 && totSubtrQ >= totChargelimit) {
			TotQChan0[width].push_back(totSubtrQ);
			//std::cout << "Trying to push back in the vector for event number ::" << evNb << ":: for pulse width = " << width << std::endl;
		}
		if (bar==1 && totSubtrQ >= totChargelimit) {TotQChan1[width].push_back(totSubtrQ);} /// Vector format : Channel(Width(Event data)).
		if (bar==2 && totSubtrQ >= totChargelimit) {TotQChan2[width].push_back(totSubtrQ);}
		if (bar==3 && totSubtrQ >= totChargelimit) {TotQChan3[width].push_back(totSubtrQ);}
		if (bar==4 && totSubtrQ >= totChargelimit) {TotQChan4[width].push_back(totSubtrQ);}
		if (bar==5 && totSubtrQ >= totChargelimit) {TotQChan5[width].push_back(totSubtrQ);}
		if (bar==6 && totSubtrQ >= totChargelimit) {TotQChan6[width].push_back(totSubtrQ);}
		if (bar==7 && totSubtrQ >= totChargelimit) {TotQChan7[width].push_back(totSubtrQ);}
		if (bar==8 && totSubtrQ >= totChargelimit) {TotQChan8[width].push_back(totSubtrQ);}
		if (bar==9 && totSubtrQ >= totChargelimit) {TotQChan9[width].push_back(totSubtrQ);}
		if (bar==10 && totSubtrQ >= totChargelimit) {TotQChan10[width].push_back(totSubtrQ);}
		if (bar==11 && totSubtrQ >= totChargelimit) {TotQChan11[width].push_back(totSubtrQ);}

		//if (bar==0) {
			//TotQChan0[width].push_back(totSubtrQ/width);
			//std::cout << "Trying to push back in the vector for event number ::" << evNb << ":: for pulse width = " << width << std::endl;
		//}
		//if (bar==1) {TotQChan1[width].push_back(totSubtrQ/width);}
		//if (bar==2) {TotQChan2[width].push_back(totSubtrQ/width);}
		//if (bar==3) {TotQChan3[width].push_back(totSubtrQ/width);}
		//if (bar==4) {TotQChan4[width].push_back(totSubtrQ/width);}
		//if (bar==5) {TotQChan5[width].push_back(totSubtrQ/width);}
		//if (bar==6) {TotQChan6[width].push_back(totSubtrQ/width);}
		//if (bar==7) {TotQChan7[width].push_back(totSubtrQ/width);}
		//if (bar==8) {TotQChan8[width].push_back(totSubtrQ/width);}
		//if (bar==9) {TotQChan9[width].push_back(totSubtrQ/width);}
		//if (bar==10) {TotQChan10[width].push_back(totSubtrQ/width);}
		//if (bar==11) {TotQChan11[width].push_back(totSubtrQ/width);}

	}
PulseWidth.push_back(12+width); 
//QChan0.clear();
}

return;
}

void TestBeamDecideWidth::onFileOpen() {
    //std::cout << "\n\n This analyzer will try to plot 1D like distributions for charge vs timesample for all channels \n\n" << std::endl;
    return;
  }

void TestBeamDecideWidth::onFileClose() {
    //std::cout << "Total number of entries in file : " << S << std::endl;
  return;
}
  
void TestBeamDecideWidth::onProcessStart() {
  std::cout << "\n\n Process starts! My analyzer should do something -- like print this \n\n" << std::endl;
  for (int i=0;i<widthlimit;i++) {
		TotQChan0.push_back({0});
		TotQChan1.push_back({0});
		TotQChan2.push_back({0});  //Initiating composite charge vectors 
		TotQChan3.push_back({0});
		TotQChan4.push_back({0});
		TotQChan5.push_back({0});
		TotQChan6.push_back({0});
		TotQChan7.push_back({0});
		TotQChan8.push_back({0});
		TotQChan9.push_back({0});
		TotQChan10.push_back({0});
		TotQChan11.push_back({0});
	}
	for(int iH=0;iH<12;iH++) {
  	ChannelWidthData[iH] = new TH1F(Form("ChannelWidthData_chan%i",iH),Form("QTot data for pulse width for channel %i; StartSample", iH),widthlimit,12,12+widthlimit);
  }
    //EventsMiss.clear();
    //TGraph *gr1 = new TGraph(30,TS0,Qmax0);
  getHistoDirectory();
  fillNb=0;
  evNb=0;
  return;
}
  

void TestBeamDecideWidth::onProcessEnd() {
	fstream file;
	file.open("/home/dhruvanshu/ldmx-sw/IntegrationWidthData.txt",ios::out); //Defining the output file to store the startsamples for integration windows for all channels.
	auto QWidth0 = AnalyseWidthData(TotQChan0);
	DrawGraph(QWidth0,0);
	auto QWidth1 = AnalyseWidthData(TotQChan1);
	DrawGraph(QWidth1,1);
	auto QWidth2 = AnalyseWidthData(TotQChan2);
	DrawGraph(QWidth2,2);
	auto QWidth3 = AnalyseWidthData(TotQChan3);
	DrawGraph(QWidth3,3);
	auto QWidth4 = AnalyseWidthData(TotQChan4);
	DrawGraph(QWidth4,4);
	auto QWidth5 = AnalyseWidthData(TotQChan5);
	DrawGraph(QWidth5,5);
	auto QWidth6 = AnalyseWidthData(TotQChan6);
	DrawGraph(QWidth6,6);
	auto QWidth7 = AnalyseWidthData(TotQChan7);
	DrawGraph(QWidth7,7);
	auto QWidth8 = AnalyseWidthData(TotQChan8);
	DrawGraph(QWidth8,8);
	auto QWidth9 = AnalyseWidthData(TotQChan9);
	DrawGraph(QWidth9,9);
	auto QWidth10 = AnalyseWidthData(TotQChan10);
	DrawGraph(QWidth10,10);
	auto QWidth11 = AnalyseWidthData(TotQChan11);
	DrawGraph(QWidth11,11);
	std::vector<std::vector<float>> HistArray[12] = {TotQChan0,TotQChan1,TotQChan2,TotQChan3,TotQChan4,TotQChan5,TotQChan6,TotQChan7,TotQChan8,TotQChan9,TotQChan10,TotQChan11};
  for(int hii=0;hii<12;hii++)
  {
  	if(HistArray[hii].size() == 0) {
  		std::cout << "Channel " << hii << " data is empty" << std::endl;
  	}
  	else{
  		std::cout << "Making final histogram for channel " << hii << std::endl;
  		float AvgOutput = MakeHistograms(HistArray[hii],ChannelWidthData[hii],hii);
  		if (hii != 8) { file << hii << "," << int(AvgOutput) << "\n"; }  //Well, channel 8 is broken so not storing its info.		
  		std::cout << "Histogram made" << std::endl;
  	}
  }
  std::cout << TotQChan0.size() << std::endl;
  file.close();
  return;
}

}

DECLARE_ANALYZER_NS(trigscint, TestBeamDecideWidth)