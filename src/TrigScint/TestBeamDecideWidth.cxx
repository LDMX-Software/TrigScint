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

int GetMaximumIndex(std::vector<float> &X)
{
	float Qmax = *std::max_element(std::begin(X),std::end(X));
	int index = std::distance(std::begin(X),std::max_element(std::begin(X),std::end(X)));
	return index+12;
}

void TestBeamDecideWidth::analyze(const framework::Event &event) {
	const auto channels{event.getCollection<trigscint::EventReadout>(inputCol_,inputPassName_)};
	int evNb = event.getEventNumber();
	for(int start=0;start<widthlimit;start++) {  //For loop which loops through different start samples. Note that this entire loop is processed once for each event, and then it would go to next event
		int nChannels = channels.size();
		for (auto chan : channels) {
			int bar = chan.getChanID();
			int FlagQ = chan.getQualityFlag();
			if (bar >= 12) { continue; }  //Excluding channels beyond 12
			float ped = peds_.at(bar); 
			int startT = 12 + start;  
			int elec = chan.getElecID(); 
			float maxQ = -999.;
			float totSubtrQ = 0.;
			std::vector<float> q = chan.getQ();
			for (int iT=startT; iT < q.size(); iT++) { //Looping over time samples. It would take charge data depending on the integration width
				float subQ = q.at(iT) - ped; //Adding pedestal correction to the charge values
				if(iT - startT < maximumpulsewidth) { //Criteria for selecting timesamples restricted by the integration width. Note that pulsewidth is preset.
					if(subQ > maxQ) { maxQ=subQ ;} 
					if(subQ > 0) {
						std::cout << "Charge contributed is " << subQ << " for startsample " << startT << " by timesample " << iT << " for bar " << bar << std::endl;
						totSubtrQ+=subQ ;
					} 
				}
				else if (subQ < 0 || q.at(iT) < 0 )
					break;
			}
			//Next, pertaining to each event, respective total charge value is stored in vectors named TotQChan. For a given event, certain channels would have charge deposit and
			//the data in vectors is pushed back successively wrt the start sample decided by the startsample loop. This is done for all available startsamples for an event at a time.
			//Only when the total charge exceeds the cutoff, the vector is filled. Cutoff value is set to approc 10k fC. A 2D histogram is also plotted for totalcharge
			// vs the respective startsample.
			if (bar==0 and totSubtrQ >= totChargelimit) {  
				TotQChan0.push_back(totSubtrQ); 
				ChargevsStartSample[0]->Fill(startT,totSubtrQ);
				//std::cout << "Charge value recorded for channel " << bar << " deposited by event " << evNb << " for start sample " << start << std::endl;
				std::cout << "Total Charge calculated as " << totSubtrQ << " for start sample " << startT << " for bar " << bar << std::endl;
 			}
			if (bar==1 and totSubtrQ >= totChargelimit) { 
				TotQChan1.push_back(totSubtrQ);
				ChargevsStartSample[1]->Fill(startT,totSubtrQ);
				std::cout << "Total Charge calculated as " << totSubtrQ << " for start sample " << startT << " for bar " << bar << std::endl;
				//std::cout << "Charge value recorded for channel " << bar << " deposited by event " << evNb << " for start sample " << start << std::endl;
			}
			if (bar==2 and totSubtrQ >= totChargelimit) { 
				TotQChan2.push_back(totSubtrQ);
				ChargevsStartSample[2]->Fill(startT,totSubtrQ);
				std::cout << "Total Charge calculated as " << totSubtrQ << " for start sample " << startT << " for bar " << bar << std::endl;
				//std::cout << "Charge value recorded for channel " << bar << " deposited by event " << evNb << " for start sample " << start << std::endl;
			}
			if (bar==3 and totSubtrQ >= totChargelimit) { 
				TotQChan3.push_back(totSubtrQ);
				ChargevsStartSample[3]->Fill(startT,totSubtrQ);
				std::cout << "Total Charge calculated as " << totSubtrQ << " for start sample " << startT << " for bar " << bar << std::endl;
				//std::cout << "Charge value recorded for channel " << bar << " deposited by event " << evNb << " for start sample " << start << std::endl;
			}
			if (bar==4 and totSubtrQ >= totChargelimit) { 
				TotQChan4.push_back(totSubtrQ);
				ChargevsStartSample[4]->Fill(startT,totSubtrQ);
				std::cout << "Total Charge calculated as " << totSubtrQ << " for start sample " << startT << " for bar " << bar << std::endl;
				//std::cout << "Charge value recorded for channel " << bar << " deposited by event " << evNb << " for start sample " << start << std::endl;
			}
			if (bar==5 and totSubtrQ >= totChargelimit) { 
				TotQChan5.push_back(totSubtrQ);
				ChargevsStartSample[5]->Fill(startT,totSubtrQ);
				std::cout << "Total Charge calculated as " << totSubtrQ << " for start sample " << startT << " for bar " << bar << std::endl;
				//std::cout << "Charge value recorded for channel " << bar << " deposited by event " << evNb << " for start sample " << start << std::endl;
			}
			if (bar==6 and totSubtrQ >= totChargelimit) { 
				TotQChan6.push_back(totSubtrQ);
				ChargevsStartSample[6]->Fill(startT,totSubtrQ);
				std::cout << "Total Charge calculated as " << totSubtrQ << " for start sample " << startT << " for bar " << bar << std::endl;
				//std::cout << "Charge value " << totSubtrQ << " recorded for channel " << bar << " deposited by event " << evNb << " for start sample " << start << std::endl;
			}
			if (bar==7 and totSubtrQ >= totChargelimit) { 
				TotQChan7.push_back(totSubtrQ);
				ChargevsStartSample[7]->Fill(startT,totSubtrQ);
				std::cout << "Total Charge calculated as " << totSubtrQ << " for start sample " << startT << " for bar " << bar << std::endl;
				//std::cout << "Charge value recorded for channel " << bar << " deposited by event " << evNb << " for start sample " << start << std::endl;
			}
			if (bar==8 and totSubtrQ >= totChargelimit) { 
				TotQChan8.push_back(totSubtrQ);
				ChargevsStartSample[8]->Fill(startT,totSubtrQ);
				std::cout << "Total Charge calculated as " << totSubtrQ << " for start sample " << startT << " for bar " << bar << std::endl;
				//std::cout << "Charge value recorded for channel " << bar << " deposited by event " << evNb << " for start sample " << start << std::endl;
			}
			if (bar==9 and totSubtrQ >= totChargelimit) { 
				TotQChan9.push_back(totSubtrQ);
				ChargevsStartSample[9]->Fill(startT,totSubtrQ);
				std::cout << "Total Charge calculated as " << totSubtrQ << " for start sample " << startT << " for bar " << bar << std::endl;
				//std::cout << "Charge value recorded for channel " << bar << " deposited by event " << evNb << " for start sample " << start << std::endl;
			}
			if (bar==10 and totSubtrQ >= totChargelimit) { 
				TotQChan10.push_back(totSubtrQ);
				ChargevsStartSample[10]->Fill(startT,totSubtrQ);
				std::cout << "Total Charge calculated as " << totSubtrQ << " for start sample " << startT << " for bar " << bar << std::endl;
				//std::cout << "Charge value recorded for channel " << bar << " deposited by event " << evNb << " for start sample " << start << std::endl;
			}
			if (bar==11 and totSubtrQ >= totChargelimit) { 
				TotQChan11.push_back(totSubtrQ);
				ChargevsStartSample[11]->Fill(startT,totSubtrQ);
				std::cout << "Total Charge calculated as " << totSubtrQ << " for start sample " << startT << " for bar " << bar << std::endl;
				//std::cout << "Charge value recorded for channel " << bar << " deposited by event " << evNb << " for start sample " << start << std::endl;
			}
	}
}//Next, now that we have respective data for all startsamples for a given event, it is now analysed to obtain the startsample favored by the event to have maximum total charge
//thereby ensuring to have the integration window enclosing the pulse. The respective start sample obtained is filled in a histogram to compare in the end.
 std::vector<float> ChargeArray[12] = {TotQChan0,TotQChan1,TotQChan2,TotQChan3,TotQChan4,TotQChan5,TotQChan6,TotQChan7,TotQChan8,TotQChan9,TotQChan10,TotQChan11};
 //for(int iq=0;iq < ChargeArray[6].size();iq++) {
 	//std::cout << "Charge deposited by event " << evNb << " for start sample " << iq + 12 << " is " << ChargeArray[6][iq] << std::endl;
 //}
 for(int k=0;k<12;k++)
 {
 	 if (ChargeArray[k].size() != 0) {
 	 	int ChargeIndex = GetMaximumIndex(ChargeArray[k]); //Function to obtain startsample for maximum total charge
 	 	ChannelWidthData[k]->Fill(ChargeIndex); //Filling respective channel histogram with the startsample obtained in previous step
 	 	std::cout << "Event " << evNb << " analysed for channel " << k  << " giving startsample of " << ChargeIndex << std::endl;
 	 }
 }
 TotQChan0.clear();
 TotQChan1.clear(); 
 TotQChan2.clear();
 TotQChan3.clear();
 TotQChan4.clear();
 TotQChan5.clear(); //Emptying the charge vectors to make them ready for next event.
 TotQChan6.clear();
 TotQChan7.clear();
 TotQChan8.clear();
 TotQChan9.clear();
 TotQChan10.clear();
 TotQChan11.clear();
 //ChargeArray.Clear();
 //std::cout << "Event " << evNb << " properly analysed. Moving to next event" << std::endl;
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
	
   
  getHistoDirectory();
  for(int iH=0;iH<12;iH++) {
  	ChannelWidthData[iH] = new TH1F(Form("ChannelWidthData_chan%i",iH),Form("QTot data for pulse width for channel %i; StartSample; Number of events", iH),widthlimit,12,12+widthlimit);
  	ChargevsStartSample[iH] = new TH2F(Form("QvsSS_chan%i",iH),Form("Charge vs start sample for channel %i; StartSample; TotalCharge[fC]",iH),10,12,22,500,totChargelimit,highchargelimit);
  	ChargevsStartSample[iH]->SetMarkerStyle(45);
		ChargevsStartSample[iH]->SetMarkerColor(kBlue);
		ChargevsStartSample[iH]->SetMarkerSize(0.3);
  }
  fillNb=0;
  evNb=0;
  return;
}
  

void TestBeamDecideWidth::onProcessEnd() {
  return;
}

}

DECLARE_ANALYZER_NS(trigscint, TestBeamDecideWidth)