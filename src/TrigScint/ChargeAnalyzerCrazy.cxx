#include "TrigScint/ChargeAnalyzerCrazy.h" 
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TStyle.h"
namespace trigscint {

  ChargeAnalyzerCrazy::ChargeAnalyzerCrazy(const std::string& name,
               framework::Process& process)
    : Analyzer(name, process) {}
  ChargeAnalyzerCrazy::~ChargeAnalyzerCrazy() {}
  
  void ChargeAnalyzerCrazy::configure(framework::config::Parameters &parameters){

    inputCol_  = parameters.getParameter< std::string >("inputCollection");
    inputPassName_  = parameters.getParameter< std::string >("inputPassName");
    peds_  = parameters.getParameter< std::vector<double> >("pedestals");
    gain_  = parameters.getParameter< std::vector<double> >("gain");
    startSample_  = parameters.getParameter< int >("startSample");
    EventDisplayNumber_ = parameters.getParameter< int >("EventDisplayNumber");

    std::cout << " [ ChargeAnalyzerCrazy ] In configure(), got parameters " 
        << "\n\t inputCollection = " << inputCol_
        << "\n\t inputPassName = " << inputPassName_
        << "\n\t startSample = " << startSample_
        << "\n\t pedestals[0] = " << peds_[0]
        << "\n\t gain[0] = " << gain_[0]
        << "\n\t EventDisplayNumber = " << EventDisplayNumber_
        << "\t." << std::endl;

    return;
  }

  std::pair<std::vector<float>,std::vector<float>> GetMaximumCharge(std::vector<std::vector<float>> X) 
  {
    std::vector<float> TS;
    std::vector<float> Q;
    for(int iT=0;iT<30;iT++) {  //Looping around timesample, the loop will pick out charge values corresponding to time sample i
      std::vector<float> QiT;
      for(int j=0;j<X.size();j++) {
        std::vector<float> TQ = X[j];         //Loop to check all elements of the array
        if(TQ[0] == iT) {
          QiT.push_back(TQ[1]);
        }                                                      
      }
      float Qmax = *std::max_element(std::begin(QiT), std::end(QiT));
      TS.push_back(iT);
      Q.push_back(Qmax);
    }
    return make_pair(TS,Q);
  }

  void DrawGraph(std::vector<std::vector<float>> X, std::vector<std::vector<float>> Y, std::vector<std::vector<float>> Z, Int_t display)  //X is Timesample, Y is charge
  {
    std::cout << "Drae" << std::endl;
    TCanvas *c1 = new TCanvas();
    gStyle->SetTitleFontSize(0.08);
    gStyle->SetTitleSize(0.06,"x");
    gStyle->SetTitleSize(0.06,"y");
    c1->Divide(2,6);
    TGraph *gr[12];
    int ChFl=0;
    for(int g=0;g<12;g++)
    { 
      c1->cd(g+1);
      //gPad->SetLogy();
      gr[g] = new TGraph();
      std::vector<float> XTS = X[g];
      std::vector<float> YQ = Y[g];
      float Qmx = *std::max_element(std::begin(YQ),std::end(YQ));
      for(unsigned int i=0;i<XTS.size();i++) {
        gr[g]->SetPoint(gr[g]->GetN(),XTS[i],YQ[i]);
      }
      for(unsigned int f=0;f<Z.size();f++) {
        std::vector<float> Ini = Z[f];
        if (Ini[0] == g) {
          ChFl = int(Ini[1]);
        }
      }
      gr[g]->SetLineWidth(4);
      gr[g]->SetFillColorAlpha(kBlue,0.35);
      gr[g]->GetXaxis()->SetLabelSize(0.1);
      gr[g]->GetYaxis()->SetLabelSize(0.1);
      gr[g]->GetXaxis()->SetTitle("TimeSamples");
      //gr[g]->GetXaxis()->SetTitleSize(0.3);
      gr[g]->GetYaxis()->SetTitle("Charge[fC]");
      //gr[g]->GetYaxis()->SetTitleSize(0.3);
      gr[g]->SetTitle(Form("Channel%d, Flag value %d",g,ChFl));
      gr[g]->GetYaxis()->SetRangeUser(-1.0,13000);
      //gr[g]->SetLogy();
      //gr[g]->SetFillStyle(3005);
      //gr[g]->SetTitleSize(0.3);
      gr[g]->Draw("AB");
      gPad->SetGrid(1,1);
      gPad->Update();
    }
    TString outs = "/home/dhruvanshu/LDMX_Analysis_Files/PeakStartSamples/Noise/FlagBarPlot_" + std::to_string(display) + ".pdf";
    c1->Print(outs);  
  }
  void ChargeAnalyzerCrazy::analyze(const framework::Event &event) {

    const auto channels{event.getCollection<trigscint::EventReadout>(inputCol_,inputPassName_)};
    const auto chunnels{event.getCollection<trigscint::TestBeamHit>("testBeamHitsUp",inputPassName_)};
  evNb = event.getEventNumber();
  //int display = EventDisplayNumber_;
  int nChannels = channels.size();

  if (evNb == EventDisplayNumber_) {
    for (auto chun : chunnels) {
      int barI = chun.getBarID();
      float FlagHit = chun.getQualityFlag();

      std::cout << "Bar ID " << barI << " has hit flag value " << FlagHit << std::endl;
      std::vector<float> ChanFlag = {barI,FlagHit};
      FlagList.push_back(ChanFlag);
    }
    for (auto chan : channels) {

      int bar = chan.getChanID();
      int elec = chan.getElecID();
      float FlagQ = chan.getQualityFlag();

      std::cout << "Bar ID " << bar << " has channel flag value " << FlagQ << std::endl;
      //std::vector<float> ChanFlag = {bar,Flag};
      //FlagList.push_back(ChanFlag);
      std::vector<float> q = chan.getQ();
      //std::vector<float> Qts; //a vector storing time sample and charge
      for(unsigned int i=0;i<q.size();i++) {
        //float j = i;
        //Qts.push_back(i);  //Time sample as first element of the vector
        //Qts.push_back(q[i]); // Charge as second element of the vector

        if (bar==0) { 
          QChan0.push_back(q[i]);
          //std::cout << bar << "=" << Qts[0] << ":" << Qts[1] << std::endl;
        }
        if (bar==1) { 
          QChan1.push_back(q[i]);  
        }
        if (bar==2) { QChan2.push_back(q[i]);}
        if (bar==3) { 
          QChan3.push_back(q[i]);
          if (i < 25 && i > 14) {
            QL+=q[i];
          }
          if (i > 25) {
            QH+=q[i];
          }

        }
        if (bar==4) { QChan4.push_back(q[i]);}
        if (bar==5) { QChan5.push_back(q[i]);}
        if (bar==6) { QChan6.push_back(q[i]);}
        if (bar==7) { 
          QChan7.push_back(q[i]);
          //std::cout << "Charge value deposited at time sample " << i << " is " << q[i] << " fC" << std::endl;
        }
        if (bar==8) { QChan8.push_back(q[i]);}
        if (bar==9) { QChan9.push_back(q[i]);}
        if (bar==10) { QChan10.push_back(q[i]);}
        if (bar==11) { QChan11.push_back(q[i]);}
        if (bar==12) { QChan12.push_back(q[i]);}
        if (bar==13) { QChan13.push_back(q[i]);}
        if (bar==14) { QChan14.push_back(q[i]);}
        if (bar==15) { QChan15.push_back(q[i]);}
        //std::cout << bar << "=" << Qts[0] << ":" << Qts[1] << std::endl;
      }
      //Qts.clear();
    }
  } 
  //std::cout << "All vectors successfully developed" << std::endl;
  return;
  }

  void ChargeAnalyzerCrazy::onFileOpen() {
    std::cout << "\n\n This analyzer will try to plot 1D like distributions for charge vs timesample for all channels \n\n" << std::endl;
    return;
  }

  void ChargeAnalyzerCrazy::onFileClose() {
    //std::cout << "Total number of entries in file : " << S << std::endl;
    return;
  }
  
  void ChargeAnalyzerCrazy::onProcessStart() {
    std::cout << "\n\n Process starts! My analyzer should do something -- like print this \n\n" << std::endl;
    //EventsMiss.clear();
    //TGraph *gr1 = new TGraph(30,TS0,Qmax0);
    getHistoDirectory();
    return;
  }
  

  void ChargeAnalyzerCrazy::onProcessEnd() {
    std::vector<float> arr[12] = {QChan0,QChan1,QChan2,QChan3,QChan4,QChan5,QChan6,QChan7,QChan8,QChan9,QChan10,QChan11};
    //TSAll.push_back({0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29})
    for(unsigned int o=0;o<12;o++)
    {
      if (arr[o].size() == 0)
      {
        std::cout << "Channel number " << o << " is ety" << std::endl;
        QAll.push_back({0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0});
        TSAll.push_back({0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29}); //This is to ensure root doesnt crash in case of an empty input. 
      }                                                                                                     //Highly plausible considering all channels arent depositing charge from all events.    
      else 
      {
        //auto Q = GetMaximumCharge(arr[o]);
        //std::cout << "Got it" << std::endl;
        //std::vector<float> TSmax = Q.first;
        //std::vector<float> Qmax = Q.second;
        QAll.push_back(arr[o]);
        TSAll.push_back({0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29});
      }
      
      //TSmax.clear();
      //Qmax.clear();
    }
    //std::cout << "Checking the charge values before plotting" << std::endl;
    //std::sort(QAll[2].begin(), QAll[2].end());
    //for(int qI=0;qI < QAll[2].size(); qI++) {
      //std::cout << "Charge deposited on channel 2 at time sample " << qI << " is " << QAll[2][qI] + 1 << std::endl;
      //std::cout << "Charge deposited on channel 3 at time sample " << qI << " is " << QAll[3][qI] - 3 << std::endl;
    //}
    for(int f=0; f<FlagList.size() ;f++) {
      std::cout << "Channel ID " << FlagList[f][0] << " has flag " << FlagList[f][1] << std::endl; 
    }
    DrawGraph(TSAll,QAll,FlagList,EventDisplayNumber_);
    //DrawGraph(TSAll,QAll);
    std::cout << "Sum of charge for TS < 25 is " << QL << " fC" << std::endl;
    std::cout << "Sum of charge for TS > 25 is " << QH << " fC" << std::endl;
    return;
  }


}

DECLARE_ANALYZER_NS(trigscint, ChargeAnalyzerCrazy)
