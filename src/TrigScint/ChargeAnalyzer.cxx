#include "TrigScint/ChargeAnalyzer.h"

namespace trigscint {

  ChargeAnalyzer::ChargeAnalyzer(const std::string& name,
               framework::Process& process)
    : Analyzer(name, process) {}
  ChargeAnalyzer::~ChargeAnalyzer() {}
  
  void ChargeAnalyzer::configure(framework::config::Parameters &parameters){

    inputCol_  = parameters.getParameter< std::string >("inputCollection");
    inputPassName_  = parameters.getParameter< std::string >("inputPassName");
    peds_  = parameters.getParameter< std::vector<double> >("pedestals");
    gain_  = parameters.getParameter< std::vector<double> >("gain");
    startSample_  = parameters.getParameter< int >("startSample");

    std::cout << " [ QIEAnalyzer ] In configure(), got parameters " 
        << "\n\t inputCollection = " << inputCol_
        << "\n\t inputPassName = " << inputPassName_
        << "\n\t startSample = " << startSample_
        << "\n\t pedestals[0] = " << peds_[0]
        << "\n\t gain[0] = " << gain_[0]
        << "\t." << std::endl;

    return;
  }

  void ChargeAnalyzer::analyze(const framework::Event &event) {

    //std::cout << "getting collection" << std::endl;
    const auto channels{event.getCollection<trigscint::EventReadout>(inputCol_, inputPassName_)};
    //std::cout << "got..." << std::endl;
  int evNb = event.getEventNumber();
  int nChannels = channels.size();
  //int S=0;
  
  for (auto chan : channels) {
    //int bar = chan.getBarID();
    int bar = chan.getChanID();
    std::vector<float> q = chan.getQ();
    float QTot = chan.getTotQ();
    float QAvg = chan.getAvgQ();
    float QMin = chan.getMinQ();
    float QMax = chan.getMaxQ();
    float QMed = chan.getMedQ();

    hQTot->Fill(QTot);
    if (QTot < Qlow_thr) {
      hQTotvschan_low->SetMarkerStyle(45);
      hQTotvschan_low->SetMarkerColor(kBlue);
      hQTotvschan_low->SetMarkerSize(0.3);
      hQTot_low[bar]->Fill(QTot);
      hQTotvschan_low->Fill(bar,QTot);
    }
    if (QTot > Qlow_thr && QTot < Qmed_thr) {
      hQTotvschan_med->SetMarkerStyle(45);
      hQTotvschan_med->SetMarkerColor(kBlue);
      hQTotvschan_med->SetMarkerSize(0.3);
      hQTot_med[bar]->Fill(QTot);
      hQTotvschan_med->Fill(bar,QTot);
    }
    if (QTot > Qmed_thr) {
      hQTotvschan_high->SetMarkerStyle(45);
      hQTotvschan_high->SetMarkerColor(kBlue);
      hQTotvschan_high->SetMarkerSize(0.3);
      hQTot_high[bar]->Fill(QTot);
      hQTotvschan_high->Fill(bar,QTot);
    }
    //hQTotvsTS[bar]->Fill()

    if (QAvg < Qlow_thr) {
      hQAvgvschan_low->SetMarkerStyle(45);
      hQAvgvschan_low->SetMarkerColor(kBlue);
      hQAvgvschan_low->SetMarkerSize(0.3);
      hQAvg_low[bar]->Fill(QAvg);
      hQAvgvschan_low->Fill(bar,QAvg);
    }
    if (QAvg > Qlow_thr && QAvg < Qmed_thr) {
      hQAvgvschan_med->SetMarkerStyle(45);
      hQAvgvschan_med->SetMarkerColor(kBlue);
      hQAvgvschan_med->SetMarkerSize(0.3);
      hQAvg_med[bar]->Fill(QAvg);
      hQAvgvschan_med->Fill(bar,QAvg);
    }
    if (QAvg > Qmed_thr) {
      hQAvgvschan_high->SetMarkerStyle(45);
      hQAvgvschan_high->SetMarkerColor(kBlue);
      hQAvgvschan_high->SetMarkerSize(0.3);
      hQAvg_high[bar]->Fill(QAvg);
      hQAvgvschan_high->Fill(bar,QAvg);
    }

    if (QMin < Qlow_thr) {
      hQMinvschan_low->SetMarkerStyle(45);
      hQMinvschan_low->SetMarkerColor(kBlue);
      hQMinvschan_low->SetMarkerSize(0.3);
      hQMin_low[bar]->Fill(QMin);
      hQMinvschan_low->Fill(bar,QMin);
    }
    if (QMin > Qlow_thr && QMin < Qmed_thr) {
      hQMinvschan_med->SetMarkerStyle(45);
      hQMinvschan_med->SetMarkerColor(kBlue);
      hQMinvschan_med->SetMarkerSize(0.3);
      hQMin_med[bar]->Fill(QMin);
      hQMinvschan_med->Fill(bar,QMin);
    }
    if (QMin > Qmed_thr) {
      hQMinvschan_high->SetMarkerStyle(45);
      hQMinvschan_high->SetMarkerColor(kBlue);
      hQMinvschan_high->SetMarkerSize(0.3);
      hQMin_high[bar]->Fill(QMin);
      hQMinvschan_high->Fill(bar,QMin);
    }

    if (QMax < Qlow_thr) {
      hQMaxvschan_low->SetMarkerStyle(45);
      hQMaxvschan_low->SetMarkerColor(kBlue);
      hQMaxvschan_low->SetMarkerSize(0.3);
      hQMax_low[bar]->Fill(QMax);
      hQMaxvschan_low->Fill(bar,QMax);
    }
    if (QMax > Qlow_thr && QMax < Qmed_thr) {
      hQMaxvschan_med->SetMarkerStyle(45);
      hQMaxvschan_med->SetMarkerColor(kBlue);
      hQMaxvschan_med->SetMarkerSize(0.3);
      hQMax_med[bar]->Fill(QMax);
      hQMaxvschan_med->Fill(bar,QMax);
    }
    if (QMax > Qmed_thr) {
      hQMaxvschan_high->SetMarkerStyle(45);
      hQMaxvschan_high->SetMarkerColor(kBlue);
      hQMaxvschan_high->SetMarkerSize(0.3);
      hQMax_high[bar]->Fill(QMax);
      hQMaxvschan_high->Fill(bar,QMax);
    }

    if (QMed < Qlow_thr) {
      hQMedvschan_low->SetMarkerStyle(45);
      hQMedvschan_low->SetMarkerColor(kBlue);
      hQMedvschan_low->SetMarkerSize(0.3);
      hQMed_low[bar]->Fill(QMed);
      hQMedvschan_low->Fill(bar,QMed);
    }
    if (QMed > Qlow_thr && QMed < Qmed_thr) {
      hQMedvschan_med->SetMarkerStyle(45);
      hQMedvschan_med->SetMarkerColor(kBlue);
      hQMedvschan_med->SetMarkerSize(0.3);
      hQMed_med[bar]->Fill(QMed);
      hQMedvschan_med->Fill(bar,QMed);
    }
    if (QMed > Qmed_thr) {
      hQMedvschan_high->SetMarkerStyle(45);
      hQMedvschan_high->SetMarkerColor(kBlue);
      hQMedvschan_high->SetMarkerSize(0.3);
      hQMed_high[bar]->Fill(QMed);
      hQMedvschan_high->Fill(bar,QMed);
    }

      for(unsigned int i=0;i<q.size();i++) {
        hQvsTS_low[bar]->SetMarkerStyle(45);
        hQvsTS_low[bar]->SetMarkerColor(kBlue);
        hQvsTS_low[bar]->SetMarkerSize(0.3);
        hQvsTS_med[bar]->SetMarkerStyle(45);
        hQvsTS_med[bar]->SetMarkerColor(kBlue);
        hQvsTS_med[bar]->SetMarkerSize(0.3);
        hQvsTS_high[bar]->SetMarkerStyle(45);
        hQvsTS_high[bar]->SetMarkerColor(kBlue);
        hQvsTS_high[bar]->SetMarkerSize(0.3);
        if (q[i] < Qlow_thr) {
          hQ_low[bar]->Fill(q[i]);
          hQvsTS_low[bar]->Fill(i,q[i]);
        }
        if (q[i] > Qlow_thr && q[i]<Qmed_thr) {
          hQ_med[bar]->Fill(q[i]);
          hQvsTS_med[bar]->Fill(i,q[i]);
        }
        if (q[i] > Qmed_thr) {
          hQ_high[bar]->Fill(q[i]);
          hQvsTS_high[bar]->Fill(i,q[i]);
        }
      }
    
  }
  return;
  }

  void ChargeAnalyzer::onFileOpen() {
    std::cout << "\n\n File is opening! My analyzer should do something -- like print this \n\n" << std::endl;
    //std::cout << "Total number of entries in file : " << S << std::endl;
    return;
  }

  void ChargeAnalyzer::onFileClose() {
    //std::cout << "Total number of entries in file : " << S << std::endl;
    return;
  }
  
  void ChargeAnalyzer::onProcessStart() {
    std::cout << "\n\n Process starts! My analyzer should do something -- like print this \n\n" << std::endl;

    getHistoDirectory();

  
  //int nTimeSamp=40;
  int PEmax=400;
  //int nPEbins=2*PEmax;
  //float Qmax=PEmax/(6250./4.e6);
  //float Qmin=-100;
  //int nQbins=(Qmax-Qmin)/4;
  int nQbins_low = Qlow_thr;
  int nQbins_med = (Qmed_thr-Qlow_thr);
  int nQbins_high = (Qhigh_thr-Qmed_thr);

  for (int iB=0; iB<nChannels; iB++) {
    hQvsTS_low[iB] = new TH2F(Form("hQvsTS_low_chan%i",iB),Form("Charge vs timesample for chan%i (Q < Qlow_thr); Timesample ; Q", iB),30,0,29,nQbins_low,-100,Qlow_thr);
    hQvsTS_med[iB] = new TH2F(Form("hQvsTS_med_chan%i",iB),Form("Charge vs timesample for chan%i; Timesample ; Q", iB),30,0,29,nQbins_med/10,Qlow_thr,Qmed_thr);
    hQvsTS_high[iB] = new TH2F(Form("hQvsTS_high_chan%i",iB),Form("Charge vs timesample for chan%i; Timesample ; Q", iB),30,0,29,nQbins_high/15,Qmed_thr,Qhigh_thr);
    hQ_low[iB] = new TH1F(Form("hQ_low_chan%i",iB), Form("Qs for chan%i (Q < Qlow_thr); Q", iB),nQbins_low,-100,Qlow_thr);
    hQ_med[iB] = new TH1F(Form("hQ_med_chan%i",iB), Form("Qs for chan%i (Qlow_thr < Q < Qmed_thr); Q", iB),nQbins_med/10,Qlow_thr,Qmed_thr);
    hQ_high[iB] = new TH1F(Form("hQ_high_chan%i",iB), Form("Qs for chan%i (Q > Qmed_thr); Q", iB),nQbins_high/15,Qmed_thr,Qhigh_thr);
    hQTot_low[iB] = new TH1F(Form("hQTot_low_chan%i",iB), Form("Total Q for chan%i (Q < Qlow_thr); QTot", iB),nQbins_low,-100,Qlow_thr);
    hQTot_med[iB] = new TH1F(Form("hQTot_med_chan%i",iB), Form("Total Q for chan%i (Qlow_thr < Q < Qmed_thr); QTot", iB),nQbins_med/10,Qlow_thr,Qmed_thr);
    hQTot_high[iB] = new TH1F(Form("hQTot_high_chan%i",iB), Form("Total Q for chan%i (Q > Qmed_thr); QTot", iB),nQbins_high/15,Qmed_thr,Qhigh_thr);
    hQAvg_low[iB] = new TH1F(Form("hQAvg_low_chan%i",iB), Form("Average Q for chan%i (Q < Qlow_thr); QAvg", iB),nQbins_low,-5.e2,Qlow_thr);
    hQAvg_med[iB] = new TH1F(Form("hQAvg_med_chan%i",iB), Form("Average Q for chan%i (Qlow_thr < Q < Qmed_thr); QAvg", iB),nQbins_med/10,Qlow_thr,Qmed_thr);
    hQAvg_high[iB] = new TH1F(Form("hQAvg_high_chan%i",iB), Form("Average Q for chan%i (Q > Qmed_thr); QAvg", iB),nQbins_high/15,Qmed_thr,Qhigh_thr);
    hQMin_low[iB] = new TH1F(Form("hQMin_low_chan%i",iB), Form("Minimum Q for chan%i (Q < Qlow_thr); QMin", iB),nQbins_low,-5.e2,Qlow_thr);
    hQMin_med[iB] = new TH1F(Form("hQMin_med_chan%i",iB), Form("Minimum Q for chan%i (Qlow_thr < Q < Qmed_thr); QMin", iB),nQbins_med/10,Qlow_thr,Qmed_thr);
    hQMin_high[iB] = new TH1F(Form("hQMin_high_chan%i",iB), Form("Minimum Q for chan%i (Q > Qmed_thr); QMin", iB),nQbins_high/15,Qmed_thr,Qhigh_thr);
    hQMax_low[iB] = new TH1F(Form("hQMax_low_chan%i",iB), Form("Maximum Q for chan%i (Q < Qlow_thr); QMax", iB),nQbins_low,-5.e2,Qlow_thr);
    hQMax_med[iB] = new TH1F(Form("hQMax_med_chan%i",iB), Form("Maximum Q for chan%i (Qlow_thr < Q < Qmed_thr); QMax", iB),nQbins_med/10,Qlow_thr,Qmed_thr);
    hQMax_high[iB] = new TH1F(Form("hQMax_high_chan%i",iB), Form("Maximum Q for chan%i (Q > Qmed_thr); QMax", iB),nQbins_high/15,Qmed_thr,Qhigh_thr);
    hQMed_low[iB] = new TH1F(Form("hQMed_low_chan%i",iB), Form("Median Q for chan%i (Q < Qlow_thr); QMed", iB),nQbins_low,-5.e2,Qlow_thr);
    hQMed_med[iB] = new TH1F(Form("hQMed_med_chan%i",iB), Form("Median Q for chan%i (Qlow_thr < Q < Qmed_thr); QMed", iB),nQbins_med/10,Qlow_thr,Qmed_thr);
    hQMed_high[iB] = new TH1F(Form("hQMed_high_chan%i",iB), Form("Median Q for chan%i (Q > Qmed_thr); QMed", iB),nQbins_high/15,Qmed_thr,Qhigh_thr);
    hQTotvschan_low = new TH2F(Form("hQTotvschan_low"),Form("QTot vs channel (Q < Qlow_thr); channel ; QTot"),16,0,15,nQbins_low,-100,Qlow_thr);
    hQTotvschan_med = new TH2F(Form("hQTotvschan_med"),Form("QTot vs channel (Qlow_thr < Q < Qmed_thr) ; channel ; QTot"),16,0,15,nQbins_med/10,Qlow_thr,Qmed_thr);
    hQTotvschan_high = new TH2F(Form("hQTotvschan_high"),Form("QTot vs channel (Q > Qmed_thr); channel ; QTot"),16,0,15,nQbins_high/15,Qmed_thr,Qhigh_thr);
    hQAvgvschan_low = new TH2F(Form("hQAvgvschan_low"),Form("QAvg vs channel (Q < Qlow_thr); channel ; QAvg"),16,0,15,nQbins_low,-5.e2,Qlow_thr);
    hQAvgvschan_med = new TH2F(Form("hQAvgvschan_med"),Form("QAvg vs channel (Qlow_thr < Q < Qmed_thr); channel ; QAvg"),16,0,15,nQbins_med/10,Qlow_thr,Qmed_thr);
    hQAvgvschan_high = new TH2F(Form("hQAvgvschan_high"),Form("QAvg vs channel (Q > Qmed_thr); channel ; QAvg"),16,0,15,nQbins_high/15,Qmed_thr,Qhigh_thr);
    hQMinvschan_low = new TH2F(Form("hQMinvschan_low"),Form("QMax vs channel (Q < Qlow_thr); channel (Q < Qlow_thr); QMax"),16,0,15,nQbins_low,-5.e2,Qlow_thr);
    hQMinvschan_med = new TH2F(Form("hQMinvschan_med"),Form("QMax vs channel (Qlow_thr < Q < Qmed_thr); channel ; QMax"),16,0,15,nQbins_med/10,Qlow_thr,Qmed_thr);
    hQMinvschan_high = new TH2F(Form("hQMinvschan_high"),Form("QMax vs channel (Q > Qmed_thr); channel ; QMax"),16,0,15,nQbins_high/15,Qmed_thr,Qhigh_thr);
    hQMaxvschan_low = new TH2F(Form("hQMaxvschan_low"),Form("QMin vs channel (Q < Qlow_thr); channel; QMin"),16,0,15,nQbins_low,-5.e2,Qlow_thr);
    hQMaxvschan_med = new TH2F(Form("hQMaxvschan_med"),Form("QMin vs channel (Qlow_thr < Q < Qmed_thr); channel ; QMin"),16,0,15,nQbins_med/10,Qlow_thr,Qmed_thr);
    hQMaxvschan_high = new TH2F(Form("hQMaxvschan_high"),Form("QMin vs channel (Q > Qmed_thr); channel ; QMin"),16,0,15,nQbins_high/15,Qmed_thr,Qhigh_thr);
    hQMedvschan_low = new TH2F(Form("hQMedvschan_low"),Form("QMed vs channel (Q < Qlow_thr) ; channel ; QMed"),16,0,15,nQbins_low,-5.e2,Qlow_thr);
    hQMedvschan_med = new TH2F(Form("hQMedvschan_med"),Form("QMed vs channel (Qlow_thr < Q < Qmed_thr); channel ; QMed"),16,0,15,nQbins_med/10,Qlow_thr,Qmed_thr);
    hQMedvschan_high = new TH2F(Form("hQMedvschan_high"),Form("QMed vs channel (Q > Qmed_thr); channel ; QMed"),16,0,15,nQbins_high/15,Qmed_thr,Qhigh_thr);
    hQTot = new TH1F(Form("hQTotal"), Form("Total charge for all channels"),500,-500,500.e3);
  }

  fillNb=0;
  evNb=0;

  
    return;
  }
  

  void ChargeAnalyzer::onProcessEnd() {

    return;
  }


}

DECLARE_ANALYZER_NS(trigscint, ChargeAnalyzer)