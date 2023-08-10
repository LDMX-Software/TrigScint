#include <stdio.h>
//#include <math.h>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <fstream>
#include <iostream>
#include <iomanip>

void setTDRStyle() {
  
  // For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(600); //Height of canvas
  gStyle->SetCanvasDefW(600); //Width of canvas
  gStyle->SetCanvasDefX(0);   //POsition on screen
  gStyle->SetCanvasDefY(0);

  // For the Pad:
  gStyle->SetPadBorderMode(0);
  // gStyle->SetPadBorderSize(Width_t size = 1);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  // For the frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(2);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);

  // For the histo:
  // gStyle->SetHistFillColor(1);
  // gStyle->SetHistFillStyle(0);
  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);
  // gStyle->SetLegoInnerR(Float_t rad = 0.5);
  // gStyle->SetNumberContours(Int_t number = 20);

  gStyle->SetEndErrorSize(2);
  //gStyle->SetErrorMarker(20);
  gStyle->SetErrorX(0.);
  
  gStyle->SetMarkerStyle(7);

  //For the fit/function:
  gStyle->SetOptFit(0);
  gStyle->SetFitFormat("5.4g");
  gStyle->SetFuncColor(2);
  gStyle->SetFuncStyle(1);
  gStyle->SetFuncWidth(1);

  gStyle->SetFuncWidth(1);

  //For the date:
  gStyle->SetOptDate(0);
  // gStyle->SetDateX(Float_t x = 0.01);
  // gStyle->SetDateY(Float_t y = 0.01);

  gStyle->SetOptStat(0);

  /*
  // For the statistics box:
  gStyle->SetOptFile(0);
  gStyle->SetOptStat("mr");
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.04);///---> gStyle->SetStatFontSize(0.025);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatH(0.1);
  gStyle->SetStatW(0.2);///---> gStyle->SetStatW(0.15);
  */

  // gStyle->SetStatStyle(Style_t style = 1001);
  // gStyle->SetStatX(Float_t x = 0);
  // gStyle->SetStatY(Float_t y = 0);

  // Margins:
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.04);

  // For the Global title:

  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);
  // gStyle->SetTitleH(0); // Set the height of the title box
  // gStyle->SetTitleW(0); // Set the width of the title box
  // gStyle->SetTitleX(0); // Set the position of the title box
  // gStyle->SetTitleY(0.985); // Set the position of the title box
  // gStyle->SetTitleStyle(Style_t style = 1001);
  // gStyle->SetTitleBorderSize(2);

  // For the axis titles:
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.06, "XYZ");
  // gStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // gStyle->SetTitleYSize(Float_t size = 0.02);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.25);
  // gStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:

  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);

  // Change for log plots:
  gStyle->SetOptLogx(0);
  gStyle->SetOptLogy(0);
  gStyle->SetOptLogz(0);

  // Postscript options:
  
  gStyle->SetPaperSize(20.,20.);
  // gStyle->SetLineScalePS(Float_t scale = 3);
  // gStyle->SetLineStyleString(Int_t i, const char* text);
  // gStyle->SetHeaderPS(const char* header);
  // gStyle->SetTitlePS(const char* pstitle);

  // gStyle->SetBarOffset(Float_t baroff = 0.5);
  // gStyle->SetBarWidth(Float_t barwidth = 0.5);
  // gStyle->SetPaintTextFormat(const char* format = "g");
  // gStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // gStyle->SetTimeOffset(Double_t toffset);
  // gStyle->SetHistMinimumZero(kTRUE);

  gStyle->SetPalette(1);

  //   gStyle->cd();

  gROOT->ForceStyle();
}

TH1* ChangedPlotterGeom(std::string filename){
    /*
     *You Should Only Have To Modify This Section
     * */
    int real=0;
    TTree *helper;
    if(real==0){
        TFile *myFile = TFile::Open(filename.c_str());
            //("unpacked_4gev_electrons_Apr02_0922_reformat_30timeSamplesFrom0_linearize_hits_clusters_dataShaped.root");
                //event.root");
        //"unpacked_4GeVelectrons_Apr1_1500_reformat_30timeSamplesFrom0_linearize_hits_clusters_dataShaped.root");
        //("unpacked_positive4GevElectrons_Mar28_1536_reformat_30timeSamplesFrom0_linearize_hits_clusters_dataShaped.root");
        //("out9nogap_dataShaped.root");
        //unpacked_ldmx_captan_out_22-04-2022_15-51-45__244_reformat_30timeSamplesFrom0_linearize_hits_clusters_dataShaped.root");
        //unpacked_ldmx_captan_out_13-04-2022_00-00-38__187_reformat_30timeSamplesFrom0_linearize_hits_clusters_dataShaped.root");
        
        
        TDirectory *dir = myFile->GetDirectory("tShaped");
        dir->GetObject("Helper",helper);
    }else{
        TFile *myFile = TFile::Open("out8.root");
        //myFile->Print("");
        myFile->GetObject("LDMX_Events",helper);
        //dangus->Print("");
        //dangus->GetObject("LDMX_Events",helper);
    }
        
    
    float IntCharge[17];
    float PENumber[17];
    float ChargeAtTimeForChan[17][30];
    float HitEfficiency[12];
    float HitEfficiencyRatio[12]; 
    //float barStatus[12]={2,1,0,0,0,0,-1,-2,3,2,4,-2};
    float barStatus[12]={2,1,0,0,0,0,0,0,0,0,0,0};
    TH1* Dist=new TH1F("PE","aoij",200,0.0,20.0);
    TH1* DivideR=new TH1F("DivR","pekeR",12.0,0.0,12.0);
    TH1* DivideL=new TH1F("DivL","pekeL",12.0,0.0,12.0);
    

    TH1* Divide=new TH1F("Div","peke",12.0,0.0,12.0);
    TH1* NumbersBotL=new TH1F("NBL","oiajw",12,0.0,12.0);
    TH1* NumbersTopL=new TH1F("Left Triggered Effficiency","Left Triggered Hit Efficiency",12,0.0,12.0);
    TH1* NumbersBotR=new TH1F("NBR","oiajw",12,0.0,12.0);
    TH1* NumbersTopR=new TH1F("Events Per Bar Above Threshold","Events Per Bar Above Threshold",12,0.0,12.0);
    float specialNum;
    
    if(real==0){
        helper->SetBranchAddress("IntCharge",IntCharge);
        //helper->SetBranchAddress("IntCharge[1]",&IntCharge[1]);
        helper->SetBranchAddress("PENumber",PENumber);
        //helper->SetBranchAddress("HitEfficiency",HitEfficiency);
    }else{
        //helper->SetBranchAddress("IntCharge",IntCharge);
        helper->SetBranchAddress("TriggerPadUpDigiHits_sim.pe_",PENumber);
        //helper->SetBranchAddress("HitEfficiency",HitEfficiency);
    }
        
    int counter=0.0;
    //Thresh is the tag probe threshold, thresh2 is the probe, and thresh3 is to ensure a neighboring cell isn't activated yielding probe activation.
    float thresh = 60.0;
    float thresh2 = 22.0;
    float thresh3 = 22.0;
    for(int i=0; i<12;i++){
        NumbersBotR->Fill((float)i);
        NumbersBotL->Fill((float)i);
    }
    float LeftB[12] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    float RightB[12] = {.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0};
    float LeftT[12] = {.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0};
    float RightT[12] = {.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0,.0};
    for(int i=0; i<helper->GetEntries(); i++){
        helper->GetEntry(i); 
        //Dist->Fill(PENumber[0]);
        for(int J=0;J<12;J++){          
            if((PENumber[J]>thresh) and (PENumber[J-1]>thresh) and (PENumber[J+1]<thresh2)){
                LeftT[J]++;
            }
            if((PENumber[J]>thresh) and (PENumber[J+1]>thresh) and (PENumber[J-1]<thresh2)){
                RightT[J]++;
            }
            if((PENumber[J]>thresh)){
                RightB[J]++;
                NumbersTopR->Fill(J);
            }
        }
    }
    for(int J =0;J<12;J++){
        std::cout<<LeftT[J]<<std::endl;
        std::cout<<LeftT[J]/RightB[J]<<std::endl;
        std::cout<<RightT[J]<<std::endl;
        std::cout<<RightT[J]/RightB[J]<<std::endl;
        std::cout<<RightB[J]<<std::endl;
        std::cout<<LeftT[J]/RightB[J]+RightT[J]/RightB[J]<<std::endl;
        std::cout<<"\n"<<std::endl;
    }
    
    //NumbersTopR->Draw("hist same e"); 
    //gPad->SaveAs("coolshite.png"); 
    NumbersTopR->Scale(4./NumbersTopR->Integral());
    return NumbersTopR;

}

TH1* ChangedPlotterMod(std::string filename,int MODE){
    /*
     *You Should Only Have To Modify This Section
     * */
    //c->cd();
    int real=0;
    TTree *helper;
    if(real==0){
        TFile *myFile = TFile::Open(filename.c_str());//("unpacked_positive4GevElectrons_Mar28_1455_reformat_30timeSamplesFrom0_linearize_hits_clusters_dataShaped.root");//unpacked_ldmx_captan_out_13-04-2022_00-00-38__187_reformat_30timeSamplesFrom0_linearize_hits_clusters_dataShaped.root");
        TDirectory *dir = myFile->GetDirectory("tShaped");
        dir->GetObject("Helper",helper);
    }else{
        TFile *myFile = TFile::Open("out9_dataShaped.root");
        //myFile->Print("");
        myFile->GetObject("LDMX_Events",helper);
        //dangus->Print("");
        //dangus->GetObject("LDMX_Events",helper);
    }
        
    
    float IntCharge[50];
    float PENumber[12];
    float ChargeAtTimeForChan[50][30];
    float HitEfficiency[12];
    float HitEfficiencyRatio[50]; 
    //float barStatus[12]={2,1,0,0,0,0,-1,-2,3,2,4,-2};
    float barStatus[12];//={2,1,0,0,0,0,0,0,0,0,0,0};
    for(int k=0;k<12;k++){
        barStatus[k]=0;
    }
    barStatus[0]=2;barStatus[1]=1;barStatus[11]=-2;barStatus[10]=-1;
    TH1* Dist=new TH1F("PE","aoij",200,0.0,20.0);
    TH1* DivideR=new TH1F("DivR","pekeR",12.0,0.0,12.0);
    TH1* DivideL=new TH1F("DivL","pekeL",12.0,0.0,12.0);
    TH1* HitsAboveThresh=new TH1F("Hits Above Threshold","Hits Above Threshold",12.0,0.0,12.0);
    TH1* DoubleHits=new TH1F("Number of Hit Pairs (i,i+2)","Number of Hit Pairs (i,i+2)",8.0,0.0,8.0);
        
    TH1* Divide=new TH1F("Div","peke",12.0,0.0,12.0);
    TH1* NumbersBotL=new TH1F("NBL","oiajw",12,0.0,12.0);
    TH1* NumbersTopL=new TH1F("Tag Probe Efficiency","Tag Probe Effficiency (Mar 30), Left Tag Red Right Blue",12,0.0,12.0);
    TH1* NumbersBotR=new TH1F("NBR","oiajw",12,0.0,12.0);
    TH1* NumbersTopR=new TH1F("Channel Efficiency","Channel Efficiency",12,0.0,12.0);
    float specialNum;
    
    if(real==0){
        helper->SetBranchAddress("IntCharge",IntCharge);
        //helper->SetBranchAddress("IntCharge[1]",&IntCharge[1]);
        helper->SetBranchAddress("PENumber",PENumber);
        //helper->SetBranchAddress("HitEfficiency",HitEfficiency);
    }else{
        //helper->SetBranchAddress("IntCharge",IntCharge);
        helper->SetBranchAddress("TriggerPadUpDigiHits_sim.pe_",PENumber);
        //helper->SetBranchAddress("HitEfficiency",HitEfficiency);
    }
        
    int counter=0.0;
    //Thresh is the tag probe threshold, thresh2 is the probe, and thresh3 is to ensure a neighboring cell isn't activated yielding probe activation.
    float thresh = 50.0;
    float thresh2 = 50.0;
    float thresh3 = 10;
    for(int i=0; i<50;i++){
        NumbersBotR->Fill((float)i);
        NumbersBotL->Fill((float)i);
    }
    for(int i=0; i<helper->GetEntries(); i++){
        helper->GetEntry(i); 
        Dist->Fill(PENumber[0]);
        float SUMS=0;
        float SUMS2=0;
        bool doSum2 = true;
        for(int J=0;J<12;J++){
            if(PENumber[J]>thresh){
                SUMS+=1.0;
                if(doSum2 and PENumber[J+1]>thresh){
                    J+=1;
                    SUMS2+=1;
                }
            }
        }
        HitsAboveThresh->Fill(SUMS);
        DoubleHits->Fill(SUMS2);
        for(int J=0;J<50;J++){
            if(barStatus[J]==3){
                //DEAD CHANNEL
                continue;
            }
            
            if(barStatus[J]==0){
                //Healthy channel surrounded by healthy channel
                if((PENumber[J-1]>thresh) and (PENumber[J-2]<thresh3) and (PENumber[J+1]<thresh3)){
                    NumbersBotL->Fill((float)J);
                    if(PENumber[J]>thresh2){
                        NumbersTopL->Fill((float)J);
                    }
                }
                if((PENumber[J+1]>thresh) and (PENumber[J+2]<thresh3) and (PENumber[J-1]<thresh3)){
                    NumbersBotR->Fill((float)J);
                    if(PENumber[J]>thresh2){
                        NumbersTopR->Fill((float)J);
                    }
                }
            }
            
            if(barStatus[J]==2){
                 //Dead channel immediattely to your left
                if((PENumber[J+1]>thresh) and (PENumber[J+2])<thresh3){
                    NumbersBotR->Fill((float)J);
                    if(PENumber[J]>thresh2){
                        NumbersTopR->Fill((float)J);
                    }
                } 
            }
            
            
            if(barStatus[J]==-2){
                //Dead channel immediattely to your right
                if((PENumber[J-1]>thresh) and (PENumber[J-2])<thresh3){
                    NumbersBotL->Fill((float)J);
                    if(PENumber[J]>thresh2){
                        NumbersTopL->Fill((float)J);
                    }
                }
            }
            
            
            if(barStatus[J]==1){
                //Dead Channel 2 away to your left
                if(PENumber[J-1]>thresh){
                    NumbersBotL->Fill((float)J);
                    if(PENumber[J]>thresh2){
                        NumbersTopL->Fill((float)J);
                    }
                } 
                if((PENumber[J+1]>thresh) and (PENumber[J+2]<thresh3)){
                    NumbersBotR->Fill((float)J);
                    if(PENumber[J]>thresh2){
                        NumbersTopR->Fill((float)J);
                    }
                }   
            }

            if(barStatus[J]==-1){
                //Dead Channel 2 away to your right
                if((PENumber[J-1]>thresh) and (PENumber[J-2])<thresh3){
                    NumbersBotL->Fill((float)J);
                    if(PENumber[J]>thresh2){
                        NumbersTopL->Fill((float)J);
                    }
                }
                if(PENumber[J+1]>thresh){
                    NumbersBotR->Fill((float)J);
                    if(PENumber[J]>thresh2){
                        NumbersTopR->Fill((float)J);
                    }
                } 
            }
            
            if(barStatus[J]==4){
                //Dead channel 2 away on your left and right.
                if(PENumber[J-1]>thresh){
                    NumbersBotL->Fill((float)J);
                    if(PENumber[J]>thresh2){
                        NumbersTopL->Fill((float)J);
                    }
                }
                if(PENumber[J+1]>thresh){
                    NumbersBotR->Fill((float)J);
                    if(PENumber[J]>thresh2){
                        NumbersTopR->Fill((float)J);
                    }
                } 
            }
        }
    }
    //Populating the scaling histogram
    float width=.125;
    float leftw=.62;
    for(int I = 0; I<100; I++){
        for(int J=0; J<50; J++){
            if(barStatus[J]==0){
                if(I<100*(1-width)){
                    DivideR->Fill(J);
                    DivideL->Fill(J);
                    Divide->Fill(J);
                }
            }
            else if(barStatus[J]==-2){
                if(I<100*(1-width)){
                    DivideR->Fill(J);
                    DivideL->Fill(J);
                    Divide->Fill(J);
                }
            }
            else if(barStatus[J]==2){
                if(I<100*(1-width)){
                    DivideR->Fill(J);
                    DivideL->Fill(J);
                    Divide->Fill(J);   
                }
            }
            
            else if(barStatus[J]==-1){
                if(I<100*(1-width)){
                    DivideL->Fill(J);
                    Divide->Fill(J);
                }
                if(I<100*(1-leftw-width)){
                    DivideR->Fill(J);
                }
            }
            else if(barStatus[J]==1){
                if(I<100*(1-width)){
                    DivideR->Fill(J);
                    Divide->Fill(J);
                }
                if(I<100*(leftw)){
                    DivideL->Fill(J);
                }
            }
            
            else if(barStatus[J]==4){
                if(I<90){
                    Divide->Fill(J);
                }
                if(I<42){
                    DivideL->Fill(J);
                }
                if(I<42){
                    DivideR->Fill(J);
                }
            }
        }
    }
    Divide->Scale(.01);
    DivideL->Scale(.01);
    //DivideL->Scale(2);
    DivideR->Scale(.01);
    //DivideR->Scale(2);
    /*
    TGraph *g = new TGraph();
    g->SetTitle("Max Efficiency for .3 mm Gaps,3 mm Bars Shifted by 1.5");
    g->GetYaxis()->SetTitle("Efficiency");
    g->GetXaxis()->SetTitle("Channel");
    float Helper[12]={.45,.9,.9,.9,.9,.9,.9,.45,0,.45,.9,.45};
    for(int J=0;J<12;J++){
        //if(not(J==8)){
        //    Helper[J]-=.1;
        //}
        g->SetPoint(2*J,float(J),Helper[J]);
        g->SetPoint(2*J+1,float(J)+1.0,Helper[J]);
    }
     */

    //TCanvas *c1 = new TCanvas("c1");
    NumbersTopR->GetYaxis()->SetTitle("Efficiency");
    NumbersTopL->GetXaxis()->SetTitle("Channel");
    std::cout<<counter<<std::endl;
    
    NumbersTopR->Divide(NumbersTopR,NumbersBotR,1.,1.,"B");
    //NumbersTopR->Divide(NumbersTopR,DivideR,1.,1.,"B");
    NumbersTopL->Divide(NumbersTopL,NumbersBotL,1.,1.,"B");
    //NumbersTopL->Divide(NumbersTopL,DivideL,1.,1.,"B");

    //NumbersTopR->Add(NumbersTopL);
    //NumbersTopR->Scale(.5);
    NumbersTopR->SetLineColor(kRed);
    NumbersTopL->SetLineColor(kBlue);
    //NumbersTopR->Divide(NumbersTopR,DivideR,1.,1.,"B");
    //NumbersTopL->Divide(NumbersTopL,DivideL,1.,1.,"B");
    //g->GetXaxis()->SetLimits(0.0,12.0);
   
    //setTDRStyle();
    TH1* Max=new TH1F("Max","pekeR",12.0,0.0,12.0);
    for (int bin = 0; bin <= NumbersTopL->GetNbinsX(); bin++) {
        if(NumbersTopR->GetBinContent(bin)>NumbersTopL->GetBinContent(bin)){
            Max->Fill(bin-1, NumbersTopR->GetBinContent(bin)); 
            Max->SetBinError(bin,NumbersTopR->GetBinError(bin));
        }else{
            Max->Fill(bin-1, NumbersTopL->GetBinContent(bin)); 
            Max->SetBinError(bin,NumbersTopL->GetBinError(bin));
        }
    }
    /*if(MODE==0){
        Max->SetLineColor(kRed);
    }else{
        Max->SetLineColor(kBlue);
    }
    Max->GetYaxis()->SetTitle("Tag and Probe Metric");
    Max->GetXaxis()->SetTitle("Bar Index");
    Max->SetAxisRange(0.0,1.0,"Y");
    Max->Draw("hist SAME e");*/
    //NumbersTopL->Draw("hist same e"); 
    //NumbersTopR->Draw("hist same e"); 
    //DoubleHits->Draw("hist e");
    //auto legend = new TLegend(0.1,0.1,.48,.3);
    //legend->AddEntry("Channel Efficiency","1 Electron Hit Efficiency");
    //legend->AddEntry("Max Efficiency","Max Efficiency");
    //legend->Draw("same");
    return Max;
}

void ChangedPlotterModOverlay(){
        //auto c = new TCanvas("c","c",600,600);
    auto h1 = ChangedPlotterMod("event3.root",0); 
    //FOR MONTE CARLO
    auto h5 = ChangedPlotterMod("Gap00_plots.root",1);
    auto h2 = ChangedPlotterGeom("event3.root");
    
    auto h3 = ChangedPlotterMod("May28.root",0);
    auto h4 = ChangedPlotterGeom("May28.root");

    h1->SetLineColor(kRed);
    h2->SetLineColor(kRed);

    h3->SetLineColor(kBlue);
    h4->SetLineColor(kBlue);

    h5->SetLineColor(kBlack);

    h1->GetYaxis()->SetTitle("Tag and Probe Metric");
    h1->GetXaxis()->SetTitle("Bar Index"); 
    setTDRStyle();
    h1->Draw("hist e");
    h2->Draw("hist e same");
    h3->Draw("hist e same");
    h4->Draw("hist e same");
    h5->Draw("hist e same");
    //c->Draw("SAME");
    gPad->SaveAs("shiftblocktop_3.png"); 
}

