using namespace std;

Double_t background(Double_t *x, Double_t *par) {
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t langaufixGsigmafun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 3000.0;      // number of convolution steps
      Double_t sc =  15.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;
      Double_t GSigma = TMath::Sqrt(par[1]);

      // MP shift correction
      mpc = par[1] - mpshift * par[0];
      //mpc=par[1];
      // Range of convolution integral
      //xlow = x[0] - sc * par[3];
      //xupp = x[0] + sc * par[3];
      xlow = x[0] - sc * GSigma;
      xupp = x[0] + sc * GSigma;

      step = (xupp-xlow) / np;

      //std::cout << "Input is " << x[0] << std::endl;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,GSigma);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,GSigma);
      }

      return (par[2] * step * sum * invsq2pi / GSigma);
}

Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 1000.0;      // number of convolution steps
      Double_t sc =  15.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0];

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      //std::cout << "Input is " << x[0] << std::endl;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}

Double_t FermiFunction(Double_t *x, Double_t *par) {
  return par[0]/(1+ TMath::Exp((x[0]-60.0)/par[1]));
}

Double_t StepFunction(Double_t *x, Double_t *par) {
  if (x[0] < par[0]) { 
    //std::cout << "Values check : " << x[0] << "::" << par[0] << std::endl;
    return par[1];
  }
  else {  
    //std::cout << "Values check : " << x[0] << "::" << par[0] << std::endl;
    return 0.0;
  }
}

Double_t LanGauPlusPoly(Double_t *x, Double_t *par) {
  //return langaufun(x,par) + background(x,&par[3]);
  double value1,value2,value;
  double langpar[3] = {par[0],par[1],par[2]};
  //double langpar[4] = {par[0],par[1],par[2],par[3]};
  double polypar[3] = {par[3],par[4],par[5]};
  //double polypar[3] = {par[4],-2.0*par[5]*par[1],par[5]};
  value1 = langaufixGsigmafun(x,langpar);
  //value2 = par[3] - 2.0*par[4]*par[1]*x[0] + par[4]*x[0]*x[0];
  value2 = background(x,polypar);
  value = value1 + value2;
  //value = par[6]*value1 + (1-par[6])*value2;
  return value;
}

Double_t LanGauPlusFermi(Double_t *x, Double_t *par) {
  //return langaufun(x,par) + FermiFunction(x,&par[2]);
  double value1,value2;
  double langpar[3] = {par[0],par[1],par[2]};
  double fermpar[2] = {par[3],par[4]};
  //double langpar[4] = {par[0],par[1],par[2],par[3]};
  //double fermpar[3] = {par[3],par[4],par[5]};
  value1 = langaufixGsigmafun(x,langpar);
  value2 = FermiFunction(x,fermpar);
  return value1+value2;
}

void setTDRstyle()
{
    // For the canvas:
  gStyle->SetCanvasBorderMode(1);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(1000); //Height of canvas
  gStyle->SetCanvasDefW(1400); //Width of canvas
  gStyle->SetCanvasDefX(10);   //POsition on screen
  gStyle->SetCanvasDefY(0);

  // For the Pad:
  gStyle->SetPadBorderMode(0);
  // gStyle->SetPadBorderSize(Width_t size = 1);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(2);
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
  //gStyle->SetErrorX(0.);
  //gStyle->SetErrorY(0.);
  
  gStyle->SetMarkerStyle(8);
  gStyle->SetMarkerSize(0.7);

  //For the fit/function:
  //gStyle->SetOptFit(0);
  gStyle->SetOptFit(1111);
  gStyle->SetFitFormat("5.4g");
  //gStyle->SetFuncColor(2);
  //gStyle->SetFuncStyle(1);
  //gStyle->SetFuncWidth(1);

  //For the date:
  gStyle->SetOptDate(0);
  // gStyle->SetDateX(Float_t x = 0.01);
  // gStyle->SetDateY(Float_t y = 0.01);

  //gStyle->SetOptStat("iou");

  
  // For the statistics box:
  //gStyle->SetOptFile(0);
  gStyle->SetOptStat(111111);
  //gStyle->SetStats(0);
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.03);///---> gStyle->SetStatFontSize(0.025);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("5.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatH(0.15);
  gStyle->SetStatW(0.15);///---> gStyle->SetStatW(0.15);
  

  // gStyle->SetStatStyle(Style_t style = 1001);
  // gStyle->SetStatX(Float_t x = 0);
  // gStyle->SetStatY(Float_t y = 0);

  // Margins:
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.04);

  // For the Global title:

  gStyle->SetOptTitle(1);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  //gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleX(0.26);
  gStyle->SetTitleH(0.07); // Set the height of the title box
  // gStyle->SetTitleW(0); // Set the width of the title box
  // gStyle->SetTitleX(0); // Set the position of the title box
  gStyle->SetTitleY(0.95); // Set the position of the title box
  //gStyle->SetTitleStyle(Style_t style = 1001);
  // gStyle->SetTitleBorderSize(2);

  // For the axis titles:

  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.08, "XYZ");
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
}

//void FitMerge(TFile *hzz, Int_t bins, float_t xmin, float_t xmax, TString &plot_title, TString &hsig1, TString &hsig2, TString &hsig3, TString &hsig4, TString &hsig5, TString &hsig6, TString &hsig7, TString &hsig8, TString &hsig9, TString &hsig10, TString &hsig11, TString &hsig12, TString &location)
//void FitMerge(TFile *hzz, Int_t bins, float_t xmin, float_t xmax, TString &plot_title, TString HSig[12], TString &location, TString &PlotIdentifier, std::vector<float> &LangauParam0, std::vector<float> &LangauParam1,std::vector<float> &LangauParam2, std::vector<float> &LangauParam3, std::vector<float> &LangauParam4, std::vector<float> &LangauParam5)
void FitMerge(TFile *hzz, Int_t bins, float_t xmin, float_t xmax, TString &plot_title, TString HSig[12], TString &location, TString &PlotIdentifier)
{
  TCanvas *c2 = new TCanvas();
  //TPaveStats *statsbox;

  setTDRstyle();

  //std::vector<TH1F*> SigList;
  TString HistNames[12] = {"Signal1","Signal2","Signal3","Signal4","Signal5","Signal6","Signal7","Signal8","Signal9","Signal10","Signal11","Signal12"};
  c2->Divide(3,4,0.005,0.01);
  //for(int si=0;si < 12; si++)
 //{
    //TH1F *sig = (TH1F*)gROOT->FindObject(HistNames[si]);
    //delete sig;
    //sig = new TH1F(HistNames[si],HSig[si],bins,xmin,xmax); //Chan8
    //sig = (TH1F*)hzz->Get(HSig[si]);
    //SigList.push_back(sig);
  //}
  
  const char *labels[12] = {"#bf{Channel 0}","#bf{Channel 1}","#bf{Channel 2}","#bf{Channel 3}","#bf{Channel 4}","#bf{Channel 5}","#bf{Channel 6}","#bf{Channel 7}","#bf{Channel 8}","#bf{Channel 9}","#bf{Channel 10}","#bf{Channel 11}"};
  for(int g=0;g<12;g++)
  {
    c2->cd(g+1);
    if (g != 8) {
    //std::cout << "Using histogram " << SigList.at(g) << " for plotting" << std::endl;
      TH1F *h = (TH1F*)gROOT->FindObject(HistNames[g]);
      delete h;
      TH1F *sig = new TH1F(HistNames[g],HSig[g],bins,xmin,xmax); //Chan8
      sig = (TH1F*)hzz->Get(HSig[g]);

      //TH1F *h = (TH1F*)gROOT->FindObject("Signal1");
      //delete h;
      //TString hsig = "hpe_ChanID0";
      //TH1F *sig = new TH1F("Signal1",hsig,bins,xmin,xmax); //Chan8
      //sig = (TH1F*)hzz->Get(hsig);
      sig->GetXaxis()->SetTitle("PE");
      sig->GetYaxis()->SetTitle("Number of events");
      TString title = labels[g];

      sig->SetLineWidth(2);
      sig->SetLineColor(kBlue);
      sig->SetTitle(title);
      sig->GetXaxis()->SetLabelSize(0.08);
      sig->GetXaxis()->SetTitleSize(0.08);
      sig->GetYaxis()->SetLabelSize(0.08);
      sig->GetYaxis()->SetTitleSize(0.08);
      sig->GetYaxis()->SetTitleOffset(0.85);
      sig->GetXaxis()->SetTitleOffset(0.75);
      sig->SetStats(0);

      sig->Draw();

      //TF1 *LF = new TF1("LF",LanGauPlusPoly,20,200,5);
      //LF->SetLineColor(kOrange);
      //LF->SetParameters(LangauParam0[g],LangauParam1[g],LangauParam2[g],LangauParam3[g],LangauParam4[g],LangauParam5[g]);
      //LF->Draw("same");
      if (g==2) {
        c2->Update();
    //gPad->Modified();
        TLatex *text = new TLatex(gPad->GetUxmax()-47,gPad->GetUymax(),"#it{LDMX Internal}");
    //text->SetNDC();
        text->SetTextSize(0.07);
        text->Draw();
      }
      gPad->Update();
    }
    else {
      TH1F *h = (TH1F*)gROOT->FindObject(HistNames[g]);
      delete h;
      TH1F *sig2 = (TH1F*)hzz->Get(HSig[0])->Clone();
      TF1 *fun = sig2->GetFunction("f_lan");
      //TF1 *fun2 = hIn->GetFunction("fGausBac");
      sig2->GetListOfFunctions()->Remove(fun);
      for(unsigned int ibin=0;ibin <= bins; ibin++) {
        sig2->SetBinContent(ibin, 0);
      }
      sig2->GetXaxis()->SetTitle("PE");
      sig2->GetYaxis()->SetTitle("Number of events");
      TString title = labels[g];

      sig2->SetLineWidth(2);
      sig2->SetLineColor(kBlue);
      sig2->SetTitle(title);
      sig2->GetXaxis()->SetLabelSize(0.08);
      sig2->GetXaxis()->SetTitleSize(0.08);
      sig2->GetYaxis()->SetLabelSize(0.08);
      sig2->GetYaxis()->SetTitleSize(0.08);
      sig2->GetYaxis()->SetTitleOffset(0.85);
      sig2->GetXaxis()->SetTitleOffset(0.75);
      sig2->SetStats(0);

      sig2->Draw();
    }
  }
  gPad->Update();
  c2->Update();
  //gPad->Modified();
  //TLatex *text = new TLatex(gPad->GetUxmax()-47,gPad->GetUymax(),"#it{LDMX Internal}");
  //text->SetNDC();
  //text->SetTextSize(0.05);
  //text->Draw();

  TString outs = location + "/" + PlotIdentifier + "_" + plot_title + "_langaupoly_fitfunction_combined_response_corrected.pdf";
  //TString outs = location + "/" + PlotIdentifier + "_" + plot_title + "_langaupoly_fitfunction_combined.pdf";
  c2->Print(outs);
}


void LDMX_SummarizeMipResponses(TString file_location, TString plot_location, TString Identifier, TString RunNumber)
{
  fstream file;
  Double_t parameters[6];
  int channelID;
  //TString sigfile = file_location + "/unpacked_reprocessed_ldmx_captan_out_" + Identifier +  "_digi_linearize_reprocessed_startsamples_pulsewidth_7_hits_gain_corrected_response_LangauPoly_FixGSigma.root";
  TString sigfile = file_location + "/unpacked_reprocessed_ldmx_captan_out_"+ Identifier +"_reprocessed_startsamples_pulsewidth_7_hits_gain_mip_corrected_response_LangauPoly_FixGSigma.root";
  //file.open(datafile,ios::in);

  TString plot_title_fit("hpe_chanID");
  TString plot_title_all("hPE_Mips_combined");

  TFile *hzz = TFile::Open(sigfile);

  //std::vector<float> LangauParams0;
  //std::vector<float> LangauParams1;
  //std::vector<float> LangauParams2;
  //std::vector<float> LangauParams3;
  //std::vector<float> LangauParams4;
  //std::vector<float> LangauParams5;

  //while(!file.eof()) {
    //file >> channelID >> parameters[0] >> parameters[1] >> parameters[2] >> parameters[3] >> parameters[4] >> parameters[5];
    //if (channelID != 17) { 
      //LangauParams0.push_back(parameters[0]);
      //LangauParams1.push_back(parameters[1]);
      //LangauParams2.push_back(parameters[2]);
      //LangauParams3.push_back(parameters[3]);
      //LangauParams4.push_back(parameters[4]);
      //LangauParams5.push_back(parameters[5]);
    //}
  //}

  TString plot_combined[12] = {"hpe_chanID0","hpe_chanID1","hpe_chanID2","hpe_chanID3","hpe_chanID4","hpe_chanID5","hpe_chanID6","hpe_chanID7","hpe_chanID8","hpe_chanID9","hpe_chanID10","hpe_chanID11"};
  FitMerge(hzz,81,20,200,plot_title_all,plot_combined,plot_location,RunNumber);

}