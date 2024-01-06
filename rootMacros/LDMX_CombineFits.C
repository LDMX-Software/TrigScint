//Root File to plot combination fit results

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
  gStyle->SetTitleFontSize(0.04);
  gStyle->SetTitleX(0.26);
  gStyle->SetTitleH(0.07); // Set the height of the title box
  // gStyle->SetTitleW(0); // Set the width of the title box
  // gStyle->SetTitleX(0); // Set the position of the title box
  gStyle->SetTitleY(0.95); // Set the position of the title box
  // gStyle->SetTitleStyle(Style_t style = 1001);
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

void PlotFitResiduals(TFile *hzz, Int_t bins, float_t xmin, float_t xmax, TString &plot_title, TString &hsig, Double_t parameters[6], Int_t lowerlimitfit, Int_t upperlimitfit, Int_t channelID, TString &location)
{
  setTDRstyle();
  gStyle->SetPadLeftMargin(0.1);
  TCanvas *c3 = new TCanvas();
  c3->Divide(1,2,0,0);
  TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.3,1.0,1.0);
  pad1->SetBottomMargin(0.01);
  TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.3);
  pad2->SetBottomMargin(0.2);
  //pad2->SetTopMargin(0);
  pad1->Draw();
  pad2->Draw();
  //TPaveStats *statsbox;

  setTDRstyle();
  //gStyle->SetOptFit(1111);

  TH1F *h = (TH1F*)gROOT->FindObject("Signal");
  delete h;
  TH1F *h2 = (TH1F*)gROOT->FindObject("SignalResiduals");
  delete h2;

  TH1F *sig = new TH1F("Signal",hsig,bins,xmin,xmax);
  sig = (TH1F*)hzz->Get(hsig);
  TH1F *sigres = new TH1F("SignalResiduals","",bins,xmin,xmax);
  TF1 *FitFunc = new TF1("SignalFit",LanGauPlusPoly,lowerlimitfit,upperlimitfit,6);
  FitFunc->SetParameters(parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5]);

  //sig->Sumw2();
  for(unsigned int ibin=1;ibin <= bins; ibin++) {
    float diff = sig->GetBinContent(ibin) - FitFunc->Eval(sig->GetBinCenter(ibin));
    float residual;
    if (sig->GetBinContent(ibin) == 0)
      residual = 0.0;
    else
      residual = diff/sig->GetBinError(ibin);
    std::cout << "Residuals data " << diff << "::" << sig->GetBinError(ibin) << "::" << sig->GetBinContent(ibin) << "::" << sig->GetBinCenter(ibin) << "::" << FitFunc->Eval(sig->GetBinCenter(ibin)) << std::endl;
    sigres->SetBinContent(ibin,residual);
    //sigres->SetBinError(ibin,0.0);
    sigres->SetBinError(ibin,1);
  }

  //c3->Divide(1,2);

  pad1->cd();
  //sig->GetXaxis()->SetTitle("PE");
  sig->GetYaxis()->SetTitle("Number of events");
  sig->GetYaxis()->SetTitleSize(0.06);
  sig->GetYaxis()->SetLabelSize(0.05);
  //sig->GetYaxis()->SetTitleSize(0.11);
  //sig->GetXaxis()->SetLabelSize(0.);
  sig->GetYaxis()->SetTitleOffset(0.65);
  sig->SetStats(0);
  TString title = "Channel " + to_string(channelID);
  if (channelID == 17) { title = "AllChannels" ;}
  //gPad->Modified();
  //gPad->Update();

  sig->SetLineWidth(2);
  sig->SetLineColor(kBlue);
  sig->SetTitle(title);
  //sig->SetTitleSize(0.06);
  //sig->SetStats(1);
  sig->Draw("E");
  c3->Update();
  TLatex *text = new TLatex(gPad->GetUxmax()-32,gPad->GetUymax(),"#it{LDMX Internal}");
  text->SetTextSize(0.05);
  text->Draw();

  //TString YTitle = "#frac{Measure #minus Fit}{MeasureUncertainty}";
  TString YTitle = "#frac{X_{measure} #minus X_{fit}}{#sigma_{measure}}";
  pad2->cd();
  sigres->GetXaxis()->SetTitle("PE");
  sigres->GetXaxis()->SetTitleOffset(0.75);
  sigres->GetYaxis()->SetTitle(YTitle);
  sigres->SetLineWidth(2);
  sigres->SetLineColor(kBlue);
  sigres->GetYaxis()->SetRangeUser(-4,4);
  sigres->GetXaxis()->SetLabelSize(0.13);
  sigres->GetXaxis()->SetTitleSize(0.13);
  sigres->GetYaxis()->SetLabelSize(0.12);
  sigres->GetYaxis()->SetTitleSize(0.13);
  sigres->GetYaxis()->SetTitleOffset(0.30);
  sigres->SetStats(0);

  sigres->Draw("E");  

  TString outs = location + plot_title + "_langaupoly_" + channelID + "_fitfunctionresiduals.pdf";
  c3->Print(outs); 
}

void FitPrint(TFile *hzz, Int_t bins, float_t xmin, float_t xmax, TString &plot_title, TString &hsig, Int_t lowerlimitfit, Int_t upperlimitfit, Int_t channelID, TString &location)
{
  TCanvas *c2 = new TCanvas();
  //TPaveStats *statsbox;

  setTDRstyle();

  TH1F *h = (TH1F*)gROOT->FindObject("Signal");
  delete h;
  TH1F *sig = new TH1F("Signal",hsig,bins,xmin,xmax);
  TF1 *f = sig->GetFunction("f_lan");
  sig->GetListOfFunctions()->Remove(f);
  sig = (TH1F*)hzz->Get(hsig);
  sig->GetXaxis()->SetTitle("PE");
  sig->GetYaxis()->SetTitle("Number of events");
  TString title = "Channel " + to_string(channelID);
  if (channelID == 17) { title = "AllChannels" ;}

  sig->SetLineWidth(2);
  sig->SetLineColor(kBlue);
  sig->SetTitle(title);
;
  sig->Draw();
  c2->Update();
  TLatex *text = new TLatex(gPad->GetUxmax()-47,gPad->GetUymax(),"#it{LDMX Internal}");
  text->SetTextSize(0.05);
  text->Draw();

  TString outs = location + plot_title + "_langaupoly_" + channelID + "_fitfunction.pdf";
  c2->Print(outs);

}

void FitSeperator(TFile *hzz, Int_t bins, float_t xmin, float_t xmax, TString &plot_title, TString &hsig, Double_t parameters[6], Int_t lowerlimitfit, Int_t upperlimitfit, Int_t channelID, TString &location, Int_t drawlegend)
{
	TCanvas *c1 = new TCanvas();

  setTDRstyle();

  TH1F *h = (TH1F*)gROOT->FindObject("Signal");
  delete h;
	TH1F *sig = new TH1F("Signal",hsig,bins,xmin,xmax);
	sig = (TH1F*)hzz->Get(hsig);
  sig->GetXaxis()->SetTitle("PE");
  sig->GetYaxis()->SetTitle("Number of events");
  TString title = "Channel " + to_string(channelID);
  if (channelID == 17) { title = "AllChannels" ;}
  //gPad->Modified();
  //gPad->Update();

	sig->SetLineWidth(2);
	sig->SetLineColor(kBlue);
	sig->SetTitle(title);
	sig->SetStats(0);

	TF1 *LGF = new TF1("LGF",langaufixGsigmafun,lowerlimitfit,upperlimitfit,3);
	LGF->SetLineColor(kMagenta);
	LGF->SetParameters(parameters[0],parameters[1],parameters[2]);

	TF1 *Poly = new TF1("Poly",background,lowerlimitfit,upperlimitfit,3);
	Poly->SetLineColor(kBlack);
	Poly->SetParameters(parameters[3],parameters[4],parameters[5]);

  TF1 *Comb = new TF1("Comb",LanGauPlusPoly,lowerlimitfit,upperlimitfit,6);
  Comb->SetLineColor(kRed);
  Comb->SetParameters(parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5]);

  TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
  legend->AddEntry(sig,"Data");
  legend->AddEntry(Comb,"Signal + Background");
  legend->AddEntry(LGF,"Langau Signal");
  //legend->AddEntry(Step,"Fermi Background");
  legend->AddEntry(Poly,"Polynomial Background");
  //legend->AddEntry(LGF_mod,"Langau contribution");
  //legend->AddEntry(Poly_mod,"Background contribution");
  if (drawlegend == 0) { 
    sig->Draw();
    c1->Update();
    std::cout << channelID << "::" << gPad->GetUymax() << "::" << gPad->GetUxmax() << std::endl;
    TLatex *text = new TLatex(gPad->GetUxmax()-47,gPad->GetUymax()+2902,"#it{LDMX Internal}");
    text->SetTextSize(0.05);
    text->Draw();
    LGF->Draw("same");
    Poly->Draw("same");
    Comb->Draw("same");
    c1->SetLogy();
    legend->Draw();
    //TLatex text;
    //text.SetTextSize(0.05);
    //text.DrawLatex(gPad->GetUxmax(),gPad->GetUymax(),"#it{LDMX Internal}");
    std::cout << "After drawing legend =" << channelID << "::" << gPad->GetUymax() << "::" << gPad->GetUxmax() << std::endl;
    TString outs = location + plot_title + "_langaupoly_" + channelID + "_log.pdf";
    c1->Print(outs);
  }
  else {
    sig->Draw();
    c1->Update();
    LGF->Draw("same");
    Poly->Draw("same");
    Comb->Draw("same");
    c1->SetLogy();
    std::cout << channelID << "::" << gPad->GetUymax() << "::" << gPad->GetUxmax() << std::endl;
    TLatex *text = new TLatex(gPad->GetUxmax()-47,gPad->GetUymax(),"#it{LDMX Internal}");
    text->SetTextSize(0.05);
    text->Draw();
    std::cout << "After drawing legend =" << channelID << "::" << gPad->GetUymax() << "::" << gPad->GetUxmax() << std::endl;
    TString outs = location + plot_title + "_langaupoly_" + channelID + "_log.pdf";
    c1->Print(outs);
  }
}

//void LDMX_CombineFits(const char *sigfile, const char *datafile, const char *identifier_functions, const char *identifier_fits, const char *identifier_residuals)
void LDMX_CombineFits(TString location, TString plots_location, TString Identifier)
{
	fstream file;
	Double_t parameters[6];
	int channelID;
  //TString sigfile = location + "/unpacked_reprocessed_ldmx_captan_out_" + Identifier +  "_digi_linearize_reprocessed_startsamples_pulsewidth_7_hits_gain_corrected_response_LangauPoly_FixGSigma.root";
	//TString datafile = location + "/unpacked_reprocessed_ldmx_captan_out_" + Identifier +  "_digi_linearize_reprocessed_startsamples_pulsewidth_7_hits_gain_corrected_response_LangauPolyParams_FixGSigma.txt";
  TString sigfile = location + "/FitRootFiles/unpacked_reprocessed_ldmx_captan_out_" + Identifier +  "_reprocessed_startsamples_pulsewidth_7_hits_gain_mip_corrected_response_LangauPoly_FixGSigma.root";
  TString datafile = location + "/FitParameters/unpacked_reprocessed_ldmx_captan_out_" + Identifier +  "_reprocessed_startsamples_pulsewidth_7_hits_gain_mip_corrected_response_LangauPolyParams_FixGSigma.txt";
  file.open(datafile,ios::in);
  //TString loc = "/sdf/group/ldmx/users/dhparmar/LDMX_RAWData/validation_plots/FitPlots/LangauPolyPlots/"
  //TString loc = "/home/dhruvanshu/LDMX_Analysis_Files/FitTests/FocusTalkRuns/MipResponseExampleFits/";
  TString location_fits = plots_location + "/Fits_ResponseCorrected/";
  TString location_functions = plots_location + "/FitFunctions_ResponseCorrected/";
  TString location_residuals = plots_location + "/FitResiduals_ResponseCorrected/";
  gSystem->mkdir(location_residuals);
  gSystem->mkdir(location_functions);
  gSystem->mkdir(location_fits);
	TString plot_title_fit("hpe_chanID");
  TString plot_title_all("hPE_all");

	TFile *hzz = TFile::Open(sigfile);

	while(!file.eof()) {
		file >> channelID >> parameters[0] >> parameters[1] >> parameters[2] >> parameters[3] >> parameters[4] >> parameters[5]; 
    if (file.eof()) break;
    if(channelID == 17) { 
      TString plot_channel = plot_title_all;
      FitPrint(hzz,81,20,200,plot_channel,plot_channel,20,200,channelID,location_fits);
      FitPrint(hzz,81,20,200,plot_channel,plot_channel,20,200,channelID,location_fits);
      FitSeperator(hzz,81,20,200,plot_channel,plot_channel,parameters,20,200,channelID,location_functions,0);
      FitSeperator(hzz,81,20,200,plot_channel,plot_channel,parameters,20,200,channelID,location_functions,0);
      PlotFitResiduals(hzz,81,20,200,plot_channel,plot_channel,parameters,20,200,channelID,location_residuals);
    }
    else { 
      std::cout << "Plotting histograms for channel " << channelID << std::endl;
      TString plot_channel = plot_title_fit + channelID;
      FitPrint(hzz,81,20,200,plot_channel,plot_channel,20,200,channelID,location_fits);
      FitSeperator(hzz,81,20,200,plot_channel,plot_channel,parameters,20,200,channelID,location_functions,1);
      PlotFitResiduals(hzz,81,20,200,plot_channel,plot_channel,parameters,20,200,channelID,location_residuals);
      PlotFitResiduals(hzz,81,20,200,plot_channel,plot_channel,parameters,20,200,channelID,location_residuals);
    }
		//FitSeperator(hzz,45,20,200,plot_channel,plot_channel,parameters,20,200,channelID);	
	}
  file.close();

}