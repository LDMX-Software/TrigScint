//Mild modifications and hacks done by Dhruvanshu Parmar to make the script running
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TF1NormSum.h>
#include <TGraphErrors.h>


#include <string>
#include <vector>
#include <iostream>
#include <fstream>


//LinearFit

Double_t Linear(Double_t *x, Double_t *par) {
  return par[0]*x[0] + par[1];
}

// Quadratic background function
Double_t background(Double_t *x, Double_t *par) {
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

// Lorentzian Peak function
Double_t lorentzianPeak(Double_t *x, Double_t *par) {
  return (0.5*par[0]*par[1]/TMath::Pi()) / TMath::Max(1.e-10,
													  (x[0]-par[2])*(x[0]-par[2])+ .25*par[1]*par[1]);
}


//some function with a broader low side
Double_t asymPeak(Double_t *x, Double_t *par) {
  Double_t expo=TMath::Exp((x[0]-par[1])/par[0]);
  par[2]=sqrt(fabs(par[1]));
  return 1/par[0]*expo*TMath::Exp( -expo );
  //f(x)=1βe−x−μβe−e−x−μβ
}

// Gaussian Peak function                                                                                                                                    
Double_t gaussianPeak(Double_t *x, Double_t *par) {
  return par[0]*TMath::Exp(  -0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2])  );
}

Double_t landaufun(Double_t *x, Double_t *par) {

  Double_t mpshift  = -0.22278298;
  Double_t mpc;

  mpc = par[1] - mpshift * par[0];

  return TMath::Landau(x[0],mpc,par[0]) / par[0];

}


// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
  //  return background(x,par) + lorentzianPeak(x,&par[3]);
  return gaussianPeak(x,par) + Linear(x,&par[2]);
  //return gaussianPeak(x,par) + background(x,&par[3]);
}

Double_t StepFunction(Double_t *x, Double_t *par) {
  if (x[0] < 40.0) { 
    //std::cout << "Values check : " << x[0] << "::" << par[0] << std::endl;
    return par[0];
  }
  else {  
    //std::cout << "Values check : " << x[0] << "::" << par[0] << std::endl;
    return 0.0;
  }
}

Double_t FermiFunction(Double_t *x, Double_t *par) {
  //return par[0]/(1+ TMath::Exp((x[0]-30.0)/par[1]));
  //return 20.0/(1+ TMath::Exp((x[0]-par[0])/par[1]));
  return par[0]/(1+ TMath::Exp((x[0]-par[1])/4.0));
  //return par[0]/(1+ TMath::Exp((x[0]-par[1])/par[2]));
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
      //Double_t GSigma = TMath::Sqrt(par[1]);

      // MP shift correction
      mpc = par[1] - mpshift * par[0];
      //mpc=par[1];
      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];
      //xlow = x[0] - sc * GSigma;
      //xupp = x[0] + sc * GSigma;

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


Double_t LanGauPlusPoly(Double_t *x, Double_t *par) {
  //return langaufun(x,par) + background(x,&par[3]);
  double value1,value2,value;
  double langpar[3] = {par[0],par[1],par[2]};
  //double langpar[4] = {par[0],par[1],par[2],par[3]};
  double polypar[3] = {par[3],par[4],par[5]};
  //double polypar[3] = {par[3],-2.0*par[4]*par[1],par[4]}; //For ax^2+bx+c, setting b = -2*a*x_0
  value1 = langaufixGsigmafun(x,langpar);
  //value2 = par[3] - 2.0*par[4]*par[1]*x[0] + par[4]*x[0]*x[0];
  value2 = background(x,polypar);
  value = value1 + value2;
  return value;

}

Double_t LanGauCombine(Double_t *x, Double_t *par) {
  //return langaufun(x,par) + StepFunction(x,&par[1]);
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
      Double_t value;
      Double_t GSigma = TMath::Sqrt(par[1]);


      // MP shift correction
      mpc = par[1] - mpshift * par[0];

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
         //sum += fland * TMath::Gaus(x[0],xx,par[3]);
         sum += fland * TMath::Gaus(x[0],xx,GSigma);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         //sum += fland * TMath::Gaus(x[0],xx,par[3]);
         sum += fland * TMath::Gaus(x[0],xx,GSigma);
      }

      //value = (par[2] * step * sum * invsq2pi / par[3]) + (par[4]/(1+ TMath::Exp((x[0]-60.0)/par[5]))); //Uncomment to include Fermi function
      value = (par[2] * step * sum * invsq2pi / GSigma) + (par[3]/(1+ TMath::Exp((x[0]-60.0)/par[4]))); //Uncomment to include Fermi function
      //value = (par[2] * step * sum * invsq2pi / par[3]) + par[4] + par[5]*x[0] + par[6]*x[0]*x[0]; //Uncomment to include Poly function
      //if (x[0] < par[4] || x[0] == par[4]) {  //Uncomment to include step function
        //value = (par[2] * step * sum * invsq2pi / par[3]) + par[5]; 
      //}
      //if (x[0] > par[4]) {
        //value = (par[2] * step * sum * invsq2pi / par[3]);
      //}
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

Double_t LanDauPlusPoly(Double_t *x, Double_t *par) {
  return landaufun(x,par) + background(x,&par[3]);
}

std::vector<double> MeasureEfficiency(std::vector<double> &parameters) //Function calculating efficiency by ratio of integral from MIP threshold and integral from 20 PE peak.
{
  TF1 *LGF = new TF1("LGF",langaufixGsigmafun,20,200,3);
  //int Threshold[3] = {40,50,60};
  int Threshold = 50;
  //double MPVI = parameters[1];
  double MPVR[3] = {parameters[1]-10.0,parameters[1],parameters[1]+10.0};
  double IntTotal=0.0;
  double IntThr=0.0;
  double Eff;
  std::vector<double> EfficiencyVector;
  //LGF->SetParameters(parameters[0],parameters[1],parameters[2],parameters[3]);
  //for(unsigned int i=0;i < 3; i++) {
  for(unsigned int j=0;j<3;j++) {
    LGF->SetParameters(parameters[0],MPVR[j],parameters[2]);
    IntTotal=LGF->Integral(20,200);
    IntThr=LGF->Integral(Threshold,200);
    Eff = IntThr/IntTotal;
    EfficiencyVector.push_back(Eff);
    std::cout << Threshold << "::" << MPVR[j] << "::" << IntThr << "::" << IntTotal << "::" << Eff << std::endl;
  }
  //}
  return EfficiencyVector;
}

std::pair<float,float> FindPeak(std::vector<double> &parameters, double mpvuncer, double widthuncer, int index, double covmw)
{
  std::vector<double> TestPeaks;
  std::vector<double> Locations;
  double PeakMaxLoc;
  double MPV_factor;
  //Double_t parameters_langau[3],parameters_langaupoly[6],parameters_langaufermi[5];
  //if (index == 3) {Double_t parameters_langau[3] = {parameters[0],parameters[1],parameters[2]};}
  //if (index == 4) {Double_t parameters_langaupoly[6] = {parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5]};}
  //if (index == 5) {Double_t parameters_langaufermi[5] = {parameters[0],parameters[1],parameters[2],parameters[3],parameters[4]};}
  Double_t StartPeak = parameters[1] - parameters[0];
  Double_t EndPeak = parameters[1] + parameters[0];
  Double_t TestPeak;
  while (StartPeak < EndPeak) {
    //TestPeak = LanGauPlusPoly(&StartPeak,parameters_list);
    if (index == 3) {
      Double_t parameters_langau[3] = {parameters[0],parameters[1],parameters[2]};
      TestPeak = langaufun(&StartPeak,parameters_langau);
      //TestPeaks.push_back(TestPeak);
      //Locations.push_back(StartPeak);
    }
    if (index == 4) { 
      Double_t parameters_langaupoly[6] = {parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5]};
      TestPeak = LanGauPlusPoly(&StartPeak,parameters_langaupoly);
      //std::cout << "TestPeak is" << TestPeak << " using mpv " << parameters_langaupoly[0] << std::endl;
      //TestPeaks.push_back(TestPeak);
      //Locations.push_back(StartPeak);
    }
    if (index == 5) { 
      Double_t parameters_langaufermi[5] = {parameters[0],parameters[1],parameters[2],parameters[3],parameters[4]};
      TestPeak = LanGauPlusFermi(&StartPeak,parameters_langaufermi);
      //TestPeaks.push_back(TestPeak);
      //Locations.push_back(StartPeak);
    }
    //std::cout << TestPeak << "==" << index << std::endl;
    TestPeaks.push_back(TestPeak);
    Locations.push_back(StartPeak);    
    StartPeak+=0.001;
  }
  double PeakMax = *std::max_element(std::begin(TestPeaks), std::end(TestPeaks));
  for(unsigned int ser=0;ser < TestPeaks.size(); ser++) {
   if (TestPeaks[ser] == PeakMax) { 
    std::cout << "Peak found at " << Locations[ser] << " with value " << TestPeaks[ser] << std::endl;;
    PeakMaxLoc = Locations[ser];
  } 
  }
  MPV_factor = (PeakMaxLoc - parameters[1])/parameters[0];
  double MPV_factor_uncer = TMath::Sqrt(TMath::Power(mpvuncer,2) + TMath::Power(MPV_factor*widthuncer,2) + 2*MPV_factor*TMath::Power(covmw,2));
  std::cout << mpvuncer << "::" << widthuncer << "::" << MPV_factor << std::endl;
  return make_pair(PeakMaxLoc,MPV_factor_uncer);
}

//std::pair<int,int> MinimizeChi(TH1F *hsig)
//Changed passname to "hits" from "TBReco" coz the hits clusters do not posses that passname in collection
void extractResponse_final(string inFile, int deadChannel=-1, string passName = "hits", int nSamples = 30, int verbosity=2, bool isSim=false, bool doClean=true) 
{
  // macro extracting the per-channel MIP response based on calibrated hits
  //

  // we need hits, that are cleaned and calibrated
  // the average PE response, fit it --> this is what we calibrate to
  // then for each channel, fit and get mean
  // store the per-channel scale factor so we can apply it in the MIP regime 
  // done
  

  bool verbosePrint = (verbosity == 0 ); //adhere to ldmx printout lavel numbering
  
  //gStyle->SetOptFit(1111);

  TFile * fIn = TFile::Open( inFile.c_str(), "READ");
  TTree * tree = (TTree*)fIn->Get("LDMX_Events");
  fIn->ls();
  int exampleEvNb = 1;

  string digiName="testBeamHitsUp";
  int nChannels = 12;
  
  vector<string> vars = {"pe"}; // this variable is totQ/nSamp
  vector <int> maxVals = {200};  // don't need to go high, focus on the single-few PE peaks 
  vector <int> minVals = {20};  // should cover most negative pedestals 
  vector <float> binFactor = {0.45}; // histogram binning: bins per integer 

  TCanvas * c1 = new TCanvas("c1", "plot canvas", 600, 500);

  vector<string> cuts;
  vector<TH1F*> v_hOut;

  string cleanStr="";
  if (doClean)
  {
    cleanStr=Form("(%s_%s.flag_ == 0 || %s_%s.flag_ == 4) &&", digiName.c_str(),passName.c_str(), digiName.c_str(),passName.c_str());
    //cleanStr=Form("(%s_%s.flag_ == 14 ) &&", digiName.c_str(),passName.c_str(), digiName.c_str(),passName.c_str());
    std::cout << "Clean" << std::endl; 
  }
	
  
  for (int iC = 0 ; iC < nChannels ; iC++)
	cuts.push_back( Form("%s %s_%s.barID_ == %i", cleanStr.c_str(), digiName.c_str(),passName.c_str(), iC));

  c1->SetRightMargin( 1.5*c1->GetRightMargin());

  for (unsigned int iV = 0; iV < vars.size(); iV++) {
	for (unsigned int iC = 0; iC < cuts.size(); iC++) {
	  string cut =cuts[iC];
	  //set up binning
	  string bins = Form("%i, %i, %i", (int)(binFactor[iV]*(maxVals[iV]-minVals[iV])), minVals[iV], maxVals[iV]);
	  if (iC == 0)
		std::cout << bins << std::endl;
	  //draw from the tree
	  string bin = Form("%i", iC);
	  string nSamp = Form("%i", nSamples);
	  std::cout<<"At channel " << bin << std::endl;
	  std::cout<<"Using cut " << cut << std::endl;
	  tree->Draw( (digiName+"_"+passName+"."+vars[iV]+"_ >> h"+bin+"("+bins+")").c_str(), cut.c_str() );
	  //get them each and keep for later
	  TH1F *hOut = (TH1F*)gDirectory->Get(Form("h%i", iC)); 
	  if (!hOut || hOut->IsZombie()) {
		std::cout << "No histogram for channel " << iC << ", skipping" << std::endl;
		continue;
	  }
	  //set up axes 	  
	  string title = ";"+vars[iV]+" "+cut+";entries";
	  hOut->SetTitle( title.c_str() );
	  hOut->GetYaxis()->SetRangeUser(0.001, 1.5*hOut->GetMaximum()) ;
	  //draw 
	  hOut->Draw();
	  hOut->SetName(Form("h%s_chanID%i", vars[iV].c_str(), iC));
	  v_hOut.push_back( (TH1F*)hOut->Clone() );
	}//over cuts
  }//over variables 


  //make a histogram with all channels

  TH1F * hAll = (TH1F*)v_hOut[8]->Clone("hPE_all");
  for (unsigned int iC = 0; iC < nChannels; iC++) { //Changed initial iC to '0' from '1' coz I think it will skip channel 0 in earlier case.
	if (iC != 8) {//the if block for some reason is not recognizing deadChannel, hence the script throws an exception while processing channel 8. 
	  //continue;//A tiny hack but replacing 'deadChannel' by '8', i.e specifically telling script to skip channel 8 works. 
	hAll->Add(v_hOut[iC]);
    }
  }
  hAll->Add(v_hOut[8],-1); //I noticed that the earlier version of the script used to clone histogram of channel 0, which in the loop used to get repeated hence got counted twice in avgResponse. I tried
//resolving it by cloning information from the dead channel because it was easier to delete chan 8 histogram in the end.
  int upperlimitfit = 200;
  int lowerlimitfit = 20;
  fstream Response;// responseFile.Data());
  TString responseFile = inFile;
  responseFile = responseFile.ReplaceAll(".root", "");
  Response.open(responseFile + "_response_LangauPolyParams_FixGSigma.txt",ios::out); //Modification done by Dhruvanshu. Apparently, just passing the stream as input was not creating a file due to good()
  //function not counting it as an acceptable stream. Now this script is able to output the data as a text file.
  if (!Response.good()) {
    std::cout << "File not ready to use" << std::endl;
  }
  TF1 *Nx = new TF1("Nx",LanGauPlusPoly,lowerlimitfit,upperlimitfit,6);
  Nx->SetParNames("Width","MP","Area","Const","CoeffX","CoeffX2"); //Parameters name are self explanatory. Look into the functions to get a detailed idea about the parameters of langau function.
  Nx->SetParameters(1.8,80.0,10000.0,0.5);			  //Uncomment while fitting a polynomial background function.	
  TFitResultPtr fPtr = hAll->Fit("Nx", "QSR");
  float avgResponse = fPtr->Parameter(1);
  float avgResponseError = fPtr->ParError(1);
  std::cout << "Got MPV for all channels combined: " << avgResponse << std::endl;
  Response << 17 << " " << fPtr->Parameter(0) << " " << fPtr->Parameter(1) << " " << fPtr->Parameter(2) << " " << fPtr->Parameter(3) << " " << fPtr->Parameter(4) << " " <<  fPtr->Parameter(5) << "\n";
  
  fstream EffFile;
  EffFile.open(responseFile + "_LangauPoly_FixGSigma_EffValues.txt",ios::out);
  std::vector<float> ChannelChiSquare; //Stores the chisquare/NDF value
  std::vector<float> Integrals; //stores the number of MIPS
  int ChannelFitIndex[12] = {3,3,3,3,3,3,3,3,0,3,3,3}; //This is here because initially I was performing a chisquare optimization for the range. We noticed that it affected the tail fitting, so
  //we are now fixed to one fit range. FitIndex here corresponds to the type of function to fit as can be seen in the following loop.
  //float fitRangeWidth=200.; // take this window to either side of mean to catch peak 
  TGraphErrors* gMPVs = new TGraphErrors(nChannels);
  for (unsigned int iH = 0;  iH < v_hOut.size(); iH++) {
	if (iH == 8 or v_hOut.at(iH)->Integral() == 0) { //Seems like deadChannel is not getting recognized, so I replaced the deadChannel variable by specifically '8' coz we know channel 8 is dead
	  ChannelChiSquare.push_back(0.0);
    Integrals.push_back(0);
    std::cout << "Checking for channel " << iH << std::endl;
    continue;
  }
	//do the actual peak fitting
  else {
    std::cout << "----> Fitting response for channel " << iH << std::endl;
    std::cout << "Checking for channel " << iH << " with number of elements " << v_hOut.at(iH)->Integral() << std::endl;
    TF1 *f_lan = new TF1("f_lan",LanGauPlusPoly,lowerlimitfit,upperlimitfit,6); //Uncomment when using LangauPlusPoly
    f_lan->SetParameters(1.8,80.0,10000.0,0.5);
    f_lan->SetParName(0,"Width");
    f_lan->SetParName(1,"MP");
    f_lan->SetParName(2,"Area");
    f_lan->SetParName(3,"Constant");
    f_lan->SetParName(4,"CoeffX");
    f_lan->SetParName(5,"CoeffX2");
    
    TFitResultPtr fResult = v_hOut.at(iH)->Fit("f_lan","QSR");
    std::cout << "For channel " << iH << ", got MPV \t" << fResult->Parameter(1) << std::endl;
    TMatrixD cov = fResult->GetCovarianceMatrix();
    cov.Print(); 
    std::cout << "Covariance =" << cov[0][1] << std::endl; 
    
    std::vector<double> FitParams = {fResult->Parameter(0),fResult->Parameter(1),fResult->Parameter(2)}; //Only signal function considered. Remember, we are setting GSigma = sqrt(MPV)
    auto EV = MeasureEfficiency(FitParams);

    gMPVs->SetPoint( iH, iH, fResult->Parameter(1));
    gMPVs->SetPointError( iH, 0, fResult->ParError(1));
    ChannelChiSquare.push_back(f_lan->GetChisquare()/f_lan->GetNDF());
    Integrals.push_back(v_hOut.at(iH)->Integral());
    
    Response << iH << " " << fResult->Parameter(0) << " " << fResult->Parameter(1) << " " << fResult->Parameter(2) << " " << fResult->Parameter(3) << " " << fResult->Parameter(4) << " " << fResult->Parameter(5) << "\n";
    EffFile << iH << "," << EV[0] << "," << EV[1] << "," << EV[2] << "\n";
  }
  }
  Response.close();
  EffFile.close();
  
  //store this to a root file for later plotting
  TString outFile = inFile;
  outFile=outFile.ReplaceAll(".root", "_response_LangauPoly_FixGSigma.root"); //Root file storing all the fits
  std::cout << "using outfile name" "" << outFile << std::endl;

  TFile * fOut = TFile::Open( outFile.Data(), "RECREATE");
  fOut->cd();
  fOut->ls();
  std::cout << "Got " << v_hOut.size() << " histograms" << std::endl;
  hAll->Write();
  int nChanToWrite = deadChannel > -1 ? nChannels - 1 : nChannels; //adjust for if we have a dead channel in this run
   for (unsigned int iH = 0;  iH < v_hOut.size(); iH++) {
	v_hOut[iH]->Write();
  }
  if (deadChannel > -1 ) { //if there in fact is one (TODO: allow for list?)
	gMPVs->RemovePoint(deadChannel);
  }
  gMPVs->SetTitle(";channel ID;MIP peak fit MPV");
  gMPVs->Write("gMIPResponseVsChanID");
  fOut->Close();
  fIn->Close();
  
  //write responses to a .txt file for easy extraction of calibration in the hit reconstruction step 
  //std::cout << "using outfile name" "" << responseFile << std::endl;
  fstream ResponseFin;// responseFile.Data());
  TString responseFinFile = inFile;
  responseFinFile = responseFinFile.ReplaceAll(".root", ""); //Store the response PEs with some other useful data.
  //TString responseFile2 = responseFile + "_response.txt";
  //std::cout << responseFile2 << std::endl;
  ResponseFin.open(responseFinFile + "_response_LangauPolyMips_FixGSigma.txt",ios::out); //Modification done by Dhruvanshu. Apparently, just passing the stream as input was not creating a file due to good()
  //function not counting it as an acceptable stream. Now this script is able to output the data as a text file.
  if (!ResponseFin.good()) {
    std::cout << "File not ready to use" << std::endl;
  }
  for (int iP = 0; iP < gMPVs->GetN(); iP++) {
	double x, MPV, MPVerror;
  if (iP == 8) {
    ResponseFin << 8 << "," << 0 << "," << 0 <<  "," << 0 << "," << 0 << "," << 0 << "," << 0 << "," << 0 << "\n";
  }
  else {
    gMPVs->GetPoint(iP, x, MPV);
    MPVerror = gMPVs->GetErrorY(x);
    //ResponseFin << x << "," << avgResponse/MPV  << "," << Integrals[iP] << "," << MPV << "," << ChannelChiSquare[iP] << "," << ChannelFitIndex[iP] << "," << ChannelFitLowLimit[iP] << "\n";  
    ResponseFin << x << "," << avgResponse << "," << avgResponseError << "," << MPV << "," << MPVerror << "," << Integrals[iP] << "," << ChannelChiSquare[iP] << "," << avgResponse/MPV << "\n";
    std::cout << x << "," << avgResponse/MPV << std::endl;
  }
  }
  std::cout << "Data written" << std::endl;
  ResponseFin.close();
    //Done.
}

