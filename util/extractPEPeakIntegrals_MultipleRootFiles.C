#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TChain.h>

#include <string>
#include <vector>
#include <iostream>
#include <fstream>


//TGraphErrors * isolateAndFitPeaks( TH1F * hIn, float width, bool verbose, bool isSim, int nSamples);
TGraphErrors * findAndFitPeaks( TH1F * hIn, float width, bool verbose, bool isSim);


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


// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
  //  return background(x,par) + lorentzianPeak(x,&par[3]);
  return gaussianPeak(x,par) + background(x,&par[3]);
}

Double_t ExponentialGaus(Double_t *x, Double_t *par) {
	double value_gaus, value_expo;
	value_gaus = par[0]*TMath::Exp(  -0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2])  );
	value_expo = par[3]*TMath::Exp( -1.0*x[0]*par[4]);
	return value_gaus + value_expo;
}

Double_t ExponentialBackground(Double_t *x, Double_t *par) {
	return par[0]*TMath::Exp( -1.0*par[1]*x[0]);
}

void DrawGaussians(TH1F *hIn, TString &plottitle, std::vector<float> &lowlimits, std::vector<float> &upperlimits, std::vector<float> &GausConst, std::vector<float> &GausMean, std::vector<float> &GausSigma, std::vector<float> &ExpConst, std::vector<float> &ExpScale, TString &plot, TString &location)
{
	TCanvas *c1 = new TCanvas();

	int colorlist[9] = {1,2,3,5,6,7,8,9,10};
	hIn->SetLineWidth(2);
	hIn->SetLineColor(kBlue);
	hIn->SetTitle(plottitle);
	//sig->SetStats(111111);
	TF1 *fun = hIn->GetFunction("fBkg");
	TF1 *fun2 = hIn->GetFunction("fGausBac");
  hIn->GetListOfFunctions()->Remove(fun);
  hIn->GetListOfFunctions()->Remove(fun2);
	hIn->Draw();
	//c1->SetLogy();
 	for(int i=0;i < lowlimits.size(); i++)
	{
		std::cout << "FitParameters for drawing " << "===" << GausConst[i] << "===" << GausMean[i] << "===" << GausSigma[i] << "===" << ExpConst[i] << "===" << ExpScale[i] << std::endl;
		//TF1 *GPF = new TF1("GPF",gaussianPeak,lowlimits[i],upperlimits[i],3);
		TF1 *GPF = new TF1("GPF",ExponentialGaus,lowlimits[i],upperlimits[i],5);
		TF1 *Gau = new TF1("Gau",gaussianPeak,lowlimits[i],upperlimits[i],3);
		GPF->SetLineColor(colorlist[i]);
		GPF->SetParameters(GausConst[i],GausMean[i],GausSigma[i],ExpConst[i],ExpScale[i]);
		Gau->SetLineColor(7);
		Gau->SetParameters(GausConst[i],GausMean[i],GausSigma[i]);
		GPF->Draw("same");
		Gau->Draw("same");
	}
	//TF1 *GausBac = new TF1("GausBac",ExponentialBackground,0,5,2);
	//GausBac->SetParameters(ExpConst[0],ExpScale[0]);
	//GausBac->SetLineColor(6);
	//GausBac->Draw("same");
	TString outs = location + "/" + plot + ".png";
	TString outs2 = location + "/" + plot + "_log.png";
	c1->Print(outs);
	c1->SetLogy();
	c1->Print(outs2);
	//hIn->Clear();
	//return PeakData;
}

TMatrixD ConstructCovarianceMatrix(TMatrixD FitCovMat)
{
  TMatrix GausCovMatrix(3,3);
  GausCovMatrix(0,0) = FitCovMat(0,0);
  GausCovMatrix(0,1) = FitCovMat(0,1);
  GausCovMatrix(0,2) = FitCovMat(0,2);
  GausCovMatrix(1,0) = FitCovMat(1,0);
  GausCovMatrix(1,1) = FitCovMat(1,1);
  GausCovMatrix(1,2) = FitCovMat(1,2);
  GausCovMatrix(2,0) = FitCovMat(2,0);
  GausCovMatrix(2,1) = FitCovMat(2,1);
  GausCovMatrix(2,2) = FitCovMat(2,2);

  return GausCovMatrix;
}

std::vector<std::vector<float>> isolateAndFitPeaks( TH1F * hIn, float width, bool verbose, bool isSim, int nSamples )
{
  float mean = hIn->GetBinCenter( hIn->GetMaximumBin() );
//make a clone where we can iteratively remove what is to the left of peak of interest 
  TH1F * hToFit = (TH1F*)hIn->Clone();
  float oldMean = hIn->GetXaxis()->GetXmin()+1.01;
  float maxVal = hIn->GetXaxis()->GetXmax();
  float minVal = oldMean;
  minVal = TMath::Max(mean-width,minVal);
  float max=TMath::Min(maxVal, mean+width); //start a bit more narrow, we don't know distance between peaks yet                                              
  TF1 * fBkg = new TF1("fBkg", background, 0,20,3);
  hIn->Fit(fBkg, "RS", "", mean, maxVal); // get the background fit by fitting from first peak to end
  std::cout << "Starting with fit interval (" << minVal << ", " << max << ")" << std::endl;
  TF1 * fGaus = new TF1("fGaus", "gaus", hIn->GetXaxis()->GetXmin(), hIn->GetXaxis()->GetXmax());
  //TF1 * fGaus = new TF1("fGaus","gaus",minVal, max);
  TF1 *fGausExp = new TF1("fGausExp", ExponentialGaus, hIn->GetXaxis()->GetXmin(), hIn->GetXaxis()->GetXmax(), 5);
  TF1 *fGausBac = new TF1("fGausBac", ExponentialBackground, hIn->GetXaxis()->GetXmin(), 10+hIn->GetXaxis()->GetXmax(), 2);
  //  TF1 * fGaus = new TF1("fGaus", lorentzianPeak, hIn->GetXaxis()->GetXmin(), hIn->GetXaxis()->GetXmax(), 3);
  //TF1 * fGaus = new TF1("fGaus", lorentzianPeak, minVal, maxVal, 3);
  //TF1 * fGaus = new TF1("fGaus", asymPeak, hIn->GetXaxis()->GetXmin(), hIn->GetXaxis()->GetXmax(), 3);
  //TF1 * fGaus = new TF1("fGaus", asymPeak, minVal, max, 3);
  vector<TFitResultPtr> vPtrs;
  TFitResultPtr ptrBack;
  int iPstart = 0;
  if (isSim)
    iPstart=1; //no real pedestal peak in MC
  
  int iP=iPstart;
  std::vector<float> Con;
  std::vector<float> Mean;
  std::vector<float> Sigma;
  std::vector<float> LowLimit;
  std::vector<float> HighLimit;
 	std::vector<float> ExpConst;
 	std::vector<float> ExpScale;
 	std::vector<float> GausIntegral;
 	std::vector<float> GausIntegralError;

  int norm = 0;
  float sigma=0;                                                                                                        
  bool hasAdjustedWidth=false;

  //fGausBac->SetParameter(0,)
  //fGausBac->SetParLimits(0,0.0,10e5);
  //fGausBac->SetParLimits(1,0.0,5.0);  //Fitting common exponential background
  //ptrBack = hIn->Fit(fGausBac, "QRS", "same", oldMean+1.0, maxVal);//mean-width, mean+width);
  std::cout << "Parameters for the fix fitting is " << oldMean << "::" << maxVal << std::endl;
  while ( max <= maxVal ) {
    if (verbose) {
      std::cout<<"Fitting for peak " << iP << ": around mean=" << mean
               << " between " << minVal << " and " << max
               <<  std::endl;
    }

    //if (iP == 0) {
    	//fGausBac->SetParLimits(0,0.0,10e5);
  		//fGausBac->SetParLimits(1,0.0,5.0);
  		//ptrBack = hIn->Fit(fGausBac, "QRS", "same", oldMean, maxVal);//mean-width, mean+width);
    //}

    //fGaus->SetParameter(0,hToFit->GetBinContent(hToFit->FindBin(mean)));
    //fGaus->SetParameter(1,mean);
    //fGaus->SetParameter(2,0.5*(max-minVal));
    //fGaus->SetParLimits(1,minVal,max);

    fGausExp->SetParameter(0,hToFit->GetBinContent(hToFit->FindBin(mean)));
    fGausExp->SetParameter(1,mean);
    fGausExp->SetParameter(2,0.5*(max-minVal));
    fGausExp->SetParLimits(0,0.0,hToFit->GetBinContent(hToFit->FindBin(mean))); 
    fGausExp->SetParLimits(1,minVal,max);
    fGausExp->SetParLimits(2,0.01,2*(max-minVal));
    fGausExp->SetParLimits(3,0.0,10e5);
    //fGausExp->FixParameter(3,ptrBack->Parameter(0)); //Using the parameters of the common exponential to fit the gaussian peaks
    //fGausExp->FixParameter(4,ptrBack->Parameter(1));
    fGausExp->SetParLimits(4,0.0,5.0);
    TFitResultPtr ptr = hToFit->Fit(fGausExp, "QRS", "same", minVal, max);//mean-width, mean+width);
    //if (iP > iPstart) ptr = hToFit->Fit(fGausExp, "QRS", "same", minVal, max);//mean-width, mean+width);
    std::cout << "Some check params :" << mean << "::" << width << "::" << maxVal << "::" << oldMean << "::" << sigma << std::endl;
    std::cout << "Checking fit params :" << ptr->Parameter(0) << "::" << ptr->Parameter(1) << "::" << ptr->Parameter(2) << std::endl;     
	if (ptr->Parameter(1) < oldMean || ptr->Parameter(2) < 0.2*sigma) {//ptr is not 0 --> not converged, try again with narrower range, given that we probably nailed the peak already                                                    
 	    std::cout << "Ok, trying to narrow down the peak until the two sigma width is smaller than the difference between max and min" << std::endl;
 	    minVal=mean-0.4*width;
      max=mean+0.4*width;
      if (verbose) {
        std::cout<<"\tFitting again for peak " << iP << ": around mean=" << mean
                 << " between " << minVal << " and " << max
                 <<  std::endl;
      }
      ptr = hToFit->Fit(fGausExp, "RS", "same",  minVal, max);// mean-0.8*width, mean+0.8*width);
      //if (iP > iPstart) ptr = hToFit->Fit(fGausExp, "RS", "same",  minVal, max);// mean-0.8*width, mean+0.8*width);                                                                    
   }
   
   if (iP > iPstart && !hasAdjustedWidth) {
	 width=(ptr->Parameter(1)-oldMean)/2.; //total width is half distance between the two peaks                                                             
	 hasAdjustedWidth=true;
	 //TF1 *fGaus = new TF1("fGaus", ExponentialGaus, hIn->GetXaxis()->GetXmin(), hIn->GetXaxis()->GetXmax(), 5);
	 if (verbose) {
	   std::cout<<"\tUpdating width to " << width
				<<  std::endl;
	 }
   }
   mean=ptr->Parameter(1);
   sigma=ptr->Parameter(2);
   LowLimit.push_back(minVal);
   HighLimit.push_back(max);
   //if (iP > iPstart) {}
  // }
  // else {
   	//ExpConst.push_back(0.0);
   	//ExpScale.push_back(0.0);
  // }
   TF1 *GausInt = new TF1("GausInt",gaussianPeak,minVal,max,3);
   TMatrixD covTot = ptr->GetCovarianceMatrix();
   TMatrixD GausCovMatrix = ConstructCovarianceMatrix(covTot);
   std::cout << "Covariance Matrix before sub-filling" << std::endl;
   covTot.Print();
   std::cout << "Covariance Matrix after sub-filling" << std::endl;
   GausCovMatrix.Print();
   GausInt->SetParameters(ptr->Parameter(0),ptr->Parameter(1),ptr->Parameter(2));
   GausInt->SetParError(0,ptr->ParError(0));
   GausInt->SetParError(1,ptr->ParError(1));
   GausInt->SetParError(2,ptr->ParError(2));
   float IntLow = mean - 3*sigma;
   float IntHigh = mean + 3*sigma;
   if (iP == 0)
    IntLow = mean - sigma;
   float PeakIntegral = GausInt->Integral(IntLow,IntHigh);
   float PeakIntegralError = GausInt->IntegralError(IntLow,IntHigh,GausInt->GetParameters(),GausCovMatrix.GetMatrixArray());
   std:;cout << "Checking some integral limits :: " << mean << ";" << sigma << ";" << minVal << ";" << max << std::endl;
   fGaus->DrawCopy("same");
   //could keep track of last bin which was reset, for some speed gain
   std::cout << "Trying to eliminate elements upto " << hToFit->FindBin( mean+TMath::Min(width,float(7*sigma))) << " using width " << width << std::endl;
   for (int iB=1; iB<hToFit->FindBin( mean+TMath::Min(width,float(7*sigma)) ); iB++) {
   	hToFit->SetBinContent(iB, 0);
    //std::cout << hToFit->GetBinContent(iB) << std::endl;
   }   
   oldMean=mean;
   mean=hToFit->GetBinCenter( hToFit->GetMaximumBin() );
   std::cout << "Mean after modification is " << mean << " with minVal " << minVal << std::endl;
   //minVal=TMath::Max(minVal,float(mean-0.5*width));
   minVal = mean-0.5*width;
   max=mean+width;
   std::cout << "Now minVal is " << minVal << " and max is " << max << std::endl;      								
   //norm+=ptr->Parameter(0);   //keep track of how much stats we have left to work with
   if (verbose) {
	 std::cout<<"\tUpdating sum of peak heights to " << norm
			  <<  std::endl;
   }
   // should probably assess fit quality somewhere here too                                                                                                   
   vPtrs.push_back( ptr );
   if (verbose) {
	 std::cout<<"\t\tFor peak " << iP << ", got mean=" << ptr->Parameter(1) <<
	   "\tand sigma=" << ptr->Parameter(2) <<  std::endl;
	 std::cout<<"\t\tUpdated mean to: " <<  mean << std::endl;
   }
   std::cout<<"\t\tSaving parameters for peak " << iP << ", got mean=" << vPtrs.at(iP)->Parameter(1) <<
	   "\tand sigma=" << vPtrs.at(iP)->Parameter(2) <<  std::endl;
   Con.push_back(vPtrs.at(iP)->Parameter(0));
   Mean.push_back(vPtrs.at(iP)->Parameter(1));
   Sigma.push_back(vPtrs.at(iP)->Parameter(2));
   ExpConst.push_back(vPtrs.at(iP)->Parameter(3));
   ExpScale.push_back(vPtrs.at(iP)->Parameter(4));
   GausIntegral.push_back(PeakIntegral);
   GausIntegralError.push_back(PeakIntegralError);
   std::cout << "ChiSquare for the fit is " << fGausExp->GetChisquare()/fGausExp->GetNDF() << " for peak " << iP << std::endl;
   //Integrals.push_back(PeakIntegral);
   //IntegralsError.push_back(PeakIntegralError);
   iP++;
   if ( iP > 10 ) //more than enough, and avoid eternal while loop                                                                                            
	 	break;
   if (!isSim) {
   	float bkgLevel=0;
   	bkgLevel=fBkg->Eval(ptr->Parameter(1));
   //   if (hIn->GetEntries()-norm < 0.05*hIn->GetEntries() || ptr->Parameter(0)-bkgLevel< 8) { //don't fit peaks with just a few entries                          
  	if (hToFit->Integral() < 0.0001*hIn->Integral() || ptr->Parameter(0)-bkgLevel< 8 || hToFit->Integral() < 10) { //don't fit peaks with just a few entries                          
			std::cout<<"\tHitting stats break factor at peak integral " << ptr->Parameter(0)
			  << " and (in data) fitted background level " << bkgLevel
			  << " and histogram entries " << hToFit->Integral()
			  << " out of total from start " << hIn->Integral()
			  <<  std::endl;
	 		break;
   	}

   }
	 //bkgLevel=fBkg->Eval(ptr->Parameter(1));
   //   if (hIn->GetEntries()-norm < 0.05*hIn->GetEntries() || ptr->Parameter(0)-bkgLevel< 8) { //don't fit peaks with just a few entries                          
   //if (hToFit->Integral() < 0.01*hIn->Integral() || ptr->Parameter(0)-bkgLevel< 8 || hToFit->Integral() < 10) { //don't fit peaks with just a few entries                          
	 //std::cout<<"\tHitting stats break factor at peak integral " << ptr->Parameter(0)
	//		  << " and (in data) fitted background level " << bkgLevel
		//	  << " and histogram entries " << hToFit->Integral()
			//  << " out of total from start " << hIn->Integral()
			 // <<  std::endl;
	 //break;
   //}
   
  }

  std::vector<std::vector<float>> ParArray;
  ParArray.push_back(LowLimit);
  ParArray.push_back(HighLimit);
  ParArray.push_back(Con);
  ParArray.push_back(Mean);
  ParArray.push_back(Sigma);
  ParArray.push_back(ExpConst);
  ParArray.push_back(ExpScale);
  ParArray.push_back(GausIntegral);
  ParArray.push_back(GausIntegralError);
  //const int nPoints = vPtrs.size();
  //double x[(const int) nPoints];
  //double ex[(const int) nPoints];
  //double y[(const int) nPoints];
  //double ey[(const int) nPoints];
  //float noise=0;
  //for ( int iP=0; iP<nPoints; iP++) {
	//x[iP] = iPstart+iP;
	//ex[iP] = 0;
	//y[iP] = vPtrs.at(iP)->Parameter(1);
	//ey[iP] = vPtrs.at(iP)->Parameter(2);
	//if ( iP > 0 ) { //extract the noise level from relative size of peaks
	  //std::cout << "Estimate Lambda(peak " << iP << ") = " << vPtrs.at(iP)->Parameter(0)/vPtrs.at(iP-1)->Parameter(0)/(iP*nSamples) << std::endl;
	  //noise+=vPtrs.at(iP)->Parameter(0)/vPtrs.at(iP-1)->Parameter(0)/(iP*nSamples);
	//}
 //}
  //std::cout<<"Average noise/time sample from all used peaks = " << noise/(nPoints-1) <<std::endl;
  //std::cout<<"Average noise/event from all used peaks = " << noise/(nPoints-1)*nSamples <<std::endl;
  
  //TGraphErrors * g= new TGraphErrors( nPoints, x, y, ex, ey);
  
  return ParArray; 
}

void extractPEPeakIntegrals_MultipleRootFiles(TString File1, TString File2, TString File3, TString location, int deadChannel=-1, string decodePassName = "hits", int nSamples = 30, int verbosity=2, bool isSim=false, bool doClean=true)
//void extractPEPeakIntegrals_MultipleRootFiles(TString File1, TString File2, TString location, int deadChannel=-1, string decodePassName = "hits", int nSamples = 30, int verbosity=2, bool isSim=false, bool doClean=true)
//void extractPEPeakIntegrals_MultipleRootFiles(string inFile, TString location, int deadChannel=-1, string decodePassName = "hits", int nSamples = 30, int verbosity=2, bool isSim=false, bool doClean=true)
{
  // macro extracting the per-channel gain and pedestal based on total Q histograms, and single PE peak fitting

  bool verbosePrint = (verbosity == 2 ); //adhere to ldmx printout lavel numbering
  
  //TFile * fIn = TFile::Open( inFile.c_str(), "READ");
  //TTree * tree = (TTree*)fIn->Get("LDMX_Events");
  //fIn->ls();

  TString fIn1 = File1 + ".root";
  TString fIn2 = File2 + ".root";
  TString fIn3 = File3 + ".root";
  TChain *chain = new TChain("LDMX_Events");
  chain->Add(fIn1);
  chain->Add(fIn2);
  chain->Add(fIn3);
  chain->ls();
  TTree * tree = chain;
  tree->ls();
  int exampleEvNb = 1;

  string digiName="testBeamHitsUp";
  int nChannels = 12;
  
  vector<string> vars = {"pe"}; // this variable is totQ/nSamp
  vector <int> maxVals = {4};  // don't need to go high, focus on the single-few PE peaks 
  vector <int> minVals = {-1};  // should cover most negative pedestals 
  vector <float> binFactor = {150.0/5}; // histogram binning 

  TCanvas * c1 = new TCanvas("c1", "plot canvas", 600, 500);

  vector<string> cuts;
  vector<TH1F*> v_hOut;

  string cleanStr="";
  if (doClean) //no need to do flag == 4, it's not set until hit reconstruction anyway
    //cleanStr=Form("%s_%s.flag_==0 &&", digiName.c_str(),decodePassName.c_str());
  	cleanStr=Form("(%s_%s.flag_ == 0 || %s_%s.flag_ == 4) &&", digiName.c_str(),decodePassName.c_str(), digiName.c_str(),decodePassName.c_str());
  
  for (int iC = 0 ; iC < nChannels ; iC++)
	cuts.push_back( Form("%s %s_%s.barID_==%i", cleanStr.c_str(), digiName.c_str(),decodePassName.c_str(), iC));

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
      //std::cout << (nSamp+"*("+digiName+"_"+decodePassName+"."+vars[iV]+"_) >> h"+bin+"("+bins+")").c_str() << std::endl;
      std::cout << cut.c_str() << std::endl;
	  //tree->Draw( (nSamp+"*("+digiName+"_"+decodePassName+"."+vars[iV]+"_) >> h"+bin+"("+bins+")").c_str(), cut.c_str() );
	  tree->Draw( (digiName+"_"+decodePassName+"."+vars[iV]+"_ >> h"+bin+"("+bins+")").c_str(), cut.c_str() );
	  //get them each and keep for later
	  TH1F *hOut = (TH1F*)gDirectory->Get(Form("h%i", iC)); 
	  if (!hOut || hOut->IsZombie()) {
		std::cout << "No histogram for channel " << iC << ", skipping" << std::endl;
		continue;
	  }
	  //set up axes 	  
	  string title = ";"+vars[iV]+" "+cut+";entries";
	  hOut->SetTitle( title.c_str() );
	  //hOut->GetYaxis()->SetRangeUser(0.001, 1.5*hOut->GetMaximum()) ;
	  //draw 
	  hOut->Draw();
	  hOut->SetName(Form("h%s_chanID%i", vars[iV].c_str(), iC));
	  v_hOut.push_back( (TH1F*)hOut->Clone() );
	}//over cuts
  }//over variables 


  // ok so now we have all the histograms. we just need to fit them...
  // it is reasonable to assume that we will have a hunch how wide each PE peak will be.
  // we can pass that param to the fitting function. 

  float fitRangeWidth=0.5; // take this window to either side of mean to catch peak 
  //vector <TGraphErrors*> v_g;
  //TGraphErrors* gGains = new TGraphErrors(nChannels);
  //TGraphErrors* gPeds = new TGraphErrors(nChannels);
  //float scaleFac=6250.; //conversion from fC to e
  TString TitleAll("FitAll_");
  TString PlotTitle("NoiseFits_");

  TString IntegralFile, IntegralErrorFile;
  TString BaseFileName = File1.ReplaceAll("_TS_0_4","_TS_0_14_combined");
  IntegralFile = BaseFileName + "_Integral.txt";
  IntegralErrorFile = BaseFileName + "_IntegralError.txt";
  ofstream fIntegral;// gainFile.Data());
  fIntegral.open( IntegralFile.Data());
  ofstream fIntegralError;
  fIntegralError.open( IntegralErrorFile.Data());
 
  for (unsigned int iH = 0;  iH < v_hOut.size(); iH++) {
	if (iH == deadChannel)
	  continue;
	//do the actual peak fitting 
	std::cout<< "----> Fitting peaks for channel " << iH << std::endl;
	//	TGraphErrors * g= findAndFitPeaks( v_hOut.at(iH), fitRangeWidth, verbosePrint, isSim); 
	TString SaveNameFit = TitleAll + iH;
	TString TitleFit = PlotTitle + iH;
	auto ParameterArray = isolateAndFitPeaks( v_hOut.at(iH), fitRangeWidth, verbosePrint, isSim, nSamples);
	DrawGaussians(v_hOut.at(iH),TitleFit,ParameterArray.at(0),ParameterArray.at(1),ParameterArray.at(2),ParameterArray.at(3),ParameterArray.at(4),ParameterArray.at(5),ParameterArray.at(6),SaveNameFit,location); 
  fIntegral << iH << "," << ParameterArray.at(7)[0] << "," << ParameterArray.at(7)[1] << "," << ParameterArray.at(7)[2] << "," << ParameterArray.at(7)[3] << "\n";
  fIntegralError << iH << "," << ParameterArray.at(8)[0] << "," << ParameterArray.at(8)[1] << "," << ParameterArray.at(8)[2] << "," << ParameterArray.at(8)[3] << "\n";
  }
  fIntegral.close();
  fIntegralError.close();
}


TGraphErrors* findAndFitPeaks( TH1F * hIn, float width, bool verbose, bool isSim )
{

  float mean = hIn->GetBinCenter( hIn->GetMaximumBin() );
  float maxVal = hIn->GetXaxis()->GetXmax();
  std::cout << "Limit on maxVal = " << maxVal << std::endl;
  float oldMean = hIn->GetXaxis()->GetXmin()+5;
  float minVal = oldMean;
  minVal = TMath::Max(mean-width,minVal);
	float max=TMath::Min(maxVal, mean+width); //start a bit more narrow, we don't know distance between peaks yet
  std::cout << "Starting with fit interval (" << minVal << ", " << max << ")" << std::endl;
  vector<TFitResultPtr> vPtrs;
  int iPstart = 0;
  if (isSim)
	iPstart=1; //no real pedestal peak in MC 
  int iP=iPstart;
  int norm = 0;
  TF1 * fGaus = new TF1("fGaus", "gaus", -50, 5000);
  TF1 * fBkg = new TF1("fBkg", background, -50, 5000,3);
  hIn->Fit(fBkg, "RS", "", mean+10, maxVal); // get the background fit by fitting from first peak to end
  //
  /*
  TF1 * fSum= new TF1("fSum", "fBkg+gaus(3)", mean, 5000);
  //*/
  TF1 *fSum = new TF1("fSum",fitFunction,mean, 5000, 6);
  fSum->FixParameter(3, 0);//fBkg->GetParameter(0));
  fSum->FixParameter(4, 0);//fBkg->GetParameter(1));
  fSum->FixParameter(5, 0);//fBkg->GetParameter(2));
  bool hasAdjustedWidth=false;
  while ( max <= maxVal ) {
	if (verbose) {
	  std::cout<<"Fitting for peak " << iP << ": around mean=" << mean
			   << " between " << minVal << " and " << max 
			   <<  std::endl;
	}
	TFitResultPtr ptr = hIn->Fit(fGaus, "QRS", "same", minVal, max);//mean-width, mean+width);
	//potentially adjust width based on fit result?
	//parameter order: norm, mean, sigma
	if (ptr || ptr->Parameter(1) < oldMean) {//ptr is not 0 --> not converged, try again with narrower range
	  minVal=mean-1.2*width;
	  max=mean+0.8*width;
	  if (verbose) {
		std::cout<<"\tFitting again for peak " << iP << ": around mean=" << mean
				 << " between " << minVal << " and " << max
				 <<  std::endl;
	  }
	  ptr = hIn->Fit(fGaus, "RS", "same",  minVal, max);// mean-0.8*width, mean+0.8*width);
	}
	mean=ptr->Parameter(1);
	if (ptr->Parameter(1) < oldMean+width and iP > 2) //we got some hits already, mean is decreasing so this is getting wonky, interrupt
	  break;
	if ( oldMean < 0) //need to do something different in the start, where we don't know the 1PE interval between peaks
	  mean+= 1.5*width; //width is not yet adjusted based on peak distance 
	else {
	  if (verbose) {
		std::cout<<"\tUpdated next mean estimate from fit result " << mean << " to: " <<  mean+ptr->Parameter(1)-oldMean << " using old mean " << oldMean<< std::endl;
	  }
	  mean+=ptr->Parameter(1)-oldMean;
	  if (!hasAdjustedWidth) {
		width=(ptr->Parameter(1)-oldMean)/4.; //total width is half distance between the two peaks 
		hasAdjustedWidth=true;
		if (verbose) {
		  std::cout<<"\tUpdating width to " << width
				   <<  std::endl;
		}
	  }
	  fGaus->DrawCopy("same");
	}
	oldMean=ptr->Parameter(1);
	minVal=TMath::Max(minVal, float(mean-1.25*width));
	max=TMath::Min(maxVal, float(mean+1.25*width)); 

	norm+=ptr->Parameter(0);//ptr->Parameter(0); //TMath::Sqrt(ptr->Parameter(2)*TMath::Pi())*ptr->Parameter(0); //keep track of how much stats we have left to work with
    if (verbose) {
	  std::cout<<"\tUpdating sum of peak heights to " << norm 
			   <<  std::endl;
	}
	// should probably assess fit quality somewhere here too 
	vPtrs.push_back( ptr );
    if (verbose) {
	  std::cout<<"\t\tFor peak " << iP << ", got mean=" << ptr->Parameter(1) <<
		"\tand sigma=" << ptr->Parameter(2) <<  std::endl;
	  std::cout<<"\t\tUpdated mean to: " <<  mean << std::endl;
	}
	iP++;
	if ( iP > 10 ) //more than enough, and avoid eternal while loop
	  break;
	float bkgLevel=0;
	if (!isSim)
	  bkgLevel=fBkg->Eval(ptr->Parameter(1));  
	if (hIn->GetEntries()-norm < 0.05*hIn->GetEntries() || ptr->Parameter(0)-bkgLevel< 8) { //don't fit peaks with just a few entries 
	  std::cout<<"\tHitting stats break factor at peak integral " << ptr->Parameter(0)
			   << " and (in data) fitted background level " << bkgLevel
             << " and histogram entries " << hIn->GetEntries()
             <<  std::endl;
	  break;
	}
  }	

    const int nPoints = vPtrs.size();
    double x[(const int) nPoints];
    double ex[(const int) nPoints];
    double y[(const int) nPoints];
    double ey[(const int) nPoints];
    for ( int iP=0; iP<nPoints; iP++) {
      x[iP] = iPstart+iP;
      ex[iP] = 0;
      y[iP] = vPtrs.at(iP)->Parameter(1);;
      ey[iP] = vPtrs.at(iP)->Parameter(2);;
    }

    TGraphErrors * g= new TGraphErrors( nPoints, x, y, ex, ey);

	return g;
}

