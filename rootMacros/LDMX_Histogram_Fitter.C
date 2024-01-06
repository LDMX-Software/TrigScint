//Root file to fit functions to LDMX histograms
using namespace std;

double HistogramFitterLDMX(TFile *hzz, Int_t bins, float_t xmin, float_t xmax, TString &hsig, TString &plottitle) 
{
	TCanvas *c2 = new TCanvas();

	TH1F *sig = new TH1F("Signal",hsig,bins,xmin,xmax);
	sig = (TH1F*)hzz->Get(hsig);

	TF1 *Nx = new TF1("Nx","gaus",50,120);
	sig->Fit("Nx","R");
	double MIP = Nx->GetParameter(1);
	//double MIP_spread = Nx->GetParameter(2);
	sig->Draw();
	TString outs1 = "/home/dhruvanshu/ldmx-sw/LDMX_Plots/GaussianFit_" + plottitle + ".png";
	c2->Print(outs1);
	sig->Clear();
	return MIP;
}

void LDMX_Histogram_Fitter(const char *sigfile)
{
	TFile *hzz = TFile::Open(sigfile);
	double MIP_peak;
	std::vector<double> MIPs;
	std::vector<double> bars;
	//std::vector<double> MIPsigma;
	TString plot_title_PEmed2("hPE_med2_bar");
	for(unsigned int iB=0;iB < 16; iB++)
	{
		TString plot_PEmed2 = "Hits/" + plot_title_PEmed2 + iB;
		TString plot_PEmed2_title = plot_title_PEmed2 + iB;
		MIP_peak = HistogramFitterLDMX(hzz,200,10,150,plot_PEmed2,plot_PEmed2_title);
		bars.push_back(iB);
		MIPs.push_back(MIP_peak);
		//MIPsigma.push_back(MIP_SD);
	}

	TCanvas *c1 = new TCanvas();
	auto gr1 = new TGraph();
	for(unsigned int i=0;i<bars.size();i++){
		gr1->SetMarkerColor(4);
		gr1->SetMarkerStyle(20);
		gr1->SetTitle("MIP peaks;Bar_ID;MIP PE");
		gr1->AddPoint(bars[i],MIPs[i]);
		//std::cout << bars[i] << ":" << MIPs[i] << std::endl;
	}
	gr1->GetYaxis()->SetRangeUser(70,100);
	gr1->Draw("");
	TString outs = "/home/dhruvanshu/ldmx-sw/LDMX_Plots/MIP_peaks_50_120.png";
	c1->Print(outs);

}

