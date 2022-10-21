//Root file to print out histograms plotted using LDMX framework
using namespace std;

void HistogramPlotterLDMX(TFile *hzz, Int_t bins, float_t xmin, float_t xmax, TString &plot_title, TString &hsig, TString &location, Int_t k)
{
	TCanvas *c1 = new TCanvas();

	gStyle->SetLabelSize(0.048,"x");
	gStyle->SetLabelSize(0.048,"y");
	gStyle->SetTitleSize(0.045,"x");
	gStyle->SetTitleSize(0.045,"y");

	TH1F *sig = new TH1F("Signal",hsig,bins,xmin,xmax);
	sig = (TH1F*)hzz->Get(hsig);

	sig->SetLineWidth(2);
	sig->SetLineColor(kBlue);
	sig->SetTitle(plot_title);
	if (k==1) {
		c1->SetLogy();
	}

	//TString text = "Number of bins : " + to_string(bins);
	//TLatex *t = new TLatex(160,150000,text);
	sig->Draw();
	//t->Draw();

	if (k==1) {
		TString outs = "/home/dhruvanshu/ldmx-sw/LDMX_Plots/" + location + plot_title + "_log.png";
		c1->Print(outs);
	}
	else {
		TString outs = "/home/dhruvanshu/ldmx-sw/LDMX_Plots/" + location + plot_title + ".png";
		c1->Print(outs);
	}
	
	hsig.Clear();

}

void HistogramPlotter2D(TFile *hzz, Int_t bins_x, float_t xmin, float_t xmax, Int_t bins_y, float_t ymin, float_t ymax, TString &plot_title, TString &hsig, TString &location, Int_t j)
{
	TCanvas *c2 = new TCanvas();

	gStyle->SetLabelSize(0.048,"x");
	gStyle->SetLabelSize(0.048,"y");
	gStyle->SetTitleSize(0.045,"x");
	gStyle->SetTitleSize(0.045,"y");
	gStyle->SetPalette(kRainBow);

	TH2F *sig = new TH2F("Signal",hsig,bins_x,xmin,xmax,bins_y,ymin,ymax);
	sig = (TH2F*)hzz->Get(hsig);

	if (j==0) {
		sig->SetTitle(plot_title);
		sig->SetMarkerStyle(45);
		sig->SetMarkerColor(kBlue);
		sig->SetMarkerSize(0.3);
		//std::cout << sig->Integral() << std::endl;		
		sig->SetStats(0);
		c2->Update();
		sig->Draw();
		//c2->SetLogy();
	}

	if (j==1) {
		sig->GetYaxis()->SetRangeUser(63.8,64.3);
		sig->SetTitle(plot_title);
		sig->SetMarkerStyle(45);
		sig->SetMarkerColor(kBlue);
		sig->SetMarkerSize(0.2);
		sig->SetStats(0);
		sig->Draw();
	}

	if (j==2) {
		sig->SetTitle(plot_title);
		sig->SetMarkerStyle(45);
		sig->SetMarkerColor(kBlue);
		sig->SetMarkerSize(0.3);
		sig->SetStats(0);
		c2->SetLogy();
		sig->Draw();
	}

	if (j==3) {
		sig->SetTitle(plot_title);
		sig->SetStats(0);
		gStyle->SetPalette(kCool);
		c2->SetLogz();
		sig->Draw("colz");
	}

	if (j==4) {
		sig->SetTitle(plot_title);
		sig->SetStats(0);
		gStyle->SetPalette(kCool);
		//c2->SetLogz();
		sig->Draw("colz");
	}

	//sig->Draw();


	if (j==2) {
		TString outs = "/home/dhruvanshu/ldmx-sw/LDMX_Plots/" + location + plot_title + "_log.png";
		c2->Print(outs);
		hsig.Clear();
	}
	if (j==3 || j==4) {
		TString outs = "/home/dhruvanshu/ldmx-sw/LDMX_Plots/" + location + plot_title + "_color.png";
		c2->Print(outs);
		hsig.Clear();
	}
	else {
		TString outs = "/home/dhruvanshu/ldmx-sw/LDMX_Plots/" + location + plot_title + ".png";
		c2->Print(outs);
		hsig.Clear();
	}
} 

void LDMX_Histogram_Plotter(const char *sigfile)
{
	TFile *hzz = TFile::Open(sigfile);
	TString location_ADC("24-04-2022_08-41-40__277/ADC/");
	TString location_Charge("24-04-2022_08-41-40__277/Charge/");
	TString location_Hits("24-04-2022_08-41-40__277/Hits/");

	int Qlow_thr = 2.e2;
    int Qmed_thr = 5.e3;
    int Qmed_thr2 = 30.e3;
    int Qhigh_thr = 10.e3;
    int Qhigh_thr2 = 64.e3;
    	
    int PE_low = 1;
    int PE_med1 = 10;
    int PE_med2 = 150;
    int PE_high = 400;

    int nQbins_low = Qlow_thr;
  	int nQbins_med = (Qmed_thr-Qlow_thr);
  	int nQbins_high = (Qhigh_thr-Qmed_thr);
  	int nPE_bins = 200;
  	int PEmax = 1200;

	TString plot_title_ADC("hADC_chan");
	TString plot_title_TDC("hTDC_chan");
	TString plot_title_ADC_total("hADC_Total");
	TString plot_title_TDC_total("hTDC_Total");
	TString plot_title_ADCvsTS("hADCvsTS_chan");
	TString plot_title_TDCvsTS("hTDCvsTS_chan");
	//TString plot_title_QvsTS("hQvsTS_chan");
	TString plot_ADC_total = "ADC/" + plot_title_ADC_total;
	TString plot_TDC_total = "ADC/" + plot_title_TDC_total;

	TString plot_title_Qlow("hQ_low_chan");
	TString plot_title_Qmed("hQ_med_chan");
	TString plot_title_Qhigh("hQ_high_chan");
	TString plot_title_QTot_low("hQTot_low_chan");
	TString plot_title_QTot_med("hQTot_med_chan");
	TString plot_title_QTot_med2("hQTot_med2_chan");
	TString plot_title_QTot_high("hQTot_high_chan");
	TString plot_title_QTotal("hQTotal_channel");
	TString plot_title_QlowvsTS("hQvsTS_low_chan");
	TString plot_title_QmedvsTS("hQvsTS_med_chan");
	TString plot_title_QhighvsTS("hQvsTS_high_chan");
	TString plot_title_QTotvschan_low("hQTotvschan_low");
	TString plot_title_QTotvschan_low_ver2("hQTotvschan_low_ver2");
	TString plot_title_QTotvschan_med("hQTotvschan_med");
	TString plot_title_QTotvschan_med2("hQTotvschan_med2");
	TString plot_title_QTotvschan_high("hQTotvschan_high");
	TString plot_QTotvschan_low = "Charge/" + plot_title_QTotvschan_low;
	TString plot_QTotvschan_low_ver2 = "Charge/" + plot_title_QTotvschan_low_ver2;
	TString plot_QTotvschan_med = "Charge/" + plot_title_QTotvschan_med;
	TString plot_QTotvschan_med2 = "Charge/" + plot_title_QTotvschan_med2;
	TString plot_QTotvschan_high = "Charge/" + plot_title_QTotvschan_high;
	//TString plot_title_QTotalvschannel("hQTotvschan");

	TString plot_title_PElow("hPE_low_bar");
	TString plot_title_PEmed1("hPE_med1_bar");
	TString plot_title_PEmed2("hPE_med2_bar");
	TString plot_title_PEhigh("hPE_high_bar");
	TString plot_title_PE("hPE");
	TString plot_title_PE2D_low("hPE2D_low");
	TString plot_title_PE2D_med1("hPE2D_med1");
	TString plot_title_PE2D_med2("hPE2D_med2");
	TString plot_title_PE2D_high("hPE2D_high");
	TString plot_title_PE2D("hPE2D");
	TString plot_PE2D_low = "Hits/" + plot_title_PE2D_low;
	TString plot_PE2D_med1 = "Hits/" + plot_title_PE2D_med1;
	TString plot_PE2D_med2 = "Hits/" + plot_title_PE2D_med2;
	TString plot_PE2D_high = "Hits/" + plot_title_PE2D_high;
	TString plot_PE2D = "Hits/" + plot_title_PE2D;
	TString plot_title_ChanvsElec("hElecvsChan");
	TString plot_title_BarvsChan("hBarvsChan");
	TString plot_ChanvsElec = "Charge/" + plot_title_ChanvsElec;
	TString plot_BarvsChan = "Charge/" + plot_title_BarvsChan;

	HistogramPlotterLDMX(hzz,80,0,240,plot_title_ADC_total,plot_ADC_total,location_ADC,0);
	HistogramPlotterLDMX(hzz,20,60,70,plot_title_TDC_total,plot_TDC_total,location_ADC,0);
	HistogramPlotter2D(hzz,15,0,15,100,-1,PE_low,plot_title_QTotvschan_low,plot_QTotvschan_low,location_Charge,0);
	HistogramPlotter2D(hzz,15,0,15,100,0,PE_low,plot_title_QTotvschan_low_ver2,plot_QTotvschan_low_ver2,location_Charge,3);
	HistogramPlotter2D(hzz,15,0,15,100,PE_low,PE_med1,plot_title_QTotvschan_med,plot_QTotvschan_med,location_Charge,0);
	HistogramPlotter2D(hzz,15,0,15,100,PE_med1,PE_med2,plot_title_QTotvschan_med2,plot_QTotvschan_med2,location_Charge,0);
	HistogramPlotter2D(hzz,15,0,15,100,PE_med2,PE_high,plot_title_QTotvschan_high,plot_QTotvschan_high,location_Charge,0);
	HistogramPlotter2D(hzz,15,0,15,nPE_bins,0,PE_low,plot_title_PE2D_low,plot_PE2D_low,location_Hits,3);
	HistogramPlotter2D(hzz,15,0,15,nPE_bins,PE_low,PE_med1,plot_title_PE2D_med1,plot_PE2D_med1,location_Hits,3);
	HistogramPlotter2D(hzz,15,0,15,nPE_bins,PE_med1,PE_med2,plot_title_PE2D_med2,plot_PE2D_med2,location_Hits,3);
	HistogramPlotter2D(hzz,15,0,15,nPE_bins,PE_med2,PEmax,plot_title_PE2D_high,plot_PE2D_high,location_Hits,3);
	HistogramPlotter2D(hzz,15,0,15,nPE_bins,0,PEmax,plot_title_PE2D,plot_PE2D,location_Hits,3);
	HistogramPlotter2D(hzz,15,0,15,15,0,15,plot_title_ChanvsElec,plot_ChanvsElec,location_Charge,4);
	HistogramPlotter2D(hzz,15,0,15,15,0,15,plot_title_BarvsChan,plot_BarvsChan,location_Charge,4);

	for(unsigned int iB=0;iB < 16; iB++){
		TString plot_ADC = "ADC/" + plot_title_ADC + iB;
		TString filename_ADC = plot_title_ADC + iB;
		TString plot_TDC = "ADC/" + plot_title_TDC + iB;
		TString filename_TDC = plot_title_TDC + iB;
		TString plot_ADC2D = "ADC/" + plot_title_ADCvsTS + iB;
		TString filename_ADC2D = plot_title_ADCvsTS + iB;
		TString plot_TDC2D = "ADC/" + plot_title_TDCvsTS + iB;
		TString filename_TDC2D = plot_title_TDCvsTS + iB;

		TString plot_Qlow = "Charge/" + plot_title_Qlow + iB;
		TString filename_Qlow = plot_title_Qlow + iB;
		TString plot_Qmed = "Charge/" + plot_title_Qmed + iB;
		TString filename_Qmed = plot_title_Qmed + iB;
		TString plot_Qhigh = "Charge/" + plot_title_Qhigh + iB;
		TString filename_Qhigh = plot_title_Qhigh + iB;
		TString plot_QTot_low = "Charge/" + plot_title_QTot_low + iB;
		TString filename_QTot_low = plot_title_QTot_low + iB;
		TString plot_QTot_med = "Charge/" + plot_title_QTot_med + iB;
		TString filename_QTot_med = plot_title_QTot_med + iB;
		TString plot_QTot_med2 = "Charge/" + plot_title_QTot_med2 + iB;
		TString filename_QTot_med2 = plot_title_QTot_med2 + iB;
		TString plot_QTot_high = "Charge/" + plot_title_QTot_high + iB;
		TString filename_QTot_high = plot_title_QTot_high + iB;
		TString plot_QTotal = "Charge/" + plot_title_QTotal + iB;
		TString filename_QTotal = plot_title_QTotal + iB;

		TString plot_PElow = "Hits/" + plot_title_PElow + iB;
		TString filename_PElow = plot_title_PElow + iB;
		TString plot_PEmed1 = "Hits/" + plot_title_PEmed1 + iB;
		TString filename_PEmed1 = plot_title_PEmed1 + iB;
		TString plot_PEmed2 = "Hits/" + plot_title_PEmed2 + iB;
		TString filename_PEmed2 = plot_title_PEmed2 + iB;
		TString plot_PEhigh = "Hits/" + plot_title_PEhigh + iB;
		TString filename_PEhigh = plot_title_PEhigh + iB;

		HistogramPlotterLDMX(hzz,80,0,240,filename_ADC,plot_ADC,location_ADC,1);
		HistogramPlotterLDMX(hzz,20,60,70,filename_TDC,plot_TDC,location_ADC,1);
		HistogramPlotter2D(hzz,29,0,29,80,0,240,filename_ADC2D,plot_ADC2D,location_ADC,0);
		//HistogramPlotter2D(hzz,30,0,29,80,0,240,filename_ADC2D,plot_ADC2D,2);
		HistogramPlotter2D(hzz,29,0,29,40,60,70,filename_TDC2D,plot_TDC2D,location_ADC,1);
		//HistogramPlotter2D(hzz,30,0,29,100,-10.e3,350.e3,filename_Q2D,plot_Q2D,0);
		HistogramPlotterLDMX(hzz,nQbins_low,-100,Qlow_thr,filename_Qlow,plot_Qlow,location_Charge,0);
		HistogramPlotterLDMX(hzz,nQbins_med/10,Qlow_thr,Qmed_thr,filename_Qmed,plot_Qmed,location_Charge,0);
		HistogramPlotterLDMX(hzz,nQbins_high/15,Qmed_thr,Qhigh_thr,filename_Qhigh,plot_Qhigh,location_Charge,0);
		HistogramPlotterLDMX(hzz,500,-1,PE_low,filename_QTot_low,plot_QTot_low,location_Charge,0);
		HistogramPlotterLDMX(hzz,100,PE_low,PE_med1,filename_QTot_med,plot_QTot_med,location_Charge,0);
		HistogramPlotterLDMX(hzz,100,PE_med1,PE_med2,filename_QTot_med2,plot_QTot_med2,location_Charge,0);
		HistogramPlotterLDMX(hzz,100,PE_med2,PE_high,filename_QTot_high,plot_QTot_high,location_Charge,0);
		HistogramPlotterLDMX(hzz,100,-10,400,filename_QTotal,plot_QTotal,location_Charge,1);
		HistogramPlotterLDMX(hzz,100,0,PE_low,filename_PElow,plot_PElow,location_Hits,1);
		HistogramPlotterLDMX(hzz,100,PE_low,PE_med1,filename_PEmed1,plot_PEmed1,location_Hits,1);
		HistogramPlotterLDMX(hzz,100,PE_med1,PE_med2,filename_PEmed2,plot_PEmed2,location_Hits,1);
		HistogramPlotterLDMX(hzz,100,PE_med2,PEmax,filename_PEhigh,plot_PEhigh,location_Hits,1);
	}

	for(unsigned int iB=0;iB<16;iB++){
		TString plot_ADC2D = "ADC/" + plot_title_ADCvsTS + iB;
		TString filename_ADC2D = plot_title_ADCvsTS + iB;
		TString plot_QlowvsTS = "Charge/" + plot_title_QlowvsTS + iB;
		TString filename_QlowvsTS = plot_title_QlowvsTS + iB;
		TString plot_QmedvsTS = "Charge/" + plot_title_QmedvsTS + iB;
		TString filename_QmedvsTS = plot_title_QmedvsTS + iB;
		TString plot_QhighvsTS = "Charge/" + plot_title_QhighvsTS + iB;
		TString filename_QhighvsTS = plot_title_QhighvsTS + iB; 
		HistogramPlotter2D(hzz,29,0,29,80,0,240,filename_ADC2D,plot_ADC2D,location_ADC,0);
		HistogramPlotter2D(hzz,29,0,29,nQbins_low,-100,Qlow_thr,filename_QlowvsTS,plot_QlowvsTS,location_Charge,3);
		HistogramPlotter2D(hzz,29,0,29,nQbins_high/15,Qmed_thr,Qhigh_thr,filename_QhighvsTS,plot_QhighvsTS,location_Charge,3);
		HistogramPlotter2D(hzz,29,0,29,nQbins_med/10,Qlow_thr,Qmed_thr,filename_QmedvsTS,plot_QmedvsTS,location_Charge,3);
	}

	
	//HistogramPlotterLDMX(hzz,10,60,70,"hTDC_Total","test/hTDC_Total");
}