//Root file to plot distributions for integration windows for all channels
using namespace std;

//void IntegrationWindowIndex(TFile *hzz1, TFile *hzz2, TFile *hzz3 ,TFile *hzz4 ,TFile *hzz5 ,TFile *hzz6 ,TFile *hzz7 ,TFile *hzz8, TFile *hzz9, TFile *hzz10, Int_t bins, float_t xmin, float_t xmax, TString &plot_title, TString &hsig, TString &channel)
void IntegrationWindowIndex(TFile *hzz1, TFile *hzz2, TFile *hzz3 ,TFile *hzz4 ,TFile *hzz5, Int_t bins, float_t xmin, float_t xmax, TString &plot_title, TString &hsig, TString &channel)
{

	TCanvas *c1 = new TCanvas();

	gStyle->SetLabelSize(0.042,"x");
	gStyle->SetLabelSize(0.042,"y");
	gStyle->SetTitleSize(0.040,"x");
	gStyle->SetTitleSize(0.040,"y");

	TH1F *h1 = (TH1F*)gROOT->FindObject("Signal1");
  	delete h1;
  	TH1F *h2 = (TH1F*)gROOT->FindObject("Signal2");
  	delete h2;
  	TH1F *h3 = (TH1F*)gROOT->FindObject("Signal3");
  	delete h3;
  	TH1F *h4 = (TH1F*)gROOT->FindObject("Signal4");
  	delete h4;
  	TH1F *h5 = (TH1F*)gROOT->FindObject("Signal5");
  	delete h5;
  	TH1F *h6 = (TH1F*)gROOT->FindObject("Signal6");
  	delete h6;
	TH1F *sig1 = new TH1F("Signal1",hsig,bins,xmin,xmax);
	sig1 = (TH1F*)hzz1->Get(hsig);
	TH1F *sig2 = new TH1F("Signal2",hsig,bins,xmin,xmax);
	sig2 = (TH1F*)hzz2->Get(hsig);
	TH1F *sig3 = new TH1F("Signal3",hsig,bins,xmin,xmax);
	sig3 = (TH1F*)hzz3->Get(hsig);
	TH1F *sig4 = new TH1F("Signal4",hsig,bins,xmin,xmax);
	sig4 = (TH1F*)hzz4->Get(hsig);
	TH1F *sig5 = new TH1F("Signal5",hsig,bins,xmin,xmax);
	sig5 = (TH1F*)hzz5->Get(hsig);
	//TH1F *sig6 = new TH1F("Signal6",hsig,bins,xmin,xmax);
	//sig6 = (TH1F*)hzz6->Get(hsig);
	//TH1F *sig7 = new TH1F("Signal7",hsig,bins,xmin,xmax);
	//sig7= (TH1F*)hzz7->Get(hsig);
	//TH1F *sig8 = new TH1F("Signal8",hsig,bins,xmin,xmax);
	//sig8 = (TH1F*)hzz8->Get(hsig);
	//TH1F *sig9 = new TH1F("Signal9",hsig,bins,xmin,xmax);
	//sig9 = (TH1F*)hzz9->Get(hsig);
	//TH1F *sig10 = new TH1F("Signal10",hsig,bins,xmin,xmax);
	//sig10 = (TH1F*)hzz10->Get(hsig);

	THStack *Win = new THStack("Win","Channel_" + channel + ";PE; No of Events");
	//TString title = "Channel_" + channel;

	sig1->SetLineWidth(2);
	sig2->SetLineWidth(2);
	sig3->SetLineWidth(2);
	sig4->SetLineWidth(2);
	sig5->SetLineWidth(2);
	//sig6->SetLineWidth(2);
	//sig7->SetLineWidth(2);
	//sig8->SetLineWidth(2);
	//sig9->SetLineWidth(2);
	//sig10->SetLineWidth(2);

	sig1->SetLineColor(kBlue);
	sig2->SetLineColor(kBlack);
	sig5->SetLineColor(kMagenta);
	sig3->SetLineColor(kRed);
	sig4->SetLineColor(kGreen);
	//sig6->SetLineColor(kYellow);
	//sig7->SetLineColor(kGreen);
	//sig8->SetLineColor(kCyan);
	//sig9->SetLineColor(kBlack);
	//sig10->SetLineColor(kMagenta);

	TLegend *legZ2 = new TLegend(0.7,0.7,0.9,0.9);
	//legZ2->SetTextSize(500);
	legZ2->AddEntry(sig1,"0-4");
	legZ2->AddEntry(sig2,"5-9");
	legZ2->AddEntry(sig3,"15-19");
	legZ2->AddEntry(sig4,"20-24");
	legZ2->AddEntry(sig5,"25-29");
	//legZ2->AddEntry(sig6,"20-23");
	//legZ2->AddEntry(sig7,"17-04-2022_21-14-11__179");
	//legZ2->AddEntry(sig8,"12-04-2022_21-52-34__182");
	//legZ2->AddEntry(sig9,"19-04-2022_11-26-03__205");
	//legZ2->AddEntry(sig10,"20-04-2022_20-40-43__217");
	Win->SetTitle(plot_title);	

	Win->Add(sig1);
	Win->Add(sig2);
	Win->Add(sig3);
	Win->Add(sig4);
	Win->Add(sig5);
	//Win->Add(sig6);
	//Win->Add(sig7);
	//Win->Add(sig8);
	//Win->Add(sig9);
	//Win->Add(sig10);
	Win->Draw("HIST NOSTACK");
	Win->GetXaxis()->SetTitle("PE");
	Win->GetYaxis()->SetTitle("No of events");
	c1->Modified();
	c1->SetLogy();
	legZ2->Draw();
	TString outs = "/home/dhruvanshu/LDMX_Analysis_Files/NoiseStudyCalculations/IntegWindows/Width5/Updated/Channel_" + channel + "_SinglePEPeaksWindowData_Width5.pdf";
	c1->Print(outs);
}

void LDMX_PlotAllWindow()
{

	TFile *hzz1 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/NoiseStudyCalculations/IntegWindows/Width5/unpacked_reprocessed_ldmx_captan_out_16-04-2022_19-22-54__160_digi_linearize_reprocessed_startsamples_pulsewidth_5_hits_gain_corrected_TS_0_4_analysed.root");
	TFile *hzz2 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/NoiseStudyCalculations/IntegWindows/Width5/unpacked_reprocessed_ldmx_captan_out_16-04-2022_19-22-54__160_digi_linearize_reprocessed_startsamples_pulsewidth_5_hits_gain_corrected_TS_5_9_analysed.root");
	//TFile *hzz3 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/NoiseStudyCalculations/IntegWindows/Width5/unpacked_reprocessed_ldmx_captan_out_16-04-2022_19-22-54__160_digi_linearize_reprocessed_startsamples_pulsewidth_5_hits_gain_corrected_TS_10_14_analysed.root");
	TFile *hzz3 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/NoiseStudyCalculations/IntegWindows/Width5/unpacked_reprocessed_ldmx_captan_out_16-04-2022_19-22-54__160_digi_linearize_reprocessed_startsamples_pulsewidth_5_hits_gain_corrected_TS_15_19_analysed.root");
	TFile *hzz4 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/NoiseStudyCalculations/IntegWindows/Width5/unpacked_reprocessed_ldmx_captan_out_16-04-2022_19-22-54__160_digi_linearize_reprocessed_startsamples_pulsewidth_5_hits_gain_corrected_TS_20_24_analysed.root");
	TFile *hzz5 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/NoiseStudyCalculations/IntegWindows/Width5/unpacked_reprocessed_ldmx_captan_out_16-04-2022_19-22-54__160_digi_linearize_reprocessed_startsamples_pulsewidth_5_hits_gain_corrected_TS_25_29_analysed.root");
	//TFile *hzz6 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/NoiseStudyCalculations/IntegWindows/Width5/unpacked_reprocessed_ldmx_captan_out_16-04-2022_19-22-54__160_digi_linearize_reprocessed_startsamples_pulsewidth_4_hits_gain_corrected_TS_20_23_analysed.root");
	//TFile *hzz6 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/unpacked_ldmx_captan_out_17-04-2022_21-14-11__175_reformat_30timeSamplesFrom0_linearize_hits_clusters_width.root");
	//TFile *hzz7 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/unpacked_ldmx_captan_out_17-04-2022_21-14-11__179_reformat_30timeSamplesFrom0_linearize_hits_clusters_width.root");
	//TFile *hzz8 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/unpacked_ldmx_captan_out_12-04-2022_21-52-34__182_reformat_30timeSamplesFrom0_linearize_hits_clusters_width.root");
	//TFile *hzz9 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/unpacked_ldmx_captan_out_19-04-2022_11-26-03__205_reformat_30timeSamplesFrom0_linearize_hits_clusters_width.root");
	//TFile *hzz10 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/unpacked_ldmx_captan_out_20-04-2022_20-40-43__217_reformat_30timeSamplesFrom0_linearize_hits_clusters_width.root");

	//int widthlimit = 10;
	TString Peaks[12] = {" (PeakTS=17,Width=5)"," (PeakTS=17,Width=5)"," (PeakTS=17,Width=5)"," (PeakTS=17,Width=5)"," (PeakTS=17,Width=5)"," (PeakTS=17,Width=5)"," (PeakTS=15,Width=5)"," (PeakTS=15,Width=5)"," (PeakTS=16,Width=5)"," (PeakTS=15,Width=5)"," (PeakTS=17,Width=5)"," (PeakTS=17,Width=5)"};
	//TString Peaks[12] = {" (PeakTS=17)"," (PeakTS=17)"," (PeakTS=17)"," (PeakTS=17)"," (PeakTS=17)"," (PeakTS=17)"," (PeakTS=15)"," (PeakTS=15)"," (PeakTS=16)"," (PeakTS=15)"," (PeakTS=17)"," (PeakTS=17)"};
	TString plot_title_sample("hPE_med1_bar_flag04_");
	for(unsigned int i=0;i<12;i++) {
		TString plot_title_chan = "Hits/" + plot_title_sample + i;
		TString plot_channel = "Channel_";
		plot_channel += to_string(i);
		plot_channel += Peaks[i];
		TString channel_id = to_string(i);
		IntegrationWindowIndex(hzz1,hzz2,hzz3,hzz4,hzz5,200,0,10,plot_channel,plot_title_chan,channel_id);		
	}

}