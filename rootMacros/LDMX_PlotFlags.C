//Root file to plot distributions for integration windows for all channels
using namespace std;

void FlagIndex(TFile *hzz, Int_t bins, float_t xmin, float_t xmax, TString &plot_title, TString &hsig1, TString &hsig2, TString &hsig3, TString &hsig4, TString &hsig5, TString &channel)
{

	TCanvas *c1 = new TCanvas();

	gStyle->SetLabelSize(0.042,"x");
	gStyle->SetLabelSize(0.042,"y");
	gStyle->SetTitleSize(0.040,"x");
	gStyle->SetTitleSize(0.040,"y");

	TH1F *sig1 = new TH1F("Signal1",hsig1,bins,xmin,xmax);
	sig1 = (TH1F*)hzz->Get(hsig1);
	TH1F *sig2 = new TH1F("Signal2",hsig2,bins,xmin,xmax);
	sig2 = (TH1F*)hzz->Get(hsig2);
	TH1F *sig3 = new TH1F("Signal3",hsig3,bins,xmin,xmax);
	sig3 = (TH1F*)hzz->Get(hsig3);
	TH1F *sig4 = new TH1F("Signal4",hsig4,bins,xmin,xmax);
	sig4 = (TH1F*)hzz->Get(hsig4);
	TH1F *sig5 = new TH1F("Signal5",hsig5,bins,xmin,xmax);
	sig5 = (TH1F*)hzz->Get(hsig5);

	THStack *Win = new THStack("Win","Channel_" + channel + ";PEs; Number of Events");

	sig1->SetLineWidth(2); //All
	sig2->SetLineWidth(2); //Only 0
	sig3->SetLineWidth(2); //Only 4
	sig4->SetLineWidth(2); //Only 4 or 0
	sig5->SetLineWidth(2); //Neither 4 nor 0

	sig1->SetLineColor(kBlue);
	sig2->SetLineColor(kGreen);
	sig3->SetLineColor(kBlack);
	sig4->SetLineColor(kRed);
	sig5->SetLineColor(kMagenta);


	TLegend *legZ2 = new TLegend(0.6,0.6,0.9,0.9);
	//legZ2->SetTextSize(500);
	legZ2->AddEntry(sig1,"All Flags");
	legZ2->AddEntry(sig2,"Only Flag 0");
	legZ2->AddEntry(sig3,"Only Flag 4");
	legZ2->AddEntry(sig4,"Only Flag 0 or 4");
	legZ2->AddEntry(sig5,"Neither Flag 0 nor 4");
	legZ2->SetTextSize(0.035);

	Win->SetTitle(plot_title);	

	Win->Add(sig1);
	Win->Add(sig2);
	Win->Add(sig3);
	Win->Add(sig4);
	Win->Add(sig5);
	//c1->SetLogy();
	Win->Draw("HIST NOSTACK");
	Win->GetXaxis()->SetTitle("PE");
	Win->GetYaxis()->SetTitle("Number of events");
	Win->GetXaxis()->SetLabelSize(0.05);
  	Win->GetXaxis()->SetTitleSize(0.05);
  	Win->GetYaxis()->SetLabelSize(0.05);
  	Win->GetYaxis()->SetTitleSize(0.05);
  	Win->GetYaxis()->SetTitleOffset(0.9);
	Win->SetTitle("");
	c1->Modified();
	c1->Update();
  	TLatex *text = new TLatex(gPad->GetUxmax()-2.4,gPad->GetUymax()+1450,"#it{LDMX Internal}");
  	text->SetTextSize(0.05);
  	text->Draw();
	legZ2->Draw();
	TString outs = "/home/dhruvanshu/ldmx-sw/Barplots/PaperPlotsRevised/Channel_" + channel + "_MIPFlagData_60_150.pdf";
	c1->Print(outs);
}

void LDMX_PlotFlags()
{

	//TFile *hzz1 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/unpacked_ldmx_captan_out_17-04-2022_21-14-11__179_reformat_30timeSamplesFrom0_linearize_hits_clusters_linearize_18_5_5_hits_full_multiple_analysed_v2.root");
	TFile *hzz1 = TFile::Open("/home/dhruvanshu/ldmx-sw/Barplots/MIPFlags/unpacked_ldmx_captan_out_17-04-2022_21-14-11__179_reformat_30timeSamplesFrom0_linearize_hits_clusters_18_5_5_hits_analysed.root");
	//TFile *hzz1 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/CollaborationMeeting_Plots/unpacked_reprocessed_ldmx_captan_out_16-04-2022_19-22-54__160_digi_linearize_reprocessed_startsamples_pulsewidth_7_hits_analysed.root");

	int widthlimit = 10;
	int PE_low = 1;
    int PE_med1 = 60;
    int PE_med2 = 150;
    int nPE_bins = 100;

	TString plot_title_sample1("hPE_med2_bar");
	TString plot_title_sample2("hPE_med2_bar_flag0_");
	TString plot_title_sample3("hPE_med2_bar_flag4_");
	TString plot_title_sample4("hPE_med2_bar_flag04_");
	TString plot_title_sample5("hPE_med2_bar_flagN04_");
	for(unsigned int i=0;i<12;i++) {
		TString plot_title1_chan = "Hits/" + plot_title_sample1 + i;
		TString plot_title2_chan = "Hits/" + plot_title_sample2 + i;
		TString plot_title3_chan = "Hits/" + plot_title_sample3 + i;
		TString plot_title4_chan = "Hits/" + plot_title_sample4 + i;
		TString plot_title5_chan = "Hits/" + plot_title_sample5 + i;
		//TString plot_channel1 = plot_title_sample1 + i;
		//TString plot_channel2 = plot_title_sample2 + i;
		//TString plot_channel3 = plot_title_sample3 + i;
		//TString plot_channel4 = plot_title_sample4 + i;
		//TString plot_channel5 = plot_title_sample5 + i;
		TString channel_id = to_string(i);
		TString plot_title = "Channel_" + to_string(i) + " PE peaks data flag comparison";
		FlagIndex(hzz1,nPE_bins,PE_med1,PE_med2,plot_title,plot_title1_chan,plot_title2_chan,plot_title3_chan,plot_title4_chan,plot_title5_chan,channel_id);		
	}

}