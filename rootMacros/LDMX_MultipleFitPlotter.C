//Root file to plot distributions for integration windows for all channels
using namespace std;

void DrawFits(TFile *hzz1, TFile *hzz2, TFile *hzz3, TFile *hzz4, TFile *hzz5, TFile *hzz6, TFile *hzz7, TFile *hzz8, TFile *hzz9,Int_t bins, float_t xmin, float_t xmax, TString &plot_title, TString &hsig, TString &channel)
{
	TCanvas *c1 = new TCanvas();

	gStyle->SetLabelSize(0.048,"x");
	gStyle->SetLabelSize(0.048,"y");
	gStyle->SetTitleSize(0.046,"x");
	gStyle->SetTitleSize(0.046,"y");
	gStyle->SetOptFit(1111);
	gStyle->SetOptStat("nei");

	c1->Divide(3,3);
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
	TH1F *sig6 = new TH1F("Signal6",hsig,bins,xmin,xmax);
	sig6 = (TH1F*)hzz6->Get(hsig);
	TH1F *sig7 = new TH1F("Signal7",hsig,bins,xmin,xmax);
	sig7= (TH1F*)hzz7->Get(hsig);
	TH1F *sig8 = new TH1F("Signal8",hsig,bins,xmin,xmax);
	sig8 = (TH1F*)hzz8->Get(hsig);
	TH1F *sig9 = new TH1F("Signal9",hsig,bins,xmin,xmax);
	sig9 = (TH1F*)hzz9->Get(hsig);
	TH1F *sig[9] = {sig1,sig2,sig3,sig4,sig5,sig6,sig7,sig8,sig9};
	const char *labels[9] = {"10k events","20k events","30k events","40k events","50k events","60k events","70k events","80k events","90k events"};
	for(int g=0;g<9;g++)
	{
		c1->cd(g+1);
		sig[g]->SetLineWidth(2);
		sig[g]->SetLineColor(kBlue);
		sig[g]->SetTitle(plot_title);
		sig[g]->SetStats(11111111);
		float ymax = sig[g]->GetMaximum();

		TLatex *text = new TLatex(40,ymax/1.2,labels[g]);
		sig[g]->Draw();
		text->SetTextSize(0.05);
		text->Draw();
		//gPad->SetGrid(1,1);
      	gPad->Update();
	}
	TString outs = "/home/dhruvanshu/ldmx-sw/LDMX_Plots/MIPCalibs/MIPCalib_Stability_Check_Channel_" + channel + "_SameParams.pdf";
    c1->Print(outs);

}

void LDMX_MultipleFitPlotter()
{
	//fstream file;
	TFile *hzz1 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/unpacked_ldmx_captan_out_17-04-2022_21-14-11__179_reformat_30timeSamplesFrom0_linearize_hits_clusters_linearize_18_5_5_hits_10k_response_fin_SameParams.root");
	TFile *hzz2 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/unpacked_ldmx_captan_out_17-04-2022_21-14-11__179_reformat_30timeSamplesFrom0_linearize_hits_clusters_linearize_18_5_5_hits_20k_response_fin_SameParams.root");
	TFile *hzz3 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/unpacked_ldmx_captan_out_17-04-2022_21-14-11__179_reformat_30timeSamplesFrom0_linearize_hits_clusters_linearize_18_5_5_hits_30k_response_fin_SameParams.root");
	TFile *hzz4 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/unpacked_ldmx_captan_out_17-04-2022_21-14-11__179_reformat_30timeSamplesFrom0_linearize_hits_clusters_linearize_18_5_5_hits_40k_response_fin_SameParams.root");
	TFile *hzz5 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/unpacked_ldmx_captan_out_17-04-2022_21-14-11__179_reformat_30timeSamplesFrom0_linearize_hits_clusters_linearize_18_5_5_hits_50k_response_fin_SameParams.root");
	TFile *hzz6 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/unpacked_ldmx_captan_out_17-04-2022_21-14-11__179_reformat_30timeSamplesFrom0_linearize_hits_clusters_linearize_18_5_5_hits_60k_response_fin_SameParams.root");
	TFile *hzz7 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/unpacked_ldmx_captan_out_17-04-2022_21-14-11__179_reformat_30timeSamplesFrom0_linearize_hits_clusters_linearize_18_5_5_hits_70k_response_fin_SameParams.root");
	TFile *hzz8 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/unpacked_ldmx_captan_out_17-04-2022_21-14-11__179_reformat_30timeSamplesFrom0_linearize_hits_clusters_linearize_18_5_5_hits_80k_response_fin_SameParams.root");
	TFile *hzz9 = TFile::Open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/unpacked_ldmx_captan_out_17-04-2022_21-14-11__179_reformat_30timeSamplesFrom0_linearize_hits_clusters_linearize_18_5_5_hits_90k_response_fin_SameParams.root");


	TString plot_title_sample("hpe_chanID");
	//TString plot_title_2D("QvsSS_chan");
	//std::vector<int> SampleIndices;
	for(unsigned int i=0;i<12;i++) {
		//TString plot_title_chan = "Charge/" + plot_title_sample + i;
		TString plot_channel = plot_title_sample + i;
		TString channel_id = to_string(i);
		//TString plot_title_chan2D = "Charge/" + plot_title_2D + i;
		//TString plot_channel2D = plot_title_2D + i;
		DrawFits(hzz1,hzz2,hzz3,hzz4,hzz5,hzz6,hzz7,hzz8,hzz9,45,20,200,plot_channel,plot_channel,channel_id);
	}
}