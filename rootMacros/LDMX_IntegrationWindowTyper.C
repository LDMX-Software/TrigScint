//Root file to plot distributions for integration windows for all channels
using namespace std;

int IntegrationWindowIndex(TFile *hzz, Int_t bins, float_t xmin, float_t xmax, TString &plot_title, TString &hsig, TString &location)
{
	TCanvas *c1 = new TCanvas();

	gStyle->SetLabelSize(0.042,"x");
	gStyle->SetLabelSize(0.042,"y");
	gStyle->SetTitleSize(0.040,"x");
	gStyle->SetTitleSize(0.040,"y");

	TH1F *sig = new TH1F("Signal",hsig,bins,xmin,xmax);
	sig = (TH1F*)hzz->Get(hsig);

	sig->SetLineWidth(2);
	sig->SetLineColor(kBlue);
	sig->SetTitle(plot_title);
	sig->SetStats(111111);

	int binmax = sig->GetMaximumBin();
	float index = sig->GetXaxis()->GetBinCenter(binmax);
	Float_t ymax = sig->GetMaximum();
	TLatex *text = new TLatex(19,ymax/2.0,"PulseWidth = 7");

	sig->Draw();
	text->SetTextSize(0.04);
	text->Draw();
	//TString outs = "/home/dhruvanshu/ldmx-sw/LDMX_Plots/LDMX_Analysis_Run5/" + location + "/" + plot_title + ".png";
	//c1->Print(outs);	

	//sig->clear();
	return int(index);
}


void LDMX_IntegrationWindowTyper(const char *sigfile, const char *identifier)
{
	fstream file;
	TFile *hzz = TFile::Open(sigfile);
	TString location = identifier;
	//file.open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/unpacked_ldmx_captan_out_" + location + "_reformat_30timeSamplesFrom0_linearize_hits_clusters_width.txt",ios::out);
	//file.open("/home/dhruvanshu/LDMX_Analysis_Files/LDMX_Analysis_Run1/instrumentation_run287_v3_rotated180_digi_linearize_width.txt",ios::out);
	file.open(location + ".txt",ios::out);
	int widthlimit = 10;
	float totChargelimit = 10000.0;    	
	float highchargelimit = 8*totChargelimit;
	TString plot_title_sample("ChannelWidthData_chan");
	TString plot_title_2D("QvsSS_chan");
	std::vector<int> SampleIndices;
	for(unsigned int i=0;i<12;i++) {
		TString plot_title_chan = "Charge/" + plot_title_sample + i;
		TString plot_channel = plot_title_sample + i;
		TString plot_title_chan2D = "Charge/" + plot_title_2D + i;
		TString plot_channel2D = plot_title_2D + i;
		int StartSample = IntegrationWindowIndex(hzz,widthlimit,12,12+widthlimit,plot_channel,plot_title_chan,location);
		file << i << "," << StartSample << "\n" ;
		std::cout << " Peak Location found at " << StartSample << " for channel " << i << std::endl;
		SampleIndices.push_back(StartSample);
	}
	TString plot_title_all("StartSampleAllChannels");
	file.close();
}