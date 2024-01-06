using namespace std;
void LDMX_Histogram_Integrator(const char *sigfile, const char *inputfile) 
{
	TFile *hzz = TFile::Open(sigfile);
	double MIPs,Events;
	TString plot_title_MIP("hQTot_med_chan_MIPcut");
	TString plot_title_Total("hQTotal_channel");
	fstream file;
	file.open(inputfile,ios::out);
	for(unsigned int iB=0;iB<16;iB++)
	{
		
		TString MIPCuts = "Charge/" + plot_title_MIP + iB;
		TString TotalEvents = "Charge/" + plot_title_Total + iB;

		TH1F *sig = new TH1F("Signal",MIPCuts,100,50,400);
		sig = (TH1F*)hzz->Get(MIPCuts);
		TH1F *sig_total = new TH1F("Total Signal",TotalEvents,100,-10,400);
		sig_total = (TH1F*)hzz->Get(TotalEvents);
		
		MIPs = sig->Integral();
		Events = sig_total->Integral();

		file << iB << "," << MIPs << "," << Events << "\n" ;
		sig->Clear();
		sig_total->Clear();
	}
	file.close();
}