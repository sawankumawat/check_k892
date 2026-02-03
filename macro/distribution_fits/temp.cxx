using namespace std;

void temp(){
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/systematic2022_new/KsKs_Channel/higher-mass-resonances";
    TFile *f1 = new TFile((path + "/mult_0-100/Spectra/spectra_.root").c_str(), "READ");
    TFile *f2 = new TFile((path + "/fits/FitParam.root").c_str(), "READ");
    TH1F *h1 = (TH1F *)f1->Get("hMass1525");
    TH1F *h2 = (TH1F *)f2->Get("Mult_0_100/hMass_1525");
    cout<<"Number of bins in h1: "<<h1->GetNbinsX()<<endl;
    cout<<"Number of bins in h2: "<<h2->GetNbinsX()<<endl;
    for(int ibin = 1; ibin <= h1->GetNbinsX(); ibin++){
        cout<<"Bin "<<ibin<<": h1 = "<<h1->GetBinContent(ibin)<<" , h2 = "<<h2->GetBinContent(ibin)<<endl;
    }
}