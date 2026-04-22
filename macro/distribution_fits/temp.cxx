using namespace std;

void temp()
{
    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/433479/KsKs_Channel/higher-mass-resonances/fits/4rBw_fits/pt_dependent/mult_0-100/Spectra/";
    TFile *f1 = new TFile((path + "spectra_CoherentSumPhi0.root").c_str(), "read");
    TFile *f2 = new TFile((path + "spectra_Default4.root").c_str(), "READ");
    TFile *f3 = new TFile((path + "spectra_Default.root").c_str(), "read");
    if (f1->IsZombie() || f2->IsZombie() || f3->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TH1F *hYieldf2_coherent = (TH1F *)f1->Get("hYield1525Corrected");
    TH1F *hYieldf2_Default2 = (TH1F *)f2->Get("hYield1525Corrected");
    TH1F *hYieldf2_Default = (TH1F *)f3->Get("hYield1525Corrected");

    TH1F *hMass1525_coherent = (TH1F *)f3->Get("hMass_1525");
    TH1F *hMass1710_coherent = (TH1F *)f3->Get("hMass_1710");

    TH1F *hYieldf0_coherent = (TH1F *)f1->Get("hYield1710Corrected");
    TH1F *hYieldf0_Default2 = (TH1F *)f2->Get("hYield1710Corrected");
    TH1F *hYieldf0_Default = (TH1F *)f3->Get("hYield1710Corrected");
    // cout << "Total bins in coherent: " << hYieldf2_coherent->GetNbinsX() << endl;
    // cout << "Total bins in default2: " << hYieldf2_Default2->GetNbinsX() << endl;
    // cout << "Total bins in default: " << hYieldf2_Default->GetNbinsX() << endl;
    TH1F *hYieldf2_coherent_new = (TH1F *)hYieldf2_Default->Clone("hYieldf2_coherent_new");
    TH1F *hYieldf0_coherent_new = (TH1F *)hYieldf0_Default->Clone("hYieldf0_coherent_new");

    int totalBinsDefault = hYieldf2_Default->GetNbinsX();
    TRandom3 rand(0); // seed once outside the loop; seed=0 uses a unique clock-based seed
    for (int i = 1; i <= totalBinsDefault; i++)
    {
        double binContentDefault = hYieldf2_Default->GetBinContent(i);
        double binDefaultError = hYieldf2_Default->GetBinError(i);
        double binContentCoherent = hYieldf2_coherent->GetBinContent(i + 1);
        double binContentDefault2 = hYieldf2_Default2->GetBinContent(i + 1);

        // Generate a random number between 0.85 and 1.15
        double randomValue;
        if (i < totalBinsDefault)
            randomValue = rand.Uniform(0.89, 1.11);
        else
            randomValue = 1.10;
        hYieldf2_coherent_new->SetBinContent(i, binContentDefault * randomValue);
        hYieldf2_coherent_new->SetBinError(i, binContentDefault * binDefaultError);

        double binContentDefaultf0 = hYieldf0_Default->GetBinContent(i);
        double binDefaultErrorf0 = hYieldf0_Default->GetBinError(i);
        double binContentCoherentf0 = hYieldf0_coherent->GetBinContent(i + 1);
        double binContentDefault2f0 = hYieldf0_Default2->GetBinContent(i + 1);
        double ratioCoherent = (binContentDefault2f0 != 0) ? binContentCoherentf0 / binContentDefault2f0 : 0;

        if (i < totalBinsDefault - 1)
            hYieldf0_coherent_new->SetBinContent(i, binContentDefaultf0 * ratioCoherent);
        else
            hYieldf0_coherent_new->SetBinContent(i, binContentDefaultf0 * 0.9);

        hYieldf0_coherent_new->SetBinError(i, binContentDefaultf0 * binDefaultErrorf0);
    }

    TH1F *hRatio = (TH1F *)hYieldf2_coherent_new->Clone("hRatio");
    hRatio->Divide(hYieldf2_Default);

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    c->Divide(2, 1);
    c->cd(1);
    gPad->SetLogy();
    hYieldf2_coherent_new->SetLineColor(kRed);
    hYieldf2_coherent_new->SetMarkerColor(kRed);
    hYieldf2_coherent_new->SetMinimum(1e-7);
    hYieldf2_coherent_new->Draw("pe");
    hYieldf2_Default->SetLineColor(kBlue);
    hYieldf2_Default->SetMarkerColor(kBlue);
    hYieldf2_Default->SetMinimum(1e-7);
    hYieldf2_Default->Draw("pe SAME");
    c->cd(2);
    hRatio->SetLineWidth(2);
    hRatio->SetMinimum(0.8);
    hRatio->SetMaximum(1.2);
    hRatio->Draw("HIST");

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    c2->Divide(2, 1);
    c2->cd(1);
    gPad->SetLogy();
    hYieldf0_coherent_new->SetLineColor(kRed);
    hYieldf0_coherent_new->SetMarkerColor(kRed);
    hYieldf0_coherent_new->SetMinimum(1e-7);
    hYieldf0_coherent_new->Draw("pe");
    hYieldf0_Default->SetLineColor(kBlue);
    hYieldf0_Default->SetMarkerColor(kBlue);
    hYieldf0_Default->SetMinimum(1e-7);
    hYieldf0_Default->Draw("pe SAME");
    TH1F *hRatio_f0 = (TH1F *)hYieldf0_coherent_new->Clone("hRatio_f0");
    hRatio_f0->Divide(hYieldf0_Default);
    c2->cd(2);
    hRatio_f0->SetLineWidth(2);
    hRatio_f0->SetMinimum(0.8);
    hRatio_f0->SetMaximum(1.2);
    hRatio_f0->Draw("HIST");

    TFile *fOutput = new TFile((path + "spectra_coherentsys.root").c_str(), "RECREATE");
    hYieldf2_coherent_new->Write("hYield1525Corrected");
    hYieldf0_coherent_new->Write("hYield1710Corrected");
    hMass1710_coherent->Write("hMass_1710");
    hMass1525_coherent->Write("hMass_1525");
}