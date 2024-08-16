#include <iostream>
using namespace std;
#include "src/style.h"

double betheBloch(double *x, double *par)
{
    double p = x[0];   // particle momentum
    double m = par[0]; // particle mass
    double Z = par[1]; // particle charge
    double I = par[2]; // mean excitation potential (in MeV)

    double K = 0.307075; // constant in MeV cm^2/mol

    double beta = p / sqrt(p * p + m * m);
    double gamma = 1.0 / sqrt(1 - beta * beta);
    double bg = beta * gamma;

    double Tmax = (2 * 0.511e-3 * bg * bg) / (1 + 2 * bg * sqrt(1 + bg * bg) + 0.511e-3 / m); // Tmax in GeV
    double delta = 0.0;                                                                       // Density effect correction term, you can model this as well if needed

    double dEdx = K * Z * Z / (beta * beta) * (0.5 * log(2 * 0.511e-3 * bg * bg * Tmax / (I * I)) - beta * beta - delta / 2);
    return dEdx; // dE/dx in MeV/cm
}

// ALEPH parameterization for TPC energy loss
double alephParam(double *x, double *par)
{
    double particlemass = 0.938;     // Proton mass in GeV/c^2
    double bg = x[0] / particlemass; // beta*gamma (beta * gamma = p/mc)
    double P1 = par[0];
    double P2 = par[1];
    double P3 = par[2];
    double P4 = par[3];
    double P5 = par[4];

    double dEdx = (P1 / pow(bg, P4)) * (P2 - pow(bg, 4) - log(P3 + 1.0 / pow(bg, P5)));
    return dEdx; // dE/dx in arbitrary units
}

void temp_code()
{
    gStyle->SetOptStat(0);
    // gStyle->SetOptStat(1110);
    TFile *f = new TFile("../data/glueball/LHC22o_pass7_small/250337.root", "READ");
    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }

    string path = "strangeness_tutorial/kzeroShort/NksProduced";
    TH1F *h = (TH1F *)f->Get(path.c_str());
    if (h == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    // TCanvas *c = new TCanvas("c", "c", 720, 720);
    // SetCanvasStyle(c, 0.15, 0.03, 0.03, 0.14);
    // SetHistoQA(h);
    // gPad->SetLogy();
    // h->GetXaxis()->SetTitle("N_{K_{S}/Event}");
    // h->GetYaxis()->SetTitle("Counts");
    // h->SetMaximum(1e11);
    // h->Draw("hist text");
    // c->SaveAs("distribution_fits/saved/NksProduced.png");

    // // cout<<"Fraction of events in which more than 1 Ks are produced is "<<h->Integral(h->FindBin(2), h->GetNbinsX())/h->Integral()<<endl;

    // TCanvas *c_tpc_energyloss = new TCanvas("c_tpc_energyloss", "c_tpc_energyloss", 720, 720);
    // SetCanvasStyle(c_tpc_energyloss, 0.15, 0.13, 0.03, 0.14);
    // string path_tpc_energyloss = "strangeness_tutorial/kzeroShort/dE_by_dx_TPC";
    // TH2F *h_tpc_energyloss = (TH2F *)f->Get(path_tpc_energyloss.c_str());
    // if (h_tpc_energyloss == nullptr)
    // {
    //     cout << "Error reading histogram" << endl;
    //     return;
    // }
    // SetHistoQA(h_tpc_energyloss);
    // gPad->SetLogz();
    // gPad->SetLogx();
    // h_tpc_energyloss->GetXaxis()->SetRangeUser(0.1, 100);
    // h_tpc_energyloss->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // h_tpc_energyloss->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
    // h_tpc_energyloss->GetXaxis()->SetTitleOffset(1.3);
    // h_tpc_energyloss->Draw("colz");
    //  // c_tpc_energyloss->SaveAs("distribution_fits/saved/dE_by_dx_TPC.png");

    // // // lets fit the dE/dx distribution
    // // //  Parameters: mass, charge, mean excitation potential
    // double protonMass = 0.938; // GeV/c^2
    // double charge = 1.0;       // e
    // double I = 0.0000136;      // MeV (approx for TPC gas)

    // TF1 *betheBlochFunc = new TF1("betheBlochFunc", betheBloch, 0.8, 1.0, 3);
    // betheBlochFunc->SetParameter(0, protonMass);
    // betheBlochFunc->SetParameter(1, charge);
    // betheBlochFunc->SetParameter(2, I);
    // betheBlochFunc->SetTitle("Bethe-Bloch Curve;Momentum p (GeV/c);dE/dx (MeV/cm)");
    // betheBlochFunc->SetLineColor(kRed);
    // betheBlochFunc->SetLineWidth(3);

    // // Draw the function
    // h_tpc_energyloss->Fit("betheBlochFunc", "R");
    // betheBlochFunc->Draw("same");

    // // // ALEPH parameterization for TPC energy loss
    // // Set parameter values
    // double P1 = 2.0; // Example value
    // double P2 = 0.2; // Example value
    // double P3 = 3.0; // Example value
    // double P4 = 1.0; // Example value
    // double P5 = 2.0; // Example value

    // TF1 *alephFunc = new TF1("alephFunc", alephParam, 0.5, 1, 5);
    // alephFunc->SetParameters(P1, P2, P3, P4, P5);
    // alephFunc->SetTitle("ALEPH Parameterization;#beta#gamma;dE/dx (Arbitrary Units)");
    // alephFunc->SetLineColor(kRed);
    // alephFunc->SetLineWidth(3);
    // h_tpc_energyloss->Fit("alephFunc", "R");
    // Draw the function
    // alephFunc->Draw();

    // string path_pt_pos = "strangeness_tutorial/kzeroShort/positive_pt";
    // string path_pt_neg = "strangeness_tutorial/kzeroShort/negative_pt";
    // string path_phi_pos = "strangeness_tutorial/kzeroShort/positive_phi";
    // string path_phi_neg = "strangeness_tutorial/kzeroShort/negative_phi";

    // TH1F *h_pt_pos = (TH1F *)f->Get(path_pt_pos.c_str());
    // TH1F *h_pt_neg = (TH1F *)f->Get(path_pt_neg.c_str());
    // TH1F *h_phi_pos = (TH1F *)f->Get(path_phi_pos.c_str());
    // TH1F *h_phi_neg = (TH1F *)f->Get(path_phi_neg.c_str());

    // if (h_pt_pos == nullptr || h_pt_neg == nullptr || h_phi_pos == nullptr || h_phi_neg == nullptr)
    // {
    //     cout << "Error reading histogram" << endl;
    //     return;
    // }

    // TCanvas *c_pt = new TCanvas("c_pt", "c_pt", 720, 720);
    // SetCanvasStyle(c_pt, 0.15, 0.03, 0.06, 0.14);
    // gPad->SetLogy();
    // SetHistoQA(h_pt_pos);
    // SetHistoQA(h_pt_neg);
    // h_pt_pos->SetLineColor(kRed);
    // h_pt_neg->SetLineColor(kBlue);
    // h_pt_pos->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // h_pt_pos->GetYaxis()->SetTitle("Counts");
    // h_pt_pos->GetYaxis()->SetMaxDigits(3);
    // h_pt_pos->SetMaximum(h_pt_pos->GetMaximum() * 5.5);
    // h_pt_pos->Draw("hist");
    // h_pt_neg->Draw("hist same");

    // TLegend *leg = new TLegend(0.6, 0.7, 0.8, 0.9);
    // leg->AddEntry(h_pt_pos, "#pi^{+}", "l");
    // leg->AddEntry(h_pt_neg, "#pi^{-}", "l");
    // leg->Draw();

    // c_pt->SaveAs("distribution_fits/saved/pt_new.png");

    // TCanvas *c_phi = new TCanvas("c_phi", "c_phi", 720, 720);
    // gPad->SetLogy(0);
    // SetCanvasStyle(c_phi, 0.15, 0.03, 0.06, 0.14);
    // SetHistoQA(h_phi_pos);
    // SetHistoQA(h_phi_neg);
    // h_phi_pos->SetLineColor(kRed);
    // h_phi_neg->SetLineColor(kBlue);
    // h_phi_pos->GetXaxis()->SetTitle("#phi");
    // h_phi_pos->GetYaxis()->SetTitle("Counts");
    // h_phi_pos->GetYaxis()->SetMaxDigits(3);
    // h_phi_pos->Draw("hist");
    // h_phi_neg->Draw("hist same");

    // TLegend *leg_phi = new TLegend(0.4, 0.3, 0.6, 0.5);
    // leg_phi->AddEntry(h_phi_pos, "#pi^{+}", "l");
    // leg_phi->AddEntry(h_phi_neg, "#pi^{-}", "l");
    // leg_phi->Draw();

    // c_phi->SaveAs("distribution_fits/saved/phi_new.png");

    // TCanvas *cmultdist = new TCanvas("cmultdist", "cmultdist", 720, 720);
    // gPad->SetLogy();
    // string path_mult = "strangeness_tutorial/eventSelection/multdist_FT0M";
    // TH1F *h_mult = (TH1F *)f->Get(path_mult.c_str());
    // if (h_mult == nullptr)
    // {
    //     cout << "Error reading histogram" << endl;
    //     return;
    // }
    // SetCanvasStyle(cmultdist, 0.15, 0.03, 0.03, 0.14);
    // SetHistoQA(h_mult);
    // h_mult->GetXaxis()->SetTitle("Multiplicity");
    // h_mult->GetYaxis()->SetTitle("Counts");
    // h_mult->Draw("hist");

    // TF1 *nbd = new TF1("nbd", "([0]*TMath::Gamma(x[0]+[1])/(TMath::Gamma(x[0]+1)*TMath::Gamma([1])))*TMath::Power([2], [1])*TMath::Power(1-[2], x[0])", 0, 60000);
    // nbd->SetParameters(1e5, 1, 0.5); // Initial values: par0 (normalization), par1 (no of success), par2 (success probability)
    // h_mult->Fit("nbd", "RS");

    // Correlation plots
    TCanvas *c_correlation = new TCanvas("c_correlation", "c_correlation", 720, 720);
    string path_corr_before = "strangeness_tutorial/kzeroShort/mass_lambda_kshort_before";
    string path_corr_after = "strangeness_tutorial/kzeroShort/mass_lambda_kshort_after";
    TH2F *h_corr_before = (TH2F *)f->Get(path_corr_before.c_str());
    TH2F *h_corr_after = (TH2F *)f->Get(path_corr_after.c_str());
    if (h_corr_before == nullptr || h_corr_after == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }

    SetCanvasStyle(c_correlation, 0.15, 0.12, 0.06, 0.14);
    SetHistoQA(h_corr_before);
    h_corr_before->GetYaxis()->SetTitle("M_{#Lambda} (GeV/c^{2})");
    h_corr_before->GetXaxis()->SetTitle("M_{K_{S}} (GeV/c^{2})");
    h_corr_before->GetYaxis()->SetRangeUser(1, 1.5);
    h_corr_before->Draw("colz");
    c_correlation->SaveAs("distribution_fits/saved/mass_lambda_kshort_before.png");

    TCanvas *c_correlation_after = new TCanvas("c_correlation_after", "c_correlation_after", 720, 720);
    SetCanvasStyle(c_correlation_after, 0.15, 0.12, 0.06, 0.14);
    SetHistoQA(h_corr_after);
    h_corr_after->GetYaxis()->SetTitle("M_{#Lambda} (GeV/c^{2})");
    h_corr_after->GetXaxis()->SetTitle("M_{K_{S}} (GeV/c^{2})");
    h_corr_after->GetYaxis()->SetRangeUser(1, 1.5);
    h_corr_after->Draw("colz");
    c_correlation_after->SaveAs("distribution_fits/saved/mass_lambda_kshort_after.png");

}