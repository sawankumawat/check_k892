#include <iostream>
#include "TDatabasePDG.h"
#include "../src/fitfunc.h"
#include "../src/common_glue.h"
#include "../src/style.h"

TDatabasePDG *pdg = new TDatabasePDG();
TF1 *fitgaus(TH1 *h, double ksmass, double kswidth);
TF1 *fitgauspol2(TH1 *h, double ksmass, double kswidth);
TF1 *CBpol2(TH1 *h, double *parameters);
TF1 *CB(TH1 *h, double ksmass, double kswidth);
void SetHistoStyle_temp(TH1 *h, Int_t MCol, Int_t MSty, double binwidth);

void gaussian_fit_Ks()
{
    double ksmass = pdg->GetParticle(310)->Mass();
    // double kswidth = pdg->GetParticle(310)->Width();
    double kswidth = 0.005;
    cout << "PDG mass: " << ksmass << " PDG width: " << kswidth << endl;
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);
    TLatex *t2 = new TLatex();
    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.045);
    t2->SetTextFont(42);

    //********************************************************************

    const string outputQAfolder_str = kSignalOutput + "/" + kchannel + "/" + kfoldername + "/QA"; // path for root file
    TFile *f = new TFile(("../" + outputQAfolder_str + "/KsInvMass.root").c_str(), "read");
    THnSparseF *hsparse = (THnSparseF *)f->Get("kshort_2dsparse");
    if (hsparse == nullptr)
    {
        cout << "THnSparse not found" << endl;
        return;
    }
    THnSparseF *hsparseClone = (THnSparseF *)hsparse->Clone("hsparseClone");
    TH1F *hpt = (TH1F *)hsparse->Projection(1);
    int ptbinlow = hsparseClone->GetAxis(1)->FindBin(0.0);
    int ptbinhigh = hsparseClone->GetAxis(1)->FindBin(30.0);
    hsparseClone->GetAxis(1)->SetRange(ptbinlow, ptbinhigh);

    TH1F *hInvMass = (TH1F *)hsparseClone->Projection(0);
    double binwidth = hInvMass->GetBinWidth(1);
    cout << "Bin width: " << binwidth * 1000 << endl;

    TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
    SetCanvasStyle(c1, 0.145, 0.05, 0.05, 0.13);
    SetHistoStyle_temp(hInvMass, 1, 20, binwidth);
    hInvMass->Draw("pe");

    TH1F *hInvMassClone1 = (TH1F *)hInvMass->Clone("hInvMassClone1");
    TH1F *hInvMassClone2 = (TH1F *)hInvMass->Clone("hInvMassClone2");
    hInvMassClone1->SetLineColor(0);
    hInvMassClone1->SetLineWidth(0);
    hInvMassClone1->SetFillColor(5);
    hInvMassClone1->GetXaxis()->SetRangeUser(ksmass - 3 * kswidth, ksmass + 3 * kswidth);
    hInvMassClone1->Draw("E3 hist same");
    hInvMassClone2->Draw("pe same");
    // TF1 *fit = fitgaus(hInvMass, ksmass, kswidth);
    // TF1 *fit = fitgauspol2(hInvMass, ksmass, kswidth);
    // TF1 *fit = CBpol2(hInvMass, ksmass, kswidth);
    TF1 *fit = CB(hInvMass, ksmass, kswidth);
    double parameters[5];
    for (int i = 0; i < 5; i++)
    {
        parameters[i] = fit->GetParameter(i);
    }
    TF1 *fit2 = CBpol2(hInvMass, parameters);

    TLegend *lp2 = DrawLegend(0.14, 0.7, 0.4, 0.9);
    lp2->SetTextSize(0.04);
    lp2->SetTextFont(42);
    lp2->SetFillStyle(0);
    lp2->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    lp2->AddEntry((TObject *)0, "FT0M, 0-100%", "");
    lp2->AddEntry((TObject *)0, "|#it{y}| < 0.5", "");
    lp2->Draw("same");

    // rotational
    TLegend *lp3 = DrawLegend(0.65, 0.3, 0.9, 0.5);
    lp3->SetFillStyle(0);
    lp3->SetTextFont(42);
    lp3->AddEntry(hInvMass, "Data", "pe");
    lp3->AddEntry(fit2, "Gaussian Fit", "l");
    lp3->AddEntry(hInvMassClone1, "Signal", "f");
    lp3->SetTextSize(0.04);
    lp3->Draw("same");

    gPad->Update();
    TPaveStats *st = (TPaveStats *)hInvMass->FindObject("stats");
    st->SetX1NDC(0.57);
    st->SetX2NDC(0.95);
    st->SetY1NDC(0.53);
    st->SetY2NDC(0.95);
    st->Draw("same");

    // // Now we will plot the Ks invariant mass distribution as a function of pT
    // const int Nptbins = 12;
    // float ptbins[Nptbins + 1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 12.0, 20.0, 30.0};
    // TH1F *hInvMassPt[Nptbins];
    // TCanvas *c2 = new TCanvas("c2", "c2", 720, 720);
    // c2->Divide(4, 3);
    // for (int ipt = 0; ipt < Nptbins; ipt++)
    // {
    //     c2->cd(ipt + 1);
    //     hsparse->GetAxis(1)->SetRangeUser(ptbins[ipt], ptbins[ipt + 1]);
    //     hInvMassPt[ipt] = (TH1F *)hsparse->Projection(0);
    //     SetHistoStyle_temp(hInvMassPt[ipt], 1, 20, binwidth);
    //     hInvMassPt[ipt]->Draw("pe same");
    //     t2->DrawLatex(0.2, 0.8, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptbins[ipt], ptbins[ipt + 1]));
    //     TF1 *fitpt = fitgaus(hInvMassPt[ipt], ksmass, kswidth);
    // }
}

void SetHistoStyle_temp(TH1 *h, Int_t MCol, Int_t MSty, double binwidth)
{
    h->SetLineWidth(2);
    h->SetTitle(0);
    h->GetXaxis()->SetNdivisions(509);
    h->GetYaxis()->SetNdivisions(509);
    h->GetXaxis()->SetLabelOffset(0.015);
    h->GetXaxis()->SetLabelFont(42);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetTitleFont(42);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetTickLength(0.04);
    h->GetYaxis()->CenterTitle(true);
    h->GetXaxis()->CenterTitle(true);
    h->GetYaxis()->SetLabelOffset(0.035);
    h->GetYaxis()->SetLabelFont(42);
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetTickLength(0.04);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleFont(42);

    // other parameters
    h->GetYaxis()->SetTitleOffset(1.6);
    h->GetXaxis()->SetTitleOffset(1.3);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.0);
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
    h->SetMaximum(1.2 * h->GetMaximum());
    h->GetXaxis()->SetRangeUser(0.45, 0.55);
    h->GetYaxis()->SetTitle(Form("Counts /(%f MeV/#it{c^{2}})", binwidth));
    h->GetXaxis()->SetTitle("#it{M}_{#pi^{+}#pi^{-}} (GeV/#it{c^{2}})");
}

TF1 *fitgaus(TH1 *h, double ksmass, double kswidth)
{
    TF1 *fit = new TF1("fit", "gaus", ksmass - 1.5 * kswidth, ksmass + 1 * kswidth);
    fit->SetParameter(1, ksmass);
    fit->SetParameter(2, kswidth);
    h->Fit(fit, "REI"); // Assuming you meant to use the passed histogram 'h' instead of 'hInvMass'
    fit->SetLineColor(1);
    fit->SetLineWidth(2);
    fit->Draw("SAME");
    return fit;
}

TF1 *fitgauspol2(TH1 *h, double ksmass, double kswidth)
{
    TF1 *fit = new TF1("fit", "gaus(0)+pol2(3)", ksmass - 10 * kswidth, ksmass + 10 * kswidth);
    fit->SetParameter(1, ksmass);
    fit->SetParameter(2, kswidth);
    fit->SetLineColor(1);
    fit->SetLineWidth(2);
    h->Fit(fit, "REI");
    fit->Draw("SAME");
    return fit;
}

TF1 *CBpol2(TH1 *h, double *parameters)
{
    double mass = parameters[1];
    double width = parameters[2];
    TF1 *fit = new TF1("fit", CrystalBallpol2, mass - 15 * width, mass + 15 * width, 8);
    fit->SetParameter(0, parameters[0]);
    fit->SetParameter(1, parameters[1]);
    fit->SetParameter(2, parameters[2]);
    fit->SetParameter(3, parameters[3]);
    fit->SetParameter(4, parameters[4]);
    fit->SetLineColor(1);
    fit->SetLineWidth(2);
    h->Fit(fit, "REBMS+");
    fit->Draw("SAME");
    return fit;
}

TF1 *CB(TH1 *h, double ksmass, double kswidth)
{
    TF1 *fit = new TF1("fit", CrystalBall, ksmass - 2.5 * kswidth, ksmass + 3.0 * kswidth, 5);
    fit->SetParameter(0, 8.1e7);
    fit->SetParameter(1, ksmass);
    fit->SetParameter(2, kswidth);
    fit->SetLineColor(1);
    fit->SetLineWidth(2);
    h->Fit(fit, "REBMS0+");
    // fit->Draw("SAME");
    return fit;
}
