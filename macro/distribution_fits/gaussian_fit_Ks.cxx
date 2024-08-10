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
TF1 *doubleCB(TH1 *h, double *parameters);
TF1 *doubleCBpol2(TH1 *h, double *parameters);
void SetHistoStyle_temp(TH1 *h, Int_t MCol, Int_t MSty, double binwidth);

void gaussian_fit_Ks()
{
    // configurables *********************
    bool saveplots = true;
    bool showpTwise_study = false;
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);
    int rebin = 1;
    // configurables *********************

    double ksmass = pdg->GetParticle(310)->Mass();
    // double kswidth = pdg->GetParticle(310)->Width();
    double kswidth = 0.005;
    cout << "PDG mass: " << ksmass << " PDG width: " << kswidth << endl;

    TLatex *t2 = new TLatex();
    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.05);
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

    TH1F *hInvMass = (TH1F *)hsparseClone->Projection(0, "E");
    hInvMass->Rebin(rebin);
    double binwidth = hInvMass->GetBinWidth(1);
    cout << "Bin width: " << binwidth << endl;

    TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
    SetCanvasStyle(c1, 0.145, 0.05, 0.05, 0.13);
    SetHistoStyle_temp(hInvMass, 1, 20, binwidth);
    hInvMass->Draw("pe");

    TH1F *hInvMassClone1 = (TH1F *)hInvMass->Clone("hInvMassClone1");
    TH1F *hInvMassClone2 = (TH1F *)hInvMass->Clone("hInvMassClone2");
    hInvMassClone1->SetLineColor(0);
    hInvMassClone1->SetLineWidth(0);
    hInvMassClone1->SetFillColor(5);
    hInvMassClone1->GetXaxis()->SetRangeUser(ksmass - 4 * kswidth, ksmass + 4 * kswidth);
    hInvMassClone1->Draw("E3 hist same");
    hInvMassClone2->Draw("pe same");
    // // TF1 *fit = fitgaus(hInvMass, ksmass, kswidth);
    // TF1 *fit3 = fitgauspol2(hInvMass, ksmass, kswidth);
    // fit3->SetNpx(1e6);
    // fit3->SetLineColor(2);
    // fit3->SetLineWidth(2);
    // fit3->Draw("SAME");
    // double *parameters_temp = fit3->GetParameters();

    TF1 *fit = CB(hInvMass, ksmass, kswidth);
    double parameters[5];
    for (int i = 0; i < 5; i++)
    {
        parameters[i] = fit->GetParameter(i);
    }
    // TF1 *fit3 = CBpol2(hInvMass, parameters);
    TF1 *fit2 = doubleCB(hInvMass, parameters);
    double parameters2[7];
    for (int i = 0; i < 7; i++)
    {
        parameters2[i] = fit2->GetParameter(i);
    }
    TF1 *fit3 = doubleCBpol2(hInvMass, parameters2);

    // lets calculate the no. of Ks pairs in the signal region in +- 2 sigma, 3 sigma, 4 sigma, and 5 sigma
    cout << "binwidth: " << binwidth << "\n";
    double allks = hInvMass->Integral();
    double Npairs2sigma = fit3->Integral(ksmass - 2 * kswidth, ksmass + 2 * kswidth) / binwidth;
    double Npairs3sigma = fit3->Integral(ksmass - 3 * kswidth, ksmass + 3 * kswidth) / binwidth;
    double Npairs4sigma = fit3->Integral(ksmass - 4 * kswidth, ksmass + 4 * kswidth) / binwidth;
    double Npairs5sigma = fit3->Integral(ksmass - 5 * kswidth, ksmass + 5 * kswidth) / binwidth;
    cout << "No. and percentage of Ks in +- 2 sigma: " << Npairs2sigma << " " << Npairs2sigma / allks << "\n"
         << "No. and percentage of Ks in +- 3 sigma: " << Npairs3sigma << " " << Npairs3sigma / allks << "\n"
         << "No. and percentage of Ks in +- 4 sigma: " << Npairs4sigma << " " << Npairs4sigma / allks << "\n"
         << "No. and percentage of Ks in +- 5 sigma: " << Npairs5sigma << " " << Npairs5sigma / allks << "\n";

    TLegend *lp2 = DrawLegend(0.14, 0.7, 0.4, 0.9);
    lp2->SetTextSize(0.04);
    lp2->SetTextFont(42);
    lp2->SetFillStyle(0);
    lp2->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    lp2->AddEntry((TObject *)0, "FT0M, 0-100%", "");
    lp2->AddEntry((TObject *)0, "|#it{y}| < 0.5", "");
    lp2->Draw("same");

    // legend
    TLegend *lp3 = DrawLegend(0.6, 0.55, 0.9, 0.75);
    lp3->SetFillStyle(0);
    lp3->SetTextFont(42);
    lp3->AddEntry(hInvMass, "Data", "pe");
    lp3->AddEntry(hInvMassClone1, "Signal", "f");
    lp3->AddEntry(fit3, "CB + pol2 fit", "l");
    lp3->SetTextSize(0.04);
    lp3->Draw("same");

    gPad->Modified();
    gPad->Update();
    TPaveStats *st = (TPaveStats *)hInvMass->FindObject("stats");
    st->SetX1NDC(0.60); // 0.57
    st->SetX2NDC(0.95);
    st->SetY1NDC(0.78); // 0.53
    st->SetY2NDC(0.95);
    st->Draw("same");

    TLegend *lp1 = DrawLegend(0.45, 0.45, 0.9, 0.55);
    lp1->SetFillStyle(0);
    lp1->SetTextFont(42);
    lp1->SetTextSize(0.03);
    lp1->AddEntry((TObject *)0, Form("Mean = %.7f #pm %.7f", fit3->GetParameter(1), fit3->GetParError(1)), "");
    lp1->AddEntry((TObject *)0, Form("Sigma = %.7f #pm %.7f", fit3->GetParameter(2), fit3->GetParError(2)), "");
    lp1->Draw("same");
    if (saveplots)
    {
        c1->SaveAs("gaussian_fit_Ks_fullpt.png");
    }
    if (showpTwise_study)
    {
        // Now we will plot the Ks invariant mass distribution as a function of pT
        const int Nptbins = 16;
        float ptbins[Nptbins + 1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 9.0, 12.0, 15.0, 20.0, 30.0};
        TH1F *hInvMassPt[Nptbins];
        TCanvas *c2 = new TCanvas("c2", "c2", 720, 720);
        c2->Divide(4, 4);
        TH1F *hpTwiseMass = new TH1F("hpTwiseMass", "hpTwiseMass", Nptbins, ptbins);
        TH1F *hpTwiseWidth = new TH1F("hpTwiseWidth", "hpTwiseWidth", Nptbins, ptbins);
        TH1F *hpTwiseNorm = new TH1F("hpTwiseNorm", "hpTwiseNorm", Nptbins, ptbins);
        for (int ipt = 0; ipt < Nptbins; ipt++)
        {
            c2->cd(ipt + 1);
            THnSparse *hsprase_clone2 = (THnSparse *)hsparse->Clone("hsprase_clone2");
            int lowptbin = hsprase_clone2->GetAxis(1)->FindBin(ptbins[ipt]);
            int highptbin = hsprase_clone2->GetAxis(1)->FindBin(ptbins[ipt + 1]);
            hsprase_clone2->GetAxis(1)->SetRange(lowptbin, highptbin);
            hInvMassPt[ipt] = (TH1F *)hsprase_clone2->Projection(0);
            hInvMassPt[ipt]->Rebin(rebin);
            SetHistoStyle_temp(hInvMassPt[ipt], 1, 20, binwidth);
            hInvMassPt[ipt]->SetTitle("");
            hInvMassPt[ipt]->Draw("pe");
            t2->DrawLatex(0.16, 0.82, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptbins[ipt], ptbins[ipt + 1]));
            TF1 *fitpt1;
            TF1 *fitpt2;
            TF1 *fitpt3;
            if (ipt < 8)
            {
                fitpt3 = doubleCBpol2(hInvMassPt[ipt], parameters2);
            }
            else if (ipt >= 8)
            {
                fitpt1 = CB(hInvMassPt[ipt], ksmass, kswidth);
                double parameters1[5];
                for (int i = 0; i < 5; i++)
                {
                    parameters1[i] = fitpt1->GetParameter(i);
                }
                fitpt2 = doubleCB(hInvMassPt[ipt], parameters1);
                double parameters3[7];
                for (int i = 0; i < 7; i++)
                {
                    parameters3[i] = fitpt2->GetParameter(i);
                }
                cout << "\nnow fitting with doubleCBpol2 with pT bin: " << ipt << "\n\n";
                fitpt3 = doubleCBpol2(hInvMassPt[ipt], parameters3);
                // double parameters3[7] = {fit3->GetParameter(0), fit3->GetParameter(1), fit3->GetParameter(2), fit3->GetParameter(3)+2.0, fit3->GetParameter(4), fit3->GetParameter(5), fit3->GetParameter(6)};
                // fitpt3 = doubleCBpol2(hInvMassPt[ipt], parameters3);
            }

            // trying gaussian fit instead of double CB
            //  TF1 *fitpt3 = new TF1("fitpt3", "gaus(0)+pol2(3)", ksmass - 10 * kswidth, ksmass + 10 * kswidth);
            //  fitpt3->SetParameter(1, ksmass);
            //  fitpt3->SetParameter(2, kswidth);
            //  // fitpt3->SetParameter(3, parameters_temp[3]);
            //  // fitpt3->SetParameter(4, parameters_temp[4]);
            //  // fitpt3->SetParameter(5, parameters_temp[5]);
            //  hInvMassPt[ipt]->Fit(fitpt3, "REBMS0+");
            //  fitpt3->Draw("same");
            //  cout << "fit mean is " << fitpt3->GetParameter(1) << " and fit width is " << fitpt3->GetParameter(2) << endl;
            hpTwiseMass->SetBinContent(ipt + 1, fitpt3->GetParameter(1));
            hpTwiseMass->SetBinError(ipt + 1, fitpt3->GetParError(1));
            hpTwiseWidth->SetBinContent(ipt + 1, abs(fitpt3->GetParameter(2)));
            hpTwiseWidth->SetBinError(ipt + 1, fitpt3->GetParError(2));
            hpTwiseNorm->SetBinContent(ipt + 1, fitpt3->GetParameter(0));
            hpTwiseNorm->SetBinError(ipt + 1, fitpt3->GetParError(0));

            TLatex lat;
            lat.SetNDC();
            lat.SetTextFont(42);
            lat.SetTextSize(0.04);
            lat.DrawLatex(0.57, 0.75, Form("Mean = %.3f #pm %.3f", fitpt3->GetParameter(1), fitpt3->GetParError(1)));
            lat.DrawLatex(0.57, 0.7, Form("Sigma = %.3f #pm %.3f", abs(fitpt3->GetParameter(2)), fitpt3->GetParError(2)));
        }
        if (saveplots)
        {
            c2->SaveAs("gaussian_fit_Ks_differential_ptbins.png");
        }

        TCanvas *c3 = new TCanvas("c3", "c3", 720, 720);
        SetCanvasStyle(c3, 0.22, 0.05, 0.05, 0.13);
        SetHistoQA(hpTwiseMass);
        hpTwiseMass->SetMarkerSize(1);
        hpTwiseMass->GetYaxis()->SetTitle("Fit Mean (GeV/#it{c^{2}})");
        hpTwiseMass->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hpTwiseMass->GetYaxis()->SetMaxDigits(3);
        hpTwiseMass->GetYaxis()->SetTitleOffset(2.3);
        hpTwiseMass->SetStats(0);
        hpTwiseMass->GetYaxis()->SetRangeUser(0.49, 0.502);
        hpTwiseMass->Draw("pe");
        TLine *line = new TLine(0.0, ksmass, 30.0, ksmass);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(kRed);
        line->Draw("l same");
        TLegend *lp4 = DrawLegend(0.6, 0.75, 0.9, 0.85);
        lp4->SetFillStyle(0);
        lp4->SetTextFont(42);
        lp4->SetTextSize(0.04);
        lp4->AddEntry(hpTwiseMass, "Fit Mean", "lpe");
        lp4->AddEntry(line, "PDG Mass", "l");
        lp4->Draw("same");
        if (saveplots)
        {
            c3->SaveAs("saved/gaussian_fit_Ks_differential_ptbins_mean.png");
        }

        TCanvas *c4 = new TCanvas("c4", "c4", 720, 720);
        SetCanvasStyle(c4, 0.13, 0.05, 0.05, 0.13);
        SetHistoQA(hpTwiseWidth);
        hpTwiseWidth->SetMarkerSize(1);
        hpTwiseWidth->GetYaxis()->SetTitle("Fit Width (GeV/#it{c^{2}})");
        hpTwiseWidth->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hpTwiseWidth->GetYaxis()->SetMaxDigits(3);
        hpTwiseWidth->GetYaxis()->SetTitleOffset(1.2);
        hpTwiseWidth->SetStats(0);
        hpTwiseWidth->GetYaxis()->SetRangeUser(0.0, 0.015);
        hpTwiseWidth->Draw("pe");
        if (saveplots)
        {
            c4->SaveAs("saved/gaussian_fit_Ks_differential_ptbins_width.png");
        }

        c4->Clear();
        gPad->SetLogy();
        SetCanvasStyle(c4, 0.13, 0.05, 0.05, 0.13);
        SetHistoQA(hpTwiseNorm);
        hpTwiseNorm->SetMarkerSize(1);
        hpTwiseNorm->GetYaxis()->SetTitle("Fit Normalization");
        hpTwiseNorm->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hpTwiseNorm->GetYaxis()->SetMaxDigits(3);
        hpTwiseNorm->GetYaxis()->SetTitleOffset(1.3);
        hpTwiseNorm->SetStats(0);
        hpTwiseNorm->GetYaxis()->SetRangeUser(50, 1.0e8);
        hpTwiseNorm->Draw("pe");
        if (saveplots)
        {
            c4->SaveAs("saved/gaussian_fit_Ks_differential_ptbins_norm.png");
        }
    }
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
    h->GetYaxis()->SetTitle(Form("Counts /(%.5f GeV/#it{c^{2}})", binwidth));
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
    fit->SetLineColor(2);
    fit->SetLineWidth(2);
    h->Fit(fit, "REBMS0+");
    // fit->Draw("SAME");
    return fit;
}

TF1 *CB(TH1 *h, double ksmass, double kswidth)
{
    TF1 *fit = new TF1("fit", CrystalBall, ksmass - 2.5 * kswidth, ksmass + 3.0 * kswidth, 5);
    fit->SetParNames("Norm", "Mean", "Sigma", "Alpha", "n");
    fit->SetParameter(0, 8.1e7);
    fit->SetParameter(1, ksmass);
    fit->SetParameter(2, kswidth);
    fit->SetLineColor(1);
    fit->SetLineWidth(2);
    h->Fit(fit, "REBMSQ0+");
    // fit->Draw("SAME");
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
    h->Fit(fit, "REBMSQ+");
    fit->Draw("SAME");
    return fit;
}

TF1 *doubleCB(TH1 *h, double *parameters)
{
    double mass = parameters[1];
    double width = parameters[2];
    TF1 *fit = new TF1("fit", DoubleCrystalBall, mass - 3.5 * width, mass + 4.0 * width, 7);
    fit->SetParNames("Norm", "Mean", "Sigma", "AlphaL", "n1", "AlphaR", "n2");
    fit->SetParameter(0, parameters[0]);
    fit->SetParameter(1, mass);
    fit->SetParameter(2, width);
    fit->FixParameter(3, parameters[3]);
    fit->FixParameter(4, parameters[4]);
    fit->SetParameter(5, parameters[3]);
    fit->SetParameter(6, parameters[4]);
    fit->SetLineColor(1);
    fit->SetLineWidth(2);
    h->Fit(fit, "REBMSQ0+");
    // fit->Draw("SAME");
    return fit;
}

TF1 *doubleCBpol2(TH1 *h, double *parameters)
{
    double mass = parameters[1];
    double width = parameters[2];
    TF1 *fit = new TF1("fit", DoubleCrystalBallpol2, mass - 15 * width, mass + 15 * width, 10);
    fit->SetParNames("Norm", "Mean", "Sigma", "AlphaL", "n1", "AlphaR", "n2", "p0", "p1", "p2");
    fit->SetParameter(0, parameters[0]);
    fit->SetParameter(1, mass);
    // fit->SetParLimits(2, 0.0, 20.0);
    fit->SetParameter(2, width);
    fit->SetParLimits(3, 0.0, 10.0);
    fit->SetParameter(3, parameters[3]);
    fit->SetParLimits(4, 0.0, 10.0);
    fit->SetParameter(4, parameters[4]);
    fit->SetParLimits(5, 0.0, 10.0);
    fit->SetParameter(5, parameters[5]);
    fit->SetParLimits(6, 0.0, 10.0);
    fit->SetParameter(6, parameters[6]);
    if (sizeof(parameters) == 10)
    {
        fit->SetParameter(7, parameters[7]);
        fit->SetParameter(8, parameters[8]);
        fit->SetParameter(9, parameters[9]);
    }
    fit->SetLineColor(2);
    fit->SetLineWidth(2);
    h->Fit(fit, "REBMS");
    fit->SetNpx(1e6);
    fit->Draw("SAME");
    return fit;
}
