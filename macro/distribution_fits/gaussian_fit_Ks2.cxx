#include <iostream>
#include "TDatabasePDG.h"
#include "../src/fitfunc.h"
#include "../src/common_glue.h"
#include "../src/style.h"

TDatabasePDG *pdg = new TDatabasePDG();
double ksmass = pdg->GetParticle(310)->Mass();
double kswidth = 0.005;

TF1 *fitgaus(TH1 *h, double ksmass, double kswidth);
TF1 *fitgauspol2(TH1 *h, double ksmass, double kswidth);
TF1 *CBpol2(TH1 *h, double *parameters, bool mainfit);
TF1 *CB(TH1 *h, double ksmass, double kswidth);
TF1 *doubleCB(TH1 *h, double *parameters, bool mainfit);
TF1 *doubleCBpol2(TH1 *h, double *parameters, bool mainfit, TLegend *leg = nullptr, float legendsize = 0.04, float rlow = 15, float rhigh = 15);
TF1 *doubleCBpol3(TH1 *h, double *parameters, bool mainfit, TLegend *leg = nullptr, float legendsize = 0.04, float rlow = 15, float rhigh = 15);
TF1 *doubleCBpol1(TH1 *h, double *parameters, bool mainfit);
void SetHistoStyle_temp(TH1 *h, Int_t MCol, Int_t MSty, double binwidth);
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);
void directfit(TF1 *&fit, TH1 *h, double *parameters, bool mainfit = false);
void perform_fit(TF1 *&fit, TH1 *h, bool mainfit = false);

void gaussian_fit_Ks2()
{
    // configurables *********************
    bool saveplots = false;
    bool showpt_study = false;
    gStyle->SetOptStat(1110);
    gStyle->SetFitFormat("7.7g"); // 6 significant digits
    gStyle->SetOptFit(1111);
    string whichpass = "pass7";

    int rebin = 3; // for rebin 3, change the normalization factor in CB function (and use 254232 output)
    // configurables *********************

    // double kswidth = pdg->GetParticle(310)->Width();
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
    double pad1Size, pad2Size;
    canvas_style(c1, pad1Size, pad2Size);
    c1->cd(1);
    SetHistoStyle_temp(hInvMass, 1, 20, binwidth);
    hInvMass->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hInvMass->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hInvMass->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hInvMass->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hInvMass->GetYaxis()->SetTitleOffset(1.2);
    hInvMass->SetMinimum(0.001);
    hInvMass->Draw("pe");
    auto noofevents = hInvMass->Integral();

    TH1F *hInvMassClone1 = (TH1F *)hInvMass->Clone("hInvMassClone1");
    TH1F *hInvMassClone2 = (TH1F *)hInvMass->Clone("hInvMassClone2");
    hInvMassClone1->SetLineColor(0);
    hInvMassClone1->SetLineWidth(0);
    hInvMassClone1->SetFillColor(5);
    hInvMassClone1->GetXaxis()->SetRangeUser(ksmass - 4.3 * kswidth, ksmass + 4 * kswidth);
    hInvMassClone1->Draw("E3 hist same");
    hInvMassClone2->Draw("pe same");

    TF1 *fit = CB(hInvMass, ksmass, kswidth); // single crystal ball fit
    double parameters[5];
    for (int i = 0; i < 5; i++)
    {
        parameters[i] = fit->GetParameter(i);
    }
    fit->SetLineColor(2);

    TF1 *fit2 = doubleCB(hInvMass, parameters, false); // double crystal ball fit
    double parameters2[7];
    for (int i = 0; i < 7; i++)
    {
        parameters2[i] = fit2->GetParameter(i);
    }
    fit2->SetLineColor(2);
    fit2->SetLineWidth(2);
    // fit2->Draw("SAME");
    // legend
    TLegend *lp3 = DrawLegend(0.17, 0.25, 0.52, 0.64);
    lp3->SetFillStyle(0);
    lp3->SetTextFont(42);
    lp3->AddEntry(hInvMass, "Data", "pe");
    lp3->AddEntry(hInvMassClone1, "Signal", "f");
    // lp3->AddEntry(fit3, "CB + pol2 fit", "l");

    TF1 *fit3 = doubleCBpol2(hInvMass, parameters2, true, lp3, 0.04, 15, 15); // double crystal ball with pol2 fit
    fit3->Draw("SAME");

    // iteration for successful fit
    float rangelow = 8;
    float rangehigh = 8;
    bool fitsuccessfully = false;
    // while (!fitsuccessfully)
    // {
    //     fit3 = doubleCBpol2(hInvMass, parameters2, true, lp3, 0.04, rangelow, rangehigh); // double crystal ball with pol2 fit
    //     TString fitStatus = gMinuit->fCstatu;
    //     if (fitStatus.Contains("CONVERGED"))
    //     {
    //         fitsuccessfully = true;
    //         cout << "Fit converged successfully" << endl;
    //         cout << "The value of range low and range high is: " << rangelow << " " << rangehigh << endl;
    //     }
    //     else
    //     {
    //         rangelow += 0.1;
    //         if (rangelow > 20)
    //         {
    //             rangehigh += 0.1;
    //             rangelow = 8;
    //         }
    //         if (rangehigh > 20)
    //         {
    //             cout << "Fit failed to converge" << endl;
    //             break;
    //         }
    //         cout << "\n\nRange low and high is " << rangelow << " " << rangehigh << "\n\n";
    //     }
    // }
    // TF1 *fit3 = doubleCBpol1(hInvMass, parameters, false); // double crystal ball with pol1 fit
    cout << "The value and error of alpha Left is: " << fit3->GetParameter(3) << " " << fit3->GetParError(3) << endl;
    cout << "The value and error of alpha right is : " << fit3->GetParameter(5) << " " << fit3->GetParError(5) << endl;

    // lets calculate the no. of Ks pairs in the signal region in +- 2 sigma, 3 sigma, 4 sigma, and 5 sigma
    cout << "binwidth: " << binwidth << "\n";
    double allks = hInvMass->Integral();
    double Npairs2sigma = fit3->Integral(ksmass - 2 * kswidth, ksmass + 2 * kswidth) / binwidth;
    double Npairs3sigma = fit3->Integral(ksmass - 3 * kswidth, ksmass + 3 * kswidth) / binwidth;
    double Npairs4sigma = fit3->Integral(ksmass - 4 * kswidth, ksmass + 4 * kswidth) / binwidth;
    double Npairs5sigma = fit3->Integral(ksmass - 5 * kswidth, ksmass + 5 * kswidth) / binwidth;
    double Npairs6sigma = fit3->Integral(ksmass - 6 * kswidth, ksmass + 6 * kswidth) / binwidth;
    cout << "No. and percentage of Ks in +- 2 sigma: " << Npairs2sigma << " " << Npairs2sigma / allks << "\n"
         << "No. and percentage of Ks in +- 3 sigma: " << Npairs3sigma << " " << Npairs3sigma / allks << "\n"
         << "No. and percentage of Ks in +- 4 sigma: " << Npairs4sigma << " " << Npairs4sigma / allks << "\n"
         << "No. and percentage of Ks in +- 5 sigma: " << Npairs5sigma << " " << Npairs5sigma / allks << "\n"
         << "No. and percentage of Ks in +- 6 sigma: " << Npairs6sigma << " " << Npairs6sigma / allks << "\n";

    // Purity of signal = S/(S+B)
    // S = signal yield from functional integration fit after subtracting the background
    // S+B = total yield from the invariant mass distribution before subtracting the background
    // purity for 4 sigma invariant mass region
    double sigbkg2 = hInvMass->Integral(hInvMass->GetXaxis()->FindBin(ksmass - 2 * kswidth), hInvMass->GetXaxis()->FindBin(ksmass + 2 * kswidth));
    double sigbkg3 = hInvMass->Integral(hInvMass->GetXaxis()->FindBin(ksmass - 3 * kswidth), hInvMass->GetXaxis()->FindBin(ksmass + 3 * kswidth));
    double sigbkg4 = hInvMass->Integral(hInvMass->GetXaxis()->FindBin(ksmass - 4 * kswidth), hInvMass->GetXaxis()->FindBin(ksmass + 4 * kswidth));
    double sigbkg5 = hInvMass->Integral(hInvMass->GetXaxis()->FindBin(ksmass - 5 * kswidth), hInvMass->GetXaxis()->FindBin(ksmass + 5 * kswidth));
    double sigbkg6 = hInvMass->Integral(hInvMass->GetXaxis()->FindBin(ksmass - 6 * kswidth), hInvMass->GetXaxis()->FindBin(ksmass + 6 * kswidth));
    double purity2 = Npairs2sigma / sigbkg2;
    double purity3 = Npairs3sigma / sigbkg3;
    double purity4 = Npairs4sigma / sigbkg4;
    double purity5 = Npairs5sigma / sigbkg5;
    double purity6 = Npairs6sigma / sigbkg6;
    cout << "Purity of signal in 2 sigma: " << purity2 << "\n"
         << "Purity of signal in 3 sigma: " << purity3 << "\n"
         << "Purity of signal in 4 sigma: " << purity4 << "\n"
         << "Purity of signal in 5 sigma: " << purity5 << "\n"
         << "Purity of signal in 6 sigma: " << purity6 << "\n";

    TLegend *lp2 = DrawLegend(0.16, 0.68, 0.42, 0.85);
    lp2->SetTextSize(0.04);
    lp2->SetTextFont(42);
    lp2->SetFillStyle(0);
    lp2->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    lp2->AddEntry((TObject *)0, "FT0M, 0-100%", "");
    lp2->AddEntry((TObject *)0, "|#it{y}| < 0.5", "");
    lp2->Draw("same");

    gPad->Modified();
    gPad->Update();
    TPaveStats *st = (TPaveStats *)hInvMass->FindObject("stats");
    st->SetX1NDC(0.6); // 0.60
    st->SetX2NDC(0.95);
    st->SetY1NDC(0.3); // 0.78
    st->SetY2NDC(0.95);
    st->Draw("same");

    TLine *linevert = new TLine(ksmass - 4 * kswidth, 0, ksmass - 4 * kswidth, hInvMass->GetMaximum());
    linevert->SetLineStyle(2);
    linevert->SetLineColor(3);
    linevert->SetLineWidth(3);
    linevert->Draw("same");
    TLine *linevert2 = new TLine(ksmass + 4 * kswidth, 0, ksmass + 4 * kswidth, hInvMass->GetMaximum());
    linevert2->SetLineStyle(2);
    linevert2->SetLineColor(3);
    linevert2->SetLineWidth(3);
    linevert2->Draw("same");

    TLegend *lp1 = DrawLegend(0.52, 0.35, 0.9, 0.45);
    lp1->SetFillStyle(0);
    lp1->SetTextFont(42);
    lp1->SetTextSize(0.03);
    lp1->AddEntry((TObject *)0, Form("Mean = %.3f #pm %.2e", fit3->GetParameter(1), fit3->GetParError(1)), "");
    lp1->AddEntry((TObject *)0, Form("Sigma = %.3f #pm %.2e", fit3->GetParameter(2), fit3->GetParError(2)), "");
    // lp1->Draw("same");

    // lets calculate the fit to the data ratio
    c1->cd(2);
    TH1F *hInvMassRatio = (TH1F *)hInvMass->Clone("hInvMassRatio");
    // for (int i = 0; i < hInvMass->GetNbinsX(); i++)
    // {
    //     hInvMassRatio->SetBinContent(i + 1, hInvMass->GetBinContent(i + 1) / fit3->Eval(hInvMass->GetBinCenter(i + 1)));
    //     hInvMassRatio->SetBinError(i + 1, hInvMass->GetBinError(i + 1) / fit3->Eval(hInvMass->GetBinCenter(i + 1)));
    // }
    hInvMassRatio->Divide(fit3);

    SetHistoStyle_temp(hInvMassRatio, 1, 20, binwidth);
    hInvMassRatio->SetMarkerSize(0.5);
    hInvMassRatio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hInvMassRatio->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    hInvMassRatio->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hInvMassRatio->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    hInvMassRatio->GetYaxis()->SetRangeUser(0.6, 1.39);
    hInvMassRatio->GetYaxis()->SetTitle("Data/Fit");
    hInvMassRatio->GetYaxis()->SetTitleOffset(0.5);
    hInvMassRatio->GetXaxis()->SetTitleOffset(1.15);
    hInvMassRatio->GetYaxis()->SetNdivisions(505);
    hInvMassRatio->SetStats(0);
    hInvMassRatio->Draw("pe");

    TFile *foutput = new TFile(Form("saved/%s/output_rebin%d.root", whichpass.c_str(), rebin), "recreate");
    if (saveplots)
    {
        c1->SaveAs(Form("saved/%s/gaussfit_Ks_rebin%d.png", whichpass.c_str(), rebin));
    }
    c1->Write("ks_fit");

    if (showpt_study)
    {
        // Now we will plot the Ks invariant mass distribution as a function of pT
        const int Nptbins = 16;
        double param_allpt[Nptbins][10];
        float ptbins[Nptbins + 1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 9.0, 12.0, 15.0, 20.0, 30.0};
        TH1F *hInvMassPt[Nptbins];
        TCanvas *c2 = new TCanvas("c2", "c2", 720, 720);
        c2->Divide(4, 4);
        TH1F *hpTwiseMass = new TH1F("hpTwiseMass", "hpTwiseMass", Nptbins, ptbins);
        TH1F *hpTwiseWidth = new TH1F("hpTwiseWidth", "hpTwiseWidth", Nptbins, ptbins);
        TH1F *hpTwise_yield_int = new TH1F("hpTwise_yield_int", "hpTwise_yield_int", Nptbins, ptbins);
        TH1F *hpTwise_yield_bin = new TH1F("hpTwise_yield_bin", "hpTwise_yield_bin", Nptbins, ptbins);
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

            // bool conditions[] = { //for pass 7 with low statistics
            //     rebin == 1 && ipt < 7,                 // Condition 0
            //     (rebin == 2) && ipt < 4,               // Condition 1
            //     (rebin == 3) && ipt < 5,               // Condition 2
            //     rebin == 1 && (ipt == 9 || ipt == 10), // Condition 3
            //     rebin == 1 && ipt == Nptbins - 2,      // Condition 4
            //     rebin == 2 && ipt == Nptbins - 2,       // Condition 5
            //     rebin == 3 && ipt == Nptbins - 1       // Condition 6
            // };

            bool conditions[] = {
                // for pass 8 with full statistics
                (rebin == 1 || rebin == 3) && ipt < 4,                     // Condition 0
                (rebin == 2) && ipt < 3,                                   // Condition 1
                rebin == 10 && (ipt == Nptbins - 2 || ipt == Nptbins - 1), // Condition 2
                rebin == 10 && (ipt == 9 || ipt == 10),                    // Condition 3
                rebin == 1 && ipt == Nptbins - 2,                          // Condition 4
                rebin == 3 && ipt == Nptbins - 1,                          // Condition 5
                (rebin == 1 || rebin == 2) && ipt == Nptbins - 1           // Condition 6
            };

            if (conditions[0])
            {
                directfit(fitpt3, hInvMassPt[ipt], parameters2, false);
            }
            else if (conditions[1])
            {
                directfit(fitpt3, hInvMassPt[ipt], parameters2, false);
            }
            else if (conditions[2])
            {
                directfit(fitpt3, hInvMassPt[ipt], parameters2, false);
            }

            else if (conditions[3])
            {
                fitpt3 = doubleCBpol2(hInvMassPt[ipt], param_allpt[5], false);
                fitpt3->Draw("SAME");
            }
            else if (conditions[4])
            {
                fitpt3 = doubleCBpol2(hInvMassPt[ipt], param_allpt[Nptbins - 4], false);
                fitpt3->Draw("SAME");
            }
            else if (conditions[5])
            {
                fitpt3 = doubleCBpol2(hInvMassPt[ipt], param_allpt[Nptbins - 4], false);
                fitpt3->Draw("SAME");
            }
            else if (conditions[6])
            {
                fitpt3 = doubleCBpol2(hInvMassPt[ipt], param_allpt[Nptbins - 4], false);
                fitpt3->Draw("SAME");
            }
            else
            {
                fitpt3 = doubleCBpol2(hInvMassPt[ipt], param_allpt[ipt - 1], false);
                fitpt3->Draw("SAME");
            }

            for (int iparam = 0; iparam < 10; iparam++)
            {
                param_allpt[ipt][iparam] = fitpt3->GetParameter(iparam);
            }

            hpTwiseMass->SetBinContent(ipt + 1, fitpt3->GetParameter(1));
            hpTwiseMass->SetBinError(ipt + 1, fitpt3->GetParError(1));
            hpTwiseWidth->SetBinContent(ipt + 1, abs(fitpt3->GetParameter(2)));
            hpTwiseWidth->SetBinError(ipt + 1, fitpt3->GetParError(2));

            // TLatex lat;
            // lat.SetNDC();
            // lat.SetTextFont(42);
            // lat.SetTextSize(0.04);
            // lat.DrawLatex(0.57, 0.75, Form("Mean = %.3f #pm %.3f", fitpt3->GetParameter(1), fitpt3->GetParError(1)));
            // lat.DrawLatex(0.57, 0.7, Form("Sigma = %.3f #pm %.3f", abs(fitpt3->GetParameter(2)), fitpt3->GetParError(2)));

            // Yield from bin counting method calculation
            TF1 *fitFcn2_plusm = new TF1("fitFcn2_plusm", DoubleCrystalBallpol2, fitpt3->GetXmin(), fitpt3->GetXmax(), 7);
            TF1 *fitFcn2_minusm = new TF1("fitFcn2_minusm", DoubleCrystalBallpol2, fitpt3->GetXmin(), fitpt3->GetXmax(), 7);

            for (int ipar = 0; ipar < 7; ipar++)
            {
                if (ipar == 1)
                {
                    fitFcn2_plusm->FixParameter(ipar, fitpt3->GetParameter(ipar) + fitpt3->GetParError(ipar));
                    fitFcn2_minusm->FixParameter(ipar, fitpt3->GetParameter(ipar) - fitpt3->GetParError(ipar));
                }
                else
                {
                    fitFcn2_plusm->FixParameter(ipar, fitpt3->GetParameter(ipar));
                    fitFcn2_minusm->FixParameter(ipar, fitpt3->GetParameter(ipar));
                }
            }

            // finding the histogram bin for the mass plus and minus its two times the width (the width is not known, so we will take the value of 10 MeV = 0.01 GeV for now)

            TF1 *double_CB_fit = new TF1("double_CB_fit", DoubleCrystalBall, fitpt3->GetXmin(), fitpt3->GetXmax(), 7);
            TF1 *pol2_fit = new TF1("pol2_fit", polynomial2, fitpt3->GetXmin(), fitpt3->GetXmax(), 3);
            for (int ipar = 0; ipar < 10; ipar++)
            {
                if (ipar < 7)
                {
                    double_CB_fit->FixParameter(ipar, fitpt3->GetParameter(ipar));
                }
                else
                {
                    pol2_fit->FixParameter(ipar - 7, fitpt3->GetParameter(ipar));
                }
            }
            auto ptbinwidth = ptbins[ipt + 1] - ptbins[ipt];

            double total_yield = double_CB_fit->Integral(0.2, 0.8) / (binwidth * ptbinwidth * noofevents);
            double total_yield_error = fitpt3->IntegralError(0.2, 0.8) / (binwidth * ptbinwidth * noofevents);
            cout << "Yield from functional integration: " << total_yield << endl;
            cout << "Error in yield from functional integration: " << total_yield_error << endl;
            if (total_yield_error > total_yield)
            {
                total_yield_error = 0;
            }
            hpTwise_yield_int->SetBinContent(ipt + 1, total_yield);
            hpTwise_yield_int->SetBinError(ipt + 1, total_yield_error);

            auto bin_min = hInvMassPt[ipt]->FindBin(ksmass - 2 * 0.01);
            auto bin_max = hInvMassPt[ipt]->FindBin(ksmass + 2 * 0.01);
            double bc_error;
            double Yield_bincount_hist = hInvMassPt[ipt]->IntegralAndError(bin_min, bin_max, bc_error);
            cout << "bc_error is " << bc_error << endl;

            double bkgvalue = pol2_fit->Integral(hInvMassPt[ipt]->GetBinLowEdge(bin_min), hInvMassPt[ipt]->GetBinLowEdge(bin_max + 1));
            double Integral_BW_withsigma = double_CB_fit->Integral(hInvMassPt[ipt]->GetBinLowEdge(bin_min), hInvMassPt[ipt]->GetBinLowEdge(bin_max + 1));

            auto fYield_BinCount = Yield_bincount_hist - (bkgvalue / binwidth);
            auto YieldIntegral_BW = double_CB_fit->Integral(0.2, 0.8) / binwidth;
            auto Yfraction_cBW = (Integral_BW_withsigma / YieldIntegral_BW);

            auto sum_tail_correction = (double_CB_fit->Integral(0.2, hInvMassPt[ipt]->GetBinLowEdge(bin_min)) + double_CB_fit->Integral(hInvMassPt[ipt]->GetBinLowEdge(bin_max + 1), 0.8)) / binwidth;
            auto Total_Ybincounting = (sum_tail_correction + fYield_BinCount) / (ptbinwidth * noofevents);
            cout << "Total_Ybincounting: " << Total_Ybincounting << endl;
            // cout<<"Error in function integration yield: "<<double_CB_fit->IntegralError(ksmass - 3 * 0.01, ksmass + 3 * 0.01) / binwidth<<endl;

            auto Tail_correction_plusm = (fitFcn2_plusm->Integral(0.2, hInvMassPt[ipt]->GetBinLowEdge(bin_min)) + (fitFcn2_plusm->Integral(hInvMassPt[ipt]->GetBinLowEdge(bin_max + 1), 0.8))) / binwidth;
            cout << "Tail_correction_plusm: " << Tail_correction_plusm << endl;
            auto Tail_correction_minusm = (fitFcn2_minusm->Integral(0.2, hInvMassPt[ipt]->GetBinLowEdge(bin_min)) + (fitFcn2_minusm->Integral(hInvMassPt[ipt]->GetBinLowEdge(bin_max + 1), 0.8))) / binwidth;
            cout << "Tail_correction_minusm: " << Tail_correction_minusm << endl;
            auto Error_2 = (Tail_correction_plusm - Tail_correction_minusm) / 2;
            cout << "Error_2: " << Error_2 << endl;
            // auto Final_pro_error = sqrt(pow(bc_error, 2) + pow(Error_2, 2)) / (ptbinwidth * noofevents);
            auto Final_pro_error = sqrt(pow(bc_error, 2)) / (ptbinwidth * noofevents);
            cout << "Final_pro_error: " << Final_pro_error << endl;

            hpTwise_yield_bin->SetBinContent(ipt + 1, Total_Ybincounting);
            hpTwise_yield_bin->SetBinError(ipt + 1, Final_pro_error);
        }
        if (saveplots)
        {
            c2->SaveAs(Form("saved/%s/gaussfit_Ks_all_ptbins_rebin%d.png", whichpass.c_str(), rebin));
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
        hpTwiseMass->GetYaxis()->SetRangeUser(0.492, 0.502);
        hpTwiseMass->Draw("pe");
        hpTwiseMass->Write("ks_mass_fit");
        TLine *line = new TLine(0.0, ksmass, 30.0, ksmass);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(kRed);
        line->Draw("l same");
        // create a band around the line of width equal to the pdg error bar (0.013 * 1e-3)
        double xvalues[2] = {0.0, 30.0};
        double yvalues[2] = {ksmass, ksmass};
        double exvalues[2] = {0.0, 0.0};
        double eyvalues[2] = {0.000013, 0.000013};
        TGraphErrors *band = new TGraphErrors(2, xvalues, yvalues, exvalues, eyvalues);
        band->SetFillColor(kBlue); // Blue color with 30% opacity
        band->SetFillStyle(3001);
        band->SetLineColor(kBlue);
        band->Draw("3"); // "3" option draws a filled band
        TLegend *lp4 = DrawLegend(0.23, 0.73, 0.9, 0.9);
        lp4->SetFillStyle(0);
        lp4->SetTextFont(42);
        lp4->SetTextSize(0.035);
        lp4->AddEntry(hpTwiseMass, "Fit Mean", "lpe");
        // lp4->AddEntry(line, "PDG Mass ", "l");
        lp4->AddEntry(band, "PDG Mass (0.4976 #pm 1.3e-5 MeV/c^{2})", "f");
        lp4->Draw("same");
        if (saveplots)
        {
            c3->SaveAs(Form("saved/%s/gaussfit_Ks_mean_rebin%d.png", whichpass.c_str(), rebin));
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
        hpTwiseWidth->Write("ks_width_fit");
        if (saveplots)
        {
            c4->SaveAs(Form("saved/%s/gaussfit_Ks_width_rebin%d.png", whichpass.c_str(), rebin));
        }

        TCanvas *c5 = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c5, 0.3, 0.03, 0.03, 0.2);
        canvas_style(c5, pad1Size, pad2Size);
        c5->cd(1);
        c5->SetLogy();
        gPad->SetLogy();
        SetHistoQA(hpTwise_yield_int);
        SetHistoQA(hpTwise_yield_bin);
        SetHistoStyle(hpTwise_yield_int, 1, 20, 1, 0.04 / pad1Size, 0.04 / pad1Size, 0.04 / pad1Size, 0.04 / pad1Size, 1.3, 1.3);
        hpTwise_yield_int->SetMarkerSize(1);
        hpTwise_yield_bin->SetMarkerSize(1);
        hpTwise_yield_int->SetMarkerStyle(21);
        hpTwise_yield_int->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
        hpTwise_yield_int->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hpTwise_yield_int->GetYaxis()->SetMaxDigits(3);
        hpTwise_yield_int->GetYaxis()->SetTitleOffset(1.15);
        hpTwise_yield_int->SetStats(0);
        hpTwise_yield_int->Draw("pe");
        hpTwise_yield_bin->SetMarkerColor(kRed);
        hpTwise_yield_bin->SetLineColor(kRed);
        hpTwise_yield_bin->Draw("pe same");
        hpTwise_yield_int->Write("ks_yield_int");
        hpTwise_yield_bin->Write("ks_yield_bin");
        TLegend *lp5 = DrawLegend(0.4, 0.75, 0.7, 0.85);
        lp5->SetFillStyle(0);
        lp5->SetTextFont(42);
        lp5->SetTextSize(0.04);
        lp5->AddEntry(hpTwise_yield_int, "Functional integration", "lpe");
        lp5->AddEntry(hpTwise_yield_bin, "Bin counting", "lpe");
        lp5->Draw("same");

        TH1F *hratios = (TH1F *)hpTwise_yield_int->Clone("hratios");
        hratios->Divide(hpTwise_yield_bin);

        c5->cd(2);
        gPad->SetLogy(0);
        SetHistoStyle(hratios, 1, 53, 1, 0.04 / pad2Size, 0.04 / pad2Size, 0.04 / pad2Size, 0.04 / pad2Size, 1.13, 1.8);
        hratios->GetYaxis()->SetTitle("FI / BC");
        hratios->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        hratios->GetYaxis()->SetTitleOffset(0.55);
        hratios->GetXaxis()->SetTitleOffset(1.1);
        hratios->GetYaxis()->SetRangeUser(0.9, 1.098);
        hratios->GetYaxis()->SetNdivisions(505);
        hratios->Draw("pe");
        TLine *line1 = new TLine(0.0, 1.0, 30.0, 1.0);
        line1->SetLineStyle(2);
        line1->SetLineWidth(2);
        line1->Draw("l same");

        if (saveplots)
        {
            c5->SaveAs(Form("saved/%s/gaussfit_Ks_yield_rebin%d.png", whichpass.c_str(), rebin));
        }
    }
}

void directfit(TF1 *&fit, TH1 *h, double *parameters, bool mainfit = false)
{
    fit = doubleCBpol2(h, parameters, mainfit);
    fit->Draw("SAME");
}

void perform_fit(TF1 *&fit, TH1 *h, bool mainfit = false)
{
    TF1 *fitpt1 = CB(h, ksmass, kswidth);
    double parameters1[5];
    for (int i = 0; i < 5; i++)
    {
        parameters1[i] = fitpt1->GetParameter(i);
    }
    TF1 *fitpt2 = doubleCB(h, parameters1, mainfit);
    double parameters3[7];
    for (int i = 0; i < 7; i++)
    {
        parameters3[i] = fitpt2->GetParameter(i);
    }
    fit = doubleCBpol2(h, parameters3, mainfit);
    fit->Draw("SAME");
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
    h->SetMarkerSize(0.8);
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
    h->Fit(fit, "REBMS0");
    // fit->Draw("SAME");
    return fit;
}

TF1 *CB(TH1 *h, double ksmass, double kswidth)
{
    TF1 *fit = new TF1("fit", CrystalBall, ksmass - 2.5 * kswidth, ksmass + 3.5 * kswidth, 5);
    fit->SetParNames("Norm", "Mean", "Sigma", "Alpha", "n");
    fit->SetParameter(0, 6.9e9); // if rebin == 3, then 5e8 else 8e7
    fit->SetParameter(1, ksmass);
    fit->SetParameter(2, kswidth);
    fit->SetParameter(3, 1.5);
    fit->SetParameter(4, 5.0);
    fit->SetLineColor(1);
    fit->SetLineWidth(2);
    h->Fit(fit, "REBMS0");
    // fit->Draw("SAME");
    return fit;
}

TF1 *CBpol2(TH1 *h, double *parameters, bool mainfit)
{
    double mass = parameters[1];
    double width = parameters[2];
    TF1 *fit = new TF1("fit", CrystalBallpol2, mass - 15 * width, mass + 15 * width, 8);
    fit->SetParameter(0, parameters[0]);
    fit->SetParameter(1, mass);
    fit->SetParameter(2, width);
    fit->SetParameter(3, parameters[3]);
    fit->SetParameter(4, parameters[4]);
    fit->SetLineColor(1);
    fit->SetLineWidth(2);
    h->Fit(fit, "REBMSQ0");
    // fit->Draw("SAME");

    // drawing separate fits
    if (mainfit)
    {
        TF1 *fitCB = new TF1("fitCB", CrystalBall, mass - 15 * width, mass + 15 * width, 5);
        TF1 *fitpol2 = new TF1("fitpol2", polynomial2, mass - 15 * width, mass + 15 * width, 3);
        for (int ipar = 0; ipar < 5; ipar++)
        {
            fitCB->SetParameter(ipar, fit->GetParameter(ipar));
            fitCB->SetLineColor(4);
            fitCB->SetLineStyle(2);
            fitCB->SetLineWidth(2);
        }
        fitpol2->SetParameter(0, fit->GetParameter(5));
        fitpol2->SetParameter(1, fit->GetParameter(6));
        fitpol2->SetParameter(2, fit->GetParameter(7));
        fitpol2->SetLineColor(28);
        fitpol2->SetLineStyle(2);
        fitpol2->SetLineWidth(2);
        fitCB->Draw("SAME");
        fitpol2->Draw("SAME");
        TLegend *leg = new TLegend(0.17, 0.75, 0.5, 0.93);
        SetLegendStyle(leg);
        leg->AddEntry(fit, "Crystal Ball + pol2 fit", "l");
        leg->AddEntry(fitCB, "Crystal Ball", "l");
        leg->AddEntry(fitpol2, "pol2", "l");
        leg->Draw("same");
    }
    return fit;
}

TF1 *doubleCB(TH1 *h, double *parameters, bool mainfit)
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
    h->Fit(fit, "REBMS0");
    // fit->Draw("SAME");

    // drawing separate fits
    if (mainfit)
    {
        TF1 *fitCBleft = new TF1("fitCBleft", CrystalBall, mass - 3.5 * width, mass + 4.0 * width, 5);
        TF1 *fitCBright = new TF1("fitCBright", CrystalBall, mass - 3.5 * width, mass + 4.0 * width, 5);
        for (int ipar = 0; ipar < 5; ipar++)
        {
            fitCBleft->SetParameter(ipar, fit->GetParameter(ipar));
            if (ipar < 3)
                fitCBright->SetParameter(ipar, fit->GetParameter(ipar));
            else
                fitCBright->SetParameter(ipar, fit->GetParameter(ipar + 2));
            fitCBleft->SetLineColor(4);
            fitCBright->SetLineColor(6);
            fitCBleft->SetLineStyle(2);
            fitCBright->SetLineStyle(2);
            fitCBleft->SetLineWidth(2);
            fitCBright->SetLineWidth(2);
        }
        fitCBleft->Draw("SAME");
        fitCBright->Draw("SAME");
        TLegend *leg = new TLegend(0.17, 0.75, 0.5, 0.93);
        SetLegendStyle(leg);
        leg->AddEntry(fit, "Double Crystal Ball fit", "l");
        leg->AddEntry(fitCBleft, "Left Crystal Ball", "l");
        leg->AddEntry(fitCBright, "Right Crystal Ball", "l");
        leg->Draw("same");
    }
    return fit;
}

TF1 *doubleCBpol2(TH1 *h, double *parameters, bool mainfit, TLegend *leg = nullptr, float legendsize = 0.04, float rlow = 15, float rhigh = 15)
{
    double mass = parameters[1];
    double width = parameters[2];
    TF1 *fit = new TF1("fit", DoubleCrystalBallpol2, mass - 15 * width, mass + 15 * width, 10);
    fit->SetParNames("Norm", "Mean", "Sigma", "AlphaL", "n1", "AlphaR", "n2", "p0", "p1", "p2");
    fit->SetParameter(0, parameters[0]);
    fit->SetParameter(1, mass);
    fit->SetParLimits(2, 0.0, 20.0);
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
    // fit->Draw("SAME");

    // drawing separate fits
    if (mainfit)
    {
        TF1 *fitCBleft = new TF1("fitCBleft", CrystalBall, mass - 30 * width, mass + 30 * width, 5);
        TF1 *fitCBright = new TF1("fitCBright", CrystalBall, mass - 30 * width, mass + 30 * width, 5);
        TF1 *fitpol2 = new TF1("fitpol2", polynomial2, mass - 30 * width, mass + 30 * width, 3);
        for (int ipar = 0; ipar < 5; ipar++)
        {
            fitCBleft->SetParameter(ipar, fit->GetParameter(ipar));
            if (ipar < 3)
                fitCBright->SetParameter(ipar, fit->GetParameter(ipar));
            else
                fitCBright->SetParameter(ipar, fit->GetParameter(ipar + 2));
            fitCBleft->SetLineColor(4);
            fitCBright->SetLineColor(6);
            fitCBleft->SetLineStyle(2);
            fitCBright->SetLineStyle(2);
            fitCBleft->SetLineWidth(2);
            fitCBright->SetLineWidth(2);
            if (ipar < 3)
            {
                fitpol2->SetParameter(ipar, fit->GetParameter(ipar + 7));
            }
            fitpol2->SetLineColor(28);
            fitpol2->SetLineStyle(2);
            fitpol2->SetLineWidth(2);
        }
        fitCBleft->Draw("SAME");
        fitCBright->Draw("SAME");
        fitpol2->Draw("SAME");

        // TLegend *leg = new TLegend(0.17, 0.75, 0.5, 0.93);
        if (leg != nullptr)
        {
            SetLegendStyle(leg);
            leg->SetTextSize(legendsize);
            leg->AddEntry(fit, "Double CB + pol2", "l");
            leg->AddEntry(fitCBleft, "Left CB", "l");
            leg->AddEntry(fitCBright, "Right CB", "l");
            leg->AddEntry(fitpol2, "Polynomial 2", "l");
            leg->Draw("same");
        }
    }
    return fit;
}

TF1 *doubleCBpol3(TH1 *h, double *parameters, bool mainfit, TLegend *leg = nullptr, float legendsize = 0.04, float rlow = 15, float rhigh = 15)
{
    double mass = parameters[1];
    double width = parameters[2];
    TF1 *fit = new TF1("fit", DoubleCrystalBallpol3, mass - rlow * width, mass + rhigh * width, 11);
    fit->SetParNames("Norm", "Mean", "Sigma", "AlphaL", "n1", "AlphaR", "n2", "p0", "p1", "p2", "p3");
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

    if (sizeof(parameters) == 11)
    {
        fit->SetParameter(7, parameters[7]);
        fit->SetParameter(8, parameters[8]);
        fit->SetParameter(9, parameters[9]);
        fit->SetParameter(10, parameters[10]);
    }
    fit->SetLineColor(2);
    fit->SetLineWidth(2);
    h->Fit(fit, "REBMS");
    fit->SetNpx(1e6);
    // fit->Draw("SAME");

    // drawing separate fits
    if (mainfit)
    {
        TF1 *fitCBleft = new TF1("fitCBleft", CrystalBall, mass - 30 * width, mass + 30 * width, 5);
        TF1 *fitCBright = new TF1("fitCBright", CrystalBall, mass - 30 * width, mass + 30 * width, 5);
        TF1 *fitpol3 = new TF1("fitpol3", polynomial3, mass - 30 * width, mass + 30 * width, 4);
        for (int ipar = 0; ipar < 5; ipar++)
        {
            fitCBleft->SetParameter(ipar, fit->GetParameter(ipar));
            if (ipar < 3)
                fitCBright->SetParameter(ipar, fit->GetParameter(ipar));
            else
                fitCBright->SetParameter(ipar, fit->GetParameter(ipar + 2));
            fitCBleft->SetLineColor(4);
            fitCBright->SetLineColor(6);
            fitCBleft->SetLineStyle(2);
            fitCBright->SetLineStyle(2);
            fitCBleft->SetLineWidth(2);
            fitCBright->SetLineWidth(2);
            if (ipar < 4)
            {
                fitpol3->SetParameter(ipar, fit->GetParameter(ipar + 7));
            }
            fitpol3->SetLineColor(28);
            fitpol3->SetLineStyle(2);
            fitpol3->SetLineWidth(2);
        }
        fitCBleft->Draw("SAME");
        fitCBright->Draw("SAME");
        fitpol3->Draw("SAME");

        // TLegend *leg = new TLegend(0.17, 0.75, 0.5, 0.93);
        if (leg != nullptr)
        {
            SetLegendStyle(leg);
            leg->SetTextSize(legendsize);
            leg->AddEntry(fit, "Double CB + pol2", "l");
            leg->AddEntry(fitCBleft, "Left CB", "l");
            leg->AddEntry(fitCBright, "Right CB", "l");
            leg->AddEntry(fitpol3, "Polynomial 3", "l");
            leg->Draw("same");
        }
    }
    return fit;
}

TF1 *doubleCBpol1(TH1 *h, double *parameters, bool mainfit)
{
    double mass = parameters[1];
    double width = parameters[2];
    TF1 *fit = new TF1("fit", DoubleCrystalBallpol1, mass - 15 * width, mass + 15 * width, 9);
    fit->SetParNames("Norm", "Mean", "Sigma", "AlphaL", "n1", "AlphaR", "n2", "p0", "p1");
    fit->SetParameter(0, parameters[0]);
    fit->SetParameter(1, mass);
    // fit->SetParLimits(2, 0.0, 20.0);
    fit->SetParameter(2, width);
    // fit->SetParLimits(3, 0.0, 10.0);
    fit->SetParameter(3, parameters[3]);
    // fit->SetParLimits(4, 0.0, 10.0);
    fit->SetParameter(4, parameters[4]);
    // fit->SetParLimits(5, 0.0, 10.0);
    fit->SetParameter(5, parameters[5]);
    // fit->SetParLimits(6, 0.0, 10.0);
    fit->SetParameter(6, parameters[6]);
    if (sizeof(parameters) == 9)
    {
        fit->SetParameter(7, parameters[7]);
        fit->SetParameter(8, parameters[8]);
    }
    fit->SetLineColor(2);
    fit->SetLineWidth(2);
    h->Fit(fit, "REBMS0");
    fit->SetNpx(1e6);
    // fit->Draw("SAME");

    // drawing separate fits
    if (mainfit)
    {
        TF1 *fitCBleft = new TF1("fitCBleft", CrystalBall, mass - 30 * width, mass + 30 * width, 5);
        TF1 *fitCBright = new TF1("fitCBright", CrystalBall, mass - 30 * width, mass + 30 * width, 5);
        TF1 *fitpol1 = new TF1("fitpol1", polynomial1, mass - 30 * width, mass + 30 * width, 2);
        for (int ipar = 0; ipar < 5; ipar++)
        {
            fitCBleft->SetParameter(ipar, fit->GetParameter(ipar));
            if (ipar < 3)
                fitCBright->SetParameter(ipar, fit->GetParameter(ipar));
            else
                fitCBright->SetParameter(ipar, fit->GetParameter(ipar + 2));
            fitCBleft->SetLineColor(4);
            fitCBright->SetLineColor(6);
            fitCBleft->SetLineStyle(2);
            fitCBright->SetLineStyle(3);
            fitCBleft->SetLineWidth(2);
            fitCBright->SetLineWidth(2);
            if (ipar < 2)
            {
                fitpol1->SetParameter(ipar, fit->GetParameter(ipar + 7));
            }
            fitpol1->SetLineColor(28);
            fitpol1->SetLineStyle(4);
            fitpol1->SetLineWidth(2);
        }
        fitCBleft->Draw("SAME");
        fitCBright->Draw("SAME");
        fitpol1->Draw("SAME");

        TLegend *leg = new TLegend(0.17, 0.75, 0.5, 0.93);
        SetLegendStyle(leg);
        leg->AddEntry(fit, "Double Crystal Ball + pol1 fit", "l");
        leg->AddEntry(fitCBleft, "Left Crystal Ball", "l");
        leg->AddEntry(fitCBright, "Right Crystal Ball", "l");
        leg->AddEntry(fitpol1, "Polynomial 1", "l");
        leg->Draw("same");
    }
    return fit;
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.02);
    pad2->SetRightMargin(0.02);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.14);
    pad2->SetLeftMargin(0.14);
    pad1->SetTopMargin(0.08);
    pad1->SetBottomMargin(0.0001);
    pad2->SetTopMargin(0.001);
}
