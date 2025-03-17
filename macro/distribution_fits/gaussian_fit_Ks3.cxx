#include <iostream>
#include "TDatabasePDG.h"
// #include "../src/fitfunc.h"
#include "../src/common_glue.h"
#include "../src/style.h"

TDatabasePDG *pdg = new TDatabasePDG();
Double_t fitgaus(Double_t *x, Double_t *par);
Double_t fitgauspol2(Double_t *x, Double_t *par);
Double_t CBpol2(Double_t *x, Double_t *par);
Double_t CBpol1(Double_t *x, Double_t *par);
Double_t CB(Double_t *x, Double_t *par);
Double_t doubleCB(Double_t *x, Double_t *par);
Double_t doubleCBpol2(Double_t *x, Double_t *par);
Double_t doubleCBpol1(Double_t *x, Double_t *par);
void SetHistoStyle_temp(TH1 *h, Int_t MCol, Int_t MSty, double binwidth);
void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void gaussian_fit_Ks3()
{

    // gStyle->SetOptStat(1110);
    gStyle->SetOptStat(0);
    // gStyle->SetOptFit(1111);
    gStyle->SetOptFit(0);

    int rebin = 1;
    // configurables *********************

    double ksmass = pdg->GetParticle(310)->Mass();
    double kswidth = 0.004;

    TLatex *t2 = new TLatex();
    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.05);
    t2->SetTextFont(42);

    //********************************************************************

    string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/351471/KsKs_Channel/higher-mass-resonances_100SigmaKs/QA/KsInvMass.root";
    TFile *f = new TFile(path.c_str(), "read");
    if (f == nullptr)
    {
        cout << "File not found" << endl;
        return;
    }
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
    SetCanvasStyle(c1, 0.14, 0.03, 0.05, 0.13);
    SetHistoStyle_temp(hInvMass, 1, 20, binwidth);
    hInvMass->GetYaxis()->SetMaxDigits(3);
    hInvMass->GetYaxis()->SetTitleOffset(1.5);
    hInvMass->SetMinimum(-100);
    hInvMass->SetMarkerStyle(20);
    hInvMass->SetMarkerSize(1.0);
    hInvMass->Draw("pe");
    auto noofevents = hInvMass->Integral();

    TF1 *gauspol2 = new TF1("gauspol2", fitgauspol2, 0.45, 0.55, 6);
    gauspol2->SetParameter(0, 1e9);
    gauspol2->SetParameter(1, ksmass);
    gauspol2->SetParameter(2, kswidth);
    gauspol2->SetParameter(3, 1);
    gauspol2->SetParameter(4, 1);
    gauspol2->SetParameter(5, 1);
    // hInvMass->Fit("gauspol2", "R");

    TF1 *crystalball = new TF1("crystalball", CB, 0.487, 0.507, 5);
    crystalball->SetParameter(0, 1e9);
    crystalball->SetParameter(1, ksmass);
    crystalball->SetParameter(2, kswidth);
    crystalball->SetParameter(3, 1);
    crystalball->SetParameter(4, 1);
    hInvMass->Fit("crystalball", "REBMS0");

    TF1 *crystalballpol1 = new TF1("crystalballpol1", CBpol1, 0.46, 0.54, 7);
    crystalballpol1->SetParameter(0, crystalball->GetParameter(0));
    crystalballpol1->SetParameter(1, crystalball->GetParameter(1));
    crystalballpol1->SetParameter(2, crystalball->GetParameter(2));
    crystalballpol1->SetParameter(3, crystalball->GetParameter(3));
    crystalballpol1->SetParameter(4, crystalball->GetParameter(4));
    crystalballpol1->SetParameter(5, 5e7);
    crystalballpol1->SetParameter(6, 6e7);
    hInvMass->Fit("crystalballpol1", "REBMS0");

    TF1 *doublecrystalBallpol1 = new TF1("doublecrystalBallpol1", doubleCBpol1, 0.46, 0.54, 9);
    doublecrystalBallpol1->SetParameter(0, crystalballpol1->GetParameter(0));
    doublecrystalBallpol1->SetParameter(1, crystalballpol1->GetParameter(1));
    doublecrystalBallpol1->SetParameter(2, crystalballpol1->GetParameter(2));
    doublecrystalBallpol1->SetParameter(3, crystalballpol1->GetParameter(3));
    doublecrystalBallpol1->SetParameter(4, crystalballpol1->GetParameter(4));
    doublecrystalBallpol1->SetParameter(5, crystalballpol1->GetParameter(3));
    doublecrystalBallpol1->SetParameter(6, crystalballpol1->GetParameter(4));
    doublecrystalBallpol1->SetParameter(7, crystalballpol1->GetParameter(5));
    doublecrystalBallpol1->SetParameter(8, crystalballpol1->GetParameter(6));
    hInvMass->Fit("doublecrystalBallpol1", "REBMS");

    double mean = doublecrystalBallpol1->GetParameter(1);
    double width = doublecrystalBallpol1->GetParameter(2);

    TH1F *hInvMassClone1 = (TH1F *)hInvMass->Clone("hInvMassClone1");
    TH1F *hInvMassClone2 = (TH1F *)hInvMass->Clone("hInvMassClone2");
    hInvMassClone1->SetLineColor(0);
    hInvMassClone1->SetLineWidth(0);
    hInvMassClone1->SetFillColor(5);
    hInvMassClone1->GetXaxis()->SetRangeUser(mean - 3.0 * width, mean + 3.0 * width);
    hInvMassClone1->Draw("E3 hist same");
    hInvMassClone2->Draw("pe same");

    cout << "mass is " << mean * 1000 << " MeV" << " #pm " << doublecrystalBallpol1->GetParError(1) * 1000 << " MeV" << endl;
    cout << "width is " << width * 1000 << " MeV" << " #pm " << doublecrystalBallpol1->GetParError(2) * 1000 << " MeV" << endl;

    TF1 *pol1 = new TF1("pol1", "[0]+x[0]*[1]", doublecrystalBallpol1->GetXmin(), doublecrystalBallpol1->GetXmax());
    pol1->SetParameter(0, doublecrystalBallpol1->GetParameter(7));
    pol1->SetParameter(1, doublecrystalBallpol1->GetParameter(8));
    pol1->SetLineColor(4);
    pol1->SetLineStyle(2);
    pol1->SetLineWidth(2);
    pol1->Draw("SAME");

    TF1 *doublecrystalBall = new TF1("doublecrystalBall", doubleCB, doublecrystalBallpol1->GetXmin(), doublecrystalBallpol1->GetXmax(), 7);
    for (int i = 0; i < 7; i++)
    {
        doublecrystalBall->SetParameter(i, doublecrystalBallpol1->GetParameter(i));
    }
    doublecrystalBall->SetLineColor(6);
    doublecrystalBall->SetLineWidth(2);
    doublecrystalBall->SetLineStyle(2);
    // doublecrystalBall->Draw("SAME");

    double areabkg = pol1->Integral(mean - 3 * width, mean + 3 * width) / binwidth;
    double areasigbkg = hInvMass->Integral(hInvMass->FindBin(mean - 3 * width), hInvMass->FindBin(mean + 3 * width));
    cout << "background area is " << areabkg << endl;
    cout << "signal+bkg area is " << areasigbkg << endl;
    double areasig = areasigbkg - areabkg;
    cout<<"Signal area: "<<areasig<<endl;
    double purity = areasig * 100 / areasigbkg;
    cout << "Purity: " << purity << endl;
    cout<<"Total counts percentage in the signal region: "<<areasigbkg*100/noofevents<<endl;

    TLegend *lp2 = DrawLegend(0.16, 0.62, 0.42, 0.91);
    lp2->SetTextSize(0.035);
    lp2->SetTextFont(42);
    lp2->SetFillStyle(0);
    lp2->AddEntry((TObject *)0, "pp, #sqrt{#it{s}} = 13.6 TeV", "");
    lp2->AddEntry((TObject *)0, "FT0M, 0-100%", "");
    lp2->AddEntry((TObject *)0, "|#it{y}| < 0.5", "");
    lp2->AddEntry((TObject *)0, "p_{T}: 0-30 GeV/#it{c}", "");
    lp2->Draw("same");

    TLegend *lp3 = DrawLegend(0.6, 0.72, 0.9, 0.9);
    lp3->SetFillStyle(0);
    lp3->SetTextSize(0.035);
    lp3->SetTextFont(42);
    lp3->AddEntry(hInvMass, "Data (K^{0}_{s} #rightarrow #pi^{+}#pi^{-})", "lpe");
    lp3->AddEntry(doublecrystalBallpol1, "Fit", "l");
    lp3->AddEntry(pol1, "Pol 1", "l");
    lp3->Draw("same");

    TLegend *lp4 = DrawLegend(0.48, 0.5, 0.9, 0.65);
    lp4->SetFillStyle(0);
    lp4->SetBorderSize(0);
    lp4->SetTextFont(42);
    lp4->SetTextSize(0.03);
    lp4->AddEntry((TObject *)0, Form("Mass = %.2f #pm %.1e MeV", doublecrystalBallpol1->GetParameter(1) * 1000, doublecrystalBallpol1->GetParError(1) * 1000), "");
    lp4->AddEntry((TObject *)0, Form("Width = %.2f #pm %.1e MeV", doublecrystalBallpol1->GetParameter(2) * 1000, doublecrystalBallpol1->GetParError(2) * 1000), "");
    lp4->Draw("same");

    c1->SaveAs("saved/gaussfit_Ks3.pdf");

    // // gPad->Update();
    // // TPaveStats *st = (TPaveStats *)hInvMass->FindObject("stats");
    // // st->SetX1NDC(0.6); // 0.60
    // // st->SetX2NDC(0.95);
    // // st->SetY1NDC(0.3); // 0.78
    // // st->SetY2NDC(0.95);
    // // st->Draw("same");
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
    // h->GetXaxis()->SetRangeUser(0.45, 0.55);
    h->GetYaxis()->SetTitle(Form("Counts /(%.3f GeV/#it{c^{2}})", binwidth));
    h->GetXaxis()->SetTitle("#it{M}_{#pi^{+}#pi^{-}} (GeV/#it{c^{2}})");
}

Double_t fitgaus(Double_t *x, Double_t *par)
{
    double gauss = par[0] * TMath::Exp(-0.5 * ((x[0] - par[1]) / par[2]) * ((x[0] - par[1]) / par[2]));
    return gauss;
}

Double_t fitgauspol2(Double_t *x, Double_t *par)
{
    double gauss = par[0] * TMath::Exp(-0.5 * ((x[0] - par[1]) / par[2]) * ((x[0] - par[1]) / par[2]));
    double pol2 = par[3] + par[4] * x[0] + par[5] * x[0] * x[0];
    return gauss + pol2;
}

Double_t CB(Double_t *x, Double_t *par)
{
    // par[0] normalization
    // par[1] mean of gaussian
    // par[2] sigma of gaussian
    // par[3] alpha
    // par[4] n

    double t = (x[0] - par[1]) / par[2];
    double absAlpha_L = fabs(par[3]);
    double n = par[4];
    double y1 = 0;

    if (t < -absAlpha_L)
    {
        double a = exp(-0.5 * absAlpha_L * absAlpha_L) * TMath::Power(n / absAlpha_L, n);
        double b = (n / absAlpha_L) - absAlpha_L;
        y1 = par[0] * (a / TMath::Power(b - t, n));
    }
    else if (t >= -absAlpha_L)
    {
        y1 = par[0] * exp(-0.5 * t * t);
    }
    return y1;
}

Double_t CBpol1(Double_t *x, Double_t *par)
{
    // par[0] normalization
    // par[1] mean of gaussian
    // par[2] sigma of gaussian
    // par[3] alpha
    // par[4] n
    // par[5] p0
    // par[6] p1
    // par[7] p2

    double t = (x[0] - par[1]) / par[2];
    double absAlpha_L = fabs(par[3]);
    double n = par[4];
    double y1 = 0;
    double pol1 = par[5] + par[6] * x[0];

    if (t < -absAlpha_L)
    {
        double a = exp(-0.5 * absAlpha_L * absAlpha_L) * TMath::Power(n / absAlpha_L, n);
        double b = (n / absAlpha_L) - absAlpha_L;
        y1 = par[0] * (a / TMath::Power(b - t, n));
    }
    else if (t >= -absAlpha_L)
    {
        y1 = par[0] * exp(-0.5 * t * t);
    }
    return y1 + pol1;
}

Double_t CBpol2(Double_t *x, Double_t *par)
{
    // par[0] normalization
    // par[1] mean of gaussian
    // par[2] sigma of gaussian
    // par[3] alpha
    // par[4] n
    // par[5] p0
    // par[6] p1
    // par[7] p2

    double t = (x[0] - par[1]) / par[2];
    double absAlpha_L = fabs(par[3]);
    double n = par[4];
    double y1 = 0;
    double pol2 = par[5] + par[6] * x[0] + par[7] * x[0] * x[0];

    if (t < -absAlpha_L)
    {
        double a = exp(-0.5 * absAlpha_L * absAlpha_L) * TMath::Power(n / absAlpha_L, n);
        double b = (n / absAlpha_L) - absAlpha_L;
        y1 = par[0] * (a / TMath::Power(b - t, n));
    }
    else if (t >= -absAlpha_L)
    {
        y1 = par[0] * exp(-0.5 * t * t);
    }
    return y1 + pol2;
}

Double_t doubleCB(Double_t *x, Double_t *par)
{
    // par[0] normalization
    // par[1] mean of gaussian
    // par[2] sigma of gaussian
    // par[3] alpha left
    // par[4] n1
    // par[5] alpha right
    // par[6] n2

    double t = (x[0] - par[1]) / par[2];
    double absAlpha_L = fabs(par[3]);
    double n1 = par[4];
    double absAlpha_R = fabs(par[5]);
    double n2 = par[6];
    double y1 = 0;

    if ((t >= -absAlpha_L) && (t < absAlpha_R))
    {
        y1 = par[0] * exp(-0.5 * t * t);
    }
    else if (t < -absAlpha_L)
    {
        double a = exp(-0.5 * absAlpha_L * absAlpha_L) * TMath::Power(n1 / absAlpha_L, n1);
        double b = n1 / absAlpha_L - absAlpha_L;
        y1 = par[0] * (a / TMath::Power(b - t, n1));
    }
    else if (t >= absAlpha_R)
    {
        double a = exp(-0.5 * absAlpha_R * absAlpha_R) * TMath::Power(n2 / absAlpha_R, n2);
        double b = n2 / absAlpha_R - absAlpha_R;
        y1 = par[0] * (a / TMath::Power(b + t, n2));
    }

    return y1;
}

Double_t doubleCBpol1(Double_t *x, Double_t *par)
{
    // par[0] normalization
    // par[1] mean of gaussian
    // par[2] sigma of gaussian
    // par[3] alpha left
    // par[4] n1
    // par[5] alpha right
    // par[6] n2
    // par[7] p0
    // par[8] p1
    // par[9] p2

    double t = (x[0] - par[1]) / par[2];
    double absAlpha_L = fabs(par[3]);
    double n1 = par[4];
    double absAlpha_R = fabs(par[5]);
    double n2 = par[6];
    double y1 = 0;
    double pol2 = par[7] + par[8] * x[0] + par[9] * x[0] * x[0];
    double pol1 = par[7] + par[8] * x[0];

    if ((t >= -absAlpha_L) && (t < absAlpha_R))
    {
        y1 = par[0] * exp(-0.5 * t * t);
    }
    else if (t < -absAlpha_L)
    {
        double a = exp(-0.5 * absAlpha_L * absAlpha_L) * TMath::Power(n1 / absAlpha_L, n1);
        double b = n1 / absAlpha_L - absAlpha_L;
        y1 = par[0] * (a / TMath::Power(b - t, n1));
    }
    else if (t >= absAlpha_R)
    {
        double a = exp(-0.5 * absAlpha_R * absAlpha_R) * TMath::Power(n2 / absAlpha_R, n2);
        double b = n2 / absAlpha_R - absAlpha_R;
        y1 = par[0] * (a / TMath::Power(b + t, n2));
    }

    return y1 + pol1;
}

Double_t doubleCBpol2(Double_t *x, Double_t *par)
{
    // par[0] normalization
    // par[1] mean of gaussian
    // par[2] sigma of gaussian
    // par[3] alpha left
    // par[4] n1
    // par[5] alpha right
    // par[6] n2
    // par[7] p0
    // par[8] p1
    // par[9] p2

    double t = (x[0] - par[1]) / par[2];
    double absAlpha_L = fabs(par[3]);
    double n1 = par[4];
    double absAlpha_R = fabs(par[5]);
    double n2 = par[6];
    double y1 = 0;
    double pol2 = par[7] + par[8] * x[0] + par[9] * x[0] * x[0];
    double pol1 = par[7] + par[8] * x[0];

    if ((t >= -absAlpha_L) && (t < absAlpha_R))
    {
        y1 = par[0] * exp(-0.5 * t * t);
    }
    else if (t < -absAlpha_L)
    {
        double a = exp(-0.5 * absAlpha_L * absAlpha_L) * TMath::Power(n1 / absAlpha_L, n1);
        double b = n1 / absAlpha_L - absAlpha_L;
        y1 = par[0] * (a / TMath::Power(b - t, n1));
    }
    else if (t >= absAlpha_R)
    {
        double a = exp(-0.5 * absAlpha_R * absAlpha_R) * TMath::Power(n2 / absAlpha_R, n2);
        double b = n2 / absAlpha_R - absAlpha_R;
        y1 = par[0] * (a / TMath::Power(b + t, n2));
    }

    return y1 + pol2;
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
