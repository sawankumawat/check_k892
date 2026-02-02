#include <iostream>
#include "../macro/src/style.h"
#include "../macro/src/fitfunc.h"
#include "../macro/src/initializations.h"
using namespace std;

void normalize_hist(TH1F *hdef, TH1F *hbkg, double norm_range_low, double norm_range_high)
{
    double integral_def = hdef->Integral(hdef->FindBin(norm_range_low), hdef->FindBin(norm_range_high));
    double integral_bkg = hbkg->Integral(hbkg->FindBin(norm_range_low), hbkg->FindBin(norm_range_high));
    double ratio = integral_def / integral_bkg;
    hbkg->Scale(ratio);
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);

void different_rot_kstar()
{
    gStyle->SetOptStat(0);
    int rebin = 2;
    TFile *fkstar = new TFile("kstar_different_rotBkg.root");
    if (fkstar->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TH1F *hmult = (TH1F *)fkstar->Get("mult_dist");
    TH1F *hkstar_unlike = (TH1F *)fkstar->Get("kstar_ptspectra");
    TH1F *hkstar_rot1 = (TH1F *)fkstar->Get("kstar_rotDef");  // angle different pi/10
    TH1F *hkstar_rot2 = (TH1F *)fkstar->Get("kstar_rotLow");  // angle different pi/15
    TH1F *hkstar_rot3 = (TH1F *)fkstar->Get("kstar_rotHigh"); // angle different pi/8
    TH1F *hkstar_like = (TH1F *)fkstar->Get("kstar_like");

    if (hmult == nullptr || hkstar_unlike == nullptr || hkstar_rot1 == nullptr || hkstar_rot2 == nullptr || hkstar_rot3 == nullptr || hkstar_like == nullptr)
    {
        cout << "Error reading histograms" << endl;
        return;
    }

    float norm_range_low = 1.2;
    float norm_range_high = 1.3;
    normalize_hist(hkstar_unlike, hkstar_rot1, norm_range_low, norm_range_high);
    normalize_hist(hkstar_unlike, hkstar_rot2, norm_range_low, norm_range_high);
    normalize_hist(hkstar_unlike, hkstar_rot3, norm_range_low, norm_range_high);

    TCanvas *c = new TCanvas("", "", 720, 720);
    double pad1Size, pad2Size;
    canvas_style(c, pad1Size, pad2Size);
    c->cd(1);
    SetHistoQA(hkstar_unlike);
    hkstar_unlike->SetMarkerStyle(20);
    hkstar_unlike->SetMarkerSize(0.5);
    hkstar_unlike->SetTitle(0);
    hkstar_unlike->GetXaxis()->SetTitle("M_{K^{#pm}#pi^{#mp}} [GeV/c^{2}]");
    hkstar_unlike->GetYaxis()->SetTitle(Form("Counts / %.2f MeV/c^{2}", hkstar_unlike->GetBinWidth(1) * 1000));
    hkstar_unlike->GetXaxis()->SetRangeUser(0.6, 1.35);
    hkstar_unlike->GetXaxis()->SetTitleSize(0.04 / pad1Size);
    hkstar_unlike->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    hkstar_unlike->GetXaxis()->SetLabelSize(0.04 / pad1Size);
    hkstar_unlike->GetYaxis()->SetLabelSize(0.04 / pad1Size);
    hkstar_unlike->GetYaxis()->SetTitleOffset(1.2);
    hkstar_unlike->SetMinimum(30e3);
    hkstar_unlike->Draw("pe");
    hkstar_rot1->SetMarkerStyle(21);
    hkstar_rot1->SetMarkerSize(0.5);
    hkstar_rot1->SetLineColor(kRed);
    hkstar_rot1->SetMarkerColor(kRed);
    hkstar_rot1->Draw("pe same");
    hkstar_rot2->SetMarkerStyle(22);
    hkstar_rot2->SetMarkerSize(0.5);
    hkstar_rot2->SetLineColor(kBlue);
    hkstar_rot2->SetMarkerColor(kBlue);
    hkstar_rot2->Draw("pe same");
    hkstar_rot3->SetMarkerStyle(23);
    hkstar_rot3->SetMarkerSize(0.5);
    hkstar_rot3->SetLineColor(kGreen);
    hkstar_rot3->SetMarkerColor(kGreen);
    hkstar_rot3->Draw("pe same");
    hkstar_like->SetMarkerStyle(24);
    hkstar_like->SetMarkerSize(0.5);
    hkstar_like->SetLineColor(kMagenta);
    hkstar_like->SetMarkerColor(kMagenta);
    hkstar_like->Draw("pe same");

    TLegend *leg = new TLegend(0.6, 0.78, 0.92, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04 / pad1Size);
    leg->AddEntry((TObject *)0, "Pythia simulations", "");
    leg->AddEntry((TObject *)0, "K*^{0} (892) spectra", "");
    leg->Draw();

    TLegend *leg2 = new TLegend(0.30, 0.20, 0.65, 0.60);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.04 / pad1Size);
    leg2->AddEntry(hkstar_unlike, "Unlike sign", "pel");
    leg2->AddEntry(hkstar_like, "Like sign", "pel");
    leg2->AddEntry(hkstar_rot1, "Rotated #pi #pm #pi/10", "pel");
    leg2->AddEntry(hkstar_rot2, "Rotated #pi #pm #pi/15", "pel");
    leg2->AddEntry(hkstar_rot3, "Rotated #pi #pm #pi/8", "pel");
    leg2->Draw();

    TH1F *hratio_rot1 = (TH1F *)hkstar_like->Clone();
    hratio_rot1->Divide(hkstar_rot1);
    TH1F *hratio_rot2 = (TH1F *)hkstar_like->Clone();
    hratio_rot2->Divide(hkstar_rot2);
    TH1F *hratio_rot3 = (TH1F *)hkstar_like->Clone();
    hratio_rot3->Divide(hkstar_rot3);

    c->cd(2);
    SetHistoQA(hratio_rot1);
    hratio_rot1->SetTitle(0);
    hratio_rot1->GetXaxis()->SetTitle("M_{K^{#pm}#pi^{#mp}} [GeV/c^{2}]");
    hratio_rot1->GetYaxis()->SetTitle("(Rot. / Like)");
    hratio_rot1->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    hratio_rot1->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    hratio_rot1->GetXaxis()->SetLabelSize(0.04 / pad2Size);
    hratio_rot1->GetYaxis()->SetLabelSize(0.04 / pad2Size);
    hratio_rot1->GetXaxis()->SetRangeUser(0.6, 1.35);
    hratio_rot1->GetYaxis()->SetRangeUser(0.93, 1.07);
    hratio_rot1->GetYaxis()->SetNdivisions(505);
    hratio_rot1->GetYaxis()->SetTitleOffset(0.5);
    hratio_rot1->Draw("pe");
    hratio_rot2->SetMarkerStyle(22);
    hratio_rot2->SetMarkerSize(0.5);
    hratio_rot2->SetLineColor(kBlue);
    hratio_rot2->SetMarkerColor(kBlue);
    hratio_rot2->Draw("pe same");
    hratio_rot3->SetMarkerStyle(23);
    hratio_rot3->SetMarkerSize(0.5);
    hratio_rot3->SetLineColor(kGreen);
    hratio_rot3->SetMarkerColor(kGreen);
    hratio_rot3->Draw("pe same");
    TLine *line = new TLine(0.6, 1.0, 1.35, 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kBlue);
    line->Draw();
    c->SaveAs("kstar_different_rotBkg.pdf");

    TH1F *invsub_like = (TH1F *)hkstar_unlike->Clone();
    invsub_like->Add(hkstar_like, -1);
    TH1F *invsub_rot1 = (TH1F *)hkstar_unlike->Clone();
    invsub_rot1->Add(hkstar_rot1, -1);
    TH1F *invsub_rot2 = (TH1F *)hkstar_unlike->Clone();
    invsub_rot2->Add(hkstar_rot2, -1);
    TH1F *invsub_rot3 = (TH1F *)hkstar_unlike->Clone();
    invsub_rot3->Add(hkstar_rot3, -1);

    // TH1F *invsub_selected = (TH1F *)invsub_like->Clone(); // selected like-sign combinatorial background
    TH1F *invsub_selected = (TH1F *)invsub_rot2->Clone(); // selected rot 1 combinatorial background

    invsub_selected->Rebin(rebin);
    invsub_selected->GetYaxis()->SetTitle(Form("Counts / %.1f MeV/c^{2}", invsub_selected->GetBinWidth(1) * 1000));
    TCanvas *c2 = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c2, 0.15, 0.005, 0.05, 0.15);
    SetHistoQA(invsub_selected);
    invsub_selected->GetYaxis()->SetMaxDigits(3);
    invsub_selected->Draw("pe");

    TF1 *fitFcn, *fitFcn1;
    fitFcn = new TF1("fitfunc", BreitWignerpoly3, 0.69, 1.3, 7);
    fitFcn1 = new TF1("fitfunc1", polynomial3, 0.69, 1.3, 4);
    TF1 *fitFcn2 = new TF1("fitFcn2", BW, 0.65, 1.3, 3); // only signal

    // fitFcn->SetParLimits(0, 0.80, 0.98); // Mass
    fitFcn->SetParameter(0, 0.895);   // Mass
    fitFcn->SetParLimits(2, 0, 10e9); // Yield
    fitFcn->SetParameter(1, 0.047);   // width
    fitFcn->SetParNames("Mass", "Width", "Yield", "A", "B", "C", "D");
    invsub_selected->Fit(fitFcn, "REBMSQ+");
    Double_t *par = fitFcn->GetParameters();
    fitFcn2->SetParameters(&par[0]);
    fitFcn1->SetParameters(&par[3]);
    fitFcn1->SetLineColor(4);
    fitFcn1->SetLineStyle(2);
    fitFcn1->SetLineWidth(4);
    fitFcn2->SetLineColor(6);
    fitFcn2->SetLineStyle(2);
    fitFcn2->SetLineWidth(4);
    fitFcn->SetLineWidth(4);
    fitFcn->SetLineWidth(2);
    fitFcn1->SetLineWidth(2);
    fitFcn2->SetLineWidth(2);
    fitFcn->Draw("same");
    fitFcn1->Draw("same");
    fitFcn2->Draw("same");

    TLegend *pag2 = new TLegend(0.2, 0.7, 0.45, 0.9);
    pag2->SetBorderSize(0);
    pag2->SetTextFont(42);
    pag2->SetTextSize(0.04);
    pag2->SetFillStyle(0);
    pag2->AddEntry(fitFcn, "BW+pol3");
    pag2->AddEntry(fitFcn1, "BW");
    pag2->AddEntry(fitFcn2, "pol3");
    pag2->Draw();

    // significance of signal
    int bmin = invsub_selected->GetXaxis()->FindBin(masspdg - 2 * widthpdg);
    int bmax = invsub_selected->GetXaxis()->FindBin(masspdg + 2 * widthpdg);
    auto binwidth_file = (hkstar_unlike->GetXaxis()->GetXmax() - hkstar_unlike->GetXaxis()->GetXmin()) * rebin / hkstar_unlike->GetXaxis()->GetNbins();

    double significance_den = TMath::Sqrt(hkstar_unlike->Integral(bmin, bmax));
    double significance_num = (fitFcn2->Integral(masspdg - 2 * widthpdg, masspdg + 2 * widthpdg)) / (binwidth_file);
    double ratio = significance_num / significance_den; // significance of signal
    double chi2ndf = (fitFcn->GetChisquare()) / (fitFcn->GetNDF());
    double fitmean = fitFcn->GetParameter(0);
    double fitmeanerror = fitFcn->GetParError(0);
    double fitwidth = fitFcn->GetParameter(1);
    double fitwidtherror = fitFcn->GetParError(1);

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.04);
    lat.SetTextFont(42);
    lat.DrawLatex(0.57, 0.85, Form("#Chi^{2}/NDF: %.1f", chi2ndf));
    // lat.DrawLatex(0.57, 0.85, Form("#Chi^{2}/NDF: %.1f/%d", fitFcn->GetChisquare(), fitFcn->GetNDF()));

    lat.DrawLatex(0.57, 0.8, Form("Significance: %.2f", ratio));
    lat.DrawLatex(0.57, 0.75, Form("Mass: %.3f #pm %.5f", fitmean, fitmeanerror));
    lat.DrawLatex(0.57, 0.7, Form("Width: %.3f #pm %.5f", fitwidth, fitwidtherror));
    lat.SetTextSize(0.035);
    lat.SetTextColor(kRed+3);
    // lat.DrawLatex(0.50, 0.6, "Unlike sign - Like sign");
    // lat.DrawLatex(0.50, 0.6, "Unlike sign - Rotated (#pi #pm #pi/10)");
    lat.DrawLatex(0.50, 0.6, "Unlike sign - Rotated (#pi #pm #pi/15)");
    // lat.DrawLatex(0.50, 0.6, "Unlike sign - Rotated (#pi #pm #pi/8)");

    c2->SaveAs("fit_rot2.png");
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    // SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.02);
    pad2->SetRightMargin(0.02);
    pad2->SetBottomMargin(0.35);
    pad1->SetLeftMargin(0.145);
    pad2->SetLeftMargin(0.145);
    pad1->SetTopMargin(0.06);
    pad1->SetBottomMargin(0.003);
    pad2->SetTopMargin(0.01);
}