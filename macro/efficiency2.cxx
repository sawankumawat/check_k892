#include <iostream>
#include "src/style.h"
#include "src/initializations.h"

using namespace std;

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
    pad1->SetRightMargin(0.009);
    pad2->SetRightMargin(0.009);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.14);
    pad2->SetLeftMargin(0.14);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.002);
    pad2->SetTopMargin(0.04);
}

void efficiency2()
{
    string sysvarname = "_bypassTOF";
    string sysvarname2 = "_bypass_TOF";
    // string sysvarname = "";
    // string sysvarname2 = "";
    TFile *fileraw = new TFile("/home/sawan/k892_postprocessing/output/LHC22o_pass6_small_INEL/common/signal_2.root", "READ");
    TFile *fileeff = new TFile("/home/sawan/k892_postprocessing/mc/LHC24b1/204671.root", "READ");

    if (fileeff->IsZombie() || fileraw->IsZombie())
    {
        cout << "Error opening files" << endl;
        return;
    }
    TH3F *h2gen = (TH3F *)fileeff->Get("lf-k892analysis_bypass_TOF/k892Gen");
    TH3F *h2genanti = (TH3F *)fileeff->Get("lf-k892analysis_bypass_TOF/k892GenAnti");
    TH2F *h2rec = (TH2F *)fileeff->Get("lf-k892analysis_bypass_TOF/k892Rec");
    TH2F *h2recanti = (TH2F *)fileeff->Get("lf-k892analysis_bypass_TOF/k892RecAnti");
    if (h2rec == nullptr || h2recanti == nullptr || h2gen == nullptr || h2genanti == nullptr)
    {
        cout << "Error reading histograms" << endl;
        return;
    }

    TH1D *h1gen = h2gen->ProjectionY("h1gen", h2gen->GetXaxis()->FindBin(0.0 - 0.001), h2gen->GetXaxis()->FindBin(0.0 + 0.001), -1, -1, "E");
    TH1D *h1genanti = h2genanti->ProjectionY("h1genanti", h2genanti->GetXaxis()->FindBin(0.0 - 0.001), h2genanti->GetXaxis()->FindBin(0.0 + 0.001), -1, -1, "E");

    TH1D *h1gensel8 = h2gen->ProjectionY("h1gensel8", h2gen->GetXaxis()->FindBin(1.0 - 0.001), h2gen->GetXaxis()->FindBin(1.0 + 0.001), -1, -1, "E");
    TH1D *h1genantisel8 = h2genanti->ProjectionY("h1genantisel8", h2genanti->GetXaxis()->FindBin(1.0 - 0.001), h2genanti->GetXaxis()->FindBin(1.0 + 0.001), -1, -1, "E");

    TH1D *h1genallev = h2gen->ProjectionY("h1genallev", h2gen->GetXaxis()->FindBin(2.0 - 0.001), h2gen->GetXaxis()->FindBin(2.0 + 0.001), -1, -1, "E");
    TH1D *h1genantiallev = h2genanti->ProjectionY("h1genantiallev", h2genanti->GetXaxis()->FindBin(2.0 - 0.001), h2genanti->GetXaxis()->FindBin(2.0 + 0.001), -1, -1, "E");

    TH1D *h1INELgt0 = h2gen->ProjectionY("h1INELgt0", h2gen->GetXaxis()->FindBin(0.0 - 0.001), h2gen->GetXaxis()->FindBin(0.0 + 0.001), h2gen->GetZaxis()->FindBin(0.0 + 1e-7), h2gen->GetZaxis()->FindBin(100 - 1e-7), "E");
    TH1D *h1INELgt0anti = h2genanti->ProjectionY("h1INELgt0anti", h2genanti->GetXaxis()->FindBin(0.0 - 0.001), h2genanti->GetXaxis()->FindBin(0.0 + 0.001), h2genanti->GetZaxis()->FindBin(0.0 + 1e-7), h2genanti->GetZaxis()->FindBin(100 - 1e-7), "E");

    TH1D *h1rec = h2rec->ProjectionX("h1rec");
    TH1D *h1recanti = h2recanti->ProjectionX("h1recanti");

    if (h1recanti == nullptr || h1rec == nullptr || h1gen == nullptr || h1genanti == nullptr || h1gensel8 == nullptr || h1genantisel8 == nullptr || h1genallev == nullptr || h1genantiallev == nullptr || h1INELgt0 == nullptr || h1INELgt0anti == nullptr)
    {
        cout << "Error reading efficiency histograms" << endl;
        return;
    }
    h1gen->Add(h1genanti, 1);
    h1rec->Add(h1recanti, 1);
    h1gensel8->Add(h1genantisel8, 1);
    h1genallev->Add(h1genantiallev, 1);
    h1INELgt0->Add(h1INELgt0anti, 1); // INEL > 0 with vtz cut

    TH1F *hyieldraw = (TH1F *)fileraw->Get(("lf-k892analysis" + sysvarname + "/K892/0/hRawYields").c_str());
    if (hyieldraw == nullptr)
    {
        cout << "Error reading yield histogram" << endl;
        return;
    }
    TH1F *heffvtz = (TH1F *)hyieldraw->Clone();
    TH1F *heffVtzSel8 = (TH1F *)hyieldraw->Clone();
    TH1F *heffallevsel = (TH1F *)hyieldraw->Clone();
    TH1F *heffINELgt0 = (TH1F *)hyieldraw->Clone();
    cout << "bins: " << heffvtz->GetNbinsX() << endl;
    TH1F *hratio1 = (TH1F *)h1gen->Clone();
    TH1F *hratio2 = (TH1F *)h1gen->Clone();
    for (int i = 0; i < heffvtz->GetNbinsX(); i++)
    {
        lowpt = pT_bins[i];
        highpt = pT_bins[i + 1];
        cout << "lowpt: " << lowpt << " highpt: " << highpt << endl;
        double efficiency = h1rec->Integral(h1rec->GetXaxis()->FindBin(lowpt), h1rec->GetXaxis()->FindBin(highpt)) / h1gen->Integral(h1gen->GetXaxis()->FindBin(lowpt), h1gen->GetXaxis()->FindBin(highpt));
        double efficiencysel8 = h1rec->Integral(h1rec->GetXaxis()->FindBin(lowpt), h1rec->GetXaxis()->FindBin(highpt)) / h1gensel8->Integral(h1gensel8->GetXaxis()->FindBin(lowpt), h1gensel8->GetXaxis()->FindBin(highpt));
        double efficiencyall = h1rec->Integral(h1rec->GetXaxis()->FindBin(lowpt), h1rec->GetXaxis()->FindBin(highpt)) / h1genallev->Integral(h1genallev->GetXaxis()->FindBin(lowpt), h1genallev->GetXaxis()->FindBin(highpt));
        double effciencyINELgt0 = h1rec->Integral(h1rec->GetXaxis()->FindBin(lowpt), h1rec->GetXaxis()->FindBin(highpt)) / h1INELgt0->Integral(h1INELgt0->GetXaxis()->FindBin(lowpt), h1INELgt0->GetXaxis()->FindBin(highpt));
        cout << "Efficiency: " << efficiency << "Efficiencysel8" << efficiencysel8 << " Efficiencyall " << efficiencyall << endl;

        heffvtz->SetBinContent(i + 1, efficiency);
        heffVtzSel8->SetBinContent(i + 1, efficiencysel8);
        heffallevsel->SetBinContent(i + 1, efficiencyall);
        heffINELgt0->SetBinContent(i + 1, effciencyINELgt0);
        hyieldraw->SetBinContent(i + 1, hyieldraw->GetBinContent(i + 1) / efficiency);

        // ratios
        double ratio1 = efficiency / efficiencysel8;
        double ratio2 = efficiency / efficiencyall;
        hratio1->SetBinContent(i + 1, ratio1);
        hratio2->SetBinContent(i + 1, ratio2);
    }
    // TCanvas *c1 = new TCanvas();
    // heffvtz->Draw();

    TCanvas *c2 = new TCanvas();
    heffVtzSel8->Draw();

    TCanvas *c3 = new TCanvas();
    heffallevsel->Draw();

    // TCanvas *c4 = new TCanvas();
    // gPad->SetLogy();
    // hyieldraw->Draw();

    // TCanvas *c6 = new TCanvas();
    // heffINELgt0->Draw();

    TFile *fileeff10 = new TFile("efficiencies.root", "RECREATE");
    heffvtz->Write("INEL10");
    heffINELgt0->Write("INEL10gt0");
    heffVtzSel8->Write("INEL10_sel8");
    heffallevsel->Write("INEL10_alleventsel");

    // // Ratios

    // TCanvas *c5 = new TCanvas("c5", "c5", 720, 720);
    // SetCanvasStyle(c5, 0.25, 0.03, 0.03, 0.15);
    // double pad1Size, pad2Size;
    // canvas_style(c5, pad1Size, pad2Size);
    // c5->cd(1);
    // SetHistoStyle(h1, 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
    // h1->GetYaxis()->SetTitleSize(0.04 / pad1Size);
    // h1->SetMaximum(h1->GetMaximum() * 2);
    // h1->GetYaxis()->SetTitleOffset(1.15);
    // h1->GetXaxis()->SetTitleOffset(1.02);
    // h1->SetMarkerStyle(22);
    // h1->SetMarkerSize(1);
    // h1->SetMinimum(7e-6);
    // h1->Draw("pe");
    // c5->SetLogy(1);
    // gPad->SetLogy(1);
    // gr->SetMarkerStyle(29);
    // gr->SetMarkerSize(1);
    // gr->SetMarkerColor(kRed);
    // gr->SetLineColor(kRed);
    // gr->Draw("psame");
    // // h900->SetMarkerStyle(23);
    // // h900->SetMarkerSize(1);
    // // h900->SetMarkerColor(kGreen);
    // // h900->SetLineColor(kGreen);
    // // h900->Draw("pesame");

    // TLegend *leg = new TLegend(0.5, 0.5, 0.8, 0.8);
    // SetLegendStyle(leg);
    // leg->AddEntry(h1, "pp 13.6 TeV (This Analysis)", "lpe");
    // // leg->AddEntry(fLevyTsallis, "Levy-Tsallis Fit (pp 13.6 TeV)", "l");
    // leg->AddEntry(gr, "pp 13 TeV (Published)", "lpe");
    // // leg->AddEntry(h900, "pp 900 GeV", "lpe");
    // leg->SetTextSize(0.05);
    // leg->Draw();

    // c5->cd(2);
    // TH1F *hdummy = (TH1F *)h1->Clone();
    // for (int i = 0; i < hdummy->GetNbinsX(); i++)
    // {
    //     hdummy->SetBinContent(i + 1, 0);
    //     hdummy->SetBinError(i + 1, 0);
    // }

    // SetgrgrStyle(gratio, 1, 53, 1, 0.05, 0.05, 0.04 / pad2Size, 0.04 / pad2Size, 1.13, 1.8);
    // gratio->GetYaxis()->SetTitleSize(0.04 / pad2Size);
    // gratio->GetXaxis()->SetTitleSize(0.04 / pad2Size);
    // gratio->SetMarkerStyle(23);
    // gratio->SetMarkerSize(1.0);
    // gratio->SetMarkerColor(kRed);
    // gratio->SetLineColor(kRed);
    // gratio->GetYaxis()->SetTitle("#frac{This Analysis}{Published}");
    // gratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // gratio->GetXaxis()->CenterTitle(1);
    // // gratio->GetYaxis()->CenterTitle(1);
    // gratio->GetYaxis()->SetTitleOffset(0.45);
    // gratio->GetYaxis()->SetNdivisions(505);
    // gratio->GetXaxis()->SetRangeUser(0, 15);
    // // hdummy->Draw();
    // gratio->Draw("ap");
}