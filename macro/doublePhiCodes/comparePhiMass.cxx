#include <iostream>
#include "src/style.h"

TFile *OpenFile(const string &path);
template <typename T>
T *GetHisto(TFile *f, const std::string &name);

// No difference between opti9 and opti10 for shifting the kaon momentum

void comparePhiMass()
{
    TString savePath = "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/PhiInvMass";
    TFile *f1 = OpenFile((savePath + "/PhiParams26.root").Data());
    // TFile *f2 = OpenFile((savePath + "/PhiParams25.root").Data());
    // TFile *f2 = OpenFile((savePath + "/PhiParams25_Shifted.root").Data());
    TFile *f2 = OpenFile((savePath + "/PhiParams25_aiam.root").Data());

    TFile *f3 = OpenFile("PhiMomentumScaleVsPt/PhiMomentumScaleVsPt.root");
    TGraphErrors *gPhiMassSBcode2025 = GetHisto<TGraphErrors>(f3, "gPhiMass2025");
    TGraphErrors *gPhiMassSBcode2026 = GetHisto<TGraphErrors>(f3, "gPhiMass2026");

    TGraphErrors *gMassVsPt1 = GetHisto<TGraphErrors>(f1, "gMassVsPt");
    TGraphErrors *gMassVsPt2 = GetHisto<TGraphErrors>(f2, "gMassVsPt");

    // TFile *f4 = OpenFile((savePath + "/PhiParams25_aiamShifted.root").Data());
    TFile *f4 = OpenFile((savePath + "/PhiParams25_aiamShifted_opti10.root").Data());
    TGraphErrors *gMassVsPt3 = GetHisto<TGraphErrors>(f4, "gMassVsPt");

    TCanvas *cMassVsPt = new TCanvas("cMassVsPt", "Mass vs Pt", 720, 720);
    SetCanvasStyle(cMassVsPt, 0.20, 0.03, 0.05, 0.15);
    gMassVsPt1->SetMarkerStyle(20);
    gMassVsPt1->SetMarkerColor(kRed);
    gMassVsPt1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    gMassVsPt1->GetYaxis()->SetTitle("M_{#Phi} (GeV/#it{c}^{2})");
    gMassVsPt1->GetYaxis()->SetRangeUser(1.0168, 1.0218);
    SetGraphErrorStyle(gMassVsPt1);
    gMassVsPt1->GetYaxis()->SetTitleOffset(2.1);
    gMassVsPt1->SetLineColor(kRed);
    gMassVsPt1->Draw("APE");
    gMassVsPt2->SetMarkerStyle(21);
    gMassVsPt2->SetMarkerColor(kBlue);
    gMassVsPt2->SetLineColor(kBlue);
    gMassVsPt2->Draw("pe same");

    // SetGraphErrorStyle(gPhiMassSBcode2026);
    // gPhiMassSBcode2026->SetMarkerColor(kGreen + 2);
    // gPhiMassSBcode2026->SetLineColor(kGreen + 2);
    // gPhiMassSBcode2026->Draw("p same");
    // SetGraphErrorStyle(gPhiMassSBcode2025);
    // gPhiMassSBcode2025->SetMarkerColor(kMagenta);
    // gPhiMassSBcode2025->SetLineColor(kMagenta);
    // gPhiMassSBcode2025->Draw("p same");

    // SetGraphErrorStyle(gMassVsPt3);
    // gMassVsPt3->SetMarkerColor(kGreen + 2);
    // gMassVsPt3->SetLineColor(kGreen + 2);
    // gMassVsPt3->Draw("pe same");

    TLine *linePDG = new TLine(0.5, 1.019460, 30, 1.019460);
    linePDG->SetLineColor(kBlack);
    linePDG->SetLineStyle(7);
    linePDG->Draw();
    TBox *boxPDG = new TBox(0.5, 1.019460 - 0.000016, 30, 1.019460 + 0.000016);
    boxPDG->SetFillColor(kGray);
    boxPDG->SetFillStyle(3001);
    boxPDG->Draw();

    TLegend *legend = new TLegend(0.22, 0.75, 0.85, 0.85);
    // legend->SetNColumns(2);
    // legend->AddEntry((TObject *)0, "Sawan code", "");
    // legend->AddEntry((TObject *)0, "Sourav Bhaiya code", "");
    // legend->AddEntry(gMassVsPt1, "LHC26_skimmed", "p");
    // legend->AddEntry(gPhiMassSBcode2026, "LHC26_skimmed", "p");
    // legend->AddEntry(gMassVsPt2, "LHC25_skimmed", "p");
    // // legend->AddEntry(linePDG, "PDG Mass", "l");
    // legend->AddEntry(gPhiMassSBcode2025, "LHC25_skimmed", "p");

    legend->AddEntry(gMassVsPt1, "LHC26_skimmed", "p");
    legend->AddEntry(gMassVsPt2, "LHC25_skimmed", "p");
    // legend->AddEntry(gMassVsPt3, "LHC25_skimmed (Kaon momentum shifted)", "p");
    legend->AddEntry(linePDG, "PDG Mass", "l");

    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->SetBorderSize(0);
    legend->Draw();
    cMassVsPt->SaveAs("/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/PhiMassShift/PhiMassVsPt_Compare.png");

    TGraphErrors *gMassShift = new TGraphErrors(gMassVsPt1->GetN());
    for (int i = 0; i < gMassVsPt1->GetN(); i++)
    {
        double x1, y1, x2, y2;
        gMassVsPt1->GetPoint(i, x1, y1);
        gMassVsPt2->GetPoint(i, x2, y2);
        gMassShift->SetPoint(i, x1, y1 - y2);
        gMassShift->SetPointError(i, 0, sqrt(pow(gMassVsPt1->GetErrorY(i), 2) + pow(gMassVsPt2->GetErrorY(i), 2)));
        cout << "Bin " << i << ": Pt = " << x1
             << ", Mass (2026): " << y1 << ", Mass (2025): " << y2 << endl;
    }
    TCanvas *cMassShift = new TCanvas("cMassShift", "Mass Shift vs Pt", 720, 720);
    SetCanvasStyle(cMassShift, 0.20, 0.03, 0.05, 0.15);
    gMassShift->SetMarkerStyle(20);
    gMassShift->SetMarkerColor(kBlack);
    gMassShift->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    gMassShift->GetYaxis()->SetTitle("#Delta M_{#Phi} (GeV/#it{c}^{2})");
    gMassShift->GetYaxis()->SetRangeUser(-0.0007, 0.0017);
    SetGraphErrorStyle(gMassShift);
    gMassShift->GetYaxis()->SetTitleOffset(2.1);
    gMassShift->Draw("AP");
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(42);
    latex->SetTextSize(0.035);
    latex->DrawLatex(0.25, 0.85, "#Delta M_{#Phi} (LHC25_skimmed - LHC26_skimmed)");
    // cMassShift->SaveAs("PhiMassShiftVsPt.png");

    // TGraphErrors *gPurity25 = GetHisto<TGraphErrors>(f2, "gPurity");
    // TGraphErrors *gPurity26 = GetHisto<TGraphErrors>(f1, "gPurity");
    // TCanvas *cPurityVsPt = new TCanvas("cPurityVsPt", "Purity vs Pt", 720, 720);
    // SetCanvasStyle(cPurityVsPt, 0.15, 0.03, 0.05, 0.15);
    // SetGraphErrorStyle(gPurity25);
    // gPurity25->SetMarkerStyle(21);
    // gPurity25->SetMarkerColor(kBlue);
    // gPurity25->SetLineColor(kBlue);
    // gPurity26->SetMarkerStyle(20);
    // gPurity26->SetMarkerColor(kRed);
    // gPurity26->SetLineColor(kRed);
    // gPurity25->Draw("APE");
    // gPurity26->Draw("PE SAME");
    // latex->SetTextSize(0.035);
    // latex->DrawLatex(0.25, 0.88, "Fit function : Voigtian + Pol2");
    // latex->DrawLatex(0.25, 0.81, "Purity window: #it{M}_{#Phi} #pm 0.005 GeV/#it{c}^{2}");

    // TLegend *legendPurity = new TLegend(0.22, 0.75, 0.85, 0.85);
    // legendPurity->AddEntry(gPurity26, "LHC26_skimmed", "p");
    // legendPurity->AddEntry(gPurity25, "LHC25_skimmed", "p");
    // legendPurity->SetTextFont(42);
    // legendPurity->SetTextSize(0.03);
    // legendPurity->SetBorderSize(0);
    // legendPurity->Draw();
    // cPurityVsPt->SaveAs("/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/PhiMassShift/PhiPurityVsPt_Compare.png");
}

//==============End of the main code==================

TFile *OpenFile(const string &path)
{
    TFile *f = new TFile(path.c_str(), "read");
    if (f->IsZombie())
    {
        cout << "Error: File not found: " << path << endl;
        return nullptr;
    }
    return f;
}

template <typename T>
T *GetHisto(TFile *f, const std::string &name)
{
    T *histo = dynamic_cast<T *>(f->Get(name.c_str()));

    if (!histo)
    {
        std::cout << "Error: histo " << name
                  << " not found in file " << f->GetName() << std::endl;
        return nullptr;
    }

    return histo;
}
