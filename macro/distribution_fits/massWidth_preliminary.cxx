#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TBox.h>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"

void massWidth_preliminary()
{
    // integrated pt range
    //  float mass_alice = 1709.1, stat_mass = 4.5, sys_low_mass = 14.1, sys_high_mass = 14.1;
    //  float width_alice = 145.6, stat_width = 15.1, sys_low_width = 41.9, sys_high_width = 41.9;

    // pt 3 to 5 GeV/c
    float mass_alice = 1691.9, stat_mass = 3.7, sys_low_mass = 9.2, sys_high_mass = 9.2;
    float width_alice = 146.6, stat_width = 13.5, sys_low_width = 22.3, sys_high_width = 22.3;

    float mass_PDG = 1733, pdg_stat_mass = 8;
    float width_PDG = 150, pdg_stat_width = 12;

    TCanvas *c = new TCanvas("c", "Mass and Width Comparison", 720, 720);
    SetCanvasStyle(c, 0.15, 0.01, 0.01, 0.14);
    c->Divide(1, 2); // Divide canvas into 2 vertical pads

    // Top panel for Mass
    c->cd(1);
    TPad *pad1 = (TPad *)gPad;
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.18);
    gPad->SetRightMargin(0.01);
    gPad->SetLeftMargin(0.01);

    TGraphErrors *g_mass = new TGraphErrors(1);
    SetGrapherrorStyle(g_mass);
    g_mass->SetPoint(0, mass_alice, 1);
    g_mass->SetPointError(0, stat_mass, 0);
    g_mass->SetMarkerStyle(20);
    g_mass->SetMarkerColor(1);
    g_mass->SetLineColor(1);
    g_mass->SetMarkerSize(1.5);
    g_mass->SetLineWidth(2);
    g_mass->GetYaxis()->SetLabelSize(0);
    g_mass->GetYaxis()->SetTickLength(0);
    g_mass->GetXaxis()->SetTitleSize(0.06);
    g_mass->GetXaxis()->SetLabelSize(0.06);
    g_mass->Draw("AP");

    if (sys_low_mass != 0)
    {
        TBox *box = new TBox(mass_alice - sys_low_mass, 0.97, mass_alice + sys_high_mass, 1.03);
        box->SetFillStyle(0);
        box->SetLineColor(1);
        box->SetLineWidth(2);
        box->Draw("same");
    }

    g_mass->GetXaxis()->SetTitle("Mass (MeV/c^{2})");
    g_mass->GetXaxis()->SetLimits(1670, 1810);
    g_mass->SetMinimum(0.5);
    g_mass->SetMaximum(1.5);

    TBox *pdg_box = new TBox(mass_PDG - pdg_stat_mass, 0.5, mass_PDG + pdg_stat_mass, 1.49);
    pdg_box->SetFillColor(kYellow + 1);
    pdg_box->SetLineColor(kYellow + 1);
    pdg_box->SetFillStyle(3002);
    pdg_box->Draw("same");

    TLine *pdg_line = new TLine(mass_PDG, 0.5, mass_PDG, 1.49);
    pdg_line->SetLineColor(kRed);
    pdg_line->SetLineWidth(2);
    pdg_line->Draw("same");

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.06);
    lat.SetTextFont(42);
    lat.DrawLatex(0.05, 0.89, "ALICE");
    lat.DrawLatex(0.05, 0.79, "pp, #sqrt{s} = 13.6 TeV");
    // lat.DrawLatex(0.05, 0.69, "FT0M (0-100%), |y|<0.5");

    TLatex lat2;
    lat2.SetNDC();
    lat2.SetTextSize(0.05);
    lat2.SetTextFont(42);
    lat2.DrawLatex(0.6, 0.29, "Uncertainties: stat.(bar), sys.(box)");

    TLegend *leg = new TLegend(0.65, 0.65, 0.95, 0.95);
    leg->SetTextFont(42);
    leg->SetTextSize(0.06);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetHeader("Mass f_{0}(1710)");
    leg->AddEntry(g_mass, "pp 13.6 TeV", "lp");
    leg->AddEntry(pdg_box, "PDG 2024", "F");
    leg->Draw();

    TLegend *legtemp = new TLegend(0.65, 0.65, 0.95, 0.95);
    legtemp->SetTextFont(42);
    legtemp->SetTextSize(0.06);
    legtemp->SetBorderSize(0);
    legtemp->SetFillStyle(0);
    legtemp->AddEntry((TObject *)0, "", "");
    legtemp->AddEntry((TObject *)0, "", "");
    legtemp->AddEntry(pdg_line, " ", "e");
    legtemp->Draw();

    // Bottom panel for Width
    c->cd(2);
    TPad *pad2 = (TPad *)gPad;
    gPad->SetTopMargin(0.0);
    gPad->SetBottomMargin(0.18);
    gPad->SetRightMargin(0.01);
    gPad->SetLeftMargin(0.01);

    TGraphErrors *g_width = new TGraphErrors(1);
    SetGrapherrorStyle(g_width);
    g_width->SetPoint(0, width_alice, 1);
    g_width->SetPointError(0, stat_width, 0);
    g_width->SetMarkerStyle(21);
    g_width->SetMarkerColor(4);
    g_width->SetLineColor(4);
    g_width->SetMarkerSize(1.5);
    g_width->SetLineWidth(2);
    g_width->GetYaxis()->SetLabelSize(0);
    g_width->GetYaxis()->SetTickLength(0);
    g_width->GetXaxis()->SetTitleSize(0.06);
    g_width->GetXaxis()->SetLabelSize(0.06);
    g_width->Draw("AP");

    if (sys_low_width != 0)
    {
        TBox *box = new TBox(width_alice - sys_low_width, 0.97, width_alice + sys_high_width, 1.03);
        box->SetFillStyle(0);
        box->SetLineColor(4);
        box->SetLineWidth(2);
        box->Draw("same");
    }

    g_width->GetXaxis()->SetTitle("Width (MeV/c^{2})");
    g_width->GetXaxis()->SetLimits(40, 289);
    g_width->SetMinimum(0.5);
    g_width->SetMaximum(1.5);

    TBox *pdg_box_width = new TBox(width_PDG - pdg_stat_width, 0.5, width_PDG + pdg_stat_width, 1.49);
    pdg_box_width->SetFillColor(kMagenta - 10);
    pdg_box_width->SetLineColor(kMagenta - 10);
    pdg_box_width->SetFillStyle(3002);
    pdg_box_width->Draw("same");

    TLine *pdg_line_width = new TLine(width_PDG, 0.5, width_PDG, 1.49);
    pdg_line_width->SetLineColor(kRed + 3);
    pdg_line_width->SetLineWidth(2);
    pdg_line_width->Draw("same");

    lat2.DrawLatex(0.6, 0.29, "Uncertainties: stat.(bar), sys.(box)");

    TLegend *leg2 = new TLegend(0.65, 0.65, 0.95, 0.95);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.06);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetHeader("Width f_{0}(1710)");
    leg2->AddEntry(g_width, "pp 13.6 TeV", "lp");
    leg2->AddEntry(pdg_box_width, "PDG 2024", "F");
    leg2->Draw();

    TLegend *legtemp2 = new TLegend(0.65, 0.65, 0.95, 0.95);
    legtemp2->SetTextFont(42);
    legtemp2->SetTextSize(0.06);
    legtemp2->SetBorderSize(0);
    legtemp2->SetFillStyle(0);
    legtemp2->AddEntry((TObject *)0, "", "");
    legtemp2->AddEntry((TObject *)0, "", "");
    legtemp2->AddEntry(pdg_line_width, " ", "e");
    legtemp2->Draw();

    TLatex lat3;
    lat3.SetNDC();
    lat3.SetTextSize(0.06);
    lat3.SetTextFont(42);
    lat3.DrawLatex(0.05, 0.89, "FT0M (0-100%), |y|<0.5");
    // lat3.DrawLatex(0.05, 0.79, "0.0 < #it{p}_{T} < 30 GeV/c");
    lat3.DrawLatex(0.05, 0.79, "3.0 < #it{p}_{T} < 5 GeV/c");

    // c->SaveAs("/home/sawan/Music/f1710mass_width.png");
    c->SaveAs("/home/sawan/Music/f1710mass_width_highpt1.png");
}
