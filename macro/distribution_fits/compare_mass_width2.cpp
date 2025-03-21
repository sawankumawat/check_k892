#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TBox.h>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"

string path = "/home/sawan/check_k892/output/glueball/LHC22o_pass7_small/358932/KsKs_Channel/higher-mass-resonances_id24937/fits/";

void compare_mass_width2()
{

    vector<float> mass_alice = {1709.1, 4.5, 14.1, 14.1};  // mass, stat. err, sys. err low, sys. err high
    vector<float> width_alice = {145.6, 15.1, 41.9, 41.9}; // width, stat. err, sys. err low, sys. err high

    vector<float> mass_BESIII = {1765, 2, 1, 1}; // e+e-
    vector<float> width_BESIII = {146, 3, 7, 1};

    vector<float> mass_HERA = {1692, 6, 3, 9}; // ep
    vector<float> width_HERA = {125, 12, 32, 19};

    vector<float> mass_WA02 = {1730, 15, 0, 0}; // pp -> KK and KsKs
    vector<float> width_WA02 = {100, 25, 0, 0};

    vector<float> mass_L3 = {1767, 14, 0, 0}; // e+e-
    vector<float> width_L3 = {187, 60, 0, 0};

    vector<float> mass_PDG = {1733, 8, 0, 0}; // PDG
    vector<float> width_PDG = {150, 12, 0, 0};

    TCanvas *c1 = new TCanvas("c1", "Mass Comparison", 720, 720);
    SetCanvasStyle(c1, 0.15, 0.01, 0.01, 0.14);
    vector<const char *> experiments = {"WA02", "L3", "HERA", "BESIII", "ALICE"};

    // Y positions for different experiments
    vector<float> y_positions = {1, 2, 3, 4, 5};

    const int N = 5; // PDG is not included as a data point
    float x[N], y_mass[N], ex_stat_mass[N], ex_sys_low_mass[N], ex_sys_high_mass[N];
    float y_width[N], ex_stat_width[N], ex_sys_low_width[N], ex_sys_high_width[N];

    vector<vector<float>> masses = {mass_WA02, mass_L3, mass_HERA, mass_BESIII, mass_alice};
    vector<vector<float>> widths = {width_WA02, width_L3, width_HERA, width_BESIII, width_alice};

    // Define different marker styles for each experiment
    int markers[] = {21, 20, 22, 23, 33}; // Different ROOT marker styles
    int colors[] = {1, 2, 4, kGreen + 3, 28};
    float markersizes[] = {1.5, 1.5, 1.5, 1.5, 1.8};
    TGraphErrors *g_stat_mass[N];
    TGraphErrors *g_stat_width[N];

    for (int i = 0; i < N; i++)
    {
        x[i] = masses[i][0];                // Mass value
        y_mass[i] = y_positions[i];         // Y-position
        ex_stat_mass[i] = masses[i][1];     // Statistical error
        ex_sys_low_mass[i] = masses[i][2];  // Systematic error (low)
        ex_sys_high_mass[i] = masses[i][3]; // Systematic error (high)

        g_stat_mass[i] = new TGraphErrors();
        SetGrapherrorStyle(g_stat_mass[i]);
        g_stat_mass[i]->SetPoint(i, x[i], y_mass[i]);
        g_stat_mass[i]->SetPointError(i, ex_stat_mass[i], 0);
        g_stat_mass[i]->SetMarkerStyle(markers[N - i - 1]); // Assign unique marker style
        g_stat_mass[i]->SetMarkerColor(colors[N - 1 - i]);
        g_stat_mass[i]->SetLineColor(colors[N - 1 - i]);
        g_stat_mass[i]->SetMarkerSize(markersizes[N - 1 - i]);
        g_stat_mass[i]->SetLineWidth(2);

        g_stat_mass[i]->SetTitle(0);
        if (i == 0)
            g_stat_mass[i]->Draw("AP");
        else
            g_stat_mass[i]->Draw("P same");
    }

    // // Systematic errors (asymmetric)
    // TGraphAsymmErrors *g_sys = new TGraphAsymmErrors(N, x, y, ex_sys_low_mass, ex_sys_high_mass, 0, 0);
    // g_sys->SetMarkerStyle(21);
    // g_sys->SetMarkerSize(1.2);
    // g_sys->SetMarkerColor(kBlack);
    // g_sys->SetFillColorAlpha(kRed, 0.3); // Red shaded error bars
    // g_sys->SetLineColor(kRed);
    // g_sys->Draw("p same");

    // Drawing systematic error bands
    for (int i = 0; i < N; i++)
    {
        if (ex_sys_low_mass[i] == 0)
            continue;
        TBox *box = new TBox(x[i] - ex_sys_low_mass[i], y_mass[i] - 0.08, x[i] + ex_sys_high_mass[i], y_mass[i] + 0.08);
        box->SetFillStyle(0);
        box->SetFillColor(0);
        box->SetLineColor(colors[N - 1 - i]);
        box->SetLineWidth(2);
        box->Draw("same");
    }

    // Formatting the axis
    g_stat_mass[0]->GetXaxis()->SetTitle("Mass (MeV/c^{2})");
    g_stat_mass[0]->GetYaxis()->SetTitle("");
    g_stat_mass[0]->GetYaxis()->SetNdivisions(6);
    g_stat_mass[0]->GetXaxis()->SetLimits(1655, 1880);
    g_stat_mass[0]->SetMinimum(0.5);
    g_stat_mass[0]->SetMaximum(5.5);
    g_stat_mass[0]->GetXaxis()->SetNdivisions(508);
    g_stat_mass[0]->GetYaxis()->SetLabelSize(0.0001);

    // Draw PDG vertical band
    float pdg_x = mass_PDG[0];
    float pdg_stat_mass = mass_PDG[1];

    TBox *pdg_box = new TBox(pdg_x - pdg_stat_mass, 0.5, pdg_x + pdg_stat_mass, 5.49);
    pdg_box->SetFillColor(kYellow + 1);
    pdg_box->SetLineColor(kYellow + 1);
    pdg_box->SetFillStyle(3002);
    pdg_box->Draw("same");

    TLine *pdg_line = new TLine(pdg_x, 0.5, pdg_x, 5.49);
    pdg_line->SetLineColor(kRed);
    pdg_line->SetLineWidth(2);
    pdg_line->Draw("same");

    // Get the minimum x-axis value
    double xMin = g_stat_mass[0]->GetXaxis()->GetXmin();

    // Define a small offset to position the text just to the left of the y-axis
    double offset = 0.15 * (g_stat_mass[0]->GetXaxis()->GetXmax() - xMin);

    // Experiment labels
    TLatex latex;
    latex.SetTextSize(0.04);
    // latex.SetNDC();
    latex.SetTextFont(42);
    for (int i = 0; i < N; i++)
    {
        latex.DrawLatex(xMin - offset, y_positions[i], experiments[i]);
    }

    // Legend
    TLegend *leg = new TLegend(0.63, 0.80, 0.95, 0.87);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(pdg_box, "PDG 2024", "F");
    leg->Draw();

    TLegend *leg2 = new TLegend(0.63, 0.80, 0.95, 0.87);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.03);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry(pdg_line, " ", "e");
    leg2->Draw();

    TLatex lat3;
    lat3.SetNDC();
    lat3.SetTextSize(0.023);
    lat3.SetTextFont(42);
    // lat3.DrawLatex(0.15, 0.89, "ALICE");
    // lat3.DrawLatex(0.15, 0.85, "pp, #sqrt{s} = 13.6 TeV");
    // lat3.DrawLatex(0.15, 0.81, "FT0M (0-100%), |y|<0.5");
    lat3.DrawLatex(0.6, 0.88, "Uncertainties: stat.(bars), sys.(boxes)");

    c1->SaveAs((path + "mass_comparison.pdf").c_str());
    c1->SaveAs((path + "mass_comparison.png").c_str());

    TCanvas *c2 = new TCanvas("c2", "Width Comparison", 720, 720);
    SetCanvasStyle(c2, 0.15, 0.01, 0.01, 0.14);

    for (int i = 0; i < N; i++)
    {
        x[i] = widths[i][0];                 // Width value
        y_width[i] = y_positions[i];         // Y-position
        ex_stat_width[i] = widths[i][1];     // Statistical error
        ex_sys_low_width[i] = widths[i][2];  // Systematic error (low)
        ex_sys_high_width[i] = widths[i][3]; // Systematic error (high)

        g_stat_width[i] = new TGraphErrors();
        SetGrapherrorStyle(g_stat_width[i]);
        g_stat_width[i]->SetPoint(i, x[i], y_width[i]);
        g_stat_width[i]->SetPointError(i, ex_stat_width[i], 0);
        g_stat_width[i]->SetMarkerStyle(markers[N - i - 1]); // Assign unique marker style
        g_stat_width[i]->SetMarkerColor(colors[N - 1 - i]);
        g_stat_width[i]->SetLineColor(colors[N - 1 - i]);
        g_stat_width[i]->SetMarkerSize(markersizes[N - 1 - i]);
        g_stat_width[i]->SetLineWidth(2);

        g_stat_width[i]->SetTitle(0);
        if (i == 0)
            g_stat_width[i]->Draw("AP");
        else
            g_stat_width[i]->Draw("P same");
    }

    // // Systematic errors (asymmetric)
    // TGraphAsymmErrors *g_sys = new TGraphAsymmErrors(N, x, y, ex_sys_low_width, ex_sys_high_width, 0, 0);
    // g_sys->SetMarkerStyle(21);
    // g_sys->SetMarkerSize(1.2);
    // g_sys->SetMarkerColor(kBlack);
    // g_sys->SetFillColorAlpha(kRed, 0.3); // Red shaded error bars
    // g_sys->SetLineColor(kRed);
    // g_sys->Draw("p same");

    // Drawing systematic error bands
    for (int i = 0; i < N; i++)
    {
        if (ex_sys_low_width[i] == 0)
            continue;

        TBox *box = new TBox(x[i] - ex_sys_low_width[i], y_width[i] - 0.08, x[i] + ex_sys_high_width[i], y_width[i] + 0.08);
        box->SetFillStyle(0);
        box->SetFillColor(0);
        box->SetLineColor(colors[N - 1 - i]);
        box->SetLineWidth(2);
        box->Draw("same");
    }

    // Formatting the axis
    g_stat_width[0]->GetXaxis()->SetTitle("Width (MeV/c^{2})");
    g_stat_width[0]->GetYaxis()->SetTitle("");
    g_stat_width[0]->GetYaxis()->SetNdivisions(6);
    g_stat_width[0]->GetXaxis()->SetLimits(50, 349);
    g_stat_width[0]->SetMinimum(0.5);
    g_stat_width[0]->SetMaximum(5.5);
    g_stat_width[0]->GetXaxis()->SetNdivisions(508);
    g_stat_width[0]->GetYaxis()->SetLabelSize(0.0001);

    // Draw PDG vertical band
    float pdg_width = width_PDG[0];
    float pdg_stat_width = width_PDG[1];

    TBox *pdg_box_width = new TBox(pdg_width - pdg_stat_width, 0.5, pdg_width + pdg_stat_width, 5.49);
    pdg_box_width->SetFillColor(kYellow + 1);
    pdg_box_width->SetLineColor(kYellow + 1);
    pdg_box_width->SetFillStyle(3002);
    pdg_box_width->Draw("same");

    TLine *pdg_line_width = new TLine(pdg_width, 0.5, pdg_width, 5.49);
    pdg_line_width->SetLineColor(kRed);
    pdg_line_width->SetLineWidth(2);
    pdg_line_width->Draw("same");

    // Get the minimum x-axis value
    double xMin_width = g_stat_width[0]->GetXaxis()->GetXmin();

    // Define a small offset to position the text just to the left of the y-axis
    double offset_width = 0.15 * (g_stat_width[0]->GetXaxis()->GetXmax() - xMin_width);

    // Experiment labels
    TLatex latex_width;
    latex_width.SetTextSize(0.04);
    // latex_width.SetNDC();
    latex_width.SetTextFont(42);
    for (int i = 0; i < N; i++)
    {
        latex_width.DrawLatex(xMin_width - offset_width, y_positions[i], experiments[i]);
    }

    // Legend
    TLegend *leg_width = new TLegend(0.63, 0.80, 0.95, 0.87);
    leg_width->SetTextFont(42);
    leg_width->SetTextSize(0.03);
    leg_width->SetBorderSize(0);
    leg_width->SetFillStyle(0);
    leg_width->AddEntry(pdg_box_width, "PDG 2024", "F");
    leg_width->Draw();

    TLegend *leg2_width = new TLegend(0.63, 0.80, 0.95, 0.87);
    leg2_width->SetTextFont(42);
    leg2_width->SetTextSize(0.03);
    leg2_width->SetBorderSize(0);
    leg2_width->SetFillStyle(0);
    leg2_width->AddEntry(pdg_line_width, " ", "e");
    leg2_width->Draw();

    TLatex lat3_width;
    lat3_width.SetNDC();
    lat3_width.SetTextSize(0.023);
    lat3_width.SetTextFont(42);
    // lat3_width.DrawLatex(0.15, 0.89, "ALICE");
    // lat3_width.DrawLatex(0.15, 0.85, "pp, #sqrt{s} = 13.6 TeV");
    // lat3_width.DrawLatex(0.15, 0.81, "FT0M (0-100%), |y|<0.5");
    lat3_width.DrawLatex(0.6, 0.88, "Uncertainties: stat.(bars), sys.(boxes)");

    c2->SaveAs((path + "width_comparison.pdf").c_str());
    c2->SaveAs((path + "width_comparison.png").c_str());
}