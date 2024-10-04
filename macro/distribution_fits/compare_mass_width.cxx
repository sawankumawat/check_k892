#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TBox.h>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"

void compare_mass_width()
{
    const string kResBkg = "ROTATED";
    // path to the output folder
    const string outputfolder_str = "../" + kSignalOutput + "/" + kchannel + "/" + kfoldername;
    const string fits_folder_str = outputfolder_str + "/fits/";

    TFile *fexpol = new TFile((fits_folder_str + "fitparams_" + kResBkg + "_expol" + ".root").c_str(), "read");
    TFile *fBoltzman = new TFile((fits_folder_str + "fitparams_" + kResBkg + "_Boltzman" + ".root").c_str(), "read");
    TFile *fboth[2] = {fexpol, fBoltzman};
    if (fexpol->IsZombie() || fBoltzman->IsZombie())
    {
        cout << "Error opening files" << endl;
        return;
    }
    TGraphErrors *r_mass[2];
    TGraphErrors *r_width[2];
    double masses
    for (int i = 0; i < 2; i++)
    {
        r_mass[i] = (TGraphErrors *)fboth[i]->Get(Form("mass_all", i));
        r_width[i] = (TGraphErrors *)fboth[i]->Get(Form("width_all", i));
        if(r_mass[i] == nullptr || r_width[i] == nullptr)
        {
            cout << "Error reading graphs" << endl;
            return;
        }
    }

    TCanvas *c = new TCanvas("c", "Mass Plot", 800, 600);
    SetCanvasStyle(c, 0.25, 0.03, 0.05, 0.14);

    const int n = 5;

    double mass[n] = {1.05, 1.15, 1.25, 1.35, 1.45};    // Mass values
    double massErr[n] = {0.02, 0.02, 0.02, 0.02, 0.02}; // Mass errors (x-errors)
    double y[n] = {1, 2, 3, 4, 5};                      // Y-values (dummy for plot)
    double yErr[n] = {0.1, 0.15, 0.12, 0.2, 0.14};      // Y-errors (dummy for plot)

    TGraphErrors *graph = new TGraphErrors(n, mass, y, massErr, yErr);
    graph->SetTitle("; Mass (GeV/c)");
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    graph->Draw("AP");

    // Add vertical line for PDG value (1.27 GeV/c)
    double pdg_mass = 1.27;
    TLine *line = new TLine(pdg_mass, 0, pdg_mass, 6); // Vertical line at 1.27 GeV
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->Draw("same");

    // Add a band around PDG value representing error
    double pdg_error = 0.03;                                                // PDG error
    TBox *box = new TBox(pdg_mass - pdg_error, 0, pdg_mass + pdg_error, 6); // Error band
    box->SetFillColorAlpha(kRed, 0.01);                                     // Set more transparency (0.1 means highly transparent)
    box->SetFillStyle(3003);
    box->Draw("same");

    // Set axis ranges
    graph->GetXaxis()->SetLimits(1.0, 2.0); // Set X-axis from 1 to 2 GeV/c
    graph->GetYaxis()->SetRangeUser(0, 6);  // Adjust Y-axis range
    graph->GetYaxis()->SetNdivisions(5);
    graph->GetYaxis()->SetBinLabel(1, "");

    TLatex latex;
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.7, 1, "Measurement A");
    latex.DrawLatex(0.7, 2, "Measurement B");
    latex.DrawLatex(0.7, 3, "Measurement C");
    latex.DrawLatex(0.7, 4, "Measurement D");
    latex.DrawLatex(0.7, 5, "Measurement E");
    c->Update();
}
