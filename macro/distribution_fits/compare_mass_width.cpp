#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TBox.h>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"

void plot(TGraphErrors *gr_mass_expol[3], TGraphErrors *gr_mass_boltzman[3], TGraphErrors *gr_HERA_mass[3], TLine *line[3], TBox *box[3], float limits[3][2], string resonance[3], string xaxistitle, double *pdg, double *pdgerr);
const string kResBkg = "ROTATED";
// path to the output folder
const string outputfolder_str = "../" + kSignalOutput + "/" + kchannel + "/" + kfoldername;
const string fits_folder_str = outputfolder_str + "/fits/";
TString save_folder = fits_folder_str + "compare_mass_width/";
const string save_folder_str = fits_folder_str + "compare_mass_width/";
void compare_mass_width()
{

    // Create the folder to save the fitted distributions
    if (gSystem->mkdir(save_folder, kTRUE))
    {
        std::cout << "Folder " << save_folder << " created successfully." << std::endl;
    }
    else
    {
        std::cout << "Creating folder " << save_folder << std::endl;
    }

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
    float mass_values[2][6];
    float mass_values_errors[2][6];
    float width_values[2][6];
    float width_values_errors[2][6];

    for (int i = 0; i < 2; i++)
    {
        r_mass[i] = (TGraphErrors *)fboth[i]->Get("mass_all");
        r_width[i] = (TGraphErrors *)fboth[i]->Get("width_all");
        if (r_mass[i] == nullptr || r_width[i] == nullptr)
        {
            cout << "Error reading graphs" << endl;
            return;
        }
        for (int j = 0; j < 6; j++)
        {
            mass_values[i][j] = r_mass[i]->GetY()[j];
            width_values[i][j] = r_width[i]->GetY()[j];
            mass_values_errors[i][j] = r_mass[i]->GetEY()[j];
            width_values_errors[i][j] = r_width[i]->GetEY()[j];
        }
    }

    // TGraphErrors *graph = new TGraphErrors(n, mass, y, massErr, yErr);
    TGraphErrors *gr_mass_expol[3];
    TGraphErrors *gr_width_expol[3];
    TGraphErrors *gr_mass_boltzman[3];
    TGraphErrors *gr_width_boltzman[3];
    TGraphErrors *gr_HERA_mass[3];
    TGraphErrors *gr_HERA_width[3];
    TLine *line[3];

    double pdgmasses[3] = {f1270Mass, f1525Mass, f1710Mass};
    double pdgmasserrors[3] = {f1270MassErr, f1525MassErr, f1710MassErr};
    double pdgwidths[3] = {f1270Width, f1525Width, f1710Width};
    double pdgwidtherrors[3] = {f1270WidthErr, f1525WidthErr, f1710WidthErr};
    float limits_mass[3][2] = {{1.26, 1.33}, {1.505, 1.53}, {1.66, 1.77}};
    float limits_width[3][2] = {{0.02, 0.38}, {0.05, 0.15}, {0.06, 0.4}};
    string resonances[3] = {"f1270", "f1525", "f1710"};

    // HERA reuslts without interference terms
    float HERA_mass[3] = {1.304, 1.523, 1.692};
    float HERA_mass_errors[3] = {0.006, 0.003, 0.006};
    float HERA_width[3] = {0.061, 0.071, 0.125};
    float HERA_width_errors[3] = {0.011, 0.005, 0.012};
    TBox *box[3];

    for (int i = 0; i < 3; i++)
    {
        gr_mass_expol[i] = new TGraphErrors(3);
        gr_width_expol[i] = new TGraphErrors(3);
        gr_mass_boltzman[i] = new TGraphErrors(3);
        gr_width_boltzman[i] = new TGraphErrors(3);
        gr_HERA_mass[i] = new TGraphErrors(3);
        gr_HERA_width[i] = new TGraphErrors(3);
        SetGrapherrorStyle(gr_mass_expol[i]);
        SetGrapherrorStyle(gr_mass_boltzman[i]);
        SetGrapherrorStyle(gr_width_expol[i]);
        SetGrapherrorStyle(gr_width_boltzman[i]);

        gr_mass_expol[i]->SetPoint(i, mass_values[0][i + 3], 1);
        gr_mass_expol[i]->SetPointError(i, mass_values_errors[0][i + 3], 0);
        gr_mass_boltzman[i]->SetPoint(i, mass_values[1][i + 3], 2);
        gr_mass_boltzman[i]->SetPointError(i, mass_values_errors[1][i + 3], 0);
        gr_width_expol[i]->SetPoint(i, width_values[0][i + 3], 1);
        gr_width_expol[i]->SetPointError(i, width_values_errors[0][i + 3], 0);
        gr_width_boltzman[i]->SetPoint(i, width_values[1][i + 3], 2);
        gr_width_boltzman[i]->SetPointError(i, width_values_errors[1][i + 3], 0);
        gr_HERA_mass[i]->SetPoint(i, HERA_mass[i], 3);
        gr_HERA_mass[i]->SetPointError(i, HERA_mass_errors[i], 0);
        gr_HERA_width[i]->SetPoint(i, HERA_width[i], 3);
        gr_HERA_width[i]->SetPointError(i, HERA_width_errors[i], 0);
    }

    plot(gr_mass_expol, gr_mass_boltzman, gr_HERA_mass, line, box, limits_mass, resonances, "Mass", pdgmasses, pdgmasserrors);
    plot(gr_width_expol, gr_width_boltzman, gr_HERA_width, line, box, limits_width, resonances, "Width", pdgwidths, pdgwidtherrors);

} // end of main function

void plot(TGraphErrors *gr_mass_expol[3], TGraphErrors *gr_mass_boltzman[3], TGraphErrors *gr_HERA_mass[3], TLine *line[3], TBox *box[3], float limits[3][2], string resonances[3], string xaxistitle, double *pdg, double *pdgerr)
{
    for (int i = 0; i < 3; i++)
    {

        TCanvas *c = new TCanvas("", "", 800, 600);
        SetCanvasStyle(c, 0.2, 0.03, 0.05, 0.14);
        gr_mass_expol[i]->SetTitle(Form("; %s (GeV/c)", xaxistitle.c_str()));
        gr_mass_expol[i]->SetMarkerStyle(20);
        gr_mass_expol[i]->SetMarkerColor(kBlue);
        gr_mass_expol[i]->SetLineColor(kBlue);
        gr_mass_boltzman[i]->SetMarkerStyle(21);
        gr_mass_boltzman[i]->SetMarkerColor(1);
        gr_mass_boltzman[i]->SetLineColor(1);
        gr_HERA_mass[i]->SetMarkerStyle(22);
        gr_HERA_mass[i]->SetMarkerColor(kGreen);
        gr_HERA_mass[i]->SetLineColor(kGreen);
        gr_mass_expol[i]->Draw("AP");
        gr_mass_boltzman[i]->Draw("P SAME");
        gr_HERA_mass[i]->Draw("P SAME");
        line[i] = new TLine(pdg[i], 0, pdg[i], 3.5);
        line[i]->SetLineColor(kRed);
        line[i]->SetLineWidth(2);
        line[i]->Draw("same");
        box[i] = new TBox(pdg[i] - pdgerr[i], 0, pdg[i] + pdgerr[i], 3.5);
        box[i]->SetFillColor(kYellow);
        box[i]->SetFillStyle(3001);
        box[i]->Draw("same");

        TLegend *leg = new TLegend(0.7, 0.2, 0.95, 0.55);
        leg->SetFillStyle(0);
        // leg->SetBorderSize(0);
        leg->AddEntry(gr_HERA_mass[i], "HERA", "lp");
        leg->AddEntry(gr_mass_boltzman[i], "Boltzman fit", "lp");
        leg->AddEntry(gr_mass_expol[i], "Expol fit", "lp");
        leg->AddEntry(box[i], Form("PDG %s", resonances[i].c_str()), "f");
        leg->Draw("same");
        TLegend *leg2 = new TLegend(0.7, 0.2, 0.95, 0.55); // this is created to show line and box in the legend
        leg2->SetFillStyle(0);
        leg2->SetBorderSize(0);
        leg2->AddEntry((TObject *)0, "", "");
        leg2->AddEntry((TObject *)0, "", "");
        leg2->AddEntry((TObject *)0, "", "");
        leg2->AddEntry(line[i], " ", "e");
        leg2->Draw("same");

        gr_mass_expol[i]->GetXaxis()->SetLimits(limits[i][0], limits[i][1]); // Set X-axis from 1 to 2 GeV/c
        gr_mass_expol[i]->GetYaxis()->SetRangeUser(0, 4);                    // Adjust Y-axis range
        gr_mass_expol[i]->GetYaxis()->SetNdivisions(3);
        gr_mass_expol[i]->GetYaxis()->SetBinLabel(1, "");

        // Get the minimum x-axis value
        double xMin = gr_mass_expol[i]->GetXaxis()->GetXmin();

        // Define a small offset to position the text just to the left of the y-axis
        double offset = 0.25 * (gr_mass_expol[i]->GetXaxis()->GetXmax() - xMin);

        TLatex latex;
        latex.SetTextSize(0.04);
        latex.DrawLatex(xMin - offset, 1, "Expol fit");
        latex.DrawLatex(xMin - offset, 2, "Boltzman fit");
        latex.DrawLatex(xMin - offset, 3, "HERA Expt.");

        c->SaveAs((save_folder_str + xaxistitle + resonances[i] + ".png").c_str());
    }
}
