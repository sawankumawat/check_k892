#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TBox.h>
#include "../src/style.h"
#include "../src/fitfunc.h"
#include "../src/common_glue.h"
#include "../src/fitting_range_glue.h"

void plot(TGraphErrors *gr_mass_expol[4], TGraphErrors *gr_mass_boltzman[4], TGraphErrors *gr_HERA_mass[4], TGraphErrors *gr_mass_coherent[4], TGraphErrors *gr_mass_L3[4], TGraphErrors *gr_mass_BES[4], TGraphErrors *gr_mass_WA02[4], TGraphErrors *gr_mass_BESIII[4], TLine *line[4], TBox *box[4], float limits[4][2], string resonance[4], string xaxistitle, double *pdg, double *pdgerr);
const string kResBkg = "ROTATED";
// path to the output folder
const string outputfolder_str = "../" + kSignalOutput + "/" + kchannel + "/" + kfoldername;
const string fits_folder_str = outputfolder_str + "/fits/";
TString save_folder = fits_folder_str + "compare_mass_width/";
const string save_folder_str = fits_folder_str + "compare_mass_width/";

// total resonances to be included in the plots
const int tot_reso = 4;
// *************************************

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

    // TFile *fexpol = new TFile((fits_folder_str + "fitparams_" + kResBkg + "_expol" + ".root").c_str(), "read");
    // TFile *fBoltzman = new TFile((fits_folder_str + "fitparams_" + kResBkg + "_Boltzman" + ".root").c_str(), "read");
    // TFile *fboth[2] = {fexpol, fBoltzman};
    // if (fexpol->IsZombie() || fBoltzman->IsZombie())
    // {
    //     cout << "Error opening files" << endl;
    //     return;
    // }

    // TGraphErrors *r_mass[2];
    // TGraphErrors *r_width[2];
    // float mass_values[2][6];
    // float mass_values_errors[2][6];
    // float width_values[2][6];
    // float width_values_errors[2][6];

    // // 3rBW + expol fit
    // float mass_expol[tot_reso] = {1.273, 1.273, 1.515, 1.694};
    // float width_expol[tot_reso] = {0.190, 0.190, 0.07795, 0.2054};
    // float mass_expol_err[tot_reso] = {0.002, 0.002, 0.001, 0.005};
    // float width_expol_err[tot_reso] = {0.010, 0.010, 0.00189, 0.0129};

    // // 3rBW + Boltzman fit
    // float mass_boltzman[tot_reso] = {1.297, 1.297, 1.514, 1.693};
    // float width_boltzman[tot_reso] = {0.163, 0.163, 0.08869, 0.1579};
    // float mass_boltzman_err[tot_reso] = {0.001, 0.001, 0.001, 0.004};
    // float width_boltzman_err[tot_reso] = {0.0001, 0.0001, 0.0019, 0.008};

    // 3rBW + expol1 fit
    float mass_boltzman[tot_reso] = {1.2799, 1.2799, 1.5136, 1.6917};
    float width_boltzman[tot_reso] = {0.1643, 0.1643, 0.0803, 0.2131};
    float mass_boltzman_err[tot_reso] = {0.0021, 0.0021, 0.0005, 0.0047};
    float width_boltzman_err[tot_reso] = {0.0085, 0.0085, 0.0019, 0.0150};

    // 4rBW + expol1 fit
    float mass_expol[tot_reso] = {1.2184, 1.3048, 1.5130, 1.6936};
    float width_expol[tot_reso] = {0.1412, 0.1043, 0.0852, 0.2027};
    float mass_expol_err[tot_reso] = {0.0048, 0.0020, 0.0006, 0.0046};
    float width_expol_err[tot_reso] = {0.0171, 0.0079, 0.0021, 0.0178};

    // Coherent BW sum fit
    //  float mass_coherent[tot_reso] = {1.19, 1.33, 1.513, 1.705};
    //  float width_coherent[tot_reso] = {0.1557, 0.09088, 0.07935, 0.1705};
    //  float mass_coherent_err[tot_reso] = {0.01, 0.004, 0.001, 0.005};
    //  float width_coherent_err[tot_reso] = {0.0229, 0.00756, 0.00182, 0.0139};

    float mass_coherent[tot_reso] = {1.1900, 1.1900, 1.513, 1.705};
    float width_coherent[tot_reso] = {0.1557, 0.1557, 0.07935, 0.1705};
    float mass_coherent_err[tot_reso] = {0.0100, 0.0100, 0.001, 0.005};
    float width_coherent_err[tot_reso] = {0.0229, 0.0229, 0.00182, 0.0139};

    // HERA reuslts without interference terms
    float HERA_mass[tot_reso] = {1.304, 1.304, 1.523, 1.692};
    float HERA_mass_errors[tot_reso] = {0.006, 0.006, 0.003, 0.006};
    float HERA_width[tot_reso] = {0.061, 0.061, 0.071, 0.125};
    float HERA_width_errors[tot_reso] = {0.011, 0.011, 0.005, 0.012};

    // L3 results without interference terms
    float L3_mass[tot_reso] = {1.239, 1.239, 1.523, 1.767};
    float L3_mass_errors[tot_reso] = {0.006, 0.006, 0.006, 0.014};
    float L3_width[tot_reso] = {0.078, 0.078, 0.1, 0.187};
    float L3_width_errors[tot_reso] = {0.019, 0.019, 0.015, 0.06};

    // BES results without interference terms
    float BES_mass[tot_reso] = {1.262, 1.262, 0, 1.765};
    float BES_mass_errors[tot_reso] = {0.002, 0.002, 0, 0.013};
    float BES_width[tot_reso] = {0.175, 0.175, 0, 0.145};
    float BES_width_errors[tot_reso] = {0.006, 0.006, 0, 0.069};

    // WA02 experiment
    float WA02_mass[tot_reso] = {1.305, 1.305, 1.515, 1.730};
    float WA02_mass_errors[tot_reso] = {0.02, 0.02, 0.015, 0.015};
    float WA02_width[tot_reso] = {0.132, 0.132, 0.070, 0.1};
    float WA02_width_errors[tot_reso] = {0.025, 0.025, 0.025, 0.025};

    // BESIII experiment
    float BESIII_mass[tot_reso] = {0, 0, 0, 1.765};
    float BESIII_mass_errors[tot_reso] = {0, 0, 0, 0.002};
    float BESIII_width[tot_reso] = {0, 0, 0, 0.146};
    float BESIII_width_errors[tot_reso] = {0, 0, 0, 0.003};

    // for (int i = 0; i < 2; i++)
    // {
    //     r_mass[i] = (TGraphErrors *)fboth[i]->Get("mass_all");
    //     r_width[i] = (TGraphErrors *)fboth[i]->Get("width_all");
    //     if (r_mass[i] == nullptr || r_width[i] == nullptr)
    //     {
    //         cout << "Error reading graphs" << endl;
    //         return;
    //     }
    //     for (int j = 0; j < 6; j++)
    //     {
    //         mass_values[i][j] = r_mass[i]->GetY()[j];
    //         width_values[i][j] = r_width[i]->GetY()[j];
    //         mass_values_errors[i][j] = r_mass[i]->GetEY()[j];
    //         width_values_errors[i][j] = r_width[i]->GetEY()[j];
    //     }
    // }

    // TGraphErrors *graph = new TGraphErrors(n, mass, y, massErr, yErr);
    TGraphErrors *gr_mass_expol[tot_reso];
    TGraphErrors *gr_width_expol[tot_reso];
    TGraphErrors *gr_mass_boltzman[tot_reso];
    TGraphErrors *gr_width_boltzman[tot_reso];
    TGraphErrors *gr_mass_coherent[tot_reso];
    TGraphErrors *gr_width_coherent[tot_reso];
    TGraphErrors *gr_HERA_mass[tot_reso];
    TGraphErrors *gr_HERA_width[tot_reso];
    TGraphErrors *gr_L3_mass[tot_reso];
    TGraphErrors *gr_L3_width[tot_reso];
    TGraphErrors *gr_BES_mass[tot_reso];
    TGraphErrors *gr_BES_width[tot_reso];
    TGraphErrors *gr_WA02_mass[tot_reso];
    TGraphErrors *gr_WA02_width[tot_reso];
    TGraphErrors *gr_BESIII_mass[tot_reso];
    TGraphErrors *gr_BESIII_width[tot_reso];
    TLine *line[tot_reso];

    double pdgmasses[tot_reso] = {f1270Mass, a1320Mass, f1525Mass, f1710Mass};
    double pdgmasserrors[tot_reso] = {f1270MassErr, a1320MassErr, f1525MassErr, f1710MassErr};
    double pdgwidths[tot_reso] = {f1270Width, a1320Width, f1525Width, f1710Width};
    double pdgwidtherrors[tot_reso] = {f1270WidthErr, a1320WidthErr, f1525WidthErr, f1710WidthErr};
    float limits_mass[tot_reso][2] = {{1.20, 1.379}, {1.241, 1.419}, {1.495, 1.548}, {1.63, 1.81}};
    float limits_width[tot_reso][2] = {{0.02, 0.38}, {0.015, 0.33}, {0.05, 0.13}, {0.06, 0.4}};
    string resonances[tot_reso] = {"f_{2}(1270)", "a_{2}(1320)^{0}", "f'_{2}(1525)", "f_{0}(1710)"};
    TBox *box[tot_reso];

    for (int i = 0; i < tot_reso; i++)
    {
        gr_mass_expol[i] = new TGraphErrors();
        gr_width_expol[i] = new TGraphErrors();
        gr_mass_boltzman[i] = new TGraphErrors();
        gr_width_boltzman[i] = new TGraphErrors();
        gr_mass_coherent[i] = new TGraphErrors();
        gr_width_coherent[i] = new TGraphErrors();
        gr_HERA_mass[i] = new TGraphErrors();
        gr_HERA_width[i] = new TGraphErrors();
        gr_L3_mass[i] = new TGraphErrors();
        gr_L3_width[i] = new TGraphErrors();
        gr_BES_mass[i] = new TGraphErrors();
        gr_BES_width[i] = new TGraphErrors();
        gr_WA02_mass[i] = new TGraphErrors();
        gr_WA02_width[i] = new TGraphErrors();
        gr_BESIII_mass[i] = new TGraphErrors();
        gr_BESIII_width[i] = new TGraphErrors();

        SetGrapherrorStyle(gr_mass_expol[i]);
        SetGrapherrorStyle(gr_width_expol[i]);
        SetGrapherrorStyle(gr_mass_boltzman[i]);
        SetGrapherrorStyle(gr_width_boltzman[i]);
        SetGrapherrorStyle(gr_mass_coherent[i]);
        SetGrapherrorStyle(gr_width_coherent[i]);
        SetGrapherrorStyle(gr_HERA_mass[i]);
        SetGrapherrorStyle(gr_HERA_width[i]);
        SetGrapherrorStyle(gr_L3_mass[i]);
        SetGrapherrorStyle(gr_L3_width[i]);
        SetGrapherrorStyle(gr_BES_mass[i]);
        SetGrapherrorStyle(gr_BES_width[i]);
        SetGrapherrorStyle(gr_WA02_mass[i]);
        SetGrapherrorStyle(gr_WA02_width[i]);
        SetGrapherrorStyle(gr_BESIII_mass[i]);
        SetGrapherrorStyle(gr_BESIII_width[i]);

        // gr_HERA_mass[i]->SetPoint(i, HERA_mass[i], 1);
        // gr_HERA_mass[i]->SetPointError(i, HERA_mass_errors[i], 0);
        // gr_HERA_width[i]->SetPoint(i, HERA_width[i], 1);
        // gr_HERA_width[i]->SetPointError(i, HERA_width_errors[i], 0);

        // gr_L3_mass[i]->SetPoint(i, L3_mass[i], 2);
        // gr_L3_mass[i]->SetPointError(i, L3_mass_errors[i], 0);
        // gr_L3_width[i]->SetPoint(i, L3_width[i], 2);
        // gr_L3_width[i]->SetPointError(i, L3_width_errors[i], 0);

        gr_mass_boltzman[i]->SetPoint(i, mass_boltzman[i], 1);
        gr_mass_boltzman[i]->SetPointError(i, mass_boltzman_err[i], 0);
        gr_width_boltzman[i]->SetPoint(i, width_boltzman[i], 1);
        gr_width_boltzman[i]->SetPointError(i, width_boltzman_err[i], 0);

        gr_mass_expol[i]->SetPoint(i, mass_expol[i], 2);
        gr_mass_expol[i]->SetPointError(i, mass_expol_err[i], 0);
        gr_width_expol[i]->SetPoint(i, width_expol[i], 2);
        gr_width_expol[i]->SetPointError(i, width_expol_err[i], 0);

        // gr_BES_mass[i]->SetPoint(i, BES_mass[i], 3);
        // gr_BES_mass[i]->SetPointError(i, BES_mass_errors[i], 0);
        // gr_BES_width[i]->SetPoint(i, BES_width[i], 3);
        // gr_BES_width[i]->SetPointError(i, BES_width_errors[i], 0);

        // gr_WA02_mass[i]->SetPoint(i, WA02_mass[i], 4);
        // gr_WA02_mass[i]->SetPointError(i, WA02_mass_errors[i], 0);
        // gr_WA02_width[i]->SetPoint(i, WA02_width[i], 4);
        // gr_WA02_width[i]->SetPointError(i, WA02_width_errors[i], 0);

        // gr_BESIII_mass[i]->SetPoint(i, BESIII_mass[i], 5);
        // gr_BESIII_mass[i]->SetPointError(i, BESIII_mass_errors[i], 0);
        // gr_BESIII_width[i]->SetPoint(i, BESIII_width[i], 5);
        // gr_BESIII_width[i]->SetPointError(i, BESIII_width_errors[i], 0);

        // gr_mass_boltzman[i]->SetPoint(i, mass_boltzman[i], 6);
        // gr_mass_boltzman[i]->SetPointError(i, mass_boltzman_err[i], 0);
        // gr_width_boltzman[i]->SetPoint(i, width_boltzman[i], 6);
        // gr_width_boltzman[i]->SetPointError(i, width_boltzman_err[i], 0);

        // gr_mass_expol[i]->SetPoint(i, mass_expol[i], 7);
        // gr_mass_expol[i]->SetPointError(i, mass_expol_err[i], 0);
        // gr_width_expol[i]->SetPoint(i, width_expol[i], 7);
        // gr_width_expol[i]->SetPointError(i, width_expol_err[i], 0);

        // gr_mass_coherent[i]->SetPoint(i, mass_coherent[i], 8);
        // gr_mass_coherent[i]->SetPointError(i, mass_coherent_err[i], 0);
        // gr_width_coherent[i]->SetPoint(i, width_coherent[i], 8);
        // gr_width_coherent[i]->SetPointError(i, width_coherent_err[i], 0);
    }

    plot(gr_mass_expol, gr_mass_boltzman, gr_HERA_mass, gr_mass_coherent, gr_L3_mass, gr_BES_mass, gr_WA02_mass, gr_BESIII_mass, line, box, limits_mass, resonances, "Mass", pdgmasses, pdgmasserrors);

    plot(gr_width_expol, gr_width_boltzman, gr_HERA_width, gr_width_coherent, gr_L3_width, gr_BES_width, gr_WA02_width, gr_BESIII_width, line, box, limits_width, resonances, "Width", pdgwidths, pdgwidtherrors);

} // end of main function

void plot(TGraphErrors *gr_mass_expol[4], TGraphErrors *gr_mass_boltzman[4], TGraphErrors *gr_HERA_mass[4], TGraphErrors *gr_mass_coherent[4], TGraphErrors *gr_mass_L3[4], TGraphErrors *gr_mass_BES[4], TGraphErrors *gr_mass_WA02[4], TGraphErrors *gr_mass_BESIII[4], TLine *line[4], TBox *box[4], float limits[4][2], string resonances[4], string xaxistitle, double *pdg, double *pdgerr)
{
    for (int i = 0; i < tot_reso; i++)
    {

        TCanvas *c = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c, 0.15, 0.03, 0.05, 0.14);
        gr_mass_expol[i]->SetTitle(Form("; %s (GeV/c)", xaxistitle.c_str()));
        gr_mass_expol[i]->SetMarkerStyle(20);
        gr_mass_expol[i]->SetMarkerColor(1);
        gr_mass_expol[i]->SetLineColor(1);

        gr_mass_boltzman[i]->SetMarkerStyle(21);
        gr_mass_boltzman[i]->SetMarkerColor(2);
        gr_mass_boltzman[i]->SetLineColor(2);

        gr_mass_coherent[i]->SetMarkerStyle(22);
        gr_mass_coherent[i]->SetMarkerColor(4);
        gr_mass_coherent[i]->SetLineColor(4);

        gr_HERA_mass[i]->SetMarkerStyle(23);
        gr_HERA_mass[i]->SetMarkerColor(6);
        gr_HERA_mass[i]->SetLineColor(6);

        gr_mass_L3[i]->SetMarkerStyle(24);
        gr_mass_L3[i]->SetMarkerColor(9);
        gr_mass_L3[i]->SetLineColor(9);

        gr_mass_BES[i]->SetMarkerStyle(25);
        gr_mass_BES[i]->SetMarkerColor(28);
        gr_mass_BES[i]->SetLineColor(28);

        gr_mass_WA02[i]->SetMarkerStyle(26);
        gr_mass_WA02[i]->SetMarkerColor(46);
        gr_mass_WA02[i]->SetLineColor(46);

        gr_mass_BESIII[i]->SetMarkerStyle(27);
        gr_mass_BESIII[i]->SetMarkerColor(42);
        gr_mass_BESIII[i]->SetLineColor(42);

        gr_mass_expol[i]->Draw("AP");
        gr_mass_boltzman[i]->Draw("P SAME");
        // gr_mass_coherent[i]->Draw("P SAME");
        // gr_HERA_mass[i]->Draw("P SAME");
        // gr_mass_L3[i]->Draw("P SAME");
        // gr_mass_BES[i]->Draw("P SAME");
        // gr_mass_WA02[i]->Draw("P SAME");
        // gr_mass_BESIII[i]->Draw("P SAME");

        line[i] = new TLine(pdg[i], 0, pdg[i], 2.5);
        line[i]->SetLineColor(kRed);
        line[i]->SetLineWidth(2);
        line[i]->Draw("same");
        box[i] = new TBox(pdg[i] - pdgerr[i], 0, pdg[i] + pdgerr[i], 2.5);
        box[i]->SetFillColor(kYellow);
        box[i]->SetFillStyle(3001);
        box[i]->Draw("same");

        TLegend *leg = new TLegend(0.60, 0.65, 0.94, 0.9);
        leg->SetFillStyle(0);
        // leg->SetBorderSize(0);
        // leg->AddEntry(gr_mass_coherent[i], "Coherent fit", "lp");
        leg->AddEntry(gr_mass_expol[i], "4rBW + Expol 1 fit", "lp");
        leg->AddEntry(gr_mass_boltzman[i], "3rBW + Expol 1 fit", "lp");
        // leg->AddEntry(gr_mass_BESIII[i], "BESIII expt.", "lp");
        // leg->AddEntry(gr_mass_WA02[i], "WA02 expt.", "lp");
        // leg->AddEntry(gr_mass_BES[i], "BES expt.", "lp");
        // leg->AddEntry(gr_mass_L3[i], "L3 expt.", "lp");
        // leg->AddEntry(gr_HERA_mass[i], "HERA expt.", "lp");
        leg->AddEntry(box[i], Form("PDG %s", resonances[i].c_str()), "f");
        leg->Draw("same");
        TLegend *leg2 = new TLegend(0.60, 0.65, 0.94, 0.9); // this is created to show line and box in the legend
        leg2->SetFillStyle(0);
        leg2->SetBorderSize(0);
        leg2->AddEntry((TObject *)0, "", "");
        leg2->AddEntry((TObject *)0, "", "");
        // leg2->AddEntry((TObject *)0, "", "");
        // leg2->AddEntry((TObject *)0, "", "");
        // leg2->AddEntry((TObject *)0, "", "");
        // leg2->AddEntry((TObject *)0, "", "");
        // leg2->AddEntry((TObject *)0, "", "");
        // leg2->AddEntry((TObject *)0, "", "");
        leg2->AddEntry(line[i], " ", "e");
        leg2->Draw("same");

        gr_mass_expol[i]->GetXaxis()->SetLimits(limits[i][0], limits[i][1]); // Set X-axis from 1 to 2 GeV/c
        gr_mass_expol[i]->GetYaxis()->SetRangeUser(0, 3);                    // Change here for y axis range
        gr_mass_expol[i]->GetYaxis()->SetNdivisions(3);
        gr_mass_expol[i]->GetYaxis()->SetBinLabel(1, "");

        // Get the minimum x-axis value
        double xMin = gr_mass_expol[i]->GetXaxis()->GetXmin();

        // Define a small offset to position the text just to the left of the y-axis
        double offset = 0.15 * (gr_mass_expol[i]->GetXaxis()->GetXmax() - xMin);

        TLatex latex;
        latex.SetTextSize(0.04);
        // latex.DrawLatex(xMin - offset, 1, "HERA Experiment");
        // latex.DrawLatex(xMin - offset, 2, "L3 Experiment");
        // latex.DrawLatex(xMin - offset, 3, "BES Expt.");
        // latex.DrawLatex(xMin - offset, 4, "WA02 Expt.");
        // latex.DrawLatex(xMin - offset, 5, "BESIII Expt.");
        // latex.DrawLatex(xMin - offset, 6, "3BW + Boltzman");
        // latex.DrawLatex(xMin - offset, 7, "3BW + Expol fit");
        // latex.DrawLatex(xMin - offset, 8, "Coherent fit");

        latex.DrawLatex(xMin - offset, 1, "3rBW");
        latex.DrawLatex(xMin - offset, 2, "4rBW");

        c->SaveAs((save_folder_str + xaxistitle + resonances[i] + ".png").c_str());
    }
}
