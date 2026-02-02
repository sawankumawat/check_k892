
#include "../src/style.h"
using namespace std;

void INEL_crossSection()
{
    //=======================Paper (PHYSICAL REVIEW C 97, 054910 (2018))===========================
    //============Get value without fitting (fit parameters given in the paper)=================
    double sqrt_s_TeV = 13.6;

    // Convert to s in GeV^2
    double s_GeV2 = pow(sqrt_s_TeV * 1e3, 2);

    double L = log(s_GeV2); // natural log

    // ---- Parameters from Table I ----

    // ln(s)
    double A1 = -3.33, B1 = 4.195, n1 = 1.0;

    // ln^2(s)
    double A2 = 25.0, B2 = 0.146, n2 = 2.0;

    // ln^n(s)
    double A3 = 29.8, B3 = 0.038, n3 = 2.43;

    // ---- Evaluate ----
    double sigma_ln = A1 + B1 * pow(L, n1);
    double sigma_ln2 = A2 + B2 * pow(L, n2);
    double sigma_lnn = A3 + B3 * pow(L, n3);

    double model_err = fabs(sigma_ln - sigma_lnn) / 2.4;

    std::cout << "ln(s) value used = " << L << std::endl;

    std::cout << "σ_NN(13.6 TeV) ln^2(s) = " << sigma_ln2 << " +/- " << model_err << " mb\n";

    /*
    //================To get from fits (not recommended since the values are taken from the fit only)========
    const int N = 18;

    // sqrt(s) in TeV
    double s[N] = {0.2, 0.9, 2.76, 5.02, 5.44, 5.5, 7.0, 8.0, 8.16, 8.8,
                   10.6, 13.0, 14.0, 17.0, 27.0, 39.0, 63.0, 100.0};

    // sigma_NN in mb
    double sigma[N] = {41.6, 52.2, 61.8, 67.6, 68.4, 68.5, 70.9, 72.3, 72.5, 73.3,
                       75.3, 77.6, 78.4, 80.6, 86.0, 90.5, 96.5, 102.6};

    // uncertainties
    double err[N] = {0.6, 1.0, 0.9, 0.6, 0.5, 0.5, 0.4, 0.5, 0.5, 0.6,
                     0.7, 1.0, 1.1, 1.5, 2.4, 3.3, 4.6, 6.0};

    double xerr[N] = {0};

    // Graph
    TGraphErrors *gr = new TGraphErrors(N, s, sigma, xerr, err);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.1);
    gr->SetTitle("NN Inelastic Cross Section; #sqrt{s} (TeV); #sigma_{NN} (mb)");

    // ---------------- FITS ----------------

    // // n = 1  →  A + B ln(s)
    // TF1 *fit_n1 = new TF1("fit_n1", "[0] + [1]*pow(log(x),1)", 0.2, 120);
    // fit_n1->SetParameters(40, 10);
    // fit_n1->SetLineColor(kRed);

    // // n = 2  →  A + B ln^2(s)
    // TF1 *fit_n2 = new TF1("fit_n2", "[0] + [1]*pow(log(x),2)", 0.2, 120);
    // fit_n2->SetParameters(40, 1);
    // fit_n2->SetLineColor(kBlue);

    // n free → A + B ln^n(s)
    TF1 *fit_free = new TF1("fit_free", "[0] + [1]*pow(log(x), [2])", 0.2, 120);
    fit_free->SetParameters(55, 9, 1.0); // initial guess
    // fit_free->SetParLimits(2, 1.0, 2.0); // limit n between 1.0 and 2.0
    fit_free->FixParameter(2, 2.0); //
    fit_free->SetParNames("A", "B", "n");
    fit_free->SetLineColor(kGreen + 2);

    // Perform fits
    // gr->Fit(fit_n1, "REBMS");
    // gr->Fit(fit_n2, "REBMS+");
    gr->Fit(fit_free, "RI");

    // ---------------- PLOT ----------------
    TCanvas *c1 = new TCanvas("c1", "sigmaNN fits", 900, 700);
    gr->Draw("AP");
    // fit_n1->Draw("same");
    // fit_n2->Draw("same");
    fit_free->Draw("same");

    // TLegend *leg = new TLegend(0.15, 0.60, 0.48, 0.85);
    // leg->AddEntry(gr, "Data", "p");
    // leg->AddEntry(fit_n1, "A + B ln(s)  (n=1)", "l");
    // leg->AddEntry(fit_n2, "A + B ln^{2}(s) (n=2)", "l");
    // leg->AddEntry(fit_free, "A + B ln^{n}(s) (n free)", "l");
    // leg->Draw();

    // // ---------------- EXTRAPOLATION ----------------
    // double s_eval = 13.6;

    // double val_n1 = fit_n1->Eval(s_eval);
    // double val_n2 = fit_n2->Eval(s_eval);
    // double val_free = fit_free->Eval(s_eval);

    // // Model uncertainties (paper method)
    // double err_n1_n2 = fabs(val_n1 - val_n2);
    // double err_n1_free = fabs(val_n1 - val_free);

    // // Choose reference = n=1 (as in literature tables)
    // double sigma_final = val_n1;

    // // Conservative uncertainty
    // double total_err = TMath::Max(err_n1_n2, err_n1_free);

    // // ---------------- OUTPUT ----------------
    // std::cout << "\n================ FIT RESULTS ================\n";
    // std::cout << "n=1  fit  σ(13.6) = " << val_n1 << " mb\n";
    // std::cout << "n=2  fit  σ(13.6) = " << val_n2 << " mb\n";
    // std::cout << "n free fit σ(13.6) = " << val_free << " mb\n";
    // std::cout << "---------------------------------------------\n";
    // std::cout << "Uncertainty |n=1 - n=2|     = " << err_n1_n2 << " mb\n";
    // std::cout << "Uncertainty |n=1 - n_free|  = " << err_n1_free << " mb\n";
    // std::cout << "---------------------------------------------\n";
    // std::cout << "FINAL RESULT:\n";
    // std::cout << "σ_NN(13.6 TeV) = " << sigma_final
    //           << " ± " << total_err << " mb\n";
    // std::cout << "=============================================\n";

    // std::cout << "\nFree-n fit result: n = "
    //           << fit_free->GetParameter(2)
    //           << " ± "
    //           << fit_free->GetParError(2) << std::endl;
    */
}
