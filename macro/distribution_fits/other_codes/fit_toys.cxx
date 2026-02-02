#include <iostream>
using namespace std;
#include "../src/style.h"

void fit_toys()
{
    // Load your toy histogram
    TFile *f = new TFile("toy_significance_NullModel.root"); // adjust filename
    if (f->IsZombie())
    {
        cerr << "Error: File could not be opened!" << endl;
        return;
    }
    TH1F *h = (TH1F *)f->Get("Normalized_toy_null"); // histogram of q0 from toys
    h->SetBinContent(h->GetXaxis()->FindBin(-0.1), -999);

    // --- Fit with chi2 distribution (1 d.o.f.)
    TF1 *fchi = new TF1("fchi", "[0]*ROOT::Math::chisquared_pdf(x,1)", 0, 70);
    fchi->SetParameter(0, 1.0);
    h->Fit(fchi, "R");

    // Canvas setup
    TCanvas *c = new TCanvas("c", "q0 distribution", 720, 720);
    SetCanvasStyle(c, 0.14, 0.06, 0.03, 0.13);
    SetHistoQA(h);
    h->SetLineColor(kBlue + 2);
    h->SetMarkerStyle(20);
    // h->SetTitle("Distribution of q_{0} from Toys; q_{0} = -2#Delta ln L; Normalized Frequency");
    h->SetTitle("; q_{0} = -2#Delta ln L; Normalized Counts");
    h->SetMinimum(-0.03);
    h->SetMaximum(0.9);
    h->SetStats(0);
    h->Draw("E");

    // Draw fit
    fchi->SetLineColor(kRed);
    fchi->SetLineWidth(2);
    fchi->Draw("SAME");

    // Observed test statistic
    double q0_obs = 1788.99;

    // Conservative significance from toys
    int Ntoys = h->GetEntries();
    int N_exceed = 0;
    for (int i = 1; i <= h->GetNbinsX(); i++)
    {
        if (h->GetBinLowEdge(i) >= q0_obs)
        {
            N_exceed += h->GetBinContent(i) * Ntoys;
        }
    }

    double pval;
    if (N_exceed == 0)
    {
        // conservative estimate if no toys exceed
        pval = 1.0 / Ntoys;
    }
    else
    {
        pval = double(N_exceed) / double(Ntoys);
    }

    double Z = ROOT::Math::gaussian_quantile_c(pval, 1); // Convert p -> Z

    // Draw observed line
    TLine *line = new TLine(q0_obs, 0, q0_obs, h->GetMaximum());
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->SetLineStyle(2);
    line->Draw("SAME");

    // Legend
    TLegend *leg = new TLegend(0.52, 0.75, 0.85, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.033);
    leg->AddEntry(h, "Toy MC (null hypothesis)", "lep");
    leg->AddEntry(fchi, "#chi^{2}(1) fit", "l");
    leg->AddEntry(line, Form("Observed q_{0} = %.2f", q0_obs), "l");
    leg->Draw();

    // Text with significance
    TLatex latex;
    latex.SetNDC();
    // latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.18, 0.88, "q_{0}^{obs} > All toys");
    // latex.DrawLatex(0.16, 0.83, Form("p-value < %.2g", pval));
    // latex.DrawLatex(0.16, 0.78, Form("Significance Z #geq %.2f #sigma", Z));

    c->SaveAs("toy_model_q0_fit.png");
}
