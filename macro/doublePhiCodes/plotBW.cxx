void plotBW()
{
    // Resonance parameters
    const double mass = 2.714;  // GeV/c^2
    const double width = 0.012; // GeV

    // Breit-Wigner function
    TF1 *bw = new TF1("bw",
                      "[2]*TMath::BreitWigner(x,[0],[1])",
                      2.65, 2.78);

    bw->SetParameters(mass, width, 1.0); // mass, width, normalization

    // Style
    bw->SetTitle(";Invariant Mass (GeV/#it{c}^{2});Arbitrary Units");
    bw->SetLineColor(kBlue + 2);
    bw->SetLineWidth(3);

    TCanvas *c = new TCanvas("c", "Breit-Wigner", 800, 600);
    bw->Draw();

    // Draw vertical line at resonance mass
    double ymax = bw->GetMaximum();

    TLine *l = new TLine(mass, 0, mass, ymax);
    l->SetLineStyle(2);
    l->SetLineColor(kRed);
    l->SetLineWidth(2);
    l->Draw("same");

    // Label
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.18, 0.85,
                    Form("M = %.3f GeV/#it{c}^{2}", mass));
    latex.DrawLatex(0.18, 0.80,
                    Form("#Gamma = %.0f MeV", width * 1000));

    c->Update();
}