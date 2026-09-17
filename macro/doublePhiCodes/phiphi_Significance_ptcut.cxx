#include "src/style.h"
using namespace std;

void phiphi_Significance_ptcut()
{
    // List of pT cuts to scan
    std::vector<double> ptCuts = {6.0, 7.0, 8.0, 9.0, 10.0};

    std::vector<double> x_pt;
    std::vector<double> y_sigTotal;
    std::vector<double> y_sigUncorr;
    std::vector<double> y_sigCorr;

    // Loop through files and extract data
    for (double pt : ptCuts)
    {
        // Construct filename matching pattern (e.g., YieldRatios_pt6.0.txt)
        std::ostringstream filename;
        filename << "/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/Template/pTvariation/YieldRatios_pt" << std::fixed << std::setprecision(1) << pt << ".txt";

        std::ifstream inFile(filename.str());
        if (!inFile.is_open())
        {
            std::cerr << "Warning: Could not open file " << filename.str() << std::endl;
            continue;
        }

        std::string line;
        double ratioTotal = -1.0;
        double ratioUncorr = -1.0;
        double ratioCorr = -1.0;

        while (std::getline(inFile, line))
        {
            if (line.find("Signal / Total Background") != std::string::npos)
            {
                std::size_t colonPos = line.find(":");
                if (colonPos != std::string::npos)
                {
                    ratioTotal = std::stod(line.substr(colonPos + 1));
                }
            }
            if (line.find("Signal / Uncorrelated Background") != std::string::npos)
            {
                std::size_t colonPos = line.find(":");
                if (colonPos != std::string::npos)
                {
                    ratioUncorr = std::stod(line.substr(colonPos + 1));
                }
            }
            if (line.find("Signal / Correlated Background") != std::string::npos)
            {
                std::size_t colonPos = line.find(":");
                if (colonPos != std::string::npos)
                {
                    ratioCorr = std::stod(line.substr(colonPos + 1));
                }
            }
        }
        inFile.close();

        if (ratioTotal > 0 && ratioUncorr > 0 && ratioCorr > 0)
        {
            x_pt.push_back(pt);
            y_sigTotal.push_back(ratioTotal);
            y_sigUncorr.push_back(ratioUncorr);
            y_sigCorr.push_back(ratioCorr);
            std::cout << "Successfully parsed " << filename.str()
                      << " -> Sig/Total: " << ratioTotal
                      << ", Sig/Uncorr: " << ratioUncorr
                      << ", Sig/Corr: " << ratioCorr << std::endl;
        }
    }

    if (x_pt.empty())
    {
        std::cerr << "Error: No valid data extracted from files!" << std::endl;
        return;
    }

    // Create Graphs
    int nPoints = x_pt.size();
    TGraph *grTotal = new TGraph(nPoints, &x_pt[0], &y_sigTotal[0]);
    TGraph *grUncorr = new TGraph(nPoints, &x_pt[0], &y_sigUncorr[0]);
    TGraph *grCorr = new TGraph(nPoints, &x_pt[0], &y_sigCorr[0]);

    // Graph styling
    SetGraphStyle(grTotal);
    grTotal->SetMarkerStyle(20); // Solid circle
    grTotal->SetMarkerSize(1.2);
    grTotal->SetMarkerColor(kRed + 1);
    grTotal->SetLineColor(kRed + 1);
    grTotal->SetLineWidth(2);

    SetGraphStyle(grUncorr);
    grUncorr->SetMarkerStyle(21); // Solid square
    grUncorr->SetMarkerSize(1.2);
    grUncorr->SetMarkerColor(kBlue + 1);
    grUncorr->SetLineColor(kBlue + 1);
    grUncorr->SetLineWidth(2);

    SetGraphStyle(grCorr);
    grCorr->SetMarkerStyle(22); // Solid triangle
    grCorr->SetMarkerSize(1.2);
    grCorr->SetMarkerColor(kGreen + 2);
    grCorr->SetLineColor(kGreen + 2);
    grCorr->SetLineWidth(2);

    // Canvas setup
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1", "Signal to Background Ratios vs pT Cut", 720, 720);
    SetCanvasStyle(c1, 0.16, 0.03, 0.05, 0.15);
    gPad->SetGrid();

    // Draw graph with frame axes
    grTotal->SetTitle("; #it{p}_{T}^{#phi#phi} Cut (GeV/ #it{c});Yield Ratio");
    grTotal->GetXaxis()->SetTitleSize(0.045);
    grTotal->GetYaxis()->SetTitleSize(0.045);
    grTotal->GetXaxis()->SetTitleOffset(1.25);
    grTotal->GetYaxis()->SetTitleOffset(1.6);
    grTotal->GetXaxis()->SetNdivisions(505);

    // // Set plot axis limits
    // double maxVal = *std::max_element(y_sigUncorr.begin(), y_sigUncorr.end()) * 1.05;
    // grTotal->GetYaxis()->SetRangeUser(0.0, maxVal);
    grTotal->GetYaxis()->SetRangeUser(0.0, 0.068); // Set a fixed range for better comparison

    grTotal->Draw("APL"); // A = Axes, P = Points, L = Line
    // grUncorr->Draw("PL SAME");
    // grCorr->Draw("PL SAME");

    // Add Legend
    TLegend *leg = new TLegend(0.35, 0.78, 0.5, 0.88);
    leg->SetBorderSize(0);
    // leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    // leg->SetHeader("2025 triggered dataset");
    leg->AddEntry(grTotal, "Signal / Total Background", "p");
    // leg->AddEntry(grUncorr, "Signal / Uncorrelated Bkg (non-SS)", "p");
    // leg->AddEntry(grCorr, "Signal / Correlated Bkg (SS)", "p");
    leg->Draw();

    c1->Update();
    c1->SaveAs("/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/Template/pTvariation/SignalToBkg_vs_pTCut.png");
    // c1->SaveAs("SignalToBkg_vs_pTCut.pdf");

    TCanvas *c2 = new TCanvas("c2", "Signal by uncorrelated background vs pT Cut", 720, 720);
    SetCanvasStyle(c2, 0.16, 0.03, 0.05, 0.15);
    gPad->SetGrid();
    grUncorr->SetTitle("; #it{p}_{T}^{#phi#phi} Cut (GeV/ #it{c});Yield Ratio");
    SetGraphStyle(grUncorr);
    grUncorr->GetXaxis()->SetTitleSize(0.045);
    grUncorr->GetYaxis()->SetTitleSize(0.045);
    grUncorr->GetXaxis()->SetTitleOffset(1.25);
    grUncorr->GetYaxis()->SetTitleOffset(1.6);
    grUncorr->GetXaxis()->SetNdivisions(505);
    grUncorr->GetYaxis()->SetRangeUser(0.0, 0.123); // Set a fixed range for better comparison
    grUncorr->Draw("APL");

    TLegend *leg2 = new TLegend(0.35, 0.8, 0.8, 0.90);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.033);
    leg2->AddEntry(grUncorr, "Signal / Uncorrelated Background", "p");
    leg2->Draw();
    c2->Update();
    c2->SaveAs("/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/Template/pTvariation/SignalToUncorrBkg_vs_pTCut.png");

    double width[] = {67.6, 55.4, 51.9, 33.7, 40.9};
    double widthErr[] = {11.4, 11.6, 10.6, 9.2, 11.6};

    TGraphErrors *grWidth = new TGraphErrors(nPoints, &x_pt[0], width, nullptr, widthErr);
    SetGraphErrorStyle(grWidth);
    grWidth->SetMarkerStyle(20);
    grWidth->SetMarkerSize(1.2);
    grWidth->SetMarkerColor(kGreen + 2);
    grWidth->SetLineColor(kGreen + 2);
    grWidth->SetLineWidth(2);

    TCanvas *cWidth = new TCanvas("cWidth", "Width vs pT Cut", 720, 720);
    SetCanvasStyle(cWidth, 0.15, 0.03, 0.05, 0.15);
    gPad->SetGrid();
    grWidth->SetTitle("; #it{p}_{T}^{#phi#phi} Cut (GeV/ #it{c});Width (MeV/#it{c}^{2})");
    grWidth->GetXaxis()->SetTitleSize(0.045);
    grWidth->GetYaxis()->SetTitleSize(0.045);
    grWidth->GetXaxis()->SetTitleOffset(1.25);
    grWidth->GetYaxis()->SetTitleOffset(1.45);
    grWidth->SetMinimum(8);
    grWidth->SetMaximum(95);
    grWidth->GetXaxis()->SetNdivisions(505);
    grWidth->Draw("AP");

    TLegend *legWidth = new TLegend(0.65, 0.8, 0.90, 0.88);
    legWidth->SetBorderSize(0);
    legWidth->SetTextSize(0.03);
    legWidth->AddEntry(grWidth, "#phi#phi signal width", "p");
    legWidth->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.03);
    lat.SetTextFont(42);
    // lat.DrawLatex(0.18, 0.85, "2026 triggered dataset");
    cWidth->Update();
    cWidth->SaveAs("/home/sawan/Storage/check_k892/output/doublePhi/LocalTests/Template/pTvariation/Width_vs_pTCut.png");
}