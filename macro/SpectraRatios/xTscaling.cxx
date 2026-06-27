#include <iostream>
#include <iomanip>
#include "../src/style.h"
using namespace std;

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

TH1D *GetHisto(TFile *f, const string &name)
{
    TH1D *histo = (TH1D *)f->Get(name.c_str());

    if (!histo || histo == nullptr)
    {
        cout << "Error: histo " << name << " not found in file " << f->GetName() << endl;
        return nullptr;
    }

    SetHistoQA(histo);
    histo->SetTitle(0);
    return histo;
}

TGraphErrors *GetGraph(TFile *f, const string &name)
{
    TGraphErrors *graph = (TGraphErrors *)f->Get(name.c_str());

    if (!graph || graph == nullptr)
    {
        cout << "Error: graph " << name << " not found in file " << f->GetName() << endl;
        return nullptr;
    }

    SetGraphErrorStyle(graph);
    graph->SetTitle(0);
    return graph;
}
TGraphErrors *HistToXTGraph(TH1D *hStat, TH1D *hRelSys, double sqrts, double sigmaINEL, double nAverage = 1.0, bool alreadyNormalizedYield = false, float AverageKstar = 1.0);
TGraphErrors *GraphToXTGraph(TGraphErrors *gInput, TH1D *hStat, TH1D *hSys, double sqrts, double sigmaINEL, double nAverage = 1.0, bool alreadyNormalizedYield = false, float AverageKstar = 1.0);
void GetGraphRange(TGraph *g, double &xmin, double &xmax);
TGraph *CalculateN(TGraphErrors *g1, TGraphErrors *g2, double sqrts1, double sqrts2);

void xTscaling()
{
    // gStyle->SetOptFit(1111);
    string path1 = "../../output/kstar/LHC22o_pass7/708297/kstarqa/hInvMass/corrected_spectra_0_120.root";
    string sysPath = "../../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/SystematicsPlots/SysUncert.root";
    TFile *fINEL = OpenFile(path1);
    TFile *fSystematics = OpenFile(sysPath);
    TH1D *hSpectraINELStat = GetHisto(fINEL, "mult_0-120/corrected_spectra_Integral_final");
    TH1D *hSpectraINELSys = (TH1D *)hSpectraINELStat->Clone("hSpectraINELSys");
    TH1D *hRelUncert = GetHisto(fSystematics, "hTotalSysSmoothed_0_100");

    string pathRun13TeV = "HEP_data/pp13TeV_INEL.root";
    TFile *fRun13TeV = OpenFile(pathRun13TeV);
    TGraphErrors *gSpectraRun13TeV = GetGraph(fRun13TeV, "Table 4/Graph1D_y1");
    TH1D *hStatError = GetHisto(fRun13TeV, "Table 4/Hist1D_y1_e1");
    TH1D *hTotalSysError = GetHisto(fRun13TeV, "Table 4/Hist1D_y1_e2");

    string pathRun7TeV = "HEP_data/HEPData_8TeV_INEL_Kstar_Phi.root";
    TFile *fRun7TeV = OpenFile(pathRun7TeV);
    TGraphErrors *gSpectraRun7TeV = GetGraph(fRun7TeV, "Table 3/Graph1D_y1");
    TH1D *hStatErrorRun7TeV = GetHisto(fRun7TeV, "Table 3/Hist1D_y1_e1");
    TH1D *hTotalSysErrorRun7TeV = GetHisto(fRun7TeV, "Table 3/Hist1D_y1_e2");

    string pathRun2p76TeV = "HEP_data/HEPData_2p76TeV.root";
    TFile *fRun2p76TeV = OpenFile(pathRun2p76TeV);
    TGraphErrors *gSpectraRun2p76TeV = GetGraph(fRun2p76TeV, "Table 1/Graph1D_y1");
    TH1D *hStatErrorRun2p76TeV = GetHisto(fRun2p76TeV, "Table 1/Hist1D_y1_e1");
    TH1D *hTotalSysErrorRun2p76TeV = GetHisto(fRun2p76TeV, "Table 1/Hist1D_y1_e2");

    const double sqrts136 = 13600.;
    const double sqrts13 = 13000.;
    const double sqrts7 = 7000.;
    const double sqrts276 = 2760.;

    const double sigma_inel136 = 77.904; // mb
    const double sigma_inel13 = 77.6;    // mb
    const double sigma_inel7 = 70.9;     // mb
    const double sigma_inel276 = 61.8;   // mb

    // I have scaled all other graphs than 13 TeV instead of scaling 13 TeV only, due to published result with K* + anit-kstar sum only.

    // For the calculation of N we do not need to scale with sqrt(s)^2
    auto gXT136_n = HistToXTGraph(hSpectraINELStat, hRelUncert, sqrts136, sigma_inel136, 1.0, false, 2.0); // Since in 13 TeV it is just sum and not average
    auto gXT13_n = GraphToXTGraph(gSpectraRun13TeV, hStatError, hTotalSysError, sqrts13, sigma_inel13, 1.0, false, 1.0);
    auto gXT7_n = GraphToXTGraph(gSpectraRun7TeV, hStatErrorRun7TeV, hTotalSysErrorRun7TeV, sqrts7, sigma_inel7, 1.0, false, 2.0);
    auto gXT276_n = GraphToXTGraph(gSpectraRun2p76TeV, hStatErrorRun2p76TeV, hTotalSysErrorRun2p76TeV, sqrts276, sigma_inel276, 1.0, true, 2.0);

    // auto gN13 = CalculateN(gXT136_n, gXT13_n, sqrts136, sqrts13); // No use, energy difference is too small
    auto gN7 = CalculateN(gXT136_n, gXT7_n, sqrts136, sqrts7);
    auto gN276 = CalculateN(gXT136_n, gXT276_n, sqrts136, sqrts276);
    auto gn7_276 = CalculateN(gXT7_n, gXT276_n, sqrts7, sqrts276);
    auto gn13_276 = CalculateN(gXT13_n, gXT276_n, sqrts13, sqrts276);
    auto gn13_7 = CalculateN(gXT13_n, gXT7_n, sqrts13, sqrts7);
    string labels[5] = {
        "n#left(#frac{Y(13.6)}{Y(7)}#right)",
        "n#left(#frac{Y(13.6)}{Y(2.76)}#right)",
        "n#left(#frac{Y(7)}{Y(2.76)}#right)",
        "n#left(#frac{Y(13)}{Y(2.76)}#right)",
        "n#left(#frac{Y(13)}{Y(7)}#right)"};

    vector<TGraph *> gNList = {gN7, gN276, gn7_276, gn13_276, gn13_7};

    TCanvas *cN = new TCanvas("cN", "cN", 1080, 720);
    SetCanvasStyle(cN, 0.13, 0.10, 0.01, 0.12);
    cN->Divide(3, 2);
    vector<double> nValues;

    vector<vector<double>> fitRanges = {
        {1.6e-3, 3.4e-3},  // 13.6 TeV / 7 TeV
        {1.8e-3, 4.5e-3},  // 13.6 TeV / 2.76 TeV
        {1.8e-3, 5.45e-3}, // 7 TeV / 2.76 TeV
        {1.6e-3, 2.8e-3},  // 13 TeV / 2.76 TeV
        {0.75e-3, 2.8e-3}   // 13 TeV / 7 TeV
    };

    // Fit all graphs and then calculate the average n(xT) value
    for (int i = 0; i < gNList.size(); ++i)
    {
        cN->cd(i + 1);
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.11);
        gPad->SetTopMargin(0.05);
        gPad->SetBottomMargin(0.11);
        TGraph *gN = gNList[i];
        gN->GetXaxis()->SetTitle("x_{T}=2p_{T}/#sqrt{s}");
        gN->GetYaxis()->SetTitle("n(x_{T})");
        gN->GetXaxis()->SetMaxDigits(3);
        // gN->GetXaxis()->SetNdivisions(505);
        gN->SetMarkerStyle(21);
        gN->SetMarkerColor(kRed);
        gN->SetLineColor(kRed);
        gN->SetMaximum(6.3);
        gN->Draw("AP");

        TF1 *pol0 = new TF1(Form("pol0_%d", i), "[0]", fitRanges[i][0], fitRanges[i][1]);
        pol0->SetLineColor(kBlue + 2);
        pol0->SetLineWidth(2);
        pol0->SetParameter(0, 4.7);
        gN->Fit(pol0, "REBMS");
        TLatex lat;
        lat.SetNDC();
        lat.SetTextSize(0.045);
        lat.SetTextFont(42);
        lat.DrawLatex(0.5, 0.4, labels[i].c_str());
        lat.DrawLatex(0.5, 0.3, Form("<n> = %.2f", pol0->GetParameter(0)));
        nValues.push_back(pol0->GetParameter(0));

        cout << "n value from fit is " << pol0->GetParameter(0) << endl;
    }
    cN->SaveAs("Plots/n_xT.png");

    double nAverage = std::accumulate(nValues.begin(), nValues.end(), 0.0) / nValues.size();
    cout << "Average n value is " << nAverage << endl;

    // nAverage = 4.53; // Run2 value

    auto gXT136 = HistToXTGraph(hSpectraINELStat, hRelUncert, sqrts136, sigma_inel136, pow(sqrts136, nAverage), false, 2.0); // Since in 13 TeV it is just sum and not average
    auto gXT13 = GraphToXTGraph(gSpectraRun13TeV, hStatError, hTotalSysError, sqrts13, sigma_inel13, pow(sqrts13, nAverage), false, 1.0);
    auto gXT7 = GraphToXTGraph(gSpectraRun7TeV, hStatErrorRun7TeV, hTotalSysErrorRun7TeV, sqrts7, sigma_inel7, pow(sqrts7, nAverage), false, 2.0);
    auto gXT276 = GraphToXTGraph(gSpectraRun2p76TeV, hStatErrorRun2p76TeV, hTotalSysErrorRun2p76TeV, sqrts276, sigma_inel276, pow(sqrts276, nAverage), true, 2.0);

    gXT7->SetMarkerSize(1.7);
    gXT276->SetMarkerSize(1.7);

    TCanvas *cXT = new TCanvas("cXT", "cXT", 720, 720);
    SetCanvasStyle(cXT, 0.16, 0.06, 0.01, 0.14);
    gPad->SetLogy();
    gPad->SetLogx();
    TH1D *hSpectraDummy = new TH1D("hSpectraDummy", "", 10000000, 1e-7, 1.0);
    SetHistoQA(hSpectraDummy);
    hSpectraDummy->GetXaxis()->SetRangeUser(9e-7, 9e-2);
    // hSpectraDummy->GetYaxis()->SetRangeUser(8e+11, 8e+23);
    hSpectraDummy->GetYaxis()->SetRangeUser(8e+9, 9.9e+21);
    hSpectraDummy->GetYaxis()->SetTitleOffset(1.6);
    hSpectraDummy->GetXaxis()->SetTitleOffset(1.4);
    hSpectraDummy->SetStats(0);
    hSpectraDummy->GetYaxis()->SetNdivisions(505);
    hSpectraDummy->GetYaxis()->SetTitle(Form("#sqrt{s}^{%.2f} d^{3}#sigma/dp^{3} (mb GeV^{-2}c^{3})", nAverage));
    hSpectraDummy->GetXaxis()->SetTitle("#it{x}_{T}");
    hSpectraDummy->Draw("pe");
    gXT136->SetMarkerStyle(20);
    gXT136->SetMarkerColor(kBlue);
    gXT136->SetLineColor(kBlue);
    gXT136->Draw("Pe same");
    gXT13->SetMarkerStyle(21);
    gXT13->SetMarkerColor(kRed);
    gXT13->SetLineColor(kRed);
    gXT13->Draw("P same");
    gXT7->SetMarkerStyle(22);
    gXT7->SetMarkerColor(kGreen + 2);
    gXT7->SetLineColor(kGreen + 2);
    gXT7->Draw("P same");
    gXT276->SetMarkerStyle(23);
    gXT276->SetMarkerColor(kOrange + 7);
    gXT276->SetLineColor(kOrange + 7);
    gXT276->Draw("P same");

    // Create a combined graph
    TGraphErrors *gXTAll = new TGraphErrors();
    int ip = 0;

    auto AddPoints = [&](TGraphErrors *g)
    {
        for (int i = 0; i < g->GetN(); ++i)
        {
            double x, y;
            g->GetPoint(i, x, y);
            gXTAll->SetPoint(ip, x, y);
            gXTAll->SetPointError(ip, g->GetErrorX(i), g->GetErrorY(i));
            ip++;
        }
    };

    AddPoints(gXT276);
    AddPoints(gXT7);
    AddPoints(gXT13);
    AddPoints(gXT136);

    TF1 *fitPL = new TF1("fitPL", "[0]*pow(x,[1])*pow(1+x,[2])", 2e-3, 9e-3);
    fitPL->SetLineColor(kBlack);
    fitPL->SetLineWidth(2);
    fitPL->SetParameters(1.8, -6.0, 5.0);
    fitPL->SetParLimits(0, -1000, 1000);
    fitPL->SetParLimits(1, -10.0, 10.0);
    fitPL->SetParLimits(2, -1000.0, 1000.0);
    gXTAll->Fit(fitPL, "REBMS");
    fitPL->Draw("same");

    TLegend *leg = new TLegend(0.65, 0.62, 0.93, 0.95);
    SetLegendStyle(leg);
    leg->SetTextSize(0.035);
    leg->AddEntry((TObject *)0, "K*^{0}", "");
    leg->AddEntry(gXT136, "13.6 TeV", "p");
    leg->AddEntry(gXT13, "13 TeV", "p");
    leg->AddEntry(gXT7, "7 TeV", "p");
    leg->AddEntry(gXT276, "2.76 TeV", "p");
    leg->AddEntry(fitPL, "Combined", "l");
    leg->AddEntry((TObject *)0, "power-law fit", "");
    leg->Draw();
    cXT->SaveAs("Plots/xT_scaling.png");
}

TGraphErrors *HistToXTGraph(TH1D *hStat, TH1D *hRelSys, double sqrts, double sigmaINEL, double nAverage = 1.0, bool alreadyNormalizedYield = false, float AverageKstar = 1.0)
{
    TGraphErrors *g = new TGraphErrors();

    for (int ibin = 1; ibin <= hStat->GetNbinsX(); ++ibin)
    {
        double pt = hStat->GetBinCenter(ibin);
        double dpt = hStat->GetBinWidth(ibin) / 2.0;

        double y = hStat->GetBinContent(ibin);
        double estat = hStat->GetBinError(ibin);

        double esys = hRelSys ? hRelSys->GetBinContent(ibin) * y : 0.0;
        double etot = sqrt(estat * estat + esys * esys);

        double invY = (alreadyNormalizedYield) ? (y * sigmaINEL * nAverage * AverageKstar) : (y * sigmaINEL * nAverage * AverageKstar / (2.0 * TMath::Pi() * pt));
        double invEY = (alreadyNormalizedYield) ? (etot * sigmaINEL * nAverage * AverageKstar) : (etot * sigmaINEL * nAverage * AverageKstar / (2.0 * TMath::Pi() * pt));

        double xT = 2.0 * pt / sqrts;
        double exT = 2.0 * dpt / sqrts;

        int p = g->GetN();

        g->SetPoint(p, xT, invY);
        g->SetPointError(p, exT, invEY);
    }
    SetGraphErrorStyle(g);
    g->SetMarkerSize(1.5);

    return g;
}

TGraphErrors *GraphToXTGraph(TGraphErrors *gInput, TH1D *hStat, TH1D *hSys, double sqrts, double sigmaINEL, double nAverage = 1.0, bool alreadyNormalizedYield = false, float AverageKstar = 1.0)
{
    TGraphErrors *g = new TGraphErrors();

    for (int i = 0; i < gInput->GetN(); ++i)
    {
        double pt, y;
        gInput->GetPoint(i, pt, y);

        double ex = gInput->GetErrorX(i);
        double stat = hStat->GetBinContent(i + 1);
        double sys = hSys->GetBinContent(i + 1);

        double total = sqrt(stat * stat + sys * sys);

        double invY = (alreadyNormalizedYield) ? (y * sigmaINEL * nAverage * AverageKstar) : (y * sigmaINEL * nAverage * AverageKstar / (2.0 * TMath::Pi() * pt));
        double invEY = (alreadyNormalizedYield) ? (total * sigmaINEL * nAverage * AverageKstar) : (total * sigmaINEL * nAverage * AverageKstar / (2.0 * TMath::Pi() * pt));

        double xT = 2.0 * pt / sqrts;
        double exT = 2.0 * ex / sqrts;

        g->SetPoint(i, xT, invY);
        g->SetPointError(i, exT, invEY);
    }
    SetGraphErrorStyle(g);
    g->SetMarkerSize(1.5);

    return g;
}

void GetGraphRange(TGraph *g, double &xmin, double &xmax)
{
    xmin = 1e30;
    xmax = -1e30;

    for (int i = 0; i < g->GetN(); i++)
    {
        double x, y;
        g->GetPoint(i, x, y);

        xmin = std::min(xmin, x);
        xmax = std::max(xmax, x);
    }
}

TGraph *CalculateN(TGraphErrors *g1, TGraphErrors *g2, double sqrts1, double sqrts2)
{
    double xmin1, xmax1;
    double xmin2, xmax2;

    GetGraphRange(g1, xmin1, xmax1);
    GetGraphRange(g2, xmin2, xmax2);

    double xmin = std::max(xmin1, xmin2);
    double xmax = std::min(xmax1, xmax2);

    TSpline3 spline("spline", g2);

    TGraph *gN = new TGraph();

    for (int i = 0; i < g1->GetN(); i++)
    {
        double xT, y1;

        g1->GetPoint(i, xT, y1);

        if (xT < xmin || xT > xmax)
            continue;

        double y2 = spline.Eval(xT);

        if (y1 <= 0 || y2 <= 0)
            continue;

        // cout << std::setprecision(6)
        //      << "xT = " << xT
        //      << "  y136 = " << y1
        //      << "  yOther = " << y2
        //      << "  ratio = " << y1 / y2
        //      << endl;

        double n = -log(y1 / y2) / log(sqrts1 / sqrts2);
        // Formula given in https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.105.062002 (formula 3)

        gN->SetPoint(gN->GetN(), xT, n);
    }
    // cout<<endl;
    SetGraphStyle(gN);
    return gN;
}