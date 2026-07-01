#include <iostream>
#include <iomanip>
#include "../src/style.h"
using namespace std;

TFile *OpenFile(const string &path);
TH1D *GetHisto(TFile *f, const string &name);
TGraphErrors *GetGraph(TFile *f, const string &name);
TH1D *RebinToMatch(const TH1D *hInput, const TH1D *hReference, const TString &newName = "");
void canvas_style(TCanvas *c, double pad1Size = 0.7, double pad2Size = 0.3, double leftMargin = 0.16, double rightMargin = 0.06, double topMargin = 0.02, double bottomMargin = 0.35);

struct ClosestMultiplicity
{
    bool useAverage = false;
    int idx1 = -1;      // first point
    int idx2 = -1;      // second point (-1 if single point)
    double value = 0.0; // selected multiplicity (or average)
};
ClosestMultiplicity FindClosestMultiplicity(TGraphErrors *g, double target);
TH1D *BuildSpectrum(TFile *file, const ClosestMultiplicity &match, const char *histPrefix = "hkstarPt_");
TH1D *BuildRatioHistogram(TH1D *hNumerator, TH1D *hDenominator, const TString &name);

void compareSpectra()
{
    TFile *fSpectras[11];
    int multClasses[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    TH1D *hSpectra[11];
    TH1D *hSpectraSys[11];
    string sysPath = "../../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/SystematicsPlots/SysUncert.root";
    TFile *fSystematics = OpenFile(sysPath);
    TH1D *hRelUncert = GetHisto(fSystematics, "hTotalSysSmoothed_0_100");

    for (int i = 0; i < 11; i++)
    {
        double multLow, multHigh;
        if (i == 0)
        {
            multLow = 0;
            multHigh = 100;
        }
        else
        {
            multLow = multClasses[i - 1];
            multHigh = multClasses[i];
        }
        fSpectras[i] = OpenFile(Form("../../output/kstar/LHC22o_pass7/679906/kstarqa/hInvMass/corrected_spectra_%.0f_%.0f.root", multLow, multHigh));
        hSpectra[i] = GetHisto(fSpectras[i], Form("mult_%.0f-%.0f/corrected_spectra_Integral_final", multLow, multHigh));
        hSpectraSys[i] = (TH1D *)hSpectra[i]->Clone(Form("hSpectraSys_%.0f-%.0f", multLow, multHigh));
        for (int bin = 1; bin <= hSpectraSys[i]->GetNbinsX(); ++bin)
        {
            double content = hSpectraSys[i]->GetBinContent(bin);
            double relUncert = hRelUncert->GetBinContent(bin);
            double sysError = content * relUncert;
            hSpectraSys[i]->SetBinError(bin, sysError);
        }
    }

    //===============================================
    // ============EPOS tested locally================
    //===============================================
    TFile *fEPOS = OpenFile("EPOS_finalQA_CorrectpTCutPiKp.root");

    //===============================================
    // ===========Pythia tested locally================
    //===============================================
    TFile *fPythiaMonash = OpenFile("../../pythia/Pythia_MonashLocal.root");
    TFile *fPythiaMonashNoCR = OpenFile("../../pythia/Pythia_MonashWoCRLocal.root");
    TFile *fPythiaShoving = OpenFile("../../pythia/Pythia_ShovingLocal.root");
    TFile *fPythiaMonashRescattering = OpenFile("../../pythia/Pythia_MonashRescatteringLocal.root");
    TFile *fPythiaRopes = OpenFile("../../pythia/Pythia_RopesLocal.root");
    TFile *fPythiaRescattering = OpenFile("../../pythia/Pythia_RescatteringLocal.root");

    enum PythiaModel
    {
        kPythiaMonash,
        kPythiaMonashNoCR,
        kPythiaShoving,
        kPythiaRopes,
        kPythiaMonashRes,
        kPythiaRes,
        kNPythiaModels
    };

    const char *modelLabelLocal[kNPythiaModels] = {
        "Pythia Monash",
        "Pythia Monash No CR",
        "Pythia Shoving",
        "Pythia Ropes",
        "Pythia Monash Rescattering",
        "Pythia Rescattering"};

    vector<TFile *> fPythiaModels = {fPythiaMonash, fPythiaMonashNoCR, fPythiaShoving, fPythiaRopes, fPythiaMonashRescattering, fPythiaRescattering};
    TGraphErrors *gPythiaYieldLocal[kNPythiaModels];
    int centrality[11] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    double dnch_detaRun3[] = {21.78, 18.48, 15.76, 13.89, 12.50, 10.86, 9.09, 7.63, 5.87, 3.69};

    TH1D *hSpectraHighestMult[kNPythiaModels];
    TH1D *hSpectraMult20to30[kNPythiaModels];
    vector<vector<double>> dNdEtaPythia(kNPythiaModels, vector<double>(11, 0.0));

    int colorsPythia[kNPythiaModels + 100] = {kMagenta, kBrown, kGreen + 2, kGray + 2, kBlue + 1, kAzure + 7, kOrange + 2};

    TH1D *hSpectraModel[10][kNPythiaModels];
    TH1D *hSpecRebinned[10][kNPythiaModels];
    TH1D *hSpectraModelEPOS[10];
    TH1D *hSpecRebinnedEPOS[10];
    vector<double> dNdEtaEPOS(10, 0.0);
    TGraphErrors *gEPOSYield = GetGraph(fEPOS, "IST9_ITY80/kstar_vs_mult");

    for (int ialice = 0; ialice < 10; ialice++)
    {
        for (int imodel = 0; imodel < kNPythiaModels; imodel++)
        {
            gPythiaYieldLocal[imodel] = GetGraph(fPythiaModels[imodel], "gMeanYield_kstar");

            int totalPoints = gPythiaYieldLocal[imodel]->GetN();
            cout << "Model: " << modelLabelLocal[imodel] << ", Total points: " << totalPoints << endl;

            for (int ipoint = 0; ipoint < totalPoints; ipoint++)
            {
                double x, y;
                gPythiaYieldLocal[imodel]->GetPoint(ipoint, x, y);
                cout << "Point " << ipoint << ": Multiplicity " << x << ", Yield: " << y << endl;
            }
            cout << endl;
            ClosestMultiplicity match = FindClosestMultiplicity(gPythiaYieldLocal[imodel], dnch_detaRun3[ialice]);
            dNdEtaPythia[imodel][ialice] = match.value;

            if (match.useAverage)
            {
                cout << modelLabelLocal[imodel]
                     << " : Average of points "
                     << match.idx1 << " and "
                     << match.idx2
                     << " (value = " << match.value << ")"
                     << endl;
            }
            else
            {
                cout << modelLabelLocal[imodel]
                     << " : Closest point "
                     << match.idx1
                     << " (value = " << match.value << ")"
                     << endl;
            }
            cout << endl;

            hSpectraModel[ialice][imodel] = BuildSpectrum(fPythiaModels[imodel], match);
            hSpecRebinned[ialice][imodel] = RebinToMatch(hSpectraModel[ialice][imodel], hSpectra[ialice], Form("hSpec_%d_%s", ialice, modelLabelLocal[imodel]));

            // Average K*892 and anti-K*892
            hSpecRebinned[ialice][imodel]->Scale(0.5);
            hSpecRebinned[ialice][imodel]->SetLineColor(colorsPythia[imodel]);
            hSpecRebinned[ialice][imodel]->SetLineStyle(2);
            hSpecRebinned[ialice][imodel]->SetLineWidth(3);
        }

        // Same for EPOS
        int totalPointsEPOS = gEPOSYield->GetN();
        cout << "EPOS: Total points: " << totalPointsEPOS << endl;
        for (int ipoint = 0; ipoint < totalPointsEPOS; ipoint++)
        {
            double x, y;
            gEPOSYield->GetPoint(ipoint, x, y);
            cout << "Point " << ipoint << ": Multiplicity " << x << ", Yield: " << y << endl;
        }
        cout << endl;
        ClosestMultiplicity matchEPOS = FindClosestMultiplicity(gEPOSYield, dnch_detaRun3[ialice]);

        if (matchEPOS.useAverage)
        {
            cout << "EPOS : Average of points "
                 << matchEPOS.idx1 << " and "
                 << matchEPOS.idx2
                 << " (value = " << matchEPOS.value << ")"
                 << endl;
        }
        else
        {
            cout << "EPOS : Closest point "
                 << matchEPOS.idx1
                 << " (value = " << matchEPOS.value << ")"
                 << endl;
        }
        cout << endl;
        dNdEtaEPOS[ialice] = matchEPOS.value;

        hSpectraModelEPOS[ialice] = BuildSpectrum(fEPOS, matchEPOS, "IST9_ITY80/hPtMB_kstar_IST9_ITY80_Cent");
        hSpecRebinnedEPOS[ialice] = RebinToMatch(hSpectraModelEPOS[ialice], hSpectra[ialice], Form("hSpecEPOS_%d", ialice));

        // Average K*892 and anti-K*892
        hSpecRebinnedEPOS[ialice]->Scale(0.5);
        hSpecRebinnedEPOS[ialice]->SetLineColor(kGreen - 2);
        hSpecRebinnedEPOS[ialice]->SetLineStyle(2);
        hSpecRebinnedEPOS[ialice]->SetLineWidth(3);
    }

    // //======================================================
    // //    ===========EPOS local model (min Bias only)===========
    // //======================================================
    // TFile *fEPOS2 = OpenFile("ModelRootFiles/EPOS_finalQA_ptCut_FinerBins.root");
    // TH1D *hEPOS_Yield = GetHisto(fEPOS2, "IST9_ITY80/hPtMB_kstar");
    // TCanvas *cEPOS = new TCanvas("cEPOS", "cEPOS", 720, 720);
    // SetCanvasStyle(cEPOS, 0.15, 0.03, 0.05, 0.15);
    // gPad->SetLogy();
    // hSpectra[0]->Draw("PE");

    // hEPOS_Yield->Scale(1.0 / hEPOS_Yield->Integral());
    // hEPOS_Yield->Scale(0.2); // Average of K* and K*bar

    // TH1D *hEPOS_rebinned = RebinToMatch(hEPOS_Yield, hSpectra[0], "hEPOS_rebinned");
    // hEPOS_rebinned->SetLineColor(kRed + 1);
    // hEPOS_rebinned->SetLineStyle(2);
    // hEPOS_rebinned->Draw("HIST SAME");

    // //Check from hyperloop epos file
    // TFile *fEPOS_hyperloop = OpenFile("ModelRootFiles/EPOS_Hydro.root");
    // TH2D *hMultSpectra = (TH2D *)fEPOS_hyperloop->Get("mc-particle-prediction/prediction/pt/FT0AC/Kstar");
    // TH1D *hSpectraHydro = hMultSpectra->ProjectionX("hSpectraHydro");
    // hSpectraHydro->Scale(1.0 / hSpectraHydro->Integral());
    // hSpectraHydro->Scale(0.18);        // Average of K* and K*bar
    // TH1D *hSpectraHydro_rebinned = RebinToMatch(hSpectraHydro, hSpectra[0], "hSpectraHydro_rebinned");
    // hSpectraHydro_rebinned->SetLineColor(kBlue + 1);
    // hSpectraHydro_rebinned->SetLineStyle(2);
    // hSpectraHydro_rebinned->Draw("HIST SAME");

    //============================================================
    //=================Pythia (Central production)========================
    //============================================================
    TFile *fMC = OpenFile("../..//mc/LHC24f3c/679945.root");
    THnSparse *hGenSparse = (THnSparse *)fMC->Get("kstarqa/hInvMass/hk892GenpT");
    if (!hGenSparse || hGenSparse == nullptr)
    {
        cout << "Error: hk892GenpTCalib1 not found in file " << fMC->GetName() << endl;
        return;
    }

    for (int iMult = 1; iMult < 11; iMult++)
    {

        // double multlow = 0;
        // double multhigh = 1;
        double multlow = multClasses[iMult - 1];
        double multhigh = multClasses[iMult];

        int lowbinMultGen = hGenSparse->GetAxis(1)->FindBin(multlow + 1e-3);
        int highbinMultGen = hGenSparse->GetAxis(1)->FindBin(multhigh - 1e-3);

        hGenSparse->GetAxis(1)->SetRange(lowbinMultGen, highbinMultGen);

        TH1D *hGenMult = GetHisto(fMC, "kstarqa/hInvMass/h1RecMult");
        hGenMult->SetName(Form("hGenMult_%.0f_%.0f", multlow, multhigh));
        int totalEvents = hGenMult->Integral(lowbinMultGen, highbinMultGen);

        TH1D *hGenSpectra = hGenSparse->Projection(0, "E");
        hGenSpectra->SetName(Form("hGenSpectra_%.0f_%.0f", multlow, multhigh));
        hGenSpectra->Scale(1.0 / totalEvents);
        hGenSpectra->Scale(0.5);        // Average of K* and K*bar
        hGenSpectra->Scale(1.0 / 0.66); // Branching ratio correction
        TH1D *hGenSpectraRebinned = RebinToMatch(hGenSpectra, hSpectra[0], "hGenSpectraRebinned");
        hGenSpectraRebinned->SetName(Form("hGenSpectraRebinned_%.0f_%.0f", multlow, multhigh));

        // 0: Min Bias
        // 1: 0-1%
        // 2: 1-5%
        // 3: 5-10%
        // 4: 10-15%
        // 5: 15-20%
        // 6: 20-30%
        // 7: 30-40%
        // 8: 40-50%
        // 9: 50-70%
        // 10: 70-100%

        // int WhichCent = 1; // Do no give 0, as I have not initiated MinBias spectra for pythia. For that I need to do average sum over all points.
        int WhichCent = iMult;

        TCanvas *cPythiaCentral = new TCanvas(Form("cPythiaCentral_%.0f_%.0f", multlow, multhigh), Form("cPythiaCentral_%.0f_%.0f", multlow, multhigh), 720, 720);
        double pad1Size = 0.7, pad2Size = 0.3;
        canvas_style(cPythiaCentral, pad1Size, pad2Size);
        SetCanvasStyle(cPythiaCentral, 0.15, 0.03, 0.05, 0.13);
        cPythiaCentral->cd(1);

        gPad->SetLogy();
        hSpectra[WhichCent]->SetMaximum(hSpectra[WhichCent]->GetMaximum() * 15);
        hSpectra[WhichCent]->SetMinimum(hSpectra[WhichCent]->GetMinimum() * 0.2);
        hSpectra[WhichCent]->GetXaxis()->SetTitleSize(0.04 / pad1Size);
        hSpectra[WhichCent]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
        hSpectra[WhichCent]->GetXaxis()->SetLabelSize(0.04 / pad1Size);
        hSpectra[WhichCent]->GetYaxis()->SetLabelSize(0.04 / pad1Size);
        hSpectra[WhichCent]->GetYaxis()->SetTitleOffset(1.8 * pad1Size);
        hSpectra[WhichCent]->Draw("PE");
        hSpectraSys[WhichCent]->SetFillStyle(0);
        hSpectraSys[WhichCent]->Draw("E2 SAME");
        hGenSpectraRebinned->SetLineColor(kRed + 1);
        hGenSpectraRebinned->SetLineStyle(2);
        hGenSpectraRebinned->SetLineWidth(3);
        hGenSpectraRebinned->Draw("HIST SAME");
        hSpecRebinned[WhichCent - 1][kPythiaMonash]->Draw("HIST SAME");
        // hSpecRebinned[WhichCent - 1][kPythiaShoving]->Draw("HIST SAME");
        hSpecRebinned[WhichCent - 1][kPythiaRopes]->SetLineStyle(1);
        hSpecRebinned[WhichCent - 1][kPythiaRopes]->Draw("HIST SAME");
        // hSpecRebinned[WhichCent - 1][kPythiaMonashRes]->Draw("HIST SAME");
        // hSpecRebinned[WhichCent - 1][kPythiaRes]->Draw("HIST SAME");
        hSpecRebinnedEPOS[WhichCent - 1]->Draw("HIST SAME");

        TLegend *legend = new TLegend(0.48, 0.68, 0.93, 0.92);
        SetLegendStyle(legend);
        legend->SetTextSize(0.03);
        legend->AddEntry(hSpectra[WhichCent], Form("pp, %d-%d%%", centrality[WhichCent - 1], centrality[WhichCent]), "pe");
        legend->AddEntry(hGenSpectraRebinned, "Pythia Monash (Central Prod.)", "l");
        legend->AddEntry(hSpecRebinned[WhichCent - 1][kPythiaMonash], modelLabelLocal[kPythiaMonash], "l");
        // legend->AddEntry(hSpecRebinned[WhichCent - 1][kPythiaShoving], modelLabelLocal[kPythiaShoving], "l");
        legend->AddEntry(hSpecRebinned[WhichCent - 1][kPythiaRopes], modelLabelLocal[kPythiaRopes], "l");
        // legend->AddEntry(hSpecRebinned[WhichCent - 1][kPythiaMonashRes], modelLabelLocal[kPythiaMonashRes], "l");
        // legend->AddEntry(hSpecRebinned[WhichCent -1][kPythiaRes], modelLabelLocal[kPythiaRes], "l");
        legend->AddEntry(hSpecRebinnedEPOS[WhichCent - 1], "EPOS", "l");
        legend->Draw();

        //============================================
        //=======Bottom panel for ratio plots========
        //============================================
        cPythiaCentral->cd(2);
        TH1D *hRatioGen = BuildRatioHistogram(hSpectra[WhichCent], hGenSpectraRebinned, Form("hRatioGen_%.0f_%.0f", multlow, multhigh));

        TH1D *hRatioMonash = BuildRatioHistogram(hSpectra[WhichCent], hSpecRebinned[WhichCent - 1][kPythiaMonash], Form("hRatioMonash_%.0f_%.0f", multlow, multhigh));

        TH1D *hRatioMonashRes = BuildRatioHistogram(hSpectra[WhichCent], hSpecRebinned[WhichCent - 1][kPythiaMonashRes], Form("hRatioMonashRes_%.0f_%.0f", multlow, multhigh));

        TH1D *hRatioShoving = BuildRatioHistogram(hSpectra[WhichCent], hSpecRebinned[WhichCent - 1][kPythiaShoving], Form("hRatioShoving_%.0f_%.0f", multlow, multhigh));

        TH1D *hRatioRopes = BuildRatioHistogram(hSpectra[WhichCent], hSpecRebinned[WhichCent - 1][kPythiaRopes], Form("hRatioRopes_%.0f_%.0f", multlow, multhigh));

        TH1D *hRatioEPOS = BuildRatioHistogram(hSpectra[WhichCent], hSpecRebinnedEPOS[WhichCent - 1], Form("hRatioEPOS_%.0f_%.0f", multlow, multhigh));

        hRatioGen->GetXaxis()->SetTitleSize(0.04 / pad2Size);
        hRatioGen->GetYaxis()->SetTitleSize(0.036 / pad2Size);
        hRatioGen->GetXaxis()->SetLabelSize(0.04 / pad2Size);
        hRatioGen->GetYaxis()->SetLabelSize(0.04 / pad2Size);
        hRatioGen->GetYaxis()->SetTitleOffset(2.0 * pad2Size);
        hRatioGen->GetYaxis()->SetTitle("Data / Model");
        hRatioGen->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hRatioGen->SetMaximum(2.3);
        hRatioGen->SetMinimum(0.0);
        hRatioGen->GetYaxis()->SetNdivisions(505);
        hRatioGen->Draw("HIST"); // First one initializes the frame

        TLine *line = new TLine(hRatioGen->GetXaxis()->GetXmin(), 1.0, hRatioGen->GetXaxis()->GetXmax(), 1.0);
        line->SetLineColor(kBlack);
        line->SetLineStyle(8);
        line->Draw("SAME");

        // --- NEW: Data Systematic Uncertainty Band Around 1.0 ---
        TH1D *hRatioSysData = (TH1D *)hSpectraSys[WhichCent]->Clone(Form("hRatioSysData_%.0f_%.0f", multlow, multhigh));
        hRatioSysData->Divide(hSpectra[WhichCent]); // Sets bin contents to 1.0, scales errors relatively
        hRatioSysData->SetFillColor(kGray);         // Color of the shaded band
        hRatioSysData->SetFillStyle(3001);          // Semi-transparent/hashed fill style
        hRatioSysData->SetMarkerStyle(0);           // Remove markers
        hRatioSysData->Draw("E2 SAME");             // Draw as an error band before model curves

        hRatioMonash->Draw("HIST SAME");
        // hRatioShoving->Draw("HIST SAME");
        // hRatioMonashRes->Draw("HIST SAME");
        hRatioRopes->Draw("HIST SAME");
        hRatioEPOS->Draw("HIST SAME");
        hRatioGen->Draw("HIST SAME"); // Draw last to ensure it's on top

        cPythiaCentral->SaveAs(Form("Plots/SpectraCompare/SpectraCompare_%.0f_%.0f.png", multlow, multhigh));

        hGenSparse->GetAxis(1)->SetRange(0, -1);
    }

    // // Lets write in table form the dNdEta values from data, EPOS and Pythia models
    // cout << "Multiplicity Class\tData dNch/dEta\tEPOS dNch/dEta\tPythia Monash dNch/dEta\tPythia Monash No CR dNch/dEta\tPythia Shoving dNch/dEta\tPythia Ropes dNch/dEta\tPythia Monash Rescattering dNch/dEta\tPythia Rescattering dNch/dEta" << endl;
    // for (int i = 0; i < 10; i++)
    // {
    //     cout << Form("%d-%d\t", centrality[i], centrality[i + 1]) << Form("%.2f\t", dnch_detaRun3[i]) << Form("%.2f\t", dNdEtaEPOS[i]);
    //     for (int imodel = 0; imodel < kNPythiaModels; imodel++)
    //     {
    //         cout << Form("%.2f\t", dNdEtaPythia[imodel][i]);
    //     }
    //     cout << endl;
    // }
}

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

TH1D *RebinToMatch(const TH1D *hInput, const TH1D *hReference, const TString &newName = "")
{
    if (!hInput || !hReference)
        return nullptr;

    TString hName = newName;
    if (hName.IsNull())
        hName = Form("%s_rebinned", hInput->GetName());

    TH1D *hRebinned = (TH1D *)hReference->Clone(hName);
    hRebinned->Reset();
    hRebinned->Sumw2();

    for (int i = 1; i <= hReference->GetNbinsX(); ++i)
    {
        double low = hReference->GetBinLowEdge(i);
        double high = hReference->GetBinLowEdge(i + 1);

        int binLow = hInput->GetXaxis()->FindBin(low + 1e-8);
        int binHigh = hInput->GetXaxis()->FindBin(high - 1e-8);

        double err = 0.0;
        double integral = hInput->IntegralAndError(binLow, binHigh, err);

        double width = high - low;
        // cout << "Bin " << i << ", pT Low " << low << " High " << high << ", width " << width << endl;

        hRebinned->SetBinContent(i, integral / width);
        hRebinned->SetBinError(i, err / width);
    }

    return hRebinned;
}

ClosestMultiplicity FindClosestMultiplicity(TGraphErrors *g, double target)
{
    ClosestMultiplicity result;

    const int n = g->GetN();
    if (n == 0)
        return result;

    // Store all multiplicities
    std::vector<double> mult(n);
    for (int i = 0; i < n; i++)
    {
        double x, y;
        g->GetPoint(i, x, y);
        mult[i] = x;
    }

    // Find closest point
    int closest = 0;
    double bestDiff = std::abs(mult[0] - target);

    for (int i = 1; i < n; i++)
    {
        double diff = std::abs(mult[i] - target);
        if (diff < bestDiff)
        {
            bestDiff = diff;
            closest = i;
        }
    }

    // Default answer = closest point
    result.useAverage = false;
    result.idx1 = closest;
    result.value = mult[closest];

    // Check neighbour averages
    auto CheckAverage = [&](int i1, int i2)
    {
        if (i1 < 0 || i2 < 0 || i1 >= n || i2 >= n)
            return;

        double avg = 0.5 * (mult[i1] + mult[i2]);
        double diff = std::abs(avg - target);

        if (diff < std::abs(result.value - target))
        {
            result.useAverage = true;
            result.idx1 = i1;
            result.idx2 = i2;
            result.value = avg;
        }
    };

    CheckAverage(closest - 1, closest);
    CheckAverage(closest, closest + 1);

    return result;
}

TH1D *BuildSpectrum(TFile *file,
                    const ClosestMultiplicity &match,
                    const char *histPrefix = "hkstarPt_")
{
    if (!file)
        return nullptr;

    // ---------- Single multiplicity ----------
    if (!match.useAverage)
    {
        TH1D *h = GetHisto(file, Form("%s%d", histPrefix, match.idx1));
        h = (TH1D *)h->Clone(Form("%s_clone_%d", histPrefix, match.idx1));
        // h->Scale(1.0 / h->Integral());
        return h;
    }

    // ---------- Average of two multiplicities ----------
    TH1D *h1 = GetHisto(file, Form("%s%d", histPrefix, match.idx1));
    TH1D *h2 = GetHisto(file, Form("%s%d", histPrefix, match.idx2));

    h1 = (TH1D *)h1->Clone(Form("%s_avg_%d_%d",
                                histPrefix,
                                match.idx1,
                                match.idx2));

    h2 = (TH1D *)h2->Clone(Form("%s_tmp_%d", histPrefix, match.idx2));

    // h1->Scale(1.0 / h1->Integral());
    // h2->Scale(1.0 / h2->Integral());

    h1->Add(h2);
    h1->Scale(0.5);

    delete h2;

    return h1;
}

void canvas_style(TCanvas *c, double pad1Size = 0.7, double pad2Size = 0.3, double leftMargin = 0.16, double rightMargin = 0.06, double topMargin = 0.02, double bottomMargin = 0.35)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    // pad2Size = 0.3; // Size of the first pad
    // pad1Size = 1 - pad2Size;

    pad1->SetPad(0, pad2Size, 1, 1); // x1, y1, x2, y2 (top pad)
    pad2->SetPad(0, 0, 1, pad2Size); // x1, y1, x2, y2 (bottom pad)
    pad1->SetRightMargin(rightMargin);
    pad2->SetRightMargin(rightMargin);
    pad1->SetLeftMargin(leftMargin);
    pad2->SetLeftMargin(leftMargin);
    pad1->SetTopMargin(topMargin);
    pad2->SetTopMargin(0.001);
    pad2->SetBottomMargin(bottomMargin);
    pad1->SetBottomMargin(0.001);

    pad1->SetTickx(1);
    pad1->SetTicky(1);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
}

TH1D *BuildRatioHistogram(TH1D *hNumerator,
                          TH1D *hDenominator,
                          const TString &name)
{
    TH1D *hRatio = (TH1D *)hNumerator->Clone(name);
    hRatio->Reset("ICES");
    hRatio->Divide(hNumerator, hDenominator, 1.0, 1.0, "B");

    hRatio->SetLineColor(hDenominator->GetLineColor());
    hRatio->SetLineStyle(hDenominator->GetLineStyle());
    hRatio->SetLineWidth(hDenominator->GetLineWidth());

    return hRatio;
}
