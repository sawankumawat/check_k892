#include <iostream>
#include <iomanip>
#include "../src/style.h"
using namespace std;

TFile *OpenFile(const string &path);
TH1D *GetHisto(TFile *f, const string &name);
TGraphErrors *GetGraph(TFile *f, const string &name);
void canvas_style(TCanvas *c, double pad1Size = 0.7, double pad2Size = 0.3, double leftMargin = 0.16, double rightMargin = 0.06, double topMargin = 0.02, double bottomMargin = 0.35);

TGraph *MakeGraph(TH1 *h);

TH1D *BuildRatioHistogram(TH1D *hNumerator, TH1D *hDenominator, const TString &name);

void compareSpectraHyperloop2()
{
    TFile *fSpectras[11];
    int multClasses[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
    TH1D *hSpectraKstar[11];
    TH1D *hSpectraKstarSys[11];
    TH1D *hSpectraPhi[11];
    TH1D *hSpectraPhiSys[11];

    TFile *fSpectraPhi = OpenFile("PhiResults/PhiAnalysisResults_pprun3.root");
    TFile *fSpectraKstar = OpenFile("/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED/KstarSpectra.root");

    string folderPhi[] = {"MB", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10"};

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

        hSpectraKstar[i] = GetHisto(fSpectraKstar, Form("SpectraStat_%.0f_%.0f", multLow, multHigh));
        hSpectraKstarSys[i] = GetHisto(fSpectraKstar, Form("SpectraSys_%.0f_%.0f", multLow, multHigh));

        hSpectraPhi[i] = GetHisto(fSpectraPhi, Form("Spectra/%s/hStat", folderPhi[i].c_str()));
        hSpectraPhiSys[i] = GetHisto(fSpectraPhi, Form("Spectra/%s/hSys", folderPhi[i].c_str()));
    }

    //===============================================
    // ===========Pythia and EPOS from hyperloop================
    //===============================================
    TFile *ModelsHyperloop = OpenFile("ModelRootFiles/ModelResults2_fixedBins.root");

    enum Model
    {
        kEPOS_Hydro,
        kPythiaCR,
        kPythiaMonash,
        kPythiaRopes,
        kPythiaShoving,
        kPythiaMonashRescattering,
        kNModels
    };

    const char *modelLabel[kNModels] = {
        "EPOS Hydro",
        "Pythia CR",
        "Pythia Monash",
        "Pythia Ropes",
        "Pythia Shoving",
        "Pythia Monash Rescattering"};

    string hyperloopModels[kNModels] = {"EPOS_Hydro", "Pythia_CR", "Pythia_Monash", "Pythia_Ropes", "Pythia_Shoving", "Pythia_Monash_Rescattering"};

    TGraphErrors *gYield[kNModels], *gYieldPhi[kNModels];
    double dnch_detaRun3[] = {21.78, 18.48, 15.76, 13.89, 12.50, 10.86, 9.09, 7.63, 5.87, 3.69};

    vector<vector<double>> dNdEtaPythia(kNModels, vector<double>(11, 0.0));

    int colorsPythia[kNModels + 100] = {kBrown, kAzure + 7, kMagenta, kGreen + 2, kBlue + 1, kOrange + 2, kGray + 2};

    TH1D *hSpecModel[11][kNModels];
    vector<double> dNdEtaEPOS(11, 0.0);

    for (int ialice = 0; ialice < 11; ialice++)
    {
        for (int imodel = 0; imodel < kNModels; imodel++)
        {
            gYield[imodel] = GetGraph(ModelsHyperloop, hyperloopModels[imodel] + "/Kstar/gMeanYield_Kstar");

            int totalPoints = gYield[imodel]->GetN();
            if (ialice == 0)
            {
                cout << "Model: " << modelLabel[imodel] << ", Total points: " << totalPoints << endl;

                for (int ipoint = 0; ipoint < totalPoints; ipoint++)
                {
                    double x, y;
                    gYield[imodel]->GetPoint(ipoint, x, y);
                    cout << "Point " << ipoint << ": Multiplicity " << x << ", Multiplicity data: " << dnch_detaRun3[ipoint] << endl;
                }
                cout << endl;
            }

            // Average K*892 and anti-K*892
            string histName = (ialice == 0) ? "MinBias" : Form("%s", to_string(ialice - 1).c_str());
            hSpecModel[ialice][imodel] = GetHisto(ModelsHyperloop, Form("%s/Kstar/hPt_Kstar_%s", hyperloopModels[imodel].c_str(), histName.c_str()));
            hSpecModel[ialice][imodel]->Scale(0.5);
            hSpecModel[ialice][imodel]->SetLineColor(colorsPythia[imodel]);
            hSpecModel[ialice][imodel]->SetLineStyle(2);
            hSpecModel[ialice][imodel]->SetLineWidth(3);
        }
    }

    //============================================================
    //=================Resutls========================
    //============================================================

    for (int iMult = 0; iMult < 11; iMult++)
    {

        double multlow, multhigh;
        if (iMult == 0)
        {
            multlow = 0;
            multhigh = 100;
        }
        else
        {
            multlow = multClasses[iMult - 1];
            multhigh = multClasses[iMult];
        }

        TCanvas *cPythiaCentral = new TCanvas(Form("cPythiaCentral_%.0f_%.0f", multlow, multhigh), Form("cPythiaCentral_%.0f_%.0f", multlow, multhigh), 720, 720);
        double pad1Size = 0.7, pad2Size = 0.3;
        canvas_style(cPythiaCentral, pad1Size, pad2Size);
        SetCanvasStyle(cPythiaCentral, 0.15, 0.03, 0.05, 0.13);
        cPythiaCentral->cd(1);

        gPad->SetLogy();
        hSpectraKstar[iMult]->SetMaximum(hSpectraKstar[iMult]->GetMaximum() * 15);
        hSpectraKstar[iMult]->SetMinimum(hSpectraKstar[iMult]->GetMinimum() * 0.2);
        hSpectraKstar[iMult]->GetXaxis()->SetTitleSize(0.04 / pad1Size);
        hSpectraKstar[iMult]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
        hSpectraKstar[iMult]->GetXaxis()->SetLabelSize(0.04 / pad1Size);
        hSpectraKstar[iMult]->GetYaxis()->SetLabelSize(0.04 / pad1Size);
        hSpectraKstar[iMult]->GetYaxis()->SetTitleOffset(1.8 * pad1Size);
        hSpectraKstar[iMult]->Draw("PE");
        hSpectraKstarSys[iMult]->SetFillStyle(0);
        hSpectraKstarSys[iMult]->Draw("E2 SAME");

        MakeGraph(hSpecModel[iMult][kPythiaMonash])->Draw("L SAME");
        hSpecModel[iMult][kPythiaRopes]->SetLineStyle(2);
        MakeGraph(hSpecModel[iMult][kPythiaRopes])->Draw("L SAME");
        hSpecModel[iMult][kPythiaShoving]->SetLineStyle(2);
        MakeGraph(hSpecModel[iMult][kPythiaShoving])->Draw("L SAME");
        // MakeGraph(hSpecModel[iMult][kEPOS_Hydro])->Draw("L SAME");
        // MakeGraph(hSpecModel[iMult][kPythiaMonashRescattering])->Draw("L SAME");

        TLegend *legend = new TLegend(0.48, 0.6, 0.93, 0.92);
        SetLegendStyle(legend);
        legend->SetTextSize(0.04);
        legend->AddEntry(hSpectraKstar[iMult], Form("pp, %.0f-%.0f%%", multlow, multhigh), "pe");
        legend->AddEntry(hSpecModel[iMult][kPythiaMonash], modelLabel[kPythiaMonash], "l");
        legend->AddEntry(hSpecModel[iMult][kPythiaShoving], modelLabel[kPythiaShoving], "l");
        legend->AddEntry(hSpecModel[iMult][kPythiaRopes], modelLabel[kPythiaRopes], "l");
        // legend->AddEntry(hSpecModel[iMult][kPythiaMonashRes], modelLabel[kPythiaMonashRes], "l");
        // legend->AddEntry(hSpecModel[iMult -1][kPythiaRes], modelLabel[kPythiaRes], "l");
        // legend->AddEntry(hSpecModel[iMult][kEPOS_Hydro], "EPOS", "l");
        // legend->AddEntry(hSpecModel[iMult][kPythiaMonashRescattering], modelLabel[kPythiaMonashRescattering], "l");
        legend->Draw();

        //============================================
        //=======Bottom panel for ratio plots========
        //============================================
        cPythiaCentral->cd(2);
        // gPad->SetLogy();

        TH1D *hRatioMonash = BuildRatioHistogram(hSpectraKstar[iMult], hSpecModel[iMult][kPythiaMonash], Form("hRatioMonash_%.0f_%.0f", multlow, multhigh));

        TH1D *hRatioMonashRes = BuildRatioHistogram(hSpectraKstar[iMult], hSpecModel[iMult][kPythiaMonashRescattering], Form("hRatioMonashRes_%.0f_%.0f", multlow, multhigh));

        TH1D *hRatioShoving = BuildRatioHistogram(hSpectraKstar[iMult], hSpecModel[iMult][kPythiaShoving], Form("hRatioShoving_%.0f_%.0f", multlow, multhigh));

        TH1D *hRatioRopes = BuildRatioHistogram(hSpectraKstar[iMult], hSpecModel[iMult][kPythiaRopes], Form("hRatioRopes_%.0f_%.0f", multlow, multhigh));

        TH1D *hRatioEPOS = BuildRatioHistogram(hSpectraKstar[iMult], hSpecModel[iMult][kEPOS_Hydro], Form("hRatioEPOS_%.0f_%.0f", multlow, multhigh));

        TH1D *hRatioMonashRescattering = BuildRatioHistogram(hSpectraKstar[iMult], hSpecModel[iMult][kPythiaMonashRescattering], Form("hRatioMonashRescattering_%.0f_%.0f", multlow, multhigh));

        // --- NEW: Data Systematic Uncertainty Band Around 1.0 ---
        TH1D *hRatioSysData = (TH1D *)hSpectraKstarSys[iMult]->Clone(Form("hRatioSysData_%.0f_%.0f", multlow, multhigh));
        hRatioSysData->Divide(hSpectraKstar[iMult]); // Sets bin contents to 1.0, scales errors relatively
        hRatioSysData->SetFillColor(kGray);          // Color of the shaded band
        hRatioSysData->SetFillStyle(3001);           // Semi-transparent/hashed fill style
        hRatioSysData->SetMarkerStyle(0);            // Remove markers
        hRatioSysData->GetXaxis()->SetTitleSize(0.04 / pad2Size);
        hRatioSysData->GetYaxis()->SetTitleSize(0.036 / pad2Size);
        hRatioSysData->GetXaxis()->SetLabelSize(0.04 / pad2Size);
        hRatioSysData->GetYaxis()->SetLabelSize(0.04 / pad2Size);
        hRatioSysData->GetYaxis()->SetTitleOffset(2.0 * pad2Size);
        hRatioSysData->GetYaxis()->SetTitle("Data / Model");
        hRatioSysData->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hRatioSysData->SetMaximum(2.05);
        // hRatioSysData->SetMaximum(9.0);
        hRatioSysData->SetMinimum(0.3);
        // hRatioSysData->GetYaxis()->SetMoreLogLabels(kTRUE);
        hRatioSysData->GetYaxis()->SetNdivisions(505);
        hRatioSysData->Draw("E2"); // Draw as an error band before model curves

        TLine *line = new TLine(0.0, 1.0, 20.0, 1.0);
        line->SetLineColor(kBlack);
        line->SetLineStyle(8);
        line->Draw("SAME");

        hRatioMonash->Draw("HIST SAME");
        hRatioShoving->Draw("HIST SAME");
        // hRatioMonashRes->Draw("HIST SAME");
        hRatioRopes->Draw("HIST SAME");
        // hRatioEPOS->GetXaxis()->SetRangeUser(0.0, 6.0);
        // hRatioEPOS->Draw("HIST SAME");
        // RestrictHistogramX(hRatioEPOS, 0.0, 6.0, Form("hRatioEPOS_%.0f_%.0f", multlow, multhigh))->Draw("HIST SAME");
        // hRatioMonashRescattering->Draw("HIST SAME");

        cPythiaCentral->SaveAs(Form("Plots/SpectraCompare/SpectraCompare_%.0f_%.0f.png", multlow, multhigh));
    }
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

TH1D *BuildRatioHistogram(TH1D *hNumerator, TH1D *hDenominator, const TString &name)
{
    TH1D *hRatio = (TH1D *)hNumerator->Clone(name);
    hRatio->Reset("ICES");
    hRatio->Divide(hNumerator, hDenominator, 1.0, 1.0, "B");

    hRatio->SetLineColor(hDenominator->GetLineColor());
    hRatio->SetLineStyle(hDenominator->GetLineStyle());
    hRatio->SetLineWidth(hDenominator->GetLineWidth());

    return hRatio;
}

TGraph *MakeGraph(TH1 *h)
{
    int n = h->GetNbinsX();
    TGraph *gr = new TGraph(n);

    for (int i = 1; i <= n; i++)
    {
        gr->SetPoint(i - 1, h->GetBinCenter(i), h->GetBinContent(i));
    }

    gr->SetLineColor(h->GetLineColor());
    gr->SetLineStyle(h->GetLineStyle());
    gr->SetLineWidth(h->GetLineWidth());

    return gr;
}
