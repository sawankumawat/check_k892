#include <iostream>
#include "src/style.h"
#include "src/fitfunc.h"
#include "src/initializations.h"

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size);
void calculateEfficiency(TFile *fileEff, TH1F *hyieldIntegral, TH1F *heff, const string &MCpath, int multlow1, int multhigh1);
int colors[] = {kBlue + 2, kRed + 1, kGreen + 2, kMagenta + 2, kCyan + 2, kOrange + 7, kViolet + 3, kPink + 1, kAzure + 7, kTeal + 7};

void mc_closure()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    bool isMinBias = true;
    bool isINEL = false;
    bool compareAfterEfficiencyCorrection = false;

    // string MCfile = "657468"; // 2024 (Old with mother id checked in kstarqa code)
    // string MCfile = "665348"; // 2024 (Mother id check commented in the kstarqa code)
    // string MCfile = "666966"; // pp reference MC
    // string MCfile = "667890"; // 2024 (MC_closure, MC_closure_INEL, MC_closure_OnlyTPC: all with TOF shift)
    // string MCfile = "669655"; // 2024 (MC_closure)
    // string MCfile = "673285"; // 2024 (TOF3: MC_closure, MC_closure_INEL, MC_closure_MID0p3)
    // string MCfile = "674418"; // 2024 (TOF3 with checks on Mother: MC_closure, MC_closure_INEL, MC_closure_MID0p3, MC_closure_MID, MC_closure_NoITSROF, MC_closure_WithoutTOFShift)
    string MCfile = "677471"; // 2024 (MC_closure, MC_closure_INEL, MC_closure_MID0p3, MC_closure_MID, MC_closure_NoITSROF, MC_closure_WithoutTOFShift, MC_closure_OnlyTPC, MC_closure_PVContributor)

    string MC_path = "MC_closure_PVContributor";

    string path1 = "../output/kstar/LHC22o_pass7/MC_closure/" + MCfile + "/kstarqa_" + MC_path + "/hInvMass"; // path for yield.root file (from rec MC)
    string path2 = "../data/kstar/LHC22o_pass7/MC_closure/";                                                  // MC file path

    TString savePath = path1 + "/MC_closure_plots";
    gSystem->mkdir(savePath, kTRUE);

    TFile *fspectra1 = new TFile((path1 + ((isINEL) ? "/yield_INEL.root" : "/yield.root")).c_str(), "read");
    TFile *fspectra2 = new TFile((path2 + MCfile + ".root").c_str(), "read");
    if (fspectra1->IsZombie() || fspectra2->IsZombie())
    {
        cout << "Error: files not found" << endl;
        return;
    }

    float mult_classes[] = {0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
    int numofmultbins = sizeof(mult_classes) / sizeof(mult_classes[0]) - 1;
    if (isMinBias || isINEL)
    {
        numofmultbins = 0; // Only one bin for MB or INEL
    }

    TH1F *hmult1[numofmultbins + 1];
    TH1F *hmult2[numofmultbins + 1];
    TH1F *heff[numofmultbins + 1];

    THnSparseF *hSparseRec;
    hSparseRec = (compareAfterEfficiencyCorrection) ? (THnSparseF *)fspectra2->Get(Form("kstarqa_%s/hInvMass/hk892GenpT", MC_path.c_str())) : (THnSparseF *)fspectra2->Get(Form("kstarqa_%s/hInvMass/h2KstarRecpt2", MC_path.c_str()));
    if (hSparseRec == nullptr)
    {
        cout << "Error reading efficiency histogram MC in the path "<< Form("kstarqa_%s/hInvMass/h2KstarRecpt2", MC_path.c_str()) << endl;
        return;
    }
    TH1D *h1recmult = (TH1D *)fspectra2->Get(Form("kstarqa_%s/hInvMass/h1RecMult", MC_path.c_str()));
    double multhigh, multlow;
    TCanvas *cEfficiency = new TCanvas("cEfficiency", "cEfficiency", 720, 720);
    SetCanvasStyle(cEfficiency, 0.15, 0.03, 0.03, 0.15);
    TFile *fOutput = new TFile((savePath + "/MCClosure.root"), "RECREATE");

    for (int imult = 0; imult < numofmultbins + 1; imult++)
    {
        fOutput->cd();
        TDirectory *dir = fOutput->mkdir(Form("mult_%.0f-%.0f", (imult == 0) ? 0 : mult_classes[imult - 1], (imult == 0) ? ((isINEL) ? 120 : 100) : mult_classes[imult]));
        dir->cd();
        if (imult == 0)
        {
            multlow = 0;
            multhigh = (isINEL) ? 120 : 100;
        }
        else
        {
            multlow = mult_classes[imult - 1];
            multhigh = mult_classes[imult];
        }

        hmult1[imult] = (TH1F *)fspectra1->Get(Form("mult_%.0f-%.0f/yield_bincount", multlow, multhigh));
        // hmult1[imult] = (TH1F *)fspectra1->Get(Form("mult_%.0f-%.0f/yield_integral", multlow, multhigh));
        if (hmult1[imult] == nullptr)
        {
            cout << "Histogram hmult1 not found" << endl;
            return;
        }
        hmult2[imult] = (TH1F *)hmult1[imult]->Clone();
        heff[imult] = (TH1F *)hmult1[imult]->Clone(Form("heff_mult_%.0f-%.0f", multlow, multhigh));

        if (compareAfterEfficiencyCorrection)
            calculateEfficiency(fspectra2, hmult1[imult], heff[imult], MC_path, multlow, multhigh);

        int lowbinMultRec = hSparseRec->GetAxis(1)->FindBin(multlow + 1e-5);
        int highbinMultRec = hSparseRec->GetAxis(1)->FindBin(multhigh - 1e-5);
        hSparseRec->GetAxis(1)->SetRange(lowbinMultRec, highbinMultRec);

        TH1D *h1rec_Original = hSparseRec->Projection(0, "E");
        TH1F *h1rec = (TH1F *)h1rec_Original->Rebin(Npt, Form("h1rec_rebinned_mult_%.0f-%.0f", multlow, multhigh), pT_bins);
        // cout << "Number of bins in h1rec after rebinning is " << h1rec->GetNbinsX() << endl;

        int entries = h1recmult->Integral(h1recmult->GetXaxis()->FindBin(multlow + 1e-5), h1recmult->GetXaxis()->FindBin(multhigh - 1e-5));
        cout << "multiplicity class " << multlow << " - " << multhigh << " : " << entries << endl;

        for (int i = 0; i < hmult1[imult]->GetNbinsX(); i++)
        {
            //// Both methods are same if we first rebin the h1rec or we integrate to take the yield value.
            // float highpt = pT_bins[i + 1];
            // float lowpt = pT_bins[i];
            // float ptbinwidth = highpt - lowpt;
            // // cout << "pT bin " << lowpt << " - " << highpt << endl;

            float BR = 0.67;
            float ptbinwidth = h1rec->GetBinWidth(i + 1);
            double yield = h1rec->GetBinContent(i + 1);
            double err = h1rec->GetBinError(i + 1) / (ptbinwidth * BR * entries);
            double nrec = yield / (ptbinwidth * BR * entries);
            cout << "ptbinwidth " << ptbinwidth << " yield " << yield << " err " << err << endl;

            // double nrec = h1rec->Integral(h1rec->GetXaxis()->FindBin(lowpt + 1e-5), h1rec->GetXaxis()->FindBin(highpt - 1e-5)) / (ptbinwidth * BR * entries);
            hmult2[imult]->SetBinContent(i + 1, nrec);
            hmult2[imult]->SetBinError(i + 1, err);
        }

        TH1F *hratio1 = (TH1F *)hmult1[imult]->Clone(Form("ratio_mult_%.0f-%.0f", multlow, multhigh));
        hratio1->Divide(hmult2[imult]);

        TCanvas *c1 = new TCanvas("c1", "c1", 720, 720);
        SetCanvasStyle(c1, 0.25, 0.03, 0.03, 0.15);
        double pad1Size, pad2Size;
        canvas_style(c1, pad1Size, pad2Size);
        c1->cd(1);
        gPad->SetLogy(1);
        SetHistoStyle(hmult1[imult], 1, 53, 1, 0.05, 0.05, 0.04 / pad1Size, 0.04 / pad1Size, 1.13, 1.8);
        hmult1[imult]->GetYaxis()->SetTitleSize(0.04 / pad1Size);
        hmult1[imult]->SetMaximum(hmult1[imult]->GetMaximum() * 10);
        hmult1[imult]->SetMinimum(hmult1[imult]->GetMinimum() * 0.8);
        hmult1[imult]->GetYaxis()->SetTitleOffset(1.30);
        hmult1[imult]->GetXaxis()->SetTitleOffset(1.02);
        hmult1[imult]->SetMarkerStyle(20);
        hmult1[imult]->SetMarkerSize(1);
        // hmult1[imult]->GetXaxis()->SetRangeUser(0, 15);
        hmult1[imult]->Draw("pe");
        hmult2[imult]->SetMarkerStyle(21);
        hmult2[imult]->SetMarkerSize(1);
        hmult2[imult]->SetMarkerColor(kBlue);
        hmult2[imult]->SetLineColor(kBlue);
        hmult2[imult]->SetLineWidth(2);
        hmult2[imult]->Draw("pe same");
        hmult1[imult]->Write("FitYield");
        hmult2[imult]->Write("ReconstructedSpectra");

        TLegend *leg = new TLegend(0.46, 0.64, 0.9, 0.91);
        SetLegendStyle(leg);
        leg->SetHeader(Form("Multiplicity: %.0f-%.0f%%", multlow, multhigh));
        leg->AddEntry(hmult2[imult], "MC Generated", "lpe");
        leg->AddEntry(hmult1[imult], "MC Rec. (AxE corrected)", "lpe");
        leg->SetTextSize(0.04);
        leg->Draw();

        c1->cd(2);
        gPad->SetGridy(1);
        TH1F *hdummy = (TH1F *)hmult1[0]->Clone();
        for (int i = 0; i < hdummy->GetNbinsX(); i++)
        {
            hdummy->SetBinContent(i + 1, 0);
            hdummy->SetBinError(i + 1, 0);
        }

        SetHistoQA(hratio1);
        hratio1->GetYaxis()->SetTitleSize(0.035 / pad2Size);
        hratio1->GetXaxis()->SetTitleSize(0.04 / pad2Size);
        hratio1->GetYaxis()->SetLabelSize(0.04 / pad2Size);
        hratio1->GetXaxis()->SetLabelSize(0.04 / pad2Size);
        hratio1->SetMarkerStyle(20);
        hratio1->SetMarkerSize(1.1);
        // hratio1->SetMarkerColor(kBlue);
        // hratio1->SetLineColor(kBlue);
        hratio1->GetYaxis()->SetTitle("#frac{MC Rec.}{MC Gen.}");
        hratio1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        hratio1->GetXaxis()->CenterTitle(1);
        hratio1->GetYaxis()->SetTitleOffset(0.6);
        hratio1->GetXaxis()->SetTitleOffset(1.1);
        hratio1->GetYaxis()->SetNdivisions(505);
        hratio1->SetMaximum(1.23);
        hratio1->SetMinimum(0.75);
        // hratio1->GetXaxis()->SetRangeUser(0, 15);
        hratio1->Draw("p");
        hratio1->Write("Ratio");

        TLine *line = new TLine(0, 1, 10, 1);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(1);
        line->Draw();
        c1->SaveAs(savePath + Form("/MCclosure_%.0f-%.0f.pdf", multlow, multhigh));

        cEfficiency->cd();
        heff[imult]->SetMarkerStyle(20);
        heff[imult]->SetMarkerSize(1.2);
        heff[imult]->SetMarkerColor(colors[imult]);
        heff[imult]->SetLineColor(colors[imult]);
        heff[imult]->SetLineWidth(2);
        heff[imult]->Draw("pe same");
    }

    cout << "\n\nSelections: " << endl;
    cout << "Is minimum bias: " << ((isMinBias) ? "✅" : "❌") << endl;
    cout << "Is INEL: " << (isINEL ? "✅" : "❌") << endl;
    cout << "Comparing after efficiency correction: " << (compareAfterEfficiencyCorrection ? "✅" : "❌") << endl;
}

void calculateEfficiency(TFile *fileEff, TH1F *hyieldIntegral, TH1F *heff, const string &MCpath, int multlow1, int multhigh1)
{

    THnSparseF *hSpraseGen = (THnSparseF *)fileEff->Get(Form("kstarqa_%s/hInvMass/hk892GenpTCalib1", MCpath.c_str()));
    THnSparseF *hSparseRec = (THnSparseF *)fileEff->Get(Form("kstarqa_%s/hInvMass/h2KstarRecptCalib1", MCpath.c_str()));
    if (hSpraseGen == nullptr || hSparseRec == nullptr)
    {
        cout << "Error reading efficiency histograms " << Form("kstarqa_%s/hInvMass/hk892GenpTCalib1", MCpath.c_str()) << endl;
        return;
    }
    TH1D *h1gen;
    TH1D *h1rec;

    int lowbinMultGen = hSpraseGen->GetAxis(1)->FindBin(multlow1 + 1e-5);
    int highbinMultGen = hSpraseGen->GetAxis(1)->FindBin(multhigh1 - 1e-5);
    hSpraseGen->GetAxis(1)->SetRange(lowbinMultGen, highbinMultGen);

    int lowbinMultRec = hSparseRec->GetAxis(1)->FindBin(multlow1 + 1e-5);
    int highbinMultRec = hSparseRec->GetAxis(1)->FindBin(multhigh1 - 1e-5);
    hSparseRec->GetAxis(1)->SetRange(lowbinMultRec, highbinMultRec);

    h1gen = hSpraseGen->Projection(0, "E");
    h1rec = hSparseRec->Projection(0, "E");

    for (int i = 0; i < heff->GetNbinsX(); i++)
    {
        double lowpt = pT_bins[i];
        double highpt = pT_bins[i + 1];

        double nrec = h1rec->Integral(h1rec->GetXaxis()->FindBin(lowpt), h1rec->GetXaxis()->FindBin(highpt));
        double ngen = h1gen->Integral(h1gen->GetXaxis()->FindBin(lowpt), h1gen->GetXaxis()->FindBin(highpt));
        double efficiency = nrec / ngen;
        double efficiencyerr = sqrt(abs(((nrec + 1) / (ngen + 2)) * ((nrec + 2) / (ngen + 3) - (nrec + 1) / (ngen + 2))));

        cout << "Efficiency: " << efficiency << " +/- " << efficiencyerr << endl;

        heff->SetBinContent(i + 1, efficiency);
        heff->SetBinError(i + 1, efficiencyerr);
        hyieldIntegral->SetBinContent(i + 1, hyieldIntegral->GetBinContent(i + 1) / (efficiency));

        double errorinyieldIntegral = hyieldIntegral->GetBinError(i + 1);
        double rawyieldvalueIntegral = hyieldIntegral->GetBinContent(i + 1);
        hyieldIntegral->SetBinError(i + 1, sqrt(pow(errorinyieldIntegral / efficiency, 2) + pow(rawyieldvalueIntegral * efficiencyerr / (efficiency * efficiency), 2)));
    }
}

void canvas_style(TCanvas *c, double &pad1Size, double &pad2Size)
{
    SetCanvasStyle(c, 0.15, 0.005, 0.05, 0.15);
    c->Divide(1, 2, 0, 0);
    TPad *pad1 = (TPad *)c->GetPad(1);
    TPad *pad2 = (TPad *)c->GetPad(2);
    pad2Size = 0.3; // Size of the first pad
    pad1Size = 1 - pad2Size;

    pad1->SetPad(0, 0.3, 1, 1); // x1, y1, x2, y2 (top pad)
    pad2->SetPad(0, 0, 1, 0.3);
    pad1->SetRightMargin(0.06);
    pad2->SetRightMargin(0.06);
    pad2->SetBottomMargin(0.33);
    pad1->SetLeftMargin(0.16);
    pad2->SetLeftMargin(0.16);
    pad1->SetTopMargin(0.02);
    pad1->SetBottomMargin(0.001);
    pad2->SetTopMargin(0.001);

    pad1->SetTicks(1, 1);
    pad2->SetTicks(1, 1);
}