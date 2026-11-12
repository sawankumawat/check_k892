#include <iostream>
#include <iomanip>
#include "../src/style.h"
using namespace std;

//--------Important info---------------------//
//=====In data:
// 1. K*0, K*+- is averaged
// 2. Pi,K,p is averaged
//=====In EPOS:
// 1. Only K*0 is averaged.
//===== In Pythia (Local):
// 1. Nothing is averaged. All as summed.

// Important: Do not scale any graph without cloning it first, otherwise the original graph will be scaled and it will affect other plots.
TGraphErrors *SmoothGraph(const TGraphErrors *gr);
TGraphErrors *MakeRatio(const TGraphErrors *numerator, const TGraphErrors *denominator, bool isModel = false, float CorrFactorDen = 1);
TGraphErrors **MakeRatioUncorr(TGraphErrors **numerator, TGraphErrors **denominator, float CorrFactorDen = 1, int nGraphs = 3);
TGraphErrors *GetGraph(TFile *f, const string &name);
TFile *OpenFile(const string &path);

void setStyle(TGraphErrors *gr, int color, int style)
{
    gr->SetLineColor(color);
    gr->SetLineStyle(style);
}

void ScaleGraph(TGraph *gr, double scale);
TGraphErrors **DivideByMult(TGraphErrors **gr, double WhichMultPoint, double tolerance = 0.5, int noOfGraphs = 3);
TGraphErrors *DivideByMultModel(TGraphErrors *gr, double WhichMultPoint, double tolerance = 0.5);
void RestrictModelXaxis(TGraphErrors *gr, double xMin, double xMax);

void checkEPOS()
{
    // TFile *fEPOS_New = OpenFile("/home/sawan/Storage/EPOS4_Resonance/FinalEPOSResults_BongFix2.root");
    // TFile *fEPOS_New2 = OpenFile("/home/sawan/Storage/EPOS4_Resonance/FinalEPOSResults_BongFix2.root");

    TFile *fEPOS_New = OpenFile("/home/sawan/Storage/EPOS4_Resonance/FinalEPOSResults_BongFix_noIST9.root");
    TFile *fEPOS_New2 = OpenFile("/home/sawan/Storage/EPOS4_Resonance/FinalEPOSResults_BongFix_noIST9.root");

    TGraphErrors *gEPOS_KstardNdy_UrQMDOff = GetGraph(fEPOS_New, "g_dNdy_opt0");
    TGraphErrors *gEPOS_KstardNdy_UrQMDOn = GetGraph(fEPOS_New, "g_dNdy_opt1");
    TGraphErrors *gEPOS_KstarMeanpT_UrQMDOff = GetGraph(fEPOS_New, "g_pT_opt0");
    TGraphErrors *gEPOS_KstarMeanpT_UrQMDOn = GetGraph(fEPOS_New, "g_pT_opt1");

    TGraphErrors *gEPOS_PhidNdy_UrQMDOff = GetGraph(fEPOS_New, "g_phi_dNdy_opt0");
    TGraphErrors *gEPOS_PhidNdy_UrQMDOn = GetGraph(fEPOS_New, "g_phi_dNdy_opt1");
    TGraphErrors *gEPOS_PhiMeanpT_UrQMDOff = GetGraph(fEPOS_New, "g_phi_pT_opt0");
    TGraphErrors *gEPOS_PhiMeanpT_UrQMDOn = GetGraph(fEPOS_New, "g_phi_pT_opt1");

    TCanvas *cdNdyKstarEPOS = new TCanvas("cdNdyKstarEPOS", "cdNdyKstarEPOS", 720, 720);
    cdNdyKstarEPOS->Divide(1, 2, 0.01, 0.01);
    cdNdyKstarEPOS->cd(1);
    gPad->SetTopMargin(0.0);
    gPad->SetBottomMargin(0.0);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.001);
    SetGraphErrorStyle(gEPOS_KstardNdy_UrQMDOff);
    gEPOS_KstardNdy_UrQMDOff->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gEPOS_KstardNdy_UrQMDOff->GetYaxis()->SetTitle("dN/dy");
    gEPOS_KstardNdy_UrQMDOff->SetLineColor(kRed);
    gEPOS_KstardNdy_UrQMDOff->SetLineStyle(2);
    gEPOS_KstardNdy_UrQMDOff->SetLineWidth(3);
    gEPOS_KstardNdy_UrQMDOff->Draw("AL");
    gEPOS_KstardNdy_UrQMDOn->SetLineColor(kBlue - 6);
    gEPOS_KstardNdy_UrQMDOn->SetLineStyle(9);
    gEPOS_KstardNdy_UrQMDOn->SetLineWidth(3);
    gEPOS_KstardNdy_UrQMDOn->Draw("L same");

    TLegend *legEPOS_Kstar = new TLegend(0.22, 0.75, 0.85, 0.85);
    legEPOS_Kstar->SetBorderSize(0);
    legEPOS_Kstar->SetFillStyle(0);
    legEPOS_Kstar->SetTextSize(0.05);
    legEPOS_Kstar->AddEntry(gEPOS_KstardNdy_UrQMDOff, "EPOS UrQMD OFF", "l");
    legEPOS_Kstar->AddEntry(gEPOS_KstardNdy_UrQMDOn, "EPOS UrQMD ON", "l");
    legEPOS_Kstar->Draw();

    cdNdyKstarEPOS->cd(2);
    gPad->SetTopMargin(0.0);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.001);
    int totalPointsEPOS_Kstar = gEPOS_KstardNdy_UrQMDOff->GetN();
    TGraphErrors *gRatioEPOS = new TGraphErrors(totalPointsEPOS_Kstar);
    for (int i = 0; i < totalPointsEPOS_Kstar; i++)
    {
        double x, yOff, yOn;
        gEPOS_KstardNdy_UrQMDOff->GetPoint(i, x, yOff);
        gEPOS_KstardNdy_UrQMDOn->GetPoint(i, x, yOn);
        double ratio = yOn / yOff;

        gRatioEPOS->SetPoint(i, x, ratio);
        gRatioEPOS->SetPointError(i, 0, 0);
    }
    SetGraphErrorStyle(gRatioEPOS);
    gRatioEPOS->SetLineColor(kBlack);
    gRatioEPOS->SetLineStyle(1);
    gRatioEPOS->SetLineWidth(2);
    gRatioEPOS->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gRatioEPOS->GetYaxis()->SetTitle("EPOS UrQMD ON / OFF");
    gRatioEPOS->Draw("AL");

    //====================================================
    //==============Kstar/Ka ratio========================
    //====================================================

    TGraphErrors *gdNdy_Kstar_Off = GetGraph(fEPOS_New2, "g_dNdy_opt0");
    TGraphErrors *gdNdy_Kstar_On = GetGraph(fEPOS_New2, "g_dNdy_opt1");
    TGraphErrors *gdNdy_K0s_Off = GetGraph(fEPOS_New2, "g_stable_dNdy_opt0_K0s");
    TGraphErrors *gdNdy_K0s_On = GetGraph(fEPOS_New2, "g_stable_dNdy_opt1_K0s");
    TGraphErrors *gdNdy_Kaon_Off = GetGraph(fEPOS_New2, "g_stable_dNdy_opt0_kaon");
    TGraphErrors *gdNdy_Kaon_On = GetGraph(fEPOS_New2, "g_stable_dNdy_opt1_kaon");
    TGraphErrors *gdNdy_Pion_Off = GetGraph(fEPOS_New2, "g_stable_dNdy_opt0_pi");
    TGraphErrors *gdNdy_Pion_On = GetGraph(fEPOS_New2, "g_stable_dNdy_opt1_pi");
    TGraphErrors *gdNdy_Phi_Off = GetGraph(fEPOS_New2, "g_phi_dNdy_opt0");
    TGraphErrors *gdNdy_Phi_On = GetGraph(fEPOS_New2, "g_phi_dNdy_opt1");

    TGraphErrors *gKstarKaRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gKstarKaRatio_EPOSNew_UrQMDOn = new TGraphErrors();
    TGraphErrors *gKstarPiRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gKstarPiRatio_EPOSNew_UrQMDOn = new TGraphErrors();
    TGraphErrors *gKstarK0sRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gKstarK0sRatio_EPOSNew_UrQMDOn = new TGraphErrors();

    TGraphErrors *gPhiKaRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gPhiKaRatio_EPOSNew_UrQMDOn = new TGraphErrors();
    TGraphErrors *gPhiPiRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gPhiPiRatio_EPOSNew_UrQMDOn = new TGraphErrors();
    TGraphErrors *gPhiK0sRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gPhiK0sRatio_EPOSNew_UrQMDOn = new TGraphErrors();

    TGraphErrors *gKstarPhiRatio_EPOSNew_UrQMDOff = new TGraphErrors();
    TGraphErrors *gKstarPhiRatio_EPOSNew_UrQMDOn = new TGraphErrors();

    int graphPointsEPOSNew = gdNdy_Kaon_On->GetN();
    for (int i = 0; i < graphPointsEPOSNew; i++)
    {
        double xK0s, yK0s, yKaon, yPion, yKstarOff, yKstarOn;
        double yK0sON, yKaonON, yPionON, yPhiOff, yPhiOn;

        gdNdy_Phi_Off->GetPoint(i, xK0s, yPhiOff);
        gdNdy_Phi_On->GetPoint(i, xK0s, yPhiOn);
        cout << "xPhi: " << xK0s << ", yPhiOff: " << yPhiOff << ", yPhiOn: " << yPhiOn << endl;

        gdNdy_Kstar_Off->GetPoint(i, xK0s, yKstarOff);
        gdNdy_Kstar_On->GetPoint(i, xK0s, yKstarOn);
        cout << "xKstar: " << xK0s << ", yKstarOff: " << yKstarOff << ", yKstarOn: " << yKstarOn << endl;

        gdNdy_K0s_Off->GetPoint(i, xK0s, yK0s);
        gdNdy_Kaon_Off->GetPoint(i, xK0s, yKaon);
        gdNdy_Pion_Off->GetPoint(i, xK0s, yPion);

        gdNdy_K0s_On->GetPoint(i, xK0s, yK0sON);
        gdNdy_Kaon_On->GetPoint(i, xK0s, yKaonON);
        gdNdy_Pion_On->GetPoint(i, xK0s, yPionON);

        gKstarKaRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yKstarOff / yKaon);
        gKstarKaRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yKstarOn / yKaon);
        gKstarPiRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yKstarOff / yPion);
        gKstarPiRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yKstarOn / yPion);
        gKstarK0sRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yKstarOff / yK0s);
        gKstarK0sRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yKstarOn / yK0s);

        gKstarPhiRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yKstarOff / yPhiOff);
        gKstarPhiRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yKstarOn / yPhiOn);
    }

    int graphPointsEPOSNewPhi = gdNdy_Phi_On->GetN();
    for (int i = 0; i < graphPointsEPOSNewPhi; i++)
    {
        double xK0s, yK0s, yKaon, yPion, yPhiOff, yPhiOn;
        double yK0sON, yKaonON, yPionON;
        gdNdy_Phi_Off->GetPoint(i, xK0s, yPhiOff);
        gdNdy_Phi_On->GetPoint(i, xK0s, yPhiOn);

        gdNdy_K0s_Off->GetPoint(i, xK0s, yK0s);
        gdNdy_Kaon_Off->GetPoint(i, xK0s, yKaon);
        gdNdy_Pion_Off->GetPoint(i, xK0s, yPion);

        gdNdy_K0s_On->GetPoint(i, xK0s, yK0sON);
        gdNdy_Kaon_On->GetPoint(i, xK0s, yKaonON);
        gdNdy_Pion_On->GetPoint(i, xK0s, yPionON);

        gPhiKaRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yPhiOff / yKaon);
        gPhiKaRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yPhiOn / yKaon);
        gPhiPiRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yPhiOff / yPion);
        gPhiPiRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yPhiOn / yPion);
        gPhiK0sRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yPhiOff / yK0s);
        gPhiK0sRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yPhiOn / yK0s);
    }

    TGraphErrors *gKstarKaRatio_UrQMDOffSmooth = SmoothGraph(gKstarKaRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gKstarKaRatio_UrQMDOnSmooth = SmoothGraph(gKstarKaRatio_EPOSNew_UrQMDOn);
    TGraphErrors *gKstarPiRatio_UrQMDOffSmooth = SmoothGraph(gKstarPiRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gKstarPiRatio_UrQMDOnSmooth = SmoothGraph(gKstarPiRatio_EPOSNew_UrQMDOn);
    TGraphErrors *gKstarK0sRatio_UrQMDOffSmooth = SmoothGraph(gKstarK0sRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gKstarK0sRatio_UrQMDOnSmooth = SmoothGraph(gKstarK0sRatio_EPOSNew_UrQMDOn);

    TGraphErrors *gPhiKaRatio_UrQMDOffSmooth = SmoothGraph(gPhiKaRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gPhiKaRatio_UrQMDOnSmooth = SmoothGraph(gPhiKaRatio_EPOSNew_UrQMDOn);
    TGraphErrors *gPhiPiRatio_UrQMDOffSmooth = SmoothGraph(gPhiPiRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gPhiPiRatio_UrQMDOnSmooth = SmoothGraph(gPhiPiRatio_EPOSNew_UrQMDOn);
    TGraphErrors *gPhiK0sRatio_UrQMDOffSmooth = SmoothGraph(gPhiK0sRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gPhiK0sRatio_UrQMDOnSmooth = SmoothGraph(gPhiK0sRatio_EPOSNew_UrQMDOn);

    TGraphErrors *gKstarPhiRatio_UrQMDOffSmooth = SmoothGraph(gKstarPhiRatio_EPOSNew_UrQMDOff);
    TGraphErrors *gKstarPhiRatio_UrQMDOnSmooth = SmoothGraph(gKstarPhiRatio_EPOSNew_UrQMDOn);

    TCanvas *cKstarKaRatioEPOS = new TCanvas("cKstarKaRatioEPOS", "cKstarKaRatioEPOS", 720, 720);
    cKstarKaRatioEPOS->Divide(1, 2, 0.01, 0.01);
    cKstarKaRatioEPOS->cd(1);
    gPad->SetTopMargin(0.0);
    gPad->SetBottomMargin(0.0);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.001);
    SetGraphErrorStyle(gKstarKaRatio_UrQMDOffSmooth);
    gKstarKaRatio_UrQMDOffSmooth->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gKstarKaRatio_UrQMDOffSmooth->GetYaxis()->SetTitle("K* / K");
    gKstarKaRatio_UrQMDOffSmooth->GetYaxis()->SetRangeUser(0.27, 0.41);
    gKstarKaRatio_UrQMDOffSmooth->SetLineColor(kRed);
    gKstarKaRatio_UrQMDOffSmooth->SetLineStyle(2);
    gKstarKaRatio_UrQMDOffSmooth->SetLineWidth(3);
    gKstarKaRatio_UrQMDOffSmooth->Draw("AL");
    gKstarKaRatio_UrQMDOnSmooth->SetLineColor(kBlue - 6);
    gKstarKaRatio_UrQMDOnSmooth->SetLineStyle(9);
    gKstarKaRatio_UrQMDOnSmooth->SetLineWidth(3);
    gKstarKaRatio_UrQMDOnSmooth->Draw("L same");

    TLegend *legKstarKaRatioEPOS = new TLegend(0.22, 0.72, 0.85, 0.85);
    legKstarKaRatioEPOS->SetBorderSize(0);
    legKstarKaRatioEPOS->SetFillStyle(0);
    legKstarKaRatioEPOS->SetTextSize(0.05);
    legKstarKaRatioEPOS->AddEntry(gKstarKaRatio_UrQMDOffSmooth, "EPOS UrQMD OFF", "l");
    legKstarKaRatioEPOS->AddEntry(gKstarKaRatio_UrQMDOnSmooth, "EPOS UrQMD ON", "l");
    legKstarKaRatioEPOS->Draw();

    cKstarKaRatioEPOS->cd(2);
    gPad->SetTopMargin(0.0);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.001);
    int totalPointsEPOS_KstarKa = gKstarKaRatio_UrQMDOffSmooth->GetN();
    TGraphErrors *gRatioKstarKaEPOS = new TGraphErrors(totalPointsEPOS_KstarKa);
    for (int i = 0; i < totalPointsEPOS_KstarKa; i++)
    {
        double x, yOff, yOn;
        gKstarKaRatio_UrQMDOffSmooth->GetPoint(i, x, yOff);
        gKstarKaRatio_UrQMDOnSmooth->GetPoint(i, x, yOn);
        double ratio = (yOff != 0) ? yOn / yOff : 0;
        gRatioKstarKaEPOS->SetPoint(i, x, ratio);
        gRatioKstarKaEPOS->SetPointError(i, 0, 0);
    }
    SetGraphErrorStyle(gRatioKstarKaEPOS);
    gRatioKstarKaEPOS->SetLineColor(kBlack);
    gRatioKstarKaEPOS->SetLineStyle(1);
    gRatioKstarKaEPOS->SetLineWidth(2);
    gRatioKstarKaEPOS->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gRatioKstarKaEPOS->GetYaxis()->SetTitle("EPOS UrQMD ON / OFF");
    gRatioKstarKaEPOS->Draw("AL");
}

TGraphErrors *MakeRatio(const TGraphErrors *numerator, const TGraphErrors *denominator, bool isModel = false, float CorrFactorDen = 1)
{
    if (numerator == nullptr || denominator == nullptr)
    {
        cout << "Error: null input graph while building ratio" << endl;
        return nullptr;
    }

    if (numerator->GetN() != denominator->GetN())
    {
        cout << "Error: ratio inputs have different number of points" << endl;
        return nullptr;
    }

    TGraphErrors *ratio = new TGraphErrors(numerator->GetN());
    for (int i = 0; i < numerator->GetN(); ++i)
    {
        double xNumerator = 0;
        double yNumerator = 0;
        double xDenominator = 0;
        double yDenominator = 0;
        numerator->GetPoint(i, xNumerator, yNumerator);
        denominator->GetPoint(i, xDenominator, yDenominator);

        // if (std::round(xNumerator * 10.0) != std::round(xDenominator * 10.0))
        // {
        //     cout << "Warning!!!!! x values do not match for ratio at point "
        //          << i << ": " << xNumerator << " vs " << xDenominator
        //          << ", for graphs " << numerator->GetName()
        //          << " and " << denominator->GetName() << endl;
        // }

        yDenominator *= CorrFactorDen; // Apply correction factor to denominator if needed

        double yRatio = (yDenominator != 0) ? yNumerator / yDenominator : 0;
        ratio->SetPoint(i, xNumerator, yRatio);

        double xErrorNumerator = numerator->GetErrorX(i);
        double numeratorError = numerator->GetErrorY(i);
        double denominatorError = denominator->GetErrorY(i);
        denominatorError *= CorrFactorDen; // Apply correction factor to denominator error if needed
        double yRatioError = 0;
        if (yDenominator != 0)
        {
            yRatioError = sqrt(pow(numeratorError / yDenominator, 2) + pow(yNumerator * denominatorError / (yDenominator * yDenominator), 2));
        }
        if (isModel)
            ratio->SetPointError(i, 0, 0);
        else
            ratio->SetPointError(i, xErrorNumerator, yRatioError);
    }

    SetGraphErrorStyle(ratio);

    return ratio;
}

TGraphErrors **MakeRatioUncorr(TGraphErrors **numerator, TGraphErrors **denominator, float CorrFactorDen = 1, int nGraphs = 3)
{
    if (!numerator || !denominator)
    {
        cout << "Error: null input graph array" << endl;
        return nullptr;
    }

    if (nGraphs < 3)
    {
        cout << "Error: nGraphs must be at least 3" << endl;
        return nullptr;
    }

    TGraphErrors **ratio = new TGraphErrors *[3];

    for (int ig = 0; ig < 3; ig++)
    {
        if (!numerator[ig] || !denominator[ig])
        {
            cout << "Error: null graph at index " << ig << endl;
            delete[] ratio;
            return nullptr;
        }

        if (numerator[ig]->GetN() != denominator[ig]->GetN())
        {
            cout << "Error: different number of points in numerator " << "and denominator graphs" << endl;

            cout << "ig = " << ig << ", numerator = " << numerator[ig]->GetN() << ", denominator = " << denominator[ig]->GetN() << endl;

            delete[] ratio;
            return nullptr;
        }

        ratio[ig] = new TGraphErrors(numerator[ig]->GetN());

        for (int i = 0; i < numerator[ig]->GetN(); i++)
        {
            double xNum, yNum;
            double xDen, yDen;

            numerator[ig]->GetPoint(i, xNum, yNum);
            denominator[ig]->GetPoint(i, xDen, yDen);

            // Apply correction factor to denominator
            yDen *= CorrFactorDen;

            double yRatio = (yDen != 0.) ? yNum / yDen : 0.;

            ratio[ig]->SetPoint(i, xNum, yRatio);

            // ---------------------------------------------------------
            // Select the errors according to the requested definition
            // ---------------------------------------------------------

            double numErr = 0.;
            double denErr = 0.;

            if (ig == 0)
            {
                // Ratio 0:
                // stat(num) / stat(den)
                numErr = numerator[0]->GetErrorY(i);
                denErr = denominator[0]->GetErrorY(i);
            }
            else if (ig == 1)
            {
                // Ratio 1:
                // uncorr(num) / total(den) (wrong)
                // total (num) / total (den) (correct)
                numErr = numerator[1]->GetErrorY(i);
                denErr = denominator[1]->GetErrorY(i);
            }
            else if (ig == 2)
            {
                // Ratio 2:
                // uncorr(num) / uncorr(den)
                numErr = numerator[2]->GetErrorY(i);
                denErr = denominator[2]->GetErrorY(i);
            }

            // Correction factor also scales the denominator uncertainty
            denErr *= CorrFactorDen;

            double yRatioErr = 0.;

            if (yDen != 0.)
            {
                yRatioErr = sqrt(
                    pow(numErr / yDen, 2) +
                    pow(yNum * denErr / (yDen * yDen), 2));
            }

            double xErr = numerator[ig]->GetErrorX(i);

            ratio[ig]->SetPointError(i, xErr, yRatioErr);
        }

        SetGraphErrorStyle(ratio[ig]);
    }

    return ratio;
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

TGraphErrors *GetGraph(TFile *f, const string &name)
{
    TGraphErrors *graph = (TGraphErrors *)f->Get(name.c_str());

    if (!graph || graph == nullptr)
    {
        cout << "Error: graph " << name
             << " not found in file " << f->GetName() << endl;
        return nullptr;
    }

    SetGraphErrorStyle(graph);
    graph->SetTitle(0);
    return graph;
}

void ScaleGraph(TGraph *gr, double scale)
{
    if (!gr)
        return;

    for (int i = 0; i < gr->GetN(); ++i)
    {
        double x, y;
        gr->GetPoint(i, x, y);
        gr->SetPoint(i, x, y * scale);
    }

    if (auto *ge = dynamic_cast<TGraphErrors *>(gr))
    {
        for (int i = 0; i < ge->GetN(); ++i)
        {
            ge->SetPointError(i, ge->GetErrorX(i), ge->GetErrorY(i) * scale);
        }
    }
}

TGraphErrors **DivideByMult(TGraphErrors **gr, double WhichMultPoint, double tolerance = 0.5, int noOfGraphs = 3)
{
    if (!gr || !gr[0] || !gr[2])
    {
        cout << "Error: null graph array or required graphs" << endl;
        return nullptr;
    }

    if (noOfGraphs < 3)
    {
        cout << "Error: noOfGraphs must be at least 3" << endl;
        return nullptr;
    }

    // Output:
    // grCopy[0] = ratio with statistical uncertainty
    // grCopy[1] = ratio with uncorrelated uncertainty
    TGraphErrors **grCopy = new TGraphErrors *[2];

    // ------------------------------------------------------------
    // Find the reference multiplicity point
    // ------------------------------------------------------------

    int refIndex = -1;

    if (WhichMultPoint >= 0)
    {
        for (int i = 0; i < gr[0]->GetN(); i++)
        {
            double x, y;
            gr[0]->GetPoint(i, x, y);

            if (fabs(x - WhichMultPoint) < tolerance)
            {
                refIndex = i;
                break;
            }
        }

        if (refIndex < 0)
        {
            cout << "Error: Could not find multiplicity point "
                 << WhichMultPoint << endl;

            delete[] grCopy;
            return nullptr;
        }
    }
    else
    {
        // Use the lowest multiplicity point
        double xFirst, yFirst;
        double xLast, yLast;

        gr[0]->GetPoint(0, xFirst, yFirst);
        gr[0]->GetPoint(gr[0]->GetN() - 1, xLast, yLast);

        refIndex = (xFirst < xLast) ? 0 : gr[0]->GetN() - 1;
    }

    // ------------------------------------------------------------
    // Reference value
    // ------------------------------------------------------------

    double xRef, yRef;
    gr[0]->GetPoint(refIndex, xRef, yRef);

    if (yRef == 0.)
    {
        cout << "Error: Reference value is zero" << endl;
        delete[] grCopy;
        return nullptr;
    }

    // Reference uncertainties
    double yRefStatErr = gr[0]->GetErrorY(refIndex);
    double yRefUncorrErr = gr[2]->GetErrorY(refIndex);

    // ------------------------------------------------------------
    // Create the two output graphs
    // ------------------------------------------------------------

    grCopy[0] = (TGraphErrors *)gr[0]->Clone();
    grCopy[1] = (TGraphErrors *)gr[2]->Clone();

    // ------------------------------------------------------------
    // Calculate ratios
    // ------------------------------------------------------------

    for (int i = 0; i < gr[0]->GetN(); i++)
    {
        double x, y;
        gr[0]->GetPoint(i, x, y);

        double xerr = gr[0]->GetErrorX(i);

        // ========================================================
        // 1. STATISTICAL RATIO
        // ========================================================

        double yStatErr = gr[0]->GetErrorY(i);

        double ratio = y / yRef;

        double ratioStatErr = sqrt(
            pow(yStatErr / yRef, 2) +
            pow(y * yRefStatErr / (yRef * yRef), 2));

        grCopy[0]->SetPoint(i, x, ratio);
        grCopy[0]->SetPointError(i, xerr, ratioStatErr);

        // ========================================================
        // 2. UNCORRELATED SYSTEMATIC RATIO
        // ========================================================

        double yUncorrErr = gr[2]->GetErrorY(i);

        double ratioUncorrErr = sqrt(pow(yUncorrErr / yRef, 2) + pow(y * yRefUncorrErr / (yRef * yRef), 2));

        grCopy[1]->SetPoint(i, x, ratio);
        grCopy[1]->SetPointError(i, xerr, ratioUncorrErr);
    }

    SetGraphErrorStyle(grCopy[0]);
    SetGraphErrorStyle(grCopy[1]);

    return grCopy;
}

void RestrictModelXaxis(TGraphErrors *gr, double xMin, double xMax)
{
    for (int i = gr->GetN() - 1; i >= 0; --i)
    {
        double x, y;
        gr->GetPoint(i, x, y);
        if (x < xMin || x > xMax)
        {
            gr->RemovePoint(i);
        }
    }
}

TGraphErrors *DivideByMultModel(TGraphErrors *gr, double WhichMultPoint, double tolerance = 0.5)
{
    TGraphErrors *grCopy = (TGraphErrors *)gr->Clone();

    double yGivenMult = 1.0;
    double yGivenMultErr = 0.0;

    for (int i = 0; i < gr->GetN(); i++)
    {
        double x, y;
        gr->GetPoint(i, x, y);

        if (fabs(x - WhichMultPoint) < tolerance)
        {
            yGivenMult = y;
            yGivenMultErr = gr->GetErrorY(i);
            break;
        }
    }

    for (int i = 0; i < gr->GetN(); i++)
    {
        double x, y;
        gr->GetPoint(i, x, y);

        double yerr = gr->GetErrorY(i);
        double xerr = gr->GetErrorX(i);

        // Normalization point
        if (fabs(x - WhichMultPoint) < tolerance)
        {
            grCopy->SetPoint(i, x, 1.0);
            grCopy->SetPointError(i, xerr, 0.0);
            continue;
        }

        double ratio = y / yGivenMult;

        double ratioErr = sqrt(pow(yerr / yGivenMult, 2) + pow(yGivenMultErr * y / (yGivenMult * yGivenMult), 2));

        grCopy->SetPoint(i, x, ratio);
        grCopy->SetPointError(i, xerr, ratioErr);
    }

    return grCopy;
}

TGraphErrors *SmoothGraph(const TGraphErrors *gr)
{
    int n = gr->GetN();

    TGraphErrors *grSmooth = new TGraphErrors(n);

    for (int i = 0; i < n; i++)
    {

        double x, y;
        gr->GetPoint(i, x, y);

        double ySmooth = 0.0;
        double errSmooth = 0.0;
        int nPoints = 0;

        // Average neighbouring points
        for (int j = i - 1; j <= i + 1; j++)
        {

            if (j < 0 || j >= n)
                continue;

            double xj, yj;
            gr->GetPoint(j, xj, yj);

            double eyj = gr->GetErrorY(j);

            ySmooth += yj;
            errSmooth += eyj * eyj;

            nPoints++;
        }

        ySmooth /= nPoints;

        // Error propagation for the average
        errSmooth = std::sqrt(errSmooth) / nPoints;

        grSmooth->SetPoint(i, x, ySmooth);
        grSmooth->SetPointError(
            i,
            gr->GetErrorX(i),
            errSmooth);
    }

    return grSmooth;
}
