#include "src/common.h"
bool useLSbkg = 1;
TH1 *GetNorBkg(TH1 *hSig, TH1 *hBkg, Int_t isnormleft, Int_t isnormright); // returns a pointer to new TH1 histogram (GetNorBkg) representing the normalized background distribution. This is just a prototype, the main function is defined below the void signal loop.

// function
Double_t BreitWignerpoly2(Double_t *x, Double_t *par)
{
    // BW = (0.5 * yield * Gamma) / ((x - Mean)^2 + 0.25 * Gamma^2)
    double BW = (0.5 * par[2] * par[1] / TMath::Pi()) / ((x[0] - par[0]) * (x[0] - par[0]) + 0.25 * par[1] * par[1]);
    // double poly = par[5] + par[4]*x[0] + par[3]*x[0]*x[0];
    double poly3 = par[6] + par[5] * x[0] + par[4] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
    return (BW + poly3);
}
Double_t BW(Double_t *x, Double_t *par)
{
    return (0.5 * par[2] * par[1] / TMath::Pi() / ((x[0] - par[0]) * (x[0] - par[0]) + 0.25 * par[1] * par[1]));
}
Double_t polynomial3(Double_t *x, Double_t *par)
{
    double poly3 = par[3] + par[2] * x[0] + par[1] * x[0] * x[0] + par[0] * x[0] * x[0] * x[0];
    return (poly3);
}
  const  int Npt = 100;
// fiiting parameters array///////////////////////
  Double_t Yield[Npt];
  Double_t Mass[Npt];
  Double_t Width[Npt];
  Double_t poly0[Npt];
  Double_t poly1[Npt];
  Double_t poly2[Npt];
  Double_t poly3[Npt];
  Double_t ErrorMass[Npt];
  Double_t ErrorWidth[Npt];
  Double_t ErrorYield[Npt];
  Double_t Chi2Ndf[Npt];

// double lowfitrange[] = {0.7, 0.75, 0.78, 0.78, 0.72, 0.72, 0.7, 0.72, 0.72, 0.72, 0.72, 0.68, 0.66, 0.74};
// double highfitrange[] = {1.1, 1.1, 1.06, 1.1, 1.06, 1.05, 1.04, 1.04, 1.06, 1.06, 1.06, 1.12, 1.05, 1.12};

// normalization
Int_t normleft = 1;
std::vector<Double_t> NormalizeRange_L = {0.65, 0.75}; // Norminal
Int_t normright = 0;
std::vector<Double_t> NormalizeRange_R = {1.7, 1.9}; // Norminal

// Fit range
vector<Double_t> rangefit = {0.8, 0.85, 0.95, 1.05}; // left edge, signal peak left edge, signal peak right edge, right
// vector<Double_t> range_fit2 = {0.77, 0.85, 0.95, 1.};

void signal2(int chosenDir = 0, int chosenAnti = -1, int chosenPart = -1, int chosenCent = -1)
{
    // gStyle->SetOptStat(0);
    // gStyle->SetOptTitle(0);
    // gStyle->SetLegendBorderSize(0);
    gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1110);
  gStyle->SetStatColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFillColor(0);
  gStyle->SetLineColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleColor(1);

    TFile inputFile(kDataFilename.data());
    TString outputFileName = kSignalOutput.data();
    TFile *outputFile = (TFile *)new TFile(outputFileName, "recreate");

    // Fitting functions
    TCanvas *fCanvas = new TCanvas("FitCheck", "Fits", kCanvasW, kCanvasH);
    auto th3_Signal = (TH3F *)inputFile.Get(Form("%s/h3k892invmassDS", kOutputName.c_str()));
    auto th3_ME = (TH3F *)inputFile.Get(Form("%s/h3k892invmassME", kOutputName.c_str()));
    auto th3_LS = (TH3F *)inputFile.Get(Form("%s/h3k892invmassLS", kOutputName.c_str()));
    outputFile->cd();

    // Preparation
    std::vector<Double_t> pt_points_e = {};
    for (int i = 0; i < (int)kpTbin.size() - 1; i++)
        pt_points_e.push_back(kpTbin.at(i + 1) - kpTbin[i]); // takes the width and add in pt_points vector
    // Arrays
    std::vector<Double_t> fitmean = {};
    std::vector<Double_t> fitmean_err = {};
    std::vector<Double_t> fitsigma = {};
    std::vector<Double_t> fitsigma_err = {};
    std::vector<Double_t> fitgamma = {};
    std::vector<Double_t> fitgamma_err = {};
    std::vector<Double_t> RawYield = {};
    std::vector<Double_t> RawYield_err = {};
    // Canvas
    TCanvas *cInvMasspT = new TCanvas("cInvMasspT", "", kCanvasW, kCanvasH);
    cInvMasspT->Draw();
    cInvMasspT->SetTickx();
    cInvMasspT->SetTicky();
    cInvMasspT->SetLeftMargin(0.13);
    cInvMasspT->SetTopMargin(0.07);
    cInvMasspT->SetRightMargin(0.01);
    TCanvas *cInvMasspTFit = new TCanvas("cInvMasspTFit", "", kCanvasW, kCanvasH);
    cInvMasspTFit->Draw();
    cInvMasspTFit->SetTickx();
    cInvMasspTFit->SetTicky();
    cInvMasspTFit->SetLeftMargin(0.13);
    cInvMasspTFit->SetTopMargin(0.07);
    cInvMasspTFit->SetRightMargin(0.01);
    // Text
    TLatex *t2 = new TLatex();
    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.04);
    // Fit
    TString fitName = "BW";
    TString bkgFormula = "pol2";
    vector<Double_t> par_fit = {0.98, 0.004, 0.0487}; // first is constant coefficient, second and third is linear and quadratic coefficients.
    vector<Int_t> fix_fit = {0, 0, 1};
    TH1 *hBkg_orig = nullptr;
    TH1 *hBkg = nullptr;
    for (Int_t iPt = 0, nPt = kpTbin.size() - 1; iPt < nPt; ++iPt)
    {
        auto hSig = th3_Signal->ProjectionZ(
            Form("hSig_%d", iPt), -1, -1, th3_Signal->GetYaxis()->FindBin(kpTbin[iPt]),
            th3_Signal->GetYaxis()->FindBin(kpTbin[iPt + 1])); // from the 3D histogram th3_signal,  1D histograms of invariant mass are made for which all multiplicity is taken (-1, -1 for x axis) and the transverse momentum is selected for difference 0.1, i.e. 0.1 to 0.2, then 0.2 to 0.3 .......
        if (useLSbkg)
        {
            hBkg_orig = th3_LS->ProjectionZ(
                Form("hBkg_%d", iPt), -1, -1, th3_LS->GetYaxis()->FindBin(kpTbin[iPt]),
                th3_LS->GetYaxis()->FindBin(kpTbin[iPt + 1]));
            hBkg = (TH1 *)hBkg_orig->Clone(Form("hBkg_%d", iPt));
        }
        else
        {
            hBkg_orig = th3_ME->ProjectionZ(
                Form("hBkg_%d", iPt), -1, -1, th3_ME->GetYaxis()->FindBin(kpTbin[iPt]),
                th3_ME->GetYaxis()->FindBin(kpTbin[iPt + 1]));
            hBkg = GetNorBkg(hSig, hBkg_orig, normleft, normright); // normalize the ME(Mixed Event) bkg
        }
        hBkg->SetLineColor(2);
        // rebin
        hSig->Rebin(kRebin);
        hBkg->Rebin(kRebin);
        hSig->Write(Form("hSignal_%i", iPt));
        hBkg->Write(Form("hBkg_%i", iPt));

        hSig->GetXaxis()->SetRangeUser(kDrawRange[0], kDrawRange[1]); // originally the range was from 0.6 to 1.5 which is now changed to 0.7 to 1.4
        cInvMasspT->cd();
        hSig->Draw("E");
        hBkg->Draw("E same");

        auto legend = new TLegend(0.65, 0.67, 0.9, 0.87);
        legend->SetFillStyle(0);
        legend->AddEntry(hSig, "data", "LE");
        legend->AddEntry(hBkg, "Bkg", "LE");
        // legend->AddEntry(temp, "Normalization Region", "F");
        legend->Draw();
        t2->DrawLatex(0.22, 0.95, "#bf{K(892)^{0} #rightarrow #pi + K}");
        t2->DrawLatex(0.75, 0.95,
                      Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", kpTbin[iPt],
                           kpTbin[iPt + 1]));
        cInvMasspT->Write(Form("cInvMasspT_%i", iPt));
        // cInvMasspT->SaveAs(Form("%s/cInvMasspT_%i.pdf", kFiguresFolder.c_str(), iPt));
        cInvMasspT->SaveAs(Form("%s/cInvMasspT_%i.png", kFiguresFolder.c_str(), iPt));

        TH1F *Sigfit = (TH1F *)hSig->Clone();
        TH1F *Bkgfit = (TH1F *)hBkg->Clone();
        Sigfit->Add(Bkgfit, -1); // signal - bkg
        // Sigfit->GetXaxis()->SetRangeUser(kDrawFitRange[0], kDrawFitRange[1]);

        //     /// FIT procedure from AliPhysics/PWGLF/RESONANCES/macros/utils/VoigtianFit.C
        //     TString name = Sigfit->GetName();
        //     std::cout << "using BWFit():\n  h=" << name << "\n  background bkgFormula=" << bkgFormula << std::endl;
        //     std::cout << "  mass=" << par_fit[0];
        //     if (fix_fit[0])
        //         std::cout << " (fixed)" << std::endl;   //????????????????
        //     else
        //         std::cout << " (free)" << std::endl;
        //     std::cout << "  resolution=" << par_fit[1];
        //     if (fix_fit[1])
        //         std::cout << " (fixed)" << std::endl;
        //     else
        //         std::cout << " (free)" << std::endl;
        //     std::cout << "  width=" << par_fit[2];
        //     if (fix_fit[2])
        //         std::cout << " (fixed)" << std::endl;
        //     else
        //         std::cout << " (free)" << std::endl;
        //     std::cout << "  range=" << range_fit[0] << " " << range_fit[1] << " " << range_fit[2] << " " << range_fit[3] << std::endl;

        //     int jLoop;

        //     // create copy of histogram h with peak removed
            TH1 *SigWithoutPeak = (TH1 *)Sigfit->Clone(Form("%s_nopeak", name.Data()));
            for (jLoop = Sigfit->GetXaxis()->FindBin(1.000001 * range_fit[1]); jLoop <= Sigfit->GetXaxis()->FindBin(0.999999 * range_fit[2]); jLoop++)
            {
                SigWithoutPeak->SetBinContent(jLoop, 0.);
                SigWithoutPeak->SetBinError(jLoop, 0.);
            }

        //     // get initial estimate of background
        //     TF1 *fBkg;
        //     // if (iPt<4)
        //     // {
        //     //     fBkg = new TF1(Form("%s_back", name.Data()), bkgFormula.Data(), range_fit2[0], range_fit2[3]);
        //     // }
        //     // else
        //     fBkg = new TF1(Form("%s_back", name.Data()), bkgFormula.Data(), range_fit[0], range_fit[3]);
        //     SigWithoutPeak->Fit(fBkg, "RQNI");

        //     // define peak fit function
        //     int vp = fBkg->GetNpar();
        //    TF1 *fSigBkg = new TF1(Form("%s_peak" , name.Data()), Form("%s+[%i]*TMath::Voigt(x-[%i],[%i],[%i])", bkgFormula.Data(), vp, vp + 1, vp + 2, vp + 3), range_fit[0], range_fit[3]); // vp is just the parameter number. the pol2 will take 3 parameters i.e 0,1,2 and voigt takes 4 parameters which are defined by using vp.

        //     // set initial parameter values, only the peak height is free.
        //     fSigBkg->SetParLimits(vp, 1e1, 1.e8);
        //     for (jLoop = 0; jLoop < vp; jLoop++)
        //     {
        //         fSigBkg->SetParameter(jLoop, fBkg->GetParameter(jLoop));  //initialize the parameters for pol2+voigt
        //         fSigBkg->FixParameter(jLoop, fBkg->GetParameter(jLoop)); //fixing the parameters of pol2 in pol2+voigt
        //     }
        //     fSigBkg->SetParameter(vp, Sigfit->GetBinContent(Sigfit->GetXaxis()->FindBin(0.5 * (range_fit[2] - range_fit[1]))) - fBkg->Eval(0.5 * (range_fit[2] - range_fit[1])));
        //     for (jLoop = 0; jLoop < 3; jLoop++)
        //     {
        //         fSigBkg->SetParameter(vp + jLoop + 1, par_fit[jLoop]);
        //         fSigBkg->FixParameter(vp + jLoop + 1, par_fit[jLoop]);
        //     }

        //     Sigfit->Fit(fSigBkg, "RQN");

        //     if (!fix_fit[2])
        //     { // release width
        //         fSigBkg->ReleaseParameter(vp + 3);
        //         fSigBkg->SetParError(vp + 3, 0.1 * fSigBkg->GetParameter(vp + 3));
        //         Sigfit->Fit(fSigBkg, "RQN");
        //     }

        //     if (!fix_fit[1])
        //     { // release resolution
        //         fSigBkg->ReleaseParameter(vp + 2);
        //         fSigBkg->SetParError(vp + 2, 0.1 * fSigBkg->GetParameter(vp + 2));
        //         Sigfit->Fit(fSigBkg, "RQN");
        //     }

        //     if (!fix_fit[0])
        //     { // release mass
        //         fSigBkg->ReleaseParameter(vp + 1);
        //         fSigBkg->SetParError(vp + 1, 0.1 * fSigBkg->GetParameter(vp + 1));
        //         Sigfit->Fit(fSigBkg, "RQN");
        //     }

        //     // release background constant parameter
        //     fSigBkg->ReleaseParameter(0);
        //     fSigBkg->SetParError(0, fBkg->GetParError(0));
        //     Sigfit->Fit(fSigBkg, "RQN");

        //     // release other background parameters
        //     for (jLoop = 1; jLoop < vp; jLoop++)
        //     {
        //         fSigBkg->ReleaseParameter(jLoop);
        //         fSigBkg->SetParError(jLoop, fBkg->GetParError(jLoop));
        //     }
        //     Sigfit->Fit(fSigBkg, "RQN");
        //     // Double_t chisq = fSigBkg->GetChisquare();
        //     // Double_t ndf = fSigBkg->GetNDF();
        //     // Double_t chindf = chisq/ndf;
        //     // cout<< " The value of chi-square is for iteration "<<iPt<<" is "<<chisq<<". The value of ndf is "<<ndf<<". The value of chi-square by ndf is "<<chindf<<endl;

        //     // Configure fit options
        //     fBkg->SetLineStyle(2);
        //     fBkg->SetLineColor(kBlue);

        //     // final fit
        //     std::cerr << "doing final fit" << std::endl;
        //     Sigfit->Fit(fSigBkg, "RQNI");

        //     /// FIT procedure DONE

        //     cInvMasspTFit->cd();
        //     Sigfit->Draw("E");
        //     fSigBkg->Draw("SAME");
        //     fBkg->Draw("SAME");
        //     t2->DrawLatex(0.22, 0.95, "#bf{K(892)^{0} #rightarrow #pi + K}");
        //     t2->DrawLatex(0.75, 0.95,
        //                   Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", kpTbin[iPt],
        //                        kpTbin[iPt + 1]));
        //     cInvMasspTFit->Write(Form("cInvMasspTFit_%i", iPt));
        //     cInvMasspTFit->SaveAs(Form("%s/cInvMasspTFit_%i.png", kFiguresFolder.c_str(), iPt));
        //     // cInvMasspTFit->SaveAs(Form("%s/cInvMasspTFit_%i.pdf", kFiguresFolder.c_str(), iPt));

        TF1 *fitFcn = new TF1("fitfunc", BreitWignerpoly2, rangefit[0], rangefit[3], 7); // sig+bkg fit function
        TF1 *fitFcn1 = new TF1("fitfunc1", polynomial3, rangefit[0], rangefit[3], 4);    // only residualbkg
        TF1 *fitFcn2 = new TF1("fitFcn2", BW, rangefit[0], rangefit[3], 3);              // only signal

        fitFcn->SetParLimits(0, 0.87, 0.91); // mass
        // fitFcn->FixParameter(1,0.047); //width
        fitFcn->SetParameter(1, 0.047); // width
        fitFcn->SetParameter(2, 1000);  // yield

        fitFcn->SetParameter(3, 100);
        fitFcn->SetParameter(4, 100);
        fitFcn->SetParameter(5, 100);
        fitFcn->SetParameter(6, 100);

        /*
        fitFcn->SetParLimits(3,-200,200);
        fitFcn->SetParLimits(4,-200,200);
        fitFcn->SetParLimits(5,-200,200);
        */

        // fitFcn->SetParNames("Mass","Width","Yield","A","B","C","D");
        fitFcn->SetParNames("Mass", "Width", "Yield", "C", "B", "A");
        Sigfit->Fit(fitFcn, "REBMS+");

        //**Extraction of fitting parameters******************************************************************************

        Double_t *par = fitFcn->GetParameters();

        Mass[iPt] = fitFcn->GetParameter(0);
        Width[iPt] = fitFcn->GetParameter(1);
        Yield[iPt] = fitFcn->GetParameter(2);
        /*poly0[iPt]=fitFcn->GetParameter(3);
        poly1[iPt]=fitFcn->GetParameter(4);
        poly2[iPt]=fitFcn->GetParameter(5);*/
        poly2[iPt] = fitFcn->GetParameter(3);
        poly1[iPt] = fitFcn->GetParameter(4);
        poly0[iPt] = fitFcn->GetParameter(5);
        poly3[iPt] = fitFcn->GetParameter(6);

        fitFcn2->SetParameters(&par[0]);
        fitFcn1->SetParameters(&par[3]);

        ErrorMass[iPt] = fitFcn->GetParError(0);
        ErrorWidth[iPt] = fitFcn->GetParError(1);
        ErrorYield[iPt] = fitFcn->GetParError(2);
        Chi2Ndf[iPt] = (fitFcn->GetChisquare()) / (fitFcn->GetNDF());

        cInvMasspTFit->cd();
            fitFcn1-> SetLineStyle(2);
            fitFcn1-> SetLineColor(8);
            fitFcn2-> SetLineStyle(2);
            fitFcn2-> SetLineColor(7);
            Sigfit->Draw("E");
            fitFcn1->Draw("SAME");
            fitFcn2->Draw("SAME");
            t2->DrawLatex(0.22, 0.95, "#bf{K(892)^{0} #rightarrow #pi + K}");
            t2->DrawLatex(0.75, 0.95,
                          Form("#bf{%.1f < #it{p}_{T} < %.1f GeV/c}", kpTbin[iPt],
                               kpTbin[iPt + 1]));
            cInvMasspTFit->Write(Form("cInvMasspTFit_%i", iPt));
            cInvMasspTFit->SaveAs(Form("%s/cInvMasspTFit_%i.png", kFiguresFolder2.c_str(), iPt));
    }

    outputFile->Close();
}

TH1 *GetNorBkg(TH1 *hSig, TH1 *hBkg, Int_t isnormleft, Int_t isnormright)
{
    Double_t normalization_data = 0.;
    Double_t normalization_mixed = 0.;
    if (isnormleft == 1)
    {
        normalization_data +=
            hSig->Integral(hSig->GetXaxis()->FindBin(NormalizeRange_L[0]),
                           hSig->GetXaxis()->FindBin(NormalizeRange_L[1]) - 1); // why -1 bin ??????
        normalization_mixed +=
            hBkg->Integral(hBkg->GetXaxis()->FindBin(NormalizeRange_L[0]),
                           hBkg->GetXaxis()->FindBin(NormalizeRange_L[1]) - 1);
    }
    if (isnormright == 1)
    {
        normalization_data +=
            hSig->Integral(hSig->GetXaxis()->FindBin(NormalizeRange_R[0]),
                           hSig->GetXaxis()->FindBin(NormalizeRange_R[1]) - 1);
        normalization_mixed +=
            hBkg->Integral(hBkg->GetXaxis()->FindBin(NormalizeRange_R[0]),
                           hBkg->GetXaxis()->FindBin(NormalizeRange_R[1]) - 1);
    }
    hBkg->Scale(normalization_data / normalization_mixed);
    hBkg->SetLineColor(kRed);
    hBkg->SetMarkerColor(kRed);

    return hBkg;
}