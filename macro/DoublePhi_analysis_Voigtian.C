// DoublePhi_analysis_Voigtian.C
//
// The single-phi fit, template smoothing fit, and direct M(phi phi) fit each
// have an independent background selector with identical numbering:
//   0 pol2, 1 pol3, 2 exp(pol2), 3 exp(pol3),
//   4 Chebyshev-2, 5 Chebyshev-3, 6 Bernstein-2, 7 Bernstein-3.
// The direct ROOT fit and RooStats local-p0 scan always use the same selected
// direct background.
//
// Optional template fit:
//   a global 2D m(phi1)-m(phi2) fit defines the SS/SB/BS/BB decomposition;
//   local posterior sums provide the SS and non-SS M(phi phi) templates.
//   Both template components are smoothed with templateBkgModel.
//   The signal window is excluded while fitting the template shapes.  The two
//   templates are then normalized in the sidebands, subtracted, and the
//   residual is fitted with a Voigtian.
//
// Public call:
//   DoublePhi_analysis_Voigtian(inputFile, sparseName, ptMin, ptMax,
//       deltaMMax, phiBkgModel, templateBkgModel, directBkgModel,
//       runTemplateFit, phiFitMin, phiFitMax, xResolution)
//
// Example: direct Chebyshev-3 fit
//   root -l -b -q 'DoublePhi_analysis_Voigtian.C+("AnalysisResults.root",\
//     "doublephimeson/SEMassDoublePhi",8,100,0.005,3,0,5,0,1.000,1.040,0.015)'
//
// Example: exp(pol3) direct fit plus Bernstein-2 template fit
//   root -l -b -q 'DoublePhi_analysis_Voigtian.C+("AnalysisResults.root",\
//     "doublephimeson/SEMassDoublePhi",8,100,0.005,3,6,3,1,1.000,1.040,0.015)'
//
// Assumed THnSparse axes:
//   0 M(phi phi), 1 pT(phi phi), 3 pair rapidity,
//   4 M(phi1), 5 M(phi2), 6 stored DeltaM_phi.
//
// All performed 2D fits, their input histograms, and the numerical fit summary
// are saved.  Ordinary projections use THnSparse::SetRange + Projection("E");
// the only explicit sparse loop is the posterior P_SS weighting.
//
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <initializer_list>
#include <limits>
#include <memory>
#include <utility>

#include "TFile.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TBox.h"
#include "TMath.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMatrixDSym.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TNamed.h"
#include "TString.h"
#include "TDirectory.h"
#include "TPaveText.h"
#include "TVirtualFitter.h"

#include "Math/MinimizerOptions.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGenericPdf.h"
#include "RooFormulaVar.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooBernstein.h"
#include "RooChebychev.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooMsgService.h"
#include "RooWorkspace.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HypoTestResult.h"

using namespace RooFit;

namespace {

const int kAxisMPair   = 0;
const int kAxisPtPair  = 1;
const int kAxisRapidity = 3;
const int kAxisM1      = 4;
const int kAxisM2      = 5;
const int kAxisDeltaM  = 6;

// Pair-rapidity selection shared by every native sparse projection and by the
// one explicit posterior loop.  It is configured once in DoublePhi_analysis.
double gRapidityMin = 0.0;
double gRapidityMax = 0.8;

double gRejectMin = -1.0;
double gRejectMax = -1.0;

// One selector is used consistently by the single-phi, template, and direct
// double-phi fits.
enum BackgroundModel {
  kPol2 = 0,
  kPol3 = 1,
  kExpPol2 = 2,
  kExpPol3 = 3,
  kChebyshev2 = 4,
  kChebyshev3 = 5,
  kBernstein2 = 6,
  kBernstein3 = 7
};

int    gContModel = kPol2;
int    gContNPar  = 4;
double gContXMin  = 0.0;
double gContXMax  = 1.0;

int    gPhiBackgroundModel = kExpPol3;
double gPhiFitMin = 1.000;
double gPhiFitMax = 1.040;

// Single-phi signal model used in every 2D m(phi1)-m(phi2) fit.
// The common Breit-Wigner mean and width are fitted globally.
double gPhiBWGammaInit = 0.0042;
double gPhiBWGammaMin = 0.0010;
double gPhiBWGammaMax = 0.0200;

bool IsValidBackgroundModel(int model)
{
  return model >= kPol2 && model <= kBernstein3;
}

int BackgroundOrder(int model)
{
  switch (model) {
    case kPol3:
    case kExpPol3:
    case kChebyshev3:
    case kBernstein3:
      return 3;
    default:
      return 2;
  }
}

const char* BackgroundLabel(int model)
{
  switch (model) {
    case kPol2:        return "pol2";
    case kPol3:        return "pol3";
    case kExpPol2:     return "exp(pol2)";
    case kExpPol3:     return "exp(pol3)";
    case kChebyshev2:  return "Chebyshev-2";
    case kChebyshev3:  return "Chebyshev-3";
    case kBernstein2:  return "Bernstein-2";
    case kBernstein3:  return "Bernstein-3";
    default:           return "unknown";
  }
}

TString BackgroundTag(int model)
{
  switch (model) {
    case kPol2:        return "pol2";
    case kPol3:        return "pol3";
    case kExpPol2:     return "expPol2";
    case kExpPol3:     return "expPol3";
    case kChebyshev2:  return "Chebyshev2";
    case kChebyshev3:  return "Chebyshev3";
    case kBernstein2:  return "Bernstein2";
    case kBernstein3:  return "Bernstein3";
    default:           return "unknown";
  }
}

double BackgroundValue(int model,
                       double x,
                       double xmin,
                       double xmax,
                       const double* p)
{
  if (!p || !(xmax > xmin)) return 0.0;

  const double u = std::max(0.0, std::min(1.0, (x-xmin)/(xmax-xmin)));
  const double t = 2.0*u-1.0;
  const double t2 = t*t;
  const double t3 = t2*t;
  const int order = BackgroundOrder(model);

  if (model == kExpPol2 || model == kExpPol3) {
    double exponent = p[0]+p[1]*t+p[2]*t2;
    if (order == 3) exponent += p[3]*t3;
    exponent = std::max(-700.0, std::min(700.0, exponent));
    return std::exp(exponent);
  }

  double shape = 0.0;
  if (model == kChebyshev2 || model == kChebyshev3) {
    shape = 1.0+p[1]*t+p[2]*(2.0*t2-1.0);
    if (order == 3) shape += p[3]*(4.0*t3-3.0*t);
  } else if (model == kBernstein2) {
    const double v = 1.0-u;
    shape = v*v+p[1]*(2.0*u*v)+p[2]*u*u;
  } else if (model == kBernstein3) {
    const double v = 1.0-u;
    shape = v*v*v+p[1]*(3.0*u*v*v)+p[2]*(3.0*u*u*v)+p[3]*u*u*u;
  } else {
    shape = 1.0+p[1]*t+p[2]*t2;
    if (order == 3) shape += p[3]*t3;
  }

  return std::max(1.0e-12, p[0]*shape);
}

void ConfigureBackgroundParameters(TF1* function,
                                   int offset,
                                   int model,
                                   double average,
                                   double maximum)
{
  if (!function) return;
  const int order = BackgroundOrder(model);
  const bool isExponential = model == kExpPol2 || model == kExpPol3;
  const bool isBernstein = model == kBernstein2 || model == kBernstein3;
  const bool isChebyshev = model == kChebyshev2 || model == kChebyshev3;

  function->SetParName(offset, isExponential ? "logB_{0}" : "B_{0}");
  if (isExponential) {
    const double logAverage = std::log(std::max(1.0e-12, average));
    function->SetParameter(offset, logAverage);
    function->SetParLimits(offset, logAverage-10.0, logAverage+10.0);
  } else {
    function->SetParameter(offset, std::max(1.0e-9, average));
    function->SetParLimits(
        offset, 1.0e-12, std::max(10.0, 100.0*maximum));
  }

  for (int i = 1; i <= 3; ++i) {
    function->SetParName(offset+i, Form("B_{%d}", i));
    if (i > order) {
      function->FixParameter(offset+i, 0.0);
      continue;
    }
    function->SetParameter(offset+i, isBernstein ? 1.0 : 0.0);
    if (isBernstein) {
      function->SetParLimits(offset+i, 1.0e-4, 100.0);
    } else if (isChebyshev) {
      function->SetParLimits(offset+i, -1.0, 1.0);
    } else if (isExponential) {
      function->SetParLimits(offset+i, -10.0, 10.0);
    } else {
      function->SetParLimits(offset+i, -20.0, 20.0);
    }
  }
}

struct Fit2DResult {
  bool ok = false;

  double mean  = 1.0194;
  double meanErr = 0.0;
  double gamma = 0.0042; // fitted Breit-Wigner width
  double gammaErr = 0.0;

  double nSS = 0.0;
  double nSB = 0.0;
  double nBS = 0.0;
  double nBB = 0.0;

  double nSSErr = 0.0;
  double nSBErr = 0.0;
  double nBSErr = 0.0;
  double nBBErr = 0.0;

  double c1m1 = 0.0, c2m1 = 0.0, c3m1 = 0.0;
  double c1m2 = 0.0, c2m2 = 0.0, c3m2 = 0.0;

  int chebOrder = 2;
  double phiMassMin = 1.000;
  double phiMassMax = 1.040;
};

struct WindowResult {
  double nData = 0.0;
  double nDataErr = 0.0;
  double nCont = 0.0;
  double nExcess = 0.0;
  double ratio = 0.0;
  double fracExcess = 0.0;
  double significance = 0.0; // nExcess / nDataErr, continuum uncertainty not included
};

void SetNiceHist(TH1* h) {
  if (!h) return;
  h->SetStats(false);
  h->SetLineWidth(2);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.9);
}


void StyleAxis1D(TH1* h)
{
  if (!h) return;
  h->SetStats(false);
  h->GetXaxis()->SetTitleSize(0.050);
  h->GetYaxis()->SetTitleSize(0.050);
  h->GetXaxis()->SetLabelSize(0.043);
  h->GetYaxis()->SetLabelSize(0.043);
  h->GetXaxis()->SetTitleOffset(1.05);
  h->GetYaxis()->SetTitleOffset(1.25);
  h->GetXaxis()->SetNdivisions(505);
  h->GetYaxis()->SetNdivisions(505);
}

void StyleAxis2D(TH2* h)
{
  if (!h) return;
  h->SetStats(false);
  h->GetXaxis()->SetTitleSize(0.047);
  h->GetYaxis()->SetTitleSize(0.047);
  h->GetZaxis()->SetTitleSize(0.043);
  h->GetXaxis()->SetLabelSize(0.040);
  h->GetYaxis()->SetLabelSize(0.040);
  h->GetZaxis()->SetLabelSize(0.037);
  h->GetXaxis()->SetTitleOffset(1.08);
  h->GetYaxis()->SetTitleOffset(1.23);
  h->GetZaxis()->SetTitleOffset(1.18);
  h->GetXaxis()->SetNdivisions(505);
  h->GetYaxis()->SetNdivisions(505);
}

double HistogramMinimum(const TH1* h)
{
  if (!h || h->GetNbinsX() <= 0) return 0.0;
  double ymin = h->GetBinContent(1);
  for (int ib = 2; ib <= h->GetNbinsX(); ++ib) {
    ymin = std::min(ymin, h->GetBinContent(ib));
  }
  return ymin;
}

double HistogramMaximum(const TH1* h)
{
  if (!h || h->GetNbinsX() <= 0) return 0.0;
  double ymax = h->GetBinContent(1);
  for (int ib = 2; ib <= h->GetNbinsX(); ++ib) {
    ymax = std::max(ymax, h->GetBinContent(ib));
  }
  return ymax;
}

void SetCountDisplayRange(TH1* h, const TH1* secondHistogram = nullptr)
{
  if (!h) return;
  double ymin = HistogramMinimum(h);
  double ymax = HistogramMaximum(h);
  if (secondHistogram) {
    ymin = std::min(ymin, HistogramMinimum(secondHistogram));
    ymax = std::max(ymax, HistogramMaximum(secondHistogram));
  }

  // Requested plotting convention for count spectra.
  h->SetMinimum(0.7 * ymin);
  h->SetMaximum(1.3 * ymax);
}

void SetCountDisplayRangeMany(
    TH1* h,
    std::initializer_list<const TH1*> additionalHistograms)
{
  if (!h) return;
  double ymin = HistogramMinimum(h);
  double ymax = HistogramMaximum(h);
  for (const TH1* other : additionalHistograms) {
    if (!other) continue;
    ymin = std::min(ymin, HistogramMinimum(other));
    ymax = std::max(ymax, HistogramMaximum(other));
  }
  h->SetMinimum(0.7 * ymin);
  h->SetMaximum(1.3 * ymax);
}

void SetSymmetricDisplayRange(TH1* h, double headroom = 1.22)
{
  if (!h) return;
  double maxAbs = 0.0;
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
    const double y = h->GetBinContent(ib);
    const double e = h->GetBinError(ib);
    maxAbs = std::max(maxAbs, std::abs(y) + e);
  }
  maxAbs = std::max(1.0, maxAbs) * headroom;
  h->SetMinimum(-maxAbs);
  h->SetMaximum(maxAbs);
}

double ContinuumValueNoReject(double x, const double* p)
{
  return BackgroundValue(gContModel, x, gContXMin, gContXMax, p);
}

double ContinuumNoReject(double* x, double* p)
{
  return ContinuumValueNoReject(x[0], p);
}

double ContinuumReject(double* x, double* p)
{
  // Used only during the continuum fit.  The signal/exotic window is rejected
  // from the fit, but the final continuum must still be evaluated in that
  // window as an interpolation.
  if (gRejectMax > gRejectMin && x[0] > gRejectMin && x[0] < gRejectMax) {
    TF1::RejectPoint();
    return 0.0;
  }
  return ContinuumValueNoReject(x[0], p);
}

const char* ContinuumModelLabel()
{
  return BackgroundLabel(gContModel);
}

int ExactNBinsFromWidth(double xmin,
                         double xmax,
                         double requestedWidth,
                         const char* label)
{
  if (!(xmax > xmin)) {
    std::cerr << "ERROR: invalid range for " << label
              << ": [" << xmin << ", " << xmax << "]" << std::endl;
    return -1;
  }
  if (!(requestedWidth > 0.0) || !std::isfinite(requestedWidth)) {
    std::cerr << "ERROR: " << label
              << " must be positive and finite, got "
              << requestedWidth << std::endl;
    return -1;
  }

  const double exact = (xmax - xmin) / requestedWidth;
  const long long rounded = std::llround(exact);
  const double tolerance = 1.0e-8 * std::max(1.0, std::abs(exact));

  if (rounded < 1 || std::abs(exact - static_cast<double>(rounded)) > tolerance) {
    std::cerr << "ERROR: Mpair range [" << xmin << ", " << xmax
              << "] is not an integer multiple of " << label
              << " = " << requestedWidth << ".  Exact number of bins would be "
              << exact << std::endl;
    return -1;
  }

  return static_cast<int>(rounded);
}

double PhiPuritySignalModel(double* x, double* p)
{
  return std::max(0.0, p[0]) * TMath::BreitWigner(
      x[0], p[1], std::max(1.0e-9, p[2]));
}

double PhiPurityBackgroundModel(double* x, double* p)
{
  return BackgroundValue(
      gPhiBackgroundModel, x[0], gPhiFitMin, gPhiFitMax, p);
}

double PhiPurityTotalModel(double* x, double* p)
{
  return PhiPuritySignalModel(x, p)
       + PhiPurityBackgroundModel(x, p+3);
}

void PhiPurityWindowYields(const TH1D* hMass,
                           const double* parameters,
                           double windowMin,
                           double windowMax,
                           double& signal,
                           double& background)
{
  signal = 0.0;
  background = 0.0;
  if (!hMass || !parameters || !(windowMax > windowMin)) return;

  for (int ib = 1; ib <= hMass->GetNbinsX(); ++ib) {
    const double binLo = hMass->GetXaxis()->GetBinLowEdge(ib);
    const double binHi = hMass->GetXaxis()->GetBinUpEdge(ib);
    const double overlapLo = std::max(binLo, windowMin);
    const double overlapHi = std::min(binHi, windowMax);
    if (!(overlapHi > overlapLo)) continue;
    const double fraction = (overlapHi-overlapLo)/(binHi-binLo);
    double x[1] = {hMass->GetXaxis()->GetBinCenter(ib)};
    signal += fraction*PhiPuritySignalModel(x, const_cast<double*>(parameters));
    background += fraction*PhiPurityBackgroundModel(
        x, const_cast<double*>(parameters+3));
  }
}

struct PhiPurityVsPtResult {
  std::unique_ptr<TGraphAsymmErrors> purityGraph;
  int successfulFits = 0;
};

PhiPurityVsPtResult DrawPhiPurityVsPt(
    TH2* hMassVsPt,
    double purityHalfWidth,
    int phiBkgModel,
    double requestedFitMin,
    double requestedFitMax,
    const char* outPrefix)
{
  PhiPurityVsPtResult out;
  if (!hMassVsPt || !(purityHalfWidth > 0.0) || !outPrefix ||
      !IsValidBackgroundModel(phiBkgModel) ||
      !(requestedFitMax > requestedFitMin)) return out;

  const std::vector<double> ptEdges = {
      0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0,
      4.0, 5.0, 6.0, 8.0, 10.0, 20.0, 50.0};
  const int nPtBins = static_cast<int>(ptEdges.size())-1;
  const double fitMin = std::max(
      requestedFitMin, hMassVsPt->GetXaxis()->GetXmin());
  const double fitMax = std::min(
      requestedFitMax, hMassVsPt->GetXaxis()->GetXmax());
  if (!(fitMax > fitMin)) {
    std::cerr << "WARNING: hPhiMassVsPt mass axis does not overlap "
              << "the requested range [" << requestedFitMin << ","
              << requestedFitMax << "] GeV/c^2." << std::endl;
    return out;
  }
  gPhiBackgroundModel = phiBkgModel;
  gPhiFitMin = fitMin;
  gPhiFitMax = fitMax;

  out.purityGraph.reset(new TGraphAsymmErrors());
  out.purityGraph->SetName("gPhiPurity_vsPt");

  std::ofstream csv(Form("%s_phiPurity_vsPt.csv", outPrefix));
  csv << "pt_low,pt_high,fit_status,entries,mean,mean_error,"
         "gamma,gamma_error,purity_window_low,purity_window_high,"
         "signal_in_window,background_in_window,purity,purity_error\n";

  auto* fitCanvas = new TCanvas(
      "cPhiPurityVsPtFits", "single-phi mass fits versus pT", 1900, 1500);
  fitCanvas->Divide(4, 4, 0.002, 0.002);

  std::vector<std::unique_ptr<TH1D>> massProjections;
  std::vector<std::unique_ptr<TF1>> totalFunctions;
  std::vector<std::unique_ptr<TF1>> signalFunctions;
  std::vector<std::unique_ptr<TF1>> backgroundFunctions;
  massProjections.reserve(nPtBins);
  totalFunctions.reserve(nPtBins);
  signalFunctions.reserve(nPtBins);
  backgroundFunctions.reserve(nPtBins);

  int graphPoint = 0;
  for (int ipt = 0; ipt < nPtBins; ++ipt) {
    const double ptLo = ptEdges[ipt];
    const double ptHi = ptEdges[ipt+1];
    const int yFirst = std::max(1, std::min(
        hMassVsPt->GetNbinsY(),
        hMassVsPt->GetYaxis()->FindBin(ptLo+1.0e-4)));
    const int yLast = std::max(1, std::min(
        hMassVsPt->GetNbinsY(),
        hMassVsPt->GetYaxis()->FindBin(ptHi-1.0e-4)));

    std::unique_ptr<TH1D> hMass(hMassVsPt->ProjectionX(
        Form("hPhiMass_ptBin_%02d_%.1f_%.1f", ipt, ptLo, ptHi),
        yFirst, yLast, "e"));
    if (!hMass) {
      csv << ptLo << "," << ptHi << ",-999,0,0,0,0,0,0,0,0,0,0,0\n";
      continue;
    }
    hMass->SetDirectory(nullptr);
    const int fitBinLo = hMass->FindBin(fitMin+1.0e-9);
    const int fitBinHi = hMass->FindBin(fitMax-1.0e-9);
    const double entries = hMass->Integral(fitBinLo, fitBinHi);
    if (!(entries > 0.0)) {
      csv << ptLo << "," << ptHi << ",-998,0,0,0,0,0,0,0,0,0,0,0\n";
      continue;
    }

    const double average = std::max(
        1.0e-6, entries/std::max(1, fitBinHi-fitBinLo+1));
    const double peak = HistogramMaximum(hMass.get());
    const double meanInit = 1.0194;
    const double gammaInit = 0.0042;
    double xAtPeak[1] = {meanInit};
    double unitSignal[3] = {1.0, meanInit, gammaInit};
    const double unitPeak = std::max(
        1.0e-12, PhiPuritySignalModel(xAtPeak, unitSignal));
    const double signalNormInit = std::max(
        1.0e-9, (peak-average)/unitPeak);

    std::unique_ptr<TF1> total(new TF1(
        Form("fPhiMassTotal_ptBin_%02d", ipt),
        PhiPurityTotalModel, fitMin, fitMax, 7));
    total->SetNpx(1000);
    total->SetParameter(0, signalNormInit);
    total->SetParameter(1, meanInit);
    total->SetParameter(2, gammaInit);
    total->SetParName(0, "N_{BW}");
    total->SetParName(1, "m_{#phi}");
    total->SetParName(2, "#Gamma_{#phi}");
    total->SetParLimits(0, 0.0, std::max(1.0, 20.0*entries));
    total->SetParLimits(1, 1.0160, 1.0220);
    total->SetParLimits(2, 0.0010, 0.0150);
    ConfigureBackgroundParameters(
        total.get(), 3, phiBkgModel, average, peak);

    // N keeps ownership of the TF1 here instead of attaching it to TH1.
    TFitResultPtr fitResult = hMass->Fit(total.get(), "ERQNBS");


auto corr = fitResult->GetCorrelationMatrix();

std::cout << "\nCorrelation matrix:\n\n";

std::cout << std::setw(15) << "";
for (int j = 0; j < fitResult->NPar(); ++j)
    std::cout << std::setw(15) << fitResult->ParName(j);
std::cout << "\n";

for (int i = 0; i < fitResult->NPar(); ++i) {
    std::cout << std::setw(15) << fitResult->ParName(i);
    for (int j = 0; j < fitResult->NPar(); ++j)
        std::cout << std::setw(15) << std::fixed << std::setprecision(3)
                  << corr(i,j);
    std::cout << "\n";
}


    int fitStatus = fitResult.Get() ? fitResult->Status()
                                    : static_cast<int>(fitResult);
    if (fitStatus != 0) {
      fitResult = hMass->Fit(total.get(), "RQMN0S");
      fitStatus = fitResult.Get() ? fitResult->Status()
                                  : static_cast<int>(fitResult);
    }

    double parameters[7] = {0.0};
    for (int ip = 0; ip < 7; ++ip) parameters[ip] = total->GetParameter(ip);
    const double windowLo = std::max(fitMin,
        parameters[1]-purityHalfWidth);
    const double windowHi = std::min(fitMax,
        parameters[1]+purityHalfWidth);
    double signal = 0.0;
    double background = 0.0;
    PhiPurityWindowYields(
        hMass.get(), parameters, windowLo, windowHi, signal, background);
    const double denominator = signal+background;
    const double purity = denominator > 0.0 ? signal/denominator : 0.0;

    double purityError = 0.0;
    if (fitResult.Get() && fitResult->CovMatrixStatus() > 0 &&
        denominator > 0.0) {
      const TMatrixDSym covariance = fitResult->GetCovarianceMatrix();
      double gradient[7] = {0.0};
      for (int ip = 0; ip < 7; ++ip) {
        double step = std::max(
            1.0e-7,
            std::max(std::abs(parameters[ip])*1.0e-5,
                     0.1*std::abs(total->GetParError(ip))));
        if ((ip == 0 || ip == 2) && parameters[ip] > 0.0) {
          step = std::min(step, 0.45*parameters[ip]);
        }
        double plus[7] = {0.0};
        double minus[7] = {0.0};
        for (int jp = 0; jp < 7; ++jp) {
          plus[jp] = parameters[jp];
          minus[jp] = parameters[jp];
        }
        plus[ip] += step;
        minus[ip] -= step;
        double sPlus = 0.0, bPlus = 0.0;
        double sMinus = 0.0, bMinus = 0.0;
        const double plusWindowLo = std::max(
            fitMin, plus[1]-purityHalfWidth);
        const double plusWindowHi = std::min(
            fitMax, plus[1]+purityHalfWidth);
        const double minusWindowLo = std::max(
            fitMin, minus[1]-purityHalfWidth);
        const double minusWindowHi = std::min(
            fitMax, minus[1]+purityHalfWidth);
        PhiPurityWindowYields(
            hMass.get(), plus, plusWindowLo, plusWindowHi, sPlus, bPlus);
        PhiPurityWindowYields(
            hMass.get(), minus, minusWindowLo, minusWindowHi,
            sMinus, bMinus);
        const double pPlus = (sPlus+bPlus) > 0.0
            ? sPlus/(sPlus+bPlus) : purity;
        const double pMinus = (sMinus+bMinus) > 0.0
            ? sMinus/(sMinus+bMinus) : purity;
        gradient[ip] = (pPlus-pMinus)/(2.0*step);
      }
      double variance = 0.0;
      for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
          variance += gradient[i]*covariance(i,j)*gradient[j];
        }
      }
      purityError = std::sqrt(std::max(0.0, variance));
    }

    if (fitStatus == 0 && std::isfinite(purity) &&
        std::isfinite(purityError)) {
      const double ptCentre = std::sqrt(ptLo*ptHi);
      out.purityGraph->SetPoint(graphPoint, ptCentre, purity);
      out.purityGraph->SetPointError(
          graphPoint, ptCentre-ptLo, ptHi-ptCentre,
          purityError, purityError);
      ++graphPoint;
      ++out.successfulFits;
    }

    std::unique_ptr<TF1> signalFunction(new TF1(
        Form("fPhiMassSignal_ptBin_%02d", ipt),
        PhiPuritySignalModel, fitMin, fitMax, 3));
    signalFunction->SetParameters(parameters[0], parameters[1], parameters[2]);
    signalFunction->SetLineColor(kMagenta+1);
    signalFunction->SetLineStyle(2);
    signalFunction->SetLineWidth(3);
    std::unique_ptr<TF1> backgroundFunction(new TF1(
        Form("fPhiMassBackground_ptBin_%02d", ipt),
        PhiPurityBackgroundModel, fitMin, fitMax, 4));
    backgroundFunction->SetParameters(
        parameters[3], parameters[4], parameters[5], parameters[6]);
    backgroundFunction->SetLineColor(kBlue+1);
    backgroundFunction->SetLineStyle(7);
    backgroundFunction->SetLineWidth(3);
    total->SetLineColor(kRed+1);
    total->SetLineWidth(3);

    fitCanvas->cd(ipt+1);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.03);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.10);
    gPad->SetTicks(1, 1);
    hMass->SetTitle(Form(
        "%.1f < p_{T}^{#phi} < %.1f GeV/#it{c};"
        "M_{K^{+}K^{-}} (GeV/#it{c}^{2});counts", ptLo, ptHi));
    hMass->SetMarkerStyle(20);
    hMass->SetMarkerSize(0.45);
    hMass->SetLineColor(kBlack);
    hMass->GetXaxis()->SetRangeUser(fitMin, fitMax);
    hMass->SetMinimum(0.7*HistogramMinimum(hMass.get()));
    hMass->SetMaximum(1.3*std::max(
        HistogramMaximum(hMass.get()), total->GetMaximum(fitMin,fitMax)));
    StyleAxis1D(hMass.get());
    hMass->Draw("E1");
    total->Draw("same");
    backgroundFunction->Draw("same");
    signalFunction->Draw("same");
    hMass->Draw("E1 same");

    TLatex fitText;
    fitText.SetNDC();
    fitText.SetTextSize(0.040);
    fitText.DrawLatex(0.17, 0.84,
        Form("P_{#phi}=%.3f #pm %.3f", purity, purityError));
    fitText.DrawLatex(0.17, 0.78,
        Form("m_{#phi}=%.5f, #Gamma=%.4f", parameters[1], parameters[2]));
    fitText.DrawLatex(0.17, 0.72, Form("fit status=%d", fitStatus));
    if (ipt == 0) {
      auto* fitLegend = new TLegend(0.54, 0.57, 0.95, 0.86);
      fitLegend->SetBorderSize(0);
      fitLegend->SetFillStyle(0);
      fitLegend->SetTextSize(0.037);
      fitLegend->AddEntry(hMass.get(), "data", "lep");
      fitLegend->AddEntry(
          total.get(), Form("BW + %s", BackgroundLabel(phiBkgModel)), "l");
      fitLegend->AddEntry(signalFunction.get(), "Breit-Wigner", "l");
      fitLegend->AddEntry(
          backgroundFunction.get(),
          Form("%s background", BackgroundLabel(phiBkgModel)), "l");
      fitLegend->Draw();
    }

    csv << ptLo << "," << ptHi << "," << fitStatus << "," << entries
        << "," << parameters[1] << "," << total->GetParError(1)
        << "," << parameters[2] << "," << total->GetParError(2)
        << "," << windowLo << "," << windowHi
        << "," << signal << "," << background
        << "," << purity << "," << purityError << "\n";

    massProjections.push_back(std::move(hMass));
    totalFunctions.push_back(std::move(total));
    signalFunctions.push_back(std::move(signalFunction));
    backgroundFunctions.push_back(std::move(backgroundFunction));
  }
  csv.close();

  fitCanvas->SaveAs(Form("%s_phiPurity_vsPt_fits.png", outPrefix));
  fitCanvas->SaveAs(Form("%s_phiPurity_vsPt_fits.pdf", outPrefix));
  delete fitCanvas;

  if (out.purityGraph->GetN() > 0) {
    auto* purityCanvas = new TCanvas(
        "cPhiPurityVsPt", "single-phi purity versus pT", 1000, 780);
    purityCanvas->SetLeftMargin(0.14);
    purityCanvas->SetRightMargin(0.04);
    purityCanvas->SetBottomMargin(0.13);
    purityCanvas->SetTopMargin(0.08);
    purityCanvas->SetTicks(1, 1);
    purityCanvas->SetGridy();
    purityCanvas->SetLogx();

    out.purityGraph->SetTitle(
        "Single-#phi purity versus transverse momentum;"
        "p_{T}^{#phi} (GeV/#it{c});P_{#phi}=S/(S+B)");
    out.purityGraph->SetLineColor(kBlue+1);
    out.purityGraph->SetMarkerColor(kBlue+1);
    out.purityGraph->SetLineWidth(3);
    out.purityGraph->SetMarkerStyle(20);
    out.purityGraph->SetMarkerSize(1.1);
    out.purityGraph->GetXaxis()->SetLimits(0.45, 55.0);
    out.purityGraph->SetMinimum(0.0);
    out.purityGraph->SetMaximum(1.05);
    out.purityGraph->GetXaxis()->SetTitleSize(0.050);
    out.purityGraph->GetYaxis()->SetTitleSize(0.050);
    out.purityGraph->GetXaxis()->SetLabelSize(0.043);
    out.purityGraph->GetYaxis()->SetLabelSize(0.043);
    out.purityGraph->GetYaxis()->SetTitleOffset(1.25);
    out.purityGraph->Draw("APZ");

    TLatex purityText;
    purityText.SetNDC();
    purityText.SetTextSize(0.034);
    purityText.DrawLatex(
        0.18, 0.88,
        Form("mass model: Breit-Wigner + %s", BackgroundLabel(phiBkgModel)));
    purityText.DrawLatex(0.18, 0.83,
        Form("purity window: |M_{K^{+}K^{-}}-m_{#phi}^{fit}| < %.3f GeV/#it{c}^{2}",
             purityHalfWidth));

    purityCanvas->SaveAs(Form("%s_Fig8_phiPurity_vsPt.png", outPrefix));
    purityCanvas->SaveAs(Form("%s_Fig8_phiPurity_vsPt.pdf", outPrefix));
    delete purityCanvas;
  }

  return out;
}


TF1* MakeContinuumInterpolationFunction(TF1* fFitReject,
                                        const char* name,
                                        double xmin,
                                        double xmax)
{
  if (!fFitReject) return nullptr;
  auto* f = new TF1(name, ContinuumNoReject, xmin, xmax, gContNPar);
  f->SetNpx(1000);
  for (int i = 0; i < gContNPar; ++i) {
    f->SetParameter(i, fFitReject->GetParameter(i));
    f->SetParError(i, fFitReject->GetParError(i));
  }
  f->SetLineColor(kRed);
  f->SetLineWidth(3);
  return f;
}

// Save and restore every THnSparse axis range.  Each native projection below
// starts from a clean (unrestricted) sparse, applies only its own cuts, and
// leaves the caller's axis ranges unchanged.
class SparseAxisRangeGuard {
 public:
  explicit SparseAxisRangeGuard(THnSparseF* h) : fSparse(h) {
    if (!fSparse) return;
    fRanges.reserve(fSparse->GetNdimensions());
    for (int i = 0; i < fSparse->GetNdimensions(); ++i) {
      TAxis* ax = fSparse->GetAxis(i);
      fRanges.push_back(
          {ax->GetFirst(), ax->GetLast(), ax->TestBit(TAxis::kAxisRange)});
      ax->SetRange(0, 0);
    }
  }

  ~SparseAxisRangeGuard() {
    if (!fSparse) return;
    for (int i = 0; i < fSparse->GetNdimensions(); ++i) {
      if (fRanges[i].active) {
        fSparse->GetAxis(i)->SetRange(fRanges[i].first, fRanges[i].last);
      } else {
        fSparse->GetAxis(i)->SetRange(0, 0);
      }
    }
  }

 private:
  struct SavedRange {
    int first = 0;
    int last = 0;
    bool active = false;
  };
  THnSparseF* fSparse = nullptr;
  std::vector<SavedRange> fRanges;
};

// Apply edge-safe THnSparse ranges.  The +/-1e-4 offsets ensure that a cut
// placed exactly on a bin edge selects the intended first and last bins.
bool SetAxisRangeByValue(TAxis* ax, double low, double high)
{
  if (!ax || !(high > low)) return false;
  const double axisMin = ax->GetXmin();
  const double axisMax = ax->GetXmax();
  const double selectedLow = std::max(low, axisMin);
  const double selectedHigh = std::min(high, axisMax);
  if (!(selectedHigh > selectedLow)) return false;
  const int first = std::max(1, std::min(
      ax->GetNbins(), ax->FindBin(selectedLow + 1.0e-4)));
  const int last = std::max(1, std::min(
      ax->GetNbins(), ax->FindBin(selectedHigh - 1.0e-4)));
  if (first > last) return false;
  ax->SetRange(first, last);
  return true;
}

bool ApplyRapiditySelection(THnSparseF* hSparse)
{
  if (!hSparse || hSparse->GetNdimensions() <= kAxisRapidity) return false;
  return SetAxisRangeByValue(
      hSparse->GetAxis(kAxisRapidity), gRapidityMin, gRapidityMax);
}

TH1D* RebinNativePairProjection(std::unique_ptr<TH1D> projected,
                                const char* name,
                                int requestedBins,
                                double requestedMin,
                                double requestedMax)
{
  if (!projected || requestedBins <= 0) return nullptr;
  projected->SetDirectory(nullptr);

  const int sourceBins = projected->GetNbinsX();
  const double sourceMin = projected->GetXaxis()->GetBinLowEdge(1);
  const double sourceMax = projected->GetXaxis()->GetBinUpEdge(sourceBins);
  const double edgeTolerance =
      1.0e-9 * std::max({1.0, std::abs(requestedMin), std::abs(requestedMax)});
  if (std::abs(sourceMin - requestedMin) > edgeTolerance ||
      std::abs(sourceMax - requestedMax) > edgeTolerance) {
    std::cerr << "ERROR: requested M(phi phi) range [" << requestedMin << ", "
              << requestedMax << "] does not coincide with THnSparse axis edges. "
              << "Native projection gives [" << sourceMin << ", " << sourceMax
              << "].\n";
    return nullptr;
  }
  if (sourceBins == requestedBins) {
    projected->SetName(name);
    return projected.release();
  }
  if (sourceBins < requestedBins || sourceBins % requestedBins != 0) {
    std::cerr << "ERROR: native M(phi phi) projection has " << sourceBins
              << " bins in the selected range; it cannot be rebinned exactly to "
              << requestedBins << " bins.  Change finalPairBinWidth to an integer "
              << "multiple of the THnSparse axis bin width.\n";
    return nullptr;
  }

  const int group = sourceBins / requestedBins;
  TH1D* rebinned = dynamic_cast<TH1D*>(projected->Rebin(group, name));
  if (!rebinned) return nullptr;
  rebinned->SetDirectory(nullptr);
  if (rebinned == projected.get()) return projected.release();
  return rebinned;
}

TH1D* ProjectPairMassNativeRanges(THnSparseF* hSparse,
                                  const char* name,
                                  double ptMin,
                                  double ptMax,
                                  double m1Min,
                                  double m1Max,
                                  double m2Min,
                                  double m2Max,
                                  bool applyDMRange,
                                  double dmMin,
                                  double dmMax,
                                  double mPairMin,
                                  double mPairMax,
                                  int nPairBins)
{
  if (!hSparse) return nullptr;
  SparseAxisRangeGuard guard(hSparse);

  if (!ApplyRapiditySelection(hSparse) ||
      !SetAxisRangeByValue(hSparse->GetAxis(kAxisMPair), mPairMin, mPairMax) ||
      !SetAxisRangeByValue(hSparse->GetAxis(kAxisPtPair), ptMin, ptMax) ||
      !SetAxisRangeByValue(hSparse->GetAxis(kAxisM1), m1Min, m1Max) ||
      !SetAxisRangeByValue(hSparse->GetAxis(kAxisM2), m2Min, m2Max) ||
      (applyDMRange &&
       !SetAxisRangeByValue(hSparse->GetAxis(kAxisDeltaM), dmMin, dmMax))) {
    return nullptr;
  }

  // ROOT performs the projection and propagates the THnSparse Sumw2 errors.
  std::unique_ptr<TH1D> projected(hSparse->Projection(kAxisMPair, "E"));
  TH1D* result = RebinNativePairProjection(
      std::move(projected), name, nPairBins, mPairMin, mPairMax);
  if (result) {
    result->SetTitle(";M_{#phi#phi} (GeV/#it{c}^{2});candidate pairs");
  }
  return result;
}

TH3D* BuildM1M2DeltaMFromSparse(THnSparseF* hSparse,
                                const char* name,
                                double ptMin,
                                double ptMax,
                                double mPairMin,
                                double mPairMax,
                                double phiMassMin,
                                double phiMassMax,
                                double dmMax = 0.050)
{
  if (!hSparse) return nullptr;

  const int nDim = hSparse->GetNdimensions();
  if (nDim <= std::max({kAxisMPair, kAxisPtPair, kAxisRapidity,
                        kAxisM1, kAxisM2, kAxisDeltaM})) {
    std::cerr << "ERROR: THnSparse has only " << nDim
              << " dimensions, but axes 0,1,3,4,5,6 are requested.\n";
    return nullptr;
  }

  SparseAxisRangeGuard guard(hSparse);
  if (!ApplyRapiditySelection(hSparse) ||
      !SetAxisRangeByValue(hSparse->GetAxis(kAxisMPair), mPairMin, mPairMax) ||
      !SetAxisRangeByValue(hSparse->GetAxis(kAxisPtPair), ptMin, ptMax) ||
      !SetAxisRangeByValue(hSparse->GetAxis(kAxisM1), phiMassMin, phiMassMax) ||
      !SetAxisRangeByValue(hSparse->GetAxis(kAxisM2), phiMassMin, phiMassMax) ||
      !SetAxisRangeByValue(hSparse->GetAxis(kAxisDeltaM), 0.0, dmMax)) {
    return nullptr;
  }

  // For the 3D overload the ROOT order is Projection(xDim,yDim,zDim).
  TH3D* h3 = hSparse->Projection(kAxisM1, kAxisM2, kAxisDeltaM, "E");
  if (!h3) return nullptr;
  h3->SetName(name);
  h3->SetDirectory(nullptr);
  h3->SetTitle(";M_{#phi,1} (GeV/#it{c}^{2});M_{#phi,2} (GeV/#it{c}^{2});stored #DeltaM_{#phi} (GeV/#it{c}^{2})");
  std::cout << "Native THnSparse 3D projection " << name
            << ": integral = " << h3->Integral() << std::endl;
  return h3;
}

RooFitResult* MinimizeExtendedNLL(RooAbsPdf& model, RooDataHist& data) {
  std::unique_ptr<RooAbsReal> nll(model.createNLL(data, Extended(true), Offset(true)));
  RooMinimizer minim(*nll);
  minim.setPrintLevel(-1);
  minim.setStrategy(1);
  minim.setEps(1e-6);
  minim.optimizeConst(2);

  int status = minim.minimize("Minuit", "Migrad");
  if (status != 0) {
    std::cout << "First Migrad status = " << status << ", retry strategy 2\n";
    minim.setStrategy(2);
    status = minim.minimize("Minuit", "Migrad");
  }
  minim.hesse();
  return minim.save("fit2DResult", "fit2DResult");
}

void FillChebCoeffList(RooArgList& list,
                       std::vector<std::unique_ptr<RooRealVar>>& store,
                       const char* prefix,
                       int order)
{
  const int ord = std::max(1, std::min(order, 3));
  for (int i = 1; i <= ord; ++i) {
    store.emplace_back(new RooRealVar(Form("%s_c%d", prefix, i),
                                      Form("%s_c%d", prefix, i), 0.0, -1.0, 1.0));
    list.add(*store.back());
  }
}

double ChebPdf(double x, double xmin, double xmax, double c1, double c2, double c3, int order);
double BWPdf(double x, double xmin, double xmax, double mean, double gamma);

Fit2DResult FitFull2DMass(TH2D* h2,
                          double phiMassMin,
                          double phiMassMax,
                          int chebOrder,
                          const char* outPrefix)
{
  Fit2DResult res;
  res.phiMassMin = phiMassMin;
  res.phiMassMax = phiMassMax;
  res.chebOrder = std::max(1, std::min(chebOrder, 3));

  if (!h2 || h2->Integral() <= 0.0) {
    std::cerr << "ERROR: empty h2 passed to FitFull2DMass.\n";
    return res;
  }

  RooRealVar m1("m1", "M_{#phi,1}", phiMassMin, phiMassMax, "GeV/c^{2}");
  RooRealVar m2("m2", "M_{#phi,2}", phiMassMin, phiMassMax, "GeV/c^{2}");
  RooDataHist data2D("data2D", "data2D", RooArgList(m1, m2), h2);

  RooRealVar mean("mean", "m_{#phi}", 1.01945, 1.0170, 1.0215);
  RooRealVar gamma("gamma", "#Gamma_{#phi}", gPhiBWGammaInit,
                   gPhiBWGammaMin, gPhiBWGammaMax);

  RooBreitWigner sig1("sig1", "BW(M1)", m1, mean, gamma);
  RooBreitWigner sig2("sig2", "BW(M2)", m2, mean, gamma);

  RooArgList coeffs1;
  RooArgList coeffs2;
  std::vector<std::unique_ptr<RooRealVar>> coeffStore1;
  std::vector<std::unique_ptr<RooRealVar>> coeffStore2;
  FillChebCoeffList(coeffs1, coeffStore1, "m1", res.chebOrder);
  FillChebCoeffList(coeffs2, coeffStore2, "m2", res.chebOrder);

  RooChebychev bkg1("bkg1", "B(M1)", m1, coeffs1);
  RooChebychev bkg2("bkg2", "B(M2)", m2, coeffs2);

  RooProdPdf pdfSS("pdfSS", "S(M1)S(M2)", RooArgSet(sig1, sig2));
  RooProdPdf pdfSB("pdfSB", "S(M1)B(M2)", RooArgSet(sig1, bkg2));
  RooProdPdf pdfBS("pdfBS", "B(M1)S(M2)", RooArgSet(bkg1, sig2));
  RooProdPdf pdfBB("pdfBB", "B(M1)B(M2)", RooArgSet(bkg1, bkg2));

  const double nTot = h2->Integral(1, h2->GetNbinsX(), 1, h2->GetNbinsY());

  RooRealVar nSS("nSS", "N_{SS}", 0.20 * nTot, 0.0, 2.0 * nTot);
  RooRealVar nSB("nSB", "N_{SB}", 0.25 * nTot, 0.0, 2.0 * nTot);
  RooRealVar nBS("nBS", "N_{BS}", 0.25 * nTot, 0.0, 2.0 * nTot);
  RooRealVar nBB("nBB", "N_{BB}", 0.30 * nTot, 0.0, 2.0 * nTot);

  RooAddPdf model2D("model2D", "SS+SB+BS+BB",
                    RooArgList(pdfSS, pdfSB, pdfBS, pdfBB),
                    RooArgList(nSS, nSB, nBS, nBB));

  std::unique_ptr<RooFitResult> fit(MinimizeExtendedNLL(model2D, data2D));

  res.ok = fit && fit->status() == 0;
  res.mean = mean.getVal();
  res.meanErr = mean.getError();
  res.gamma = gamma.getVal();
  res.gammaErr = gamma.getError();
  res.nSS = nSS.getVal();
  res.nSB = nSB.getVal();
  res.nBS = nBS.getVal();
  res.nBB = nBB.getVal();
  res.nSSErr = nSS.getError();
  res.nSBErr = nSB.getError();
  res.nBSErr = nBS.getError();
  res.nBBErr = nBB.getError();

  if (coeffStore1.size() > 0) res.c1m1 = coeffStore1[0]->getVal();
  if (coeffStore1.size() > 1) res.c2m1 = coeffStore1[1]->getVal();
  if (coeffStore1.size() > 2) res.c3m1 = coeffStore1[2]->getVal();
  if (coeffStore2.size() > 0) res.c1m2 = coeffStore2[0]->getVal();
  if (coeffStore2.size() > 1) res.c2m2 = coeffStore2[1]->getVal();
  if (coeffStore2.size() > 2) res.c3m2 = coeffStore2[2]->getVal();

  std::cout << "\n========== Full 2D phi-mass fit result ==========" << std::endl;
  std::cout << "fit status = " << (fit ? fit->status() : -999) << std::endl;
  std::cout << "mean                  = " << res.mean
            << " +/- " << res.meanErr << " GeV/c^2" << std::endl;
  std::cout << "Breit-Wigner Gamma_phi (full width) = "
            << res.gamma << " +/- " << res.gammaErr << " GeV/c^2 = "
            << 1000.0*res.gamma << " +/- " << 1000.0*res.gammaErr
            << " MeV/c^2" << std::endl;
  std::cout << "N_SS  = " << res.nSS << " +/- " << res.nSSErr << std::endl;
  std::cout << "N_SB  = " << res.nSB << " +/- " << res.nSBErr << std::endl;
  std::cout << "N_BS  = " << res.nBS << " +/- " << res.nBSErr << std::endl;
  std::cout << "N_BB  = " << res.nBB << " +/- " << res.nBBErr << std::endl;
  std::cout << "=================================================\n" << std::endl;

  RooPlot* fr1 = m1.frame(Bins(h2->GetNbinsX()));
  data2D.plotOn(fr1, Name("data_m1"));
  model2D.plotOn(fr1, Name("model_m1"), LineWidth(3));
  model2D.plotOn(fr1, Components(pdfSS), Name("SS_m1"), LineStyle(kDashDotted), LineWidth(2));
  model2D.plotOn(fr1, Components(RooArgSet(pdfSB, pdfBS)), Name("oneS_m1"), LineStyle(kDotted), LineWidth(2));
  model2D.plotOn(fr1, Components(pdfBB), Name("BB_m1"), LineStyle(kDashed), LineWidth(2));

  RooPlot* fr2 = m2.frame(Bins(h2->GetNbinsY()));
  data2D.plotOn(fr2, Name("data_m2"));
  model2D.plotOn(fr2, Name("model_m2"), LineWidth(3));
  model2D.plotOn(fr2, Components(pdfSS), Name("SS_m2"), LineStyle(kDashDotted), LineWidth(2));
  model2D.plotOn(fr2, Components(RooArgSet(pdfSB, pdfBS)), Name("oneS_m2"), LineStyle(kDotted), LineWidth(2));
  model2D.plotOn(fr2, Components(pdfBB), Name("BB_m2"), LineStyle(kDashed), LineWidth(2));

  auto* cproj = new TCanvas("c2Dfit_projection_checks", "c2Dfit_projection_checks", 1800, 800);
  cproj->Divide(2, 1);

  cproj->cd(1);
  gPad->SetLeftMargin(0.13);
  fr1->SetTitle("2D fit projection on M_{#phi,1}");
  fr1->GetYaxis()->SetTitle("Counts");
  fr1->Draw();
  auto* leg1 = new TLegend(0.55, 0.65, 0.88, 0.88);
  leg1->SetBorderSize(0); leg1->SetFillStyle(0);
  leg1->AddEntry(fr1->findObject("data_m1"), "Data", "lep");
  leg1->AddEntry(fr1->findObject("model_m1"), "Total fit", "l");
  leg1->AddEntry(fr1->findObject("SS_m1"), "SS", "l");
  leg1->AddEntry(fr1->findObject("oneS_m1"), "SB+BS", "l");
  leg1->AddEntry(fr1->findObject("BB_m1"), "BB", "l");
  leg1->Draw();

  cproj->cd(2);
  gPad->SetLeftMargin(0.13);
  fr2->SetTitle("2D fit projection on M_{#phi,2}");
  fr2->GetYaxis()->SetTitle("Counts");
  fr2->Draw();
  auto* leg2 = new TLegend(0.55, 0.65, 0.88, 0.88);
  leg2->SetBorderSize(0); leg2->SetFillStyle(0);
  leg2->AddEntry(fr2->findObject("data_m2"), "Data", "lep");
  leg2->AddEntry(fr2->findObject("model_m2"), "Total fit", "l");
  leg2->AddEntry(fr2->findObject("SS_m2"), "SS", "l");
  leg2->AddEntry(fr2->findObject("oneS_m2"), "SB+BS", "l");
  leg2->AddEntry(fr2->findObject("BB_m2"), "BB", "l");
  leg2->Draw();

  cproj->SaveAs(Form("%s_2Dfit_projection_checks.png", outPrefix));
  cproj->SaveAs(Form("%s_2Dfit_projection_checks.pdf", outPrefix));

  delete cproj;
  delete fr1;
  delete fr2;

  // Save the actual 2D invariant-mass distribution used by the fit.
  auto* c2dData = new TCanvas(
    "c2Dfit_data_global",
    "global 2D phi-mass data",
    900,
    800);
  c2dData->SetRightMargin(0.14);
  h2->SetStats(false);
  h2->Draw("COLZ");
  c2dData->SaveAs(Form("%s_2Dfit_data.png", outPrefix));
  c2dData->SaveAs(Form("%s_2Dfit_data.pdf", outPrefix));
  delete c2dData;

  // Save the numerical global-fit result next to the plots.
  std::ofstream fitSummary(Form("%s_2Dfit_result.csv", outPrefix));
  fitSummary
    << "status,mean,mean_err,phi_bw_width,phi_bw_width_err,"
    << "N_SS,N_SS_err,N_SB,N_SB_err,N_BS,N_BS_err,N_BB,N_BB_err,"
    << "c1_m1,c2_m1,c3_m1,c1_m2,c2_m2,c3_m2\n";
  fitSummary
    << (fit ? fit->status() : -999) << ","
    << res.mean << ","
    << res.meanErr << ","
    << res.gamma << ","
    << res.gammaErr << ","
    << res.nSS << "," << res.nSSErr << ","
    << res.nSB << "," << res.nSBErr << ","
    << res.nBS << "," << res.nBSErr << ","
    << res.nBB << "," << res.nBBErr << ","
    << res.c1m1 << "," << res.c2m1 << "," << res.c3m1 << ","
    << res.c1m2 << "," << res.c2m2 << "," << res.c3m2 << "\n";
  fitSummary.close();

  return res;
}

double ChebRaw(double x, double xmin, double xmax,
               double c1, double c2, double c3, int order)
{
  const double t = 2.0 * (x - xmin) / (xmax - xmin) - 1.0;
  double val = 1.0;
  if (order >= 1) val += c1 * t;
  if (order >= 2) val += c2 * (2.0 * t * t - 1.0);
  if (order >= 3) val += c3 * (4.0 * t * t * t - 3.0 * t);
  return std::max(0.0, val);
}

double ChebNormIntegral(double xmin, double xmax, double c2, int order)
{
  const double width = xmax - xmin;
  double integ = width;
  if (order >= 2) integ *= (1.0 - c2 / 3.0);
  if (integ <= 0.0 || !std::isfinite(integ)) integ = width;
  return integ;
}

double ChebPdf(double x, double xmin, double xmax,
               double c1, double c2, double c3, int order)
{
  return ChebRaw(x, xmin, xmax, c1, c2, c3, order) /
         ChebNormIntegral(xmin, xmax, c2, order);
}

double BWRaw(double x, double mean, double gamma)
{
  const double half = 0.5 * gamma;
  return 1.0 / ((x - mean) * (x - mean) + half * half);
}

double BWNormIntegral(double xmin, double xmax, double mean, double gamma)
{
  const double half = 0.5 * gamma;
  if (half <= 0.0) return 1.0;
  return (std::atan((xmax - mean) / half) -
          std::atan((xmin - mean) / half)) / half;
}

double BWPdf(double x, double xmin, double xmax, double mean, double gamma)
{
  const double integ = BWNormIntegral(xmin, xmax, mean, gamma);
  if (integ <= 0.0 || !std::isfinite(integ)) return 0.0;
  return BWRaw(x, mean, gamma) / integ;
}

double PosteriorSSProbability(double m1,
                              double m2,
                              const Fit2DResult& fit)
{
  if (!fit.ok) return -1.0;
  const double s1 = BWPdf(m1, fit.phiMassMin, fit.phiMassMax,
                          fit.mean, fit.gamma);
  const double s2 = BWPdf(m2, fit.phiMassMin, fit.phiMassMax,
                          fit.mean, fit.gamma);
  const double b1 = ChebPdf(m1, fit.phiMassMin, fit.phiMassMax,
                           fit.c1m1, fit.c2m1, fit.c3m1,
                           fit.chebOrder);
  const double b2 = ChebPdf(m2, fit.phiMassMin, fit.phiMassMax,
                           fit.c1m2, fit.c2m2, fit.c3m2,
                           fit.chebOrder);

  const double aSS = std::max(0.0, fit.nSS) * s1 * s2;
  const double total = aSS
      + std::max(0.0, fit.nSB) * s1 * b2
      + std::max(0.0, fit.nBS) * b1 * s2
      + std::max(0.0, fit.nBB) * b1 * b2;
  if (!(total > 0.0) || !std::isfinite(total)) return -1.0;
  return std::max(0.0, std::min(1.0, aSS / total));
}

struct DeltaMSignificanceScanResult {
  std::unique_ptr<TGraph> significanceGraph;
  std::unique_ptr<TGraph> purityGraph;
  std::unique_ptr<TGraph> purityWeightedSignificanceGraph;
  std::unique_ptr<TGraph> signalGraph;
  std::unique_ptr<TGraph> backgroundGraph;
  double bestCut = 0.0;
  double bestSignificance = 0.0;
  double bestPurityWeightedCut = 0.0;
  double bestPurityWeightedSignificance = 0.0;
  double significanceAtSelection = 0.0;
  double purityWeightedAtSelection = 0.0;
};

DeltaMSignificanceScanResult DrawDeltaMTruePairSignificance(
    TH2D* h2Mass,
    const Fit2DResult& fit,
    double scanMin,
    double scanMax,
    double scanStep,
    double selectionCut,
    double ptMin,
    double ptMax,
    double rapidityMin,
    double rapidityMax,
    double mPairMin,
    double mPairMax,
    const char* outPrefix)
{
  DeltaMSignificanceScanResult out;
  if (!h2Mass || !fit.ok || !(scanStep > 0.0) || !(scanMax >= scanMin)) {
    return out;
  }

  scanMin = std::max(0.0, scanMin);
  const int nPoints = std::max(
      1, static_cast<int>(std::floor((scanMax-scanMin)/scanStep + 0.5)) + 1);
  out.significanceGraph.reset(new TGraph());
  out.purityGraph.reset(new TGraph());
  out.purityWeightedSignificanceGraph.reset(new TGraph());
  out.signalGraph.reset(new TGraph());
  out.backgroundGraph.reset(new TGraph());
  out.significanceGraph->SetName("gDeltaM_truePhiPhi_significance");
  out.purityGraph->SetName("gDeltaM_truePhiPhi_purity_SoverB");
  out.purityWeightedSignificanceGraph->SetName(
      "gDeltaM_truePhiPhi_purityWeightedSignificanceProxy");
  out.signalGraph->SetName("gDeltaM_truePhiPhi_signalYield");
  out.backgroundGraph->SetName("gDeltaM_nonSS_backgroundYield");

  std::ofstream csv(Form("%s_deltaM_pairSignificance.csv", outPrefix));
  csv << "DeltaM_max,S_truePhiPhi,B_nonSS,S_over_sqrt_SplusB,S_over_B,"
         "S_over_SplusB,purity_weighted_significance_proxy\n";

  double closestDistance = std::numeric_limits<double>::infinity();
  double maximumPurity = 0.0;
  for (int ip = 0; ip < nPoints; ++ip) {
    const double cut = std::min(scanMax, scanMin + ip*scanStep);
    double signal = 0.0;
    double background = 0.0;

    for (int ix = 1; ix <= h2Mass->GetNbinsX(); ++ix) {
      const double m1 = h2Mass->GetXaxis()->GetBinCenter(ix);
      for (int iy = 1; iy <= h2Mass->GetNbinsY(); ++iy) {
        const double count = h2Mass->GetBinContent(ix, iy);
        if (!(count > 0.0)) continue;
        const double m2 = h2Mass->GetYaxis()->GetBinCenter(iy);
        const double radialDeltaM = std::hypot(m1-fit.mean, m2-fit.mean);
        if (radialDeltaM >= cut) continue;
        const double pSS = PosteriorSSProbability(m1, m2, fit);
        if (!(pSS >= 0.0)) continue;
        signal += count*pSS;
        background += count*(1.0-pSS);
      }
    }

    const double denominator = signal + background;
    const double significance = denominator > 0.0
        ? signal/std::sqrt(denominator) : 0.0;
    const double signalToBackground = background > 0.0
        ? signal/background : 0.0;
    const double purity = denominator > 0.0 ? signal/denominator : 0.0;
    const double purityWeightedSignificance = significance*purity;
    out.significanceGraph->SetPoint(ip, cut, significance);
    out.purityGraph->SetPoint(ip, cut, signalToBackground);
    out.purityWeightedSignificanceGraph->SetPoint(
        ip, cut, purityWeightedSignificance);
    out.signalGraph->SetPoint(ip, cut, signal);
    out.backgroundGraph->SetPoint(ip, cut, background);
    csv << cut << "," << signal << "," << background << ","
        << significance << "," << signalToBackground << ","
        << purity << "," << purityWeightedSignificance << "\n";
    maximumPurity = std::max(maximumPurity, signalToBackground);

    if (significance > out.bestSignificance) {
      out.bestSignificance = significance;
      out.bestCut = cut;
    }
    if (purityWeightedSignificance > out.bestPurityWeightedSignificance) {
      out.bestPurityWeightedSignificance = purityWeightedSignificance;
      out.bestPurityWeightedCut = cut;
    }
    const double distance = std::abs(cut-selectionCut);
    if (distance < closestDistance) {
      closestDistance = distance;
      out.significanceAtSelection = significance;
      out.purityWeightedAtSelection = purityWeightedSignificance;
    }
  }
  csv.close();

  auto* c = new TCanvas("cDeltaMTruePairSignificance",
                        "true phi-phi significance and purity versus DeltaM",
                        1050, 780);
  c->SetLeftMargin(0.14);
  c->SetRightMargin(0.14);
  c->SetBottomMargin(0.13);
  c->SetTopMargin(0.08);
  c->SetTicks(1, 1);
  c->SetGridy();

  out.significanceGraph->SetTitle(
      "True #phi#phi-pair significance and purity versus radial #DeltaM_{#phi};"
      "#DeltaM_{#phi}^{max} (GeV/#it{c}^{2});"
      "S_{#phi#phi}/#sqrt{S_{#phi#phi}+B_{non-SS}}");
  out.significanceGraph->SetLineColor(kBlue+1);
  out.significanceGraph->SetMarkerColor(kBlue+1);
  out.significanceGraph->SetLineWidth(3);
  out.significanceGraph->SetMarkerStyle(20);
  out.significanceGraph->SetMarkerSize(0.9);
  out.significanceGraph->GetXaxis()->SetTitleSize(0.050);
  out.significanceGraph->GetYaxis()->SetTitleSize(0.050);
  out.significanceGraph->GetXaxis()->SetLabelSize(0.043);
  out.significanceGraph->GetYaxis()->SetLabelSize(0.043);
  out.significanceGraph->GetYaxis()->SetTitleOffset(1.25);
  const double significanceDisplayMax =
      std::max(1.0, 1.22*out.bestSignificance);
  out.significanceGraph->SetMinimum(0.0);
  out.significanceGraph->SetMaximum(significanceDisplayMax);
  out.significanceGraph->Draw("ALP");

  // Draw S/B on the same canvas with an independent right-hand scale.  The
  // stored purity graph keeps the physical S/B values; only this display copy
  // is rescaled to the left-axis coordinate system.
  const double purityDisplayMax = std::max(1.0e-6, 1.22*maximumPurity);
  std::unique_ptr<TGraph> purityDisplay(new TGraph(*out.purityGraph));
  for (int ip = 0; ip < purityDisplay->GetN(); ++ip) {
    double x = 0.0;
    double purity = 0.0;
    purityDisplay->GetPoint(ip, x, purity);
    purityDisplay->SetPoint(
        ip, x, purity*significanceDisplayMax/purityDisplayMax);
  }
  purityDisplay->SetLineColor(kMagenta+2);
  purityDisplay->SetMarkerColor(kMagenta+2);
  purityDisplay->SetLineWidth(3);
  purityDisplay->SetLineStyle(7);
  purityDisplay->SetMarkerStyle(24);
  purityDisplay->SetMarkerSize(0.9);
  purityDisplay->Draw("LP same");

  gPad->Update();
  auto* purityAxis = new TGaxis(
      gPad->GetUxmax(), 0.0, gPad->GetUxmax(), significanceDisplayMax,
      0.0, purityDisplayMax, 510, "+L");
  purityAxis->SetTitle("purity  S_{#phi#phi}/B_{non-SS}");
  purityAxis->SetTitleColor(kMagenta+2);
  purityAxis->SetLabelColor(kMagenta+2);
  purityAxis->SetLineColor(kMagenta+2);
  purityAxis->SetTitleSize(0.050);
  purityAxis->SetLabelSize(0.043);
  purityAxis->SetTitleOffset(1.15);
  purityAxis->Draw();

  auto* currentLine = new TLine(
      selectionCut, 0.0, selectionCut,
      significanceDisplayMax);
  currentLine->SetLineColor(kRed+1);
  currentLine->SetLineStyle(2);
  currentLine->SetLineWidth(3);
  //  currentLine->Draw("same");
  auto* bestSignificanceLine = new TLine(
      out.bestCut, 0.0, out.bestCut,
      significanceDisplayMax);
  bestSignificanceLine->SetLineColor(kBlue+1);
  bestSignificanceLine->SetLineStyle(3);
  bestSignificanceLine->SetLineWidth(2);
  //bestSignificanceLine->Draw("same");

  auto* leg = new TLegend(0.43, 0.63, 0.89, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.032);
  leg->AddEntry(out.significanceGraph.get(),
                "significance  S/#sqrt{S+B}", "lp");
  leg->AddEntry(purityDisplay.get(), "purity  S/B", "lp");
  leg->AddEntry(currentLine,
                Form("analysis cut %.3f: Z=%.2f", selectionCut,
                     out.significanceAtSelection), "l");
  leg->AddEntry(bestSignificanceLine,
                Form("maximum Z %.2f at %.3f", out.bestSignificance,
                     out.bestCut), "l");
  leg->Draw();

  TLatex text;
  text.SetNDC();
  text.SetTextSize(0.033);
  text.DrawLatex(0.17, 0.88,
                 Form("%.1f < p_{T}^{#phi#phi} < %.1f GeV/#it{c}",
                      ptMin, ptMax));
  text.DrawLatex(0.17, 0.83,
                 Form("%.2f < y_{#phi#phi} < %.2f", rapidityMin,
                      rapidityMax));
  text.DrawLatex(0.17, 0.78,
                 Form("%.3f < M_{#phi#phi} < %.3f GeV/#it{c}^{2}",
                      mPairMin, mPairMax));
  text.DrawLatex(0.17, 0.73,
                 "#DeltaM_{#phi}=#sqrt{(m_{1}-m_{#phi})^{2}+(m_{2}-m_{#phi})^{2}}");

  c->SaveAs(Form("%s_Fig6_deltaM_significancePurity.png", outPrefix));
  c->SaveAs(Form("%s_Fig6_deltaM_significancePurity.pdf", outPrefix));
  delete c;

  // Keep the combined optimization metric on its own figure.  This is a
  // purity-weighted significance proxy, not a formal discovery significance.
  auto* cPurityWeighted = new TCanvas(
      "cDeltaMPurityWeightedSignificance",
      "significance times purity versus DeltaM", 1000, 780);
  cPurityWeighted->SetLeftMargin(0.14);
  cPurityWeighted->SetRightMargin(0.04);
  cPurityWeighted->SetBottomMargin(0.13);
  cPurityWeighted->SetTopMargin(0.08);
  cPurityWeighted->SetTicks(1, 1);
  cPurityWeighted->SetGridy();

  out.purityWeightedSignificanceGraph->SetTitle(
      "Significance #times purity versus radial #DeltaM_{#phi};"
      "#DeltaM_{#phi}^{max} (GeV/#it{c}^{2});"
      "[S/#sqrt{S+B}] #times [S/(S+B)]");
  out.purityWeightedSignificanceGraph->SetLineColor(kGreen+2);
  out.purityWeightedSignificanceGraph->SetMarkerColor(kGreen+2);
  out.purityWeightedSignificanceGraph->SetLineWidth(3);
  out.purityWeightedSignificanceGraph->SetMarkerStyle(21);
  out.purityWeightedSignificanceGraph->SetMarkerSize(0.9);
  out.purityWeightedSignificanceGraph->GetXaxis()->SetTitleSize(0.050);
  out.purityWeightedSignificanceGraph->GetYaxis()->SetTitleSize(0.050);
  out.purityWeightedSignificanceGraph->GetXaxis()->SetLabelSize(0.043);
  out.purityWeightedSignificanceGraph->GetYaxis()->SetLabelSize(0.043);
  out.purityWeightedSignificanceGraph->GetYaxis()->SetTitleOffset(1.25);
  const double purityWeightedDisplayMax = std::max(
      1.0, 1.22*out.bestPurityWeightedSignificance);
  out.purityWeightedSignificanceGraph->SetMinimum(0.0);
  out.purityWeightedSignificanceGraph->SetMaximum(purityWeightedDisplayMax);
  out.purityWeightedSignificanceGraph->Draw("ALP");

  auto* currentPurityWeightedLine = new TLine(
      selectionCut, 0.0, selectionCut, purityWeightedDisplayMax);
  currentPurityWeightedLine->SetLineColor(kRed+1);
  currentPurityWeightedLine->SetLineStyle(2);
  currentPurityWeightedLine->SetLineWidth(3);
  // currentPurityWeightedLine->Draw("same");

  auto* bestPurityWeightedLine = new TLine(
      out.bestPurityWeightedCut, 0.0,
      out.bestPurityWeightedCut, purityWeightedDisplayMax);
  bestPurityWeightedLine->SetLineColor(kGreen+3);
  bestPurityWeightedLine->SetLineStyle(7);
  bestPurityWeightedLine->SetLineWidth(3);
  // bestPurityWeightedLine->Draw("same");

  auto* purityWeightedLegend = new TLegend(0.43, 0.68, 0.91, 0.90);
  purityWeightedLegend->SetBorderSize(0);
  purityWeightedLegend->SetFillStyle(0);
  purityWeightedLegend->SetTextSize(0.032);
  purityWeightedLegend->AddEntry(
      out.purityWeightedSignificanceGraph.get(),
      "purity-weighted significance proxy", "lp");
  purityWeightedLegend->AddEntry(
      currentPurityWeightedLine,
      Form("analysis cut %.3f: %.2f", selectionCut,
           out.purityWeightedAtSelection), "l");
 
  purityWeightedLegend->Draw();

  TLatex purityWeightedText;
  purityWeightedText.SetNDC();
  purityWeightedText.SetTextSize(0.033);
  purityWeightedText.DrawLatex(
      0.17, 0.88,
      Form("%.1f < p_{T}^{#phi#phi} < %.1f GeV/#it{c}", ptMin, ptMax));

  purityWeightedText.DrawLatex(
      0.17, 0.78,
      Form("full range: %.3f < M_{#phi#phi} < %.3f GeV/#it{c}^{2}",
           mPairMin, mPairMax));
  
  cPurityWeighted->SaveAs(Form(
      "%s_Fig7_deltaM_significanceTimesPurity.png", outPrefix));
  cPurityWeighted->SaveAs(Form(
      "%s_Fig7_deltaM_significanceTimesPurity.pdf", outPrefix));
  delete cPurityWeighted;
  return out;
}


TH2D* BuildM1M2ForMPairBinFromSparse(THnSparseF* hSparse,
                                     const char* name,
                                     double ptMin,
                                     double ptMax,
                                     double phiMassMin,
                                     double phiMassMax,
                                     double mPairLo,
                                     double mPairHi)
{
  if (!hSparse) return nullptr;

  const int nDim = hSparse->GetNdimensions();
  if (nDim <= std::max({kAxisMPair, kAxisPtPair, kAxisRapidity,
                        kAxisM1, kAxisM2})) {
    std::cerr << "ERROR: THnSparse has only " << nDim
              << " dimensions, but axes 0,1,3,4,5 are requested.\n";
    return nullptr;
  }

  SparseAxisRangeGuard guard(hSparse);
  // Deliberately do not set any range on kAxisDeltaM here: every 2D
  // m(phi1)-m(phi2) fit must use the full DeltaM distribution.
  if (!ApplyRapiditySelection(hSparse) ||
      !SetAxisRangeByValue(hSparse->GetAxis(kAxisMPair), mPairLo, mPairHi) ||
      !SetAxisRangeByValue(hSparse->GetAxis(kAxisPtPair), ptMin, ptMax) ||
      !SetAxisRangeByValue(hSparse->GetAxis(kAxisM1), phiMassMin, phiMassMax) ||
      !SetAxisRangeByValue(hSparse->GetAxis(kAxisM2), phiMassMin, phiMassMax)) {
    return nullptr;
  }

  // ROOT's TH2 overload is Projection(yDim,xDim): first argument is y.
  TH2D* h2 = hSparse->Projection(kAxisM2, kAxisM1, "E");
  if (!h2) return nullptr;
  h2->SetName(name);
  h2->SetDirectory(nullptr);
  h2->SetTitle(";M_{#phi,1} (GeV/#it{c}^{2});M_{#phi,2} (GeV/#it{c}^{2})");
  return h2;
}

Fit2DResult Fit2DMassOneMPairBin(TH2D* h2,
                                 const Fit2DResult& globalFit,
                                 double phiMassMin,
                                 double phiMassMax,
                                 int chebOrder,
                                 int yieldMode,
                                 const char* savePrefix = "")
{
  // yieldMode = 1: fix all shape parameters from the global 2D fit; float only SS/SB/BS/BB yields.
  // yieldMode = 2: float shape parameters as well. This is less stable and needs high statistics per bin.
  Fit2DResult res;
  res.phiMassMin = phiMassMin;
  res.phiMassMax = phiMassMax;
  res.chebOrder = std::max(1, std::min(chebOrder, 3));

  if (!h2 || h2->Integral() <= 0.0) return res;

  const bool fixShape = (yieldMode == 1);
  const double nTot = h2->Integral(1, h2->GetNbinsX(), 1, h2->GetNbinsY());
  if (nTot <= 0.0) return res;

  RooRealVar m1("m1", "M_{#phi,1}", phiMassMin, phiMassMax, "GeV/c^{2}");
  RooRealVar m2("m2", "M_{#phi,2}", phiMassMin, phiMassMax, "GeV/c^{2}");
  RooDataHist data2D("data2D_bin", "data2D_bin", RooArgList(m1, m2), h2);

  RooRealVar mean("mean_bin", "m_{#phi}", globalFit.mean, 1.0170, 1.0215);
  RooRealVar gamma("gamma_bin", "#Gamma_{#phi}", globalFit.gamma,
                   gPhiBWGammaMin, gPhiBWGammaMax);
  if (fixShape) {
    mean.setConstant(true);
    gamma.setConstant(true);
  }

  RooBreitWigner sig1("sig1_bin", "BW(M1)", m1, mean, gamma);
  RooBreitWigner sig2("sig2_bin", "BW(M2)", m2, mean, gamma);

  RooRealVar c1m1("m1_c1_bin", "m1_c1_bin", globalFit.c1m1, -1.0, 1.0);
  RooRealVar c2m1("m1_c2_bin", "m1_c2_bin", globalFit.c2m1, -1.0, 1.0);
  RooRealVar c3m1("m1_c3_bin", "m1_c3_bin", globalFit.c3m1, -1.0, 1.0);
  RooRealVar c1m2("m2_c1_bin", "m2_c1_bin", globalFit.c1m2, -1.0, 1.0);
  RooRealVar c2m2("m2_c2_bin", "m2_c2_bin", globalFit.c2m2, -1.0, 1.0);
  RooRealVar c3m2("m2_c3_bin", "m2_c3_bin", globalFit.c3m2, -1.0, 1.0);

  if (fixShape) {
    c1m1.setConstant(true); c2m1.setConstant(true); c3m1.setConstant(true);
    c1m2.setConstant(true); c2m2.setConstant(true); c3m2.setConstant(true);
  }

  RooArgList coeffs1;
  RooArgList coeffs2;
  if (res.chebOrder >= 1) { coeffs1.add(c1m1); coeffs2.add(c1m2); }
  if (res.chebOrder >= 2) { coeffs1.add(c2m1); coeffs2.add(c2m2); }
  if (res.chebOrder >= 3) { coeffs1.add(c3m1); coeffs2.add(c3m2); }

  RooChebychev bkg1("bkg1_bin", "B(M1)", m1, coeffs1);
  RooChebychev bkg2("bkg2_bin", "B(M2)", m2, coeffs2);

  RooProdPdf pdfSS("pdfSS_bin", "S(M1)S(M2)", RooArgSet(sig1, sig2));
  RooProdPdf pdfSB("pdfSB_bin", "S(M1)B(M2)", RooArgSet(sig1, bkg2));
  RooProdPdf pdfBS("pdfBS_bin", "B(M1)S(M2)", RooArgSet(bkg1, sig2));
  RooProdPdf pdfBB("pdfBB_bin", "B(M1)B(M2)", RooArgSet(bkg1, bkg2));

  const double gTot = std::max(1.0, globalFit.nSS + globalFit.nSB + globalFit.nBS + globalFit.nBB);
  const double fSS = std::max(0.02, globalFit.nSS / gTot);
  const double fSB = std::max(0.02, globalFit.nSB / gTot);
  const double fBS = std::max(0.02, globalFit.nBS / gTot);
  const double fBB = std::max(0.02, globalFit.nBB / gTot);

  RooRealVar nSS("nSS_bin", "N_{SS}", fSS * nTot, 0.0, 2.0 * nTot);
  RooRealVar nSB("nSB_bin", "N_{SB}", fSB * nTot, 0.0, 2.0 * nTot);
  RooRealVar nBS("nBS_bin", "N_{BS}", fBS * nTot, 0.0, 2.0 * nTot);
  RooRealVar nBB("nBB_bin", "N_{BB}", fBB * nTot, 0.0, 2.0 * nTot);

  RooAddPdf model2D("model2D_bin", "SS+SB+BS+BB",
                    RooArgList(pdfSS, pdfSB, pdfBS, pdfBB),
                    RooArgList(nSS, nSB, nBS, nBB));

  std::unique_ptr<RooFitResult> fit(MinimizeExtendedNLL(model2D, data2D));

  res.ok = fit && fit->status() == 0;
  res.mean = mean.getVal();
  res.meanErr = mean.getError();
  res.gamma = gamma.getVal();
  res.gammaErr = gamma.getError();
  res.nSS = nSS.getVal();
  res.nSB = nSB.getVal();
  res.nBS = nBS.getVal();
  res.nBB = nBB.getVal();
  res.nSSErr = nSS.getError();
  res.nSBErr = nSB.getError();
  res.nBSErr = nBS.getError();
  res.nBBErr = nBB.getError();
  res.c1m1 = c1m1.getVal(); res.c2m1 = c2m1.getVal(); res.c3m1 = c3m1.getVal();
  res.c1m2 = c1m2.getVal(); res.c2m2 = c2m2.getVal(); res.c3m2 = c3m2.getVal();

  if (savePrefix && TString(savePrefix).Length() > 0) {
    RooPlot* fr1 = m1.frame(Bins(h2->GetNbinsX()));
    data2D.plotOn(fr1, Name("data_m1_bin"));
    model2D.plotOn(fr1, Name("model_m1_bin"), LineColor(kRed+1), LineWidth(3));
    model2D.plotOn(fr1, Components(pdfSS), Name("SS_m1_bin"), LineColor(kBlue+1), LineStyle(kDashed), LineWidth(2));
    model2D.plotOn(fr1, Components(RooArgSet(pdfSB, pdfBS)), Name("SBBS_m1_bin"), LineColor(kGreen+2), LineStyle(kDotted), LineWidth(2));
    model2D.plotOn(fr1, Components(pdfBB), Name("BB_m1_bin"), LineColor(kMagenta+1), LineStyle(kDashDotted), LineWidth(2));

    RooPlot* fr2 = m2.frame(Bins(h2->GetNbinsY()));
    data2D.plotOn(fr2, Name("data_m2_bin"));
    model2D.plotOn(fr2, Name("model_m2_bin"), LineColor(kRed+1), LineWidth(3));
    model2D.plotOn(fr2, Components(pdfSS), Name("SS_m2_bin"), LineColor(kBlue+1), LineStyle(kDashed), LineWidth(2));
    model2D.plotOn(fr2, Components(RooArgSet(pdfSB, pdfBS)), Name("SBBS_m2_bin"), LineColor(kGreen+2), LineStyle(kDotted), LineWidth(2));
    model2D.plotOn(fr2, Components(pdfBB), Name("BB_m2_bin"), LineColor(kMagenta+1), LineStyle(kDashDotted), LineWidth(2));

    auto* cfit = new TCanvas(Form("c_%s_proj", gSystem->BaseName(savePrefix)),
                             Form("%s projections", savePrefix), 1600, 800);
    cfit->Divide(2,1);
    cfit->cd(1); gPad->SetLeftMargin(0.13); fr1->Draw();
    cfit->cd(2); gPad->SetLeftMargin(0.13); fr2->Draw();
    cfit->SaveAs(Form("%s_projections.png", savePrefix));
    cfit->SaveAs(Form("%s_projections.pdf", savePrefix));
    delete cfit;
    delete fr1;
    delete fr2;

    auto* c2d = new TCanvas(Form("c_%s_2D", gSystem->BaseName(savePrefix)),
                            Form("%s 2D", savePrefix), 900, 800);
    gPad->SetRightMargin(0.14);
    h2->SetStats(false);
    h2->Draw("COLZ");
    c2d->SaveAs(Form("%s_h2data.png", savePrefix));
    c2d->SaveAs(Form("%s_h2data.pdf", savePrefix));
    delete c2d;
  }

  return res;
}


const char* SelectionLabel(int selectionMode);

bool PassSelectedRegion(double m1,
                        double m2,
                        double dm,
                        double mRef,
                        double dmSigMax,
                        int selectionMode,
                        double rectHalfWidth);

struct SelectedPosteriorYields {
  bool ok = false;
  double raw = 0.0;
  double rawVar = 0.0;
  double nSS = 0.0;
  double nNonSS = 0.0;
  double statVarSS = 0.0;
  double statVarNonSS = 0.0;
  double fitVarDiag = 0.0;
  double invalidModelWeight = 0.0;
  Long64_t usedSparseBins = 0;
};

SelectedPosteriorYields ExtractSelectedPosteriorYields(
    THnSparseF* hSparse,
    const Fit2DResult& fit,
    double ptMin,
    double ptMax,
    double phiMassMin,
    double phiMassMax,
    double dmCut,
    double mRef,
    int selectionMode,
    double rectHalfWidth,
    double mPairLo,
    double mPairHi)
{
  // The full-range 2D fit determines the component model in this M(phi phi)
  // bin.  The exact selected sample is then decomposed entry by entry (more
  // precisely, occupied sparse-bin by occupied sparse-bin) using the local
  // posterior probability
  //
  //   P_SS(m1,m2) = N_SS S1 S2 /
  //                  (N_SS S1 S2 + N_SB S1 B2
  //                   + N_BS B1 S2 + N_BB B1 B2).
  //
  // Crucially, PassSelectedRegion uses the ACTUAL stored DeltaM axis for
  // selectionMode==0.  No global DeltaM survival efficiency is assumed.
  SelectedPosteriorYields out;
  if (!hSparse || !fit.ok) return out;

  const int nDim = hSparse->GetNdimensions();
  if (nDim <= std::max({kAxisMPair, kAxisPtPair, kAxisRapidity,
                        kAxisM1, kAxisM2, kAxisDeltaM})) {
    std::cerr << "ERROR: THnSparse has too few dimensions for posterior selected-yield extraction.\n";
    return out;
  }

  TAxis* axMPair = hSparse->GetAxis(kAxisMPair);
  TAxis* axPt    = hSparse->GetAxis(kAxisPtPair);
  TAxis* axY     = hSparse->GetAxis(kAxisRapidity);
  TAxis* axM1    = hSparse->GetAxis(kAxisM1);
  TAxis* axM2    = hSparse->GetAxis(kAxisM2);
  TAxis* axDM    = hSparse->GetAxis(kAxisDeltaM);
  if (!axMPair || !axPt || !axY || !axM1 || !axM2 || !axDM) return out;

  const double nSS = std::max(0.0, fit.nSS);
  const double nSB = std::max(0.0, fit.nSB);
  const double nBS = std::max(0.0, fit.nBS);
  const double nBB = std::max(0.0, fit.nBB);

  // Derivatives of the selected N_SS with respect to the four fitted yields.
  // They are used only for a diagonal approximation to the fit-yield error.
  double dSel_dNSS = 0.0;
  double dSel_dNSB = 0.0;
  double dSel_dNBS = 0.0;
  double dSel_dNBB = 0.0;

  std::vector<int> idx(nDim, 0);
  for (Long64_t ib = 0; ib < hSparse->GetNbins(); ++ib) {
    const double c = hSparse->GetBinContent(ib, idx.data());
    if (c <= 0.0) continue;
    const double cVar = hSparse->GetBinError2(ib);

    const double mPair = axMPair->GetBinCenter(idx[kAxisMPair]);
    const double pt    = axPt->GetBinCenter(idx[kAxisPtPair]);
    const double rapidity = axY->GetBinCenter(idx[kAxisRapidity]);
    const double m1    = axM1->GetBinCenter(idx[kAxisM1]);
    const double m2    = axM2->GetBinCenter(idx[kAxisM2]);
    const double dm    = axDM->GetBinCenter(idx[kAxisDeltaM]);

    if (mPair < mPairLo || mPair >= mPairHi) continue;
    if (pt < ptMin || pt >= ptMax) continue;
    if (rapidity < gRapidityMin || rapidity >= gRapidityMax) continue;
    if (m1 < phiMassMin || m1 >= phiMassMax) continue;
    if (m2 < phiMassMin || m2 >= phiMassMax) continue;
    if (!PassSelectedRegion(m1, m2, dm, mRef, dmCut,
                            selectionMode, rectHalfWidth)) continue;

    const double s1 = BWPdf(m1, fit.phiMassMin, fit.phiMassMax,
                             fit.mean, fit.gamma);
    const double s2 = BWPdf(m2, fit.phiMassMin, fit.phiMassMax,
                             fit.mean, fit.gamma);
    const double b1 = ChebPdf(m1, fit.phiMassMin, fit.phiMassMax,
                              fit.c1m1, fit.c2m1, fit.c3m1,
                              fit.chebOrder);
    const double b2 = ChebPdf(m2, fit.phiMassMin, fit.phiMassMax,
                              fit.c1m2, fit.c2m2, fit.c3m2,
                              fit.chebOrder);

    const double qSS = s1 * s2;
    const double qSB = s1 * b2;
    const double qBS = b1 * s2;
    const double qBB = b1 * b2;

    out.raw += c;
    out.rawVar += cVar;

    const double aSS = nSS * qSS;
    const double aSB = nSB * qSB;
    const double aBS = nBS * qBS;
    const double aBB = nBB * qBB;
    const double total = aSS + aSB + aBS + aBB;
    if (!(total > 0.0) || !std::isfinite(total)) {
      out.invalidModelWeight += c;
      continue;
    }

    const double pSS = std::max(0.0, std::min(1.0, aSS / total));
    const double pNonSS = 1.0 - pSS;

    out.nSS += c * pSS;
    out.nNonSS += c * pNonSS;
    out.statVarSS += cVar * pSS * pSS;
    out.statVarNonSS += cVar * pNonSS * pNonSS;
    out.usedSparseBins++;

    // Analytic derivatives of pSS with respect to the fitted component yields.
    const double invT2 = 1.0 / (total * total);
    dSel_dNSS += c * qSS * (total - aSS) * invT2;
    dSel_dNSB += c * (-aSS * qSB) * invT2;
    dSel_dNBS += c * (-aSS * qBS) * invT2;
    dSel_dNBB += c * (-aSS * qBB) * invT2;
  }

  // Diagonal approximation only: RooFit yield correlations are not retained in
  // Fit2DResult.  This is still more local than propagating one global cut
  // efficiency.  The central values do not depend on this approximation.
  const double eNSS = (fit.nSSErr > 0.0 ? fit.nSSErr : 0.0);
  const double eNSB = (fit.nSBErr > 0.0 ? fit.nSBErr : 0.0);
  const double eNBS = (fit.nBSErr > 0.0 ? fit.nBSErr : 0.0);
  const double eNBB = (fit.nBBErr > 0.0 ? fit.nBBErr : 0.0);
  out.fitVarDiag = dSel_dNSS * dSel_dNSS * eNSS * eNSS
                 + dSel_dNSB * dSel_dNSB * eNSB * eNSB
                 + dSel_dNBS * dSel_dNBS * eNBS * eNBS
                 + dSel_dNBB * dSel_dNBB * eNBB * eNBB;

  const double invalidFrac = (out.raw > 0.0 ? out.invalidModelWeight / out.raw : 1.0);
  out.ok = (out.raw > 0.0 && invalidFrac < 1.0e-9 &&
            std::isfinite(out.nSS) && std::isfinite(out.nNonSS));
  return out;
}

TH1D* BuildSSAndNonSSYieldsVsMPair_PosteriorSelected(
                                                       THnSparseF* hSparse,
                                                       const Fit2DResult& globalFit,
                                                       const char* nameSS,
                                                       double ptMin,
                                                       double ptMax,
                                                       double phiMassMin,
                                                       double phiMassMax,
                                                       int chebOrder,
                                                       double dmCut,
                                                       double mPairMin,
                                                       double mPairMax,
                                                       int nPairBins,
                                                       double minRawFor2DFit,
                                                       TH1D*& hNonSSOut,
                                                       TH1D*& hRawOut,
                                                       TH1D*& hFitStatusOut,
                                                       int selectionMode = 0,
                                                       double rectHalfWidth = 0.010,
                                                       int pairMassFitMode = 1,
                                                       const char* saveFitDir = "")
{
  // Selected-sample strategy:
  //   pairMassFitMode=0: use the global 2D-fit component fractions in every
  //                      M(phi phi) bin; no local 2D fit is performed.
  //   pairMassFitMode=1: fit the full m1-m2 plane separately in every M(phi phi)
  //                      bin, with global shape parameters fixed and four yields free.
  //   pairMassFitMode=2: fit the full m1-m2 plane separately in every M(phi phi)
  //                      bin, allowing the mass-shape parameters to float as well.
  //   Then loop over the actual THnSparse entries in that M(phi phi) bin.
  //   3) Apply the exact requested selection.  For selectionMode==0 this is the
  //      exact STORED DeltaM-axis cut, including its pT-dependent construction.
  //   4) For every selected entry, sum the local posterior P_SS(m1,m2).
  //
  // Hence
  //   N_SS(selected)    = sum_selected c * P_SS,
  //   N_nonSS(selected) = sum_selected c * (1-P_SS),
  // and closure to the selected raw count is exact up to floating precision.
  // No assumption is made that a DeltaM efficiency is constant versus M(phi phi).

  auto* hSS = new TH1D(nameSS,
                       ";M_{#phi#phi} (GeV/#it{c}^{2});N_{SS}^{selected} from posterior sum",
                       nPairBins, mPairMin, mPairMax);
  hSS->SetDirectory(nullptr);
  hSS->Sumw2();

  hNonSSOut = new TH1D(Form("%s_nonSS", nameSS),
                       ";M_{#phi#phi} (GeV/#it{c}^{2});N_{nonSS}^{selected} from posterior sum",
                       nPairBins, mPairMin, mPairMax);
  hNonSSOut->SetDirectory(nullptr);
  hNonSSOut->Sumw2();

  hRawOut = new TH1D(Form("%s_raw", nameSS),
                     ";M_{#phi#phi} (GeV/#it{c}^{2});raw candidate pairs in selected region",
                     nPairBins, mPairMin, mPairMax);
  hRawOut->SetDirectory(nullptr);
  hRawOut->Sumw2();

  hFitStatusOut = new TH1D(Form("%s_fitStatus", nameSS),
                           ";M_{#phi#phi} (GeV/#it{c}^{2});status flag",
                           nPairBins, mPairMin, mPairMax);
  hFitStatusOut->SetDirectory(nullptr);
  hFitStatusOut->Sumw2();

  const double mRef = globalFit.mean; // used only by rectangular selections

  std::cout << "\n========== Per-Mpair full-range 2D fit + exact selected posterior sum ==========" << std::endl;
  std::cout << "Extraction region = " << SelectionLabel(selectionMode)
            << ", dmCut = " << dmCut
            << ", rectHalfWidth = " << rectHalfWidth << std::endl;
  std::cout << "pairMassFitMode = " << pairMassFitMode
            << " (0=global fractions, 1=local yields/global shapes, "
            << "2=local yields and local shapes)" << std::endl;
  if (pairMassFitMode == 0) {
    std::cout << "No local 2D mass fit is performed; global component fractions are used." << std::endl;
  } else {
    std::cout << "Each Mpair bin is fitted without a DeltaM cut; then the exact stored selection is applied while summing P_SS and 1-P_SS." << std::endl;
  }
  std::cout << "No global component survival efficiency is used.\n" << std::endl;

  std::unique_ptr<TFile> all2DInputFile;
  std::unique_ptr<std::ofstream> all2DFitSummary;

  if (saveFitDir && TString(saveFitDir).Length() > 0) {
    all2DInputFile.reset(
      TFile::Open(
        Form("%s/all_2D_invariant_mass_inputs.root", saveFitDir),
        "RECREATE"));

    all2DFitSummary.reset(
      new std::ofstream(
        Form("%s/all_2D_invariant_mass_fit_results.csv", saveFitDir)));

    (*all2DFitSummary)
      << "bin,M_low,M_high,raw_full,status,"
      << "mean,phi_bw_width,"
      << "N_SS,N_SS_err,N_SB,N_SB_err,N_BS,N_BS_err,N_BB,N_BB_err,"
      << "selected_raw,selected_SS,selected_SS_err,"
      << "selected_nonSS,selected_nonSS_err,closure\n";
  }

  for (int ib = 1; ib <= nPairBins; ++ib) {
    const double lo = hSS->GetXaxis()->GetBinLowEdge(ib);
    const double hi = hSS->GetXaxis()->GetBinUpEdge(ib);

    std::unique_ptr<TH2D> h2Full(BuildM1M2ForMPairBinFromSparse(
        hSparse,
        Form("h2_fullDM_M1M2_MpairBin_%03d", ib),
        ptMin, ptMax,
        phiMassMin, phiMassMax,
        lo, hi));

    const double rawFull = h2Full
      ? h2Full->Integral(1, h2Full->GetNbinsX(), 1, h2Full->GetNbinsY())
      : 0.0;

    // Preserve every coarse-bin 2D invariant-mass distribution, including
    // low-statistics bins that are skipped by the fitter.
    if (all2DInputFile && !all2DInputFile->IsZombie() && h2Full) {
      all2DInputFile->cd();
      h2Full->Write(
        Form("h2_MpairBin_%03d_M%.4f_%.4f", ib, lo, hi));
    }

    const bool needsLocalFit = (pairMassFitMode != 0);
    if (!h2Full || rawFull <= 0.0 ||
        (needsLocalFit && rawFull < minRawFor2DFit)) {
      hFitStatusOut->SetBinContent(ib, -2);
      if (all2DFitSummary) {
        (*all2DFitSummary)
          << ib << "," << lo << "," << hi << "," << rawFull << ",-2,"
          << "0,0,0,0,0,0,0,0,0,0,0,"
          << "0,0,0,0,0,0\n";
      }
      std::cout << Form("bin %02d  M=[%.4f,%.4f): rawFull=%.0f  SKIP low stat",
                        ib, lo, hi, rawFull) << std::endl;
      continue;
    }

    Fit2DResult rFull;
    if (pairMassFitMode == 0) {
      // P_SS depends only on relative component normalizations.  Reusing the
      // global result therefore implements the global-fraction option.
      rFull = globalFit;
      rFull.ok = globalFit.ok;
    } else {
      TString savePrefix = "";
      if (saveFitDir && TString(saveFitDir).Length() > 0) {
        savePrefix = Form("%s/bin%03d_M%.3f_%.3f_fullDM",
                          saveFitDir, ib, lo, hi);
      }

      const int localYieldMode = (pairMassFitMode == 2 ? 2 : 1);
      rFull = Fit2DMassOneMPairBin(
          h2Full.get(), globalFit,
          phiMassMin, phiMassMax,
          chebOrder,
          localYieldMode,
          savePrefix.Data());
    }

    if (!rFull.ok) {
      hFitStatusOut->SetBinContent(ib, 1);
      if (all2DFitSummary) {
        (*all2DFitSummary)
          << ib << "," << lo << "," << hi << "," << rawFull << ",1,"
          << rFull.mean << "," << rFull.gamma << ","
          << rFull.nSS << "," << rFull.nSSErr << ","
          << rFull.nSB << "," << rFull.nSBErr << ","
          << rFull.nBS << "," << rFull.nBSErr << ","
          << rFull.nBB << "," << rFull.nBBErr << ","
          << "0,0,0,0,0,0\n";
      }
      std::cout << Form("bin %02d  M=[%.4f,%.4f): rawFull=%.0f  2D component model failed",
                        ib, lo, hi, rawFull) << std::endl;
      continue;
    }

    const SelectedPosteriorYields sel = ExtractSelectedPosteriorYields(
        hSparse, rFull,
        ptMin, ptMax,
        phiMassMin, phiMassMax,
        dmCut, mRef,
        selectionMode, rectHalfWidth,
        lo, hi);

    if (!sel.ok) {
      hFitStatusOut->SetBinContent(ib, 2);
      if (all2DFitSummary) {
        (*all2DFitSummary)
          << ib << "," << lo << "," << hi << "," << rawFull << ",2,"
          << rFull.mean << "," << rFull.gamma << ","
          << rFull.nSS << "," << rFull.nSSErr << ","
          << rFull.nSB << "," << rFull.nSBErr << ","
          << rFull.nBS << "," << rFull.nBSErr << ","
          << rFull.nBB << "," << rFull.nBBErr << ","
          << "0,0,0,0,0,0\n";
      }
      std::cout << Form("bin %02d  M=[%.4f,%.4f): rawFull=%.0f  selected posterior extraction failed/empty",
                        ib, lo, hi, rawFull) << std::endl;
      continue;
    }

    const double closure = sel.raw - sel.nSS - sel.nNonSS;
    const double eSS = std::sqrt(std::max(0.0, sel.statVarSS + sel.fitVarDiag));
    const double eNonSS = std::sqrt(std::max(0.0, sel.statVarNonSS + sel.fitVarDiag));

    hRawOut->SetBinContent(ib, sel.raw);
    hRawOut->SetBinError(ib, std::sqrt(std::max(0.0, sel.rawVar)));
    hSS->SetBinContent(ib, sel.nSS);
    hSS->SetBinError(ib, eSS);
    hNonSSOut->SetBinContent(ib, sel.nNonSS);
    hNonSSOut->SetBinError(ib, eNonSS);
    hFitStatusOut->SetBinContent(ib, 0);

    if (all2DFitSummary) {
      (*all2DFitSummary)
        << ib << "," << lo << "," << hi << "," << rawFull << ",0,"
        << rFull.mean << "," << rFull.gamma << ","
        << rFull.nSS << "," << rFull.nSSErr << ","
        << rFull.nSB << "," << rFull.nSBErr << ","
        << rFull.nBS << "," << rFull.nBSErr << ","
        << rFull.nBB << "," << rFull.nBBErr << ","
        << sel.raw << ","
        << sel.nSS << "," << eSS << ","
        << sel.nNonSS << "," << eNonSS << ","
        << closure << "\n";
    }

    std::cout << Form(
      "bin %02d  M=[%.4f,%.4f): rawFull=%.0f rawSelected=%.1f  "
      "fullYields=(SS %.1f,SB %.1f,BS %.1f,BB %.1f)  "
      "N_SS(selected)=%.2f +/- %.2f  N_nonSS(selected)=%.2f +/- %.2f  "
      "closure=%.3g  invalidModelWeight=%.3g  sparseBins=%lld",
      ib, lo, hi, rawFull, sel.raw,
      rFull.nSS, rFull.nSB, rFull.nBS, rFull.nBB,
      sel.nSS, eSS, sel.nNonSS, eNonSS,
      closure, sel.invalidModelWeight,
      static_cast<long long>(sel.usedSparseBins)) << std::endl;
  }

  if (all2DFitSummary) {
    all2DFitSummary->close();
  }
  if (all2DInputFile && !all2DInputFile->IsZombie()) {
    all2DInputFile->Write();
    all2DInputFile->Close();
  }

  std::cout << "===============================================================================\n" << std::endl;
  return hSS;
}


} // namespace


namespace {

void ScaleHistogramToUnitArea(TH1D* h, double floorFrac = 1.0e-9)
{
  if (!h) return;
  h->SetDirectory(nullptr);
  h->Sumw2(false);

  double maxVal = 0.0;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    maxVal = std::max(maxVal, h->GetBinContent(i));
  }
  const double floor = std::max(1.0e-30, floorFrac * std::max(1.0, maxVal));

  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    double y = h->GetBinContent(i);
    if (!std::isfinite(y) || y < floor) y = floor;
    h->SetBinContent(i, y);
    h->SetBinError(i, 0.0);
  }

  const double integral = h->Integral(1, h->GetNbinsX());
  if (integral > 0.0) h->Scale(1.0 / integral);
}


const char* SelectionLabel(int selectionMode)
{
  if (selectionMode == 1) return "rectSig";
  if (selectionMode == 2) return "rectSide";
  if (selectionMode == 3) return "fullMass";
  return "dmSig";
}

bool PassSelectedRegion(double m1,
                        double m2,
                        double dm,
                        double mRef,
                        double dmSigMax,
                        int selectionMode,
                        double rectHalfWidth)
{
  if (selectionMode == 1) {
    return (std::abs(m1 - mRef) < rectHalfWidth &&
            std::abs(m2 - mRef) < rectHalfWidth);
  }
  if (selectionMode == 2) {
    return !(std::abs(m1 - mRef) < rectHalfWidth &&
             std::abs(m2 - mRef) < rectHalfWidth);
  }
  if (selectionMode == 3) {
    // Full selected phi-candidate mass square: no DeltaM or rectangular cut.
    return true;
  }
  return (dm >= 0.0 && dm < dmSigMax);
}

TH1D* ProjectPairMassNativeSelection(THnSparseF* hSparse,
                                     const char* name,
                                     double ptMin,
                                     double ptMax,
                                     double phiMassMin,
                                     double phiMassMax,
                                     double dmSigMax,
                                     double mRef,
                                     int selectionMode,
                                     double rectHalfWidth,
                                     double mPairMin,
                                     double mPairMax,
                                     int nPairBins)
{
  if (!hSparse) return nullptr;

  const int nDim = hSparse->GetNdimensions();
  if (nDim <= std::max({kAxisMPair, kAxisPtPair, kAxisRapidity,
                        kAxisM1, kAxisM2, kAxisDeltaM})) {
    std::cerr << "ERROR: THnSparse has too few dimensions for ProjectPairMassNativeSelection.\n";
    return nullptr;
  }

  if (selectionMode == 0) {
    return ProjectPairMassNativeRanges(
        hSparse, name, ptMin, ptMax,
        phiMassMin, phiMassMax, phiMassMin, phiMassMax,
        true, 0.0, dmSigMax,
        mPairMin, mPairMax, nPairBins);
  }
  if (selectionMode == 1) {
    return ProjectPairMassNativeRanges(
        hSparse, name, ptMin, ptMax,
        mRef - rectHalfWidth, mRef + rectHalfWidth,
        mRef - rectHalfWidth, mRef + rectHalfWidth,
        false, 0.0, 0.0,
        mPairMin, mPairMax, nPairBins);
  }
  if (selectionMode == 3) {
    return ProjectPairMassNativeRanges(
        hSparse, name, ptMin, ptMax,
        phiMassMin, phiMassMax, phiMassMin, phiMassMax,
        false, 0.0, 0.0,
        mPairMin, mPairMax, nPairBins);
  }
  if (selectionMode != 2) {
    std::cerr << "ERROR: invalid selectionMode=" << selectionMode << std::endl;
    return nullptr;
  }

  // Rectangle complement, constructed as four disjoint native projections.
  // Using full-minus-central would overestimate the errors because the two
  // histograms are correlated.
  const double rectLo = std::max(phiMassMin, mRef - rectHalfWidth);
  const double rectHi = std::min(phiMassMax, mRef + rectHalfWidth);
  std::vector<std::unique_ptr<TH1D>> pieces;
  auto addPiece = [&](double m1Lo, double m1Hi,
                      double m2Lo, double m2Hi,
                      const char* suffix) {
    if (!(m1Hi > m1Lo) || !(m2Hi > m2Lo)) return;
    pieces.emplace_back(ProjectPairMassNativeRanges(
        hSparse, Form("%s_%s", name, suffix), ptMin, ptMax,
        m1Lo, m1Hi, m2Lo, m2Hi,
        false, 0.0, 0.0,
        mPairMin, mPairMax, nPairBins));
  };

  addPiece(phiMassMin, rectLo, phiMassMin, phiMassMax, "m1Low");
  addPiece(rectHi, phiMassMax, phiMassMin, phiMassMax, "m1High");
  addPiece(rectLo, rectHi, phiMassMin, rectLo, "m2Low");
  addPiece(rectLo, rectHi, rectHi, phiMassMax, "m2High");

  TH1D* result = nullptr;
  for (auto& piece : pieces) {
    if (!piece) continue;
    if (!result) {
      result = dynamic_cast<TH1D*>(piece->Clone(name));
      if (result) result->SetDirectory(nullptr);
    } else {
      result->Add(piece.get());
    }
  }
  return result;
}

TH1D* ProjectPairMassNative(THnSparseF* hSparse,
                            const char* name,
                            double ptMin,
                            double ptMax,
                            double phiMassMin,
                            double phiMassMax,
                            double dmMin,
                            double dmMax,
                            double mPairMin,
                            double mPairMax,
                            int nPairBins)
{
  if (!hSparse) return nullptr;

  const int nDim = hSparse->GetNdimensions();
  if (nDim <= std::max({kAxisMPair, kAxisPtPair, kAxisRapidity,
                        kAxisM1, kAxisM2, kAxisDeltaM})) {
    std::cerr << "ERROR: THnSparse has only " << nDim
              << " dimensions, but axes 0,1,3,4,5,6 are requested.\n";
    return nullptr;
  }

  return ProjectPairMassNativeRanges(
      hSparse, name, ptMin, ptMax,
      phiMassMin, phiMassMax, phiMassMin, phiMassMax,
      true, dmMin, dmMax,
      mPairMin, mPairMax, nPairBins);
}

TF1* FitSmoothTruePhiPhiContinuum(TH1D* hTrueYield,
                                  const char* name,
                                  double mPairMin,
                                  double mPairMax,
                                  double xMass,
                                  double xHalfWindow,
                                  int templateBkgModel)
{
  if (!hTrueYield) return nullptr;

  gContModel = templateBkgModel;
  gContNPar = 4; // maximum background parameter count; order-2 fixes p3
  gContXMin = mPairMin;
  gContXMax = mPairMax;

  if (xMass > 0.0 && xHalfWindow > 0.0) {
    gRejectMin = xMass - xHalfWindow;
    gRejectMax = xMass + xHalfWindow;
  } else {
    gRejectMin = -1.0;
    gRejectMax = -1.0;
  }

  auto* fReject = new TF1(Form("%s_reject", name), ContinuumReject,
                          mPairMin, mPairMax, gContNPar);
  fReject->SetNpx(1000);

  const double yMean = hTrueYield->Integral() / std::max(1, hTrueYield->GetNbinsX());
  ConfigureBackgroundParameters(
      fReject, 0, templateBkgModel, yMean, HistogramMaximum(hTrueYield));

  hTrueYield->Fit(fReject, "RQ0");

  auto* fNoReject = MakeContinuumInterpolationFunction(fReject, name, mPairMin, mPairMax);
  delete fReject;
  return fNoReject;
}

TH1D* BuildSmoothTemplateHistFromTF1(TH1D* hRef,
                                     TF1* f,
                                     const char* name,
                                     const char* componentLabel)
{
  if (!hRef || !f) return nullptr;
  auto* h = dynamic_cast<TH1D*>(hRef->Clone(name));
  if (!h) return nullptr;
  h->Reset("ICESM");
  h->SetDirectory(nullptr);
  h->Sumw2(false);
  h->SetTitle(Form("%s (%s);M_{#phi#phi} (GeV/#it{c}^{2});unit-area template",
                   componentLabel ? componentLabel : "background template",
                   ContinuumModelLabel()));

  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double lo = h->GetXaxis()->GetBinLowEdge(i);
    const double hi = h->GetXaxis()->GetBinUpEdge(i);
    double y = f->Integral(lo, hi);
    if (!std::isfinite(y) || y < 0.0) y = 0.0;
    h->SetBinContent(i, y);
    h->SetBinError(i, 0.0);
  }
  ScaleHistogramToUnitArea(h);
  return h;
}

double VoigtRaw(double x, double mean, double gamma, double resolution)
{
  if (!(gamma > 0.0) || !(resolution > 0.0)) return 0.0;
  // TMath::Voigt uses Gaussian sigma and the Lorentzian full width (Gamma),
  // matching RooVoigtian's parameter convention.
  return TMath::Voigt(x - mean, resolution, gamma, 4);
}

// Simple ROOT Voigtian used by both X spectrum fits.
// Parameters are always:
//   p0 = total signal yield, p1 = Lorentzian full width Gamma_X,
//   p2 = M_X, p3 = fixed Gaussian detector sigma.
// TMath::Voigt is a normalized density, therefore multiplication by the
// histogram bin width converts p0 into the expected counts per bin.
double gResidualVoigtBinWidth = 0.008;
double gDirectVoigtBinWidth = 0.008;
double gDirectFitMin = 2.5;
double gDirectFitMax = 2.9;
int gDirectBackgroundModel = 0;

double RootVoigtYieldResidual(double* x, double* p)
{
  return gResidualVoigtBinWidth * p[0] *
      TMath::Voigt(x[0] - p[2], p[3], p[1], 4);
}

double RootVoigtYieldDirect(double* x, double* p)
{
  return gDirectVoigtBinWidth * p[0] *
      TMath::Voigt(x[0] - p[2], p[3], p[1], 4);
}

double DirectBackgroundTF1(double* x, double* p)
{
  return BackgroundValue(
      gDirectBackgroundModel, x[0], gDirectFitMin, gDirectFitMax, p);
}

double DirectVoigtBackgroundTF1(double* x, double* p)
{
  return RootVoigtYieldDirect(x, p) + DirectBackgroundTF1(x, p+4);
}

double VoigtRawIntegral(double xmin,
                        double xmax,
                        double mean,
                        double gamma,
                        double resolution)
{
  if (!(xmax > xmin) || !(gamma > 0.0) || !(resolution > 0.0)) return 0.0;

  // Fixed-order Simpson integration is deterministic and sufficiently precise
  // for the narrow X signal over the 2.5--2.9 GeV/c^2 analysis range.
  constexpr int nSteps = 400;
  const double step = (xmax - xmin) / static_cast<double>(nSteps);
  double sum = VoigtRaw(xmin, mean, gamma, resolution)
             + VoigtRaw(xmax, mean, gamma, resolution);
  for (int i = 1; i < nSteps; ++i) {
    const double x = xmin + i * step;
    sum += (i % 2 ? 4.0 : 2.0) * VoigtRaw(x, mean, gamma, resolution);
  }
  return sum * step / 3.0;
}

double RootVoigtIntegralFraction(double xlo,
                                 double xhi,
                                 double mean,
                                 double gamma,
                                 double resolution)
{
  if (!(xhi > xlo)) return 0.0;
  const double integral = VoigtRawIntegral(
      xlo, xhi, mean, gamma, resolution);
  if (!std::isfinite(integral)) return 0.0;
  return std::max(0.0, std::min(1.0, integral));
}

TH1D* MakeExpectedFromTemplates(TH1D* hData,
                                TH1D* hTrueTpl,
                                TH1D* hOtherTpl,
                                double nTrue,
                                double nOther,
                                double nSig,
                                double sigMean,
                                double sigGamma,
                                double sigResolution,
                                double fitMin,
                                double fitMax,
                                const char* name,
                                int component)
{
  // component: 0 total, 1 truePhiPhi, 2 other, 3 signal
  if (!hData) return nullptr;
  auto* h = dynamic_cast<TH1D*>(hData->Clone(name));
  if (!h) return nullptr;
  h->Reset("ICESM");
  h->SetDirectory(nullptr);
  h->Sumw2(false);

  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double x = h->GetXaxis()->GetBinCenter(i);
    const double dx = h->GetXaxis()->GetBinWidth(i);
    const double yTrue  = (hTrueTpl  ? nTrue  * hTrueTpl->GetBinContent(i)  : 0.0);
    const double yOther = (hOtherTpl ? nOther * hOtherTpl->GetBinContent(i) : 0.0);
    // Use the same simple ROOT density as the fitted signal.  nSig is the
    // total Voigt yield; dx converts the normalized density to counts/bin.
    const double ySig = nSig * VoigtRaw(
        x, sigMean, sigGamma, sigResolution) * dx;

    double y = yTrue + yOther + ySig;
    if (component == 1) y = yTrue;
    if (component == 2) y = yOther;
    if (component == 3) y = ySig;

    h->SetBinContent(i, y);
    h->SetBinError(i, 0.0);
  }
  return h;
}

TH1D* MakeResidualHist(TH1D* hData, TH1D* hModel, const char* name, bool pull)
{
  if (!hData || !hModel) return nullptr;
  auto* h = dynamic_cast<TH1D*>(hData->Clone(name));
  if (!h) return nullptr;
  h->Reset("ICESM");
  h->SetDirectory(nullptr);
  h->SetTitle(pull ? ";M_{#phi#phi};(data-fit)/#sigma_{data}" : ";M_{#phi#phi};data-fit");

  for (int i = 1; i <= hData->GetNbinsX(); ++i) {
    const double d = hData->GetBinContent(i);
    const double ed = hData->GetBinError(i) > 0.0 ? hData->GetBinError(i) : std::sqrt(std::max(0.0, d));
    const double m = hModel->GetBinContent(i);
    h->SetBinContent(i, pull ? ((ed > 0.0) ? (d - m) / ed : 0.0) : (d - m));
    h->SetBinError(i, pull ? 1.0 : ed);
  }
  return h;
}

TGraph* MakeTF1Graph(TF1* f, const char* name, double xmin, double xmax, int n = 200)
{
  if (!f || xmax <= xmin) return nullptr;
  auto* g = new TGraph(n);
  g->SetName(name);
  for (int i = 0; i < n; ++i) {
    const double x = xmin + (xmax - xmin) * (i + 0.5) / n;
    g->SetPoint(i, x, f->Eval(x));
  }
  return g;
}

TH1D* MakeResidualVsTF1(TH1D* hData, TF1* f, const char* name)
{
  if (!hData || !f) return nullptr;
  auto* h = dynamic_cast<TH1D*>(hData->Clone(name));
  if (!h) return nullptr;
  h->Reset("ICESM");
  h->SetDirectory(nullptr);
  h->SetTitle(";M_{#phi#phi} (GeV/#it{c}^{2});data - fit");
  h->GetYaxis()->SetTitleOffset(0.8);
  for (int i = 1; i <= hData->GetNbinsX(); ++i) {
    const double y = hData->GetBinContent(i);
    const double e = hData->GetBinError(i) > 0.0 ? hData->GetBinError(i) : std::sqrt(std::max(0.0, y));
    const double x = hData->GetXaxis()->GetBinCenter(i);
    h->SetBinContent(i, y - f->Eval(x));
    h->SetBinError(i, e);
  }
  return h;
}

double FitNormToUnitTemplate(TH1D* hData, TH1D* hUnitTpl, double& normErr)
{
  normErr = 0.0;
  if (!hData || !hUnitTpl) return 0.0;
  double num = 0.0;
  double den = 0.0;
  for (int i = 1; i <= hData->GetNbinsX(); ++i) {
    const double d = hData->GetBinContent(i);
    const double e = hData->GetBinError(i) > 0.0 ? hData->GetBinError(i) : std::sqrt(std::max(1.0, d));
    const double t = hUnitTpl->GetBinContent(i);
    if (e <= 0.0 || t <= 0.0) continue;
    num += d * t / (e * e);
    den += t * t / (e * e);
  }
  if (den <= 0.0) return 0.0;
  normErr = std::sqrt(1.0 / den);
  return num / den;
}


void DrawFig1TrueContinuum(TH1D* hTrueYield,
                           TF1* fTrue,
                           double rejectLo,
                           double rejectHi,
                           const char* outPrefix)
{
  if (!hTrueYield || !fTrue) return;
  std::unique_ptr<TH1D> hRes(MakeResidualVsTF1(hTrueYield, fTrue, "hFig1_true_data_minus_fit"));
  auto* c = new TCanvas("cFig1_trueContinuum", "Fig.1 true continuum", 1500, 650);
  c->Divide(2,1);

  c->cd(1);
  gPad->SetLeftMargin(0.13);
  hTrueYield->SetTitle("Fig.1: true-#phi#phi continuum;M_{#phi#phi} (GeV/#it{c}^{2});Y_{SS}");
  hTrueYield->SetMarkerStyle(20);
  SetCountDisplayRange(hTrueYield);
  hTrueYield->Draw("E");
  const double xmin = hTrueYield->GetXaxis()->GetXmin();
  const double xmax = hTrueYield->GetXaxis()->GetXmax();
  std::unique_ptr<TGraph> gL(MakeTF1Graph(fTrue, "gTrue_fit_left", xmin, rejectLo));
  std::unique_ptr<TGraph> gR(MakeTF1Graph(fTrue, "gTrue_fit_right", rejectHi, xmax));
  std::unique_ptr<TGraph> gI(MakeTF1Graph(fTrue, "gTrue_interp", rejectLo, rejectHi));
  for (auto* g : {gL.get(), gR.get()}) {
    if (!g) continue;
    g->SetLineColor(kRed + 1); g->SetLineWidth(3); g->SetLineStyle(1); g->Draw("L same");
  }
  if (gI) { gI->SetLineColor(kRed + 1); gI->SetLineWidth(3); gI->SetLineStyle(2); gI->Draw("L same"); }
  auto* leg = new TLegend(0.52, 0.70, 0.88, 0.88);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->AddEntry(hTrueYield, "Y_{SS}", "lep");
  leg->AddEntry(gL.get(), Form("%s fit region", ContinuumModelLabel()), "l");
  leg->AddEntry(gI.get(), "interpolation", "l");
  leg->Draw();
  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.035);
  lat.DrawLatex(0.16, 0.84, Form("Rejected: %.2f < M < %.2f", rejectLo, rejectHi));

  c->cd(2);
  gPad->SetLeftMargin(0.13);
  gPad->SetGridy();
  hRes->SetTitle("Fig.1: true-#phi#phi residual;M_{#phi#phi} (GeV/#it{c}^{2});Y_{SS} - fit");
  hRes->SetMarkerStyle(20);
  hRes->Draw("E");
  auto* l0 = new TLine(xmin, 0.0, xmax, 0.0);
  l0->SetLineStyle(2); l0->Draw("same");

  c->SaveAs(Form("%s_Fig1_trueContinuum.png", outPrefix));
  c->SaveAs(Form("%s_Fig1_trueContinuum.pdf", outPrefix));
  delete c;
}

void DrawFig2BkgAndControl(TH1D* hNonSS,
                           TF1* fOther,
                           TH1D* hControl,
                           TH1D* hOtherTpl,
                           const char* outPrefix,
                           int selectionMode,
                           double ctrlMin,
                           double ctrlMax)
{
  if (!hNonSS || !fOther) return;
  auto* c = new TCanvas("cFig2_bkgControl", "Fig.2 bkg/control", 1500, 650);
  c->Divide(2,1);

  c->cd(1);
  gPad->SetLeftMargin(0.13);
  hNonSS->SetTitle("#phi K K + 4K background construction;M_{#phi#phi} (GeV/#it{c}^{2});N_{nonSS}");
  hNonSS->SetMarkerStyle(20);
  SetCountDisplayRange(hNonSS);
  hNonSS->Draw("E");
  std::unique_ptr<TGraph> gB(MakeTF1Graph(fOther, "gNonSS_fit", hNonSS->GetXaxis()->GetXmin(), hNonSS->GetXaxis()->GetXmax()));
  if (gB) { gB->SetLineColor(kRed+1); gB->SetLineWidth(3); gB->Draw("L same"); }
  auto* leg1 = new TLegend(0.55, 0.75, 0.88, 0.88);
  leg1->SetBorderSize(0); leg1->SetFillStyle(0);
  leg1->AddEntry(hNonSS, "#phi K K + 4K yield", "lep");
  leg1->AddEntry(gB.get(), Form("%s fit", ContinuumModelLabel()), "l");
  leg1->Draw();

  c->cd(2);
  gPad->SetLeftMargin(0.13);
  if (hControl && hOtherTpl) {
    if (selectionMode == 0)
      hControl->SetTitle(Form("Fig.2: control region %.3f < #DeltaM < %.3f fitted by non-SS template;M_{#phi#phi} (GeV/#it{c}^{2});counts", ctrlMin, ctrlMax));
    else
      hControl->SetTitle("Fig.2: rectangular sideband/complement fitted by non-SS template;M_{#phi#phi} (GeV/#it{c}^{2});counts");
    hControl->SetMarkerStyle(20);
    SetCountDisplayRange(hControl);
    hControl->Draw("E");
    double normErr = 0.0;
    const double norm = FitNormToUnitTemplate(hControl, hOtherTpl, normErr);
    std::unique_ptr<TH1D> hScaled(dynamic_cast<TH1D*>(hOtherTpl->Clone("hFig2_nonSS_template_scaled_to_control")));
    double chi2 = 0.0;
    int ndf = 0;
    if (hScaled) {
      hScaled->SetDirectory(nullptr);
      hScaled->Scale(norm);
      hScaled->SetLineColor(kRed+1);
      hScaled->SetLineWidth(3);
      SetCountDisplayRange(hControl, hScaled.get());
      hControl->Draw("E");
      hScaled->Draw("HIST same");
      hControl->Draw("E same");

      for (int ib = 1; ib <= hControl->GetNbinsX(); ++ib) {
        const double d = hControl->GetBinContent(ib);
        const double e = hControl->GetBinError(ib) > 0.0 ? hControl->GetBinError(ib) : std::sqrt(std::max(1.0, d));
        const double m = hScaled->GetBinContent(ib);
        if (e <= 0.0) continue;
        chi2 += (d - m) * (d - m) / (e * e);
        ndf++;
      }
      ndf = std::max(0, ndf - 1); // one fitted parameter: template normalization
    }
    auto* leg2 = new TLegend(0.48, 0.70, 0.88, 0.88);
    leg2->SetBorderSize(0); leg2->SetFillStyle(0);
    leg2->AddEntry(hControl, "control data", "lep");
    leg2->AddEntry(hScaled.get(), Form("non-SS template #times N, N=%.0f #pm %.0f", norm, normErr), "l");
    leg2->AddEntry((TObject*)nullptr, Form("#chi^{2}/ndf = %.1f/%d", chi2, ndf), "");
    leg2->Draw();

    std::cout << Form("Fig.2 control-template fit: norm = %.3f +/- %.3f, chi2/ndf = %.3f/%d",
                      norm, normErr, chi2, ndf) << std::endl;
  }
  c->SaveAs(Form("%s_Fig2_nonSS_control.png", outPrefix));
  c->SaveAs(Form("%s_Fig2_nonSS_control.pdf", outPrefix));
  delete c;
}


void DrawTemplateConstructionTwoPanel(TH1D* hTrueYield,
                                      TF1* fTrueSmooth,
                                      TH1D* hNonSSYield,
                                      TF1* fNonSSSmooth,
                                      double rejectLo,
                                      double rejectHi,
                                      bool rejectNonSS,
                                      const char* outPrefix)
{
  if (!hTrueYield || !fTrueSmooth || !hNonSSYield || !fNonSSSmooth) return;

  auto* c = new TCanvas("cTemplateConstructionTwoPanel",
                        "SS and non-SS template construction", 1900, 760);
  c->Divide(2, 1, 0.002, 0.002);

  // Keep every graph alive until after SaveAs. ROOT pads store pointers to the
  // drawn objects; deleting a temporary TGraph inside the panel lambda removes
  // the visible fit curve from the saved canvas.
  std::vector<std::unique_ptr<TGraph>> fitGraphs;

  auto keepGraph = [&](TGraph* graph) -> TGraph* {
    if (!graph) return nullptr;
    fitGraphs.emplace_back(graph);
    return fitGraphs.back().get();
  };

  auto drawPanel = [&](int ipad, TH1D* h, TF1* f, const char* tag,
                       const char* title, const char* dataLabel, bool reject) {
    c->cd(ipad);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.035);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.075);
    gPad->SetTicks(1, 1);

    const double xmin = h->GetXaxis()->GetXmin();
    const double xmax = h->GetXaxis()->GetXmax();
    const double rLo = std::max(xmin, rejectLo);
    const double rHi = std::min(xmax, rejectHi);
    const bool hasReject = reject && rHi > rLo;

    h->SetTitle(title);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.78);
    h->SetLineColor(kBlack);
    StyleAxis1D(h);
    SetCountDisplayRange(h);
    h->Draw("E1");

    TGraph* gLegendFit = nullptr;
    TGraph* gInterpolation = nullptr;
    if (hasReject) {
      auto* gLeft = keepGraph(MakeTF1Graph(f, Form("g%s_left", tag), xmin, rLo));
      auto* gRight = keepGraph(MakeTF1Graph(f, Form("g%s_right", tag), rHi, xmax));
      gInterpolation = keepGraph(MakeTF1Graph(f, Form("g%s_interp", tag), rLo, rHi));
      for (auto* graph : {gLeft, gRight}) {
        if (!graph) continue;
        graph->SetLineColor(kRed + 1);
        graph->SetLineWidth(3);
        graph->Draw("L same");
        if (!gLegendFit) gLegendFit = graph;
      }
      if (gInterpolation) {
        gInterpolation->SetLineColor(kRed + 1);
        gInterpolation->SetLineWidth(3);
        gInterpolation->SetLineStyle(2);
        gInterpolation->Draw("L same");
      }

      auto* lLo = new TLine(rLo, h->GetMinimum(), rLo, h->GetMaximum());
      auto* lHi = new TLine(rHi, h->GetMinimum(), rHi, h->GetMaximum());
      for (auto* line : {lLo, lHi}) {
        line->SetLineColor(kGray + 2);
        line->SetLineStyle(3);
        line->SetLineWidth(2);
        line->Draw("same");
      }
    } else {
      auto* gFull = keepGraph(MakeTF1Graph(f, Form("g%s_full", tag), xmin, xmax));
      if (gFull) {
        gFull->SetLineColor(kRed + 1);
        gFull->SetLineWidth(3);
        gFull->Draw("L same");
        gLegendFit = gFull;
      }
    }

    h->Draw("E1 same");

    auto* leg = new TLegend(0.43, 0.68, 0.92, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(h, dataLabel, "lep");
    if (hasReject) {
      if (gLegendFit)
        leg->AddEntry(gLegendFit, Form("%s fit outside window", ContinuumModelLabel()), "l");
      if (gInterpolation)
        leg->AddEntry(gInterpolation, "interpolation through excluded window", "l");
    } else if (gLegendFit) {
      leg->AddEntry(gLegendFit, Form("%s full-range fit", ContinuumModelLabel()), "l");
    }
    leg->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.034);
    lat.DrawLatex(0.17, 0.86,
                  hasReject ? Form("Excluded: %.3f < M_{#phi#phi} < %.3f", rejectLo, rejectHi)
                            : "No signal-window exclusion");
  };

  drawPanel(1, hTrueYield, fTrueSmooth, "SS",
            Form("SS template (%s);M_{#phi#phi} (GeV/#it{c}^{2});selected N_{SS}",
                 ContinuumModelLabel()),
            "extracted N_{SS}", true);
  drawPanel(2, hNonSSYield, fNonSSSmooth, "NonSS",
            Form("Non-SS template (%s);M_{#phi#phi} (GeV/#it{c}^{2});selected N_{nonSS}",
                 ContinuumModelLabel()),
            "extracted N_{nonSS}", rejectNonSS);

  c->SaveAs(Form("%s_Fig1_templateConstruction_twoPanel.png", outPrefix));
  c->SaveAs(Form("%s_Fig1_templateConstruction_twoPanel.pdf", outPrefix));
  delete c;
}


struct TwoTemplateSidebandFitResult {
  bool ok = false;

  // Full-range yields multiplying the unit-area templates.
  double nTrue = 0.0;
  double nOther = 0.0;
  double nTrueErr = 0.0;
  double nOtherErr = 0.0;
  double covTrueOther = 0.0;

  // Scale factors relative to the nominal yields extracted from the
  // coarse-bin 2D mass decomposition.  finalTemplateNormMode=1 fixes
  // both of these scale factors to exactly one.
  double nominalTrue = 0.0;
  double nominalOther = 0.0;
  double scaleTrue = 0.0;
  double scaleOther = 0.0;
  double scaleTrueErr = 0.0;
  double scaleOtherErr = 0.0;
  int normalizationMode = 0;

  double chi2 = 0.0;
  int ndf = 0;
  int nPoints = 0;
  int nPars = 0;
};

TwoTemplateSidebandFitResult FitTwoTemplateNormsSidebandWLS(
    TH1D* hData,
    TH1D* hTrueTpl,
    TH1D* hOtherTpl,
    double rejectLo,
    double rejectHi,
    int normalizationMode,
    double nominalTrueYield,
    double nominalOtherYield)
{
  // Unit-area template convention:
  //
  //   model_i = N_true*T_true_i + N_other*T_nonSS_i.
  //
  // The nominal yields are the integrals of the selected SS and non-SS
  // distributions obtained from the coarse-bin 2D mass decomposition.
  //
  // normalizationMode=0:
  //   N_true and N_other are fitted freely in the sidebands.  The reported
  //   scale factors are N/nominal.
  //
  // normalizationMode=1:
  //   scaleTrue=scaleOther=1, therefore
  //   N_true=nominalTrueYield and N_other=nominalOtherYield.
  TwoTemplateSidebandFitResult r;
  if (!hData || !hTrueTpl || !hOtherTpl) return r;

  r.normalizationMode = normalizationMode;
  r.nominalTrue = std::max(0.0, nominalTrueYield);
  r.nominalOther = std::max(0.0, nominalOtherYield);

  double a00 = 0.0;
  double a01 = 0.0;
  double a11 = 0.0;
  double b0 = 0.0;
  double b1 = 0.0;

  for (int i = 1; i <= hData->GetNbinsX(); ++i) {
    const double x = hData->GetXaxis()->GetBinCenter(i);
    if (rejectHi > rejectLo && x > rejectLo && x < rejectHi) continue;

    const double d = hData->GetBinContent(i);
    const double e =
      hData->GetBinError(i) > 0.0
        ? hData->GetBinError(i)
        : std::sqrt(std::max(1.0, d));

    if (e <= 0.0) continue;

    const double t0 = hTrueTpl->GetBinContent(i);
    const double t1 = hOtherTpl->GetBinContent(i);
    if (t0 <= 0.0 && t1 <= 0.0) continue;

    const double w = 1.0 / (e * e);
    a00 += w * t0 * t0;
    a01 += w * t0 * t1;
    a11 += w * t1 * t1;
    b0 += w * d * t0;
    b1 += w * d * t1;
    r.nPoints++;
  }

  if (normalizationMode == 1) {
    // Fix the dimensionless scale factors to one.
    r.nTrue = r.nominalTrue;
    r.nOther = r.nominalOther;
    r.nTrueErr = 0.0;
    r.nOtherErr = 0.0;
    r.covTrueOther = 0.0;
    r.scaleTrue = 1.0;
    r.scaleOther = 1.0;
    r.scaleTrueErr = 0.0;
    r.scaleOtherErr = 0.0;
    r.nPars = 0;
  } else {
    const double det = a00 * a11 - a01 * a01;
    bool twoPar = false;

    if (det > 1.0e-30 && std::isfinite(det)) {
      r.nTrue = (b0 * a11 - b1 * a01) / det;
      r.nOther = (a00 * b1 - a01 * b0) / det;

      twoPar =
        std::isfinite(r.nTrue) &&
        std::isfinite(r.nOther) &&
        r.nTrue >= 0.0 &&
        r.nOther >= 0.0;

      if (twoPar) {
        r.nTrueErr =
          std::sqrt(std::max(0.0, a11 / det));
        r.nOtherErr =
          std::sqrt(std::max(0.0, a00 / det));
        r.covTrueOther = -a01 / det;
        r.nPars = 2;
      }
    }

    // If the unconstrained two-template solution is negative, retain only
    // the better positive one-template solution.
    if (!twoPar) {
      if (a00 > 0.0 && a11 > 0.0) {
        const double n0Only = b0 / a00;
        const double n1Only = b1 / a11;
        double chi0 = 0.0;
        double chi1 = 0.0;

        for (int i = 1; i <= hData->GetNbinsX(); ++i) {
          const double x =
            hData->GetXaxis()->GetBinCenter(i);

          if (rejectHi > rejectLo &&
              x > rejectLo &&
              x < rejectHi) {
            continue;
          }

          const double d =
            hData->GetBinContent(i);

          const double e =
            hData->GetBinError(i) > 0.0
              ? hData->GetBinError(i)
              : std::sqrt(std::max(1.0, d));

          if (e <= 0.0) continue;

          const double t0 =
            hTrueTpl->GetBinContent(i);

          const double t1 =
            hOtherTpl->GetBinContent(i);

          const double m0 =
            std::max(0.0, n0Only) * t0;

          const double m1 =
            std::max(0.0, n1Only) * t1;

          chi0 += (d - m0) * (d - m0) / (e * e);
          chi1 += (d - m1) * (d - m1) / (e * e);
        }

        if (chi0 < chi1) {
          r.nTrue = std::max(0.0, n0Only);
          r.nOther = 0.0;
          r.nTrueErr = std::sqrt(1.0 / a00);
          r.nOtherErr = 0.0;
        } else {
          r.nTrue = 0.0;
          r.nOther = std::max(0.0, n1Only);
          r.nTrueErr = 0.0;
          r.nOtherErr = std::sqrt(1.0 / a11);
        }

        r.covTrueOther = 0.0;
        r.nPars = 1;
      }
    }

    r.scaleTrue =
      r.nominalTrue > 0.0
        ? r.nTrue / r.nominalTrue
        : 0.0;

    r.scaleOther =
      r.nominalOther > 0.0
        ? r.nOther / r.nominalOther
        : 0.0;

    r.scaleTrueErr =
      r.nominalTrue > 0.0
        ? r.nTrueErr / r.nominalTrue
        : 0.0;

    r.scaleOtherErr =
      r.nominalOther > 0.0
        ? r.nOtherErr / r.nominalOther
        : 0.0;
  }

  r.chi2 = 0.0;
  for (int i = 1; i <= hData->GetNbinsX(); ++i) {
    const double x =
      hData->GetXaxis()->GetBinCenter(i);

    if (rejectHi > rejectLo &&
        x > rejectLo &&
        x < rejectHi) {
      continue;
    }

    const double d =
      hData->GetBinContent(i);

    const double e =
      hData->GetBinError(i) > 0.0
        ? hData->GetBinError(i)
        : std::sqrt(std::max(1.0, d));

    if (e <= 0.0) continue;

    const double model =
      r.nTrue * hTrueTpl->GetBinContent(i) +
      r.nOther * hOtherTpl->GetBinContent(i);

    r.chi2 +=
      (d - model) * (d - model) / (e * e);
  }

  r.ndf = std::max(0, r.nPoints - r.nPars);
  r.ok =
    r.nPoints > r.nPars &&
    std::isfinite(r.chi2);

  return r;
}

TH1D* MakeFunctionHistOnDataBins(TH1D* hRef, TF1* f, const char* name)
{
  if (!hRef || !f) return nullptr;
  auto* h = dynamic_cast<TH1D*>(hRef->Clone(name));
  if (!h) return nullptr;
  h->Reset("ICESM");
  h->SetDirectory(nullptr);
  h->Sumw2(false);
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double x = h->GetXaxis()->GetBinCenter(i);
    h->SetBinContent(i, f->Eval(x));
    h->SetBinError(i, 0.0);
  }
  return h;
}

TF1* FitResidualVoigtian(TH1D* hResidual,
                         double residualFitLo,
                         double residualFitHi,
                         double fullFitMin,
                         double fullFitMax,
                         double meanInit,
                         double meanMin,
                         double meanMax,
                         double gammaInit,
                         double gammaMin,
                         double gammaMax,
                         double resolution)
{
  if (!hResidual || residualFitHi <= residualFitLo) return nullptr;
  gResidualVoigtBinWidth = hResidual->GetXaxis()->GetBinWidth(1);

  double positiveExcess = 0.0;
  double absResidual = 0.0;
  for (int i = 1; i <= hResidual->GetNbinsX(); ++i) {
    const double x = hResidual->GetXaxis()->GetBinCenter(i);
    if (x <= residualFitLo || x >= residualFitHi) continue;
    const double y = hResidual->GetBinContent(i);
    if (y > 0.0) positiveExcess += y;
    absResidual += std::abs(y);
  }
  const double nInit = std::max(1.0, positiveExcess);
  const double nMax = std::max(10.0, 3.0 * std::max(positiveExcess, absResidual));

  auto* fVoigt = new TF1("fResidualVoigt", RootVoigtYieldResidual,
                         residualFitLo, residualFitHi, 4);
  fVoigt->SetNpx(1000);
  fVoigt->SetParNames("N_{X}", "#Gamma_{X}", "M_{X}", "#sigma_{res}");
  fVoigt->SetParameters(nInit, gammaInit, meanInit, resolution);
  fVoigt->SetParLimits(0, 0.0, nMax);
  fVoigt->SetParLimits(1, 0.001, 0.2);
  fVoigt->SetParLimits(2, meanMin, meanMax);
  fVoigt->FixParameter(3, resolution);
  fVoigt->SetLineColor(kMagenta + 1);
  fVoigt->SetLineWidth(3);

  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter::SetMaxIterations(100000);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit", "Migrad");
  hResidual->Fit(fVoigt, "RQN0");
  hResidual->Fit(fVoigt, "RQSN0");
  return fVoigt;
}

void SaveBkgSubtractedFitCSV(TH1D* hData,
                             TH1D* hTrueBkg,
                             TH1D* hOtherBkg,
                             TH1D* hTotalBkg,
                             TH1D* hResidual,
                             TF1* fResidualVoigt,
                             const char* outPrefix)
{
  std::ofstream out(Form("%s_fit_components.csv", outPrefix));
  out << "M_center,M_low,M_high,data,data_err,truePhiPhi_bkg,nonSS_bkg,total_bkg,residual,residual_err,Voigtian_residual_fit\n";
  if (!hData) return;

  for (int i = 1; i <= hData->GetNbinsX(); ++i) {
    const double x = hData->GetXaxis()->GetBinCenter(i);
    out << x << ","
        << hData->GetXaxis()->GetBinLowEdge(i) << ","
        << hData->GetXaxis()->GetBinUpEdge(i) << ","
        << hData->GetBinContent(i) << ","
        << hData->GetBinError(i) << ","
        << (hTrueBkg ? hTrueBkg->GetBinContent(i) : 0.0) << ","
        << (hOtherBkg ? hOtherBkg->GetBinContent(i) : 0.0) << ","
        << (hTotalBkg ? hTotalBkg->GetBinContent(i) : 0.0) << ","
        << (hResidual ? hResidual->GetBinContent(i) : 0.0) << ","
        << (hResidual ? hResidual->GetBinError(i) : 0.0) << ","
        << (fResidualVoigt ? fResidualVoigt->Eval(x) : 0.0) << "\n";
  }
  out.close();
}

struct SHMDataObservableResult {
  bool ok = false;

  double obsMassMin = 0.0;
  double obsMassMax = 0.0;

  double nXTotal = 0.0;
  double nXTotalErr = 0.0;
  double nXWindow = 0.0;
  double nXWindowErr = 0.0;
  double signalWindowFraction = 0.0;

  double nTruePhiPhiWindow = 0.0;
  double nTruePhiPhiWindowErr = 0.0;
  double trueTemplateFraction = 0.0;

  double nNonSSWindow = 0.0;
  double nTotalBkgWindow = 0.0;
  double nDataWindow = 0.0;
  double nResidualWindow = 0.0;

  double rData = 0.0;
  double rDataErr = 0.0;
  double rSHM = -1.0;
  double dataOverSHM = -1.0;
  double dataOverSHMErr = 0.0;
};

double IntegralHistRangeFractional(const TH1D* h, double xmin, double xmax)
{
  if (!h || xmax <= xmin) return 0.0;
  double sum = 0.0;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double lo = h->GetXaxis()->GetBinLowEdge(i);
    const double hi = h->GetXaxis()->GetBinUpEdge(i);
    const double ovLo = std::max(lo, xmin);
    const double ovHi = std::min(hi, xmax);
    if (ovHi <= ovLo) continue;
    const double frac = (ovHi - ovLo) / (hi - lo);
    sum += frac * h->GetBinContent(i);
  }
  return sum;
}

double IntegralHistErrorRangeFractional(const TH1D* h, double xmin, double xmax)
{
  if (!h || xmax <= xmin) return 0.0;
  double e2 = 0.0;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double lo = h->GetXaxis()->GetBinLowEdge(i);
    const double hi = h->GetXaxis()->GetBinUpEdge(i);
    const double ovLo = std::max(lo, xmin);
    const double ovHi = std::min(hi, xmax);
    if (ovHi <= ovLo) continue;
    const double frac = (ovHi - ovLo) / (hi - lo);
    const double e = h->GetBinError(i);
    e2 += frac * frac * e * e;
  }
  return std::sqrt(std::max(0.0, e2));
}

SHMDataObservableResult ComputeSHMDataObservable(TH1D* hData,
                                                 TH1D* hTrueBkg,
                                                 TH1D* hOtherBkg,
                                                 TH1D* hTotalBkg,
                                                 TH1D* hResidual,
                                                 TH1D* hTrueUnitTemplate,
                                                 TF1* fResidualVoigt,
                                                 const TwoTemplateSidebandFitResult& bkgFit,
                                                 const char* outPrefix,
                                                 double obsLo,
                                                 double obsHi,
                                                 double fullMassLo,
                                                 double fullMassHi,
                                                 double ptMin,
                                                 double ptMax,
                                                 double rapidityMin,
                                                 double rapidityMax,
                                                 int selectionMode,
                                                 double selectionCut,
                                                 double rSHMFiducial)
{
  SHMDataObservableResult r;
  r.obsMassMin = obsLo;
  r.obsMassMax = obsHi;
  r.rSHM = rSHMFiducial;

  if (!hData || !hTrueBkg || !hOtherBkg || !hTotalBkg || !hResidual ||
      !fResidualVoigt || obsHi <= obsLo) {
    std::cerr << "WARNING: cannot compute SHM data observable; missing histograms/function or bad observable window." << std::endl;
    return r;
  }

  const double nXTotal = fResidualVoigt->GetParameter(0);
  const double eNXTotal = fResidualVoigt->GetParError(0);
  const double gX = fResidualVoigt->GetParameter(1);
  const double mX = fResidualVoigt->GetParameter(2);
  const double resolution = fResidualVoigt->GetParameter(3);
  const double signalFrac = RootVoigtIntegralFraction(
      obsLo, obsHi, mX, gX, resolution);

  r.nXTotal = nXTotal;
  r.nXTotalErr = eNXTotal;
  r.signalWindowFraction = signalFrac;
  r.nXWindow = nXTotal * signalFrac;
  r.nXWindowErr = eNXTotal * signalFrac; // mean/width covariance not included

  r.nTruePhiPhiWindow = IntegralHistRangeFractional(hTrueBkg, obsLo, obsHi);
  if (hTrueUnitTemplate) {
    r.trueTemplateFraction = IntegralHistRangeFractional(hTrueUnitTemplate, obsLo, obsHi);
    r.nTruePhiPhiWindowErr = r.trueTemplateFraction * bkgFit.nTrueErr; // template-shape uncertainty not included
  } else {
    r.nTruePhiPhiWindowErr = IntegralHistErrorRangeFractional(hTrueBkg, obsLo, obsHi);
  }

  r.nNonSSWindow = IntegralHistRangeFractional(hOtherBkg, obsLo, obsHi);
  r.nTotalBkgWindow = IntegralHistRangeFractional(hTotalBkg, obsLo, obsHi);
  r.nDataWindow = IntegralHistRangeFractional(hData, obsLo, obsHi);
  r.nResidualWindow = IntegralHistRangeFractional(hResidual, obsLo, obsHi);

  if (r.nTruePhiPhiWindow > 0.0) {
    r.rData = r.nXWindow / r.nTruePhiPhiWindow;
    const double relNX = (r.nXWindow > 0.0 ? r.nXWindowErr / r.nXWindow : 0.0);
    const double relDen = (r.nTruePhiPhiWindowErr > 0.0 ? r.nTruePhiPhiWindowErr / r.nTruePhiPhiWindow : 0.0);
    r.rDataErr = r.rData * std::sqrt(relNX * relNX + relDen * relDen);
    r.ok = true;
  }

  if (rSHMFiducial > 0.0 && r.rData > 0.0) {
    r.dataOverSHM = r.rData / rSHMFiducial;
    r.dataOverSHMErr = r.rDataErr / rSHMFiducial;
  }

  std::ofstream out(Form("%s_SHM_data_observable.csv", outPrefix));
  out << "observable,value,error,comment\n";
  out << "M_window_low," << r.obsMassMin << ",0,GeV/c2\n";
  out << "M_window_high," << r.obsMassMax << ",0,GeV/c2\n";
  out << "pt_pair_min," << ptMin << ",0,GeV/c\n";
  out << "pt_pair_max," << ptMax << ",0,GeV/c\n";
  out << "rapidity_pair_min," << rapidityMin << ",0,sparse axis 3\n";
  out << "rapidity_pair_max," << rapidityMax << ",0,sparse axis 3; upper edge excluded\n";
  out << "selectionMode," << selectionMode << ",0,0=storedDeltaM 1=rectangle\n";
  out << "selectionCut," << selectionCut << ",0,DeltaM or rectangular half width\n";
  out << "finalTemplateNormMode," << bkgFit.normalizationMode << ",0,0=free scale factors 1=fixed to one\n";
  out << "true_template_nominal_yield," << bkgFit.nominalTrue << ",0,nominal yield before final sideband rescaling\n";
  out << "true_template_scale," << bkgFit.scaleTrue << "," << bkgFit.scaleTrueErr << ",final true-phi-phi template scale\n";
  out << "nonSS_template_nominal_yield," << bkgFit.nominalOther << ",0,nominal yield before final sideband rescaling\n";
  out << "nonSS_template_scale," << bkgFit.scaleOther << "," << bkgFit.scaleOtherErr << ",final non-SS template scale\n";
  out << "X_resolution_sigma," << resolution << ",0,fixed Gaussian sigma in GeV/c2\n";
  out << "N_X_total_Voigtian," << r.nXTotal << "," << r.nXTotalErr << ",total yield parameter of normalized TMath::Voigt\n";
  out << "Voigtian_fraction_in_window," << r.signalWindowFraction << ",0,TMath::Voigt integral fraction in observable mass window\n";
  out << "N_X_window," << r.nXWindow << "," << r.nXWindowErr << ",X yield in observable mass window\n";
  out << "N_truePhiPhi_window," << r.nTruePhiPhiWindow << "," << r.nTruePhiPhiWindowErr << ",true phi-phi continuum/background in same mass window\n";
  out << "truePhiPhi_template_fraction_window," << r.trueTemplateFraction << ",0,unit-area true phi-phi template integral in window\n";
  out << "N_nonSS_window," << r.nNonSSWindow << ",0,non-SS background in same mass window\n";
  out << "N_total_bkg_window," << r.nTotalBkgWindow << ",0,total background in same mass window\n";
  out << "N_data_window," << r.nDataWindow << "," << IntegralHistErrorRangeFractional(hData, obsLo, obsHi) << ",raw selected data in same mass window\n";
  out << "N_residual_window," << r.nResidualWindow << "," << IntegralHistErrorRangeFractional(hResidual, obsLo, obsHi) << ",data minus background integral in same mass window\n";
  out << "R_data_X_over_truePhiPhi," << r.rData << "," << r.rDataErr << ",N_X_window / N_truePhiPhi_window\n";
  out << "R_SHM_fiducial_input," << r.rSHM << ",0,optional input from SHM folding macro\n";
  out << "data_over_SHM," << r.dataOverSHM << "," << r.dataOverSHMErr << ",(R_data)/(R_SHM_fiducial_input)\n";
  out.close();

  std::cout << "\n========== Data observable for SHM comparison ==========" << std::endl;
  std::cout << "Observable mass window         = [" << obsLo << ", " << obsHi << "] GeV/c^2" << std::endl;
  std::cout << "Pair rapidity window           = [" << rapidityMin << ", "
            << rapidityMax << ")" << std::endl;
  std::cout << "N_X window                    = " << r.nXWindow << " +/- " << r.nXWindowErr
            << "  (Voigtian fraction=" << r.signalWindowFraction << ")" << std::endl;
  std::cout << "N_truePhiPhi continuum window = " << r.nTruePhiPhiWindow << " +/- " << r.nTruePhiPhiWindowErr << std::endl;
  std::cout << "N_nonSS window                = " << r.nNonSSWindow << std::endl;
  std::cout << "R_data = N_X/N_truePhiPhi     = " << r.rData << " +/- " << r.rDataErr << std::endl;
  if (rSHMFiducial > 0.0) {
    std::cout << "R_SHM_fiducial input          = " << rSHMFiducial << std::endl;
    std::cout << "Data/SHM                     = " << r.dataOverSHM << " +/- " << r.dataOverSHMErr << std::endl;
  }
  std::cout << "Saved: " << outPrefix << "_SHM_data_observable.csv" << std::endl;
  std::cout << "=======================================================\n" << std::endl;

  return r;
}

void DrawFig3BkgSubtractedFit(TH1D* hData,
                              TH1D* hTrueBkg,
                              TH1D* hOtherBkg,
                              TH1D* hTotalBkg,
                              TH1D* hResidual,
                              TF1* fResidualVoigt,
                              const TwoTemplateSidebandFitResult& bkgFit,
                              const char* outPrefix,
                              double rejectLo,
                              double rejectHi)
{
  if (!hData || !hTotalBkg || !hResidual) return;

  auto* c = new TCanvas("cFig3_bkgSubtractedFit",
                        "template fit and background-subtracted spectrum",
                        1800, 780);
  c->Divide(2, 1);

  c->cd(1);
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.12);
  hData->SetTitle(Form(
      "Template background fit (%s shapes);M_{#phi#phi} (GeV/#it{c}^{2});counts",
      ContinuumModelLabel()));
  hData->SetMarkerStyle(20);
  hData->SetMarkerSize(0.8);
  hData->SetLineColor(kBlack);
  StyleAxis1D(hData);
  // Include every plotted component in the displayed range.  Using the data
  // alone hid the lower SS and non-SS templates in this panel.
  SetCountDisplayRangeMany(
      hData, {hTrueBkg, hOtherBkg, hTotalBkg});
  hData->Draw("E1");

  if (hTrueBkg) {
    hTrueBkg->SetLineColor(kBlue + 1);
    hTrueBkg->SetLineStyle(2);
    hTrueBkg->SetLineWidth(3);
    hTrueBkg->Draw("HIST same");
  }
  if (hOtherBkg) {
    hOtherBkg->SetLineColor(kGreen + 2);
    hOtherBkg->SetLineStyle(7);
    hOtherBkg->SetLineWidth(3);
    hOtherBkg->Draw("HIST same");
  }
  hTotalBkg->SetLineColor(kRed + 1);
  hTotalBkg->SetLineWidth(3);
  hTotalBkg->Draw("HIST same");
  hData->Draw("E same");

  gPad->Update();
  auto* lLo = new TLine(rejectLo, gPad->GetUymin(), rejectLo, gPad->GetUymax());
  auto* lHi = new TLine(rejectHi, gPad->GetUymin(), rejectHi, gPad->GetUymax());
  lLo->SetLineStyle(3); lHi->SetLineStyle(3);
  lLo->Draw("same"); lHi->Draw("same");

  auto* leg = new TLegend(0.49, 0.55, 0.88, 0.88);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->AddEntry(hData, "data", "lep");
  leg->AddEntry(hTotalBkg, "total background", "l");
  leg->AddEntry(hTrueBkg, "SS continuum", "l");
  leg->AddEntry(hOtherBkg, "non-SS background", "l");
  leg->AddEntry((TObject*)nullptr,
                Form("fit excludes %.3f < M < %.3f", rejectLo, rejectHi), "");
  if (bkgFit.ndf > 0)
    leg->AddEntry((TObject*)nullptr,
                  Form("#chi^{2}/ndf = %.1f/%d", bkgFit.chi2, bkgFit.ndf), "");
  leg->Draw();

  c->cd(2);
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.12);
  gPad->SetGridy();
  hResidual->SetTitle("Background-subtracted pair mass;M_{#phi#phi} (GeV/#it{c}^{2});data - background");
  hResidual->SetMarkerStyle(20);
  hResidual->SetMarkerSize(0.8);
  StyleAxis1D(hResidual);
  SetSymmetricDisplayRange(hResidual, 1.20);
  hResidual->Draw("E1");
  auto* l0 = new TLine(hResidual->GetXaxis()->GetXmin(), 0.0,
                       hResidual->GetXaxis()->GetXmax(), 0.0);
  l0->SetLineStyle(2); l0->Draw("same");

  if (fResidualVoigt) {
    fResidualVoigt->SetLineColor(kMagenta + 1);
    fResidualVoigt->SetLineWidth(3);
    fResidualVoigt->Draw("same");
    TLatex lat;
    lat.SetNDC(); lat.SetTextSize(0.035);
    const double nX = fResidualVoigt->GetParameter(0);
    const double eNX = fResidualVoigt->GetParError(0);
    lat.DrawLatex(0.16, 0.88, Form("N_{X}=%.1f #pm %.1f", nX, eNX));
    lat.DrawLatex(0.16, 0.83,
                  Form("M_{X}=%.4f #pm %.4f",
                       fResidualVoigt->GetParameter(2), fResidualVoigt->GetParError(2)));
    lat.DrawLatex(0.16, 0.78,
                  Form("#Gamma_{X}=%.4f #pm %.4f",
                       fResidualVoigt->GetParameter(1), fResidualVoigt->GetParError(1)));
    lat.DrawLatex(0.16, 0.73,
                  Form("yield/error=%.2f (not local Z)", eNX > 0.0 ? nX/eNX : 0.0));
    lat.DrawLatex(0.16, 0.68,
                  Form("#sigma_{res}=%.1f MeV/#it{c}^{2} (fixed)",
                       1000.0 * fResidualVoigt->GetParameter(3)));
  }

  c->SaveAs(Form("%s_Fig3_finalTemplateFit.png", outPrefix));
  c->SaveAs(Form("%s_Fig3_finalTemplateFit.pdf", outPrefix));
  delete c;
}


struct AsymptoticPoint {
  bool ok = false;
  int fitStatus = -999;
  double mass = 0.0;
  double gamma = 0.0;
  double nSigHat = 0.0;
  double p0 = 0.5;
  double z = 0.0;
};

AsymptoticPoint ComputeAsymptoticLocalP0(TH1D* hMass,
                                          double massHypothesis,
                                          double gammaFixed,
                                          double resolutionFixed,
                                          double fitMin,
                                          double fitMax,
                                          int pointIndex,
                                          int directBkgModel)
{
  AsymptoticPoint out;
  out.mass = massHypothesis;
  out.gamma = gammaFixed;
  if (!hMass || hMass->Integral() <= 0.0 || !(fitMax > fitMin) ||
      !(gammaFixed > 0.0) || !(resolutionFixed > 0.0)) return out;

  const TString tag = pointIndex < 0 ? "p0_best" : Form("p0_%04d", pointIndex);
  RooRealVar mass(Form("mass_%s", tag.Data()), "M_{#phi#phi}",
                  fitMin, fitMax, "GeV/c^{2}");
  mass.setRange("fitRange", fitMin, fitMax);
  RooDataHist data(Form("data_%s", tag.Data()), "binned pair-mass data",
                   RooArgList(mass), hMass);

  const int bMin = hMass->FindBin(fitMin + 1e-9);
  const int bMax = hMass->FindBin(fitMax - 1e-9);
  const double nTot = hMass->Integral(bMin, bMax);
  if (!(nTot > 0.0)) return out;

  RooRealVar mean(Form("mean_%s", tag.Data()), "M_{X}", massHypothesis);
  RooRealVar gamma(Form("gamma_%s", tag.Data()), "#Gamma_{X}", gammaFixed);
  RooRealVar resolution(Form("resolution_%s", tag.Data()), "#sigma_{res}",
                        resolutionFixed);
  mean.setConstant(true);
  gamma.setConstant(true);
  resolution.setConstant(true);
  RooVoigtian signal(Form("signal_%s", tag.Data()), "Voigtian signal",
                     mass, mean, gamma, resolution);

  const int bkgOrder = BackgroundOrder(directBkgModel);
  const bool usesThirdOrder = bkgOrder == 3;
  const bool isBernstein = directBkgModel == kBernstein2 ||
                           directBkgModel == kBernstein3;
  const bool isChebyshev = directBkgModel == kChebyshev2 ||
                           directBkgModel == kChebyshev3;
  const bool isExponential = directBkgModel == kExpPol2 ||
                             directBkgModel == kExpPol3;
  const double bInit = isBernstein ? 1.0 : 0.0;
  const double shapeMin = isBernstein ? 1.0e-4
      : (isChebyshev ? -1.0 : (isExponential ? -10.0 : -20.0));
  const double shapeMax = isBernstein ? 100.0
      : (isChebyshev ? 1.0 : (isExponential ? 10.0 : 20.0));
  RooRealVar b0(Form("b0_%s", tag.Data()), "b_{0}", 1.0);
  b0.setConstant(true);
  RooRealVar b1(Form("b1_%s", tag.Data()), "b_{1}", bInit, shapeMin, shapeMax);
  RooRealVar b2(Form("b2_%s", tag.Data()), "b_{2}", bInit, shapeMin, shapeMax);
  RooRealVar b3(Form("b3_%s", tag.Data()), "b_{3}", bInit, shapeMin, shapeMax);
  RooArgList bkgCoeffs(b1, b2, b3);
  RooArgList order2Coeffs(b1, b2);
  RooArgList bernstein2Coeffs(b0, b1, b2);
  RooArgList bernstein3Coeffs(b0, b1, b2, b3);
  RooArgList formulaArgs2(mass, b1, b2);
  RooArgList formulaArgs3(mass, b1, b2, b3);
  std::unique_ptr<RooAbsPdf> background;
  const TString scaledMass = Form(
      "((2.0*(@0-%.17g)/%.17g)-1.0)", fitMin, fitMax-fitMin);
  const TString pol2Formula = Form(
      "1+@1*(%s)+@2*(%s)*(%s)",
      scaledMass.Data(), scaledMass.Data(), scaledMass.Data());
  const TString pol3Formula = Form(
      "%s+@3*(%s)*(%s)*(%s)", pol2Formula.Data(),
      scaledMass.Data(), scaledMass.Data(), scaledMass.Data());
  const TString expPol2Formula = Form(
      "exp(@1*(%s)+@2*(%s)*(%s))",
      scaledMass.Data(), scaledMass.Data(), scaledMass.Data());
  const TString expPol3Formula = Form(
      "exp(@1*(%s)+@2*(%s)*(%s)+@3*(%s)*(%s)*(%s))",
      scaledMass.Data(), scaledMass.Data(), scaledMass.Data(),
      scaledMass.Data(), scaledMass.Data(), scaledMass.Data());
  switch (directBkgModel) {
    case kPol2:
      background.reset(new RooGenericPdf(
          Form("background_%s", tag.Data()), "pol2",
          pol2Formula.Data(), formulaArgs2));
      break;
    case kPol3:
      background.reset(new RooGenericPdf(
          Form("background_%s", tag.Data()), "pol3",
          pol3Formula.Data(), formulaArgs3));
      break;
    case kExpPol2:
      background.reset(new RooGenericPdf(
          Form("background_%s", tag.Data()), "exp(pol2)",
          expPol2Formula.Data(), formulaArgs2));
      break;
    case kExpPol3:
      background.reset(new RooGenericPdf(
          Form("background_%s", tag.Data()), "exp(pol3)",
          expPol3Formula.Data(), formulaArgs3));
      break;
    case kBernstein2:
      background.reset(new RooBernstein(
          Form("background_%s", tag.Data()), "Bernstein-2",
          mass, bernstein2Coeffs));
      break;
    case kBernstein3:
      background.reset(new RooBernstein(
          Form("background_%s", tag.Data()), "Bernstein-3",
          mass, bernstein3Coeffs));
      break;
    case kChebyshev2:
      background.reset(new RooChebychev(
          Form("background_%s", tag.Data()), "Chebyshev-2",
          mass, order2Coeffs));
      break;
    case kChebyshev3:
      background.reset(new RooChebychev(
          Form("background_%s", tag.Data()), "Chebyshev-3",
          mass, bkgCoeffs));
      break;
    default:
      background.reset(new RooGenericPdf(
          Form("background_%s", tag.Data()), "exp(pol3)",
          expPol3Formula.Data(), formulaArgs3));
      break;
  }

  RooRealVar nSig(Form("nSig_%s", tag.Data()), "N_{sig}",
                  0.02 * nTot, 0.0, 3.0 * nTot);
  RooRealVar nBkg(Form("nBkg_%s", tag.Data()), "N_{bkg}",
                  0.98 * nTot, 0.0, 3.0 * nTot);
  RooAddPdf model(Form("model_%s", tag.Data()),
                  Form("Voigtian + %s", BackgroundLabel(directBkgModel)),
                  RooArgList(signal, *background), RooArgList(nSig, nBkg));

  auto& msg = RooMsgService::instance();
  const auto oldLevel = msg.globalKillBelow();
  msg.setGlobalKillBelow(RooFit::ERROR);
  std::unique_ptr<RooFitResult> seedFit(
      model.fitTo(data, Save(true), Extended(true), Range("fitRange"),
                  Minimizer("Minuit", "Migrad"),
                  Strategy(1), PrintLevel(-1), Warnings(false), Verbose(false)));
  msg.setGlobalKillBelow(oldLevel);

  out.fitStatus = seedFit ? seedFit->status() : -999;
  out.nSigHat = std::max(0.0, nSig.getVal());

  RooWorkspace workspace(Form("workspace_%s", tag.Data()));
  workspace.import(model);
  workspace.import(data);

  auto* wPdf = workspace.pdf(model.GetName());
  auto* wData = workspace.data(data.GetName());
  auto* wMass = workspace.var(mass.GetName());
  auto* wNSig = workspace.var(nSig.GetName());
  auto* wNBkg = workspace.var(nBkg.GetName());
  auto* wB1 = workspace.var(b1.GetName());
  auto* wB2 = workspace.var(b2.GetName());
  auto* wB3 = workspace.var(b3.GetName());
  if (!wPdf || !wData || !wMass || !wNSig || !wNBkg ||
      !wB1 || !wB2 || (usesThirdOrder && !wB3)) {
    std::cerr << "WARNING: failed to build RooStats workspace at M="
              << massHypothesis << std::endl;
    return out;
  }

  RooArgSet observables(*wMass);
  RooArgSet poi(*wNSig);
  RooArgSet nuisances(*wNBkg, *wB1, *wB2);
  if (usesThirdOrder) nuisances.add(*wB3);

  RooStats::ModelConfig sbModel(Form("SPlusB_%s", tag.Data()), &workspace);
  sbModel.SetPdf(*wPdf);
  sbModel.SetObservables(observables);
  sbModel.SetParametersOfInterest(poi);
  sbModel.SetNuisanceParameters(nuisances);
  wNSig->setVal(std::max(1.0e-6, out.nSigHat));
  sbModel.SetSnapshot(poi);

  RooStats::ModelConfig bModel(Form("BOnly_%s", tag.Data()), &workspace);
  bModel.SetPdf(*wPdf);
  bModel.SetObservables(observables);
  bModel.SetParametersOfInterest(poi);
  bModel.SetNuisanceParameters(nuisances);
  wNSig->setVal(0.0);
  bModel.SetSnapshot(poi);

  RooStats::AsymptoticCalculator calculator(*wData, sbModel, bModel);
  calculator.SetOneSidedDiscovery(true);
  RooStats::AsymptoticCalculator::SetPrintLevel(-1);
  std::unique_ptr<RooStats::HypoTestResult> result(calculator.GetHypoTest());
  if (!result) return out;

  out.p0 = result->NullPValue();
  out.z = result->Significance();
  out.ok = std::isfinite(out.p0) && std::isfinite(out.z) &&
           out.p0 >= 0.0 && out.p0 <= 1.0;
  if (!out.ok) {
    out.p0 = 0.5;
    out.z = 0.0;
  }
  return out;
}

struct DirectPairFitResult {
  bool ok = false;
  int status = -999;
  double nSig = 0.0, nSigErr = 0.0;
  double nBkg = 0.0, nBkgErr = 0.0;
  double mean = 0.0, meanErr = 0.0;
  double gamma = 0.0, gammaErr = 0.0;
  double resolution = 0.0;
  double chi2Ndf = 0.0;
  double localP0 = 0.5;
  double localZ = 0.0;
  TF1* totalFunction = nullptr;
  TF1* signalFunction = nullptr;
  TF1* backgroundFunction = nullptr;
};

DirectPairFitResult DrawDirectPairMassFit(TH1D* hMass,
                                          double dmCut,
                                          double fitMin,
                                          double fitMax,
                                          double meanMin,
                                          double meanMax,
                                          double gammaMin,
                                          double gammaMax,
                                          double xResolution,
                                          int directBkgModel)
{
  DirectPairFitResult out;
  if (!hMass || hMass->Integral() <= 0.0) return out;

  gDirectVoigtBinWidth = hMass->GetXaxis()->GetBinWidth(1);
  gDirectFitMin = fitMin;
  gDirectFitMax = fitMax;
  gDirectBackgroundModel = directBkgModel;

  const int firstBin = hMass->FindBin(fitMin + 1.e-9);
  const int lastBin = hMass->FindBin(fitMax - 1.e-9);
  const double nTot = hMass->Integral(firstBin, lastBin);
  const int nBins = std::max(1, lastBin-firstBin+1);
  const double average = std::max(1.0, nTot/nBins);
  const double meanInit = 0.5*(meanMin+meanMax);
  const double gammaInit = 0.5*(gammaMin+gammaMax);
  const double fitFraction = std::max(
      1.0e-3, VoigtRawIntegral(fitMin, fitMax, meanInit,
                               gammaInit, xResolution));
  const double signalYieldInit = 0.03*nTot/fitFraction;

  auto* total = new TF1("fDirectVoigtBackground", DirectVoigtBackgroundTF1,
                        fitMin, fitMax, 8);
  total->SetNpx(2000);
  total->SetParameter(0, signalYieldInit);
  total->SetParameter(1, gammaInit);
  total->SetParameter(2, meanInit);
  total->SetParameter(3, xResolution);
  total->SetParName(0, "N_{X}");
  total->SetParName(1, "#Gamma_{X}");
  total->SetParName(2, "M_{X}");
  total->SetParName(3, "#sigma_{res}");
  total->SetParLimits(0, 0.0, std::max(10.0, 10.0*nTot));
  total->SetParLimits(1, gammaMin, gammaMax);
  total->SetParLimits(2, meanMin, meanMax);
  total->FixParameter(3, xResolution);
  ConfigureBackgroundParameters(
      total, 4, directBkgModel, average, HistogramMaximum(hMass));

  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter::SetMaxIterations(100000);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit", "Migrad");
  hMass->Fit(total, "RQN0");
TFitResultPtr fit = hMass->Fit(total, "RQSN0");

// =====================================================
// Print covariance/correlation information for M(phi phi)
// =====================================================
if (fit.Get()) {

  std::cout << "\n============================================\n";
  std::cout << " Direct M(phi phi) fit parameter correlations\n";
  std::cout << "============================================\n";

  std::cout << "Fit status       = " << fit->Status() << "\n";
  std::cout << "CovMatrixStatus  = " << fit->CovMatrixStatus() << "\n";

  const int npar = fit->NPar();

  // Parameter values
  std::cout << "\nFit parameters:\n";
  for (int i = 0; i < npar; ++i) {
    std::cout << std::setw(15) << fit->ParName(i)
              << " = "
              << std::setw(12) << fit->Parameter(i)
              << " +/- "
              << fit->ParError(i)
              << "\n";
  }

  // Correlation matrix
  TMatrixDSym corr = fit->GetCorrelationMatrix();

  std::cout << "\nCorrelation matrix:\n\n";

  std::cout << std::setw(15) << "";
  for (int j = 0; j < npar; ++j)
    std::cout << std::setw(15) << fit->ParName(j);
  std::cout << "\n";

  for (int i = 0; i < npar; ++i) {
    std::cout << std::setw(15) << fit->ParName(i);

    for (int j = 0; j < npar; ++j) {
      std::cout << std::setw(15)
                << std::fixed << std::setprecision(3)
                << corr(i,j);
    }

    std::cout << "\n";
  }

  std::cout << "============================================\n\n";
}

out.status = static_cast<int>(fit);
  out.ok = fit.Get() && out.status == 0 && fit->CovMatrixStatus() >= 2;
  out.nSig = total->GetParameter(0);
  out.nSigErr = total->GetParError(0);
  out.mean = total->GetParameter(2);
  out.meanErr = total->GetParError(2);
  out.gamma = total->GetParameter(1);
  out.gammaErr = total->GetParError(1);
  out.resolution = xResolution;
  out.chi2Ndf = total->GetNDF() > 0
      ? total->GetChisquare()/total->GetNDF() : 0.0;

  auto* signal = new TF1("fDirectVoigtSignal", RootVoigtYieldDirect,
                         fitMin, fitMax, 4);
  signal->SetParameters(total->GetParameter(0), total->GetParameter(1),
                        total->GetParameter(2), total->GetParameter(3));
  signal->SetNpx(2000);
  auto* background = new TF1("fDirectBackground", DirectBackgroundTF1,
                             fitMin, fitMax, 4);
  background->SetParameters(total->GetParameter(4), total->GetParameter(5),
                            total->GetParameter(6), total->GetParameter(7));
  background->SetNpx(2000);

  out.nBkg = background->Integral(fitMin, fitMax)/gDirectVoigtBinWidth;
  if (fit.Get()) {
    double bkgParameters[4];
    double bkgCovariance[16];
    for (int i = 0; i < 4; ++i) {
      bkgParameters[i] = total->GetParameter(i+4);
      for (int j = 0; j < 4; ++j)
        bkgCovariance[4*i+j] = fit->CovMatrix(i+4,j+4);
    }
    out.nBkgErr = background->IntegralError(
        fitMin, fitMax, bkgParameters, bkgCovariance) /
        gDirectVoigtBinWidth;
  }

  out.totalFunction = total;
  out.signalFunction = signal;
  out.backgroundFunction = background;

  if (out.ok) {
    const AsymptoticPoint local = ComputeAsymptoticLocalP0(
        hMass, out.mean, out.gamma, xResolution,
        fitMin, fitMax, -1, directBkgModel);
    if (local.ok) {
      out.localP0 = local.p0;
      out.localZ = local.z;
    }
  }

  hMass->SetTitle(Form("Direct ROOT Voigtian + %s fit;"
                       "M_{#phi#phi} (GeV/#it{c}^{2});Counts",
                       BackgroundLabel(directBkgModel)));
  hMass->SetMarkerStyle(20);
  hMass->SetMarkerSize(0.62);
  hMass->SetLineColor(kBlack);
  StyleAxis1D(hMass);
  SetCountDisplayRange(hMass);
  hMass->Draw("E1");
  total->SetLineColor(kRed + 1);
  total->SetLineWidth(3);
  background->SetLineColor(kBlue + 1);
  background->SetLineStyle(kDashed);
  background->SetLineWidth(3);
  signal->SetLineColor(kMagenta + 1);
  signal->SetLineStyle(kDotted);
  signal->SetLineWidth(3);
  total->Draw("same");
  background->Draw("same");
  signal->Draw("same");
  hMass->Draw("E1 same");

  auto* leg = new TLegend(0.48, 0.67, 0.91, 0.90);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.033);
  leg->AddEntry(hMass, Form("#DeltaM_{#phi}<%.3f data", dmCut), "lep");
  leg->AddEntry(total,
                Form("Voigtian + %s", BackgroundLabel(directBkgModel)), "l");
  leg->AddEntry(signal, "TMath::Voigt", "l");
  leg->AddEntry(background,
                Form("%s background", BackgroundLabel(directBkgModel)), "l");
  leg->Draw();

  TLatex lat;
  lat.SetNDC(); lat.SetTextSize(0.031);
  lat.DrawLatex(0.17, 0.88, Form("M_{X}=%.4f #pm %.4f", out.mean, out.meanErr));
  lat.DrawLatex(0.17, 0.83, Form("#Gamma_{X}=%.4f #pm %.4f", out.gamma, out.gammaErr));
  lat.DrawLatex(0.17, 0.78, Form("N_{sig}=%.1f #pm %.1f", out.nSig, out.nSigErr));
  lat.DrawLatex(0.17, 0.73, Form("#chi^{2}/ndf=%.2f", out.chi2Ndf));
  lat.DrawLatex(0.17, 0.68, Form("local p_{0}=%.3g  (Z=%.2f)", out.localP0, out.localZ));
  lat.DrawLatex(0.17, 0.63,
                Form("#sigma_{res}=%.1f MeV/#it{c}^{2} (fixed)",
                     1000.0 * xResolution));
  return out;
}

DirectPairFitResult DrawDirectFitThreePanel(TH2D* h2Mass,
                                            TH1D* hSignal,
                                            TH1D* hOutside,
                                            double phiMean,
                                            double ptMin,
                                            double ptMax,
                                            double rapidityMin,
                                            double rapidityMax,
                                            double dmCut,
                                            double fitMin,
                                            double fitMax,
                                            double meanMin,
                                            double meanMax,
                                            double gammaMin,
                                            double gammaMax,
                                            double xResolution,
                                            int directBkgModel,
                                            const char* outPrefix)
{
  DirectPairFitResult out;
  if (!h2Mass || !hSignal || !hOutside) return out;

  std::unique_ptr<TH1D> hOutsideScaled(
      static_cast<TH1D*>(hOutside->Clone("hPairMass_outsideDeltaM_scaled")));
  hOutsideScaled->SetDirectory(nullptr);
  const double normLo = std::max(fitMin, fitMax - 0.06);
  const int bLo = hSignal->FindBin(normLo + 1e-9);
  const int bHi = hSignal->FindBin(fitMax - 1e-9);
  const double intSig = hSignal->Integral(bLo, bHi);
  const double intOut = hOutsideScaled->Integral(bLo, bHi);
  const double outsideScale = intOut > 0.0 ? intSig / intOut : 1.0;
  hOutsideScaled->Scale(outsideScale);

  auto* c = new TCanvas("cDirectPairFitThreePanel",
                        "2D masses, DeltaM comparison, direct pair-mass fit",
                        2550, 780);
  c->Divide(3, 1, 0.002, 0.002);

  c->cd(1);
  gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.17);
  gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.075); gPad->SetTicks(1, 1);
  h2Mass->SetTitle("Two-dimensional #phi-candidate masses;M_{#phi,1} (GeV/#it{c}^{2});M_{#phi,2} (GeV/#it{c}^{2})");
  h2Mass->GetZaxis()->SetTitle("Counts");
  StyleAxis2D(h2Mass);
  h2Mass->Draw("COLZ");
  auto* signalCircle = new TEllipse(phiMean, phiMean, dmCut, dmCut);
  signalCircle->SetFillStyle(0);
  signalCircle->SetLineColor(kRed + 1);
  signalCircle->SetLineWidth(4);
  signalCircle->Draw("same");
  TLatex lat;
  lat.SetNDC(); lat.SetTextSize(0.033);
  lat.DrawLatex(0.17, 0.91,
                Form("%.1f < p_{T}^{#phi#phi} < %.1f GeV/c", ptMin, ptMax));
  lat.DrawLatex(0.17, 0.86,
                Form("%.2f < y_{#phi#phi} < %.2f", rapidityMin,
                     rapidityMax));
  lat.SetTextColor(kRed + 1);
  lat.DrawLatex(0.17, 0.81,
                Form("red circle: #DeltaM_{#phi}<%.3f", dmCut));
  lat.SetTextColor(kBlack);

  c->cd(2);
  gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.035);
  gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.075); gPad->SetTicks(1, 1);
  hSignal->SetTitle("Selected and outside-region pair mass;M_{#phi#phi} (GeV/#it{c}^{2});Counts");
  hSignal->SetMarkerStyle(20); hSignal->SetMarkerSize(0.68);
  hSignal->SetLineColor(kRed + 1); hSignal->SetMarkerColor(kRed + 1);
  hOutsideScaled->SetLineColor(kBlue + 1); hOutsideScaled->SetLineWidth(3);
  hOutsideScaled->SetMarkerColor(kBlue + 1); hOutsideScaled->SetMarkerStyle(24);
  StyleAxis1D(hSignal);
  SetCountDisplayRange(hSignal, hOutsideScaled.get());
  hSignal->Draw("E1");
  hOutsideScaled->Draw("HIST same");
  hOutsideScaled->Draw("E1 same");
  hSignal->Draw("E1 same");
  auto* leg = new TLegend(0.41, 0.70, 0.92, 0.91);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
  leg->AddEntry(hSignal, Form("#DeltaM_{#phi}<%.3f", dmCut), "lep");
  leg->AddEntry(hOutsideScaled.get(),
                Form("#DeltaM_{#phi}#geq%.3f, scaled #times %.3g", dmCut, outsideScale),
                "lep");
  leg->AddEntry((TObject*)nullptr,
                Form("normalized in %.3f < M < %.3f", normLo, fitMax), "");
  leg->Draw();

  c->cd(3);
  gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.035);
  gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.075); gPad->SetTicks(1, 1);
  out = DrawDirectPairMassFit(hSignal, dmCut, fitMin, fitMax,
                              meanMin, meanMax, gammaMin, gammaMax,
                              xResolution, directBkgModel);

  c->SaveAs(Form("%s_Fig4_directPairMass_threePanel.png", outPrefix));
  c->SaveAs(Form("%s_Fig4_directPairMass_threePanel.pdf", outPrefix));
  delete c;
  delete out.totalFunction;
  delete out.signalFunction;
  delete out.backgroundFunction;
  out.totalFunction = nullptr;
  out.signalFunction = nullptr;
  out.backgroundFunction = nullptr;

  std::ofstream csv(Form("%s_directPairMassFit.csv", outPrefix));
  csv << "background_model,status,Nsig,Nsig_err,Nbkg,Nbkg_err,mean,mean_err,gamma,gamma_err,resolution_sigma,chi2_ndf,local_p0,local_Z,outside_scale,pt_min,pt_max,rapidity_min,rapidity_max\n";
  csv << BackgroundLabel(directBkgModel) << ","
      << out.status << "," << out.nSig << "," << out.nSigErr << ","
      << out.nBkg << "," << out.nBkgErr << "," << out.mean << ","
      << out.meanErr << "," << out.gamma << "," << out.gammaErr << ","
      << out.resolution << ","
      << out.chi2Ndf << "," << out.localP0 << "," << out.localZ << ","
      << outsideScale << "," << ptMin << "," << ptMax << ","
      << rapidityMin << "," << rapidityMax << "\n";
  return out;
}

SHMDataObservableResult ComputeDirectSHMDataObservable(
    const DirectPairFitResult& xFit,
    const SelectedPosteriorYields& selectedPhiPhi,
    const Fit2DResult& phiPairFit,
    const char* outPrefix,
    double obsLo,
    double obsHi,
    double fullMassLo,
    double fullMassHi,
    double ptMin,
    double ptMax,
    double rapidityMin,
    double rapidityMax,
    double deltaMCut,
    int directBkgModel,
    double rSHMFiducial)
{
  SHMDataObservableResult r;
  r.obsMassMin = obsLo;
  r.obsMassMax = obsHi;
  r.rSHM = rSHMFiducial;

  if (!xFit.ok || !selectedPhiPhi.ok || !(obsHi > obsLo)) {
    std::cerr << "WARNING: default non-template R extraction failed: "
                 "invalid direct X fit or selected 2D phi-pair result."
              << std::endl;
    return r;
  }

  r.nXTotal = xFit.nSig;
  r.nXTotalErr = xFit.nSigErr;
  r.signalWindowFraction = RootVoigtIntegralFraction(
      obsLo, obsHi, xFit.mean, xFit.gamma, xFit.resolution);
  r.nXWindow = r.nXTotal*r.signalWindowFraction;

  // Propagate the direct-fit normalization and the change of the normalized
  // TMath::Voigt fraction under +/-1 sigma variations of M_X and Gamma_X.
  // Covariances among these three quantities are not retained by the compact
  // DirectPairFitResult, so this is a transparent diagonal approximation.
  const double fracMeanUp = RootVoigtIntegralFraction(
      obsLo, obsHi, xFit.mean+xFit.meanErr,
      xFit.gamma, xFit.resolution);
  const double fracMeanDown = RootVoigtIntegralFraction(
      obsLo, obsHi, xFit.mean-xFit.meanErr,
      xFit.gamma, xFit.resolution);
  const double fracGammaUp = RootVoigtIntegralFraction(
      obsLo, obsHi, xFit.mean,
      xFit.gamma+xFit.gammaErr, xFit.resolution);
  const double fracGammaDown = RootVoigtIntegralFraction(
      obsLo, obsHi, xFit.mean,
      std::max(1.0e-9,xFit.gamma-xFit.gammaErr), xFit.resolution);
  const double fracMeanErr = 0.5*std::abs(fracMeanUp-fracMeanDown);
  const double fracGammaErr = 0.5*std::abs(fracGammaUp-fracGammaDown);
  const double fracErr = std::hypot(fracMeanErr,fracGammaErr);
  r.nXWindowErr = std::hypot(
      r.signalWindowFraction*r.nXTotalErr,
      r.nXTotal*fracErr);

  // The full, uncut 2D mass square is fitted with BWxBW, BWxB, BxBW and
  // BxB.  The true-pair denominator is then the posterior SS sum over the
  // ACTUAL stored DeltaM-selected sparse entries in the same Mphiphi window.
  // This avoids fitting a BWxBW PDF to a distribution truncated by a circle.
  const double nSelectedSSInclusive = selectedPhiPhi.nSS;
  const double nSelectedSSInclusiveErr = std::sqrt(std::max(
      0.0,selectedPhiPhi.statVarSS+selectedPhiPhi.fitVarDiag));
  // The 2D SS component contains both independent phi-phi continuum and any
  // X->phi phi signal, because both are genuine SS pairs.  The SHM denominator
  // is the independent-pair contribution, so remove the fitted X yield in the
  // same fiducial window.  The covariance is not retained; errors are combined
  // diagonally and both inclusive and subtracted yields are written below.
  r.nTruePhiPhiWindow = nSelectedSSInclusive-r.nXWindow;
  r.nTruePhiPhiWindowErr = std::hypot(
      nSelectedSSInclusiveErr,r.nXWindowErr);
  r.nNonSSWindow = selectedPhiPhi.nNonSS;
  r.nDataWindow = selectedPhiPhi.raw;

  if (r.nTruePhiPhiWindow > 0.0) {
    r.rData = r.nXWindow/r.nTruePhiPhiWindow;
    const double relX = r.nXWindow > 0.0
      ? r.nXWindowErr/r.nXWindow : 0.0;
    const double relPair = r.nTruePhiPhiWindowErr > 0.0
      ? r.nTruePhiPhiWindowErr/r.nTruePhiPhiWindow : 0.0;
    r.rDataErr = r.rData*std::hypot(relX,relPair);
    r.ok = true;
  }
  if (rSHMFiducial > 0.0 && r.rData > 0.0) {
    r.dataOverSHM = r.rData/rSHMFiducial;
    r.dataOverSHMErr = r.rDataErr/rSHMFiducial;
  }

  std::ofstream out(Form("%s_SHM_data_observable_direct.csv",outPrefix));
  out << "observable,value,error,comment\n";
  out << "fit_mode,0,0,direct Voigtian+"
      << BackgroundLabel(directBkgModel)
      << " fit; no template fit\n";
  out << "M_window_low," << obsLo << ",0,GeV/c2\n";
  out << "M_window_high," << obsHi << ",0,GeV/c2\n";
  out << "pt_pair_min," << ptMin << ",0,GeV/c\n";
  out << "pt_pair_max," << ptMax << ",0,GeV/c\n";
  out << "rapidity_pair_min," << rapidityMin << ",0,sparse axis 3\n";
  out << "rapidity_pair_max," << rapidityMax << ",0,sparse axis 3; upper edge excluded\n";
  out << "DeltaM_max," << deltaMCut << ",0,stored DeltaM cut\n";
  out << "phi_BW_mean," << phiPairFit.mean << ","
      << phiPairFit.meanErr << ",GeV/c2\n";
  out << "phi_BW_Gamma_full_width," << phiPairFit.gamma << ","
      << phiPairFit.gammaErr << ",GeV/c2\n";
  out << "N_truePhiPhi_full2D_before_DeltaM," << phiPairFit.nSS << ","
      << phiPairFit.nSSErr << ",extended BWxBW component yield\n";
  out << "X_resolution_sigma," << xFit.resolution
      << ",0,fixed Gaussian sigma in GeV/c2\n";
  out << "N_X_total_Voigtian," << r.nXTotal << "," << r.nXTotalErr
      << ",total yield parameter of normalized TMath::Voigt\n";
  out << "Voigtian_fraction_in_window," << r.signalWindowFraction << ","
      << fracErr << ",TMath::Voigt integral fraction\n";
  out << "N_X_window," << r.nXWindow << "," << r.nXWindowErr
      << ",X yield inside the common Mphiphi window\n";
  out << "N_truePhiPhi_selected_inclusive_X," << nSelectedSSInclusive << ","
      << nSelectedSSInclusiveErr
      << ",posterior SS yield after DeltaM; contains continuum plus X\n";
  out << "N_truePhiPhi_continuum_window," << r.nTruePhiPhiWindow << ","
      << r.nTruePhiPhiWindowErr
      << ",selected SS minus fitted X; denominator of R_data\n";
  out << "N_nonSS_window," << r.nNonSSWindow << ",0,posterior non-SS yield after cut\n";
  out << "N_selected_data_window," << r.nDataWindow << ","
      << std::sqrt(std::max(0.0,selectedPhiPhi.rawVar))
      << ",raw entries after cut\n";
  out << "R_data_fiducial," << r.rData << "," << r.rDataErr
      << ",N_X_window/N_truePhiPhi_window\n";
  out << "R_SHM_fiducial_input," << r.rSHM
      << ",0,optional model input\n";
  out << "data_over_SHM," << r.dataOverSHM << ","
      << r.dataOverSHMErr << ",R_data/R_SHM\n";
  out.close();

  std::cout << "\n========== DEFAULT NON-TEMPLATE R EXTRACTION ==========" << std::endl;
  std::cout << "Signal and denominator window = [" << obsLo << ", "
            << obsHi << "] GeV/c^2" << std::endl;
  std::cout << "pT pair window                = [" << ptMin << ", "
            << ptMax << "] GeV/c" << std::endl;
  std::cout << "rapidity pair window          = [" << rapidityMin << ", "
            << rapidityMax << ")" << std::endl;
  std::cout << "stored DeltaM cut             < " << deltaMCut
            << " GeV/c^2" << std::endl;
  std::cout << "phi model                     = Breit-Wigner (never Voigtian)"
            << std::endl;
  std::cout << "Gamma_phi (full BW width)     = " << phiPairFit.gamma
            << " +/- " << phiPairFit.gammaErr << " GeV/c^2 = "
            << 1000.0*phiPairFit.gamma << " +/- "
            << 1000.0*phiPairFit.gammaErr << " MeV/c^2" << std::endl;
  std::cout << "X model                       = Voigtian, sigma_res = "
            << 1000.0*xFit.resolution << " MeV/c^2 (fixed)" << std::endl;
  std::cout << "N_X total from direct fit     = " << r.nXTotal
            << " +/- " << r.nXTotalErr << std::endl;
  std::cout << "N_X in common mass window     = " << r.nXWindow
            << " +/- " << r.nXWindowErr
            << " (Voigtian fraction=" << r.signalWindowFraction << ")" << std::endl;
  std::cout << "N_true phi-phi before DeltaM  = " << phiPairFit.nSS
            << " +/- " << phiPairFit.nSSErr
            << " (extended 2D BWxBW yield)" << std::endl;
  std::cout << "N_true phi-phi after DeltaM   = " << nSelectedSSInclusive
            << " +/- " << nSelectedSSInclusiveErr
            << " (inclusive of X)" << std::endl;
  std::cout << "N_independent phi-phi cont.   = " << r.nTruePhiPhiWindow
            << " +/- " << r.nTruePhiPhiWindowErr << std::endl;
  std::cout << "R_data^fid = N_X/N_phiphi     = " << r.rData
            << " +/- " << r.rDataErr << std::endl;
  if (rSHMFiducial > 0.0) {
    std::cout << "R_SHM^fid input               = " << rSHMFiducial << std::endl;
    std::cout << "R_data/R_SHM                 = " << r.dataOverSHM
              << " +/- " << r.dataOverSHMErr << std::endl;
  }
  std::cout << "Saved: " << outPrefix
            << "_SHM_data_observable_direct.csv" << std::endl;
  std::cout << "=======================================================\n" << std::endl;
  return r;
}

struct LocalP0ScanResult {
  std::unique_ptr<TGraph> p0Graph;
  std::unique_ptr<TGraph> zGraph;
  std::unique_ptr<TGraphAsymmErrors> p0Band;
  std::unique_ptr<TGraphAsymmErrors> zBand;
  double bestMass = 0.0;
  double bestP0 = 0.5;
  double bestZ = 0.0;
  double gamma = 0.0;
  double gammaErr = 0.0;
  double fittedMean = 0.0;
  double fittedMeanErr = 0.0;
  int bandMode = 0;
};

LocalP0ScanResult DrawAsymptoticLocalP0Scan(TH1D* hMass,
                                            double fitMin,
                                            double fitMax,
                                            double scanMin,
                                            double scanMax,
                                            double scanStep,
                                            double gammaFixed,
                                            double gammaError,
                                            double resolutionFixed,
                                            double fittedMean,
                                            double fittedMeanError,
                                            int bandMode,
                                            int directBkgModel,
                                            const char* outPrefix)
{
  // bandMode:
  //   0 = central p0/Z curves only.
  //   1 = envelope from Gamma_X = Gamma_hat +/- sigma_Gamma.
  //   2 = envelope from the 3x3 grid
  //       Gamma_X = Gamma_hat, Gamma_hat +/- sigma_Gamma and
  //       signal-location shift = 0, +/- sigma_M.
  //
  // This is a fit-parameter-variation/sensitivity envelope.  It is not a
  // frequentist confidence interval on p0.  At every tested mass, Nsig, Nbkg
  // and the selected background coefficients are still profiled by RooStats.  The fitted
  // mean itself is the scan coordinate, so sigma_M is additionally displayed
  // as a vertical band around the best-fit mass.
  LocalP0ScanResult out;
  out.gamma = gammaFixed;
  out.gammaErr = std::max(0.0, gammaError);
  out.fittedMean = fittedMean;
  out.fittedMeanErr = std::max(0.0, fittedMeanError);
  out.bandMode = std::max(0, std::min(bandMode, 2));

  if (!hMass || !(fitMax > fitMin) || !(scanMax > scanMin) ||
      !(scanStep > 0.0) || !(gammaFixed > 0.0) ||
      !(resolutionFixed > 0.0)) return out;

  scanMin = std::max(scanMin, fitMin);
  scanMax = std::min(scanMax, fitMax);
  if (!(scanMax > scanMin)) return out;

  out.p0Graph.reset(new TGraph());
  out.zGraph.reset(new TGraph());
  out.p0Graph->SetName("gLocalP0_asymptotic");
  out.zGraph->SetName("gLocalZ_asymptotic");

  const bool makeVariationBand =
      out.bandMode > 0 &&
      ((out.gammaErr > 0.0) ||
       (out.bandMode >= 2 && out.fittedMeanErr > 0.0));

  if (makeVariationBand) {
    out.p0Band.reset(new TGraphAsymmErrors());
    out.zBand.reset(new TGraphAsymmErrors());
    out.p0Band->SetName("gLocalP0_parameterVariationBand");
    out.zBand->SetName("gLocalZ_parameterVariationBand");
  }

  std::ofstream csv(Form("%s_localP0_asymptotic.csv", outPrefix));
  csv << "background_model,mass,gamma_central,gamma_error,resolution_sigma,mean_error,p0,Z,"
      << "p0_low,p0_high,Z_low,Z_high,fit_status,nsig_hat,n_valid_variations\n";

  int graphPoint = 0;
  int scanIndex = 0;
  double minimumPositiveP0 = 0.5;
  double maximumEnvelopeZ = 0.0;

  for (double massHyp = scanMin;
       massHyp <= scanMax + 0.5 * scanStep;
       massHyp += scanStep, ++scanIndex) {
    const int centralIndex = 100 * scanIndex;
    const AsymptoticPoint point = ComputeAsymptoticLocalP0(
        hMass, massHyp, gammaFixed, resolutionFixed,
        fitMin, fitMax, centralIndex,
        directBkgModel);
    if (!point.ok) {
      csv << BackgroundLabel(directBkgModel) << ","
          << massHyp << "," << gammaFixed << "," << out.gammaErr << ","
          << resolutionFixed << ","
          << out.fittedMeanErr << "," << point.p0 << "," << point.z << ","
          << point.p0 << "," << point.p0 << "," << point.z << "," << point.z
          << "," << point.fitStatus << "," << point.nSigHat << ",0\n";
      continue;
    }

    double p0Low = point.p0;
    double p0High = point.p0;
    double zLow = point.z;
    double zHigh = point.z;
    int nValidVariations = 1;

    if (makeVariationBand) {
      std::vector<double> gammaValues{gammaFixed};
      if (out.gammaErr > 0.0) {
        gammaValues.push_back(std::max(1.0e-6, gammaFixed - out.gammaErr));
        gammaValues.push_back(gammaFixed + out.gammaErr);
      }

      std::vector<double> massShifts{0.0};
      if (out.bandMode >= 2 && out.fittedMeanErr > 0.0) {
        massShifts.push_back(-out.fittedMeanErr);
        massShifts.push_back(+out.fittedMeanErr);
      }

      int variationIndex = 1;
      for (double massShift : massShifts) {
        for (double gammaValue : gammaValues) {
          const bool isCentral =
              std::abs(massShift) < 1.0e-15 &&
              std::abs(gammaValue - gammaFixed) <
                1.0e-15 * std::max(1.0, std::abs(gammaFixed));
          if (isCentral) continue;

          const double variedMass = massHyp + massShift;
          if (variedMass <= fitMin || variedMass >= fitMax) {
            ++variationIndex;
            continue;
          }

          const AsymptoticPoint varied = ComputeAsymptoticLocalP0(
              hMass, variedMass, gammaValue, resolutionFixed, fitMin, fitMax,
              centralIndex + variationIndex, directBkgModel);
          ++variationIndex;
          if (!varied.ok) continue;

          p0Low = std::min(p0Low, varied.p0);
          p0High = std::max(p0High, varied.p0);
          zLow = std::min(zLow, varied.z);
          zHigh = std::max(zHigh, varied.z);
          ++nValidVariations;
        }
      }
    }

    const double plottedP0 = std::max(1.0e-300, point.p0);
    const double plottedP0Low = std::max(1.0e-300, p0Low);
    const double plottedP0High = std::max(plottedP0, p0High);

    out.p0Graph->SetPoint(graphPoint, massHyp, plottedP0);
    out.zGraph->SetPoint(graphPoint, massHyp, point.z);

    if (out.p0Band && out.zBand) {
      out.p0Band->SetPoint(graphPoint, massHyp, plottedP0);
      out.p0Band->SetPointError(
          graphPoint, 0.0, 0.0,
          std::max(0.0, plottedP0 - plottedP0Low),
          std::max(0.0, plottedP0High - plottedP0));
      out.zBand->SetPoint(graphPoint, massHyp, point.z);
      out.zBand->SetPointError(
          graphPoint, 0.0, 0.0,
          std::max(0.0, point.z - zLow),
          std::max(0.0, zHigh - point.z));
    }

    csv << BackgroundLabel(directBkgModel) << ","
        << massHyp << "," << gammaFixed << "," << out.gammaErr << ","
        << resolutionFixed << ","
        << out.fittedMeanErr << "," << point.p0 << "," << point.z << ","
        << p0Low << "," << p0High << "," << zLow << "," << zHigh << ","
        << point.fitStatus << "," << point.nSigHat << ","
        << nValidVariations << "\n";

    ++graphPoint;
    minimumPositiveP0 = std::min(minimumPositiveP0, plottedP0Low);
    maximumEnvelopeZ = std::max(maximumEnvelopeZ, zHigh);
    if (graphPoint == 1 || point.z > out.bestZ) {
      out.bestZ = point.z;
      out.bestP0 = point.p0;
      out.bestMass = massHyp;
    }
  }
  csv.close();

  if (graphPoint == 0) {
    out.p0Graph.reset();
    out.zGraph.reset();
    out.p0Band.reset();
    out.zBand.reset();
    return out;
  }

  out.p0Graph->SetLineColor(kBlue + 1);
  out.p0Graph->SetMarkerColor(kBlue + 1);
  out.p0Graph->SetLineWidth(3);
  out.p0Graph->SetMarkerStyle(20);
  out.p0Graph->SetMarkerSize(0.72);
  out.zGraph->SetLineColor(kRed + 1);
  out.zGraph->SetMarkerColor(kRed + 1);
  out.zGraph->SetLineWidth(3);
  out.zGraph->SetMarkerStyle(20);
  out.zGraph->SetMarkerSize(0.72);

  if (out.p0Band && out.zBand) {
    out.p0Band->SetFillColorAlpha(kAzure - 9, 0.45);
    out.p0Band->SetLineColor(kAzure - 9);
    out.zBand->SetFillColorAlpha(kRed - 9, 0.40);
    out.zBand->SetLineColor(kRed - 9);
  }

  auto* c = new TCanvas("cLocalP0Asymptotic", "local p0 and significance", 1900, 760);
  c->Divide(2, 1, 0.002, 0.002);

  const double p0Min =
      std::max(1.0e-16, std::min(1.0e-2, 0.40 * minimumPositiveP0));
  const double zMax =
      1.18 * std::max(1.0, std::max(out.bestZ, maximumEnvelopeZ));

  c->cd(1);
  gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.035);
  gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.075);
  gPad->SetTicks(1, 1); gPad->SetLogy(); gPad->SetGridy();
  out.p0Graph->SetMinimum(p0Min);
  out.p0Graph->SetMaximum(0.65);
  out.p0Graph->SetTitle(
      "Asymptotic local p_{0};M_{X} hypothesis (GeV/#it{c}^{2});local p_{0}");
  out.p0Graph->Draw("AL");
  out.p0Graph->GetXaxis()->SetLimits(scanMin, scanMax);
  out.p0Graph->GetXaxis()->SetTitleSize(0.050);
  out.p0Graph->GetYaxis()->SetTitleSize(0.050);
  out.p0Graph->GetXaxis()->SetLabelSize(0.043);
  out.p0Graph->GetYaxis()->SetLabelSize(0.043);
  out.p0Graph->GetYaxis()->SetTitleOffset(1.20);
  out.p0Graph->GetXaxis()->SetNdivisions(505);

  TBox* meanBandP0 = nullptr;
  if (out.fittedMeanErr > 0.0 && out.fittedMean > fitMin &&
      out.fittedMean < fitMax) {
    meanBandP0 = new TBox(
        std::max(scanMin, out.fittedMean - out.fittedMeanErr), p0Min,
        std::min(scanMax, out.fittedMean + out.fittedMeanErr), 0.65);
    meanBandP0->SetFillColorAlpha(kGray + 1, 0.20);
    meanBandP0->SetLineColor(0);
    meanBandP0->Draw("same");
  }
  if (out.p0Band) out.p0Band->Draw("3 same");
  out.p0Graph->Draw("LP same");

  TLatex sigmaLabel;
  sigmaLabel.SetTextSize(0.029);
  sigmaLabel.SetTextColor(kGray + 2);
  for (int iz = 1; iz <= 5; ++iz) {
    const double pSigma = 0.5 * TMath::Erfc(iz / TMath::Sqrt2());
    if (pSigma < p0Min || pSigma > 0.65) continue;
    auto* line = new TLine(scanMin, pSigma, scanMax, pSigma);
    line->SetLineColor(kGray + 1);
    line->SetLineStyle(3);
    line->Draw("same");
    sigmaLabel.DrawLatex(scanMax - 0.045 * (scanMax - scanMin),
                         1.08 * pSigma, Form("%d#sigma", iz));
  }

  auto* legP0 = new TLegend(0.54, 0.72, 0.91, 0.91);
  legP0->SetBorderSize(0); legP0->SetFillStyle(0); legP0->SetTextSize(0.030);
  legP0->AddEntry(out.p0Graph.get(), "central fixed-shape scan", "lp");
  if (out.p0Band) {
    legP0->AddEntry(
        out.p0Band.get(),
        out.bandMode >= 2
          ? "fit-parameter #pm1#sigma envelope"
          : "#Gamma_{X} #pm 1#sigma envelope",
        "f");
  }
  if (meanBandP0) {
    legP0->AddEntry(meanBandP0, "direct-fit M_{X} #pm 1#sigma", "f");
  }
  legP0->Draw();

  TLatex info;
  info.SetNDC(); info.SetTextSize(0.030);
  info.DrawLatex(
      0.17, 0.88,
      Form("fixed #Gamma_{X}=%.4f #pm %.4f GeV/#it{c}^{2}",
           gammaFixed, out.gammaErr));
  info.DrawLatex(
      0.17, 0.83,
      Form("minimum central p_{0}=%.3g at %.4f", out.bestP0, out.bestMass));
  info.DrawLatex(
      0.17, 0.78,
      Form("background: %s", BackgroundLabel(directBkgModel)));
  info.DrawLatex(
      0.17, 0.73,
      Form("Voigtian #sigma_{res}=%.1f MeV/#it{c}^{2} (fixed)",
           1000.0 * resolutionFixed));

  c->cd(2);
  gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.035);
  gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.075);
  gPad->SetTicks(1, 1); gPad->SetGridy();
  out.zGraph->SetMinimum(0.0);
  out.zGraph->SetMaximum(zMax);
  out.zGraph->SetTitle(
      "Asymptotic local significance;M_{X} hypothesis (GeV/#it{c}^{2});local Z");
  out.zGraph->Draw("AL");
  out.zGraph->GetXaxis()->SetLimits(scanMin, scanMax);
  out.zGraph->GetXaxis()->SetTitleSize(0.050);
  out.zGraph->GetYaxis()->SetTitleSize(0.050);
  out.zGraph->GetXaxis()->SetLabelSize(0.043);
  out.zGraph->GetYaxis()->SetLabelSize(0.043);
  out.zGraph->GetYaxis()->SetTitleOffset(1.20);
  out.zGraph->GetXaxis()->SetNdivisions(505);

  TBox* meanBandZ = nullptr;
  if (out.fittedMeanErr > 0.0 && out.fittedMean > fitMin &&
      out.fittedMean < fitMax) {
    meanBandZ = new TBox(
        std::max(scanMin, out.fittedMean - out.fittedMeanErr), 0.0,
        std::min(scanMax, out.fittedMean + out.fittedMeanErr), zMax);
    meanBandZ->SetFillColorAlpha(kGray + 1, 0.20);
    meanBandZ->SetLineColor(0);
    meanBandZ->Draw("same");
  }
  if (out.zBand) out.zBand->Draw("3 same");
  out.zGraph->Draw("LP same");

  auto* legZ = new TLegend(0.54, 0.72, 0.91, 0.91);
  legZ->SetBorderSize(0); legZ->SetFillStyle(0); legZ->SetTextSize(0.030);
  legZ->AddEntry(out.zGraph.get(), "central fixed-shape scan", "lp");
  if (out.zBand) {
    legZ->AddEntry(
        out.zBand.get(),
        out.bandMode >= 2
          ? "fit-parameter #pm1#sigma envelope"
          : "#Gamma_{X} #pm 1#sigma envelope",
        "f");
  }
  if (meanBandZ) {
    legZ->AddEntry(meanBandZ, "direct-fit M_{X} #pm 1#sigma", "f");
  }
  legZ->Draw();

  info.DrawLatex(
      0.17, 0.88,
      Form("maximum central Z=%.2f at %.4f GeV/#it{c}^{2}",
           out.bestZ, out.bestMass));
  info.DrawLatex(0.17, 0.83, "one-sided discovery test");
  info.DrawLatex(0.17, 0.78, "band: parameter-variation envelope, not CI");

  c->SaveAs(Form("%s_Fig5_localP0_asymptotic.png", outPrefix));
  c->SaveAs(Form("%s_Fig5_localP0_asymptotic.pdf", outPrefix));
  delete c;
  return out;
}


} // namespace

namespace {

void RunDoublePhiAnalysis(
				 //const char* inputFile = "New_approach/doublephi26_new.root",
				  const char* inputFile = "AnalysisResults_25scale.root",
                                  const char* sparseName = "doublephimeson/SEMassDoublePhi",
                                  double ptMin = 8.0,
                                  double ptMax = 100.0,
                                  double phiMassMin = 1.0,
                                  double phiMassMax = 1.039,
                                  int chebOrder = 3,
                                  double mPairMin = 2.4,
                                  double mPairMax = 2.9,
                                  double templatePairBinWidth = 0.008,
                                  double finalPairBinWidth = 0.008,
                                  int selectionMode = 0,
                                  double selectionCut = 0.005,
				  double localScanStep = 0.005,
				  double normalFitDeltaM = 0.006,
				  int templateBkgModel = kPol2,
                                  double trueRejectMin = 2.65,
                                  double trueRejectMax = 2.73,
                                  double xMassMin = 2.63,
                                  double xMassMax = 2.75,
                                  double gammaMin = 0.001,
                                  double gammaMax = 0.050,
                                  int save2DFits = 1,
                                  double ctrlMin = 0.010,
                                  double ctrlMax = 0.20,
                                  int pairSignalMode = 0,
                                  int makeDebugPlots = 0,
                                  double rSHMFiducial = -1.0,
                                  int finalTemplateNormMode = 1,
                                  int pairMassFitMode = 1,
                                  double phiBWGammaInit = 0.0042,
                                  double phiBWGammaMin = 0.0010,
                                  double phiBWGammaMax = 0.0200,
                                  int templateRejectMode = 1,
                                 
                                  int makeLocalP0Scan = 1,
                                  double localScanMin = 2.55,
                                  double localScanMax = 2.85,
                                 
                                  double localScanGamma = -1.0,
                                  int localScanBandMode = 1,
                                  double localScanGammaError = -1.0,
                                  int useTemplateFit = 0,
                                  int directBkgModel = kExpPol3,
                                  double rapidityMin = 0.0,
                                  double rapidityMax = 0.8,
                                  int makeDeltaMSignificanceScan = 1,
                                  double deltaMScanMin = 0.001,
                                  double deltaMScanMax = 0.030,
                                  double deltaMScanStep = 0.001,
                                  int makePhiPurityVsPt = 1,
                                  const char* phiMassVsPtName = "doublephimeson/hPhiMassVsPt",
                                  int phiBkgModel = kExpPol3,
                                  double phiFitMin = 1.000,
                                  double phiFitMax = 1.040,
                                  double phiPurityHalfWidth = 0.007,
                                  double xResolution = 0.015)
{
  // The pair-mass binning is deliberately split into two independent controls:
  //   templatePairBinWidth for the per-Mpair 2D fits/template extraction,
  //   finalPairBinWidth    for the final selected spectrum and template/Voigtian fit.
  // selectionMode = 0: final fitted data uses stored DeltaM axis with DeltaM < selectionCut.
  // selectionMode = 1: final fitted data uses rectangular mass window |m1-mphi|<selectionCut and |m2-mphi|<selectionCut.
  // pairSignalMode controls which entries are included in the posterior sum used to extract selected Y_SS and nonSS template shapes. The 2D mass fit itself is always done in the full phi-mass range.
  const double dmSigMax = selectionCut;      // used only when selectionMode==0
  const double rectHalfWidth = selectionCut; // used only when selectionMode==1

  // Used only for Fig. 2 right-panel diagnostic/control comparison.
  // It is NOT used to build the final nonSS template or final selected-spectrum fit.
  // For selectionMode==0 this is the stored-DeltaM control region: ctrlMin < DeltaM < ctrlMax.
  // For selectionMode==1 the right panel uses the rectangular sideband/complement region instead.
  const double dmCtrlMin = ctrlMin;
  const double dmCtrlMax = ctrlMax;

  if (!(ptMax > ptMin) || !(phiMassMax > phiMassMin) ||
      !(mPairMax > mPairMin) || !(rapidityMax > rapidityMin)) {
    std::cerr << "ERROR: invalid pT, phi-mass, pair-mass, or rapidity range."
              << std::endl;
    return;
  }
  if (!(selectionCut > 0.0) || !(ctrlMax > ctrlMin) || ctrlMin < 0.0) {
    std::cerr << "ERROR: require selectionCut>0 and ctrlMax>ctrlMin>=0."
              << std::endl;
    return;
  }
  if (!(xMassMax > xMassMin) || xMassMin < mPairMin ||
      xMassMax > mPairMax || !(trueRejectMax > trueRejectMin) ||
      trueRejectMin < mPairMin || trueRejectMax > mPairMax) {
    std::cerr << "ERROR: X observable and template-rejection mass windows "
                 "must be ordered and contained inside the pair-mass range."
              << std::endl;
    return;
  }
  if (!(xResolution > 0.0)) {
    std::cerr << "ERROR: xResolution must be positive; it is the fixed "
                 "Gaussian sigma of the X Voigtian in GeV/c^2."
              << std::endl;
    return;
  }
  if (selectionMode != 0 && selectionMode != 1) {
    std::cerr << "ERROR: selectionMode must be 0 (stored DeltaM) or 1 "
                 "(rectangular phi-mass window)." << std::endl;
    return;
  }
  if (chebOrder < 1 || chebOrder > 3) {
    std::cout << "WARNING: chebOrder=" << chebOrder
              << " is outside [1,3]; clamping it." << std::endl;
    chebOrder = std::max(1, std::min(chebOrder, 3));
  }
  if (makeDeltaMSignificanceScan != 0 &&
      makeDeltaMSignificanceScan != 1) {
    std::cout << "WARNING: makeDeltaMSignificanceScan must be 0 or 1; "
                 "using 1." << std::endl;
    makeDeltaMSignificanceScan = 1;
  }
  if (!(deltaMScanStep > 0.0) || !(deltaMScanMax >= deltaMScanMin) ||
      deltaMScanMin < 0.0) {
    std::cerr << "ERROR: invalid DeltaM significance scan range or step."
              << std::endl;
    return;
  }
  if (makePhiPurityVsPt != 0 && makePhiPurityVsPt != 1) {
    std::cout << "WARNING: makePhiPurityVsPt must be 0 or 1; using 1."
              << std::endl;
    makePhiPurityVsPt = 1;
  }
  if (makePhiPurityVsPt &&
      (!(phiPurityHalfWidth > 0.0) || !phiMassVsPtName ||
       TString(phiMassVsPtName).IsNull() || !(phiFitMax > phiFitMin))) {
    std::cerr << "ERROR: phi-purity plotting requires a histogram name and "
                 "phiPurityHalfWidth>0 and phiFitMax>phiFitMin." << std::endl;
    return;
  }
  if (!IsValidBackgroundModel(phiBkgModel)) {
    std::cout << "WARNING: invalid phiBkgModel. Using exp(pol3)."
              << std::endl;
    phiBkgModel = kExpPol3;
  }
  gRapidityMin = rapidityMin;
  gRapidityMax = rapidityMax;

  // pairSignalMode controls ONLY the template signal/nonSS extraction, not the final data spectrum:
  //   0 = use the same selection as the final fit (selectionMode)
  //   1 = force stored DeltaM < selectionCut
  //   2 = force rectangular |m1-mphi|<selectionCut and |m2-mphi|<selectionCut
  //   3 = use the full phi-candidate mass square [phiMassMin,phiMassMax]^2 (no DeltaM/rect cut)
  int pairSelectionMode = selectionMode;
  if (pairSignalMode == 1) pairSelectionMode = 0;
  else if (pairSignalMode == 2) pairSelectionMode = 1;
  else if (pairSignalMode == 3) pairSelectionMode = 3;
  else if (pairSignalMode != 0) {
    std::cout << "WARNING: invalid pairSignalMode=" << pairSignalMode
              << ". Using pairSignalMode=0, same as final selection." << std::endl;
    pairSignalMode = 0;
    pairSelectionMode = selectionMode;
  }

  const double minRawFor2DFit = 50.0;

  const int nTemplatePairBins = ExactNBinsFromWidth(
    mPairMin, mPairMax, templatePairBinWidth, "templatePairBinWidth");
  const int nFinalPairBins = ExactNBinsFromWidth(
    mPairMin, mPairMax, finalPairBinWidth, "finalPairBinWidth");

  if (nTemplatePairBins <= 0 || nFinalPairBins <= 0) {
    return;
  }

  const double actualTemplatePairBinWidth =
    (mPairMax - mPairMin) / static_cast<double>(nTemplatePairBins);
  const double actualFinalPairBinWidth =
    (mPairMax - mPairMin) / static_cast<double>(nFinalPairBins);

  // Final pair-invariant-mass signal is a Voigtian with fixed resolution.
  const double xMassInit = 0.5 * (xMassMin + xMassMax);
  const double bwGammaInit = 0.5 * (gammaMin + gammaMax);

  // Breit-Wigner signal model for each single-phi mass in the 2D fits.
  gPhiBWGammaInit = phiBWGammaInit;
  gPhiBWGammaMin = phiBWGammaMin;
  gPhiBWGammaMax = phiBWGammaMax;

  gStyle->SetOptStat(0);
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter::SetMaxIterations(100000);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit", "Migrad");
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  if (!IsValidBackgroundModel(templateBkgModel)) {
    std::cout << "WARNING: invalid templateBkgModel. Using pol2."
              << std::endl;
    templateBkgModel = kPol2;
  }
  if (finalTemplateNormMode != 0 && finalTemplateNormMode != 1) {
    std::cout << "WARNING: invalid finalTemplateNormMode="
              << finalTemplateNormMode
              << ". Using 1 (fixed template scale factors)." << std::endl;
    finalTemplateNormMode = 1;
  }

  if (pairMassFitMode < 0 || pairMassFitMode > 2) {
    std::cout << "WARNING: invalid pairMassFitMode="
              << pairMassFitMode
              << ". Using 1 (separate Mpair-bin fits with global shapes fixed)."
              << std::endl;
    pairMassFitMode = 1;
  }

  if (!(gammaMax > gammaMin) || !(gammaMin > 0.0)) {
    std::cerr << "ERROR: require gammaMax > gammaMin > 0 for the intrinsic "
                 "Breit-Wigner width of the final X Voigtian."
              << std::endl;
    return;
  }
  if (!(phiBWGammaMax > phiBWGammaMin) ||
      !(phiBWGammaMin > 0.0) ||
      !(phiBWGammaInit >= phiBWGammaMin) ||
      !(phiBWGammaInit <= phiBWGammaMax)) {
    std::cerr << "ERROR: require phiBWGammaMax > phiBWGammaMin > 0 and "
              << "phiBWGammaMin <= phiBWGammaInit <= phiBWGammaMax." << std::endl;
    return;
  }
  if (templateRejectMode != 0 && templateRejectMode != 1) {
    std::cout << "WARNING: templateRejectMode must be 0 (SS only) or 1 "
                 "(SS and non-SS). Using 1.\n";
    templateRejectMode = 1;
  }
  if (!(normalFitDeltaM > 0.0)) normalFitDeltaM = selectionCut;
  if (useTemplateFit != 0 && useTemplateFit != 1) {
    std::cout << "WARNING: useTemplateFit must be 0 or 1. Using default 1.\n";
    useTemplateFit = 1;
  }
  if (!IsValidBackgroundModel(directBkgModel)) {
    std::cout << "WARNING: invalid directBkgModel. Using exp(pol3).\n";
    directBkgModel = kExpPol3;
  }
  std::cout << "Template background: "
            << BackgroundLabel(templateBkgModel) << std::endl;

  // True-phi-phi continuum-template fit excludes this fixed mass window.
  // This is intentionally decoupled from the final X mass/width limits.
  if (!(trueRejectMax > trueRejectMin)) {
    std::cerr << "ERROR: trueRejectMax must be larger than trueRejectMin." << std::endl;
    return;
  }
  const double rejectLo = trueRejectMin;
  const double rejectHi = trueRejectMax;
  const double xHalfWindow = 0.5 * (rejectHi - rejectLo);
  const double xRejectCenter = 0.5 * (rejectHi + rejectLo);

  // R_data is always evaluated in exactly [xMassMin,xMassMax].
  const double rDataMassMin = xMassMin;
  const double rDataMassMax = xMassMax;

  const TString templateBkgTag = BackgroundTag(templateBkgModel);
  const TString directBkgTag = BackgroundTag(directBkgModel);
  const TString phiBkgTag = BackgroundTag(phiBkgModel);
  TString outPrefixStr = Form(
    "DoublePhi_v7_%s_%s_cut%.3f_pair%s_fit%d_norm%d_rej%d_templateBkg%s_tpl%.4f_final%.4f_phiBW%.4f_phiBkg%s_phiRange%.3f_%.3f_xRes%.4f_directBkg%s_y%.2f_%.2f",
    useTemplateFit ? "template" : "defaultVoigt",
    SelectionLabel(selectionMode), selectionCut, SelectionLabel(pairSelectionMode),
    pairMassFitMode, finalTemplateNormMode, templateRejectMode,
    templateBkgTag.Data(),
    actualTemplatePairBinWidth, actualFinalPairBinWidth, phiBWGammaInit,
    phiBkgTag.Data(), phiFitMin, phiFitMax,
    xResolution,
    directBkgTag.Data(), rapidityMin, rapidityMax);
  const char* outPrefix = outPrefixStr.Data();
  TString fitDir = Form("%s_2DmassFits", outPrefix);
  if (save2DFits) gSystem->mkdir(fitDir.Data(), kTRUE);

  std::cout << "\n========== Requested model ==========" << std::endl;
  std::cout << "Primary extraction mode: "
            << (useTemplateFit
                  ? "template background subtraction"
                  : Form("DEFAULT direct TMath::Voigt+%s Minuit fit plus integrated 2D phi-pair fit",
                         BackgroundLabel(directBkgModel)))
            << std::endl;
  std::cout << "Single-phi signals use Breit-Wigner; X spectrum fits use TMath::Voigt with Minuit."
            << std::endl;
  std::cout << "Configurable interface: pt=[" << ptMin << "," << ptMax
            << "], rapidity=[" << rapidityMin << "," << rapidityMax << ")"
            << ", selectionMode=" << selectionMode
            << ", selectionCut=" << selectionCut
            << ", ctrl=[" << dmCtrlMin << "," << dmCtrlMax << "]"
            << ", pairSignalMode=" << pairSignalMode
            << " (extraction region=" << SelectionLabel(pairSelectionMode) << ")"
            << ", pairMassFitMode=" << pairMassFitMode
            << ", finalTemplateNormMode=" << finalTemplateNormMode
            << ", templateBkgModel=" << templateBkgModel
            << " (" << BackgroundLabel(templateBkgModel) << ")"
            << ", directBkgModel=" << directBkgModel
            << " (" << BackgroundLabel(directBkgModel) << ")"
            << ", phiBkgModel=" << phiBkgModel
            << " (" << BackgroundLabel(phiBkgModel) << ")"
            << std::endl;
  std::cout << "Single-phi mass fit/display range: ["
            << phiFitMin << "," << phiFitMax << "] GeV/c^2" << std::endl;
  std::cout << "Mass settings: phiMass=[" << phiMassMin << "," << phiMassMax
            << "], chebOrder=" << chebOrder
            << ", Mpair=[" << mPairMin << "," << mPairMax << "]" << std::endl;
  std::cout << "Template-extraction binning: requested width="
            << templatePairBinWidth
            << ", actual width=" << actualTemplatePairBinWidth
            << ", bins=" << nTemplatePairBins << std::endl;
  std::cout << "Final pair-spectrum binning: requested width="
            << finalPairBinWidth
            << ", actual width=" << actualFinalPairBinWidth
            << ", bins=" << nFinalPairBins << std::endl;
  std::cout << "True-template rejection window = [" << trueRejectMin << "," << trueRejectMax << "]" << std::endl;
  std::cout << "2D single-phi signal: Breit-Wigner, common fitted width initialized at "
            << phiBWGammaInit << " in [" << phiBWGammaMin << ","
            << phiBWGammaMax << "]" << std::endl;
  std::cout << "Template rejection mode = " << templateRejectMode
            << (templateRejectMode == 0 ? " (SS only)" : " (SS and non-SS)")
            << std::endl;
  std::cout << "Direct pair-mass fit uses stored DeltaM < " << normalFitDeltaM << std::endl;
  std::cout << "Local scan parameter band mode = " << localScanBandMode
            << " (0=none, 1=Gamma +/-1sigma, 2=Gamma and mass-location +/-1sigma)"
            << std::endl;
  std::cout << "Save all performed 2D invariant-mass fits = "
            << save2DFits
            << " (default is 1)." << std::endl;
  std::cout << "Final pair-mass signal: Voigtian with M_X=["
            << xMassMin << "," << xMassMax
            << "], Gamma_X=[" << gammaMin << "," << gammaMax << "]"
            << ", fixed sigma_res=" << xResolution << " GeV/c^2 ("
            << 1000.0*xResolution << " MeV/c^2)" << std::endl;
  std::cout << "R_data mass window is fixed to xMassMin--xMassMax: M=["
            << rDataMassMin << "," << rDataMassMax
            << "], optional R_SHM_fiducial=" << rSHMFiducial << std::endl;
  std::cout << "Final data selectionMode = " << selectionMode << " (0=stored DeltaM, 1=rectangle around fitted phi mass), rectHalfWidth=" << rectHalfWidth << std::endl;
  std::cout << "Pair-signal extraction mode = " << pairSignalMode
            << " (0=same as final, 1=DeltaM, 2=rectangle, 3=full phi-mass square), effective extraction region="
            << SelectionLabel(pairSelectionMode) << std::endl;
  std::cout << "Template 1: N_SS(Mpair) from exact selected-entry posterior sums; pairMassFitMode="
            << pairMassFitMode << ", then analytic continuum fit" << std::endl;
  std::cout << "            local 2D fits/template points use "
            << nTemplatePairBins << " bins of width "
            << actualTemplatePairBinWidth << " GeV/c^2" << std::endl;
  std::cout << "            true-phi-phi fit excludes [" << rejectLo << ", " << rejectHi << "]" << std::endl;
  std::cout << "Template 2: N_nonSS(Mpair) from the complementary exact selected-entry posterior sums; "
            << (templateRejectMode == 1
                  ? Form("fit excludes [%.3f, %.3f]", rejectLo, rejectHi)
                  : "fit uses the full pair-mass range")
            << std::endl;
  std::cout << "Final fit: data(" << SelectionLabel(selectionMode)
            << ", cut=" << selectionCut
            << ") = N_true*T_trueFixed + N_bkg*T_nonSSFixed; "
            << "template scale factors are "
            << (finalTemplateNormMode == 0 ? "free" : "fixed to 1")
            << ". Residual signal = Voigtian with free mass and intrinsic width,"
            << " fixed detector resolution."
            << std::endl;
  std::cout << "           selected data and final template/Voigtian fit use "
            << nFinalPairBins << " bins of width "
            << actualFinalPairBinWidth << " GeV/c^2" << std::endl;
  std::cout << "Fig.2 right panel: fit control histogram with fixed nonSS template shape and one free normalization. ";
  if (selectionMode == 0) std::cout << "Control = stored DeltaM in [" << dmCtrlMin << "," << dmCtrlMax << "]." << std::endl;
  else std::cout << "Control = rectangular sideband/complement region." << std::endl;
  std::cout << "=====================================\n" << std::endl;

  TFile* f = TFile::Open(inputFile);
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: cannot open input file: " << inputFile << std::endl;
    return;
  }

  auto* hSparse = dynamic_cast<THnSparseF*>(f->Get(sparseName));
  if (!hSparse) {
    std::cerr << "ERROR: cannot find THnSparseF: " << sparseName << std::endl;
    f->Close();
    return;
  }
  if (hSparse->GetNdimensions() <= kAxisDeltaM ||
      !hSparse->GetAxis(kAxisRapidity)) {
    std::cerr << "ERROR: the input THnSparse does not contain the expected "
                 "axes 0,1,3,4,5,6." << std::endl;
    f->Close();
    return;
  }
  const TAxis* rapidityAxis = hSparse->GetAxis(kAxisRapidity);
  if (rapidityMax <= rapidityAxis->GetXmin() ||
      rapidityMin >= rapidityAxis->GetXmax()) {
    std::cerr << "ERROR: requested rapidity range [" << rapidityMin << ","
              << rapidityMax << ") does not overlap sparse axis 3 ["
              << rapidityAxis->GetXmin() << "," << rapidityAxis->GetXmax()
              << "]." << std::endl;
    f->Close();
    return;
  }
  std::cout << "Applied pair-rapidity selection on sparse axis 3: "
            << rapidityMin << " <= y_{phiphi} < " << rapidityMax
            << std::endl;

  PhiPurityVsPtResult phiPurityVsPt;
  if (makePhiPurityVsPt) {
    auto* hPhiMassVsPt = dynamic_cast<TH2*>(f->Get(phiMassVsPtName));
    if (!hPhiMassVsPt) {
      std::cerr << "WARNING: cannot find 2D phi mass-versus-pT histogram: "
                << phiMassVsPtName
                << ". The phi-purity plot will be skipped." << std::endl;
    } else {
      phiPurityVsPt = DrawPhiPurityVsPt(
          hPhiMassVsPt, phiPurityHalfWidth,
          phiBkgModel, phiFitMin, phiFitMax, outPrefix);
      std::cout << "Single-phi purity versus pT: "
                << phiPurityVsPt.successfulFits << "/13 fits succeeded."
                << std::endl;
    }
  }

  // ------------------------------------------------------------
  // 1) Global 2D phi-mass fit to define stable signal/background shapes.
  // ------------------------------------------------------------
  // Global phi-mass shape fit: use all available DeltaM_phi with a common Breit-Wigner mean and width.
  // The global fit is only used to fix the shapes and non-SS composition; yields below are extracted with DeltaM_phi < dmSigMax.
  const double dmMaxForGlobalShape = hSparse->GetAxis(kAxisDeltaM)
                                   ? hSparse->GetAxis(kAxisDeltaM)->GetXmax()
                                   : std::max(dmCtrlMax, dmSigMax);
  std::unique_ptr<TH3D> h3(BuildM1M2DeltaMFromSparse(hSparse,
                                                     "h3_global_M1_M2_storedDeltaM_noDMCut",
                                                     ptMin,
                                                     ptMax,
                                                     mPairMin,
                                                     mPairMax,
                                                     phiMassMin,
                                                     phiMassMax,
                                                     dmMaxForGlobalShape));
  if (!h3 || h3->Integral() <= 0.0) {
    std::cerr << "ERROR: failed to build global M1-M2-DeltaM histogram." << std::endl;
    f->Close();
    return;
  }

  std::unique_ptr<TH2D> h2Fit(BuildM1M2ForMPairBinFromSparse(
      hSparse,
      "h2_global_M1M2_for_shape_fit",
      ptMin, ptMax,
      phiMassMin, phiMassMax,
      mPairMin, mPairMax));
  if (!h2Fit || h2Fit->Integral() <= 0.0) {
    std::cerr << "ERROR: failed to build global M1-M2 histogram." << std::endl;
    f->Close();
    return;
  }

  TString global2DPrefix = save2DFits ? Form("%s/global2D", fitDir.Data()) : TString(outPrefix);
  Fit2DResult globalFit = FitFull2DMass(h2Fit.get(), phiMassMin, phiMassMax, chebOrder, global2DPrefix.Data());
  if (!globalFit.ok) {
    std::cout << "WARNING: global 2D shape fit did not return status 0. Continue, but inspect plots." << std::endl;
  }

  // ------------------------------------------------------------
  // 2) Fine-binned data histogram in the signal DeltaM region.
  //    This histogram defines the binning used by the final template fit and Voigtian fit.
  //    A DeltaM-control histogram is saved only as a diagnostic.
  // ------------------------------------------------------------
  std::unique_ptr<TH1D> hData(ProjectPairMassNativeSelection(hSparse,
                                                    Form("hData_%s_signalRegion", SelectionLabel(selectionMode)),
                                                    ptMin, ptMax,
                                                    phiMassMin, phiMassMax,
                                                    dmSigMax,
                                                    globalFit.mean,
                                                    selectionMode,
                                                    rectHalfWidth,
                                                    mPairMin, mPairMax,
                                                    nFinalPairBins));

  std::unique_ptr<TH1D> hOtherRaw;
  if (selectionMode == 0) {
    hOtherRaw.reset(ProjectPairMassNative(hSparse,
                                          "hOtherBkg_dmControlRegion_raw",
                                          ptMin, ptMax,
                                          phiMassMin, phiMassMax,
                                          dmCtrlMin, dmCtrlMax,
                                          mPairMin, mPairMax,
                                          nFinalPairBins));
  } else {
    hOtherRaw.reset(ProjectPairMassNativeSelection(hSparse,
                                          "hOtherBkg_rectSideband_raw",
                                          ptMin, ptMax,
                                          phiMassMin, phiMassMax,
                                          dmSigMax,
                                          globalFit.mean,
                                          2,
                                          rectHalfWidth,
                                          mPairMin, mPairMax,
                                          nFinalPairBins));
  }

  if (!hData || hData->Integral() <= 0.0) {
    std::cerr << "ERROR: empty DeltaM signal-region data histogram." << std::endl;
    f->Close();
    return;
  }

  SetNiceHist(hData.get());
  if (hOtherRaw) SetNiceHist(hOtherRaw.get());

  // Independent direct pair-mass fit requested for the fixed stored-DeltaM signal region.
  std::unique_ptr<TH1D> hDirectSignal(ProjectPairMassNative(
      hSparse, "hPairMass_direct_deltaM_signal",
      ptMin, ptMax, phiMassMin, phiMassMax,
      0.0, normalFitDeltaM, mPairMin, mPairMax, nFinalPairBins));
  std::unique_ptr<TH1D> hDirectOutside(ProjectPairMassNative(
      hSparse, "hPairMass_direct_deltaM_outside",
      ptMin, ptMax, phiMassMin, phiMassMax,
      normalFitDeltaM, dmMaxForGlobalShape,
      mPairMin, mPairMax, nFinalPairBins));

  // Default R denominator: in the same Mphiphi and pT window used for R_data,
  // fit the complete 2D phi-candidate mass square before imposing DeltaM.
  // Applying the circular cut before an ordinary rectangular BWxBW fit would
  // truncate all four component PDFs and bias their extended yields.
  std::unique_ptr<TH2D> h2DirectPhiPairWindow(
      BuildM1M2ForMPairBinFromSparse(
          hSparse,
          "h2_directR_fullDeltaM_M1M2_inXWindow",
          ptMin, ptMax,
          phiMassMin, phiMassMax,
          rDataMassMin, rDataMassMax));
  Fit2DResult directPhiPairFit;
  SelectedPosteriorYields directSelectedPhiPair;
  if (h2DirectPhiPairWindow && h2DirectPhiPairWindow->Integral() > 0.0) {
    const TString direct2DPrefix = Form("%s_directDenominator2D",outPrefix);
    directPhiPairFit = FitFull2DMass(
        h2DirectPhiPairWindow.get(),phiMassMin,phiMassMax,chebOrder,
        direct2DPrefix.Data());
    directSelectedPhiPair = ExtractSelectedPosteriorYields(
        hSparse,directPhiPairFit,
        ptMin,ptMax,
        phiMassMin,phiMassMax,
        normalFitDeltaM,directPhiPairFit.mean,
        0,normalFitDeltaM,
        rDataMassMin,rDataMassMax);
    if (!directSelectedPhiPair.ok) {
      std::cout << "WARNING: exact DeltaM-selected posterior true-pair "
                   "extraction failed."
                << std::endl;
    }
  } else {
    std::cout << "WARNING: the 2D phi-pair denominator histogram in the "
                 "R_data mass window is empty."
              << std::endl;
  }

  DeltaMSignificanceScanResult deltaMSignificanceScan;
  // Optimize the phi-pair DeltaM selection over the complete analysed pair-
  // mass range.  Do not tune this pair-quality requirement in the X window.
  if (makeDeltaMSignificanceScan && h2Fit && globalFit.ok) {
    deltaMSignificanceScan = DrawDeltaMTruePairSignificance(
        h2Fit.get(), globalFit,
        deltaMScanMin, deltaMScanMax, deltaMScanStep,
        normalFitDeltaM,
        ptMin, ptMax, rapidityMin, rapidityMax,
        mPairMin, mPairMax, outPrefix);
    std::cout << "Full-Mphiphi-range DeltaM scan: maximum S/sqrt(S+B)="
              << deltaMSignificanceScan.bestSignificance
              << " at radial DeltaM<"
              << deltaMSignificanceScan.bestCut << " GeV/c^2; at the analysis "
              << "cut " << normalFitDeltaM << " it is "
              << deltaMSignificanceScan.significanceAtSelection
              << ". Maximum significance-times-purity proxy="
              << deltaMSignificanceScan.bestPurityWeightedSignificance
              << " at radial DeltaM<"
              << deltaMSignificanceScan.bestPurityWeightedCut << " GeV/c^2"
              << std::endl;
  }

  DirectPairFitResult directFit;
  if (hDirectSignal && hDirectSignal->Integral() > 0.0 &&
      hDirectOutside && hDirectOutside->Integral() > 0.0) {
    directFit = DrawDirectFitThreePanel(
        h2DirectPhiPairWindow ? h2DirectPhiPairWindow.get() : h2Fit.get(),
        hDirectSignal.get(), hDirectOutside.get(),
        directPhiPairFit.ok ? directPhiPairFit.mean : globalFit.mean,
        ptMin, ptMax, rapidityMin, rapidityMax, normalFitDeltaM,
        mPairMin, mPairMax, xMassMin, xMassMax, gammaMin, gammaMax,
        xResolution, directBkgModel, outPrefix);
  } else {
    std::cout << "WARNING: direct pair-mass fit or outside-region comparison is empty.\n";
  }

  LocalP0ScanResult localP0Scan;
  if (makeLocalP0Scan && directFit.ok && hDirectSignal &&
      hDirectSignal->Integral() > 0.0) {
    localScanMin = std::max(localScanMin, mPairMin);
    localScanMax = std::min(localScanMax, mPairMax);
    if (!(localScanMax > localScanMin)) {
      localScanMin = mPairMin;
      localScanMax = mPairMax;
    }
    if (!(localScanStep > 0.0)) localScanStep = 0.005;
    const double scanGamma = localScanGamma > 0.0
                               ? localScanGamma
                               : (directFit.gamma > 0.0 ? directFit.gamma
                                                        : 0.5 * (gammaMin + gammaMax));
    const double scanGammaError =
        localScanGammaError >= 0.0
          ? localScanGammaError
          : ((localScanGamma <= 0.0 && directFit.gammaErr > 0.0)
               ? directFit.gammaErr
               : 0.0);
    localP0Scan = DrawAsymptoticLocalP0Scan(
        hDirectSignal.get(), mPairMin, mPairMax,
        localScanMin, localScanMax, localScanStep,
        scanGamma, scanGammaError, xResolution,
        directFit.mean, directFit.meanErr,
        localScanBandMode, directBkgModel, outPrefix);
    if (localP0Scan.p0Graph) {
      std::cout << "Asymptotic local scan: minimum p0=" << localP0Scan.bestP0
                << ", maximum Z=" << localP0Scan.bestZ
                << " at M=" << localP0Scan.bestMass
                << " GeV/c^2, Gamma=" << localP0Scan.gamma
                << " +/- " << localP0Scan.gammaErr
                << " GeV/c^2, bandMode=" << localP0Scan.bandMode
                << std::endl;
    }
  } else if (makeLocalP0Scan && !directFit.ok) {
    std::cout << "WARNING: local-p0 scan skipped because the direct pair-mass "
                 "fit did not converge." << std::endl;
  }

  // This is the requested default observable.  It is independent of the
  // optional template-background construction below.
  const SHMDataObservableResult directSHMDataObs =
      ComputeDirectSHMDataObservable(
          directFit,directSelectedPhiPair,directPhiPairFit,
          outPrefix,
          rDataMassMin,rDataMassMax,
          mPairMin,mPairMax,
          ptMin,ptMax,
          rapidityMin,rapidityMax,
          normalFitDeltaM,directBkgModel,rSHMFiducial);

  if (!useTemplateFit) {
    TFile fout(Form("%s_output.root",outPrefix),"RECREATE");
    h3->Write();
    h2Fit->Write();
    hData->Write();
    if (hOtherRaw) hOtherRaw->Write();
    if (hDirectSignal) hDirectSignal->Write();
    if (hDirectOutside) hDirectOutside->Write();
    if (h2DirectPhiPairWindow) h2DirectPhiPairWindow->Write();
    if (phiPurityVsPt.purityGraph && phiPurityVsPt.successfulFits > 0)
      phiPurityVsPt.purityGraph->Write();
    if (deltaMSignificanceScan.significanceGraph)
      deltaMSignificanceScan.significanceGraph->Write();
    if (deltaMSignificanceScan.purityGraph)
      deltaMSignificanceScan.purityGraph->Write();
    if (deltaMSignificanceScan.purityWeightedSignificanceGraph)
      deltaMSignificanceScan.purityWeightedSignificanceGraph->Write();
    if (deltaMSignificanceScan.signalGraph)
      deltaMSignificanceScan.signalGraph->Write();
    if (deltaMSignificanceScan.backgroundGraph)
      deltaMSignificanceScan.backgroundGraph->Write();
    if (localP0Scan.p0Graph)
      localP0Scan.p0Graph->Write("gLocalP0_asymptotic");
    if (localP0Scan.zGraph)
      localP0Scan.zGraph->Write("gLocalZ_asymptotic");
    if (localP0Scan.p0Band)
      localP0Scan.p0Band->Write("gLocalP0_parameterVariationBand");
    if (localP0Scan.zBand)
      localP0Scan.zBand->Write("gLocalZ_parameterVariationBand");
    TNamed directDefinition(
      "default_non_template_R_definition",
      Form(
        "No template fit. X: direct TMath::Voigt+%s Minuit fit to stored DeltaM<%.9g; "
        "NX in M=[%.9g,%.9g] is obtained from the normalized Voigt integral. "
        "Denominator: full-square 2D BW/B decomposition in the same M and pT "
        "window followed by exact posterior SS sum over stored DeltaM<%.9g; "
        "the fitted X contribution is subtracted from this inclusive SS yield "
        "to obtain the independent phi-phi continuum denominator. "
        "Gamma_phi=%.9g +/- %.9g GeV is the full RooBreitWigner width. "
        "R_data_fid=%.9g +/- %.9g; optional R_SHM_input=%.9g; "
        "R_data/R_SHM=%.9g +/- %.9g. X resolution sigma=%.9g GeV/c2 (fixed).",
        BackgroundLabel(directBkgModel),
        normalFitDeltaM,rDataMassMin,rDataMassMax,normalFitDeltaM,
        directPhiPairFit.gamma,directPhiPairFit.gammaErr,
        directSHMDataObs.rData,directSHMDataObs.rDataErr,
        directSHMDataObs.rSHM,directSHMDataObs.dataOverSHM,
        directSHMDataObs.dataOverSHMErr,xResolution));
    directDefinition.Write();
    TNamed selectionDefinition(
      "analysis_selection",
      Form("pT=[%.9g,%.9g), rapidity(axis3)=[%.9g,%.9g), "
           "Mphiphi=[%.9g,%.9g), storedDeltaM<%.9g",
           ptMin,ptMax,rapidityMin,rapidityMax,
           mPairMin,mPairMax,normalFitDeltaM));
    selectionDefinition.Write();
    if (phiPurityVsPt.purityGraph && phiPurityVsPt.successfulFits > 0) {
      TNamed phiPurityDefinition(
        "phi_purity_vs_pt_definition",
        Form("Input=%s; pT bins=[0.5,0.8,1.2,1.6,2.0,2.5,3.0,4.0,"
             "5.0,6.0,8.0,10.0,20.0,50.0] GeV/c; model=Breit-Wigner+"
             "%s; fit/display range=[%.9g,%.9g] GeV/c2; "
             "purity=S/(S+B) in fitted mean +/- %.9g GeV/c2; "
             "successful fits=%d/13.",
             phiMassVsPtName,BackgroundLabel(phiBkgModel),
             phiFitMin,phiFitMax,phiPurityHalfWidth,
             phiPurityVsPt.successfulFits));
      phiPurityDefinition.Write();
    }
    if (deltaMSignificanceScan.significanceGraph) {
      TNamed deltaMScanDefinition(
        "deltaM_pair_significance_definition",
        Form("Radial DeltaM in the fitted full-Mphiphi 2D mass plane "
             "Mphiphi=[%.9g,%.9g); "
             "S=sum(Nbin*P_SS), B=sum(Nbin*(1-P_SS)), "
             "significance=S/sqrt(S+B), purity=S/B, "
             "purity-weighted significance proxy="
             "[S/sqrt(S+B)]*[S/(S+B)]. "
             "Best significance cut=%.9g with Z=%.9g; "
             "best proxy cut=%.9g with proxy=%.9g.",
             mPairMin,mPairMax,
             deltaMSignificanceScan.bestCut,
             deltaMSignificanceScan.bestSignificance,
             deltaMSignificanceScan.bestPurityWeightedCut,
             deltaMSignificanceScan.bestPurityWeightedSignificance));
      deltaMScanDefinition.Write();
    }
    fout.Close();

    std::cout << "Saved default non-template outputs with prefix: "
              << outPrefix << std::endl;
    std::cout << "  " << outPrefix
              << "_Fig4_directPairMass_threePanel.png/pdf\n";
    std::cout << "  " << outPrefix
              << "_directDenominator2D_2Dfit_projection_checks.png/pdf\n";
    std::cout << "  " << outPrefix
              << "_directDenominator2D_2Dfit_data.png/pdf\n";
    if (makeDeltaMSignificanceScan) {
      std::cout << "  " << outPrefix
                << "_Fig6_deltaM_significancePurity.png/pdf\n";
      std::cout << "  " << outPrefix
                << "_Fig7_deltaM_significanceTimesPurity.png/pdf\n";
      std::cout << "  " << outPrefix
                << "_deltaM_pairSignificance.csv\n";
    }
    if (phiPurityVsPt.successfulFits > 0) {
      std::cout << "  " << outPrefix << "_Fig8_phiPurity_vsPt.png/pdf\n";
      std::cout << "  " << outPrefix << "_phiPurity_vsPt_fits.png/pdf\n";
      std::cout << "  " << outPrefix << "_phiPurity_vsPt.csv\n";
    }
    std::cout << "  " << outPrefix
              << "_SHM_data_observable_direct.csv\n";
    std::cout << "  " << outPrefix << "_output.root\n";
    f->Close();
    return;
  }

  // ------------------------------------------------------------
  // 3) Coarse-binned pair-yield template distributions.
  //    One full 2D phi-mass fit is performed in every coarse Mpair bin.
  //    These can use the final selection, forced DeltaM selection, forced rectangular
  //    selection, or the full phi-candidate mass square, controlled by pairSignalMode.
  //    The final data spectrum itself still uses selectionMode.
  // ------------------------------------------------------------
  std::cout << "\nTemplate signal/nonSS extraction: using " << SelectionLabel(pairSelectionMode)
            << " with selectionCut=" << selectionCut << "\n"
            << "The final fitted data spectrum still uses " << SelectionLabel(selectionMode) << ".\n"
            << "The 2D mass shape parameters are fixed from the global full-DeltaM fit.\n"
            << "After each per-Mpair full-range fit, the exact selected sparse entries are decomposed using local P_SS and 1-P_SS weights.\n"
            << "For pairSignalMode=3, the posterior sum is performed over the full phi-candidate mass square.\n"
            << std::endl;

  TH1D* hNonSSRaw = nullptr;
  TH1D* hRawForSS = nullptr;
  TH1D* hFitStatus = nullptr;
  std::unique_ptr<TH1D> hTrueYield(BuildSSAndNonSSYieldsVsMPair_PosteriorSelected(
                                                     hSparse,
                                                     globalFit,
                                                     "hTruePhiPhiYield_pairTemplate_posteriorSelected",
                                                     ptMin, ptMax,
                                                     phiMassMin, phiMassMax,
                                                     chebOrder,
                                                     dmSigMax,
                                                     mPairMin, mPairMax,
                                                     nTemplatePairBins,
                                                     minRawFor2DFit,
                                                     hNonSSRaw,
                                                     hRawForSS,
                                                     hFitStatus,
                                                     pairSelectionMode,
                                                     rectHalfWidth,
                                                     pairMassFitMode,
                                                     save2DFits ? fitDir.Data() : ""));

  std::unique_ptr<TH1D> hNonSSYield(hNonSSRaw);
  std::unique_ptr<TH1D> hRawForSSPtr(hRawForSS);
  std::unique_ptr<TH1D> hFitStatusPtr(hFitStatus);

  if (!hTrueYield || hTrueYield->Integral() <= 0.0) {
    std::cerr << "ERROR: empty true-phi-phi yield histogram." << std::endl;
    f->Close();
    return;
  }
  if (!hNonSSYield || hNonSSYield->Integral() <= 0.0) {
    std::cerr << "WARNING: nonSS yield histogram is empty or zero. The second background will be poorly constrained." << std::endl;
  }

  SetNiceHist(hTrueYield.get());
  SetNiceHist(hNonSSYield.get());

  // Fit the true-phi-phi continuum with the selected analytic model, excluding the signal window.
  // The default rejected interval is 2.650--2.750 GeV/c^2.
  const TString trueContinuumName = Form(
      "fTruePhiPhiContinuum_%s_rejectSignal", templateBkgTag.Data());
  const TString nonSSContinuumName = Form(
      "fNonSSBkg_%s_%s", templateBkgTag.Data(),
      templateRejectMode == 1 ? "rejectSignal" : "fullRange");
  std::unique_ptr<TF1> fTrueSmooth(FitSmoothTruePhiPhiContinuum(hTrueYield.get(),
                                                                trueContinuumName.Data(),
                                                                mPairMin,
                                                                mPairMax,
                                                                xRejectCenter,
                                                                xHalfWindow,
                                                                templateBkgModel));
  if (!fTrueSmooth) {
    std::cerr << "ERROR: true-phi-phi continuum smoothing failed." << std::endl;
    f->Close();
    return;
  }

  // templateRejectMode=0: reject the signal window only from the SS fit.
  // templateRejectMode=1: reject it from both SS and non-SS fits.
  const bool rejectNonSS = (templateRejectMode == 1);
  std::unique_ptr<TF1> fOtherSmooth(FitSmoothTruePhiPhiContinuum(
      hNonSSYield.get(),
      nonSSContinuumName.Data(),
      mPairMin, mPairMax,
      rejectNonSS ? xRejectCenter : -1.0,
      rejectNonSS ? xHalfWindow : 0.0,
      templateBkgModel));
  if (!fOtherSmooth) {
    std::cerr << "ERROR: nonSS background smoothing failed." << std::endl;
    f->Close();
    return;
  }

  const TString trueTemplateName = Form(
      "hTemplate_truePhiPhi_%s_unitArea", templateBkgTag.Data());
  const TString nonSSTemplateName = Form(
      "hTemplate_nonSSBkg_%s_unitArea", templateBkgTag.Data());
  std::unique_ptr<TH1D> hTrueTpl(BuildSmoothTemplateHistFromTF1(
      hData.get(), fTrueSmooth.get(), trueTemplateName.Data(),
      "true #phi#phi continuum"));
  std::unique_ptr<TH1D> hOtherTpl(BuildSmoothTemplateHistFromTF1(
      hData.get(), fOtherSmooth.get(), nonSSTemplateName.Data(),
      "non-SS background"));
  if (!hTrueTpl || !hOtherTpl) {
    std::cerr << "ERROR: failed to build unit-area analytic templates." << std::endl;
    f->Close();
    return;
  }

  // Compact two-panel template-construction plot requested for the analysis note/talk.
  // It is always saved because it documents exactly how the two fixed background templates are obtained.
  DrawTemplateConstructionTwoPanel(hTrueYield.get(),
                                   fTrueSmooth.get(),
                                   hNonSSYield.get(),
                                   fOtherSmooth.get(),
                                   rejectLo,
                                   rejectHi,
                                   rejectNonSS,
                                   outPrefix);

  if (makeDebugPlots) {
    DrawFig1TrueContinuum(hTrueYield.get(), fTrueSmooth.get(), rejectLo, rejectHi, outPrefix);
    DrawFig2BkgAndControl(hNonSSYield.get(), fOtherSmooth.get(), hOtherRaw.get(), hOtherTpl.get(), outPrefix, selectionMode, dmCtrlMin, dmCtrlMax);
  }

  // Nominal template normalizations before the final sideband rescaling.
  // The continuum functions were fitted to counts per coarse Mpair bin.
  // Therefore integral(function dM)/coarseBinWidth gives the corresponding
  // full-range yield.  For the true-phi-phi template this excludes the exotic
  // contribution because the signal window was rejected in the smooth fit.
  const double nominalTrueYield =
    fTrueSmooth
      ? fTrueSmooth->Integral(mPairMin, mPairMax) /
        actualTemplatePairBinWidth
      : 0.0;

  const double nominalOtherYield =
    fOtherSmooth
      ? fOtherSmooth->Integral(mPairMin, mPairMax) /
        actualTemplatePairBinWidth
      : 0.0;

  // ------------------------------------------------------------
  // 4) Requested final procedure:
  //    First fit the selected Mphiphi data with only the two fixed background
  //    templates away from the signal/exotic window.  Then subtract this total
  //    background and fit the residual pair invariant mass with a Voigtian.
  //    Breit-Wigner shapes are used for each single-phi
  //    signal in the global and per-Mpair 2D mass fits.
  // ------------------------------------------------------------
  const double nTot = hData->Integral(1, hData->GetNbinsX());
  const TwoTemplateSidebandFitResult bkgFit =
    FitTwoTemplateNormsSidebandWLS(
      hData.get(),
      hTrueTpl.get(),
      hOtherTpl.get(),
      rejectLo,
      rejectHi,
      finalTemplateNormMode,
      nominalTrueYield,
      nominalOtherYield);

  const double nTrueVal  = bkgFit.nTrue;
  const double nTrueErr  = bkgFit.nTrueErr;
  const double nOtherVal = bkgFit.nOther;
  const double nOtherErr = bkgFit.nOtherErr;

  std::cout << "\n========== Background-only sideband template fit ==========" << std::endl;
  std::cout << "Fit region           = full M range excluding " << rejectLo << " < M < " << rejectHi << std::endl;
  std::cout << "Data yield           = " << nTot << std::endl;
  std::cout << "Normalization mode   = " << finalTemplateNormMode
            << (finalTemplateNormMode == 0
                  ? " (free scale factors)"
                  : " (scale factors fixed to 1)")
            << std::endl;
  std::cout << "Nominal true yield   = " << bkgFit.nominalTrue
            << ", scale = " << bkgFit.scaleTrue
            << " +/- " << bkgFit.scaleTrueErr << std::endl;
  std::cout << "Nominal nonSS yield  = " << bkgFit.nominalOther
            << ", scale = " << bkgFit.scaleOther
            << " +/- " << bkgFit.scaleOtherErr << std::endl;
  std::cout << "N_truePhiPhi bkg     = " << nTrueVal << " +/- " << nTrueErr << std::endl;
  std::cout << "N_nonSS bkg          = " << nOtherVal << " +/- " << nOtherErr << std::endl;
  std::cout << "chi2/ndf             = " << bkgFit.chi2 << "/" << bkgFit.ndf << std::endl;
  if (!bkgFit.ok) {
    std::cout << "WARNING: sideband-only background-template fit may be underconstrained. "
              << "Check the output plot and try a wider sideband/rejection definition." << std::endl;
  }
  std::cout << "==========================================================\n" << std::endl;

  // ------------------------------------------------------------
  // 5) Build background components, subtract them, and fit the final
  //    pair-invariant-mass residual with a Voigtian.
  // ------------------------------------------------------------
  std::unique_ptr<TH1D> hTrueScaled(MakeExpectedFromTemplates(hData.get(), hTrueTpl.get(), hOtherTpl.get(),
                                                              nTrueVal, nOtherVal, 0.0,
                                                              xMassInit, bwGammaInit,
                                                              xResolution,
                                                              mPairMin, mPairMax,
                                                              "hFitComponent_truePhiPhi_bkgOnly", 1));
  std::unique_ptr<TH1D> hOtherScaled(MakeExpectedFromTemplates(hData.get(), hTrueTpl.get(), hOtherTpl.get(),
                                                               nTrueVal, nOtherVal, 0.0,
                                                               xMassInit, bwGammaInit,
                                                               xResolution,
                                                               mPairMin, mPairMax,
                                                               "hFitComponent_nonSS_bkgOnly", 2));
  std::unique_ptr<TH1D> hTotalBkg(MakeExpectedFromTemplates(hData.get(), hTrueTpl.get(), hOtherTpl.get(),
                                                            nTrueVal, nOtherVal, 0.0,
                                                            xMassInit, bwGammaInit,
                                                            xResolution,
                                                            mPairMin, mPairMax,
                                                            "hFitComponent_totalBkg_sidebandFit", 0));
  std::unique_ptr<TH1D> hResidual(MakeResidualHist(hData.get(), hTotalBkg.get(),
                                                   "hResidual_data_minus_sidebandBkg", false));

  const double residualFitLo = mPairMin;
  const double residualFitHi = mPairMax;
  std::unique_ptr<TF1> fResidualVoigt(FitResidualVoigtian(
      hResidual.get(), residualFitLo, residualFitHi,
      mPairMin, mPairMax, xMassInit, xMassMin, xMassMax,
      bwGammaInit, gammaMin, gammaMax, xResolution));
  std::unique_ptr<TH1D> hResidualVoigt(MakeFunctionHistOnDataBins(
      hData.get(), fResidualVoigt.get(),
      "hResidualVoigt_fitFunction_onBins"));

  if (fResidualVoigt) {
    const double nSigVal   = fResidualVoigt->GetParameter(0);
    const double nSigErr   = fResidualVoigt->GetParError(0);
    const double gammaXVal = fResidualVoigt->GetParameter(1);
    const double gammaXErr = fResidualVoigt->GetParError(1);
    const double meanXVal  = fResidualVoigt->GetParameter(2);
    const double meanXErr  = fResidualVoigt->GetParError(2);
    const double yieldOverError = (nSigErr > 0.0 ? nSigVal / nSigErr : 0.0);

    std::cout << "\n========== Residual Voigtian fit ===============" << std::endl;
    std::cout << "Residual fit range   = " << residualFitLo << " < M < " << residualFitHi << std::endl;
    std::cout << "N_X                  = " << nSigVal << " +/- " << nSigErr << std::endl;
    std::cout << "M_X                  = " << meanXVal << " +/- " << meanXErr << std::endl;
    std::cout << "Gamma_X              = " << gammaXVal << " +/- " << gammaXErr << std::endl;
    std::cout << "sigma_res (fixed)     = " << xResolution << " GeV/c^2 = "
              << 1000.0*xResolution << " MeV/c^2" << std::endl;
    std::cout << "Yield/error (not Z)  = " << yieldOverError << std::endl;
    std::cout << "==============================================\n" << std::endl;
  } else {
    std::cout << "WARNING: residual Voigtian fit was not produced." << std::endl;
  }

  DrawFig3BkgSubtractedFit(hData.get(),
                           hTrueScaled.get(),
                           hOtherScaled.get(),
                           hTotalBkg.get(),
                           hResidual.get(),
                           fResidualVoigt.get(),
                           bkgFit,
                           outPrefix,
                           rejectLo,
                           rejectHi);

  SaveBkgSubtractedFitCSV(hData.get(), hTrueScaled.get(), hOtherScaled.get(),
                          hTotalBkg.get(), hResidual.get(), fResidualVoigt.get(), outPrefix);

  // R_data is always evaluated in exactly [xMassMin,xMassMax].
  const SHMDataObservableResult shmDataObs = ComputeSHMDataObservable(hData.get(),
                                                                      hTrueScaled.get(),
                                                                      hOtherScaled.get(),
                                                                      hTotalBkg.get(),
                                                                      hResidual.get(),
                                                                      hTrueTpl.get(),
                                                                      fResidualVoigt.get(),
                                                                      bkgFit,
                                                                      outPrefix,
                                                                      xMassMin,
                                                                      xMassMax,
                                                                      mPairMin,
                                                                      mPairMax,
                                                                      ptMin,
                                                                      ptMax,
                                                                      rapidityMin,
                                                                      rapidityMax,
                                                                      selectionMode,
                                                                      selectionCut,
                                                                      rSHMFiducial);

  // Final template-fit canvas is drawn left/right: background fit and residual Voigtian.

  TFile fout(Form("%s_output.root", outPrefix), "RECREATE");
  h3->Write();
  h2Fit->Write();
  hData->Write();
  if (hOtherRaw) hOtherRaw->Write();
  if (hDirectSignal) hDirectSignal->Write();
  if (hDirectOutside) hDirectOutside->Write();
  if (h2DirectPhiPairWindow) h2DirectPhiPairWindow->Write();
  if (phiPurityVsPt.purityGraph && phiPurityVsPt.successfulFits > 0)
    phiPurityVsPt.purityGraph->Write();
  if (deltaMSignificanceScan.significanceGraph)
    deltaMSignificanceScan.significanceGraph->Write();
  if (deltaMSignificanceScan.purityGraph)
    deltaMSignificanceScan.purityGraph->Write();
  if (deltaMSignificanceScan.purityWeightedSignificanceGraph)
    deltaMSignificanceScan.purityWeightedSignificanceGraph->Write();
  if (deltaMSignificanceScan.signalGraph)
    deltaMSignificanceScan.signalGraph->Write();
  if (deltaMSignificanceScan.backgroundGraph)
    deltaMSignificanceScan.backgroundGraph->Write();
  hTrueYield->Write();
  if (hRawForSSPtr) hRawForSSPtr->Write();
  if (hFitStatusPtr) hFitStatusPtr->Write();
  hTrueTpl->Write();
  hOtherTpl->Write();
  hTrueScaled->Write();
  hOtherScaled->Write();
  hTotalBkg->Write();
  hResidual->Write();
  if (hResidualVoigt) hResidualVoigt->Write();
  if (fResidualVoigt)
    fResidualVoigt->Write("fResidualVoigt_fit_dataMinusBkg");
  if (localP0Scan.p0Graph) localP0Scan.p0Graph->Write("gLocalP0_asymptotic");
  if (localP0Scan.zGraph) localP0Scan.zGraph->Write("gLocalZ_asymptotic");
  if (localP0Scan.p0Band) {
    localP0Scan.p0Band->Write("gLocalP0_parameterVariationBand");
  }
  if (localP0Scan.zBand) {
    localP0Scan.zBand->Write("gLocalZ_parameterVariationBand");
  }
  TNamed localP0Def(
    "local_p0_definition",
    Form("Direct DeltaM<%.9g Voigtian+%s spectrum; one-sided RooStats AsymptoticCalculator. "
         "Direct-fit point: p0=%.9g, Z=%.9g, M=%.9g, Gamma=%.9g. "
         "The Gaussian resolution sigma is fixed to %.9g GeV/c2. "
         "Scan: best central p0=%.9g, best central Z=%.9g at M=%.9g; "
         "Gamma=%.9g +/- %.9g, fitted M=%.9g +/- %.9g, bandMode=%d. "
         "The band is a +/-1sigma fit-parameter-variation envelope, not a confidence interval on p0.",
         normalFitDeltaM, BackgroundLabel(directBkgModel),
         directFit.localP0, directFit.localZ,
         directFit.mean, directFit.gamma, xResolution,
         localP0Scan.bestP0, localP0Scan.bestZ,
         localP0Scan.bestMass, localP0Scan.gamma, localP0Scan.gammaErr,
         localP0Scan.fittedMean, localP0Scan.fittedMeanErr,
         localP0Scan.bandMode));
  localP0Def.Write();
  if (deltaMSignificanceScan.significanceGraph) {
    TNamed deltaMScanDef(
      "deltaM_pair_significance_definition",
      Form("From the full 2D mass histogram in Mphiphi=[%.9g,%.9g), "
           "pT=[%.9g,%.9g), rapidity(axis3)=[%.9g,%.9g). "
           "Each bin is decomposed with the fitted posterior P_SS; "
           "S=sum(Nbin*P_SS), B=sum(Nbin*(1-P_SS)), "
           "significance=S/sqrt(S+B), purity=S/B, "
           "purity-weighted significance proxy="
           "[S/sqrt(S+B)]*[S/(S+B)]. "
           "Best significance cut=%.9g with Z=%.9g; "
           "best proxy cut=%.9g with proxy=%.9g.",
           mPairMin,mPairMax,ptMin,ptMax,rapidityMin,rapidityMax,
           deltaMSignificanceScan.bestCut,
           deltaMSignificanceScan.bestSignificance,
           deltaMSignificanceScan.bestPurityWeightedCut,
           deltaMSignificanceScan.bestPurityWeightedSignificance));
    deltaMScanDef.Write();
  }
  TNamed shmObsDef(
    "SHM_data_observable",
    Form(
      "M=[%.6g,%.6g], NX_window=%.9g +/- %.9g, "
      "NtruePhiPhi_window=%.9g +/- %.9g, R_data=%.9g +/- %.9g, "
      "GammaX=%.9g +/- %.9g, XresolutionSigma=%.9g, phiBWGammaGlobal=%.9g, "
      "templateScales=(%.9g,%.9g), R_SHM_input=%.9g, "
      "data_over_SHM=%.9g +/- %.9g, rapidity(axis3)=[%.9g,%.9g)",
      shmDataObs.obsMassMin, shmDataObs.obsMassMax,
      shmDataObs.nXWindow, shmDataObs.nXWindowErr,
      shmDataObs.nTruePhiPhiWindow, shmDataObs.nTruePhiPhiWindowErr,
      shmDataObs.rData, shmDataObs.rDataErr,
      fResidualVoigt ? fResidualVoigt->GetParameter(1) : 0.0,
      fResidualVoigt ? fResidualVoigt->GetParError(1) : 0.0,
      xResolution, globalFit.gamma, bkgFit.scaleTrue, bkgFit.scaleOther,
      shmDataObs.rSHM, shmDataObs.dataOverSHM, shmDataObs.dataOverSHMErr,
      rapidityMin,rapidityMax));
  shmObsDef.Write();
  if (phiPurityVsPt.purityGraph && phiPurityVsPt.successfulFits > 0) {
    TNamed phiPurityDefinition(
      "phi_purity_vs_pt_definition",
      Form("Input=%s; pT bins=[0.5,0.8,1.2,1.6,2.0,2.5,3.0,4.0,"
           "5.0,6.0,8.0,10.0,20.0,50.0] GeV/c; model=Breit-Wigner+"
           "%s; fit/display range=[%.9g,%.9g] GeV/c2; "
           "purity=S/(S+B) in fitted mean +/- %.9g GeV/c2; "
           "successful fits=%d/13.",
           phiMassVsPtName,BackgroundLabel(phiBkgModel),
           phiFitMin,phiFitMax,phiPurityHalfWidth,
           phiPurityVsPt.successfulFits));
    phiPurityDefinition.Write();
  }
  fTrueSmooth->Write(trueContinuumName.Data());
  fOtherSmooth->Write(nonSSContinuumName.Data());
  if (hNonSSYield) hNonSSYield->Write();
  TNamed modelDef(
    "model_definition",
    Form(
      "The selected Mphiphi data are described by unit-area SS and nonSS templates "
      "built with %s. Pair rapidity axis 3 is restricted to [%.9g,%.9g). "
      "templateRejectMode=%d (0 rejects the signal window only from SS; 1 rejects it from both). "
      "finalTemplateNormMode=%d (0 floats template scales; 1 fixes them to one). "
      "Each single-phi mass in the 2D decomposition uses a Breit-Wigner with global fitted mean %.9g and width %.9g. "
      "pairMassFitMode=%d (0 global fractions, 1 local yields with fixed global shapes, 2 local mean/width/background shapes free). "
      "The residual pair mass and the independent DeltaM<%.9g direct fit both use binWidth*N_X*TMath::Voigt with Minuit and fixed Gaussian sigma %.9g GeV/c2. "
      "The direct spectrum also has a one-sided RooStats AsymptoticCalculator local-p0 scan. "
      "R_data is computed in [%.9g,%.9g].",
      ContinuumModelLabel(),rapidityMin,rapidityMax,
      templateRejectMode, finalTemplateNormMode, globalFit.mean, globalFit.gamma,
      pairMassFitMode, normalFitDeltaM, xResolution, xMassMin, xMassMax));
  modelDef.Write();
  fout.Close();

  std::cout << "Saved outputs with prefix: " << outPrefix << std::endl;
  std::cout << "  " << outPrefix << "_Fig1_templateConstruction_twoPanel.png/pdf\n";
  if (makeDebugPlots) {
    std::cout << "  " << outPrefix << "_Fig1_trueContinuum.png/pdf\n";
    std::cout << "  " << outPrefix << "_Fig2_nonSS_control.png/pdf\n";
  }
  std::cout << "  " << outPrefix << "_Fig3_finalTemplateFit.png/pdf\n";
  std::cout << "  " << outPrefix << "_Fig4_directPairMass_threePanel.png/pdf\n";
  std::cout << "  " << outPrefix << "_directPairMassFit.csv\n";
  if (makeLocalP0Scan) {
    std::cout << "  " << outPrefix << "_Fig5_localP0_asymptotic.png/pdf\n";
    std::cout << "  " << outPrefix << "_localP0_asymptotic.csv\n";
  }
  if (makeDeltaMSignificanceScan) {
    std::cout << "  " << outPrefix
              << "_Fig6_deltaM_significancePurity.png/pdf\n";
    std::cout << "  " << outPrefix
              << "_Fig7_deltaM_significanceTimesPurity.png/pdf\n";
    std::cout << "  " << outPrefix
              << "_deltaM_pairSignificance.csv\n";
  }
  if (phiPurityVsPt.successfulFits > 0) {
    std::cout << "  " << outPrefix << "_Fig8_phiPurity_vsPt.png/pdf\n";
    std::cout << "  " << outPrefix << "_phiPurity_vsPt_fits.png/pdf\n";
    std::cout << "  " << outPrefix << "_phiPurity_vsPt.csv\n";
  }
  std::cout << "  " << outPrefix << "_fit_components.csv\n";
  std::cout << "  " << outPrefix << "_SHM_data_observable.csv\n";
  std::cout << "  " << outPrefix << "_output.root\n";
  if (save2DFits) {
    std::cout << "  all performed 2D mass-fit plots in folder: " << fitDir << "\n";
    std::cout << "  " << fitDir << "/all_2D_invariant_mass_inputs.root\n";
    std::cout << "  " << fitDir << "/all_2D_invariant_mass_fit_results.csv\n";
  }

  f->Close();
}

} // namespace

// Compact public entry point.
//
// phiBkgModel, templateBkgModel, and directBkgModel use the same numbering:
//   0 pol2, 1 pol3, 2 exp(pol2), 3 exp(pol3),
//   4 Chebyshev-2, 5 Chebyshev-3, 6 Bernstein-2, 7 Bernstein-3.
// runTemplateFit:
//   0 = direct M(phi phi) fit only, 1 = also run the template analysis.
void DoublePhi_analysis_Voigtian(
    const char* inputFile = "../AnalysisResults_pid2003.root",
    const char* sparseName = "doublephimeson/SEMassDoublePhi",
    double ptMin = 6.0,
    double ptMax = 100.0,
    double deltaMMax = 0.004,
    int phiBkgModel = 0,
    int templateBkgModel = 5,
    int directBkgModel = 5,
    int runTemplateFit = 1,
    double phiFitMin = 1.005,
    double phiFitMax = 1.035,
    double xResolution = 0.015)
{
  RunDoublePhiAnalysis(
      inputFile,
      sparseName,
      ptMin,
      ptMax,
      1.005,   // phiMassMin
      1.035,   // phiMassMax
      3,       // Chebyshev order in the 2D single-phi background
      2.5,     // mPairMin
      2.9,     // mPairMax
      0.01,   // templatePairBinWidth
      0.01,   // finalPairBinWidth
      0,       // stored-DeltaM selection
      deltaMMax,
      0.005,   // localScanStep
      deltaMMax,
      templateBkgModel,
      2.65,    // template rejection window
      2.73,
      2.63,    // X-yield/mean window
      2.75,
      0.001,   // Gamma_X limits
      0.10,
      1,       // save all 2D fits
      0.010,   // DeltaM control region
      0.20,
      0,       // template uses the same selected region as data
      0,       // extra debug plots off
      -1.0,    // optional R_SHM input disabled
      1,       // fixed nominal template normalizations
      1,       // local 2D yields, global shapes fixed
      0.0042,  // single-phi BW width initialization and limits
      0.0010,
      0.0200,
      1,       // reject signal window from both templates
      1,       // local-p0 scan on
      2.55,
      2.85,
      -1.0,    // scan uses fitted Gamma_X
      1,       // Gamma_X variation band
      -1.0,    // scan uses fitted Gamma_X uncertainty
      runTemplateFit,
      directBkgModel,
      0.0,     // pair rapidity range
      0.8,
      1,       // DeltaM significance/purity scan on
      0.001,
      0.030,
      0.001,
      1,       // phi purity versus pT on
      "doublephimeson/hPhiMassVsPt",
      phiBkgModel,
      phiFitMin,
      phiFitMax,
      0.007,
      xResolution);
}
