

#include <TAxis.h>
#include <TCanvas.h>
#include <TClass.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TVirtualFitter.h>

#include <Math/MinimizerOptions.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace PhiPtCalibration
{

  constexpr double kMassK = 0.493677;      // GeV/c^2
  constexpr double kMassPhiPDG = 1.019461; // GeV/c^2
  // TMath::BreitWigner(..., Gamma) and TMath::Voigt(..., sigma, Gamma, 4)
  // use the same full Lorentzian width Gamma.  Keep this definition identical in
  // data fits, single-phi closure toys and the X-resolution propagation.
  constexpr double kWidthPhi = 0.0042500; // Lorentzian FWHM, GeV/c^2

  // The projected mass histograms are fitted without rebinning.  This value is
  // updated from the native histogram axis before every fit.  Multiplication by
  // the bin width makes par[0] the fitted signal yield, as in the reference code.
  double gPhiMassBinWidth = 0.001;

  Double_t PhiVoigtPol2(Double_t *x, Double_t *par)
  {
    const double signal = gPhiMassBinWidth * par[0] * TMath::Voigt(x[0] - par[2], par[3], par[1], 4);
    const double background = par[4] + par[5] * x[0] + par[6] * x[0] * x[0];
    return signal + background;
  }

  Double_t PhiVoigtOnly(Double_t *x, Double_t *par)
  {
    return gPhiMassBinWidth * par[0] * TMath::Voigt(x[0] - par[2], par[3], par[1], 4);
  }

  Double_t PhiPol2(Double_t *x, Double_t *par)
  {
    return par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
  }

  struct FitPoint
  {
    bool valid = false;
    int status = -999;
    int covQual = -1;
    double ptLow = 0.;
    double ptHigh = 0.;
    double mean = 0.;
    double meanError = 0.;
    double sigma = 0.;
    double sigmaError = 0.;
    double nSignal = 0.;
    double nSignalError = 0.;
    double nBackground = 0.;
    double chi2 = 0.;
    int ndf = 0;
  };

  struct ToyPoint
  {
    bool valid = false;
    long long accepted = 0;
    double slope = 0.; // GeV per unit epsilon
    double slopeError = 0.;
    double targetShift = 0.; // m_phi(2026)-m_phi(2025), GeV
    double targetShiftError = 0.;
    double epsilonCorrection = 0.; // applied directly to 2025 kaons
    double epsilonError = 0.;
    double reproducedShift = 0.;
    double reproducedShiftError = 0.;
  };

  std::string SafeName(const std::string &input)
  {
    std::string output = input;
    for (char &c : output)
    {
      if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_'))
        c = '_';
    }
    return output;
  }

  std::vector<double> ParseEdges(const char *text)
  {
    std::vector<double> edges;
    if (!text)
      return edges;
    std::stringstream stream(text);
    std::string token;
    while (std::getline(stream, token, ','))
    {
      std::stringstream valueStream(token);
      double value = 0.;
      if (valueStream >> value)
        edges.push_back(value);
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
  }

  TObject *FindObjectRecursive(TDirectory *directory, const std::string &name)
  {
    if (!directory)
      return nullptr;
    if (TObject *direct = directory->Get(name.c_str()))
      return direct;

    TIter next(directory->GetListOfKeys());
    while (TKey *key = dynamic_cast<TKey *>(next()))
    {
      TClass *objectClass = gROOT->GetClass(key->GetClassName());
      if (!objectClass)
        continue;
      if (name == key->GetName())
        return key->ReadObj();
      if (!objectClass->InheritsFrom(TDirectory::Class()))
        continue;
      TObject *object = key->ReadObj();
      TDirectory *subdirectory = dynamic_cast<TDirectory *>(object);
      if (!subdirectory)
      {
        delete object;
        continue;
      }
      if (TObject *found = FindObjectRecursive(subdirectory, name))
        return found;
    }
    return nullptr;
  }

  TH2 *GetMassVsPt(TFile *file, const std::string &pathOrName)
  {
    if (!file || file->IsZombie())
      return nullptr;
    TObject *object = file->Get(pathOrName.c_str());
    if (!object)
    {
      const std::size_t slash = pathOrName.find_last_of('/');
      const std::string shortName = slash == std::string::npos
                                        ? pathOrName
                                        : pathOrName.substr(slash + 1);
      object = FindObjectRecursive(file, shortName);
    }
    TH2 *histogram = dynamic_cast<TH2 *>(object);
    if (!histogram)
    {
      std::cerr << "ERROR: TH2 '" << pathOrName << "' was not found in "
                << file->GetName() << std::endl;
      return nullptr;
    }
    if (!(histogram->GetXaxis()->GetXmin() < kMassPhiPDG &&
          histogram->GetXaxis()->GetXmax() > kMassPhiPDG))
    {
      std::cerr << "ERROR: x axis of " << histogram->GetName()
                << " does not contain the phi mass. Expected x=mass, y=pT."
                << std::endl;
      return nullptr;
    }
    return histogram;
  }

  TH1D *ProjectMass(TH2 *histogram, double ptLow, double ptHigh, const std::string &name)
  {
    if (!histogram || ptHigh <= ptLow)
      return nullptr;
    TAxis *ptAxis = histogram->GetYaxis();
    int first = ptAxis->FindBin(ptLow + 1.e-9);
    int last = ptAxis->FindBin(ptHigh - 1.e-9);
    first = std::max(1, std::min(ptAxis->GetNbins(), first));
    last = std::max(first, std::min(ptAxis->GetNbins(), last));
    TH1D *projection = histogram->ProjectionX(name.c_str(), first, last, "e");
    if (projection)
    {
      projection->SetDirectory(nullptr);
      projection->SetTitle(Form("%.3g < p_{T,#phi} < %.3g GeV/c", ptLow, ptHigh));
      projection->GetXaxis()->SetTitle("m_{K^{+}K^{-}} (GeV/c^{2})");
      projection->GetYaxis()->SetTitle("Candidates");
    }
    return projection;
  }

  double IntegralInRange(TH1D *histogram, double low, double high)
  {
    if (!histogram)
      return 0.;
    const int first = std::max(1, histogram->GetXaxis()->FindBin(low + 1.e-9));
    const int last = std::min(histogram->GetNbinsX(), histogram->GetXaxis()->FindBin(high - 1.e-9));
    return histogram->Integral(first, last);
  }

  std::vector<double> HistogramEdgesInRange(TH1D *histogram, double low, double high)
  {
    std::vector<double> edges;
    if (!histogram || high <= low)
      return edges;
    edges.push_back(low);
    for (int b = 1; b <= histogram->GetNbinsX(); ++b)
    {
      const double edge = histogram->GetXaxis()->GetBinUpEdge(b);
      if (edge > low + 1.e-10 && edge < high - 1.e-10)
        edges.push_back(edge);
    }
    edges.push_back(high);
    return edges;
  }

  FitPoint FitOnePtBin(TH1D *histogram, const std::string &label,
                       int binIndex, double ptLow, double ptHigh,
                       double fitLow, double fitHigh,
                       TDirectory *outputDirectory, const std::string &pdfName)
  {
    FitPoint result;
    result.ptLow = ptLow;
    result.ptHigh = ptHigh;
    if (!histogram || IntegralInRange(histogram, fitLow, fitHigh) < 500.)
    {
      std::cerr << "WARNING: insufficient entries for " << label << " pT bin "
                << ptLow << "-" << ptHigh << std::endl;
      return result;
    }

    const std::string suffix = SafeName(label) + Form("_pt%02d", binIndex);
    const double total = IntegralInRange(histogram, fitLow, fitHigh);
    const int firstBin = histogram->GetXaxis()->FindBin(fitLow + 1.e-9);
    const int lastBin = histogram->GetXaxis()->FindBin(fitHigh - 1.e-9);
    const double binWidth = histogram->GetXaxis()->GetBinWidth(
        histogram->GetXaxis()->FindBin(kMassPhiPDG));
    gPhiMassBinWidth = binWidth;
    const int nEdgeBins = std::max(1, (lastBin - firstBin + 1) / 10);
    double leftBackground = 0.;
    double rightBackground = 0.;
    for (int i = 0; i < nEdgeBins; ++i)
    {
      leftBackground += histogram->GetBinContent(firstBin + i);
      rightBackground += histogram->GetBinContent(lastBin - i);
    }
    leftBackground /= nEdgeBins;
    rightBackground /= nEdgeBins;
    const double backgroundAtCenter = 0.5 * (leftBackground + rightBackground);
    const double backgroundSlope = (rightBackground - leftBackground) /
                                   std::max(1.e-6, fitHigh - fitLow);
    const double seedSigma = 0.0018;
    const int signalFirstBin = histogram->GetXaxis()->FindBin(1.010 + 1.e-9);
    const int signalLastBin = histogram->GetXaxis()->FindBin(1.030 - 1.e-9);
    const int numberOfSignalBins = signalLastBin - signalFirstBin + 1;
    const double yieldSeed = std::max(
        100., histogram->Integral(signalFirstBin, signalLastBin) -
                  backgroundAtCenter * numberOfSignalBins);

    // Same plain-ROOT parameterization as the supplied working example:
    // binWidth * Yield * TMath::Voigt(x-Mass, Resolution, Width, 4) + pol2.
    TF1 model(("fVoigtPol2_" + suffix).c_str(), PhiVoigtPol2,
              fitLow, fitHigh, 7);
    model.SetParNames("Yield", "Gamma", "mPhi", "Sigma",
                      "Bkg0", "Bkg1", "Bkg2");
    model.SetParameters(yieldSeed, kWidthPhi, kMassPhiPDG, seedSigma,
                        backgroundAtCenter - backgroundSlope * kMassPhiPDG,
                        backgroundSlope, 0.);
    model.SetParLimits(0, 0., 2. * total);
    model.FixParameter(1, kWidthPhi);
    model.SetParLimits(2, 1.010, 1.030);
    model.SetParLimits(3, 0.0001, 0.0100);
    model.SetNpx(10000);

    TVirtualFitter::SetDefaultFitter("Minuit2");
    TVirtualFitter::SetMaxIterations(10000);
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
    TFitResultPtr fitResult = histogram->Fit(&model, "Q0ERMS");

    result.status = static_cast<int>(fitResult);
    result.covQual = fitResult.Get() ? fitResult->CovMatrixStatus() : -1;
    result.mean = model.GetParameter(2);
    result.meanError = model.GetParError(2);
    result.sigma = model.GetParameter(3);
    result.sigmaError = model.GetParError(3);
    result.chi2 = model.GetChisquare();
    result.ndf = model.GetNDF();

    TF1 signalFunction(("fSignal_" + suffix).c_str(), PhiVoigtOnly,
                       fitLow, fitHigh, 4);
    signalFunction.SetParameters(model.GetParameter(0), kWidthPhi,
                                 model.GetParameter(2), model.GetParameter(3));
    signalFunction.SetNpx(10000);
    TF1 backgroundFunction(("fBackground_" + suffix).c_str(), PhiPol2,
                           fitLow, fitHigh, 3);
    backgroundFunction.SetParameters(model.GetParameter(4), model.GetParameter(5),
                                     model.GetParameter(6));
    backgroundFunction.SetNpx(10000);
    result.nSignal = model.GetParameter(0);
    result.nSignalError = model.GetParError(0);
    result.nBackground = backgroundFunction.Integral(fitLow, fitHigh) / binWidth;
    result.valid = fitResult.Get() &&
                   std::isfinite(result.mean) && result.meanError > 0. &&
                   std::isfinite(result.sigma) && result.sigmaError > 0. &&
                   result.sigma > 0.0001 && result.sigma < 0.0100;

    TH1D *hModel = dynamic_cast<TH1D *>(histogram->Clone(("hModel_" + suffix).c_str()));
    TH1D *hSignal = dynamic_cast<TH1D *>(histogram->Clone(("hSignal_" + suffix).c_str()));
    TH1D *hBackground = dynamic_cast<TH1D *>(histogram->Clone(("hBackground_" + suffix).c_str()));
    for (TH1D *h : {hModel, hSignal, hBackground})
    {
      h->Reset("ICES");
      h->SetDirectory(nullptr);
    }

    for (int b = 1; b <= histogram->GetNbinsX(); ++b)
    {
      const double x = histogram->GetXaxis()->GetBinCenter(b);
      if (x < fitLow || x > fitHigh)
        continue;
      const double expectedSignal = signalFunction.Eval(x);
      const double expectedBackground = backgroundFunction.Eval(x);
      hSignal->SetBinContent(b, expectedSignal);
      hBackground->SetBinContent(b, expectedBackground);
      hModel->SetBinContent(b, expectedSignal + expectedBackground);
    }

    histogram->SetMarkerStyle(20);
    histogram->SetMarkerSize(0.75);
    histogram->SetLineColor(kBlack);
    hModel->SetLineColor(kRed + 1);
    hModel->SetLineWidth(3);
    hSignal->SetLineColor(kBlue + 1);
    hSignal->SetLineWidth(2);
    hSignal->SetLineStyle(2);
    hBackground->SetLineColor(kGreen + 2);
    hBackground->SetLineWidth(2);
    hBackground->SetLineStyle(7);

    TCanvas canvas(("c_" + suffix).c_str(), "phi fit", 900, 750);
    canvas.SetLeftMargin(0.13);
    canvas.SetBottomMargin(0.12);
    histogram->GetXaxis()->SetRangeUser(fitLow, fitHigh);
    histogram->SetMaximum(1.16 * histogram->GetMaximum());
    histogram->Draw("E");
    hModel->Draw("HIST SAME");
    hSignal->Draw("HIST SAME");
    hBackground->Draw("HIST SAME");
    histogram->Draw("E SAME");

    TLegend legend(0.58, 0.66, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(histogram, (label + " data").c_str(), "lep");
    legend.AddEntry(hModel, "Voigtian + pol2", "l");
    legend.AddEntry(hSignal, "Voigtian signal", "l");
    legend.AddEntry(hBackground, "pol2 background", "l");
    legend.Draw();

    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.034);
    text.DrawLatex(0.16, 0.88, Form("%s: %.3g < p_{T,#phi} < %.3g GeV/c", label.c_str(), ptLow, ptHigh));
    text.DrawLatex(0.16, 0.83, Form("m_{#phi}=%.3f #pm %.3f MeV/c^{2}", 1000. * result.mean, 1000. * result.meanError));
    text.DrawLatex(0.16, 0.78, Form("#sigma_{#phi}=%.3f #pm %.3f MeV/c^{2}", 1000. * result.sigma, 1000. * result.sigmaError));
    text.DrawLatex(0.16, 0.73,
                   Form("#Gamma_{#phi}=%.2f MeV/c^{2} (fixed)",
                        1000. * kWidthPhi));
    text.DrawLatex(0.16, 0.68, Form("status=%d, covQual=%d, #chi^{2}/ndf=%.1f/%d", result.status, result.covQual, result.chi2, result.ndf));
    canvas.Print(pdfName.c_str());

    if (outputDirectory)
    {
      outputDirectory->cd();
      histogram->Write();
      hModel->Write();
      hSignal->Write();
      hBackground->Write();
      model.Write();
      signalFunction.Write();
      backgroundFunction.Write();
      if (fitResult.Get())
        fitResult->Write(("fitResult_" + suffix).c_str());
    }
    delete hModel;
    delete hSignal;
    delete hBackground;
    return result;
  }

  TH1D *BuildSignalWeightedPt(TH2 *h2025, TH2 *h2026)
  {
    if (!h2025 || !h2026)
      return nullptr;
    TH1D *combined = nullptr;
    int sourceIndex = 0;
    for (TH2 *source : {h2025, h2026})
    {
      const int signalFirst = source->GetXaxis()->FindBin(1.010 + 1.e-9);
      const int signalLast = source->GetXaxis()->FindBin(1.030 - 1.e-9);
      const int leftFirst = source->GetXaxis()->FindBin(1.000 + 1.e-9);
      const int leftLast = source->GetXaxis()->FindBin(1.010 - 1.e-9);
      const int rightFirst = source->GetXaxis()->FindBin(1.030 + 1.e-9);
      const int rightLast = source->GetXaxis()->FindBin(1.040 - 1.e-9);
      std::unique_ptr<TH1D> signal(source->ProjectionY(
          Form("tmpPtSignal_%d", sourceIndex), signalFirst, signalLast, "e"));
      std::unique_ptr<TH1D> left(source->ProjectionY(
          Form("tmpPtLeft_%d", sourceIndex), leftFirst, leftLast, "e"));
      std::unique_ptr<TH1D> right(source->ProjectionY(
          Form("tmpPtRight_%d", sourceIndex), rightFirst, rightLast, "e"));
      signal->SetDirectory(nullptr);
      left->SetDirectory(nullptr);
      right->SetDirectory(nullptr);
      signal->Add(left.get(), -1.0);
      signal->Add(right.get(), -1.0);
      for (int b = 1; b <= signal->GetNbinsX(); ++b)
      {
        if (signal->GetBinContent(b) < 0.)
          signal->SetBinContent(b, 0.);
      }
      if (!combined)
      {
        combined = dynamic_cast<TH1D *>(signal->Clone("hPhiPtForToy"));
        combined->SetDirectory(nullptr);
      }
      else
      {
        combined->Add(signal.get());
      }
      ++sourceIndex;
    }
    if (combined)
    {
      combined->SetTitle("Sideband-subtracted #phi p_{T} spectrum used by toy");
      combined->GetXaxis()->SetTitle("p_{T,#phi} (GeV/c)");
      combined->GetYaxis()->SetTitle("Signal-weighted candidates");
    }
    return combined;
  }

  double SamplePt(TH1D *spectrum, double ptLow, double ptHigh, TRandom3 &random)
  {
    if (!spectrum || ptHigh <= ptLow)
      return ptLow;
    std::vector<double> lowEdges;
    std::vector<double> highEdges;
    std::vector<double> cumulative;
    double sum = 0.;
    for (int b = 1; b <= spectrum->GetNbinsX(); ++b)
    {
      const double low = std::max(ptLow, spectrum->GetXaxis()->GetBinLowEdge(b));
      const double high = std::min(ptHigh, spectrum->GetXaxis()->GetBinUpEdge(b));
      if (high <= low)
        continue;
      const double weight = std::max(0., spectrum->GetBinContent(b)) *
                            (high - low) / spectrum->GetXaxis()->GetBinWidth(b);
      if (weight <= 0.)
        continue;
      sum += weight;
      lowEdges.push_back(low);
      highEdges.push_back(high);
      cumulative.push_back(sum);
    }
    if (sum <= 0.)
      return random.Uniform(ptLow, ptHigh);
    const double selected = random.Uniform(0., sum);
    const std::size_t index = std::lower_bound(cumulative.begin(), cumulative.end(),
                                               selected) -
                              cumulative.begin();
    return random.Uniform(lowEdges[index], highEdges[index]);
  }

  bool GeneratePhiDecay(double generatedPhiMass, double phiPt, double yMax,
                        TRandom3 &random,
                        TLorentzVector &kaonPlus, TLorentzVector &kaonMinus)
  {
    const double phiY = random.Uniform(-yMax, yMax);
    const double phiAzimuth = random.Uniform(0., TMath::TwoPi());
    const double transverseMass = std::sqrt(generatedPhiMass * generatedPhiMass +
                                            phiPt * phiPt);
    TLorentzVector phi;
    // Build the exact rapidity four-vector explicitly; this is also stable at
    // pT=0, unlike converting rapidity to pseudorapidity.
    const double pz = transverseMass * std::sinh(phiY);
    const double energy = transverseMass * std::cosh(phiY);
    phi.SetPxPyPzE(phiPt * std::cos(phiAzimuth), phiPt * std::sin(phiAzimuth),
                   pz, energy);

    const double q2 = 0.25 * generatedPhiMass * generatedPhiMass -
                      kMassK * kMassK;
    if (q2 <= 0.)
      return false;
    const double q = std::sqrt(q2);
    const double cosTheta = random.Uniform(-1., 1.);
    const double sinTheta = std::sqrt(std::max(0., 1. - cosTheta * cosTheta));
    const double decayAzimuth = random.Uniform(0., TMath::TwoPi());
    const TVector3 momentum(q * sinTheta * std::cos(decayAzimuth),
                            q * sinTheta * std::sin(decayAzimuth),
                            q * cosTheta);
    const double kaonEnergy = std::sqrt(q * q + kMassK * kMassK);
    kaonPlus.SetPxPyPzE(momentum.X(), momentum.Y(), momentum.Z(), kaonEnergy);
    kaonMinus.SetPxPyPzE(-momentum.X(), -momentum.Y(), -momentum.Z(), kaonEnergy);
    kaonPlus.Boost(phi.BoostVector());
    kaonMinus.Boost(phi.BoostVector());
    return true;
  }

  TLorentzVector ScaleMomentum(const TLorentzVector &input, double epsilon)
  {
    const double factor = 1. + epsilon;
    const TVector3 momentum = factor * input.Vect();
    TLorentzVector output;
    output.SetVectM(momentum, kMassK);
    return output;
  }

  struct AverageResult
  {
    long long accepted = 0;
    double mean = 0.;
    double error = 0.;
  };

  struct ToyPeakFit
  {
    bool valid = false;
    int status = -999;
    int covQual = -1;
    double mean = 0.;
    double meanError = 0.;
    double sigma = 0.;
    double sigmaError = 0.;
  };

  ToyPeakFit FitToyPhiPeak(TH1D *histogram, double fitLow = 1.000,
                           double fitHigh = 1.03899,
                           double seedSigma = 0.0018)
  {
    ToyPeakFit result;
    if (!histogram || histogram->Integral() < 1000. || fitHigh <= fitLow)
      return result;
    static int fitCounter = 0;
    TF1 function(Form("fToyPhi_%d", fitCounter++),
                 "[0]*TMath::Voigt(x-[1],[2],[3],4)", fitLow, fitHigh);
    function.SetParNames("normalization", "mPhi", "sigmaPhi", "GammaPhi");
    const double peakDensity = TMath::Voigt(0., seedSigma, kWidthPhi, 4);
    function.SetParameters(histogram->GetMaximum() /
                               std::max(peakDensity, 1.e-12),
                           kMassPhiPDG, seedSigma,
                           kWidthPhi);
    function.SetParLimits(0, 0., 100. * histogram->GetMaximum());
    function.SetParLimits(1, 1.010, 1.030);
    function.SetParLimits(2, 1.e-6, 0.020);
    function.FixParameter(3, kWidthPhi);
    TVirtualFitter::SetDefaultFitter("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
    TFitResultPtr fit = histogram->Fit(&function, "Q0RSN");
    result.status = static_cast<int>(fit);
    result.covQual = fit.Get() ? fit->CovMatrixStatus() : -1;
    result.mean = function.GetParameter(1);
    result.meanError = function.GetParError(1);
    result.sigma = std::abs(function.GetParameter(2));
    result.sigmaError = function.GetParError(2);
    result.valid = fit.Get() && result.status == 0 && result.covQual >= 2 &&
                   std::isfinite(result.mean) &&
                   std::isfinite(result.meanError) && result.meanError > 0. &&
                   std::isfinite(result.sigma) && result.sigma > 0. &&
                   std::isfinite(result.sigmaError) && result.sigmaError > 0.;
    return result;
  }

  AverageResult AveragePhiMassShift(TH1D *ptSpectrum, double ptLow, double ptHigh,
                                    double epsilon, long long requestedAccepted,
                                    unsigned int seed, double yMax,
                                    double kaonPtMin, double kaonEtaMax)
  {
    AverageResult result;
    if (requestedAccepted <= 0)
      return result;
    TRandom3 random(seed);
    long long generated = 0;
    const long long maximumGenerated = std::max(1000LL, 100LL * requestedAccepted);
    long double sum = 0.;
    long double sum2 = 0.;
    static int toyCallCounter = 0;
    const int callIndex = toyCallCounter++;
    TH1D nominalHistogram(Form("hToyPhiNominalTemporary_%d", callIndex), "",
                          240, 0.990, 1.050);
    TH1D correctedHistogram(Form("hToyPhiCorrectedTemporary_%d", callIndex), "",
                            240, 0.990, 1.050);
    nominalHistogram.SetDirectory(nullptr);
    correctedHistogram.SetDirectory(nullptr);
    while (result.accepted < requestedAccepted && generated < maximumGenerated)
    {
      ++generated;
      const double phiPt = SamplePt(ptSpectrum, ptLow, ptHigh, random);
      double generatedPhiMass = 0.;
      for (int trial = 0; trial < 1000; ++trial)
      {
        generatedPhiMass = random.BreitWigner(kMassPhiPDG, kWidthPhi);
        if (generatedPhiMass > 2. * kMassK && generatedPhiMass >= 0.990 &&
            generatedPhiMass <= 1.050)
          break;
      }
      if (generatedPhiMass <= 2. * kMassK || generatedPhiMass < 0.990 ||
          generatedPhiMass > 1.050)
        continue;
      TLorentzVector kaonPlus;
      TLorentzVector kaonMinus;
      if (!GeneratePhiDecay(generatedPhiMass, phiPt, yMax, random,
                            kaonPlus, kaonMinus))
        continue;
      // Select the uncorrected 2025-like candidate. The analysis should reapply
      // its selections after the correction when producing final spectra.
      if (kaonPlus.Pt() < kaonPtMin || kaonMinus.Pt() < kaonPtMin)
        continue;
      if (std::abs(kaonPlus.Eta()) > kaonEtaMax ||
          std::abs(kaonMinus.Eta()) > kaonEtaMax)
        continue;

      const double nominalMass = (kaonPlus + kaonMinus).M();
      const TLorentzVector correctedPlus = ScaleMomentum(kaonPlus, epsilon);
      const TLorentzVector correctedMinus = ScaleMomentum(kaonMinus, epsilon);
      const double correctedMass = (correctedPlus + correctedMinus).M();
      const double shift = correctedMass - nominalMass;
      sum += shift;
      sum2 += shift * shift;
      nominalHistogram.Fill(nominalMass);
      correctedHistogram.Fill(correctedMass);
      ++result.accepted;
    }
    if (result.accepted > 0)
      result.mean = static_cast<double>(sum / result.accepted);
    if (result.accepted > 1)
    {
      const long double variance = std::max(0.L, sum2 / result.accepted -
                                                     (sum / result.accepted) *
                                                         (sum / result.accepted));
      result.error = std::sqrt(static_cast<double>(variance / result.accepted));
    }
    // The calibration observable is the fitted peak displacement, matching the
    // data procedure. Keep the event-average shift above only as a robust
    // fallback if a toy fit fails.
    const ToyPeakFit nominalFit = FitToyPhiPeak(&nominalHistogram);
    const ToyPeakFit correctedFit = FitToyPhiPeak(&correctedHistogram);
    if (nominalFit.valid && correctedFit.valid)
    {
      result.mean = correctedFit.mean - nominalFit.mean;
      result.error = std::hypot(correctedFit.meanError, nominalFit.meanError);
    }
    return result;
  }

  ToyPoint DetermineCorrection(TH1D *ptSpectrum, const FitPoint &fit2025,
                               const FitPoint &fit2026, int binIndex,
                               long long toyEvents, double yMax,
                               double kaonPtMin, double kaonEtaMax,
                               TDirectory *outputDirectory,
                               const std::string &responsePdf)
  {
    ToyPoint result;
    if (!fit2025.valid || !fit2026.valid)
      return result;
    result.targetShift = fit2026.mean - fit2025.mean;
    result.targetShiftError = std::hypot(fit2026.meanError, fit2025.meanError);

    const std::vector<double> epsilonScan = {-0.020, -0.010, -0.005,
                                             0.005, 0.010, 0.020};
    TGraphErrors response(static_cast<int>(epsilonScan.size()));
    response.SetName(Form("gPhiMassResponse_pt%02d", binIndex));
    response.SetTitle(Form("#phi mass response, %.3g < p_{T,#phi} < %.3g GeV/c;"
                           "#epsilon;#Delta m_{#phi}^{MC} (GeV/c^{2})",
                           fit2025.ptLow, fit2025.ptHigh));
    for (std::size_t i = 0; i < epsilonScan.size(); ++i)
    {
      const AverageResult point = AveragePhiMassShift(
          ptSpectrum, fit2025.ptLow, fit2025.ptHigh, epsilonScan[i], toyEvents,
          41711u + 1009u * binIndex, yMax, kaonPtMin, kaonEtaMax);
      response.SetPoint(static_cast<int>(i), epsilonScan[i], point.mean);
      response.SetPointError(static_cast<int>(i), 0., point.error);
      result.accepted = point.accepted;
    }

    TF1 linear(Form("fPhiMassResponse_pt%02d", binIndex), "[0]*x", -0.025, 0.025);
    linear.SetParameter(0, 0.05);
    linear.SetParName(0, "dmPhi_dEpsilon");
    TFitResultPtr linearFit = response.Fit(&linear, "Q0RS");
    result.slope = linear.GetParameter(0);
    result.slopeError = linear.GetParError(0);
    if (static_cast<int>(linearFit) != 0 || !std::isfinite(result.slope) ||
        std::abs(result.slope) < 1.e-8)
      return result;

    result.epsilonCorrection = result.targetShift / result.slope;
    const double relativeSlopeError = result.slopeError / result.slope;
    result.epsilonError = std::sqrt(
        std::pow(result.targetShiftError / result.slope, 2) +
        std::pow(result.epsilonCorrection * relativeSlopeError, 2));

    const AverageResult closure = AveragePhiMassShift(
        ptSpectrum, fit2025.ptLow, fit2025.ptHigh, result.epsilonCorrection,
        toyEvents, 41711u + 1009u * binIndex, yMax, kaonPtMin, kaonEtaMax);
    result.reproducedShift = closure.mean;
    result.reproducedShiftError = closure.error;
    result.valid = std::isfinite(result.epsilonCorrection) &&
                   std::abs(result.epsilonCorrection) < 0.10;

    TCanvas canvas(Form("cPhiMassResponse_pt%02d", binIndex),
                   "phi mass response", 900, 750);
    canvas.SetLeftMargin(0.13);
    response.SetMarkerStyle(20);
    response.SetMarkerColor(kBlue + 1);
    response.SetLineColor(kBlue + 1);
    response.Draw("AP");
    linear.SetLineColor(kRed + 1);
    linear.SetLineWidth(2);
    linear.Draw("SAME");
    TLine targetLine(-0.025, result.targetShift, 0.025, result.targetShift);
    targetLine.SetLineColor(kGreen + 2);
    targetLine.SetLineStyle(2);
    targetLine.Draw();
    TLine epsilonLine(result.epsilonCorrection, response.GetYaxis()->GetXmin(),
                      result.epsilonCorrection, response.GetYaxis()->GetXmax());
    epsilonLine.SetLineColor(kMagenta + 2);
    epsilonLine.SetLineStyle(2);
    epsilonLine.Draw();
    TLatex annotation;
    annotation.SetNDC();
    annotation.SetTextSize(0.034);
    annotation.DrawLatex(0.16, 0.88,
                         Form("target m_{#phi}^{2026}-m_{#phi}^{2025}=%.3f MeV/c^{2}",
                              1000. * result.targetShift));
    annotation.DrawLatex(0.16, 0.83,
                         Form("#epsilon_{corr}=%.4f #pm %.4f %%",
                              100. * result.epsilonCorrection, 100. * result.epsilonError));
    canvas.Print(responsePdf.c_str());

    if (outputDirectory)
    {
      outputDirectory->cd();
      response.Write();
      linear.Write();
      canvas.Write();
    }
    return result;
  }

  TGraphErrors *MakeGraph(const char *name, const char *title,
                          const std::vector<FitPoint> &points,
                          bool useSigma, double scale)
  {
    TGraphErrors *graph = new TGraphErrors();
    graph->SetName(name);
    graph->SetTitle(title);
    int n = 0;
    for (const FitPoint &point : points)
    {
      if (!point.valid)
        continue;
      const double x = 0.5 * (point.ptLow + point.ptHigh);
      const double ex = 0.5 * (point.ptHigh - point.ptLow);
      const double y = scale * (useSigma ? point.sigma : point.mean);
      const double ey = scale * (useSigma ? point.sigmaError : point.meanError);
      graph->SetPoint(n, x, y);
      graph->SetPointError(n, ex, ey);
      ++n;
    }
    return graph;
  }

  void StyleGraph(TGraphErrors *graph, int color, int marker)
  {
    if (!graph)
      return;
    graph->SetLineColor(color);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(marker);
    graph->SetMarkerSize(1.0);
    graph->SetLineWidth(2);
  }

  struct XResolutionPoint
  {
    double pt = 0.;
    double ptError = 0.;
    double sigmaPhi = 0.;
    double sigmaPhiError = 0.;
    double relativeKaonResolution = 0.;
    double relativeKaonResolutionError = 0.;
    double closureSigmaPhi = 0.;
    double closureSigmaPhiError = 0.;
  };

  struct XResolutionResult
  {
    bool valid = false;
    int fitStatus = -999;
    int covQual = -1;
    long long accepted = 0;
    long long sampledPt = 0;
    double meanMeV = 0.;
    double meanErrorMeV = 0.;
    double sigmaMeV = 0.;
    double sigmaErrorMeV = 0.;
    double rmsMeV = 0.;
    int observedFitStatus = -999;
    int observedCovQual = -1;
    double observedMeanMeV = 0.;
    double observedMeanErrorMeV = 0.;
    double observedSigmaMeV = 0.;
    double observedSigmaErrorMeV = 0.;
    double observedFwhmMeV = 0.;
  };

  double SamplePhysicalPhiMass(TRandom3 &random)
  {
    for (int trial = 0; trial < 1000; ++trial)
    {
      const double mass = random.BreitWigner(kMassPhiPDG, kWidthPhi);
      if (mass > 2. * kMassK && mass < kMassPhiPDG + 10. * kWidthPhi)
        return mass;
    }
    return kMassPhiPDG;
  }

  double SamplePhysicalXMass(TRandom3 &random, double poleMass,
                             double intrinsicWidth, double minimumMass,
                             double maximumMass)
  {
    if (maximumMass <= minimumMass)
      return 0.;
    if (intrinsicWidth <= 0.)
      return poleMass > minimumMass && poleMass < maximumMass ? poleMass : 0.;
    const double low = std::max(minimumMass,
                                poleMass - 10. * intrinsicWidth);
    const double high = std::min(maximumMass,
                                 poleMass + 10. * intrinsicWidth);
    if (high <= low)
      return 0.;
    for (int trial = 0; trial < 10000; ++trial)
    {
      const double mass = random.BreitWigner(poleMass, intrinsicWidth);
      if (mass > low && mass < high)
        return mass;
    }
    return poleMass > low && poleMass < high ? poleMass : 0.;
  }

  TLorentzVector FourVectorFromPtY(double pt, double rapidity, double azimuth,
                                   double mass)
  {
    const double transverseMass = std::sqrt(mass * mass + pt * pt);
    TLorentzVector vector;
    vector.SetPxPyPzE(pt * std::cos(azimuth), pt * std::sin(azimuth),
                      transverseMass * std::sinh(rapidity),
                      transverseMass * std::cosh(rapidity));
    return vector;
  }

  bool TwoBodyDecay(const TLorentzVector &parent, double daughterMass1,
                    double daughterMass2, TRandom3 &random,
                    TLorentzVector &daughter1, TLorentzVector &daughter2)
  {
    const double parentMass = parent.M();
    if (parentMass <= daughterMass1 + daughterMass2)
      return false;
    const double first = parentMass * parentMass -
                         std::pow(daughterMass1 + daughterMass2, 2);
    const double second = parentMass * parentMass -
                          std::pow(daughterMass1 - daughterMass2, 2);
    if (first <= 0. || second <= 0.)
      return false;
    const double momentum = std::sqrt(first * second) / (2. * parentMass);
    const double cosTheta = random.Uniform(-1., 1.);
    const double sinTheta = std::sqrt(std::max(0., 1. - cosTheta * cosTheta));
    const double azimuth = random.Uniform(0., TMath::TwoPi());
    const TVector3 p(momentum * sinTheta * std::cos(azimuth),
                     momentum * sinTheta * std::sin(azimuth),
                     momentum * cosTheta);
    daughter1.SetVectM(p, daughterMass1);
    daughter2.SetVectM(-p, daughterMass2);
    daughter1.Boost(parent.BoostVector());
    daughter2.Boost(parent.BoostVector());
    return true;
  }

  TLorentzVector SmearKaonMomentum(const TLorentzVector &truth,
                                   double relativeResolution,
                                   double gaussianPull)
  {
    const double scale = std::max(1.e-6,
                                  1. + relativeResolution * gaussianPull);
    TLorentzVector reconstructed;
    reconstructed.SetVectM(scale * truth.Vect(), kMassK);
    return reconstructed;
  }

  double SimulatePhiDetectorSigma(double phiPt, double relativeKaonResolution,
                                  long long requestedAccepted, double yMax,
                                  double kaonPtMin, double kaonEtaMax,
                                  unsigned int seed, double fitLow,
                                  double fitHigh, double &sigmaError)
  {
    sigmaError = 0.;
    TRandom3 random(seed);
    long long accepted = 0;
    long long generated = 0;
    static int histogramCounter = 0;
    TH1D reconstructedMass(
        Form("hPhiResolutionToyTemporary_%d", histogramCounter++), "",
        240, 0.990, 1.050);
    reconstructedMass.SetDirectory(nullptr);
    const long long maximumGenerated = std::max(
        requestedAccepted + 1000LL, 100LL * requestedAccepted);
    while (accepted < requestedAccepted && generated < maximumGenerated)
    {
      ++generated;
      const double generatedMass = SamplePhysicalPhiMass(random);
      TLorentzVector kaonPlus;
      TLorentzVector kaonMinus;
      if (!GeneratePhiDecay(generatedMass, phiPt, yMax, random,
                            kaonPlus, kaonMinus))
        continue;
      const TLorentzVector kaonPlusReco = SmearKaonMomentum(
          kaonPlus, relativeKaonResolution, random.Gaus(0., 1.));
      const TLorentzVector kaonMinusReco = SmearKaonMomentum(
          kaonMinus, relativeKaonResolution, random.Gaus(0., 1.));
      if (kaonPlusReco.Pt() < kaonPtMin ||
          kaonMinusReco.Pt() < kaonPtMin ||
          std::abs(kaonPlusReco.Eta()) > kaonEtaMax ||
          std::abs(kaonMinusReco.Eta()) > kaonEtaMax)
        continue;
      const TLorentzVector reconstructed = kaonPlusReco + kaonMinusReco;
      reconstructedMass.Fill(reconstructed.M());
      ++accepted;
    }
    if (accepted < 1000)
      return -1.;

    // This is the same observable used in data: fit M(KK)_reco with a Voigtian
    // whose Lorentzian FWHM is fixed, and return its Gaussian sigma.  Do not use
    // RMS[M(KK)_reco-Mphi_true], which is a different resolution definition.
    const ToyPeakFit peak = FitToyPhiPeak(
        &reconstructedMass, fitLow, fitHigh, 0.0018);
    if (!peak.valid)
      return -1.;
    sigmaError = peak.sigmaError;
    return peak.sigma;
  }

  double InferKaonResolution(double phiPt, double targetSigmaPhi,
                             long long toyEvents, double yMax,
                             double kaonPtMin, double kaonEtaMax,
                             unsigned int seed, double fitLow, double fitHigh,
                             double &closureSigmaPhi,
                             double &closureSigmaPhiError)
  {
    if (targetSigmaPhi <= 0.)
    {
      closureSigmaPhi = 0.;
      closureSigmaPhiError = 0.;
      return 0.;
    }
    double low = 0.;
    double high = 0.01;
    double temporaryError = 0.;
    double sigmaHigh = SimulatePhiDetectorSigma(
        phiPt, high, toyEvents, yMax, kaonPtMin, kaonEtaMax, seed,
        fitLow, fitHigh, temporaryError);
    while ((sigmaHigh <= 0. || sigmaHigh < targetSigmaPhi) && high < 0.50)
    {
      high *= 2.;
      sigmaHigh = SimulatePhiDetectorSigma(
          phiPt, high, toyEvents, yMax, kaonPtMin, kaonEtaMax, seed,
          fitLow, fitHigh, temporaryError);
    }
    if (sigmaHigh <= 0.)
    {
      closureSigmaPhi = 0.;
      closureSigmaPhiError = 0.;
      return 0.;
    }
    for (int iteration = 0; iteration < 18; ++iteration)
    {
      const double middle = 0.5 * (low + high);
      const double sigmaMiddle = SimulatePhiDetectorSigma(
          phiPt, middle, toyEvents, yMax, kaonPtMin, kaonEtaMax, seed,
          fitLow, fitHigh, temporaryError);
      if (sigmaMiddle > 0. && sigmaMiddle < targetSigmaPhi)
        low = middle;
      else
        high = middle;
    }
    const double result = 0.5 * (low + high);
    closureSigmaPhi = SimulatePhiDetectorSigma(
        phiPt, result, 2 * toyEvents, yMax, kaonPtMin, kaonEtaMax,
        seed + 100003U, fitLow, fitHigh, closureSigmaPhiError);
    return result;
  }

  double InterpolateKaonResolution(const std::vector<XResolutionPoint> &points,
                                   double phiPt)
  {
    if (points.empty())
      return 0.;
    if (phiPt <= points.front().pt)
      return points.front().relativeKaonResolution;
    if (phiPt >= points.back().pt)
      return points.back().relativeKaonResolution;
    for (std::size_t i = 1; i < points.size(); ++i)
    {
      if (phiPt > points[i].pt)
        continue;
      const double fraction = (phiPt - points[i - 1].pt) /
                              (points[i].pt - points[i - 1].pt);
      return points[i - 1].relativeKaonResolution + fraction *
                                                        (points[i].relativeKaonResolution -
                                                         points[i - 1].relativeKaonResolution);
    }
    return points.back().relativeKaonResolution;
  }

  THnSparse *GetPairSparse(TFile *file, const std::string &pathOrName)
  {
    if (!file || file->IsZombie())
      return nullptr;
    TObject *object = file->Get(pathOrName.c_str());
    if (!object)
    {
      const std::size_t slash = pathOrName.find_last_of('/');
      const std::string shortName = slash == std::string::npos
                                        ? pathOrName
                                        : pathOrName.substr(slash + 1);
      object = FindObjectRecursive(file, shortName);
    }
    THnSparse *sparse = dynamic_cast<THnSparse *>(object);
    if (!sparse)
      std::cerr << "ERROR: THnSparse '" << pathOrName << "' was not found in "
                << file->GetName() << std::endl;
    return sparse;
  }

  TH1D *ProjectSelectedXPt(THnSparse *sparse, int massAxis, int ptAxis,
                           double massLow, double massHigh, double ptMin,
                           const std::string &name)
  {
    if (!sparse || massAxis < 0 || ptAxis < 0 ||
        massAxis >= sparse->GetNdimensions() ||
        ptAxis >= sparse->GetNdimensions() || massAxis == ptAxis)
      return nullptr;
    TAxis *mAxis = sparse->GetAxis(massAxis);
    TAxis *pAxis = sparse->GetAxis(ptAxis);
    const int oldMassFirst = mAxis->GetFirst();
    const int oldMassLast = mAxis->GetLast();
    const int oldPtFirst = pAxis->GetFirst();
    const int oldPtLast = pAxis->GetLast();
    mAxis->SetRangeUser(massLow + 1.e-9, massHigh - 1.e-9);
    pAxis->SetRangeUser(ptMin + 1.e-9, pAxis->GetXmax() - 1.e-9);
    TH1D *result = dynamic_cast<TH1D *>(sparse->Projection(ptAxis, "E"));
    mAxis->SetRange(oldMassFirst, oldMassLast);
    pAxis->SetRange(oldPtFirst, oldPtLast);
    if (!result)
      return nullptr;
    result->SetDirectory(nullptr);
    result->SetName(name.c_str());
    result->SetTitle("Reconstructed selected X p_{T} shape;"
                     "p_{T,X} (GeV/c);Candidates");
    return result;
  }

  XResolutionResult EstimateXResolution(
      const std::string &label, TFile *inputFile,
      const std::string &sparsePathOrName, int sparseMassAxis, int sparsePtAxis,
      double shapeMassLow, double shapeMassHigh, double xPtMin, double xMass,
      double xIntrinsicWidth, double deltaMassMax, double phiFitLow,
      double phiFitHigh,
      const std::vector<FitPoint> &phiFits,
      long long phiCalibrationEvents,
      long long xAcceptedEvents, double phiPtMin, double phiPtMax,
      double phiRapidityMax, double xRapidityMax, double kaonPtMin,
      double kaonEtaMax, unsigned int seed, TDirectory *outputDirectory,
      const std::string &pdfName)
  {
    XResolutionResult result;
    std::vector<XResolutionPoint> resolutionPoints;
    for (const FitPoint &fit : phiFits)
    {
      if (!fit.valid || fit.sigma <= 0.)
        continue;
      XResolutionPoint point;
      point.pt = 0.5 * (fit.ptLow + fit.ptHigh);
      point.ptError = 0.5 * (fit.ptHigh - fit.ptLow);
      point.sigmaPhi = fit.sigma;
      point.sigmaPhiError = fit.sigmaError;
      point.relativeKaonResolution = InferKaonResolution(
          point.pt, point.sigmaPhi, phiCalibrationEvents, phiRapidityMax,
          kaonPtMin, kaonEtaMax,
          seed + 1009U * static_cast<unsigned int>(resolutionPoints.size()),
          phiFitLow, phiFitHigh, point.closureSigmaPhi,
          point.closureSigmaPhiError);
      if (point.relativeKaonResolution <= 0. ||
          point.closureSigmaPhi <= 0.)
      {
        std::cerr << "WARNING: phi-resolution inversion failed for " << label
                  << " at pT=" << point.pt << " GeV/c" << std::endl;
        continue;
      }
      point.relativeKaonResolutionError = point.sigmaPhi > 0.
                                              ? point.relativeKaonResolution * point.sigmaPhiError / point.sigmaPhi
                                              : 0.;
      resolutionPoints.push_back(point);
    }
    if (resolutionPoints.size() < 2)
    {
      std::cerr << "ERROR: insufficient phi-resolution points for X toy "
                << label << std::endl;
      return result;
    }

    THnSparse *sparse = GetPairSparse(inputFile, sparsePathOrName);
    if (!sparse)
      return result;
    std::unique_ptr<TH1D> xPtShape(ProjectSelectedXPt(
        sparse, sparseMassAxis, sparsePtAxis, shapeMassLow, shapeMassHigh,
        xPtMin, "hSelectedXPtShape_" + SafeName(label)));
    if (!xPtShape || xPtShape->Integral() <= 0.)
    {
      std::cerr << "ERROR: empty selected X-pT projection for " << label
                << std::endl;
      return result;
    }

    TGraphErrors effectiveKaonResolution;
    effectiveKaonResolution.SetName(
        ("gEffectiveKaonResolution_" + SafeName(label)).c_str());
    effectiveKaonResolution.SetTitle(
        (label + ";p_{T,#phi} (GeV/c);#sigma(p_{K})/p_{K} (%)").c_str());
    for (std::size_t i = 0; i < resolutionPoints.size(); ++i)
    {
      const XResolutionPoint &point = resolutionPoints[i];
      effectiveKaonResolution.SetPoint(
          i, point.pt, 100. * point.relativeKaonResolution);
      effectiveKaonResolution.SetPointError(
          i, point.ptError, 100. * point.relativeKaonResolutionError);
    }
    StyleGraph(&effectiveKaonResolution, kMagenta + 2, 20);

    TGraphErrors phiResolutionTarget;
    TGraphErrors phiResolutionClosure;
    phiResolutionTarget.SetName(
        ("gPhiResolutionTarget_" + SafeName(label)).c_str());
    phiResolutionClosure.SetName(
        ("gPhiResolutionClosure_" + SafeName(label)).c_str());
    phiResolutionTarget.SetTitle(
        (label + ";p_{T,#phi} (GeV/c);#sigma_{#phi} (MeV/c^{2})").c_str());
    phiResolutionClosure.SetTitle(phiResolutionTarget.GetTitle());
    for (std::size_t i = 0; i < resolutionPoints.size(); ++i)
    {
      const XResolutionPoint &point = resolutionPoints[i];
      phiResolutionTarget.SetPoint(i, point.pt, 1000. * point.sigmaPhi);
      phiResolutionTarget.SetPointError(
          i, point.ptError, 1000. * point.sigmaPhiError);
      phiResolutionClosure.SetPoint(
          i, point.pt, 1000. * point.closureSigmaPhi);
      phiResolutionClosure.SetPointError(
          i, point.ptError, 1000. * point.closureSigmaPhiError);
      std::cout << "  " << label << " phi closure pT=" << point.pt
                << " GeV/c: target sigma=" << 1000. * point.sigmaPhi
                << " MeV, toy-fit sigma=" << 1000. * point.closureSigmaPhi
                << " MeV, sigma(pK)/pK="
                << 100. * point.relativeKaonResolution << "%" << std::endl;
    }
    StyleGraph(&phiResolutionTarget, kBlack, 20);
    StyleGraph(&phiResolutionClosure, kBlue + 1, 24);

    TH1D residual(("hXMassResidual_" + SafeName(label)).c_str(),
                  (label + ";M_{4K}^{reco}-M_{X}^{true} (MeV/c^{2});Events").c_str(),
                  500, -100., 100.);
    residual.SetDirectory(nullptr);
    TH1D generatedMass(("hXGeneratedMass_" + SafeName(label)).c_str(),
                       (label + ";M_{X}^{true} (GeV/c^{2});Events").c_str(),
                       400, shapeMassLow, shapeMassHigh);
    TH1D reconstructedMass(
        ("hXReconstructedMass_" + SafeName(label)).c_str(),
        (label + ";M_{4K}^{reco} (GeV/c^{2});Events").c_str(),
        400, shapeMassLow, shapeMassHigh);
    generatedMass.SetDirectory(nullptr);
    reconstructedMass.SetDirectory(nullptr);
    TRandom3 random(seed + 700001U);
    const double ptMaximum = xPtShape->GetXaxis()->GetXmax();
    while (result.accepted < xAcceptedEvents &&
           result.sampledPt < 20LL * xAcceptedEvents)
    {
      ++result.sampledPt;
      // This pT is already a reconstructed, fully selected distribution. Hold
      // the draw fixed while finding an accepted phase-space decay, so the cuts
      // below do not apply an additional acceptance weight to the input pT shape.
      const double xPt = SamplePt(xPtShape.get(), xPtMin, ptMaximum, random);
      bool acceptedCandidate = false;
      double xResidualMeV = 0.;
      double xTrueMass = 0.;
      double xRecoMass = 0.;
      for (int topologyTrial = 0; topologyTrial < 5000; ++topologyTrial)
      {
        const double phiMass1 = SamplePhysicalPhiMass(random);
        const double phiMass2 = SamplePhysicalPhiMass(random);
        xTrueMass = SamplePhysicalXMass(
            random, xMass, xIntrinsicWidth,
            std::max(shapeMassLow, phiMass1 + phiMass2 + 1.e-9),
            shapeMassHigh);
        if (xTrueMass <= 0.)
          continue;
        const TLorentzVector x = FourVectorFromPtY(
            xPt, random.Uniform(-xRapidityMax, xRapidityMax),
            random.Uniform(0., TMath::TwoPi()), xTrueMass);
        TLorentzVector phi1;
        TLorentzVector phi2;
        if (!TwoBodyDecay(x, phiMass1, phiMass2, random, phi1, phi2))
          continue;
        TLorentzVector k1Plus;
        TLorentzVector k1Minus;
        TLorentzVector k2Plus;
        TLorentzVector k2Minus;
        if (!TwoBodyDecay(phi1, kMassK, kMassK, random, k1Plus, k1Minus) ||
            !TwoBodyDecay(phi2, kMassK, kMassK, random, k2Plus, k2Minus))
          continue;
        const double resolution1 =
            InterpolateKaonResolution(resolutionPoints, phi1.Pt());
        const double resolution2 =
            InterpolateKaonResolution(resolutionPoints, phi2.Pt());
        const TLorentzVector k1PlusReco = SmearKaonMomentum(
            k1Plus, resolution1, random.Gaus(0., 1.));
        const TLorentzVector k1MinusReco = SmearKaonMomentum(
            k1Minus, resolution1, random.Gaus(0., 1.));
        const TLorentzVector k2PlusReco = SmearKaonMomentum(
            k2Plus, resolution2, random.Gaus(0., 1.));
        const TLorentzVector k2MinusReco = SmearKaonMomentum(
            k2Minus, resolution2, random.Gaus(0., 1.));
        const TLorentzVector phi1Reco = k1PlusReco + k1MinusReco;
        const TLorentzVector phi2Reco = k2PlusReco + k2MinusReco;

        // Apply all daughter selections to reconstructed objects, after detector
        // smearing and before filling the X response.  The sampled pT(X) value is
        // held fixed while regenerating the topology, so its reconstructed,
        // already-selected THnSparse shape is not acceptance-weighted again.
        if (phi1Reco.Pt() < phiPtMin || phi1Reco.Pt() > phiPtMax ||
            phi2Reco.Pt() < phiPtMin || phi2Reco.Pt() > phiPtMax ||
            std::abs(phi1Reco.Rapidity()) > phiRapidityMax ||
            std::abs(phi2Reco.Rapidity()) > phiRapidityMax)
          continue;
        if (k1PlusReco.Pt() < kaonPtMin || k1MinusReco.Pt() < kaonPtMin ||
            k2PlusReco.Pt() < kaonPtMin || k2MinusReco.Pt() < kaonPtMin ||
            std::abs(k1PlusReco.Eta()) > kaonEtaMax ||
            std::abs(k1MinusReco.Eta()) > kaonEtaMax ||
            std::abs(k2PlusReco.Eta()) > kaonEtaMax ||
            std::abs(k2MinusReco.Eta()) > kaonEtaMax)
          continue;

        // Same circular double-phi mass selection as in data.
        const double deltaMass = std::hypot(phi1Reco.M() - kMassPhiPDG,
                                            phi2Reco.M() - kMassPhiPDG);
        if (deltaMass >= deltaMassMax)
          continue;

        xRecoMass = (phi1Reco + phi2Reco).M();
        xResidualMeV = 1000. * (xRecoMass - xTrueMass);
        acceptedCandidate = true;
        break;
      }
      if (!acceptedCandidate)
        continue;
      residual.Fill(xResidualMeV);
      generatedMass.Fill(xTrueMass);
      reconstructedMass.Fill(xRecoMass);
      ++result.accepted;
    }

    result.rmsMeV = residual.GetRMS();
    const double fitHalfWidth = std::max(2., 2.5 * result.rmsMeV);
    TF1 gaussian(("fXResolution_" + SafeName(label)).c_str(), "gaus",
                 residual.GetMean() - fitHalfWidth,
                 residual.GetMean() + fitHalfWidth);
    gaussian.SetParameters(residual.GetMaximum(), residual.GetMean(),
                           std::max(0.5, result.rmsMeV));
    gaussian.SetLineColor(kRed + 1);
    gaussian.SetLineWidth(3);
    TVirtualFitter::SetDefaultFitter("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
    TFitResultPtr fitResult = residual.Fit(&gaussian, "Q0RSN");
    result.fitStatus = static_cast<int>(fitResult);
    result.covQual = fitResult.Get() ? fitResult->CovMatrixStatus() : -1;
    result.meanMeV = gaussian.GetParameter(1);
    result.meanErrorMeV = gaussian.GetParError(1);
    result.sigmaMeV = std::abs(gaussian.GetParameter(2));
    result.sigmaErrorMeV = gaussian.GetParError(2);
    result.valid = fitResult.Get() && result.fitStatus == 0 &&
                   result.covQual >= 2 && result.sigmaMeV > 0.;

    // Fit the observed reconstructed peak separately.  Gamma_X is the
    // Lorentzian FWHM and is fixed; the fitted Gaussian sigma is the detector
    // term after integrating over the generated line shape and selections.
    const double observedFitLow = std::max(shapeMassLow,
                                           xMass - 5. * xIntrinsicWidth);
    const double observedFitHigh = std::min(shapeMassHigh,
                                            xMass + 5. * xIntrinsicWidth);
    TF1 observedModel(("fXObservedVoigt_" + SafeName(label)).c_str(),
                      "[0]*TMath::Voigt(x-[1],[2],[3],4)",
                      observedFitLow, observedFitHigh);
    observedModel.SetParNames("normalization", "mX", "sigmaX", "GammaX");
    const double observedSeedSigma = std::max(0.0005,
                                              result.sigmaMeV / 1000.);
    const double observedPeakDensity = TMath::Voigt(
        0., observedSeedSigma, xIntrinsicWidth, 4);
    observedModel.SetParameters(
        reconstructedMass.GetMaximum() /
            std::max(observedPeakDensity, 1.e-12),
        xMass, observedSeedSigma, xIntrinsicWidth);
    observedModel.SetParLimits(0, 0., 100. * reconstructedMass.GetMaximum());
    observedModel.SetParLimits(1, observedFitLow, observedFitHigh);
    observedModel.SetParLimits(2, 1.e-5, 0.050);
    observedModel.FixParameter(3, xIntrinsicWidth);
    observedModel.SetLineColor(kRed + 1);
    observedModel.SetLineWidth(3);
    TFitResultPtr observedFitResult = reconstructedMass.Fit(
        &observedModel, "Q0RSN");
    result.observedFitStatus = static_cast<int>(observedFitResult);
    result.observedCovQual = observedFitResult.Get()
                                 ? observedFitResult->CovMatrixStatus()
                                 : -1;
    result.observedMeanMeV = 1000. * observedModel.GetParameter(1);
    result.observedMeanErrorMeV = 1000. * observedModel.GetParError(1);
    result.observedSigmaMeV = 1000. * std::abs(observedModel.GetParameter(2));
    result.observedSigmaErrorMeV = 1000. * observedModel.GetParError(2);
    const double gammaMeV = 1000. * xIntrinsicWidth;
    const double gaussianFwhmMeV = 2.354820045 * result.observedSigmaMeV;
    result.observedFwhmMeV = 0.5346 * gammaMeV +
                             std::sqrt(0.2166 * gammaMeV * gammaMeV +
                                       gaussianFwhmMeV * gaussianFwhmMeV);

    TCanvas canvas(("cXResolution_" + SafeName(label)).c_str(),
                   "X detector response", 2200, 600);
    canvas.Divide(4, 1);
    canvas.cd(1);
    gPad->SetLeftMargin(0.14);
    phiResolutionTarget.Draw("AP");
    phiResolutionClosure.Draw("P SAME");
    TLegend closureLegend(0.48, 0.75, 0.88, 0.88);
    closureLegend.SetBorderSize(0);
    closureLegend.AddEntry(&phiResolutionTarget, "data Voigt #sigma", "lep");
    closureLegend.AddEntry(&phiResolutionClosure, "toy Voigt #sigma", "lep");
    closureLegend.Draw();
    canvas.cd(2);
    gPad->SetLeftMargin(0.14);
    effectiveKaonResolution.Draw("AP");
    canvas.cd(3);
    gPad->SetLeftMargin(0.14);
    residual.SetMarkerStyle(20);
    residual.SetMarkerSize(0.65);
    residual.Draw("E");
    gaussian.Draw("L SAME");
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.036);
    text.DrawLatex(0.16, 0.88,
                   Form("%s: #sigma_{X}=%.3f #pm %.3f MeV/c^{2}",
                        label.c_str(), result.sigmaMeV,
                        result.sigmaErrorMeV));
    text.DrawLatex(0.16, 0.82,
                   Form("m_{X}=%.3f GeV, #Gamma_{X}=%.1f MeV",
                        xMass, 1000. * xIntrinsicWidth));
    text.DrawLatex(0.16, 0.76,
                   "p_{T,X} input: reconstructed and selected THnSparse shape");
    text.DrawLatex(0.16, 0.70,
                   Form("#DeltaM_{#phi}<%.1f MeV/c^{2} after smearing",
                        1000. * deltaMassMax));
    canvas.cd(4);
    gPad->SetLeftMargin(0.14);
    reconstructedMass.SetMarkerStyle(20);
    reconstructedMass.SetMarkerSize(0.55);
    reconstructedMass.SetLineColor(kBlack);
    generatedMass.SetLineColor(kBlue + 1);
    generatedMass.SetLineWidth(2);
    reconstructedMass.Draw("E");
    generatedMass.Draw("HIST SAME");
    observedModel.Draw("L SAME");
    TLegend observedLegend(0.49, 0.72, 0.88, 0.88);
    observedLegend.SetBorderSize(0);
    observedLegend.AddEntry(&reconstructedMass, "reconstructed", "lep");
    observedLegend.AddEntry(&generatedMass, "generated BW", "l");
    observedLegend.AddEntry(&observedModel, "Voigt fit", "l");
    observedLegend.Draw();
    TLatex observedText;
    observedText.SetNDC();
    observedText.SetTextSize(0.032);
    observedText.DrawLatex(
        0.16, 0.66,
        Form("#sigma_{det}^{Voigt}=%.2f #pm %.2f MeV/c^{2}",
             result.observedSigmaMeV, result.observedSigmaErrorMeV));
    observedText.DrawLatex(
        0.16, 0.60,
        Form("Voigt FWHM #approx %.2f MeV/c^{2}",
             result.observedFwhmMeV));
    canvas.Print(pdfName.c_str());

    if (outputDirectory)
    {
      outputDirectory->cd();
      xPtShape->Write();
      phiResolutionTarget.Write();
      phiResolutionClosure.Write();
      effectiveKaonResolution.Write();
      residual.Write();
      generatedMass.Write();
      reconstructedMass.Write();
      gaussian.Write();
      observedModel.Write();
      canvas.Write();
      if (fitResult.Get())
        fitResult->Write(("fitResultXResolution_" + SafeName(label)).c_str());
      if (observedFitResult.Get())
        observedFitResult->Write(
            ("fitResultXObservedVoigt_" + SafeName(label)).c_str());
    }
    return result;
  }

  void WriteCorrectionHeader(const std::string &path,
                             const std::vector<double> &ptEdges,
                             const std::vector<ToyPoint> &toyPoints)
  {
    std::ofstream output(path);
    output << "#ifndef PHI_MOMENTUM_SCALE_2025_TO_2026_H\n";
    output << "#define PHI_MOMENTUM_SCALE_2025_TO_2026_H\n\n";
    output << "// Generated by CalibratePhiMomentumScaleVsPt.C\n";
    output << "// epsilon is applied to both kaon three-momenta of a phi candidate:\n";
    output << "//   p_corr = (1 + epsilon) * p_2025\n";
    output << "// Recalculate each kaon energy using the fixed kaon mass afterwards.\n\n";
    output << "namespace PhiMomentumScale2025To2026 {\n";
    output << "constexpr int kNBins = " << toyPoints.size() << ";\n";
    output << std::setprecision(12);
    output << "constexpr double kPtEdges[kNBins + 1] = {";
    for (std::size_t i = 0; i < ptEdges.size(); ++i)
    {
      if (i)
        output << ", ";
      output << ptEdges[i];
    }
    output << "};\n";
    output << "constexpr double kEpsilon[kNBins] = {";
    for (std::size_t i = 0; i < toyPoints.size(); ++i)
    {
      if (i)
        output << ", ";
      output << (toyPoints[i].valid ? toyPoints[i].epsilonCorrection : 0.0);
    }
    output << "};\n";
    output << "constexpr bool kValid[kNBins] = {";
    for (std::size_t i = 0; i < toyPoints.size(); ++i)
    {
      if (i)
        output << ", ";
      output << (toyPoints[i].valid ? "true" : "false");
    }
    output << "};\n\n";
    output << "inline double Epsilon(double reconstructedPhiPt)\n";
    output << "{\n";
    output << "  for (int i = 0; i < kNBins; ++i) {\n";
    output << "    if (reconstructedPhiPt >= kPtEdges[i] &&\n";
    output << "        reconstructedPhiPt < kPtEdges[i + 1])\n";
    output << "      return kValid[i] ? kEpsilon[i] : 0.0;\n";
    output << "  }\n";
    output << "  return 0.0; // no extrapolation outside the calibrated range\n";
    output << "}\n\n";
    output << "inline double ScaleFactor(double reconstructedPhiPt)\n";
    output << "{\n";
    output << "  return 1.0 + Epsilon(reconstructedPhiPt);\n";
    output << "}\n";
    output << "} // namespace PhiMomentumScale2025To2026\n\n";
    output << "#endif\n";
  }

} // namespace PhiPtCalibration

void CalibratePhiMomentumScaleVsPt(
    const char *file2026Name = "/home/sawan/alice/practice/OutputDoublePhi/New/processopti5/AnalysisResults_PID2003.root",
    const char *file2025Name = "/home/sawan/alice/practice/OutputDoublePhi/New/processopti5/LHC25/AnalysisResults_LHC25_PID2003.root",
    // const char *file2026Name = "/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/AnalysisResults_LHC26_PID2003.root",
    // const char *file2025Name = "/home/sawan/alice/practice/OutputDoublePhi/New/processopti8/LHC25/AnalysisResults_LHC25_PID2003.root",
    const char *histogramPathOrName = "hPhiMassVsPt",
    const char *outputDirectoryName = "PhiMomentumScaleVsPt",
    const char *ptBins = "0.8,1,1.2,1.6,2,2.5,3,4,5,6,7,8,9,10,12,14,16,20",
    long long toyAcceptedEventsPerBin = 150000,
    double fitLow = 1.000,
    double fitHigh = 1.03899,
    double yMax = 0.8,
    double kaonPtMin = 0.15,
    double kaonEtaMax = 0.8,
    bool runXResolutionToy = true,
    const char *xSparsePathOrName = "SEMassDoublePhi",
    int xSparseMassAxis = 0,
    int xSparsePtAxis = 1,
    double xShapeMassLow = 2.6,
    double xShapeMassHigh = 2.8,
    double xPtMin = 9.0,
    double xMass = 2.690,
    long long phiResolutionToyEventsPerPoint = 50000,
    long long xResolutionAcceptedEvents = 300000,
    double xRapidityGenerationMax = 0.8,
    unsigned int xResolutionSeed = 735193U,
    double xDeltaMassMax = 0.005,
    double xIntrinsicWidth = 0.012)
{
  using namespace PhiPtCalibration;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  const std::vector<double> ptEdges = ParseEdges(ptBins);
  if (ptEdges.size() < 2)
  {
    std::cerr << "ERROR: at least two pT edges are required." << std::endl;
    return;
  }
  if (fitHigh <= fitLow)
  {
    std::cerr << "ERROR: invalid mass fit range." << std::endl;
    return;
  }

  std::unique_ptr<TFile> file2025(TFile::Open(file2025Name, "READ"));
  std::unique_ptr<TFile> file2026(TFile::Open(file2026Name, "READ"));
  if (!file2025 || file2025->IsZombie() || !file2026 || file2026->IsZombie())
  {
    std::cerr << "ERROR: could not open one or both input files." << std::endl;
    return;
  }
  TH2 *h2_2025 = GetMassVsPt(file2025.get(), histogramPathOrName);
  TH2 *h2_2026 = GetMassVsPt(file2026.get(), histogramPathOrName);
  if (!h2_2025 || !h2_2026)
    return;

  gSystem->mkdir(outputDirectoryName, true);
  const std::string outputDirectory(outputDirectoryName);
  const std::string outputRoot = outputDirectory + "/PhiMomentumScaleVsPt.root";
  std::unique_ptr<TFile> output(TFile::Open(outputRoot.c_str(), "RECREATE"));
  if (!output || output->IsZombie())
  {
    std::cerr << "ERROR: could not create " << outputRoot << std::endl;
    return;
  }
  TDirectory *directory2025 = output->mkdir("fits_2025");
  TDirectory *directory2026 = output->mkdir("fits_2026");
  TDirectory *directoryMC = output->mkdir("momentum_scale_MC");
  TDirectory *directoryX2025 = output->mkdir("x_resolution_2025");
  TDirectory *directoryX2026 = output->mkdir("x_resolution_2026");

  output->cd();
  TH2 *input2025 = dynamic_cast<TH2 *>(h2_2025->Clone("hPhiMassVsPt_2025"));
  TH2 *input2026 = dynamic_cast<TH2 *>(h2_2026->Clone("hPhiMassVsPt_2026"));
  input2025->Write();
  input2026->Write();

  const std::string pdf2025 = outputDirectory + "/01_phi_fits_2025.pdf";
  const std::string pdf2026 = outputDirectory + "/02_phi_fits_2026.pdf";
  {
    TCanvas opener("opener2025", "", 10, 10);
    opener.Print((pdf2025 + "[").c_str());
  }
  {
    TCanvas opener("opener2026", "", 10, 10);
    opener.Print((pdf2026 + "[").c_str());
  }

  const std::size_t numberOfPtBins = ptEdges.size() - 1;
  std::vector<FitPoint> fits2025(numberOfPtBins);
  std::vector<FitPoint> fits2026(numberOfPtBins);
  for (std::size_t i = 0; i < numberOfPtBins; ++i)
  {
    TH1D *mass2025 = ProjectMass(h2_2025, ptEdges[i], ptEdges[i + 1],
                                 Form("hMass_2025_pt%02zu", i));
    TH1D *mass2026 = ProjectMass(h2_2026, ptEdges[i], ptEdges[i + 1],
                                 Form("hMass_2026_pt%02zu", i));
    fits2025[i] = FitOnePtBin(mass2025, "2025", static_cast<int>(i),
                              ptEdges[i], ptEdges[i + 1], fitLow, fitHigh,
                              directory2025, pdf2025);
    fits2026[i] = FitOnePtBin(mass2026, "2026", static_cast<int>(i),
                              ptEdges[i], ptEdges[i + 1], fitLow, fitHigh,
                              directory2026, pdf2026);
    delete mass2025;
    delete mass2026;
  }
  {
    TCanvas closer("closer2025", "", 10, 10);
    closer.Print((pdf2025 + "]").c_str());
  }
  {
    TCanvas closer("closer2026", "", 10, 10);
    closer.Print((pdf2026 + "]").c_str());
  }

  std::unique_ptr<TH1D> phiPtForToy(BuildSignalWeightedPt(h2_2025, h2_2026));
  if (!phiPtForToy || phiPtForToy->Integral() <= 0.)
  {
    std::cerr << "ERROR: could not build the phi pT spectrum for the toy."
              << std::endl;
    return;
  }
  output->cd();
  phiPtForToy->Write();

  const std::string pdfMC = outputDirectory + "/07_mc_phi_mass_response.pdf";
  {
    TCanvas opener("openerMC", "", 10, 10);
    opener.Print((pdfMC + "[").c_str());
  }
  std::vector<ToyPoint> toyPoints(numberOfPtBins);
  for (std::size_t i = 0; i < numberOfPtBins; ++i)
  {
    toyPoints[i] = DetermineCorrection(phiPtForToy.get(), fits2025[i], fits2026[i],
                                       static_cast<int>(i),
                                       toyAcceptedEventsPerBin, yMax,
                                       kaonPtMin, kaonEtaMax, directoryMC,
                                       pdfMC);
    if (toyPoints[i].valid)
    {
      std::cout << std::fixed << std::setprecision(4)
                << "pT " << ptEdges[i] << "-" << ptEdges[i + 1]
                << " GeV/c: m25=" << 1000. * fits2025[i].mean
                << " MeV, m26=" << 1000. * fits2026[i].mean
                << " MeV, epsilon_corr="
                << 100. * toyPoints[i].epsilonCorrection << " +/- "
                << 100. * toyPoints[i].epsilonError << " %" << std::endl;
    }
  }
  {
    TCanvas closer("closerMC", "", 10, 10);
    closer.Print((pdfMC + "]").c_str());
  }

  std::unique_ptr<TGraphErrors> gMass2025(MakeGraph(
      "gPhiMass2025", "#phi mass peak versus p_{T};p_{T,#phi} (GeV/c);"
                      "m_{#phi} (GeV/c^{2})",
      fits2025, false, 1.0));
  std::unique_ptr<TGraphErrors> gMass2026(MakeGraph(
      "gPhiMass2026", "#phi mass peak versus p_{T};p_{T,#phi} (GeV/c);"
                      "m_{#phi} (GeV/c^{2})",
      fits2026, false, 1.0));
  std::unique_ptr<TGraphErrors> gSigma2025(MakeGraph(
      "gPhiSigma2025", "#phi Gaussian resolution versus p_{T};"
                       "p_{T,#phi} (GeV/c);#sigma_{#phi} (MeV/c^{2})",
      fits2025, true, 1000.));
  std::unique_ptr<TGraphErrors> gSigma2026(MakeGraph(
      "gPhiSigma2026", "#phi Gaussian resolution versus p_{T};"
                       "p_{T,#phi} (GeV/c);#sigma_{#phi} (MeV/c^{2})",
      fits2026, true, 1000.));
  StyleGraph(gMass2025.get(), kRed + 1, 20);
  StyleGraph(gMass2026.get(), kBlue + 1, 21);
  StyleGraph(gSigma2025.get(), kRed + 1, 20);
  StyleGraph(gSigma2026.get(), kBlue + 1, 21);

  TGraphErrors gDeltaMass;
  gDeltaMass.SetName("gPhiMassDifference2025Minus2026");
  gDeltaMass.SetTitle("#phi mass difference; p_{T,#phi} (GeV/c);"
                      "m_{#phi}^{2025}-m_{#phi}^{2026} (MeV/c^{2})");
  TGraphErrors gEpsilon;
  gEpsilon.SetName("gEpsilonCorrection2025To2026");
  gEpsilon.SetTitle("2025 kaon momentum correction; p_{T,#phi} (GeV/c);"
                    "#epsilon_{corr} (%)");
  TGraphErrors gMass2025Corrected;
  gMass2025Corrected.SetName("gPhiMass2025CorrectedPrediction");
  gMass2025Corrected.SetTitle("Predicted corrected 2025 #phi mass;"
                              "p_{T,#phi} (GeV/c);m_{#phi} (GeV/c^{2})");
  int massDifferenceIndex = 0;
  int epsilonIndex = 0;
  int correctedIndex = 0;
  for (std::size_t i = 0; i < numberOfPtBins; ++i)
  {
    if (!fits2025[i].valid || !fits2026[i].valid)
      continue;
    const double x = 0.5 * (ptEdges[i] + ptEdges[i + 1]);
    const double ex = 0.5 * (ptEdges[i + 1] - ptEdges[i]);
    gDeltaMass.SetPoint(massDifferenceIndex, x,
                        1000. * (fits2025[i].mean - fits2026[i].mean));
    gDeltaMass.SetPointError(massDifferenceIndex, ex,
                             1000. * std::hypot(fits2025[i].meanError,
                                                fits2026[i].meanError));
    ++massDifferenceIndex;
    if (toyPoints[i].valid)
    {
      gEpsilon.SetPoint(epsilonIndex, x, 100. * toyPoints[i].epsilonCorrection);
      gEpsilon.SetPointError(epsilonIndex, ex, 100. * toyPoints[i].epsilonError);
      ++epsilonIndex;
      const double correctedMass = fits2025[i].mean +
                                   toyPoints[i].reproducedShift;
      const double correctedError = std::sqrt(
          fits2025[i].meanError * fits2025[i].meanError +
          toyPoints[i].reproducedShiftError * toyPoints[i].reproducedShiftError);
      gMass2025Corrected.SetPoint(correctedIndex, x, correctedMass);
      gMass2025Corrected.SetPointError(correctedIndex, ex, correctedError);
      ++correctedIndex;
    }
  }
  StyleGraph(&gDeltaMass, kBlack, 20);
  StyleGraph(&gEpsilon, kMagenta + 2, 20);
  StyleGraph(&gMass2025Corrected, kMagenta + 2, 24);

  // Parametrize the correction points with the same function implemented in
  // processopti5.  The graph y values and p0, p1, p4 and p5 are in percent.
  TF1 correctionParameterization(
      "fMomentumCorrectionSigmoidQuadratic",
      "[0]+[1]/(1.0+exp(-(x-[2])/[3]))+[4]*x+[5]*x*x",
      ptEdges.front(), ptEdges.back());
  correctionParameterization.SetParNames("p0_percent", "p1_percent",
                                         "p2_GeV", "p3_GeV",
                                         "p4_percent_per_GeV",
                                         "p5_percent_per_GeV2");
  correctionParameterization.SetParameters(
      0.749950, 0.421824, 1.614360, 0.290147, 0.0100, 0.0015);
  correctionParameterization.SetParLimits(0, 0.0, 2.0);
  correctionParameterization.SetParLimits(1, 0.0, 2.0);
  correctionParameterization.SetParLimits(2, 0.8, 3.0);
  correctionParameterization.SetParLimits(3, 0.05, 1.5);
  correctionParameterization.SetParLimits(4, -0.20, 0.20);
  correctionParameterization.SetParLimits(5, 0.0, 0.02);
  correctionParameterization.SetLineColor(kBlue + 1);
  correctionParameterization.SetLineWidth(3);
  correctionParameterization.SetNpx(10000);

  TFitResultPtr correctionFitResult(-1);
  int correctionFitStatus = -999;
  int correctionFitCovQual = -1;
  bool correctionFitValid = false;
  bool correctionFitDrawable = false;
  double correctionAICc = std::numeric_limits<double>::quiet_NaN();
  if (gEpsilon.GetN() >= 8)
  {
    TVirtualFitter::SetDefaultFitter("Minuit2");
    TVirtualFitter::SetMaxIterations(100000);
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
    // The x errors are pT-bin widths, not measurement uncertainties. Fit an
    // otherwise identical graph with ex=0 so they do not enter the chi2.
    TGraphErrors correctionFitGraph;
    for (int point = 0; point < gEpsilon.GetN(); ++point)
    {
      double x = 0.;
      double y = 0.;
      gEpsilon.GetPoint(point, x, y);
      correctionFitGraph.SetPoint(point, x, y);
      correctionFitGraph.SetPointError(point, 0., gEpsilon.GetErrorY(point));
    }
    // First Migrad pass moves away from the seeds. The second pass stores the
    // result and covariance. Do not request Minos/Improve here: those options
    // produced the misleading encoded status 4000 although covQual was 3.
    correctionFitGraph.Fit(&correctionParameterization, "Q0RN");
    correctionFitResult = correctionFitGraph.Fit(
        &correctionParameterization, "Q0RSN");
    correctionFitStatus = static_cast<int>(correctionFitResult);
    correctionFitCovQual = correctionFitResult.Get()
                               ? correctionFitResult->CovMatrixStatus()
                               : -1;
    correctionFitDrawable = correctionFitResult.Get();
    for (int parameter = 0; parameter < 6; ++parameter)
      correctionFitDrawable = correctionFitDrawable && std::isfinite(
                                                           correctionParameterization.GetParameter(parameter));
    correctionFitValid = correctionFitResult.Get() &&
                         correctionFitStatus == 0 && correctionFitCovQual >= 2;
    const int numberOfParameters = 6;
    const int numberOfPoints = gEpsilon.GetN();
    if (numberOfPoints > numberOfParameters + 1)
    {
      correctionAICc = correctionParameterization.GetChisquare() +
                       2. * numberOfParameters +
                       2. * numberOfParameters * (numberOfParameters + 1.) /
                           (numberOfPoints - numberOfParameters - 1.);
    }
  }
  std::cout << "\nMomentum-correction parameterization: status/covQual = "
            << correctionFitStatus << "/" << correctionFitCovQual
            << ", chi2/ndf = " << correctionParameterization.GetChisquare()
            << "/" << correctionParameterization.GetNDF()
            << ", AICc = " << correctionAICc << '\n';
  for (int parameter = 0; parameter < 6; ++parameter)
  {
    std::cout << "  p" << parameter << " = "
              << correctionParameterization.GetParameter(parameter)
              << " +/- " << correctionParameterization.GetParError(parameter)
              << '\n';
  }

  TCanvas massCanvas("cPhiMassVsPt", "phi mass", 950, 750);
  massCanvas.SetLeftMargin(0.13);
  gMass2025->SetMinimum(1.01880);
  gMass2025->SetMaximum(1.02);
  gMass2025->Draw("AP");
  gMass2025->GetXaxis()->SetLimits(ptEdges.front(), ptEdges.back());
  gMass2026->Draw("P SAME");
  TLine pdgLine(ptEdges.front(), kMassPhiPDG, ptEdges.back(), kMassPhiPDG);
  pdgLine.SetLineStyle(2);
  pdgLine.SetLineColor(kGray + 2);
  pdgLine.Draw();
  TLegend massLegend(0.67, 0.73, 0.88, 0.88);
  massLegend.SetBorderSize(0);
  massLegend.AddEntry(gMass2025.get(), "2025", "lep");
  massLegend.AddEntry(gMass2026.get(), "2026 reference", "lep");
  massLegend.AddEntry(&pdgLine, "PDG mass", "l");
  massLegend.Draw();
  massCanvas.Print((outputDirectory + "/03_phi_mass_vs_pt.pdf").c_str());

  TCanvas sigmaCanvas("cPhiSigmaVsPt", "phi resolution", 950, 750);
  sigmaCanvas.SetLeftMargin(0.13);
  gSigma2025->Draw("AP");
  gSigma2025->GetXaxis()->SetLimits(ptEdges.front(), ptEdges.back());
  gSigma2026->Draw("P SAME");
  TLegend sigmaLegend(0.67, 0.76, 0.88, 0.88);
  sigmaLegend.SetBorderSize(0);
  sigmaLegend.AddEntry(gSigma2025.get(), "2025", "lep");
  sigmaLegend.AddEntry(gSigma2026.get(), "2026", "lep");
  sigmaLegend.Draw();
  sigmaCanvas.Print((outputDirectory + "/04_phi_resolution_vs_pt.pdf").c_str());

  TCanvas differenceCanvas("cPhiDifferenceVsPt", "phi difference", 950, 750);
  differenceCanvas.SetLeftMargin(0.13);
  gDeltaMass.Draw("AP");
  gDeltaMass.GetXaxis()->SetLimits(ptEdges.front(), ptEdges.back());
  TLine zeroMass(ptEdges.front(), 0., ptEdges.back(), 0.);
  zeroMass.SetLineStyle(2);
  zeroMass.Draw();
  differenceCanvas.Print((outputDirectory + "/05_phi_mass_difference_vs_pt.pdf").c_str());

  TCanvas epsilonCanvas("cEpsilonVsPt", "momentum correction", 950, 750);
  epsilonCanvas.SetLeftMargin(0.13);
  gEpsilon.SetMinimum(0.6);
  gEpsilon.SetMaximum(2.4);
  gEpsilon.Draw("AP");
  gEpsilon.GetXaxis()->SetLimits(ptEdges.front(), ptEdges.back());
  if (correctionFitDrawable)
    correctionParameterization.Draw("L SAME");
  TLine zeroEpsilon(ptEdges.front(), 0., ptEdges.back(), 0.);
  zeroEpsilon.SetLineStyle(2);
  zeroEpsilon.Draw();
  TLatex convention;
  convention.SetNDC();
  convention.SetTextSize(0.035);
  convention.DrawLatex(0.16, 0.88,
                       "#vec{p}_{K}^{corr}=(1+#epsilon_{corr})#vec{p}_{K}^{2025}");
  if (correctionFitDrawable)
  {
    convention.SetTextColor(kBlue + 1);
    convention.SetTextSize(0.030);
    convention.DrawLatex(0.16, 0.83, "sigmoid + quadratic tail");
    convention.DrawLatex(
        0.16, 0.78,
        Form("#chi^{2}/ndf=%.1f/%d, AICc=%.1f",
             correctionParameterization.GetChisquare(),
             correctionParameterization.GetNDF(), correctionAICc));
    convention.DrawLatex(
        0.16, 0.73,
        Form("status/covQual=%d/%d", correctionFitStatus,
             correctionFitCovQual));
    convention.SetTextColor(kBlack);
  }
  epsilonCanvas.Print((outputDirectory + "/06_momentum_correction_vs_pt.pdf").c_str());

  TCanvas closureCanvas("cCorrectedPhiMassVsPt", "correction closure", 950, 750);
  closureCanvas.SetLeftMargin(0.13);
  gMass2026->Draw("AP");
  gMass2026->GetXaxis()->SetLimits(ptEdges.front(), ptEdges.back());
  gMass2025Corrected.Draw("P SAME");
  TLegend closureLegend(0.58, 0.75, 0.88, 0.88);
  closureLegend.SetBorderSize(0);
  closureLegend.AddEntry(gMass2026.get(), "2026 reference", "lep");
  closureLegend.AddEntry(&gMass2025Corrected, "2025 corrected (MC prediction)", "lep");
  closureLegend.Draw();
  closureCanvas.Print((outputDirectory + "/08_predicted_phi_mass_closure.pdf").c_str());

  XResolutionResult xResolution2025;
  XResolutionResult xResolution2026;
  if (runXResolutionToy)
  {
    if (xShapeMassHigh <= xShapeMassLow || xPtMin < 0. ||
        xMass <= 2. * kMassPhiPDG || xIntrinsicWidth <= 0. ||
        xDeltaMassMax <= 0. ||
        phiResolutionToyEventsPerPoint <= 0 ||
        xResolutionAcceptedEvents <= 0 || xRapidityGenerationMax <= 0.)
    {
      std::cerr << "ERROR: invalid X-resolution toy configuration; "
                   "skipping X-resolution estimation."
                << std::endl;
    }
    else
    {
      std::cout << "\nEstimating X detector resolution. The pT(X) inputs "
                   "are reconstructed, already-selected THnSparse shapes."
                << std::endl;
      xResolution2025 = EstimateXResolution(
          "2025", file2025.get(), xSparsePathOrName,
          xSparseMassAxis, xSparsePtAxis, xShapeMassLow, xShapeMassHigh,
          xPtMin, xMass, xIntrinsicWidth, xDeltaMassMax,
          fitLow, fitHigh, fits2025,
          phiResolutionToyEventsPerPoint,
          xResolutionAcceptedEvents, ptEdges.front(), ptEdges.back(), yMax,
          xRapidityGenerationMax, kaonPtMin, kaonEtaMax, xResolutionSeed,
          directoryX2025, outputDirectory + "/09_x_resolution_2025.pdf");
      xResolution2026 = EstimateXResolution(
          "2026", file2026.get(), xSparsePathOrName,
          xSparseMassAxis, xSparsePtAxis, xShapeMassLow, xShapeMassHigh,
          xPtMin, xMass, xIntrinsicWidth, xDeltaMassMax,
          fitLow, fitHigh, fits2026,
          phiResolutionToyEventsPerPoint,
          xResolutionAcceptedEvents, ptEdges.front(), ptEdges.back(), yMax,
          xRapidityGenerationMax, kaonPtMin, kaonEtaMax,
          xResolutionSeed + 900001U, directoryX2026,
          outputDirectory + "/10_x_resolution_2026.pdf");
      std::cout << "  sigma_X(2025) = " << xResolution2025.sigmaMeV
                << " +/- " << xResolution2025.sigmaErrorMeV
                << " MeV/c^2\n"
                << "  sigma_X(2026) = " << xResolution2026.sigmaMeV
                << " +/- " << xResolution2026.sigmaErrorMeV
                << " MeV/c^2\n"
                << "  observed Voigt FWHM(2025) = "
                << xResolution2025.observedFwhmMeV << " MeV/c^2\n"
                << "  observed Voigt FWHM(2026) = "
                << xResolution2026.observedFwhmMeV << " MeV/c^2"
                << std::endl;
    }
  }

  output->cd();
  gMass2025->Write();
  gMass2026->Write();
  gSigma2025->Write();
  gSigma2026->Write();
  gDeltaMass.Write();
  gEpsilon.Write();
  correctionParameterization.Write();
  if (correctionFitResult.Get())
    correctionFitResult->Write("fitResultMomentumCorrectionSigmoidQuadratic");
  gMass2025Corrected.Write();
  massCanvas.Write();
  sigmaCanvas.Write();
  differenceCanvas.Write();
  epsilonCanvas.Write();
  closureCanvas.Write();

  const std::string csvPath = outputDirectory + "/phi_momentum_scale_vs_pt.csv";
  std::ofstream csv(csvPath);
  csv << "pt_low_GeV,pt_high_GeV,fit2025_valid,mphi2025_GeV,mphi2025_err_GeV,"
         "sigma2025_GeV,sigma2025_err_GeV,fit2026_valid,mphi2026_GeV,"
         "mphi2026_err_GeV,sigma2026_GeV,sigma2026_err_GeV,"
         "delta_m_2025_minus_2026_GeV,mc_slope_GeV_per_epsilon,"
         "epsilon_correction_2025_to_2026,epsilon_error,closure_shift_GeV\n";
  csv << std::setprecision(12);
  for (std::size_t i = 0; i < numberOfPtBins; ++i)
  {
    csv << ptEdges[i] << ',' << ptEdges[i + 1] << ','
        << fits2025[i].valid << ',' << fits2025[i].mean << ','
        << fits2025[i].meanError << ',' << fits2025[i].sigma << ','
        << fits2025[i].sigmaError << ',' << fits2026[i].valid << ','
        << fits2026[i].mean << ',' << fits2026[i].meanError << ','
        << fits2026[i].sigma << ',' << fits2026[i].sigmaError << ','
        << (fits2025[i].mean - fits2026[i].mean) << ','
        << toyPoints[i].slope << ',' << toyPoints[i].epsilonCorrection << ','
        << toyPoints[i].epsilonError << ',' << toyPoints[i].reproducedShift << '\n';
  }
  csv.close();

  const std::string configPath = outputDirectory +
                                 "/momentum_correction_config.json";
  std::ofstream config(configPath);
  config << std::setprecision(10);
  config << "{\n  \"doublephimeson\": {\n";
  config << "    \"cfgApplyKaonMomentumCorrection\": "
         << (correctionFitValid ? "true" : "false") << ",\n";
  config << "    \"cfgMomCorrP0Percent\": "
         << correctionParameterization.GetParameter(0) << ",\n";
  config << "    \"cfgMomCorrP1Percent\": "
         << correctionParameterization.GetParameter(1) << ",\n";
  config << "    \"cfgMomCorrP2GeV\": "
         << correctionParameterization.GetParameter(2) << ",\n";
  config << "    \"cfgMomCorrP3GeV\": "
         << correctionParameterization.GetParameter(3) << ",\n";
  config << "    \"cfgMomCorrP4PercentPerGeV\": "
         << correctionParameterization.GetParameter(4) << ",\n";
  config << "    \"cfgMomCorrP5PercentPerGeV2\": "
         << correctionParameterization.GetParameter(5) << ",\n";
  config << "    \"cfgMomCorrPtMin\": " << ptEdges.front() << ",\n";
  config << "    \"cfgMomCorrPtMax\": " << ptEdges.back() << "\n";
  config << "  }\n}\n";
  config.close();

  const std::string summaryPath = outputDirectory + "/summary.txt";
  std::ofstream summary(summaryPath);
  summary << "Phi momentum-scale calibration: 2025 -> 2026\n\n";
  summary << "2025 file: " << file2025Name << "\n";
  summary << "2026 reference file: " << file2026Name << "\n";
  summary << "Histogram: " << histogramPathOrName << " (x=mass, y=pT)\n";
  summary << "Fit: binWidth*[0]*TMath::Voigt(x-[2],[3],[1],4) + pol2, Minuit2/Migrad\n";
  summary << "Invariant-mass projection: native binning retained (no Rebin)\n";
  summary << "Fixed Lorentzian width: " << 1000. * kWidthPhi << " MeV/c^2\n";
  summary << "Fit range: " << fitLow << "-" << fitHigh << " GeV/c^2\n";
  summary << "Toy: isotropic phi->K+K-, flat |y_phi|<" << yMax
          << ", kaon pT>" << kaonPtMin << " GeV/c, |eta_K|<" << kaonEtaMax
          << "\n";
  summary << "Toy events: " << toyAcceptedEventsPerBin << " accepted per pT bin\n\n";
  summary << "Definition: pK_corr = (1 + epsilon_corr) pK_2025.\n";
  summary << "A positive epsilon_corr increases the 2025 kaon momenta.\n\n";
  summary << "Correction parameterization:\n";
  summary << "epsilon_corr(%) = p0 + p1/[1+exp(-(pT_phi-p2)/p3)] "
             "+ p4*pT_phi + p5*pT_phi^2\n";
  summary << "status=" << correctionFitStatus
          << ", covQual=" << correctionFitCovQual
          << ", valid=" << correctionFitValid << "\n";
  summary << "chi2/ndf=" << correctionParameterization.GetChisquare()
          << "/" << correctionParameterization.GetNDF()
          << ", AICc=" << correctionAICc << "\n";
  for (int parameter = 0; parameter < 6; ++parameter)
  {
    summary << "p" << parameter << "="
            << correctionParameterization.GetParameter(parameter)
            << "+/-" << correctionParameterization.GetParError(parameter)
            << "\n";
  }
  summary << "\n";
  for (std::size_t i = 0; i < numberOfPtBins; ++i)
  {
    summary << "pT " << ptEdges[i] << "-" << ptEdges[i + 1] << " GeV/c: ";
    if (!toyPoints[i].valid)
    {
      summary << "INVALID\n";
      continue;
    }
    summary << "m25=" << 1000. * fits2025[i].mean << "+/-"
            << 1000. * fits2025[i].meanError << " MeV, m26="
            << 1000. * fits2026[i].mean << "+/-"
            << 1000. * fits2026[i].meanError << " MeV, sigma25="
            << 1000. * fits2025[i].sigma << "+/-"
            << 1000. * fits2025[i].sigmaError << " MeV, sigma26="
            << 1000. * fits2026[i].sigma << "+/-"
            << 1000. * fits2026[i].sigmaError << " MeV, epsilon_corr="
            << 100. * toyPoints[i].epsilonCorrection << "+/-"
            << 100. * toyPoints[i].epsilonError << "%\n";
  }
  summary << "\nX detector-resolution propagation:\n";
  summary << "enabled=" << runXResolutionToy << "\n";
  if (runXResolutionToy)
  {
    summary << "THnSparse=" << xSparsePathOrName
            << ", massAxis=" << xSparseMassAxis
            << ", ptAxis=" << xSparsePtAxis << "\n";
    summary << "Selected reconstructed shape: " << xShapeMassLow
            << " < M(phi phi) < " << xShapeMassHigh
            << " GeV/c^2, pT(X) > " << xPtMin << " GeV/c\n";
    summary << "The THnSparse pT(X) projection is already reconstructed and "
               "fully selected. It is not multiplied by an additional "
               "acceptance. For every sampled pT(X), the phase-space decay "
               "is generated conditionally inside the daughter acceptance, "
               "so the input pT shape is retained.\n";
    summary << "Reconstructed candidate cut before filling X response: "
            << "DeltaM_phi=sqrt[(m_phi1-m_PDG)^2+(m_phi2-m_PDG)^2] < "
            << xDeltaMassMax << " GeV/c^2. The pT(X) draw is held fixed "
                                "while regenerating until this cut is satisfied.\n";
    summary << "X pole mass=" << xMass
            << " GeV/c^2, generated with Breit-Wigner Gamma_X="
            << 1000. * xIntrinsicWidth << " MeV/c^2; Gamma_phi="
            << 1000. * kWidthPhi << " MeV/c^2.\n";
    summary << "For each pT(phi) point, the toy generates a Breit-Wigner phi "
               "with this same Gamma_phi, smears the two kaon momenta, fills "
               "M(KK)_reco, and fits it with TMath::Voigt using Gamma_phi "
               "fixed. The inferred sigma(pK)/pK is accepted only through "
               "this Voigt-sigma closure; Gamma_phi is therefore not counted "
               "as detector resolution.\n";
    summary << "2025: valid=" << xResolution2025.valid
            << ", accepted=" << xResolution2025.accepted
            << ", status/covQual=" << xResolution2025.fitStatus << "/"
            << xResolution2025.covQual
            << ", mean=" << xResolution2025.meanMeV << "+/-"
            << xResolution2025.meanErrorMeV
            << " MeV/c^2, sigma=" << xResolution2025.sigmaMeV << "+/-"
            << xResolution2025.sigmaErrorMeV
            << " MeV/c^2, RMS=" << xResolution2025.rmsMeV
            << " MeV/c^2, observedMean="
            << xResolution2025.observedMeanMeV << "+/-"
            << xResolution2025.observedMeanErrorMeV
            << " MeV/c^2, observedVoigtSigma="
            << xResolution2025.observedSigmaMeV << "+/-"
            << xResolution2025.observedSigmaErrorMeV
            << " MeV/c^2, observedVoigtFWHM="
            << xResolution2025.observedFwhmMeV << " MeV/c^2\n";
    summary << "2026: valid=" << xResolution2026.valid
            << ", accepted=" << xResolution2026.accepted
            << ", status/covQual=" << xResolution2026.fitStatus << "/"
            << xResolution2026.covQual
            << ", mean=" << xResolution2026.meanMeV << "+/-"
            << xResolution2026.meanErrorMeV
            << " MeV/c^2, sigma=" << xResolution2026.sigmaMeV << "+/-"
            << xResolution2026.sigmaErrorMeV
            << " MeV/c^2, RMS=" << xResolution2026.rmsMeV
            << " MeV/c^2, observedMean="
            << xResolution2026.observedMeanMeV << "+/-"
            << xResolution2026.observedMeanErrorMeV
            << " MeV/c^2, observedVoigtSigma="
            << xResolution2026.observedSigmaMeV << "+/-"
            << xResolution2026.observedSigmaErrorMeV
            << " MeV/c^2, observedVoigtFWHM="
            << xResolution2026.observedFwhmMeV << " MeV/c^2\n";
  }
  summary.close();

  WriteCorrectionHeader(outputDirectory + "/PhiMomentumScale2025To2026.h",
                        ptEdges, toyPoints);
  output->Write();
  output->Close();

  std::cout << "\nFinished. Main outputs:\n"
            << "  " << outputRoot << "\n"
            << "  " << csvPath << "\n"
            << "  " << configPath << "\n"
            << "  " << outputDirectory + "/PhiMomentumScale2025To2026.h" << "\n"
            << "  " << outputDirectory + "/09_x_resolution_2025.pdf" << "\n"
            << "  " << outputDirectory + "/10_x_resolution_2026.pdf" << "\n"
            << "  " << summaryPath << std::endl;
}
