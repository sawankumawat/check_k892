import numpy as np
import ROOT
from pathlib import Path
## Configuration ##
kTestRun = 0
kCentTest = 0
kIsINEL = 0

kRootExtension = ".root"

kBaseInputDir = "../../"
# dataset = "LHC22cde" # or LHC22m
# dataset = "LHC22m" # or LHC22m
dataset = "LHC22o" # or LHC22m

kBaseOutputDir = ""
if dataset == "LHC22cde":
    kBaseOutputDir = "../../output/LHC22cde_pass4"
if dataset == "LHC22m":
    kBaseOutputDir = "../../output/LHC22m_pass4"
if dataset == "LHC22o":
    kBaseOutputDir = "../../output/LHC22o_pass4"

if kIsINEL > 0:
    kBaseOutputDir += "_INEL"

kBaseOutputDir += "/"

kFiguresFolder = kBaseOutputDir + "images/"
kOutputName = "lf-k892analysis"

kDataFilename = ""
if dataset == "LHC22cde":
    kDataFilename = kBaseInputDir + "data/LHC22cde_pass4/107888_89_90.root"
if dataset == "LHC22m":
    kDataFilename = kBaseInputDir + "data/LHC22m_pass4_small/107887.root"

kMCfilename = ""
if dataset == "LHC22cde":
    kMCfilename = kBaseInputDir +"mc/LHC22h1b2/AnalysisResults_gen900.root"
    kMCfilename = kBaseInputDir + "mc/LHC23f3/AnalysisResults_pp13.6TeV.root"

kInitOutput = kBaseOutputDir + "common/init" + kRootExtension
kSignalOutput = kBaseOutputDir + "common/signal" + kRootExtension
kEfficiencyOutput = kBaseOutputDir + "common/efficiency" + kRootExtension
kSpectraOutput = kBaseOutputDir + "common/spectra" + kRootExtension
kSystematicsOutput = kBaseOutputDir + "common/systematics" + kRootExtension

kInitOutputPart = kInitOutput[:-len(kRootExtension)]
kSignalOutputPart = kSignalOutput[:-len(kRootExtension)]
kEfficiencyOutputPart = kEfficiencyOutput[:-len(kRootExtension)]
kSpectraOutputPart = kSpectraOutput[:-len(kRootExtension)]
kYieldMeanOutputPart = kBaseOutputDir + "yieldmean"

kSystematicsDetailOutput = kBaseOutputDir + "systematic_details/systematics" + kRootExtension

kFilterListResoTaskNames = "lf-reso2initializer"
kFilterListNames = "lf-k892analysis"
kFilterListNamesMC = "lf-k892analysis"

masspdg = 0.895
widthpdg = 0.047
branchingratio = 0.66
kMaterialBudgetError = 0.01

kDoMCQA = 0
kParticleType = ["K892", "AntiK892"]

kAntiBinSize = 1
kRebin = 5

kMultiplicityBins = []
if dataset == "LHC22cde":
    kMultiplicityBins = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
if dataset == "LHC22m":
    kMultiplicityBins = [0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]

kDrawRange = [0.65, 1.4]
kDrawFitRange = [0.65, 1.3]

kpTbin = []
if dataset == "LHC22cde":
    kpTbin = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 5.0]
if dataset == "LHC22m":
    kpTbin =[0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0]
    
kLastpTbin = kpTbin[-1]

kNPtBins = len(kpTbin)
kNMultiplicityBins = len(kMultiplicityBins) - 1


## Drawing configuration
kCanvasW = 960
kCanvasH = 720
kDrawRangePlot = np.array([1.26,1.8])
kPtBinsPlot = np.array([0.7, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 3.2, 4.0, 4.8, 6])
kMaterialColors = [ROOT.kOrange + 10, ROOT.kOrange - 3, ROOT.kOrange, ROOT.kSpring - 6, ROOT.kTeal + 2, ROOT.kAzure + 7,  ROOT.kAzure - 1, ROOT.kViolet -6, ROOT.kViolet - 1,  ROOT.kPink + 8, ROOT.kMagenta + 3]

SystematicBins = [
    "", 
    "_id5195", "_id5196", "_id5197", "_id5198", "_id5199", "_id5200", "_id5201", "_id5202", "_id5203", "_id5204"]
SystematicBinsName = [
    "",    "", "_PVTrk", "_DCAxy_loose", "_DCAz_loose", "_DCA_loose", "_PID25", "_PID35", "_PID20", "_PID40", "_Trk_loose"]

# Plotter
kCanvasW = 720
kCanvasH = 720

def GetCanvas(cname, w, h):
    c = ROOT.TCanvas(cname,cname,w,h)
    c.SetTickx()
    c.SetTicky()
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.13)
    c.SetBottomMargin(0.12)
    c.SetRightMargin(0.01)
    c.SetFillStyle(0)
    return c

def SaveCanvas(c, name = "temp", path = "figs/", ftype = "pdf"):
    # Save canvas with path. if path is not there, make a folder.
    Path(path).mkdir(parents=True, exist_ok=True)
    c.SaveAs(path+name+"."+ftype)

# Default Latex
t0 = ROOT.TLatex()
t0R = ROOT.TLatex()
t00 = ROOT.TLatex()
t00R = ROOT.TLatex()
t = ROOT.TLatex()
tR = ROOT.TLatex()
t1 = ROOT.TLatex()
t1R = ROOT.TLatex()
t2 = ROOT.TLatex()
t2R = ROOT.TLatex()
t3 = ROOT.TLatex()
t3R = ROOT.TLatex()
t4 = ROOT.TLatex()
t4R = ROOT.TLatex()

# Common latex
# for memo, super small
t00.SetNDC()
t00.SetTextSize(0.025)
t00R.SetNDC()
t00R.SetTextSize(0.025)
t00R.SetTextAlign(33)
# for memo, super small
t0.SetNDC()
t0.SetTextSize(0.030)
t0R.SetNDC()
t0R.SetTextSize(0.030)
t0R.SetTextAlign(33)
# for memo, small
t.SetNDC()
t.SetTextSize(0.035)
tR.SetNDC()
tR.SetTextSize(0.035)
tR.SetTextAlign(33)
# for intermidiate
t1.SetNDC()
t1.SetTextSize(0.05)
t1R.SetNDC()
t1R.SetTextSize(0.05)
t1R.SetTextAlign(33)
# for warning, big
t2.SetNDC()
t2.SetTextSize(0.04)
t2R.SetNDC()
t2R.SetTextSize(0.04)
t2R.SetTextAlign(33)
# for small pad, huge
t3.SetNDC()
t3.SetTextSize(0.06)
t3R.SetNDC()
t3R.SetTextSize(0.06)
t3R.SetTextAlign(33)
# 
t4.SetNDC()
t4.SetTextSize(0.08)
t4R.SetNDC()
t4R.SetTextSize(0.08)
t4R.SetTextAlign(33)

if __name__ == "__main__":
    print("Common variables loaded")