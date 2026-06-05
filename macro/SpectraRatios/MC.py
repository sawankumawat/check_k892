#!/usr/bin/env python3


from ROOT import TFile, TGraph, TGraphErrors
from common import *


def target_mult(multiplicity):
    target_multiplicity = {"I": [35.8, 0.5],
                           "II": [32.2, 0.4],
                           "III": [30.1, 0.4]}
    return target_multiplicity[multiplicity]


def normalization_h(hv0, hv0nch, multiplicity):
    gmean = TGraph()
    bins_of_interest = []
    for i in range(1, hv0nch.GetNbinsX()+1):
        if i > 200:
            continue
        hpp = hv0nch.ProjectionY("hpp", i, i)
        mean = hpp.GetMean()
        del hpp
        gmean.AddPoint(hv0nch.GetXaxis().GetBinCenter(i), mean)
        if mean > multiplicity[0]-multiplicity[1] and mean < multiplicity[0]+multiplicity[1]:
            bins_of_interest.append(i)
    v0mrange = [bins_of_interest[0], bins_of_interest[-1]]
    if 0:
        draw_nice_canvas(hv0nch.GetName())
        hv0nch.Draw("COL")
        print("Bins of interest", bins_of_interest)
        gmin = TGraph()
        gmin.AddPoint(0, multiplicity[0]-multiplicity[1])
        gmin.AddPoint(300, multiplicity[0]-multiplicity[1])
        gmax = TGraph()
        gmax.AddPoint(0, multiplicity[0]+multiplicity[1])
        gmax.AddPoint(300, multiplicity[0]+multiplicity[1])
        gmin.Draw("L")
        gmax.Draw("L")
        gmean.Draw("PL")

        draw_nice_canvas("AverageNch")
        hv0nc_p = hv0nch.ProjectionY("hv0nc_p", *v0mrange)
        hv0nc_p.Draw()
        update_all_canvases()

        draw_nice_canvas("Normalization")
        hv0.Draw()
        gmax_norm = TGraph()
        gmax_norm.AddPoint(v0mrange[1], 0)
        gmax_norm.AddPoint(v0mrange[1], hv0.GetMaximum())
        gmax_norm.Draw("PL")
        gmin_norm = TGraph()
        gmin_norm.AddPoint(v0mrange[0], 0)
        gmin_norm.AddPoint(v0mrange[0], hv0.GetMaximum())
        gmin_norm.Draw("PL")

        input("Press enter to continue")
        print(norm, v0mrange)
    norm = hv0.Integral(v0mrange[0], v0mrange[1])
    return norm, v0mrange


def normalization(d, multiplicity):
    # d.ls()
    hv0 = d.FindObject("fHistV0MMult")
    hv0nch = d.FindObject("fHistNchVsV0MMult")
    return normalization_h(hv0, hv0nch, multiplicity)


def get_spectra_predictions(particle="XiMinus",
                            multiplicity="I",
                            monash=True,
                            draw=False):
    multiplicity = target_mult(multiplicity)

    if monash:
        f = "MC/TrainOutput/AnalysisResults_3065_20240416-1524.root"
        f = "MC/TrainOutput/MergedMonash.root"
    else:
        f = "MC/TrainOutput/AnalysisResults_3066_20240416-1524.root"
        f = "MC/TrainOutput/MergedRopes.root"
    f = TFile(f, "READ")
    # f.ls()
    d = f.Get("PWGLF_MCPredictions/cListminbias")
    # d.ls()
    number_of_events, v0mrange = normalization(d, multiplicity)     #it will return the number of events in I multiplicity class and v0mrange means the starting and the end bin of the multiplicity for v0M 
    h = "fHistDDPt"
    h = "fHistDDYield"
    h = "fHistPtVsV0MMult"
    h = d.FindObject(h+"_"+particle)
    h = h.ProjectionY(h.GetName()+f"_{v0mrange}", *v0mrange)
    h.Scale(1/number_of_events, "width")
    # h.Scale(1, "width")
    if draw:
        if monash:
            draw_nice_canvas("Monash")
        else:
            draw_nice_canvas("Ropes")
        h.Draw()
        update_all_canvases()
    h.SetDirectory(0)
    f.Close()
    return h


def make_yields(particle="XiMinus",
                monash=True,
                draw=False):
    if monash:
        f = "MC/TrainOutput/AnalysisResults_3065_20240416-1524.root"
        f = "MC/TrainOutput/MergedMonash.root"
    else:
        f = "MC/TrainOutput/AnalysisResults_3066_20240416-1524.root"
        f = "MC/TrainOutput/MergedRopes.root"
    f = TFile(f, "READ")
    # f.ls()
    d = f.Get("PWGLF_MCPredictions/cListminbias")
    # d.ls()
    h = "fHistPtVsV0MMult"
    h = d.FindObject(h+"_"+particle)
    hp = h.ProjectionX(h.GetName()+"_yield")
    hv0nch = d.FindObject("fHistNchVsV0MMult")
    hv0 = d.FindObject("fHistV0MMult")
    prediction = TGraphErrors()
    prediction.SetName(("Monash" if monash else "Ropes") + particle)
    for i in range(1, hp.GetNbinsX()+1):
        x = hp.GetXaxis().GetBinCenter(i)
        if x > 150:
            continue
        xi = hv0nch.GetXaxis().FindBin(x)
        hpp = hv0nch.ProjectionY("hpp", xi, xi)
        y = hpp.GetMean()
        prediction.AddPoint(i,y, hp.GetBinContent(i)/hv0.GetBinContent(i))
#        prediction.SetPointError(prediction.GetN()-1, hpp.GetRMS(), hp.GetBinError(i)/hv0.GetBinContent(i))
    if draw:
        draw_nice_canvas("Corr")
        h.Draw("COL")
        draw_nice_canvas("Yield")
        hp.Draw()
        draw_nice_canvas("Yieldvsmult")
        prediction.Draw("ALP")
        update_all_canvases()
        input("Press enter to continue")
    return prediction


def make_meanpt(particle="XiMinus",
                monash=True,
                draw=False):
    if monash:
        f = "MC/TrainOutput/AnalysisResults_3065_20240416-1524.root"
        f = "MC/TrainOutput/MergedMonash.root"
    else:
        f = "MC/TrainOutput/AnalysisResults_3066_20240416-1524.root"
        f = "MC/TrainOutput/MergedRopes.root"
    f = TFile(f, "READ")
    f.ls()
    d = f.Get("PWGLF_MCPredictions/cListminbias")
    # d.ls()
    h = "fHistPtVsV0MMult"
    h = d.FindObject(h+"_"+particle)
    if 0:
        draw_nice_canvas("fHistPtVsV0MMult", replace=False)
        h.Draw("COL")
        update_all_canvases()
        input("Press enter to continue")

    hv0nch = d.FindObject("fHistNchVsV0MMult")
    prediction = TGraphErrors()
    prediction.SetName(("Monash" if monash else "Ropes") + particle)
    for i in range(1, h.GetXaxis().GetNbins()+1):
        x = h.GetXaxis().GetBinCenter(i)
        if x > 150:
            continue
        xi = hv0nch.GetXaxis().FindBin(x)
        hpp = hv0nch.ProjectionY("hpp", xi, xi)
        if 0:
            draw_nice_canvas("MeanPtExtraction", replace=False)
            hpp.Draw()
            update_all_canvases()
            input("Press enter to continue")

        y = hpp.GetMean()
        # print("Multiplicity at forward", x, "corresponds to multiplicity at midrapidity", y)
        meanpt = h.ProjectionY("extra_pt", i, i)
        if 0:
            draw_nice_canvas("MeanPtExtraction", replace=False)
            meanpt.Draw()
            update_all_canvases()
            input("Press enter to continue")
        prediction.AddPoint(y, meanpt.GetMean())
        prediction.SetPointError(prediction.GetN()-1, hpp.GetRMS(), meanpt.GetMeanError())
        # prediction.AddPoint(y, hp.GetBinContent(i)/hv0.GetBinContent(i))
        # prediction.SetPointError(prediction.GetN()-1, hpp.GetMeanError(), hp.GetBinError(i)/hv0.GetBinContent(i))
    if draw :
        draw_nice_canvas("Corr")
        h.Draw("COL")
        draw_nice_canvas("Ptvsmult")
        prediction.Draw("ALP")
        update_all_canvases()
        input("Press enter to continue")
    return prediction


def getpred(paritcle="XiMinus", multiplicity="I", monash=True):
    f = f"MC/{'Monash' if monash else 'Ropes'}_{paritcle}_{multiplicity}.root"
    f = f"MC/Merged{'Monash' if monash else 'Ropes'}_{paritcle}_{multiplicity}.root"
    f = TFile(f, "READ")
    h = f.Get(f.GetListOfKeys().At(0).GetName())
    h.SetDirectory(0)
    return h


def getpredyields(paritcle="XiMinus", monash=True):
    f = f"MC/{'Monash' if monash else 'Ropes'}_Yields_{paritcle}.root"
    f = f"MC/Merged{'Monash' if monash else 'Ropes'}_Yields_{paritcle}.root"
    f = TFile(f, "READ")
    h = f.Get(f.GetListOfKeys().At(0).GetName())
    # print(f.GetName(), "HEy", h)
    return h


def make_spectra(multiplicities=["III"], particles=["XiMinus", "XiPlus"]):
    for i in multiplicities:
        for j in particles:
            get_spectra_predictions(j, i, monash=True, draw=False).SaveAs(f"MC/MergedMonash_{j}_{i}.root")
            get_spectra_predictions(j, i, monash=False, draw=False).SaveAs(f"MC/MergedRopes_{j}_{i}.root")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--multiplicities", "-m", nargs=1, choices=["I", "II", "III"], default=None)
    parser.add_argument("--particles", "-p", nargs="+", choices=["XiMinus", "XiPlus",
                                                                 "OmegaMinus", "OmegaPlus",
                                                                 "PiMinus", "PiPlus",
                                                                 "KMinus", "KPlus",
                                                                 "AntiProton", "Proton"])
    parser.add_argument("--mode", "-M", default="pt", choices=["spectra", "yields", "pt"])
    args = parser.parse_args()
    if args.mode == "spectra":
        make_spectra(particles=args.particles, multiplicities=args.multiplicities)
    elif args.mode == "pt":
        for i in args.particles:
            make_meanpt(i, monash=True).SaveAs(f"MC/MergedMonash_MeanPt_{i}.root")
            make_meanpt(i, monash=False).SaveAs(f"MC/MergedRopes_MeanPt_{i}.root")
    elif args.mode == "yields":
        for i in args.particles:
            make_yields(i, monash=True).SaveAs(f"MC/MergedMonash_Yields_{i}.root")
            make_yields(i, monash=False).SaveAs(f"MC/MergedRopes_Yields_{i}.root")
    else:
        raise ValueError("Unknown mode")
