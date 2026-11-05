#include <iostream>
#include "src/style.h"

void makeQAplotsData()
{
    gStyle->SetOptStat(0);
    string outputtype = "pdf";
    int colors[] = {kBlue + 2, kRed + 1, kGreen + 2, kMagenta + 2, kCyan + 2, kOrange + 7, kViolet + 3, kPink + 1, kAzure + 7, kTeal + 7, kBlue};

    TString data_path = "/home/sawan/Storage/check_k892/output/kstar/LHC22o_pass7/749276/kstarqa/hInvMass/ROTATED";
    TFile *fileeff = new TFile("/home/sawan/Storage/check_k892/data/kstar/LHC22o_pass7/794581.root", "READ");

    string QaPath = "kstarqa";

    TString output_QA_folder = data_path + "/QA";
    TString output_QA_folder_before = data_path + "/QA/before";
    gSystem->mkdir(output_QA_folder, kTRUE);
    gSystem->mkdir(output_QA_folder_before, kTRUE);

    TH3F *hDCAxy3D = (TH3F *)fileeff->Get(Form("%s/eventSelection/hDcaxy_cent_pt", QaPath.c_str()));
    TH3F *hDCAz3D = (TH3F *)fileeff->Get(Form("%s/eventSelection/hDcaz_cent_pt", QaPath.c_str()));
    TH1F *hMult = (TH1F *)fileeff->Get(Form("%s/eventSelection/hMultiplicity", QaPath.c_str()));
    TH1F *hvz = (TH1F *)fileeff->Get(Form("%s/eventSelection/hVertexZRec", QaPath.c_str()));
    TH1F *hoccupancy = (TH1F *)fileeff->Get(Form("%s/eventSelection/hOccupancy", QaPath.c_str()));
    TH1F *hEventCut = (TH1F *)fileeff->Get(Form("%s/eventSelection/hEventCut", QaPath.c_str()));
    TH1F *htracksRecMC = (TH1F *)fileeff->Get(Form("%s/eventSelection/recMCparticles", QaPath.c_str()));
    TH1D *hDCAxy = hDCAxy3D->ProjectionX("hDCAxy", -1, -1, -1, -1);
    TH1D *hDCAz = hDCAz3D->ProjectionX("hDCAz", -1, -1, -1, -1);
    if (hDCAxy == nullptr || hDCAz == nullptr || hvz == nullptr || hoccupancy == nullptr || hEventCut == nullptr || htracksRecMC == nullptr)
    {
        cerr << "Event selection histograms not found!!!!!!!!!!!!" << endl;
        return;
    }

    TH3F *hNsigmaTPCTOFKaon = (TH3F *)fileeff->Get(Form("%s/hPID/After/hNsigma_TPC_TOF_Ka_after", QaPath.c_str()));
    TH3F *hNsigmaTPCTOFPion = (TH3F *)fileeff->Get(Form("%s/hPID/After/hNsigma_TPC_TOF_Pi_after", QaPath.c_str()));
    TH3F *hNsigmaTPCKaon = (TH3F *)fileeff->Get(Form("%s/hPID/After/hTPCnsigKa_mult_pt", QaPath.c_str()));
    TH3F *hNsigmaTPCPion = (TH3F *)fileeff->Get(Form("%s/hPID/After/hTPCnsigPi_mult_pt", QaPath.c_str()));
    TH3F *hNsigmaTOFKaon = (TH3F *)fileeff->Get(Form("%s/hPID/After/hTOFnsigKa_mult_pt", QaPath.c_str()));
    TH3F *hNsigmaTOFPion = (TH3F *)fileeff->Get(Form("%s/hPID/After/hTOFnsigPi_mult_pt", QaPath.c_str()));
    if (hNsigmaTPCTOFKaon == nullptr || hNsigmaTPCTOFPion == nullptr || hNsigmaTPCKaon == nullptr || hNsigmaTPCPion == nullptr || hNsigmaTOFKaon == nullptr || hNsigmaTOFPion == nullptr)
    {
        cerr << "PID histograms after selection not found!!!!!!!!!!!!" << endl;
        return;
    }

    TH3F *hNsigmaTPCTOFKaon_before = (TH3F *)fileeff->Get(Form("%s/hPID/Before/hNsigma_TPC_TOF_Ka_before", QaPath.c_str()));
    TH3F *hNsigmaTPCTOFPion_before = (TH3F *)fileeff->Get(Form("%s/hPID/Before/hNsigma_TPC_TOF_Pi_before", QaPath.c_str()));
    TH3F *hNsigmaTPCKaon_before = (TH3F *)fileeff->Get(Form("%s/hPID/Before/hTPCnsigKa_mult_pt", QaPath.c_str()));
    TH3F *hNsigmaTPCPion_before = (TH3F *)fileeff->Get(Form("%s/hPID/Before/hTPCnsigPi_mult_pt", QaPath.c_str()));
    TH3F *hNsigmaTOFKaon_before = (TH3F *)fileeff->Get(Form("%s/hPID/Before/hTOFnsigKa_mult_pt", QaPath.c_str()));
    TH3F *hNsigmaTOFPion_before = (TH3F *)fileeff->Get(Form("%s/hPID/Before/hTOFnsigPi_mult_pt", QaPath.c_str()));
    if (hNsigmaTPCTOFKaon_before == nullptr || hNsigmaTPCTOFPion_before == nullptr || hNsigmaTPCKaon_before == nullptr || hNsigmaTPCPion_before == nullptr || hNsigmaTOFKaon_before == nullptr || hNsigmaTOFPion_before == nullptr)
    {
        cerr << "PID histograms before selection not found!!!!!!!!!!!!" << endl;
        return;
    }

    // outputtype = "png";
    TCanvas *cMult = new TCanvas("cMult", "Multiplicity", 720, 720);
    SetCanvasStyle(cMult, 0.14, 0.03, 0.06, 0.14);
    SetHistoQA(hMult);
    // hMult->GetXaxis()->SetTitle("Multiplicity (%)");
    hMult->GetXaxis()->SetTitle("Multiplicity (%)");
    hMult->GetYaxis()->SetTitle("Counts");
    hMult->GetXaxis()->SetRangeUser(0, 109.5);
    hMult->SetMaximum(0.13e9);
    hMult->Draw();
    cMult->SaveAs(output_QA_folder + ("/Multiplicity." + outputtype).c_str());

    TCanvas *cvz = new TCanvas("cvz", "Vertex Z", 720, 720);
    SetCanvasStyle(cvz, 0.14, 0.03, 0.06, 0.14);
    SetHistoQA(hvz);
    hvz->GetXaxis()->SetTitle("Vertex Z (cm)");
    hvz->GetYaxis()->SetTitle("Counts");
    hvz->Draw();
    cvz->SaveAs(output_QA_folder + ("/VertexZ." + outputtype).c_str());

    TCanvas *cTracksData = new TCanvas("cTracksData", "Tracks Data", 1440, 720);
    SetCanvasStyle(cTracksData, 0.1, 0.06, 0.06, 0.16);
    SetHistoQA(htracksRecMC);
    // gPad->SetLogy(1);
    gPad->SetGrid(1, 0);
    gPad->SetLogy(1);
    htracksRecMC->SetTitle("Tracks Rec MC");
    htracksRecMC->GetYaxis()->SetTitle("Counts");
    htracksRecMC->GetYaxis()->SetTitleOffset(0.8);
    htracksRecMC->GetXaxis()->SetRangeUser(0, 17);
    // htracksRecMC->SetMaximum(htracksRecMC->GetMaximum() * 10);
    htracksRecMC->Draw();
    cTracksData->SaveAs(output_QA_folder + ("/TracksRecMC." + outputtype).c_str());

    TCanvas *cEventCut = new TCanvas("cEventCut", "Event Cut", 720, 720);
    SetCanvasStyle(cEventCut, 0.13, 0.05, 0.06, 0.17);
    SetHistoQA(hEventCut);
    gPad->SetGrid(1, 0);
    // gPad->SetLogy(1);
    hEventCut->GetXaxis()->SetBinLabel(4, "INEL > 0");
    hEventCut->SetBinContent(4, hEventCut->GetBinContent(9));
    hEventCut->SetTitle("Event selections");
    hEventCut->GetYaxis()->SetTitle("Counts");
    hEventCut->GetXaxis()->SetRangeUser(0, 4);
    hEventCut->SetMinimum(0);
    hEventCut->SetMaximum(hEventCut->GetBinContent(1) * 1.2);
    hEventCut->Draw();
    double firstBin = hEventCut->GetBinContent(1);

    TLatex latex;
    latex.SetTextSize(0.03);
    latex.SetTextAlign(22); // centered

    for (int i = 1; i <= 4; i++)
    {

        double content = hEventCut->GetBinContent(i);

        double percent = 0;
        if (firstBin > 0)
            percent = (content / firstBin) * 100.0;

        TString label = Form("%.1f%%", percent);

        double x = hEventCut->GetBinCenter(i);
        double y = content * 1.05; // slightly above bin

        latex.DrawLatex(x, y, label);
    }
    cEventCut->SaveAs(output_QA_folder + ("/EventCut." + outputtype).c_str());

    //===============PID plots==========================
    // Plots to make. 2D plots TPCKa, TPCPi, TOFKa, TOFPi as a function of pT. Then 2D plot of TPC vs TOF of Ka and Pi at different pT ranges.

    double pTrangesForTPCandTOF[] = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0};
    vector<vector<int>> multRangesForTPCandTOF = {{0, 1}, {1, 5}, {10, 20}, {30, 40}, {50, 70}, {70, 100}};
    TH2F *h2DTPCKaon = (TH2F *)hNsigmaTPCKaon->Project3D("xz");
    TH2F *h2DTPCPion = (TH2F *)hNsigmaTPCPion->Project3D("xz");
    TH2F *h2DTOFKaon = (TH2F *)hNsigmaTOFKaon->Project3D("xz");
    TH2F *h2DTOFPion = (TH2F *)hNsigmaTOFPion->Project3D("xz");
    TH2 *h2DNsigmaTPCTOFKaon[6];
    TH2 *h2DNsigmaTPCTOFPion[6];
    TH1F *h1DNsigmaTPCKaon[6];
    TH1F *h1DNsigmaTPCPion[6];
    TH1F *h1DNsigmaTOFKaon[6];
    TH1F *h1DNsigmaTOFPion[6];

    TH1F *h1DNsigmaTPCKaon_before[6];
    TH1F *h1DNsigmaTPCPion_before[6];
    TH1F *h1DNsigmaTOFKaon_before[6];
    TH1F *h1DNsigmaTOFPion_before[6];

    TH1F *h1DNsigmaTPCKaon_pt[6];
    TH1F *h1DNsigmaTPCPion_pt[6];
    TH1F *h1DNsigmaTOFKaon_pt[6];
    TH1F *h1DNsigmaTOFPion_pt[6];

    for (int i = 0; i < 6; i++)
    {
        int lowBinpT_z = hNsigmaTPCTOFKaon->GetZaxis()->FindBin(pTrangesForTPCandTOF[i] + 0.001);
        int highBinpT_z = hNsigmaTPCTOFKaon->GetZaxis()->FindBin(pTrangesForTPCandTOF[i + 1] - 0.001);

        hNsigmaTPCTOFKaon->GetZaxis()->SetRange(lowBinpT_z, highBinpT_z);
        hNsigmaTPCTOFPion->GetZaxis()->SetRange(lowBinpT_z, highBinpT_z);

        h2DNsigmaTPCTOFKaon[i] = (TH2 *)hNsigmaTPCTOFKaon->Project3D("yx"); // the axis are switched during projection, so we get TPC on y and TOF on x
        h2DNsigmaTPCTOFPion[i] = (TH2 *)hNsigmaTPCTOFPion->Project3D("yx");

        h2DNsigmaTPCTOFKaon[i]->SetName(Form("xy_Ka_%d", i));
        h2DNsigmaTPCTOFPion[i]->SetName(Form("xy_Pi_%d", i));

        // Now lets make 1D projections with limit on multiplicity ranges
        int lowBinMult = hNsigmaTPCKaon->GetYaxis()->FindBin(multRangesForTPCandTOF[i][0] + 0.001);
        int highBinMult = hNsigmaTPCKaon->GetYaxis()->FindBin(multRangesForTPCandTOF[i][1] - 0.001);

        h1DNsigmaTPCKaon[i] = (TH1F *)hNsigmaTPCKaon->ProjectionX(Form("h1D_TPC_Ka_%d", i), lowBinMult, highBinMult, -1, -1);
        h1DNsigmaTPCPion[i] = (TH1F *)hNsigmaTPCPion->ProjectionX(Form("h1D_TPC_Pi_%d", i), lowBinMult, highBinMult, -1, -1);
        h1DNsigmaTOFKaon[i] = (TH1F *)hNsigmaTOFKaon->ProjectionX(Form("h1D_TOF_Ka_%d", i), lowBinMult, highBinMult, -1, -1);
        h1DNsigmaTOFPion[i] = (TH1F *)hNsigmaTOFPion->ProjectionX(Form("h1D_TOF_Pi_%d", i), lowBinMult, highBinMult, -1, -1);

        h1DNsigmaTPCKaon_before[i] = (TH1F *)hNsigmaTPCKaon_before->ProjectionX(Form("h1D_TPC_Ka_before_%d", i), lowBinMult, highBinMult, -1, -1);
        h1DNsigmaTPCPion_before[i] = (TH1F *)hNsigmaTPCPion_before->ProjectionX(Form("h1D_TPC_Pi_before_%d", i), lowBinMult, highBinMult, -1, -1);
        h1DNsigmaTOFKaon_before[i] = (TH1F *)hNsigmaTOFKaon_before->ProjectionX(Form("h1D_TOF_Ka_before_%d", i), lowBinMult, highBinMult, -1, -1);
        h1DNsigmaTOFPion_before[i] = (TH1F *)hNsigmaTOFPion_before->ProjectionX(Form("h1D_TOF_Pi_before_%d", i), lowBinMult, highBinMult, -1, -1);

        int lowBinpT = hNsigmaTPCKaon->GetZaxis()->FindBin(pTrangesForTPCandTOF[i] + 0.001);
        int highBinpT = hNsigmaTPCKaon->GetZaxis()->FindBin(pTrangesForTPCandTOF[i + 1] - 0.001);

        h1DNsigmaTPCKaon_pt[i] = (TH1F *)hNsigmaTPCKaon->ProjectionX(Form("h1D_TPC_Ka_pt_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT);
        h1DNsigmaTPCPion_pt[i] = (TH1F *)hNsigmaTPCPion->ProjectionX(Form("h1D_TPC_Pi_pt_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT);
        h1DNsigmaTOFKaon_pt[i] = (TH1F *)hNsigmaTOFKaon->ProjectionX(Form("h1D_TOF_Ka_pt_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT);
        h1DNsigmaTOFPion_pt[i] = (TH1F *)hNsigmaTOFPion->ProjectionX(Form("h1D_TOF_Pi_pt_%d", i), lowBinMult, highBinMult, lowBinpT, highBinpT);
    }

    TCanvas *cNsigmaTPCKaon = new TCanvas("cNsigmaTPCKaon", "Nsigma TPC Kaon", 720, 720);
    SetCanvasStyle(cNsigmaTPCKaon, 0.14, 0.15, 0.06, 0.14);
    SetHistoQA2D(h2DTPCKaon);
    // gPad->SetLogz(1);
    h2DTPCKaon->GetYaxis()->SetTitle("n#sigma_{TPC}");
    h2DTPCKaon->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h2DTPCKaon->GetYaxis()->SetRangeUser(-3, 3);
    h2DTPCKaon->GetXaxis()->SetRangeUser(0, 5);
    h2DTPCKaon->Draw("colz");
    cNsigmaTPCKaon->SaveAs(output_QA_folder + ("/NsigmaTPCKaon2D." + outputtype).c_str());

    TCanvas *cNsigmaTPCPion = new TCanvas("cNsigmaTPCPion", "Nsigma TPC Pion", 720, 720);
    SetCanvasStyle(cNsigmaTPCPion, 0.14, 0.15, 0.06, 0.14);
    SetHistoQA2D(h2DTPCPion);
    // gPad->SetLogz(1);
    h2DTPCPion->GetYaxis()->SetTitle("n#sigma_{TPC}");
    h2DTPCPion->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h2DTPCPion->GetYaxis()->SetRangeUser(-3, 3);
    h2DTPCPion->GetXaxis()->SetRangeUser(0, 5);
    h2DTPCPion->Draw("colz");
    cNsigmaTPCPion->SaveAs(output_QA_folder + ("/NsigmaTPCPion2D." + outputtype).c_str());

    TCanvas *cNsigmaTOFKaon = new TCanvas("cNsigmaTOFKaon", "Nsigma TOF Kaon", 720, 720);
    SetCanvasStyle(cNsigmaTOFKaon, 0.14, 0.15, 0.06, 0.14);
    SetHistoQA2D(h2DTOFKaon);
    // gPad->SetLogz(1);
    h2DTOFKaon->GetYaxis()->SetTitle("n#sigma_{TOF}");
    h2DTOFKaon->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h2DTOFKaon->GetYaxis()->SetRangeUser(-3, 3);
    h2DTOFKaon->GetXaxis()->SetRangeUser(0, 5);
    h2DTOFKaon->Draw("colz");
    cNsigmaTOFKaon->SaveAs(output_QA_folder + ("/NsigmaTOFKaon2D." + outputtype).c_str());

    TCanvas *cNsigmaTOFPion = new TCanvas("cNsigmaTOFPion", "Nsigma TOF Pion", 720, 720);
    SetCanvasStyle(cNsigmaTOFPion, 0.14, 0.15, 0.06, 0.14);
    SetHistoQA2D(h2DTOFPion);
    // gPad->SetLogz(1);
    h2DTOFPion->GetYaxis()->SetTitle("n#sigma_{TOF}");
    h2DTOFPion->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    h2DTOFPion->GetYaxis()->SetRangeUser(-3, 3);
    h2DTOFPion->GetXaxis()->SetRangeUser(0, 5);
    h2DTOFPion->Draw("colz");
    cNsigmaTOFPion->SaveAs(output_QA_folder + ("/NsigmaTOFPion2D." + outputtype).c_str());

    TCanvas *cNsigmaTPCTOFKaon = new TCanvas("cNsigmaTPCTOFKaon", "Nsigma TPC vs TOF Kaon", 1080, 720);
    TCanvas *cNsigmaTPCTOFPion = new TCanvas("cNsigmaTPCTOFPion", "Nsigma TPC vs TOF Pion", 1080, 720);
    SetCanvasStyle(cNsigmaTPCTOFKaon, 0.14, 0.15, 0.06, 0.14);
    SetCanvasStyle(cNsigmaTPCTOFPion, 0.14, 0.15, 0.06, 0.14);
    cNsigmaTPCTOFKaon->Divide(3, 2);
    cNsigmaTPCTOFPion->Divide(3, 2);

    TLatex latPID;
    latPID.SetNDC();
    latPID.SetTextFont(42);
    latPID.SetTextSize(0.06);
    latPID.SetTextColor(kCyan);

    for (int i = 0; i < 6; i++)
    {
        cNsigmaTPCTOFKaon->cd(i + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        gPad->SetBottomMargin(0.13);
        gPad->SetTopMargin(0.06);
        SetCanvasStyle(cNsigmaTPCTOFKaon, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(h2DNsigmaTPCTOFKaon[i]);
        h2DNsigmaTPCTOFKaon[i]->GetXaxis()->SetTitle("n#sigma_{TPC}");
        h2DNsigmaTPCTOFKaon[i]->GetYaxis()->SetTitle("n#sigma_{TOF}");
        h2DNsigmaTPCTOFKaon[i]->GetYaxis()->SetRangeUser(-3, 3);
        h2DNsigmaTPCTOFKaon[i]->GetXaxis()->SetRangeUser(-3, 3);
        h2DNsigmaTPCTOFKaon[i]->Draw("colz");
        latPID.DrawLatex(0.2, 0.85, Form("p_{T}: %.1f-%.1f GeV/c", pTrangesForTPCandTOF[i], pTrangesForTPCandTOF[i + 1]));
        // cout<<"pT range: "<<pTrangesForTPCandTOF[i]<<"-"<<pTrangesForTPCandTOF[i + 1]<<" GeV/c, mean TPC nSigma: "<<h2DNsigmaTPCTOFKaon[i]->GetMean(1)<<", mean TOF nSigma: "<<h2DNsigmaTPCTOFKaon[i]->GetMean(2)<<endl;

        cNsigmaTPCTOFPion->cd(i + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        gPad->SetBottomMargin(0.13);
        gPad->SetTopMargin(0.06);
        SetCanvasStyle(cNsigmaTPCTOFPion, 0.14, 0.15, 0.06, 0.14);
        SetHistoQA2D(h2DNsigmaTPCTOFPion[i]);
        h2DNsigmaTPCTOFPion[i]->GetXaxis()->SetTitle("n#sigma_{TPC}");
        h2DNsigmaTPCTOFPion[i]->GetYaxis()->SetTitle("n#sigma_{TOF}");
        h2DNsigmaTPCTOFPion[i]->GetYaxis()->SetRangeUser(-3, 3);
        h2DNsigmaTPCTOFPion[i]->GetXaxis()->SetRangeUser(-3, 3);
        h2DNsigmaTPCTOFPion[i]->Draw("colz");
        latPID.DrawLatex(0.2, 0.85, Form("p_{T}: %.1f-%.1f GeV/c", pTrangesForTPCandTOF[i], pTrangesForTPCandTOF[i + 1]));
    }

    cNsigmaTPCTOFKaon->SaveAs(output_QA_folder + ("/NsigmaTPCTOFKaon_2D." + outputtype).c_str());
    cNsigmaTPCTOFPion->SaveAs(output_QA_folder + ("/NsigmaTPCTOFPion_2D." + outputtype).c_str());

    TCanvas *cNsigmaTPCKaon_1D = new TCanvas("cNsigmaTPCKaon_1D", "Nsigma TPC Kaon 1D", 720, 720);
    SetCanvasStyle(cNsigmaTPCKaon_1D, 0.14, 0.05, 0.06, 0.14);
    cNsigmaTPCKaon_1D->SetGrid(1, 1);

    TCanvas *cNsigmaTPCKaon_1D_before = new TCanvas("cNsigmaTPCKaon_1D_before", "Nsigma TPC Kaon 1D", 720, 720);
    SetCanvasStyle(cNsigmaTPCKaon_1D_before, 0.14, 0.05, 0.06, 0.14);
    cNsigmaTPCKaon_1D_before->SetGrid(1, 1);

    TCanvas *cNsigmaTPCKaon_1D_pt = new TCanvas("cNsigmaTPCKaon_1D_pt", "Nsigma TPC Kaon 1D", 720, 720);
    SetCanvasStyle(cNsigmaTPCKaon_1D_pt, 0.14, 0.05, 0.06, 0.14);
    cNsigmaTPCKaon_1D_pt->SetGrid(1, 1);

    TLegend *leg1DPID = new TLegend(0.17, 0.82, 0.6, 0.92);
    leg1DPID->SetNColumns(3);
    leg1DPID->SetTextSize(0.028);
    leg1DPID->SetFillStyle(0);
    leg1DPID->SetBorderSize(0);
    leg1DPID->SetTextFont(42);

    TLegend *leg1DPID_pt = new TLegend(0.17, 0.82, 0.8, 0.92);
    leg1DPID_pt->SetNColumns(3);
    leg1DPID_pt->SetTextSize(0.028);
    leg1DPID_pt->SetFillStyle(0);
    leg1DPID_pt->SetBorderSize(0);
    leg1DPID_pt->SetTextFont(42);
    leg1DPID_pt->SetHeader("#it{p}_{T} ranges");

    for (int i = 0; i < 6; i++)
    {
        cNsigmaTPCKaon_1D->cd();
        SetHistoQA(h1DNsigmaTPCKaon[i]);
        h1DNsigmaTPCKaon[i]->SetMarkerColor(colors[i]);
        h1DNsigmaTPCKaon[i]->GetXaxis()->SetTitle("n#sigma_{TPC} K^{#pm}");
        h1DNsigmaTPCKaon[i]->GetYaxis()->SetTitle("Counts");
        h1DNsigmaTPCKaon[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        h1DNsigmaTPCKaon[i]->SetMaximum(h1DNsigmaTPCKaon[i]->GetMaximum() * 15);
        h1DNsigmaTPCKaon[i]->Draw("p same");
        leg1DPID->AddEntry(h1DNsigmaTPCKaon[i], Form("%d-%d%%", multRangesForTPCandTOF[i][0], multRangesForTPCandTOF[i][1]), "p");
        leg1DPID_pt->AddEntry(h1DNsigmaTPCKaon_pt[i], Form("%.1f-%.1f", pTrangesForTPCandTOF[i], pTrangesForTPCandTOF[i + 1]), "p");

        cNsigmaTPCKaon_1D_before->cd();
        SetHistoQA(h1DNsigmaTPCKaon_before[i]);
        h1DNsigmaTPCKaon_before[i]->SetMarkerColor(colors[i]);
        h1DNsigmaTPCKaon_before[i]->GetXaxis()->SetTitle("n#sigma_{TPC} K^{#pm}");
        h1DNsigmaTPCKaon_before[i]->GetYaxis()->SetTitle("Counts");
        h1DNsigmaTPCKaon_before[i]->GetXaxis()->SetRangeUser(-5, 5);
        h1DNsigmaTPCKaon_before[i]->SetMaximum(h1DNsigmaTPCKaon_before[i]->GetMaximum() * 15);
        h1DNsigmaTPCKaon_before[i]->Draw("p same");

        cNsigmaTPCKaon_1D_pt->cd();
        SetHistoQA(h1DNsigmaTPCKaon_pt[i]);
        h1DNsigmaTPCKaon_pt[i]->SetMarkerColor(colors[i]);
        h1DNsigmaTPCKaon_pt[i]->GetXaxis()->SetTitle("n#sigma_{TPC} K^{#pm}");
        h1DNsigmaTPCKaon_pt[i]->GetYaxis()->SetTitle("Counts");
        h1DNsigmaTPCKaon_pt[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        h1DNsigmaTPCKaon_pt[i]->SetMaximum(h1DNsigmaTPCKaon_pt[i]->GetMaximum() * 15);
        h1DNsigmaTPCKaon_pt[i]->Draw("p same");
    }
    cNsigmaTPCKaon_1D->cd();
    leg1DPID->Draw();
    TLine *lineTPCKaon = new TLine(0, 0, 0, h1DNsigmaTPCKaon[0]->GetMaximum());
    lineTPCKaon->SetLineStyle(2);
    lineTPCKaon->SetLineColor(2);
    lineTPCKaon->Draw();
    cNsigmaTPCKaon_1D->SaveAs(output_QA_folder + ("/NsigmaTPCKaon_1D." + outputtype).c_str());

    cNsigmaTPCKaon_1D_before->cd();
    leg1DPID->Draw();
    lineTPCKaon->SetY2(h1DNsigmaTPCKaon_before[0]->GetMaximum());
    lineTPCKaon->Draw();
    cNsigmaTPCKaon_1D_before->SaveAs(output_QA_folder_before + ("/NsigmaTPCKaon_1D_before." + outputtype).c_str());

    cNsigmaTPCKaon_1D_pt->cd();
    leg1DPID_pt->Draw();
    lineTPCKaon->SetY2(h1DNsigmaTPCKaon_pt[0]->GetMaximum());
    lineTPCKaon->Draw();
    cNsigmaTPCKaon_1D_pt->SaveAs(output_QA_folder + ("/NsigmaTPCKaon_1D_pt." + outputtype).c_str());

    TCanvas *cNsigmaTPCPion_1D = new TCanvas("cNsigmaTPCPion_1D", "Nsigma TPC Pion 1D", 720, 720);
    SetCanvasStyle(cNsigmaTPCPion_1D, 0.14, 0.05, 0.06, 0.14);
    cNsigmaTPCPion_1D->SetGrid(1, 1);

    TCanvas *cNsigmaTPCPion_1D_before = new TCanvas("cNsigmaTPCPion_1D_before", "Nsigma TPC Pion 1D", 720, 720);
    SetCanvasStyle(cNsigmaTPCPion_1D_before, 0.14, 0.05, 0.06, 0.14);
    cNsigmaTPCPion_1D_before->SetGrid(1, 1);

    TCanvas *cNsigmaTPCPion_1D_pt = new TCanvas("cNsigmaTPCPion_1D_pt", "Nsigma TPC Pion 1D", 720, 720);
    SetCanvasStyle(cNsigmaTPCPion_1D_pt, 0.14, 0.05, 0.06, 0.14);
    cNsigmaTPCPion_1D_pt->SetGrid(1, 1);

    for (int i = 0; i < 6; i++)
    {
        cNsigmaTPCPion_1D->cd();
        SetHistoQA(h1DNsigmaTPCPion[i]);
        h1DNsigmaTPCPion[i]->SetMarkerColor(colors[i]);
        h1DNsigmaTPCPion[i]->GetXaxis()->SetTitle("n#sigma_{TPC} #pi^{#pm}");
        h1DNsigmaTPCPion[i]->GetYaxis()->SetTitle("Counts");
        h1DNsigmaTPCPion[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        h1DNsigmaTPCPion[i]->SetMaximum(h1DNsigmaTPCPion[i]->GetMaximum() * 15);
        h1DNsigmaTPCPion[i]->Draw("p same");

        cNsigmaTPCPion_1D_before->cd();
        SetHistoQA(h1DNsigmaTPCPion_before[i]);
        h1DNsigmaTPCPion_before[i]->SetMarkerColor(colors[i]);
        h1DNsigmaTPCPion_before[i]->GetXaxis()->SetTitle("n#sigma_{TPC} #pi^{#pm}");
        h1DNsigmaTPCPion_before[i]->GetYaxis()->SetTitle("Counts");
        h1DNsigmaTPCPion_before[i]->GetXaxis()->SetRangeUser(-5, 5);
        h1DNsigmaTPCPion_before[i]->SetMaximum(h1DNsigmaTPCPion_before[i]->GetMaximum() * 15);
        h1DNsigmaTPCPion_before[i]->Draw("p same");

        cNsigmaTPCPion_1D_pt->cd();
        SetHistoQA(h1DNsigmaTPCPion_pt[i]);
        h1DNsigmaTPCPion_pt[i]->SetMarkerColor(colors[i]);
        h1DNsigmaTPCPion_pt[i]->GetXaxis()->SetTitle("n#sigma_{TPC} #pi^{#pm}");
        h1DNsigmaTPCPion_pt[i]->GetYaxis()->SetTitle("Counts");
        h1DNsigmaTPCPion_pt[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        h1DNsigmaTPCPion_pt[i]->SetMaximum(h1DNsigmaTPCPion_pt[i]->GetMaximum() * 12);
        h1DNsigmaTPCPion_pt[i]->Draw("p same");
    }
    cNsigmaTPCPion_1D->cd();
    leg1DPID->Draw();
    TLine *lineTPCPion = new TLine(0, 0, 0, h1DNsigmaTPCPion[0]->GetMaximum());
    lineTPCPion->SetLineStyle(2);
    lineTPCPion->SetLineColor(2);
    lineTPCPion->Draw();
    cNsigmaTPCPion_1D->SaveAs(output_QA_folder + ("/NsigmaTPCPion_1D." + outputtype).c_str());

    cNsigmaTPCPion_1D_before->cd();
    leg1DPID->Draw();
    lineTPCPion->SetY2(h1DNsigmaTPCPion_before[0]->GetMaximum());
    lineTPCPion->Draw();
    cNsigmaTPCPion_1D_before->SaveAs(output_QA_folder_before + ("/NsigmaTPCPion_1D_before." + outputtype).c_str());

    cNsigmaTPCPion_1D_pt->cd();
    leg1DPID_pt->Draw();
    lineTPCPion->SetY2(h1DNsigmaTPCPion_pt[0]->GetMaximum());
    lineTPCPion->Draw();
    cNsigmaTPCPion_1D_pt->SaveAs(output_QA_folder + ("/NsigmaTPCPion_1D_pt." + outputtype).c_str());

    TCanvas *cNsigmaTOFKaon_1D = new TCanvas("cNsigmaTOFKaon_1D", "Nsigma TOF Kaon 1D", 720, 720);
    SetCanvasStyle(cNsigmaTOFKaon_1D, 0.14, 0.05, 0.06, 0.14);
    cNsigmaTOFKaon_1D->SetGrid(1, 1);

    TCanvas *cNsigmaTOFKaon_1D_before = new TCanvas("cNsigmaTOFKaon_1D_before", "Nsigma TOF Kaon 1D", 720, 720);
    SetCanvasStyle(cNsigmaTOFKaon_1D_before, 0.14, 0.05, 0.06, 0.14);
    cNsigmaTOFKaon_1D_before->SetGrid(1, 1);

    TCanvas *cNsigmaTOFKaon_1D_pt = new TCanvas("cNsigmaTOFKaon_1D_pt", "Nsigma TOF Kaon 1D", 720, 720);
    SetCanvasStyle(cNsigmaTOFKaon_1D_pt, 0.14, 0.05, 0.06, 0.14);
    cNsigmaTOFKaon_1D_pt->SetGrid(1, 1);

    TF1 *fitGaus = new TF1("fitGaus", "gaus", -2, 2);
    fitGaus->SetParameter(0, 1);
    fitGaus->SetParameter(1, -0.8);
    fitGaus->SetParameter(2, 1);
    for (int i = 0; i < 6; i++)
    {
        cNsigmaTOFKaon_1D->cd();
        SetHistoQA(h1DNsigmaTOFKaon[i]);
        h1DNsigmaTOFKaon[i]->SetMarkerColor(colors[i]);
        h1DNsigmaTOFKaon[i]->GetXaxis()->SetTitle("n#sigma_{TOF} K^{#pm}");
        h1DNsigmaTOFKaon[i]->GetYaxis()->SetTitle("Counts");
        h1DNsigmaTOFKaon[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        h1DNsigmaTOFKaon[i]->SetMaximum(h1DNsigmaTOFKaon[i]->GetMaximum() * 15);
        h1DNsigmaTOFKaon[i]->Draw("p same");

        cNsigmaTOFKaon_1D_before->cd();
        SetHistoQA(h1DNsigmaTOFKaon_before[i]);
        h1DNsigmaTOFKaon_before[i]->SetMarkerColor(colors[i]);
        h1DNsigmaTOFKaon_before[i]->GetXaxis()->SetTitle("n#sigma_{TOF} K^{#pm}");
        h1DNsigmaTOFKaon_before[i]->GetYaxis()->SetTitle("Counts");
        h1DNsigmaTOFKaon_before[i]->GetXaxis()->SetRangeUser(-5, 5);
        h1DNsigmaTOFKaon_before[i]->SetMaximum(h1DNsigmaTOFKaon_before[i]->GetMaximum() * 15);
        h1DNsigmaTOFKaon_before[i]->Draw("p same");

        cNsigmaTOFKaon_1D_pt->cd();
        SetHistoQA(h1DNsigmaTOFKaon_pt[i]);
        h1DNsigmaTOFKaon_pt[i]->SetMarkerColor(colors[i]);
        h1DNsigmaTOFKaon_pt[i]->GetXaxis()->SetTitle("n#sigma_{TOF} K^{#pm}");
        h1DNsigmaTOFKaon_pt[i]->GetYaxis()->SetTitle("Counts");
        h1DNsigmaTOFKaon_pt[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        h1DNsigmaTOFKaon_pt[i]->SetMaximum(h1DNsigmaTOFKaon_pt[i]->GetMaximum() * 100);
        h1DNsigmaTOFKaon_pt[i]->Draw("p same");
    }
    cNsigmaTOFKaon_1D->cd();
    // h1DNsigmaTOFKaon[3]->Fit(fitGaus, "REBMS");
    h1DNsigmaTOFKaon[3]->SetStats(1);
    leg1DPID->Draw();
    TLine *lineTOFKaon = new TLine(0, 0, 0, h1DNsigmaTOFKaon[0]->GetMaximum());
    lineTOFKaon->SetLineStyle(2);
    lineTOFKaon->SetLineColor(2);
    lineTOFKaon->Draw();
    cNsigmaTOFKaon_1D->SaveAs(output_QA_folder + ("/NsigmaTOFKaon_1D." + outputtype).c_str());

    cNsigmaTOFKaon_1D_before->cd();
    leg1DPID->Draw();
    lineTOFKaon->SetY2(h1DNsigmaTOFKaon_before[0]->GetMaximum());
    lineTOFKaon->Draw();
    cNsigmaTOFKaon_1D_before->SaveAs(output_QA_folder_before + ("/NsigmaTOFKaon_1D_before." + outputtype).c_str());

    cNsigmaTOFKaon_1D_pt->cd();
    leg1DPID_pt->Draw();
    lineTOFKaon->SetY2(h1DNsigmaTOFKaon_pt[0]->GetMaximum());
    lineTOFKaon->Draw();
    cNsigmaTOFKaon_1D_pt->SaveAs(output_QA_folder + ("/NsigmaTOFKaon_1D_pt." + outputtype).c_str());

    TCanvas *cNsigmaTOFPion_1D = new TCanvas("cNsigmaTOFPion_1D", "Nsigma TOF Pion 1D", 720, 720);
    SetCanvasStyle(cNsigmaTOFPion_1D, 0.14, 0.05, 0.06, 0.14);
    cNsigmaTOFPion_1D->SetGrid(1, 1);

    TCanvas *cNsigmaTOFPion_1D_before = new TCanvas("cNsigmaTOFPion_1D_before", "Nsigma TOF Pion 1D", 720, 720);
    SetCanvasStyle(cNsigmaTOFPion_1D_before, 0.14, 0.05, 0.06, 0.14);
    cNsigmaTOFPion_1D_before->SetGrid(1, 1);

    TCanvas *cNsigmaTOFPion_1D_pt = new TCanvas("cNsigmaTOFPion_1D_pt", "Nsigma TOF Pion 1D", 720, 720);
    SetCanvasStyle(cNsigmaTOFPion_1D_pt, 0.14, 0.05, 0.06, 0.14);
    cNsigmaTOFPion_1D_pt->SetGrid(1, 1);

    for (int i = 0; i < 6; i++)
    {
        cNsigmaTOFPion_1D->cd();
        SetHistoQA(h1DNsigmaTOFPion[i]);
        h1DNsigmaTOFPion[i]->SetMarkerColor(colors[i]);
        h1DNsigmaTOFPion[i]->GetXaxis()->SetTitle("n#sigma_{TOF} #pi^{#pm}");
        h1DNsigmaTOFPion[i]->GetYaxis()->SetTitle("Counts");
        h1DNsigmaTOFPion[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        h1DNsigmaTOFPion[i]->SetMaximum(h1DNsigmaTOFPion[i]->GetMaximum() * 15);
        h1DNsigmaTOFPion[i]->Draw("p same");

        cNsigmaTOFPion_1D_before->cd();
        SetHistoQA(h1DNsigmaTOFPion_before[i]);
        h1DNsigmaTOFPion_before[i]->SetMarkerColor(colors[i]);
        h1DNsigmaTOFPion_before[i]->GetXaxis()->SetTitle("n#sigma_{TOF} #pi^{#pm}");
        h1DNsigmaTOFPion_before[i]->GetYaxis()->SetTitle("Counts");
        h1DNsigmaTOFPion_before[i]->GetXaxis()->SetRangeUser(-5, 5);
        h1DNsigmaTOFPion_before[i]->SetMaximum(h1DNsigmaTOFPion_before[i]->GetMaximum() * 15);
        h1DNsigmaTOFPion_before[i]->Draw("p same");

        cNsigmaTOFPion_1D_pt->cd();
        SetHistoQA(h1DNsigmaTOFPion_pt[i]);
        h1DNsigmaTOFPion_pt[i]->SetMarkerColor(colors[i]);
        h1DNsigmaTOFPion_pt[i]->GetXaxis()->SetTitle("n#sigma_{TOF} #pi^{#pm}");
        h1DNsigmaTOFPion_pt[i]->GetYaxis()->SetTitle("Counts");
        h1DNsigmaTOFPion_pt[i]->GetXaxis()->SetRangeUser(-2.5, 2.5);
        h1DNsigmaTOFPion_pt[i]->SetMaximum(h1DNsigmaTOFPion_pt[i]->GetMaximum() * 25);
        h1DNsigmaTOFPion_pt[i]->Draw("p same");
    }
    cNsigmaTOFPion_1D->cd();
    // h1DNsigmaTOFPion[5]->Fit(fitGaus, "REBMS");
    leg1DPID->Draw();
    TLine *lineTOFPion = new TLine(0, 0, 0, h1DNsigmaTOFPion[0]->GetMaximum());
    lineTOFPion->SetLineStyle(2);
    lineTOFPion->SetLineColor(2);
    lineTOFPion->Draw();
    cNsigmaTOFPion_1D->SaveAs(output_QA_folder + ("/NsigmaTOFPion_1D." + outputtype).c_str());

    cNsigmaTOFPion_1D_before->cd();
    leg1DPID->Draw();
    lineTOFPion->SetY2(h1DNsigmaTOFPion_before[0]->GetMaximum());
    lineTOFPion->Draw();
    cNsigmaTOFPion_1D_before->SaveAs(output_QA_folder_before + ("/NsigmaTOFPion_1D_before." + outputtype).c_str());

    cNsigmaTOFPion_1D_pt->cd();
    leg1DPID_pt->Draw();
    lineTOFPion->SetY2(h1DNsigmaTOFPion_pt[0]->GetMaximum());
    lineTOFPion->Draw();
    cNsigmaTOFPion_1D_pt->SaveAs(output_QA_folder + ("/NsigmaTOFPion_1D_pt." + outputtype).c_str());
}
