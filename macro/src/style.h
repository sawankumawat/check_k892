#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

void SetHistoStyle(TH1 *h, Int_t mcolor, Int_t mstyle, Float_t msize, Float_t Tsizex, Float_t Tsizey, Float_t Lsizex, Float_t Lsizey, Float_t Offsetx, Float_t Offsety)
{
    h->SetMarkerColor(mcolor);
    h->SetMarkerStyle(mstyle);
    h->SetMarkerSize(msize);
    h->GetXaxis()->SetTitleSize(Tsizex);
    h->GetYaxis()->SetTitleSize(Tsizey);
    h->GetXaxis()->SetLabelSize(Lsizex);
    h->GetYaxis()->SetLabelSize(Lsizey);
    h->GetXaxis()->SetTitleOffset(Offsetx);
    h->GetYaxis()->SetTitleOffset(Offsety);
    h->GetXaxis()->CenterTitle(true);
    h->GetYaxis()->CenterTitle(true);
}

void SetgrStyle(TGraph *h, Int_t mcolor, Int_t mstyle, Float_t msize, Float_t Tsizex, Float_t Tsizey, Float_t Lsizex, Float_t Lsizey, Float_t Offsetx, Float_t Offsety)
{
    h->SetMarkerColor(mcolor);
    h->SetMarkerStyle(mstyle);
    h->SetMarkerSize(msize);
    h->GetXaxis()->SetTitleSize(Tsizex);
    h->GetYaxis()->SetTitleSize(Tsizey);
    h->GetXaxis()->SetLabelSize(Lsizex);
    h->GetYaxis()->SetLabelSize(Lsizey);
    h->GetXaxis()->SetTitleOffset(Offsetx);
    h->GetYaxis()->SetTitleOffset(Offsety);
    h->GetXaxis()->CenterTitle(true);
    h->GetYaxis()->CenterTitle(true);
}

void SetHistoQA(TH1 *h)
{
    h->SetMarkerStyle(8);
    h->SetMarkerSize(1);
    h->SetMarkerColor(1);
    h->SetLineColor(1);
    h->SetLineWidth(2);
    h->SetTitle(0);
    // h->GetXaxis()->SetNdivisions(506);
    // h->GetYaxis()->SetNdivisions(505);
    h->GetXaxis()->SetLabelOffset(0.015);
    h->GetXaxis()->SetLabelFont(42);
    h->GetXaxis()->SetTitleFont(42);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTickLength(0.02);
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetXaxis()->CenterTitle(true);
    h->GetYaxis()->CenterTitle(true);
    h->GetYaxis()->SetTitleOffset(1.4);
    h->GetYaxis()->SetDecimals(false);
    h->GetYaxis()->SetLabelOffset(0.015);
    h->GetYaxis()->SetLabelFont(42);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTickLength(0.02);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleFont(42);
    h->GetYaxis()->SetMaxDigits(3);
}

void SetHistoQA2D(TH2 *h)
{
    h->SetMarkerStyle(8);
    h->SetMarkerSize(0.5);
    h->SetMarkerColor(1);
    h->SetLineColor(1);
    h->SetLineWidth(2);
    h->SetTitle(0);
    // h->GetXaxis()->SetNdivisions(506);
    // h->GetYaxis()->SetNdivisions(505);
    h->GetXaxis()->SetLabelOffset(0.015);
    h->GetXaxis()->SetLabelFont(42);
    h->GetXaxis()->SetTitleFont(42);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTickLength(0.04);
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetXaxis()->CenterTitle(true);

    h->GetYaxis()->CenterTitle(true);
    h->GetYaxis()->SetTitleOffset(1.4);
    h->GetYaxis()->SetDecimals(false);
    h->GetYaxis()->SetLabelOffset(0.015);
    h->GetYaxis()->SetLabelFont(42);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTickLength(0.04);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleFont(42);

    h->GetZaxis()->CenterTitle(true);
    h->GetZaxis()->SetTitleOffset(1.4);
    h->GetZaxis()->SetDecimals(false);
    h->GetZaxis()->SetLabelOffset(0.015);
    h->GetZaxis()->SetLabelFont(42);
    h->GetZaxis()->SetLabelSize(0.045);
    h->GetZaxis()->SetTickLength(0.04);
    h->GetZaxis()->SetTitleSize(0.05);
    h->GetZaxis()->SetTitleFont(42);
    h->GetZaxis()->SetMaxDigits(3);
}

void SetGrapherrorStyle(TGraphErrors *gr)
{
    gr->SetTitle(0);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.0);
    gr->SetLineWidth(2);
    gr->GetXaxis()->CenterTitle(true);
    gr->GetXaxis()->SetLabelOffset(0.015);
    gr->GetXaxis()->SetLabelFont(42);
    gr->GetXaxis()->SetTitleFont(42);
    gr->GetXaxis()->SetLabelSize(0.045);
    gr->GetXaxis()->SetTitleSize(0.05);
    gr->GetXaxis()->SetTickLength(0.04);
    gr->GetXaxis()->SetTitleOffset(1.4);
    gr->GetYaxis()->SetTitleOffset(1.5);
    gr->GetYaxis()->CenterTitle(true);
    gr->GetYaxis()->SetLabelOffset(0.015);
    gr->GetYaxis()->SetLabelFont(42);
    gr->GetYaxis()->SetLabelSize(0.045);
    gr->GetYaxis()->SetTickLength(0.04);
    gr->GetYaxis()->SetTitleSize(0.05);
    gr->GetYaxis()->SetTitleFont(42);
}

void SetGraphStyle(TGraph *gr, int markercolor, int linecolor)
{
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(markercolor);
    gr->SetLineColor(linecolor);
    gr->SetMarkerSize(1.4);
    gr->SetLineWidth(2);
    gr->GetXaxis()->CenterTitle(true);
    gr->GetXaxis()->SetNdivisions(506);
    gr->GetYaxis()->SetNdivisions(505);
    gr->GetXaxis()->SetLabelOffset(0.015);
    gr->GetXaxis()->SetLabelFont(42);
    gr->GetXaxis()->SetTitleFont(42);
    gr->GetXaxis()->SetLabelSize(0.04);
    gr->GetXaxis()->SetTitleSize(0.045);
    gr->GetXaxis()->SetTickLength(0.04);
    gr->GetXaxis()->SetTitleOffset(1.2);
    gr->GetYaxis()->SetTitleOffset(1.3);
    gr->GetYaxis()->CenterTitle(true);
    gr->GetYaxis()->SetDecimals(false);
    gr->GetYaxis()->SetLabelOffset(0.015);
    gr->GetYaxis()->SetLabelFont(42);
    gr->GetYaxis()->SetLabelSize(0.04);
    gr->GetYaxis()->SetTickLength(0.04);
    gr->GetYaxis()->SetTitleSize(0.045);
    gr->GetYaxis()->SetTitleFont(42);
}

void SetCanvasStyle(TCanvas *c, float leftmargin, float rightmargin, float topmargin, float bottommargin)
{
    c->Range(0, 0, 1, 1);
    // c->SetBorderSize(2);
    c->SetBorderMode(0);
    c->SetFillColor(10);
    c->SetFrameFillColor(10);
    c->SetFrameLineWidth(2);
    // c->SetLeftMargin(0.2);
    // c->SetRightMargin(0.05);
    // c->SetTopMargin(0.08);
    // c->SetBottomMargin(0.2);
    c->SetLeftMargin(leftmargin);
    c->SetRightMargin(rightmargin);
    c->SetTopMargin(topmargin);
    c->SetBottomMargin(bottommargin);
    c->SetTicks(1, 1);
}

void SetLegendStyle(TLegend *l)
{
    l->SetBorderSize(0);
    l->SetTextFont(42);
    l->SetFillStyle(0);
}

void SetLineStyle(TLine *line, int linecolor)
{
    line->SetLineStyle(2);
    line->SetLineColor(linecolor);
    line->SetLineWidth(3);
    line->Draw();
}

//*****TLegend Class******************************************************************

TLegend *DrawLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{

    TLegend *legend = new TLegend(x1, y1, x2, y2);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->SetLineColor(0);
    legend->SetShadowColor(0);
    return legend;
}

void gstyle()
{
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
    // gStyle->SetErrorX(0);
}
