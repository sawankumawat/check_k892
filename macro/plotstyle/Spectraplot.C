#include<iostream>
#include<cmath>

using namespace std;


TLegend *DrawLegend(Double_t x1,Double_t y1,Double_t x2,Double_t y2)
{

  TLegend *legend = new TLegend(x1,y1,x2,y2);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->SetLineColor(0);
  legend->SetShadowColor(0);
  return legend;
}


 


void Spectraplot()

{

  MyStyle();

  TFile *f0 = new TFile("Correctedspec0_100_rebin1stbin.root","READ");  
  TH1D  *h0 = (TH1D*)f0->Get("hYield_kstar_stat");

  
  TH1D *h_0=Plothisto(h0, kRed, 24); //v2default
  
 
 TH1D *h_dummy=new TH1D("h_dummy","h_dummy",1,0.0,15.0);
 h_dummy=Plothisto(h_dummy,1,1);
 //h_dummy->SetMinimum(0.6);
 //h_dummy->SetMaximum(1.6);
 h_dummy->SetMinimum(0.0000000010);
 h_dummy->SetMaximum(0.00010);

 

 TCanvas *c11 = DrawCanvas("c11");
 c11->cd();
 TPad *pad = MyPad(0.01,0.01,0.99, 0.99, 0.15,0.06,0.02,0.13, "pad");
 pad->cd();
 h_dummy->Draw("X0");

 h_0->SetFillStyle(1);
 h_0->Draw("PSAME"); 


 TLegend *lp1=DrawLegend(0.6,0.85,0.7,0.9);
lp1->SetTextFont(22);
lp1->SetBorderSize(0);
lp1->SetFillStyle(0);
lp1->SetFillColor(0);
lp1->SetTextSize(0.04);
lp1->SetShadowColor(1);
lp1->AddEntry((TObject*)0,"ALICE","");

TLegend *lp2=DrawLegend(0.6,0.65,0.7,0.82);

lp2->SetTextFont(42);
lp2->SetBorderSize(0);
lp2->SetFillStyle(0);
lp2->SetFillColor(0);
lp2->SetTextSize(0.04);
lp2->SetShadowColor(1);
lp2->AddEntry((TObject*)0,"f_{1}(1285)","");
lp2->AddEntry((TObject*)0,"pp #sqrt{#it{s}} = 13 TeV","");
 lp2->AddEntry((TObject*)0,"0-100 %","");


 

   TLegend *lp = DrawLegend(0.45,0.18,0.6,0.22);
  lp->SetTextSize(.04);
  lp->SetHeader("Uncertainties: stat.(bars)");
  

  lp->Draw();
  lp1->Draw();
  lp2->Draw();
  
  c11->SaveAs("spectra_rebin1stbin.eps");


}

 
  
  
  

TCanvas * DrawCanvas(TString opt="c")
{
  TCanvas *c1 = new TCanvas(opt.Data(),opt.Data(),10,10,600,600);
  c1->cd(1);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.15);
  return c1;
}
TPad *MyPad(Double_t x1 = 0.01,Double_t y1 = 0.01,Double_t x2 = 0.49, Double_t y2= 0.49, Float_t lm  =0.2,Float_t rm = 0.05,Float_t tm =0.05,Float_t bm = 0.2,TString name = "pad")
{
  TPad *c1_1 = new TPad(name.Data(), name.Data(),x1,y1,x2,y2);
  c1_1->Draw();
  c1_1->cd();
  c1_1->cd()->SetLogy(); 
  c1_1->Range(0,0,1,1);
  c1_1->SetBorderSize(2);
  c1_1->SetBorderMode(0);
  c1_1->SetFillColor(10);
  c1_1->SetFrameFillColor(10);
  c1_1->SetFrameLineWidth(2);
  c1_1->SetLeftMargin(lm);
  c1_1->SetRightMargin(rm);
  c1_1->SetTopMargin(tm);
  c1_1->SetBottomMargin(bm);
  c1_1->SetTicks(1,1);
  return c1_1;
} 





TH1D * Plothisto(TH1D *h1, Int_t MarkerColor = 1, Int_t MarkerStyle = 20)
{
  
  h1->SetTitle("");
  h1->SetMarkerStyle(MarkerStyle);
  h1->SetMarkerColor(MarkerColor);
  h1->SetLineColor(MarkerColor);
  
  h1->SetMarkerSize(1.4);
  h1->SetLineWidth(2);
  h1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  h1->GetXaxis()->CenterTitle(false);
  h1->GetXaxis()->SetNdivisions(506);
  h1->GetYaxis()->SetNdivisions(505);
  h1->GetXaxis()->SetLabelOffset(0.015);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetLabelSize(0.04);
  h1->GetXaxis()->SetTitleSize(0.04);
  h1->GetXaxis()->SetTickLength(0.04);
  h1->GetXaxis()->SetTitleOffset(1.2);
  h1->GetYaxis()->SetTitleOffset(1.7);
  h1->GetYaxis()->CenterTitle(true);
  h1->GetYaxis()->SetDecimals(false);
  //h1->GetYaxis()->SetNdivisions(310);
  h1->GetYaxis()->SetLabelOffset(0.015);
  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetLabelSize(0.04);
  h1->GetYaxis()->SetTickLength(0.04);
  h1->GetYaxis()->SetTitleSize(0.04);
  h1->GetYaxis()->SetTitle("(1/#it{N}_{evt})(d^{2}#it{N}/d#it{p}_{T}d#it{y}) (GeV/#it{c})^{-1}");
  h1->GetYaxis()->SetTitleFont(42);
  return h1;
}




TGraphAsymmErrors * PlotGraph(Int_t NdataPoint, Int_t MarkerColor = 1, Int_t MarkerStyle = 20, Double_t *X, Double_t *ErXL, Double_t *ErXH, Double_t *Y, Double_t *ErYL, Double_t *ErYH)
{
  TGraphAsymmErrors *gr = new TGraphAsymmErrors(NdataPoint, X ,Y , ErXL, ErXH, ErYL, ErYH);
  gr->SetTitle("");
  gr->SetMarkerStyle(MarkerStyle);
  gr->SetMarkerColor(MarkerColor);
  gr->SetMarkerSize(1.8);
  gr->SetLineWidth(2);
  gr->SetLineColor(MarkerColor);
  gr->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  gr->GetXaxis()->CenterTitle(false);
  gr->GetXaxis()->SetNdivisions(506);
  gr->GetYaxis()->SetNdivisions(505);
  gr->GetXaxis()->SetLabelOffset(0.015);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetTickLength(0.04);
  gr->GetXaxis()->SetTitleOffset(1.25);
  gr->GetYaxis()->SetTitleOffset(1.45);
  gr->GetYaxis()->CenterTitle(false);
  gr->GetYaxis()->SetDecimals(false);
  //gr->GetYaxis()->SetNdivisions(310);
  gr->GetYaxis()->SetLabelOffset(0.015);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetTickLength(0.04);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitle("#it{#rho}_{00}"); 
  gr->GetYaxis()->SetTitleFont(42);
  return gr;
}




//==================================
TLegend *DrawLegend(Double_t x1,Double_t y1,Double_t x2,Double_t y2)
{
  //TLegend *legend = new TLegend(0.5,0.65,0.88,0.85);
  TLegend *legend = new TLegend(x1,y1,x2,y2);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->SetLineColor(0);
  legend->SetShadowColor(0);
  //legend->AddEntry(gr1,"(0 - 100) %","p");
  //legend->AddEntry(func1,"p_{0}[ 1 + 2 v_{2}^{obs} cos(2(#Phi - #Psi))]","l");
  return legend;
}
//=====================================
TLatex *DrawText(Double_t x = 0, Double_t y = 0,Int_t tColor = 2,TString name)
{
  TLatex* tex = new TLatex(x,y,name.Data());
  tex->SetTextSize(0.04);
  tex->SetTextColor(tColor);
  tex->SetTextFont(42);
  //tex->Draw();
  return tex;
}
//==================================
void MyStyle()
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
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
  gStyle->SetLineWidth(2);
  //gStyle->SetErrorX(0.2);

}

TGraphErrors *DrawFrame(TGraphErrors *h, Int_t MCol,Int_t MSty, Bool_t mrk = 0)
{
  //h->GetXaxis()->SetTitle("d#it{N}_{ch}/d#eta");
  h->GetXaxis()->SetTitle("<#it{N}_{part}>");
  if(mrk){
    h->SetMarkerColor(MCol);
    h->SetMarkerStyle(MSty);
    h->SetMarkerSize(1.4);
  }
  h->SetLineColor(MCol);
  h->SetLineWidth(2);
  h->GetXaxis()->CenterTitle(false);
  h->GetXaxis()->SetNdivisions(505);
  h->GetYaxis()->SetNdivisions(505);
  h->GetXaxis()->SetLabelOffset(0.017);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTickLength(0.06);
  h->GetXaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetYaxis()->CenterTitle(false);
  h->GetYaxis()->SetDecimals(false);
  //h->GetYaxis()->SetNdivisions(310);
  h->GetYaxis()->SetLabelOffset(0.015);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTickLength(0.04);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitle("#it{#rho_{00}}"); 
  h->GetYaxis()->SetTitleFont(42);
  //h->Draw("AP");
  return h;
}







