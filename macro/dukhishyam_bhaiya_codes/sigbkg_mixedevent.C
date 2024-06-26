Int_t fReBinSpectra = 1;
void sigbkg_mixedevent()
{
  MyStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
   gStyle->SetOptTitle(0);

   //   TFile *fInputFile2 = new TFile("Output_K0sK0s_0to100_ptbin0to20.root");

    TFile *fInputFile2 = new TFile("Output_K0sK0s_0to100_ptbin1to10.root");
    fInputFile2->ls();
    Int_t STCol = 1;
    Int_t STSty = 24;
    Int_t STCol1 = 2;
    Int_t STSty1 = 24;
   
    Int_t CBCol = 4;
    Int_t CBSty = 24;

    int maxs=170E+7;
    int mins=0.01E+9;
    
    Int_t pTBin1 = 7;
    
    //TString ptbn ="1.2 #leq #it{p}_{T} < 1.4 GeV/#it{c}";
    
    TString ptbn ="1 #leq #it{p}_{T} < 10 GeV/#it{c}";
  //TString ptbn ="0.6 #leq #it{p}_{T} < 0.8 GeV/#it{c}";
  //TString cos ="0.6 #leq cos#it{#theta*} < 0.8";

   TLegend *l2 = DrawLegend(0.28,0.5,0.4,0.7);
   TGaxis::SetMaxDigits(3);
   //  TH1D * hScCBkg  = (TH1D*)fInputFile2->Get("hbkg");
   TH1D * hScCBkg  = (TH1D*)fInputFile2->Get("hbkgmix");

  //hScCBkg->Draw();
  TH1D * hSigTot1  = (TH1D*)fInputFile2->Get("hsig");
    
  
  ///hSigTot1->Draw();

  // hScCBkg->Draw("");
  //return 0;


    // hSigTot1->SetTitle("");
  
  TCanvas *c1 = DrawCanvas("c11");
  //c11->Range(0,0,1,1);
  c11->cd();
  //TPad *p3 = MyPad(0.01, 0.01, 0.99, 0.99, 0.17, 0.02, 0.06, 0.12, "pad3");
  //TH1F *hf1 = DrawFrame(hSigTot1,STCol,STSty,0);
  //TH1F *hf2 = DrawFrame(hScCBkg,STCol1,STSty1,0);
  //hf1->GetXaxis()->SetRangeUser(0.6,1.5);
  // hf1->SetMaximum(maxs);
  //hf1->SetMinimum(mins);

  // hSigTot1->GetXaxis()->SetRangeUser(0.6,1.2);
  // hSigTot1->Draw("");
  // hScCBkg->Draw("");
  // hf1->Draw();
  // hf2->Draw("same");
   // hSigTot1->Draw();
   //hScCBkg->Draw("same");

TPad *p3 = MyPad(0.01, 0.01, 0.99, 0.99, 0.15, 0.06, 0.06, 0.15, "pad3");
  // TH1F *hf1 = DrawFrame(hSigTot1,STCol,STSty,1);
  // TH1F *hf2 = DrawFrame(hScCBkg,STCol1,STSty1,1);
  // hf1->GetXaxis()->SetRangeUser(lowInv,highInv);
   hSigTot1 = DrawFrame(hSigTot1,STCol,STSty,1);
   hScCBkg = DrawFrame(hScCBkg,STCol1,STSty1,1);
   // hf1->SetMaximum(maxs);
   // hf1->SetMinimum(mins);
   //hf1->Draw();
   // hf2->Draw("same");
   //   hSigTot1->SetMaximum(maxs);
   // hSigTot1->SetMinimum(mins);

    hSigTot1->GetXaxis()->SetRangeUser(0.95,2.5);
    //    hSigTot1->GetXaxis()->SetTitle("");
    //    hSigTot1->GetYaxis()->SetRangeUser(0.1e+5,0.7e+5);

    // hSigTot1->GetYaxis()->SetRangeUser(0.1e+6,0.28e+6);

    //    hScCBkg->Scale(1.0/10.0);
    
    hScCBkg->SetMarkerSize(0.9);
    hSigTot1->SetMarkerSize(0.9);
    
    hSigTot1->Draw("");
    hScCBkg->Draw("same");
    
   TLegend *lp2 = DrawLegend(0.51,0.62,0.7,0.92);
  lp2->SetTextSize(.05);
  lp2->AddEntry((TObject*)0,"#bf{ALICE Performance}","");
  lp2->AddEntry((TObject*)0,"pp, #sqrt{#it{s}} = 13 TeV ","");
  lp2->AddEntry((TObject*)0,"V0M (0-100%)","");
  lp2->AddEntry((TObject*)0,"|#it{y}| < 0.5","");
  lp2->AddEntry((TObject*)0,Form("%s",ptbn.Data()),"");
  lp2->Draw("same");

  TLegend* lp3=DrawLegend(0.4,0.2,0.64,0.34);
  lp3->AddEntry(hSigTot1,"Same-event pairs","p");
  lp3->AddEntry(hScCBkg,"Mixed-event pairs","p");
  lp3->SetTextSize(.05);
  lp3->Draw("same");

  p3->Modified();
  c11->cd();

  //  c11->SaveAs("Mixlow_InvMassNormIntermediatePt_pPb5TeV_bkg.gif");
  // c11->Print("Mixlow_InvMassNormIntermediatePt_pPb5TeV_bkg.eps");
  // c11->Print("Mixlow_InvMassNormIntermediatePt_pPb5TeV_bkg.root");
  // c11->Print("Mixlow_sigbkg_pPb5TeV_bkg.pdf");

  /*
    c11->SaveAs("Rotational_low_InvMassNormIntermediatePt_pPb5TeV_bkg.gif");
  c11->Print("Rotational_low_InvMassNormIntermediatePt_pPb5TeV_bkg.eps");
  c11->Print("Rotational_low_InvMassNormIntermediatePt_pPb5TeV_bkg.root");
  c11->Print("Rotational_low_sigbkg_pPb5TeV_bkg.pdf");
  */

  // c11->SaveAs("Preliminary_Invmass_K0sK0s_0to100_pt1to10_mixedevent.png");
  //  c11->Print("Preliminary_Invmass_K0sK0s_0to100_pt1to10_mixedevent.eps");


  c11->SaveAs("Performance_Invmass_K0sK0s_0to100_pt1to10_mixedevent.png");
  c11->Print("Performance_Invmass_K0sK0s_0to100_pt1to10_mixedevent.eps");

  
  //  c11->Print("Mi.root");


}

TCanvas * DrawCanvas(TString opt="c")
{
  TCanvas *c1 = new TCanvas(opt.Data(),opt.Data(),10,10,700,600);
  c1->cd(1);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.2);
  return c1;
}
TPad *MyPad(Double_t x1 = 0.01,Double_t y1 = 0.01,Double_t x2 = 0.49, Double_t y2= 0.49, Float_t lm  =0.2,Float_t rm = 0.05,Float_t tm =0.05,Float_t bm = 0.2,TString name = "pad")
{
  TPad *c1_1 = new TPad(name.Data(), name.Data(),x1,y1,x2,y2);
  c1_1->Draw();
  c1_1->cd();
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
TGraphErrors * PlotGraph(Int_t NdataPoint, Int_t MarkerColor = 1, Int_t MarkerStyle = 20, Double_t *X, Double_t *ErX,Double_t *Y, Double_t *ErY)
{
  TGraphErrors *gr = new TGraphErrors(NdataPoint, X ,Y , ErX, ErY);
  gr->SetTitle("");
  gr->SetMarkerStyle(MarkerStyle);
  gr->SetMarkerColor(MarkerColor);
  gr->SetMarkerSize(1.5);
  //gr->GetXaxis()->SetTitle("m_{T} - m_{0} (GeV/c^{2})");
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetTitleOffset(1.25);
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetXaxis()->CenterTitle(true);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetNdivisions(509);
  //gr->GetYaxis()->SetTitle("#frac{1}{2#pi m_{T}N_{Ev}}#frac{dN^{2}}{dydm_{T}} (GeV/c)^{-1}");
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetYaxis()->CenterTitle(true);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleOffset(1.7);
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetNdivisions(509);
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
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
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
}
Double_t VoigtPoly2(Double_t *x, Double_t *par)
{
  return 0.001*fReBinSpectra*par[0]*TMath::Voigt(x[0] - par[2], par[3], par[1], 4) + (par[4] + par[5]*x[0] + par[6]*x[0]*x[0]);
}
Double_t Voigt(Double_t *x, Double_t *par)
{
  return 0.001*fReBinSpectra*par[0]*TMath::Voigt(x[0] - par[2], par[3], par[1], 4) ;
}
Double_t Poly2(Double_t *x, Double_t *par)
{
  return  par[0] + par[1]*x[0]+par[2]*x[0]*x[0];
}
Double_t BWPoly2(Double_t *x, Double_t *par)
{
  return (par[0]*par[1]*0.01*fReBinSpectra)/(2*3.14159)/((x[0]-par[2])*(x[0]-par[2])+(par[1]*par[1])/4.)+ par[3]+ par[4]*x[0]+par[5]*x[0]*x[0];
}
Double_t BW(Double_t *x, Double_t *par)
{
  return (1.0*0.01*fReBinSpectra/(2*3.14159))*((par[0]*par[1])/((x[0]-par[2])*(x[0]-par[2]) +(par[1]/2)*(par[1]/2)));
}
void *DrawFrame(TH1 *h, Int_t MCol,Int_t MSty, Bool_t mrk = 0)
{
  h->GetXaxis()->SetTitle("#it{M}_{K^{0}_{S}K^{0}_{S}} (GeV/#it{c^{2}})");
  if(mrk){
    h->SetMarkerColor(MCol);
    h->SetMarkerStyle(MSty);
    h->SetMarkerSize(1.);
  }
  h->SetLineColor(MCol);
  h->SetLineWidth(2);
  //h->GetXaxis()->CenterTitle(true);
  h->GetXaxis()->SetNdivisions(509);
  h->GetYaxis()->SetNdivisions(508);
  h->GetXaxis()->SetLabelOffset(0.015);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitleSize(0.06);
  //  h->GetXaxis()->SetTickLength(0.03);
  h->GetXaxis()->SetTitleOffset(1.15);
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetXaxis()->CenterTitle(true);
  h->GetYaxis()->CenterTitle(true);
  h->GetYaxis()->SetDecimals(false);
  //h->GetYaxis()->SetNdivisions(310);
  h->GetYaxis()->SetLabelOffset(0.015);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTickLength(0.02);
  h->GetXaxis()->SetTickLength(0.02);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitle("Counts / (0.01 GeV/#it{c^{2}})"); 
  h->GetYaxis()->SetTitleFont(42);
  h->Draw("");
  return h;
}
