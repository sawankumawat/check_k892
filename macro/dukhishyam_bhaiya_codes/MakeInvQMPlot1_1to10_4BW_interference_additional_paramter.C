Int_t fReBinSpectra = 1;
void MakeInvQMPlot1_1to10_4BW_interference_additional_paramter()
{
  MyStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
 //TFile *fInputFile1 = new TFile("Output_K0sK0s_0to100_ptbin0p5to10.root");

  //  TFile *fInputFile1 = new TFile("K0sK0s_signal_pTrange_0p5to10_0to100.root");
  TFile *fInputFile1 = new TFile("K0sK0s_signal_0to100_1to10GeV.root");

  fInputFile1->ls();

  //  return 0;
  //  TH1D * hfsigbkg =  (TH1D*)fInputFile1->Get("fHistNum0");//hItBgOneBw  hMassPt1TC
  TH1D * hfsigbkg =  (TH1D*)fInputFile1->Get("fHistNum0");//hItBgOneBw  hMassPt1TC

  TFitResultPtr r;
   
  //TFile *fInputFile1 = new TFile("AnalysisResultsKStar_KstarPM_Rotational_Invmass_plot_signal_rot.root");
  fInputFile1->ls();
   //   return 0;

  Int_t STCol = 1;
  Int_t STSty = 21;
  
  Int_t CBCol = 1;
  Int_t CBSty = 21;

  Float_t lowInv = 1.1;
  Float_t highInv = 2.0;

  Float_t maxst =820E3;
  Float_t minst = .05E6; 

  Float_t maxs = 2.9E+6;
  Float_t mins = 0.4E+6; 

  //Double_t fitlow = 0.66;
  //Double_t fithi = 1.1;
  
  //  Double_t fitlow = 1.23;
  // Double_t fithi = 1.95;
  
  Double_t fitlow = 1.26;
  Double_t fithi = 1.96;

  // Double_t fitlow = 1.38;
  // Double_t fithi = 1.88;
  

  //  Int_t pTBin1 = 12;
  //TString ptbn ="1.2 #leq #it{p}_{T} < 1.8 GeV/#it{c}";
  TString ptbn ="1.0 #leq #it{p}_{T} < 10 GeV/#it{c}";
  
  //TString ptbn ="4.0 #leq #it{p}_{T} < 4.5 GeV/#it{c}";
  //TString ptbn ="2.5 #leq #it{p}_{T} < 3.0 GeV/#it{c}";
  //TString ptbn ="2.0 #leq #it{p}_{T} < 2.5 GeV/#it{c}";

  //TString ptbn ="0.6 #leq #it{p}_{T} < 0.8 GeV/#it{c}";
  //TString ptbn ="0.6 #leq #it{p}_{T} < 0.8 GeV/#it{c}";
  //TString cos ="0.6 #leq cos#it{#theta*} < 0.8";
  
  TLegend *l2 = DrawLegend(0.25,0.05,0.4,0.15);
  
  //TGaxis::SetMaxDigits(3);
   
  TGaxis::SetMaxDigits(2);
   
  /*  TF1 *fBgOnly1 = new TF1("fBgOnly1",Poly2,fitlow,fithi,3);
      TF1 *fSigBg1 = new TF1("fSigBg1",BWPoly2,fitlow,fithi,6);
      TF1 *fSigOnly1 = new TF1("fSigOnly1",BW,fitlow,fithi,3);
  */
   
   //   TF1 *fBgOnly1 = new TF1("fBgOnly1",Expol,fitlow,fithi,4);
   // TF1 *fSigBg1 = new TF1("fSigBg1",BreitWignerExpol,fitlow,fithi,7);
   // TF1 *fSigOnly1 = new TF1("fSigOnly1",BW,fitlow,fithi,3);

  //   TF1 *fSigBg1 = new TF1("fitfunc",rBreitWigner,fitlow,fithi,17);
   TF1 *fSigBg1 = new TF1("fitfunc",rBreitWigner,fitlow,fithi,20);
   TF1 *fBgOnly1 = new TF1("fitfunc1",polynomial2,fitlow,fithi,3);

   //TF1 *fSigOnly1 = new TF1("fitFcn2",rBW,fitlow,fithi,14);

   TF1 *fSigOnly1 = new TF1("fitFcn2",rBW,fitlow,fithi,17);



   
   
   
   // TF1 *fSig1 = new TF1("fsig1",rBW1,1.1,1.4,3);
   //TF1 *fSig2 = new TF1("fsig2",rBW1,1.4,1.6,3);
   // TF1 *fSig3 = new TF1("fsig3",rBW1,1.6,1.8,3);
   
   
   //    TF1 *fSig1 = new TF1("fsig1",BreitWigner,1.1,1.4,3);
   // TF1 *fSig2 = new TF1("fsig2",BreitWigner,1.4,1.6,3);
   //TF1 *fSig3 = new TF1("fsig3",BreitWigner,1.6,1.8,3);

   
   TF1 *fSig1 = new TF1("fsig1",rBW1,1.1,1.4,3);
   TF1 *fSig2 = new TF1("fsig2",rBW1,1.4,1.6,3);
   TF1 *fSig3 = new TF1("fsig3",rBW1,1.6,1.8,3);
    

  
   
   fSig1->SetLineColor(1);
   fSig1->SetLineStyle(1);
   fSig1->SetLineWidth(2);

   fSig2->SetLineColor(4);
   fSig2->SetLineStyle(2);
   fSig2->SetLineWidth(2);
   
   fSig3->SetLineColor(3);
   fSig3->SetLineStyle(1);
   fSig3->SetLineWidth(2);
   
   
   
   fBgOnly1->SetLineColor(4);
   fBgOnly1->SetLineStyle(2);
   fBgOnly1->SetLineWidth(2);
   
   fSigBg1->SetLineColor(2);
   fSigBg1->SetLineStyle(1);
   fSigBg1->SetLineWidth(2);

   fSigOnly1->SetLineColor(6);
   fSigOnly1->SetLineStyle(2);
   fSigOnly1->SetLineWidth(2);
   
   fSigBg1->SetParameter(0,1000); //yield
   //   fSigBg1->SetParLimits(1,1.15,1.3); //inv. mass peak range
   //fSigBg1->SetParLimits(2,0.03,0.1); //width-fixed
   // fSigBg1->SetParameter(1,1.270); //inv. mass peak range
   // fSigBg1->SetParameter(2,0.186); //width-fixed

    fSigBg1->SetParLimits(1,1.1,1.4); //inv. mass peak range                                                                                                          
    fSigBg1->SetParLimits(2,0.12,0.16); //width-fixed                                                                                                                 

   fSigBg1->SetParameter(3,1000); //yield
   //fSigBg1->SetParLimits(4,1.15,1.4); //inv. mass peak range
   // fSigBg1->SetParLimits(5,0.07,0.1); //width-fixed
   
   fSigBg1->SetParameter(4,1.320); //inv. mass peak range
   fSigBg1->SetParameter(5,0.1); //width-fixed
   
  
   fSigBg1->SetParameter(6,10000); //yield
   fSigBg1->SetParLimits(7,1.45,1.65); //inv. mass peak range
   fSigBg1->SetParLimits(8,0.03,0.16); //width-fixed

   //   fSigBg1->SetParameter(7,1.525); //inv. mass peak range
   // fSigBg1->SetParameter(8,0.073); //width-fixed
   
   fSigBg1->SetParameter(9,1000); //yield
   fSigBg1->SetParLimits(10,1.6,1.8); //inv. mass peak range
   fSigBg1->SetParLimits(11,0.08,0.171); //width-fixed

   //   fSigBg1->SetParameter(10,1.710); //inv. mass peak range
   // fSigBg1->SetParameter(11,0.139); //width-fixed

   //   fSigBg1->SetParameter(12,0.08,); //width-fixed
   // fSigBg1->SetParameter(13,0.08,0.17); //width-fixed
  
  fSigBg1->SetParNames("N_{1270}","M_{1270}","#Gamma_{1270}"," N_{1320}","M_{1320}","#Gamma_{1320}","N_{1525}","M_{1525}","#Gamma_{1525}");
  
  TH1D * hSigTemp1 =  (TH1D*)hfsigbkg->Clone("hfsig_f0");//hItBgOneBw  hMassPt1TC

  TH1D * hSigTemp2 =  (TH1D*)hfsigbkg->Clone("hfsig_f2");//hItBgOneBw  hMassPt1TC                                                        

  TH1D * hSigFit1  =  (TH1D*)hSigTemp1->Clone("hSigFit1");
  hSigFit1->Draw();

  //  return 0;

  
  
  hSigFit1->SetTitle("");
  hSigFit1->SetMarkerSize(1.0);
  fSigBg1->SetNpx(10000);
  
  r=hSigFit1->Fit(fSigBg1,"REBMS+"); //signal after like subtraction likesign  bkg 
  
  //hSigFit1->Fit(fSigBg1,"REBMS","",fitlow,fithi);

  /*
  Int_t fitres1;
  Int_t trials1 = 0;
  trials1 = 0;
  do{
    // fitres1 = hSigFit1->Fit(fSigBg1, "REBMS");
    Double_t *Partemp   = fSigBg1->GetParameters();
    fSigBg1->SetParameters(&Partemp[0]);
    Printf("Trial: %d %d", trials1++,fitres1);
    if(trials1 > 1) {
      Printf("FIT DOES NOT CONVERGE IN LINE %d",__LINE__);
      break;
    }
  }
  while (fitres1 != 0);
  */
  
  Double_t *Par1   = fSigBg1->GetParameters();
  Double_t *Erpar1 = fSigBg1->GetParErrors();

  Double_t yield, mass, width, yielderr, masserr,widtherr,mass2, width2, yielderr, masserr2,widtherr2,mass3, width3, yielderr3, masserr3,widtherr3;
  yield = fSigBg1->GetParameter(0);
  yielderr = fSigBg1->GetParError(0);

  mass = fSigBg1->GetParameter(1);
  masserr = fSigBg1->GetParError(1);


  mass2 = fSigBg1->GetParameter(4);
  masserr2 = fSigBg1->GetParError(4);
  
  mass3 = fSigBg1->GetParameter(7);
  masserr3 = fSigBg1->GetParError(7);
  

  width = 1000*fSigBg1->GetParameter(2);
  widtherr = 1000*fSigBg1->GetParError(2);
  
  width2 = 1000*fSigBg1->GetParameter(5);
  widtherr2 = 1000*fSigBg1->GetParError(5);

    
  width3 = 1000*fSigBg1->GetParameter(8);
  widtherr3 = 1000*fSigBg1->GetParError(8);
  
  cout << mass2 << "\t-----width" << width3 << endl;
  
  //  cout << "Yield--before residual subtraction--" << fSigBg1->GetParameter(0) << endl; 

  // fSigOnly1->Draw();
  
  
  fSigOnly1->SetParameters(&Par1[0]);
  fSigOnly1->SetParErrors(&Erpar1[0]);

  //cout << "Yield--after residual subtraction--" << fSigOnly1->GetParameter(0) << endl; 
  
 
  //return 0;
  
  cout << "\t-----significance_num---" << (fSigOnly1->Integral(1.2,1.4)/hSigFit1->GetBinWidth(1)) << "\t" << (fSigOnly1->Integral(1.3,1.75)/hSigFit1->GetBinWidth(1)) << "\t" <<(fSigOnly1->Integral(1.6,2.0)/hSigFit1->GetBinWidth(1))<<endl;

  
  //Set fiiting function range
   
  
   

  
  /*
  fSigOnly1->SetParameter(0,yield);
  fSigOnly1->SetParError(0,yielderr);

  fSigOnly1->SetParameter(0,mass);
  fSigOnly1->SetParError(0,masserr);
  
  fSigOnly1->SetParameter(0,width);
  fSigOnly1->SetParError(0,widtherr);
  */
  
  //  fSigOnly1->Draw();
  // return 0;
  
  fBgOnly1->SetParameters(&Par1[17]);
  fBgOnly1->SetParErrors(&Erpar1[17]);

  // fBgOnly1->Draw();
  // return 0;

  
  TCanvas *c1 = DrawCanvas("c11");
  //c11->Range(0,0,1,1);
  c11->cd();
  c11->SetLeftMargin(0.5);

  //hSigFit1->Draw("");
  //fSigBg1->Draw("same");
  //fBgOnly1->Draw("same");


  /*
  
  TH1D * hSigFit1_SResBkg  =  (TH1D*)hSigTemp1->Clone("hSigFit1");

  TH1D * hSig=  (TH1D*)hSigTemp1->Clone("hSigFit1_new");

  hSigFit1_SResBkg->Reset();

  for(Int_t ip =0; ip < hSigFit1_SResBkg->GetNbinsX(); ip++)
    {
      float x=hSigFit1_SResBkg->GetBinCenter(ip);
      float y=fBgOnly1->Eval(x);
      // value = hSigTemp1->GetBinContent(ip+1) - (fBgOnly1->Integral(ip+1,ip+2)/hSigFit1_SResBkg->GetBinWidth(2)); 
      
      cout << "signal-----" <<  y << endl;
      
      hSigFit1_SResBkg->SetBinContent(ip+1,y);


      
      //     hSigFit1_SResBkg->SetBinError(ip+1,value_err);
    }

  //hSig->Draw();
  // return 0;

  
  //hSig->Add(,-1);
  TCanvas *c1 = new TCanvas("c1","",500,500);
  c1->cd();
  
 hSig->Draw("p");
 hSigFit1_SResBkg->Draw("p same");
  return 0;
  */

  
  //return 0;
  TPad *p3 = MyPad(0.01, 0.01, 0.99, 0.99, 0.14, 0.02, 0.06, 0.12, "pad3");
  p3->SetLeftMargin(0.16);
  /*  TH1 *hf1 = DrawFrame(hSigFit1,STCol,STSty,0);
   hf1->GetXaxis()->SetRangeUser(lowInv,highInv);
  hf1->SetMaximum(1.3);
  hf1->SetMinimum(0.65);
  hf1->Draw("");
  */
  hSigFit1 = DrawFrame(hSigFit1,STCol,STSty,1);
  hSigFit1->GetXaxis()->SetRangeUser(1.1,2.2);
  hSigFit1->SetMarkerSize(0.8);

  // hSigFit1->GetXaxis()->SetRangeUser(1.1,1.9);
  hSigFit1->GetYaxis()->SetRangeUser(-1000.0,0.18e+6);

  //hSigFit1->SetMinimum(-5000.00);
  //hSigFit1->SetMaximum(0.073e+6);
  //hSigFit1->SetMaximum(0.046e+6);
  //fSigOnly->SetFillColor(2);

  
  hSigFit1->Draw("ep");
  fSigOnly1->Draw("same e3");
  hSigFit1->Draw("same l"); 
  fSigBg1->Draw("same l");
  fBgOnly1->Draw("same l");
  
  /*  fSig1->SetParameter(0,fSigOnly1->GetParameter(0));
  fSig1->SetParError(0,fSigOnly1->GetParameter(0));
  fSig1->SetParameter(1,fSigOnly1->GetParameter(1));
  fSig1->SetParError(1,fSigOnly1->GetParameter(1));
  fSig1->SetParameter(2,fSigOnly1->GetParameter(2));
  fSig1->SetParError(2,fSigOnly1->GetParameter(2));
  
  
  fSig2->SetParameter(0,fSigOnly1->GetParameter(6));
  fSig2->SetParameter(1,fSigOnly1->GetParameter(7));
  fSig2->SetParameter(2,fSigOnly1->GetParameter(8));
  
  fSig3->SetParameter(0,fSigOnly1->GetParameter(9));
  fSig3->SetParameter(1,fSigOnly1->GetParameter(10));
  fSig3->SetParameter(2,fSigOnly1->GetParameter(11));
  */

  // fSig1->Draw("same l");
  // fSig2->Draw("same l");
  // fSig3->Draw("same l");
  
  // hSigFit1->SetMarkerColor(1);
  // hSigFit1->SetMarkerStyle(20);
   
  // TLegend *lp2 = DrawLegend(0.38,0.67,0.56,0.9);

  // TLegend *lp2 = DrawLegend(0.3,0.67,0.5,0.9);
  //rotational 
 TLegend *lp2 = DrawLegend(0.28,0.8,0.52,0.9);
 lp2->SetTextSize(.03);
 // lp2->AddEntry((TObject*)0,"#bf{ALICE Preliminary}","");
  //  lp2->AddEntry((TObject*)0,"V0A Multiplicity Event Classes","");
  lp2->AddEntry((TObject*)0,"pp #sqrt{#it{s}} = 13 TeV (0-100%)","");
  lp2->AddEntry((TObject*)0,"|#it{y}| < 0.5","");

  //lp2->AddEntry((TObject*)0,"K^{*0} (#bar{K^{*0}}) #rightarrow K^{+}#pi^{-} (K^{-}#pi^{+})","");
  //lp2->AddEntry((TObject*)0,"K^{*0} (#bar{K^{*0}}) #rightarrow K^{+}#pi^{-} (K^{-}#pi^{+})","");
  //lp2->AddEntry((TObject*)0,"K*^{0} (#bar{K}*^{0}) #rightarrow K^{+}#pi^{-} (K^{-}#pi^{+})","");
  //lp2->AddEntry((TObject*)0,"K*^{#pm} #rightarrow K^{0}_{S} #pi^{#pm}","");

   //lp2->AddEntry((TObject*)0,"|#it{y}| < 0.5","");
  // lp2->AddEntry((TObject*)0,"-0.5 < #it{y} < 0","");
   lp2->Draw("same");

   //  TLegend* lp3=DrawLegend(0.19,0.17,0.40,0.42);
   //  TLegend* lp3=DrawLegend(0.19,0.17,0.40,0.42);
   //  TLegend* lp3=DrawLegend(0.25,0.2,0.40,0.35);

   //rotational 
   TLegend* lp3=DrawLegend(0.5,0.6,0.8,0.8);
  lp3->AddEntry((TObject*)0,Form("%s",ptbn.Data()),"");
  //  lp3->AddEntry((TObject*)0,Form("%s",cos.Data()),"");
  lp3->AddEntry(hSigFit1,"Data","p");
  lp3->AddEntry(fSigBg1,"4rBW + Residual BG","l");
  lp3->AddEntry(fBgOnly1,"Residual BG","l");
  lp3->SetTextSize(.04);

  // lp3->AddEntry(hSigFit1,"Data (stat. uncert.)","p");
  //lp3->AddEntry(fSigOnly1,"Breit-Wigner Peak Fit","l");
  //lp3->AddEntry(fBgOnly1,"Residual BG","l");
  //lp3->SetTextSize(.04);
  lp3->Draw("same");

  p3->Modified();
  c11->cd();

  //  TCanvas *c22 = new TCanvas("c22","c22",500,500);
  //c22->cd();

  TCanvas *c22 = DrawCanvas("c22");
  c22->cd();
  c22->SetLeftMargin(0.1);
  MyStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TH1D * hSigFit1_SResBkg  =  (TH1D*)hSigTemp1->Clone("hSigFit1new");
  TH1D * hSig=  (TH1D*)hSigTemp2->Clone("hSigFit1_new");                                                                                                          
    hSigFit1_SResBkg->Reset();                                                                                                                                       
                                                                                                                                                                      
  for(Int_t ip =0; ip < hSigFit1_SResBkg->GetNbinsX(); ip++)                                                                                                          
    {                                                                                                                                                                 
      float x=hSigFit1_SResBkg->GetBinCenter(ip+1);
      // float x=hSigFit1_SResBkg->GetBinLowEdge(ip+1);
      float y=fBgOnly1->Eval(x);                                                                                                                                            //value = hSigTemp1->GetBinContent(ip+1) - (fBgOnly1->Integral(ip+1,ip+2)/hSigFit1_SResBkg->GetBinWidth(2));
      cout << "\----ip " <<  ip <<  "---position__x" << x  <<  "bkg-----" <<  y << "\t signal------" << hSig->GetBinContent(ip+1) << endl;
      hSigFit1_SResBkg->SetBinContent(ip+1,y);        
    }

  //  TF1 *fSigOnly2 = new TF1("fitFcn2",rBW,fitlow,fithi,14);

  //  TF1 *fSigOnly2 = new TF1("fitFcn2",rBW,1.24,1.88,14);
  
  TF1 *fSigOnly2 = new TF1("fitFcn2",rBW,1.26,1.97,17);
  fSigOnly2->SetParameter(0,1000); //yield
  //  fSigOnly2->SetParameter(1,1.270); //inv. mass peak range
  // fSigOnly2->SetParameter(2,0.186); //width-fixed
   fSigOnly2->SetParLimits(1,1.1,1.4); //inv. mass peak range
   fSigOnly2->SetParLimits(2,0.12,0.16); //width-fixed
   
   fSigOnly2->SetParameter(3,1000); //yield
   fSigOnly2->SetParameter(4,1.320); //inv. mass peak range
   fSigOnly2->SetParameter(5,0.1); //width-fixed
   
  
   fSigOnly2->SetParameter(6,10000); //yield
   fSigOnly2->SetParLimits(7,1.45,1.65); //inv. mass peak range
   fSigOnly2->SetParLimits(8,0.03,0.16); //width-fixed

   fSigOnly2->SetParameter(9,1000); //yield
   fSigOnly2->SetParLimits(10,1.6,1.8); //inv. mass peak range
   fSigOnly2->SetParLimits(11,0.08,0.171); //width-fixed
  
   hSig->Add(hSigFit1_SResBkg,-1);
   TFile *fout = new TFile("K0sK0s_signal_pp13TeV_ptrange_1to10_4rBW_new.root","recreate");
   fout->cd();
   hSig->Write();
   
   r=hSig->Fit(fSigOnly2,"REBMS+");

   hSig = DrawFrame(hSig,STCol,STSty,1);
   hSig->GetXaxis()->SetRangeUser(1.1,2.2);
   hSig->SetMarkerSize(0.8);
   hSig->GetYaxis()->SetRangeUser(-100.0,20000.0);
   //hSigFit1_SResBkg->Draw();
   hSig->Draw("p");
   hSig->SetTitle("");
   hSig->SetMarkerSize(1.0);
   // fSig->SetNpx(10000);
   hSig->GetXaxis()->SetRangeUser(1.2,2.0);
   //   fSigOnly1->Draw("p same"); 

 TLegend *lp2 = DrawLegend(0.28,0.6,0.52,0.8);
 lp2->SetTextSize(.03);
 // lp2->AddEntry((TObject*)0,"#bf{ALICE Preliminary}","");
  //  lp2->AddEntry((TObject*)0,"V0A Multiplicity Event Classes","");
  lp2->AddEntry((TObject*)0,"pp #sqrt{#it{s}} = 13 TeV (0-100%)","");
  lp2->AddEntry((TObject*)0,"|#it{y}| < 0.5","");

  //lp2->AddEntry((TObject*)0,"K^{*0} (#bar{K^{*0}}) #rightarrow K^{+}#pi^{-} (K^{-}#pi^{+})","");
  //lp2->AddEntry((TObject*)0,"K^{*0} (#bar{K^{*0}}) #rightarrow K^{+}#pi^{-} (K^{-}#pi^{+})","");
  //lp2->AddEntry((TObject*)0,"K*^{0} (#bar{K}*^{0}) #rightarrow K^{+}#pi^{-} (K^{-}#pi^{+})","");
  //lp2->AddEntry((TObject*)0,"K*^{#pm} #rightarrow K^{0}_{S} #pi^{#pm}","");

   //lp2->AddEntry((TObject*)0,"|#it{y}| < 0.5","");
  // lp2->AddEntry((TObject*)0,"-0.5 < #it{y} < 0","");
   lp2->Draw("same");

   //  TLegend* lp3=DrawLegend(0.19,0.17,0.40,0.42);
   //  TLegend* lp3=DrawLegend(0.19,0.17,0.40,0.42);
   //  TLegend* lp3=DrawLegend(0.25,0.2,0.40,0.35);

   //rotational 
   TLegend* lp3=DrawLegend(0.5,0.6,0.8,0.8);
  lp3->AddEntry((TObject*)0,Form("%s",ptbn.Data()),"");
  //  lp3->AddEntry((TObject*)0,Form("%s",cos.Data()),"");
  lp3->AddEntry(hSig,"Data","p");
  lp3->AddEntry(fSigBg1,"4rBW","l");
  lp3->SetTextSize(.04);
  lp3->Draw("same");

   
  
  //    c11->SaveAs("Mix_lowpt_pPb5TeV_signal_new.gif");
  // c11->Print("Mix_lowpt_pPb5TeV_signal_new.eps");
   // c11->Print("InvMassFitIntermediatePt_pPb8TeV_kstar.root");
  //  c11->Print("Mix_lowpt_pPb5TeV_signal_new.pdf");


  /*   c11->SaveAs("Mix_highpt_pPb5TeV_signal_new.gif");
   c11->Print("Mix_highpt_pPb5TeV_signal_new.eps");
   // c11->Print("InvMassFitIntermediatePt_pPb8TeV_kstar.root");
  c11->Print("Mix_highpt_pPb5TeV_signal_new.pdf");
  */

   //c11->SaveAs("Mix_highpt_pPb5TeV_signal_high.gif");
   //c11->Print("Mix_highpt_pPb5TeV_signal_high.eps");
  //c11->Print("InvMassFitIntermediatePt_pPb8TeV_kstar.root");
  //c11->Print("Mix_highpt_pPb5TeV_signal_high.pdf");
     
  //   c11->SaveAs("Invmass_Preliminary_signal_k0sk0s_pp13TeV_0to20.png");
  // c11->Print("Invmass_Preliminary_signal_k0sk0s_pp13TeV_0to20.eps");
  // c11->Print("Invmass_Preliminary_signal_k0sk0s_pp13TeV_0to20.root");
  // c11->Print("Invmass_Preliminary_signal_k0sk0s_pp13TeV_0to20.pdf");

  //  Mass and width comparsion plot:
  /*
  TCanvas *c2 = DrawCanvas("c22");
  //c11->Range(0,0,1,1);
  c22->cd();
  c22->SetLeftMargin(0.15);
  c22->SetRightMargin(0.15);
  c22->SetTopMargin(0.15);
  c22->SetBottomMargin(0.15);
  
  //  Double_t Mass[16]= {1.2,1.25,1.30,1.35,1.40,1.45,1.50,1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.90,2.0};
  // Double_t Mass_value[16]= {1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1};
  // Double_t MassError[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  //  Double_t Mass[16]= {1.2,1.25,1.30,1.35,1.40,1.45,1.50,1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.90,2.0};
  // Double_t Mass_value[16]= {1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1};
  //  Double_t MassError[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  // Double_t Mass[2]= {1.2};
  // Double_t Mass_value[16]= {1};
  // Double_t MassError[16]={0};

  //  TGraphErrors *gr = new TGraphErrors(1,Mass,Mass_value,MassError,MassError);
  // gr->SetMarkerColor(1);
  // gr->SetMarkerSize(1);
  // gr->SetMarkerStyle(20);      
  // gr->Draw();
  //  Double_t Mass[1] = {1.275};
  // Double_t Mass_error[1] = {0.1};
  
  Double_t Mass[1] = {1275.5};
  Double_t Mass_error[1] = {0.1};

  Int_t const n = 1;
  Float_t x[n]  = {1.2755};
  Float_t y[n]  = {1};


  Float_t pdg2x[n]  = {1.525};
  Float_t pdg2y[n]  = {1};

  Float_t pdg3x[n]  = {1.723};
  Float_t pdg3y[n]  = {1};
  

  Float_t x1[n]  = {mass};
  Float_t y1[n]  = {1};

  Float_t x2[n]  = {mass2};
  Float_t y2[n]  = {1};
  
  Float_t x3[n]  = {mass3};
  Float_t y3[n]  = {1};



  
  // create the error arrays
  //  Float_t ex[n] = {.1867/2.0};
  Float_t ey[n] = {0.0};
  Float_t ex[n] = {.0008};

  Float_t pdg2ey[n] = {0.0};
  Float_t pdg2ex[n] = {.005};

  Float_t pdg3ey[n] = {0.0};
  Float_t pdg3ex[n] = {.006};

  

  Float_t ey1[n] = {0.0};
  Float_t ex1[n] = {masserr};

  Float_t ey2[n] = {0.0};
  Float_t ex2[n] = {masserr2};

  Float_t ey3[n] = {0.0};
  Float_t ex3[n] = {masserr3};

  //  TH1F* hf1=new TH1F("frame1","",1,1240,1320);
  TH1F* hf1=new TH1F("frame1","",1,1.2,2.0);
  hf1->SetXTitle("Mass (GeV/c^{2})");
  hf1->GetYaxis()->SetTitleOffset(0.9);
  hf1->GetYaxis()->SetNdivisions(10);
  hf1->GetYaxis()->SetTitleOffset(1.1);
  hf1->GetYaxis()->SetTitleSize(0.06);
  hf1->GetXaxis()->SetTitleSize(0.06);

  // create the TGraphErrors and draw it
  gr = new TGraphErrors(n,x,y,ex,ey);
  //   gr->SetTitle("TGraphErrors Example");
  gr->SetMarkerColor(1);
  gr->SetMarkerStyle(24);


  grpdg2 = new TGraphErrors(n,pdg2x,pdg2y,pdg2ex,pdg2ey);
  grpdg2->SetMarkerColor(2);
  grpdg2->SetMarkerStyle(24);

  grpdg3 = new TGraphErrors(n,pdg3x,pdg3y,pdg3ex,pdg3ey);
  //   gr->SetTitle("TGraphErrors Example");
  grpdg3->SetMarkerColor(4);
  grpdg3->SetMarkerStyle(24);
  
  // gr->SetLimits(1200,1400);

  gr1 = new TGraphErrors(n,x1,y1,ex1,ey1);
  //   gr->SetTitle("TGraphErrors Example");
  gr1->SetMarkerColor(1);
  gr1->SetMarkerSize(1.2);
  gr1->SetMarkerStyle(25);


  gr2 = new TGraphErrors(n,x2,y2,ex2,ey2);
  //   gr->SetTitle("TGraphErrors Example");
  gr2->SetMarkerColor(2);
  gr2->SetMarkerSize(1.2);
  gr2->SetMarkerStyle(25);


  gr3 = new TGraphErrors(n,x3,y3,ex3,ey3);
  //   gr->SetTitle("TGraphErrors Example");
  gr3->SetMarkerColor(4);
  gr3->SetMarkerSize(1.2);
  gr3->SetMarkerStyle(25);

    
  hf1->Draw("");
  hf1->GetYaxis()->SetRangeUser(0,2);
  gr->Draw("P same");
  grpdg2->Draw("P same");
  grpdg3->Draw("P same");
  gr1->Draw("P same");
  gr2->Draw("P same");
  gr3->Draw("P same");

  TLine *line = new TLine(1.275-0.0008,0,1.275-0.0008,1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  //  line->Draw();

  TLine *line1 = new TLine(1.275+0.0008,0,1.275+0.0008,1);
  line1->SetLineColor(1);
  line1->SetLineStyle(2);
  // line1->Draw();


  TLegend* lp3=DrawLegend(0.2,0.7,0.4,0.8);
  //  lp3->AddEntry((TObject*)0,Form("%s",ptbn.Data()),"");
  lp3->AddEntry(gr,"PDG (2018) (square marker)","");
  lp3->AddEntry(gr1,"This studies (open circle)","");
  lp3->SetTextSize(.03);
  lp3->Draw("same");
  */
  
  /**********************************************************/
  // width 
  /* TCanvas *c3 = DrawCanvas("c33");
  //c11->Range(0,0,1,1);
  c33->cd();
  c33->SetLeftMargin(0.15);
  c33->SetRightMargin(0.15);
  c33->SetTopMargin(0.15);
  c33->SetBottomMargin(0.15);
  


  Double_t Mass[1] = {1275.5};
  Double_t Mass_error[1] = {0.1};

  Int_t const n = 1;
  Float_t x[n]  = {1.275};
  Float_t y[n]  = {186.7};

  Float_t pdg2x[n]  = {1.525};
  Float_t pdg2y[n]  = {73.0};

  Float_t pdg3x[n]  = {1.713};
  Float_t pdg3y[n]  = {139.0};
  

  Float_t x1[n]  = {mass};
  Float_t y1[n]  = {width};

  Float_t x2[n]  = {mass2};
  Float_t y2[n]  = {width2};
  
  Float_t x3[n]  = {mass3};
  Float_t y3[n]  = {width3};


  // create the error arrays
  //  Float_t ex[n] = {.1867/2.0};
  Float_t ey[n] = {2};
  Float_t ex[n] = {0.0008};

  Float_t pdg2ey[n] = {6};
  Float_t pdg2ex[n] = {0.005};

  Float_t pdg3ey[n] = {8};
  Float_t pdg3ex[n] = {0.006};

  

  Float_t ey1[n] = {widtherr};
  Float_t ex1[n] = {masserr};

  Float_t ey2[n] = {widtherr2};
  Float_t ex2[n] = {masserr2};

  Float_t ey3[n] = {widtherr3};
  Float_t ex3[n] = {masserr3};

  //  TH1F* hf1=new TH1F("frame1","",1,1240,1320);
  TH1F* hf1=new TH1F("frame1","",1,1.2,2.0);
  hf1->SetYTitle("Width (MeV/c^{2})");
  hf1->SetXTitle("Mass (GeV/c^{2})");
  hf1->GetYaxis()->SetTitleOffset(0.9);
  hf1->GetYaxis()->SetNdivisions(10);
  hf1->GetYaxis()->SetTitleOffset(1.1);
  hf1->GetYaxis()->SetTitleSize(0.06);
  hf1->GetXaxis()->SetTitleSize(0.06);

  // create the TGraphErrors and draw it
  gr = new TGraphErrors(n,x,y,ex,ey);
  //   gr->SetTitle("TGraphErrors Example");
  gr->SetMarkerColor(1);
  gr->SetMarkerStyle(25);


  grpdg2 = new TGraphErrors(n,pdg2x,pdg2y,pdg2ex,pdg2ey);
  grpdg2->SetMarkerColor(2);
  grpdg2->SetMarkerStyle(25);

  grpdg3 = new TGraphErrors(n,pdg3x,pdg3y,pdg3ex,pdg3ey);
  //   gr->SetTitle("TGraphErrors Example");
  grpdg3->SetMarkerColor(4);
  grpdg3->SetMarkerStyle(25);
  
  // gr->SetLimits(1200,1400);

  gr1 = new TGraphErrors(n,x1,y1,ex1,ey1);
  //   gr->SetTitle("TGraphErrors Example");
  gr1->SetMarkerColor(1);
  gr1->SetMarkerSize(1.2);
  gr1->SetMarkerStyle(24);


  gr2 = new TGraphErrors(n,x2,y2,ex2,ey2);
  //   gr->SetTitle("TGraphErrors Example");
  gr2->SetMarkerColor(2);
  gr2->SetMarkerSize(1.2);
  gr2->SetMarkerStyle(24);


  gr3 = new TGraphErrors(n,x3,y3,ex3,ey3);
  //   gr->SetTitle("TGraphErrors Example");
  gr3->SetMarkerColor(4);
  gr3->SetMarkerSize(1.2);
  gr3->SetMarkerStyle(24);

    
  hf1->Draw("");
   hf1->GetYaxis()->SetRangeUser(0,300);
    gr->Draw("P same");
   grpdg2->Draw("P same");
    grpdg3->Draw("P same");
     gr1->Draw("P same");
    gr2->Draw("P same");
    gr3->Draw("P same");

  TLine *line = new TLine(1.275-0.0008,0,1.275-0.0008,1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  //  line->Draw();

  TLine *line1 = new TLine(1.275+0.0008,0,1.275+0.0008,1);
  line1->SetLineColor(1);
  line1->SetLineStyle(2);
  // line1->Draw();


  TLegend* lp3=DrawLegend(0.2,0.7,0.4,0.8);
  //  lp3->AddEntry((TObject*)0,Form("%s",ptbn.Data()),"");
  lp3->AddEntry(gr,"PDG (2018) (square marker)","");
  lp3->AddEntry(gr1,"This studies (open circle)","");
  lp3->SetTextSize(.03);
  lp3->Draw("same");
  */
    
}

TCanvas * DrawCanvas(TString opt="c")
{
  TCanvas *c1 = new TCanvas(opt.Data(),opt.Data(),10,10,600,600);
  c1->cd(1);
  c1->SetLeftMargin(0.25);
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
  gr->GetXaxis()->SetTitleOffset(1.35);
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


Double_t rBreitWigner(Double_t *x, Double_t *par)
{

  ///double rBW1 = par[0]*par[2]/(2*3.14159)/((x[0]-par[1])**2+par[2]**2/4.);
  //double rBW2 = par[3]*par[5]/(2*3.14159)/((x[0]-par[4])**2+par[5]**2/4.);
  // double rBW3 = par[6]*par[8]/(2*3.14159)/((x[0]-par[7])**2+par[9]**2/4.);
  // double rBW4 = par[9]*par[11]/(2*3.14159)/((x[0]-par[10])**2+par[11]**2/4.);

  double  npart1 = x[0]*x[0]- 4*(0.4976*0.4976);
  double  dpart1 = par[1]*par[1] - 4*(0.4976*0.4976);
  double  dpart2 = par[4]*par[4] - 4*(0.4976*0.4976);
  double  dpart3 = par[7]*par[7] - 4*(0.4976*0.4976);
  double  dpart4 = par[10]*par[10] - 4*(0.4976*0.4976);

  double Gamma = par[2];
  double Gamma1 = par[5];
  double Gamma2 = par[8];
  double Gamma3 = par[11];

  double rBW1 = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);
  double rBW2 = par[3]*x[0]*par[4]*Gamma1/(TMath::Power((x[0]*x[0]-par[4]*par[4]),2.0)+par[4]*par[4]*Gamma1*Gamma1);
  double rBW3 = par[6]*x[0]*par[7]*Gamma2/(TMath::Power((x[0]*x[0]-par[7]*par[7]),2.0)+par[7]*par[7]*Gamma2*Gamma2);
  double rBW4 = par[9]*x[0]*par[10]*Gamma2/(TMath::Power((x[0]*x[0]-par[10]*par[10]),2.0)+par[10]*par[10]*Gamma2*Gamma2);

  // double sigfun =par[12]*(TMath::Power(((5*rBW1)-(3*rBW2)+2*rBW3),2))+par[13]*(TMath::Power(rBW4,2.0));

   double sigfun =par[12]*(TMath::Power(((par[13]*rBW1)+(par[14]*rBW2)+par[15]*rBW3),2))+par[16]*(TMath::Power(rBW4,2.0));

  double poly2 = par[17]*TMath::Power(x[0],par[18])*TMath::Exp(-x[0]*par[19]);


  return(sigfun+poly2);
    
}


Double_t rBW(Double_t *x, Double_t *par)
{
  //  double  npart1 = x[0]*x[0]- 4*(0.4976*0.4976);
  // double  dpart1 = par[1]*par[1] - 4*(0.4976*0.4976);

  //double Gamma = par[2]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);
  //double rBW = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);
  // double BF = TMath::Sqrt((x[0]*x[0]+par[3]*par[3]));
  // double PS =  x[0]*TMath::Exp(-BF/par[4])/(BF);

  // double rBW = ;

  double  npart1 = x[0]*x[0]- 4*(0.4976*0.4976);
  double  dpart1 = par[1]*par[1] - 4*(0.4976*0.4976);
  double  dpart2 = par[4]*par[4] - 4*(0.4976*0.4976);
  double  dpart3 = par[7]*par[7] - 4*(0.4976*0.4976);
  double  dpart4 = par[10]*par[10] - 4*(0.4976*0.4976);
  
  //  double Gamma = par[2]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);
  // double Gamma1 = par[5]*(TMath::Power(par[4]/x[0],1.0))*TMath::Power((npart1)/(dpart2),1.5);
  // double Gamma2 = par[8]*(TMath::Power(par[7]/x[0],1.0))*TMath::Power((npart1)/(dpart3),1.5);
  // double Gamma3 = par[11]*(TMath::Power(par[10]/x[0],1.0))*TMath::Power((npart1)/(dpart4),1.5);


  
  //double rBW1 = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);
  //  double rBW2 = par[3]*x[0]*par[4]*Gamma1/(TMath::Power((x[0]*x[0]-par[4]*par[4]),2.0)+par[4]*par[4]*Gamma1*Gamma1);
  // double rBW3 = par[6]*x[0]*par[7]*Gamma2/(TMath::Power((x[0]*x[0]-par[7]*par[7]),2.0)+par[7]*par[7]*Gamma2*Gamma2);
  // double rBW4 = par[9]*x[0]*par[10]*Gamma3/(TMath::Power((x[0]*x[0]-par[10]*par[10]),2.0)+par[10]*par[10]*Gamma3*Gamma3);
  // double rBW1 = par[0]*par[1]*par[1]*Gamma/(TMath::Power(( -x[0]*x[0]+par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);

  double Gamma = par[2];
  double Gamma1 = par[5];
  double Gamma2 = par[8];
  double Gamma3 = par[11];
  
  /*
  double  npart1 = x[0]*x[0]- 4*(0.4976*0.4976);
  double  dpart1 = par[1]*par[1] - 4*(0.4976*0.4976);
  double  dpart2 = par[4]*par[4] - 4*(0.4976*0.4976);
  
  double Gamma = par[2]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);
  
  double Gamma1 = par[5]*(TMath::Power(par[4]/x[0],1.0))*TMath::Power((npart1)/(dpart2),1.5);
  
  double rBW1 = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);

  double rBW2 = par[3]*x[0]*par[4]*Gamma1/(TMath::Power((x[0]*x[0]-par[4]*par[4]),2.0)+par[4]*par[4]*Gamma1*Gamma1);
  */


  /* double  npart1 = x[0]*x[0]- 4*(0.4976*0.4976);
  double  dpart1 = par[1]*par[1] - 4*(0.4976*0.4976);

  double Gamma = par[2]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);
  
  double Gamma1 = par[5]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);
  
  double rBW1 = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);

  double rBW2 = par[3]*x[0]*par[4]*Gamma/(TMath::Power((x[0]*x[0]-par[4]*par[4]),2.0)+par[4]*par[4]*Gamma1*Gamma1);
  */

  //  return(rBW1+rBW2+rBW3+rBW4);
    //  return(rBW1);

   
    double rBW1 = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);
    double rBW2 = par[3]*x[0]*par[4]*Gamma1/(TMath::Power((x[0]*x[0]-par[4]*par[4]),2.0)+par[4]*par[4]*Gamma1*Gamma1);
    double rBW3 = par[6]*x[0]*par[7]*Gamma2/(TMath::Power((x[0]*x[0]-par[7]*par[7]),2.0)+par[7]*par[7]*Gamma2*Gamma2);
    double rBW4 = par[9]*x[0]*par[10]*Gamma2/(TMath::Power((x[0]*x[0]-par[10]*par[10]),2.0)+par[10]*par[10]*Gamma2*Gamma2);

    //  double sigfun =par[12]*(TMath::Power((5*rBW1 -3*rBW2 +2*rBW3),2))+par[13]*(TMath::Power(rBW4,2.0));

        double sigfun =par[12]*(TMath::Power((par[13]*rBW1 +par[14]*rBW2 +par[15]*rBW3),2))+par[16]*(TMath::Power(rBW4,2.0));
      
      //      return(rBW1+rBW2+rBW3+rBW4);
      return(sigfun);


}


Double_t rBW1(Double_t *x, Double_t *par)
{
  //  double  npart1 = x[0]*x[0]- 4*(0.4976*0.4976);
  // double  dpart1 = par[1]*par[1] - 4*(0.4976*0.4976);

  //double Gamma = par[2]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);
  //double rBW = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);
  // double BF = TMath::Sqrt((x[0]*x[0]+par[3]*par[3]));
  // double PS =  x[0]*TMath::Exp(-BF/par[4])/(BF);

  // double rBW = ;

  double  npart1 = x[0]*x[0]- 4*(0.4976*0.4976);
  double  dpart1 = par[1]*par[1] - 4*(0.4976*0.4976);
  double  dpart2 = par[4]*par[4] - 4*(0.4976*0.4976);
  double  dpart3 = par[7]*par[7] - 4*(0.4976*0.4976);
  double  dpart4 = par[10]*par[10] - 4*(0.4976*0.4976);
  
  //  double Gamma = par[2]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);
  // double Gamma1 = par[5]*(TMath::Power(par[4]/x[0],1.0))*TMath::Power((npart1)/(dpart2),1.5);
  // double Gamma2 = par[8]*(TMath::Power(par[7]/x[0],1.0))*TMath::Power((npart1)/(dpart3),1.5);
  // double Gamma3 = par[11]*(TMath::Power(par[10]/x[0],1.0))*TMath::Power((npart1)/(dpart4),1.5);


  
  //double rBW1 = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);
  //  double rBW2 = par[3]*x[0]*par[4]*Gamma1/(TMath::Power((x[0]*x[0]-par[4]*par[4]),2.0)+par[4]*par[4]*Gamma1*Gamma1);
  // double rBW3 = par[6]*x[0]*par[7]*Gamma2/(TMath::Power((x[0]*x[0]-par[7]*par[7]),2.0)+par[7]*par[7]*Gamma2*Gamma2);
  // double rBW4 = par[9]*x[0]*par[10]*Gamma3/(TMath::Power((x[0]*x[0]-par[10]*par[10]),2.0)+par[10]*par[10]*Gamma3*Gamma3);
  // double rBW1 = par[0]*par[1]*par[1]*Gamma/(TMath::Power(( -x[0]*x[0]+par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);

  double Gamma = par[2];
  double Gamma1 = par[5];
  double Gamma2 = par[8];
  double Gamma3 = par[11];
  
  /*
  double  npart1 = x[0]*x[0]- 4*(0.4976*0.4976);
  double  dpart1 = par[1]*par[1] - 4*(0.4976*0.4976);
  double  dpart2 = par[4]*par[4] - 4*(0.4976*0.4976);
  
  double Gamma = par[2]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);
  
  double Gamma1 = par[5]*(TMath::Power(par[4]/x[0],1.0))*TMath::Power((npart1)/(dpart2),1.5);
  
  double rBW1 = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);

  double rBW2 = par[3]*x[0]*par[4]*Gamma1/(TMath::Power((x[0]*x[0]-par[4]*par[4]),2.0)+par[4]*par[4]*Gamma1*Gamma1);
  */


  /* double  npart1 = x[0]*x[0]- 4*(0.4976*0.4976);
  double  dpart1 = par[1]*par[1] - 4*(0.4976*0.4976);

  double Gamma = par[2]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);
  
  double Gamma1 = par[5]*(TMath::Power(par[1]/x[0],1.0))*TMath::Power((npart1)/(dpart1),1.5);
  
  double rBW1 = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);

  double rBW2 = par[3]*x[0]*par[4]*Gamma/(TMath::Power((x[0]*x[0]-par[4]*par[4]),2.0)+par[4]*par[4]*Gamma1*Gamma1);
  */

  //  return(rBW1+rBW2+rBW3+rBW4);
    //  return(rBW1);

    double rBW1 = par[0]*x[0]*par[1]*Gamma/(TMath::Power((x[0]*x[0]-par[1]*par[1]),2.0)+par[1]*par[1]*Gamma*Gamma);
    double rBW2 = par[3]*x[0]*par[4]*Gamma1/(TMath::Power((x[0]*x[0]-par[4]*par[4]),2.0)+par[4]*par[4]*Gamma1*Gamma1);
    double rBW3 = par[6]*x[0]*par[7]*Gamma2/(TMath::Power((x[0]*x[0]-par[7]*par[7]),2.0)+par[7]*par[7]*Gamma2*Gamma2);
    return(rBW1);
}

Double_t BreitWigner(Double_t *x, Double_t *par) 
{  
  return (par[0]*par[2]/(2*3.14159)/((x[0]-par[1])**2+par[2]**2/4.)+ par[3]+x[0]*par[4]+ x[0]*x[0]*par[5]);
  
}

Double_t polynomial2(Double_t *x, Double_t *par) 
{ 
  // double poly2 = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
  
  double poly2 = par[0]*TMath::Power(x[0],par[1])*TMath::Exp(-x[0]*par[2]);
  
  //  double poly2 = par[3]*TMath::Power(x[0],par[4])*TMath::Exp(-x[0]*par[5]);
  //  double poly2 = par[0]*TMath::Sqrt(TMath::Power(x[0] -(2*0.497),par[1]))*TMath::Power(par[2],1.5)*TMath::Exp(-par[2]*(TMath::Power(x[0]-(2*0.497)),par[1]));
  // double poly2 = par[3]*exp(x[0]*par[4]);

  //  double poly2 = par[0]*exp(x[0]*par[1])+par[2]*exp(x[0]*par[3]*x[0])+par[4]*exp(x[0]*x[0]*x[0]*par[5]);
 
  return (poly2);
  // return (bg);

  //  return ((x[0]- 0.63718)**par[3])*exp(par[0] + x[0]*par[1] + x[0]*x[0]*par[2]);
}

Double_t BW(Double_t *x, Double_t *par) 
{
  return  par[0]*par[2]/(2*3.14159)/((x[0]-par[1])**2+par[2]**2/4.);
  //return (0.5*par[0]*par[1]/TMath::Pi() /((x[0]-par[2])*(x[0]-par[2]) + 0.25*par[1]*par[1]));
}

void *DrawFrame(TH1 *h, Int_t MCol,Int_t MSty, Bool_t mrk = 0)
{
  //  h->GetXaxis()->SetTitle("#it{M}_{K^{0}_{S}#pi^{#pm}} (GeV/#it{c^{2}})");
  
   h->GetXaxis()->SetTitle("#it{M}_{K^{0}_{S}K^{0}_{S}} (GeV/#it{c^{2}})");
   if(mrk){
    h->SetMarkerColor(MCol);
    h->SetMarkerStyle(MSty);
    h->SetMarkerSize(0.5);
  }
  h->SetLineColor(MCol);
  h->SetLineWidth(2);
  //h->GetXaxis()->CenterTitle(true);
  h->GetXaxis()->SetNdivisions(509);
  h->GetYaxis()->SetNdivisions(509);
  h->GetXaxis()->SetLabelOffset(0.015);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetTitleFont(42);
  // h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetTickLength(0.01);
  h->GetYaxis()->SetTickLength(0.01);
  h->GetXaxis()->SetTitleOffset(1.4);
  h->GetYaxis()->SetTitleOffset(1.65);
  //h->GetYaxis()->CenterTitle(true);
  h->GetYaxis()->SetDecimals(true);
  //h->GetYaxis()->SetNdivisions(310);
  h->GetYaxis()->SetLabelOffset(0.015);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTickLength(0.04);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitle("Counts / (0.02 GeV/#it{c^{2}})"); 
  h->GetYaxis()->SetTitleFont(42);
  h->Draw("");
  return h;
}
