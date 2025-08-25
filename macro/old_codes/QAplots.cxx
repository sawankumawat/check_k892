#include<iostream>
using namespace std;
#include <iostream>
#include <cmath>
#include "TArrow.h"
#include "src/style.h"
#include "src/common.h"

void QAplots(){
    TFile *fInputFile = new TFile(kDataFilename.c_str(), "Read");
    string foldername = "lf-reso2initializer/Event/";
    TH1F *hvz = (TH1F *)fInputFile->Get((foldername + "posZ").c_str());
    TH1F *mult = (TH1F *)fInputFile->Get((foldername + "CentFT0C").c_str());
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1000);
    SetCanvasStyle2(c1, 0.15, 0.05, 0.05, 0.15);
    SetHistoQA(hvz);
    hvz->GetYaxis()->SetTitle("Counts");
    hvz->GetYaxis()->SetTitleOffset(1.5);
    hvz->GetXaxis()->SetTitle("V_{Z} (cm)");
    hvz->Draw();
    c1->SaveAs((kSignalOutput + "/vz.png").c_str());
    c1->Clear();
    SetHistoQA(mult);
    mult->GetYaxis()->SetTitle("Counts");
    mult->GetYaxis()->SetTitleOffset(1.5);
    mult->GetXaxis()->SetTitle("Centrality (%)");
    mult->Draw();
    cout<<"The number of bins in the centrality histogram is: "<<mult->GetNbinsX()<<endl;
    cout<<"The low x-axis limit is: "<<mult->GetXaxis()->GetBinLowEdge(1)<<endl;
    cout<<"The high x-axis limit is: "<<mult->GetXaxis()->GetBinUpEdge(mult->GetNbinsX())<<endl;
    c1->SaveAs((kSignalOutput + "/centrality.png").c_str());

}
