#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include "../src/common_glue.h"
#include "../src/style.h"
using namespace std;

void analyze_mc_13TeV()
{
    TFile *fmc = new TFile("/home/sawan/check_k892/mc/AnalysisResultsMC_KsKspp_13TeV_nopileupmother.root", "read");
    if (fmc->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    // print the names of all lists in the root file
    TIter next(fmc->GetListOfKeys());
    TKey *key;
    cout << "The folders in the root file are: \n";
    while ((key = (TKey *)next()))
    {
        cout << key->GetName() << endl;
    }

    TList *list = (TList *)fmc->Get("RsnOut_output2023_3.0_3.0_0.0_0.00_12.000_0_1_1_1.0_1.00_0.0_0.000_0.0");
    if (list == nullptr)
    {
        cout << "Error opening list" << endl;
        return;
    }

    // Iterate in the list to see all objects
    TIter iterator2(list);
    TObject *obj;

    while ((obj = iterator2()))
    {
        cout << "Object Name: " << obj->GetName()
             << ", Object Type: " << obj->ClassName() << endl;
    }

    THnSparseF *hmc1525 = (THnSparseF *)list->FindObject("RsnTaskF0_f0_MotherMCf1525_MIX");
    if(hmc1525 == nullptr)
    {
        cout << "Error opening histogram" << endl;
        return;
    }
    TH1F *hmc1525_proj = (TH1F *)hmc1525->Projection(0);
    if(hmc1525_proj == nullptr)
    {
        cout << "Error opening histogram" << endl;
        return;
    }
    TCanvas *c = new TCanvas("", "", 720, 720);
    SetCanvasStyle(c, 0.14, 0.03, 0.05, 0.14);
    hmc1525_proj->GetXaxis()->SetRangeUser(1,3);
    SetHistoQA(hmc1525_proj);
    hmc1525_proj->Draw();
}