#include <iostream>
using namespace std;
#include "style.h"

void plots()
{
    TFile *f = new TFile("f1710.root", "read");
    // TFile *f = new TFile("f1525.root", "read");
    if (f->IsZombie())
    {
        cout << "Error opening file" << endl;
        return;
    }
    TH1F *GenpT = (TH1F *)f->Get("strangeness_tutorial/MChists/Genf1710");
    TH1F *GenMass = (TH1F *)f->Get("strangeness_tutorial/MChists/Genf1710_mass");
    TH1F *MC_mult = (TH1F *)f->Get("strangeness_tutorial/MChists/MC_mult");
    TH1F *MC_mult_sel = (TH1F *)f->Get("strangeness_tutorial/MChists/MC_mult_after_event_sel");
    TH1F *recpt1 = (TH1F *)f->Get("strangeness_tutorial/MChists/Recf1710_pt1");
    TH1F *recpt2 = (TH1F *)f->Get("strangeness_tutorial/MChists/Recf1710_pt2");
    TH1F *recMass = (TH1F *)f->Get("strangeness_tutorial/MChists/Recf1710_mass");
    if (GenpT == nullptr || GenMass == nullptr || MC_mult == nullptr || MC_mult_sel == nullptr || recpt1 == nullptr || recpt2 == nullptr || recMass == nullptr)
    {
        cout << "Error reading histogram" << endl;
        return;
    }
    vector<string> histnames = {"Gen p_{T} f_{0}(1710)", "Gen Mass f_{0}(1710)", "Multiplicity", "Multiplicity after event selection", "Rec f_{0}(1710) p_{T}", "Rec f_{0}(1710) p_{T}", "Rec f_{0}(1710) Mass"};
    // vector<string> histnames = {"Gen p_{T} f_{2}(1525)", "Gen Mass f_{2}(1525)", "Multiplicity", "Multiplicity after event selection", "Rec f_{2}(1525) p_{T}", "Rec f_{2}(1525) p_{T}", "Rec f_{2}(1525) Mass"};
    vector<TH1F *> hists = {GenpT, GenMass, MC_mult, MC_mult_sel, recpt1, recpt2, recMass};
    int counter = 0;
    for (auto hist : hists)
    {
        SetHistostyle2(hist);
        hist->SetTitle(histnames[counter].c_str());
        hist->GetYaxis()->SetTitle("Counts");
        TCanvas *c = new TCanvas("", "", 720, 720);
        SetCanvasStyle(c, 0.15, 0.05, 0.08, 0.15);
        hist->Draw();

        std::string filename = Form("plots/%s.png", histnames[counter].c_str());
        c->SaveAs(filename.c_str());
        counter++;
    }
}