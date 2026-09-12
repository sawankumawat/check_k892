#ifdef __CLING__
#pragma cling optimize(0)
#endif
void c1()
{
//=========Macro generated from canvas: c1/c1
//=========  (Tue Sep  8 23:45:14 2026) by ROOT version 6.36.000
   TCanvas *c1 = new TCanvas("c1", "c1", 1650, 87, 720, 720);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(1);
   TColor::SetPalette(57, nullptr);
   c1->Range(0,0,1,1);
   c1->SetFillColor(10);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.1754386);
   c1->SetRightMargin(0.006027281);
   c1->SetTopMargin(0.02586207);
   c1->SetBottomMargin(0.08333334);
   c1->SetFrameLineWidth(2);
   c1->SetFrameBorderMode(0);
   
// ------------>Primitives in pad: c1_1
   TPad *c1_1__0 = new TPad("c1_1", "c1_1", 0, 0.3, 1, 1);
   c1_1__0->Draw();
   c1_1__0->cd();
   c1_1__0->Range(-2.051282,-5.826774,10.76923,0.2322343);
   c1_1__0->SetFillColor(0);
   c1_1__0->SetBorderMode(0);
   c1_1__0->SetBorderSize(2);
   c1_1__0->SetLogy();
   c1_1__0->SetTickx(1);
   c1_1__0->SetTicky(1);
   c1_1__0->SetLeftMargin(0.16);
   c1_1__0->SetRightMargin(0.06);
   c1_1__0->SetTopMargin(0.02);
   c1_1__0->SetBottomMargin(0.001);
   c1_1__0->SetFrameBorderMode(0);
   c1_1__0->SetFrameBorderMode(0);
   
   std::vector<Double_t> hmultClone0__1_x_vect0{
      0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8,
      2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8,
      10, 12, 15, 20
   };
   TH1F *hmultClone0__1 = new TH1F("hmultClone0__1", "", 23, hmultClone0__1_x_vect0.data());
   hmultClone0__1->SetBinContent(1,0.03358076140284538);
   hmultClone0__1->SetBinContent(2,0.07438167184591293);
   hmultClone0__1->SetBinContent(3,0.08301561325788498);
   hmultClone0__1->SetBinContent(4,0.0860920175909996);
   hmultClone0__1->SetBinContent(5,0.07832219451665878);
   hmultClone0__1->SetBinContent(6,0.0675927922129631);
   hmultClone0__1->SetBinContent(7,0.05259978398680687);
   hmultClone0__1->SetBinContent(8,0.04025240987539291);
   hmultClone0__1->SetBinContent(9,0.03160271048545837);
   hmultClone0__1->SetBinContent(10,0.0242327693849802);
   hmultClone0__1->SetBinContent(11,0.01559113245457411);
   hmultClone0__1->SetBinContent(12,0.008395189419388771);
   hmultClone0__1->SetBinContent(13,0.004507861565798521);
   hmultClone0__1->SetBinContent(14,0.002540447516366839);
   hmultClone0__1->SetBinContent(15,0.001491296570748091);
   hmultClone0__1->SetBinContent(16,0.000916872697416693);
   hmultClone0__1->SetBinContent(17,0.0004666532040573657);
   hmultClone0__1->SetBinContent(18,0.0001960183144547045);
   hmultClone0__1->SetBinContent(19,0.0001012189459288493);
   hmultClone0__1->SetBinContent(20,4.396432268549688e-05);
   hmultClone0__1->SetBinContent(21,1.597803566255607e-05);
   hmultClone0__1->SetBinContent(22,5.759656232839916e-06);
   hmultClone0__1->SetBinContent(23,1.678967123552866e-06);
   hmultClone0__1->SetBinError(1,0.0002365586315876048);
   hmultClone0__1->SetBinError(2,0.0004135218288015013);
   hmultClone0__1->SetBinError(3,0.0003202998006379551);
   hmultClone0__1->SetBinError(4,0.0002602874913860922);
   hmultClone0__1->SetBinError(5,0.0002122213625431193);
   hmultClone0__1->SetBinError(6,0.0001718440467903207);
   hmultClone0__1->SetBinError(7,0.000125343630964814);
   hmultClone0__1->SetBinError(8,7.797725042151744e-05);
   hmultClone0__1->SetBinError(9,6.075786550090196e-05);
   hmultClone0__1->SetBinError(10,4.593250762964202e-05);
   hmultClone0__1->SetBinError(11,2.045140049703351e-05);
   hmultClone0__1->SetBinError(12,1.538091790546306e-05);
   hmultClone0__1->SetBinError(13,1.033313151028091e-05);
   hmultClone0__1->SetBinError(14,7.312722306347435e-06);
   hmultClone0__1->SetBinError(15,4.949205495980177e-06);
   hmultClone0__1->SetBinError(16,4.073827423270718e-06);
   hmultClone0__1->SetBinError(17,2.105917988638474e-06);
   hmultClone0__1->SetBinError(18,1.247127244056787e-06);
   hmultClone0__1->SetBinError(19,8.439556610260486e-07);
   hmultClone0__1->SetBinError(20,4.18394737439702e-07);
   hmultClone0__1->SetBinError(21,2.095111225911694e-07);
   hmultClone0__1->SetBinError(22,1.060650911098701e-07);
   hmultClone0__1->SetBinError(23,5.115761266051584e-08);
   hmultClone0__1->SetMinimum(1.511070411197579e-06);
   hmultClone0__1->SetMaximum(1.291380263864994);
   hmultClone0__1->SetEntries(733864.5426864941);
   hmultClone0__1->SetStats(0);
   hmultClone0__1->SetLineColor(TColor::GetColor("#0000ff"));
   hmultClone0__1->SetLineWidth(2);
   hmultClone0__1->SetMarkerColor(TColor::GetColor("#0000ff"));
   hmultClone0__1->SetMarkerStyle(20);
   hmultClone0__1->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
   hmultClone0__1->GetXaxis()->SetRange(1, 20);
   hmultClone0__1->GetXaxis()->CenterTitle(true);
   hmultClone0__1->GetXaxis()->SetLabelFont(42);
   hmultClone0__1->GetXaxis()->SetLabelOffset(0.01499999966472387);
   hmultClone0__1->GetXaxis()->SetLabelSize(0.05714285746216774);
   hmultClone0__1->GetXaxis()->SetTitleSize(0.05000000074505806);
   hmultClone0__1->GetXaxis()->SetTickLength(0.01999999955296516);
   hmultClone0__1->GetXaxis()->SetTitleOffset(1.019999980926514);
   hmultClone0__1->GetXaxis()->SetTitleFont(42);
   hmultClone0__1->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
   hmultClone0__1->GetYaxis()->CenterTitle(true);
   hmultClone0__1->GetYaxis()->SetNdivisions(3000510);
   hmultClone0__1->GetYaxis()->SetLabelFont(42);
   hmultClone0__1->GetYaxis()->SetLabelOffset(0.01499999966472387);
   hmultClone0__1->GetYaxis()->SetLabelSize(0.05714285746216774);
   hmultClone0__1->GetYaxis()->SetTitleSize(0.05714285746216774);
   hmultClone0__1->GetYaxis()->SetTickLength(0.01999999955296516);
   hmultClone0__1->GetYaxis()->SetTitleOffset(1.299999952316284);
   hmultClone0__1->GetYaxis()->SetTitleFont(42);
   hmultClone0__1->GetZaxis()->SetLabelFont(42);
   hmultClone0__1->GetZaxis()->SetTitleOffset(1);
   hmultClone0__1->GetZaxis()->SetTitleFont(42);
   hmultClone0__1->Draw("pe");
   
   std::vector<Double_t> h2__2_x_vect1{
      0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8,
      2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8,
      10, 12, 15, 20
   };
   TH1F *h2__2 = new TH1F("h2__2", "", 23, h2__2_x_vect1.data());
   h2__2->SetBinContent(1,0.03358076140284538);
   h2__2->SetBinContent(2,0.07438167184591293);
   h2__2->SetBinContent(3,0.08301561325788498);
   h2__2->SetBinContent(4,0.0860920175909996);
   h2__2->SetBinContent(5,0.07832219451665878);
   h2__2->SetBinContent(6,0.0675927922129631);
   h2__2->SetBinContent(7,0.05259978398680687);
   h2__2->SetBinContent(8,0.04025240987539291);
   h2__2->SetBinContent(9,0.03160271048545837);
   h2__2->SetBinContent(10,0.0242327693849802);
   h2__2->SetBinContent(11,0.01559113245457411);
   h2__2->SetBinContent(12,0.008395189419388771);
   h2__2->SetBinContent(13,0.004507861565798521);
   h2__2->SetBinContent(14,0.002540447516366839);
   h2__2->SetBinContent(15,0.001491296570748091);
   h2__2->SetBinContent(16,0.000916872697416693);
   h2__2->SetBinContent(17,0.0004666532040573657);
   h2__2->SetBinContent(18,0.0001960183144547045);
   h2__2->SetBinContent(19,0.0001012189459288493);
   h2__2->SetBinContent(20,4.396432268549688e-05);
   h2__2->SetBinContent(21,1.597803566255607e-05);
   h2__2->SetBinContent(22,5.759656232839916e-06);
   h2__2->SetBinContent(23,1.678967123552866e-06);
   h2__2->SetBinError(1,0.002025152595429375);
   h2__2->SetBinError(2,0.004689689178826617);
   h2__2->SetBinError(3,0.005777416074350172);
   h2__2->SetBinError(4,0.005372744380611333);
   h2__2->SetBinError(5,0.004144582562449234);
   h2__2->SetBinError(6,0.003256747065805754);
   h2__2->SetBinError(7,0.002728951164428464);
   h2__2->SetBinError(8,0.001960468148478417);
   h2__2->SetBinError(9,0.001343179725399879);
   h2__2->SetBinError(10,0.001007421770823719);
   h2__2->SetBinError(11,0.0007225114879758632);
   h2__2->SetBinError(12,0.0004618067527392897);
   h2__2->SetBinError(13,0.0002535138043233236);
   h2__2->SetBinError(14,0.000149614837634939);
   h2__2->SetBinError(15,8.6881395631698e-05);
   h2__2->SetBinError(16,5.469074270637473e-05);
   h2__2->SetBinError(17,3.336933894271124e-05);
   h2__2->SetBinError(18,1.449406033434981e-05);
   h2__2->SetBinError(19,7.548205258514737e-06);
   h2__2->SetBinError(20,2.822836375946779e-06);
   h2__2->SetBinError(21,8.574618346409736e-07);
   h2__2->SetBinError(22,3.811118586089456e-07);
   h2__2->SetBinError(23,1.400465920455585e-07);
   h2__2->SetEntries(733864.5426864941);
   h2__2->SetStats(0);
   h2__2->SetFillStyle(0);
   h2__2->SetLineColor(TColor::GetColor("#0000ff"));
   h2__2->SetMarkerColor(TColor::GetColor("#0000ff"));
   h2__2->SetMarkerStyle(20);
   h2__2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
   h2__2->GetXaxis()->CenterTitle(true);
   h2__2->GetXaxis()->SetLabelFont(42);
   h2__2->GetXaxis()->SetLabelOffset(0.01499999966472387);
   h2__2->GetXaxis()->SetLabelSize(0.04500000178813934);
   h2__2->GetXaxis()->SetTitleSize(0.05000000074505806);
   h2__2->GetXaxis()->SetTickLength(0.01999999955296516);
   h2__2->GetXaxis()->SetTitleOffset(1.200000047683716);
   h2__2->GetXaxis()->SetTitleFont(42);
   h2__2->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
   h2__2->GetYaxis()->CenterTitle(true);
   h2__2->GetYaxis()->SetNdivisions(3000510);
   h2__2->GetYaxis()->SetLabelFont(42);
   h2__2->GetYaxis()->SetLabelOffset(0.01499999966472387);
   h2__2->GetYaxis()->SetLabelSize(0.04500000178813934);
   h2__2->GetYaxis()->SetTitleSize(0.05000000074505806);
   h2__2->GetYaxis()->SetTickLength(0.01999999955296516);
   h2__2->GetYaxis()->SetTitleOffset(1.799999952316284);
   h2__2->GetYaxis()->SetTitleFont(42);
   h2__2->GetZaxis()->SetLabelFont(42);
   h2__2->GetZaxis()->SetTitleOffset(1);
   h2__2->GetZaxis()->SetTitleFont(42);
   h2__2->Draw("e2 same");
   
   TF1 *fitfunc1 = new TF1("*fitfunc", 0,15,4);
   // The original function : fitfunc had originally been created by:
   // TF1 *fitfunc = new TF1("fitfunc", "fitfunc", 0,15,4, 1, TF1::EAddToList::kDefault);
   fitfunc1->SetRange(0,15);
   fitfunc1->SetName("fitfunc");
   fitfunc1->SetTitle("fitfunc");
   std::vector<Double_t> fitfunc1_vect2{
      0, 0.04197899872326613, 0.07466581759638749, 0.0933464109972502, 0.0987055864699822, 0.09450490786420605, 0.08495681560232733, 0.0733271805525505, 0.06165577523029103, 0.0509982156499778,
      0.04177153130235833, 0.03403420961595233, 0.0276702234109682, 0.02249585297810795, 0.01831575127153973, 0.01494911707242883, 0.01223961037624844, 0.01005708314580731, 0.00829559908259103, 0.006870083476728011,
      0.005712750468158127, 0.004769816402105292, 0.003998681654012478, 0.003365606341458299, 0.002843837082981633, 0.002412118315955574, 0.002053519189370989, 0.001754513580775746, 0.001504260321745057, 0.001294040429887852,
      0.001116816855179274, 0.0009668895807825993, 0.0008396248743798091, 0.0007312422198320226, 0.0006386461707236933, 0.0005592932530677102, 0.0004910862764924541, 0.0004322901347884139, 0.0003814645029131055, 0.0003374098591857132,
      0.0002991240490447113, 0.000265767214864965, 0.0002366333867508843, 0.0002111273939440281, 0.0001887460399684536, 0.0001690627055754231, 0.000151714716229359, 0.0001363929462351568, 0.000122833238026801, 0.0001108092990583778,
      0.00010012680511861, 9.061849155481679e-05, 8.214005580184845e-05, 7.456672806253558e-05, 6.779039376502261e-05, 6.171717292276712e-05, 5.626537883533503e-05, 5.136379254792845e-05, 4.695020080795897e-05, 4.2970154449086e-05,
      3.937591161727644e-05, 3.612553636312427e-05, 3.318212812529687e-05, 3.051316173330453e-05, 2.808992093324822e-05, 2.588701122370278e-05, 2.388194008975692e-05, 2.205475462961451e-05, 2.0387728151326e-05, 1.886508863493755e-05,
      1.747278305453432e-05, 1.619827247354327e-05, 1.503035359644557e-05, 1.395900310626036e-05, 1.297524166073896e-05, 1.20710148783584e-05, 1.123908903210461e-05, 1.047295949639036e-05, 9.766770269941942e-06, 9.115243133139354e-06,
      8.513615198772746e-06, 7.957583786032513e-06, 7.443257693409888e-06, 6.967114070915363e-06, 6.525960198860323e-06, 6.116899572109028e-06, 5.737301767472394e-06, 5.384775639704185e-06, 5.057145449990511e-06, 4.75242958126395e-06,
      4.46882153828041e-06, 4.2046729681493e-06, 3.958478469739202e-06, 3.728861988798208e-06, 3.514564620333279e-06, 3.314433661299358e-06, 3.127412775396316e-06, 2.952533148134018e-06, 2.78890552462517e-06, 2.635713035077404e-06,
      2.492204723917293e-06, 0, 15
   };
   for (int n = 0; n < 103; n++)
      fitfunc1->SetSavedPoint(n, fitfunc1_vect2[n]);
   fitfunc1->SetFillColor(19);
   fitfunc1->SetMarkerColor(1);
   fitfunc1->SetMarkerStyle(1);
   fitfunc1->SetMarkerSize(1);
   fitfunc1->SetLineColor(TColor::GetColor("#0000ff"));
   fitfunc1->SetLineStyle(2);
   fitfunc1->SetLineWidth(2);
   fitfunc1->SetChisquare(23.80237);
   fitfunc1->SetNDF(17);
   fitfunc1->SetParameter(0, 7.168667);
   fitfunc1->SetParError(0, 0.1545016);
   fitfunc1->SetParLimits(0, 0, 0);
   fitfunc1->SetParameter(1, 0.1370742);
   fitfunc1->SetParError(1, 0.00223751);
   fitfunc1->SetParLimits(1, 0, 0);
   fitfunc1->SetParameter(2, 0.895);
   fitfunc1->SetParError(2, 0);
   fitfunc1->SetParLimits(2, 0.895, 0.895);
   fitfunc1->SetParameter(3, 0.3065515);
   fitfunc1->SetParError(3, 0.005988623);
   fitfunc1->SetParLimits(3, 0, 0);
   fitfunc1->Draw("same");
   
   std::vector<Double_t> grae_fx_vect3{ 0.2, 0.6000000000000001, 1, 1.4, 1.8, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 8.5 };
   std::vector<Double_t> grae_fy_vect4{ 0.05300698047011514, 0.09170837875551903, 0.07298830592353071, 0.04413417656610546, 0.02604805934654906, 0.01460476869853881, 0.008463470846354076, 0.004434809360399388, 0.002410225924757044, 0.00132216941487977, 0.0008237187862852278, 0.0004132322483240314, 0.0002050113771512048, 5.606951660117619e-05 };
   std::vector<Double_t> grae_fexl_vect5{ 0.2, 0.2000000000000001, 0.2, 0.2, 0.2, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 1.5 };
   std::vector<Double_t> grae_fexh_vect6{ 0.2, 0.2, 0.2, 0.2000000000000002, 0.2, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 1.5 };
   std::vector<Double_t> grae_feyl_vect7{ 0.006149467583799429, 0.006020093266295797, 0.003983369413511632, 0.002237867333534318, 0.001238462327814617, 0.0005759077884314467, 0.0003019087182659686, 0.0001583341363383146, 9.267662893275737e-05, 6.512166578221106e-05, 4.823807065969782e-05, 2.181059998188954e-05, 1.314740192587113e-05, 3.535362896139518e-06 };
   std::vector<Double_t> grae_feyh_vect8{ 0.006149467583799429, 0.006020093266295797, 0.003983369413511632, 0.002237867333534318, 0.001238462327814617, 0.0005759077884314467, 0.0003019087182659686, 0.0001583341363383146, 9.267662893275737e-05, 6.512166578221106e-05, 4.823807065969782e-05, 2.181059998188954e-05, 1.314740192587113e-05, 3.535362896139518e-06 };
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(14, grae_fx_vect3.data(), grae_fy_vect4.data(), grae_fexl_vect5.data(), grae_fexh_vect6.data(), grae_feyl_vect7.data(), grae_feyh_vect8.data());
   grae->SetName("Graph1D_y1");
   grae->SetTitle("doi:10.17182/hepdata.96957.v1/t8");
   grae->SetFillStyle(1000);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(22);
   
   TH1F *Graph_histogram1 = new TH1F("Graph_histogram1", "doi:10.17182/hepdata.96957.v1/t8", 100, 0, 11);
   Graph_histogram1->SetMinimum(4.728073833453301e-05);
   Graph_histogram1->SetMaximum(0.1074960658086258);
   Graph_histogram1->SetDirectory(nullptr);
   Graph_histogram1->SetStats(0);
   Graph_histogram1->SetLineColor(TColor::GetColor("#000099"));
   Graph_histogram1->GetXaxis()->SetTitle("$p_{\\mathrm{T}}$ [GeV/$c$]");
   Graph_histogram1->GetXaxis()->SetLabelFont(42);
   Graph_histogram1->GetXaxis()->SetTitleOffset(1);
   Graph_histogram1->GetXaxis()->SetTitleFont(42);
   Graph_histogram1->GetYaxis()->SetTitle("(1/N$_{\\mathrm{ev}}$)*d$^{2}N$/(d$p_{\\mathrm{T}}$d$y$) [$(\\mathrm{GeV}/c)^{-1}$]");
   Graph_histogram1->GetYaxis()->SetLabelFont(42);
   Graph_histogram1->GetYaxis()->SetTitleFont(42);
   Graph_histogram1->GetZaxis()->SetLabelFont(42);
   Graph_histogram1->GetZaxis()->SetTitleOffset(1);
   Graph_histogram1->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_histogram1);
   
   grae->Draw("pe ");
   
   TLegend *leg = new TLegend(0.233983, 0.066092, 0.614206, 0.265189, nullptr, "brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *legentry = leg->AddEntry("hmultClone0","pp 13.6 TeV","p");
   legentry->SetMarkerColor(TColor::GetColor("#0000ff"));
   legentry->SetMarkerStyle(20);
   legentry->SetTextFont(42);
   legentry = leg->AddEntry("Graph1D_y1","pp 13 TeV (#it{PLB 807 (2020) 135501)}","p");
   legentry->SetMarkerStyle(22);
   legentry->SetTextFont(42);
   legentry = leg->AddEntry("fitfunc","Levy-Tsallis","l");
   legentry->SetLineColor(TColor::GetColor("#0000ff"));
   legentry->SetLineStyle(2);
   legentry->SetLineWidth(2);
   legentry->SetTextFont(42);
   leg->Draw();
   TLatex *tex = new TLatex(0.23, 0.9, "ALICE");
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->SetNDC();
   tex->Draw();
   tex = new TLatex(0.23, 0.85, "|y| < 0.5, INEL > 0");
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->SetNDC();
   tex->Draw();
   tex = new TLatex(0.23, 0.8, "K*(892)^{0}");
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->SetNDC();
   tex->Draw();
   c1_1__0->Modified();
   c1->cd();
   
// ------------>Primitives in pad: c1_2
   TPad *c1_2__1 = new TPad("c1_2", "c1_2", 0, 0, 1, 0.3);
   c1_2__1->Draw();
   c1_2__1->cd();
   c1_2__1->Range(-2.053333,0.2553811,10.78,1.451196);
   c1_2__1->SetFillColor(0);
   c1_2__1->SetBorderMode(0);
   c1_2__1->SetBorderSize(2);
   c1_2__1->SetTickx(1);
   c1_2__1->SetTicky(1);
   c1_2__1->SetLeftMargin(0.16);
   c1_2__1->SetRightMargin(0.06);
   c1_2__1->SetTopMargin(0.001);
   c1_2__1->SetBottomMargin(0.33);
   c1_2__1->SetFrameBorderMode(0);
   c1_2__1->SetFrameBorderMode(0);
   
   std::vector<Double_t> gre_fx_vect9{ 0.2, 0.6000000000000001, 1, 1.4, 1.8, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 8.5 };
   std::vector<Double_t> gre_fy_vect10{ 1.023659549664156, 1.076298456143432, 1.058875162920545, 1.082053945227479, 1.062275812675236, 1.02357780400346, 1.000048037360852, 1.038935013131281, 1.053458053642057, 1.081644532148229, 1.068055388890625, 1.091244776399866, 0.9918584749726774, 0.9441945427622813 };
   std::vector<Double_t> gre_fex_vect11{ 0.2, 0.2, 0.2, 0.2000000000000001, 0.2, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 1.5 };
   std::vector<Double_t> gre_fey_vect12{ 0.1187572120063601, 0.07065240031804293, 0.05778858521697841, 0.05486662186887125, 0.05050620310112394, 0.04036259947411062, 0.03877574818492845, 0.04031810924605693, 0.04263891225464626, 0.05327493808659782, 0.06254674796252242, 0.05759643250717676, 0.06360799193319672, 0.05953449495494961 };
   TGraphErrors *gre = new TGraphErrors(14, gre_fx_vect9.data(), gre_fy_vect10.data(), gre_fex_vect11.data(), gre_fey_vect12.data());
   gre->SetName("");
   gre->SetTitle("");
   gre->SetFillStyle(1000);
   gre->SetLineColor(TColor::GetColor("#0000ff"));
   gre->SetLineWidth(2);
   gre->SetMarkerColor(TColor::GetColor("#0000ff"));
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_histogram2 = new TH1F("Graph_histogram2", "", 100, 0, 11);
   Graph_histogram2->SetMinimum(0.65);
   Graph_histogram2->SetMaximum(1.45);
   Graph_histogram2->SetDirectory(nullptr);
   Graph_histogram2->SetStats(0);
   Graph_histogram2->SetLineColor(TColor::GetColor("#000099"));
   Graph_histogram2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
   Graph_histogram2->GetXaxis()->SetRange(1, 91);
   Graph_histogram2->GetXaxis()->CenterTitle(true);
   Graph_histogram2->GetXaxis()->SetLabelFont(42);
   Graph_histogram2->GetXaxis()->SetLabelOffset(0.01499999966472387);
   Graph_histogram2->GetXaxis()->SetLabelSize(0.1333333402872086);
   Graph_histogram2->GetXaxis()->SetTitleSize(0.1333333402872086);
   Graph_histogram2->GetXaxis()->SetTickLength(0.03999999910593033);
   Graph_histogram2->GetXaxis()->SetTitleOffset(1.100000023841858);
   Graph_histogram2->GetXaxis()->SetTitleFont(42);
   Graph_histogram2->GetYaxis()->SetTitle("#frac{This Analysis}{Published}");
   Graph_histogram2->GetYaxis()->CenterTitle(true);
   Graph_histogram2->GetYaxis()->SetNdivisions(506);
   Graph_histogram2->GetYaxis()->SetLabelFont(42);
   Graph_histogram2->GetYaxis()->SetLabelOffset(0.01499999966472387);
   Graph_histogram2->GetYaxis()->SetLabelSize(0.1333333402872086);
   Graph_histogram2->GetYaxis()->SetTitleSize(0.116666667163372);
   Graph_histogram2->GetYaxis()->SetTickLength(0.03999999910593033);
   Graph_histogram2->GetYaxis()->SetTitleOffset(0.6000000238418579);
   Graph_histogram2->GetYaxis()->SetTitleFont(42);
   Graph_histogram2->GetZaxis()->SetLabelFont(42);
   Graph_histogram2->GetZaxis()->SetTitleOffset(1);
   Graph_histogram2->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_histogram2);
   
   gre->Draw("ap");
   
   std::vector<Double_t> gre_fx_vect13{ 0.2, 0.6000000000000001, 1, 1.4, 1.8, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.5, 6.5, 8.5 };
   std::vector<Double_t> gre_fy_vect14{ 1.023659549664156, 1.076298456143432, 1.058875162920545, 1.082053945227479, 1.062275812675236, 1.02357780400346, 1.000048037360852, 1.038935013131281, 1.053458053642057, 1.081644532148229, 1.068055388890625, 1.091244776399866, 0.9918584749726774, 0.9441945427622813 };
   std::vector<Double_t> gre_fex_vect15{ 0.2, 0.2, 0.2, 0.2000000000000001, 0.2, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 1.5 };
   std::vector<Double_t> gre_fey_vect16{ 0.135162040766427, 0.09748526324808093, 0.07708707118534282, 0.07607703674775471, 0.06709046455711423, 0.0622824630676312, 0.06730367094096386, 0.07098842177086627, 0.07528089974915546, 0.08251770240380818, 0.08927984454276469, 0.09698670175366095, 0.09708134301147194, 0.08496859996734024 };
   gre = new TGraphErrors(14, gre_fx_vect13.data(), gre_fy_vect14.data(), gre_fex_vect15.data(), gre_fey_vect16.data());
   gre->SetName("");
   gre->SetTitle("");
   gre->SetFillColor(TColor::GetColor("#0000ff"));
   gre->SetFillStyle(0);
   gre->SetLineColor(TColor::GetColor("#0000ff"));
   gre->SetLineWidth(2);
   gre->SetMarkerColor(TColor::GetColor("#0000ff"));
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_histogram3 = new TH1F("Graph_histogram3", "", 100, 0, 11);
   Graph_histogram3->SetMinimum(0.8263253892590825);
   Graph_histogram3->SetMaximum(1.221132031689386);
   Graph_histogram3->SetDirectory(nullptr);
   Graph_histogram3->SetStats(0);
   Graph_histogram3->SetLineColor(TColor::GetColor("#000099"));
   Graph_histogram3->GetXaxis()->SetLabelFont(42);
   Graph_histogram3->GetXaxis()->SetTitleOffset(1);
   Graph_histogram3->GetXaxis()->SetTitleFont(42);
   Graph_histogram3->GetYaxis()->SetLabelFont(42);
   Graph_histogram3->GetYaxis()->SetTitleFont(42);
   Graph_histogram3->GetZaxis()->SetLabelFont(42);
   Graph_histogram3->GetZaxis()->SetTitleOffset(1);
   Graph_histogram3->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_histogram3);
   
   gre->Draw("e2 ");
   TLine *line = new TLine(0, 1, 10, 1);
   line->SetLineStyle(2);
   line->SetLineWidth(2);
   line->Draw();
   
   TBox *box = new TBox(0, 0.9, 10, 1.1);
   box->SetFillColor(TColor::GetColor("#666666"));
   box->SetFillStyle(3003);
   box->Draw("same");
   c1_2__1->Modified();
   c1->cd();
   c1->Modified();
   c1->SetSelected(c1);
}
