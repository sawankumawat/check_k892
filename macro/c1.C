#ifdef __CLING__
#pragma cling optimize(0)
#endif
void c1()
{
//=========Macro generated from canvas: c1/c1
//=========  (Tue May 16 16:55:34 2023) by ROOT version 6.26/10
   TCanvas *c1 = new TCanvas("c1", "c1",10,64,700,500);
   c1->Range(-2.0625,-0.02161483,18.5625,0.1945927);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   Double_t Graph1D_y1_fx3001[27] = {
   0.05,
   0.15,
   0.25,
   0.35,
   0.45,
   0.55,
   0.65,
   0.75,
   0.85,
   0.95,
   1.1,
   1.3,
   1.5,
   1.7,
   1.9,
   2.2,
   2.6,
   3,
   3.4,
   3.8,
   4.5,
   5.5,
   6.5,
   7.5,
   9,
   11,
   13.5};
   Double_t Graph1D_y1_fy3001[27] = {
   0.02569625,
   0.07687868,
   0.1148528,
   0.1331602,
   0.1429375,
   0.1420053,
   0.1443418,
   0.1438145,
   0.1369708,
   0.1266951,
   0.1037723,
   0.07988674,
   0.06020794,
   0.04592486,
   0.03413558,
   0.02416249,
   0.01583164,
   0.009337281,
   0.005522419,
   0.003434073,
   0.001625859,
   0.0006570263,
   0.0003157573,
   0.0001547131,
   6.176241e-05,
   2.169473e-05,
   8.492805e-06};
   Double_t Graph1D_y1_felx3001[27] = {
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.2,
   0.2,
   0.2,
   0.2,
   0.2,
   0.5,
   0.5,
   0.5,
   0.5,
   1,
   1,
   1.5};
   Double_t Graph1D_y1_fely3001[27] = {
   0.004822243,
   0.01268301,
   0.01666997,
   0.01629568,
   0.0143071,
   0.01292793,
   0.012906,
   0.01328889,
   0.01271818,
   0.01175787,
   0.009041361,
   0.007287144,
   0.005977865,
   0.004839623,
   0.003500099,
   0.002110949,
   0.001212539,
   0.0006962834,
   0.0004371703,
   0.0003009391,
   0.0001446117,
   6.483693e-05,
   3.532806e-05,
   2.086699e-05,
   9.060963e-06,
   4.051565e-06,
   1.905413e-06};
   Double_t Graph1D_y1_fehx3001[27] = {
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.2,
   0.2,
   0.2,
   0.2,
   0.2,
   0.5,
   0.5,
   0.5,
   0.5,
   1,
   1,
   1.5};
   Double_t Graph1D_y1_fehy3001[27] = {
   0.004822243,
   0.01268301,
   0.01666997,
   0.01629568,
   0.0143071,
   0.01292793,
   0.012906,
   0.01328889,
   0.01271818,
   0.01175787,
   0.009041361,
   0.007287144,
   0.005977865,
   0.004839623,
   0.003500099,
   0.002110949,
   0.001212539,
   0.0006962834,
   0.0004371703,
   0.0003009391,
   0.0001446117,
   6.483693e-05,
   3.532806e-05,
   2.086699e-05,
   9.060963e-06,
   4.051565e-06,
   1.905413e-06};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(27,Graph1D_y1_fx3001,Graph1D_y1_fy3001,Graph1D_y1_felx3001,Graph1D_y1_fehx3001,Graph1D_y1_fely3001,Graph1D_y1_fehy3001);
   grae->SetName("Graph1D_y1");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   
   TH1F *Graph_Graph1D_y13001 = new TH1F("Graph_Graph1D_y13001","",100,0,16.5);
   Graph_Graph1D_y13001->SetMinimum(5.928653e-06);
   Graph_Graph1D_y13001->SetMaximum(0.172972);
   Graph_Graph1D_y13001->SetDirectory(0);
   Graph_Graph1D_y13001->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph1D_y13001->SetLineColor(ci);
   Graph_Graph1D_y13001->GetXaxis()->SetTitle("$p_{T}$ [$GeV/c$]");
   Graph_Graph1D_y13001->GetXaxis()->SetLabelFont(42);
   Graph_Graph1D_y13001->GetXaxis()->SetTitleOffset(1);
   Graph_Graph1D_y13001->GetXaxis()->SetTitleFont(42);
   Graph_Graph1D_y13001->GetYaxis()->SetTitle("(1/N$_{\\rm ev}$)*d$^{2}N$/(d$p_{\\rm T}$d$y$) [$(GeV/c)^{-1}$]");
   Graph_Graph1D_y13001->GetYaxis()->SetLabelFont(42);
   Graph_Graph1D_y13001->GetYaxis()->SetTitleFont(42);
   Graph_Graph1D_y13001->GetZaxis()->SetLabelFont(42);
   Graph_Graph1D_y13001->GetZaxis()->SetTitleOffset(1);
   Graph_Graph1D_y13001->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph1D_y13001);
   
   grae->Draw("alp");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
