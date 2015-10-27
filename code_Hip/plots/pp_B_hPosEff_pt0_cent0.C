{
//=========Macro generated from canvas: c1/c1
//=========  (Fri Aug  8 20:57:29 2014) by ROOT version5.34/19
   TCanvas *c1 = new TCanvas("c1", "c1",82,52,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-3.512658,-303.7975,2.816456,1594.937);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.16);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.16);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   Double_t xAxis3[12] = {-2.5, -2.1, -1.55, -1.05, -0.6, -0.1, 0.1, 0.6, 1.05, 1.55, 2.1, 2.5}; 
   
   TH1D *pp_B_hPosEff_pt0_cent0_total__3 = new TH1D("pp_B_hPosEff_pt0_cent0_total__3","mu^{+},eff,all eta (total)",11, xAxis3);
   pp_B_hPosEff_pt0_cent0_total__3->SetBinContent(1,896);
   pp_B_hPosEff_pt0_cent0_total__3->SetBinContent(2,1286);
   pp_B_hPosEff_pt0_cent0_total__3->SetBinContent(3,1144);
   pp_B_hPosEff_pt0_cent0_total__3->SetBinContent(4,1005);
   pp_B_hPosEff_pt0_cent0_total__3->SetBinContent(5,1113);
   pp_B_hPosEff_pt0_cent0_total__3->SetBinContent(6,502);
   pp_B_hPosEff_pt0_cent0_total__3->SetBinContent(7,1166);
   pp_B_hPosEff_pt0_cent0_total__3->SetBinContent(8,1012);
   pp_B_hPosEff_pt0_cent0_total__3->SetBinContent(9,1114);
   pp_B_hPosEff_pt0_cent0_total__3->SetBinContent(10,1027);
   pp_B_hPosEff_pt0_cent0_total__3->SetBinContent(11,484);
   pp_B_hPosEff_pt0_cent0_total__3->SetMinimum(0);
   pp_B_hPosEff_pt0_cent0_total__3->SetMaximum(1500);
   pp_B_hPosEff_pt0_cent0_total__3->SetEntries(10749);
   pp_B_hPosEff_pt0_cent0_total__3->SetDirectory(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   pp_B_hPosEff_pt0_cent0_total__3->SetLineColor(ci);
   pp_B_hPosEff_pt0_cent0_total__3->GetXaxis()->SetLabelFont(42);
   pp_B_hPosEff_pt0_cent0_total__3->GetXaxis()->SetLabelSize(0.035);
   pp_B_hPosEff_pt0_cent0_total__3->GetXaxis()->SetTitleSize(0.035);
   pp_B_hPosEff_pt0_cent0_total__3->GetXaxis()->SetTitleFont(42);
   pp_B_hPosEff_pt0_cent0_total__3->GetYaxis()->SetLabelFont(42);
   pp_B_hPosEff_pt0_cent0_total__3->GetYaxis()->SetLabelSize(0.035);
   pp_B_hPosEff_pt0_cent0_total__3->GetYaxis()->SetTitleSize(0.035);
   pp_B_hPosEff_pt0_cent0_total__3->GetYaxis()->SetTitleFont(42);
   pp_B_hPosEff_pt0_cent0_total__3->GetZaxis()->SetLabelFont(42);
   pp_B_hPosEff_pt0_cent0_total__3->GetZaxis()->SetLabelSize(0.035);
   pp_B_hPosEff_pt0_cent0_total__3->GetZaxis()->SetTitleSize(0.035);
   pp_B_hPosEff_pt0_cent0_total__3->GetZaxis()->SetTitleFont(42);
   pp_B_hPosEff_pt0_cent0_total__3->Draw("");
   Double_t xAxis4[12] = {-2.5, -2.1, -1.55, -1.05, -0.6, -0.1, 0.1, 0.6, 1.05, 1.55, 2.1, 2.5}; 
   
   TH1D *pp_B_hPosEff_pt0_cent0_passed__4 = new TH1D("pp_B_hPosEff_pt0_cent0_passed__4","mu^{+},eff,all eta (passed)",11, xAxis4);
   pp_B_hPosEff_pt0_cent0_passed__4->SetBinContent(1,638);
   pp_B_hPosEff_pt0_cent0_passed__4->SetBinContent(2,1013);
   pp_B_hPosEff_pt0_cent0_passed__4->SetBinContent(3,893);
   pp_B_hPosEff_pt0_cent0_passed__4->SetBinContent(4,647);
   pp_B_hPosEff_pt0_cent0_passed__4->SetBinContent(5,875);
   pp_B_hPosEff_pt0_cent0_passed__4->SetBinContent(6,175);
   pp_B_hPosEff_pt0_cent0_passed__4->SetBinContent(7,873);
   pp_B_hPosEff_pt0_cent0_passed__4->SetBinContent(8,623);
   pp_B_hPosEff_pt0_cent0_passed__4->SetBinContent(9,884);
   pp_B_hPosEff_pt0_cent0_passed__4->SetBinContent(10,873);
   pp_B_hPosEff_pt0_cent0_passed__4->SetBinContent(11,349);
   pp_B_hPosEff_pt0_cent0_passed__4->SetEntries(7843);
   pp_B_hPosEff_pt0_cent0_passed__4->SetDirectory(0);
   pp_B_hPosEff_pt0_cent0_passed__4->SetLineColor(2);
   pp_B_hPosEff_pt0_cent0_passed__4->GetXaxis()->SetLabelFont(42);
   pp_B_hPosEff_pt0_cent0_passed__4->GetXaxis()->SetLabelSize(0.035);
   pp_B_hPosEff_pt0_cent0_passed__4->GetXaxis()->SetTitleSize(0.035);
   pp_B_hPosEff_pt0_cent0_passed__4->GetXaxis()->SetTitleFont(42);
   pp_B_hPosEff_pt0_cent0_passed__4->GetYaxis()->SetLabelFont(42);
   pp_B_hPosEff_pt0_cent0_passed__4->GetYaxis()->SetLabelSize(0.035);
   pp_B_hPosEff_pt0_cent0_passed__4->GetYaxis()->SetTitleSize(0.035);
   pp_B_hPosEff_pt0_cent0_passed__4->GetYaxis()->SetTitleFont(42);
   pp_B_hPosEff_pt0_cent0_passed__4->GetZaxis()->SetLabelFont(42);
   pp_B_hPosEff_pt0_cent0_passed__4->GetZaxis()->SetLabelSize(0.035);
   pp_B_hPosEff_pt0_cent0_passed__4->GetZaxis()->SetTitleSize(0.035);
   pp_B_hPosEff_pt0_cent0_passed__4->GetZaxis()->SetTitleFont(42);
   pp_B_hPosEff_pt0_cent0_passed__4->Draw("same");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
