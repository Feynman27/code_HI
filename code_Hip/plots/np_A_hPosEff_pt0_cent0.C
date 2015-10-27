{
//=========Macro generated from canvas: c1/c1
//=========  (Fri Aug  8 21:31:20 2014) by ROOT version5.34/19
   TCanvas *c1 = new TCanvas("c1", "c1",82,52,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-3.512658,-364.5569,2.816456,1913.924);
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
   
   TH1D *np_A_hPosEff_pt0_cent0_total__3 = new TH1D("np_A_hPosEff_pt0_cent0_total__3","mu^{+},eff,all eta (total)",11, xAxis3);
   np_A_hPosEff_pt0_cent0_total__3->SetBinContent(1,655);
   np_A_hPosEff_pt0_cent0_total__3->SetBinContent(2,1467);
   np_A_hPosEff_pt0_cent0_total__3->SetBinContent(3,1575);
   np_A_hPosEff_pt0_cent0_total__3->SetBinContent(4,1483);
   np_A_hPosEff_pt0_cent0_total__3->SetBinContent(5,1577);
   np_A_hPosEff_pt0_cent0_total__3->SetBinContent(6,602);
   np_A_hPosEff_pt0_cent0_total__3->SetBinContent(7,1400);
   np_A_hPosEff_pt0_cent0_total__3->SetBinContent(8,1220);
   np_A_hPosEff_pt0_cent0_total__3->SetBinContent(9,1311);
   np_A_hPosEff_pt0_cent0_total__3->SetBinContent(10,1307);
   np_A_hPosEff_pt0_cent0_total__3->SetBinContent(11,772);
   np_A_hPosEff_pt0_cent0_total__3->SetMinimum(0);
   np_A_hPosEff_pt0_cent0_total__3->SetMaximum(1800);
   np_A_hPosEff_pt0_cent0_total__3->SetEntries(13369);
   np_A_hPosEff_pt0_cent0_total__3->SetDirectory(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   np_A_hPosEff_pt0_cent0_total__3->SetLineColor(ci);
   np_A_hPosEff_pt0_cent0_total__3->GetXaxis()->SetLabelFont(42);
   np_A_hPosEff_pt0_cent0_total__3->GetXaxis()->SetLabelSize(0.035);
   np_A_hPosEff_pt0_cent0_total__3->GetXaxis()->SetTitleSize(0.035);
   np_A_hPosEff_pt0_cent0_total__3->GetXaxis()->SetTitleFont(42);
   np_A_hPosEff_pt0_cent0_total__3->GetYaxis()->SetLabelFont(42);
   np_A_hPosEff_pt0_cent0_total__3->GetYaxis()->SetLabelSize(0.035);
   np_A_hPosEff_pt0_cent0_total__3->GetYaxis()->SetTitleSize(0.035);
   np_A_hPosEff_pt0_cent0_total__3->GetYaxis()->SetTitleFont(42);
   np_A_hPosEff_pt0_cent0_total__3->GetZaxis()->SetLabelFont(42);
   np_A_hPosEff_pt0_cent0_total__3->GetZaxis()->SetLabelSize(0.035);
   np_A_hPosEff_pt0_cent0_total__3->GetZaxis()->SetTitleSize(0.035);
   np_A_hPosEff_pt0_cent0_total__3->GetZaxis()->SetTitleFont(42);
   np_A_hPosEff_pt0_cent0_total__3->Draw("");
   Double_t xAxis4[12] = {-2.5, -2.1, -1.55, -1.05, -0.6, -0.1, 0.1, 0.6, 1.05, 1.55, 2.1, 2.5}; 
   
   TH1D *np_A_hPosEff_pt0_cent0_passed__4 = new TH1D("np_A_hPosEff_pt0_cent0_passed__4","mu^{+},eff,all eta (passed)",11, xAxis4);
   np_A_hPosEff_pt0_cent0_passed__4->SetBinContent(1,468);
   np_A_hPosEff_pt0_cent0_passed__4->SetBinContent(2,1154);
   np_A_hPosEff_pt0_cent0_passed__4->SetBinContent(3,1200);
   np_A_hPosEff_pt0_cent0_passed__4->SetBinContent(4,984);
   np_A_hPosEff_pt0_cent0_passed__4->SetBinContent(5,1151);
   np_A_hPosEff_pt0_cent0_passed__4->SetBinContent(6,191);
   np_A_hPosEff_pt0_cent0_passed__4->SetBinContent(7,1063);
   np_A_hPosEff_pt0_cent0_passed__4->SetBinContent(8,852);
   np_A_hPosEff_pt0_cent0_passed__4->SetBinContent(9,1093);
   np_A_hPosEff_pt0_cent0_passed__4->SetBinContent(10,1088);
   np_A_hPosEff_pt0_cent0_passed__4->SetBinContent(11,561);
   np_A_hPosEff_pt0_cent0_passed__4->SetEntries(9805);
   np_A_hPosEff_pt0_cent0_passed__4->SetDirectory(0);
   np_A_hPosEff_pt0_cent0_passed__4->SetLineColor(2);
   np_A_hPosEff_pt0_cent0_passed__4->GetXaxis()->SetLabelFont(42);
   np_A_hPosEff_pt0_cent0_passed__4->GetXaxis()->SetLabelSize(0.035);
   np_A_hPosEff_pt0_cent0_passed__4->GetXaxis()->SetTitleSize(0.035);
   np_A_hPosEff_pt0_cent0_passed__4->GetXaxis()->SetTitleFont(42);
   np_A_hPosEff_pt0_cent0_passed__4->GetYaxis()->SetLabelFont(42);
   np_A_hPosEff_pt0_cent0_passed__4->GetYaxis()->SetLabelSize(0.035);
   np_A_hPosEff_pt0_cent0_passed__4->GetYaxis()->SetTitleSize(0.035);
   np_A_hPosEff_pt0_cent0_passed__4->GetYaxis()->SetTitleFont(42);
   np_A_hPosEff_pt0_cent0_passed__4->GetZaxis()->SetLabelFont(42);
   np_A_hPosEff_pt0_cent0_passed__4->GetZaxis()->SetLabelSize(0.035);
   np_A_hPosEff_pt0_cent0_passed__4->GetZaxis()->SetTitleSize(0.035);
   np_A_hPosEff_pt0_cent0_passed__4->GetZaxis()->SetTitleFont(42);
   np_A_hPosEff_pt0_cent0_passed__4->Draw("same");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
