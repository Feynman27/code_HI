{
//=========Macro generated from canvas: c1/c1
//=========  (Fri Aug  8 20:54:24 2014) by ROOT version5.34/19
   TCanvas *c1 = new TCanvas("c1", "c1",82,52,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-3.512658,-162.0253,2.816456,850.6329);
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
   Double_t xAxis1[12] = {-2.5, -2.1, -1.55, -1.05, -0.6, -0.1, 0.1, 0.6, 1.05, 1.55, 2.1, 2.5}; 
   
   TH1D *pn_B_hPosEff_pt0_cent0_total__1 = new TH1D("pn_B_hPosEff_pt0_cent0_total__1","mu^{+},eff,all eta (total)",11, xAxis1);
   pn_B_hPosEff_pt0_cent0_total__1->SetBinContent(1,389);
   pn_B_hPosEff_pt0_cent0_total__1->SetBinContent(2,614);
   pn_B_hPosEff_pt0_cent0_total__1->SetBinContent(3,612);
   pn_B_hPosEff_pt0_cent0_total__1->SetBinContent(4,566);
   pn_B_hPosEff_pt0_cent0_total__1->SetBinContent(5,619);
   pn_B_hPosEff_pt0_cent0_total__1->SetBinContent(6,260);
   pn_B_hPosEff_pt0_cent0_total__1->SetBinContent(7,679);
   pn_B_hPosEff_pt0_cent0_total__1->SetBinContent(8,596);
   pn_B_hPosEff_pt0_cent0_total__1->SetBinContent(9,685);
   pn_B_hPosEff_pt0_cent0_total__1->SetBinContent(10,625);
   pn_B_hPosEff_pt0_cent0_total__1->SetBinContent(11,298);
   pn_B_hPosEff_pt0_cent0_total__1->SetMinimum(0);
   pn_B_hPosEff_pt0_cent0_total__1->SetMaximum(800);
   pn_B_hPosEff_pt0_cent0_total__1->SetEntries(5943);
   pn_B_hPosEff_pt0_cent0_total__1->SetDirectory(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   pn_B_hPosEff_pt0_cent0_total__1->SetLineColor(ci);
   pn_B_hPosEff_pt0_cent0_total__1->GetXaxis()->SetTitle("eta");
   pn_B_hPosEff_pt0_cent0_total__1->GetXaxis()->SetLabelFont(42);
   pn_B_hPosEff_pt0_cent0_total__1->GetXaxis()->SetLabelSize(0.035);
   pn_B_hPosEff_pt0_cent0_total__1->GetXaxis()->SetTitleSize(0.035);
   pn_B_hPosEff_pt0_cent0_total__1->GetXaxis()->SetTitleFont(42);
   pn_B_hPosEff_pt0_cent0_total__1->GetYaxis()->SetLabelFont(42);
   pn_B_hPosEff_pt0_cent0_total__1->GetYaxis()->SetLabelSize(0.035);
   pn_B_hPosEff_pt0_cent0_total__1->GetYaxis()->SetTitleSize(0.035);
   pn_B_hPosEff_pt0_cent0_total__1->GetYaxis()->SetTitleFont(42);
   pn_B_hPosEff_pt0_cent0_total__1->GetZaxis()->SetLabelFont(42);
   pn_B_hPosEff_pt0_cent0_total__1->GetZaxis()->SetLabelSize(0.035);
   pn_B_hPosEff_pt0_cent0_total__1->GetZaxis()->SetTitleSize(0.035);
   pn_B_hPosEff_pt0_cent0_total__1->GetZaxis()->SetTitleFont(42);
   pn_B_hPosEff_pt0_cent0_total__1->Draw("");
   Double_t xAxis2[12] = {-2.5, -2.1, -1.55, -1.05, -0.6, -0.1, 0.1, 0.6, 1.05, 1.55, 2.1, 2.5}; 
   
   TH1D *pn_B_hPosEff_pt0_cent0_passed__2 = new TH1D("pn_B_hPosEff_pt0_cent0_passed__2","mu^{+},eff,all eta (passed)",11, xAxis2);
   pn_B_hPosEff_pt0_cent0_passed__2->SetBinContent(1,291);
   pn_B_hPosEff_pt0_cent0_passed__2->SetBinContent(2,479);
   pn_B_hPosEff_pt0_cent0_passed__2->SetBinContent(3,464);
   pn_B_hPosEff_pt0_cent0_passed__2->SetBinContent(4,377);
   pn_B_hPosEff_pt0_cent0_passed__2->SetBinContent(5,451);
   pn_B_hPosEff_pt0_cent0_passed__2->SetBinContent(6,93);
   pn_B_hPosEff_pt0_cent0_passed__2->SetBinContent(7,502);
   pn_B_hPosEff_pt0_cent0_passed__2->SetBinContent(8,387);
   pn_B_hPosEff_pt0_cent0_passed__2->SetBinContent(9,551);
   pn_B_hPosEff_pt0_cent0_passed__2->SetBinContent(10,521);
   pn_B_hPosEff_pt0_cent0_passed__2->SetBinContent(11,207);
   pn_B_hPosEff_pt0_cent0_passed__2->SetEntries(4323);
   pn_B_hPosEff_pt0_cent0_passed__2->SetDirectory(0);
   pn_B_hPosEff_pt0_cent0_passed__2->SetLineColor(2);
   pn_B_hPosEff_pt0_cent0_passed__2->GetXaxis()->SetLabelFont(42);
   pn_B_hPosEff_pt0_cent0_passed__2->GetXaxis()->SetLabelSize(0.035);
   pn_B_hPosEff_pt0_cent0_passed__2->GetXaxis()->SetTitleSize(0.035);
   pn_B_hPosEff_pt0_cent0_passed__2->GetXaxis()->SetTitleFont(42);
   pn_B_hPosEff_pt0_cent0_passed__2->GetYaxis()->SetLabelFont(42);
   pn_B_hPosEff_pt0_cent0_passed__2->GetYaxis()->SetLabelSize(0.035);
   pn_B_hPosEff_pt0_cent0_passed__2->GetYaxis()->SetTitleSize(0.035);
   pn_B_hPosEff_pt0_cent0_passed__2->GetYaxis()->SetTitleFont(42);
   pn_B_hPosEff_pt0_cent0_passed__2->GetZaxis()->SetLabelFont(42);
   pn_B_hPosEff_pt0_cent0_passed__2->GetZaxis()->SetLabelSize(0.035);
   pn_B_hPosEff_pt0_cent0_passed__2->GetZaxis()->SetTitleSize(0.035);
   pn_B_hPosEff_pt0_cent0_passed__2->GetZaxis()->SetTitleFont(42);
   pn_B_hPosEff_pt0_cent0_passed__2->Draw("same");
   TLatex *   tex = new TLatex(-1.639386,706.8869,"#mu^{+},pn,PeriodB");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
