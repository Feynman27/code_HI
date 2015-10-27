{
//=========Macro generated from canvas: EffectiveSignalCent0TrkPt0/EffectiveSignalCent0TrkPt0
//=========  (Thu Jun 27 09:42:23 2013) by ROOT version5.32/00
   TCanvas *EffectiveSignalCent0TrkPt0 = new TCanvas("EffectiveSignalCent0TrkPt0", "EffectiveSignalCent0TrkPt0",1682,52,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   EffectiveSignalCent0TrkPt0->Range(-0.204962,-324.0506,1.076051,1701.266);
   EffectiveSignalCent0TrkPt0->SetFillColor(0);
   EffectiveSignalCent0TrkPt0->SetBorderMode(0);
   EffectiveSignalCent0TrkPt0->SetBorderSize(2);
   EffectiveSignalCent0TrkPt0->SetTickx(1);
   EffectiveSignalCent0TrkPt0->SetTicky(1);
   EffectiveSignalCent0TrkPt0->SetLeftMargin(0.16);
   EffectiveSignalCent0TrkPt0->SetRightMargin(0.05);
   EffectiveSignalCent0TrkPt0->SetTopMargin(0.05);
   EffectiveSignalCent0TrkPt0->SetBottomMargin(0.16);
   EffectiveSignalCent0TrkPt0->SetFrameBorderMode(0);
   EffectiveSignalCent0TrkPt0->SetFrameBorderMode(0);
   
   TH1F *hDummy2__1 = new TH1F("hDummy2__1","hDummy2",50,0,1.1);
   hDummy2__1->SetMinimum(0);
   hDummy2__1->SetMaximum(1600);
   hDummy2__1->SetDirectory(0);
   hDummy2__1->SetStats(0);
   hDummy2__1->SetLineWidth(2);
   hDummy2__1->SetMarkerStyle(20);
   hDummy2__1->SetMarkerSize(1.2);
   hDummy2__1->GetXaxis()->SetTitle("i_{#mu} = #frac{#Sigma p_{T}^{trk}(#Delta R, p_{T}^{trk})}{p_{T}^{#mu}}");
   hDummy2__1->GetXaxis()->SetRange(1,46);
   hDummy2__1->GetXaxis()->SetLabelFont(42);
   hDummy2__1->GetXaxis()->SetTitleSize(0.03);
   hDummy2__1->GetXaxis()->SetTitleOffset(1.8);
   hDummy2__1->GetXaxis()->SetTitleFont(42);
   hDummy2__1->GetYaxis()->SetTitle("N_{eff}");
   hDummy2__1->GetYaxis()->SetLabelFont(42);
   hDummy2__1->GetYaxis()->SetTitleSize(0.05);
   hDummy2__1->GetYaxis()->SetTitleOffset(1.4);
   hDummy2__1->GetYaxis()->SetTitleFont(42);
   hDummy2__1->GetZaxis()->SetLabelFont(42);
   hDummy2__1->GetZaxis()->SetLabelSize(0.05);
   hDummy2__1->GetZaxis()->SetTitleSize(0.05);
   hDummy2__1->GetZaxis()->SetTitleFont(42);
   hDummy2__1->Draw("");
   TLatex *   tex = new TLatex(0.67,0.41,"p_{T}^{trk}>0.5GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.75,0.48,"0-5%");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
   EffectiveSignalCent0TrkPt0->Modified();
   EffectiveSignalCent0TrkPt0->cd();
   EffectiveSignalCent0TrkPt0->SetSelected(EffectiveSignalCent0TrkPt0);
}
