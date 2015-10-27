{
//=========Macro generated from canvas: cMinus0/cMinus0
//=========  (Thu Dec 19 13:50:33 2013) by ROOT version5.34/13
   TCanvas *cMinus0 = new TCanvas("cMinus0", "cMinus0",550,80,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cMinus0->Range(-0.5326582,-0.2025316,2.796456,1.063291);
   cMinus0->SetFillColor(0);
   cMinus0->SetBorderMode(0);
   cMinus0->SetBorderSize(2);
   cMinus0->SetTickx(1);
   cMinus0->SetTicky(1);
   cMinus0->SetLeftMargin(0.16);
   cMinus0->SetRightMargin(0.05);
   cMinus0->SetTopMargin(0.05);
   cMinus0->SetBottomMargin(0.16);
   cMinus0->SetFrameBorderMode(0);
   cMinus0->SetFrameBorderMode(0);
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(9);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,0.225,0.5560454);
   grae->SetPointError(0,0.125,0.125,0.02384805,0.02354019);
   grae->SetPoint(1,0.475,0.5433753);
   grae->SetPointError(1,0.125,0.125,0.02365522,0.02340175);
   grae->SetPoint(2,0.7,0.5846025);
   grae->SetPointError(2,0.1,0.1,0.0259739,0.02549075);
   grae->SetPoint(3,0.925,0.5505488);
   grae->SetPointError(3,0.125,0.125,0.02360179,0.02334843);
   grae->SetPoint(4,1.21,0.432867);
   grae->SetPointError(4,0.16,0.16,0.02103812,0.02111763);
   grae->SetPoint(5,1.445,0.5122262);
   grae->SetPointError(5,0.075,0.075,0.03042791,0.03007654);
   grae->SetPoint(6,1.63,0.5171619);
   grae->SetPointError(6,0.11,0.11,0.02533743,0.02504905);
   grae->SetPoint(7,1.92,0.4519408);
   grae->SetPointError(7,0.18,0.18,0.02022128,0.0202332);
   grae->SetPoint(8,2.25,0.2620026);
   grae->SetPointError(8,0.15,0.15,0.01912502,0.01957449);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,0,2.63);
   Graph_Graph1->SetMinimum(0);
   Graph_Graph1->SetMaximum(1);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   Graph_Graph1->SetLineWidth(2);
   Graph_Graph1->SetMarkerStyle(20);
   Graph_Graph1->SetMarkerSize(1.2);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph1);
   
   grae->Draw("ape");
   
   grae = new TGraphAsymmErrors(9);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(24);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,0.225,0.5553962);
   grae->SetPointError(0,0.125,0.125,0.02291429,0.02263467);
   grae->SetPoint(1,0.475,0.543115);
   grae->SetPointError(1,0.125,0.125,0.02261526,0.02238637);
   grae->SetPoint(2,0.7,0.5836567);
   grae->SetPointError(2,0.1,0.1,0.02518965,0.02474209);
   grae->SetPoint(3,0.925,0.5517042);
   grae->SetPointError(3,0.125,0.125,0.02344163,0.02318886);
   grae->SetPoint(4,1.21,0.431133);
   grae->SetPointError(4,0.16,0.16,0.02158312,0.02167202);
   grae->SetPoint(5,1.445,0.511878);
   grae->SetPointError(5,0.075,0.075,0.03195326,0.0315657);
   grae->SetPoint(6,1.63,0.5166727);
   grae->SetPointError(6,0.11,0.11,0.02677489,0.02645364);
   grae->SetPoint(7,1.92,0.4506768);
   grae->SetPointError(7,0.18,0.18,0.0210381,0.02105452);
   grae->SetPoint(8,2.25,0.2469554);
   grae->SetPointError(8,0.15,0.15,0.01776631,0.01819857);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,0,2.63);
   Graph_Graph2->SetMinimum(0.1912682);
   Graph_Graph2->SetMaximum(0.6463197);
   Graph_Graph2->SetDirectory(0);
   Graph_Graph2->SetStats(0);
   Graph_Graph2->SetLineWidth(2);
   Graph_Graph2->SetMarkerStyle(20);
   Graph_Graph2->SetMarkerSize(1.2);
   Graph_Graph2->GetXaxis()->SetLabelFont(42);
   Graph_Graph2->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph2->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph2->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph2->GetXaxis()->SetTitleFont(42);
   Graph_Graph2->GetYaxis()->SetLabelFont(42);
   Graph_Graph2->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph2->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph2->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph2->GetYaxis()->SetTitleFont(42);
   Graph_Graph2->GetZaxis()->SetLabelFont(42);
   Graph_Graph2->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph2->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph2->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph2);
   
   grae->Draw("pe");
   
   TLegend *leg = new TLegend(0.1845638,0.2045455,0.454698,0.3548951,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph0","Nominal","pe");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","ReWeighted","pe");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   leg->Draw();
   TLatex *   tex = new TLatex(0.5302013,0.2954545,"#mu^{-}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.1153846);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.5285235,0.222028,"0-5%");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
   cMinus0->Modified();
   cMinus0->cd();
   cMinus0->SetSelected(cMinus0);
}
