{
//=========Macro generated from canvas: cEtaPlus/cEtaPlus
//=========  (Sun Aug 11 21:54:01 2013) by ROOT version5.32/00
   TCanvas *cEtaPlus = new TCanvas("cEtaPlus", "cEtaPlus",769,52,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cEtaPlus->Range(-0.3951658,-2.025316,2.567746,10.63291);
   cEtaPlus->SetFillColor(0);
   cEtaPlus->SetBorderMode(0);
   cEtaPlus->SetBorderSize(2);
   cEtaPlus->SetTickx(1);
   cEtaPlus->SetTicky(1);
   cEtaPlus->SetLeftMargin(0.16);
   cEtaPlus->SetRightMargin(0.05);
   cEtaPlus->SetTopMargin(0.05);
   cEtaPlus->SetBottomMargin(0.16);
   cEtaPlus->SetFrameBorderMode(0);
   cEtaPlus->SetFrameBorderMode(0);
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(9);
   grae->SetName("grWpc");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#cc0000");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#cc0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(23);
   grae->SetMarkerSize(1.5);
   grae->SetPoint(0,0.225,8.701805);
   grae->SetPointError(0,0.125,0.125,0.3783174,0.3783174);
   grae->SetPoint(1,0.475,7.649763);
   grae->SetPointError(1,0.125,0.125,0.3463731,0.3463731);
   grae->SetPoint(2,0.7,7.896176);
   grae->SetPointError(2,0.1,0.1,0.3840557,0.3840557);
   grae->SetPoint(3,0.925,7.436997);
   grae->SetPointError(3,0.125,0.125,0.3345559,0.3345559);
   grae->SetPoint(4,1.175,6.334746);
   grae->SetPointError(4,0.125,0.125,0.3451581,0.3451581);
   grae->SetPoint(5,1.425,5.989713);
   grae->SetPointError(5,0.125,0.125,0.3222833,0.3222833);
   grae->SetPoint(6,1.7,4.934416);
   grae->SetPointError(6,0.15,0.15,0.2775243,0.2775243);
   grae->SetPoint(7,1.975,4.234815);
   grae->SetPointError(7,0.125,0.125,0.2985153,0.2985153);
   grae->SetPoint(8,2.25,3.747698);
   grae->SetPointError(8,0.15,0.15,0.287349,0.287349);
   
   TH1F *Graph_grWpc1 = new TH1F("Graph_grWpc1","Graph",100,0,2.63);
   Graph_grWpc1->SetMinimum(0);
   Graph_grWpc1->SetMaximum(10);
   Graph_grWpc1->SetDirectory(0);
   Graph_grWpc1->SetStats(0);
   Graph_grWpc1->SetLineWidth(2);
   Graph_grWpc1->SetMarkerStyle(20);
   Graph_grWpc1->SetMarkerSize(1.2);
   Graph_grWpc1->GetXaxis()->SetTitle("|#eta|");
   Graph_grWpc1->GetXaxis()->SetRange(4,92);
   Graph_grWpc1->GetXaxis()->SetLabelFont(42);
   Graph_grWpc1->GetXaxis()->SetTitleSize(0.05);
   Graph_grWpc1->GetXaxis()->SetTitleOffset(1.4);
   Graph_grWpc1->GetXaxis()->SetTitleFont(42);
   Graph_grWpc1->GetYaxis()->SetTitle("#frac{dN}{d#eta} #frac{1}{N_{events} } #frac{10^{9}}{#LT N_{coll} #GT}");
   Graph_grWpc1->GetYaxis()->SetLabelFont(42);
   Graph_grWpc1->GetYaxis()->SetTitleOffset(1.4);
   Graph_grWpc1->GetYaxis()->SetTitleFont(42);
   Graph_grWpc1->GetZaxis()->SetLabelFont(42);
   Graph_grWpc1->GetZaxis()->SetLabelSize(0.05);
   Graph_grWpc1->GetZaxis()->SetTitleSize(0.05);
   Graph_grWpc1->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_grWpc1);
   
   grae->Draw("ape");
   
   grae = new TGraphAsymmErrors(9);
   grae->SetName("grWpUncorrelatedSystc");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#cc9999");
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#cc0000");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#cc0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(23);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,0.225,8.701805);
   grae->SetPointError(0,0.125,0.125,0.493279,0.493279);
   grae->SetPoint(1,0.475,7.649763);
   grae->SetPointError(1,0.125,0.125,0.4414774,0.4414774);
   grae->SetPoint(2,0.7,7.896176);
   grae->SetPointError(2,0.1,0.1,0.4638913,0.4638913);
   grae->SetPoint(3,0.925,7.436997);
   grae->SetPointError(3,0.125,0.125,0.4172951,0.4172951);
   grae->SetPoint(4,1.175,6.334746);
   grae->SetPointError(4,0.125,0.125,0.4047249,0.4047249);
   grae->SetPoint(5,1.425,5.989713);
   grae->SetPointError(5,0.125,0.125,0.3718353,0.3718353);
   grae->SetPoint(6,1.7,4.934416);
   grae->SetPointError(6,0.15,0.15,0.301129,0.301129);
   grae->SetPoint(7,1.975,4.234815);
   grae->SetPointError(7,0.125,0.125,0.3141863,0.3141863);
   grae->SetPoint(8,2.25,3.747698);
   grae->SetPointError(8,0.15,0.15,0.3027305,0.3027305);
   
   TH1F *Graph_grWpUncorrelatedSystc2 = new TH1F("Graph_grWpUncorrelatedSystc2","Graph",100,0,2.63);
   Graph_grWpUncorrelatedSystc2->SetMinimum(2.869956);
   Graph_grWpUncorrelatedSystc2->SetMaximum(9.770096);
   Graph_grWpUncorrelatedSystc2->SetDirectory(0);
   Graph_grWpUncorrelatedSystc2->SetStats(0);
   Graph_grWpUncorrelatedSystc2->SetLineWidth(2);
   Graph_grWpUncorrelatedSystc2->SetMarkerStyle(20);
   Graph_grWpUncorrelatedSystc2->SetMarkerSize(1.2);
   Graph_grWpUncorrelatedSystc2->GetXaxis()->SetLabelFont(42);
   Graph_grWpUncorrelatedSystc2->GetXaxis()->SetLabelSize(0.05);
   Graph_grWpUncorrelatedSystc2->GetXaxis()->SetTitleSize(0.05);
   Graph_grWpUncorrelatedSystc2->GetXaxis()->SetTitleOffset(1.4);
   Graph_grWpUncorrelatedSystc2->GetXaxis()->SetTitleFont(42);
   Graph_grWpUncorrelatedSystc2->GetYaxis()->SetLabelFont(42);
   Graph_grWpUncorrelatedSystc2->GetYaxis()->SetLabelSize(0.05);
   Graph_grWpUncorrelatedSystc2->GetYaxis()->SetTitleSize(0.05);
   Graph_grWpUncorrelatedSystc2->GetYaxis()->SetTitleOffset(1.4);
   Graph_grWpUncorrelatedSystc2->GetYaxis()->SetTitleFont(42);
   Graph_grWpUncorrelatedSystc2->GetZaxis()->SetLabelFont(42);
   Graph_grWpUncorrelatedSystc2->GetZaxis()->SetLabelSize(0.05);
   Graph_grWpUncorrelatedSystc2->GetZaxis()->SetTitleSize(0.05);
   Graph_grWpUncorrelatedSystc2->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_grWpUncorrelatedSystc2);
   
   grae->Draw("e2");
   
   grae = new TGraphAsymmErrors(9);
   grae->SetName("grWpCorrelatedSystc");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ff0000");
   grae->SetFillColor(ci);
   grae->SetFillStyle(3004);

   ci = TColor::GetColor("#cc0000");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,0.225,8.701805);
   grae->SetPointError(0,0.125,0.125,0.2594396,0.2594396);
   grae->SetPoint(1,0.475,7.649763);
   grae->SetPointError(1,0.125,0.125,0.2714556,0.2714556);
   grae->SetPoint(2,0.7,7.896176);
   grae->SetPointError(2,0.1,0.1,0.2244766,0.2244766);
   grae->SetPoint(3,0.925,7.436997);
   grae->SetPointError(3,0.125,0.125,0.2403047,0.2403047);
   grae->SetPoint(4,1.175,6.334746);
   grae->SetPointError(4,0.125,0.125,0.2743154,0.2743154);
   grae->SetPoint(5,1.425,5.989713);
   grae->SetPointError(5,0.125,0.125,0.2085459,0.2085459);
   grae->SetPoint(6,1.7,4.934416);
   grae->SetPointError(6,0.15,0.15,0.2222889,0.2222889);
   grae->SetPoint(7,1.975,4.234815);
   grae->SetPointError(7,0.125,0.125,0.1486127,0.1486127);
   grae->SetPoint(8,2.25,3.747698);
   grae->SetPointError(8,0.15,0.15,0.1801444,0.1801444);
   
   TH1F *Graph_grWpCorrelatedSystc3 = new TH1F("Graph_grWpCorrelatedSystc3","Graph",100,0,2.63);
   Graph_grWpCorrelatedSystc3->SetMinimum(3.028184);
   Graph_grWpCorrelatedSystc3->SetMaximum(9.500614);
   Graph_grWpCorrelatedSystc3->SetDirectory(0);
   Graph_grWpCorrelatedSystc3->SetStats(0);
   Graph_grWpCorrelatedSystc3->SetLineWidth(2);
   Graph_grWpCorrelatedSystc3->SetMarkerStyle(20);
   Graph_grWpCorrelatedSystc3->SetMarkerSize(1.2);
   Graph_grWpCorrelatedSystc3->GetXaxis()->SetLabelFont(42);
   Graph_grWpCorrelatedSystc3->GetXaxis()->SetLabelSize(0.05);
   Graph_grWpCorrelatedSystc3->GetXaxis()->SetTitleSize(0.05);
   Graph_grWpCorrelatedSystc3->GetXaxis()->SetTitleOffset(1.4);
   Graph_grWpCorrelatedSystc3->GetXaxis()->SetTitleFont(42);
   Graph_grWpCorrelatedSystc3->GetYaxis()->SetLabelFont(42);
   Graph_grWpCorrelatedSystc3->GetYaxis()->SetLabelSize(0.05);
   Graph_grWpCorrelatedSystc3->GetYaxis()->SetTitleSize(0.05);
   Graph_grWpCorrelatedSystc3->GetYaxis()->SetTitleOffset(1.4);
   Graph_grWpCorrelatedSystc3->GetYaxis()->SetTitleFont(42);
   Graph_grWpCorrelatedSystc3->GetZaxis()->SetLabelFont(42);
   Graph_grWpCorrelatedSystc3->GetZaxis()->SetLabelSize(0.05);
   Graph_grWpCorrelatedSystc3->GetZaxis()->SetTitleSize(0.05);
   Graph_grWpCorrelatedSystc3->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_grWpCorrelatedSystc3);
   
   grae->Draw("e2");
   
   grae = new TGraphAsymmErrors(9);
   grae->SetName("grWpc");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#cc0000");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#cc0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(23);
   grae->SetMarkerSize(1.5);
   grae->SetPoint(0,0.225,8.701805);
   grae->SetPointError(0,0.125,0.125,0.3783174,0.3783174);
   grae->SetPoint(1,0.475,7.649763);
   grae->SetPointError(1,0.125,0.125,0.3463731,0.3463731);
   grae->SetPoint(2,0.7,7.896176);
   grae->SetPointError(2,0.1,0.1,0.3840557,0.3840557);
   grae->SetPoint(3,0.925,7.436997);
   grae->SetPointError(3,0.125,0.125,0.3345559,0.3345559);
   grae->SetPoint(4,1.175,6.334746);
   grae->SetPointError(4,0.125,0.125,0.3451581,0.3451581);
   grae->SetPoint(5,1.425,5.989713);
   grae->SetPointError(5,0.125,0.125,0.3222833,0.3222833);
   grae->SetPoint(6,1.7,4.934416);
   grae->SetPointError(6,0.15,0.15,0.2775243,0.2775243);
   grae->SetPoint(7,1.975,4.234815);
   grae->SetPointError(7,0.125,0.125,0.2985153,0.2985153);
   grae->SetPoint(8,2.25,3.747698);
   grae->SetPointError(8,0.15,0.15,0.287349,0.287349);
   
   TH1F *Graph_Graph_grWpc14 = new TH1F("Graph_Graph_grWpc14","Graph",100,0,2.63);
   Graph_Graph_grWpc14->SetMinimum(0);
   Graph_Graph_grWpc14->SetMaximum(10);
   Graph_Graph_grWpc14->SetDirectory(0);
   Graph_Graph_grWpc14->SetStats(0);
   Graph_Graph_grWpc14->SetLineWidth(2);
   Graph_Graph_grWpc14->SetMarkerStyle(20);
   Graph_Graph_grWpc14->SetMarkerSize(1.2);
   Graph_Graph_grWpc14->GetXaxis()->SetTitle("|#eta|");
   Graph_Graph_grWpc14->GetXaxis()->SetRange(4,92);
   Graph_Graph_grWpc14->GetXaxis()->SetLabelFont(42);
   Graph_Graph_grWpc14->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph_grWpc14->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph_grWpc14->GetXaxis()->SetTitleFont(42);
   Graph_Graph_grWpc14->GetYaxis()->SetTitle("#frac{dN}{d#eta} #frac{1}{N_{events} } #frac{10^{9}}{#LT N_{coll} #GT}");
   Graph_Graph_grWpc14->GetYaxis()->SetLabelFont(42);
   Graph_Graph_grWpc14->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph_grWpc14->GetYaxis()->SetTitleFont(42);
   Graph_Graph_grWpc14->GetZaxis()->SetLabelFont(42);
   Graph_Graph_grWpc14->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_grWpc14->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph_grWpc14->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_grWpc14);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(9);
   grae->SetName("grWpc");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#cc0000");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#cc0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(23);
   grae->SetMarkerSize(1.5);
   grae->SetPoint(0,0.225,8.701805);
   grae->SetPointError(0,0.125,0.125,0.3783174,0.3783174);
   grae->SetPoint(1,0.475,7.649763);
   grae->SetPointError(1,0.125,0.125,0.3463731,0.3463731);
   grae->SetPoint(2,0.7,7.896176);
   grae->SetPointError(2,0.1,0.1,0.3840557,0.3840557);
   grae->SetPoint(3,0.925,7.436997);
   grae->SetPointError(3,0.125,0.125,0.3345559,0.3345559);
   grae->SetPoint(4,1.175,6.334746);
   grae->SetPointError(4,0.125,0.125,0.3451581,0.3451581);
   grae->SetPoint(5,1.425,5.989713);
   grae->SetPointError(5,0.125,0.125,0.3222833,0.3222833);
   grae->SetPoint(6,1.7,4.934416);
   grae->SetPointError(6,0.15,0.15,0.2775243,0.2775243);
   grae->SetPoint(7,1.975,4.234815);
   grae->SetPointError(7,0.125,0.125,0.2985153,0.2985153);
   grae->SetPoint(8,2.25,3.747698);
   grae->SetPointError(8,0.15,0.15,0.287349,0.287349);
   
   TH1F *Graph_Graph_Graph_grWpc145 = new TH1F("Graph_Graph_Graph_grWpc145","Graph",100,0,2.63);
   Graph_Graph_Graph_grWpc145->SetMinimum(0);
   Graph_Graph_Graph_grWpc145->SetMaximum(10);
   Graph_Graph_Graph_grWpc145->SetDirectory(0);
   Graph_Graph_Graph_grWpc145->SetStats(0);
   Graph_Graph_Graph_grWpc145->SetLineWidth(2);
   Graph_Graph_Graph_grWpc145->SetMarkerStyle(20);
   Graph_Graph_Graph_grWpc145->SetMarkerSize(1.2);
   Graph_Graph_Graph_grWpc145->GetXaxis()->SetTitle("|#eta|");
   Graph_Graph_Graph_grWpc145->GetXaxis()->SetRange(4,92);
   Graph_Graph_Graph_grWpc145->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph_grWpc145->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph_Graph_grWpc145->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph_Graph_grWpc145->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph_grWpc145->GetYaxis()->SetTitle("#frac{dN}{d#eta} #frac{1}{N_{events} } #frac{10^{9}}{#LT N_{coll} #GT}");
   Graph_Graph_Graph_grWpc145->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph_grWpc145->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph_Graph_grWpc145->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph_grWpc145->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph_grWpc145->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_Graph_grWpc145->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph_Graph_grWpc145->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph_grWpc145);
   
   grae->Draw("pe");
   
   TLegend *leg = new TLegend(0.192953,0.2027972,0.5067114,0.3793706,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04545455);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("grWpc","Data","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#cc0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.5);
   entry=leg->AddEntry("NULL","nn","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#00cccc");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(2);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("NULL","np(pn)","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#cc00ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(2);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("NULL","pp","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#00ff00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(2);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("NULL","All","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(27);
   entry->SetMarkerSize(1.7);
   leg->Draw();
   TLatex *   tex = new TLatex(0.4244966,0.2377622,"W^{-}#rightarrow#mu^{-}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.06993007);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.4211409,0.3129371,"POW+PY8 NLO");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
   cEtaPlus->Modified();
   cEtaPlus->cd();
   cEtaPlus->SetSelected(cEtaPlus);
}
