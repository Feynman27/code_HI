{
//=========Macro generated from canvas: cEtaMinus/cEtaMinus
//=========  (Sun Aug 11 21:46:17 2013) by ROOT version5.32/00
   TCanvas *cEtaMinus = new TCanvas("cEtaMinus", "cEtaMinus",733,64,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cEtaMinus->Range(-0.3951658,-2.025316,2.567746,10.63291);
   cEtaMinus->SetFillColor(0);
   cEtaMinus->SetBorderMode(0);
   cEtaMinus->SetBorderSize(2);
   cEtaMinus->SetTickx(1);
   cEtaMinus->SetTicky(1);
   cEtaMinus->SetLeftMargin(0.16);
   cEtaMinus->SetRightMargin(0.05);
   cEtaMinus->SetTopMargin(0.05);
   cEtaMinus->SetBottomMargin(0.16);
   cEtaMinus->SetFrameBorderMode(0);
   cEtaMinus->SetFrameBorderMode(0);
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(9);
   grae->SetName("grWmc");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#0000ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(22);
   grae->SetMarkerSize(1.5);
   grae->SetPoint(0,0.225,7.223949);
   grae->SetPointError(0,0.125,0.125,0.3539886,0.3539886);
   grae->SetPoint(1,0.475,7.266781);
   grae->SetPointError(1,0.125,0.125,0.3453692,0.3453692);
   grae->SetPoint(2,0.7,6.404811);
   grae->SetPointError(2,0.1,0.1,0.3535848,0.3535848);
   grae->SetPoint(3,0.925,6.357371);
   grae->SetPointError(3,0.125,0.125,0.3125256,0.3125256);
   grae->SetPoint(4,1.175,5.362904);
   grae->SetPointError(4,0.125,0.125,0.3182656,0.3182656);
   grae->SetPoint(5,1.425,6.333934);
   grae->SetPointError(5,0.125,0.125,0.3362998,0.3362998);
   grae->SetPoint(6,1.7,5.164461);
   grae->SetPointError(6,0.15,0.15,0.2815426,0.2815426);
   grae->SetPoint(7,1.975,5.520734);
   grae->SetPointError(7,0.125,0.125,0.3472644,0.3472644);
   grae->SetPoint(8,2.25,7.064405);
   grae->SetPointError(8,0.15,0.15,0.3996699,0.3996699);
   
   TH1F *Graph_grWmc6 = new TH1F("Graph_grWmc6","Graph",100,0,2.63);
   Graph_grWmc6->SetMinimum(0);
   Graph_grWmc6->SetMaximum(10);
   Graph_grWmc6->SetDirectory(0);
   Graph_grWmc6->SetStats(0);
   Graph_grWmc6->SetLineWidth(2);
   Graph_grWmc6->SetMarkerStyle(20);
   Graph_grWmc6->SetMarkerSize(1.2);
   Graph_grWmc6->GetXaxis()->SetTitle("|#eta|");
   Graph_grWmc6->GetXaxis()->SetRange(4,92);
   Graph_grWmc6->GetXaxis()->SetLabelFont(42);
   Graph_grWmc6->GetXaxis()->SetTitleSize(0.05);
   Graph_grWmc6->GetXaxis()->SetTitleOffset(1.4);
   Graph_grWmc6->GetXaxis()->SetTitleFont(42);
   Graph_grWmc6->GetYaxis()->SetTitle("#frac{dN}{d#eta} #frac{1}{N_{events} } #frac{10^{9}}{#LT N_{coll} #GT}");
   Graph_grWmc6->GetYaxis()->SetLabelFont(42);
   Graph_grWmc6->GetYaxis()->SetTitleOffset(1.4);
   Graph_grWmc6->GetYaxis()->SetTitleFont(42);
   Graph_grWmc6->GetZaxis()->SetLabelFont(42);
   Graph_grWmc6->GetZaxis()->SetLabelSize(0.05);
   Graph_grWmc6->GetZaxis()->SetTitleSize(0.05);
   Graph_grWmc6->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_grWmc6);
   
   grae->Draw("ape");
   
   grae = new TGraphAsymmErrors(9);
   grae->SetName("grWmUncorrelatedSystc");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#9999cc");
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#0000ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(22);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,0.225,7.223949);
   grae->SetPointError(0,0.125,0.125,0.4391351,0.4391351);
   grae->SetPoint(1,0.475,7.266781);
   grae->SetPointError(1,0.125,0.125,0.4398886,0.4398886);
   grae->SetPoint(2,0.7,6.404811);
   grae->SetPointError(2,0.1,0.1,0.4079118,0.4079118);
   grae->SetPoint(3,0.925,6.357371);
   grae->SetPointError(3,0.125,0.125,0.3829125,0.3829125);
   grae->SetPoint(4,1.175,5.362904);
   grae->SetPointError(4,0.125,0.125,0.3631035,0.3631035);
   grae->SetPoint(5,1.425,6.333934);
   grae->SetPointError(5,0.125,0.125,0.3908296,0.3908296);
   grae->SetPoint(6,1.7,5.164461);
   grae->SetPointError(6,0.15,0.15,0.3078254,0.3078254);
   grae->SetPoint(7,1.975,5.520734);
   grae->SetPointError(7,0.125,0.125,0.3698181,0.3698181);
   grae->SetPoint(8,2.25,7.064405);
   grae->SetPointError(8,0.15,0.15,0.4453209,0.4453209);
   
   TH1F *Graph_grWmUncorrelatedSystc7 = new TH1F("Graph_grWmUncorrelatedSystc7","Graph",100,0,2.63);
   Graph_grWmUncorrelatedSystc7->SetMinimum(4.571633);
   Graph_grWmUncorrelatedSystc7->SetMaximum(7.991674);
   Graph_grWmUncorrelatedSystc7->SetDirectory(0);
   Graph_grWmUncorrelatedSystc7->SetStats(0);
   Graph_grWmUncorrelatedSystc7->SetLineWidth(2);
   Graph_grWmUncorrelatedSystc7->SetMarkerStyle(20);
   Graph_grWmUncorrelatedSystc7->SetMarkerSize(1.2);
   Graph_grWmUncorrelatedSystc7->GetXaxis()->SetLabelFont(42);
   Graph_grWmUncorrelatedSystc7->GetXaxis()->SetLabelSize(0.05);
   Graph_grWmUncorrelatedSystc7->GetXaxis()->SetTitleSize(0.05);
   Graph_grWmUncorrelatedSystc7->GetXaxis()->SetTitleOffset(1.4);
   Graph_grWmUncorrelatedSystc7->GetXaxis()->SetTitleFont(42);
   Graph_grWmUncorrelatedSystc7->GetYaxis()->SetLabelFont(42);
   Graph_grWmUncorrelatedSystc7->GetYaxis()->SetLabelSize(0.05);
   Graph_grWmUncorrelatedSystc7->GetYaxis()->SetTitleSize(0.05);
   Graph_grWmUncorrelatedSystc7->GetYaxis()->SetTitleOffset(1.4);
   Graph_grWmUncorrelatedSystc7->GetYaxis()->SetTitleFont(42);
   Graph_grWmUncorrelatedSystc7->GetZaxis()->SetLabelFont(42);
   Graph_grWmUncorrelatedSystc7->GetZaxis()->SetLabelSize(0.05);
   Graph_grWmUncorrelatedSystc7->GetZaxis()->SetTitleSize(0.05);
   Graph_grWmUncorrelatedSystc7->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_grWmUncorrelatedSystc7);
   
   grae->Draw("e2");
   
   grae = new TGraphAsymmErrors(9);
   grae->SetName("grWmCorrelatedSystc");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#0000ff");
   grae->SetFillColor(ci);
   grae->SetFillStyle(3005);

   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,0.225,7.223949);
   grae->SetPointError(0,0.125,0.125,0.2324875,0.2324875);
   grae->SetPoint(1,0.475,7.266781);
   grae->SetPointError(1,0.125,0.125,0.2795504,0.2795504);
   grae->SetPoint(2,0.7,6.404811);
   grae->SetPointError(2,0.1,0.1,0.2464683,0.2464683);
   grae->SetPoint(3,0.925,6.357371);
   grae->SetPointError(3,0.125,0.125,0.2583279,0.2583279);
   grae->SetPoint(4,1.175,5.362904);
   grae->SetPointError(4,0.125,0.125,0.2457286,0.2457286);
   grae->SetPoint(5,1.425,6.333934);
   grae->SetPointError(5,0.125,0.125,0.2480184,0.2480184);
   grae->SetPoint(6,1.7,5.164461);
   grae->SetPointError(6,0.15,0.15,0.1801373,0.1801373);
   grae->SetPoint(7,1.975,5.520734);
   grae->SetPointError(7,0.125,0.125,0.1559647,0.1559647);
   grae->SetPoint(8,2.25,7.064405);
   grae->SetPointError(8,0.15,0.15,0.2182971,0.2182971);
   
   TH1F *Graph_grWmCorrelatedSystc8 = new TH1F("Graph_grWmCorrelatedSystc8","Graph",100,0,2.63);
   Graph_grWmCorrelatedSystc8->SetMinimum(4.728123);
   Graph_grWmCorrelatedSystc8->SetMaximum(7.802533);
   Graph_grWmCorrelatedSystc8->SetDirectory(0);
   Graph_grWmCorrelatedSystc8->SetStats(0);
   Graph_grWmCorrelatedSystc8->SetLineWidth(2);
   Graph_grWmCorrelatedSystc8->SetMarkerStyle(20);
   Graph_grWmCorrelatedSystc8->SetMarkerSize(1.2);
   Graph_grWmCorrelatedSystc8->GetXaxis()->SetLabelFont(42);
   Graph_grWmCorrelatedSystc8->GetXaxis()->SetLabelSize(0.05);
   Graph_grWmCorrelatedSystc8->GetXaxis()->SetTitleSize(0.05);
   Graph_grWmCorrelatedSystc8->GetXaxis()->SetTitleOffset(1.4);
   Graph_grWmCorrelatedSystc8->GetXaxis()->SetTitleFont(42);
   Graph_grWmCorrelatedSystc8->GetYaxis()->SetLabelFont(42);
   Graph_grWmCorrelatedSystc8->GetYaxis()->SetLabelSize(0.05);
   Graph_grWmCorrelatedSystc8->GetYaxis()->SetTitleSize(0.05);
   Graph_grWmCorrelatedSystc8->GetYaxis()->SetTitleOffset(1.4);
   Graph_grWmCorrelatedSystc8->GetYaxis()->SetTitleFont(42);
   Graph_grWmCorrelatedSystc8->GetZaxis()->SetLabelFont(42);
   Graph_grWmCorrelatedSystc8->GetZaxis()->SetLabelSize(0.05);
   Graph_grWmCorrelatedSystc8->GetZaxis()->SetTitleSize(0.05);
   Graph_grWmCorrelatedSystc8->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_grWmCorrelatedSystc8);
   
   grae->Draw("e2");
   
   grae = new TGraphAsymmErrors(9);
   grae->SetName("grWmc");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#0000ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(22);
   grae->SetMarkerSize(1.5);
   grae->SetPoint(0,0.225,7.223949);
   grae->SetPointError(0,0.125,0.125,0.3539886,0.3539886);
   grae->SetPoint(1,0.475,7.266781);
   grae->SetPointError(1,0.125,0.125,0.3453692,0.3453692);
   grae->SetPoint(2,0.7,6.404811);
   grae->SetPointError(2,0.1,0.1,0.3535848,0.3535848);
   grae->SetPoint(3,0.925,6.357371);
   grae->SetPointError(3,0.125,0.125,0.3125256,0.3125256);
   grae->SetPoint(4,1.175,5.362904);
   grae->SetPointError(4,0.125,0.125,0.3182656,0.3182656);
   grae->SetPoint(5,1.425,6.333934);
   grae->SetPointError(5,0.125,0.125,0.3362998,0.3362998);
   grae->SetPoint(6,1.7,5.164461);
   grae->SetPointError(6,0.15,0.15,0.2815426,0.2815426);
   grae->SetPoint(7,1.975,5.520734);
   grae->SetPointError(7,0.125,0.125,0.3472644,0.3472644);
   grae->SetPoint(8,2.25,7.064405);
   grae->SetPointError(8,0.15,0.15,0.3996699,0.3996699);
   
   TH1F *Graph_Graph_grWmc69 = new TH1F("Graph_Graph_grWmc69","Graph",100,0,2.63);
   Graph_Graph_grWmc69->SetMinimum(0);
   Graph_Graph_grWmc69->SetMaximum(10);
   Graph_Graph_grWmc69->SetDirectory(0);
   Graph_Graph_grWmc69->SetStats(0);
   Graph_Graph_grWmc69->SetLineWidth(2);
   Graph_Graph_grWmc69->SetMarkerStyle(20);
   Graph_Graph_grWmc69->SetMarkerSize(1.2);
   Graph_Graph_grWmc69->GetXaxis()->SetTitle("|#eta|");
   Graph_Graph_grWmc69->GetXaxis()->SetRange(4,92);
   Graph_Graph_grWmc69->GetXaxis()->SetLabelFont(42);
   Graph_Graph_grWmc69->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph_grWmc69->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph_grWmc69->GetXaxis()->SetTitleFont(42);
   Graph_Graph_grWmc69->GetYaxis()->SetTitle("#frac{dN}{d#eta} #frac{1}{N_{events} } #frac{10^{9}}{#LT N_{coll} #GT}");
   Graph_Graph_grWmc69->GetYaxis()->SetLabelFont(42);
   Graph_Graph_grWmc69->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph_grWmc69->GetYaxis()->SetTitleFont(42);
   Graph_Graph_grWmc69->GetZaxis()->SetLabelFont(42);
   Graph_Graph_grWmc69->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_grWmc69->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph_grWmc69->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_grWmc69);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(9);
   grae->SetName("grWmc");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#0000ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(22);
   grae->SetMarkerSize(1.5);
   grae->SetPoint(0,0.225,7.223949);
   grae->SetPointError(0,0.125,0.125,0.3539886,0.3539886);
   grae->SetPoint(1,0.475,7.266781);
   grae->SetPointError(1,0.125,0.125,0.3453692,0.3453692);
   grae->SetPoint(2,0.7,6.404811);
   grae->SetPointError(2,0.1,0.1,0.3535848,0.3535848);
   grae->SetPoint(3,0.925,6.357371);
   grae->SetPointError(3,0.125,0.125,0.3125256,0.3125256);
   grae->SetPoint(4,1.175,5.362904);
   grae->SetPointError(4,0.125,0.125,0.3182656,0.3182656);
   grae->SetPoint(5,1.425,6.333934);
   grae->SetPointError(5,0.125,0.125,0.3362998,0.3362998);
   grae->SetPoint(6,1.7,5.164461);
   grae->SetPointError(6,0.15,0.15,0.2815426,0.2815426);
   grae->SetPoint(7,1.975,5.520734);
   grae->SetPointError(7,0.125,0.125,0.3472644,0.3472644);
   grae->SetPoint(8,2.25,7.064405);
   grae->SetPointError(8,0.15,0.15,0.3996699,0.3996699);
   
   TH1F *Graph_Graph_Graph_grWmc6910 = new TH1F("Graph_Graph_Graph_grWmc6910","Graph",100,0,2.63);
   Graph_Graph_Graph_grWmc6910->SetMinimum(0);
   Graph_Graph_Graph_grWmc6910->SetMaximum(10);
   Graph_Graph_Graph_grWmc6910->SetDirectory(0);
   Graph_Graph_Graph_grWmc6910->SetStats(0);
   Graph_Graph_Graph_grWmc6910->SetLineWidth(2);
   Graph_Graph_Graph_grWmc6910->SetMarkerStyle(20);
   Graph_Graph_Graph_grWmc6910->SetMarkerSize(1.2);
   Graph_Graph_Graph_grWmc6910->GetXaxis()->SetTitle("|#eta|");
   Graph_Graph_Graph_grWmc6910->GetXaxis()->SetRange(4,92);
   Graph_Graph_Graph_grWmc6910->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph_grWmc6910->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph_Graph_grWmc6910->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph_Graph_grWmc6910->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph_grWmc6910->GetYaxis()->SetTitle("#frac{dN}{d#eta} #frac{1}{N_{events} } #frac{10^{9}}{#LT N_{coll} #GT}");
   Graph_Graph_Graph_grWmc6910->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph_grWmc6910->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph_Graph_grWmc6910->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph_grWmc6910->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph_grWmc6910->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_Graph_grWmc6910->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph_Graph_grWmc6910->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph_grWmc6910);
   
   grae->Draw("pe");
   
   TLegend *leg = new TLegend(0.192953,0.2027972,0.5067114,0.3793706,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04545455);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("grWmc","Data","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
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

   ci = TColor::GetColor("#0000ff");
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
   cEtaMinus->Modified();
   cEtaMinus->cd();
   cEtaMinus->SetSelected(cEtaMinus);
}
