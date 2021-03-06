{
//=========Macro generated from canvas: c3/c3
//=========  (Thu Oct  3 19:25:05 2013) by ROOT version5.34/09
   TCanvas *c3 = new TCanvas("c3", "c3",2379,133,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c3->Range(-0.3639066,-0.5241593,2.546497,0.2477876);
   c3->SetFillColor(0);
   c3->SetBorderMode(0);
   c3->SetBorderSize(2);
   c3->SetTickx(1);
   c3->SetTicky(1);
   c3->SetLeftMargin(0.159396);
   c3->SetRightMargin(0.05033557);
   c3->SetTopMargin(0.04895105);
   c3->SetBottomMargin(0.1608392);
   c3->SetFrameBorderMode(0);
   c3->SetFrameBorderMode(0);
   Double_t xAxis1[10] = {0.1, 0.35, 0.6, 0.8, 1.05, 1.3, 1.55, 1.85, 2.1, 2.4}; 
   
   TH1D *hassFit_cent0__1 = new TH1D("hassFit_cent0__1","hassFit_cent0",9, xAxis1);
   hassFit_cent0__1->SetMinimum(-0.4);
   hassFit_cent0__1->SetMaximum(0.21);
   hassFit_cent0__1->SetDirectory(0);
   hassFit_cent0__1->SetStats(0);
   hassFit_cent0__1->SetLineWidth(2);
   hassFit_cent0__1->SetMarkerStyle(20);
   hassFit_cent0__1->SetMarkerSize(1.2);
   hassFit_cent0__1->GetXaxis()->SetTitle("|#eta_{#mu}|");
   hassFit_cent0__1->GetXaxis()->SetRange(1,9);
   hassFit_cent0__1->GetXaxis()->SetLabelFont(42);
   hassFit_cent0__1->GetXaxis()->SetTitleSize(0.05);
   hassFit_cent0__1->GetXaxis()->SetTitleOffset(1.4);
   hassFit_cent0__1->GetXaxis()->SetTitleFont(42);
   hassFit_cent0__1->GetYaxis()->SetTitle("A_{#mu}");
   hassFit_cent0__1->GetYaxis()->SetNoExponent();
   hassFit_cent0__1->GetYaxis()->SetLabelFont(42);
   hassFit_cent0__1->GetYaxis()->SetTitleSize(0.05);
   hassFit_cent0__1->GetYaxis()->SetTitleOffset(1.4);
   hassFit_cent0__1->GetYaxis()->SetTitleFont(42);
   hassFit_cent0__1->GetZaxis()->SetLabelFont(42);
   hassFit_cent0__1->GetZaxis()->SetLabelSize(0.05);
   hassFit_cent0__1->GetZaxis()->SetTitleSize(0.05);
   hassFit_cent0__1->GetZaxis()->SetTitleFont(42);
   hassFit_cent0__1->Draw("");
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(9);
   grae->SetName("Graph");
   grae->SetTitle("Graph");
   grae->SetFillColor(38);
   grae->SetFillStyle(3004);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,0.225,0.1202851);
   grae->SetPointError(0,0.125,0.125,0.04458126,0.04458126);
   grae->SetPoint(1,0.475,0.03544523);
   grae->SetPointError(1,0.125,0.125,0.0445474,0.0445474);
   grae->SetPoint(2,0.7,0.1333614);
   grae->SetPointError(2,0.1,0.1,0.04393111,0.04393111);
   grae->SetPoint(3,0.925,0.09142858);
   grae->SetPointError(3,0.125,0.125,0.04272781,0.04272781);
   grae->SetPoint(4,1.175,0.1013146);
   grae->SetPointError(4,0.125,0.125,0.04815646,0.04815646);
   grae->SetPoint(5,1.425,0.003587187);
   grae->SetPointError(5,0.125,0.125,0.04458778,0.04458778);
   grae->SetPoint(6,1.7,-0.01331378);
   grae->SetPointError(6,0.15,0.15,0.0520316,0.0520316);
   grae->SetPoint(7,1.975,-0.1063848);
   grae->SetPointError(7,0.125,0.125,0.05049361,0.05049361);
   grae->SetPoint(8,2.25,-0.3011886);
   grae->SetPointError(8,0.15,0.15,0.0537708,0.0537708);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,0,2.63);
   Graph_Graph1->SetMinimum(-0.4081846);
   Graph_Graph1->SetMaximum(0.2305177);
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
   
   grae->Draw("2");
   Double_t xAxis2[10] = {0.1, 0.35, 0.6, 0.8, 1.05, 1.3, 1.55, 1.85, 2.1, 2.4}; 
   
   TH1F *hEtaPlusMcAll__2 = new TH1F("hEtaPlusMcAll__2","hEtaPlusMcAll",9, xAxis2);
   hEtaPlusMcAll__2->SetBinContent(1,0.0860365);
   hEtaPlusMcAll__2->SetBinContent(2,0.09523783);
   hEtaPlusMcAll__2->SetBinContent(3,0.09096931);
   hEtaPlusMcAll__2->SetBinContent(4,0.07248296);
   hEtaPlusMcAll__2->SetBinContent(5,0.05074221);
   hEtaPlusMcAll__2->SetBinContent(6,0.02617276);
   hEtaPlusMcAll__2->SetBinContent(7,-0.04210219);
   hEtaPlusMcAll__2->SetBinContent(8,-0.1277015);
   hEtaPlusMcAll__2->SetBinContent(9,-0.2514366);
   hEtaPlusMcAll__2->SetBinError(1,0.005035776);
   hEtaPlusMcAll__2->SetBinError(2,0.005049663);
   hEtaPlusMcAll__2->SetBinError(3,0.00570708);
   hEtaPlusMcAll__2->SetBinError(4,0.005178113);
   hEtaPlusMcAll__2->SetBinError(5,0.005270762);
   hEtaPlusMcAll__2->SetBinError(6,0.005443933);
   hEtaPlusMcAll__2->SetBinError(7,0.005198336);
   hEtaPlusMcAll__2->SetBinError(8,0.00607069);
   hEtaPlusMcAll__2->SetBinError(9,0.006008842);
   hEtaPlusMcAll__2->SetEntries(0.0005840431);
   hEtaPlusMcAll__2->SetDirectory(0);
   hEtaPlusMcAll__2->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#0000ff");
   hEtaPlusMcAll__2->SetFillColor(ci);
   hEtaPlusMcAll__2->SetFillStyle(3005);

   ci = TColor::GetColor("#0000ff");
   hEtaPlusMcAll__2->SetMarkerColor(ci);
   hEtaPlusMcAll__2->SetMarkerStyle(27);
   hEtaPlusMcAll__2->SetMarkerSize(1.9);
   hEtaPlusMcAll__2->GetXaxis()->SetLabelFont(42);
   hEtaPlusMcAll__2->GetXaxis()->SetLabelSize(0.05);
   hEtaPlusMcAll__2->GetXaxis()->SetTitleSize(0.05);
   hEtaPlusMcAll__2->GetXaxis()->SetTitleOffset(1.4);
   hEtaPlusMcAll__2->GetXaxis()->SetTitleFont(42);
   hEtaPlusMcAll__2->GetYaxis()->SetLabelFont(42);
   hEtaPlusMcAll__2->GetYaxis()->SetLabelSize(0.05);
   hEtaPlusMcAll__2->GetYaxis()->SetTitleSize(0.05);
   hEtaPlusMcAll__2->GetYaxis()->SetTitleOffset(1.4);
   hEtaPlusMcAll__2->GetYaxis()->SetTitleFont(42);
   hEtaPlusMcAll__2->GetZaxis()->SetLabelFont(42);
   hEtaPlusMcAll__2->GetZaxis()->SetLabelSize(0.05);
   hEtaPlusMcAll__2->GetZaxis()->SetTitleSize(0.05);
   hEtaPlusMcAll__2->GetZaxis()->SetTitleFont(42);
   hEtaPlusMcAll__2->Draw("pe2same");
   Double_t xAxis3[10] = {0.1, 0.35, 0.6, 0.8, 1.05, 1.3, 1.55, 1.85, 2.1, 2.4}; 
   
   TH1F *hEtaPlusMSTW__3 = new TH1F("hEtaPlusMSTW__3","hEtaPlusMSTW",9, xAxis3);
   hEtaPlusMSTW__3->SetBinContent(1,0.1018884);
   hEtaPlusMSTW__3->SetBinContent(2,0.09499385);
   hEtaPlusMSTW__3->SetBinContent(3,0.09110907);
   hEtaPlusMSTW__3->SetBinContent(4,0.08350987);
   hEtaPlusMSTW__3->SetBinContent(5,0.06339273);
   hEtaPlusMSTW__3->SetBinContent(6,0.03476641);
   hEtaPlusMSTW__3->SetBinContent(7,-0.02327944);
   hEtaPlusMSTW__3->SetBinContent(8,-0.1224817);
   hEtaPlusMSTW__3->SetBinContent(9,-0.2472788);
   hEtaPlusMSTW__3->SetBinError(1,0.00175888);
   hEtaPlusMSTW__3->SetBinError(2,0.001761619);
   hEtaPlusMSTW__3->SetBinError(3,0.00198385);
   hEtaPlusMSTW__3->SetBinError(4,0.00178604);
   hEtaPlusMSTW__3->SetBinError(5,0.001816671);
   hEtaPlusMSTW__3->SetBinError(6,0.001863759);
   hEtaPlusMSTW__3->SetBinError(7,0.001775201);
   hEtaPlusMSTW__3->SetBinError(8,0.002057001);
   hEtaPlusMSTW__3->SetBinError(9,0.002017936);
   hEtaPlusMSTW__3->SetEntries(180.387);
   hEtaPlusMSTW__3->SetDirectory(0);
   hEtaPlusMSTW__3->SetStats(0);

   ci = TColor::GetColor("#ff0000");
   hEtaPlusMSTW__3->SetFillColor(ci);
   hEtaPlusMSTW__3->SetFillStyle(3006);

   ci = TColor::GetColor("#ff0000");
   hEtaPlusMSTW__3->SetMarkerColor(ci);
   hEtaPlusMSTW__3->SetMarkerStyle(33);
   hEtaPlusMSTW__3->SetMarkerSize(1.9);
   hEtaPlusMSTW__3->GetXaxis()->SetLabelFont(42);
   hEtaPlusMSTW__3->GetXaxis()->SetLabelSize(0.05);
   hEtaPlusMSTW__3->GetXaxis()->SetTitleSize(0.05);
   hEtaPlusMSTW__3->GetXaxis()->SetTitleOffset(1.4);
   hEtaPlusMSTW__3->GetXaxis()->SetTitleFont(42);
   hEtaPlusMSTW__3->GetYaxis()->SetLabelFont(42);
   hEtaPlusMSTW__3->GetYaxis()->SetLabelSize(0.05);
   hEtaPlusMSTW__3->GetYaxis()->SetTitleSize(0.05);
   hEtaPlusMSTW__3->GetYaxis()->SetTitleOffset(1.4);
   hEtaPlusMSTW__3->GetYaxis()->SetTitleFont(42);
   hEtaPlusMSTW__3->GetZaxis()->SetLabelFont(42);
   hEtaPlusMSTW__3->GetZaxis()->SetLabelSize(0.05);
   hEtaPlusMSTW__3->GetZaxis()->SetTitleSize(0.05);
   hEtaPlusMSTW__3->GetZaxis()->SetTitleFont(42);
   hEtaPlusMSTW__3->Draw("pe2same");
   
   grae = new TGraphAsymmErrors(9);
   grae->SetName("Graph");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetLineWidth(3);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,0.225,0.1202851);
   grae->SetPointError(0,0.125,0.125,0.03274769,0.03274769);
   grae->SetPoint(1,0.475,0.03544523);
   grae->SetPointError(1,0.125,0.125,0.03280151,0.03280151);
   grae->SetPoint(2,0.7,0.1333614);
   grae->SetPointError(2,0.1,0.1,0.03684502,0.03684502);
   grae->SetPoint(3,0.925,0.09142858);
   grae->SetPointError(3,0.125,0.125,0.03328536,0.03328536);
   grae->SetPoint(4,1.175,0.1013146);
   grae->SetPointError(4,0.125,0.125,0.04029875,0.04029875);
   grae->SetPoint(5,1.425,0.003587187);
   grae->SetPointError(5,0.125,0.125,0.0377735,0.0377735);
   grae->SetPoint(6,1.7,-0.01331378);
   grae->SetPointError(6,0.15,0.15,0.03912068,0.03912068);
   grae->SetPoint(7,1.975,-0.1063848);
   grae->SetPointError(7,0.125,0.125,0.04713353,0.04713353);
   grae->SetPoint(8,2.25,-0.3011886);
   grae->SetPointError(8,0.15,0.15,0.04774354,0.04774354);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,0,2.63);
   Graph_Graph2->SetMinimum(-0.400846);
   Graph_Graph2->SetMaximum(0.2221203);
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
   Double_t xAxis4[10] = {0.1, 0.35, 0.6, 0.8, 1.05, 1.3, 1.55, 1.85, 2.1, 2.4}; 
   
   TH1D *hassFit_cent0__4 = new TH1D("hassFit_cent0__4","hassFit_cent0",9, xAxis4);
   hassFit_cent0__4->SetMinimum(-0.4);
   hassFit_cent0__4->SetMaximum(0.21);
   hassFit_cent0__4->SetDirectory(0);
   hassFit_cent0__4->SetStats(0);
   hassFit_cent0__4->SetLineWidth(2);
   hassFit_cent0__4->SetMarkerStyle(20);
   hassFit_cent0__4->SetMarkerSize(1.2);
   hassFit_cent0__4->GetXaxis()->SetTitle("|#eta_{#mu}|");
   hassFit_cent0__4->GetXaxis()->SetRange(0,9);
   hassFit_cent0__4->GetXaxis()->SetLabelFont(42);
   hassFit_cent0__4->GetXaxis()->SetTitleSize(0.05);
   hassFit_cent0__4->GetXaxis()->SetTitleOffset(1.4);
   hassFit_cent0__4->GetXaxis()->SetTitleFont(42);
   hassFit_cent0__4->GetYaxis()->SetTitle("A_{#mu}");
   hassFit_cent0__4->GetYaxis()->SetNoExponent();
   hassFit_cent0__4->GetYaxis()->SetLabelFont(42);
   hassFit_cent0__4->GetYaxis()->SetTitleSize(0.05);
   hassFit_cent0__4->GetYaxis()->SetTitleOffset(1.4);
   hassFit_cent0__4->GetYaxis()->SetTitleFont(42);
   hassFit_cent0__4->GetZaxis()->SetLabelFont(42);
   hassFit_cent0__4->GetZaxis()->SetLabelSize(0.05);
   hassFit_cent0__4->GetZaxis()->SetTitleSize(0.05);
   hassFit_cent0__4->GetZaxis()->SetTitleFont(42);
   hassFit_cent0__4->Draw("sameaxig");
   
   TLegend *leg = new TLegend(0.1661074,0.1713287,0.5352349,0.4073427,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03846154);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph","Data 2011","PEF");
   entry->SetFillColor(38);
   entry->SetFillStyle(3004);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(62);
   entry=leg->AddEntry("hEtaPlusMcAll","MRST LO*","pef");

   ci = TColor::GetColor("#0000ff");
   entry->SetFillColor(ci);
   entry->SetFillStyle(3005);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(27);
   entry->SetMarkerSize(1.9);
   entry->SetTextFont(62);
   entry=leg->AddEntry("hEtaPlusMSTW","MSTW NLO","pef");

   ci = TColor::GetColor("#ff0000");
   entry->SetFillColor(ci);
   entry->SetFillStyle(3006);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(1.9);
   entry->SetTextFont(62);
   leg->Draw();
   
   leg = new TLegend(0.4563758,0.2062937,0.7466443,0.3811189,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03846154);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("hAsymmDummyStat","Stat. uncertainty","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(2);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph","Total uncertainty","f");
   entry->SetFillColor(38);
   entry->SetFillStyle(3004);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   leg->Draw();
   TLatex *   tex = new TLatex(0.1828859,0.4335664,"#int Ldt #approx 0.14 nb^{-1}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03846154);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.4731544,0.4370629,"#sqrt{s_{NN}} = 2.76 TeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04195804);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.6073826,0.8566434,"ATLAS");
tex->SetNDC();
   tex->SetTextFont(72);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.7684564,0.8566434,"Internal");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
   c3->Modified();
   c3->cd();
   c3->SetSelected(c3);
}
