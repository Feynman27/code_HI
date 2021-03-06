{
//=========Macro generated from canvas: cRcp/cRcp
//=========  (Sun Nov  3 14:26:50 2013) by ROOT version5.34/10
   TCanvas *cRcp = new TCanvas("cRcp", "cRcp",784,52,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cRcp->Range(0,0,1,1);
   cRcp->SetFillColor(0);
   cRcp->SetBorderMode(0);
   cRcp->SetBorderSize(2);
   cRcp->SetTickx(1);
   cRcp->SetTicky(1);
   cRcp->SetLeftMargin(0.16);
   cRcp->SetRightMargin(0.05);
   cRcp->SetTopMargin(0.05);
   cRcp->SetBottomMargin(0.16);
   cRcp->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: pad2
   TPad *pad2 = new TPad("pad2", "",0,0,1,1);
   pad2->Draw();
   pad2->cd();
   pad2->Range(0,0,1,1);
   pad2->SetFillColor(0);
   pad2->SetBorderMode(0);
   pad2->SetBorderSize(2);
   pad2->SetTickx(1);
   pad2->SetTicky(1);
   pad2->SetLeftMargin(0.16);
   pad2->SetRightMargin(0.05);
   pad2->SetTopMargin(0.05);
   pad2->SetBottomMargin(0.16);
   pad2->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: pad1
   pad1 = new TPad("pad1", "",0,0,1,1);
   pad1->Draw();
   pad1->cd();
   pad1->Range(-85.06329,-7.493671,446.5823,39.34177);
   pad1->SetFillColor(0);
   pad1->SetFillStyle(4000);
   pad1->SetBorderMode(0);
   pad1->SetBorderSize(2);
   pad1->SetTickx(1);
   pad1->SetTicky(1);
   pad1->SetLeftMargin(0.16);
   pad1->SetRightMargin(0.05);
   pad1->SetTopMargin(0.05);
   pad1->SetBottomMargin(0.16);
   pad1->SetFrameBorderMode(0);
   pad1->SetFrameBorderMode(0);
   
   TH1F *hDummy__1 = new TH1F("hDummy__1","hDummy",100,0,420);
   hDummy__1->SetMinimum(0);
   hDummy__1->SetMaximum(37);
   hDummy__1->SetDirectory(0);
   hDummy__1->SetStats(0);
   hDummy__1->SetLineWidth(2);
   hDummy__1->SetMarkerStyle(20);
   hDummy__1->SetMarkerSize(1.2);
   hDummy__1->GetXaxis()->SetTitle("#LT N_{part} #GT");
   hDummy__1->GetXaxis()->SetLabelFont(42);
   hDummy__1->GetXaxis()->SetLabelSize(0.036);
   hDummy__1->GetXaxis()->SetTitleSize(0.05);
   hDummy__1->GetXaxis()->SetTitleOffset(1.4);
   hDummy__1->GetXaxis()->SetTitleFont(42);
   hDummy__1->GetYaxis()->SetTitle("#frac{10^{9}}{#LT N_{coll} #GT}#frac{N^{W#rightarrow#mu,fiducial}}{N_{events}}");
   hDummy__1->GetYaxis()->SetLabelFont(42);
   hDummy__1->GetYaxis()->SetLabelSize(0.03);
   hDummy__1->GetYaxis()->SetTitleOffset(1.4);
   hDummy__1->GetYaxis()->SetTitleFont(42);
   hDummy__1->GetZaxis()->SetLabelFont(42);
   hDummy__1->GetZaxis()->SetLabelSize(0.05);
   hDummy__1->GetZaxis()->SetTitleSize(0.05);
   hDummy__1->GetZaxis()->SetTitleFont(42);
   hDummy__1->Draw("0");
   
   TF1 *funcPythia = new TF1("funcPythia","[0]",0,420);
   funcPythia->SetFillColor(19);
   funcPythia->SetFillStyle(0);
   funcPythia->SetMarkerStyle(20);
   funcPythia->SetMarkerSize(1.2);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#666666");
   funcPythia->SetLineColor(ci);
   funcPythia->SetLineWidth(3);
   funcPythia->SetLineStyle(2);
   funcPythia->GetXaxis()->SetLabelFont(42);
   funcPythia->GetXaxis()->SetLabelSize(0.05);
   funcPythia->GetXaxis()->SetTitleSize(0.05);
   funcPythia->GetXaxis()->SetTitleOffset(1.4);
   funcPythia->GetXaxis()->SetTitleFont(42);
   funcPythia->GetYaxis()->SetLabelFont(42);
   funcPythia->GetYaxis()->SetLabelSize(0.05);
   funcPythia->GetYaxis()->SetTitleSize(0.05);
   funcPythia->GetYaxis()->SetTitleOffset(1.4);
   funcPythia->GetYaxis()->SetTitleFont(42);
   funcPythia->SetParameter(0,23.79885);
   funcPythia->SetParError(0,0);
   funcPythia->SetParLimits(0,23.79885,23.79885);
   funcPythia->Draw("same");
   
   TF1 *funcMSTWnlo = new TF1("funcMSTWnlo","[0]",0,420);
   funcMSTWnlo->SetFillColor(19);
   funcMSTWnlo->SetFillStyle(0);
   funcMSTWnlo->SetMarkerStyle(20);
   funcMSTWnlo->SetMarkerSize(1.2);

   ci = TColor::GetColor("#666666");
   funcMSTWnlo->SetLineColor(ci);
   funcMSTWnlo->SetLineWidth(3);
   funcMSTWnlo->SetLineStyle(9);
   funcMSTWnlo->GetXaxis()->SetLabelFont(42);
   funcMSTWnlo->GetXaxis()->SetLabelSize(0.05);
   funcMSTWnlo->GetXaxis()->SetTitleSize(0.05);
   funcMSTWnlo->GetXaxis()->SetTitleOffset(1.4);
   funcMSTWnlo->GetXaxis()->SetTitleFont(42);
   funcMSTWnlo->GetYaxis()->SetLabelFont(42);
   funcMSTWnlo->GetYaxis()->SetLabelSize(0.05);
   funcMSTWnlo->GetYaxis()->SetTitleSize(0.05);
   funcMSTWnlo->GetYaxis()->SetTitleOffset(1.4);
   funcMSTWnlo->GetYaxis()->SetTitleFont(42);
   funcMSTWnlo->SetParameter(0,27.71782);
   funcMSTWnlo->SetParError(0,0);
   funcMSTWnlo->SetParLimits(0,27.71782,27.71782);
   funcMSTWnlo->Draw("same");
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(6);
   grae->SetName("graphWSumc");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetMarkerStyle(33);
   grae->SetMarkerSize(1.8);
   grae->SetPoint(0,395.16,27.70695);
   grae->SetPointError(0,13.9108,13.9108,0.7727571,0.7727571);
   grae->SetPoint(1,343.26,27.28461);
   grae->SetPointError(1,14.97234,14.97234,0.8589446,0.8589446);
   grae->SetPoint(2,294.88,29.56996);
   grae->SetPointError(2,15.66444,15.66444,0.9792449,0.9792449);
   grae->SetPoint(3,252.52,30.57403);
   grae->SetPointError(3,15.83232,15.83232,1.12038,1.12038);
   grae->SetPoint(4,170.83,30.54106);
   grae->SetPointError(4,16.10358,16.10358,0.7212261,0.7212261);
   grae->SetPoint(5,58.93,27.13788);
   grae->SetPointError(5,14.34243,14.34243,1.038295,1.038295);
   
   TH1F *Graph_graphWSumc1 = new TH1F("Graph_graphWSumc1","Graph",100,8.139247,445.5191);
   Graph_graphWSumc1->SetMinimum(25.5401);
   Graph_graphWSumc1->SetMaximum(32.25389);
   Graph_graphWSumc1->SetDirectory(0);
   Graph_graphWSumc1->SetStats(0);
   Graph_graphWSumc1->SetLineWidth(2);
   Graph_graphWSumc1->SetMarkerStyle(20);
   Graph_graphWSumc1->SetMarkerSize(1.2);
   Graph_graphWSumc1->GetXaxis()->SetLabelFont(42);
   Graph_graphWSumc1->GetXaxis()->SetLabelSize(0.05);
   Graph_graphWSumc1->GetXaxis()->SetTitleSize(0.05);
   Graph_graphWSumc1->GetXaxis()->SetTitleOffset(1.4);
   Graph_graphWSumc1->GetXaxis()->SetTitleFont(42);
   Graph_graphWSumc1->GetYaxis()->SetLabelFont(42);
   Graph_graphWSumc1->GetYaxis()->SetLabelSize(0.05);
   Graph_graphWSumc1->GetYaxis()->SetTitleSize(0.05);
   Graph_graphWSumc1->GetYaxis()->SetTitleOffset(1.4);
   Graph_graphWSumc1->GetYaxis()->SetTitleFont(42);
   Graph_graphWSumc1->GetZaxis()->SetLabelFont(42);
   Graph_graphWSumc1->GetZaxis()->SetLabelSize(0.05);
   Graph_graphWSumc1->GetZaxis()->SetTitleSize(0.05);
   Graph_graphWSumc1->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_graphWSumc1);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("graphWSumSystc");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#999999");
   grae->SetFillColor(ci);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,395.16,27.70695);
   grae->SetPointError(0,13.9108,13.9108,0.9338258,0.9338258);
   grae->SetPoint(1,343.26,27.28461);
   grae->SetPointError(1,14.97234,14.97234,0.9748864,0.9748864);
   grae->SetPoint(2,294.88,29.56996);
   grae->SetPointError(2,15.66444,15.66444,1.164303,1.164303);
   grae->SetPoint(3,252.52,30.57403);
   grae->SetPointError(3,15.83232,15.83232,1.244725,1.244725);
   grae->SetPoint(4,170.83,30.54106);
   grae->SetPointError(4,16.10358,16.10358,0.7878689,0.7878689);
   grae->SetPoint(5,58.93,27.13788);
   grae->SetPointError(5,14.34243,14.34243,1.073423,1.073423);
   
   TH1F *Graph_graphWSumSystc2 = new TH1F("Graph_graphWSumSystc2","Graph",100,8.139247,445.5191);
   Graph_graphWSumSystc2->SetMinimum(25.48903);
   Graph_graphWSumSystc2->SetMaximum(32.39419);
   Graph_graphWSumSystc2->SetDirectory(0);
   Graph_graphWSumSystc2->SetStats(0);
   Graph_graphWSumSystc2->SetLineWidth(2);
   Graph_graphWSumSystc2->SetMarkerStyle(20);
   Graph_graphWSumSystc2->SetMarkerSize(1.2);
   Graph_graphWSumSystc2->GetXaxis()->SetLabelFont(42);
   Graph_graphWSumSystc2->GetXaxis()->SetLabelSize(0.05);
   Graph_graphWSumSystc2->GetXaxis()->SetTitleSize(0.05);
   Graph_graphWSumSystc2->GetXaxis()->SetTitleOffset(1.4);
   Graph_graphWSumSystc2->GetXaxis()->SetTitleFont(42);
   Graph_graphWSumSystc2->GetYaxis()->SetLabelFont(42);
   Graph_graphWSumSystc2->GetYaxis()->SetLabelSize(0.05);
   Graph_graphWSumSystc2->GetYaxis()->SetTitleSize(0.05);
   Graph_graphWSumSystc2->GetYaxis()->SetTitleOffset(1.4);
   Graph_graphWSumSystc2->GetYaxis()->SetTitleFont(42);
   Graph_graphWSumSystc2->GetZaxis()->SetLabelFont(42);
   Graph_graphWSumSystc2->GetZaxis()->SetLabelSize(0.05);
   Graph_graphWSumSystc2->GetZaxis()->SetTitleSize(0.05);
   Graph_graphWSumSystc2->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_graphWSumSystc2);
   
   grae->Draw("e2 ");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("graphWPlusSumSystc");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ff6666");
   grae->SetFillColor(ci);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,395.16,14.00369);
   grae->SetPointError(0,13.9108,13.9108,0.659828,0.659828);
   grae->SetPoint(1,343.26,13.20188);
   grae->SetPointError(1,14.97234,14.97234,0.6842623,0.6842623);
   grae->SetPoint(2,294.88,14.93654);
   grae->SetPointError(2,15.66444,15.66444,0.8378239,0.8378239);
   grae->SetPoint(3,252.52,16.3443);
   grae->SetPointError(3,15.83232,15.83232,0.9303988,0.9303988);
   grae->SetPoint(4,170.83,15.9511);
   grae->SetPointError(4,16.10358,16.10358,0.572779,0.572779);
   grae->SetPoint(5,58.93,14.06611);
   grae->SetPointError(5,14.34243,14.34243,0.7731905,0.7731905);
   
   TH1F *Graph_graphWPlusSumSystc3 = new TH1F("Graph_graphWPlusSumSystc3","Graph",100,8.139247,445.5191);
   Graph_graphWPlusSumSystc3->SetMinimum(12.04191);
   Graph_graphWPlusSumSystc3->SetMaximum(17.75041);
   Graph_graphWPlusSumSystc3->SetDirectory(0);
   Graph_graphWPlusSumSystc3->SetStats(0);
   Graph_graphWPlusSumSystc3->SetLineWidth(2);
   Graph_graphWPlusSumSystc3->SetMarkerStyle(20);
   Graph_graphWPlusSumSystc3->SetMarkerSize(1.2);
   Graph_graphWPlusSumSystc3->GetXaxis()->SetLabelFont(42);
   Graph_graphWPlusSumSystc3->GetXaxis()->SetLabelSize(0.05);
   Graph_graphWPlusSumSystc3->GetXaxis()->SetTitleSize(0.05);
   Graph_graphWPlusSumSystc3->GetXaxis()->SetTitleOffset(1.4);
   Graph_graphWPlusSumSystc3->GetXaxis()->SetTitleFont(42);
   Graph_graphWPlusSumSystc3->GetYaxis()->SetLabelFont(42);
   Graph_graphWPlusSumSystc3->GetYaxis()->SetLabelSize(0.05);
   Graph_graphWPlusSumSystc3->GetYaxis()->SetTitleSize(0.05);
   Graph_graphWPlusSumSystc3->GetYaxis()->SetTitleOffset(1.4);
   Graph_graphWPlusSumSystc3->GetYaxis()->SetTitleFont(42);
   Graph_graphWPlusSumSystc3->GetZaxis()->SetLabelFont(42);
   Graph_graphWPlusSumSystc3->GetZaxis()->SetLabelSize(0.05);
   Graph_graphWPlusSumSystc3->GetZaxis()->SetTitleSize(0.05);
   Graph_graphWPlusSumSystc3->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_graphWPlusSumSystc3);
   
   grae->Draw("e2 ");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("graphWMinusSumSystc");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#6666ff");
   grae->SetFillColor(ci);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,395.16,13.70327);
   grae->SetPointError(0,13.9108,13.9108,0.6608007,0.6608007);
   grae->SetPoint(1,343.26,14.08272);
   grae->SetPointError(1,14.97234,14.97234,0.6943979,0.6943979);
   grae->SetPoint(2,294.88,14.63342);
   grae->SetPointError(2,15.66444,15.66444,0.8084883,0.8084883);
   grae->SetPoint(3,252.52,14.22973);
   grae->SetPointError(3,15.83232,15.83232,0.82686,0.82686);
   grae->SetPoint(4,170.83,14.58996);
   grae->SetPointError(4,16.10358,16.10358,0.5409822,0.5409822);
   grae->SetPoint(5,58.93,13.07177);
   grae->SetPointError(5,14.34243,14.34243,0.7445897,0.7445897);
   
   TH1F *Graph_graphWMinusSumSystc4 = new TH1F("Graph_graphWMinusSumSystc4","Graph",100,8.139247,445.5191);
   Graph_graphWMinusSumSystc4->SetMinimum(12.01571);
   Graph_graphWMinusSumSystc4->SetMaximum(15.75339);
   Graph_graphWMinusSumSystc4->SetDirectory(0);
   Graph_graphWMinusSumSystc4->SetStats(0);
   Graph_graphWMinusSumSystc4->SetLineWidth(2);
   Graph_graphWMinusSumSystc4->SetMarkerStyle(20);
   Graph_graphWMinusSumSystc4->SetMarkerSize(1.2);
   Graph_graphWMinusSumSystc4->GetXaxis()->SetLabelFont(42);
   Graph_graphWMinusSumSystc4->GetXaxis()->SetLabelSize(0.05);
   Graph_graphWMinusSumSystc4->GetXaxis()->SetTitleSize(0.05);
   Graph_graphWMinusSumSystc4->GetXaxis()->SetTitleOffset(1.4);
   Graph_graphWMinusSumSystc4->GetXaxis()->SetTitleFont(42);
   Graph_graphWMinusSumSystc4->GetYaxis()->SetLabelFont(42);
   Graph_graphWMinusSumSystc4->GetYaxis()->SetLabelSize(0.05);
   Graph_graphWMinusSumSystc4->GetYaxis()->SetTitleSize(0.05);
   Graph_graphWMinusSumSystc4->GetYaxis()->SetTitleOffset(1.4);
   Graph_graphWMinusSumSystc4->GetYaxis()->SetTitleFont(42);
   Graph_graphWMinusSumSystc4->GetZaxis()->SetLabelFont(42);
   Graph_graphWMinusSumSystc4->GetZaxis()->SetLabelSize(0.05);
   Graph_graphWMinusSumSystc4->GetZaxis()->SetTitleSize(0.05);
   Graph_graphWMinusSumSystc4->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_graphWMinusSumSystc4);
   
   grae->Draw("e2 ");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("graphWSumSyst_correlatedc");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#999999");
   grae->SetFillColor(ci);
   grae->SetFillStyle(3244);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,382.16,27.70695);
   grae->SetPointError(0,6.9108,6.9108,2.202847,2.202847);
   grae->SetPoint(1,330.26,27.28461);
   grae->SetPointError(1,7.97234,7.97234,2.123296,2.123296);
   grae->SetPoint(2,281.88,29.56996);
   grae->SetPointError(2,8.66444,8.66444,2.287682,2.287682);
   grae->SetPoint(3,239.52,30.57403);
   grae->SetPointError(3,8.83232,8.83232,2.375929,2.375929);
   grae->SetPoint(4,157.83,30.54106);
   grae->SetPointError(4,9.10358,9.10358,2.296407,2.296407);
   grae->SetPoint(5,45.93,27.13788);
   grae->SetPointError(5,7.34243,7.34243,3.207776,3.207776);
   
   TH1F *Graph_graphWSumSyst_correlatedc5 = new TH1F("Graph_graphWSumSyst_correlatedc5","Graph",100,3.539247,424.1191);
   Graph_graphWSumSyst_correlatedc5->SetMinimum(23.02812);
   Graph_graphWSumSyst_correlatedc5->SetMaximum(33.85195);
   Graph_graphWSumSyst_correlatedc5->SetDirectory(0);
   Graph_graphWSumSyst_correlatedc5->SetStats(0);
   Graph_graphWSumSyst_correlatedc5->SetLineWidth(2);
   Graph_graphWSumSyst_correlatedc5->SetMarkerStyle(20);
   Graph_graphWSumSyst_correlatedc5->SetMarkerSize(1.2);
   Graph_graphWSumSyst_correlatedc5->GetXaxis()->SetLabelFont(42);
   Graph_graphWSumSyst_correlatedc5->GetXaxis()->SetLabelSize(0.05);
   Graph_graphWSumSyst_correlatedc5->GetXaxis()->SetTitleSize(0.05);
   Graph_graphWSumSyst_correlatedc5->GetXaxis()->SetTitleOffset(1.4);
   Graph_graphWSumSyst_correlatedc5->GetXaxis()->SetTitleFont(42);
   Graph_graphWSumSyst_correlatedc5->GetYaxis()->SetLabelFont(42);
   Graph_graphWSumSyst_correlatedc5->GetYaxis()->SetLabelSize(0.05);
   Graph_graphWSumSyst_correlatedc5->GetYaxis()->SetTitleSize(0.05);
   Graph_graphWSumSyst_correlatedc5->GetYaxis()->SetTitleOffset(1.4);
   Graph_graphWSumSyst_correlatedc5->GetYaxis()->SetTitleFont(42);
   Graph_graphWSumSyst_correlatedc5->GetZaxis()->SetLabelFont(42);
   Graph_graphWSumSyst_correlatedc5->GetZaxis()->SetLabelSize(0.05);
   Graph_graphWSumSyst_correlatedc5->GetZaxis()->SetTitleSize(0.05);
   Graph_graphWSumSyst_correlatedc5->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_graphWSumSyst_correlatedc5);
   
   grae->Draw("e2");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("graphWPlusSumSyst_correlatedc");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ff0000");
   grae->SetFillColor(ci);
   grae->SetFillStyle(3004);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,382.16,14.00369);
   grae->SetPointError(0,6.9108,6.9108,1.142452,1.142452);
   grae->SetPoint(1,330.26,13.20188);
   grae->SetPointError(1,7.97234,7.97234,1.071709,1.071709);
   grae->SetPoint(2,281.88,14.93654);
   grae->SetPointError(2,8.66444,8.66444,1.211872,1.211872);
   grae->SetPoint(3,239.52,16.3443);
   grae->SetPointError(3,8.83232,8.83232,1.323056,1.323056);
   grae->SetPoint(4,157.83,15.9511);
   grae->SetPointError(4,9.10358,9.10358,1.231042,1.231042);
   grae->SetPoint(5,45.93,14.06611);
   grae->SetPointError(5,7.34243,7.34243,1.682356,1.682356);
   
   TH1F *Graph_graphWPlusSumSyst_correlatedc6 = new TH1F("Graph_graphWPlusSumSyst_correlatedc6","Graph",100,3.539247,424.1191);
   Graph_graphWPlusSumSyst_correlatedc6->SetMinimum(11.57646);
   Graph_graphWPlusSumSyst_correlatedc6->SetMaximum(18.22107);
   Graph_graphWPlusSumSyst_correlatedc6->SetDirectory(0);
   Graph_graphWPlusSumSyst_correlatedc6->SetStats(0);
   Graph_graphWPlusSumSyst_correlatedc6->SetLineWidth(2);
   Graph_graphWPlusSumSyst_correlatedc6->SetMarkerStyle(20);
   Graph_graphWPlusSumSyst_correlatedc6->SetMarkerSize(1.2);
   Graph_graphWPlusSumSyst_correlatedc6->GetXaxis()->SetLabelFont(42);
   Graph_graphWPlusSumSyst_correlatedc6->GetXaxis()->SetLabelSize(0.05);
   Graph_graphWPlusSumSyst_correlatedc6->GetXaxis()->SetTitleSize(0.05);
   Graph_graphWPlusSumSyst_correlatedc6->GetXaxis()->SetTitleOffset(1.4);
   Graph_graphWPlusSumSyst_correlatedc6->GetXaxis()->SetTitleFont(42);
   Graph_graphWPlusSumSyst_correlatedc6->GetYaxis()->SetLabelFont(42);
   Graph_graphWPlusSumSyst_correlatedc6->GetYaxis()->SetLabelSize(0.05);
   Graph_graphWPlusSumSyst_correlatedc6->GetYaxis()->SetTitleSize(0.05);
   Graph_graphWPlusSumSyst_correlatedc6->GetYaxis()->SetTitleOffset(1.4);
   Graph_graphWPlusSumSyst_correlatedc6->GetYaxis()->SetTitleFont(42);
   Graph_graphWPlusSumSyst_correlatedc6->GetZaxis()->SetLabelFont(42);
   Graph_graphWPlusSumSyst_correlatedc6->GetZaxis()->SetLabelSize(0.05);
   Graph_graphWPlusSumSyst_correlatedc6->GetZaxis()->SetTitleSize(0.05);
   Graph_graphWPlusSumSyst_correlatedc6->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_graphWPlusSumSyst_correlatedc6);
   
   grae->Draw("e2 ");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("graphWMinusSumSyst_correlatedc");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#0000ff");
   grae->SetFillColor(ci);
   grae->SetFillStyle(3005);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,382.16,13.70327);
   grae->SetPointError(0,6.9108,6.9108,1.127756,1.127756);
   grae->SetPoint(1,330.26,14.08272);
   grae->SetPointError(1,7.97234,7.97234,1.126162,1.126162);
   grae->SetPoint(2,281.88,14.63342);
   grae->SetPointError(2,8.66444,8.66444,1.170918,1.170918);
   grae->SetPoint(3,239.52,14.22973);
   grae->SetPointError(3,8.83232,8.83232,1.160777,1.160777);
   grae->SetPoint(4,157.83,14.58996);
   grae->SetPointError(4,9.10358,9.10358,1.130321,1.130321);
   grae->SetPoint(5,45.93,13.07177);
   grae->SetPointError(5,7.34243,7.34243,1.584683,1.584683);
   
   TH1F *Graph_graphWMinusSumSyst_correlatedc7 = new TH1F("Graph_graphWMinusSumSyst_correlatedc7","Graph",100,3.539247,424.1191);
   Graph_graphWMinusSumSyst_correlatedc7->SetMinimum(11.05537);
   Graph_graphWMinusSumSyst_correlatedc7->SetMaximum(16.23607);
   Graph_graphWMinusSumSyst_correlatedc7->SetDirectory(0);
   Graph_graphWMinusSumSyst_correlatedc7->SetStats(0);
   Graph_graphWMinusSumSyst_correlatedc7->SetLineWidth(2);
   Graph_graphWMinusSumSyst_correlatedc7->SetMarkerStyle(20);
   Graph_graphWMinusSumSyst_correlatedc7->SetMarkerSize(1.2);
   Graph_graphWMinusSumSyst_correlatedc7->GetXaxis()->SetLabelFont(42);
   Graph_graphWMinusSumSyst_correlatedc7->GetXaxis()->SetLabelSize(0.05);
   Graph_graphWMinusSumSyst_correlatedc7->GetXaxis()->SetTitleSize(0.05);
   Graph_graphWMinusSumSyst_correlatedc7->GetXaxis()->SetTitleOffset(1.4);
   Graph_graphWMinusSumSyst_correlatedc7->GetXaxis()->SetTitleFont(42);
   Graph_graphWMinusSumSyst_correlatedc7->GetYaxis()->SetLabelFont(42);
   Graph_graphWMinusSumSyst_correlatedc7->GetYaxis()->SetLabelSize(0.05);
   Graph_graphWMinusSumSyst_correlatedc7->GetYaxis()->SetTitleSize(0.05);
   Graph_graphWMinusSumSyst_correlatedc7->GetYaxis()->SetTitleOffset(1.4);
   Graph_graphWMinusSumSyst_correlatedc7->GetYaxis()->SetTitleFont(42);
   Graph_graphWMinusSumSyst_correlatedc7->GetZaxis()->SetLabelFont(42);
   Graph_graphWMinusSumSyst_correlatedc7->GetZaxis()->SetLabelSize(0.05);
   Graph_graphWMinusSumSyst_correlatedc7->GetZaxis()->SetTitleSize(0.05);
   Graph_graphWMinusSumSyst_correlatedc7->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_graphWMinusSumSyst_correlatedc7);
   
   grae->Draw("e2 ");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("graphWSumDummy");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetMarkerColor(0);
   grae->SetMarkerStyle(33);
   grae->SetMarkerSize(2.2);
   grae->SetPoint(0,395.16,27.70695);
   grae->SetPointError(0,13.9108,13.9108,0.7727571,0.7727571);
   grae->SetPoint(1,343.26,27.28461);
   grae->SetPointError(1,14.97234,14.97234,0.8589446,0.8589446);
   grae->SetPoint(2,294.88,29.56996);
   grae->SetPointError(2,15.66444,15.66444,0.9792449,0.9792449);
   grae->SetPoint(3,252.52,30.57403);
   grae->SetPointError(3,15.83232,15.83232,1.12038,1.12038);
   grae->SetPoint(4,170.83,30.54106);
   grae->SetPointError(4,16.10358,16.10358,0.7212261,0.7212261);
   grae->SetPoint(5,58.93,27.13788);
   grae->SetPointError(5,14.34243,14.34243,1.038295,1.038295);
   
   TH1F *Graph_graphWSumDummy8 = new TH1F("Graph_graphWSumDummy8","Graph",100,8.139247,445.5191);
   Graph_graphWSumDummy8->SetMinimum(25.5401);
   Graph_graphWSumDummy8->SetMaximum(32.25389);
   Graph_graphWSumDummy8->SetDirectory(0);
   Graph_graphWSumDummy8->SetStats(0);
   Graph_graphWSumDummy8->SetLineWidth(2);
   Graph_graphWSumDummy8->SetMarkerStyle(20);
   Graph_graphWSumDummy8->SetMarkerSize(1.2);
   Graph_graphWSumDummy8->GetXaxis()->SetLabelFont(42);
   Graph_graphWSumDummy8->GetXaxis()->SetLabelSize(0.05);
   Graph_graphWSumDummy8->GetXaxis()->SetTitleSize(0.05);
   Graph_graphWSumDummy8->GetXaxis()->SetTitleOffset(1.4);
   Graph_graphWSumDummy8->GetXaxis()->SetTitleFont(42);
   Graph_graphWSumDummy8->GetYaxis()->SetLabelFont(42);
   Graph_graphWSumDummy8->GetYaxis()->SetLabelSize(0.05);
   Graph_graphWSumDummy8->GetYaxis()->SetTitleSize(0.05);
   Graph_graphWSumDummy8->GetYaxis()->SetTitleOffset(1.4);
   Graph_graphWSumDummy8->GetYaxis()->SetTitleFont(42);
   Graph_graphWSumDummy8->GetZaxis()->SetLabelFont(42);
   Graph_graphWSumDummy8->GetZaxis()->SetLabelSize(0.05);
   Graph_graphWSumDummy8->GetZaxis()->SetTitleSize(0.05);
   Graph_graphWSumDummy8->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_graphWSumDummy8);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("graphWPlusSumDummy");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetMarkerColor(0);
   grae->SetMarkerStyle(23);
   grae->SetMarkerSize(1.9);
   grae->SetPoint(0,395.16,14.00369);
   grae->SetPointError(0,13.9108,13.9108,0.5443481,0.5443481);
   grae->SetPoint(1,343.26,13.20188);
   grae->SetPointError(1,14.97234,14.97234,0.599948,0.599948);
   grae->SetPoint(2,294.88,14.93654);
   grae->SetPointError(2,15.66444,15.66444,0.6923321,0.6923321);
   grae->SetPoint(3,252.52,16.3443);
   grae->SetPointError(3,15.83232,15.83232,0.8229011,0.8229011);
   grae->SetPoint(4,170.83,15.9511);
   grae->SetPointError(4,16.10358,16.10358,0.5221116,0.5221116);
   grae->SetPoint(5,58.93,14.06611);
   grae->SetPointError(5,14.34243,14.34243,0.7438274,0.7438274);
   
   TH1F *Graph_graphWPlusSumDummy9 = new TH1F("Graph_graphWPlusSumDummy9","Graph",100,8.139247,445.5191);
   Graph_graphWPlusSumDummy9->SetMinimum(12.14541);
   Graph_graphWPlusSumDummy9->SetMaximum(17.62373);
   Graph_graphWPlusSumDummy9->SetDirectory(0);
   Graph_graphWPlusSumDummy9->SetStats(0);
   Graph_graphWPlusSumDummy9->SetLineWidth(2);
   Graph_graphWPlusSumDummy9->SetMarkerStyle(20);
   Graph_graphWPlusSumDummy9->SetMarkerSize(1.2);
   Graph_graphWPlusSumDummy9->GetXaxis()->SetLabelFont(42);
   Graph_graphWPlusSumDummy9->GetXaxis()->SetLabelSize(0.05);
   Graph_graphWPlusSumDummy9->GetXaxis()->SetTitleSize(0.05);
   Graph_graphWPlusSumDummy9->GetXaxis()->SetTitleOffset(1.4);
   Graph_graphWPlusSumDummy9->GetXaxis()->SetTitleFont(42);
   Graph_graphWPlusSumDummy9->GetYaxis()->SetLabelFont(42);
   Graph_graphWPlusSumDummy9->GetYaxis()->SetLabelSize(0.05);
   Graph_graphWPlusSumDummy9->GetYaxis()->SetTitleSize(0.05);
   Graph_graphWPlusSumDummy9->GetYaxis()->SetTitleOffset(1.4);
   Graph_graphWPlusSumDummy9->GetYaxis()->SetTitleFont(42);
   Graph_graphWPlusSumDummy9->GetZaxis()->SetLabelFont(42);
   Graph_graphWPlusSumDummy9->GetZaxis()->SetLabelSize(0.05);
   Graph_graphWPlusSumDummy9->GetZaxis()->SetTitleSize(0.05);
   Graph_graphWPlusSumDummy9->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_graphWPlusSumDummy9);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("graphWMinusSumDummy");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetMarkerColor(0);
   grae->SetMarkerStyle(22);
   grae->SetMarkerSize(1.9);
   grae->SetPoint(0,395.16,13.70327);
   grae->SetPointError(0,13.9108,13.9108,0.5484876,0.5484876);
   grae->SetPoint(1,343.26,14.08272);
   grae->SetPointError(1,14.97234,14.97234,0.6146936,0.6146936);
   grae->SetPoint(2,294.88,14.63342);
   grae->SetPointError(2,15.66444,15.66444,0.6925292,0.6925292);
   grae->SetPoint(3,252.52,14.22973);
   grae->SetPointError(3,15.83232,15.83232,0.7603189,0.7603189);
   grae->SetPoint(4,170.83,14.58996);
   grae->SetPointError(4,16.10358,16.10358,0.4975606,0.4975606);
   grae->SetPoint(5,58.93,13.07177);
   grae->SetPointError(5,14.34243,14.34243,0.7244151,0.7244151);
   
   TH1F *Graph_graphWMinusSumDummy10 = new TH1F("Graph_graphWMinusSumDummy10","Graph",100,8.139247,445.5191);
   Graph_graphWMinusSumDummy10->SetMinimum(12.0495);
   Graph_graphWMinusSumDummy10->SetMaximum(15.62381);
   Graph_graphWMinusSumDummy10->SetDirectory(0);
   Graph_graphWMinusSumDummy10->SetStats(0);
   Graph_graphWMinusSumDummy10->SetLineWidth(2);
   Graph_graphWMinusSumDummy10->SetMarkerStyle(20);
   Graph_graphWMinusSumDummy10->SetMarkerSize(1.2);
   Graph_graphWMinusSumDummy10->GetXaxis()->SetLabelFont(42);
   Graph_graphWMinusSumDummy10->GetXaxis()->SetLabelSize(0.05);
   Graph_graphWMinusSumDummy10->GetXaxis()->SetTitleSize(0.05);
   Graph_graphWMinusSumDummy10->GetXaxis()->SetTitleOffset(1.4);
   Graph_graphWMinusSumDummy10->GetXaxis()->SetTitleFont(42);
   Graph_graphWMinusSumDummy10->GetYaxis()->SetLabelFont(42);
   Graph_graphWMinusSumDummy10->GetYaxis()->SetLabelSize(0.05);
   Graph_graphWMinusSumDummy10->GetYaxis()->SetTitleSize(0.05);
   Graph_graphWMinusSumDummy10->GetYaxis()->SetTitleOffset(1.4);
   Graph_graphWMinusSumDummy10->GetYaxis()->SetTitleFont(42);
   Graph_graphWMinusSumDummy10->GetZaxis()->SetLabelFont(42);
   Graph_graphWMinusSumDummy10->GetZaxis()->SetLabelSize(0.05);
   Graph_graphWMinusSumDummy10->GetZaxis()->SetTitleSize(0.05);
   Graph_graphWMinusSumDummy10->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_graphWMinusSumDummy10);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("graphWSumc");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetMarkerStyle(33);
   grae->SetMarkerSize(1.8);
   grae->SetPoint(0,395.16,27.70695);
   grae->SetPointError(0,13.9108,13.9108,0.7727571,0.7727571);
   grae->SetPoint(1,343.26,27.28461);
   grae->SetPointError(1,14.97234,14.97234,0.8589446,0.8589446);
   grae->SetPoint(2,294.88,29.56996);
   grae->SetPointError(2,15.66444,15.66444,0.9792449,0.9792449);
   grae->SetPoint(3,252.52,30.57403);
   grae->SetPointError(3,15.83232,15.83232,1.12038,1.12038);
   grae->SetPoint(4,170.83,30.54106);
   grae->SetPointError(4,16.10358,16.10358,0.7212261,0.7212261);
   grae->SetPoint(5,58.93,27.13788);
   grae->SetPointError(5,14.34243,14.34243,1.038295,1.038295);
   
   TH1F *Graph_Graph_graphWSumc111 = new TH1F("Graph_Graph_graphWSumc111","Graph",100,8.139247,445.5191);
   Graph_Graph_graphWSumc111->SetMinimum(25.5401);
   Graph_Graph_graphWSumc111->SetMaximum(32.25389);
   Graph_Graph_graphWSumc111->SetDirectory(0);
   Graph_Graph_graphWSumc111->SetStats(0);
   Graph_Graph_graphWSumc111->SetLineWidth(2);
   Graph_Graph_graphWSumc111->SetMarkerStyle(20);
   Graph_Graph_graphWSumc111->SetMarkerSize(1.2);
   Graph_Graph_graphWSumc111->GetXaxis()->SetLabelFont(42);
   Graph_Graph_graphWSumc111->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_graphWSumc111->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph_graphWSumc111->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph_graphWSumc111->GetXaxis()->SetTitleFont(42);
   Graph_Graph_graphWSumc111->GetYaxis()->SetLabelFont(42);
   Graph_Graph_graphWSumc111->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_graphWSumc111->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph_graphWSumc111->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph_graphWSumc111->GetYaxis()->SetTitleFont(42);
   Graph_Graph_graphWSumc111->GetZaxis()->SetLabelFont(42);
   Graph_Graph_graphWSumc111->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_graphWSumc111->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph_graphWSumc111->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_graphWSumc111);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("graphWPlusSumc");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(23);
   grae->SetMarkerSize(1.8);
   grae->SetPoint(0,395.16,14.00369);
   grae->SetPointError(0,13.9108,13.9108,0.5443481,0.5443481);
   grae->SetPoint(1,343.26,13.20188);
   grae->SetPointError(1,14.97234,14.97234,0.599948,0.599948);
   grae->SetPoint(2,294.88,14.93654);
   grae->SetPointError(2,15.66444,15.66444,0.6923321,0.6923321);
   grae->SetPoint(3,252.52,16.3443);
   grae->SetPointError(3,15.83232,15.83232,0.8229011,0.8229011);
   grae->SetPoint(4,170.83,15.9511);
   grae->SetPointError(4,16.10358,16.10358,0.5221116,0.5221116);
   grae->SetPoint(5,58.93,14.06611);
   grae->SetPointError(5,14.34243,14.34243,0.7438274,0.7438274);
   
   TH1F *Graph_graphWPlusSumc12 = new TH1F("Graph_graphWPlusSumc12","Graph",100,8.139247,445.5191);
   Graph_graphWPlusSumc12->SetMinimum(12.14541);
   Graph_graphWPlusSumc12->SetMaximum(17.62373);
   Graph_graphWPlusSumc12->SetDirectory(0);
   Graph_graphWPlusSumc12->SetStats(0);
   Graph_graphWPlusSumc12->SetLineWidth(2);
   Graph_graphWPlusSumc12->SetMarkerStyle(20);
   Graph_graphWPlusSumc12->SetMarkerSize(1.2);
   Graph_graphWPlusSumc12->GetXaxis()->SetLabelFont(42);
   Graph_graphWPlusSumc12->GetXaxis()->SetLabelSize(0.05);
   Graph_graphWPlusSumc12->GetXaxis()->SetTitleSize(0.05);
   Graph_graphWPlusSumc12->GetXaxis()->SetTitleOffset(1.4);
   Graph_graphWPlusSumc12->GetXaxis()->SetTitleFont(42);
   Graph_graphWPlusSumc12->GetYaxis()->SetLabelFont(42);
   Graph_graphWPlusSumc12->GetYaxis()->SetLabelSize(0.05);
   Graph_graphWPlusSumc12->GetYaxis()->SetTitleSize(0.05);
   Graph_graphWPlusSumc12->GetYaxis()->SetTitleOffset(1.4);
   Graph_graphWPlusSumc12->GetYaxis()->SetTitleFont(42);
   Graph_graphWPlusSumc12->GetZaxis()->SetLabelFont(42);
   Graph_graphWPlusSumc12->GetZaxis()->SetLabelSize(0.05);
   Graph_graphWPlusSumc12->GetZaxis()->SetTitleSize(0.05);
   Graph_graphWPlusSumc12->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_graphWPlusSumc12);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("graphWMinusSumc");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(22);
   grae->SetMarkerSize(1.8);
   grae->SetPoint(0,395.16,13.70327);
   grae->SetPointError(0,13.9108,13.9108,0.5484876,0.5484876);
   grae->SetPoint(1,343.26,14.08272);
   grae->SetPointError(1,14.97234,14.97234,0.6146936,0.6146936);
   grae->SetPoint(2,294.88,14.63342);
   grae->SetPointError(2,15.66444,15.66444,0.6925292,0.6925292);
   grae->SetPoint(3,252.52,14.22973);
   grae->SetPointError(3,15.83232,15.83232,0.7603189,0.7603189);
   grae->SetPoint(4,170.83,14.58996);
   grae->SetPointError(4,16.10358,16.10358,0.4975606,0.4975606);
   grae->SetPoint(5,58.93,13.07177);
   grae->SetPointError(5,14.34243,14.34243,0.7244151,0.7244151);
   
   TH1F *Graph_graphWMinusSumc13 = new TH1F("Graph_graphWMinusSumc13","Graph",100,8.139247,445.5191);
   Graph_graphWMinusSumc13->SetMinimum(12.0495);
   Graph_graphWMinusSumc13->SetMaximum(15.62381);
   Graph_graphWMinusSumc13->SetDirectory(0);
   Graph_graphWMinusSumc13->SetStats(0);
   Graph_graphWMinusSumc13->SetLineWidth(2);
   Graph_graphWMinusSumc13->SetMarkerStyle(20);
   Graph_graphWMinusSumc13->SetMarkerSize(1.2);
   Graph_graphWMinusSumc13->GetXaxis()->SetLabelFont(42);
   Graph_graphWMinusSumc13->GetXaxis()->SetLabelSize(0.05);
   Graph_graphWMinusSumc13->GetXaxis()->SetTitleSize(0.05);
   Graph_graphWMinusSumc13->GetXaxis()->SetTitleOffset(1.4);
   Graph_graphWMinusSumc13->GetXaxis()->SetTitleFont(42);
   Graph_graphWMinusSumc13->GetYaxis()->SetLabelFont(42);
   Graph_graphWMinusSumc13->GetYaxis()->SetLabelSize(0.05);
   Graph_graphWMinusSumc13->GetYaxis()->SetTitleSize(0.05);
   Graph_graphWMinusSumc13->GetYaxis()->SetTitleOffset(1.4);
   Graph_graphWMinusSumc13->GetYaxis()->SetTitleFont(42);
   Graph_graphWMinusSumc13->GetZaxis()->SetLabelFont(42);
   Graph_graphWMinusSumc13->GetZaxis()->SetLabelSize(0.05);
   Graph_graphWMinusSumc13->GetZaxis()->SetTitleSize(0.05);
   Graph_graphWMinusSumc13->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_graphWMinusSumc13);
   
   grae->Draw("pe");
   
   TH1F *hDummy__2 = new TH1F("hDummy__2","hDummy",100,0,420);
   hDummy__2->SetMinimum(0);
   hDummy__2->SetMaximum(37);
   hDummy__2->SetDirectory(0);
   hDummy__2->SetStats(0);
   hDummy__2->SetLineWidth(2);
   hDummy__2->SetMarkerStyle(20);
   hDummy__2->SetMarkerSize(1.2);
   hDummy__2->GetXaxis()->SetTitle("#LT N_{part} #GT");
   hDummy__2->GetXaxis()->SetLabelFont(42);
   hDummy__2->GetXaxis()->SetLabelSize(0.036);
   hDummy__2->GetXaxis()->SetTitleSize(0.05);
   hDummy__2->GetXaxis()->SetTitleOffset(1.4);
   hDummy__2->GetXaxis()->SetTitleFont(42);
   hDummy__2->GetYaxis()->SetTitle("#frac{10^{9}}{#LT N_{coll} #GT}#frac{N^{W#rightarrow#mu,fiducial}}{N_{events}}");
   hDummy__2->GetYaxis()->SetLabelFont(42);
   hDummy__2->GetYaxis()->SetLabelSize(0.03);
   hDummy__2->GetYaxis()->SetTitleOffset(1.4);
   hDummy__2->GetYaxis()->SetTitleFont(42);
   hDummy__2->GetZaxis()->SetLabelFont(42);
   hDummy__2->GetZaxis()->SetLabelSize(0.05);
   hDummy__2->GetZaxis()->SetTitleSize(0.05);
   hDummy__2->GetZaxis()->SetTitleFont(42);
   hDummy__2->Draw("sameaxis");
   
   TLegend *leg = new TLegend(0.1694631,0.1730769,0.4848993,0.3741259,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03321678);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(0);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("NULL","W^{#pm}(Data)","p0f");
   entry->SetFillColor(15);
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("funcPythia","W^{#pm}(PYTHIA LO*)","l");

   ci = TColor::GetColor("#666666");
   entry->SetLineColor(ci);
   entry->SetLineStyle(2);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("funcMSTWnlo","W^{#pm}(POWHEG NLO)","l");

   ci = TColor::GetColor("#666666");
   entry->SetLineColor(ci);
   entry->SetLineStyle(9);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   
   leg = new TLegend(0.3909396,0.3094406,0.692953,0.3741259,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03321678);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(0);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("NULL","W^{+}(Data)","pf");

   ci = TColor::GetColor("#ff6666");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   leg->Draw();
   
   leg = new TLegend(0.6157718,0.3111888,0.9278523,0.3758741,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03321678);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(0);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("NULL","W^{-}(Data)","pf");

   ci = TColor::GetColor("#6666ff");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetMarkerColor(4);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   leg->Draw();
   TLatex *   tex = new TLatex(0.5872483,0.1975524,"#int Ldt #approx 0.14 nb^{-1}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03321678);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.545302,0.2692308,"Pb+Pb #sqrt{s_{NN}} = 2.76 TeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03321678);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(19.30338,33.77392,"ATLAS Preliminary");
   tex->SetTextFont(72);
   tex->SetLineWidth(2);
   tex->Draw();
   pad1->Modified();
   pad2->cd();
   pad2->Modified();
   cRcp->cd();
   cRcp->Modified();
   cRcp->cd();
   cRcp->SetSelected(cRcp);
}
