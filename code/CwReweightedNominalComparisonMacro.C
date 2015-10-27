{
//=========Macro generated from canvas: c1/c1
//=========  (Thu Dec 12 18:36:33 2013) by ROOT version5.34/13
   TCanvas *c1 = new TCanvas("c1", "c1",761,98,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-0.5326582,-0.2025316,2.796456,1.063291);
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
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(9);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");
   grae->SetLineColor(4);
   grae->SetLineWidth(2);
   grae->SetMarkerColor(4);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,0.225,0.5011887);
   grae->SetPointError(0,0.125,0.125,0.02007899,0.01997956);
   grae->SetPoint(1,0.475,0.5082583);
   grae->SetPointError(1,0.125,0.125,0.0198909,0.01980094);
   grae->SetPoint(2,0.7,0.5478534);
   grae->SetPointError(2,0.1,0.1,0.02281384,0.02255052);
   grae->SetPoint(3,0.925,0.5318238);
   grae->SetPointError(3,0.125,0.125,0.02105879,0.02089744);
   grae->SetPoint(4,1.21,0.4359913);
   grae->SetPointError(4,0.16,0.16,0.01945566,0.01951054);
   grae->SetPoint(5,1.445,0.5486026);
   grae->SetPointError(5,0.075,0.075,0.02834075,0.02781515);
   grae->SetPoint(6,1.63,0.5486248);
   grae->SetPointError(6,0.11,0.11,0.02493297,0.02456364);
   grae->SetPoint(7,1.92,0.4863851);
   grae->SetPointError(7,0.18,0.18,0.0225764,0.02246937);
   grae->SetPoint(8,2.25,0.2928353);
   grae->SetPointError(8,0.15,0.15,0.02393082,0.02451655);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,0,2.63);
   Graph_Graph1->SetMinimum(0);
   Graph_Graph1->SetMaximum(1);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   Graph_Graph1->SetLineWidth(2);
   Graph_Graph1->SetMarkerStyle(20);
   Graph_Graph1->SetMarkerSize(1.2);
   Graph_Graph1->GetXaxis()->SetTitle("|#eta|");
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetTitle("C_{W}");
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

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#6666cc");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#6666cc");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,0.225,0.5450349);
   grae->SetPointError(0,0.125,0.125,0.02081169,0.02059277);
   grae->SetPoint(1,0.475,0.5516796);
   grae->SetPointError(1,0.125,0.125,0.02059607,0.02039519);
   grae->SetPoint(2,0.7,0.5798752);
   grae->SetPointError(2,0.1,0.1,0.02322634,0.0228457);
   grae->SetPoint(3,0.925,0.5424375);
   grae->SetPointError(3,0.125,0.125,0.02122409,0.02103273);
   grae->SetPoint(4,1.21,0.4224368);
   grae->SetPointError(4,0.16,0.16,0.01908285,0.01916566);
   grae->SetPoint(5,1.445,0.5118938);
   grae->SetPointError(5,0.075,0.075,0.02764536,0.02733881);
   grae->SetPoint(6,1.63,0.5019744);
   grae->SetPointError(6,0.11,0.11,0.02410683,0.02393855);
   grae->SetPoint(7,1.92,0.4445123);
   grae->SetPointError(7,0.18,0.18,0.02154353,0.0215691);
   grae->SetPoint(8,2.25,0.2849845);
   grae->SetPointError(8,0.15,0.15,0.02341076,0.02401026);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,0,2.63);
   Graph_Graph2->SetMinimum(0.227459);
   Graph_Graph2->SetMaximum(0.6368356);
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
   TLatex *   tex = new TLatex(0.2600575,0.7415254,"#mu^{+}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.1313559);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TLegend *leg = new TLegend(0.659,0.713,0.9295,0.8636,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("NULL","Nominal","pe");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("NULL","Reweighted","pe");

   ci = TColor::GetColor("#6666cc");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#6666cc");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
