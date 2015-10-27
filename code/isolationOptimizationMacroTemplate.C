{
//=========Macro generated from canvas: c1_n2/c1_n2
//=========  (Mon Jun 24 18:45:17 2013) by ROOT version5.32/00
   TCanvas *c1_n2 = new TCanvas("c1_n2", "c1_n2",763,68,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n2->Range(-0.2227747,0.8144912,1.169567,1.002898);
   c1_n2->SetFillColor(0);
   c1_n2->SetBorderMode(0);
   c1_n2->SetBorderSize(2);
   c1_n2->SetTickx(1);
   c1_n2->SetTicky(1);
   c1_n2->SetLeftMargin(0.16);
   c1_n2->SetRightMargin(0.05);
   c1_n2->SetTopMargin(0.05);
   c1_n2->SetBottomMargin(0.16);
   c1_n2->SetFrameBorderMode(0);
   c1_n2->SetFrameBorderMode(0);
   
   TGraph *graph = new TGraph(11);
   graph->SetName("Graph");
   graph->SetTitle("Graph");
   graph->SetLineWidth(2);
   graph->SetMarkerColor(2);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.4);
   graph->SetPoint(0,0.05,0.9489273935);
   graph->SetPoint(1,0.1,0.9514858904);
   graph->SetPoint(2,0.2,0.9394405983);
   graph->SetPoint(3,0.3,0.9263055118);
   graph->SetPoint(4,0.4,0.9133561179);
   graph->SetPoint(5,0.5,0.9049750572);
   graph->SetPoint(6,0.7,0.8929591784);
   graph->SetPoint(7,1,0.8887938907);
   graph->SetPoint(8,5,0.8854540364);
   graph->SetPoint(9,10,0.8854540364);
   graph->SetPoint(10,100,0.8854540364);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,0,109.995);
   Graph_Graph1->SetMinimum(0.8446363);
   Graph_Graph1->SetMaximum(0.9934778);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   Graph_Graph1->SetLineWidth(2);
   Graph_Graph1->SetMarkerStyle(20);
   Graph_Graph1->SetMarkerSize(1.2);
   Graph_Graph1->GetXaxis()->SetTitle("i_{#mu} = #frac{#Sigma p_{T}^{trk}(#Delta R, p_{T}^{trk})}{p_{T}^{#mu}}");
   Graph_Graph1->GetXaxis()->SetRange(1,1);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.03);
   Graph_Graph1->GetXaxis()->SetTitleOffset(1.89);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetTitle("#epsilon_{eff}");
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph1);
   
   graph->Draw("ape");
   c1_n2->Modified();
   c1_n2->cd();
   c1_n2->SetSelected(c1_n2);
}
