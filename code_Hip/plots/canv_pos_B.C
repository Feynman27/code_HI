{
//=========Macro generated from canvas: canv_pos_B/#mu^{+},PeriodB
//=========  (Tue Sep 30 17:57:20 2014) by ROOT version5.34/21
   TCanvas *canv_pos_B = new TCanvas("canv_pos_B", "#mu^{+},PeriodB",0,0,700,500);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   canv_pos_B->SetHighLightColor(2);
   canv_pos_B->Range(0,0,1,1);
   canv_pos_B->SetFillColor(0);
   canv_pos_B->SetBorderMode(0);
   canv_pos_B->SetBorderSize(2);
   canv_pos_B->SetTickx(1);
   canv_pos_B->SetTicky(1);
   canv_pos_B->SetLeftMargin(0.16);
   canv_pos_B->SetRightMargin(0.05);
   canv_pos_B->SetTopMargin(0.05);
   canv_pos_B->SetBottomMargin(0.16);
   canv_pos_B->SetFrameBorderMode(0);
   
   TH1F *h = new TH1F("h","h",100,-3,3);
   h->SetMinimum(0);
   h->SetMaximum(1.1);
   h->SetStats(0);
   h->SetLineWidth(2);
   h->SetMarkerStyle(20);
   h->SetMarkerSize(1.2);
   h->GetXaxis()->SetTitle("#eta^{#mu}_{lab}");
   h->GetXaxis()->SetLabelFont(42);
   h->GetXaxis()->SetLabelSize(0.05);
   h->GetXaxis()->SetTitleSize(0.05);
   h->GetXaxis()->SetTitleOffset(1.4);
   h->GetXaxis()->SetTitleFont(42);
   h->GetYaxis()->SetTitle("C_{W}");
   h->GetYaxis()->SetLabelFont(42);
   h->GetYaxis()->SetLabelSize(0.05);
   h->GetYaxis()->SetTitleSize(0.05);
   h->GetYaxis()->SetTitleOffset(1.4);
   h->GetYaxis()->SetTitleFont(42);
   h->GetZaxis()->SetLabelFont(42);
   h->GetZaxis()->SetLabelSize(0.05);
   h->GetZaxis()->SetTitleSize(0.05);
   h->GetZaxis()->SetTitleFont(42);
   h->Draw("");
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(11);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,-2.3,0.7064613);
   grae->SetPointError(0,0.2,0.2,0.01805477,0.01744622);
   grae->SetPoint(1,-1.825,0.7864518);
   grae->SetPointError(1,0.275,0.275,0.01084963,0.01047274);
   grae->SetPoint(2,-1.3,0.7606383);
   grae->SetPointError(2,0.25,0.25,0.01079392,0.01047795);
   grae->SetPoint(3,-0.825,0.6708523);
   grae->SetPointError(3,0.225,0.225,0.01216385,0.01194543);
   grae->SetPoint(4,-0.35,0.7384891);
   grae->SetPointError(4,0.25,0.25,0.01099688,0.01071428);
   grae->SetPoint(5,0,0.3232594);
   grae->SetPointError(5,0.1,0.1,0.01870362,0.01926791);
   grae->SetPoint(6,0.35,0.7630633);
   grae->SetPointError(6,0.25,0.25,0.01121566,0.01087028);
   grae->SetPoint(7,0.825,0.6797);
   grae->SetPointError(7,0.225,0.225,0.01300753,0.01273292);
   grae->SetPoint(8,1.3,0.8294746);
   grae->SetPointError(8,0.25,0.25,0.01018483,0.009731414);
   grae->SetPoint(9,1.825,0.8398886);
   grae->SetPointError(9,0.275,0.275,0.009782679,0.009334369);
   grae->SetPoint(10,2.3,0.7004052);
   grae->SetPointError(10,0.2,0.2,0.01518892,0.01473588);
   grae->Draw("pe");
   TLatex *   tex = new TLatex(0.3461371,0.8628191,"Period B");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.1933654,0.8628191,"W^{+}#rightarrow#mu^{+}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TLegend *leg = new TLegend(0.6623563,0.2457627,0.9324713,0.4872881,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph0","0-5%","pl");

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","5-10%","pl");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","10-20%","pl");

   ci = TColor::GetColor("#ff00ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff00ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph3","20-30%","pl");

   ci = TColor::GetColor("#00ff00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ff00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph4","30-40%","pl");

   ci = TColor::GetColor("#ffcc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph5","40-60%","pl");

   ci = TColor::GetColor("#6600ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#6600ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(29);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph6","60-90%","pl");

   ci = TColor::GetColor("#00ffff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ffff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(34);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   leg->Draw();
   
   grae = new TGraphAsymmErrors(11);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(21);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,-2.3,0.7050613);
   grae->SetPointError(0,0.2,0.2,0.01818718,0.01758343);
   grae->SetPoint(1,-1.825,0.8116697);
   grae->SetPointError(1,0.275,0.275,0.01088478,0.01043359);
   grae->SetPoint(2,-1.3,0.7742095);
   grae->SetPointError(2,0.25,0.25,0.01118233,0.01081088);
   grae->SetPoint(3,-0.825,0.6651895);
   grae->SetPointError(3,0.225,0.225,0.01328978,0.01303885);
   grae->SetPoint(4,-0.35,0.749622);
   grae->SetPointError(4,0.25,0.25,0.01145509,0.0111245);
   grae->SetPoint(5,0,0.3335516);
   grae->SetPointError(5,0.1,0.1,0.01941475,0.01997362);
   grae->SetPoint(6,0.35,0.768763);
   grae->SetPointError(6,0.25,0.25,0.01139127,0.01102107);
   grae->SetPoint(7,0.825,0.6489745);
   grae->SetPointError(7,0.225,0.225,0.01364338,0.01340752);
   grae->SetPoint(8,1.3,0.8176664);
   grae->SetPointError(8,0.25,0.25,0.01090159,0.01043039);
   grae->SetPoint(9,1.825,0.8632268);
   grae->SetPointError(9,0.275,0.275,0.009564289,0.00904075);
   grae->SetPoint(10,2.3,0.688085);
   grae->SetPointError(10,0.2,0.2,0.01543355,0.01503296);
   grae->Draw("pe");
   
   leg = new TLegend(0.6623563,0.2457627,0.9324713,0.4872881,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("Graph0","0-5%","pl");

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","5-10%","pl");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","10-20%","pl");

   ci = TColor::GetColor("#ff00ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff00ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph3","20-30%","pl");

   ci = TColor::GetColor("#00ff00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ff00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph4","30-40%","pl");

   ci = TColor::GetColor("#ffcc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph5","40-60%","pl");

   ci = TColor::GetColor("#6600ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#6600ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(29);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph6","60-90%","pl");

   ci = TColor::GetColor("#00ffff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ffff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(34);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   leg->Draw();
   
   grae = new TGraphAsymmErrors(11);
   grae->SetName("Graph2");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ff00ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#ff00ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(22);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,-2.3,0.7479138);
   grae->SetPointError(0,0.2,0.2,0.01262962,0.01223049);
   grae->SetPoint(1,-1.825,0.8027638);
   grae->SetPointError(1,0.275,0.275,0.00793343,0.007704539);
   grae->SetPoint(2,-1.3,0.7540489);
   grae->SetPointError(2,0.25,0.25,0.008229275,0.008051329);
   grae->SetPoint(3,-0.825,0.6586915);
   grae->SetPointError(3,0.225,0.225,0.009346525,0.009226818);
   grae->SetPoint(4,-0.35,0.7445834);
   grae->SetPointError(4,0.25,0.25,0.008188667,0.008023496);
   grae->SetPoint(5,0,0.3351461);
   grae->SetPointError(5,0.1,0.1,0.01403855,0.01432753);
   grae->SetPoint(6,0.35,0.7707397);
   grae->SetPointError(6,0.25,0.25,0.00826671,0.008065018);
   grae->SetPoint(7,0.825,0.6541911);
   grae->SetPointError(7,0.225,0.225,0.009929655,0.009800863);
   grae->SetPoint(8,1.3,0.8474151);
   grae->SetPointError(8,0.25,0.25,0.00726441,0.006997448);
   grae->SetPoint(9,1.825,0.8653461);
   grae->SetPointError(9,0.275,0.275,0.006766315,0.006493879);
   grae->SetPoint(10,2.3,0.6937235);
   grae->SetPointError(10,0.2,0.2,0.01141402,0.01118102);
   grae->Draw("pe");
   
   leg = new TLegend(0.6623563,0.2457627,0.9324713,0.4872881,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("Graph0","0-5%","pl");

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","5-10%","pl");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","10-20%","pl");

   ci = TColor::GetColor("#ff00ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff00ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph3","20-30%","pl");

   ci = TColor::GetColor("#00ff00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ff00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph4","30-40%","pl");

   ci = TColor::GetColor("#ffcc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph5","40-60%","pl");

   ci = TColor::GetColor("#6600ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#6600ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(29);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph6","60-90%","pl");

   ci = TColor::GetColor("#00ffff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ffff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(34);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   leg->Draw();
   
   grae = new TGraphAsymmErrors(11);
   grae->SetName("Graph3");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#00ff00");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#00ff00");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(23);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,-2.3,0.7528025);
   grae->SetPointError(0,0.2,0.2,0.01299851,0.01256581);
   grae->SetPoint(1,-1.825,0.8074182);
   grae->SetPointError(1,0.275,0.275,0.008119812,0.007872655);
   grae->SetPoint(2,-1.3,0.7613537);
   grae->SetPointError(2,0.25,0.25,0.008242576,0.008054891);
   grae->SetPoint(3,-0.825,0.6734255);
   grae->SetPointError(3,0.225,0.225,0.009467404,0.009331477);
   grae->SetPoint(4,-0.35,0.7697957);
   grae->SetPointError(4,0.25,0.25,0.008329518,0.00812773);
   grae->SetPoint(5,0,0.3453466);
   grae->SetPointError(5,0.1,0.1,0.01464658,0.01493532);
   grae->SetPoint(6,0.35,0.7872195);
   grae->SetPointError(6,0.25,0.25,0.008262593,0.008036077);
   grae->SetPoint(7,0.825,0.6447764);
   grae->SetPointError(7,0.225,0.225,0.01027517,0.0101453);
   grae->SetPoint(8,1.3,0.851812);
   grae->SetPointError(8,0.25,0.25,0.007539878,0.007236763);
   grae->SetPoint(9,1.825,0.8886343);
   grae->SetPointError(9,0.275,0.275,0.006583177,0.006263415);
   grae->SetPoint(10,2.3,0.6950256);
   grae->SetPointError(10,0.2,0.2,0.01150564,0.01127054);
   grae->Draw("pe");
   
   leg = new TLegend(0.6623563,0.2457627,0.9324713,0.4872881,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("Graph0","0-5%","pl");

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","5-10%","pl");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","10-20%","pl");

   ci = TColor::GetColor("#ff00ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff00ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph3","20-30%","pl");

   ci = TColor::GetColor("#00ff00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ff00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph4","30-40%","pl");

   ci = TColor::GetColor("#ffcc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph5","40-60%","pl");

   ci = TColor::GetColor("#6600ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#6600ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(29);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph6","60-90%","pl");

   ci = TColor::GetColor("#00ffff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ffff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(34);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   leg->Draw();
   
   grae = new TGraphAsymmErrors(11);
   grae->SetName("Graph4");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ffcc00");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(33);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,-2.3,0.7580458);
   grae->SetPointError(0,0.2,0.2,0.01296666,0.01252295);
   grae->SetPoint(1,-1.825,0.814377);
   grae->SetPointError(1,0.275,0.275,0.007978068,0.007727213);
   grae->SetPoint(2,-1.3,0.7635003);
   grae->SetPointError(2,0.25,0.25,0.008453097,0.008252984);
   grae->SetPoint(3,-0.825,0.6735117);
   grae->SetPointError(3,0.225,0.225,0.00955159,0.009411738);
   grae->SetPoint(4,-0.35,0.7594297);
   grae->SetPointError(4,0.25,0.25,0.008511,0.008313397);
   grae->SetPoint(5,0,0.3236631);
   grae->SetPointError(5,0.1,0.1,0.01443686,0.01477449);
   grae->SetPoint(6,0.35,0.7813594);
   grae->SetPointError(6,0.25,0.25,0.008377603,0.008154708);
   grae->SetPoint(7,0.825,0.6841816);
   grae->SetPointError(7,0.225,0.225,0.0101672,0.009997618);
   grae->SetPoint(8,1.3,0.8423043);
   grae->SetPointError(8,0.25,0.25,0.007652578,0.007365775);
   grae->SetPoint(9,1.825,0.880887);
   grae->SetPointError(9,0.275,0.275,0.006708532,0.006396946);
   grae->SetPoint(10,2.3,0.6908424);
   grae->SetPointError(10,0.2,0.2,0.01175464,0.011515);
   grae->Draw("pe");
   
   leg = new TLegend(0.6623563,0.2457627,0.9324713,0.4872881,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("Graph0","0-5%","pl");

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","5-10%","pl");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","10-20%","pl");

   ci = TColor::GetColor("#ff00ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff00ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph3","20-30%","pl");

   ci = TColor::GetColor("#00ff00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ff00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph4","30-40%","pl");

   ci = TColor::GetColor("#ffcc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph5","40-60%","pl");

   ci = TColor::GetColor("#6600ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#6600ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(29);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph6","60-90%","pl");

   ci = TColor::GetColor("#00ffff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ffff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(34);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   leg->Draw();
   
   grae = new TGraphAsymmErrors(11);
   grae->SetName("Graph5");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#6600ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#6600ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(29);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,-2.3,0.7582852);
   grae->SetPointError(0,0.2,0.2,0.008895831,0.008682361);
   grae->SetPoint(1,-1.825,0.8347787);
   grae->SetPointError(1,0.275,0.275,0.005394462,0.005259454);
   grae->SetPoint(2,-1.3,0.7755228);
   grae->SetPointError(2,0.25,0.25,0.005796281,0.005693161);
   grae->SetPoint(3,-0.825,0.6868464);
   grae->SetPointError(3,0.225,0.225,0.006731075,0.006653946);
   grae->SetPoint(4,-0.35,0.7627871);
   grae->SetPointError(4,0.25,0.25,0.005978529,0.005877508);
   grae->SetPoint(5,0,0.326969);
   grae->SetPointError(5,0.1,0.1,0.0104546,0.01062476);
   grae->SetPoint(6,0.35,0.7828283);
   grae->SetPointError(6,0.25,0.25,0.0059175,0.005804441);
   grae->SetPoint(7,0.825,0.669238);
   grae->SetPointError(7,0.225,0.225,0.007179414,0.007102363);
   grae->SetPoint(8,1.3,0.8487102);
   grae->SetPointError(8,0.25,0.25,0.005320057,0.005172569);
   grae->SetPoint(9,1.825,0.8891524);
   grae->SetPointError(9,0.275,0.275,0.004614904,0.004453872);
   grae->SetPoint(10,2.3,0.7136069);
   grae->SetPointError(10,0.2,0.2,0.008187349,0.008047507);
   grae->Draw("pe");
   
   leg = new TLegend(0.6623563,0.2457627,0.9324713,0.4872881,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("Graph0","0-5%","pl");

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","5-10%","pl");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","10-20%","pl");

   ci = TColor::GetColor("#ff00ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff00ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph3","20-30%","pl");

   ci = TColor::GetColor("#00ff00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ff00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph4","30-40%","pl");

   ci = TColor::GetColor("#ffcc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph5","40-60%","pl");

   ci = TColor::GetColor("#6600ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#6600ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(29);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph6","60-90%","pl");

   ci = TColor::GetColor("#00ffff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ffff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(34);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   leg->Draw();
   
   grae = new TGraphAsymmErrors(11);
   grae->SetName("Graph6");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#00ffff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#00ffff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(34);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,-2.3,0.787884);
   grae->SetPointError(0,0.2,0.2,0.008903937,0.00864506);
   grae->SetPoint(1,-1.825,0.8381262);
   grae->SetPointError(1,0.275,0.275,0.005729475,0.005573841);
   grae->SetPoint(2,-1.3,0.7863949);
   grae->SetPointError(2,0.25,0.25,0.006302256,0.00617175);
   grae->SetPoint(3,-0.825,0.6855781);
   grae->SetPointError(3,0.225,0.225,0.007421662,0.007329406);
   grae->SetPoint(4,-0.35,0.7736274);
   grae->SetPointError(4,0.25,0.25,0.006541371,0.006412258);
   grae->SetPoint(5,0,0.3262151);
   grae->SetPointError(5,0.1,0.1,0.01175154,0.01196942);
   grae->SetPoint(6,0.35,0.7763314);
   grae->SetPointError(6,0.25,0.25,0.006778686,0.006636504);
   grae->SetPoint(7,0.825,0.6633825);
   grae->SetPointError(7,0.225,0.225,0.008138111,0.008042738);
   grae->SetPoint(8,1.3,0.8607855);
   grae->SetPointError(8,0.25,0.25,0.005844168,0.005646542);
   grae->SetPoint(9,1.825,0.899964);
   grae->SetPointError(9,0.275,0.275,0.005007471,0.004795267);
   grae->SetPoint(10,2.3,0.7226713);
   grae->SetPointError(10,0.2,0.2,0.009085527,0.008901589);
   grae->Draw("pe");
   
   leg = new TLegend(0.6623563,0.2457627,0.9324713,0.4872881,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("Graph0","0-5%","pl");

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","5-10%","pl");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","10-20%","pl");

   ci = TColor::GetColor("#ff00ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff00ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph3","20-30%","pl");

   ci = TColor::GetColor("#00ff00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ff00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph4","30-40%","pl");

   ci = TColor::GetColor("#ffcc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph5","40-60%","pl");

   ci = TColor::GetColor("#6600ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#6600ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(29);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph6","60-90%","pl");

   ci = TColor::GetColor("#00ffff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ffff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(34);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   leg->Draw();
   canv_pos_B->Modified();
   canv_pos_B->cd();
   canv_pos_B->SetSelected(canv_pos_B);
}
