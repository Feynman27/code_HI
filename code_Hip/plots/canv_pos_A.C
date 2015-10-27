{
//=========Macro generated from canvas: canv_pos_A/#mu^{+},PeriodA
//=========  (Tue Sep 30 17:57:19 2014) by ROOT version5.34/21
   TCanvas *canv_pos_A = new TCanvas("canv_pos_A", "#mu^{+},PeriodA",0,0,700,500);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   canv_pos_A->SetHighLightColor(2);
   canv_pos_A->Range(0,0,1,1);
   canv_pos_A->SetFillColor(0);
   canv_pos_A->SetBorderMode(0);
   canv_pos_A->SetBorderSize(2);
   canv_pos_A->SetTickx(1);
   canv_pos_A->SetTicky(1);
   canv_pos_A->SetLeftMargin(0.16);
   canv_pos_A->SetRightMargin(0.05);
   canv_pos_A->SetTopMargin(0.05);
   canv_pos_A->SetBottomMargin(0.16);
   canv_pos_A->SetFrameBorderMode(0);
   
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
   grae->SetPoint(0,-2.3,0.7338724);
   grae->SetPointError(0,0.2,0.2,0.02078201,0.01977221);
   grae->SetPoint(1,-1.825,0.78312);
   grae->SetPointError(1,0.275,0.275,0.01563978,0.01487834);
   grae->SetPoint(2,-1.3,0.7670104);
   grae->SetPointError(2,0.25,0.25,0.01633438,0.01559951);
   grae->SetPoint(3,-0.825,0.6572877);
   grae->SetPointError(3,0.225,0.225,0.01892484,0.01843746);
   grae->SetPoint(4,-0.35,0.75129);
   grae->SetPointError(4,0.25,0.25,0.01666759,0.0159904);
   grae->SetPoint(5,0,0.35411);
   grae->SetPointError(5,0.1,0.1,0.02742733,0.02838026);
   grae->SetPoint(6,0.35,0.7430248);
   grae->SetPointError(6,0.25,0.25,0.01612956,0.0155105);
   grae->SetPoint(7,0.825,0.6360369);
   grae->SetPointError(7,0.225,0.225,0.01876311,0.01835393);
   grae->SetPoint(8,1.3,0.800105);
   grae->SetPointError(8,0.25,0.25,0.01496451,0.0141854);
   grae->SetPoint(9,1.825,0.8400846);
   grae->SetPointError(9,0.275,0.275,0.01456845,0.01360737);
   grae->SetPoint(10,2.3,0.7050557);
   grae->SetPointError(10,0.2,0.2,0.02619523,0.02497804);
   grae->Draw("pe");
   TLatex *   tex = new TLatex(0.3461371,0.8628191,"Period A");
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
   entry=leg->AddEntry("Graph","5-10%","pl");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","10-20%","pl");

   ci = TColor::GetColor("#ff00ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff00ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","20-30%","pl");

   ci = TColor::GetColor("#00ff00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ff00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","30-40%","pl");

   ci = TColor::GetColor("#ffcc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","40-60%","pl");

   ci = TColor::GetColor("#6600ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#6600ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(29);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","60-90%","pl");

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
   grae->SetPoint(0,-2.3,0.7445737);
   grae->SetPointError(0,0.2,0.2,0.02255429,0.02124538);
   grae->SetPoint(1,-1.825,0.8234985);
   grae->SetPointError(1,0.275,0.275,0.01507288,0.01412165);
   grae->SetPoint(2,-1.3,0.7976583);
   grae->SetPointError(2,0.25,0.25,0.01685077,0.01589522);
   grae->SetPoint(3,-0.825,0.6783343);
   grae->SetPointError(3,0.225,0.225,0.01960893,0.01902141);
   grae->SetPoint(4,-0.35,0.770616);
   grae->SetPointError(4,0.25,0.25,0.01685941,0.01605383);
   grae->SetPoint(5,0,0.3240906);
   grae->SetPointError(5,0.1,0.1,0.02869158,0.03015045);
   grae->SetPoint(6,0.35,0.7660625);
   grae->SetPointError(6,0.25,0.25,0.01641786,0.01567637);
   grae->SetPoint(7,0.825,0.6477441);
   grae->SetPointError(7,0.225,0.225,0.0192404,0.01878534);
   grae->SetPoint(8,1.3,0.8200669);
   grae->SetPointError(8,0.25,0.25,0.0151754,0.01426676);
   grae->SetPoint(9,1.825,0.8455096);
   grae->SetPointError(9,0.275,0.275,0.01494669,0.01388478);
   grae->SetPoint(10,2.3,0.6805313);
   grae->SetPointError(10,0.2,0.2,0.02838293,0.0271678);
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
   entry=leg->AddEntry("Graph","10-20%","pl");

   ci = TColor::GetColor("#ff00ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff00ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","20-30%","pl");

   ci = TColor::GetColor("#00ff00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ff00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","30-40%","pl");

   ci = TColor::GetColor("#ffcc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","40-60%","pl");

   ci = TColor::GetColor("#6600ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#6600ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(29);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","60-90%","pl");

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
   grae->SetPoint(0,-2.3,0.7417097);
   grae->SetPointError(0,0.2,0.2,0.01570483,0.01509933);
   grae->SetPoint(1,-1.825,0.8180825);
   grae->SetPointError(1,0.275,0.275,0.01118174,0.0106805);
   grae->SetPoint(2,-1.3,0.7486617);
   grae->SetPointError(2,0.25,0.25,0.01289681,0.01248622);
   grae->SetPoint(3,-0.825,0.652039);
   grae->SetPointError(3,0.225,0.225,0.01464529,0.01437469);
   grae->SetPoint(4,-0.35,0.7644852);
   grae->SetPointError(4,0.25,0.25,0.01211255,0.01170585);
   grae->SetPoint(5,0,0.3477831);
   grae->SetPointError(5,0.1,0.1,0.02029024,0.0208351);
   grae->SetPoint(6,0.35,0.7460678);
   grae->SetPointError(6,0.25,0.25,0.01206409,0.01171121);
   grae->SetPoint(7,0.825,0.6437533);
   grae->SetPointError(7,0.225,0.225,0.01371493,0.0134874);
   grae->SetPoint(8,1.3,0.837157);
   grae->SetPointError(8,0.25,0.25,0.01041476,0.009914259);
   grae->SetPoint(9,1.825,0.8626405);
   grae->SetPointError(9,0.275,0.275,0.01019496,0.009609458);
   grae->SetPoint(10,2.3,0.7212859);
   grae->SetPointError(10,0.2,0.2,0.01884575,0.01811884);
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
   entry=leg->AddEntry("Graph","20-30%","pl");

   ci = TColor::GetColor("#00ff00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#00ff00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","30-40%","pl");

   ci = TColor::GetColor("#ffcc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","40-60%","pl");

   ci = TColor::GetColor("#6600ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#6600ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(29);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","60-90%","pl");

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
   grae->SetPoint(0,-2.3,0.7452234);
   grae->SetPointError(0,0.2,0.2,0.01606023,0.01540239);
   grae->SetPoint(1,-1.825,0.8101752);
   grae->SetPointError(1,0.275,0.275,0.01181252,0.01126993);
   grae->SetPoint(2,-1.3,0.7914621);
   grae->SetPointError(2,0.25,0.25,0.01215584,0.01165571);
   grae->SetPoint(3,-0.825,0.6645838);
   grae->SetPointError(3,0.225,0.225,0.01486415,0.01454133);
   grae->SetPoint(4,-0.35,0.7664673);
   grae->SetPointError(4,0.25,0.25,0.01221267,0.01178454);
   grae->SetPoint(5,0,0.333553);
   grae->SetPointError(5,0.1,0.1,0.02089993,0.02156611);
   grae->SetPoint(6,0.35,0.7673523);
   grae->SetPointError(6,0.25,0.25,0.01204078,0.01162973);
   grae->SetPoint(7,0.825,0.6538261);
   grae->SetPointError(7,0.225,0.225,0.01396294,0.01370473);
   grae->SetPoint(8,1.3,0.8403975);
   grae->SetPointError(8,0.25,0.25,0.01041126,0.00989859);
   grae->SetPoint(9,1.825,0.8661245);
   grae->SetPointError(9,0.275,0.275,0.0102447,0.00963543);
   grae->SetPoint(10,2.3,0.7336143);
   grae->SetPointError(10,0.2,0.2,0.0193,0.01848091);
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
   entry=leg->AddEntry("Graph","30-40%","pl");

   ci = TColor::GetColor("#ffcc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","40-60%","pl");

   ci = TColor::GetColor("#6600ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#6600ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(29);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","60-90%","pl");

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
   grae->SetPoint(0,-2.3,0.7382606);
   grae->SetPointError(0,0.2,0.2,0.01661321,0.01594594);
   grae->SetPoint(1,-1.825,0.8374589);
   grae->SetPointError(1,0.275,0.275,0.01140739,0.01078026);
   grae->SetPoint(2,-1.3,0.7754246);
   grae->SetPointError(2,0.25,0.25,0.01296533,0.01245468);
   grae->SetPoint(3,-0.825,0.6697277);
   grae->SetPointError(3,0.225,0.225,0.01496235,0.01461723);
   grae->SetPoint(4,-0.35,0.7762572);
   grae->SetPointError(4,0.25,0.25,0.01226435,0.01181333);
   grae->SetPoint(5,0,0.3588464);
   grae->SetPointError(5,0.1,0.1,0.02151722,0.02207704);
   grae->SetPoint(6,0.35,0.7894161);
   grae->SetPointError(6,0.25,0.25,0.01165197,0.01120545);
   grae->SetPoint(7,0.825,0.6734752);
   grae->SetPointError(7,0.225,0.225,0.01401291,0.01370748);
   grae->SetPoint(8,1.3,0.8435046);
   grae->SetPointError(8,0.25,0.25,0.01077415,0.0102114);
   grae->SetPoint(9,1.825,0.8728562);
   grae->SetPointError(9,0.275,0.275,0.01002635,0.009401074);
   grae->SetPoint(10,2.3,0.7340224);
   grae->SetPointError(10,0.2,0.2,0.01932027,0.01849707);
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
   entry=leg->AddEntry("Graph","40-60%","pl");

   ci = TColor::GetColor("#6600ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#6600ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(29);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","60-90%","pl");

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
   grae->SetPoint(0,-2.3,0.7512715);
   grae->SetPointError(0,0.2,0.2,0.01147539,0.01112571);
   grae->SetPoint(1,-1.825,0.832926);
   grae->SetPointError(1,0.275,0.275,0.008083231,0.007782179);
   grae->SetPoint(2,-1.3,0.7902497);
   grae->SetPointError(2,0.25,0.25,0.008912077,0.008637137);
   grae->SetPoint(3,-0.825,0.6726595);
   grae->SetPointError(3,0.225,0.225,0.01040134,0.01023507);
   grae->SetPoint(4,-0.35,0.7632313);
   grae->SetPointError(4,0.25,0.25,0.008942629,0.008719074);
   grae->SetPoint(5,0,0.3594284);
   grae->SetPointError(5,0.1,0.1,0.01543786,0.0157268);
   grae->SetPoint(6,0.35,0.7830871);
   grae->SetPointError(6,0.25,0.25,0.008469904,0.008240237);
   grae->SetPoint(7,0.825,0.6644495);
   grae->SetPointError(7,0.225,0.225,0.009881744,0.009742245);
   grae->SetPoint(8,1.3,0.8356544);
   grae->SetPointError(8,0.25,0.25,0.007460811,0.007198762);
   grae->SetPoint(9,1.825,0.88395);
   grae->SetPointError(9,0.275,0.275,0.006907987,0.006574958);
   grae->SetPoint(10,2.3,0.7255046);
   grae->SetPointError(10,0.2,0.2,0.01345123,0.0130677);
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
   entry=leg->AddEntry("Graph","60-90%","pl");

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
   grae->SetPoint(0,-2.3,0.7680076);
   grae->SetPointError(0,0.2,0.2,0.01281596,0.01233384);
   grae->SetPoint(1,-1.825,0.8475736);
   grae->SetPointError(1,0.275,0.275,0.008843597,0.00844201);
   grae->SetPoint(2,-1.3,0.7898425);
   grae->SetPointError(2,0.25,0.25,0.01011995,0.009777782);
   grae->SetPoint(3,-0.825,0.6855963);
   grae->SetPointError(3,0.225,0.225,0.01173476,0.01150427);
   grae->SetPoint(4,-0.35,0.7863845);
   grae->SetPointError(4,0.25,0.25,0.009529602,0.009236515);
   grae->SetPoint(5,0,0.3305519);
   grae->SetPointError(5,0.1,0.1,0.01665587,0.01707972);
   grae->SetPoint(6,0.35,0.7760344);
   grae->SetPointError(6,0.25,0.25,0.009534391,0.009255037);
   grae->SetPoint(7,0.825,0.6728338);
   grae->SetPointError(7,0.225,0.225,0.01110998,0.01092435);
   grae->SetPoint(8,1.3,0.8493777);
   grae->SetPointError(8,0.25,0.25,0.007975429,0.007647647);
   grae->SetPoint(9,1.825,0.9013661);
   grae->SetPointError(9,0.275,0.275,0.006737378,0.006361568);
   grae->SetPoint(10,2.3,0.7422409);
   grae->SetPointError(10,0.2,0.2,0.01397177,0.0135019);
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
   canv_pos_A->Modified();
   canv_pos_A->cd();
   canv_pos_A->SetSelected(canv_pos_A);
}
