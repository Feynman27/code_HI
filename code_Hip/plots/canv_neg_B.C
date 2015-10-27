{
//=========Macro generated from canvas: canv_neg_B/#mu^{-},PeriodB
//=========  (Tue Sep 30 17:57:18 2014) by ROOT version5.34/21
   TCanvas *canv_neg_B = new TCanvas("canv_neg_B", "#mu^{-},PeriodB",0,0,700,500);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   canv_neg_B->SetHighLightColor(2);
   canv_neg_B->Range(0,0,1,1);
   canv_neg_B->SetFillColor(0);
   canv_neg_B->SetBorderMode(0);
   canv_neg_B->SetBorderSize(2);
   canv_neg_B->SetTickx(1);
   canv_neg_B->SetTicky(1);
   canv_neg_B->SetLeftMargin(0.16);
   canv_neg_B->SetRightMargin(0.05);
   canv_neg_B->SetTopMargin(0.05);
   canv_neg_B->SetBottomMargin(0.16);
   canv_neg_B->SetFrameBorderMode(0);
   
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
   grae->SetPoint(0,-2.3,0.5988004);
   grae->SetPointError(0,0.2,0.2,0.0207611,0.02037103);
   grae->SetPoint(1,-1.825,0.8099484);
   grae->SetPointError(1,0.275,0.275,0.01253683,0.01192557);
   grae->SetPoint(2,-1.3,0.7595955);
   grae->SetPointError(2,0.25,0.25,0.01314718,0.01267216);
   grae->SetPoint(3,-0.825,0.6530935);
   grae->SetPointError(3,0.225,0.225,0.01456966,0.01430226);
   grae->SetPoint(4,-0.35,0.7683859);
   grae->SetPointError(4,0.25,0.25,0.01175876,0.01134982);
   grae->SetPoint(5,0,0.3178856);
   grae->SetPointError(5,0.1,0.1,0.02012644,0.02083628);
   grae->SetPoint(6,0.35,0.7637087);
   grae->SetPointError(6,0.25,0.25,0.01152614,0.01115085);
   grae->SetPoint(7,0.825,0.643191);
   grae->SetPointError(7,0.225,0.225,0.01355005,0.01332744);
   grae->SetPoint(8,1.3,0.7334059);
   grae->SetPointError(8,0.25,0.25,0.01181411,0.01149252);
   grae->SetPoint(9,1.825,0.7873247);
   grae->SetPointError(9,0.275,0.275,0.01012915,0.009787611);
   grae->SetPoint(10,2.3,0.6548644);
   grae->SetPointError(10,0.2,0.2,0.01400741,0.01374912);
   grae->Draw("pe");
   TLatex *   tex = new TLatex(0.3461371,0.8628191,"Period B");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.1933654,0.8628191,"W^{-}#rightarrow#mu^{-}");
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
   TLegendEntry *entry=leg->AddEntry("Graph","0-5%","pl");

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
   grae->SetPoint(0,-2.3,0.5435643);
   grae->SetPointError(0,0.2,0.2,0.02103036,0.02091175);
   grae->SetPoint(1,-1.825,0.805379);
   grae->SetPointError(1,0.275,0.275,0.01341838,0.01276626);
   grae->SetPoint(2,-1.3,0.7863929);
   grae->SetPointError(2,0.25,0.25,0.01308033,0.01252998);
   grae->SetPoint(3,-0.825,0.6403663);
   grae->SetPointError(3,0.225,0.225,0.01531912,0.01504354);
   grae->SetPoint(4,-0.35,0.766607);
   grae->SetPointError(4,0.25,0.25,0.01236557,0.01193452);
   grae->SetPoint(5,0,0.2858414);
   grae->SetPointError(5,0.1,0.1,0.01953483,0.02039252);
   grae->SetPoint(6,0.35,0.7653341);
   grae->SetPointError(6,0.25,0.25,0.01169257,0.01129448);
   grae->SetPoint(7,0.825,0.660874);
   grae->SetPointError(7,0.225,0.225,0.01395742,0.01368887);
   grae->SetPoint(8,1.3,0.7508987);
   grae->SetPointError(8,0.25,0.25,0.01210908,0.01173145);
   grae->SetPoint(9,1.825,0.7980105);
   grae->SetPointError(9,0.275,0.275,0.01058234,0.01019268);
   grae->SetPoint(10,2.3,0.6550094);
   grae->SetPointError(10,0.2,0.2,0.01486449,0.0145737);
   grae->Draw("pe");
   
   leg = new TLegend(0.6623563,0.2457627,0.9324713,0.4872881,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("Graph","0-5%","pl");

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
   grae->SetName("Graph2");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ff00ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#ff00ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(22);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,-2.3,0.5813392);
   grae->SetPointError(0,0.2,0.2,0.01529887,0.01514717);
   grae->SetPoint(1,-1.825,0.8036244);
   grae->SetPointError(1,0.275,0.275,0.00934129,0.009018516);
   grae->SetPoint(2,-1.3,0.7819917);
   grae->SetPointError(2,0.25,0.25,0.009435232,0.009149747);
   grae->SetPoint(3,-0.825,0.6597975);
   grae->SetPointError(3,0.225,0.225,0.01092548,0.01075941);
   grae->SetPoint(4,-0.35,0.7711427);
   grae->SetPointError(4,0.25,0.25,0.008938552,0.00870203);
   grae->SetPoint(5,0,0.343092);
   grae->SetPointError(5,0.1,0.1,0.01544742,0.0157745);
   grae->SetPoint(6,0.35,0.7650394);
   grae->SetPointError(6,0.25,0.25,0.008670894,0.008455175);
   grae->SetPoint(7,0.825,0.6633985);
   grae->SetPointError(7,0.225,0.225,0.0099954,0.009848672);
   grae->SetPoint(8,1.3,0.7433723);
   grae->SetPointError(8,0.25,0.25,0.008707008,0.00851786);
   grae->SetPoint(9,1.825,0.7930476);
   grae->SetPointError(9,0.275,0.275,0.007791988,0.007583052);
   grae->SetPoint(10,2.3,0.6713242);
   grae->SetPointError(10,0.2,0.2,0.01054844,0.01038259);
   grae->Draw("pe");
   
   leg = new TLegend(0.6623563,0.2457627,0.9324713,0.4872881,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("Graph","0-5%","pl");

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
   grae->SetName("Graph3");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#00ff00");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#00ff00");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(23);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,-2.3,0.5893786);
   grae->SetPointError(0,0.2,0.2,0.01548524,0.01529692);
   grae->SetPoint(1,-1.825,0.8270615);
   grae->SetPointError(1,0.275,0.275,0.009253898,0.00887233);
   grae->SetPoint(2,-1.3,0.7779732);
   grae->SetPointError(2,0.25,0.25,0.0100191,0.009705285);
   grae->SetPoint(3,-0.825,0.6678529);
   grae->SetPointError(3,0.225,0.225,0.01113099,0.01094618);
   grae->SetPoint(4,-0.35,0.7673263);
   grae->SetPointError(4,0.25,0.25,0.009050667,0.008813146);
   grae->SetPoint(5,0,0.3456331);
   grae->SetPointError(5,0.1,0.1,0.01584054,0.01619005);
   grae->SetPoint(6,0.35,0.7895587);
   grae->SetPointError(6,0.25,0.25,0.008586949,0.00833316);
   grae->SetPoint(7,0.825,0.6768296);
   grae->SetPointError(7,0.225,0.225,0.01030179,0.01013411);
   grae->SetPoint(8,1.3,0.7597908);
   grae->SetPointError(8,0.25,0.25,0.008624715,0.008417261);
   grae->SetPoint(9,1.825,0.8104483);
   grae->SetPointError(9,0.275,0.275,0.007847101,0.007610393);
   grae->SetPoint(10,2.3,0.6739392);
   grae->SetPointError(10,0.2,0.2,0.01092327,0.0107435);
   grae->Draw("pe");
   
   leg = new TLegend(0.6623563,0.2457627,0.9324713,0.4872881,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("Graph","0-5%","pl");

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
   grae->SetName("Graph4");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ffcc00");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#ffcc00");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(33);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,-2.3,0.5711606);
   grae->SetPointError(0,0.2,0.2,0.01581327,0.01567714);
   grae->SetPoint(1,-1.825,0.8255928);
   grae->SetPointError(1,0.275,0.275,0.009297352,0.008923039);
   grae->SetPoint(2,-1.3,0.7792372);
   grae->SetPointError(2,0.25,0.25,0.009984224,0.009669677);
   grae->SetPoint(3,-0.825,0.6585031);
   grae->SetPointError(3,0.225,0.225,0.01106212,0.01089507);
   grae->SetPoint(4,-0.35,0.7784441);
   grae->SetPointError(4,0.25,0.25,0.009032256,0.008776394);
   grae->SetPoint(5,0,0.3515774);
   grae->SetPointError(5,0.1,0.1,0.01614827,0.01649492);
   grae->SetPoint(6,0.35,0.7748555);
   grae->SetPointError(6,0.25,0.25,0.008729823,0.008498179);
   grae->SetPoint(7,0.825,0.6778575);
   grae->SetPointError(7,0.225,0.225,0.01003408,0.00986833);
   grae->SetPoint(8,1.3,0.751042);
   grae->SetPointError(8,0.25,0.25,0.009047271,0.008834678);
   grae->SetPoint(9,1.825,0.815195);
   grae->SetPointError(9,0.275,0.275,0.007628061,0.007394812);
   grae->SetPoint(10,2.3,0.6817187);
   grae->SetPointError(10,0.2,0.2,0.01088118,0.01069044);
   grae->Draw("pe");
   
   leg = new TLegend(0.6623563,0.2457627,0.9324713,0.4872881,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("Graph","0-5%","pl");

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
   grae->SetName("Graph5");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#6600ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#6600ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(29);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,-2.3,0.5854239);
   grae->SetPointError(0,0.2,0.2,0.0112319,0.01114758);
   grae->SetPoint(1,-1.825,0.8409121);
   grae->SetPointError(1,0.275,0.275,0.006348166,0.006148736);
   grae->SetPoint(2,-1.3,0.7818124);
   grae->SetPointError(2,0.25,0.25,0.006865982,0.006711512);
   grae->SetPoint(3,-0.825,0.664576);
   grae->SetPointError(3,0.225,0.225,0.007799431,0.007710208);
   grae->SetPoint(4,-0.35,0.7885771);
   grae->SetPointError(4,0.25,0.25,0.006295377,0.006159763);
   grae->SetPoint(5,0,0.3309308);
   grae->SetPointError(5,0.1,0.1,0.01096286,0.01114956);
   grae->SetPoint(6,0.35,0.7775087);
   grae->SetPointError(6,0.25,0.25,0.006178109,0.006058166);
   grae->SetPoint(7,0.825,0.6737833);
   grae->SetPointError(7,0.225,0.225,0.007179256,0.007098495);
   grae->SetPoint(8,1.3,0.7550098);
   grae->SetPointError(8,0.25,0.25,0.006322622,0.006213797);
   grae->SetPoint(9,1.825,0.8162555);
   grae->SetPointError(9,0.275,0.275,0.005428962,0.00530864);
   grae->SetPoint(10,2.3,0.6931052);
   grae->SetPointError(10,0.2,0.2,0.00764223,0.007538309);
   grae->Draw("pe");
   
   leg = new TLegend(0.6623563,0.2457627,0.9324713,0.4872881,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("Graph","0-5%","pl");

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
   grae->SetName("Graph6");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#00ffff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#00ffff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(34);
   grae->SetMarkerSize(1.2);
   grae->SetPoint(0,-2.3,0.581763);
   grae->SetPointError(0,0.2,0.2,0.01197954,0.01188706);
   grae->SetPoint(1,-1.825,0.8427266);
   grae->SetPointError(1,0.275,0.275,0.006933413,0.006690657);
   grae->SetPoint(2,-1.3,0.7807759);
   grae->SetPointError(2,0.25,0.25,0.007770149,0.007579823);
   grae->SetPoint(3,-0.825,0.6822521);
   grae->SetPointError(3,0.225,0.225,0.008721109,0.008592157);
   grae->SetPoint(4,-0.35,0.7965423);
   grae->SetPointError(4,0.25,0.25,0.006930113,0.006755977);
   grae->SetPoint(5,0,0.3460855);
   grae->SetPointError(5,0.1,0.1,0.01264523,0.01287103);
   grae->SetPoint(6,0.35,0.7826117);
   grae->SetPointError(6,0.25,0.25,0.006899063,0.006743264);
   grae->SetPoint(7,0.825,0.6775007);
   grae->SetPointError(7,0.225,0.225,0.00818572,0.008077995);
   grae->SetPoint(8,1.3,0.7616897);
   grae->SetPointError(8,0.25,0.25,0.006990362,0.006852089);
   grae->SetPoint(9,1.825,0.8297038);
   grae->SetPointError(9,0.275,0.275,0.005982006,0.00582089);
   grae->SetPoint(10,2.3,0.6906176);
   grae->SetPointError(10,0.2,0.2,0.008724526,0.008593284);
   grae->Draw("pe");
   
   leg = new TLegend(0.6623563,0.2457627,0.9324713,0.4872881,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("Graph","0-5%","pl");

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
   canv_neg_B->Modified();
   canv_neg_B->cd();
   canv_neg_B->SetSelected(canv_neg_B);
}
