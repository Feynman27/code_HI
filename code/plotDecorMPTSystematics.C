void plotDecorMPTSystematics(){

    TFile* fIn = TFile::Open("mptSystematicErrorEtaDistros.root");

    const int _nCentralityBins = 6;
    const int _nEtaBins = 9;

    std::vector<int> centralityBins;
    centralityBins.push_back(0);
    centralityBins.push_back(5);
    centralityBins.push_back(10);
    centralityBins.push_back(15);
    centralityBins.push_back(20);
    centralityBins.push_back(40);
    centralityBins.push_back(80);

    for(int icent = 0; icent < _nCentralityBins; icent++){

        TString sCent = ""; sCent+=centralityBins[icent]; sCent+="-"; sCent+=centralityBins[icent+1]; sCent+="%";

        TString sNameDown = "grMptDownPercentErrorEtaDistroCent";
        TString sNameUp = "grMptUpPercentErrorEtaDistroCent";

        sNameDown +=icent; sNameDown += "Charge";
        sNameUp+=icent; sNameUp+="Charge";

        TString sNameDownPlus = sNameDown; sNameDownPlus+=102;
        TString sNameUpPlus = sNameUp; sNameUpPlus+=102;

        TString sNameDownMinus = sNameDown; sNameDownMinus+=103;
        TString sNameUpMinus = sNameUp; sNameUpMinus+=103;
        
        grDownPlus = (TGraphErrors*)fIn->Get(sNameDownPlus);
        grDownMinus = (TGraphErrors*)fIn->Get(sNameDownMinus);

        grUpPlus = (TGraphErrors*)fIn->Get(sNameUpPlus);
        grUpMinus = (TGraphErrors*)fIn->Get(sNameUpMinus);

        grUpPlus->SetLineColor(2); grUpPlus->SetLineStyle(4); grUpPlus->SetLineWidth(3);
        grUpMinus->SetLineColor(2); grUpMinus->SetLineStyle(4); grUpMinus->SetLineWidth(3);

        grDownPlus->SetLineColor(4); grDownPlus->SetLineStyle(6); grDownPlus->SetLineWidth(3);
        grDownMinus->SetLineColor(4); grDownMinus->SetLineStyle(6); grDownMinus->SetLineWidth(3);

   TH1F *hDummy = new TH1F("hDummy","hDummy",100,0,2.63);
   hDummy->SetMinimum(-45.57846);
   hDummy->SetMaximum(45.57846);
   hDummy->SetDirectory(0);
   hDummy->SetStats(0);
   hDummy->SetLineWidth(2);
   hDummy->SetMarkerStyle(20);
   hDummy->SetMarkerSize(1.2);
   hDummy->GetXaxis()->SetTitle("|#eta|");
   hDummy->GetXaxis()->SetLabelFont(42);
   hDummy->GetXaxis()->SetLabelSize(0.05);
   hDummy->GetXaxis()->SetTitleSize(0.05);
   hDummy->GetXaxis()->SetTitleOffset(1.4);
   hDummy->GetXaxis()->SetTitleFont(42);
   hDummy->GetYaxis()->SetTitle("#slash{p_{T}} Relative Uncertainty [%]");
   hDummy->GetYaxis()->SetLabelFont(42);
   hDummy->GetYaxis()->SetLabelSize(0.05);
   hDummy->GetYaxis()->SetTitleSize(0.05);
   hDummy->GetYaxis()->SetTitleOffset(1.4);
   hDummy->GetYaxis()->SetTitleFont(42);

   TCanvas *c1 = new TCanvas("c1", "c1",362,134,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-0.5326582,-64.04062,2.796456,51.34789);
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
   
   hDummy->Draw();
   grDownPlus->Draw("pe same");
   grUpPlus->Draw("pe same");
   TLatex *   tex = new TLatex(0.7155172,0.8648305,"#mu^{+}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.08050848);
   tex->SetLineWidth(2);
   tex->Draw();
 
   TLegend *leg = new TLegend(0.1925287,0.7330508,0.4626437,0.8834746,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry(grUpPlus,"+1#sigma","l");
   entry->SetLineColor(2);
   entry->SetLineStyle(4);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry(grDownPlus,"-1#sigma","l");
   entry->SetLineColor(4);
   entry->SetLineStyle(6);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   tex = new TLatex(0.8045977,0.8690678,sCent);
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);


   TString sSaveNamePlus = "mptSystematicsMuPlusCent"; sSaveNamePlus+=icent;
   c1->Print(sSaveNamePlus+".pdf");
   c1->Print(sSaveNamePlus+".png");

   ///mu -
   TCanvas *c2 = new TCanvas("c2", "c2",362,134,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c2->Range(-0.5326582,-64.04062,2.796456,51.34789);
   c2->SetFillColor(0);
   c2->SetBorderMode(0);
   c2->SetBorderSize(2);
   c2->SetTickx(1);
   c2->SetTicky(1);
   c2->SetLeftMargin(0.16);
   c2->SetRightMargin(0.05);
   c2->SetTopMargin(0.05);
   c2->SetBottomMargin(0.16);
   c2->SetFrameBorderMode(0);
   c2->SetFrameBorderMode(0);
   
   hDummy->Draw();
   grDownMinus->Draw("pe same");
   grUpMinus->Draw("pe same");
   TLatex *   tex2 = new TLatex(0.7155172,0.8648305,"#mu^{-}");
   tex2->SetNDC();
   tex2->SetTextFont(42);
   tex2->SetTextSize(0.08050848);
   tex2->SetLineWidth(2);
   tex2->Draw();
 
   TLegend *leg = new TLegend(0.1925287,0.7330508,0.4626437,0.8834746,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry(grUpMinus,"+1#sigma","l");
   entry->SetLineColor(2);
   entry->SetLineStyle(4);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry(grDownMinus,"-1#sigma","l");
   entry->SetLineColor(4);
   entry->SetLineStyle(6);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   tex2 = new TLatex(0.8045977,0.8690678,sCent);
   tex2->SetNDC();
   tex2->SetTextFont(42);
   tex2->SetLineWidth(2);
   tex2->Draw();
   c2->Modified();
   c2->cd();
   c2->SetSelected(c2);


   TString sSaveNameMinus = "mptSystematicsMuMinusCent"; sSaveNameMinus+=icent;
   c2->Print(sSaveNameMinus+".pdf");
   c2->Print(sSaveNameMinus+".png");



    }
}
