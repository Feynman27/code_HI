{
    TFile *_file0 = TFile::Open("mptSystematicErrorEtaDistros.root");
    gROOT->LoadMacro("AtlasUtils.C");
    TString date = ".06.16.2013";
    TCanvas *cPlus[6];
    TCanvas *cMinus[6];
    TString arrCent[6]={"0-5%","5-10%","10-15%","15-20%","20-40%","40-80%"};
    char cNameMinus[50], cNamePlus[50];
    for(int icent=0; icent<6; ++icent){

        sprintf(cNamePlus,"canvasPlus%i",icent);
        sprintf(cNameMinus,"canvasMinus%i",icent);
        cPlus[icent] = new TCanvas(cNamePlus,cNamePlus,600,600);
        cMinus[icent] = new TCanvas(cNameMinus,cNameMinus,600,600);

        TString mptRMSBase = "grMptRMSPercentErrorEtaDistroCent";
        TString mptUpBase = "grMptUpPercentErrorEtaDistroCent";
        TString mptDownBase = "grMptDownPercentErrorEtaDistroCent";

        mptRMSBase+=icent; mptRMSBase+="Charge";
        mptUpBase+=icent; mptUpBase+="Charge";
        mptDownBase+=icent; mptDownBase+="Charge";

        TGraphErrors* grMptRMSPlus = (TGraphErrors*)_file0->Get(mptRMSBase+"102");
        TGraphErrors* grMptRMSMinus = (TGraphErrors*)_file0->Get(mptRMSBase+"103");
        TGraphErrors* grMptUpPlus = (TGraphErrors*)_file0->Get(mptUpBase+"102");
        TGraphErrors* grMptUpMinus = (TGraphErrors*)_file0->Get(mptUpBase+"103");
        TGraphErrors* grMptDownPlus = (TGraphErrors*)_file0->Get(mptDownBase+"102");
        TGraphErrors* grMptDownMinus = (TGraphErrors*)_file0->Get(mptDownBase+"103");

        TH1F* hdummy = new TH1F("hdummy","hdummy",9,0.0,2.5);
        cPlus[icent]->cd();
        hdummy->GetYaxis()->SetRangeUser(-45.0,45.0);
        hdummy->GetXaxis()->SetTitle("|#eta|");
        hdummy->GetYaxis()->SetTitle("#slash{p_{T}} Relative Uncertainty [%]");
  	    
        TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
  	    leg->SetTextFont(gStyle->GetTextFont());
	    leg->SetTextSize(gStyle->GetTextSize());
	    leg->SetBorderSize(0);
	    leg->SetFillColor(0);

	    leg->AddEntry(grMptRMSPlus,"RMS","l");
	    leg->AddEntry(grMptUpPlus,"+1#sigma","l");
	    leg->AddEntry(grMptDownPlus,"-1#sigma","l");

        hdummy->Draw();
        leg->Draw();
        myText(0.2,0.86,kBlack,arrCent[icent]);
        myText(0.2,0.78,kBlack,"#mu^{+}");
        grMptRMSPlus->SetLineColor(kBlack);
        grMptRMSPlus->Draw("pesame");
        grMptUpPlus->SetLineColor(kRed);grMptUpPlus->SetLineStyle(kDashed);
        grMptUpPlus->Draw("pe same");
        grMptDownPlus->SetLineColor(kGreen);grMptDownPlus->SetLineStyle(kDashed);
        grMptDownPlus->Draw("pe same");
        TString sNamePlus = "mptSystematicSummaryMuPlusEtaDistrosCent"; sNamePlus+=icent;
        cPlus[icent]->Print(sNamePlus+date+".pdf");
        

        cMinus[icent]->cd();
        hdummy->Draw();
        myText(0.2,0.86,kBlack,arrCent[icent]);
        myText(0.2,0.78,kBlack,"#mu^{-}");
        leg->Draw();
        grMptRMSMinus->SetLineColor(kBlack);
        grMptRMSMinus->Draw("pesame");
        grMptUpMinus->SetLineColor(kRed);grMptUpMinus->SetLineStyle(kDashed);
        grMptUpMinus->Draw("pe same");
        grMptDownMinus->SetLineColor(kGreen);grMptDownMinus->SetLineStyle(kDashed);
        grMptDownMinus->Draw("pe same");
        TString sNameMinus = "mptSystematicSummaryMuMinusEtaDistrosCent"; sNameMinus+=icent;
        cMinus[icent]->Print(sNameMinus+date+".pdf");
    }
}
