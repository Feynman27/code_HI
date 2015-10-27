{
    gROOT->LoadMacro("AtlasUtils.C");
    TFile *_file0 = TFile::Open("FitResults/signalChargeEtaDistributions_periodA.04.26.2013.root");
    TFile *_file1 = TFile::Open("FitResults/signalChargeEtaDistributions_periodB.04.26.2013.root");
    TFile *_file2 = TFile::Open("FitResults/signalChargeEtaDistributions_periodC.04.26.2013.root");

    int markerColor[] = {kRed-7,kOrange+2,kBlue,kRed,kMagenta,kAzure+10,kSpring, kGreen+2};


    _file0->cd();
    TGraphErrors* grEtaP0 = _file0->Get("grWpSystc");
    TGraphErrors* grEtaM0 = _file0->Get("grWmSystc");
    grEtaP0->SetMarkerColor(kRed);
    grEtaM0->SetMarkerColor(kRed);

    _file1->cd();
    TGraphErrors* grEtaP1 = _file1->Get("grWpSystc");
    TGraphErrors* grEtaM1 = _file1->Get("grWmSystc");
    grEtaP1->SetMarkerColor(kBlue);
    grEtaM1->SetMarkerColor(kBlue);

    _file2->cd();
    TGraphErrors* grEtaP2 = _file2->Get("grWpSystc");
    TGraphErrors* grEtaM2 = _file2->Get("grWmSystc");
    grEtaP2->SetMarkerColor(kGreen);
    grEtaM2->SetMarkerColor(kGreen);

    TCanvas* cP = new TCanvas("cP","cP",600,600);
    TH1D* hdummy = new TH1D("hdummy","hdummy",5,0.0,2.5);
    hdummy->GetYaxis()->SetRangeUser(0,1500.0);
    hdummy->GetXaxis()->SetTitle("|#eta_{#mu}|") ;
    hdummy->GetYaxis()->SetTitle("dN^{#mu}/d #eta");

    hdummy->Draw();
    grEtaP0->Draw("pesame");
    grEtaP1->Draw("pesame");
    grEtaP2->Draw("pesame");

    TLegend* leg = new TLegend(0.29,0.79,0.54,0.94);
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->AddEntry(grEtaP0,"Runs 193211-193641(L_{int}=47.4 #mub^{-1})","pe");
    leg->AddEntry(grEtaP1,"Runs 193655-193890(L_{int}=45.6 #mub^{-1})","pe");
    leg->AddEntry(grEtaP2,"Runs 194017-194382(L_{int}=47.5 #mub^{-1})","pe");
    leg->Draw();
    myText(0.29,0.48,kBlack,(char*)"#mu^{+}");
    myText(0.25,0.43,kBlack,(char*)"0-80%");
    myText(0.21,0.38,kBlack,(char*)"0.1<|#eta|<2.4");

    TCanvas* cM = new TCanvas("cM","cM",600,600);
    hdummy->Draw();
    grEtaM0->Draw("pesame");
    grEtaM1->Draw("pesame");
    grEtaM2->Draw("pesame");

    leg->Draw();

    myText(0.29,0.43,kBlack,(char*)"#mu^{-}");
    myText(0.25,0.38,kBlack,(char*)"0-80%");
    myText(0.21,0.33,kBlack,(char*)"0.1<|#eta|<2.4");

} 
