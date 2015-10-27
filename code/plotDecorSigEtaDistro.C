{

    TString date = "05.01.2013";
    gROOT->LoadMacro("AtlasUtils.C");
    TFile *_file0 = TFile::Open("signalChargeEtaDistributions_Raw.04.17.2013.root");
    TFile *_file1 = TFile::Open("signalChargeEtaDistributions_cuts1.04.21.2013.root");
    TFile *_file2 = TFile::Open("signalChargeEtaDistributions_cuts2.04.21.2013.root");
    TFile *_file3 = TFile::Open("signalChargeEtaDistributions_cuts3.04.21.2013.root");
    TFile *_file4 = TFile::Open("signalChargeEtaDistributions_cuts4.04.21.2013.root");
    TFile *_file5 = TFile::Open("signalChargeEtaDistributions_cuts5.04.21.2013.root");
    TFile *_file6 = TFile::Open("signalChargeEtaDistributions_cuts6.04.21.2013.root");

    //TFile *_file1 = TFile::Open("signalChargeEtaDistributions_CorrectedNoBkgSub.04.17.2013.root");

    int markerStyle[] = {20,29};
    int markerColor[] = {kRed-7,kOrange+2,kBlue,kRed,kMagenta,kAzure+10,kSpring, kGreen+2};

    ///plot raw distro
    _file0->cd();

    TH1D* dataEtaDistroPlusCent0c = (TH1D*)dataEtaDistroPlusCent0->Clone("dataEtaDistroPlusCent0c");
    dataEtaDistroPlusCent0c->Scale(1.0,"width");
    dataEtaDistroPlusCent0c->SetMarkerColor(markerColor[0]);
    TH1D* dataEtaDistroMinusCent0c = (TH1D*)dataEtaDistroMinusCent0->Clone("dataEtaDistroMinusCent0c");
    dataEtaDistroMinusCent0c->Scale(1.0,"width");
    dataEtaDistroMinusCent0c->SetMarkerColor(markerColor[0]);

    _file1->cd();

    TH1D* dataEtaDistroPlusCent0Corr1 = (TH1D*)dataEtaDistroPlusCent0->Clone("dataEtaDistroPlusCent0Corr1");
    dataEtaDistroPlusCent0Corr1->Scale(1.0,"width");     
    dataEtaDistroPlusCent0Corr1->SetMarkerColor(markerColor[1]);
    TH1D* dataEtaDistroMinusCent0Corr1 = (TH1D*)dataEtaDistroMinusCent0->Clone("dataEtaDistroMinusCent0Corr1");
    dataEtaDistroMinusCent0Corr1->Scale(1.0,"width");
    dataEtaDistroMinusCent0Corr1->SetMarkerColor(markerColor[1]);

    _file2->cd();

    TH1D* dataEtaDistroPlusCent0Corr2 = (TH1D*)dataEtaDistroPlusCent0->Clone("dataEtaDistroPlusCent0Corr2");
    dataEtaDistroPlusCent0Corr2->Scale(1.0,"width");     
    dataEtaDistroPlusCent0Corr2->SetMarkerColor(markerColor[2]);
    TH1D* dataEtaDistroMinusCent0Corr2 = (TH1D*)dataEtaDistroMinusCent0->Clone("dataEtaDistroMinusCent0Corr2");
    dataEtaDistroMinusCent0Corr2->Scale(1.0,"width");
    dataEtaDistroMinusCent0Corr2->SetMarkerColor(markerColor[2]);

    _file3->cd();

    TH1D* dataEtaDistroPlusCent0Corr3 = (TH1D*)dataEtaDistroPlusCent0->Clone("dataEtaDistroPlusCent0Corr3");
    dataEtaDistroPlusCent0Corr3->Scale(1.0,"width");     
    dataEtaDistroPlusCent0Corr3->SetMarkerColor(markerColor[3]);
    TH1D* dataEtaDistroMinusCent0Corr3 = (TH1D*)dataEtaDistroMinusCent0->Clone("dataEtaDistroMinusCent0Corr3");
    dataEtaDistroMinusCent0Corr3->Scale(1.0,"width");
    dataEtaDistroMinusCent0Corr3->SetMarkerColor(markerColor[3]);

    _file4->cd();

    TH1D* dataEtaDistroPlusCent0Corr4 = (TH1D*)dataEtaDistroPlusCent0->Clone("dataEtaDistroPlusCent0Corr4");
    dataEtaDistroPlusCent0Corr4->Scale(1.0,"width");     
    dataEtaDistroPlusCent0Corr4->SetMarkerColor(markerColor[4]);
    TH1D* dataEtaDistroMinusCent0Corr4 = (TH1D*)dataEtaDistroMinusCent0->Clone("dataEtaDistroMinusCent0Corr4");
    dataEtaDistroMinusCent0Corr4->Scale(1.0,"width");
    dataEtaDistroMinusCent0Corr4->SetMarkerColor(markerColor[4]);


    _file5->cd();

    TH1D* dataEtaDistroPlusCent0Corr5 = (TH1D*)dataEtaDistroPlusCent0->Clone("dataEtaDistroPlusCent0Corr5");
    dataEtaDistroPlusCent0Corr5->Scale(1.0,"width");     
    dataEtaDistroPlusCent0Corr5->SetMarkerColor(markerColor[5]);
    TH1D* dataEtaDistroMinusCent0Corr5 = (TH1D*)dataEtaDistroMinusCent0->Clone("dataEtaDistroMinusCent0Corr5");
    dataEtaDistroMinusCent0Corr5->Scale(1.0,"width");
    dataEtaDistroMinusCent0Corr5->SetMarkerColor(markerColor[5]);

    _file6->cd();

    TH1D* dataEtaDistroPlusCent0Corr6 = (TH1D*)dataEtaDistroPlusCent0->Clone("dataEtaDistroPlusCent0Corr6");
    dataEtaDistroPlusCent0Corr6->Scale(1.0,"width");     
    dataEtaDistroPlusCent0Corr6->SetMarkerColor(markerColor[6]);
    dataEtaDistroPlusCent0Corr6->SetLineColor(markerColor[6]);
    TH1D* dataEtaDistroMinusCent0Corr6 = (TH1D*)dataEtaDistroMinusCent0->Clone("dataEtaDistroMinusCent0Corr6");
    dataEtaDistroMinusCent0Corr6->Scale(1.0,"width");
    dataEtaDistroMinusCent0Corr6->SetMarkerColor(markerColor[6]);
    dataEtaDistroMinusCent0Corr6->SetLineColor(markerColor[6]);

    TH1D* mcEtaDistroPlusCent0c = (TH1D*)mcTruthEtaDistroPlusCent0->Clone("mcEtaDistroPlusCent0c");
    mcEtaDistroPlusCent0c->Scale(1.0,"width");
    mcEtaDistroPlusCent0c->Scale(dataEtaDistroPlusCent0Corr6->Integral()/mcEtaDistroPlusCent0c->Integral());
    mcEtaDistroPlusCent0c->SetMarkerColor(markerColor[7]);
    mcEtaDistroPlusCent0c->SetMarkerStyle(markerStyle[1]);
    mcEtaDistroPlusCent0c->SetMarkerSize(2.5);

    TH1D* mcEtaDistroMinusCent0c = (TH1D*)mcTruthEtaDistroMinusCent0->Clone("mcEtaDistroMinusCent0c");
    mcEtaDistroMinusCent0c->Scale(1.0,"width");
    mcEtaDistroMinusCent0c->Scale(dataEtaDistroMinusCent0Corr6->Integral()/mcEtaDistroMinusCent0c->Integral());
    mcEtaDistroMinusCent0c->SetMarkerColor(markerColor[7]);
    mcEtaDistroMinusCent0c->SetMarkerStyle(markerStyle[1]);
    mcEtaDistroMinusCent0c->SetMarkerSize(2.5);

    TCanvas* c0 = new TCanvas();
    dataEtaDistroPlusCent0Corr1->GetYaxis()->SetRangeUser(0,1250.0);
    dataEtaDistroPlusCent0Corr1->GetXaxis()->SetTitle("|#eta_{#mu}|") ;
    dataEtaDistroPlusCent0Corr1->GetYaxis()->SetTitle("dN^{#mu}/d #eta");

    dataEtaDistroPlusCent0Corr1->Draw("pe");
    dataEtaDistroPlusCent0Corr2->Draw("pesame");
    dataEtaDistroPlusCent0Corr3->Draw("pesame");
    dataEtaDistroPlusCent0Corr4->Draw("pesame");
    dataEtaDistroPlusCent0Corr5->Draw("pesame");
    dataEtaDistroPlusCent0Corr6->Draw("pesame");
    dataEtaDistroPlusCent0c->Draw("pesame");
    mcEtaDistroPlusCent0c->Draw("pesame");

    //myText(0.3,0.67,kBlack,(char*)"#epsilon_{W#rightarrow#mu} = #frac{N_{k}}{N_{gen}#left{p_{T}^{#mu,gen}>25,p_{T}^{#nu,gen}>25,m_{T}^{#mu#nu}>40#right}}");
    TLatex l;
    l.SetNDC();
    l.SetTextSize(0.033);
    l.DrawLatex(0.51,0.69,"#epsilon_{W#rightarrow#mu} =
    #frac{N_{k}}{N_{gen}#left{p_{T}^{#mu,gen}>25,p_{T}^{#nu,gen}>25,m_{T}^{#mu#nu}>40#right}}");
    TLegend* leg = new TLegend(0.17,0.74,0.44,0.93);
    //TLegend* leg = new TLegend(0.2,0.76,0.49,0.9);
    leg->SetTextFont(gStyle->GetTextFont());
    //leg->SetTextSize(gStyle->GetTextSize());
    leg->SetTextSize(0.032);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    leg->AddEntry(mcEtaDistroPlusCent0c,"PYTHIA","pe");
    leg->AddEntry(dataEtaDistroPlusCent0c,"Data (Raw)","pe");
    leg->AddEntry(dataEtaDistroPlusCent0Corr1,"k_{0}:|#eta|<2.4,p_{T}^{#mu,rec}>25,#slash{p_{T}}>25,m_{T}>40","pe");

    TLegend* leg2 = new TLegend(0.55,0.74,0.82,0.92);
    //TLegend* leg2 = new TLegend(0.52,0.75,0.79,0.93);
    leg2->SetTextFont(gStyle->GetTextFont());
    leg2->SetTextSize(gStyle->GetTextSize());
    leg2->SetTextSize(0.032);
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);

    //leg2->AddEntry(dataEtaDistroPlusCent0Corr2,"k_{1}:k_{0}+PS","pe");
    //leg2->AddEntry(dataEtaDistroPlusCent0Corr3,"k_{2}:k_{1}+Z veto","pe");
    leg2->AddEntry(dataEtaDistroPlusCent0Corr4,"k_{1}:k_{0}+PS","pe");
    leg2->AddEntry(dataEtaDistroPlusCent0Corr5,"k_{2}:k_{1}+Z veto","pe");
    leg2->AddEntry(dataEtaDistroPlusCent0Corr6,"k_{3}:k_{2}+isolation","pe");
    leg->Draw();
    leg2->Draw();



    myText(0.22,0.27,kBlack,(char*)"#mu^{+}");
    myText(0.28,0.28,kBlack,(char*)"0-5%");

    c0->Print("dataEtaDistroPlusCent0."+date+".pdf");

    TCanvas* c1 = new TCanvas();
    dataEtaDistroMinusCent0Corr1->GetYaxis()->SetRangeUser(0,1250.0);
    dataEtaDistroMinusCent0Corr1->GetXaxis()->SetTitle("|#eta_{#mu}|") ;
    dataEtaDistroMinusCent0Corr1->GetYaxis()->SetTitle("dN^{#mu}/d #eta");

    dataEtaDistroMinusCent0Corr1->Draw("pe");
    dataEtaDistroMinusCent0Corr2->Draw("pesame");
    dataEtaDistroMinusCent0Corr3->Draw("pesame");
    dataEtaDistroMinusCent0Corr4->Draw("pesame");
    dataEtaDistroMinusCent0Corr5->Draw("pesame");
    dataEtaDistroMinusCent0Corr6->Draw("pesame");
    dataEtaDistroMinusCent0c->Draw("pesame");
    mcEtaDistroMinusCent0c->Draw("pesame");


    l.SetTextSize(0.033);
    l.DrawLatex(0.51,0.69,"#epsilon_{W#rightarrow#mu} =
    #frac{N_{k}}{N_{gen}#left{p_{T}^{#mu,gen}>25,p_{T}^{#nu,gen}>25,m_{T}^{#mu#nu}>40#right}}");
    leg->Draw();
    leg2->Draw();
    myText(0.22,0.27,kBlack,(char*)"#mu^{-}");
    myText(0.28,0.28,kBlack,(char*)"0-5%");

    c1->Print("dataEtaDistroMinusCent0."+date+".pdf");
   ///40-80%
    ///plot raw distro
    _file0->cd();

    TH1D* dataEtaDistroPlusCent5c = (TH1D*)dataEtaDistroPlusCent5->Clone("dataEtaDistroPlusCent5c");
    dataEtaDistroPlusCent5c->Scale(1.0,"width");
    dataEtaDistroPlusCent5c->SetMarkerColor(markerColor[0]);
    TH1D* dataEtaDistroMinusCent5c = (TH1D*)dataEtaDistroMinusCent5->Clone("dataEtaDistroMinusCent5c");
    dataEtaDistroMinusCent5c->Scale(1.0,"width");
    dataEtaDistroMinusCent5c->SetMarkerColor(markerColor[0]);

    _file1->cd();

    TH1D* dataEtaDistroPlusCent5Corr1 = (TH1D*)dataEtaDistroPlusCent5->Clone("dataEtaDistroPlusCent5Corr1");
    dataEtaDistroPlusCent5Corr1->Scale(1.0,"width");     
    dataEtaDistroPlusCent5Corr1->SetMarkerColor(markerColor[1]);
    TH1D* dataEtaDistroMinusCent5Corr1 = (TH1D*)dataEtaDistroMinusCent5->Clone("dataEtaDistroMinusCent5Corr1");
    dataEtaDistroMinusCent5Corr1->Scale(1.0,"width");
    dataEtaDistroMinusCent5Corr1->SetMarkerColor(markerColor[1]);

    _file2->cd();

    TH1D* dataEtaDistroPlusCent5Corr2 = (TH1D*)dataEtaDistroPlusCent5->Clone("dataEtaDistroPlusCent5Corr2");
    dataEtaDistroPlusCent5Corr2->Scale(1.0,"width");     
    dataEtaDistroPlusCent5Corr2->SetMarkerColor(markerColor[2]);
    TH1D* dataEtaDistroMinusCent5Corr2 = (TH1D*)dataEtaDistroMinusCent5->Clone("dataEtaDistroMinusCent5Corr2");
    dataEtaDistroMinusCent5Corr2->Scale(1.0,"width");
    dataEtaDistroMinusCent5Corr2->SetMarkerColor(markerColor[2]);

    _file3->cd();

    TH1D* dataEtaDistroPlusCent5Corr3 = (TH1D*)dataEtaDistroPlusCent5->Clone("dataEtaDistroPlusCent5Corr3");
    dataEtaDistroPlusCent5Corr3->Scale(1.0,"width");     
    dataEtaDistroPlusCent5Corr3->SetMarkerColor(markerColor[3]);
    TH1D* dataEtaDistroMinusCent5Corr3 = (TH1D*)dataEtaDistroMinusCent5->Clone("dataEtaDistroMinusCent5Corr3");
    dataEtaDistroMinusCent5Corr3->Scale(1.0,"width");
    dataEtaDistroMinusCent5Corr3->SetMarkerColor(markerColor[3]);

    _file4->cd();

    TH1D* dataEtaDistroPlusCent5Corr4 = (TH1D*)dataEtaDistroPlusCent5->Clone("dataEtaDistroPlusCent5Corr4");
    dataEtaDistroPlusCent5Corr4->Scale(1.0,"width");     
    dataEtaDistroPlusCent5Corr4->SetMarkerColor(markerColor[4]);
    TH1D* dataEtaDistroMinusCent5Corr4 = (TH1D*)dataEtaDistroMinusCent5->Clone("dataEtaDistroMinusCent5Corr4");
    dataEtaDistroMinusCent5Corr4->Scale(1.0,"width");
    dataEtaDistroMinusCent5Corr4->SetMarkerColor(markerColor[4]);


    _file5->cd();

    TH1D* dataEtaDistroPlusCent5Corr5 = (TH1D*)dataEtaDistroPlusCent5->Clone("dataEtaDistroPlusCent5Corr5");
    dataEtaDistroPlusCent5Corr5->Scale(1.0,"width");     
    dataEtaDistroPlusCent5Corr5->SetMarkerColor(markerColor[5]);
    TH1D* dataEtaDistroMinusCent5Corr5 = (TH1D*)dataEtaDistroMinusCent5->Clone("dataEtaDistroMinusCent5Corr5");
    dataEtaDistroMinusCent5Corr5->Scale(1.0,"width");
    dataEtaDistroMinusCent5Corr5->SetMarkerColor(markerColor[5]);

    _file6->cd();

    TH1D* dataEtaDistroPlusCent5Corr6 = (TH1D*)dataEtaDistroPlusCent5->Clone("dataEtaDistroPlusCent5Corr6");
    dataEtaDistroPlusCent5Corr6->Scale(1.0,"width");     
    dataEtaDistroPlusCent5Corr6->SetMarkerColor(markerColor[6]);
    dataEtaDistroPlusCent5Corr6->SetLineColor(markerColor[6]);
    TH1D* dataEtaDistroMinusCent5Corr6 = (TH1D*)dataEtaDistroMinusCent5->Clone("dataEtaDistroMinusCent5Corr6");
    dataEtaDistroMinusCent5Corr6->Scale(1.0,"width");
    dataEtaDistroMinusCent5Corr6->SetMarkerColor(markerColor[6]);
    dataEtaDistroMinusCent5Corr6->SetLineColor(markerColor[6]);

    TH1D* mcEtaDistroPlusCent5c = (TH1D*)mcTruthEtaDistroPlusCent5->Clone("mcEtaDistroPlusCent5c");
    mcEtaDistroPlusCent5c->Scale(1.0,"width");
    mcEtaDistroPlusCent5c->Scale(dataEtaDistroPlusCent5Corr6->Integral()/mcEtaDistroPlusCent5c->Integral());
    mcEtaDistroPlusCent5c->SetMarkerColor(markerColor[7]);
    mcEtaDistroPlusCent5c->SetMarkerStyle(markerStyle[1]);
    mcEtaDistroPlusCent5c->SetMarkerSize(2.5);

    TH1D* mcEtaDistroMinusCent5c = (TH1D*)mcTruthEtaDistroMinusCent5->Clone("mcEtaDistroMinusCent5c");
    mcEtaDistroMinusCent5c->Scale(1.0,"width");
    mcEtaDistroMinusCent5c->Scale(dataEtaDistroMinusCent5Corr6->Integral()/mcEtaDistroMinusCent5c->Integral());
    mcEtaDistroMinusCent5c->SetMarkerColor(markerColor[7]);
    mcEtaDistroMinusCent5c->SetMarkerStyle(markerStyle[1]);
    mcEtaDistroMinusCent5c->SetMarkerSize(2.5);

    TCanvas* c2 = new TCanvas();
    dataEtaDistroPlusCent5Corr1->GetYaxis()->SetRangeUser(0,550.0);
    dataEtaDistroPlusCent5Corr1->GetXaxis()->SetTitle("|#eta_{#mu}|") ;
    dataEtaDistroPlusCent5Corr1->GetYaxis()->SetTitle("dN^{#mu}/d #eta");

    dataEtaDistroPlusCent5Corr1->Draw("pe");
    dataEtaDistroPlusCent5Corr2->Draw("pesame");
    dataEtaDistroPlusCent5Corr3->Draw("pesame");
    dataEtaDistroPlusCent5Corr4->Draw("pesame");
    dataEtaDistroPlusCent5Corr5->Draw("pesame");
    dataEtaDistroPlusCent5Corr6->Draw("pesame");
    dataEtaDistroPlusCent5c->Draw("pesame");
    mcEtaDistroPlusCent5c->Draw("pesame");

    //TLegend* leg = new TLegend(0.178,0.769,0.448,0.919);
    /*TLegend* leg = new TLegend(0.393,0.775,0.664,0.926);
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->AddEntry(mcEtaDistroPlusCent5c,"W#rightarrow#mu Truth","pe");
    leg->AddEntry(dataEtaDistroPlusCent5c,"Data (Raw)","pe");
    leg->AddEntry(dataEtaDistroPlusCent5Corr1,"Data (Cuts 1)","pe");
    leg->AddEntry(dataEtaDistroPlusCent5Corr2,"Data (Cuts 2)","pe");
*/
    //TLegend* leg2 = new TLegend(0.451,0.756,0.721,0.907);
/*    TLegend* leg2 = new TLegend(0.656,0.769,0.927,0.922);
    leg2->SetTextFont(gStyle->GetTextFont());
    leg2->SetTextSize(gStyle->GetTextSize());
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);

    leg2->AddEntry(dataEtaDistroPlusCent5Corr3,"Data (Cuts 3)","pe");
    leg2->AddEntry(dataEtaDistroPlusCent5Corr4,"Data (Cuts 4)","pe");
    leg2->AddEntry(dataEtaDistroPlusCent5Corr5,"Data (Cuts 5)","pe");
    leg2->AddEntry(dataEtaDistroPlusCent5Corr6,"Data (Cuts 6)","pe");
*/

    l.SetTextSize(0.033);
    l.DrawLatex(0.51,0.69,"#epsilon_{W#rightarrow#mu} =
    #frac{N_{k}}{N_{gen}#left{p_{T}^{#mu,gen}>25,p_{T}^{#nu,gen}>25,m_{T}^{#mu#nu}>40#right}}");
    leg->Draw();
    leg2->Draw();

    myText(0.22,0.27,kBlack,(char*)"#mu^{+}");
    myText(0.28,0.28,kBlack,(char*)"40-80%");

    c2->Print("dataEtaDistroPlusCent5."+date+".pdf");

    TCanvas* c3 = new TCanvas();


    dataEtaDistroMinusCent5Corr1->GetYaxis()->SetRangeUser(0,550.0);
    dataEtaDistroMinusCent5Corr1->GetXaxis()->SetTitle("|#eta_{#mu}|") ;
    dataEtaDistroMinusCent5Corr1->GetYaxis()->SetTitle("dN^{#mu}/d #eta");

    dataEtaDistroMinusCent5Corr1->Draw("pe");
    dataEtaDistroMinusCent5Corr2->Draw("pesame");
    dataEtaDistroMinusCent5Corr3->Draw("pesame");
    dataEtaDistroMinusCent5Corr4->Draw("pesame");
    dataEtaDistroMinusCent5Corr5->Draw("pesame");
    dataEtaDistroMinusCent5Corr6->Draw("pesame");
    dataEtaDistroMinusCent5c->Draw("pesame");
    mcEtaDistroMinusCent5c->Draw("pesame");

    l.SetTextSize(0.033);
    l.DrawLatex(0.51,0.69,"#epsilon_{W#rightarrow#mu} =
    #frac{N_{k}}{N_{gen}#left{p_{T}^{#mu,gen}>25,p_{T}^{#nu,gen}>25,m_{T}^{#mu#nu}>40#right}}");
    leg->Draw();
    leg2->Draw();
    myText(0.22,0.27,kBlack,(char*)"#mu^{-}");
    myText(0.28,0.28,kBlack,(char*)"40-80%");
    c3->Print("dataEtaDistroMinusCent5."+date+".pdf");


}
