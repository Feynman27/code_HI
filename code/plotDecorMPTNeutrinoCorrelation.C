{

    TFile* _file0 = TFile::Open("../MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.04.12.2013.root");

    gStyle->SetPalette(1); gROOT->LoadMacro("AtlasUtils.C");

    //TString cuts = "val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&pt>25.0&&ptcone20ID3/pt<0.1&&mt>40.0&&mt<9000.0&&((EF_mu10_MSonly_EFFS_L1ZDC&&EF_mu10_MSonly_EFFS_L1ZDC_Matched20)||(EF_mu10_MSonly_EFFS_L1TE10&&EF_mu10_MSonly_EFFS_L1TE10_Matched20)||(EF_mu10_MSonly_EFFS_L1TE20&&EF_mu10_MSonly_EFFS_L1TE20_Matched20))&&ZDY==0&&abs(eta)<2.4&&abs(mc_nu_gen_eta)<2.4";
    TString cuts = "val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&pt>25.0&&ptcone20ID3/pt<0.1&&mt>40.0&&mt<9000.0&&ZDY==0&&abs(eta)>0.1&&abs(eta)<2.4&&prompt==24&&ZDY==0";
    TString cutsCentral = cuts + "&&centrality>0.0&&centrality<0.05";
    TString cutsPeripheral = cuts + "&&centrality>0.4&&centrality<0.8";

    TCanvas *cCentral = new TCanvas("cCentral","cCentral",650,600);

    TH2D* hc = new TH2D("hc","hc",100,25.0,200.0,100,25.0,200.0);
    tree->Draw("nu_pt:mc_nu_gen_pt>>hc",cutsCentral,"colz");
    hc->GetYaxis()->SetRangeUser(25.0,115.0);
    hc->GetYaxis()->SetTitle("#slash{p_{T}}[GeV]");
    hc->GetXaxis()->SetRangeUser(25.,115.0);
    hc->GetXaxis()->SetTitle("p_{T}^{#nu}[GeV]");
    cCentral->SetLogz(1);

    TProfile *hc_p = (TProfile*) hc->ProfileX("hc_p");
    hc_p->Draw("pesame"); 
    gPad->Modified();
    gPad->Update();

    myText(0.215,0.833,kBlack,"0-5%");
    myText(0.17,0.77,kBlack,"p_{T}^{#mu}>25.0GeV,m_{T}>40GeV");
    myText(0.19,0.70,kBlack,"W#rightarrow#mu#nu MC11");

    TCanvas *cPeripheral = new TCanvas("cPeripheral","cPeripheral",650,600);
    TH2D* hp = new TH2D("hp","hp",100,25.,200.0,100,25.,200.0);
    tree->Draw("nu_pt:mc_nu_gen_pt>>hp",cutsPeripheral,"colz");
    hp->GetYaxis()->SetRangeUser(25.,115.0);
    hp->GetYaxis()->SetTitle("#slash{p_{T}}[GeV]");
    hp->GetXaxis()->SetRangeUser(25.0,115.0);
    hp->GetXaxis()->SetTitle("p_{T}^{#nu}[GeV]");
    cPeripheral->SetLogz(1);

    TProfile *hp_p = (TProfile*) hp->ProfileX("hp_p");
    hp_p->Draw("pesame"); 
    gPad->Modified();
    gPad->Update();

    myText(0.215,0.833,kBlack,"40-80%");
    myText(0.17,0.77,kBlack,"p_{T}^{#mu}>25.0GeV,m_{T}>40GeV");
    myText(0.19,0.70,kBlack,"W#rightarrow#mu#nu MC11");

}

