
{
    

    TFile* fIn = new TFile("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wtaunu_mu.NTUP_TRUTH.07.02.2013.root","read");
    TFile* outFile = new TFile("smearedMuonPtFromTau.root","recreate");
    TTree* tree = (TTree*)fIn->Get("truth");


    int mc_n;
    std::vector<int>* mc_pdgId = new std::vector<int>();
    TH1I* hMu = new TH1I("hMu","hMu",50,0,10);
    TH1I* hTau = new TH1I("hTau","hTau",500,0,10);
    TH2I* h2D = new TH2I("h2D","h2D",50,0,10,50,0,10);

    tree->SetBranchAddress("mc_pdgId",&mc_pdgId);
    tree->SetBranchAddress("mc_n",&mc_n);
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("mc_pdgId",1);
    tree->SetBranchStatus("mc_n",1);
    for(int iev=0; iev<tree->GetEntries(); ++iev){
        tree->GetEntry(iev);

        int nmu = 0, ntau = 0;
        for(int igen=0; igen<mc_n; ++igen){
            if(fabs(mc_pdgId->at(igen))==13) ++nmu;
            if(fabs(mc_pdgId->at(igen))==15) ++ntau;
        }
        hMu->Fill(nmu);
        hTau->Fill(ntau);
        h2D->Fill(ntau,nmu);
    }

    TCanvas* cmu = new TCanvas("cmu","cmu",600,600);
    hMu->Draw("hist");
    TCanvas* ctau = new TCanvas("ctau","ctau",600,600);
    hTau->Draw("hist");
    TCanvas* ccorr = new TCanvas("ccorr","ccorr",600,600);
    gStyle->SetPalette(1);
    h2D->Draw("colz");
}
