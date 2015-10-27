#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include <cmath>

double findMaxVariation(double err1, double err2, double err3, double err4){

    string sMax = "none";
    std::cout << " Variations: " << err1 << " " << err2 << " " << err3 << " " << err4 << std::endl;
    double max = -1.0;
    if(err1>max) max = err1;
    if(err2>max) max = err2;
    if(err3>max) max = err3;
    if(err4>max) max = err4;
    if(max==err1) sMax = "ud";
    if(max==err2) sMax = "du";
    if(max==err3) sMax = "uu";
    if(max==err4) sMax = "dd";
    std::cout << "Maximum: " << max << ", " << sMax << std::endl;
    return max;
}

void getAsymmetry(TH1F* hAsymm, TH1F* h1,TH1F* h2){

    //hAsymm = (TH1F*)h1->Clone();
    //hAsymm->Sumw2();
    TH1F* top = (TH1F*)hAsymm->Clone();
    TH1F* bottom = (TH1F*)hAsymm->Clone();
    top->Add(h1,h2,1.0,-1.0);
    bottom->Add(h1,h2,1.0,1.0);
    hAsymm->Divide(top,bottom);
    Int_t xmax = hAsymm->GetNbinsX();
    Double_t a,b,bot,error,da,db;
    for(Int_t i = 1; i<=xmax; ++i){
        
        a = h1->GetBinContent(i);
        b = h2->GetBinContent(i);
        bot = bottom->GetBinContent(i);
        if(bot<1e-6){}
        else{
            
            da = h1->GetBinError(i);
            db = h2->GetBinError(i);
            error = 2*TMath::Sqrt(a*a*db*db + b*b*da*da )/(bot*bot);
            hAsymm->SetBinError(i,error);
        }
    }

    delete top;
    delete bottom;
    
}

void fillHistos(TTree* t, TH1F* h, double cs){

    int nMuonsNt;
    float muonPtNt[30], muonEtaNt[30],muonMtNt[30];
    float neutrinoPtNt[30];
    int muonMotherNt[30], neutrinoMotherNt[30];

    // set branch addresses
    t->SetBranchAddress("nMuons", &nMuonsNt);
    t->SetBranchAddress("muonPt", &muonPtNt);
    t->SetBranchAddress("neutrinoPt", &neutrinoPtNt);
    t->SetBranchAddress("muonEta", &muonEtaNt);
    t->SetBranchAddress("mt", &muonMtNt);
    t->SetBranchAddress("muonMother", &muonMotherNt);
    t->SetBranchAddress("neutrinoMother", &neutrinoMotherNt);

    // turn on branches
    t->SetBranchStatus("*",0);
    t->SetBranchStatus("nMuons", 1);
    t->SetBranchStatus("muonPt", 1);
    t->SetBranchStatus("neutrinoPt", 1);
    t->SetBranchStatus("muonEta", 1);
    t->SetBranchStatus("mt", 1);
    t->SetBranchStatus("muonMother", 1);
    t->SetBranchStatus("neutrinoMother", 1);

    int nFidWs = 0;
    for(int iev=0; iev<t->GetEntries(); ++iev){

        t->GetEntry(iev);
        if(iev%10000==0) std::cout << "Event: " << iev << std::endl;

        for(int imu=0; imu<nMuonsNt; ++imu){
            if(fabs(muonMotherNt[imu]) == 24) ++nFidWs;
            if( fabs(muonMotherNt[imu]) == 24 &&
                fabs(neutrinoMotherNt[imu]) == 24 &&
                muonPtNt[imu]>25.0&&
                fabs(muonEtaNt[imu])<=2.5 &&
                neutrinoPtNt[imu]>25.0 &&
                muonMtNt[imu]>40.0
              ){
                h->Fill(fabs(muonEtaNt[imu])) ;
            } //selection
        }//imu
    }//iev

    h->Scale(1./nFidWs*cs);
    h->Scale(1,"width");

}


void CT10_plotter(){

    TString sDir = "/usatlas/u/tbales/scratch/";
    // central values
    TFile* _fPos = new TFile(sDir+"PowhegPythia8_CT10NLO_W+_07.09.2014.root","READ");
    TFile* _fNeg = new TFile(sDir+"PowhegPythia8_CT10NLO_W-_07.09.2014.root","READ");
    // scale uncertainties (muR,muF)
    TString sScaleDir = sDir/*+"scaleUncert/"*/;
    TFile* _fPosUD = new TFile(sScaleDir+"PowhegPythia8_CT10NLO_muRup_muFdown_W+_07.09.2014.root","READ");
    TFile* _fNegUD = new TFile(sScaleDir+"PowhegPythia8_CT10NLO_muRup_muFdown_W-_07.09.2014.root","READ");
    TFile* _fPosDU = new TFile(sScaleDir+"PowhegPythia8_CT10NLO_muRdown_muFup_W+_07.09.2014.root","READ");
    TFile* _fNegDU = new TFile(sScaleDir+"PowhegPythia8_CT10NLO_muRdown_muFup_W-_07.09.2014.root","READ");
    TFile* _fPosUU = new TFile(sScaleDir+"PowhegPythia8_CT10NLO_muRup_muFup_W+_07.09.2014.root","READ");
    TFile* _fNegUU = new TFile(sScaleDir+"PowhegPythia8_CT10NLO_muRup_muFup_W-_07.09.2014.root","READ");
    TFile* _fPosDD = new TFile(sScaleDir+"PowhegPythia8_CT10NLO_muRdown_muFdown_W+_07.09.2014.root","READ");
    TFile* _fNegDD = new TFile(sScaleDir+"PowhegPythia8_CT10NLO_muRdown_muFdown_W-_07.09.2014.root","READ");

    TFile* _fOut = new TFile("CT10_scaling.root","RECREATE");

    TTree* _tPos = (TTree*)_fPos->Get("T");
    TTree* _tNeg = (TTree*)_fNeg->Get("T");
    TTree* _tPosUD = (TTree*)_fPosUD->Get("T");
    TTree* _tNegUD = (TTree*)_fNegUD->Get("T");
    TTree* _tPosDU = (TTree*)_fPosDU->Get("T");
    TTree* _tNegDU = (TTree*)_fNegDU->Get("T");
    TTree* _tPosUU = (TTree*)_fPosUU->Get("T");
    TTree* _tNegUU = (TTree*)_fNegUU->Get("T");
    TTree* _tPosDD = (TTree*)_fPosDD->Get("T");
    TTree* _tNegDD = (TTree*)_fNegDD->Get("T");

    TH1F* hPos=NULL, *hNeg=NULL, *hAsymm=NULL;
    TH1F* hPosUD=NULL, *hNegUD=NULL, *hAsymmUD=NULL;
    TH1F* hPosDU=NULL, *hNegDU=NULL, *hAsymmDU=NULL;
    TH1F* hPosUU=NULL, *hNegUU=NULL, *hAsymmUU=NULL;
    TH1F* hPosDD=NULL, *hNegDD=NULL, *hAsymmDD=NULL;

    double WCrossSection_Pos = 2.181e-9/64.0e-3*1.0e9;
    double WCrossSection_Neg= 1.284e-9/64.0e-3*1.0e9;
    double WCrossSection_PosUD = 2.109e-9/64.0e-3*1.0e9;
    double WCrossSection_NegUD= 1.242e-9/64.0e-3*1.0e9;
    double WCrossSection_PosDU = 2.253e-9/64.0e-3*1.0e9;
    double WCrossSection_NegDU = 1.327e-9/64.0e-3*1.0e9;
    double WCrossSection_PosUU = 2.1798e-9/64.0e-3*1.0e9;
    double WCrossSection_NegUU = 1.2855e-9/64.0e-3*1.0e9;
    double WCrossSection_PosDD = 2.1929e-9/64.0e-3*1.0e9;
    double WCrossSection_NegDD = 1.2898e-9/64.0e-3*1.0e9;
    double wt_pp = 0.155;
    std::vector<float> vEtaBins;
    vEtaBins.push_back(0.0);
    vEtaBins.push_back(0.1);
    vEtaBins.push_back(0.35);
    vEtaBins.push_back(0.6);
    vEtaBins.push_back(0.8);
    vEtaBins.push_back(1.05);
    vEtaBins.push_back(1.37);
    vEtaBins.push_back(1.52);
    vEtaBins.push_back(1.85);
    vEtaBins.push_back(2.01);
    vEtaBins.push_back(2.5);
    const int nEtaBins = vEtaBins.size()-1;

    float arrEtaBins[nEtaBins+1];
    for(int ieta=0; ieta<vEtaBins.size(); ++ieta){

        arrEtaBins[ieta] = vEtaBins[ieta];
    }//ieta

    hPos = new TH1F("hPos","mu^{+},CT10",nEtaBins,arrEtaBins);
    hPos->Sumw2();
    hNeg = new TH1F("hNeg","mu^{-},CT10",nEtaBins,arrEtaBins);
    hNeg->Sumw2();
    hAsymm = new TH1F("hAsymm","Asymmetry,CT10",nEtaBins,arrEtaBins);
    hAsymm->Sumw2();

    hPosUD = new TH1F("hPosUD","mu^{+},CT10",nEtaBins,arrEtaBins);
    hPosUD->Sumw2();
    hNegUD = new TH1F("hNegUD","mu^{-},CT10",nEtaBins,arrEtaBins);
    hNegUD->Sumw2();
    hAsymmUD = new TH1F("hAsymmUD","AsymmetryUD,CT10",nEtaBins,arrEtaBins);
    hAsymmUD->Sumw2();

    hPosDU = new TH1F("hPosDU","mu^{+},CT10",nEtaBins,arrEtaBins);
    hPosDU->Sumw2();
    hNegDU = new TH1F("hNegDU","mu^{-},CT10",nEtaBins,arrEtaBins);
    hNegDU->Sumw2();
    hAsymmDU = new TH1F("hAsymmDU","AsymmetryDU,CT10",nEtaBins,arrEtaBins);
    hAsymmDU->Sumw2();

    hPosUU = new TH1F("hPosUU","mu^{+},CT10",nEtaBins,arrEtaBins);
    hPosUU->Sumw2();
    hNegUU = new TH1F("hNegUU","mu^{-},CT10",nEtaBins,arrEtaBins);
    hNegUU->Sumw2();
    hAsymmUU = new TH1F("hAsymmUU","AsymmetryUU,CT10",nEtaBins,arrEtaBins);
    hAsymmUU->Sumw2();

    hPosDD = new TH1F("hPosDD","mu^{+},CT10",nEtaBins,arrEtaBins);
    hPosDD->Sumw2();
    hNegDD = new TH1F("hNegDD","mu^{-},CT10",nEtaBins,arrEtaBins);
    hNegDD->Sumw2();
    hAsymmDD = new TH1F("hAsymmDD","AsymmetryDD,CT10",nEtaBins,arrEtaBins);
    hAsymmDD->Sumw2();

    fillHistos(_tPos,hPos,WCrossSection_Pos);
    fillHistos(_tNeg,hNeg,WCrossSection_Neg);
    fillHistos(_tPosUD,hPosUD,WCrossSection_PosUD);
    fillHistos(_tNegUD,hNegUD,WCrossSection_NegUD);
    fillHistos(_tPosDU,hPosDU,WCrossSection_PosDU);
    fillHistos(_tNegDU,hNegDU,WCrossSection_NegDU);
    fillHistos(_tPosUU,hPosUU,WCrossSection_PosUU);
    fillHistos(_tNegUU,hNegUU,WCrossSection_NegUU);
    fillHistos(_tPosDD,hPosDD,WCrossSection_PosDD);
    fillHistos(_tNegDD,hNegDD,WCrossSection_NegDD);

    getAsymmetry(hAsymm,hPos,hNeg);
    getAsymmetry(hAsymmUD,hPosUD,hNegUD);
    getAsymmetry(hAsymmDU,hPosDU,hNegDU);
    getAsymmetry(hAsymmUU,hPosUU,hNegUU);
    getAsymmetry(hAsymmDD,hPosDD,hNegDD);
    // set errors 
    for(int ibin=1; ibin<=hAsymm->GetNbinsX(); ++ibin){

        double errorPos = 0.0;
        double errorPosUD = fabs(hPosUD->GetBinContent(ibin)-hPos->GetBinContent(ibin));
        double errorPosDU = fabs(hPosDU->GetBinContent(ibin)-hPos->GetBinContent(ibin));
        double errorPosUU = fabs(hPosUU->GetBinContent(ibin)-hPos->GetBinContent(ibin));
        double errorPosDD = fabs(hPosDD->GetBinContent(ibin)-hPos->GetBinContent(ibin));
        errorPos = findMaxVariation(0.0,0.0,errorPosUU,errorPosDD);
        hPos->SetBinError(ibin,errorPos);

        double errorNeg = 0.0;
        double errorNegUD = fabs(hNegUD->GetBinContent(ibin)-hNeg->GetBinContent(ibin));
        double errorNegDU = fabs(hNegDU->GetBinContent(ibin)-hNeg->GetBinContent(ibin));
        double errorNegUU = fabs(hNegUU->GetBinContent(ibin)-hNeg->GetBinContent(ibin));
        double errorNegDD = fabs(hNegDD->GetBinContent(ibin)-hNeg->GetBinContent(ibin));
        errorNeg = findMaxVariation(0.0,0.0,errorNegUU,errorNegDD);
        hNeg->SetBinError(ibin,errorNeg);

        double errorAsymm = 0.0;
        double errorAsymmUD = fabs(hAsymmUD->GetBinContent(ibin)-hAsymm->GetBinContent(ibin));
        double errorAsymmDU = fabs(hAsymmDU->GetBinContent(ibin)-hAsymm->GetBinContent(ibin));
        double errorAsymmUU = fabs(hAsymmUU->GetBinContent(ibin)-hAsymm->GetBinContent(ibin));
        double errorAsymmDD = fabs(hAsymmDD->GetBinContent(ibin)-hAsymm->GetBinContent(ibin));
        errorAsymm = findMaxVariation(errorAsymmUD,errorAsymmDU,errorAsymmUU,errorAsymmDD);
        hAsymm->SetBinError(ibin,errorAsymm);
    }//ibin

    _fOut->cd();
    //hPos->Write();
    //hNeg->Write();
    //hAsymm->Write();
    
    _fOut->Write();
}
