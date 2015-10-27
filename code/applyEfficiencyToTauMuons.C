#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include <iostream>

void applyEfficiencyToTauMuons(){

    TFile* _fPtSm = new TFile("smearedMuonPtFromTau.07.09.2013.root","read");
    TFile* _fEff = new TFile("mcWEffForTauMuonStudy_07_08_2013.root","read");
    TFile* _fPtSmEff = new TFile("smearedEffAppliedMuonPtFromTau.root","recreate");

    _fPtSm->cd();
    TH1F* hPtSm = (TH1F*)_fPtSm->Get("hSmearedMuonPt");
    _fEff->cd();
    TGraphAsymmErrors* grEff = (TGraphAsymmErrors*)_fEff->Get("pEffWAccCuts");

    double* yEff = grEff->GetY();  
    TH1F*  hnew = (TH1F*)hPtSm->Clone("hnew"); 
    ///Apply efficiency bin by bin
    for(int igr=0; igr<grEff->GetN(); ++igr){
        
        double ptNominal = hPtSm->GetBinContent(igr+1);
        std::cout << "pt nominal " << ptNominal << " eff " << yEff[igr] << " = " << yEff[igr]*ptNominal << std::endl;
        hnew->SetBinContent(igr+1,yEff[igr]*ptNominal);
    }

    _fPtSmEff->cd();
    hnew->Write();
    _fPtSmEff->Write();
    _fPtSmEff->Close();
}
