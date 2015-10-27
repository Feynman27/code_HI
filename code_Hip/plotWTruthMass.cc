#include<iostream>
#include "TString.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TFile.h"

void fillMassHisto(TString sFileName,TH1F* h){
    
    TChain* _tGen = new TChain("treeGen","treeGen");
    _tGen->Add(sFileName);
    std::cout << "Number of entries in generator tree: " << _tGen->GetEntries() << std::endl;

    float mass;
    _tGen->SetBranchAddress("W_mass",&mass);
    for(int i=0;i<_tGen->GetEntries();++i){
        _tGen->GetEntry(i);
        h->Fill(mass);
    }//i

    delete _tGen;
}

int main(){

    TString sFileName = "/usatlas/u/tbales/scratch/Data2013Hip/";
    TFile* _fOut = new TFile("W_truth_mass.root","recreate");

    const int nBins = 140;
    const int massLow = 50.0;
    const int massUpp = 120.0;

    TH1F* hplus_pp_A = new TH1F("hplus_pp_A","hplus_pp_A",nBins,massLow,massUpp);
    TH1F* hplus_pp_B = new TH1F("hplus_pp_B","hplus_pp_B",nBins,massLow,massUpp);
    TH1F* hplus_pn_A = new TH1F("hplus_pn_A","hplus_pn_A",nBins,massLow,massUpp);
    TH1F* hplus_pn_B = new TH1F("hplus_pn_B","hplus_pn_B",nBins,massLow,massUpp);
    TH1F* hplus_wtd_A = new TH1F("hplus_wtd_A","hplus_wtd_A",nBins,massLow,massUpp);
    TH1F* hplus_wtd_B = new TH1F("hplus_wtd_B","hplus_wtd_B",nBins,massLow,massUpp);

    TH1F* hminus_pp_A = new TH1F("hminus_pp_A","hminus_pp_A",nBins,massLow,massUpp);
    TH1F* hminus_pp_B = new TH1F("hminus_pp_B","hminus_pp_B",nBins,massLow,massUpp);
    TH1F* hminus_pn_A = new TH1F("hminus_pn_A","hminus_pn_A",nBins,massLow,massUpp);
    TH1F* hminus_pn_B = new TH1F("hminus_pn_B","hminus_pn_B",nBins,massLow,massUpp);
    TH1F* hminus_wtd_A = new TH1F("hminus_wtd_A","hminus_wtd_A",nBins,massLow,massUpp);
    TH1F* hminus_wtd_B = new TH1F("hminus_wtd_B","hminus_wtd_B",nBins,massLow,massUpp);

    double wt_pp = 82.0/208;
    double wt_np = (208.0-82.0)/208.0;

    fillMassHisto(sFileName+"ntupleWplusmunu_pp.A.PowhegPythia8_AU2CT10.09.24.2014.root",hplus_pp_A);
    fillMassHisto(sFileName+"ntupleWplusmunu_pp.B.PowhegPythia8_AU2CT10.09.24.2014.root",hplus_pp_B);
    fillMassHisto(sFileName+"ntupleWplusmunu_pn.A.PowhegPythia8_AU2CT10.09.24.2014.root",hplus_pn_A);
    fillMassHisto(sFileName+"ntupleWplusmunu_pn.B.PowhegPythia8_AU2CT10.09.24.2014.root",hplus_pn_B);
    hplus_wtd_A->Add(hplus_pp_A,hplus_pn_A,wt_pp,wt_np);
    hplus_wtd_B->Add(hplus_pp_B,hplus_pn_B,wt_pp,wt_np);

    fillMassHisto(sFileName+"ntupleWminmunu_pp.A.PowhegPythia8_AU2CT10.09.24.2014.root",hminus_pp_A);
    fillMassHisto(sFileName+"ntupleWminmunu_pp.B.PowhegPythia8_AU2CT10.09.24.2014.root",hminus_pp_B);
    fillMassHisto(sFileName+"ntupleWminmunu_pn.A.PowhegPythia8_AU2CT10.09.24.2014.root",hminus_pn_A);
    fillMassHisto(sFileName+"ntupleWminmunu_pn.B.PowhegPythia8_AU2CT10.09.24.2014.root",hminus_pn_B);
    hminus_wtd_A->Add(hminus_pp_A,hminus_pn_A,wt_pp,wt_np);
    hminus_wtd_B->Add(hminus_pp_B,hminus_pn_B,wt_pp,wt_np);


    _fOut->cd();
    hplus_pp_A->Write() ;
    hplus_pp_B->Write();
    hplus_pn_A->Write();
    hplus_pn_B->Write();
    hplus_wtd_A->Write(); 
    hplus_wtd_B->Write();
            
    hminus_pp_A->Write();
    hminus_pp_B->Write();
    hminus_pn_A->Write();
    hminus_pn_B->Write();
    hminus_wtd_A->Write();
    hminus_wtd_B->Write();

    delete      hplus_pp_A ;
    delete      hplus_pp_B ;
    delete      hplus_pn_A ;
    delete      hplus_pn_B ;
    delete      hplus_wtd_A ; 
    delete      hplus_wtd_B ;
                       
    delete      hminus_pp_A ;
    delete      hminus_pp_B ;
    delete      hminus_pn_A ;
    delete      hminus_pn_B ;
    delete      hminus_wtd_A;
    delete      hminus_wtd_B;

    delete _fOut;
}
