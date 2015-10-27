#include "TH1.h"
#include "TVector3.h"
#include "TUnixSystem.h"
#include "TStopwatch.h"
#include "TH2.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"

#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>



void smearVar(const float nominalVar, const float smear, float& smearedVar){

    smearedVar = nominalVar+smear;

}

int getCentralityBin(float centrality){

    std::vector<double> centralityBins;
    centralityBins.push_back(0.0);
    centralityBins.push_back(0.05);
    centralityBins.push_back(0.10);
    centralityBins.push_back(0.15);
    centralityBins.push_back(0.20);
    centralityBins.push_back(0.40);
    centralityBins.push_back(0.80);
    
    const int nCentralityBins = centralityBins.size()-1;
    for(int icent=0; icent<nCentralityBins; ++icent){
        if(centrality>centralityBins[icent]&&centrality<centralityBins[icent+1]) return icent;
    }
    return nCentralityBins;
}

int getEtaBin(float eta){

    std::vector<double> etaBins;
    etaBins.push_back(0.1);
    etaBins.push_back(2.4);
    const int nEtaBins = etaBins.size()-1;

    for(int ieta=0; ieta<nEtaBins; ++ieta){ 
        if(eta>etaBins[ieta]&&eta<etaBins[ieta]) return ieta;
     }
     return nEtaBins;
}

void addModifiedMptToTree(){
    ///Start timer
    TStopwatch timer;
    timer.Start();

    TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
    TFile *fSyst = new TFile("mptHistosForSystematics.07.19.2013.root","read");
//    TString sFileIn = baseString+"HardProbesFiles/HISingleMuonHardProbesData.07.19.2013.root";
//    TString sFileIn = baseString+"MonteCarloFiles/QCD/HISingleMuonJ1_1muPYTHIADataOverlay.07.19.2013.root" ;
//    TString sFileIn = baseString+ "MonteCarloFiles/QCD/HISingleMuonJ2_1muPYTHIADataOverlay.07.19.2013.root";
//    TString sFileIn = baseString+ "MonteCarloFiles/QCD/HISingleMuonJ3_1muPYTHIADataOverlay.07.19.2013.root";
//    TString sFileIn = baseString+ "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_nn.07.19.2013.root";
//    TString sFileIn = baseString+ "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_np.07.19.2013.root";
//    TString sFileIn = baseString+ "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_pn.07.19.2013.root";
//    TString sFileIn = baseString+ "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_pp.07.19.2013.root";
//    TString sFileIn = baseString+ "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.07.19.2013.root";
//    TString sFileIn = baseString+ "MonteCarloFiles/Zmumu/HISingleMuonZmumuPYTHIAHIJINGOverlay.07.19.2013.root";

    std::cout << "Updating file : " << sFileIn << std::endl;
    TFile *fIn = new TFile(sFileIn,"update");
    if(fIn!=0) std::cout << "Updating input file." << std::endl;
    else exit(0);

    TTree* tree = (TTree*)fIn->Get("tree");

    float centralityNt;
    float mpt,mptPhi,mptMod,mptPhiMod,mtMod[50],pt[50],phi[50],eta[50];
    int nmu;

    ///Add modified nu_phi and nu_pt to the existing tree
    //TBranch* bMptMod = tree->Branch("nu_ptMod",&mptMod,"nu_ptMod/F");
    //TBranch* bMptPhiMod = tree->Branch("nu_phiMod",&mptPhiMod,"nu_phiMod/F");
    TBranch *b_nu_ptMod = tree->Branch("nu_ptMod",&mptMod,"nu_ptMod/F");
    TBranch *b_nu_phiMod = tree->Branch("nu_phiMod",&mptPhiMod,"nu_phiMod/F");
    TBranch *b_mtMod = tree->Branch("mtMod",&mtMod,"mtMod[mu_muid_n]/F");

    tree->SetBranchAddress("centrality", &centralityNt);
    tree->SetBranchAddress("nu_phi", &mptPhi);
    tree->SetBranchAddress("nu_pt", &mpt);
    tree->SetBranchAddress("pt", &pt);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("eta", &eta);
    tree->SetBranchAddress("mu_muid_n", &nmu);
    tree->SetBranchStatus("*",1) ;
/*    tree->SetBranchStatus("centrality", 1);
    tree->SetBranchStatus("nu_pt", 1);
    tree->SetBranchStatus("nu_phi", 1);
    tree->SetBranchStatus("nu_ptMod", 1);
    tree->SetBranchStatus("nu_phiMod", 1);
*/
    std::vector<double> centralityBins;
    centralityBins.push_back(0.0);
    centralityBins.push_back(0.05);
    centralityBins.push_back(0.10);
    centralityBins.push_back(0.15);
    centralityBins.push_back(0.20);
    centralityBins.push_back(0.40);
    centralityBins.push_back(0.80);
    
    const int nCentralityBins = centralityBins.size()-1;

    const int nEntries = tree->GetEntries();
    std::cout << "Number of entries: " << nEntries << std::endl;
    for(int iev=0; iev<nEntries; ++iev){

        //tree->LoadTree(iev);
        tree->GetEntry(iev);

        ///Get bin of current centrality class
        int icent = getCentralityBin(centralityNt);
        int ieta = 0;

        Double_t sampledDelPhi=0.0,sampledDelPt=0.0;

        ///retrieve 2D histo of DelPhi vs. DelPt of mpt wrt truth neutrino
        ///corresponding to this centrality class
        TString sHistName ="h2DPtPhiNeutrino_Eta"; sHistName+=ieta; sHistName+="_Cent"; sHistName+=icent;  
        TH2* hSyst = NULL;
        ///sample the 2D distribution
        if(icent<nCentralityBins) {
            hSyst = (TH2*)fSyst->Get(sHistName);
            if(hSyst==0) {std::cout << "ERROR: Cannot locate input histogram." << std::endl; exit(0);}
            hSyst->GetRandom2(sampledDelPt,sampledDelPhi);
        }
        

        ///Swing phi and modulate magnitude of mpt vector
        ///based on sampled working point in DelPhi,DelPt space
        smearVar(mptPhi,sampledDelPhi,mptPhiMod);
        smearVar(mpt,sampledDelPt,mptMod);

        ///Re-sample if mpt is negative (OK for phi)
        while(mptMod<0.0&&mptMod>-9000.&&icent<nCentralityBins) {
            if(hSyst==0) {std::cout << "ERROR: Cannot locate input histogram." << std::endl; exit(0);}
            hSyst->GetRandom2(sampledDelPt,sampledDelPhi);
            //std::cout << "sampledDelPhi: " << sampledDelPhi << " sampledDelPt: " << sampledDelPt << std::endl;
            smearVar(mpt,sampledDelPt,mptMod);
            //std::cout << "Nominal mpt: " << mpt << " smeared mpt: " << mptMod << std::endl;
        }

        for(int imu=0; imu<nmu; ++imu){
           //TVector3 vMu; vMu.SetPtEtaPhi(pt[imu],eta[imu],phi[imu]);
           //float phiMu = vMu.Phi();
           float dPhi = mptPhiMod - phi[imu];
           //float dPhi = mptPhiMod - phiMu;
           //std::cout << "Consistency check: " << dPhi1 << "=?" << dPhi << std::endl;
           if(dPhi >TMath::Pi() ) {dPhi -= TMath::TwoPi();} if(dPhi<-1.*TMath::Pi()){dPhi+=TMath::TwoPi();} //fold between -Pi and +Pi
           mtMod[imu] = (fabs(mpt) < 9000.) ? TMath::Sqrt(2.*fabs(pt[imu])*mptMod*(1.-TMath::Cos(dPhi))) : -9999. ;
        }
        if(iev%10000==0) {
            std::cout << "Event: " << iev << std::endl;
            if(sampledDelPhi!=0.0&&sampledDelPt!=0.0) {
            //if(icent==nCentralityBins) {
              std::cout << "sampledDelPhi: " << sampledDelPhi << " sampledDelPt: " << sampledDelPt << std::endl;
              std::cout << "Nominal mpt: " << mpt << " smeared mpt: " << mptMod << std::endl;
              std::cout << "Nominal phi of mpt vector: " << mptPhi << " smeared phi: " << mptPhiMod << std::endl;
             }
        }

        ///Fill new branches in pre-existing tree
        b_nu_ptMod->Fill();
        b_nu_phiMod->Fill();
        b_mtMod->Fill();

        /*if(sampledDelPhi!=0.0&&sampledDelPt!=0.0) {
            std::cout << "sampledDelPhi: " << sampledDelPhi << " sampledDelPt: " << sampledDelPt << std::endl;
            std::cout << "Nominal mpt: " << mpt << " smeared mpt: " << mptMod << std::endl;
            std::cout << "Nominal phi of mpt vector: " << mptPhi << " smeared phi: " << mptPhiMod << std::endl;
        }*/

//        tree->Fill();

    }//iev

//    tree->Print();
    tree->Write();
    delete fIn;
    delete fSyst;

    // stop timer and print results
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();
    printf("\nRealTime=%f seconds, CpuTime=%f seconds\n", rtime, ctime);
}
