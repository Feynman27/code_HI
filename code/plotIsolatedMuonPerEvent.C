#include "TH1I.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"

#include <iostream>
#include <cmath>

void plotIsolatedMuonPerEvent(){

      TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
      TString fileNameIn = baseString+"HardProbesFiles/HISingleMuonHardProbesData.04.17.2013.root";  

      float eLossNt[50];
      float scatNt[50];
      float compNt[50];
      float ptNt[50];
      float mtNt[50];
      float etaNt[50];
      float phiNt[50];
      float chargeNt[50];
      int promptNt[50];
      float centralityNt;
      float nu_ptNt;
      float ptconeNt[50];
      int valNt[50], truthMatchedNt[50],ZDYNt[50], matched1[50], matched2[50], matched3[50];
      int nmu,trig1,trig2,trig3,trig4,trig5;


      TChain* tree = new TChain("tree","tree");
      tree->Add(fileNameIn);

      std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;


     
      tree->SetBranchAddress("pt",&ptNt);
      tree->SetBranchAddress("ptcone20ID3",&ptconeNt);
      tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC",&trig1);
      tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10",&trig2);
      tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20",&trig3);
      tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC_Matched20",&matched1);
      tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10_Matched20",&matched2);
      tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20_Matched20",&matched3);

      tree->SetBranchAddress("eLoss", &eLossNt);
      tree->SetBranchAddress("ptcone20ID3", &ptconeNt);
      tree->SetBranchAddress("scat", &scatNt);
      tree->SetBranchAddress("comp", &compNt);
      tree->SetBranchAddress("pt", &ptNt);
      tree->SetBranchAddress("mt", &mtNt);
      tree->SetBranchAddress("eta", &etaNt);
      tree->SetBranchAddress("phi", &phiNt);
      tree->SetBranchAddress("charge", &chargeNt);
      tree->SetBranchAddress("val", &valNt); 
      tree->SetBranchAddress("ZDY", &ZDYNt); 
      tree->SetBranchAddress("centrality", &centralityNt);
      tree->SetBranchAddress("nu_pt", &nu_ptNt);
      tree->SetBranchAddress("mu_muid_n", &nmu);

       // --- Set branch status ---
      tree->SetBranchStatus("*",0) ;
      tree->SetBranchStatus("mu_muid_n", 1);
      tree->SetBranchStatus("eLoss", 1);
      tree->SetBranchStatus("ptcone20ID3", 1);
      tree->SetBranchStatus("scat", 1);
      tree->SetBranchStatus("comp", 1);
      tree->SetBranchStatus("pt", 1);
      tree->SetBranchStatus("mt", 1);
      tree->SetBranchStatus("eta", 1);
      tree->SetBranchStatus("phi", 1);
      tree->SetBranchStatus("charge", 1);
      tree->SetBranchStatus("val", 1); 
      tree->SetBranchStatus("ZDY", 1); 
      tree->SetBranchStatus("centrality", 1);
      tree->SetBranchStatus("nu_pt", 1);
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC",1);
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10",1);
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20",1);
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC_Matched20",1);
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10_Matched20",1);
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20_Matched20",1);

    TH1I* hIsolatedMuons = new TH1I("hIsolatedMuons","hIsolatedMuons",80,0,4);
    TH1I* hWCandidates = new TH1I("hWCandidates","hWCandidates",80,0,4);
    for(int iev = 0; iev < tree->GetEntries(); ++iev){
 
      tree->LoadTree(iev);
      tree->GetEntry(iev);
      
      int nIsolatedMuons = 0;
      int nWCandidates = 0;
      int hiMuonPtFlag = 0;
      for(int imu=0; imu<nmu; ++imu){
            
        if(ptNt[imu]>10.0){
            hiMuonPtFlag = 1;
        if(
            valNt[imu]>11 
            && abs(scatNt[imu])<4.0
            && abs(eLossNt[imu])<0.5
            && abs(etaNt[imu])>0.1 
            && abs(etaNt[imu])<2.4 
            && centralityNt>=0. 
            && centralityNt<=0.8
            && ( ( (trig1&&matched1[imu]) ||
                (trig2&&matched2[imu])
                ||(trig3&&matched3[imu]) ) )
            &&ptconeNt[imu]/ptNt[imu]<0.1
            && ZDYNt[imu]==0
            ){
                ++nIsolatedMuons;

                ///W Candidates
                if(
                    ptNt[imu]>25.0
                    && nu_ptNt >25.0
                    //&& ZDYNt[imu]==0
                    && mtNt[imu]>40.0
                  ){

                    ++nWCandidates;
                }
        }
      }
     } //imu
      
      ///Fill with number of events with n isolated/W candidates
      ///and go to next event; fill only for high pT muons (>10GeV)
      if(hiMuonPtFlag==1) 
      {
         hIsolatedMuons->Fill(nIsolatedMuons);
         hWCandidates->Fill(nWCandidates);
      }
    } //iev

    TCanvas *cIsol = new TCanvas("cIsol","cIsol",600,600);
    hIsolatedMuons->Draw("hist");
    TCanvas *cW = new TCanvas("cW","cW",600,600);
    hWCandidates->Draw("hist");
}
