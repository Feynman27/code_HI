#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLatex.h"

#include <iostream>
#include <cmath>
#include <vector>



void fillDelPhi(TString fileNameIn,TH1F* hDelPhi,double centralityLo,double centralityUpp,double etaLo,double etaUpp){

      ///reco level variables
      float phiNt[50],etaNt[50];
      float mptPhi;
      int promptNt[50];
      float centralityNt;
      float mpt;
      int nmu;

      ///generator level variables
      int mother[50],daughter[50],ngen;
      float neutrinoPhi[50],neutrinoPt[50];

      TChain* tree = new TChain("tree","tree");
      tree->Add(fileNameIn);

      std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;
        
      tree->SetBranchAddress("mc_nu_gen_mothertype", &mother);
      tree->SetBranchAddress("mc_nu_gen_phi", &neutrinoPhi);
      tree->SetBranchAddress("mc_nu_gen_type", &daughter);
      tree->SetBranchAddress("mc_nu_gen_pt", &neutrinoPt);
      tree->SetBranchAddress("mc_mu_n", &ngen);

      tree->SetBranchAddress("centrality", &centralityNt);
      tree->SetBranchAddress("nu_phi", &mptPhi);
      tree->SetBranchAddress("nu_pt", &mpt);
      tree->SetBranchAddress("mu_muid_n", &nmu);
      tree->SetBranchAddress("prompt",&promptNt);
      tree->SetBranchAddress("eta",&etaNt);

      // --- Set branch status ---
      tree->SetBranchStatus("*",0) ;
      tree->SetBranchStatus("mc_mu_n", 1);
      tree->SetBranchStatus("mc_nu_gen_mothertype", 1);
      tree->SetBranchStatus("mc_nu_gen_type", 1);
      tree->SetBranchStatus("mc_nu_gen_phi", 1);
      tree->SetBranchStatus("mc_nu_gen_pt", 1);
      tree->SetBranchStatus("mu_muid_n", 1);
      tree->SetBranchStatus("eLoss", 1);
      tree->SetBranchStatus("ptcone20ID3", 1);
      tree->SetBranchStatus("scat", 1);
      tree->SetBranchStatus("comp", 1);
      tree->SetBranchStatus("prompt",1);
      tree->SetBranchStatus("pt", 1);
      tree->SetBranchStatus("mt", 1);
      tree->SetBranchStatus("eta", 1);
      tree->SetBranchStatus("phi", 1);
      tree->SetBranchStatus("nu_phi", 1);
      tree->SetBranchStatus("val", 1); 
      tree->SetBranchStatus("ZDY", 1); 
      tree->SetBranchStatus("centrality", 1);
      tree->SetBranchStatus("nu_pt", 1);
   
   for(int iev = 0; iev < tree->GetEntries(); ++iev){
 
    tree->LoadTree(iev);
    tree->GetEntry(iev);
    if(!(centralityNt>centralityLo&&centralityNt<centralityUpp)) continue;

    ///Loop over generated muons in the sample
    ///and fill histo with dPhi of reco mpt and
    ///neutrino pt
    for(int igen=0; igen<ngen; ++igen){
 
        ///consider only neutrinos from Ws
       if(!(fabs(mother[igen])==24&&fabs(daughter[igen])==14)) continue;
         {

            for(int imu=0; imu<nmu; ++imu){
               ///consider only reconstructed muons from Ws
               if(!(promptNt[imu]==24&&fabs(etaNt[imu])>etaLo&&fabs(etaNt[imu])<etaUpp)) continue; 
                 ///Calculate dPhi btwn mpt phi and neutrino phi
                 ///if this event has a reco muon from a W
                 float dPhi = mptPhi-neutrinoPhi[igen];
                 ///Fold dPhi in[-PI,PI]
                 if(dPhi>TMath::Pi()) dPhi-=2.*TMath::Pi();
                 else if(dPhi<-1.0*TMath::Pi()) dPhi+=2.*TMath::Pi();
                 hDelPhi->Fill(dPhi);
             }//imu 
         }
    }//igen
   } //iev

}


void plotMptNeutrinoDPhi(){

   TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
   TString fileNameIn = baseString+ "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.04.12.2013.root";

   std::vector<double> centralityBins;
   centralityBins.push_back(0.0);
   centralityBins.push_back(0.05);
   centralityBins.push_back(0.10);
   centralityBins.push_back(0.15);
   centralityBins.push_back(0.20);
   centralityBins.push_back(0.40);
   centralityBins.push_back(0.80);

   const int nCentralityBins = centralityBins.size()-1;

   std::vector<double> etaBins;
   etaBins.push_back(0.1);
   etaBins.push_back(2.4);
   const int nEtaBins = etaBins.size()-1;

   std::cout << "Running over " << nCentralityBins << " centrality classes and " << nEtaBins << " eta windows." << std::endl;

   TCanvas *cPhi = new TCanvas("cPhi","cPhi",600,600);
   ///array of delta phi histos
   TH1F* hDelPhi[nEtaBins][nCentralityBins];
   char sName[50];
   for(int ieta=0; ieta<nEtaBins; ++ieta){
    for(int icent=0; icent<nCentralityBins; ++icent){

        sprintf(sName,"dPhiMptNeutrino_Eta%i_Cent%i",ieta,icent);
        hDelPhi[ieta][icent] = new TH1F(sName,sName,50,-3.2,+3.2);
        fillDelPhi(fileNameIn,hDelPhi[ieta][icent],centralityBins[icent],centralityBins[icent+1],etaBins[ieta],etaBins[ieta+1]);
        if(icent==0) {
            hDelPhi[ieta][icent]->Draw("pe");
            hDelPhi[ieta][icent]->GetXaxis()->SetTitle("#Delta#phi(#nu_{#mu}-#slash{p_{T}})");
        }
        else hDelPhi[ieta][icent]->Draw("pesame");
        //delete hDelPhi[ieta][icent];
    }//icent

    TLatex* tex = new TLatex(1.506346,1770.478,"W#rightarrow#mu#nu MC11");
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();

   }//ieta

//   outFile->Close();
}
