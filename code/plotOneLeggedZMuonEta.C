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


double foldZeroToPi(double& dPhi){

    ///Fold between -pi to pi

    if(dPhi>TMath::Pi()) dPhi -= TMath::TwoPi();
    if(dPhi<-1.*TMath::Pi()) dPhi += TMath::TwoPi();

    ///If -pi/2 or -pi, add pi
    if(dPhi==-1.0*TMath::PiOver2()||dPhi==-1.0*TMath::Pi()){

        dPhi+=TMath::Pi();
    }
    ///Quadrant IV
    else if(dPhi>-1.0*TMath::PiOver2()&&dPhi<0.0){
        dPhi+=TMath::PiOver2();
    }
    ///Quadrant III
    else if(dPhi>-1.0*TMath::Pi()&&dPhi<-1.0*TMath::PiOver2()){

        dPhi+=3.0*TMath::PiOver2();
    }
    if(dPhi<0.0||dPhi>TMath::Pi()){

        std::cout << "WARNING: Delta Phi outside [0,pi]. Expect wrong results." << std::endl;
        return -9999.0;
    }
    else return dPhi;
}

void plotOneLeggedZMuonEta(){

      TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
	  TString fileNameIn = baseString+ "MonteCarloFiles/Zmumu/HISingleMuonZmumuPYTHIAHIJINGOverlay.04.13.2013.root";

      ///reco level variables
      float eLossNt[50];
      float scatNt[50];
      float compNt[50];
      float ptNt[50];
      float mtNt[50];
      float etaNt[50];
      float phiNt[50];
      float mptPhi;
      float chargeNt[50];
      int promptNt[50];
      float centralityNt;
      float nu_ptNt;
      float ptconeNt[50];
      int valNt[50], truthMatchedNt[50],ZDYNt[50], matched1[50], matched2[50], matched3[50];
      int nmu,trig1,trig2,trig3,trig4,trig5;

      ///generator level variables
      float phiGen[50],etaGen[50],chargeGen[50];
      int mother[50],daughter[50],ngen;

      TChain* tree = new TChain("tree","tree");
      tree->Add(fileNameIn);

      std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;
        
      tree->SetBranchAddress("mc_mu_gen_mothertype", &mother);
      tree->SetBranchAddress("mc_mu_gen_type", &daughter);
      tree->SetBranchAddress("mc_mu_charge", &chargeGen);
      tree->SetBranchAddress("mc_mu_gen_eta", &etaGen);
      tree->SetBranchAddress("mc_mu_gen_phi", &phiGen);
      tree->SetBranchAddress("mc_mu_n", &ngen);


      tree->SetBranchAddress("pt",&ptNt);
      tree->SetBranchAddress("prompt",&promptNt);
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
      tree->SetBranchAddress("nu_phi", &mptPhi);
      tree->SetBranchAddress("charge", &chargeNt);
      tree->SetBranchAddress("val", &valNt); 
      tree->SetBranchAddress("ZDY", &ZDYNt); 
      tree->SetBranchAddress("centrality", &centralityNt);
      tree->SetBranchAddress("nu_pt", &nu_ptNt);
      tree->SetBranchAddress("mu_muid_n", &nmu);

       // --- Set branch status ---
      tree->SetBranchStatus("*",0) ;

      tree->SetBranchStatus("mc_mu_n", 1);
      tree->SetBranchStatus("mc_mu_gen_mothertype", 1);
      tree->SetBranchStatus("mc_mu_gen_type", 1);
      tree->SetBranchStatus("mc_mu_charge", 1);
      tree->SetBranchStatus("mc_mu_gen_eta", 1);
      tree->SetBranchStatus("mc_mu_gen_phi", 1);
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

   TH1F* hEtaLost = new TH1F("hEtaLost","hEtaLost",50,-5.0,5.0);
   TH1F* hDeltaPhiMptLost = new TH1F("hDeltaPhiMptLost","hDeltaPhiMptLost",50,0.0,3.2);

   for(int iev = 0; iev < tree->GetEntries(); ++iev){
 
    tree->LoadTree(iev);
    tree->GetEntry(iev);

    bool matchedToRec = false;

    ///Loop over generated muons in the sample
    for(int igen=0; igen<ngen; ++igen){

      ///if we've already matched a generated muon
      ///from a Z to a reconstructed muon has passed
      ///W selection, quit the generated loop and 
      ///move on to the next event
      if(matchedToRec) break;

      ///Make sure the generated muon is from a Z
      if(!(fabs(mother[igen])==23&&fabs(daughter[igen])==13)) continue;

      ///Loop over all reco muons and find 
      ///a match to the one that passes W selection
      for(int imu=0; imu<nmu; ++imu){
            
            ///Calculte eta,phi dR of gen and rec muon
            double dPhi = phiGen[igen] - phiNt[imu];
            if(dPhi>TMath::Pi()) dPhi -= TMath::TwoPi();
            if(dPhi<-1.*TMath::Pi()) dPhi += TMath::TwoPi();

            double dEta = etaGen[igen] - etaNt[imu]; 

            double dR = TMath::Sqrt( TMath::Power(dPhi,2)+TMath::Power(dEta,2) );

            ///Match generated muon to rec muon passing 
            ///W selection with dR<0.1
            if(
                promptNt[imu] == 23
                && dR < 0.1
                && valNt[imu]>11
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
                &&ptNt[imu]>25.0
                &&nu_ptNt >25.0
                &&mtNt[imu]>40.0
            ){
                ///once we've found a match, break
                ///from the rec loop and find the other
                ///generated muon
                matchedToRec = true;
                break;
            }
      } //imu

    ///Now that we've matched a generated muon from a Z with a 
    ///reconstructed muon from a Z, let's find the "lost muon"
    for(int ilost=0; ilost<ngen; ++ilost){

        ///make sure not to use the same generated muon
        if(
            !(ilost==igen)
            && matchedToRec == true 
            && fabs(mother[ilost])==23
            && fabs(daughter[ilost])==13
            && (chargeGen[igen]*chargeGen[ilost]) < 0
        )
        {
          ///Delta phi of lost muon and mpt vector
          float dPhiMptZ = phiGen[ilost] - mptPhi;
          ///Fold in [0.0,PI]
          //dPhiMptZ = foldZeroToPi(dPhiMptZ);
          ///Fill histogram of eta distribution for lost muon
          hEtaLost->Fill(etaGen[ilost]);
          hDeltaPhiMptLost->Fill(dPhiMptZ);
          break;
        }
    }//ilost

   }//igen

  } //iev

  TCanvas *cLost = new TCanvas("cLost","cLost",600,600);
  hEtaLost->Draw("HIST");
  hEtaLost->GetXaxis()->SetTitle("#eta_{Lost #mu}");
  hEtaLost->GetYaxis()->SetTitle("Muons");
  TLatex *   tex = new TLatex(-0.9085889,2930.899,"Z#rightarrow#mu#mu");
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(-1.524509,2608.524,"Monte Carlo");
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->Draw();

  TCanvas *cDelPhi = new TCanvas("cDelPhi","cDelPhi",600,600);
  hDeltaPhiMptLost->Draw("HIST");

  hDeltaPhiMptLost->GetXaxis()->SetTitle("#Delta#phi(Lost #mu-#slash{p_{T}})");
  hDeltaPhiMptLost->GetYaxis()->SetTitle("Muons");
  tex = new TLatex(1.506346,1770.478,"Z#rightarrow#mu#mu");
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(1.309252,1485.197,"Monte Carlo");
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->Draw();

}
