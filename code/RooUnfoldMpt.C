#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "/mnt/Lustre/cgrp/atlas_hi/tbalestri/RooUnfold-1.1.1/src/RooUnfold.h"
#include "/mnt/Lustre/cgrp/atlas_hi/tbalestri/RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "/mnt/Lustre/cgrp/atlas_hi/tbalestri/RooUnfold-1.1.1/src/RooUnfoldResponse.h"

#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TChain.h"
#include "TH1F.h"
#include "TMath.h"
#include "TLatex.h"
#include "TTree.h"
#include "TColor.h"
#include "TLegendEntry.h"
#include "TPad.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath> 
#include <sstream>

void RooUnfoldMpt(){
#ifdef __CINT__
    gSystem->Load("/mnt/Lustre/cgrp/atlas_hi/tbalestri/RooUnfold-1.1.1/libRooUnfold");
#endif
    double xLo = 0.0, xUpp = 120.0;
    int nBins = 24; 

    ///supply the response matrix object RooUnfoldResponse
    ///Same binning as used in final plots
    RooUnfoldResponse response(nBins,xLo,xUpp);
    TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
    TString fileNameIn = baseString+"/MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.04.12.2013.root";  

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
      float nuGenPtNt[50];
      int motherNt[50], ngen;

      TChain* tree = new TChain("tree","tree");
      tree->Add(fileNameIn);

      std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;

      tree->SetBranchAddress("mc_mu_n", &ngen);
      tree->SetBranchAddress("mc_nu_gen_pt", &nuGenPtNt);
      tree->SetBranchAddress("mc_mu_gen_mothertype", &motherNt);
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
      tree->SetBranchAddress("prompt", &promptNt);

       // --- Set branch status ---
      tree->SetBranchStatus("*",0) ;
      tree->SetBranchStatus("mc_mu_n", 1);
      tree->SetBranchStatus("mc_nu_gen_pt", 1);
      tree->SetBranchStatus("mc_mu_gen_mothertype", 1);
      tree->SetBranchStatus("mu_muid_n", 1);
      tree->SetBranchStatus("prompt", 1);
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
cout << "==================================== TRAIN ====================================" << endl;
    for(int iev = 0; iev < tree->GetEntries(); ++iev){
 
      tree->LoadTree(iev);
      tree->GetEntry(iev);

      for(int igen=0; igen<ngen; ++igen){
        for(int imu=0; imu<nmu; ++imu){

          if (
                (valNt[imu]>11) 
                && abs(scatNt[imu])<4.0
                && abs(eLossNt[imu])<0.5
                && abs(etaNt[imu])>0.1 
                && abs(etaNt[imu])<2.4 
                && ptNt[imu]>=25.0
                && centralityNt>=0. 
                && centralityNt<=0.8
                && ZDYNt[imu]==0
                && ( ptconeNt[imu]/ptNt[imu] < 0.1)
                &&mtNt[imu]>40.
                && ( promptNt[imu]==24)
		        ){
                    response.Fill(nu_ptNt,nuGenPtNt[igen]);
                }
            else response.Miss(nuGenPtNt[igen]);
         } //imu

      } //igen
   }
cout << "==================================== TEST =====================================" << endl;
   TH1D* hTrue= new TH1D("hTrue", "Test Truth", nBins, xLo, xUpp);
   TH1D* hMeas= new TH1D("hMeas", "Test Measured", nBins, xLo, xUpp);
    for(int iev = 0; iev < tree->GetEntries(); ++iev){
 
      tree->LoadTree(iev);
      tree->GetEntry(iev);

      for(int igen=0; igen<ngen; ++igen){
        ///Fill truth neutrino
        if(fabs(motherNt[igen])==24) {
            hTrue->Fill(nuGenPtNt[igen]);
        }
      }

      for(int imu=0; imu<nmu; ++imu){
            

          if (
                (valNt[imu]>11) 
                && abs(scatNt[imu])<4.0
                && abs(eLossNt[imu])<0.5
                && abs(etaNt[imu])>0.1 
                && abs(etaNt[imu])<2.4 
                && ptNt[imu]>=25.0
                && centralityNt>=0. 
                && centralityNt<=0.8
                && ZDYNt[imu]==0
                && ( ptconeNt[imu]/ptNt[imu] < 0.1)
                &&mtNt[imu]>40.
                && ( promptNt[imu]==24)
		        ){
                    hMeas->Fill(nu_ptNt);
                }
      }
    } //iev

    std::cout << "==================================== UNFOLD ===================================" << std::endl;
    RooUnfoldBayes unfold(&response, hMeas, 100);
    TH1D* hReco= (TH1D*) unfold.Hreco();

    unfold.PrintTable (std::cout, hTrue);
    hReco->Draw();
    hMeas->Draw("SAME");
    hTrue->SetLineColor(8);
    hTrue->Draw("SAME");
    
}
