#include "TH1F.h"
#include "TVector3.h"
#include "TProfile.h"
#include "TDatime.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLegend.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <string>


int getCentralityBin(float centrality){
    int index = 0;
    if(centrality>0.0&&centrality<0.05) return index; else ++index;
    if(centrality >0.05&&centrality<0.1) return index; else ++index;
    if(centrality >0.1&&centrality<0.15) return index; else ++index;
    if(centrality >0.15&&centrality<0.2) return index; else ++index;
    if(centrality >0.20&&centrality<0.40) return index; else ++index;
    if(centrality >0.40&&centrality<0.80) return index; else ++index;
    return -1;
}

float getMt(float pt, float eta, float phi, float mptPhi, float mpt){

    TVector3 lvMu;
    lvMu.SetPtEtaPhi(pt,eta,phi) ;
    float phiMuTemp = lvMu.Phi();
    if(fabs(phiMuTemp)!=fabs(phi)) {
        std::cout << phiMuTemp << "!=" << phi << std::endl;
        std::cout << "Warning: inconsistency in phi." << std::endl;
        //exit(0);
    }
    float phi_munu = phiMuTemp-mptPhi;
    if(phi_munu> TMath::Pi()) {phi_munu -= TMath::TwoPi();}
    if(phi_munu<-1*TMath::Pi()) {phi_munu+=TMath::TwoPi();}
    float mt = (fabs(mpt) < 9999.) ? TMath::Sqrt(2.*fabs(pt)*fabs(mpt)*(1.0-TMath::Cos(phi_munu))) : -9999. ;
    return mt;
}

void plotMeanMpt234GeV(){

      TDatime* time = new TDatime();
      TString sDate = "";
      sDate+=time->GetMonth(); sDate+="_"; sDate+=time->GetDay(); sDate+="_"; sDate+=time->GetYear();
      std::cout << "Today is " << sDate << std::endl;

      TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
      TString fileNameIn = baseString; 
      fileNameIn+="HardProbesFiles/HISingleMuonHardProbesData.07.13.2013.root";

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
      float mpt2,mpt3,mpt4;
      float mptPhi2,mptPhi4;
      float ptconeNt[50];
      int valNt[50], truthMatchedNt[50],ZDYNt[50], matched1[50], matched2[50], matched3[50];
      int nmu,trig1,trig2,trig3,trig4,trig5;

      TChain* tree = new TChain("tree","tree");
      tree->Add(fileNameIn);
      std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;

      tree->SetBranchAddress("pt",&ptNt);
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
      tree->SetBranchAddress("mu_muid_n", &nmu);
      tree->SetBranchAddress("nu_pt", &mpt3);
      tree->SetBranchAddress("nu_pt2000Nominal", &mpt2);
      tree->SetBranchAddress("nu_pt4000Nominal", &mpt4);
      tree->SetBranchAddress("nu_phi2000Nominal", &mptPhi2);
      tree->SetBranchAddress("nu_phi4000Nominal", &mptPhi4);
      tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC",&trig1);
      tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10",&trig2);
      tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20",&trig3);
      tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC_Matched20",&matched1);
      tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10_Matched20",&matched2);
      tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20_Matched20",&matched3);

      tree->SetBranchStatus("*",0);
      tree->SetBranchStatus("pt",1);
      tree->SetBranchStatus("eLoss", 1);
      tree->SetBranchStatus("ptcone20ID3", 1);
      tree->SetBranchStatus("scat", 1);
      tree->SetBranchStatus("comp", 1);
      tree->SetBranchStatus("pt", 1);
      tree->SetBranchStatus("mt", 1);
      tree->SetBranchStatus("eta", 1);
      tree->SetBranchStatus("phi", 1);
      tree->SetBranchStatus("nu_phi", 1);
      tree->SetBranchStatus("charge", 1);
      tree->SetBranchStatus("val", 1); 
      tree->SetBranchStatus("ZDY", 1); 
      tree->SetBranchStatus("centrality", 1);
      tree->SetBranchStatus("mu_muid_n", 1);
      tree->SetBranchStatus("nu_pt", 1);
      tree->SetBranchStatus("nu_pt2000Nominal", 1);
      tree->SetBranchStatus("nu_pt4000Nominal", 1);
      tree->SetBranchStatus("nu_pt", 1);
      tree->SetBranchStatus("nu_phi2000Nominal", 1);
      tree->SetBranchStatus("nu_phi4000Nominal", 1);
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC",1);
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10",1);
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20",1);
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC_Matched20",1);
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10_Matched20",1);
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20_Matched20",1);



      char hName2[50],hName3[50],hName4[50];
      const int nCentralityBins = 6;
      TString sCent[] = {"0-5%","5-10%","10-15%","15-20%","20-40%","40-80%"};
      float arrCent[] = {0.0,5.0,10.0,15.0,20.0,40.0,80.0};
        
      TH1F* hMpt2[nCentralityBins],*hMpt3[nCentralityBins],*hMpt4[nCentralityBins];
      TH1F* hMean2 = new TH1F("hMeanMpt2GeV","hMeanMpt2GeV",nCentralityBins,arrCent);
      TH1F* hMean3 = new TH1F("hMeanMpt3GeV","hMeanMpt3GeV",nCentralityBins,arrCent);
      TH1F* hMean4 = new TH1F("hMeanMpt4GeV","hMeanMpt4GeV",nCentralityBins,arrCent);

      for(int icent=0; icent<nCentralityBins; ++icent){
        sprintf(hName2,"hMpt2GeVCent%i",icent);
        sprintf(hName3,"hMpt3GeVCent%i",icent);
        sprintf(hName4,"hMpt4GeVCent%i",icent);
        hMpt2[icent] = new TH1F(hName2,hName2,30,0.0,120.0);
        hMpt3[icent] = new TH1F(hName3,hName3,30,0.0,120.0);
        hMpt4[icent] = new TH1F(hName4,hName4,30,0.0,120.0);
      }//icent


      for(int iev = 0; iev < tree->GetEntries(); ++iev){
 
        tree->LoadTree(iev);
        tree->GetEntry(iev);

        //if(iev==10000) break; //hack
        if(iev%10000==0) std::cout << "Event: " << iev << std::endl;

        int icent = getCentralityBin(centralityNt);
        if(icent<0) continue;
        for(int imu=0; imu<nmu; ++imu){
        
            ///Calculate new mt
            float mt2 = getMt(ptNt[imu],etaNt[imu],phiNt[imu],mptPhi2,mpt2);
            float mt4 = getMt(ptNt[imu],etaNt[imu],phiNt[imu],mptPhi4,mpt4);

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
                &&ptNt[imu]>25.0
                &&mpt2>=0.0
                &&mpt2<9000.
                &&mpt3>=0.0
                &&mpt3<9000.
                &&mpt4>=0.0
                &&mpt4<9000.
                &&mtNt[imu]>00.
                &&mt2>00.
                &&mt4>00.
              ){
               hMpt2[icent]->Fill(mpt2);
               hMpt3[icent]->Fill(mpt3);
               hMpt4[icent]->Fill(mpt4);
            }
        }//imu


    }//iev

    ///plot mean as fcn of centrality
    TCanvas* cMpt[nCentralityBins];
    char cName[30];
    for(int icent=0;icent<nCentralityBins; ++icent){

        sprintf(cName,"mptCent%i",icent);
        cMpt[icent] = new TCanvas(cName,cName,600,600);
        double mean2 = hMpt2[icent]->GetMean();
        double mean3 = hMpt3[icent]->GetMean();
        double mean4 = hMpt4[icent]->GetMean();

        hMean2->SetBinContent(icent+1,mean2);
        hMean3->SetBinContent(icent+1,mean3);
        hMean4->SetBinContent(icent+1,mean4);

        cMpt[icent]->cd();
        hMpt4[icent]->SetMarkerColor(kBlue);
        hMpt3[icent]->SetMarkerColor(kGreen);
        hMpt2[icent]->SetMarkerColor(kRed);
        hMpt4[icent]->Draw("pe");
        hMpt3[icent]->Draw("pesame");
        hMpt2[icent]->Draw("pesame");
    }//icent
    
    TCanvas *cMean = new TCanvas("cMean","cMean",600,600);
    hMean2->SetMarkerColor(kRed);
    hMean3->SetMarkerColor(kGreen);
    hMean4->SetMarkerColor(kBlue);
    hMean2->Draw("p");
    hMean3->Draw("psame");
    hMean4->Draw("psame");
}
