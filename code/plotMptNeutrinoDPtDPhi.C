#include "TH1F.h"
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

void assignMean(double& meanData, double& meanMc, const TH1F* hData, const TH1F* hMc){

   meanData = hData->GetMean(); 
   meanMc = hMc->GetMean();
   std::cout << "Data mean : " << meanData << " MC mean : " << meanMc << std::endl;
}

void assignSigma(double& sigmaData, double& sigmaMc, const TH1F* hData, const TH1F* hMc){

    sigmaData = hData->GetRMS();
    sigmaMc = hMc->GetRMS();

    std::cout << "Data sigma : " << sigmaData << " MC sigma: " << sigmaMc << std::endl;
}

void fill1DHisto(TString fileNameIn, TH1F* h, double centralityLo,double centralityUpp,double etaLo,double etaUpp, bool isMc=false){
       ///reco level variables
      std::cout << centralityLo*100. << "<Centrality<"<<centralityUpp*100. << std::endl;
      std::cout << etaLo << "<|eta|<" << etaUpp << std::endl;
      std::cout << std::endl;

      float eLossNt[50];
      float scatNt[50]; 
      float compNt[50];
      float ptNt[50];
      float mtNt[50];
      float etaNt[50];
      float phiNt[50];
      float mpt;
      float chargeNt[50];
      int promptNt[50];
      float centralityNt;
      float ptconeNt[50];
      int valNt[50], recGenMatched[50],ZDYNt[50], matched1[50], matched2[50], matched3[50];
      int nmu,trig1,trig2,trig3;

      TChain* tree = new TChain("tree","tree");
      tree->Add(fileNameIn);

      std::cout << "Filling 1D histogram. " << std::endl;
      std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;
      tree->SetBranchAddress("centrality", &centralityNt);
      tree->SetBranchAddress("mu_muid_n", &nmu);
      tree->SetBranchAddress("nu_pt", &mpt);
      if(isMc)tree->SetBranchAddress("prompt",&promptNt);
      tree->SetBranchAddress("eta",&etaNt);
      tree->SetBranchAddress("pt",&ptNt);
      if(!isMc){
        tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC",&trig1);
        tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10",&trig2);
        tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20",&trig3);
        tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC_Matched20",&matched1);
        tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10_Matched20",&matched2);
        tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20_Matched20",&matched3);
      }
      tree->SetBranchAddress("eLoss", &eLossNt);
      tree->SetBranchAddress("ptcone20ID3", &ptconeNt);
      tree->SetBranchAddress("scat", &scatNt);
      tree->SetBranchAddress("mt", &mtNt);
      tree->SetBranchAddress("val", &valNt); 
      tree->SetBranchAddress("ZDY", &ZDYNt); 
      tree->SetBranchAddress("centrality", &centralityNt);
      if(isMc)tree->SetBranchAddress("truthMatched_muid", &recGenMatched);

      // --- Set branch status ---
      tree->SetBranchStatus("*",0) ;
      tree->SetBranchStatus("centrality", 1);
      tree->SetBranchStatus("mu_muid_n", 1);
      if(isMc) {
          tree->SetBranchStatus("truthMatched_muid", 1);
          tree->SetBranchStatus("prompt",1);
      }
      tree->SetBranchStatus("nu_pt", 1);
      tree->SetBranchStatus("eta",1);
      tree->SetBranchStatus("pt",1);
      if(!isMc){
        tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC",1);
        tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10",1);
        tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20",1);
        tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC_Matched20",1);
        tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10_Matched20",1);
        tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20_Matched20",1);
      }
      tree->SetBranchStatus("eLoss", 1);
      tree->SetBranchStatus("ptcone20ID3", 1);
      tree->SetBranchStatus("scat", 1);
      tree->SetBranchStatus("mt", 1);
      tree->SetBranchStatus("val", 1); 
      tree->SetBranchStatus("ZDY", 1); 
      tree->SetBranchStatus("centrality", 1);

      for(int iev = 0; iev < tree->GetEntries(); ++iev){
 
        tree->LoadTree(iev);
        tree->GetEntry(iev);
        //if(iev==10000) break; //temp hack

        if(!(centralityNt>centralityLo&&centralityNt<centralityUpp)) continue;
        for(int imu=0; imu<nmu; ++imu){

           if(
                valNt[imu]>11
                &&fabs(scatNt[imu])<4.0
                &&fabs(eLossNt[imu])<0.5
                &&fabs(etaNt[imu])>0.1
                &&fabs(etaNt[imu])<2.4
                &&centralityNt>0.
                &&centralityNt<0.8
                &&ptconeNt[imu]/ptNt[imu]<0.1
                &&ZDYNt[imu]==0
                &&ptNt[imu]>25.0
                &&fabs(mpt)<9999.0
                &&mtNt[imu]>40.0
                && ( (isMc) || ( (trig1&&matched1[imu]) ||
                (trig2&&matched2[imu])
                ||(trig3&&matched3[imu]) ) )
                && ((!isMc)||(promptNt[imu]==24&&recGenMatched[imu]==1))
                
             ){

              h->Fill(mpt);
            }

        }//imu

      }//iev

}

void fillHistos(TString fileNameIn,TH1F* hDelPhi,TH1F* hDelPt,TH2F* hPtPhi,double shift, double smear, double centralityLo,double centralityUpp,double etaLo,double etaUpp){

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
               //if(!(promptNt[imu]==24&&fabs(etaNt[imu])>etaLo&&fabs(etaNt[imu])<etaUpp)) continue; 
               if(!(promptNt[imu]==24)) continue; 
                 ///Calculate dPhi btwn mpt phi and neutrino phi
                 ///if this event has a reco muon from a W
                 float dPhi = mptPhi-neutrinoPhi[igen];
                 ///Fold dPhi in[-PI,PI]
                 if(dPhi>TMath::Pi()) dPhi-=2.*TMath::Pi();
                 else if(dPhi<-1.0*TMath::Pi()) dPhi+=2.*TMath::Pi();
                 hDelPhi->Fill(dPhi,1/(centralityUpp-centralityLo));

                 ///Calculate dPt btwn the reco mpt and true neutrino pt
                 float dPt = (mpt+shift-neutrinoPt[igen])*smear; 
                 hDelPt->Fill(dPt,1/(centralityUpp-centralityLo));

                 ///Plot in dPhi-dPt space
                 hPtPhi->Fill(dPt,dPhi);
             }//imu 
         }
    }//igen
   } //iev

}


void plotMptNeutrinoDPtDPhi(){


   TFile* outFile = new TFile("mptHistosForSystematics.root","recreate");

   TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
   TString fileNameMCIn = baseString+ "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.04.12.2013.root";
   TString fileNameDataIn = baseString+"HardProbesFiles/HISingleMuonHardProbesData.04.17.2013.root";

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

   int marker[] = {kGreen,kBlue,kRed,kCyan,kViolet,kAzure-9};
   TString sCent[] = {"0-5%","5-10%","10-15%","15-20%","20-40%","40-80%"};
   std::cout << "Running over " << nCentralityBins << " centrality classes and " << nEtaBins << " eta windows." << std::endl;

   TCanvas *cPhi = new TCanvas("cPhi","cPhi",600,600);
   TCanvas *cPt = new TCanvas("cPt","cPt",600,600);
   ///array of delta phi histos
   TH1F* hDelPhi[nEtaBins][nCentralityBins],*hDelPt[nEtaBins][nCentralityBins];
   TH1F *hMptData[nEtaBins][nCentralityBins],*hMptMC[nEtaBins][nCentralityBins];
   TH2F* hPtPhi[nEtaBins][nCentralityBins];
   TCanvas* cPtPhi[nEtaBins][nCentralityBins];
   TPaletteAxis *palette[nEtaBins][nCentralityBins];

   char sDelPhi[50],sDelPt[50],sPtPhi[50],cName[50],sMcMpt[50],sDataMpt[50];
   for(int ieta=0; ieta<nEtaBins; ++ieta){

        TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);

    for(int icent=0; icent<nCentralityBins; ++icent){

        sprintf(sDelPhi,"dPhiMptNeutrino_Eta%i_Cent%i",ieta,icent);
        sprintf(sDelPt,"dPtMptNeutrino_Eta%i_Cent%i",ieta,icent);
        sprintf(sMcMpt,"mcMpt_Eta%i_Cent%i",ieta,icent);
        sprintf(sDataMpt,"dataMpt_Eta%i_Cent%i",ieta,icent);
        sprintf(sPtPhi,"h2DPtPhiNeutrino_Eta%i_Cent%i",ieta,icent);
        sprintf(cName,"c2D_Eta%i_Cent%i",ieta,icent);


        hDelPhi[ieta][icent] = new TH1F(sDelPhi,sDelPhi,50,-3.2,+3.2);
        hDelPt[ieta][icent] = new TH1F(sDelPt,sDelPt,100,-100.0,+100.0); 
        hMptData[ieta][icent]= new TH1F(sDataMpt,sDataMpt,100,-100.0,+100.0);
        hMptMC[ieta][icent]= new TH1F(sMcMpt,sMcMpt,100,-100.0,+100.0);
        hPtPhi[ieta][icent] = new TH2F(sPtPhi,sPtPhi,100,-100.0,+100.0,50,-3.2,+3.2);

        hPtPhi[ieta][icent]->GetXaxis()->SetLabelSize(0.04);
        hPtPhi[ieta][icent]->GetYaxis()->SetLabelSize(0.04);
        hPtPhi[ieta][icent]->GetZaxis()->SetLabelSize(0.04);

        cPtPhi[ieta][icent] = new TCanvas(cName,cName,700,700);
        cPtPhi[ieta][icent]->Range(-145.3061,-4.501695,138.7755,3.59774);
        cPtPhi[ieta][icent]->SetFillColor(0);
        cPtPhi[ieta][icent]->SetBorderMode(0);
        cPtPhi[ieta][icent]->SetBorderSize(2);
        cPtPhi[ieta][icent]->SetLogz();
        cPtPhi[ieta][icent]->SetTickx(1);
        cPtPhi[ieta][icent]->SetTicky(1);
        cPtPhi[ieta][icent]->SetLeftMargin(0.1594828);
        cPtPhi[ieta][icent]->SetRightMargin(0.1293103);
        cPtPhi[ieta][icent]->SetTopMargin(0.04910714);
        cPtPhi[ieta][icent]->SetBottomMargin(0.1607143);
        cPtPhi[ieta][icent]->SetFrameBorderMode(0);
        cPtPhi[ieta][icent]->SetFrameBorderMode(0);
        cPtPhi[ieta][icent]->SetLogz();

        palette[ieta][icent] = new TPaletteAxis(100.4082,-3.212053,113.0612,3.187947,hPtPhi[ieta][icent]);
        palette[ieta][icent]->SetLabelColor(1);
        palette[ieta][icent]->SetLabelFont(42);
        palette[ieta][icent]->SetLabelOffset(0.005);
        palette[ieta][icent]->SetLabelSize(0.04);
        palette[ieta][icent]->SetTitleOffset(1);
        palette[ieta][icent]->SetTitleSize(0.05);
        palette[ieta][icent]->SetFillColor(100);
        palette[ieta][icent]->SetFillStyle(1001);

        fill1DHisto(fileNameDataIn,hMptData[ieta][icent],centralityBins[icent],centralityBins[icent+1],etaBins[ieta],etaBins[ieta+1]);
        fill1DHisto(fileNameMCIn,hMptMC[ieta][icent],centralityBins[icent],centralityBins[icent+1],etaBins[ieta],etaBins[ieta+1],true);

        ///Normalize MC to data
        double sf = hMptData[ieta][icent]->Integral()/hMptMC[ieta][icent]->Integral();
        hMptMC[ieta][icent]->Scale(sf);

        double meanMptData=-9999.,meanMptMC=+9999.,sigmaMptData=-9999.,sigmaMptMC=0.;
        assignMean(meanMptData,meanMptMC,hMptData[ieta][icent],hMptMC[ieta][icent]);
        assignSigma(sigmaMptData,sigmaMptMC,hMptData[ieta][icent],hMptMC[ieta][icent]);

        double shiftParameter = meanMptData-meanMptMC;
        double scaleFactor = sigmaMptData/sigmaMptMC;
        std::cout << "Mpt mean shifted by " << shiftParameter << " and scaled by " << scaleFactor << std::endl;
        fillHistos(fileNameMCIn,hDelPhi[ieta][icent],hDelPt[ieta][icent],hPtPhi[ieta][icent],shiftParameter,scaleFactor,
                    centralityBins[icent],centralityBins[icent+1],etaBins[ieta],etaBins[ieta+1]);

        leg->SetTextFont(gStyle->GetTextFont());
        leg->SetTextSize(gStyle->GetTextSize());
        leg->SetBorderSize(0);
        leg->SetFillColor(0);
        leg->AddEntry(hDelPhi[ieta][icent],sCent[icent],"pe");

        if(icent==0) {
            cPhi->cd();
            leg->Draw();

            hDelPhi[ieta][icent]->SetMarkerColor(marker[icent]);
            hDelPhi[ieta][icent]->GetYaxis()->SetRangeUser(0.0,1.0e5);
            hDelPhi[ieta][icent]->Draw("pe");
            hDelPhi[ieta][icent]->GetXaxis()->SetTitle("#Delta#phi(#slash{p_{T}}-#nu_{#mu})");
            cPt->cd();
            hDelPt[ieta][icent]->SetMarkerColor(marker[icent]);
            hDelPt[ieta][icent]->GetYaxis()->SetRangeUser(0.0,2.6e4);
            hDelPt[ieta][icent]->Draw("pe");
            hDelPt[ieta][icent]->GetXaxis()->SetTitle("#Deltap_{T}(#slash{p_{T}}-#nu_{#mu})[GeV]");
        }
        else {
            cPhi->cd();
            leg->Draw();
            hDelPhi[ieta][icent]->SetMarkerColor(marker[icent]);
            hDelPhi[ieta][icent]->Draw("pesame");
            cPt->cd();
            hDelPt[ieta][icent]->SetMarkerColor(marker[icent]);
            hDelPt[ieta][icent]->Draw("pesame");
        }
        cPtPhi[ieta][icent]->cd();
        gStyle->SetPalette(1);
        TH2F* hPtPhic = (TH2F*)hPtPhi[ieta][icent]->DrawNormalized("colz",1.0);
        hPtPhic->GetYaxis()->SetTitle("#Delta#phi(#vec{#slash{p_{T}}}-#vec{#nu_{#mu}})");
        hPtPhic->GetXaxis()->SetTitle("#Deltap_{T} = (|#slash{p_{T}}|+#Delta_{Data,MC}-|#nu_{#mu}|) #frac{#sigma_{Data}}{#sigma_{MC}}[GeV]");

        TString centLabelTemp = "";
        centLabelTemp += centralityBins[icent]*100; centLabelTemp+="-";centLabelTemp+=centralityBins[icent+1]*100; centLabelTemp+="%";
        TLatex *   tex = new TLatex(0.4123563,0.3735119,centLabelTemp);
        tex->SetNDC();
        tex->SetTextFont(42);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.4267241,0.3020833,"MC11");
        tex->SetNDC();
        tex->SetTextFont(42);
        tex->SetLineWidth(2);
        tex->Draw();

        stringstream ss;
        string sSave;
        ss<<sPtPhi; ss>>sSave; sSave+="_07_20_2013.pdf";
        cPtPhi[ieta][icent]->Print((TString)sSave);
        outFile->cd();
        hPtPhic->Write();

        //delete hDelPhi[ieta][icent];
    }//icent

   }//ieta

   outFile->Close();
}
