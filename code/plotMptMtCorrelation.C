#include "TH1F.h"
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

void plotMptMtCorrelation(){

      TDatime* time = new TDatime();
      TString sDate = "";
      sDate+=time->GetMonth(); sDate+="_"; sDate+=time->GetDay(); sDate+="_"; sDate+=time->GetYear();
      std::cout << "Today is " << sDate << std::endl;
      
      bool doMc = false;

      TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
      TString fileNameIn = baseString; 
      if(doMc) fileNameIn+="MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.04.12.2013.root";
      else fileNameIn+="HardProbesFiles/HISingleMuonHardProbesData.04.17.2013.root";

      TFile* outFile = new TFile("mptMtCorrelations.root","recreate");

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
      float mpt;
      float ptconeNt[50];
      int valNt[50], truthMatchedNt[50],ZDYNt[50], matched1[50], matched2[50], matched3[50];
      int nmu,trig1,trig2,trig3,trig4,trig5;

      TChain* tree = new TChain("tree","tree");
      tree->Add(fileNameIn);
      std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;

      tree->SetBranchAddress("pt",&ptNt);
      if(doMc)tree->SetBranchAddress("prompt",&promptNt);
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
      tree->SetBranchAddress("nu_pt", &mpt);
      tree->SetBranchAddress("mu_muid_n", &nmu);


      // --- Set branch status ---
      tree->SetBranchStatus("*",0) ;
      tree->SetBranchStatus("mu_muid_n", 1);
      tree->SetBranchStatus("eLoss", 1);
      tree->SetBranchStatus("ptcone20ID3", 1);
      tree->SetBranchStatus("scat", 1);
      tree->SetBranchStatus("comp", 1);
      if(doMc)tree->SetBranchStatus("prompt",1);
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

      char hName[50];
      const int nCentralityBins = 6;
      TString sCent[] = {"0-5%","5-10%","10-15%","15-20%","20-40%","40-80%"};
      TCanvas* _c[nCentralityBins];
      TPaletteAxis *palette[nCentralityBins];

      TH2F* h2D[nCentralityBins];
      for(int icent=0; icent<nCentralityBins; ++icent){
        sprintf(hName,"h2DMptMtCorrelationCent%i",icent);
        h2D[icent] = new TH2F(hName,hName,70.0,0.0,210.0,70.0,0.0,210.0);
      }//icent

      for(int iev = 0; iev < tree->GetEntries(); ++iev){
 
        tree->LoadTree(iev);
        tree->GetEntry(iev);

        int icent = getCentralityBin(centralityNt);
        if(icent<0) continue;

        for(int imu=0; imu<nmu; ++imu){

            if(
                (promptNt[imu]==24||!doMc) 
                && valNt[imu]>11  
                && abs(scatNt[imu])<4.0
                && abs(eLossNt[imu])<0.5
                && abs(etaNt[imu])>0.1
                && abs(etaNt[imu])<2.4
                &&ptconeNt[imu]/ptNt[imu]<0.1
                && ZDYNt[imu]==0
                &&ptNt[imu]>25.0
                &&mpt>0.0
                &&mpt<9000.
                &&mtNt[imu]>0.0

              )
              h2D[icent]->Fill(mpt,mtNt[imu]);
        }//imu
       }//iev

       gStyle->SetPalette(1);
       char cName[50];
       for(int icent=0; icent<nCentralityBins; ++icent){
        
            sprintf(cName,"cCent%i",icent);
            _c[icent] = new TCanvas(cName,cName,600,600);
            _c[icent]->Range(-145.3061,-4.501695,138.7755,3.59774);
            _c[icent]->SetFillColor(0);
            _c[icent]->SetBorderMode(0);
            _c[icent]->SetBorderSize(2);
            _c[icent]->SetLogz();
            _c[icent]->SetTickx(1);
            _c[icent]->SetTicky(1);
            _c[icent]->SetLeftMargin(0.1594828);
            _c[icent]->SetRightMargin(0.1293103);
            _c[icent]->SetTopMargin(0.04910714);
            _c[icent]->SetBottomMargin(0.1607143);
            _c[icent]->SetFrameBorderMode(0);
            _c[icent]->SetFrameBorderMode(0);
            _c[icent]->SetLogz();

            palette[icent] = new TPaletteAxis(100.4082,-3.212053,113.0612,3.187947,h2D[icent]);
            palette[icent]->SetLabelColor(1);
            palette[icent]->SetLabelFont(42);
            palette[icent]->SetLabelOffset(0.005);
            palette[icent]->SetLabelSize(0.04);
            palette[icent]->SetTitleOffset(1);
            palette[icent]->SetTitleSize(0.05);
            palette[icent]->SetFillColor(100);
            palette[icent]->SetFillStyle(1001);

            h2D[icent]->GetXaxis()->SetLabelSize(0.04);
            h2D[icent]->GetYaxis()->SetLabelSize(0.04);
            h2D[icent]->GetZaxis()->SetLabelSize(0.04);
            h2D[icent]->Draw("colz");
            h2D[icent]->GetYaxis()->SetTitle("m_{T}[GeV]");
            h2D[icent]->GetXaxis()->SetTitle("#slash{p_{T}}[GeV]");
            TProfile* _pfx = h2D[icent]->ProfileX((TString)h2D[icent]->GetName()+"_pfx");
            _pfx->Draw("pesame");

            TLatex *   tex = new TLatex(0.4123563,0.3735119,sCent[icent]);
            tex->SetNDC();
            tex->SetTextFont(42);
            tex->SetLineWidth(2);
            tex->Draw();
            if(doMc)tex = new TLatex(0.4267241,0.3020833,"MC11");
            else tex = new TLatex(0.4267241,0.3020833,"Data11");
            tex->SetNDC();
            tex->SetTextFont(42);
            tex->SetLineWidth(2);
            tex->Draw();

            outFile->cd();
            h2D[icent]->Write();
            _pfx->Write();
            stringstream ss;
            string sSave;
            ss<<(h2D[icent]->GetName()); ss>>sSave; 
            if(doMc)sSave+="_MC_"; 
            else sSave+="_Data_";
            sSave+=sDate;
            _c[icent]->Print((TString)sSave+".pdf");

       }
}
