#include "TDirectory.h"
#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath> 

#include "AtlasUtils.C"
#include "AtlasStyle.C"



void binMigrationCalculator(){

    bool doEta = true;
    bool doCentrality = false;

    TFile *outFile = new TFile("binMigrationFraction.root","RECREATE");
    TString fileNameIn = "HISingleMuonWmunuPYTHIADataOverlay.04.12.2013.root";
    TFile* fIn = new TFile("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/"+fileNameIn,"READ"); 
    TTree* tree = (TTree*)fIn->Get("tree");


    std::vector <double> etaBins;
    etaBins.push_back(0.1);
    if(doEta){
    
    etaBins.push_back(0.35);
    etaBins.push_back(0.6);
    etaBins.push_back(0.8);
    etaBins.push_back(1.05);
    etaBins.push_back(1.3);
    etaBins.push_back(1.55);
    etaBins.push_back(1.85);
    etaBins.push_back(2.1);

	}
    etaBins.push_back(2.4);

    const int nEtaBins = etaBins.size()-1;
    std::cout << "nEtaBins : " << nEtaBins << std::endl;

    //centrality
    std::vector <double> centBins;
     centBins.push_back(0.0);
    if(doCentrality){
        centBins.push_back(0.05);
        centBins.push_back(0.1);
        centBins.push_back(0.15);
        centBins.push_back(0.2);
        centBins.push_back(0.4);
        
    }
    centBins.push_back(0.8);

    const int nCentralityBins = centBins.size()-1;
    std::cout << "nCentralityBins : " << nCentralityBins << std::endl;

    TString psCuts = "val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&truthMatched_muid==1";
    TString sigCuts = "&&prompt==24&&ZDY==0&&ptcone20ID3/pt<0.1&&pt>25.0&&nu_pt>25.0&&nu_pt<9000.0&&mt>40.0";

    TH1D* hRec[nEtaBins][nCentralityBins] ;
    TH1D* hRecGenInRange[nEtaBins][nCentralityBins];
    TH1D* hRecGenUp[nEtaBins][nCentralityBins];
    TH1D* hRecGenLow[nEtaBins][nCentralityBins];
    TGraphErrors* grPurity = new TGraphErrors(nEtaBins);

    double xBins[] = {0.1,0.35,0.6,0.8,1.05,1.3,1.55,1.85,2.1,2.4}; 
    int  binnum = sizeof(xBins)/sizeof(double) - 1;


    for(int icent=0; icent<nCentralityBins; ++icent){

        char cPurity[50];
        sprintf(cPurity,"hPurityEtaDepCentrality%i",icent);

        for(int ieta=0; ieta<nEtaBins; ++ieta){

          char cRec[50];
          char cRecGen[50]; 
          char cRecGenUp[50]; 
          char cRecGenLow[50];

           fIn->cd();
           sprintf(cRec,"hRecEta%iCentrality%i",ieta,icent);
           sprintf(cRecGen,"hRecGenEta%iCentrality%i",ieta,icent);
           sprintf(cRecGenUp,"hRecGenUpEta%iCentrality%i",ieta,icent);
           sprintf(cRecGenLow,"hRecGenLowEta%iCentrality%i",ieta,icent);

           hRec[ieta][icent] = new TH1D(cRec,cRec,binnum,xBins); 
           hRecGenInRange[ieta][icent] = new TH1D(cRecGen,cRecGen,binnum,xBins); 
           hRecGenUp[ieta][icent] = new TH1D(cRecGenUp,cRecGenUp,binnum,xBins); 
           hRecGenLow[ieta][icent] = new TH1D(cRecGenLow,cRecGenLow,binnum,xBins); 

           TString recCuts = psCuts+sigCuts;
           recCuts+="&&abs(eta)>"; recCuts+=etaBins[ieta];
           recCuts+="&&abs(eta)<="; recCuts+=etaBins[ieta+1];

           recCuts+="&&centrality>"; recCuts+=centBins[icent];
           recCuts+="&&centrality<"; recCuts+=centBins[icent+1];

           std::cout << "cuts:" << recCuts << std::endl;
           TString recCutsGenInRange = recCuts+"&&abs(mc_mu_gen_eta)>"; recCutsGenInRange+=etaBins[ieta];
           recCutsGenInRange+="&&abs(mc_mu_gen_eta)<="; recCutsGenInRange+=etaBins[ieta+1];
           std::cout << recCutsGenInRange << std::endl;

           TString recCutsGenLow = recCuts+"&&abs(mc_mu_gen_eta)<="; recCutsGenLow+=etaBins[ieta];
           std::cout << recCutsGenLow << std::endl;

           TString recCutsGenUpp = recCuts+"&&abs(mc_mu_gen_eta)>"; recCutsGenUpp+=etaBins[ieta+1];
           std::cout << recCutsGenUpp << std::endl;
            
           tree->Draw((TString("eta>>")+cRec).Data(),recCuts);
           tree->Draw((TString("eta>>")+cRecGen).Data(),recCutsGenInRange);
           tree->Draw((TString("eta>>")+cRecGenLow).Data(),recCutsGenLow);
           tree->Draw((TString("eta>>")+cRecGenUp).Data(),recCutsGenUpp);

           double nRec = hRec[ieta][icent]->GetEntries();
           double nRecGenInRange = hRecGenInRange[ieta][icent]->GetEntries();
           double nRecGenLow = hRecGenLow[ieta][icent]->GetEntries();
           double nRecGenUpp = hRecGenUp[ieta][icent]->GetEntries();

           std::cout << std::endl;
           std::cout << "Percentage of muons with generator level muons in range for " << etaBins[ieta] << "-" <<
            etaBins[ieta+1] << " = " << nRecGenInRange/nRec*100. << "%" << std::endl;             
           std::cout << "Percentage of muons with generator level muons with eta < " << etaBins[ieta] << " = " <<
                nRecGenLow/nRec*100. << "%" << std::endl;             
           std::cout << "Percentage of muons with generator level muons with eta > " << etaBins[ieta+1] << " = "
            << nRecGenUpp/nRec*100. << "%" << std::endl;             

           double stat = TMath::Sqrt(1./nRecGenInRange+1./nRec);
           grPurity->SetPoint(ieta, etaBins[ieta]+1./2.*(etaBins[ieta+1]-etaBins[ieta]),nRecGenInRange/nRec*100.0);
           grPurity->SetPointError(ieta,1./2.*(etaBins[ieta+1]-etaBins[ieta]),stat);

           std::cout << std::endl;

           outFile->cd();
           hRec[ieta][icent]->Write(cRec);
           hRecGenInRange[ieta][icent]->Write(cRecGen);
           hRecGenUp[ieta][icent]->Write(cRecGenUp);
           hRecGenLow[ieta][icent]->Write(cRecGenLow);

           delete hRec[ieta][icent];
           delete hRecGenInRange[ieta][icent];
           delete hRecGenUp[ieta][icent];
           delete hRecGenLow[ieta][icent];
        }
           outFile->cd();
           grPurity->Write(cPurity);


          delete grPurity;
    }

}
