#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TRandom.h"

#include <iostream>
#include <vector>
#include <cmath>

void smearGeneratorPt(){

    float smear = 0.1;
   
    TFile* fIn = new TFile("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/HIWtaumuNtuple.07.09.2013.root","read");
    TFile* outFile = new TFile("smearedMuonPtFromTau.root","recreate");
    TTree* tree = (TTree*)fIn->Get("truth");
	double xBins[85] ;
	
    double ptLow = 0.0;
	float ptmax = 200.0;
	double xlo = ptLow; double xhi = 60.0;  
	double nbins = 60;
	double binw = (xhi-xlo)/nbins;
	int ib = 0;
	for(double i=xlo ; i<xhi ; i+=binw){
		
			xBins[ib] = i;
			ib++;	
	}
	xlo = 60.0; xhi = 90.0;  
	nbins = 15;
	binw = (xhi-xlo)/nbins;
	for(double i=xlo ; i<xhi ; i+=binw){
		
			xBins[ib] = i;
			ib++;	
	}
	xlo = 90.0; xhi = 110.0;  
	nbins = 4;
	binw = (xhi-xlo)/nbins;
	for(double i=xlo ; i<xhi ; i+=binw){
		
			xBins[ib] = i;
			ib++;	
	}
	xlo = 110.0; xhi = 150.0;  
	nbins = 4;
	binw = (xhi-xlo)/nbins;
	for(double i=xlo ; i<xhi ; i+=binw){
		
			xBins[ib] = i;
			ib++;	
	}
	xlo = 150.0; xhi = ptmax; 
	nbins = 1;
	binw = (xhi-xlo)/nbins;
	for(double i=xlo ; i<=xhi ; i+=binw){

		xBins[ib] = i;
		ib++;
	}	

	int binnum = sizeof(xBins)/sizeof(double)-1 ;
    std::cout << "Number of bins : " << binnum << std::endl;

    TH1F* hSmearedMuonPt = new TH1F("hSmearedMuonPt","hSmearedMuonPt",binnum,xBins);

//    std::vector<int>* mc_pdgId = new std::vector<int>();
    int mc_n;
//    std::vector<float>* mc_pt = new std::vector<float>();
//    std::vector<float>* mc_eta = new std::vector<float>();
    float mc_pt,mc_eta;
    float mc_pdgId,mc_mptPt,mc_mt,centrality;
    tree->SetBranchAddress("mc_pdgId",&mc_pdgId);
    tree->SetBranchAddress("mc_pt",&mc_pt);
    tree->SetBranchAddress("mc_eta",&mc_eta);
    tree->SetBranchAddress("mc_mptPt",&mc_mptPt);
    tree->SetBranchAddress("mc_mt",&mc_mt);
    tree->SetBranchAddress("centrality",&centrality);
//    tree->SetBranchAddress("mc_n",&mc_n);
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("mc_pdgId",1);
    tree->SetBranchStatus("mc_pt",1);
    tree->SetBranchStatus("mc_eta",1);
    tree->SetBranchStatus("mc_mptPt",1);
    tree->SetBranchStatus("mc_mt",1);
    tree->SetBranchStatus("centrality",1);
//    tree->SetBranchStatus("mc_n",1);

    for(int igen=0; igen<tree->GetEntries(); ++igen){
        tree->GetEntry(igen);

//        for(int igen=0; igen<mc_n; ++igen){

//            std::cout << "muon pt " << mc_pt->at(igen)/1000. << std::endl;
            //if(abs(mc_pdgId->at(igen))==13&&abs(mc_eta->at(igen))>0.1&&abs(mc_eta->at(igen))<2.4){
            if(fabs(mc_pdgId)==13
                && fabs(mc_eta)>0.1
                &&fabs(mc_eta)<2.4
                &&mc_mptPt>25.0
                &&mc_mt>40.0
                &&centrality<0.8
                ){
              float mcPtSmeared = mc_pt;
              ///smear the generated pt
              mcPtSmeared +=gRandom->Gaus()*smear;

              hSmearedMuonPt->Fill(mcPtSmeared);
            }
//        }

    }

/*    delete mc_pdgId;
    delete mc_pt;
    delete mc_eta;
    */
    outFile->Write();
    outFile->Close();
}
