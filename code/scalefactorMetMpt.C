///////////////////
//This macro studies the ratio 
//of the mean calo MET and muon-subtracted
//MPT from ID tracks. This is done 
//in order to recover neutral tracks
//in the MPT calculation
//@date Nov.11,2012
//////////////////////


#include "TDirectory.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TCanvas.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>

void scalefactorMetMpt(){

TString fileName = "/tmp/tbalestr/HISingleMuonWmunuMCDataOverlay.2012.11.11.root";
TChain* tree = new TChain("tree","tree");
int nFiles = tree->Add(fileName);
std::cout <<"Filling the DataSet for "<< fileName << "...Number of files: " << nFiles << std::endl;


//_file0->cd();
#ifdef __CINT__
  gROOT->LoadMacro("AtlasUtils.C");
#endif


TCanvas c0 = TCanvas("c0","c0",600,600);

TString cutsGen ="abs(mc_mu_charge)==1";

cutsGen += "&&abs(mc_mu_gen_eta)<2.5";
cutsGen += "&&abs(charge)==1&&centrality<=0.8&&truthMatched_muid==1";
cutsGen+="&&efMatched&&(EF_mu4_MSonly_L1TE50||EF_mu4_L1VTE50)";
cutsGen+="&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0";
cutsGen+="&&nu_pt>25&&abs(nu_pt).<9999.";
cutsGen+="&&mt>40.0&&mt<300.";
unsigned int entries = tree->GetEntries();
std::cout << " reading entries " << entries << std::endl;

tree->SetBranchStatus("*",0);
tree->SetBranchStatus("mu_muid_n",1);
tree->SetBranchStatus("prompt",1);
tree->SetBranchStatus("val",1);
tree->SetBranchStatus("pt",1);
tree->SetBranchStatus("mt",1);
tree->SetBranchStatus("eta",1);
tree->SetBranchStatus("phi",1);
tree->SetBranchStatus("eLoss",1);
tree->SetBranchStatus("scat",1);
tree->SetBranchStatus("charge",1);
tree->SetBranchStatus("centrality",1);
tree->SetBranchStatus("ZDY",1);
tree->SetBranchStatus("caloMET",1);
tree->SetBranchStatus("subtractedMPT",1);
tree->SetBranchStatus("nu_pt",1);
tree->SetBranchStatus("EF_mu4_MSonly_L1TE50",1);
tree->SetBranchStatus("EF_mu4_L1VTE50",1);
tree->SetBranchStatus("efMatched",1);
tree->SetBranchStatus("truthMatched_muid",1);
tree->SetBranchStatus("mc_mu_gen_pt",1);
tree->SetBranchStatus("mc_mu_gen_eta",1);
tree->SetBranchStatus("mc_mu_charge",1);

int mu_muid_n,prompt[30],val[30],ZDY[30],EF_mu4_MSonly_L1TE50,EF_mu4_L1VTE50,efMatched[30],truthMatched_muid[30];
float pt[30],mt[30],eta[30],phi[30],eLoss[30],scat[30],charge[30],centrality[30],caloMET[30],subtractedMPT[30],nu_pt[30],mc_mu_gen_pt[30],mc_mu_gen_eta[30],mc_mu_charge[30];

tree->SetBranchAddress("mu_muid_n",&mu_muid_n);
tree->SetBranchAddress("prompt",&prompt);
tree->SetBranchAddress("val",&val);
tree->SetBranchAddress("pt",&pt);
tree->SetBranchAddress("mt",&mt);
tree->SetBranchAddress("eta",&eta);
tree->SetBranchAddress("phi",&phi);
tree->SetBranchAddress("eLoss",&eLoss);
tree->SetBranchAddress("scat",&scat);
tree->SetBranchAddress("charge",&charge);
tree->SetBranchAddress("centrality",&centrality);
tree->SetBranchAddress("ZDY",&ZDY);
tree->SetBranchAddress("caloMET",&caloMET);
tree->SetBranchAddress("subtractedMPT",&subtractedMPT);
tree->SetBranchAddress("nu_pt",&nu_pt);
tree->SetBranchAddress("EF_mu4_MSonly_L1TE50",&EF_mu4_MSonly_L1TE50);
tree->SetBranchAddress("EF_mu4_L1VTE50",&EF_mu4_L1VTE50);
tree->SetBranchAddress("efMatched",&efMatched);
tree->SetBranchAddress("truthMatched_muid",&truthMatched_muid);
tree->SetBranchAddress("mc_mu_gen_pt",&mc_mu_gen_pt);
tree->SetBranchAddress("mc_mu_gen_eta",&mc_mu_gen_eta);
tree->SetBranchAddress("mc_mu_charge",&mc_mu_charge);

TH1F* hmpt = new TH1F("hmpt","hmpt",50,0.0,100.0);
TH1F* hmet = new TH1F("hmet","hmet",50,0.0,100.0);
TH1F* hscale = new TH1F("hscale","hscale",50,0.0,100.0);

tree->Draw("subtractedMPT>>hmpt",cutsGen);
tree->Draw("caloMET>>hmet",cutsGen,"same");

float meanMPT = hmpt->GetMean();
std::cout << "Subtracted MPT mean: " << meanMPT << std::endl;
float meanMET = hmet->GetMean();
std::cout << "MET mean: " << meanMET << std::endl;
float meanRatio = meanMET/meanMPT;
std::cout << "<MET>/<SubMPT>: " << meanRatio << std::endl;


for(int i=0; i<entries; i++){
	
	tree->LoadTree(i);
	tree->GetEntry(i);

	if( i%10000==0){
     	std::cout << " event " << i << std::endl;
     	//if (i==1000000) break; //temp hack
    	}

	for(int imu=0;imu<mu_muid_n;imu++){
	  	
		if(fabs(subtractedMPT)<9999.){
		float scalefactor = meanRatio*(subtractedMPT)+pt[imu];
		hscale->Fill(scalefactor);
		}
		else continue;
	}


}
TCanvas c1 = TCanvas("c1","c1",600,600);
hscale->Draw();
}
