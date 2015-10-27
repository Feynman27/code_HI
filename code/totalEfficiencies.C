#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooDataSet.h"
#include "RooBinning.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooFormulaVar.h"
#include "RooArgList.h"
#include "RooCategory.h"

#include "TDirectory.h"
#include "TInterpreter.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath> 
#include <vector> 

#include "AtlasUtils.C"
#include "AtlasStyle.C"

#include "TriggerEfficiencies.C"

using namespace RooFit;

///////////////////////////////////////////////////////////////////////////////
// selectGenPtEta
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectGenPtEta( RooDataSet* dataSet, double ptLow, double ptUpp, double etaLow, double etaUpp, bool doMirrorEta = false )
{
  TString cut = "muonGenPt>";
  cut += ptLow;
  cut += "&&muonGenPt<";
  cut += ptUpp;
  cut += "&&((etaGen>";
  cut += etaLow;
  cut += "&&etaGen<";
  cut += etaUpp;
  cut += ")";
  if ( doMirrorEta ) {
    cut += "||(etaGen>";
    cut += -etaUpp;
    cut += "&&etaGen<";
    cut += -etaLow;
    cut += ")";
  }
  cut += ")";

  dataSet = (RooDataSet*) dataSet->reduce(cut);
  return dataSet;
}


///////////////////////////////////////////////////////////////////////////////
// fillHIMuonRecSet
///////////////////////////////////////////////////////////////////////////////
RooDataSet* fillHIMuonRecSet(const TString& pathName, const TString& fileName, RooArgSet& muonArgSet)
{
  RooDataSet* set = new RooDataSet(fileName,fileName,muonArgSet);
  
  float eLossNt[50];
  float scatNt[50];
  float compNt[50];
  float ptNt[50];
  float mtNt[50];
  float etaNt[50];
  float chargeNt[50];
  float centralityNt;
  float nu_ptNt;
  float ptconeNt[50];
  int valNt[50], ZDYNt[50], truthMatchedNt[50], promptNt[50];
  int nmu;

  // --- Load the MC tree ---
  TChain* tree = new TChain("tree","tree");
  int nFiles = tree->Add(pathName+fileName);
  
  std::cout << "Filling the RooDataSet for " << fileName << "... Number of files: " << nFiles << std::endl;
 
  // --- Set branch adresses ---
  tree->SetBranchAddress("eLoss", &eLossNt);
  //dR<0.3,pTid>3GeV
  tree->SetBranchAddress("ptcone20ID3", &ptconeNt);
  //vary ptcone variable for systematics
//  tree->SetBranchAddress("ptcone30ID3", &ptconeNt);
  tree->SetBranchAddress("scat", &scatNt);
  tree->SetBranchAddress("comp", &compNt);
  tree->SetBranchAddress("pt", &ptNt);
  tree->SetBranchAddress("mt", &mtNt);
  tree->SetBranchAddress("eta", &etaNt);
  tree->SetBranchAddress("charge", &chargeNt);
  tree->SetBranchAddress("prompt", &promptNt);
  tree->SetBranchAddress("val", &valNt); 
  tree->SetBranchAddress("ZDY", &ZDYNt); 
  tree->SetBranchAddress("centrality", &centralityNt);
  tree->SetBranchAddress("nu_pt", &nu_ptNt);
  tree->SetBranchAddress("truthMatched_muid", &truthMatchedNt);
  tree->SetBranchAddress("mu_muid_n", &nmu);

   // --- Set branch status ---
  tree->SetBranchStatus("*",0) ;
  tree->SetBranchStatus("mu_muid_n", 1);
  tree->SetBranchStatus("truthMatched_muid", 1);
  tree->SetBranchStatus("eLoss", 1);
  tree->SetBranchStatus("ptcone20ID3", 1);
//  tree->SetBranchStatus("ptcone30ID3", 1);
  tree->SetBranchStatus("scat", 1);
  tree->SetBranchStatus("comp", 1);
  tree->SetBranchStatus("pt", 1);
  tree->SetBranchStatus("mt", 1);
  tree->SetBranchStatus("eta", 1);
  tree->SetBranchStatus("charge", 1);
  tree->SetBranchStatus("val", 1); 
  tree->SetBranchStatus("ZDY", 1); 
  tree->SetBranchStatus("centrality", 1);
  tree->SetBranchStatus("nu_pt", 1);
  tree->SetBranchStatus("prompt", 1);

  std::cout << "Number of entries: " << tree->GetEntries() << std::endl; 
  for ( int i = 0; i < tree->GetEntries(); i++ ) {
//     if (i>100000) return set; // temp hack
    tree->LoadTree(i);
    tree->GetEntry(i);

    muonArgSet.setRealValue("missPt",nu_ptNt);
    muonArgSet.setRealValue("centrality",centralityNt);

    for (int imu = 0; imu<nmu;imu++){

      if ( promptNt[imu] == 23) muonArgSet.setCatLabel("muonCategory","Z");
      else if ( promptNt[imu] == 24) muonArgSet.setCatLabel("muonCategory","W");

      muonArgSet.setRealValue("muonELoss",eLossNt[imu]);
      muonArgSet.setRealValue("muonScat",scatNt[imu]);
      muonArgSet.setRealValue("muonQuality",valNt[imu]);

      double isolationTemp = -9999.0;
      isolationTemp = ptconeNt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolation",isolationTemp);
      muonArgSet.setRealValue("muonPt",ptNt[imu]);
      muonArgSet.setRealValue("muonMt",mtNt[imu]);
      muonArgSet.setRealValue("muonEta",etaNt[imu]);
      muonArgSet.setRealValue("ZDY",ZDYNt[imu]);
      muonArgSet.setRealValue("muonCharge",chargeNt[imu]);
      muonArgSet.setRealValue("muonCategory",promptNt[imu]);
      muonArgSet.setRealValue("muonGenRecMatched",truthMatchedNt[imu]);
      if ( chargeNt[imu] > 0 ) muonArgSet.setCatLabel("chargeCategory","muPlus");
      else if ( chargeNt[imu] < 0) muonArgSet.setCatLabel("chargeCategory","muMinus");
      set->add(muonArgSet);    
   }
  }

  return set;
}


RooDataSet* fillHIMuonGenSet(const TString& pathName, const TString& fileName, RooArgSet& muonGenSet)
{
  RooDataSet* set = new RooDataSet(fileName,fileName,muonGenSet);
  
  float nuGenPtNt[50], nuPhiGenNt[50];
  float muGenPtNt[50], muPhiGenNt[50];
  float mtGenNt[50];
  int motherNt[50];
  int daughterNt[50];
  float etaGenNt[50];
  float centralityNt;
  float chargeGenNt[50];
  int nmu;

  // --- Load the MC tree ---
  TChain* tree = new TChain("tree","tree");
  int nFiles = tree->Add(pathName+fileName);
  
  std::cout << "Filling the RooDataSet for " << fileName << "... Number of files: " << nFiles << std::endl;
 
  // --- Set branch adresses ---
  tree->SetBranchAddress("mc_mu_gen_pt", &muGenPtNt);
  tree->SetBranchAddress("mc_nu_gen_pt", &nuGenPtNt);
  tree->SetBranchAddress("mc_mu_gen_mothertype", &motherNt);
  tree->SetBranchAddress("mc_mu_gen_type", &daughterNt);
  tree->SetBranchAddress("mc_mu_charge", &chargeGenNt);
  tree->SetBranchAddress("mc_mu_gen_eta", &etaGenNt);
  tree->SetBranchAddress("mc_mu_gen_phi", &muPhiGenNt);
  tree->SetBranchAddress("mc_nu_gen_phi", &nuPhiGenNt);
  tree->SetBranchAddress("mc_mu_n", &nmu);
  tree->SetBranchAddress("centrality", &centralityNt);

   // --- Set branch status ---
  tree->SetBranchStatus("*",0) ;
  tree->SetBranchStatus("mc_mu_n", 1);
  tree->SetBranchStatus("centrality", 1);
  tree->SetBranchStatus("mc_mu_gen_pt", 1);
  tree->SetBranchStatus("mc_nu_gen_pt", 1);
  tree->SetBranchStatus("mc_mu_gen_mothertype", 1);
  tree->SetBranchStatus("mc_mu_gen_type", 1);
  tree->SetBranchStatus("mc_mu_charge", 1);
  tree->SetBranchStatus("mc_mu_gen_eta", 1);
  tree->SetBranchStatus("mc_mu_gen_phi", 1);
  tree->SetBranchStatus("mc_nu_gen_phi", 1);

  std::cout << "Number of entries: " << tree->GetEntries() << std::endl; 
  for ( int i = 0; i < tree->GetEntries(); i++ ) {
//     if (i>100000) return set; // temp hack
    tree->LoadTree(i);
    tree->GetEntry(i);

    muonGenSet.setRealValue("centrality",centralityNt);

    for (int imu = 0; imu<nmu;imu++){

      float dPhi = nuPhiGenNt[imu]-muPhiGenNt[imu];
      if(dPhi<-1*TMath::Pi()) dPhi += TMath::TwoPi(); if(dPhi>TMath::Pi()) dPhi -= TMath::TwoPi(); //fold btwn [-pi,pi]
      mtGenNt[imu] = TMath::Sqrt(2.0*muGenPtNt[imu]*nuGenPtNt[imu]*(1.0-TMath::Cos(dPhi)));
      muonGenSet.setRealValue("munuGenMt",mtGenNt[imu]);

      muonGenSet.setRealValue("muonGenPt",muGenPtNt[imu]);
      muonGenSet.setRealValue("nuGenPt",nuGenPtNt[imu]);
      muonGenSet.setRealValue("mother",motherNt[imu]);
      muonGenSet.setRealValue("daughter",daughterNt[imu]);
      muonGenSet.setRealValue("chargeGen",chargeGenNt[imu]);
      muonGenSet.setRealValue("etaGen",etaGenNt[imu]);

      if ( chargeGenNt[imu] > 0 ) muonGenSet.setCatLabel("chargeGenCategory","muPlus");
      else if ( chargeGenNt[imu] < 0) muonGenSet.setCatLabel("chargeGenCategory","muMinus");

      set->add(muonGenSet);    
   }
  }

  return set;
}

void totalEfficiencies(){


bool doEta = false;
bool doCentrality = false;
bool doPt = false;

///generate dictionaries
gInterpreter->GenerateDictionary("vector<RooFormulaVar>","RooFormulaVar.h;vector");
gInterpreter->GenerateDictionary("vector<TString>","TString.h;vector");
gInterpreter->GenerateDictionary("vector<RooArgList>","RooArgList.h;vector");
gInterpreter->GenerateDictionary("vector<RooDataSet*>","RooFormulaVar.h;vector");

/// --- output file ---
TString fileNameDataOut;

TString fileNameIn = "HISingleMuonWmunuPYTHIADataOverlay.04.12.2013";

///IMPORTANT:IF CHANGING ETA BINNING, YOU MUST CHANGE
///THE .txt FILE HOLDING THE TRIGGER EFFICIENCIES ACCORDINGLY
std::vector <double> etaBins;
etaBins.push_back(0.0);
//etaBins.push_back(0.1);
//etaBins.push_back(-2.4);
if(doEta){

/*	etaBins.push_back(+0.25);
	etaBins.push_back(+0.50);
	etaBins.push_back(+0.75);
	etaBins.push_back(+1.00);
	etaBins.push_back(+1.25);
	etaBins.push_back(+1.50);
	etaBins.push_back(+1.75);
	etaBins.push_back(+2.00);
	etaBins.push_back(+2.25);
*/

    
/*	etaBins.push_back(-2.35);
	etaBins.push_back(-2.3);
	etaBins.push_back(-2.25);
	etaBins.push_back(-2.2);
	etaBins.push_back(-2.15);
	etaBins.push_back(-2.1);
	etaBins.push_back(-2.05);
	etaBins.push_back(-2.0);
	etaBins.push_back(-1.975);
	etaBins.push_back(-1.950);
	etaBins.push_back(-1.925);
	etaBins.push_back(-1.9);
	etaBins.push_back(-1.8);
	etaBins.push_back(-1.7);
	etaBins.push_back(-1.6);
	etaBins.push_back(-1.5);
	etaBins.push_back(-1.4);
	etaBins.push_back(-1.35);
	etaBins.push_back(-1.3);
	etaBins.push_back(-1.25);
	etaBins.push_back(-1.2);
	etaBins.push_back(-1.15);
	etaBins.push_back(-1.1);
	etaBins.push_back(-1.05);
	etaBins.push_back(-1.0);
	etaBins.push_back(-0.9);
	etaBins.push_back(-0.8);
	etaBins.push_back(-0.7);
	etaBins.push_back(-0.6);
	etaBins.push_back(-0.5);
	etaBins.push_back(-0.4);
	etaBins.push_back(-0.3);
	etaBins.push_back(-0.2);
	etaBins.push_back(-0.1);
	etaBins.push_back(-0.075);
	etaBins.push_back(-0.05);
    etaBins.push_back(-0.025);

	etaBins.push_back(0.0);
*/

/*    etaBins.push_back(0.35);
    etaBins.push_back(0.6);
    etaBins.push_back(0.8);
    etaBins.push_back(1.05);
    etaBins.push_back(1.3);
    etaBins.push_back(1.55);
//        etaBins.push_back(1.73);
    etaBins.push_back(1.85);
    etaBins.push_back(2.1);

*/


	etaBins.push_back(0.025);
	etaBins.push_back(0.05);
	etaBins.push_back(0.075);
	etaBins.push_back(0.1);
	etaBins.push_back(0.2);
	etaBins.push_back(0.3);
	etaBins.push_back(0.4);
	etaBins.push_back(0.5);
	etaBins.push_back(0.6);
	etaBins.push_back(0.7);
	etaBins.push_back(0.8);
	etaBins.push_back(0.9);
	etaBins.push_back(1.0);
	etaBins.push_back(1.05);
	etaBins.push_back(1.1);
	etaBins.push_back(1.15);
	etaBins.push_back(1.2);
	etaBins.push_back(1.25);
	etaBins.push_back(1.3);
	etaBins.push_back(1.35);
	etaBins.push_back(1.4);
	etaBins.push_back(1.5);
	etaBins.push_back(1.6);
	etaBins.push_back(1.7);
	etaBins.push_back(1.8);
	etaBins.push_back(1.9);
	etaBins.push_back(1.925);
	etaBins.push_back(1.950);
	etaBins.push_back(1.975);
	etaBins.push_back(2.0);
	etaBins.push_back(2.05);
	etaBins.push_back(2.1);
	etaBins.push_back(2.15);
	etaBins.push_back(2.2);
	etaBins.push_back(2.25);
	etaBins.push_back(2.3);
	etaBins.push_back(2.35);
    
    etaBins.push_back(2.4);
}
etaBins.push_back(2.7);
etaBins.push_back(10.0);

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

//Npart
std::vector <double> npartBins;
if(doCentrality){
	npartBins.push_back(382.16);//0-5
	npartBins.push_back(330.26);//5-10
	npartBins.push_back(281.88);//10-15
	npartBins.push_back(239.52);//15-20
	npartBins.push_back(157.83);//20-40
	npartBins.push_back(45.93);//40-80
   
}
else npartBins.push_back(111.63); //0-80 (weighted by bin width)

const int npartNBins = npartBins.size();

if(nEtaBins==9&&nCentralityBins==6) readInputFile(nCentralityBins,nEtaBins,"triggerEffZNote_v06.txt");
else if (nEtaBins==38&&nCentralityBins==6) readInputFile(nCentralityBins,nEtaBins,"triggerEffZNote_v04.txt");
else std::cout << "WARNING: Trigger bins not initiliazed." << std::endl;
///if not mirroring eta
//readInputFile(nCentralityBins,nEtaBins,"triggerEffWNote_NoAbsEta_v01.txt");

// --- declare variables at generator level --- //
  RooRealVar  muonGenPt("muonGenPt","p_{T}",0.0,250.0,"GeV");
  RooRealVar  nuGenPt("nuGenPt","p_{T}^{#nu}",0.0,250.0,"GeV");
  RooRealVar  munuGenMt("munuGenMt","m_{T}",0.0,300.0,"GeV");
  RooRealVar  mother("mother","mother",-30.0,30.0);
  RooRealVar  daughter("daughter","daughter",-20.0,20.0);
  RooRealVar  chargeGen("chargeGen","chargeGen",-2.0,2.0);
  RooRealVar  etaGen("etaGen","etaGen",-3.0,3.0);
// --- declare cut variables at reco level --- //
  RooRealVar  muonPt("muonPt","p_{T}",0.0,250.0,"GeV");
  RooRealVar  missPt("missPt","p_{T}^{miss}",0.0,250.0,"GeV");
  RooRealVar  muonMt("muonMt","m_{T}",0.0,250.0,"GeV");
  RooRealVar  muonCharge("muonCharge","charge",-2.0,+2.0);
  RooRealVar  isolation("isolation","isolation",0.0,10.0);
  RooRealVar  centrality("centrality","centrality",0.,1.0);
  RooRealVar  muonEta("muonEta","muonEta",-3.0,+3.0);
  RooRealVar  ZDY("ZDY","ZDY",0.0,1.0);
  RooRealVar  muonQuality("muonQuality","muonQuality",0.0,15.0);
  RooRealVar  muonGenRecMatched("muonGenRecMatched","muonGenRecMatched",0.0,1.0);
  RooRealVar  muonELoss("muonELoss","muonELoss",-10.0,+10.0);
  RooRealVar  muonScat("muonScat","muonScat",-10.0,+10.0);

  TString sCutsGen = "abs(mother)==24&&abs(daughter)==13";
  RooArgList muonGenArgList(mother,daughter);
  RooFormulaVar cutsGen("cutsGen", "cutsGen", sCutsGen, muonGenArgList);

  TString sCutsFid = "abs(mother)==24&&abs(daughter)==13";
  RooFormulaVar cutsFid("cutsFid", "cutsFid", sCutsFid, RooArgList(mother,daughter));

  //reconstruction level cuts
  const unsigned int nCuts = 15;
  std::vector<TString>* vecCuts = new std::vector<TString>();

  TString sRecCutFlow0 = "muonGenRecMatched==1&&muonQuality>0"; vecCuts->push_back(sRecCutFlow0); 
  TString sRecCutFlow1 = sRecCutFlow0+"&&abs(muonEta)<2.4" ; vecCuts->push_back(sRecCutFlow1);
  TString sRecCutFlow2 = sRecCutFlow1+"&&muonQuality>1"; vecCuts->push_back(sRecCutFlow2);
  TString sRecCutFlow3 = sRecCutFlow2+"&&muonQuality>5"; vecCuts->push_back(sRecCutFlow3);
  TString sRecCutFlow4 = sRecCutFlow3+"&&muonQuality>6"; vecCuts->push_back(sRecCutFlow4);
  TString sRecCutFlow5 = sRecCutFlow4+"&&muonQuality>7"; vecCuts->push_back(sRecCutFlow5);
  TString sRecCutFlow6 = sRecCutFlow5+"&&muonQuality>9"; vecCuts->push_back(sRecCutFlow6);
  TString sRecCutFlow7 = sRecCutFlow6+"&&muonQuality>11"; vecCuts->push_back(sRecCutFlow7);
  TString sRecCutFlow8 = sRecCutFlow7+"&&abs(muonScat)<4.0"; vecCuts->push_back(sRecCutFlow8);
  TString sRecCutFlow9 = sRecCutFlow8+"&&abs(muonELoss)<0.5"; vecCuts->push_back(sRecCutFlow9);
  TString sRecCutFlow10 = sRecCutFlow9+"&&muonPt>25.0"; vecCuts->push_back(sRecCutFlow10);
  TString sRecCutFlow11 = sRecCutFlow10+"&&isolation<0.1"; vecCuts->push_back(sRecCutFlow11);
  TString sRecCutFlow12 = sRecCutFlow11+"&&ZDY==0"; vecCuts->push_back(sRecCutFlow12);
  TString sRecCutFlow13 = sRecCutFlow12+"&&missPt>25.0&&missPt<9000.0"; vecCuts->push_back(sRecCutFlow13);
  TString sRecCutFlow14 = sRecCutFlow13+"&&muonMt>40.0"; vecCuts->push_back(sRecCutFlow14);

  std::vector <RooArgList>* vecArgList  = new std::vector<RooArgList>();

  RooArgList muonRecArgList0(muonGenRecMatched,muonQuality); 
  vecArgList->push_back(muonRecArgList0);

  //RooArgList muonRecArgList1 = muonRecArgList0; muonRecArgList1.add(muonEta);
  RooArgList muonRecArgList1(muonQuality,muonGenRecMatched,muonEta); /*muonRecArgList1.add(muonEta);*/
  vecArgList->push_back(muonRecArgList1);

  RooArgList muonRecArgList2 = muonRecArgList1; 
  vecArgList->push_back(muonRecArgList2);

  RooArgList muonRecArgList3 = muonRecArgList1; 
  vecArgList->push_back(muonRecArgList3);

  RooArgList muonRecArgList4 = muonRecArgList1; 
  vecArgList->push_back(muonRecArgList4);

  RooArgList muonRecArgList5 = muonRecArgList1; 
  vecArgList->push_back(muonRecArgList5);

  RooArgList muonRecArgList6 = muonRecArgList1; 
  vecArgList->push_back(muonRecArgList6);

  RooArgList muonRecArgList7 = muonRecArgList1; 
  vecArgList->push_back(muonRecArgList7);

  RooArgList muonRecArgList8 = muonRecArgList1; muonRecArgList8.add(muonScat); 
  vecArgList->push_back(muonRecArgList8);

  RooArgList muonRecArgList9 = muonRecArgList8; muonRecArgList9.add(muonELoss); 
  vecArgList->push_back(muonRecArgList9);

  RooArgList muonRecArgList10 = muonRecArgList9; muonRecArgList10.add(muonPt); 
  vecArgList->push_back(muonRecArgList10);

  RooArgList muonRecArgList11 = muonRecArgList10; muonRecArgList11.add(isolation); 
  vecArgList->push_back(muonRecArgList11);

  RooArgList muonRecArgList12 = muonRecArgList11; muonRecArgList12.add(ZDY); 
  vecArgList->push_back(muonRecArgList12);

  RooArgList muonRecArgList13 = muonRecArgList12; muonRecArgList13.add(missPt); 
  vecArgList->push_back(muonRecArgList13);

  RooArgList muonRecArgList14 = muonRecArgList13; muonRecArgList14.add(muonMt); 
  vecArgList->push_back(muonRecArgList14);

  ///cuts on datasets at RECONSTRUCTION level
  std::vector <RooFormulaVar>* vecFormList = new std::vector<RooFormulaVar>();

  char buffer[50];
  for (int i=0; i<nCuts; ++i){
    sprintf(buffer,"cutsRecCutFlow%i",i);
    TString sTemp = (*vecCuts)[i];
    std::cout << "Cuts at index " << i << ": " << sTemp << std::endl;
    RooArgList tempRooArgList = (*vecArgList)[i];
    tempRooArgList.Print();
    RooFormulaVar tempRooFormulaVar(TString::Format(buffer),TString::Format(buffer), (*vecCuts)[i], (*vecArgList)[i]);
    vecFormList->push_back(tempRooFormulaVar);
        
  }

  RooCategory muonCategory("muonCategory","muonCategory");
  muonCategory.defineType("Z",23);
  muonCategory.defineType("W",24);

  RooCategory chargeCategory("chargeCategory","chargeCategory") ;
  chargeCategory.defineType("muMinus",-1) ;
  chargeCategory.defineType("muPlus",1) ;

  RooCategory chargeGenCategory("chargeGenCategory","chargeGenCategory") ;
  chargeGenCategory.defineType("muMinus",-1) ;
  chargeGenCategory.defineType("muPlus",1) ;

  RooArgSet muonGenArgSet(muonGenPt,nuGenPt,munuGenMt,mother,daughter,chargeGenCategory,etaGen,centrality);

  RooArgSet muonRecArgSet(muonEta,centrality,ZDY,muonCategory,chargeCategory);
  muonRecArgSet.add(muonPt);
  muonRecArgSet.add(missPt);
  muonRecArgSet.add(muonMt);
  muonRecArgSet.add(muonCharge);
  muonRecArgSet.add(isolation);
  muonRecArgSet.add(muonQuality);
  muonRecArgSet.add(muonGenRecMatched);
  muonRecArgSet.add(muonELoss);
  muonRecArgSet.add(muonScat);


  double ptmax = 300.0;
  // --- Set pt and eta bins ---
  std::vector<double> ptBins;
  ptBins.push_back(0.0);
  if(doPt){
  	ptBins.push_back(50.0);
  	//ptBins.push_back(75.0);
  }
  ptBins.push_back(ptmax);
  const int nPtBins = ptBins.size()-1;

  RooDataSet* mcWGenSet = fillHIMuonGenSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonGenArgSet); mcWGenSet->Print();
  mcWGenSet = (RooDataSet*)mcWGenSet->reduce(Cut(cutsGen)); 
  mcWGenSet = (RooDataSet*)mcWGenSet->reduce(Cut("centrality<0.8")); 
  std::cout << "Number of Wmunu evts at generator level : " << mcWGenSet->numEntries() << std::endl;

  RooDataSet* mcWGenSetPlus = (RooDataSet*)mcWGenSet->reduce("chargeGenCategory==chargeGenCategory::muPlus"); 
  RooDataSet* mcWGenSetMinus = (RooDataSet*)mcWGenSet->reduce("chargeGenCategory==chargeGenCategory::muMinus"); 

  RooDataSet* mcWFidSet = fillHIMuonGenSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonGenArgSet); 
  mcWFidSet = (RooDataSet*)mcWFidSet->reduce(Cut(cutsFid)); 
  mcWFidSet = (RooDataSet*)mcWFidSet->reduce(Cut("centrality<0.8")); 
  std::cout << "Number of Wmunu evts at generator level : " << mcWFidSet->numEntries() << std::endl;

  //construct charged sets

  RooDataSet* mcWFidSetPlus = (RooDataSet*)mcWFidSet->reduce("chargeGenCategory==chargeGenCategory::muPlus"); 
  RooDataSet* mcWFidSetMinus = (RooDataSet*)mcWFidSet->reduce("chargeGenCategory==chargeGenCategory::muMinus"); 

  double fidPlus = mcWFidSetPlus->numEntries(); double genPlus = mcWGenSetPlus->numEntries();
  double fidMinus = mcWFidSetMinus->numEntries(); double genMinus = mcWGenSetMinus->numEntries();

  std::cout << " Aw+ = " << mcWFidSetPlus->numEntries() << "/" << mcWGenSetPlus->numEntries() << " = " << fidPlus/genPlus << std::endl;
  std::cout << " Aw- = " << mcWFidSetMinus->numEntries() << "/" << mcWGenSetMinus->numEntries() << " = " << fidMinus/genMinus << std::endl;

  RooBinning b = RooBinning(0.0,300.0);
  b.addUniform(90,0.0,90.0);
  b.addUniform(53,90.0,300.0);

  std::vector <RooDataSet*> *vecDataSets = new std::vector<RooDataSet*>();

   for(int i=0; i<nCuts; ++i){

        RooDataSet* tempDataSet = fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root", muonRecArgSet); 
        RooFormulaVar tempRooFormulaVar = (*vecFormList)[i];
        tempDataSet = (RooDataSet*)tempDataSet->reduce( Cut(tempRooFormulaVar) );
        tempDataSet = (RooDataSet*)tempDataSet->reduce(Cut("muonCategory==muonCategory::W"));
        tempDataSet = (RooDataSet*)tempDataSet->reduce(Cut("centrality<0.8")); 
        vecDataSets->push_back(tempDataSet);
        //std::cout << "Number of Wmunu evts reconstructed for reco set " << i << " = " << tempDataSet->numEntries() << std::endl; 
        RooDataSet* pDataSets = (*vecDataSets)[i];
        std::cout << "Number of Wmunu evts reconstructed for reco set " << i << " = " << pDataSets->numEntries() << std::endl; 
   }
   ///cout abs and rel eff
   for(int i=0; i<nCuts; ++i){
        double denom = mcWFidSet->numEntries();
        RooDataSet* pDataSets = (*vecDataSets)[i];
        double events1 = pDataSets->numEntries();
        std::cout << "Absolute efficiency of reco set " << i << " = " << events1/denom*100 << "%"<< std::endl; 
        RooDataSet* pDataSets2 = 0;
        pDataSets2 =(*vecDataSets)[i+1];
        double events2 = pDataSets2->numEntries();
        std::cout << "Relative efficiency of reco set " << i+1 << " relative to " << i << " = " << events2/events1*100 << "%"<< std::endl;
   }

    delete vecCuts;
    delete vecArgList;
    delete vecFormList;
    delete vecDataSets;
    
} 

