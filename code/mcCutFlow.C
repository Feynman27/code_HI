#include "TCanvas.h"
#include "TDirectory.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooAddPdf.h"
#include "RooBinning.h"
#include "RooFormulaVar.h"

//#include "WAnalysisHIDep.C"
#include "AtlasUtils.C"
#include "AtlasStyle.C"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath> 

#include "correctionFactorsDep.C"

using namespace RooFit;

///////////////////////////////////////////////////////////////////////////////
//weightDS 
///////////////////////////////////////////////////////////////////////////////
RooDataSet* weightDS(RooDataSet* set, RooRealVar& w){

  RooDataSet* wSet = new RooDataSet(set->GetName(),set->GetTitle(),set,*set->get(),0,w.GetName()) ;
  std::cout << "Returning weighted d.s." << std::endl;
  std::cout << std::endl;
  wSet->Print();
  return wSet;
    
}


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

  return (RooDataSet*) dataSet->reduce(cut);
}


///////////////////////////////////////////////////////////////////////////////
// fillHIEventDataSet(used only for counting number of events in each centrality class)
///////////////////////////////////////////////////////////////////////////////
RooDataSet* fillHIEventDataSet(const TString& pathName, const TString& fileName, RooArgSet& centralityArgSet)
{
  RooDataSet* set = new RooDataSet(fileName,fileName,centralityArgSet);
  
  float centralityNt;

  // --- Load the MC tree ---
  TChain* tree = new TChain("tree","tree");
  int nFiles = tree->Add(pathName+fileName);
  
  std::cout << "Filling the Jx event RooDataSet for " << fileName << "... Number of files: " << nFiles << std::endl;
 
  // --- Set branch adresses ---
  tree->SetBranchAddress("centrality", &centralityNt);
  tree->SetBranchStatus("centrality", 1);
  
  std::cout << "Number of entries: " << tree->GetEntries() << std::endl; 
  for ( int i = 0; i < tree->GetEntries(); i++ ) {
//     if (i>100000) break; // temp hack
    tree->LoadTree(i);
    tree->GetEntry(i);

    centralityArgSet.setRealValue("centrality",centralityNt);
    set->add(centralityArgSet);
  }

  return set;
}


///////////////////////////////////////////////////////////////////////////////
// selectCentrality
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectCentrality( RooDataSet* dataSet, double centralityLow, double centralityHigh )
{
  TString cut = "centrality>";
  cut += centralityLow;
  cut += "&&centrality<";
  cut += centralityHigh;
  
  return (RooDataSet*) dataSet->reduce(cut);
}


///////////////////////////////////////////////////////////////////////////////
// selectPtEta
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectPtEta( RooDataSet* dataSet, double ptLow, double ptUpp, double etaLow, double etaUpp, bool doMirrorEta = false )
{
  TString cut = "muonPt>";
  cut += ptLow;
  cut += "&&muonPt<";
  cut += ptUpp;
  cut += "&&((muonEta>";
  cut += etaLow;
  cut += "&&muonEta<";
  cut += etaUpp;
  cut += ")";
  if ( doMirrorEta ) {
    cut += "||(muonEta>";
    cut += -etaUpp;
    cut += "&&muonEta<";
    cut += -etaLow;
    cut += ")";
  }
  cut += ")";

  return (RooDataSet*) dataSet->reduce(cut);
}


///////////////////////////////////////////////////////////////////////////////
//selectPtEtaCentrality 
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectPtEtaCentrality( RooDataSet* dataSet,double ptLow, double ptUpp, double etaLow, double etaUpp, double centralityLow, double centralityUpp, 
	bool doMirrorEta=false,bool isGen=false){

  //cut on generated eta if this is a generator level cut
  if(isGen) dataSet= selectGenPtEta( dataSet,ptLow,ptUpp, etaLow, etaUpp, doMirrorEta);
  else dataSet = selectPtEta( dataSet,ptLow,ptUpp, etaLow, etaUpp, doMirrorEta);
  TString cut = "centrality>";
  cut += centralityLow;
  cut += "&&centrality<";
  cut+= centralityUpp;

  dataSet = (RooDataSet*)dataSet->reduce(cut);
  return dataSet;
}

///////////////////////////////////////////
//return number of collision in Jx sample
//////////////////////////////////////////
double getJxEvents(RooDataSet* mcJxEvents, double centralityLow, double centralityUpp){
	//event counting dataset
    RooDataSet* mcJxEventsTemp = selectCentrality(mcJxEvents,centralityLow,centralityUpp); mcJxEventsTemp->Print();
	return mcJxEventsTemp->numEntries();
}

///////////////////////////////////////////
//return weighted Jx histogram for variable-binned input(currently only J1-J3) 
//////////////////////////////////////////
TH1F* getWeightedJxVariableBinnedHisto(double xBins[], int arrLength, double centralityLow,double centralityUpp, double ncoll, TH1F* hmcJ1Set, TH1F* hmcJ2Set, 
                                        TH1F* hmcJ3Set, RooDataSet* mcJ1Events,RooDataSet* mcJ2Events, RooDataSet* mcJ3Events, 
                                        int nBins, double xLo, double varMax){

	std::cout << "Weighting Jx samples..." << std::endl;
	std::cout << "centrality: " << centralityLow << "-" << centralityUpp << " ncoll: " << ncoll << " nBins: " << nBins << " varMax " << varMax << std::endl;
	//return the number of MC collisions in this centrality bin
	double nMCJ1 = getJxEvents(mcJ1Events, centralityLow, centralityUpp);
	double nMCJ2 = getJxEvents(mcJ2Events, centralityLow, centralityUpp);
	double nMCJ3 = getJxEvents(mcJ3Events, centralityLow, centralityUpp);

	TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",arrLength,xBins);

	//weight the Jx samples according to cross-sections; takes into account prob of find muon+jet (AMI)
	double wtJ1 = 6.6523e-03*1.8770e-4/64.0e-3; double wtJ2 = 1.4941e-02*8.2788e-6/64.0e-3; double wtJ3 = 2.4284e-02*2.9426e-7/64.0e-3;
	double arrCentWidth = centralityUpp-centralityLow;
	double evData = 1.0e9; //number of hp events
	double scaleFactor = arrCentWidth*ncoll*evData;

    std::cout << "Bins must be equal. Let's check: " << hmcQCDSet->GetNbinsX() << " =? " << hmcJ1Set->GetNbinsX() << " =? "
        << hmcJ2Set->GetNbinsX() << " =? " << hmcJ3Set->GetNbinsX() << " =? " << std::endl;
	hmcQCDSet->Add(hmcJ1Set,hmcJ2Set,(wtJ1)/(nMCJ1)*scaleFactor, (wtJ2)/(nMCJ2)*scaleFactor); 
	hmcQCDSet->Add(hmcJ3Set,(wtJ3)/(nMCJ3)*scaleFactor);

	std::cout << "J1 coeff : " << (wtJ1)/(nMCJ1)*scaleFactor << " J2 coeff : " << (wtJ2)/(nMCJ2)*scaleFactor << " J2 coeff : " 
        << (wtJ3)/(nMCJ3)*scaleFactor << std::endl;
	std::cout << "MCevents: " << nMCJ1 << " " << nMCJ2 << " " << nMCJ3 << std::endl;

	return hmcQCDSet;
}

///////////////////////////////////////////////////////////////////////////////
// selectGenPtEta
///////////////////////////////////////////////////////////////////////////////
/*RooDataSet* selectGenPtEta( RooDataSet* dataSet, double ptLow, double ptUpp, double etaLow, double etaUpp, bool doMirrorEta = false )
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
*/

///////////////////////////////////////////////////////////////////////////////
//selectPtEtaCentrality 
///////////////////////////////////////////////////////////////////////////////
/*RooDataSet* selectPtEtaCentrality( RooDataSet* dataSet,double ptLow, double ptUpp, double etaLow, double etaUpp, double centralityLow, double centralityUpp, 
	bool isGen=false, bool doMirrorEta=false){

  //cut on generated eta if this is a generator level cut
  if(isGen) dataSet= selectGenPtEta( dataSet,ptLow,ptUpp, etaLow, etaUpp, doMirrorEta);
  else dataSet = selectPtEta( dataSet,ptLow,ptUpp, etaLow, etaUpp, doMirrorEta);
  TString cut = "centrality>";
  cut += centralityLow;
  cut += "&&centrality<";
  cut+= centralityUpp;

  dataSet = (RooDataSet*)dataSet->reduce(cut);
  return dataSet;
}
*/

///////////////////////////////////////////////////////////////////////////////
// fillHIWTauDataSet
///////////////////////////////////////////////////////////////////////////////
RooDataSet* fillHIWTauDataSet(const TString& pathName, const TString& fileName, RooArgSet& muonArgSet, double wt=1.0){
    
  RooDataSet* set = new RooDataSet(fileName,fileName,muonArgSet);
  // --- Load the MC tree ---
  TChain* tree = new TChain("tree","tree");
  int nFiles = tree->Add(pathName+fileName);
  
  std::cout << "Filling the RooDataSet for " << fileName << "... Number of files: " << nFiles << std::endl;
  int nmu;
  float mc_mt[50], mc_mptPhi,mc_mptPx,mc_mptPy,mc_mptPz,mc_mptPt,centrality;
  float  mc_pt[50],mc_ptGen[50],  mc_ptNominal[50],mc_eta[50],mc_phi[50],mc_pdgId[50],mc_charge[50];
  float mc_E[50],mc_M[50];

  tree->SetBranchAddress("nTauMu",&nmu);
  tree->SetBranchAddress("mc_mt",&mc_mt);
  tree->SetBranchAddress("mc_mptPhi",&mc_mptPhi);
  tree->SetBranchAddress("mc_mptPt",&mc_mptPt);
  tree->SetBranchAddress("centrality",&centrality);
  tree->SetBranchAddress("mc_pt",&mc_pt);
  tree->SetBranchAddress("mc_ptGen",&mc_ptGen);
  tree->SetBranchAddress("mc_eta",&mc_eta);
  tree->SetBranchAddress("mc_phi",&mc_phi);
  tree->SetBranchAddress("mc_pdgId",&mc_pdgId);
  tree->SetBranchAddress("mc_charge",&mc_charge);

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("nTauMu",1);
  tree->SetBranchStatus("mc_mt",1);
  tree->SetBranchStatus("mc_mptPhi",1);
  tree->SetBranchStatus("mc_mptPt",1);
  tree->SetBranchStatus("centrality",1);
  tree->SetBranchStatus("mc_pt",1);
  tree->SetBranchStatus("mc_ptGen",1);
  tree->SetBranchStatus("mc_eta",1);
  tree->SetBranchStatus("mc_phi",1);
  tree->SetBranchStatus("mc_pdgId",1);
  tree->SetBranchStatus("mc_charge",1);

  std::cout << "Number of entries: " << tree->GetEntries() << std::endl; 
  for ( int i = 0; i < tree->GetEntries(); i++ ) {
//     if (i>100000) break; // temp hack
    tree->LoadTree(i);
    tree->GetEntry(i);

    muonArgSet.setRealValue("missPt",mc_mptPt);
    muonArgSet.setRealValue("centrality",centrality);

    for (int imu = 0; imu<nmu; ++imu){

      muonArgSet.setRealValue("w",wt);
      muonArgSet.setRealValue("muonPt",mc_pt[imu]);
      muonArgSet.setRealValue("muonGenPt",mc_ptGen[imu]);
      muonArgSet.setRealValue("muonMt",mc_mt[imu]);
      muonArgSet.setRealValue("muonEta",mc_eta[imu]);
      muonArgSet.setRealValue("muonPhi",mc_phi[imu]);

      if (mc_charge[imu] > 0. ) muonArgSet.setCatLabel("chargeCategory","muPlus");
      else if (mc_charge[imu] < 0.) muonArgSet.setCatLabel("chargeCategory","muMinus");

      set->add(muonArgSet);    
    }//imu
   }//i
   
  return set;
}


///////////////////////////////////////////////////////////////////////////////
// fillHIMuonRecSet
///////////////////////////////////////////////////////////////////////////////
RooDataSet* fillHIMuonRecSet(const TString& pathName, const TString& fileName, RooArgSet& muonArgSet, double wt=1.0)
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

  std::cout << "Number of events in tree: " << tree->GetEntries() << std::endl; 
  for ( int i = 0; i < tree->GetEntries(); i++ ) {
//     if (i>100000) return set; // temp hack
    tree->LoadTree(i);
    tree->GetEntry(i);

    muonArgSet.setRealValue("missPt",nu_ptNt);
    muonArgSet.setRealValue("centrality",centralityNt);

    for (int imu = 0; imu<nmu;imu++){


      muonArgSet.setRealValue("muonELoss",eLossNt[imu]);
      muonArgSet.setRealValue("muonScat",scatNt[imu]);
      muonArgSet.setRealValue("muonQuality",valNt[imu]);

      double isolationTemp = -9999.0;
      isolationTemp = ptconeNt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolation",isolationTemp);
      muonArgSet.setRealValue("muonPt",ptNt[imu]);
      muonArgSet.setRealValue("w",wt);
      muonArgSet.setRealValue("muonMt",mtNt[imu]);
      muonArgSet.setRealValue("muonEta",etaNt[imu]);
      muonArgSet.setRealValue("ZDY",ZDYNt[imu]);
      muonArgSet.setRealValue("muonCharge",chargeNt[imu]);
      //muonArgSet.setRealValue("muonCategory",promptNt[imu]);
      muonArgSet.setRealValue("motherRec",promptNt[imu]);
      muonArgSet.setRealValue("muonGenRecMatched",truthMatchedNt[imu]);
      //if(promptNt[imu]==24||promptNt[imu]==23) std::cout << "Prompt = " <<  promptNt[imu] << std::endl;
      if ( promptNt[imu] == 23) muonArgSet.setCatLabel("muonCategory","Z",kTRUE);
      else if ( promptNt[imu] == 24) {
          muonArgSet.setCatLabel("muonCategory","W",kTRUE);
      }
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

void Write(TFile* const outFile, TObject* const gr, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing TGraph to root file..." << std::endl;
	  gr->Write(sFile);
	  std::cout << "Done." << std::endl;
   	  gDirectory = dir;
}

TEfficiency* writeEfficiency(TFile* fEff, const TH1& hPassed, const TH1& hTotal, TString sEff){
	if(TEfficiency::CheckConsistency(hPassed,hTotal,"w")) {

		TEfficiency* pEff = 0;
		fEff->cd();
		pEff = new TEfficiency(hPassed,hTotal);
		pEff->Write(sEff);
		return pEff;
		//pEff->Draw("AP");
	}
    else return NULL;
}

TGraphAsymmErrors* calculateEfficiency(const TH1* hPassed, const TH1* hTotal, int nBins){

   TGraphAsymmErrors* grTemp = new TGraphAsymmErrors(nBins);
   grTemp->BayesDivide(hPassed,hTotal);
   return grTemp;
}

void mcCutFlow()
{	
    //gROOT->LoadMacro("AtlasUtils.C");
    ///Switches 
	bool doCharge = false ;
	bool doCentrality = true ;
	bool doEta = true; 
//    if(doCharge||doEta||doCentrality){ std::cout << "Cannot bin in eta,centrality,and or charge. " << std::endl; exit(0)}
    ///Muon quality

	float mtmax = 300.0;
	float ptmax = 200.0;

	SetAtlasStyle();
	TFile* pFile = new TFile("mcCutFlowRatios.root","recreate");
	TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";

    ///Open W,Z,Jx files

    ///Wmunu
	TString fileNameMCWIn = "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.04.12.2013";
    
    ///Zmumu
	TString fileNameMCZIn = "MonteCarloFiles/Zmumu/HISingleMuonZmumuPYTHIAHIJINGOverlay.04.13.2013";

	//J1 1 muon-filter 
	TString fileNameMCJ1In = "MonteCarloFiles/QCD/HISingleMuonJ1_1muPYTHIADataOverlay.04.13.2013";
	//J2 1 muon-filter 
	TString fileNameMCJ2In = "MonteCarloFiles/QCD/HISingleMuonJ2_1muPYTHIADataOverlay.04.13.2013";
	//J3 1 muon-filter 
	TString fileNameMCJ3In = "MonteCarloFiles/QCD/HISingleMuonJ3_1muPYTHIADataOverlay.04.13.2013";

    ///W-->tau-->mu Generated sample
    TString fileNameMCTauIn = "MonteCarloFiles/Wtaumu/HIWtaumuNtuple.07.25.2013";

    ///Set up variables to be used
    // --- declare variables at generator level --- //
    RooRealVar  muonGenPt("muonGenPt","p_{T}",0.0,350.0,"GeV");
    RooRealVar  nuGenPt("nuGenPt","p_{T}^{#nu}",0.0,250.0,"GeV");
    RooRealVar  munuGenMt("munuGenMt","m_{T}",0.0,300.0,"GeV");
    RooRealVar  mother("mother","mother",-30.0,30.0);
    RooRealVar  daughter("daughter","daughter",-20.0,20.0);
    RooRealVar  chargeGen("chargeGen","chargeGen",-2.0,2.0);
    RooRealVar  etaGen("etaGen","etaGen",-10.0,10.0);
    // --- declare cut variables at reco level --- //
    RooRealVar  muonPt("muonPt","p_{T}",0.0,350.0,"GeV");
    RooRealVar  motherRec("motherRec","motherRec",0.0,50.0,"GeV");
    RooRealVar  missPt("missPt","p_{T}^{miss}",0.0,250.0,"GeV");
    RooRealVar  muonMt("muonMt","m_{T}",0.0,250.0,"GeV");
    RooRealVar  muonCharge("muonCharge","charge",-2.0,+2.0);
    RooRealVar  isolation("isolation","isolation",0.0,10.0);
    RooRealVar  centrality("centrality","centrality",0.,1.0);
    RooRealVar  muonEta("muonEta","muonEta",-3.0,+3.0);
    RooRealVar  ZDY("ZDY","ZDY",0.0,1.0);
    RooRealVar  muonQuality("muonQuality","muonQuality",0.0,15.0);
    RooRealVar  muonGenRecMatched("muonGenRecMatched","muonGenRecMatched",0.0,2.0);
    RooRealVar  muonELoss("muonELoss","muonELoss",-10.0,+10.0);
    RooRealVar  muonScat("muonScat","muonScat",-10.0,+10.0);
    ///RooRealVar for weighting W np,pn,pp,nn samples
    RooRealVar  w("w","w",0.0,10.0);

    ///////////////////////////////
    //W generator level cuts////////
    ///////////////////////////////

    RooArgList muonGenArgList(mother,daughter); muonGenArgList.add(centrality);
    ///Generator level cuts for W
    TString sCutsWGen = "abs(mother)==24&&abs(daughter)==13&&centrality>0.0&&centrality<0.8";
    RooFormulaVar cutsWGen("cutsWGen", "cutsWGen", sCutsWGen, muonGenArgList);

    ///Geometric fiducial cut at generator level (i.e. |eta|)
    muonGenArgList.add(etaGen);
    TString sCutsWGeomFid = sCutsWGen + "&&abs(etaGen)>0.1&&abs(etaGen)<2.4"; 
    RooFormulaVar cutsWGeomFid("cutsWGeomFid", "cutsWGeomFid", sCutsWGeomFid, muonGenArgList);

    ///////////////////////////////
    //Tau generator level cuts////////
    ///////////////////////////////
    RooArgList muonTauArgList(muonEta,centrality,muonPt,missPt,muonMt,muonGenPt,w);
    TString sCutsTauGen = "";
    sCutsTauGen = "centrality>0.0&&centrality<0.8";
	RooFormulaVar cutsTauGen("cutsTauGen", "cutsTauGen", sCutsTauGen, muonTauArgList);
    TString sCutsTauGeomFid = "abs(muonEta)>0.1&&abs(muonEta)<2.4&&centrality>0.0&&centrality<0.8";
	RooFormulaVar cutsTauGeomFid("cutsTauGeomFid", "cutsTauGeomFid", sCutsTauGeomFid, muonTauArgList);

    ///////////////////////////////
    //Z generator level cuts////////
    ///////////////////////////////

    ///Generator level cuts for Z
    TString sCutsZGen = "abs(mother)==23&&abs(daughter)==13&&centrality>0.0&&centrality<0.8";
    RooFormulaVar cutsZGen("cutsZGen", "cutsZGen", sCutsZGen, muonGenArgList);

    ///Geometric fiducial cut at generator level (i.e. |eta|)
    TString sCutsZGeomFid = sCutsZGen + "&&abs(etaGen)>0.1&&abs(etaGen)<2.4"; 
    RooFormulaVar cutsZGeomFid("cutsZGeomFid", "cutsZGeomFid", sCutsZGeomFid, muonGenArgList);

    ///Note: No defined generator level cut for Jx sample

    //////////////////////////////
    ///Reconstructed level cuts///
    /////////////////////////////

    ///Reco muon matched to generated
    RooArgList muonRecArgList0(muonGenRecMatched); muonRecArgList0.add(centrality);muonRecArgList0.add(muonEta);
    TString sCutsRec0 = "muonGenRecMatched==1&&centrality>0.0&&centrality<0.8&&abs(muonEta)>0.1&&abs(muonEta)<2.4";
    RooFormulaVar cutsRec0("cutsRec0","cutsRec0",sCutsRec0,muonRecArgList0); 

    ///Reco muon quality cuts
    muonRecArgList0.add(muonQuality); muonRecArgList0.add(muonELoss); muonRecArgList0.add(muonScat); 
    TString sCutsRec1 = sCutsRec0+"&&muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0";
    RooFormulaVar cutsRec1("cutsRec1","cutsRec1",sCutsRec1,muonRecArgList0); 

    ///Z veto
/*    muonRecArgList0.add(ZDY);
    TString sCutsRec2 = sCutsRec1+"&&ZDY==0";
    RooFormulaVar cutsRec2("cutsRec2","cutsRec2",sCutsRec2,muonRecArgList0); 
*/
    ///isolation cut
    muonRecArgList0.add(isolation);
    TString sCutsRec2 = sCutsRec1+"&&isolation<0.1";
    RooFormulaVar cutsRec2("cutsRec2","cutsRec2",sCutsRec2,muonRecArgList0); 

    muonRecArgList0.add(ZDY);
    TString sCutsRec3 = sCutsRec2+"&&ZDY==0";
    RooFormulaVar cutsRec3("cutsRec3","cutsRec3",sCutsRec3,muonRecArgList0); 

    ///missing pT
    muonRecArgList0.add(missPt);
    TString sCutsRec4 = sCutsRec3+"&&missPt>25.0&&missPt<9000.0";
    RooFormulaVar cutsRec4("cutsRec4","cutsRec4",sCutsRec4,muonRecArgList0); 

    ///mT
    muonRecArgList0.add(muonMt);
    TString sCutsRec5 = sCutsRec4+"&&muonMt>40.0";
    RooFormulaVar cutsRec5("cutsRec5","cutsRec5",sCutsRec5,muonRecArgList0); 
    
    ///Cuts for generator-level tau sample
    TString sCutsTauMpt = "abs(muonEta)>0.1&&abs(muonEta)<2.4&&centrality>0.0&&centrality<0.8&&muonPt>0.0&&missPt>25.0&&missPt<9000.0";
    TString sCutsTauMt = "abs(muonEta)>0.1&&abs(muonEta)<2.4&&centrality>0.0&&centrality<0.8&&muonPt>0.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0";
	RooFormulaVar cutsTauMpt("cutsTauMpt", "cutsTauMpt", sCutsTauMpt, muonTauArgList);
	RooFormulaVar cutsTauMt("cutsTauMt", "cutsTauMt", sCutsTauMt, muonTauArgList);



    ///Reco muon mothers
    RooCategory muonCategory("muonCategory","muonCategory");
    muonCategory.defineType("Z",23);
    muonCategory.defineType("W",24);
    muonCategory.Print();

    RooCategory chargeCategory("chargeCategory","chargeCategory") ;
    chargeCategory.defineType("muMinus",-1) ;
    chargeCategory.defineType("muPlus",1) ;
    chargeCategory.Print();

    RooCategory chargeGenCategory("chargeGenCategory","chargeGenCategory") ;
    chargeGenCategory.defineType("muMinus",-1) ;
    chargeGenCategory.defineType("muPlus",1) ;

    ///pT binning
    double ptLow = 0.0;
    RooBinning b = RooBinning(ptLow,ptmax);
    b.addUniform(60,ptLow,60.0);
    b.addUniform(15,60.0,90.0);
    b.addUniform(4,90.0,110.0);
    b.addUniform(4,110.0,150.0);
    b.addUniform(1,150.0,ptmax);
    //b.addUniform(1,200.0,ptmax);

	double xBins[85] ;
	
	double xlo = ptLow; double xhi = 60.0;  
	double nbins = 60;
	double binw = (xhi-xlo)/nbins;
	int ib = 0;
	for(double i=xlo ; i<xhi ; i+=binw){
		
			xBins[ib] = i;
            std::cout << "xBins at position " << ib << " = " << xBins[ib] << std::endl;
			ib++;	
	}
	xlo = 60.0; xhi = 90.0;  
	nbins = 15;
	binw = (xhi-xlo)/nbins;
	for(double i=xlo ; i<xhi ; i+=binw){
		
			xBins[ib] = i;
            std::cout << "xBins at position " << ib << " = " << xBins[ib] << std::endl;
			ib++;	
	}
	xlo = 90.0; xhi = 110.0;  
	nbins = 4;
	binw = (xhi-xlo)/nbins;
	for(double i=xlo ; i<xhi ; i+=binw){
		
			xBins[ib] = i;
            std::cout << "xBins at position " << ib << " = " << xBins[ib] << std::endl;
			ib++;	
	}
	xlo = 110.0; xhi = 150.0;  
	nbins = 4;
	binw = (xhi-xlo)/nbins;
	for(double i=xlo ; i<xhi ; i+=binw){
		
			xBins[ib] = i;
            std::cout << "xBins at position " << ib << " = " << xBins[ib] << std::endl;
			ib++;	
	}
	xlo = 150.0; xhi = ptmax; 
	nbins = 1;
	binw = (xhi-xlo)/nbins;
	for(double i=xlo ; i<=xhi ; i+=binw){

		xBins[ib] = i;
        std::cout << "xBins at position " << ib << " = " << xBins[ib] << std::endl;
		ib++;
	}	

/*    xlo = 200.0; xhi = ptmax; 
	nbins = 1;
	binw = (xhi-xlo)/nbins;
	for(double i=xlo ; i<=xhi ; i+=binw){

		xBins[ib] = i;
        std::cout << "xBins at position " << ib << " = " << xBins[ib] << std::endl;
		ib++;
	}	
*/
	int binnum = sizeof(xBins)/sizeof(double)-1 ;
    std::cout << "Number of bins : " << binnum << std::endl;

    ///Variable sets to be filled with values
    RooArgSet muonGenArgSet(muonGenPt,nuGenPt,munuGenMt,mother,daughter,chargeGenCategory,etaGen,centrality);

    RooArgSet muonTauArgSet(muonEta,centrality,muonPt,missPt,muonMt,muonGenPt,w);

    RooArgSet muonRecArgSet(muonEta,centrality,ZDY,muonCategory,chargeCategory);
    muonRecArgSet.add(muonPt);
    muonRecArgSet.add(motherRec);
    muonRecArgSet.add(missPt);
    muonRecArgSet.add(muonMt);
    muonRecArgSet.add(muonCharge);
    muonRecArgSet.add(isolation);
    muonRecArgSet.add(muonQuality);
    muonRecArgSet.add(muonGenRecMatched);
    muonRecArgSet.add(muonELoss);
    muonRecArgSet.add(muonScat);

    ///Fill base datasets at generator level for W and Z
    RooDataSet* mcWGenSet = fillHIMuonGenSet(baseString, fileNameMCWIn+".root",muonGenArgSet); 
    mcWGenSet = (RooDataSet*)mcWGenSet->reduce(Cut(cutsWGen));
    std::cout << "Number of entries in mcWGenSet = " << mcWGenSet->numEntries() << std::endl;

    RooDataSet* mcTauGenSet = fillHIWTauDataSet(baseString, fileNameMCTauIn+".root",muonTauArgSet);  
    mcTauGenSet = (RooDataSet*)mcTauGenSet->reduce(Cut(cutsTauGen)); mcTauGenSet->Print();
    std::cout << "Number of entries in mcTauGenSet = " << mcTauGenSet->numEntries() << std::endl;

    RooDataSet* mcZGenSet = fillHIMuonGenSet(baseString, fileNameMCZIn+".root",muonGenArgSet); mcZGenSet->Print();
    mcZGenSet = (RooDataSet*)mcZGenSet->reduce(Cut(cutsZGen));
    ///Make geometric fiducial generator level subsets
    RooDataSet* mcWGeomFidSet = (RooDataSet*)mcWGenSet->reduce(Cut(cutsWGeomFid));
    std::cout << "Number of entries in mcWGenSet after geometric fiducial cut = " << mcWGeomFidSet->numEntries() << std::endl;
    RooDataSet* mcZGeomFidSet = (RooDataSet*)mcZGenSet->reduce(Cut(cutsZGeomFid)); 
    RooDataSet* mcTauGeomFidSet = (RooDataSet*)mcTauGenSet->reduce(Cut(cutsTauGeomFid)); mcTauGeomFidSet->Print();

    ///Fill datasets at reconstruction level
    ///W boson
    ///Base reco dataset
    RooDataSet* mcWSet = fillHIMuonRecSet(baseString,fileNameMCWIn+".root", muonRecArgSet);  
    std::cout << "Number of entries in mcWSet = " << mcWSet->numEntries() << std::endl;
    //mcWSet = (RooDataSet*)mcWSet->reduce(muonRecArgSet,"muonCategory==muonCategory::W"); 
    mcWSet = (RooDataSet*)mcWSet->reduce(Cut("motherRec==24")); 
    std::cout << "Number of entries in mcWSet after W mother selection = " << mcWSet->numEntries() << std::endl;
    ///Apply the cut flow
    std::cout << "Number of entries in mcWSet = " << mcWSet->numEntries() << std::endl;
    RooDataSet* mcWRecSet0 = (RooDataSet*)mcWSet->reduce(Cut(cutsRec0)); ;
    std::cout << "Number of entries in mcWRecSet0 = " << mcWRecSet0->numEntries() << std::endl;
    RooDataSet* mcWRecSet1 = (RooDataSet*)mcWSet->reduce(Cut(cutsRec1)); 
    std::cout << "Number of entries in mcWRecSet1 = " << mcWRecSet1->numEntries() << std::endl;
    RooDataSet* mcWRecSet2 = (RooDataSet*)mcWSet->reduce(Cut(cutsRec2));
    std::cout << "Number of entries in mcWRecSet2 = " << mcWRecSet2->numEntries() << std::endl;
    RooDataSet* mcWRecSet3 = (RooDataSet*)mcWSet->reduce(Cut(cutsRec3)); 
    std::cout << "Number of entries in mcWRecSet3 = " << mcWRecSet3->numEntries() << std::endl;
    RooDataSet* mcWRecSet4 = (RooDataSet*)mcWSet->reduce(Cut(cutsRec4)); 
    std::cout << "Number of entries in mcWRecSet4 = " << mcWRecSet4->numEntries() << std::endl;
    RooDataSet* mcWRecSet5 = (RooDataSet*)mcWSet->reduce(Cut(cutsRec5)); 
    std::cout << "Number of entries in mcWRecSet5 = " << mcWRecSet5->numEntries() << std::endl;
    //exit(0);
    ///Z boson
    ///Base reco dataset
    RooDataSet* mcZSet = fillHIMuonRecSet(baseString,fileNameMCZIn+".root", muonRecArgSet); mcZSet->Print(); 
    //mcZSet = (RooDataSet*)mcZSet->reduce(Cut("muonCategory==muonCategory::Z")); mcZSet->Print();
    mcZSet = (RooDataSet*)mcZSet->reduce(Cut("motherRec==23")); 
    ///Apply the cut flow
    RooDataSet* mcZRecSet0 = (RooDataSet*)mcZSet->reduce(Cut(cutsRec0));
    RooDataSet* mcZRecSet1 = (RooDataSet*)mcZSet->reduce(Cut(cutsRec1));
    RooDataSet* mcZRecSet2 = (RooDataSet*)mcZSet->reduce(Cut(cutsRec2));
    RooDataSet* mcZRecSet3 = (RooDataSet*)mcZSet->reduce(Cut(cutsRec3));
    RooDataSet* mcZRecSet4 = (RooDataSet*)mcZSet->reduce(Cut(cutsRec4));
    RooDataSet* mcZRecSet5 = (RooDataSet*)mcZSet->reduce(Cut(cutsRec5));

    ///Wtau
    ///Fill weighted d.s. with wt=N(geom,PS,Zveto,Isolation)/N(gen) taken from Wmunu MC
    //RooDataSet* mcTauRecSet0 = fillHIWTauDataSet(baseString,fileNameMCTauIn+".root",muonTauArgSet,0.619);
    //RooDataSet* mcTauRecSet0 = &RooDataSet(*mcTauGeomFidSet);
    RooDataSet* mcTauRecSet0 = mcTauGeomFidSet;
    //mcTauRecSet0 = weightDS(mcTauRecSet0,w); 
    mcTauRecSet0->Print();
    ///Cut on the weighted d.s.
    RooDataSet* mcTauRecSet1 = (RooDataSet*)mcTauRecSet0->reduce(Cut(cutsTauMpt)); mcTauRecSet1->Print();
    RooDataSet* mcTauRecSet2 = (RooDataSet*)mcTauRecSet1->reduce(Cut(cutsTauMt)); mcTauRecSet2->Print();

    ///J1 
    ///Base reco dataset
    RooDataSet* mcJ1Set = fillHIMuonRecSet(baseString,fileNameMCJ1In+".root", muonRecArgSet); mcJ1Set->Print(); 
    ///Apply the cut flow
    RooDataSet* mcJ1RecSet0 = (RooDataSet*)mcJ1Set->reduce(Cut(cutsRec0));
    RooDataSet* mcJ1RecSet1 = (RooDataSet*)mcJ1Set->reduce(Cut(cutsRec1));
    RooDataSet* mcJ1RecSet2 = (RooDataSet*)mcJ1Set->reduce(Cut(cutsRec2));
    RooDataSet* mcJ1RecSet3 = (RooDataSet*)mcJ1Set->reduce(Cut(cutsRec3));
    RooDataSet* mcJ1RecSet4 = (RooDataSet*)mcJ1Set->reduce(Cut(cutsRec4));
    RooDataSet* mcJ1RecSet5 = (RooDataSet*)mcJ1Set->reduce(Cut(cutsRec5));

    ///J2 
    ///Base reco dataset
    RooDataSet* mcJ2Set = fillHIMuonRecSet(baseString,fileNameMCJ2In+".root", muonRecArgSet); mcJ2Set->Print(); 
    ///Apply the cut flow
    RooDataSet* mcJ2RecSet0 = (RooDataSet*)mcJ2Set->reduce(Cut(cutsRec0));
    RooDataSet* mcJ2RecSet1 = (RooDataSet*)mcJ2Set->reduce(Cut(cutsRec1));
    RooDataSet* mcJ2RecSet2 = (RooDataSet*)mcJ2Set->reduce(Cut(cutsRec2));
    RooDataSet* mcJ2RecSet3 = (RooDataSet*)mcJ2Set->reduce(Cut(cutsRec3));
    RooDataSet* mcJ2RecSet4 = (RooDataSet*)mcJ2Set->reduce(Cut(cutsRec4));
    RooDataSet* mcJ2RecSet5 = (RooDataSet*)mcJ2Set->reduce(Cut(cutsRec5));

    ///J3 
    ///Base reco dataset
    RooDataSet* mcJ3Set = fillHIMuonRecSet(baseString,fileNameMCJ3In+".root", muonRecArgSet); mcJ3Set->Print(); 
    ///Apply the cut flow
    RooDataSet* mcJ3RecSet0 = (RooDataSet*)mcJ3Set->reduce(Cut(cutsRec0));
    RooDataSet* mcJ3RecSet1 = (RooDataSet*)mcJ3Set->reduce(Cut(cutsRec1));
    RooDataSet* mcJ3RecSet2 = (RooDataSet*)mcJ3Set->reduce(Cut(cutsRec2));
    RooDataSet* mcJ3RecSet3 = (RooDataSet*)mcJ3Set->reduce(Cut(cutsRec3));
    RooDataSet* mcJ3RecSet4 = (RooDataSet*)mcJ3Set->reduce(Cut(cutsRec4));
    RooDataSet* mcJ3RecSet5 = (RooDataSet*)mcJ3Set->reduce(Cut(cutsRec5));

    ///Create the histograms
    TH1F* hMcWGenSet = (TH1F*)mcWGenSet->createHistogram("hMcWGenSet",muonGenPt,Binning(b));
    TH1F* hMcZGenSet = (TH1F*)mcZGenSet->createHistogram("hMcZGenSet",muonGenPt,Binning(b));
    TH1F* hMcTauGenSet = (TH1F*)mcTauGenSet->createHistogram("hMcTauGenSet",muonPt,Binning(b));
    ///Make geometric fiducial generator level histos 
	TH1F* hMcWGeomFidSet = (TH1F*)mcWGeomFidSet->createHistogram("hMcWGeomFidSet",muonGenPt,Binning(b));
	TH1F* hMcZGeomFidSet = (TH1F*)mcZGeomFidSet->createHistogram("hMcZGeomFidSet",muonGenPt,Binning(b));
	TH1F* hMcTauGeomFidSet = (TH1F*)mcTauGeomFidSet->createHistogram("hMcTauGeomFidSet",muonPt,Binning(b));

    ///Fill datasets at reconstruction level
    ///W boson
    TH1F* hMcWRecSet0 = (TH1F*)mcWRecSet0->createHistogram("hMcWRecSet0",muonPt,Binning(b));
    TH1F* hMcWRecSet1 = (TH1F*)mcWRecSet1->createHistogram("hMcWRecSet1",muonPt,Binning(b));
    TH1F* hMcWRecSet2 = (TH1F*)mcWRecSet2->createHistogram("hMcWRecSet2",muonPt,Binning(b));
    TH1F* hMcWRecSet3 = (TH1F*)mcWRecSet3->createHistogram("hMcWRecSet3",muonPt,Binning(b));
    TH1F* hMcWRecSet4 = (TH1F*)mcWRecSet4->createHistogram("hMcWRecSet4",muonPt,Binning(b));
    TH1F* hMcWRecSet5 = (TH1F*)mcWRecSet5->createHistogram("hMcWRecSet5",muonPt,Binning(b));
   
    TList _effForTau;

    std::vector<double> etaBins;
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

    ///centrality
    std::vector <double> centralityBins;
    centralityBins.push_back(0.0);
    if(doCentrality){
        centralityBins.push_back(0.05);
        centralityBins.push_back(0.1);
        centralityBins.push_back(0.15);
        centralityBins.push_back(0.2);
        centralityBins.push_back(0.4);
        
    }
    centralityBins.push_back(0.8);

    const int nCentralityBins = centralityBins.size()-1;
    std::cout << "nCentralityBins : " << nCentralityBins << std::endl;

    RooDataSet* mcWGeomFidSubSet[nEtaBins][nCentralityBins],*mcWRecSubSet3[nEtaBins][nCentralityBins];
    RooDataSet* mcWRecSubSet5[nEtaBins][nCentralityBins];
    ///Bin in eta and centrality for tau study
    for(int icent=0; icent<nCentralityBins; ++icent){
        _effForTau.Add( new TGraphAsymmErrors(nEtaBins));
        for(int ieta=0; ieta<nEtaBins; ++ieta){
            
            ///Number of muons from Ws with pT>25GeV generated in
            ///this eta and centrality bin
            mcWGeomFidSubSet[ieta][icent] = selectPtEtaCentrality(mcWGeomFidSet,25.0,350.0,etaBins[ieta],etaBins[ieta+1],centralityBins[icent],centralityBins[icent+1],
                                                                    true,true);
            ///Number of muons from W with pT>25GeV reconstructed
            ///in this eta/centraliby bin and passed
            ///PS, isolation, and Z veto
            mcWRecSubSet3[ieta][icent] = 
                selectPtEtaCentrality(mcWRecSet3,25.0,350.0,etaBins[ieta],etaBins[ieta+1],centralityBins[icent],centralityBins[icent+1],true);
            ///use this ds for systematic study
            mcWRecSubSet5[ieta][icent] = 
                selectPtEtaCentrality(mcWRecSet5,25.0,350.0,etaBins[ieta],etaBins[ieta+1],centralityBins[icent],centralityBins[icent+1],true);

            int nGen = (int)mcWGeomFidSubSet[ieta][icent]->numEntries();
            std::cout << "Number of muons from Ws with pT>25GeV generated in eta: " << ieta << " centrality: " << icent << " = " << nGen << std::endl;
            int nPassed = (int)mcWRecSubSet3[ieta][icent]->numEntries();
            ///use for systematic study
            //std::cout << "WARNING! WARNING! WARNING!" << std::endl;
            //std::cout << "Using different cut efficiency for systematic study." << std::endl;
            //int nPassed = (int)mcWRecSubSet5[ieta][icent]->numEntries();
            std::cout << "Number of muons from W with pT>25GeV reconstructed that passed preselection, isolation, and Z veto in eta: " << 
                    ieta << " centrality: " << icent << " = " << nPassed << std::endl;
            std::cout << "That gives an efficiency in this bin of " << (double)nPassed/nGen << std::endl;
            std::cout << std::endl;
            ///Calculate efficiency using Bayesian errors 
            ///with a 68.3% C.I.
            double recEff0, errRec0Hi, errRec0Lo;
            Efficiency(nPassed,nGen,0.683,recEff0,errRec0Lo,errRec0Hi);
            errRec0Hi = errRec0Hi - recEff0;
            errRec0Lo = recEff0 - errRec0Lo;
            double binW = etaBins[ieta+1]-etaBins[ieta];
            double xEta = etaBins[ieta]+binW/2.0;
            ( (TGraphAsymmErrors*)_effForTau.At(icent))->SetPoint(ieta,xEta,recEff0); 
            ( (TGraphAsymmErrors*)_effForTau.At(icent))->SetPointError(ieta,binW/2.,binW/2.,errRec0Lo,errRec0Hi); 
        }//ieta
        ///Save TGraphs for each centrality bin
        TString sName = "effForTauStudy_cent"; sName+=icent;
        pFile->cd();
        ( (TGraphAsymmErrors*)_effForTau.At(icent))->Write(sName);
    }//icent
    
    for(int j=0; j<_effForTau.GetEntries(); ++j){
        delete _effForTau.At(j);
    }
    ///Z boson
    TH1F* hMcZRecSet0 = (TH1F*)mcZRecSet0->createHistogram("hMcZRecSet0",muonPt,Binning(b));
    TH1F* hMcZRecSet1 = (TH1F*)mcZRecSet1->createHistogram("hMcZRecSet1",muonPt,Binning(b));
    TH1F* hMcZRecSet2 = (TH1F*)mcZRecSet2->createHistogram("hMcZRecSet2",muonPt,Binning(b));
    TH1F* hMcZRecSet3 = (TH1F*)mcZRecSet3->createHistogram("hMcZRecSet3",muonPt,Binning(b));
    TH1F* hMcZRecSet4 = (TH1F*)mcZRecSet4->createHistogram("hMcZRecSet4",muonPt,Binning(b));
    TH1F* hMcZRecSet5 = (TH1F*)mcZRecSet5->createHistogram("hMcZRecSet5",muonPt,Binning(b));

    ///Wtau
    TH1F* hMcTauRecSet0 = (TH1F*)mcTauRecSet0->createHistogram("hMcTauRecSet0",muonPt,Binning(b));
    ///Create subset histos and scale by cut efficiency (N(up to Z veto)/N(in eta window))
//    TH1F* hMcTauRecSet0 = (TH1F*)hMcTauGeomFidSet->Clone("hMcTauRecSet0");
 //   hMcTauRecSet0->Scale(0.73681);
    TH1F* hMcTauRecSet1 = (TH1F*)mcTauRecSet1->createHistogram("hMcTauRecSet1",muonPt,Binning(b));
//    hMcTauRecSet1->Scale(0.73681);
    TH1F* hMcTauRecSet2 = (TH1F*)mcTauRecSet2->createHistogram("hMcTauRecSet2",muonPt,Binning(b));
//    hMcTauRecSet2->Scale(0.73681);

    ///J1 
    ///Apply the cut flow
    TH1F* hMcJ1RecSet0 = (TH1F*)mcJ1RecSet0->createHistogram("hMcJ1RecSet0",muonPt,Binning(b));
    TH1F* hMcJ1RecSet1 = (TH1F*)mcJ1RecSet1->createHistogram("hMcJ1RecSet1",muonPt,Binning(b));
    TH1F* hMcJ1RecSet2 = (TH1F*)mcJ1RecSet2->createHistogram("hMcJ1RecSet2",muonPt,Binning(b));
    TH1F* hMcJ1RecSet3 = (TH1F*)mcJ1RecSet3->createHistogram("hMcJ1RecSet3",muonPt,Binning(b));
    TH1F* hMcJ1RecSet4 = (TH1F*)mcJ1RecSet4->createHistogram("hMcJ1RecSet4",muonPt,Binning(b));
    TH1F* hMcJ1RecSet5 = (TH1F*)mcJ1RecSet5->createHistogram("hMcJ1RecSet5",muonPt,Binning(b));

    ///J2 
    TH1F* hMcJ2RecSet0 = (TH1F*)mcJ2RecSet0->createHistogram("hMcJ2RecSet0",muonPt,Binning(b));
    TH1F* hMcJ2RecSet1 = (TH1F*)mcJ2RecSet1->createHistogram("hMcJ2RecSet1",muonPt,Binning(b));
    TH1F* hMcJ2RecSet2 = (TH1F*)mcJ2RecSet2->createHistogram("hMcJ2RecSet2",muonPt,Binning(b));
    TH1F* hMcJ2RecSet3 = (TH1F*)mcJ2RecSet3->createHistogram("hMcJ2RecSet3",muonPt,Binning(b));
    TH1F* hMcJ2RecSet4 = (TH1F*)mcJ2RecSet4->createHistogram("hMcJ2RecSet4",muonPt,Binning(b));
    TH1F* hMcJ2RecSet5 = (TH1F*)mcJ2RecSet5->createHistogram("hMcJ2RecSet5",muonPt,Binning(b));

    ///J3 
    TH1F* hMcJ3RecSet0 = (TH1F*)mcJ3RecSet0->createHistogram("hMcJ3RecSet0",muonPt,Binning(b));
    TH1F* hMcJ3RecSet1 = (TH1F*)mcJ3RecSet1->createHistogram("hMcJ3RecSet1",muonPt,Binning(b));
    TH1F* hMcJ3RecSet2 = (TH1F*)mcJ3RecSet2->createHistogram("hMcJ3RecSet2",muonPt,Binning(b));
    TH1F* hMcJ3RecSet3 = (TH1F*)mcJ3RecSet3->createHistogram("hMcJ3RecSet3",muonPt,Binning(b));
    TH1F* hMcJ3RecSet4 = (TH1F*)mcJ3RecSet4->createHistogram("hMcJ3RecSet4",muonPt,Binning(b));
    TH1F* hMcJ3RecSet5 = (TH1F*)mcJ3RecSet5->createHistogram("hMcJ3RecSet5",muonPt,Binning(b));

    ///Weight the J1,2,3 histograms
    double ncoll = 452.0; //0-80%
    double centralityLow=0.0, centralityUpp=0.8;
	RooArgSet centralityArgSet(centrality);  
    ////Dataset with number of PbPb coll in 0-80% from each MC sple
    RooDataSet* mcJ1Events = fillHIEventDataSet(baseString,fileNameMCJ1In+".root",centralityArgSet );
    RooDataSet* mcJ2Events = fillHIEventDataSet(baseString,fileNameMCJ2In+".root",centralityArgSet );
    RooDataSet* mcJ3Events = fillHIEventDataSet(baseString,fileNameMCJ3In+".root",centralityArgSet );

	TH1F* hMcQCDRecSet0 = new TH1F("hMcQCDRecSet0","hMcQCDRecSet0",binnum, xBins);
    hMcQCDRecSet0 = getWeightedJxVariableBinnedHisto(xBins,hMcQCDRecSet0->GetNbinsX(),centralityLow, centralityUpp, ncoll,
                        hMcJ1RecSet0, hMcJ2RecSet0, hMcJ3RecSet0, 
                        mcJ1Events, mcJ2Events,mcJ3Events, hMcQCDRecSet0->GetNbinsX(), ptLow, ptmax);

	TH1F* hMcQCDRecSet1 = new TH1F("hMcQCDRecSet1","hMcQCDRecSet1",binnum, xBins);
    hMcQCDRecSet1 = getWeightedJxVariableBinnedHisto(xBins,hMcQCDRecSet1->GetNbinsX(),centralityLow, centralityUpp, ncoll,
                        hMcJ1RecSet1, hMcJ2RecSet1, hMcJ3RecSet1, 
                        mcJ1Events, mcJ2Events,mcJ3Events, hMcQCDRecSet1->GetNbinsX(), ptLow, ptmax);

	TH1F* hMcQCDRecSet2 = new TH1F("hMcQCDRecSet2","hMcQCDRecSet2",binnum, xBins);
    hMcQCDRecSet2 = getWeightedJxVariableBinnedHisto(xBins,hMcQCDRecSet2->GetNbinsX(),centralityLow, centralityUpp, ncoll,
                        hMcJ1RecSet2, hMcJ2RecSet2, hMcJ3RecSet2, 
                        mcJ1Events, mcJ2Events,mcJ3Events, hMcQCDRecSet2->GetNbinsX(), ptLow, ptmax);

	TH1F* hMcQCDRecSet3 = new TH1F("hMcQCDRecSet3","hMcQCDRecSet3",binnum, xBins);
    hMcQCDRecSet3 = getWeightedJxVariableBinnedHisto(xBins,hMcQCDRecSet3->GetNbinsX(),centralityLow, centralityUpp, ncoll,
                        hMcJ1RecSet3, hMcJ2RecSet3, hMcJ3RecSet3, 
                        mcJ1Events, mcJ2Events,mcJ3Events, hMcQCDRecSet3->GetNbinsX(), ptLow, ptmax);

	TH1F* hMcQCDRecSet4 = new TH1F("hMcQCDRecSet4","hMcQCDRecSet4",binnum, xBins);
    hMcQCDRecSet4 = getWeightedJxVariableBinnedHisto(xBins,hMcQCDRecSet4->GetNbinsX(),centralityLow, centralityUpp, ncoll,
                        hMcJ1RecSet4, hMcJ2RecSet4, hMcJ3RecSet4, 
                        mcJ1Events, mcJ2Events,mcJ3Events, hMcQCDRecSet4->GetNbinsX(), ptLow, ptmax);

	TH1F* hMcQCDRecSet5 = new TH1F("hMcQCDRecSet5","hMcQCDRecSet5",binnum, xBins);
    hMcQCDRecSet5 = getWeightedJxVariableBinnedHisto(xBins,hMcQCDRecSet5->GetNbinsX(),centralityLow, centralityUpp, ncoll,
                        hMcJ1RecSet5, hMcJ2RecSet5, hMcJ3RecSet5, 
                        mcJ1Events, mcJ2Events,mcJ3Events, hMcQCDRecSet5->GetNbinsX(), ptLow, ptmax);


    ///Scale by binwidth 
    ///W
    TH1F* hMcWGenSetc = (TH1F*)hMcWGenSet->Clone("hMcWGenSetc");
    hMcWGenSetc->Scale(1.0,"WIDTH");
    hMcWGenSetc->SetLineColor(kBlack);

	TH1F* hMcWGeomFidSetc = (TH1F*)hMcWGeomFidSet->Clone("hMcWGeomFidSetc");
    hMcWGeomFidSetc->Scale(1.0,"WIDTH");
    hMcWGeomFidSetc->SetLineColor(kRed);

    TH1F* hMcWRecSet0c = (TH1F*)hMcWRecSet0->Clone("hMcWRecSet0c");
    hMcWRecSet0c->Scale(1.0,"WIDTH");
    hMcWRecSet0c->SetLineColor(kMagenta);

    TH1F* hMcWRecSet1c = (TH1F*)hMcWRecSet1->Clone("hMcWRecSet1c");
    hMcWRecSet1c->Scale(1.0,"WIDTH");
    hMcWRecSet1c->SetLineColor(kGreen);

    TH1F* hMcWRecSet2c = (TH1F*)hMcWRecSet2->Clone("hMcWRecSet2c");
    hMcWRecSet2c->Scale(1.0,"WIDTH");
    hMcWRecSet2c->SetLineColor(kBlue);

    TH1F* hMcWRecSet3c = (TH1F*)hMcWRecSet3->Clone("hMcWRecSet3c");
    hMcWRecSet3c->Scale(1.0,"WIDTH");
    hMcWRecSet3c->SetLineColor(kYellow);

    TH1F* hMcWRecSet4c = (TH1F*)hMcWRecSet4->Clone("hMcWRecSet4c");
    hMcWRecSet4c->Scale(1.0,"WIDTH");
    hMcWRecSet4c->SetLineColor(kViolet+5);

    TH1F* hMcWRecSet5c = (TH1F*)hMcWRecSet5->Clone("hMcWRecSet5c");
    hMcWRecSet5c->Scale(1.0,"WIDTH");
    hMcWRecSet5c->SetLineColor(kCyan);

    ///Wtau
    TH1F* hMcTauGenSetc = (TH1F*)hMcTauGenSet->Clone("hMcTauGenSetc");
    hMcTauGenSetc->Scale(1.0,"WIDTH");
    hMcTauGenSetc->SetLineColor(kBlack);

	TH1F* hMcTauGeomFidSetc = (TH1F*)hMcTauGeomFidSet->Clone("hMcTauGeomFidSetc");
    hMcTauGeomFidSetc->Scale(1.0,"WIDTH");
    hMcTauGeomFidSetc->SetLineColor(kRed);

    TH1F* hMcTauRecSet0c = (TH1F*)hMcTauRecSet0->Clone("hMcTauRecSet0c");
    hMcTauRecSet0c->Scale(0.73681);
    hMcTauRecSet0c->Scale(1.0,"WIDTH");
    hMcTauRecSet0c->SetLineColor(kGreen);

    TH1F* hMcTauRecSet1c = (TH1F*)hMcTauRecSet1->Clone("hMcTauRecSet1c");
    hMcTauRecSet1c->Scale(0.73681);
    hMcTauRecSet1c->Scale(1.0,"WIDTH");
    hMcTauRecSet1c->SetLineColor(kViolet+5);

    TH1F* hMcTauRecSet2c = (TH1F*)hMcTauRecSet2->Clone("hMcTauRecSet2c");
    hMcTauRecSet2c->Scale(0.73681);
    hMcTauRecSet2c->Scale(1.0,"WIDTH");
    hMcTauRecSet2c->SetLineColor(kCyan);

    ///Z
    TH1F* hMcZGenSetc = (TH1F*)hMcZGenSet->Clone("hMcZGenSetc");
    hMcZGenSetc->Scale(1.0,"WIDTH");
    hMcZGenSetc->SetLineColor(kBlack);

	TH1F* hMcZGeomFidSetc = (TH1F*)hMcZGeomFidSet->Clone("hMcZGeomFidSetc");
    hMcZGeomFidSetc->Scale(1.0,"WIDTH");
    hMcZGeomFidSetc->SetLineColor(kRed);

    TH1F* hMcZRecSet0c = (TH1F*)hMcZRecSet0->Clone("hMcZRecSet0c");
    hMcZRecSet0c->Scale(1.0,"WIDTH");
    hMcZRecSet0c->SetLineColor(kMagenta);

    TH1F* hMcZRecSet1c = (TH1F*)hMcZRecSet1->Clone("hMcZRecSet1c");
    hMcZRecSet1c->Scale(1.0,"WIDTH");
    hMcZRecSet1c->SetLineColor(kGreen);

    TH1F* hMcZRecSet2c = (TH1F*)hMcZRecSet2->Clone("hMcZRecSet2c");
    hMcZRecSet2c->Scale(1.0,"WIDTH");
    hMcZRecSet2c->SetLineColor(kBlue);

    TH1F* hMcZRecSet3c = (TH1F*)hMcZRecSet3->Clone("hMcZRecSet3c");
    hMcZRecSet3c->Scale(1.0,"WIDTH");
    hMcZRecSet3c->SetLineColor(kYellow);

    TH1F* hMcZRecSet4c = (TH1F*)hMcZRecSet4->Clone("hMcZRecSet4c");
    hMcZRecSet4c->Scale(1.0,"WIDTH");
    hMcZRecSet4c->SetLineColor(kViolet+5);

    TH1F* hMcZRecSet5c = (TH1F*)hMcZRecSet5->Clone("hMcZRecSet5c");
    hMcZRecSet5c->Scale(1.0,"WIDTH");
    hMcZRecSet5c->SetLineColor(kCyan);

    ///QCD
    TH1F* hMcQCDRecSet0c = (TH1F*)hMcQCDRecSet0->Clone("hMcQCDRecSet0c");
    hMcQCDRecSet0c->Scale(1.0,"WIDTH");
    hMcQCDRecSet0c->SetLineColor(kMagenta);

    TH1F* hMcQCDRecSet1c = (TH1F*)hMcQCDRecSet1->Clone("hMcQCDRecSet1c");
    hMcQCDRecSet1c->Scale(1.0,"WIDTH");
    hMcQCDRecSet1c->SetLineColor(kGreen);

    TH1F* hMcQCDRecSet2c = (TH1F*)hMcQCDRecSet2->Clone("hMcQCDRecSet2c");
    hMcQCDRecSet2c->Scale(1.0,"WIDTH");
    hMcQCDRecSet2c->SetLineColor(kBlue);

    TH1F* hMcQCDRecSet3c = (TH1F*)hMcQCDRecSet3->Clone("hMcQCDRecSet3c");
    hMcQCDRecSet3c->Scale(1.0,"WIDTH");
    hMcQCDRecSet3c->SetLineColor(kYellow);

    TH1F* hMcQCDRecSet4c = (TH1F*)hMcQCDRecSet4->Clone("hMcQCDRecSet4c");
    hMcQCDRecSet4c->Scale(1.0,"WIDTH");
    hMcQCDRecSet4c->SetLineColor(kViolet+5);

    TH1F* hMcQCDRecSet5c = (TH1F*)hMcQCDRecSet5->Clone("hMcQCDRecSet5c");
    hMcQCDRecSet5c->Scale(1.0,"WIDTH");
    hMcQCDRecSet5c->SetLineColor(kCyan);


    ///Take the ratio relative to generated muons
    std::cout << "Calculating W ratios.." << std::endl;
	TGraphAsymmErrors* pEffW0 =calculateEfficiency(hMcWGeomFidSetc,hMcWGenSetc,binnum);
	TGraphAsymmErrors* pEffW1 =calculateEfficiency(hMcWRecSet0c,hMcWGenSetc,binnum);
	TGraphAsymmErrors* pEffW2 =calculateEfficiency(hMcWRecSet1c,hMcWGenSetc,binnum);
	TGraphAsymmErrors* pEffW3 =calculateEfficiency(hMcWRecSet2c,hMcWGenSetc,binnum);
	TGraphAsymmErrors* pEffW4 =calculateEfficiency(hMcWRecSet3c,hMcWGenSetc,binnum);
	TGraphAsymmErrors* pEffW5 =calculateEfficiency(hMcWRecSet4c,hMcWGenSetc,binnum);
	TGraphAsymmErrors* pEffW6 =calculateEfficiency(hMcWRecSet5c,hMcWGenSetc,binnum);

    std::cout << "Calculating acc+cut W ratio for tau study..." << std::endl;
	TGraphAsymmErrors* pEffWAccCuts =calculateEfficiency(hMcWRecSet3c,hMcWGeomFidSetc,binnum);
    std::cout << "Writing acc+cut W ratio for tau study..." << std::endl;
    Write(pFile,pEffWAccCuts,"pEffWCutsMt");

    ///Take the ratio relative to generated muons
    std::cout << "Calculating Tau ratios.." << std::endl;
	TGraphAsymmErrors* pEffTau0 =calculateEfficiency(hMcTauGeomFidSetc,hMcTauGenSetc,binnum);
	TGraphAsymmErrors* pEffTau1 =calculateEfficiency(hMcTauRecSet0c,hMcTauGenSetc,binnum);
	TGraphAsymmErrors* pEffTau2 =calculateEfficiency(hMcTauRecSet1c,hMcTauGenSetc,binnum);
	TGraphAsymmErrors* pEffTau3 =calculateEfficiency(hMcTauRecSet2c,hMcTauGenSetc,binnum);

    std::cout << "Calculating Z ratios.." << std::endl;
	TGraphAsymmErrors* pEffZ0 =calculateEfficiency(hMcZGeomFidSetc,hMcZGenSetc,binnum);
	TGraphAsymmErrors* pEffZ1 =calculateEfficiency(hMcZRecSet0c,hMcZGenSetc,binnum);
	TGraphAsymmErrors* pEffZ2 =calculateEfficiency(hMcZRecSet1c,hMcZGenSetc,binnum);
	TGraphAsymmErrors* pEffZ3 =calculateEfficiency(hMcZRecSet2c,hMcZGenSetc,binnum);
	TGraphAsymmErrors* pEffZ4 =calculateEfficiency(hMcZRecSet3c,hMcZGenSetc,binnum);
	TGraphAsymmErrors* pEffZ5 =calculateEfficiency(hMcZRecSet4c,hMcZGenSetc,binnum);
	TGraphAsymmErrors* pEffZ6 =calculateEfficiency(hMcZRecSet5c,hMcZGenSetc,binnum);

    std::cout << "Calculating QCD ratios.." << std::endl;
	TGraphAsymmErrors* pEffQCD0 =calculateEfficiency(hMcQCDRecSet1c,hMcQCDRecSet0c,binnum);
	TGraphAsymmErrors* pEffQCD1 =calculateEfficiency(hMcQCDRecSet2c,hMcQCDRecSet0c,binnum);
	TGraphAsymmErrors* pEffQCD2 =calculateEfficiency(hMcQCDRecSet3c,hMcQCDRecSet0c,binnum);
	TGraphAsymmErrors* pEffQCD3 =calculateEfficiency(hMcQCDRecSet4c,hMcQCDRecSet0c,binnum);
	TGraphAsymmErrors* pEffQCD4 =calculateEfficiency(hMcQCDRecSet5c,hMcQCDRecSet0c,binnum);

    ///Decorations
	TCanvas *cptW = new TCanvas("cptW","cptW",600,700);
	//TPad *pad1 = new TPad("pad1","pad1",0,0.393,0.997,0.982);
	//TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
	TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,0.7);
	pad1->SetBottomMargin(0);
	pad1->Draw();
	pad1->cd();
	hMcWGenSetc->GetXaxis()->SetTitle("p_{T}^{#mu,truth}[GeV]"); 
	hMcWGenSetc->GetYaxis()->SetRangeUser(0.15,1.0e6); 
	hMcWGenSetc->GetXaxis()->SetRangeUser(0.0,140.0); 
	hMcWGenSetc->Draw("hist");
	hMcWGeomFidSetc->Draw("hist same");
	hMcWRecSet0c->Draw("hist same");
	hMcWRecSet1c->Draw("hist same");
	hMcWRecSet2c->Draw("hist same");
	hMcWRecSet3c->Draw("hist same");
	hMcWRecSet4c->Draw("hist same");
	hMcWRecSet5c->Draw("hist same");
	hMcWGenSetc->Draw("sameaxis");
	cptW->SetLogy(); cptW->Update();
	
	TLegend* leg0 = new TLegend(0.619, 0.559, 0.876, 0.939);
	leg0->SetTextFont(gStyle->GetTextFont());
	leg0->SetTextSize(0.04);
	leg0->SetBorderSize(0);
	leg0->SetFillColor(0);
	leg0->AddEntry("hMcWGenSetc","MC truth","l");
	leg0->AddEntry("hMcWGeomFidSetc","0.1<|#eta_{gen}|<2.4","l");
	leg0->AddEntry("hMcWRecSet0c","rec #mu","l");
	leg0->AddEntry("hMcWRecSet1c","Pre-Selected","l");
	leg0->AddEntry("hMcWRecSet2c","i_{#mu}<0.1","l");
	leg0->AddEntry("hMcWRecSet3c","m_{#mu#mu}<66 GeV","l");
//	leg0->AddEntry("hMcWRecSet3c","i_{#mu}<0.1","l");
	leg0->AddEntry("hMcWRecSet4c","#slash{p_{T}}>25GeV","l");
	leg0->AddEntry("hMcWRecSet5c","m_{T}>40GeV","l");
	leg0->Draw();

	double ptCutLine = 25.0;
    TLine *line0 = new TLine(ptCutLine,0.01,ptCutLine,881201);
    line0->SetLineColor(kBlack);
    line0->SetLineStyle(kDashed);
    line0->SetLineWidth(2);
    line0->Draw();

    TString nameW = "W#rightarrow#mu#nu";
	myText(0.196,0.862,kBlack,(char*)"W#rightarrow#mu#nu");
	
	cptW->cd();
    cptW->SetLogy();
    cptW->Update();
	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3958);
	pad2->SetTopMargin(0);
	pad2->Draw();
	pad2->cd();

	TH1F* hDummy = new TH1F("hDummy","hDummy",binnum, xBins);
    hDummy->GetYaxis()->SetRangeUser(0.0,1.0);

	pEffW0->SetMarkerColor(kRed);
	pEffW1->SetMarkerColor(kMagenta);
	pEffW2->SetMarkerColor(kGreen);
	pEffW3->SetMarkerColor(kBlue);
	pEffW4->SetMarkerColor(kYellow);
	pEffW5->SetMarkerColor(kViolet+5);
	pEffW6->SetMarkerColor(kCyan);

    hDummy->Draw();
	pEffW0->Draw("pesame");
	pEffW1->Draw("pesame");
	pEffW2->Draw("pesame");
	pEffW3->Draw("pesame");
	pEffW4->Draw("pesame");
	pEffW5->Draw("pesame");
	pEffW6->Draw("pesame");
//	pEffW0->Draw("sameaxis");

    TString label = "mcCutFlow";
    cptW->Update(); cptW->Print(label+"W_log.pdf");
	cptW->Print(label+"W_log.root");

    std::cout << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "W S u m m a r y" << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "Number of generated events                                     = " << hMcWGenSetc->Integral() << std::endl; 
    std::cout << "Number of events in 0.1<|eta|<2.4                              = " << hMcWGeomFidSetc->Integral() << std::endl; 
    std::cout << "Number of events reconstructed                                 = " << hMcWRecSet0c->Integral() << std::endl; 
    std::cout << "Number of preselected events                                   = " << hMcWRecSet1c->Integral() << std::endl; 
    std::cout << "Number of events after isolation                               = " << hMcWRecSet2c->Integral() << std::endl; 
    std::cout << "Number of events after Z veto                                  = " << hMcWRecSet3c->Integral() << std::endl; 
    std::cout << "Number of events at pt>25GeV                                   = " << hMcWRecSet3c->Integral(26,binnum) << std::endl; 
    std::cout << "Number of events at mpt>25GeV                                  = " << hMcWRecSet4c->Integral(26,binnum) << std::endl; 
    std::cout << "Number of events at mt>40GeV                                   = " << hMcWRecSet5c->Integral(26,binnum) << std::endl; 


	TCanvas *cptTau = new TCanvas("cptTau","cptTau",600,700);
	TPad *padTau1 = new TPad("padTau1","padTau1",0,0.3,1,1);
	padTau1->SetBottomMargin(0);
	padTau1->Draw();
	padTau1->cd();
	hMcTauGenSetc->GetXaxis()->SetTitle("p_{T}^{#mu,truth}[GeV]"); 
	hMcTauGenSetc->GetYaxis()->SetRangeUser(0.15,1.0e6); 
	hMcTauGenSetc->GetXaxis()->SetRangeUser(0.0,140.0); 
	hMcTauGenSetc->Draw("hist");
	hMcTauGeomFidSetc->Draw("hist same");
	hMcTauRecSet0c->Draw("hist same");
	hMcTauRecSet1c->Draw("hist same");
	hMcTauRecSet2c->Draw("hist same");
	hMcTauGenSetc->Draw("sameaxis");
	cptTau->SetLogy(); cptTau->Update();
	
	TLegend* legTau = new TLegend(0.619, 0.559, 0.876, 0.939);
	legTau->SetTextFont(gStyle->GetTextFont());
	legTau->SetTextSize(0.04);
	legTau->SetBorderSize(0);
	legTau->SetFillColor(0);
	legTau->AddEntry("hMcTauGenSetc","MC truth","l");
	legTau->AddEntry("hMcTauGeomFidSetc","0.1<|#eta_{gen}|<2.4","l");
	legTau->AddEntry("hMcTauRecSet0c","PS,i_{#mu}<0.1,m_{#mu#mu}<66 GeV","l");
	legTau->AddEntry("hMcTauRecSet1c","#slash{p_{T}}>25GeV","l");
	legTau->AddEntry("hMcTauRecSet2c","m_{T}>40GeV","l");
	legTau->Draw();

    line0->SetLineColor(kBlack);
    line0->SetLineStyle(kDashed);
    line0->SetLineWidth(2);
    line0->Draw();

    TString nameTau = "W#rightarrow#tau";
	myText(0.196,0.862,kBlack,(char*)"W#rightarrow#tau");
	
	cptTau->cd();
    cptTau->SetLogy();
    cptTau->Update();
	TPad *padTau2 = new TPad("padTau2","padTau2",0,0,1,0.3958);
	padTau2->SetTopMargin(0);
	padTau2->Draw();
	padTau2->cd();

	pEffTau0->SetMarkerColor(kRed);
	pEffTau1->SetMarkerColor(kGreen);
	pEffTau2->SetMarkerColor(kViolet+5);
	pEffTau3->SetMarkerColor(kCyan);

    hDummy->Draw();
    pEffTau0->Draw("pesame");
    pEffTau1->Draw("pesame");
    pEffTau2->Draw("pesame");
    pEffTau3->Draw("pesame");

    label = "mcCutFlow";
    cptTau->Update(); cptTau->Print(label+"Tau_log.pdf");
	cptTau->Print(label+"Tau_log.root");
    std::cout << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "Tau S u m m a r y" << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "Number of generated events                                     = " << hMcTauGenSetc->Integral() << std::endl; 
    std::cout << "Number of events in 0.1<|eta|<2.4                              = " << hMcTauGeomFidSetc->Integral() << std::endl; 
    std::cout << "Number of events after cut efficiency                          = " << hMcTauRecSet0c->Integral() << std::endl; 
    std::cout << "Number of events at pt>25GeV                                   = " << hMcTauRecSet0c->Integral(26,binnum) << std::endl; 
    std::cout << "Number of events at mpt>25GeV                                  = " << hMcTauRecSet1c->Integral(26,binnum) << std::endl; 
    std::cout << "Number of events at mt>40GeV                                   = " << hMcTauRecSet2c->Integral(26,binnum) << std::endl; 

	TCanvas *cptZ = new TCanvas("cptZ","cptZ",600,700);
	TPad *padZ1 = new TPad("padZ1","padZ1",0,0.3,1,1);
	padZ1->SetBottomMargin(0);
	padZ1->Draw();
	padZ1->cd();
	hMcZGenSetc->GetXaxis()->SetTitle("p_{T}^{#mu,truth}[GeV]"); 
	hMcZGenSetc->GetYaxis()->SetRangeUser(0.15,1.0e6); 
	hMcZGenSetc->GetXaxis()->SetRangeUser(0.0,140.0); 
	hMcZGenSetc->Draw("hist");
	hMcZGeomFidSetc->Draw("hist same");
	hMcZRecSet0c->Draw("hist same");
	hMcZRecSet1c->Draw("hist same");
	hMcZRecSet2c->Draw("hist same");
	hMcZRecSet3c->Draw("hist same");
	hMcZRecSet4c->Draw("hist same");
	hMcZRecSet5c->Draw("hist same");
	hMcZGenSetc->Draw("sameaxis");
	cptZ->SetLogy(); cptZ->Update();
	
	TLegend* leg1 = new TLegend(0.619, 0.559, 0.876, 0.939);
	leg1->SetTextFont(gStyle->GetTextFont());
	leg1->SetTextSize(0.04);
	leg1->SetBorderSize(0);
	leg1->SetFillColor(0);
	leg1->AddEntry("hMcZGenSetc","MC truth","l");
	leg1->AddEntry("hMcZGeomFidSetc","0.1<|#eta_{gen}|<2.4","l");
	leg1->AddEntry("hMcZRecSet0c","rec #mu","l");
	leg1->AddEntry("hMcZRecSet1c","Pre-Selected","l");
	leg1->AddEntry("hMcZRecSet2c","i_{#mu}<0.1","l");
	leg1->AddEntry("hMcZRecSet3c","m_{#mu#mu}<66 GeV","l");
	leg1->AddEntry("hMcZRecSet4c","#slash{p_{T}}>25GeV","l");
	leg1->AddEntry("hMcZRecSet5c","m_{T}>40GeV","l");
	leg1->Draw();

    //TLine *line0 = new TLine(ptCutLine,0.01,ptCutLine,881201);
    line0->SetLineColor(kBlack);
    line0->SetLineStyle(kDashed);
    line0->SetLineWidth(2);
    line0->Draw();

    TString nameZ = "Z#rightarrow#mu#mu";
	myText(0.196,0.862,kBlack,(char*)"Z#rightarrow#mu#mu");
	
	cptZ->cd();
    cptZ->SetLogy();
    cptZ->Update();
	TPad *padZ2 = new TPad("padZ2","padZ2",0,0,1,0.3958);
	padZ2->SetTopMargin(0);
	padZ2->Draw();
	padZ2->cd();

	pEffZ0->SetMarkerColor(kRed);
	pEffZ1->SetMarkerColor(kMagenta);
	pEffZ2->SetMarkerColor(kGreen);
	pEffZ3->SetMarkerColor(kBlue);
	pEffZ4->SetMarkerColor(kYellow);
	pEffZ5->SetMarkerColor(kViolet+5);
	pEffZ6->SetMarkerColor(kCyan);

    hDummy->Draw();
	pEffZ0->Draw("pesame");
	pEffZ1->Draw("pesame");
	pEffZ2->Draw("pesame");
	pEffZ3->Draw("pesame");
	pEffZ4->Draw("pesame");
	pEffZ5->Draw("pesame");
	pEffZ6->Draw("pesame");
//	pEffZ0->Draw("sameaxis");

    label = "mcCutFlow";
    cptZ->Update(); cptZ->Print(label+"Z_log.pdf");
	cptZ->Print(label+"Z_log.root");
    std::cout << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "Z S u m m a r y" << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "Number of generated events                                     = " << hMcZGenSetc->Integral() << std::endl; 
    std::cout << "Number of events in 0.1<|eta|<2.4                              = " << hMcZGeomFidSetc->Integral() << std::endl; 
    std::cout << "Number of events reconstructed                                 = " << hMcZRecSet0c->Integral() << std::endl; 
    std::cout << "Number of preselected events                                   = " << hMcZRecSet1c->Integral() << std::endl; 
    std::cout << "Number of events after isolation                               = " << hMcZRecSet2c->Integral() << std::endl; 
    std::cout << "Number of events after Z veto                                  = " << hMcZRecSet3c->Integral() << std::endl; 
    std::cout << "Number of events at pt>25GeV                                   = " << hMcZRecSet3c->Integral(26,binnum) << std::endl; 
    std::cout << "Number of events at mpt>25GeV                                  = " << hMcZRecSet4c->Integral(26,binnum) << std::endl; 
    std::cout << "Number of events at mt>40GeV                                   = " << hMcZRecSet5c->Integral(26,binnum) << std::endl; 


	TCanvas *cptQCD = new TCanvas("cptQCD","cptQCD",600,700);
	TPad *padQCD1 = new TPad("padQCD1","padQCD1",0,0.3,1,1);
	padQCD1->SetBottomMargin(0);
	padQCD1->Draw();
	padQCD1->cd();
	hMcQCDRecSet0c->GetXaxis()->SetTitle("p_{T}^{#mu,truth}[GeV]"); 
	hMcQCDRecSet0c->GetYaxis()->SetRangeUser(0.15,1.0e6); 
	hMcQCDRecSet0c->GetXaxis()->SetRangeUser(0.0,140.0); 
	hMcQCDRecSet0c->Draw("hist");
	hMcQCDRecSet1c->Draw("hist same");
	hMcQCDRecSet2c->Draw("hist same");
	hMcQCDRecSet3c->Draw("hist same");
	hMcQCDRecSet4c->Draw("hist same");
	hMcQCDRecSet5c->Draw("hist same");
	hMcQCDRecSet0c->Draw("sameaxis");
	cptQCD->SetLogy(); cptQCD->Update();
	
	leg1 = new TLegend(0.619, 0.559, 0.876, 0.939);
	leg1->SetTextFont(gStyle->GetTextFont());
	leg1->SetTextSize(0.04);
	leg1->SetBorderSize(0);
	leg1->SetFillColor(0);
	leg1->AddEntry("hMcQCDRecSet0c","rec #mu","l");
	leg1->AddEntry("hMcQCDRecSet1c","Pre-Selected","l");
	leg1->AddEntry("hMcQCDRecSet2c","i_{#mu}<0.1","l");
	leg1->AddEntry("hMcQCDRecSet3c","m_{#mu#mu}<66 GeV","l");
	leg1->AddEntry("hMcQCDRecSet4c","#slash{p_{T}}>25GeV","l");
	leg1->AddEntry("hMcQCDRecSet5c","m_{T}>40GeV","l");
	leg1->Draw();

    //TLine *line0 = new TLine(ptCutLine,0.01,ptCutLine,881201);
    line0->SetLineColor(kBlack);
    line0->SetLineStyle(kDashed);
    line0->SetLineWidth(2);
    line0->Draw();

    TString nameQCD = "Jx #mu";
	myText(0.196,0.862,kBlack,(char*)"Jx #mu");
	
	cptQCD->cd();
    cptQCD->SetLogy();
    cptQCD->Update();
	TPad *padQCD2 = new TPad("padQCD2","padQCD2",0,0,1,0.3958);
	padQCD2->SetTopMargin(0);
	padQCD2->Draw();
	padQCD2->cd();
	pEffQCD0->SetMarkerColor(kGreen);
	pEffQCD1->SetMarkerColor(kBlue);
	pEffQCD2->SetMarkerColor(kYellow);
	pEffQCD3->SetMarkerColor(kViolet+5);
	pEffQCD4->SetMarkerColor(kCyan);

    hDummy->Draw();
	pEffQCD0->Draw("pesame");
	pEffQCD1->Draw("pesame");
	pEffQCD2->Draw("pesame");
	pEffQCD3->Draw("pesame");
	pEffQCD4->Draw("pesame");
//	pEffQCD0->Draw("sameaxis");

    label = "mcCutFlow";
    cptQCD->Update(); cptQCD->Print(label+"QCD_log.pdf");
	cptQCD->Print(label+"QCD_log.root");
    std::cout << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << "Q C D  S u m m a r y" << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << "Number of reconstructed events                                 = " << hMcQCDRecSet0c->Integral() << std::endl; 
    std::cout << "Number of preselected events                                   = " << hMcQCDRecSet1c->Integral() << std::endl; 
    std::cout << "Number of events after isolation                               = " << hMcQCDRecSet2c->Integral() << std::endl; 
    std::cout << "Number of events after Z veto                                  = " << hMcQCDRecSet3c->Integral() << std::endl; 
    std::cout << "Number of events at pt>25GeV                                   = " << hMcQCDRecSet3c->Integral(26,binnum) << std::endl; 
    std::cout << "Number of events at mpt>25GeV                                  = " << hMcQCDRecSet4c->Integral(26,binnum) << std::endl; 
    std::cout << "Number of events at mt>40GeV                                   = " << hMcQCDRecSet5c->Integral(26,binnum) << std::endl; 

}
