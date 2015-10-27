#include "TCanvas.h"
#include "TDirectory.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2DErrors.h"
#include "TList.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TChain.h"
#include "THStack.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLatex.h"
#include "TClonesArray.h"
#include "TPad.h"

#include "RooDataSet.h"
#include "RooBinning.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHist.h"

#include "WAnalysisHIDep.C"
#include "AtlasUtils.C"
#include "AtlasStyle.C"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath> 
#include <sstream>

#include "EfficiencyCorrection.C"
#include "Systematics.C"
#include "correctionFactorsDep.C"
#include "WPlotterHelper.C"

void Write(TFile* outFile, TObject* grIso, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing TGraph to root file..." << std::endl;
	  grIso->Write(sFile);
	  std::cout << "Done." << std::endl;
   	  gDirectory = dir;
}

///////////////////////////////////////////////////////////////////////////////
// fillDataSet
///////////////////////////////////////////////////////////////////////////////
RooDataSet* fillDataSet(const TString& pathName, const TString& fileName, RooArgSet& muonArgSet, int cutValue = 11)
{
  RooDataSet* set = new RooDataSet(fileName,fileName,muonArgSet);
  
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
  float ptcone10ID05Nt[50],ptcone15ID05Nt[50],ptcone20ID05Nt[50],ptcone30ID05Nt[50],ptcone40ID05Nt[50],ptcone50ID05Nt[50];
  float ptcone10ID075Nt[50],ptcone15ID075Nt[50],ptcone20ID075Nt[50],ptcone30ID075Nt[50],ptcone40ID075Nt[50],ptcone50ID075Nt[50];
  float ptcone10ID1Nt[50],ptcone15ID1Nt[50],ptcone20ID1Nt[50],ptcone30ID1Nt[50],ptcone40ID1Nt[50],ptcone50ID1Nt[50];
  float ptcone10ID2Nt[50],ptcone15ID2Nt[50],ptcone20ID2Nt[50],ptcone30ID2Nt[50],ptcone40ID2Nt[50],ptcone50ID2Nt[50];
  float ptcone10ID3Nt[50],ptcone15ID3Nt[50],ptcone20ID3Nt[50],ptcone30ID3Nt[50],ptcone40ID3Nt[50],ptcone50ID3Nt[50];
  float ptcone10ID4Nt[50],ptcone15ID4Nt[50],ptcone20ID4Nt[50],ptcone30ID4Nt[50],ptcone40ID4Nt[50],ptcone50ID4Nt[50];
  float ptcone10ID5Nt[50],ptcone15ID5Nt[50],ptcone20ID5Nt[50],ptcone30ID5Nt[50],ptcone40ID5Nt[50],ptcone50ID5Nt[50];
  float ptcone10ID6Nt[50],ptcone15ID6Nt[50],ptcone20ID6Nt[50],ptcone30ID6Nt[50],ptcone40ID6Nt[50],ptcone50ID6Nt[50];
  int valNt[50], ZDYNt[50];
  int nmu;
  

  // --- Load the MC tree ---
  TChain* tree = new TChain("tree","tree");
  int nFiles = tree->Add(pathName+fileName);
  
  std::cout << "Filling the RooDataSet for " << fileName << "... Number of files: " << nFiles << std::endl;
 
  // --- Set branch adresses ---
  tree->SetBranchAddress("mu_muid_n", &nmu);
  tree->SetBranchAddress("eLoss", &eLossNt);
  tree->SetBranchAddress("scat", &scatNt);
  tree->SetBranchAddress("comp", &compNt);
  tree->SetBranchAddress("pt", &ptNt);
  tree->SetBranchAddress("mt", &mtNt);
  tree->SetBranchAddress("eta", &etaNt);
  tree->SetBranchAddress("phi", &phiNt);
  tree->SetBranchAddress("charge", &chargeNt);
  tree->SetBranchAddress("prompt", &promptNt);
  tree->SetBranchAddress("val", &valNt); 
  tree->SetBranchAddress("ZDY", &ZDYNt); 
  tree->SetBranchAddress("centrality", &centralityNt);
  tree->SetBranchAddress("nu_pt", &nu_ptNt);

  tree->SetBranchAddress("ptcone10ID05", &ptcone10ID05Nt);
  tree->SetBranchAddress("ptcone15ID05", &ptcone15ID05Nt);
  tree->SetBranchAddress("ptcone20ID05", &ptcone20ID05Nt);
  tree->SetBranchAddress("ptcone30ID05", &ptcone30ID05Nt);
  tree->SetBranchAddress("ptcone40ID05", &ptcone40ID05Nt);
  tree->SetBranchAddress("ptcone50ID05", &ptcone50ID05Nt);

  tree->SetBranchAddress("ptcone10ID075", &ptcone10ID075Nt);
  tree->SetBranchAddress("ptcone15ID075", &ptcone15ID075Nt);
  tree->SetBranchAddress("ptcone20ID075", &ptcone20ID075Nt);
  tree->SetBranchAddress("ptcone30ID075", &ptcone30ID075Nt);
  tree->SetBranchAddress("ptcone40ID075", &ptcone40ID075Nt);
  tree->SetBranchAddress("ptcone50ID075", &ptcone50ID075Nt);
  tree->SetBranchAddress("ptcone10ID1", &ptcone10ID1Nt);
  tree->SetBranchAddress("ptcone15ID1", &ptcone15ID1Nt);
  tree->SetBranchAddress("ptcone20ID1", &ptcone20ID1Nt);
  tree->SetBranchAddress("ptcone30ID1", &ptcone30ID1Nt);
  tree->SetBranchAddress("ptcone40ID1", &ptcone40ID1Nt);
  tree->SetBranchAddress("ptcone50ID1", &ptcone50ID1Nt);

  tree->SetBranchAddress("ptcone10ID2", &ptcone10ID2Nt);
  tree->SetBranchAddress("ptcone15ID2", &ptcone15ID2Nt);
  tree->SetBranchAddress("ptcone20ID2", &ptcone20ID2Nt);
  tree->SetBranchAddress("ptcone30ID2", &ptcone30ID2Nt);
  tree->SetBranchAddress("ptcone40ID2", &ptcone40ID2Nt);
  tree->SetBranchAddress("ptcone50ID2", &ptcone50ID2Nt);

  tree->SetBranchAddress("ptcone10ID3", &ptcone10ID3Nt);
  tree->SetBranchAddress("ptcone15ID3", &ptcone15ID3Nt);
  tree->SetBranchAddress("ptcone20ID3", &ptcone20ID3Nt);
  tree->SetBranchAddress("ptcone30ID3", &ptcone30ID3Nt);
  tree->SetBranchAddress("ptcone40ID3", &ptcone40ID3Nt);
  tree->SetBranchAddress("ptcone50ID3", &ptcone50ID3Nt);

  tree->SetBranchAddress("ptcone10ID4", &ptcone10ID4Nt);
  tree->SetBranchAddress("ptcone15ID4", &ptcone15ID4Nt);
  tree->SetBranchAddress("ptcone20ID4", &ptcone20ID4Nt);
  tree->SetBranchAddress("ptcone30ID4", &ptcone30ID4Nt);
  tree->SetBranchAddress("ptcone40ID4", &ptcone40ID4Nt);
  tree->SetBranchAddress("ptcone50ID4", &ptcone50ID4Nt);

  tree->SetBranchAddress("ptcone10ID5", &ptcone10ID5Nt);
  tree->SetBranchAddress("ptcone15ID5", &ptcone15ID5Nt);
  tree->SetBranchAddress("ptcone20ID5", &ptcone20ID5Nt);
  tree->SetBranchAddress("ptcone30ID5", &ptcone30ID5Nt);
  tree->SetBranchAddress("ptcone40ID5", &ptcone40ID5Nt);
  tree->SetBranchAddress("ptcone50ID5", &ptcone50ID5Nt);

  tree->SetBranchAddress("ptcone10ID6", &ptcone10ID6Nt);
  tree->SetBranchAddress("ptcone15ID6", &ptcone15ID6Nt);
  tree->SetBranchAddress("ptcone20ID6", &ptcone20ID6Nt);
  tree->SetBranchAddress("ptcone30ID6", &ptcone30ID6Nt);
  tree->SetBranchAddress("ptcone40ID6", &ptcone40ID6Nt);
  tree->SetBranchAddress("ptcone50ID6", &ptcone50ID6Nt);


   // --- Set branch status ---
  tree->SetBranchStatus("*",0) ;
  tree->SetBranchStatus("mu_muid_n", 1);
  tree->SetBranchStatus("eLoss", 1);
  tree->SetBranchStatus("scat", 1);
  tree->SetBranchStatus("comp", 1);
  tree->SetBranchStatus("pt", 1);
  tree->SetBranchStatus("mt", 1);
  tree->SetBranchStatus("eta", 1);
  tree->SetBranchStatus("phi", 1);
  tree->SetBranchStatus("charge", 1);
  tree->SetBranchStatus("prompt", 1);
  tree->SetBranchStatus("val", 1); 
  tree->SetBranchStatus("ZDY", 1); 
  tree->SetBranchStatus("centrality", 1);
  tree->SetBranchStatus("nu_pt", 1);
  
  tree->SetBranchStatus("ptcone10ID05", 1);
  tree->SetBranchStatus("ptcone15ID05", 1);
  tree->SetBranchStatus("ptcone20ID05", 1);
  tree->SetBranchStatus("ptcone30ID05", 1);
  tree->SetBranchStatus("ptcone40ID05", 1);
  tree->SetBranchStatus("ptcone50ID05", 1);

  tree->SetBranchStatus("ptcone10ID075", 1);
  tree->SetBranchStatus("ptcone15ID075", 1);
  tree->SetBranchStatus("ptcone20ID075", 1);
  tree->SetBranchStatus("ptcone30ID075", 1);
  tree->SetBranchStatus("ptcone40ID075", 1);
  tree->SetBranchStatus("ptcone50ID075", 1);
  tree->SetBranchStatus("ptcone10ID1", 1);
  tree->SetBranchStatus("ptcone15ID1", 1);
  tree->SetBranchStatus("ptcone20ID1", 1);
  tree->SetBranchStatus("ptcone30ID1", 1);
  tree->SetBranchStatus("ptcone40ID1", 1);
  tree->SetBranchStatus("ptcone50ID1", 1);

  tree->SetBranchStatus("ptcone10ID2", 1);
  tree->SetBranchStatus("ptcone15ID2", 1);
  tree->SetBranchStatus("ptcone20ID2", 1);
  tree->SetBranchStatus("ptcone30ID2", 1);
  tree->SetBranchStatus("ptcone40ID2", 1);
  tree->SetBranchStatus("ptcone50ID2", 1);

  tree->SetBranchStatus("ptcone10ID3", 1);
  tree->SetBranchStatus("ptcone15ID3", 1);
  tree->SetBranchStatus("ptcone20ID3", 1);
  tree->SetBranchStatus("ptcone30ID3", 1);
  tree->SetBranchStatus("ptcone40ID3", 1);
  tree->SetBranchStatus("ptcone50ID3", 1);

  tree->SetBranchStatus("ptcone10ID4", 1);
  tree->SetBranchStatus("ptcone15ID4", 1);
  tree->SetBranchStatus("ptcone20ID4", 1);
  tree->SetBranchStatus("ptcone30ID4", 1);
  tree->SetBranchStatus("ptcone40ID4", 1);
  tree->SetBranchStatus("ptcone50ID4", 1);

  tree->SetBranchStatus("ptcone10ID5", 1);
  tree->SetBranchStatus("ptcone15ID5", 1);
  tree->SetBranchStatus("ptcone20ID5", 1);
  tree->SetBranchStatus("ptcone30ID5", 1);
  tree->SetBranchStatus("ptcone40ID5", 1);
  tree->SetBranchStatus("ptcone50ID5", 1);

  tree->SetBranchStatus("ptcone10ID6", 1);
  tree->SetBranchStatus("ptcone15ID6", 1);
  tree->SetBranchStatus("ptcone20ID6", 1);
  tree->SetBranchStatus("ptcone30ID6", 1);
  tree->SetBranchStatus("ptcone40ID6", 1);
  tree->SetBranchStatus("ptcone50ID6", 1);


  
  std::cout << "Number of entries: " << tree->GetEntries() << std::endl; 
  for ( int i = 0; i < tree->GetEntries(); i++ ) {
     //if (i>10000) return set; // temp hack
    tree->LoadTree(i);
    tree->GetEntry(i);

    muonArgSet.setRealValue("missPt",nu_ptNt);
    muonArgSet.setRealValue("centrality",centralityNt);
    for (int imu = 0; imu<nmu;imu++){
    if(valNt[imu]>cutValue&&fabs(eLossNt[imu])<0.5&&fabs(scatNt[imu])<4.0){

        if ( promptNt[imu] == 0 ) muonArgSet.setCatLabel("muonCategory","Unknown");
        else if ( promptNt[imu] == 1 ) muonArgSet.setCatLabel("muonCategory","HeavyFlavour");
        else if ( promptNt[imu] == 2 ) muonArgSet.setCatLabel("muonCategory","Tau");
        else if ( promptNt[imu] == 3 ) muonArgSet.setCatLabel("muonCategory","LightResonance");
        else if ( promptNt[imu] == 4 ) muonArgSet.setCatLabel("muonCategory","PionCaloDecay");
        else if ( promptNt[imu] == 5 ) muonArgSet.setCatLabel("muonCategory","PionEarlyDecay");
        else if ( promptNt[imu] == 6 ) muonArgSet.setCatLabel("muonCategory","PionLateDecay");
        else if ( promptNt[imu] == 7 ) muonArgSet.setCatLabel("muonCategory","KaonCaloDecay");
        else if ( promptNt[imu] == 8 ) muonArgSet.setCatLabel("muonCategory","KaonEarlyDecay");
        else if ( promptNt[imu] == 9 ) muonArgSet.setCatLabel("muonCategory","KaonLateDecay");
        else if ( promptNt[imu] == 10) muonArgSet.setCatLabel("muonCategory","Fake");
        else if ( promptNt[imu] == 23) muonArgSet.setCatLabel("muonCategory","Z");
        else if ( promptNt[imu] == 24) muonArgSet.setCatLabel("muonCategory","W");

      muonArgSet.setRealValue("muonELoss",eLossNt[imu]);
      muonArgSet.setRealValue("muonScat",scatNt[imu]);
      muonArgSet.setRealValue("motherRec",promptNt[imu]);

      double isolationTemp = 9999.0;
      isolationTemp = ptcone10ID05Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID05Cone10",isolationTemp);
      isolationTemp = ptcone15ID05Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID05Cone15",isolationTemp);
      isolationTemp = ptcone20ID05Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID05Cone20",isolationTemp);
      isolationTemp = ptcone30ID05Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID05Cone30",isolationTemp);
      isolationTemp = ptcone40ID05Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID05Cone40",isolationTemp);
      isolationTemp = ptcone50ID05Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID05Cone50",isolationTemp);

      isolationTemp = ptcone10ID075Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID075Cone10",isolationTemp);
      isolationTemp = ptcone15ID075Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID075Cone15",isolationTemp);
      isolationTemp = ptcone20ID075Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID075Cone20",isolationTemp);
      isolationTemp = ptcone30ID075Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID075Cone30",isolationTemp);
      isolationTemp = ptcone40ID075Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID075Cone40",isolationTemp);
      isolationTemp = ptcone50ID075Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID075Cone50",isolationTemp);

      isolationTemp = ptcone10ID1Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID1Cone10",isolationTemp);
      isolationTemp = ptcone15ID1Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID1Cone15",isolationTemp);
      isolationTemp = ptcone20ID1Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID1Cone20",isolationTemp);
      isolationTemp = ptcone30ID1Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID1Cone30",isolationTemp);
      isolationTemp = ptcone40ID1Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID1Cone40",isolationTemp);
      isolationTemp = ptcone50ID1Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID1Cone50",isolationTemp);

      isolationTemp = ptcone10ID2Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID2Cone10",isolationTemp);
      isolationTemp = ptcone15ID2Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID2Cone15",isolationTemp);
      isolationTemp = ptcone20ID2Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID2Cone20",isolationTemp);
      isolationTemp = ptcone30ID2Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID2Cone30",isolationTemp);
      isolationTemp = ptcone40ID2Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID2Cone40",isolationTemp);
      isolationTemp = ptcone50ID2Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID2Cone50",isolationTemp);

      isolationTemp = ptcone10ID3Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID3Cone10",isolationTemp);
      isolationTemp = ptcone15ID3Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID3Cone15",isolationTemp);
      isolationTemp = ptcone20ID3Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID3Cone20",isolationTemp);
      isolationTemp = ptcone30ID3Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID3Cone30",isolationTemp);
      isolationTemp = ptcone40ID3Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID3Cone40",isolationTemp);
      isolationTemp = ptcone50ID3Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID3Cone50",isolationTemp);

      isolationTemp = ptcone10ID4Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID4Cone10",isolationTemp);
      isolationTemp = ptcone15ID4Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID4Cone15",isolationTemp);
      isolationTemp = ptcone20ID4Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID4Cone20",isolationTemp);
      isolationTemp = ptcone30ID4Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID4Cone30",isolationTemp);
      isolationTemp = ptcone40ID4Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID4Cone40",isolationTemp);
      isolationTemp = ptcone50ID4Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID4Cone50",isolationTemp);

      isolationTemp = ptcone10ID5Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID5Cone10",isolationTemp);
      isolationTemp = ptcone15ID5Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID5Cone15",isolationTemp);
      isolationTemp = ptcone20ID5Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID5Cone20",isolationTemp);
      isolationTemp = ptcone30ID5Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID5Cone30",isolationTemp);
      isolationTemp = ptcone40ID5Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID5Cone40",isolationTemp);
      isolationTemp = ptcone50ID5Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID5Cone50",isolationTemp);

      isolationTemp = ptcone10ID6Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID6Cone10",isolationTemp);
      isolationTemp = ptcone15ID6Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID6Cone15",isolationTemp);
      isolationTemp = ptcone20ID6Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID6Cone20",isolationTemp);
      isolationTemp = ptcone30ID6Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID6Cone30",isolationTemp);
      isolationTemp = ptcone40ID6Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID6Cone40",isolationTemp);
      isolationTemp = ptcone50ID6Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolationID6Cone50",isolationTemp);

      muonArgSet.setRealValue("muonPt",ptNt[imu]);
      muonArgSet.setRealValue("muonMt",mtNt[imu]);
      muonArgSet.setRealValue("muonEta",etaNt[imu]);
      //muonArgSet.setRealValue("muonPhi",phiNt[imu]);
      muonArgSet.setRealValue("ZDY",ZDYNt[imu]);
      muonArgSet.setRealValue("muonCharge",chargeNt[imu]);
      muonArgSet.setRealValue("muonCategory",promptNt[imu]);
      if ( chargeNt[imu] > 0 ) muonArgSet.setCatLabel("chargeCategory","muPlus");
      else if ( chargeNt[imu] < 0) muonArgSet.setCatLabel("chargeCategory","muMinus");
      set->add(muonArgSet);    
    }
   }
  }

  return set;
}

RooDataSet* selectIDConeIsol(RooDataSet* dataSet, std::string trkPt,
            int coneRadius, double isolCut ){

                TString sIsolCut = "isolationID"; 
                sIsolCut += trkPt; 
                sIsolCut +="Cone"; 
                sIsolCut += coneRadius; 
                sIsolCut += "<"; 
                sIsolCut += isolCut;  
                std::cout << trkPt << ":" << coneRadius << ":" <<
                    isolCut << '\n' << sIsolCut << std::endl; 
                return (RooDataSet*)dataSet->reduce(Cut(sIsolCut));

}


void isoEfficiency(){

    int nbins = 400;
	double ptmax = 400.0;
    double cutValue = 11.0;

	//define bins for signal region
    //integrate from pt ptSigLo to ptmax
    double ptSigLo = 4.0;
	int ptSigBinLo = (nbins/ptmax)*ptSigLo+1; int ptSigBinUpp = (nbins/ptmax)*ptmax+1;  

    SetAtlasStyle();
	TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
	
	/// --- Open output file ---
	TDirectory *dir = gDirectory;
	TFile *outFile = new TFile("isoEfficiencyGraphs.root","RECREATE");
	gDirectory = dir;

	//data overlay
	TString fileNameSigIn = "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.04.12.2013";
  	// --- W set ---

	TString fileNameMCJ1In = "MonteCarloFiles/QCD/HISingleMuonJ1_1muPYTHIADataOverlay.04.13.2013";
	TString fileNameMCJ2In = "MonteCarloFiles/QCD/HISingleMuonJ2_1muPYTHIADataOverlay.04.13.2013";
	TString fileNameMCJ3In = "MonteCarloFiles/QCD/HISingleMuonJ3_1muPYTHIADataOverlay.04.13.2013";

    // --- declare cut variables --- //
    RooRealVar  muonPt("muonPt","p_{T}",0.0,400.0,"GeV");
    RooRealVar  missPt("missPt","p_{T}^{miss}",0.0,400.0,"GeV");
    RooRealVar  muonMt("muonMt","m_{T}",0.0,400.0,"GeV");
    RooRealVar  muonCharge("muonCharge","charge",-2.0,+2.0);
    RooRealVar  centrality("centrality","centrality",0.,1.0);
    RooRealVar  muonEta("muonEta","muonEta",-3.0,+3.0);
//    RooRealVar  muonPhi("muonPhi","muonPhi",-3.5,+3.5);
    RooRealVar  motherRec("motherRec","motherRec",0.0,250.0);
    RooRealVar  ZDY("ZDY","ZDY",0.0,1.0);

    ///List of isolation variables with different dR and pT^trk thresholds
    RooRealVar  isolationID05Cone10("isolationID05Cone10","isolationID05Cone10",0.0,100.0);
    RooRealVar  isolationID05Cone15("isolationID05Cone15","isolationID05Cone15",0.0,100.0);
    RooRealVar  isolationID05Cone20("isolationID05Cone20","isolationID05Cone20",0.0,100.0);
    RooRealVar  isolationID05Cone30("isolationID05Cone30","isolationID05Cone30",0.0,100.0);
    RooRealVar  isolationID05Cone40("isolationID05Cone40","isolationID05Cone40",0.0,100.0);
    RooRealVar  isolationID05Cone50("isolationID05Cone50","isolationID05Cone50",0.0,100.0);

    RooRealVar  isolationID075Cone10("isolationID075Cone10","isolationID075Cone10",0.0,100.0);
    RooRealVar  isolationID075Cone15("isolationID075Cone15","isolationID075Cone15",0.0,100.0);
    RooRealVar  isolationID075Cone20("isolationID075Cone20","isolationID075Cone20",0.0,100.0);
    RooRealVar  isolationID075Cone30("isolationID075Cone30","isolationID075Cone30",0.0,100.0);
    RooRealVar  isolationID075Cone40("isolationID075Cone40","isolationID075Cone40",0.0,100.0);
    RooRealVar  isolationID075Cone50("isolationID075Cone50","isolationID075Cone50",0.0,100.0);

    RooRealVar  isolationID1Cone10("isolationID1Cone10","isolationID1Cone10",0.0,100.0);
    RooRealVar  isolationID1Cone15("isolationID1Cone15","isolationID1Cone15",0.0,100.0);
    RooRealVar  isolationID1Cone20("isolationID1Cone20","isolationID1Cone20",0.0,100.0);
    RooRealVar  isolationID1Cone30("isolationID1Cone30","isolationID1Cone30",0.0,100.0);
    RooRealVar  isolationID1Cone40("isolationID1Cone40","isolationID1Cone40",0.0,100.0);
    RooRealVar  isolationID1Cone50("isolationID1Cone50","isolationID1Cone50",0.0,100.0);

    RooRealVar  isolationID2Cone10("isolationID2Cone10","isolationID2Cone10",0.0,100.0);
    RooRealVar  isolationID2Cone15("isolationID2Cone15","isolationID2Cone15",0.0,100.0);
    RooRealVar  isolationID2Cone20("isolationID2Cone20","isolationID2Cone20",0.0,100.0);
    RooRealVar  isolationID2Cone30("isolationID2Cone30","isolationID2Cone30",0.0,100.0);
    RooRealVar  isolationID2Cone40("isolationID2Cone40","isolationID2Cone40",0.0,100.0);
    RooRealVar  isolationID2Cone50("isolationID2Cone50","isolationID2Cone50",0.0,100.0);

    RooRealVar  isolationID3Cone10("isolationID3Cone10","isolationID3Cone10",0.0,100.0);
    RooRealVar  isolationID3Cone15("isolationID3Cone15","isolationID3Cone15",0.0,100.0);
    RooRealVar  isolationID3Cone20("isolationID3Cone20","isolationID3Cone20",0.0,100.0);
    RooRealVar  isolationID3Cone30("isolationID3Cone30","isolationID3Cone30",0.0,100.0);
    RooRealVar  isolationID3Cone40("isolationID3Cone40","isolationID3Cone40",0.0,100.0);
    RooRealVar  isolationID3Cone50("isolationID3Cone50","isolationID3Cone50",0.0,100.0);

    RooRealVar  isolationID4Cone10("isolationID4Cone10","isolationID4Cone10",0.0,100.0);
    RooRealVar  isolationID4Cone15("isolationID4Cone15","isolationID4Cone15",0.0,100.0);
    RooRealVar  isolationID4Cone20("isolationID4Cone20","isolationID4Cone20",0.0,100.0);
    RooRealVar  isolationID4Cone30("isolationID4Cone30","isolationID4Cone30",0.0,100.0);
    RooRealVar  isolationID4Cone40("isolationID4Cone40","isolationID4Cone40",0.0,100.0);
    RooRealVar  isolationID4Cone50("isolationID4Cone50","isolationID4Cone50",0.0,100.0);

    RooRealVar  isolationID5Cone10("isolationID5Cone10","isolationID5Cone10",0.0,100.0);
    RooRealVar  isolationID5Cone15("isolationID5Cone15","isolationID5Cone15",0.0,100.0);
    RooRealVar  isolationID5Cone20("isolationID5Cone20","isolationID5Cone20",0.0,100.0);
    RooRealVar  isolationID5Cone30("isolationID5Cone30","isolationID5Cone30",0.0,100.0);
    RooRealVar  isolationID5Cone40("isolationID5Cone40","isolationID5Cone40",0.0,100.0);
    RooRealVar  isolationID5Cone50("isolationID5Cone50","isolationID5Cone50",0.0,100.0);

    RooRealVar  isolationID6Cone10("isolationID6Cone10","isolationID6Cone10",0.0,100.0);
    RooRealVar  isolationID6Cone15("isolationID6Cone15","isolationID6Cone15",0.0,100.0);
    RooRealVar  isolationID6Cone20("isolationID6Cone20","isolationID6Cone20",0.0,100.0);
    RooRealVar  isolationID6Cone30("isolationID6Cone30","isolationID6Cone30",0.0,100.0);
    RooRealVar  isolationID6Cone40("isolationID6Cone40","isolationID6Cone40",0.0,100.0);
    RooRealVar  isolationID6Cone50("isolationID6Cone50","isolationID6Cone50",0.0,100.0);

    RooCategory muonCategory("muonCategory","muonCategory");
    muonCategory.defineType("Unknown",0);
    muonCategory.defineType("HeavyFlavour",1);
    muonCategory.defineType("Tau",2);
    muonCategory.defineType("LightResonance",3); 
    muonCategory.defineType("PionCaloDecay",4);
    muonCategory.defineType("PionEarlyDecay",5);
    muonCategory.defineType("PionLateDecay",6);
    muonCategory.defineType("KaonCaloDecay",7);
    muonCategory.defineType("KaonEarlyDecay",8);
    muonCategory.defineType("KaonLateDecay",9);
    muonCategory.defineType("Fake",10);
    muonCategory.defineType("Z",23);
    muonCategory.defineType("W",24);

    RooCategory chargeCategory("chargeCategory","chargeCategory") ;
    chargeCategory.defineType("muMinus",-1) ;
    chargeCategory.defineType("muPlus",1) ;


	TString sCutsSig =
    ///Before W sel
    //"muonPt>4.0&&abs(muonEta)>0.1&&abs(muonEta)<2.4&&ZDY==0";
    ///After W sel (minus the isolation cut)
    "muonPt>25.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0&&muonMt<400.0&&ZDY==0&&abs(muonEta)>0.1&&abs(muonEta)<2.4";
	RooFormulaVar cutsSig("cutsSig", "cutsSig", sCutsSig, RooArgList(muonPt,missPt,muonMt,ZDY,muonEta));

  	RooArgSet muonArgSet(muonPt,missPt,muonMt,muonEta,centrality,ZDY,muonCategory,chargeCategory);
    muonArgSet.add(motherRec);
    muonArgSet.add(isolationID05Cone10);
    muonArgSet.add(isolationID05Cone15);
    muonArgSet.add(isolationID05Cone20);
    muonArgSet.add(isolationID05Cone30);
    muonArgSet.add(isolationID05Cone40);
    muonArgSet.add(isolationID05Cone50);
    muonArgSet.add(isolationID075Cone10);
    muonArgSet.add(isolationID075Cone15);
    muonArgSet.add(isolationID075Cone20);
    muonArgSet.add(isolationID075Cone30);
    muonArgSet.add(isolationID075Cone40);
    muonArgSet.add(isolationID075Cone50);
    muonArgSet.add(isolationID1Cone10);
    muonArgSet.add(isolationID1Cone15);
    muonArgSet.add(isolationID1Cone20);
    muonArgSet.add(isolationID1Cone30);
    muonArgSet.add(isolationID1Cone40);
    muonArgSet.add(isolationID1Cone50);
    muonArgSet.add(isolationID2Cone10);
    muonArgSet.add(isolationID2Cone15);
    muonArgSet.add(isolationID2Cone20);
    muonArgSet.add(isolationID2Cone30);
    muonArgSet.add(isolationID2Cone40);
    muonArgSet.add(isolationID2Cone50);
    muonArgSet.add(isolationID3Cone10);
    muonArgSet.add(isolationID3Cone15);
    muonArgSet.add(isolationID3Cone20);
    muonArgSet.add(isolationID3Cone30);
    muonArgSet.add(isolationID3Cone40);
    muonArgSet.add(isolationID3Cone50);
    muonArgSet.add(isolationID4Cone10);
    muonArgSet.add(isolationID4Cone15);
    muonArgSet.add(isolationID4Cone20);
    muonArgSet.add(isolationID4Cone30);
    muonArgSet.add(isolationID4Cone40);
    muonArgSet.add(isolationID4Cone50);
    muonArgSet.add(isolationID5Cone10);
    muonArgSet.add(isolationID5Cone15);
    muonArgSet.add(isolationID5Cone20);
    muonArgSet.add(isolationID5Cone30);
    muonArgSet.add(isolationID5Cone40);
    muonArgSet.add(isolationID5Cone50);
    muonArgSet.add(isolationID6Cone10);
    muonArgSet.add(isolationID6Cone15);
    muonArgSet.add(isolationID6Cone20);
    muonArgSet.add(isolationID6Cone30);
    muonArgSet.add(isolationID6Cone40);
    muonArgSet.add(isolationID6Cone50);

	TList _grEffRad; 
	TList _grPurity; 

	std::vector <float> centralityBins;
    centralityBins.push_back(0.0);
    centralityBins.push_back(0.05);
    centralityBins.push_back(0.10);
    centralityBins.push_back(0.15);
    centralityBins.push_back(0.20);
    centralityBins.push_back(0.40);
    centralityBins.push_back(0.80);

    //ncoll
    std::vector<double> ncoll;
    ncoll.push_back(1683.3); //0-5
    ncoll.push_back(1318.0); //5-10
    ncoll.push_back(1035.4); //10-15
    ncoll.push_back(811.2); //15-20
    ncoll.push_back(440.6); //20-40
    ncoll.push_back(77.8); //40-80

	const int nCentralityBins = centralityBins.size()-1;

	std::vector <int> vecConeRad; 
    vecConeRad.push_back(10);
    vecConeRad.push_back(15);
    vecConeRad.push_back(20);
    vecConeRad.push_back(30);
    vecConeRad.push_back(40);
    vecConeRad.push_back(50);
    const int nConeRadii = vecConeRad.size();

    std::vector <std::string> vecIDTrkPt;
    vecIDTrkPt.push_back("05");
    vecIDTrkPt.push_back("075");
    vecIDTrkPt.push_back("1");
    vecIDTrkPt.push_back("2");
    vecIDTrkPt.push_back("3");
    vecIDTrkPt.push_back("4");
    vecIDTrkPt.push_back("5");
    vecIDTrkPt.push_back("6");
    const int nTrkPtBins = vecIDTrkPt.size();


	TGraph* _grBkg = new TGraph(nConeRadii);
	TGraph* _grSig = new TGraph(nConeRadii);

    std::vector <double> vecIsolCut;
    vecIsolCut.push_back(0.05);
    vecIsolCut.push_back(0.1);
    vecIsolCut.push_back(0.2);
    vecIsolCut.push_back(0.3);
    vecIsolCut.push_back(0.4);
    vecIsolCut.push_back(0.5);
    vecIsolCut.push_back(0.7);
    vecIsolCut.push_back(1.0);
    vecIsolCut.push_back(5.0);
    vecIsolCut.push_back(10.0);
    vecIsolCut.push_back(100.0);
	const int nIsolCuts = vecIsolCut.size();

    ///Number of signal counts in data before
    ///correction, after background subtraction, and
    ///after W selection
    std::vector<double> _ns;
    _ns.push_back(1338);
    _ns.push_back(1034);
    _ns.push_back(909);
    _ns.push_back(758);
    _ns.push_back(1828);
    _ns.push_back(658);

    ///Number of background counts in data before
    ///correction, after background subtraction, and
    ///after W selection
    std::vector<double> _nb;
    _nb.push_back(139);
    _nb.push_back(107);
    _nb.push_back(100);
    _nb.push_back(68);
    _nb.push_back(167);
    _nb.push_back(84);

	RooArgSet centralityArgSet(centrality);  
  	// --- Fill primary data sets ---
    ///Signal
    RooDataSet* sigSet = fillDataSet(baseString,fileNameSigIn+".root",muonArgSet, cutValue); sigSet->Print(); 	
	sigSet = (RooDataSet*)sigSet->reduce(Cut(cutsSig));  sigSet = (RooDataSet*)sigSet->reduce(Cut("motherRec==24"));
    sigSet->Print(); 
    std::cout << "Signal events b4 isolation cut: " << sigSet->numEntries() << std::endl;

    ///Background
    RooDataSet* mcJ1Set = fillDataSet(baseString,fileNameMCJ1In+".root",muonArgSet, cutValue); mcJ1Set->Print();
	mcJ1Set = (RooDataSet*)mcJ1Set->reduce(Cut(cutsSig)); mcJ1Set->Print();
    RooDataSet* mcJ1Events = fillHIEventDataSet(baseString,fileNameMCJ1In+".root",centralityArgSet );

    RooDataSet* mcJ2Set = fillDataSet(baseString,fileNameMCJ2In+".root",muonArgSet, cutValue); mcJ2Set->Print();
	mcJ2Set = (RooDataSet*)mcJ2Set->reduce(Cut(cutsSig)); mcJ2Set->Print();
    RooDataSet* mcJ2Events = fillHIEventDataSet(baseString,fileNameMCJ2In+".root",centralityArgSet );

    RooDataSet* mcJ3Set = fillDataSet(baseString,fileNameMCJ3In+".root",muonArgSet, cutValue); mcJ3Set->Print();
	mcJ3Set = (RooDataSet*)mcJ3Set->reduce(Cut(cutsSig)); mcJ3Set->Print();
    RooDataSet* mcJ3Events = fillHIEventDataSet(baseString,fileNameMCJ3In+".root",centralityArgSet );

    // --- Subdivide in bins ---

	RooDataSet* mcJ1EventSubSet[nCentralityBins];
	RooDataSet* mcJ2EventSubSet[nCentralityBins];
	RooDataSet* mcJ3EventSubSet[nCentralityBins];

    RooDataSet* sigSubSet[nCentralityBins];
    RooDataSet* mcJ1SubSet[nCentralityBins];
    RooDataSet* mcJ2SubSet[nCentralityBins];
    RooDataSet* mcJ3SubSet[nCentralityBins];

    //centrality,trkpt,coneradius,isolation cut
    RooDataSet* sigSubSetIsol[nCentralityBins][nTrkPtBins][nConeRadii][nIsolCuts];
    RooDataSet* mcJ1SubSetIsol[nCentralityBins][nTrkPtBins][nConeRadii][nIsolCuts];
    RooDataSet* mcJ2SubSetIsol[nCentralityBins][nTrkPtBins][nConeRadii][nIsolCuts];
    RooDataSet* mcJ3SubSetIsol[nCentralityBins][nTrkPtBins][nConeRadii][nIsolCuts];
  		
    RooBinning b = RooBinning(nbins,0.0,ptmax); // 1 GeV per bin

    for(int icent=0; icent<nCentralityBins; icent++){

        std::cout << "Centrality " << centralityBins[icent] << "-" << centralityBins[icent+1] << std::endl;
		//number of events in the Jx sample in this centrality bin 
		mcJ1EventSubSet[icent] = selectCentrality(mcJ1Events,centralityBins[icent], centralityBins[icent+1]);
		mcJ2EventSubSet[icent] = selectCentrality(mcJ2Events,centralityBins[icent], centralityBins[icent+1]);
		mcJ3EventSubSet[icent] = selectCentrality(mcJ3Events,centralityBins[icent], centralityBins[icent+1]);

        mcJ1SubSet[icent] = selectCentrality(mcJ1Set,centralityBins[icent], centralityBins[icent+1]);
        std::cout << "Number of J1 events before isolation cut: " <<
            mcJ1SubSet[icent]->numEntries() << std::endl;

        mcJ2SubSet[icent] = selectCentrality(mcJ2Set,centralityBins[icent], centralityBins[icent+1]);
        std::cout << "Number of J2 events before isolation cut: " <<
            mcJ2SubSet[icent]->numEntries() << std::endl;

        mcJ3SubSet[icent] = selectCentrality(mcJ3Set,centralityBins[icent], centralityBins[icent+1]);
        std::cout << "Number of J3 events before isolation cut: " <<
            mcJ3SubSet[icent]->numEntries() << std::endl;

        TH1F* hmcJ1Set = (TH1F*)mcJ1SubSet[icent]->createHistogram("hmcJ1Set",muonPt,Binning(b));
        TH1F* hmcJ2Set = (TH1F*)mcJ2SubSet[icent]->createHistogram("hmcJ2Set",muonPt,Binning(b));
        TH1F* hmcJ3Set = (TH1F*)mcJ3SubSet[icent]->createHistogram("hmcJ3Set",muonPt,Binning(b));

        //return weighted sample by Jx cross section
		TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",nbins,0.0,ptmax);
		hmcQCDSet = getWeightedJxHisto(centralityBins[icent], centralityBins[icent+1],
            ncoll[icent],hmcJ1Set,hmcJ2Set,hmcJ3Set, mcJ1EventSubSet[icent],
            mcJ2EventSubSet[icent],mcJ3EventSubSet[icent], nbins,0.0, ptmax);
        //delete mcJ1EventSubSet[icent]; delete mcJ2EventSubSet[icent]; delete mcJ3EventSubSet[icent]; 

        //subset WITHOUT isolation cut
        sigSubSet[icent] = selectCentrality(sigSet,centralityBins[icent], centralityBins[icent+1]);
        std::cout << "Number of signal events before isolation cut: " << 
            sigSubSet[icent]->numEntries() << std::endl;

        TH1F* hsigSubSet = (TH1F*)sigSubSet[icent]->createHistogram("hsigSubSet",muonPt,Binning(b));
        std::cout << "Sanity check: " << hsigSubSet->Integral() << "?=" << sigSubSet[icent]->numEntries() << std::endl;

        ///N_S in centrality bin icent before isolation
        //double nSignal = hsigSubSet->Integral();
        ///N_B in centrality bin icent before isolation
        //double nBackground = hmcQCDSet->Integral();

        ///Loop over all lower track thresholds
        for(int iTrk=0; iTrk<nTrkPtBins; iTrk++){
	
            ///Loop over all cone radii
            for(int iCone=0; iCone<nConeRadii; iCone++){

                _grEffRad.Add( new TGraphAsymmErrors(nIsolCuts));

                _grPurity.Add( new TGraph(nIsolCuts));
	
                //name to save graph as
                TString sGr = "grIsoEffTrkPt"; sGr+=vecIDTrkPt[iTrk]; 
                sGr+="cent";  sGr+=icent; 
                sGr+= "ConeRadius"; sGr+=vecConeRad[iCone];

                TString sGrRatio = sGr + "SigBkgPurity"; 

                //unique index for each graph with given centrality,
                //lower track pt, and cone radius
                int index = (icent*nTrkPtBins+iTrk)*nConeRadii+iCone;

                ///Loop over all isolation upper limits
                for(int iIsol=0; iIsol<nIsolCuts; iIsol++){
		            
                    //subset WITH isolation cut
                    sigSubSetIsol[icent][iTrk][iCone][iIsol] = selectIDConeIsol(sigSubSet[icent],
                        vecIDTrkPt[iTrk],vecConeRad[iCone] ,vecIsolCut[iIsol]);

                    std::cout << "Number of signal events after isolation cut: " <<
                        sigSubSetIsol[icent][iTrk][iCone][iIsol]->numEntries() << std::endl;
                
                    TH1F* hsigSubSetIsol = (TH1F*)sigSubSetIsol[icent][iTrk][iCone][iIsol]->createHistogram("hsigSubSetIsol",muonPt,Binning(b));

                    // --- QCD set ---
                    mcJ1SubSetIsol[icent][iTrk][iCone][iIsol] = selectIDConeIsol(mcJ1SubSet[icent],
                        vecIDTrkPt[iTrk],vecConeRad[iCone] ,vecIsolCut[iIsol]);
                    std::cout << "Number of J1 events after isolation cut: " <<
                        mcJ1SubSetIsol[icent][iTrk][iCone][iIsol]->numEntries() << std::endl;

                    mcJ2SubSetIsol[icent][iTrk][iCone][iIsol] = selectIDConeIsol(mcJ2SubSet[icent],
                         vecIDTrkPt[iTrk],vecConeRad[iCone] ,vecIsolCut[iIsol]);
                    std::cout << "Number of J2 events after isolation cut: " <<
                        mcJ2SubSetIsol[icent][iTrk][iCone][iIsol]->numEntries() << std::endl;

                    mcJ3SubSetIsol[icent][iTrk][iCone][iIsol] = selectIDConeIsol(mcJ3SubSet[icent],
                        vecIDTrkPt[iTrk],vecConeRad[iCone] ,vecIsolCut[iIsol]);
                    std::cout << "Number of J3 events after isolation cut: " <<
                        mcJ3SubSetIsol[icent][iTrk][iCone][iIsol]->numEntries() << std::endl;

                TH1F* hmcJ1SetIsol = (TH1F*)mcJ1SubSetIsol[icent][iTrk][iCone][iIsol]->createHistogram("hmcJ1SetIsol",muonPt,Binning(b));
                TH1F* hmcJ2SetIsol = (TH1F*)mcJ2SubSetIsol[icent][iTrk][iCone][iIsol]->createHistogram("hmcJ2SetIsol",muonPt,Binning(b));
                TH1F* hmcJ3SetIsol = (TH1F*)mcJ3SubSetIsol[icent][iTrk][iCone][iIsol]->createHistogram("hmcJ3SetIsol",muonPt,Binning(b));

                ///Get Bayesian errors by calculating efficiency
                ///of the added J123mu sample before and after isolation cut

                double bkgEffJx,bkgBayesLoJx,bkgBayesUpJx;
                double nJxMuNonIsol = hmcJ1Set->Integral(ptSigBinLo,ptSigBinUpp) + hmcJ2Set->Integral(ptSigBinLo,ptSigBinUpp) + hmcJ3Set->Integral(ptSigBinLo,ptSigBinUpp);
                double nJxMuIsol = hmcJ1SetIsol->Integral(ptSigBinLo,ptSigBinUpp) + 
                                    hmcJ2SetIsol->Integral(ptSigBinLo,ptSigBinUpp) + hmcJ3SetIsol->Integral(ptSigBinLo,ptSigBinUpp);

                ///Calculate efficiency of J123mu sample 
                Efficiency((int)nJxMuIsol,(int)nJxMuNonIsol,0.683,bkgEffJx,bkgBayesLoJx,bkgBayesUpJx);
                bkgBayesUpJx = bkgBayesUpJx - bkgEffJx;
                bkgBayesLoJx = bkgEffJx - bkgBayesLoJx;
                std::cout << "Background efficiency for Jx = " << bkgEffJx << "-" << bkgBayesLoJx << " +" <<
                    bkgBayesUpJx << std::endl; 

                double relBkgEffJxUp = bkgBayesUpJx/bkgEffJx; 
                double relBkgEffJxLo =  bkgBayesLoJx/bkgEffJx;
                std::cout << "Background efficiency for J1,2,3mu = " << " + " << relBkgEffJxUp*100. << "%" << " -" << relBkgEffJxLo*100. << "%" << std::endl; 

                ///Now we weight the sample
                TH1F* hmcQCDSetIsol = new TH1F("hmcQCDSetIsol","hmcQCDSetIsol",nbins,0.0,ptmax);
		        hmcQCDSetIsol = getWeightedJxHisto(centralityBins[icent], centralityBins[icent+1],
                        ncoll[icent],hmcJ1SetIsol,hmcJ2SetIsol,hmcJ3SetIsol, mcJ1EventSubSet[icent],
                        mcJ2EventSubSet[icent],mcJ3EventSubSet[icent], nbins,0.0, ptmax);

                std::cout << "Efficiency defined in pT range " << hmcQCDSetIsol->GetBinLowEdge(ptSigBinLo) << "-" << hmcQCDSetIsol->GetBinLowEdge(ptSigBinUpp);

                double bkgCounts = hmcQCDSet->Integral(ptSigBinLo,ptSigBinUpp);
                double bkgCountsPassed = hmcQCDSetIsol->Integral(ptSigBinLo,ptSigBinUpp);
                double sigCounts = hsigSubSet->Integral(ptSigBinLo,ptSigBinUpp);
                double sigCountsPassed = hsigSubSetIsol->Integral(ptSigBinLo,ptSigBinUpp);
                double bkgEff, sigEff, bkgBayesUp, bkgBayesLo, sigBayesUp,sigBayesLo;

                bkgEff = bkgCountsPassed/bkgCounts;

                ///Calculate Bayesian errors
                Efficiency((int)sigCountsPassed,(int)sigCounts,0.683,sigEff,sigBayesLo,sigBayesUp);
                sigBayesUp = sigBayesUp - sigEff;
                sigBayesLo = sigEff - sigBayesLo;
                std::cout << "Signal efficiency = " << sigEff << " +" << sigBayesUp << " -" << sigBayesLo << std::endl;
                std::cout << "Background efficiency = " << bkgEff << " +" << relBkgEffJxUp*bkgEff << " -" << relBkgEffJxLo*bkgEff << std::endl;

                ///Calculate effective signal
                double a = sigEff*_ns[icent];
                double b = bkgEff*_nb[icent];
                double purity = a/(a+b);

                double nEff = (purity*a)/(2-purity);

                ((TGraphAsymmErrors*)_grEffRad.At(index))->SetPoint(iIsol,sigEff,bkgEff);
                ((TGraphAsymmErrors*)_grEffRad.At(index))->SetPointError(iIsol,sigBayesLo,sigBayesUp,relBkgEffJxLo*bkgEff,relBkgEffJxUp*bkgEff);

                std::cout << "That gives an effective signal of " << nEff <<std::endl;
                ((TGraph*)_grPurity.At(index))->SetPoint(iIsol,vecIsolCut[iIsol],nEff);

               }//iIsol
                //after looping through all the upper cuts, save the graphs
	            Write(outFile,_grEffRad.At(index), sGr);
	            Write(outFile,_grPurity.At(index), sGrRatio);
            }//iCone
        }//iTrk

    }//icent

   outFile->Close();
    ///Loop over all cone radii
    for(int index=0; index<_grEffRad.GetEntries(); ++index){
        delete _grEffRad.At(index);
        delete _grPurity.At(index);
   }
   delete _grBkg;
   delete _grSig;
   delete outFile;

}

int main(){
    isoEfficiency();
}
