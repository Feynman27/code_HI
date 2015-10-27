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

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooAbsArg.h"
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

#include "correctionFactorsDep.C"

using namespace RooFit ;
using namespace RooStats ;

///////////////////////////////////////////////////////////////////////////////
//weightDS 
///////////////////////////////////////////////////////////////////////////////
RooDataSet* weightDS(RooDataSet* set, double wgt){

  RooRealVar wt("wt","wt",wgt);
  set->addColumn(wt);

  RooDataSet* wSet = new RooDataSet(set->GetName(),set->GetTitle(),set,*set->get(),0,wt.GetName()) ;
  std::cout << "Returning weighted d.s." << std::endl;
  std::cout << std::endl;
  wSet->Print();
  return wSet;
    
}


void plotBkgFraction(TH1F* hCutSet, TH1F* hGenSet, float CwAw, TGraphAsymmErrors* grFraction,float ieta,double etaLo,double etaUpp){
    
    int nPassed = (int)hCutSet->Integral(); 
    ///Generator filter efficiency from AMI
    float effFilter = 0.13754; 
    ///Get total number of muons from tau at any pT
    ///(i.e. not just above 17GeV where filter was set)
    int nTotal = (int)hGenSet->Integral()/effFilter;
    std::cout << "Number of muons from tau decays passed : " << nPassed << " total: " << nTotal << std::endl;
    double bkgFraction, errHi, errLo;
    Efficiency(nPassed,nTotal,0.683,bkgFraction,errLo,errHi);
    errHi = errHi-bkgFraction;
    errLo = bkgFraction-errLo;
    std::cout << std::endl;
    std::cout << "That gives a survival probability of " << bkgFraction << " +" << errHi << " -" << errLo << std::endl;
    std::cout << "Branching ratio of W-->tau-->mu:      17.4%(From PDG)" << std::endl;
    std::cout << "Ratio of W-->tau/W-->mu:              11.25/10.57 = " << 11.25/10.57 << "(From PDG)" << std::endl;
    std::cout << "[Cw*Aw]^{-1} =                       "<< 1.0/CwAw << std::endl;
    float branchRatio = 0.174, decayRatio=11.25/10.57;
    bkgFraction *=branchRatio; bkgFraction *= decayRatio; bkgFraction /= CwAw; 
    errHi *=branchRatio; errHi *= decayRatio; errHi /= CwAw; 
    errLo *=branchRatio; errLo *= decayRatio; errLo /= CwAw; 
    std::cout << "That gives a fraction of tau events per signal event (W-->tau-->mu/W-->mu) of " << bkgFraction*100. << "%" 
        << " +" << errHi*100. << " -" << errLo*100. << std::endl;
    std::cout << std::endl;
    double xEta = etaLo+(etaLo-etaUpp)/2.;
    grFraction->SetPoint(ieta,xEta,bkgFraction);
    grFraction->SetPointError(ieta,(etaLo-etaUpp)/2.,(etaLo-etaUpp)/2.,errLo,errHi);
}

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
int getEtaBin(float eta){
    int index = 0;
    if(fabs(eta)>0.1&&fabs(eta)<0.35) return index; else ++index;
    if(fabs(eta)>0.35&&fabs(eta)<0.6) return index; else ++index;
    if(fabs(eta)>0.6&&fabs(eta)<0.8) return index; else ++index;
    if(fabs(eta)>0.8&&fabs(eta)<1.05) return index; else ++index;
    if(fabs(eta)>1.05&&fabs(eta)<1.3) return index; else ++index;
    if(fabs(eta)>1.3&&fabs(eta)<1.55) return index; else ++index;
    if(fabs(eta)>1.55&&fabs(eta)<1.85) return index; else ++index;
    if(fabs(eta)>1.85&&fabs(eta)<2.1) return index; else ++index;
    if(fabs(eta)>2.1&&fabs(eta)<2.4) return index; else ++index;
    return -1;
}

/*RooDataSet* fillHIWTauDataSet(const TString& pathName, const TString& fileName, RooArgSet& muonArgSet, double wt=1.0){
    
  RooDataSet* set = new RooDataSet(fileName,fileName,muonArgSet);
  // --- Load the MC tree ---
  TChain* tree = new TChain("tree","tree");
  int nFiles = tree->Add(pathName+fileName);
  
  std::cout << "Filling the RooDataSet for " << fileName << "... Number of files: " << nFiles << std::endl;
  int nmu;
  float mc_mt[50], mc_mptPhi,mc_mptPx,mc_mptPy,mc_mptPz,mc_mptPt,centrality;
  float  mc_pt[50],  mc_ptNominal[50],mc_eta[50],mc_phi[50],mc_pdgId[50],mc_charge[50];
  float mc_M[50];

  tree->SetBranchAddress("nTauMu",&nmu);
  tree->SetBranchAddress("mc_mt",&mc_mt);
  tree->SetBranchAddress("mc_mptPhi",&mc_mptPhi);
  tree->SetBranchAddress("mc_mptPt",&mc_mptPt);
  tree->SetBranchAddress("centrality",&centrality);
  tree->SetBranchAddress("mc_pt",&mc_pt);
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
  tree->SetBranchStatus("mc_eta",1);
  tree->SetBranchStatus("mc_phi",1);
  tree->SetBranchStatus("mc_pdgId",1);
  tree->SetBranchStatus("mc_charge",1);

  std::cout << "Number of entries: " << tree->GetEntries() << std::endl; 
  for ( int i = 0; i < tree->GetEntries(); i++ ) {
//     if (i>10000) break; // temp hack
    tree->LoadTree(i);
    tree->GetEntry(i);

    muonArgSet.setRealValue("missPt",mc_mptPt);
    muonArgSet.setRealValue("centrality",centrality);

    for (int imu = 0; imu<nmu; ++imu){

      muonArgSet.setRealValue("w",wt);
      muonArgSet.setRealValue("muonPt",mc_pt[imu]);
      muonArgSet.setRealValue("muonMt",mc_mt[imu]);
      muonArgSet.setRealValue("muonEta",mc_eta[imu]);

      if (mc_charge[imu] > 0. ) muonArgSet.setCatLabel("chargeCategory","muPlus");
      else if (mc_charge[imu] < 0.) muonArgSet.setCatLabel("chargeCategory","muMinus");

      set->add(muonArgSet);    
    }//imu
   }//i
   
  return set;
}
*/
void bkgWtauPlotter(){

    ///Switch for binning
    bool doEta = true;
    bool doCentrality = true;
    bool doSystematic = false;

    ///Open file with product of CwAw as fcn
    ///of eta/centrality
    //TFile* _fCwAw = new TFile("CorrectionFactorFiles/correctionFactorsW_binbybin_9EtaBins6CentBins.07.24.2013.root","read");
    TFile* _fCwAw = new TFile("CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges.PowPy8.11.27.2013.root","read");
    //cross-check w/o mpt cut
    //TFile* _fCwAw = new TFile("crossChecks/correctionFactorsW_noMptCut.12.10.2013.root","read");
    ///Open file with efficiencies for muons passing
    ///PS,Z veto, isolation, and pT>25
//    TFile* _fCutEff = new TFile("background/mcWEffForTauMuonStudy_07_23_2013.root","read");
//
    ///Use this file for systematic study
    TFile* _fCutEff = NULL;
    if(doSystematic){
        std::cout << "WARNING! WARNING! WARNING! " << std::endl;
        std::cout << "Using different efficiencies for systematic study." << std::endl;
        _fCutEff = new TFile("systematics/mcWEffForTauMuonStudy_Systematics_07_28_2013.root","read");
    }
    else _fCutEff = new TFile("background/mcWEffForTauMuonStudy_07_23_2013.root","read");
    ///Open custom embedded Wtau ntuple
//    TFile* _fWtau = new TFile("","read");

    TFile* outFile = new TFile("fractionTauEtaCent_9etaBins6centBins.root","recreate");
	TString baseString = "/usatlas/u/tbales/scratch/";
    TString fileNameIn = "MonteCarloFiles/Wtaumu/HIWtaumuNtuple.07.30.2013";

    std::vector<double> centralityBins;
    centralityBins.push_back(0.0);
    if(doCentrality){
        centralityBins.push_back(0.05);
        centralityBins.push_back(0.10);
        centralityBins.push_back(0.15);
        centralityBins.push_back(0.20);
        centralityBins.push_back(0.40);
    }
    centralityBins.push_back(0.80);

    const int nCentralityBins = centralityBins.size()-1;

    std::vector<double> ptBins;
    ptBins.push_back(0.0);
    ptBins.push_back(300.0);
    const int nPtBins = ptBins.size()-1;

    std::vector<double> etaBins;
    etaBins.push_back(0.1);
    if(doEta){
        etaBins.push_back(0.35);
        etaBins.push_back(0.6);
        etaBins.push_back(0.8);
        etaBins.push_back(1.05);
        etaBins.push_back(1.37);
        etaBins.push_back(1.52);
        etaBins.push_back(1.74);
        etaBins.push_back(2.1);
 
    }
    etaBins.push_back(2.4);
    const int nEtaBins = etaBins.size()-1;

    // --- declare cut variables --- //
    RooRealVar muonPt("muonPt","p_{T}",0.0,350.0,"GeV");
    RooRealVar  missPt("missPt","p_{T}^{miss}",0.0,350.0,"GeV");
    RooRealVar  muonMt("muonMt","m_{T}",0.0,350.0,"GeV");
    RooRealVar  muonCharge("muonCharge","charge",-2.0,+2.0);
    RooRealVar  centrality("centrality","centrality",0.,1.0);
    RooRealVar  muonEta("muonEta","muonEta",-3.0,+3.0);
    RooRealVar  w("w","w",0.0,10.0);

    TString sCutsSig = "";
    sCutsSig = "abs(muonEta)>0.1&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0&&centrality>0.&&centrality<0.8";
    //cross-check w/o mpt cut
    //sCutsSig = "abs(muonEta)>0.1&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>0.0&&missPt<9000.0&&muonMt>40.0&&centrality>0.&&centrality<0.8";
    if(doSystematic) sCutsSig = "muonPt>25.0&&abs(muonEta)>0.1&&abs(muonEta)<2.4&&centrality>0.&&centrality<0.8";
    std::cout << " Signals cuts: "<< sCutsSig << std::endl;

    RooArgList muonTauArgList(muonEta,muonPt,missPt,muonMt,centrality);

	RooFormulaVar cutsSig("cutsSig", "cutsSig", sCutsSig, muonTauArgList);

    RooCategory chargeCategory("chargeCategory","chargeCategory") ;
    chargeCategory.defineType("muMinus",-1) ;
    chargeCategory.defineType("muPlus",1) ;
    chargeCategory.Print();

    RooArgSet muonArgSet(muonPt,missPt,muonMt,muonEta,centrality,chargeCategory,w);

    ///Fill dataset
    RooDataSet* mcTauSet = fillHIWTauDataSet(baseString,fileNameIn+".root",muonArgSet); mcTauSet->Print(); 
    ///apply selection cuts
	RooDataSet* mcTauCutSet = (RooDataSet*)mcTauSet->reduce(Cut(cutsSig)); mcTauCutSet->Print();
    ///Create subsets
    RooDataSet* mcTauSubSet[nPtBins][nEtaBins][nCentralityBins];
    TH1F* hMcTauSubSet[nPtBins][nEtaBins][nCentralityBins];
    RooDataSet* mcTauCutSubSet[nPtBins][nEtaBins][nCentralityBins];
    TH1F* hMcTauCutSubSet[nPtBins][nEtaBins][nCentralityBins];

  	RooBinning b = RooBinning(170,0.0,510.0); // 3 GeV per bin
    std::cout << "Creating subsets..." << std::endl;
    char cTauGen[50],cTauCut[50];
    float cutEff;
	for ( int i = 0; i < nPtBins; i++ ) {
	  for ( int j = 0; j < nEtaBins; j++ ) {
	    for ( int k = 0; k < nCentralityBins; k++ ){

            sprintf(cTauGen,"hTauGen_Eta%i_Cent_%i",j,k);
            sprintf(cTauCut,"hTauCut_Eta%i_Cent_%i",j,k);

            TString sCutEff = "effForTauStudy_cent"; sCutEff+= k;
            TGraphAsymmErrors* grCutEff = (TGraphAsymmErrors*)_fCutEff->Get(sCutEff);
            ///Differential in centrality and eta
            if(doCentrality&&doEta) cutEff = grCutEff->GetY()[j];
            ///Averaged over all centrality and eta
            else cutEff =0.738509 ;
            ///bin in pt,eta,and centrality
            mcTauSubSet[i][j][k] = selectPtEtaCentrality(mcTauSet,ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1],centralityBins[k], centralityBins[k+1], true);
	        hMcTauSubSet[i][j][k] = (TH1F*)mcTauSubSet[i][j][k]->createHistogram(cTauGen,muonMt,Binning(b));
            std::cout << "Number of entries at generator level after selection for eta and centrality: " <<
                hMcTauSubSet[i][j][k]->Integral() << std::endl;
            std::cout << "Consistency check: " << mcTauSubSet[i][j][k]->numEntries() << "=?" <<
                hMcTauSubSet[i][j][k]->Integral() << std::endl;

            mcTauCutSubSet[i][j][k] = selectPtEtaCentrality(mcTauCutSet,ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1],centralityBins[k], centralityBins[k+1], true);
            std::cout << "Number of entries in cut set before applying efficiency: " << mcTauCutSubSet[i][j][k]->numEntries() << std::endl;
            w.setVal(cutEff);
            std::cout << "Weighting d.s. by cut efficiency: " << w.getVal() << std::endl;
            ///apply cut efficiency
            mcTauCutSubSet[i][j][k] = weightDS(mcTauCutSubSet[i][j][k],cutEff); 
	        hMcTauCutSubSet[i][j][k] = (TH1F*)mcTauCutSubSet[i][j][k]->createHistogram(cTauCut,muonMt,Binning(b));
            std::cout << "Number of events in cut subset after weighting: " << hMcTauCutSubSet[i][j][k]->Integral() << std::endl;
            //mcTauSubSet[i][j][k]->Print();
        }//k
      }//j
    }//i

    ///TGraph of tau bkg fraction
    TList _fraction;
    float CwAw;
	for ( int i = 0; i < nPtBins; i++ ) {
	  for ( int icent = 0; icent < nCentralityBins; icent++ ){
         
         _fraction.Add(new TGraphAsymmErrors(nEtaBins));
	    for ( int ieta = 0; ieta < nEtaBins; ieta++ ) {

            std::cout << " bin "<<i<<":"<<ieta<<":"<<icent<<std::endl;
            TString sCwAw = "grCwAwEtaDistroCent"; sCwAw+=icent;
            TGraphAsymmErrors* grCwAw = (TGraphAsymmErrors*)_fCwAw->Get(sCwAw);
            ///get CwAw product
            ///Differential in centrality and eta
            if(doCentrality&&doEta) CwAw = grCwAw->GetY()[ieta];
            ///Averaged over all centrality and eta
            else CwAw = 0.359;

            //plotBkgFraction(mcTauCutSubSet[i][ieta][icent],mcTauSubSet[i][ieta][icent],CwAw, (TGraphAsymmErrors*)_fraction.At(icent), ieta, etaBins[ieta], etaBins[ieta+1]);
            plotBkgFraction(hMcTauCutSubSet[i][ieta][icent],hMcTauSubSet[i][ieta][icent],CwAw, (TGraphAsymmErrors*)_fraction.At(icent), ieta, etaBins[ieta], etaBins[ieta+1]);

        }//k
        TString sName = "tauBkgFractionCent"; sName+=icent;
        ((TGraphAsymmErrors*)_fraction.At(icent))->SetName(sName);
        outFile->cd();
        ((TGraphAsymmErrors*)_fraction.At(icent))->Write();
      }//j
    }//k
        
    std::cout << "Clean up." << std::endl;
    for(int i=0; i<_fraction.GetEntries(); ++i){

        delete _fraction.At(i);
    }
}
