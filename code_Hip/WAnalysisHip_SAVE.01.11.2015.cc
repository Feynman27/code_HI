#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath> 
#include <ctime> 
#include <sstream>
#include <vector>
#include <string>

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

#include "WAnalysisTool.h"
#include "WBookKeeper.cc"
#include "AtlasUtils.C"
#include "AtlasStyle.C"

using namespace RooFit;

RooDataSet* fillMuonDataSet(const TString& pathName, const TString& fileName, RooArgSet& muonArgSet, bool doPeriodA = false, bool doPeriodB=false, bool isMc=false){

	TChain* tree = new TChain("treeW","treeW");
	int nFiles = tree->Add(pathName+fileName);
    std::cout << "Filling the RooDataSet for " << fileName << "... Number of files: " << nFiles << std::endl;
    std::cout << "Number of entries in data tree: " << tree->GetEntries() << std::endl;

	RooDataSet* set = new RooDataSet(fileName,fileName,muonArgSet);
	float _pt       = 0;
	float _met      = 0;
	float _mt       = 0;
	float _charge   = 0;
	float _ptcone30 = 0;
	float _eta      = 0;
	float _centrality = 0;
    bool  _isPeriodA = false;
    bool  _isPeriodB = false;

	tree->SetBranchAddress("pt",&_pt);
	tree->SetBranchAddress("met",&_met);
	tree->SetBranchAddress("mt",&_mt);
	tree->SetBranchAddress("charge",&_charge);
	tree->SetBranchAddress("ptcone30",&_ptcone30);
	tree->SetBranchAddress("eta",&_eta);
	tree->SetBranchAddress("centrality",&_centrality);
	tree->SetBranchAddress("isPeriodA",&_isPeriodA);
	tree->SetBranchAddress("isPeriodB",&_isPeriodB);
	// --- Set branch status ---
  	tree->SetBranchStatus("*",0) ;
	tree->SetBranchStatus("pt",1);
	tree->SetBranchStatus("met",1);
	tree->SetBranchStatus("mt",1);
	tree->SetBranchStatus("charge",1);
	tree->SetBranchStatus("ptcone30",1);
	tree->SetBranchStatus("eta",1);
	tree->SetBranchStatus("centrality",1);
	tree->SetBranchStatus("isPeriodA",1);
	tree->SetBranchStatus("isPeriodB",1);
	
    int nmu = tree->GetEntries();
    for(int imu =0; imu < nmu; ++imu){

         tree->LoadTree(imu);
         tree->GetEntry(imu);

	      float isoTemp = _ptcone30/_pt;
          if(isoTemp>0.2) continue;
          muonArgSet.setRealValue("isolationMu",isoTemp);
	      muonArgSet.setRealValue("centrality" ,_centrality);
          // If combining periods, use PeriodB convention
          // i.e. proton is toward positive rapidity
          if( !(doPeriodA) && (_isPeriodA) ){
            float etaConv = -1.0*_eta;
	        muonArgSet.setRealValue("muonEta"  ,etaConv);
          }
          else muonArgSet.setRealValue("muonEta"    ,_eta);
	      muonArgSet.setRealValue("muonPt"     ,_pt);
	      muonArgSet.setRealValue("missEt"     ,_met);
	      muonArgSet.setRealValue("muonMt"     ,_mt);
	      muonArgSet.setRealValue("isPeriodA"     ,_isPeriodA);
	      muonArgSet.setRealValue("isPeriodB"     ,_isPeriodB);
          if(_charge<0.0) muonArgSet.setCatLabel("chargeCategory","muMinus");
	      else if(_charge>0.0) muonArgSet.setCatLabel("chargeCategory","muPlus"); 

	      set->add(muonArgSet);
    }//imu
	return set;
}

// For temporty use with private Mc @ generator level
RooDataSet* fillMuonMcSet(const TString& pathName, const TString& fileName, RooArgSet& muonArgSet){

	TChain* T = new TChain("T","T");
	int nFiles = T->Add(pathName+fileName);
        std::cout << "Filling the RooDataSet for " << fileName << "... Number of files: " << nFiles << std::endl;
	std::cout << "Number of entries in Mc tree: " << T->GetEntries() << std::endl;

	RooDataSet *set = new RooDataSet(fileName,fileName,muonArgSet);
	int muonMother[50], nmu;
	float muonEta[50],muonPt[50],muonPhi[50],neutrinoPt[50],mt[50];
	T->SetBranchAddress("muonMother",&muonMother);
	T->SetBranchAddress("muonEta",&muonEta);
	T->SetBranchAddress("muonPt",&muonPt);
	T->SetBranchAddress("neutrinoPt",&neutrinoPt);
	T->SetBranchAddress("muonPhi",&muonPhi);
	T->SetBranchAddress("mt",&mt);
	T->SetBranchAddress("nMuons",&nmu);
	
	T->SetBranchStatus("*",0);
	T->SetBranchStatus("nMuons",1);
	T->SetBranchStatus("muonMother",1);
	T->SetBranchStatus("muonEta",1);
	T->SetBranchStatus("muonPt",1);
	T->SetBranchStatus("neutrinoPt",1);
	T->SetBranchStatus("muonPhi",1);
	T->SetBranchStatus("mt",1);

	for(int iev = 0; iev < T->GetEntries(); ++iev){
           if(iev%10000==0) std::cout << "Event: " << iev << std::endl;
	   //if(iev==10000) break; //hack
           T->LoadTree(iev);
           T->GetEntry(iev);
	   for (int imu=0; imu<nmu; ++imu){
	     if( fabs(muonMother[imu]) == 24){
	      muonArgSet.setRealValue("missEt",neutrinoPt[imu]);
	      muonArgSet.setRealValue("muonPt",muonPt[imu]);
	      muonArgSet.setRealValue("muonEta",muonEta[imu]);
	      muonArgSet.setRealValue("muonMt",mt[imu]); 
		
	      set->add(muonArgSet);
	     }
	   }//imu
	}//iev

	return set;
}



void plot(TH1F* hCumulative, TString sEta, TString sCent, TString sCharge, int charge){

	TString sName = "WAnalysisHip_Mu"; sName+=charge;  
    TString cName = "c_"; 
    if(charge==1) cName+="plus";
    else if(charge==-1) cName+="minus";
    TString cTitle = "c,"+ sEta;
    cTitle+=","; cTitle+=sCent; cTitle+=","; cTitle+=sCharge;
	TCanvas *canvas = new TCanvas(cName,cTitle,600,600);	
	hCumulative->Draw("pe");
	TLatex* texEta = new TLatex(0.66,0.86,sEta);
	TLatex* texCent = new TLatex(0.72,0.79,sCent);
	texEta->SetNDC();
	texCent->SetNDC();
	texEta->Draw();
	texCent->Draw();
	canvas->Print(sName+".root");
}

void plot(TH1F* hCumulative, WAnalysisTool& toolHandle, RooRealVar& var, RooBinning b){
	
    double cutLow = b.lowBound();
    double cutUpp = b.highBound();

	toolHandle.FillHistogram(b,var,cutLow,cutUpp);
	TH1F* hDataSet = toolHandle.GetHistogram();
	//hDataSet->Scale(1.0,"width");
    hCumulative->Add(hDataSet);
	std::cout << "Current integral of cumulative histogram: " 
			<< hCumulative->Integral() << std::endl;
}

TString format(float value) {
  std::stringstream svalue;
  svalue  << std::setprecision(2) << value;
  return svalue.str();
}

void WAnalysisHip(const std::string& period,bool doBackgroundSub=false, bool doCorrectYields=false){

	bool doPeriodA = false;
	bool doPeriodB = false; 
	if(strcmp(period.c_str(),"A") == 0) {
		doPeriodA = true;
		std::cout << "*****************************" << std::endl;
		std::cout << "**Running on Period A data.**" << std::endl;
		std::cout << "*****************************" << std::endl;
	}
	else if(strcmp(period.c_str(),"B") == 0) {
		doPeriodB = true;
		std::cout << "*****************************" << std::endl;
		std::cout << "**Running on Period B data.**" << std::endl;
		std::cout << "*****************************" << std::endl;
	}
	else if(strcmp(period.c_str(),"AB") == 0) {
		std::cout << "********************************************************************" << std::endl;
		std::cout << "**                   Running on Period A and B data.              **" << std::endl;
		std::cout << "**Using Period B convention i.e. proton towards positive rapidity.**" << std::endl;
		std::cout << "********************************************************************" << std::endl;
	}
    else {
        std::cout << "ERROR: Invalid argument in WAnalysisHip(). Candidates are void WAnalysisHip(const std::string& period)" << std::endl;
        std::cout << "where period = A, B, or AB." << std::endl; return;
    }

	bool doPt = false;
	bool doEta = true;
	bool doCentrality = true;
	bool doCharge = true;

	// Declare variables to cut on
	RooRealVar  muonPt("muonPt","p_{T}",0.0,400.0,"GeV");
	RooRealVar  missEt("missEt","E_{T}^{miss}",0.0,400.0,"GeV");
	RooRealVar  muonMt("muonMt","m_{T}",0.0,400.0,"GeV");
	RooRealVar  muonCharge("muonCharge","charge",-2.0,+2.0);
	RooRealVar  isolationMu("isolationMu","isolationMu",0.0,10.0);
	RooRealVar  centrality("centrality","centrality",0.,1.0);
	RooRealVar  muonEta("muonEta","#eta^{#mu}",-3.0,+3.0);
    RooRealVar  isPeriodA("isPeriodA","A",-0.1,1.1);
    RooRealVar  isPeriodB("isPeriodB","A",-0.1,1.1);
	//RooRealVar  weight("weight","weight",0.0,10.0);
	
	TString sCutsSig = "centrality>0.0";
	RooArgList muonArgList(centrality);

	RooFormulaVar cutsSig("cutsSig", "cutsSig", sCutsSig, muonArgList);

  	RooCategory chargeCategory("chargeCategory","chargeCategory") ;
  	chargeCategory.defineType("muMinus",-1) ;
  	chargeCategory.defineType("muPlus",1) ;

	// Pt bins
	std::vector<double> ptBins;
	ptBins.push_back(25.0);
	ptBins.push_back(400.0);
	const int nPtBins = ptBins.size()-1;

	// Eta bins
	std::vector<double> etaBins;
	etaBins.push_back(-2.5);
	if(doEta){
	   etaBins.push_back(-2.1);
           //etaBins.push_back(-1.85);
           etaBins.push_back(-1.55);
           //etaBins.push_back(-1.3);
           etaBins.push_back(-1.05);
           //etaBins.push_back(-0.8);
           etaBins.push_back(-0.6);
           //etaBins.push_back(-0.35);
           etaBins.push_back(-0.1);
           etaBins.push_back(0.1);
           //etaBins.push_back(0.35);
           etaBins.push_back(0.6);
           //etaBins.push_back(0.8);
           etaBins.push_back(1.05);
           //etaBins.push_back(1.3);
           etaBins.push_back(1.55);
           //etaBins.push_back(1.85);
           etaBins.push_back(2.1);
	}
	etaBins.push_back(2.5);
	const int nEtaBins = etaBins.size()-1;
	std::cout << "No. of eta intervals: " << nEtaBins << std::endl;

	// Centrality bins
	std::vector<double> centralityBins;
	centralityBins.push_back(0.0);
	if(doCentrality){
	  centralityBins.push_back(0.05);
	  centralityBins.push_back(0.10);
	  centralityBins.push_back(0.20);
	  centralityBins.push_back(0.30); // Added to match Z
	  centralityBins.push_back(0.40);
	  centralityBins.push_back(0.60); // Added to match Z
	}
	centralityBins.push_back(0.9);
	const int nCentralityBins = centralityBins.size()-1;
	std::cout << "No. of centrality classes: " << nCentralityBins << std::endl;

    // Output file storing yields in TGraph
    // as a function of pt for each eta/centrality bin
	TString nSigPdfString;

	if(doCentrality&&doCharge&&doEta){
	      nSigPdfString = "CentChrgEta";
	}
	else if(doCentrality&&doCharge&&!doEta){
	      nSigPdfString = "CentChrg";
	}
	else if (doEta&&doCentrality&&!doCharge){
	      nSigPdfString = "CentEta";
	}
	else if (doEta&&doCharge&&!doCentrality){
	      nSigPdfString = "EtaChrg";
	}
	else if (!doEta&&doCharge&&!doCentrality){
	      nSigPdfString = "Charge";
	}
	else{ nSigPdfString = "1000000";}

	TString fileNameOut = "WAnalysis_fitResult"; fileNameOut+=nSigPdfString; 
	if(doPeriodA) fileNameOut+="_PeriodA.root";
	else if(doPeriodB) fileNameOut+="_PeriodB.root";
	else {
        fileNameOut+="_PeriodAB.root";
	}

	std::cout << "W yields will be written to file " << fileNameOut << std::endl;

    TString baseString = "/usatlas/u/tbales/scratch/Data2013Hip/";

	TString fileNameDataIn = "Period";
	if(doPeriodA) fileNameDataIn += "A/HIpSingleMuonHardProbesDataPeriodA.09.18.2013";
	else if(doPeriodB) fileNameDataIn += "B/HIpSingleMuonHardProbesDataPeriodB.09.18.2013";
	else {
        fileNameDataIn = "ntupleW.07.24.2014"; 
        // mu- Period A
/*        fileNameWmin_pp_A = "ntupleWminmunu_pp.A.PowhegPythia8_AU2CT10.09.24.2014";
        fileNameWmin_pn_A = "ntupleWminmunu_pn.A.PowhegPythia8_AU2CT10.09.24.2014";
        // mu- Period B
        fileNameWmin_pp_B = "ntupleWminmunu_pp.B.PowhegPythia8_AU2CT10.09.24.2014";
        fileNameWmin_pn_B = "ntupleWminmunu_pn.B.PowhegPythia8_AU2CT10.09.24.2014";
*/
	}

    TDirectory *dir = gDirectory;
	TFile outFile(fileNameOut,"recreate");
	gDirectory = dir;
	TString resultName = "WPtFit";
	WBookKeeper result(resultName, ptBins, etaBins, centralityBins);

	RooArgSet muonArgSet(muonEta,muonPt,missEt,muonMt,isolationMu,centrality,chargeCategory);
    muonArgSet.add(isPeriodA);
    muonArgSet.add(isPeriodB);
	//muonArgSet.add(weight);

    // Fill (integrated) data set
    RooDataSet* dataSet = fillMuonDataSet(baseString,fileNameDataIn+".root",muonArgSet,doPeriodA,doPeriodB);     
    dataSet->Print(); 	

    // Used for plotting integrated kinematics
	RooBinning bEta = RooBinning(etaBins[0],etaBins[nEtaBins]);
    for(int ieta=0; ieta<nEtaBins; ++ieta){
        bEta.addUniform(1,etaBins[ieta], etaBins[ieta+1]);
    }//ieta
	TH1F* hEtaCumulativePlus=NULL, *hEtaCumulativeMinus=NULL;
    hEtaCumulativePlus  = new TH1F("hEtaCumulativePlus","mu^{+}",bEta.numBoundaries()-1,bEta.array());
    hEtaCumulativeMinus = new TH1F("hEtaCumulativeMinus","mu^{-}",bEta.numBoundaries()-1,bEta.array());

	// Uncomment for weighting ds
	const float nProtons = 82.0;
	const float nNeutrons = 126.0;
	const float nNucleons = nProtons+nNeutrons;
    
    // Open file holding Cw values
    TFile* _fCw = NULL;
    TFile* _fEWBkg = NULL;
    if(doCorrectYields){ 
        // 5 centrality bins
        //_fCw = TFile::Open("/usatlas/u/tbales/code_Hip/CorrectionFactors/efficiencies.09.12.2014.root");
        // 7 centrality bins (matches Z analysis)
        // flipped periods in MC
        //_fCw = TFile::Open("/usatlas/u/tbales/code_Hip/CorrectionFactors/efficiencies.ZBinning.09.15.2014.root");
        // MC has correct period
        _fCw = TFile::Open("/usatlas/u/tbales/code_Hip/CorrectionFactors/efficiencies.09.24.2014.root");
        // fcal in MC re-weighted (note, errors on Cw are normal b/c weighting not allowed in TEfficiency object)
        //_fCw = TFile::Open("/usatlas/u/tbales/code_Hip/CorrectionFactors/efficiencies.fcalReWtd.09.24.2014.root");
        if(!_fCw->IsOpen()) {
            std::cout << "ERROR: Cannot open Cw file. " << std::endl; 
            exit(0);
        }
    }
    if(doBackgroundSub) {
        _fEWBkg = TFile::Open("/usatlas/u/tbales/code_Hip/BackgroundFiles/bkgEW.10.08.2014.root");
        if(!_fEWBkg->IsOpen()) {
            std::cout << "ERROR: Cannot open EW background file. " << std::endl; 
            exit(0);
        }
    }

    // Differential yields
	for(int ipt=0; ipt<nPtBins; ++ipt){
	  for(int ieta=0; ieta<nEtaBins; ++ieta){
	    for(int icent=0; icent<nCentralityBins; ++icent){

            // W^+ --> mu^+
            WAnalysisTool muPlusToolHandle(dataSet);
            muPlusToolHandle.ApplySelectionCuts(cutsSig,+1);
            muPlusToolHandle.SelectPtEtaCentrality(ptBins[ipt],ptBins[ipt+1],etaBins[ieta],etaBins[ieta+1],centralityBins[icent],centralityBins[icent+1]);
            // Save number of signal candidates to outfile
            float sigCandPlus = muPlusToolHandle.numSignalCandidates();
            float sigCandPlusStatErr = TMath::Sqrt(sigCandPlus);

            // W^- --> mu^-
            WAnalysisTool muMinusToolHandle(dataSet);
            muMinusToolHandle.ApplySelectionCuts(cutsSig,-1);
            muMinusToolHandle.SelectPtEtaCentrality(ptBins[ipt],ptBins[ipt+1],etaBins[ieta],etaBins[ieta+1],centralityBins[icent],centralityBins[icent+1]);
            float sigCandMinus = muMinusToolHandle.numSignalCandidates();
            float sigCandMinusStatErr = TMath::Sqrt(sigCandMinus);
            
            // Background subtraction and correction
            if(doBackgroundSub) {
                TGraphAsymmErrors* grBkgHandle = new TGraphAsymmErrors(nEtaBins);
                float NbkgZ_A=0.0,NbkgZ_B=0.0,NbkgTau_A=0.0,NbkgTau_B=0.0;

                std::stringstream sBkgName;
                sBkgName << "grA_muPlus_Z_pt" << ipt << "_cent" << icent;
                grBkgHandle = (TGraphAsymmErrors*)_fEWBkg->Get((sBkgName.str()).c_str());
                NbkgZ_A = grBkgHandle->GetY()[nEtaBins-1-ieta];
                sBkgName.str(std::string());
                sBkgName << "grB_muPlus_Z_pt" << ipt << "_cent" << icent;
                grBkgHandle = (TGraphAsymmErrors*)_fEWBkg->Get((sBkgName.str()).c_str());
                NbkgZ_B = grBkgHandle->GetY()[ieta];

                muPlusToolHandle.setZBackground(NbkgZ_A, NbkgZ_B);

                sBkgName.str(std::string());
                sBkgName << "grA_muMinus_Z_pt" << ipt << "_cent" << icent;
                grBkgHandle = (TGraphAsymmErrors*)_fEWBkg->Get((sBkgName.str()).c_str());
                NbkgZ_A = grBkgHandle->GetY()[nEtaBins-1-ieta];
                sBkgName.str(std::string());
                sBkgName << "grB_muMinus_Z_pt" << ipt << "_cent" << icent;
                grBkgHandle = (TGraphAsymmErrors*)_fEWBkg->Get((sBkgName.str()).c_str());
                NbkgZ_B = grBkgHandle->GetY()[ieta];

                muMinusToolHandle.setZBackground(NbkgZ_A, NbkgZ_B);

                sBkgName.str(std::string());
                sBkgName << "grA_muPlus_tau_pt" << ipt << "_cent" << icent;
                grBkgHandle = (TGraphAsymmErrors*)_fEWBkg->Get((sBkgName.str()).c_str());
                NbkgTau_A = grBkgHandle->GetY()[nEtaBins-1-ieta];
                sBkgName.str(std::string());
                sBkgName << "grB_muPlus_tau_pt" << ipt << "_cent" << icent;
                grBkgHandle = (TGraphAsymmErrors*)_fEWBkg->Get((sBkgName.str()).c_str());
                NbkgTau_B = grBkgHandle->GetY()[ieta];

                muPlusToolHandle.setTauBackground(NbkgTau_A, NbkgTau_B);

                sBkgName.str(std::string());
                sBkgName << "grA_muMinus_tau_pt" << ipt << "_cent" << icent;
                grBkgHandle = (TGraphAsymmErrors*)_fEWBkg->Get((sBkgName.str()).c_str());
                NbkgTau_A = grBkgHandle->GetY()[nEtaBins-1-ieta];
                sBkgName.str(std::string());
                sBkgName << "grB_muMinus_tau_pt" << ipt << "_cent" << icent;
                grBkgHandle = (TGraphAsymmErrors*)_fEWBkg->Get((sBkgName.str()).c_str());
                NbkgTau_B = grBkgHandle->GetY()[ieta];

                muMinusToolHandle.setTauBackground(NbkgTau_A, NbkgTau_B);

                sigCandPlus  = muPlusToolHandle.getBkgSubtractedSignalCandidates();
                sigCandMinus = muMinusToolHandle.getBkgSubtractedSignalCandidates();

                delete grBkgHandle;
            }
            // @ToDo: Add stat errors from Cw as systematic
            float Cw_A=1.0, Cw_B=1.0;
            if(doCorrectYields) {
                TGraphAsymmErrors* grCwHandle = new TGraphAsymmErrors(nEtaBins);
                // mu-
                // Open the file containing Cw factors
                std::stringstream sstrCwName;
                sstrCwName << "effPos_A_wtd_pt" << ipt << "_cent" << icent ; 
                grCwHandle = (TGraphAsymmErrors*)_fCw->Get((sstrCwName.str()).c_str());
                // For period A, eta sign is flipped
                Cw_A = grCwHandle->GetY()[nEtaBins-1-ieta];
                // Clear sstream
                sstrCwName.str(std::string());
                sstrCwName << "effPos_B_wtd_pt" << ipt << "_cent" << icent ; 
                grCwHandle = (TGraphAsymmErrors*)_fCw->Get((sstrCwName.str()).c_str());
                Cw_B = grCwHandle->GetY()[ieta];

                muPlusToolHandle.correctSignalCandidates(Cw_A,Cw_B);

                sigCandPlus = muPlusToolHandle.getNumCorrectedSignalCandidates();
                sigCandPlusStatErr = muPlusToolHandle.getStatErrorCorrectedSignalCandidates();

                // mu-
                // Clear sstream
                sstrCwName.str(std::string());
                sstrCwName << "effNeg_A_wtd_pt" << ipt << "_cent" << icent ; 
                grCwHandle = (TGraphAsymmErrors*)_fCw->Get((sstrCwName.str()).c_str());
                // For period A, eta sign is flipped
                Cw_A = grCwHandle->GetY()[nEtaBins-1-ieta];
                // Clear sstream
                sstrCwName.str(std::string());
                sstrCwName << "effNeg_B_wtd_pt" << ipt << "_cent" << icent ; 
                grCwHandle = (TGraphAsymmErrors*)_fCw->Get((sstrCwName.str()).c_str());
                Cw_B = grCwHandle->GetY()[ieta];

                muMinusToolHandle.correctSignalCandidates(Cw_A,Cw_B);

                sigCandMinus = muMinusToolHandle.getNumCorrectedSignalCandidates();
                sigCandMinusStatErr = muMinusToolHandle.getStatErrorCorrectedSignalCandidates();
                
                delete grCwHandle;
            }

            // mu+
            result.setSig(1,ipt,ieta,icent,sigCandPlus,sigCandPlusStatErr,sigCandPlusStatErr);
            plot(hEtaCumulativePlus,muPlusToolHandle,muonEta,bEta);
            // mu-
            result.setSig(-1,ipt,ieta,icent,sigCandMinus,sigCandMinusStatErr,sigCandMinusStatErr);
            plot(hEtaCumulativeMinus,muMinusToolHandle,muonEta,bEta);

            std::cout << "======================================================================  " << std::endl;
            std::cout << "=============== SUMMARY FOR PT : " << ipt << " ETA: " << ieta << " CENT: " << icent << " ==================" << std::endl;
            std::cout << "======================================================================  " << std::endl;
            std::cout << "================================ mu+ =================================  " << std::endl;
            std::cout << "Signal candidates b/4 correction:                            " << muPlusToolHandle.numSignalCandidates() << std::endl;
            std::cout << "Signal candidates after bkg subtraction:                     " << muPlusToolHandle.getBkgSubtractedSignalCandidates() << std::endl;
            std::cout << "Signal candidates after correction:                          " << sigCandPlus << std::endl;
            std::cout << "Stat. error:                                                 " << sigCandPlusStatErr/sigCandPlus*100 << " [%]" << std::endl;
            std::cout << "================================ mu- =================================  " << std::endl;
            std::cout << "Signal candidates b/4 correction:                            " << muMinusToolHandle.numSignalCandidates() << std::endl;
            std::cout << "Signal candidates after bkg subtraction:                     " << muMinusToolHandle.getBkgSubtractedSignalCandidates() << std::endl;
            std::cout << "Signal candidates after correction:                          " << sigCandMinus << std::endl;
            std::cout << "Stat. error:                                                 " << sigCandMinusStatErr/sigCandMinus*100 << " [%]" << std::endl;

        }//icent
	  }//ieta
	}//ipt

    if(_fCw->IsOpen()){
        _fCw->Close();
    }
	
	TString sEta = format(etaBins[0]); sEta+="<"; sEta+="#eta";  sEta+="<";  sEta+=format(etaBins[nEtaBins]);
	TString sCent = format(centralityBins[0]*100.0); sCent+="-"; sCent+=format(centralityBins[nCentralityBins]*100.0); sCent+="[%]";
	plot(hEtaCumulativePlus,sEta,sCent,"#mu^{+}",+1);
	plot(hEtaCumulativeMinus,sEta,sCent,"#mu^{-}",-1);

	// Save yields to outFile
	result.write(outFile,1);
	result.write(outFile,-1);

}

int main(int argc, char* argv[]){

	// Added to remove funny errors about 
	// dictionary generations
	//gROOT->ProcessLine("#include <vector>");
    std::cout << "==================================================================      " << std::endl;
    std::cout << "=============== STARTING W --> mu + nu ANALYZER    ===============      " << std::endl;
    std::cout << "==================================================================      " << std::endl;
    std::cout << "  This package will count the number of muons in the signal space       " << std::endl;
    std::cout << "define by the user and save these muons to an output file for           " << std::endl;
    std::cout << "later use. Control kinematic histograms -- provided by the user --      " << std::endl;
    std::cout << "are also generated and stored in the output file                        " << std::endl;
    std::cout << "                                                                        " << std::endl;
    std::cout << "  Provide as an arguments to this macro:                                " << std::endl;
    std::cout << "1) Period (A,B,AB) 2) Subtract Background (Y/N) 3) Correct Yields (Y/N) " << std::endl;
    std::cout << "==================================================================" << std::endl;
    std::time_t now = time(0);
    char* dt = std::ctime(&now);
    std::cout << dt << std::endl; 

    bool doBackgroundSub=false,doCorrectYields=false;
    std::string sPeriod = "";
    if(argc==2){
        if(strcmp(argv[1],"A")==0) sPeriod="A";
        else if(strcmp(argv[1],"B")==0) sPeriod="B";
        if(strcmp(argv[1],"AB")==0) sPeriod="AB";
	    WAnalysisHip(sPeriod);
    }
    else if(argc==3){
        if(strcmp(argv[2],"Y")==0) doBackgroundSub=true;
        WAnalysisHip(argv[1],doBackgroundSub);
    }
    else if(argc==4) {
        if(strcmp(argv[2],"Y")==0) doBackgroundSub=true;
        if(strcmp(argv[3],"Y")==0) doCorrectYields=true;
        WAnalysisHip(argv[1],doBackgroundSub,doCorrectYields);
    }
    else std::cout << "Please provide A, B, or AB as argument." << std::endl;
}
