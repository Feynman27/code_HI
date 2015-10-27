#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath> 
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
#include "FitResult.cc"
#include "AtlasUtils.C"
#include "AtlasStyle.C"

using namespace RooFit;

RooDataSet* fillMuonDataSet(const TString& pathName, const TString& fileName, RooArgSet& muonArgSet, bool isMc){

	TChain* tree = new TChain("data","data");
	int nFiles = tree->Add(pathName+fileName);
        std::cout << "Filling the RooDataSet for " << fileName << "... Number of files: " << nFiles << std::endl;
	std::cout << "Number of entries in data tree: " << tree->GetEntries() << std::endl;

	RooDataSet* set = new RooDataSet(fileName,fileName,muonArgSet);
	std::vector<float>* vmu_pt= 0; 
	std::vector<float>* vmu_mt= 0;
	std::vector<float>* vmu_eta= 0;
	std::vector<float>* vmu_ptcone20= 0;
	std::vector<float>* vmu_charge = 0;
	std::vector<bool>* vEF_mu8_matched = 0; 
	std::vector<bool>* visZCandidate = 0;
	std::vector<int>* vmu_quality = 0;
	float centrality, MET_RefFinal_et;

	tree->SetBranchAddress("mu_pt",&vmu_pt);
	tree->SetBranchAddress("mu_mt",&vmu_mt);
	tree->SetBranchAddress("mu_eta",&vmu_eta);
	tree->SetBranchAddress("mu_charge",&vmu_charge);
	tree->SetBranchAddress("mu_ptcone20",&vmu_ptcone20);
	tree->SetBranchAddress("EF_mu8_matched",&vEF_mu8_matched);
	tree->SetBranchAddress("isZCandidate",&visZCandidate);
	tree->SetBranchAddress("mu_quality",&vmu_quality);
	tree->SetBranchAddress("centrality",&centrality);
	tree->SetBranchAddress("MET_RefFinal_et",&MET_RefFinal_et);
	// --- Set branch status ---
  	tree->SetBranchStatus("*",0) ;
	tree->SetBranchStatus("mu_pt",1);
	tree->SetBranchStatus("mu_mt",1);
	tree->SetBranchStatus("mu_eta",1);
	tree->SetBranchStatus("mu_charge",1);
	tree->SetBranchStatus("mu_ptcone20",1);
	tree->SetBranchStatus("EF_mu8_matched",1);
	tree->SetBranchStatus("isZCandidate",1);
	tree->SetBranchStatus("mu_quality",1);
	tree->SetBranchStatus("centrality",1);
	tree->SetBranchStatus("MET_RefFinal_et",1);
	
	for(int iev = 0; iev < tree->GetEntries(); ++iev){
           if(iev%10000==0) std::cout << "Event: " << iev << std::endl;
           tree->LoadTree(iev);
           tree->GetEntry(iev);

	  int nmu = vmu_charge->size();
	  
	  for(int imu =0; imu < nmu; ++imu){
	     if(vEF_mu8_matched->at(imu)!=1||visZCandidate->at(imu)==1||vmu_quality->at(imu)<11) continue;	
	      muonArgSet.setRealValue("centrality",centrality);
	      muonArgSet.setRealValue("missEt",MET_RefFinal_et);
	      muonArgSet.setRealValue("muonPt",vmu_pt->at(imu));
	      muonArgSet.setRealValue("muonMt",vmu_mt->at(imu));
              if(vmu_charge->at(imu)<0.0) muonArgSet.setCatLabel("chargeCategory","muMinus");
	      else if(vmu_charge->at(imu)>0.0) muonArgSet.setCatLabel("chargeCategory","muPlus"); 
	      float isoTemp = vmu_ptcone20->at(imu)/vmu_pt->at(imu);
	      muonArgSet.setRealValue("isolationMu",isoTemp);
	      muonArgSet.setRealValue("muonEta",vmu_eta->at(imu));
//	      muonArgSet.setRealValue("w",wt);
	      set->add(muonArgSet);
	  }//imu
	}//iev
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



void plot(WAnalysisTool& globalToolHandle, TString sEta, TString sCent, TString sCharge, int charge){

	TH1F* hDataSet = globalToolHandle.GetHistogram();
	TCanvas *canvas = new TCanvas("canvas","canvas",600,600);	
	hDataSet->Draw("pe");
	TLatex* texEta = new TLatex(0.66,0.86,sEta);
	TLatex* texCent = new TLatex(0.72,0.79,sCent);
	texEta->SetNDC();
	texCent->SetNDC();
	texEta->Draw();
	texCent->Draw();
	TString sName = "WAnalysisHip_Mu"; sName+=charge;  
	canvas->Print(sName+".root");
}

void plot(WAnalysisTool& globalToolHandle, WAnalysisTool& toolHandle, RooRealVar& var, double cutLow, double cutUpp){
	
	RooBinning bEta = RooBinning(cutLow,cutUpp);
	bEta.addUniform(1,cutLow,-2.1);	
	bEta.addUniform(1,-2.1,-1.55);	
	bEta.addUniform(1,-1.55,-1.05);	
	bEta.addUniform(1,-1.05,-0.6);	
	bEta.addUniform(1,-0.6,-0.1);	
	bEta.addUniform(1,-0.1,0.1);	
	bEta.addUniform(1,0.1,0.6);	
	bEta.addUniform(1,0.6,1.05);	
	bEta.addUniform(1,1.05,1.55);	
	bEta.addUniform(1,1.55,2.1);	
	bEta.addUniform(1,2.1,cutUpp);	

	toolHandle.FillHistogram(bEta,var,cutLow,cutUpp);
	TH1F* hDataSet = toolHandle.GetHistogram();
	//hDataSet->Scale(1.0,"width");
	globalToolHandle.AddToStack(hDataSet);
	//TCanvas *canvas = toolHandle.plotOnCanvas(var);
	//canvas->Print("test.root");
}

TString format(float value) {
  std::stringstream svalue;
  svalue  << std::setprecision(2) << value;
  return svalue.str();
}

void WAnalysisHip(const std::string& period){

	bool doPeriodA = false;
	bool doPeriodB = false; 
	if(period.compare("A") == 0) {
		doPeriodA = true;
		std::cout << "*****************************" << std::endl;
		std::cout << "**Running on Period A data.**" << std::endl;
		std::cout << "*****************************" << std::endl;
	}
	else if(period.compare("B") == 0) {
		doPeriodB = true;
		std::cout << "*****************************" << std::endl;
		std::cout << "**Running on Period B data.**" << std::endl;
		std::cout << "*****************************" << std::endl;
	}
	else {
		std::cout << "Please choose Period A or B. " << std::endl;
		exit(EXIT_FAILURE);
	}
	bool doPt = false;
	bool doEta = true;
	bool doCentrality = true;
	bool doCharge = true;
	// Only need to build templates once
	bool doBuildTemplates = false;

	// Declare variables to cut on
	RooRealVar  muonPt("muonPt","p_{T}",0.0,400.0,"GeV");
	RooRealVar  missEt("missEt","E_{T}^{miss}",0.0,400.0,"GeV");
	RooRealVar  muonMt("muonMt","m_{T}",0.0,400.0,"GeV");
	RooRealVar  muonCharge("muonCharge","charge",-2.0,+2.0);
	RooRealVar  isolationMu("isolationMu","isolationMu",0.0,10.0);
	RooRealVar  centrality("centrality","centrality",0.,1.0);
	RooRealVar  muonEta("muonEta","#eta^{#mu}",-3.0,+3.0);
//	RooRealVar  weight("weight","weight",0.0,10.0);
	
	TString sCutsSig = "abs(muonEta)<2.5&&muonPt>25.0&&missEt>25.0&&isolationMu<0.2&&muonMt>40.0&&centrality>0.0";
	TString sMcWCutsSig = "abs(muonEta)<2.5&&muonPt>25.0&&missEt>25.0&&muonMt>40.0";
	// Test
	//TString sCutsSig = "abs(muonEta)<2.5&&muonPt>25.0&&missEt>0.0&&isolationMu<0.2&&muonMt>40.0&&centrality>0.0";
	RooArgList muonArgList(muonEta,muonPt);
	muonArgList.add(missEt);
	muonArgList.add(isolationMu);
	muonArgList.add(muonMt);
	muonArgList.add(centrality);

	RooArgList muonMcWArgList(muonEta,muonPt,missEt,muonMt);

	RooFormulaVar cutsSig("cutsSig", "cutsSig", sCutsSig, muonArgList);
	RooFormulaVar cutsMcWSig("cutsMcWSig", "cutsMcWSig", sMcWCutsSig, muonMcWArgList);

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
	std::cout << "nEtaBins:" << nEtaBins << std::endl;
	// Centrality bins
	std::vector<double> centralityBins;
	centralityBins.push_back(0.0);
	if(doCentrality){
	  centralityBins.push_back(0.05);
	  centralityBins.push_back(0.10);
	  centralityBins.push_back(0.20);
	  centralityBins.push_back(0.40);
	}
	centralityBins.push_back(0.9);
	const int nCentralityBins = centralityBins.size()-1;

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
	  std::cout << "Must choose between Period A or B." << std::endl; 
	  exit(EXIT_FAILURE);
	}
	std::cout << "W yields will be written to file " << fileNameOut << std::endl;

	TString fileNamePdf = "WAnalysisHip_templates"; 
        TString baseString = "/usatlas/u/tbales/scratch/Data2013Hip/";

	TString fileNameDataIn = "Period";
	if(doPeriodA) fileNameDataIn += "A/HIpSingleMuonHardProbesDataPeriodA.09.18.2013";
	else if(doPeriodB) fileNameDataIn += "B/HIpSingleMuonHardProbesDataPeriodB.09.18.2013";
	else {
	  std::cout << "Must choose between Period A or B." << std::endl; 
	  exit(EXIT_FAILURE);
	}

	TString baseString2 = "/usatlas/u/tbales/MonteCarloGenerators/pythia8176/RunArea/";
	// W Mc
	TString fileNameMcWPlusppIn  = "pp/W+/PowhegPythia8_MSTW2008NLO_W+_pp.08.29.2013";	
	TString fileNameMcWMinusppIn = "pp/W-/PowhegPythia8_MSTW2008NLO_W-_pp.08.29.2013";	
	TString fileNameMcWPlusnpIn  = "np/W+/PowhegPythia8_MSTW2008NLO_W+_np.08.29.2013"; // Period A	
	TString fileNameMcWMinusnpIn = "np/W-/PowhegPythia8_MSTW2008NLO_W-_np.08.29.2013"; // Period A	
	TString fileNameMcWPluspnIn  = "pn/W+/PowhegPythia8_MSTW2008NLO_W+_pn.08.29.2013"; // Period B	
	TString fileNameMcWMinuspnIn = "pn/W-/PowhegPythia8_MSTW2008NLO_W-_pn.08.29.2013"; // Period B	

        TDirectory *dir = gDirectory;
	TFile outFile(fileNameOut,"recreate");
	gDirectory = dir;
	TString resultName = "WPtFit";
	FitResult result(resultName, ptBins, etaBins, centralityBins);

	RooArgSet muonArgSet(muonPt,missEt,isolationMu,muonMt,muonEta,centrality,chargeCategory);
	RooArgSet muonMcWArgSet(muonPt,missEt,muonMt,muonEta);
	//muonArgSet.add(weight);

        // Fill data set
        RooDataSet* dataSet       = fillMuonDataSet(baseString,fileNameDataIn+".root",muonArgSet,false); dataSet->Print(); 	
        RooDataSet* mcWPlusppSet  = fillMuonMcSet(baseString2,fileNameMcWPlusppIn+".root",muonMcWArgSet);   mcWPlusppSet->Print(); 	
        RooDataSet* mcWPlusnpSet  = fillMuonMcSet(baseString2,fileNameMcWPlusnpIn+".root",muonMcWArgSet);   mcWPlusnpSet->Print(); 	
        RooDataSet* mcWPluspnSet  = fillMuonMcSet(baseString2,fileNameMcWPluspnIn+".root",muonMcWArgSet);   mcWPluspnSet->Print(); 	
        RooDataSet* mcWMinusppSet = fillMuonMcSet(baseString2,fileNameMcWMinusppIn+".root",muonMcWArgSet);   mcWMinusppSet->Print(); 	
        RooDataSet* mcWMinusnpSet = fillMuonMcSet(baseString2,fileNameMcWMinusnpIn+".root",muonMcWArgSet);   mcWMinusnpSet->Print(); 	
        RooDataSet* mcWMinuspnSet = fillMuonMcSet(baseString2,fileNameMcWMinuspnIn+".root",muonMcWArgSet);   mcWMinuspnSet->Print(); 	

	WAnalysisTool muPlusGlobalToolHandle(dataSet);
	WAnalysisTool muMinusGlobalToolHandle(dataSet);
	// W MC (pp, np, pn)
	WAnalysisTool muMcWppPlusGlobalToolHandle(mcWPlusppSet);
	WAnalysisTool muMcWnpPlusGlobalToolHandle(mcWPlusnpSet);
	WAnalysisTool muMcWpnPlusGlobalToolHandle(mcWPluspnSet);

	WAnalysisTool muMcWppMinusGlobalToolHandle(mcWMinusppSet);
	WAnalysisTool muMcWnpMinusGlobalToolHandle(mcWMinusnpSet);
	WAnalysisTool muMcWpnMinusGlobalToolHandle(mcWMinuspnSet);

	// Uncomment for weighting ds
	const float nProtons = 82.0;
	const float nNeutrons = 126.0;
	const float nNucleons = nProtons+nNeutrons;

	//weight d.s.
	// CHECK WEIGHTS
	muMcWppPlusGlobalToolHandle.weightDS(nProtons/nNucleons);
	muMcWnpPlusGlobalToolHandle.weightDS(nNeutrons/nNucleons);
	muMcWpnPlusGlobalToolHandle.weightDS(nNeutrons/nNucleons);
	muMcWppMinusGlobalToolHandle.weightDS(nProtons/nNucleons);
	muMcWnpMinusGlobalToolHandle.weightDS(nNeutrons/nNucleons);
	muMcWpnMinusGlobalToolHandle.weightDS(nNeutrons/nNucleons);

	// Test with Period A
	RooDataSet* mcWPlusSet = muMcWppPlusGlobalToolHandle.GetDataSet(); 
	RooDataSet* mcWMinusSet = muMcWppMinusGlobalToolHandle.GetDataSet(); 
	if(doPeriodA) {
		mcWPlusSet->append(*muMcWnpPlusGlobalToolHandle.GetDataSet());
	        mcWMinusSet->append(*muMcWnpMinusGlobalToolHandle.GetDataSet());
	}
	else if(doPeriodB){
		mcWPlusSet->append(*muMcWpnPlusGlobalToolHandle.GetDataSet());
	        mcWMinusSet->append(*muMcWpnMinusGlobalToolHandle.GetDataSet());
	}
	else{
	  std::cout << "Must choose Period A or B MC dataset." << std::endl;	
	  exit(EXIT_FAILURE); 
	}

	// Release intermediate Mc datasets
	delete mcWPlusppSet;
	delete mcWPlusnpSet;
	delete mcWPluspnSet;
	delete mcWMinusppSet;
	delete mcWMinusnpSet;
	delete mcWMinuspnSet;

	// merge pp,np(pn) 
	if (doBuildTemplates){
/*
          WAnalysisTool mcWPlusToolHandle(mcWPlusSet);
	  WAnalysisTool mcWMinusToolHandle(mcWMinusSet);
	  mcWPlusToolHandle.ApplySelectionCuts(cutsMcWSig);
	  mcWMinusToolHandle.ApplySelectionCuts(cutsMcWSig);
          mcWPlusToolHandle.buildTemplate(muonEta,fileNamePdf+"_WPlusEta.root");   
          mcWMinusToolHandle.buildTemplate(muonEta,fileNamePdf+"_WMinusEta.root");   
*/
	  // Move dataset to signal phase space 
	  muMcWppPlusGlobalToolHandle.ApplySelectionCuts(cutsMcWSig);
	  muMcWnpPlusGlobalToolHandle.ApplySelectionCuts(cutsMcWSig);
	  muMcWpnPlusGlobalToolHandle.ApplySelectionCuts(cutsMcWSig);

	  muMcWppMinusGlobalToolHandle.ApplySelectionCuts(cutsMcWSig);
	  muMcWnpMinusGlobalToolHandle.ApplySelectionCuts(cutsMcWSig);
	  muMcWpnMinusGlobalToolHandle.ApplySelectionCuts(cutsMcWSig);

          // Build Pdfs
          muMcWppPlusGlobalToolHandle.buildTemplate(muonEta,fileNamePdf+"_pp_WPlusEta.root");
          muMcWnpPlusGlobalToolHandle.buildTemplate(muonEta,fileNamePdf+"_np_WPlusEta.root");
          muMcWpnPlusGlobalToolHandle.buildTemplate(muonEta,fileNamePdf+"_pn_WPlusEta.root");

          muMcWppMinusGlobalToolHandle.buildTemplate(muonEta,fileNamePdf+"_pp_WMinusEta.root");
          muMcWnpMinusGlobalToolHandle.buildTemplate(muonEta,fileNamePdf+"_np_WMinusEta.root");
          muMcWpnMinusGlobalToolHandle.buildTemplate(muonEta,fileNamePdf+"_pn_WMinusEta.root");

	  std::cout << "RooKeysPdfs built. Please re-run macro with doBuildTemplates=FALSE." << std::endl;
	  exit(EXIT_SUCCESS); 
	 
	}

/*	bEta.addUniform(1,etaBins[0],-2.1);	
	bEta.addUniform(1,-2.1,-1.55);	
	bEta.addUniform(1,-1.55,-1.05);	
	bEta.addUniform(1,-1.05,-0.6);	
	bEta.addUniform(1,-0.6,-0.1);	
	bEta.addUniform(1,-0.1,+0.1);	
	bEta.addUniform(1,0.1,0.6);	
	bEta.addUniform(1,0.6,1.55);	
	bEta.addUniform(1,1.55,2.1);	
	bEta.addUniform(1,2.1,etaBins[nEtaBins]);	
*/
	//double xEtaBins[] = {etaBins[0],-2.1,-1.55,-1.05,-0.6,-0.1,0.1,0.6,1.05,1.55,2.1,etaBins[nEtaBins]};
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
		result.setSig(1,ipt,ieta,icent,sigCandPlus,sigCandPlusStatErr,sigCandPlusStatErr);
		plot(muPlusGlobalToolHandle,muPlusToolHandle,muonEta,-2.50,2.50);
	
		// W^- --> mu^-
		WAnalysisTool muMinusToolHandle(dataSet);
		muMinusToolHandle.ApplySelectionCuts(cutsSig,-1);
		muMinusToolHandle.SelectPtEtaCentrality(ptBins[ipt],ptBins[ipt+1],etaBins[ieta],etaBins[ieta+1],centralityBins[icent],centralityBins[icent+1]);
		float sigCandMinus = muMinusToolHandle.numSignalCandidates();
		float sigCandMinusStatErr = TMath::Sqrt(sigCandMinus);
		result.setSig(-1,ipt,ieta,icent,sigCandMinus,sigCandMinusStatErr,sigCandMinusStatErr);
		plot(muMinusGlobalToolHandle,muMinusToolHandle,muonEta,-2.50,2.50);

            }//icent
	  }//ieta
	}//ipt
	
	TString sEta = format(etaBins[0]); sEta+="<"; sEta+="#eta";  sEta+="<";  sEta+=format(etaBins[nEtaBins]);
	TString sCent = format(centralityBins[0]*100.0); sCent+="-"; sCent+=format(centralityBins[nCentralityBins]*100.0); sCent+="[%]";
	plot(muPlusGlobalToolHandle,sEta,sCent,"#mu^{+}",+1);
	plot(muMinusGlobalToolHandle,sEta,sCent,"#mu^{-}",-1);

	// Save yields to outFile
	result.write(outFile,1);
	result.write(outFile,-1);

}

int main(){

	// Added to remove funny errors about 
	// dictionary generations
	gROOT->ProcessLine("#include <vector>");
	WAnalysisHip("A");
	WAnalysisHip("B");
}
