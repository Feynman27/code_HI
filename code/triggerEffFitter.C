#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooPolyVar.h"
#include "RooEfficiency.h"
#include "RooPolynomial.h"
#include "RooCategory.h"
#include "RooArgSet.h"
#include "RooBinning.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"

#include "TCanvas.h"
#include "TAxis.h"
#include "TString.h"
#include "TTree.h"
#include "TMath.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TFile.h"
#include "TH1F.h"
#include "TEfficiency.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <sstream>
#include <string>
//#ifdef __MAKECINT__
//#pragma link C++ class std::vector<std::vector<int> >+;
//#pragma link C++ class std::vector<std::vector<float> >+;
//#endif
using namespace RooFit ;

void Write(TFile* const outFile, TH1F* h,char name[]){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing efficiencies to root file..." << std::endl;
	  h->Write(name);
	  std::cout << "Done." << std::endl;
      gDirectory = dir;
}

TEfficiency* writeEfficiency(TFile* fEff, const TH1& hPassed, const TH1& hTotal, char sEff[]){
	if(TEfficiency::CheckConsistency(hPassed,hTotal,"w")) {

		TEfficiency* pEff = 0;
		fEff->cd();
		pEff = new TEfficiency(hPassed,hTotal);
		pEff->Write(sEff);
		return pEff;
		//pEff->Draw("AP");
	}
    else return 0;
}



TString format(float value) {
  std::stringstream svalue;
  svalue  << std::setprecision(3) << value;
  return svalue.str();
}

////////////////////////////////////////////////////////
//selectPtEta
////////////////////////////////////////////////////////
RooDataSet* selectPtEta( RooDataSet* data, double ptLow, double ptHigh, double etaLow, double etaHigh, bool doMirrorEta = false)
{
  TString cut = "muonPt>";
  cut += ptLow;
  cut += "&&muonPt<";
  cut += ptHigh;
  cut += "&&((muonEta>";
  cut += etaLow;
  cut += "&&muonEta<";
  cut += etaHigh;
  cut += ")";
  if ( doMirrorEta){
    cut += "||(muonEta>"; 
    cut += -etaHigh;
    cut += "&&muonEta<";
    cut += -etaLow;
    cut += ")";
  }
  cut += ")";

  return (RooDataSet*) data->reduce(cut);
}
/////////////////////////////////////////////////////////
//selectPtEtaCentrality
/////////////////////////////////////////////////////////
RooDataSet* selectPtEtaCentrality( RooDataSet* data, double ptLow, double ptHigh, double etaLow, double etaHigh, double centralityLow, double centralityHigh,bool doMirrorEta = false)
{
  data = selectPtEta( data, ptLow, ptHigh, etaLow, etaHigh, doMirrorEta); 
  TString cut = "centrality>";
  cut += centralityLow;
  cut += "&&centrality<";
  cut += centralityHigh;

  return (RooDataSet*) data->reduce(cut);
}

///////////////////////////////////////////////////////////
//fillDataSet
///////////////////////////////////////////////////////////
RooDataSet* fillDataSet(const TString& fileName, RooArgSet& muonArgSet, int cutValue = 11, bool isMC = false)
{
  RooDataSet* data = new RooDataSet(fileName,fileName,muonArgSet) ;
  
  float eLossNt[50];
  float scatNt[50];
  float compNt[50];
  double ptSmNt[50];
  float ptNt[50],mtNt[50];
  float etaNt[50];
  float chargeNt[50];
  int promptNt[50];
  float centralityNt;
  float ptcone20Nt[50];
  int valNt[50], ZDYNt[50], efMatched[50];
  int nmu,trig1,trig2,trig3,trig4,trig5,L1_MU0,L1Matched[50],mbTrig1,mbTrig2;
  int matched1[50], matched2[50], matched3[50];
  float nu_ptNt;

  TChain* tree = new TChain("tree","tree");
  int nFiles = tree->Add(fileName);
  std::cout <<"Filling the RooDataSet for "<< fileName << "...Number of files: " << nFiles << std::endl;

    // --- Set branch adresses ---
  // --- Set branch adresses ---
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC",&trig1);
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10",&trig2);
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20",&trig3);
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC_Matched20",&matched1);
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10_Matched20",&matched2);
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20_Matched20",&matched3);
  tree->SetBranchAddress("EF_mbZdc_a_c_L1VTE50_trk",&mbTrig1);
  tree->SetBranchAddress("EF_L1TE50_NoAlg",&mbTrig2);
  tree->SetBranchAddress("nu_pt", &nu_ptNt);
//  tree->SetBranchAddress("EF_mu4_MSonly_L1TE50",&trig4);
//  tree->SetBranchAddress("EF_mu4_L1VTE50",&trig5);
  //tree->SetBranchAddress("EF_mu8",&trig1);
  //tree->SetBranchAddress("EF_mu8_Matched20",&trigMatched1);
  //tree->SetBranchAddress("L1_MU0",&L1_MU0);
  //tree->SetBranchAddress("L1Matched",&L1Matched);
  tree->SetBranchAddress("eLoss", &eLossNt);
  //tree->SetBranchAddress("mujetdR", &dRNt);
  tree->SetBranchAddress("ptcone20", &ptcone20Nt);
  tree->SetBranchAddress("scat", &scatNt);
  tree->SetBranchAddress("comp", &compNt);
  // --- MC_pt_smearing --- //
  //if(isMC) tree->SetBranchAddress("pTCB_smeared", &ptSmNt);
  tree->SetBranchAddress("pt", &ptNt);
  tree->SetBranchAddress("mt", &mtNt);
  tree->SetBranchAddress("eta", &etaNt);
  tree->SetBranchAddress("charge", &chargeNt);
  if(isMC) tree->SetBranchAddress("prompt", &promptNt);
  tree->SetBranchAddress("val", &valNt); 
  tree->SetBranchAddress("ZDY", &ZDYNt); 
  tree->SetBranchAddress("centrality", &centralityNt);
  tree->SetBranchAddress("mu_muid_n", &nmu);

   // --- Set branch status ---
  tree->SetBranchStatus("*",0) ;
  tree->SetBranchStatus("mu_muid_n", 1);
  tree->SetBranchStatus("eLoss", 1);
  tree->SetBranchStatus("mujetdR", 1);
  tree->SetBranchStatus("ptcone20", 1);
  tree->SetBranchStatus("nu_pt", 1);
  tree->SetBranchStatus("scat", 1);
  tree->SetBranchStatus("comp", 1);
  tree->SetBranchStatus("pt", 1);
  tree->SetBranchStatus("mt", 1);
  tree->SetBranchStatus("eta", 1);
  tree->SetBranchStatus("charge", 1);
  if(isMC) tree->SetBranchStatus("prompt", 1);
  tree->SetBranchStatus("val", 1); 
  tree->SetBranchStatus("ZDY", 1); 
  tree->SetBranchStatus("centrality", 1);
  //tree->SetBranchStatus("EF_mu10",1);
  //tree->SetBranchStatus("EF_mu10_Matched20",1);
  tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC", 1); 
  tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10", 1); 
  tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20", 1); 
  tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC_Matched20", 1); 
  tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10_Matched20", 1); 
  tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20_Matched20", 1); 
  tree->SetBranchStatus("EF_mbZdc_a_c_L1VTE50_trk",1);
  tree->SetBranchStatus("EF_L1TE50_NoAlg",1);
  
  std::cout << "Number of entries: " << tree->GetEntries() << std::endl;

  //Add to the dataset high quality muon pT values
  //accept muons that fired one of the specified triggers
  for (int i = 0; i < tree->GetEntries(); i++){
     //if(i==10000000) break; //temp hack
     tree->GetEntry(i); 
    if( !(mbTrig1)&&!(mbTrig2) ) continue;
     muonArgSet.setRealValue("centrality",centralityNt) ;
    for (int imu = 0; imu<nmu;imu++){
     if (
           valNt[imu]>11
           && abs(scatNt[imu])<4.0
           && abs(eLossNt[imu])<0.5
           && abs(etaNt[imu])>0.1
           && abs(etaNt[imu])<2.4
           && centralityNt>=0.
           && centralityNt<=0.8
           &&ptNt[imu]>0.0
           /*&&ptcone20Nt[imu]/ptNt[imu]<0.1
           && ZDYNt[imu]==0
           &&nu_ptNt >25.0
           &&ptNt[imu]>25.0
           &&mtNt[imu]>40.0
            */
        ) {

       muonArgSet.setRealValue("muonPt",ptNt[imu]); 
       muonArgSet.setRealValue("muonEta",etaNt[imu]); 
       muonArgSet.setRealValue("muonCharge",chargeNt[imu]);
       muonArgSet.setRealValue("muonELoss",eLossNt[imu]);
       muonArgSet.setRealValue("muonScat",scatNt[imu]);
       muonArgSet.setRealValue("composite",compNt[imu]);
       if(isMC)muonArgSet.setRealValue("muonCategory",promptNt[imu]);
       //muonArgSet.setRealValue("muonEFTrig",trigEFNt[imu]) ;
       //muonArgSet.setRealValue("muonEFTrigMatch",matchEFNt[imu]) ;
       
       if ( chargeNt[imu] > 0 ) muonArgSet.setCatLabel("chargeCategory","muPlus") ;
       else if ( chargeNt[imu] < 0 ) muonArgSet.setCatLabel("chargeCategory","muMinus") ;
       if ( (trig1&&matched1[imu])||(trig2&&matched2[imu])||(trig3&&matched3[imu])  ) muonArgSet.setCatLabel("cut","accept") ;
       else 
         {
           muonArgSet.setCatLabel("cut","reject")  ;
         }

       data->add(muonArgSet) ;
     }
    } //imu

  } //events
    return data ;
}

void triggerEffFitter(){
  //set the quality and variables for the trigger fitting
  int cutValue = 11;
  double ptmax = 100.0;
  double ptcutLow = 0.0;
  //one cut at a time at the moment
  bool doEta = false;
  bool doCentrality = false ;
  bool doCharge = true;
  bool doSystematics = false;

  RooRealVar muonPt("muonPt","muonPt",0.,300.) ;
  RooRealVar centrality("centrality","centrality",0.0,1.0) ;
  RooRealVar muonEta("muonEta","muonEta",-3.0,+3.0);
  RooRealVar  muonCharge("muonCharge","charge",-2.0,+2.0);
  RooRealVar  muonELoss("muonELoss","muonELoss",-1.5,+1.5);
  RooRealVar  muonScat("muonScat","muonScat",-5.0,+5.0);

/*    b.addUniform(1,kinVarLow,-2.1);
    b.addUniform(1,-2.1,-1.85);
    b.addUniform(1,-1.85,-1.55);
    b.addUniform(1,-1.55,-1.3);
    b.addUniform(1,-1.3,-1.05);
    b.addUniform(1,-1.05,-0.8);
    b.addUniform(1,-0.8,-0.6);
    b.addUniform(1,-0.6,-0.35);
    b.addUniform(1,-0.35,-0.1);
    b.addUniform(1,-0.1,0.1);
    b.addUniform(1,0.1,0.35);
*/
    RooBinning bEta(0.1,2.4) ;  
    bEta.addUniform(1,0.1,0.6);
    bEta.addUniform(1,0.6,1.05);
    bEta.addUniform(1,1.05,1.55);
    bEta.addUniform(1,1.55,2.4);

  //variable binning in pT regions
  RooBinning b0(0,100) ;  //save
  RooBinning b1(0,100) ;  
  b1.addUniform(100,0,100) ;
  b0.addUniform(18,0,9) ;  //save
  b0.addUniform(11,9,20) ; 
  b0.addUniform(10,20,40);
  b0.addUniform(4,40,60); //save
  b0.addUniform(1,60,100); //save

  RooRealVar a1("a1","a1",0.48,0.3,0.5) ;
  RooRealVar a2("a2","a2",12.,5.0,17.0) ;
  RooRealVar a3("a3","a3",10.0,0.,15.) ;
  RooRealVar a4("a4","a4",0.0,-5.0,5.0) ;
  RooRealVar a5("a5","a5",1.0,1.0e-2,4.0) ;

  //RooFormulaVar effFunc("effFunc","a4+(a1-a4)/(TMath::Exp((a2-muonPt)/a3)+1)",RooArgList(a4,a1,a2,muonPt,a3)) ;
  RooFormulaVar effFunc("effFunc","a1*( 1.0+(TMath::Erf( (muonPt-a2)/a5 )) )",RooArgList(a1,muonPt,a2,a5)) ;
  //RooFormulaVar effFunc("effFunc","(a1)/(TMath::Exp((a2-muonPt)/a3)+1)",RooArgList(a1,a2,muonPt,a3)) ;

  TFile* outFile = new TFile("triggerEfficiencies.root","recreate");
  ///Data2011 from HIDPD ntuples
  //TString fileName = "/tmp/tbalestr/HISingleMuonMB.2012.04.10.v02.root";
  /// Weizmann
  //TString fileName = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/MinimumBiasFiles/HISingleMuonMinBiasData.07.10.2013.root";
  /// BNL
  TString fileName = "/usatlas/u/tbales/scratch/data_files/HISingleMuonMinBiasData.07.10.2013.root";
  //TString fileName = "root://eosatlas//eos/atlas/user/t/tbalestr/HISingleMuon_DPD_MB_30May.root";
  //TString fileName = "root://eosatlas//eos/atlas/user/t/tbalestr/HISingleMuon_MB_10May.root";
  //TString fileName = "root://eosatlas//eos/atlas/user/t/tbalestr/HISingleMuon_notrig_11May.root";
  ///MC2011 Wmunu PYTHIA+HIJING
  //TString fileName = "root://eosatlas//eos/atlas/user/t/tbalestr/HISingleMuon_AllBeams_2Apr.root";
  //TString fileName = "root://eosatlas//eos/atlas/user/t/tbalestr/HISingleMuon_DPD_MB_30May.root";
  //TString fileName = "root://eosatlas//eos/atlas/user/t/tbalestr/HISingleMuon_MB_10May.root";
  //TString fileName = "root://eosatlas//eos/atlas/user/t/tbalestr/HISingleMuon_notrig_11May.root";
  ///MC2011 Wmunu PYTHIA+HIJING
  //TString fileName = "root://eosatlas//eos/atlas/user/t/tbalestr/HISingleMuon_AllBeams_2Apr.root";

  //category for separating those that fired trigger and those that did not
  RooCategory cut("cut","cut") ;
  cut.defineType("accept",1) ;
  cut.defineType("reject",0) ;

  RooCategory chargeCategory("chargeCategory","chargeCategory") ;
  chargeCategory.defineType("muMinus",-1) ;
  chargeCategory.defineType("muPlus",+1) ;

  RooArgSet muonArgSet(muonPt,muonEta,centrality,muonCharge,cut,chargeCategory) ;
  RooEfficiency effPdf("effPdf","effPdf",effFunc,cut,"accept") ;

  //create the dataset
  RooDataSet* data = fillDataSet(fileName,muonArgSet,cutValue, false);

  // --- Set pT, eta, and centrality bins -- //
  std::vector<double> ptBins;
  ptBins.push_back(0.0);
  ptBins.push_back(ptmax);
  const int nPtBins = ptBins.size()-1;
  std::cout << "Number of pT bins: " << nPtBins << std::endl;

  std::vector<double> etaBins;
  etaBins.push_back(0.10);
  if (doEta){
        //etaBins.push_back(0.35);
        etaBins.push_back(0.6);
        //etaBins.push_back(0.8);
        etaBins.push_back(1.05);
        //etaBins.push_back(1.3);
        etaBins.push_back(1.55);
        //etaBins.push_back(1.85);
        //etaBins.push_back(2.1);
  }
  etaBins.push_back(+2.4);
  const int nEtaBins = etaBins.size()-1;
  std::cout << "Number of eta bins: " << nEtaBins << std::endl;

  std::vector<double> centralityBins;
  centralityBins.push_back(0.00);
  if (doCentrality){
    //centralityBins.push_back(0.05);
    //centralityBins.push_back(0.10);
    //centralityBins.push_back(0.15);
    centralityBins.push_back(0.20);
    centralityBins.push_back(0.40);
  }
  centralityBins.push_back(0.80);
  const int nCentralityBins = centralityBins.size()-1;
  std::cout << "Number of centrality bins: " << nCentralityBins << std::endl;

  // ---sub-division of dataset --- //
  RooDataSet* dataSubSet[nPtBins][nEtaBins][nCentralityBins];
  RooPlot* frameEff[nPtBins][nEtaBins][nCentralityBins];
  RooPlot* frameEff_Plus[nPtBins][nEtaBins][nCentralityBins];
  RooPlot* frameEff_Minus[nPtBins][nEtaBins][nCentralityBins];
  RooPlot* frameEffEta[nPtBins][nEtaBins][nCentralityBins];
  TCanvas* cTrig = new TCanvas("cTrig","cTrig",600,600);
  TCanvas* cTrig_plus = new TCanvas("cTrig_plus","cTrig_plus",600,600);
  TCanvas* cTrig_minus = new TCanvas("cTrig_minus","cTrig_minus",600,600);
  TCanvas* cTrigEta = new TCanvas("cTrigEta","cTrigEta",600,600);

  TLatex l;
  l.SetNDC();

   for (int i=0; i<nPtBins; i++){
     for (int j=0; j<nEtaBins; j++){
       for (int k=0; k<nCentralityBins; k++){
         if(!doCharge){
         std::cout << "Subset " << i << ":"<< j << ":" << k << std::endl;
	     dataSubSet[i][j][k] = selectPtEtaCentrality( data, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true);
	     dataSubSet[i][j][k]->Print();

         //effPdf.fitTo(*dataSubSet[i][j][k],ConditionalObservables(muonPt),Range(ptcutLow,ptmax)) ;

  	     frameEff[i][j][k] = muonPt.frame();
      	 dataSubSet[i][j][k]->plotOn(frameEff[i][j][k],Efficiency(cut),Binning(b0)) ;

	     effFunc.plotOn(frameEff[i][j][k],LineColor(kRed)) ;
	     frameEff[i][j][k]->SetMaximum(1.2); frameEff[i][j][k]->SetMinimum(-0.1) ;
	     frameEff[i][j][k]->SetXTitle("#font[52]{p}_{T}^{#mu} [GeV]"); frameEff[i][j][k]->SetYTitle("#epsilon_{trigger}");
	     frameEff[i][j][k]->Draw();

         TString sCentLow = ""; 
         TString sCentHigh = ""; 
         TString sEtaLow = ""; 
         TString sEtaHigh = "";

	     TString sCentLimits = sCentLow; sCentLimits+="-"; sCentLimits+=sCentHigh; sCentLimits+="%";
	     TString sEtaLimits = sEtaLow; sEtaLimits+="<|#eta|<"; sEtaLimits+=sEtaHigh;
         l.DrawLatex(0.3,0.45, sCentLimits );
         l.DrawLatex(0.3,0.35, sEtaLimits );

         //TString plotNameLog = "effPt_"; plotNameLog+=sCentLimits; plotNameLog+=","; plotNameLog+=sEtaLimits; 
         std::stringstream plotNameLog; 
         plotNameLog << "effPt_" << sCentLimits << "," << sEtaLimits;
         std::string sPlotName = plotNameLog.str();
         std::cout << "Efficiency plateau for bin " <<  i << ":"<< j << ":" << k << a1.getVal() << std::endl;
	     cTrig->Print((sPlotName+".png").c_str());
	     //cTrig->Print(plotNameLog+".pdf");
	     //cTrig->Print(plotNameLog+".eps");
	     //cTrig->Print(sPlotName.c_str()plotNameLog+".root");
	     //cTrig->Print(plotNameLog+".C");
        }
        else{
         std::cout << "Subset " << i << ":"<< j << ":" << k << std::endl;
	     dataSubSet[i][j][k] = selectPtEtaCentrality( data, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true);
	     dataSubSet[i][j][k]->Print();

         RooDataSet* dataMinus = (RooDataSet*) dataSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
         dataMinus->Print();
         RooDataSet* dataMinusPassed = (RooDataSet*)dataMinus->reduce(Cut("cut==cut::accept"));
         dataMinusPassed->Print();
         RooDataSet* dataPlus = (RooDataSet*) dataSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
         dataPlus->Print();
         RooDataSet* dataPlusPassed = (RooDataSet*)dataPlus->reduce(Cut("cut==cut::accept"));
         dataPlusPassed->Print();

  	     frameEffEta[i][j][k] = muonEta.frame();
  	     frameEff_Plus[i][j][k] = muonPt.frame();
  	     frameEff_Minus[i][j][k] = muonPt.frame();

         effPdf.fitTo(*dataPlus,ConditionalObservables(muonPt),Range(0.0,50.0)) ;
      	 dataPlus->plotOn(frameEff_Plus[i][j][k],Efficiency(cut),Binning(b0)) ;
	     effFunc.plotOn(frameEff_Plus[i][j][k],LineColor(kRed)) ;
      	 dataPlus->plotOn(frameEffEta[i][j][k],Efficiency(cut),Binning(bEta)) ;

         effPdf.fitTo(*dataMinus,ConditionalObservables(muonPt),Range(0.0,50.0)) ;
      	 dataMinus->plotOn(frameEff_Minus[i][j][k],Efficiency(cut),Binning(b0)) ;
	     effFunc.plotOn(frameEff_Minus[i][j][k],LineColor(kRed)) ;
      	 dataMinus->plotOn(frameEffEta[i][j][k],Efficiency(cut),Binning(bEta)) ;

         char hNameMinus[50],hNamePlus[50] ;
         sprintf(hNameMinus,"hTrigEffMinusEtaCent%i",k);
         sprintf(hNamePlus,"hTrigEffPlusEtaCent%i",k);
         frameEffEta[i][j][k]->Print();
         TH1F* hdataEtaPlus = (TH1F*)dataPlus->createHistogram("hdataEtaPlus",muonEta,Binning(bEta));
         TH1F* hdataEtaPlusPassed = (TH1F*)dataPlusPassed->createHistogram("hdataEtaPlusPassed",muonEta,Binning(bEta));
         writeEfficiency(outFile,(*hdataEtaPlusPassed),(*hdataEtaPlus),hNamePlus);


         TH1F* hdataEtaMinus = (TH1F*)dataMinus->createHistogram("hdataEtaMinus",muonEta,Binning(bEta));
         TH1F* hdataEtaMinusPassed = (TH1F*)dataMinusPassed->createHistogram("hdataEtaMinusPassed",muonEta,Binning(bEta));
         writeEfficiency(outFile,(*hdataEtaMinusPassed),(*hdataEtaMinus),hNameMinus);

         // mu+
	     frameEff_Plus[i][j][k]->SetMaximum(1.2); frameEff_Plus[i][j][k]->SetMinimum(-0.1) ;
	     frameEffEta[i][j][k]->SetMaximum(1.2); frameEffEta[i][j][k]->SetMinimum(-0.1) ;
	     frameEff_Plus[i][j][k]->SetXTitle("#font[52]{p}_{T}^{#mu} [GeV]"); frameEff_Plus[i][j][k]->SetYTitle("#epsilon_{trigger}");
	     frameEffEta[i][j][k]->SetXTitle("#eta_{#mu^{+}} [GeV]"); frameEffEta[i][j][k]->SetYTitle("#epsilon_{trigger}");
         cTrig_plus->cd();
	     frameEff_Plus[i][j][k]->Draw();

         TString sCentLow = ""; 
         TString sCentHigh = ""; 
         TString sEtaLow = ""; 
         TString sEtaHigh = "";
 		 sCentLow += format(100*centralityBins[k]); //sCentLow.Remove(3);
         sCentHigh += format(100*centralityBins[k+1]); //sCentUp.Remove(3);

 		 sCentLow += format(100*centralityBins[k]); //sCentLow.Remove(3);
         sCentHigh += format(100*centralityBins[k+1]); //sCentUp.Remove(3);

		 sEtaLow += format(etaBins[j]);
		 sEtaHigh += format(etaBins[j+1]);

	     TString sCentLimits = sCentLow; sCentLimits+="-"; sCentLimits+=sCentHigh; sCentLimits+="%";
	     TString sEtaLimits = sEtaLow; sEtaLimits+="<|#eta|<"; sEtaLimits+=sEtaHigh;
         l.DrawLatex(0.3,0.45, sCentLimits );
         l.DrawLatex(0.3,0.35, sEtaLimits );

         //cTrigEta->cd();
	     //frameEffEta[i][j][k]->Draw();

         //l.DrawLatex(0.3,0.45, sCentLimits );

/*         TString plotNameLog = "effPt_"; plotNameLog+=sCentLimits; plotNameLog+=","; plotNameLog+=sEtaLimits; 
	     cTrig->Print(plotNameLog+".png");
	     cTrig->Print(plotNameLog+".pdf");
	     cTrig->Print(plotNameLog+".eps");
	     cTrig->Print(plotNameLog+".root");

         plotNameLog = "effEta_"; plotNameLog+=sCentLimits; plotNameLog+=","; plotNameLog+=sEtaLimits; 
         cTrigEta->Print(plotNameLog+".root");
*/
         std::stringstream plotNameLog; 
         plotNameLog << "effPt_" << sCentLimits << "_muplus_" << sEtaLimits;
         std::string sPlotName = plotNameLog.str();
         std::cout << "Efficiency plateau for bin " <<  i << ":"<< j << ":" << k << a1.getVal() << std::endl;
	     cTrig_plus->Print((sPlotName+".png").c_str());
	     cTrig_plus->Print((sPlotName+".pdf").c_str());
	     cTrig_plus->Print((sPlotName+".C").c_str());
	     cTrig_plus->Print("effPt_muplus.root");
	     //cTrig_plus->Print((sPlotName+".root").c_str());

         // mu-
	     frameEff_Minus[i][j][k]->SetMaximum(1.2); frameEff_Minus[i][j][k]->SetMinimum(-0.1) ;
	     frameEff_Minus[i][j][k]->SetXTitle("#font[52]{p}_{T}^{#mu} [GeV]"); frameEff_Minus[i][j][k]->SetYTitle("#epsilon_{trigger}");
         cTrig_minus->cd();
	     frameEff_Minus[i][j][k]->Draw();

         sCentLow = ""; 
         sCentHigh = ""; 
         sEtaLow = ""; 
         sEtaHigh = "";
 		 sCentLow += format(100*centralityBins[k]); //sCentLow.Remove(3);
         sCentHigh += format(100*centralityBins[k+1]); //sCentUp.Remove(3);

 		 sCentLow += format(100*centralityBins[k]); //sCentLow.Remove(3);
         sCentHigh += format(100*centralityBins[k+1]); //sCentUp.Remove(3);

		 sEtaLow += format(etaBins[j]);
		 sEtaHigh += format(etaBins[j+1]);

	     sCentLimits = sCentLow; sCentLimits+="-"; sCentLimits+=sCentHigh; sCentLimits+="%";
	     sEtaLimits = sEtaLow; sEtaLimits+="<|#eta|<"; sEtaLimits+=sEtaHigh;
         l.DrawLatex(0.3,0.45, sCentLimits );
         l.DrawLatex(0.3,0.35, sEtaLimits );

         plotNameLog << "effPt_" << sCentLimits << "_muminus_" << sEtaLimits;
         sPlotName = plotNameLog.str();
         std::cout << "Efficiency plateau for bin " <<  i << ":"<< j << ":" << k << a1.getVal() << std::endl;
	     cTrig_minus->Print((sPlotName+".png").c_str());
	     cTrig_minus->Print((sPlotName+".pdf").c_str());
	     cTrig_minus->Print((sPlotName+".C").c_str());
	     //cTrig_minus->Print((sPlotName+".root").c_str());
	     cTrig_minus->Print("effPt_muminus.root");

        }

       }
     }
   }


    //outFile->Write();
    outFile->Close();
/*  
  //Eta and centrality binning
  std::cout << "Fitting efficiency over all centrality and eta..." << std::endl;
  ///use for unbinned ML fit
  effPdf.fitTo(*data,ConditionalObservables(muonPt),Range(ptcutLow,ptmax), Minos(kTRUE));
  RooPlot* frame1 = muonPt.frame(Title("Data (all, accepted)"));
  //plot all data
  data->plotOn(frame1,Binning(b1), DataError(RooAbsData::SumW2));
  //plot accepted data
  data->plotOn(frame1,Cut("cut==cut::accept"),Binning(b1),MarkerColor(kRed),LineColor(kRed), DataError(RooAbsData::SumW2)) ;

  RooPlot* frame2 = muonPt.frame(Title("Fitted efficiency")) ;

  frame1->SetXTitle("#font[52]{p}_{T}^{#mu} [GeV]") ; frame2->SetXTitle("#font[52]{p}_{T}^{#mu} [GeV]");

  frame1->SetYTitle("Muons/[GeV]"); frame2->SetYTitle("#epsilon_{trigger}");
  
  std::cout << "plotting... " << std::endl;
  //data->plotOn(frame2,Efficiency(cut),Binning(b0)) ;
  data->plotOn(frame2,Efficiency(cut),Binning(b0)) ;
  frame2->SetMarkerSize(0.6);
  //data->plotOn(frame2,Efficiency(cut),Binning(b0)) ;
  //data->plotOn(frame2,Efficiency(cut),Binning(b0),DataError(RooAbsData::SumW2)) ;
  ///uncomment to plot eff curve
  effFunc.plotOn(frame2,LineColor(kRed)) ;
  std::cout << "done."<<std::endl;
  //std::cout << "chi^2 " << frame2->chiSquare() << endl;
  frame2->SetMaximum(1.2); frame2->SetMinimum(-0.1) ;

*/
 }
