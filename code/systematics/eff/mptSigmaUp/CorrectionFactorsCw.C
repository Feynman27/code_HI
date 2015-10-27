#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooDataSet.h"
#include "RooBinning.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooFormulaVar.h"

#include "TDirectory.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath> 

#include "AtlasUtils.C"
#include "AtlasStyle.C"

#include "TriggerEfficiencies.C"

using namespace RooFit;


///////////////////////////////////////////////////////////////////////////////
//plot difference in charge Cw 
///////////////////////////////////////////////////////////////////////////////
plotChargeDiffCw(TGraph* grDiff, TGraphErrors* grPlus, TGraphErrors* grMinus, int ieta, int icent, double npart){

	double* yPlus = grPlus->GetY();
	double* yMinus = grMinus->GetY();

	double yTempPlus = yPlus[icent]; 
	double yTempMinus = yMinus[icent]; 

	double cWdiff = fabs(yTempPlus-yTempMinus);
	std::cout << "Difference in Cw at Eta Bin " << ieta << " and Centrality Bin " << icent << " = " << cWdiff << std::endl;
	grDiff->SetPoint(icent,npart,cWdiff);

}
///////////////////////////////////////////////////////////////////////////////
//plot difference in charge A0,A1,and A2
///////////////////////////////////////////////////////////////////////////////
plotChargeDiffAx(TGraph* grDiff, TGraphErrors* grPlus, TGraphErrors* grMinus, int ieta, double xEta){

	double* yPlus = grPlus->GetY();
	double* yMinus = grMinus->GetY();

	double yTempPlus = yPlus[ieta]; 
	double yTempMinus = yMinus[ieta]; 

	double Axdiff = fabs(yTempPlus-yTempMinus);
	std::cout << "Difference in Cw at Eta Bin " << ieta << " = " << Axdiff << std::endl;
	grDiff->SetPoint(ieta,xEta,Axdiff);

}
///////////////////////////////////////////////////////////////////////////////
//index ipt,ieta,icent 
///////////////////////////////////////////////////////////////////////////////
int indexIJK(int ich, int i, int j, int k, int nCentrality, int nEta){

	int increment = 0;	
	if(ich==102) increment = 0;
	else if(ich==103) increment = 1;
	else if(ich==104) increment = 2;
	else {
		std::cout << "WARNING: Unable to find charge identifier. Will return -9999." << std::endl;
		return -9999;
	}
	int index =  2*(nCentrality*nEta*i + nCentrality*j + k) + increment ;
	return index;
}
///////////////////////////////////////////////////////////////////////////////
//getCw 
///////////////////////////////////////////////////////////////////////////////
void getCw(RooDataSet* mcWGenSet, RooDataSet* mcWFidSet, RooDataSet* mcWRecSet, int ich, int ipt, int ieta, int icent,int nCentralityBins, int nEtaBins, 
		double npart, TGraphErrors* gr, int index, double xEta, TGraph2DErrors* gr2D){

	//number of Wmu in this eta and centrality class
	double ptGen = (double)mcWGenSet->numEntries();

	//number of muons at generator level in kinematic fiducial region within this eta and centrality class
	double ptGenInFiducial = (double)mcWFidSet->numEntries();

	//number of reconstructed signal muons after final selection within the fiducial region in this eta and centrality class
	double ptRecTemp = (double)mcWRecSet->numEntries() ;

	//fraction of signal muons in the fiducial region that were reconstructed
	double recoTemp = ptRecTemp/ptGenInFiducial;
	double recoStatTempErr = sqrt(TMath::Power(sqrt(ptRecTemp)/ptRecTemp*100,2) + TMath::Power(sqrt(ptGenInFiducial)/ptGenInFiducial*100,2))*0.01*recoTemp;

	//trigger efficiency estimate in the eta and centrality class (from Z note)
	double trigTemp = trigEfficiency(0,ieta,icent);
	double trigStatTempErr = trigEfficiencyErr(0,ieta,icent);

	double cWTemp = trigTemp*recoTemp;
	double CwStatTempErr = sqrt(TMath::Power(recoStatTempErr/recoTemp*100,2) + TMath::Power(trigStatTempErr*100,2))*0.01*cWTemp;

	//Cw centrality dependence 
	gr->SetPoint(icent,npart,cWTemp);
	gr->SetPointError(icent,0.0,CwStatTempErr);

	gr2D->SetPoint(index, npart, xEta, cWTemp); 
	gr2D->SetPointError(index, 0.0, 0.0, CwStatTempErr); 
	int index2 = indexIJK(ich,ipt,ieta,icent,nCentralityBins,nEtaBins);
/*	std::cout << " index "  << " : " <<  " charge "<< ":" <<" eta "<<":" <<" centrality " << ":" << " Generated Evts "<<"\t"<< " Fidicial " << "\t" << 
		"Trigger Eff " << "\t"<< " Reconstructed "<< "\t"<< " Cw "<< " " << "+-" << std::endl ;

	std::cout<< index2 << ":" << ich<< ":" <<ipt << ":" <<ieta<< ":" <<icent << "\t" << ptGen << "\t" << ptGenInFiducial << "\t" << 
		trigTemp << "\t" << ptRecTemp << "\t" << cWTemp << " " << CwStatTempErr << std::endl;
*/
	std::cout << cWTemp << " " << CwStatTempErr << std::endl;

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

  dataSet = (RooDataSet*) dataSet->reduce(cut);
  return dataSet;
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

  dataSet = (RooDataSet*) dataSet->reduce(cut);
  return dataSet;
}

///////////////////////////////////////////////////////////////////////////////
//selectPtEtaCentrality 
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectPtEtaCentrality( RooDataSet* dataSet,double ptLow, double ptUpp, double etaLow, double etaUpp, double centralityLow, double centralityUpp, 
	bool doMirrorEta=false, bool isGen=false){

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
  float ptcone20Nt[50];
  int valNt[50], ZDYNt[50], truthMatchedNt[50], promptNt[50];
  int nmu;

  // --- Load the MC tree ---
  TChain* tree = new TChain("tree","tree");
  int nFiles = tree->Add(pathName+fileName);
  
  std::cout << "Filling the RooDataSet for " << fileName << "... Number of files: " << nFiles << std::endl;
 
  // --- Set branch adresses ---
  tree->SetBranchAddress("eLoss", &eLossNt);
  tree->SetBranchAddress("ptcone20", &ptcone20Nt);
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
  tree->SetBranchStatus("ptcone20", 1);
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

      double isolation20Temp = -9999.0;
      isolation20Temp = ptcone20Nt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolation20",isolation20Temp);
      muonArgSet.setRealValue("muonPt",ptNt[imu]);
      muonArgSet.setRealValue("muonMt",mtNt[imu]);
      muonArgSet.setRealValue("muonEta",etaNt[imu]);
      muonArgSet.setRealValue("ZDY",ZDYNt[imu]);
      muonArgSet.setRealValue("muonCharge",chargeNt[imu]);
      muonArgSet.setRealValue("muonCategory",promptNt[imu]);
//      muonArgSet.setRealValue("muonGenRecMatched",truthMatchedNt[imu]);
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

void Write(TFile* outFile, TObject* gr, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing TGraph to root file..." << std::endl;
	  gr->Write(sFile);
	  std::cout << "Done." << std::endl;
    	  gDirectory = dir;
}



void fitParam(TGraphErrors* gr, int ich, TString sPar){

	  std::cout << "Fitting to analytic function for " << sPar << ":" << ich << std::endl;
	  TF1 f0a = TF1("f0a","[0]+[1]*x+[2]*x*x",0.0,2.7);	
	  //TF1 f0a = TF1("f0a","TMath::Sqrt([0]*x)",0.0,500.0);	
	  gr->Fit("f0a");
	  f0a.Draw("same");
	  cout << "chi2 = " << f0a.GetChisquare() << "/" << f0a.GetNDF() << ", p = " << f0a.GetProb() << endl;


}

void SaveGraph(TGraph* gr, TString sY, TString sBin, TString sLabel, bool doSystematic=false){

	TCanvas* c = new TCanvas("c","c",600,600);
	gr->Draw("ape"); 
	gr->GetYaxis()->SetTitle(sY);
	if(doSystematic) gr->GetYaxis()->SetRangeUser(0.0,0.1);
	else gr->GetYaxis()->SetRangeUser(0.0,1.0);
	gr->GetXaxis()->SetTitle("#LT N_{part} #GT");
	myText(0.33,0.89, (Color_t)kBlack, (char*)(sBin));
	c->Update();
	c->Print(sLabel+".pdf");	

}

void fit(TGraphErrors* gr, TGraphErrors* grA0, TGraphErrors* grA1, TGraphErrors* grA2, int ich, int ieta, double xEta, bool fixA2, bool fixA1){
	  std::cout << "Fitting to analytic function for " << ich << ":" << ieta << std::endl;
	  //a2 coeffs equal for mu+,mu- within 1sigma
	  //double a2 = (2.84778e-7+4.32554e-7)/2.0;
	  //double a2 = (4.56612e-07+7.19152e-07)/2.0;
	  //mu+-
	  double a2 = -1.01081e-07;
	  TF1* f0a; /*= new TF1("f0a","[0]+[1]*x+[2]*x*x",0.0,500.0);*/	
	  //f0a->FixParameter(2,a2);
 
 	//using average of both charges
/*	 double b0 = (-1.22747e-03+-1.15866e-03)/2.0;
	 double b1 = (-5.47588e-05+-2.42778e-05)/2.0;
	 double b2 = (1.24201e-04+1.09511e-04)/2.0;
	 double b3 = (1.15664e-04+8.42216e-05)/2.0;
	 double a1 = b0+b1*TMath::Cos(xEta)+b2*TMath::Cos(2.0*xEta)+b3*TMath::Cos(3.0*xEta);
*/

	 /*double b0 = -1.17977e-03;
	 double b1 = 7.08965e-05;
	 double b2 = 6.04041e-05;
	 double a1 = b0+b1*TMath::Cos(3.0*xEta)+b2*TMath::Cos(4.0*xEta);
	 */
	 double b0 = -5.38156e-04;
	 double b1 = -2.42858e-05;
	 double a1 = b0 + b1*xEta;
	 //temp hack 
	 //double a1 = -0.00116983;
	 //mu+ index
	 if(ich==102) { 
		//double four = -6.934e-4-6.0196e-5*TMath::Cos(xEta)+1.406e-5*TMath::Cos(2.*xEta)-9.778e-6*TMath::Cos(3.*xEta)+5.0098e-5*TMath::Cos(xEta*4.);
		//double a1 = -7.668e-4+1.219e-5*xEta*xEta;
		f0a = new TF1("f0a","[0]+[1]*x+[2]*x*x",0.0,500.0);	
		
	 	if(fixA1) f0a->FixParameter(1,a1);
	 	if(fixA2) f0a->FixParameter(2,a2);
		
 	 }
	  //mu- index
	  else if(ich==103) { 
		//double four = -7.7963e-4-1.0117e-4*TMath::Cos(xEta)-4.2672e-5*TMath::Cos(2.*xEta);
		//double a1 = -8.437e-4+3.879e-5*xEta*xEta;
		f0a = new TF1("f0a","[0]+[1]*x+[2]*x*x",0.0,500.0);	
		
	 	if(fixA1) f0a->FixParameter(1,a1);
	 	if(fixA2) f0a->FixParameter(2,a2);
		
		
	  }

	  //mu+- index
	  else {
	  	f0a = new TF1("f0a","[0]+[1]*x+[2]*x*x",0.0,500.0);	

		
	 	if(fixA1) f0a->FixParameter(1,a1);
	 	if(fixA2) f0a->FixParameter(2,a2);
		
	  }
	  //TF1 f0a = TF1("f0a","TMath::Sqrt([0]*x)",0.0,500.0);	
	  gr->Fit("f0a");
	  f0a->Draw("same");
	  cout << "chi2 = " << f0a->GetChisquare() << "/" << f0a->GetNDF() << ", p = " << f0a->GetProb() << endl;

	  //fill graph pts with parameter value to study eta dependence
	  double a0 = f0a->GetParameter(0);
	  double a0err = f0a->GetParError(0);
	  if(!fixA1) a1 = f0a->GetParameter(1);
	  double a1err = f0a->GetParError(1);
	  if(!fixA2) a2 = f0a->GetParameter(2);
	  double a2err = f0a->GetParError(2);

	  grA0->SetPoint(ieta,xEta,a0);
	  grA0->SetPointError(ieta,0.0,a0err);
	  grA1->SetPoint(ieta,xEta,a1);
	  grA1->SetPointError(ieta,0.0,a1err);
	  grA2->SetPoint(ieta,xEta,a2);
	  grA2->SetPointError(ieta,0.0,a2err);
}

void plotFactor(TGraph2DErrors* gr, TString sCh){

	std::cout << "Plotting factors..." << std::endl;
	TLatex l;
	TCanvas* c = new TCanvas("c","c",600,600);
	l.SetNDC();
	l.DrawLatex(0.39,0.83,sCh);
	gr->Draw("ape");	
	gr->GetXaxis()->SetTitle("#GT N_{part} #LT");	
	gr->GetYaxis()->SetTitle("|#eta|");
	gr->GetZaxis()->SetTitle(sCh);
	gr->GetXaxis()->SetRangeUser(0.0,100.0);	
	gr->GetYaxis()->SetRangeUser(0.0,2.7);
	gr->GetZaxis()->SetRangeUser(0.0,1.0);

	c->Update();
	c->Print(sCh+".pdf"); c->Print(sCh+".root");
	std::cout << "Done." << std::endl;

}
	
void CorrectionFactorsCw(){

bool doCharge = true;
bool doCentrality = true;
bool doEta = true;
bool doPt = false;
bool doFit = true;

/// --- output file ---
TDirectory *dir = gDirectory;
TFile *outFile = new TFile("correctionFactorsW.root","RECREATE");
gDirectory = dir;

//Correction factors for fitting and mT methods
const int chargeBins = 2;
const int globalBin  = 1;

//read in trigger efficiencies
readInputFile();

TString fileNameIn = "HISingleMuonMCWmunu.12.25.2012";

std::vector <double> etaBins;
etaBins.push_back(0.0);
if(doEta){
	etaBins.push_back(0.25);
	etaBins.push_back(0.5);
	etaBins.push_back(0.75);
	etaBins.push_back(1.0);
	etaBins.push_back(1.25);
	etaBins.push_back(1.5);
	etaBins.push_back(1.75);
	etaBins.push_back(2.0);
	etaBins.push_back(2.25);
	}
etaBins.push_back(2.5);

const int nEtaBins = etaBins.size()-1;

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

std::cout << "Initialiazing TGraph with " << nCentralityBins*nEtaBins << " points. " << std::endl;
TGraph2DErrors* gr2DPlus = new TGraph2DErrors(nCentralityBins*nEtaBins);
TGraph2DErrors* gr2DMinus = new TGraph2DErrors(nCentralityBins*nEtaBins);
TGraph2DErrors* gr2D = new TGraph2DErrors(nCentralityBins*nEtaBins);

TGraphErrors* grCentPlus = new TGraphErrors(nCentralityBins);
TGraphErrors* grCentMinus = new TGraphErrors(nCentralityBins);
TGraph* grCentDiff = new TGraphErrors(nCentralityBins);
TGraphErrors* grCent = new TGraphErrors(nCentralityBins);

//graphs for fit parameter eta dependence
TGraphErrors* grA0Plus = new TGraphErrors(nEtaBins);
TGraphErrors* grA1Plus = new TGraphErrors(nEtaBins);
TGraphErrors* grA2Plus = new TGraphErrors(nEtaBins);

TGraphErrors* grA0Minus = new TGraphErrors(nEtaBins);
TGraphErrors* grA1Minus = new TGraphErrors(nEtaBins);
TGraphErrors* grA2Minus = new TGraphErrors(nEtaBins);

//difference in plus/minus for systematics
TGraph* grA0Diff = new TGraphErrors(nEtaBins);
TGraph* grA1Diff = new TGraphErrors(nEtaBins);
TGraph* grA2Diff = new TGraphErrors(nEtaBins);

TGraphErrors* grA0 = new TGraphErrors(nEtaBins);
TGraphErrors* grA1 = new TGraphErrors(nEtaBins);
TGraphErrors* grA2 = new TGraphErrors(nEtaBins);

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
  RooRealVar  isolation20("isolation20","isolation20",0.0,10.0);
  RooRealVar  centrality("centrality","centrality",0.,1.0);
  RooRealVar  muonEta("muonEta","muonEta",-3.0,+3.0);
  RooRealVar  ZDY("ZDY","ZDY",0.0,1.0);
  RooRealVar  muonQuality("muonQuality","muonQuality",0.0,15.0);
//  RooRealVar  muonGenRecMatched("muonGenRecMatched","muonGenRecMatched",0.0,2.0);
  RooRealVar  muonELoss("muonELoss","muonELoss",-10.0,+10.0);
  RooRealVar  muonScat("muonScat","muonScat",-10.0,+10.0);

  TString sCutsGen = "abs(mother)==24&&abs(daughter)==13";
  RooArgList muonGenArgList(mother,daughter);
  RooFormulaVar cutsGen("cutsGen", "cutsGen", sCutsGen, muonGenArgList);
  //fiducial cuts
  //TString sCutsFid = "muonGenPt>25.0&&abs(etaGen)<2.5&&nuGenPt>25.0&&munuGenMt>40.0&&abs(mother)==24&&abs(daughter)==13";

  //systematics for +1sigma pt and mpt
  std::cout << "Doing +1 sigma systematics for MPT..." << std::endl;
  TString sCutsFid = "muonGenPt>25.0&&abs(etaGen)<2.5&&nuGenPt>36.05&&munuGenMt>40.0&&abs(mother)==24&&abs(daughter)==13";
  RooFormulaVar cutsFid("cutsFid", "cutsFid", sCutsFid, RooArgList(muonGenPt,etaGen,nuGenPt,munuGenMt,mother,daughter));

  //reconstruction level cuts
  //systematics
  TString sCutsRec = "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.5&&muonPt>25.0&&missPt>36.05&&isolation20<0.3&&muonMt>40.0&&ZDY==0";
  RooArgList muonRecArgList(muonQuality,muonELoss,muonScat,muonEta,muonPt,missPt,isolation20,muonMt);
  muonRecArgList.add(ZDY);
  RooFormulaVar cutsRec("cutsRec", "cutsRec", sCutsRec, muonRecArgList);

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
  muonRecArgSet.add(isolation20);
  muonRecArgSet.add(muonQuality);
//  muonRecArgSet.add(muonGenRecMatched);
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

  RooDataSet* mcWGenSet = fillHIMuonGenSet("/tmp/tbalestr/",fileNameIn+".root",muonGenArgSet); mcWGenSet->Print();
  mcWGenSet = (RooDataSet*)mcWGenSet->reduce(Cut(cutsGen)); 
  std::cout << "Number of Wmunu evts at generator level : " << mcWGenSet->numEntries() << std::endl;

  RooDataSet* mcWGenSetPlus = (RooDataSet*)mcWGenSet->reduce("chargeGenCategory==chargeGenCategory::muPlus"); 
  RooDataSet* mcWGenSetMinus = (RooDataSet*)mcWGenSet->reduce("chargeGenCategory==chargeGenCategory::muMinus"); 

  RooDataSet* mcWFidSet = fillHIMuonGenSet("/tmp/tbalestr/",fileNameIn+".root",muonGenArgSet); 
  mcWFidSet = (RooDataSet*)mcWFidSet->reduce(Cut(cutsFid)); 
  std::cout << "Number of Wmunu evts in fiducial region at generator level : " << mcWFidSet->numEntries() << std::endl;

  //construct charged sets
  RooDataSet* mcWFidSetPlus = (RooDataSet*)mcWFidSet->reduce("chargeGenCategory==chargeGenCategory::muPlus"); 
  RooDataSet* mcWFidSetMinus = (RooDataSet*)mcWFidSet->reduce("chargeGenCategory==chargeGenCategory::muMinus"); 

  double fidPlus = mcWFidSetPlus->numEntries(); double genPlus = mcWGenSetPlus->numEntries();
  double fidMinus = mcWFidSetMinus->numEntries(); double genMinus = mcWGenSetMinus->numEntries();

  std::cout << " Aw+ = " << mcWFidSetPlus->numEntries() << "/" << mcWGenSetPlus->numEntries() << " = " << fidPlus/genPlus << std::endl;
  std::cout << " Aw- = " << mcWFidSetMinus->numEntries() << "/" << mcWGenSetMinus->numEntries() << " = " << fidMinus/genMinus << std::endl;

  RooDataSet* mcWRecSet = fillHIMuonRecSet("/tmp/tbalestr/",fileNameIn+".root",muonRecArgSet); 
  mcWRecSet = (RooDataSet*)mcWRecSet->reduce(Cut(cutsRec)); 
  std::cout << "Number of Wmunu evts reconstructed in fiducial region : " << mcWRecSet->numEntries() << std::endl;
  // --- Subdivide in bins ---
  RooDataSet* mcWGenSubSet[nPtBins][nEtaBins][nCentralityBins];
  RooDataSet* mcWFidSubSet[nPtBins][nEtaBins][nCentralityBins];
  RooDataSet* mcWRecSubSet[nPtBins][nEtaBins][nCentralityBins];
  for ( int i = 0; i < nPtBins; i++ ) {
    for ( int j = 0; j < nEtaBins; j++ ) {
      for ( int k = 0; k < nCentralityBins; k++ ){
	mcWGenSubSet[i][j][k] = selectPtEtaCentrality(mcWGenSet ,ptBins[i],ptBins[i+1], etaBins[j], etaBins[j+1], centBins[k], centBins[k+1],true,true); 
	mcWFidSubSet[i][j][k] = selectPtEtaCentrality(mcWFidSet ,ptBins[i],ptBins[i+1], etaBins[j], etaBins[j+1], centBins[k], centBins[k+1],true,true); 
	mcWRecSubSet[i][j][k] = selectPtEtaCentrality(mcWRecSet ,ptBins[i],ptBins[i+1], etaBins[j], etaBins[j+1], centBins[k], centBins[k+1],true,false); 
    }
  }
}

	  for(int ipt=0; ipt<nPtBins; ipt++){

		for(int ieta=0; ieta<nEtaBins; ieta++){

		  double etabinLo = etaBins.at(ieta);
		  double etabinUp = etaBins.at(ieta+1); 
		  double xEta = etabinLo+(etabinUp-etabinLo)/2.0; 

		  for(int icent=0; icent<nCentralityBins; icent++){

			double centbinLo = centBins.at(icent);
			double centbinUp = centBins.at(icent+1); 

			int index = ieta*nCentralityBins+icent; 
			double npart = npartBins.at(icent); 
			
			//if(!doCharge){
				std::cout << "mu^{#pm}" << std::endl;
				getCw(mcWGenSubSet[ipt][ieta][icent],mcWFidSubSet[ipt][ieta][icent],mcWRecSubSet[ipt][ieta][icent],104,
					ipt,ieta, icent, nCentralityBins,nEtaBins,npart, grCent,index,xEta,gr2D);
			//}
			if(doCharge) {
			
		   		RooDataSet* mcWGenSetPlus  = (RooDataSet*) mcWGenSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muPlus"));
		   		RooDataSet* mcWGenSetMinus  = (RooDataSet*) mcWGenSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muMinus"));

		   		RooDataSet* mcWFidSetPlus  = (RooDataSet*) mcWFidSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muPlus"));
		   		RooDataSet* mcWFidSetMinus  = (RooDataSet*) mcWFidSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muMinus"));

		   		RooDataSet* mcWRecSetPlus  = (RooDataSet*) mcWRecSubSet[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   		RooDataSet* mcWRecSetMinus  = (RooDataSet*) mcWRecSubSet[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

				std::cout << "mu^{+}" << std::endl;
				getCw(mcWGenSetPlus,mcWFidSetPlus,mcWRecSetPlus,102,ipt,ieta, icent, nCentralityBins,nEtaBins,npart, grCentPlus,index,xEta,gr2DPlus);

				std::cout << "mu^{-}" << std::endl;
				getCw(mcWGenSetMinus,mcWFidSetMinus,mcWRecSetMinus,103,ipt,ieta, icent, nCentralityBins,nEtaBins,npart, grCentMinus,index,xEta,gr2DMinus);

				std::cout << "Plotting difference of Cw+ and Cw- for systematics." << std::endl;
				plotChargeDiffCw(grCentDiff,grCentPlus,grCentMinus,ieta, icent,npart);

			}

		  } //icent

		TString sNamePlus = "grCent";sNamePlus+=102;sNamePlus+="eta";sNamePlus+=ieta;
		TString sNameMinus = "grCent";sNameMinus+=103;sNameMinus+="eta";sNameMinus+=ieta;
		TString sNameDiff = "grCentDiff"; sNameDiff+="eta";sNameDiff+=ieta;
		TString sName = "grCent";sName+=104;sName+="eta";sName+=ieta;

		TString sEtaTemp = ""; sEtaTemp+=etabinLo; sEtaTemp+="<";sEtaTemp+="|#eta|"; sEtaTemp+="<"; sEtaTemp+=etabinUp; 


		if(doFit){
			fit(grCentPlus,grA0Plus,grA1Plus,grA2Plus,102,ieta,xEta,true,true);
			fit(grCentMinus,grA0Minus,grA1Minus,grA2Minus,103,ieta,xEta,true,true);

			std::cout << "Plotting charge difference of a0,a1,and a2 for systematics." << std::endl;
			plotChargeDiffAx(grA0Diff,grA0Plus,grA0Minus,ieta,xEta);
			plotChargeDiffAx(grA1Diff,grA1Plus,grA1Minus,ieta,xEta);
			plotChargeDiffAx(grA2Diff,grA2Plus,grA2Minus,ieta,xEta);

			fit(grCent,grA0,grA1,grA2,104,ieta,xEta,true,true);

			SaveGraph(grCentPlus,"C_{W^{+}}",sEtaTemp,sNamePlus);
			SaveGraph(grCentMinus,"C_{W^{-}}",sEtaTemp,sNameMinus);
			SaveGraph(grCentDiff,"|C_{W^{+}}-C_{W^{-}}|",sEtaTemp,sNameDiff,true);
			SaveGraph(grCent,"C_{W}",sEtaTemp,sName);

			Write(outFile,grCentPlus,sNamePlus);
			Write(outFile,grCentMinus,sNameMinus);
			Write(outFile,grCentDiff,sNameDiff);
			Write(outFile,grCent,sName);
		}
	   } //ieta
	} //ipt

	if(doFit){
		fitParam(grA0Plus,102,"a0");
		fitParam(grA1Plus,102,"a1");
		fitParam(grA2Plus,102,"a2");

		fitParam(grA0Minus,103,"a0");
		fitParam(grA1Minus,103,"a1");
		fitParam(grA2Minus,103,"a2");

		fitParam(grA0,104,"a0");
		fitParam(grA1,104,"a1");
		fitParam(grA2,104,"a2");
		
		TString sNamePlus = "gr";sNamePlus+=102;
		TString sNamePlusA0 = sNamePlus + "a0";
		TString sNamePlusA1 = sNamePlus + "a1";
		TString sNamePlusA2 = sNamePlus + "a2";

		TString sNameMinus = "gr";sNameMinus+=103;
		TString sNameMinusA0 = sNameMinus + "a0";
		TString sNameMinusA1 = sNameMinus + "a1";
		TString sNameMinusA2 = sNameMinus + "a2";

		TString sNameDiffA0 = "grA0Diff";
		TString sNameDiffA1 = "grA1Diff";
		TString sNameDiffA2 = "grA2Diff";

		TString sName = "gr";sName+=104;
		TString sNameA0 = sName + "a0";
		TString sNameA1 = sName + "a1";
		TString sNameA2 = sName + "a2";

		Write(outFile,grA0Plus,sNamePlusA0);
		Write(outFile,grA1Plus,sNamePlusA1);
		Write(outFile,grA2Plus,sNamePlusA2);

		Write(outFile,grA0Minus,sNameMinusA0);
		Write(outFile,grA1Minus,sNameMinusA1);
		Write(outFile,grA2Minus,sNameMinusA2);

		Write(outFile,grA0Diff,sNameDiffA0);
		Write(outFile,grA1Diff,sNameDiffA1);
		Write(outFile,grA2Diff,sNameDiffA2);

		Write(outFile,grA0,sNameA0);
		Write(outFile,grA1,sNameA1);
		Write(outFile,grA2,sNameA2);

		//plot 2D distribution of Cw for mu+-
		plotFactor(gr2DPlus,"C_{W^{+}}");
		plotFactor(gr2DMinus,"C_{W^{-}}");
		plotFactor(gr2D,"C_{W}");

		sNamePlus = "gr2D";sNamePlus+=102;
		sNameMinus = "gr2D";sNameMinus+=103;
		sName = "gr2D";sName+=104;

		Write(outFile,gr2DPlus,sNamePlus);
		Write(outFile,gr2DMinus,sNameMinus);
		Write(outFile,gr2D,sName);
	}
	
} //eof

int main(){
	CorrectionFactorsCw();
}
