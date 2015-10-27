#include "TDirectory.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath> 

#include "TriggerEfficiencies.C"

	
void CorrectionFactors(){

//Correction factors for fitting and mT methods
bool doMtCorrFact = true;
bool doFitCorrFact = false;

const int chargeBins = 2;
const int globalBin  = 1;
int totalBins;
int totalChargeBins;
int totalAwBins;
int totalAwChargeBins;

float arrAw[100] = {0};
float arrAwSysErrSq[100] = {0};
float arrAwStatErrSq[100] = {0};
float AwSysErr = 0;
float AwStatErr = 0;
float AwSysSumSq = 0;
float AwStatSumSq = 0;
float AwSum  = 0;
float AwpSum  = 0;
float AwmSum  = 0;
float AwMean = 0;
float AwpMean = 0;
float AwmMean = 0;

float arrCw[100] = {0};
float arrCwSysErrSq[100] = {0};
float arrCwStatErrSq[100] = {0};
float CwSysErr = 0;
float CwStatErr = 0;
float CwSysSumSq = 0;
float CwStatSumSq = 0;
float CwSum  = 0;
float CwpSum  = 0;
float CwmSum  = 0;
float CwMean = 0;
float CwpMean = 0;
float CwmMean = 0;


float arrTrig[100] = {0};
float arrTrigSysErrSq[100] = {0};
float arrTrigStatErrSq[100] = {0};
float TrigSysErr = 0;
float TrigStatErr = 0;
float TrigSysSumSq = 0;
float TrigStatSumSq = 0;
float TrigSum  = 0;
float TrigMean = 0;

float arrReco[100] = {0};
float arrRecoSysErrSq[100] = {0};
float arrRecoStatErrSq[100] = {0};
float RecoSysErr = 0;
float RecoStatErr = 0;
float RecoSysSumSq = 0;
float RecoStatSumSq = 0;
float RecoSum  = 0;
float RecoMean = 0;



//read in trigger efficiencies
readInputFile();

//TFile *_file0 = TFile::Open("/tmp/tbalestr/HISingleMuon_mcWmunuDataOverlay_aug28_v2.root");
//TString fileNameIn = "/tmp/tbalestr/HISingleMuon_mcWmunuDataOverlay_aug28_v2";
//TString fileNameIn = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/HISingleMuonWmunuMCDataOverlay.2012.11.21";
TString fileNameIn = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/HISingleMuonWmunuMCDataOverlay.2012.11.13";
TFile* fIn = new TFile(fileNameIn+".root", "READ");

if ( !fIn->IsOpen() ) {
    std::cout << fIn << " not found!" << std::endl;
    exit(0);
  }

TTree *tree = (TTree*)fIn->Get("tree");

//generator level
TH1F* hPtGenP = new TH1F("hPtGenP","hPtGenP",100,0.,100.);
TH1F* hPtGenM = new TH1F("hPtGenM","hPtGenM",100,0.,100.);
TH1F* hPtGenPM = new TH1F("hPtGenPM","hPtGenPM",100,0.,100.);

//generator level muon tracks in fiducial region
TH1F* hPtPlusGenCut = new TH1F("hPtPlusGenCut","hPtPlusGenCut",100,0.0,100.);
TH1F* hPtMinusGenCut = new TH1F("hPtMinusGenCut","hPtMinusGenCut",100,0.0,100.);
TH1F* hPtPlusMinusGenCut = new TH1F("hPtPlusMinusGenCut","hPtPlusMinusGenCut",100,0.0,100.);

//reconstructed muon tracks in fiducial region
TH1F* hPtPlusRec = new TH1F("hPtPlusRec","hPtPlusRec",100,0.,100.);
TH1F* hPtMinusRec = new TH1F("hPtMinusRec","hPtMinusRec",100,0.,100.);
TH1F* hPtPlusMinusRec = new TH1F("hPtPlusMinusRec","hPtPlusMinusRec",100,0.,100.);


//total Wmu candidates at generator level
tree->Draw("mc_mu_gen_pt>>hPtGenP","mc_mu_charge==+1");
tree->Draw("mc_mu_gen_pt>>hPtGenM","mc_mu_charge==-1");
tree->Draw("mc_mu_gen_pt>>hPtGenPM","abs(mc_mu_charge)==1");


//generator level muon tracks in fiducial region
TString cutsGlobal = "abs(mc_mu_gen_eta)<2.5";
TString cutsPlusGlobal = cutsGlobal + "&&mc_mu_charge==1"; 
TString cutsMinusGlobal = cutsGlobal + "&&mc_mu_charge==-1"; 


//Wmunu candidates at generator level after fiducial cuts
tree->Draw("mc_mu_gen_pt>>hPtPlusGenCut",cutsPlusGlobal);
tree->Draw("mc_mu_gen_pt>>hPtMinusGenCut",cutsMinusGlobal);
tree->Draw("mc_mu_gen_pt>>hPtPlusMinusGenCut",cutsGlobal);

//reconstructed Wmunu candidates
TString cutsGlobalRec; TString cutsPlusGlobalRec; TString cutsMinusGlobalRec;
if(doMtCorrFact){
	std::cout << "Calculating Cw,Aw for mT method." << std::endl;
	cutsGlobalRec = cutsGlobal + "&&centrality<=0.8&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&nu_pt>25.&&pt>25.&&mt>40.0&&mt<200.&&ptcone20/pt<0.3";
	cutsPlusGlobalRec = cutsPlusGlobal + "&&charge==+1&&centrality<=0.8&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&pt>25.&&nu_pt>25.&&mt>40.0&&mt<200.&&ptcone20/pt<0.3";
	cutsMinusGlobalRec = cutsMinusGlobal + "&&charge==-1&&centrality<=0.8&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&pt>25.&&nu_pt>25.&&mt>40.0&&mt<200.&&ptcone20/pt<0.3";
}
else if(doFitCorrFact){
	std::cout << "Calculating Cw,Aw for fitting method." << std::endl;
	cutsGlobalRec = cutsGlobal + "&&centrality<=0.8&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&pt>7.0";
	cutsPlusGlobalRec = cutsPlusGlobal + "&&charge==+1&&centrality<=0.8&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&pt>7.0";
	cutsMinusGlobalRec = cutsMinusGlobal + "&&charge==-1&&centrality<=0.8&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&pt>7.0";
}

tree->Draw("mc_mu_gen_pt>>hPtPlusRec",cutsPlusGlobalRec);
tree->Draw("mc_mu_gen_pt>>hPtMinusRec",cutsMinusGlobalRec);
tree->Draw("mc_mu_gen_pt>>hPtPlusMinusRec",cutsGlobalRec);

//calculate Aw correction factor
float ptPlusGenCut = hPtPlusGenCut->Integral();
float ptPlusGen = hPtGenP->Integral();
float aWPlus = ptPlusGenCut/ptPlusGen;
float AwPlusStatErr = sqrt(sqrt(ptPlusGen)/ptPlusGen*100 + sqrt(ptPlusGenCut)/ptPlusGenCut*100)*0.01*aWPlus;

AwSum+=aWPlus;
AwpSum+=aWPlus;
float ptMinusGenCut = hPtMinusGenCut->Integral();
float ptMinusGen = hPtGenM->Integral();
float aWMinus = ptMinusGenCut/ptMinusGen;
float AwMinusStatErr = sqrt(sqrt(ptMinusGen)/ptMinusGen*100 + sqrt(ptMinusGenCut)/ptMinusGenCut*100)*0.01*aWMinus;

AwSum+=aWMinus;
AwmSum+=aWMinus;

float ptPlusMinusGenCut = hPtPlusMinusGenCut->Integral();
float ptPlusMinusGen = hPtGenPM->Integral();
float aWPlusMinus = ptPlusMinusGenCut/ptPlusMinusGen;
float AwPlusMinusStatErr = sqrt(sqrt(ptPlusMinusGen)/ptPlusMinusGen*100 + sqrt(ptPlusMinusGenCut)/ptPlusMinusGenCut*100)*0.01*aWPlusMinus;


AwSum+=aWPlusMinus;


//calculate Cw from reco efficiency
//and trigger efficiency;
//Note:charge trig eff is integrated
//over all eta and centrality only

float ptPlusRec = hPtPlusRec->Integral();
//float ptPlusGenCut = hPtPlusGenCut->Integral();
float recoPlus = ptPlusRec/ptPlusGenCut;
float RecoPlusStatErr = sqrt(sqrt(ptPlusRec)/ptPlusRec*100 + sqrt(ptPlusGenCut)/ptPlusGenCut*100)*0.01*recoPlus;
RecoSum+=recoPlus;
float trigPlus = trigEfficiency(1,0,0);
//float trigPlus = 1.00;
float TrigPlusStatErr = trigEfficiencyErr(1,0,0);
TrigSum+=trigPlus;


float ptMinusRec = hPtMinusRec->Integral();
//float ptMinusGenCut = hPtMinusGenCut->Integral();
float recoMinus = ptMinusRec/ptMinusGenCut;
float RecoMinusStatErr = sqrt(sqrt(ptMinusRec)/ptMinusRec*100 + sqrt(ptMinusGenCut)/ptMinusGenCut*100)*0.01*recoMinus;
RecoSum+=recoMinus;
float trigMinus =trigEfficiency(-1,0,0);
//float trigMinus = 1.00;
float TrigMinusStatErr =trigEfficiencyErr(-1,0,0);
TrigSum+=trigMinus;

float ptPlusMinusRec = hPtPlusMinusRec->Integral();
//float ptPlusMinusGenCut = hPtPlusMinusGenCut->Integral();
float recoPlusMinus = ptPlusMinusRec/ptPlusMinusGenCut;
float RecoPlusMinusStatErr = sqrt(sqrt(ptPlusMinusRec)/ptPlusMinusRec*100 + sqrt(ptPlusMinusGenCut)/ptPlusMinusGenCut*100)*0.01*recoPlusMinus;
RecoSum+=recoPlusMinus;
float trigPlusMinus = trigEfficiency(0,0,0);
float TrigPlusMinusStatErr = trigEfficiencyErr(0,0,0);
TrigSum+=trigPlusMinus;

float cWPlus = trigPlus*recoPlus;
float cWMinus = trigMinus*recoMinus;
float cWPlusMinus = trigPlusMinus*recoPlusMinus;
float CwPlusStatErr = sqrt(RecoPlusStatErr/recoPlus*100 + TrigPlusStatErr*100)*0.01*cWPlus;
float CwMinusStatErr = sqrt(RecoMinusStatErr/recoMinus*100 + TrigMinusStatErr*100)*0.01*cWMinus;
float CwPlusMinusStatErr = sqrt(RecoPlusMinusStatErr/recoPlusMinus*100 + TrigPlusMinusStatErr*100)*0.01*cWPlusMinus;

CwSum+=cWPlus;
CwpSum+=cWPlus;
CwSum+=cWMinus;
CwmSum+=cWMinus;
CwSum+=cWPlusMinus;

//print results
//std::cout<<setw(82)<<"**********************************************************************************************************\n";
std::cout<<setw(82)<<"centrality"<<setw(10)<<"eta"<<setw(10)<<"charge"<<setw(14)<<"WSelection"<<setw(10)<<"trigger"<<setw(10)<<"Cw"<<setw(10)<<"Aw"<<setw(10)<<"\n";
std::cout<<setw(82)<<"0-0.8"<<setw(10)<<"0-2.5"<<setw(10)<<"+"<<setw(10)<<recoPlus<<setw(10)<<trigPlus<<setw(10)<<cWPlus<<setw(10)<<aWPlus<<"\n";
std::cout<<setw(82)<<"0-0.8"<<setw(10)<<"0-2.5"<<setw(10)<<"-"<<setw(10)<<recoMinus<<setw(10)<<trigMinus<<setw(10)<<cWMinus<<setw(10)<<aWMinus<<"\n";
std::cout<<setw(82)<<"0-0.8"<<setw(10)<<"0-2.5"<<setw(10)<<"+/-"<<setw(10)<<recoPlusMinus<<setw(10)<<trigPlusMinus<<setw(10)<<cWPlusMinus<<setw(10)<<aWPlusMinus<<"\n";

std::cout<<setw(82)<<"centrality"<<setw(10)<<"eta"<<setw(10)<<"charge"<<setw(10)<<"GenEvts"<<setw(10)<<"EvtsInAcc"<<setw(14)<<"EvtsAfterWSel"<<setw(10)<<"\n";
std::cout<<setw(82)<<"0-0.8"<<setw(10)<<"0-2.5"<<setw(10)<<"+"<<setw(10)<<ptPlusGen<<setw(10)<<ptPlusGenCut<<setw(10)<<ptPlusRec<<"\n";
std::cout<<setw(82)<<"0-0.8"<<setw(10)<<"0-2.5"<<setw(10)<<"-"<<setw(10)<<ptMinusGen<<setw(10)<<ptMinusGenCut<<setw(10)<<ptMinusRec<<"\n";
std::cout<<setw(82)<<"0-0.8"<<setw(10)<<"0-2.5"<<setw(10)<<"-"<<setw(10)<<ptPlusMinusGen<<setw(10)<<ptPlusMinusGenCut<<setw(10)<<ptPlusMinusRec<<"\n";

int icharge = 0;

//positive 
arrAw[icharge] = aWPlus;
arrAwStatErrSq[icharge] = pow(AwPlusStatErr,2);
arrTrig[icharge] = trigPlus;
arrTrigStatErrSq[icharge] = pow(TrigPlusStatErr,2);
arrReco[icharge] = recoPlus;
arrRecoStatErrSq[icharge] = pow(RecoPlusStatErr,2);
arrCw[icharge] = cWPlus;
arrCwStatErrSq[icharge] = pow(CwPlusStatErr,2);

//negative
icharge++;
arrAw[icharge] = aWMinus;
arrAwStatErrSq[icharge] = pow(AwMinusStatErr,2);
arrTrig[icharge] = trigMinus;
arrTrigStatErrSq[icharge] = pow(TrigMinusStatErr,2);
arrReco[icharge] = recoMinus;
arrRecoStatErrSq[icharge] = pow(RecoMinusStatErr,2);
arrCw[icharge] = cWMinus;
arrCwStatErrSq[icharge] = pow(CwMinusStatErr,2);

//both charges
icharge++;
arrAw[icharge] = aWPlusMinus;
arrAwStatErrSq[icharge] = pow(AwPlusMinusStatErr,2);
arrTrig[icharge] = trigPlusMinus;
arrTrigStatErrSq[icharge] = pow(TrigPlusMinusStatErr,2);
arrReco[icharge] = recoPlusMinus;
arrRecoStatErrSq[icharge] = pow(RecoPlusMinusStatErr,2);
arrCw[icharge] = cWPlusMinus;
arrCwStatErrSq[icharge] = pow(CwPlusMinusStatErr,2);

icharge = 0;

//centrality
std::vector <float> centBins;
centBins.push_back(0.0);
centBins.push_back(0.05);
centBins.push_back(0.1);
centBins.push_back(0.15);
centBins.push_back(0.2);
centBins.push_back(0.4);
centBins.push_back(0.8);

const int centralityBins = centBins.size()-1;

TH1F* hPtMinusGenCent[centralityBins]  ;
TH1F* hPtPlusGenCent[centralityBins] ;
TH1F* hPtPlusMinusGenCent[centralityBins] ;

TH1F* hPtMinusGenCutCent[centralityBins] ;
TH1F* hPtPlusGenCutCent[centralityBins] ;
TH1F* hPtPlusMinusGenCutCent[centralityBins] ;

TH1F* hPtMinusRecCent[centralityBins] ;
TH1F* hPtPlusRecCent[centralityBins] ;
TH1F* hPtPlusMinusRecCent[centralityBins] ;

	char hPlusGen[200];
	char hPlusGenCut[200];
	char hPlusRec[200];
	char hMinusGen[200];
	char hMinusGenCut[200];
	char hMinusRec[200];
	char hPlusMinusGen[200];
	char hPlusMinusGenCut[200];
	char hPlusMinusRec[200];

//fill centrality binned histos
for(int icent=0; icent<centralityBins; icent++){

	icharge=0;
	
/*	char hPlusGen[200];
	char hPlusGenCut[200];
	char hPlusRec[200];
	char hMinusGen[200];
	char hMinusGenCut[200];
	char hMinusRec[200];
	char hPlusMinusGen[200];
	char hPlusMinusGenCut[200];
	char hPlusMinusRec[200];
*/

	sprintf(hPlusGen,"hPtPlusGenCent%i",icent);
	sprintf(hMinusGen,"hPtMinusGenCent%i",icent);
	sprintf(hPlusMinusGen,"hPtPlusMinusGenCent%i",icent);

	sprintf(hPlusGenCut,"hPtPlusGenCutCent%i",icent);
	sprintf(hMinusGenCut,"hPtMinusGenCutCent%i",icent);
	sprintf(hPlusMinusGenCut,"hPtPlusMinusGenCutCent%i",icent);

	sprintf(hPlusRec,"hPtPlusRecCent%i",icent);
	sprintf(hMinusRec,"hPtMinusRecCent%i",icent);
	sprintf(hPlusMinusRec,"hPtPlusMinusRecCent%i",icent);

	hPtPlusGenCent[icent] = new TH1F(hPlusGen,hPlusGen,100,0.,100.);
	hPtMinusGenCent[icent] = new TH1F(hMinusGen,hMinusGen,100,0.,100.);
	hPtPlusMinusGenCent[icent] = new TH1F(hPlusMinusGen,hPlusMinusGen,100,0.,100.);

	hPtPlusGenCutCent[icent] = new TH1F(hPlusGenCut,hPlusGenCut,100,0.,100.);
	hPtMinusGenCutCent[icent] = new TH1F(hMinusGenCut,hMinusGenCut,100,0.,100.);
	hPtPlusMinusGenCutCent[icent] = new TH1F(hPlusMinusGenCut,hPlusMinusGenCut,100,0.,100.);

	hPtPlusRecCent[icent] = new TH1F(hPlusRec,hPlusRec,100,0.,100.);
	hPtMinusRecCent[icent] = new TH1F(hMinusRec,hMinusRec,100,0.,100.);
	hPtPlusMinusRecCent[icent] = new TH1F(hPlusMinusRec,hPlusMinusRec,100,0.,100.);


	float binLo = centBins.at(icent);
	float binUp = centBins.at(icent+1); 
	//std::cout << binUp << std::endl;

	TString binCut = "&&centrality>=";
	binCut+=binLo;
	binCut+= "&&centrality<";
	binCut+= binUp;

	//generator level muons
	tree->Draw((TString("mc_mu_gen_pt>>")+hPlusGen).Data(),"mc_mu_charge==+1"+binCut);
	tree->Draw((TString("mc_mu_gen_pt>>")+hMinusGen).Data(),"mc_mu_charge==-1"+binCut);
	tree->Draw((TString("mc_mu_gen_pt>>")+hPlusMinusGen).Data(),"abs(mc_mu_charge)==1"+binCut);

	//generator level muon tracks in fiducial region
        TString cuts = "abs(mc_mu_gen_eta)<2.5";

	cuts+=binCut;
	TString cutsPlus = cuts + "&&mc_mu_charge==1"; 
	TString cutsMinus = cuts + "&&mc_mu_charge==-1"; 

	tree->Draw((TString("mc_mu_gen_pt>>")+hPlusGenCut).Data(),cutsPlus);
	//tree->Draw(hPlusGenCut,cutsPlus);

	tree->Draw((TString("mc_mu_gen_pt>>")+hMinusGenCut).Data(),cutsMinus);

	tree->Draw((TString("mc_mu_gen_pt>>")+hPlusMinusGenCut).Data(),cuts);

	//reconstruced and select W candidates in fiducial region
	if(doMtCorrFact){
		cuts += "&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&nu_pt>25.&&pt>25.&&mt>40.0&&mt<200.&&ptcone20/pt<0.3";
		cutsPlus += "&&charge==+1&&centrality<=0.8&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&pt>25.&&nu_pt>25.&&mt>40.0&&mt<200.&&ptcone20/pt<0.3";
		cutsMinus += "&&charge==-1&&centrality<=0.8&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&pt>25.&&nu_pt>25.&&mt>40.0&&mt<200.&&ptcone20/pt<0.3";
	}
	else if(doFitCorrFact){
		cuts += "&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&pt>7.0";
		cutsPlus += "&&charge==+1&&centrality<=0.8&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&pt>7.0";
		cutsMinus += "&&charge==-1&&centrality<=0.8&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&pt>7.0";
	}

	tree->Draw((TString("mc_mu_gen_pt>>")+hPlusRec).Data(),cutsPlus);

	tree->Draw((TString("mc_mu_gen_pt>>")+hMinusRec).Data(),cutsMinus);

	tree->Draw((TString("mc_mu_gen_pt>>")+hPlusMinusRec).Data(),cuts);


	//calculate Aw correction factor
	float ptPlusGenCutTemp = hPtPlusGenCutCent[icent]->Integral();
	//std::cout << "Integral\t" << ptPlusGenCutTemp <<std::endl; exit(0);
	float ptPlusGenTemp = hPtPlusGenCent[icent]->Integral();
	float aWPlusTemp = ptPlusGenCutTemp/ptPlusGenTemp;
	float AwPlusStatErrTemp = sqrt(sqrt(ptPlusGenTemp)/ptPlusGenTemp*100 + sqrt(ptPlusGenCutTemp)/ptPlusGenCutTemp*100)*0.01*aWPlusTemp;


	float ptMinusGenCutTemp = hPtMinusGenCutCent[icent]->Integral();
	float ptMinusGenTemp = hPtMinusGenCent[icent]->Integral();
	float aWMinusTemp = ptMinusGenCutTemp/ptMinusGenTemp;
	float AwMinusStatErrTemp = sqrt(sqrt(ptMinusGenTemp)/ptMinusGenTemp*100 + sqrt(ptMinusGenCutTemp)/ptMinusGenCutTemp*100)*0.01*aWMinusTemp;

	float ptPlusMinusGenCutTemp = hPtPlusMinusGenCutCent[icent]->Integral();
	float ptPlusMinusGenTemp = hPtPlusMinusGenCent[icent]->Integral();
	float aWPlusMinusTemp = ptPlusMinusGenCutTemp/ptPlusMinusGenTemp;
	float AwPlusMinusStatErrTemp = sqrt(sqrt(ptPlusMinusGenTemp)/ptPlusMinusGenTemp*100 + sqrt(ptPlusMinusGenCutTemp)/ptPlusMinusGenCutTemp*100)*0.01*aWPlusMinusTemp;


	//calculate Cw from reco efficiency
	//and trigger efficiency

	float ptPlusRecTemp = hPtPlusRecCent[icent]->Integral();
//	float ptPlusGenCutTemp = hPtPlusGenCutCent[icent]->Integral();
	float recoPlusTemp = ptPlusRecTemp/ptPlusGenCutTemp;
	float recoPlusStatTempErr = sqrt(sqrt(ptPlusRecTemp)/ptPlusRecTemp*100 + sqrt(ptPlusGenCutTemp)/ptPlusGenCutTemp*100)*0.01*recoPlusTemp;
	float trigPlusTemp = trigEfficiency(0,0,icent+1);
	float trigPlusStatTempErr = trigEfficiencyErr(0,0,icent+1);

	float ptMinusRecTemp = hPtMinusRecCent[icent]->Integral();
	//float ptMinusGenCutTemp = hPtMinusGenCent[icent]->Integral();
	float recoMinusTemp = ptMinusRecTemp/ptMinusGenCutTemp;
	float recoMinusStatTempErr = sqrt(sqrt(ptMinusRecTemp)/ptMinusRecTemp*100 + sqrt(ptMinusGenCutTemp)/ptMinusGenCutTemp*100)*0.01*recoMinusTemp;
	float trigMinusTemp = trigEfficiency(0,0,icent+1);
	float trigMinusStatTempErr = trigEfficiencyErr(0,0,icent+1);

	float ptPlusMinusRecTemp = hPtPlusMinusRecCent[icent]->Integral();
//	float ptPlusMinusGenCutTemp = hPtPlusMinusGenCent[icent]->Integral();
	float recoPlusMinusTemp = ptPlusMinusRecTemp/ptPlusMinusGenCutTemp;
	float recoPlusMinusStatTempErr = sqrt(sqrt(ptPlusMinusRecTemp)/ptPlusMinusRecTemp*100 + sqrt(ptPlusMinusGenCutTemp)/ptPlusMinusGenCutTemp*100)*0.01*recoPlusMinusTemp;
	float trigPlusMinusTemp = trigEfficiency(0,0,icent+1);
	float trigPlusMinusStatTempErr = trigEfficiencyErr(0,0,icent+1);

	float cWPlusTemp = trigPlusTemp*recoPlusTemp;
	float cWMinusTemp = trigMinusTemp*recoMinusTemp;
	float cWPlusMinusTemp = trigPlusMinusTemp*recoPlusMinusTemp;
	float CwPlusStatTempErr = sqrt(recoPlusStatTempErr/recoPlusTemp*100 + trigPlusStatTempErr*100)*0.01*cWPlusTemp;
	float CwMinusStatTempErr = sqrt(recoMinusStatTempErr/recoMinusTemp*100 + trigMinusStatTempErr*100)*0.01*cWMinusTemp;
	float CwPlusMinusStatTempErr = sqrt(recoPlusMinusStatTempErr/recoPlusMinusTemp*100 + trigPlusMinusStatTempErr*100)*0.01*cWPlusMinusTemp;

	AwSum+=aWPlusTemp;
	AwpSum+=aWPlusTemp;
	AwSum+=aWMinusTemp;
	AwmSum+=aWMinusTemp;
	AwSum+=aWPlusMinusTemp;

	TrigSum+=trigPlusTemp;
	TrigSum+=trigMinusTemp;
	TrigSum+=trigPlusMinusTemp;

	RecoSum+=recoPlusTemp;
	RecoSum+=recoMinusTemp;
	RecoSum+=recoPlusMinusTemp;


	CwSum+=cWPlusTemp;
	CwpSum+=cWPlusTemp;
	CwSum+=cWMinusTemp;
	CwmSum+=cWMinusTemp;
	CwSum+=cWPlusMinusTemp;



	std::cout<<setw(82)<<"centrality"<<setw(10)<<"eta"<<setw(10)<<"charge"<<setw(14)<<"WSelection"<<setw(10)<<"trigger"<<setw(10)<<"Cw"<<setw(10)<<"Aw"<<setw(10)<<"\n";
	std::cout<<setw(82)<<binLo<<"-"<<binUp<<setw(10)<<"0-2.5"<<setw(10)<<"+"<<setw(10)<<recoPlusTemp<<setw(10)<<trigPlusTemp<<setw(10)<<cWPlusTemp<<setw(10)<<aWPlusTemp<<"\n";
	std::cout<<setw(82)<<binLo<<"-"<<binUp<<setw(10)<<"0-2.5"<<setw(10)<<"-"<<setw(10)<<recoMinusTemp<<setw(10)<<trigMinusTemp<<setw(10)<<cWMinusTemp<<setw(10)<<aWMinusTemp<<"\n";
	std::cout<<setw(82)<<binLo<<"-"<<binUp<<setw(10)<<"0-2.5"<<setw(10)<<"+/-"<<setw(10)<<recoPlusMinusTemp<<setw(10)<<trigPlusMinusTemp<<setw(10)<<cWPlusMinusTemp<<setw(10)<<aWPlusMinusTemp<<"\n";

	std::cout<<setw(82)<<"centrality"<<setw(10)<<"eta"<<setw(10)<<"charge"<<setw(10)<<"GenEvts"<<setw(10)<<"EvtsInAcc"<<setw(14)<<"EvtsAfterWSel"<<setw(10)<<"\n";
	std::cout<<setw(82)<<binLo<<"-"<<binUp<<setw(10)<<"0-2.5"<<setw(10)<<"+"<<setw(10)<<ptPlusGenTemp<<setw(10)<<ptPlusGenCutTemp<<setw(10)<<ptPlusRecTemp<<"\n";
	std::cout<<setw(82)<<binLo<<"-"<<binUp<<setw(10)<<"0-2.5"<<setw(10)<<"-"<<setw(10)<<ptMinusGenTemp<<setw(10)<<ptMinusGenCutTemp<<setw(10)<<ptMinusRecTemp<<"\n";
	std::cout<<setw(82)<<binLo<<"-"<<binUp<<setw(10)<<"0-2.5"<<setw(10)<<"+/-"<<setw(10)<<ptPlusMinusGenTemp<<setw(10)<<ptPlusMinusGenCutTemp<<setw(10)<<ptPlusMinusRecTemp<<"\n";


	arrAw[(icent+1)*3+icharge] = aWPlusTemp;
	arrAwStatErrSq[(icent+1)*3+icharge] = pow(AwPlusStatErrTemp,2);
	arrTrig[(icent+1)*3+icharge] = trigPlusTemp;
	arrTrigStatErrSq[(icent+1)*3+icharge] = pow(trigPlusStatTempErr,2);
	arrReco[(icent+1)*3+icharge] = recoPlusTemp;
	arrRecoStatErrSq[(icent+1)*3+icharge] = pow(recoPlusStatTempErr,2);
	arrCw[(icent+1)*3+icharge] = cWPlusTemp;
	arrCwStatErrSq[(icent+1)*3+icharge] = pow(CwPlusStatTempErr,2);

	icharge++;
	arrAw[(icent+1)*3+icharge] = aWMinusTemp;
	arrAwStatErrSq[(icent+1)*3+icharge] = pow(AwMinusStatErrTemp,2);
	arrTrig[(icent+1)*3+icharge] = trigMinusTemp;
	arrTrigStatErrSq[(icent+1)*3+icharge] = pow(trigMinusStatTempErr,2);
	arrReco[(icent+1)*3+icharge] = recoMinusTemp;
	arrRecoStatErrSq[(icent+1)*3+icharge] = pow(recoMinusStatTempErr,2);
	arrCw[(icent+1)*3+icharge] = cWMinusTemp;
	arrCwStatErrSq[(icent+1)*3+icharge] = pow(CwMinusStatTempErr,2);

	icharge++;
	arrAw[(icent+1)*3+icharge] = aWPlusMinusTemp;
	arrAwStatErrSq[(icent+1)*3+icharge] = pow(AwPlusMinusStatErrTemp,2);
	arrTrig[(icent+1)*3+icharge] = trigPlusMinusTemp;
	arrTrigStatErrSq[(icent+1)*3+icharge] = pow(trigPlusMinusStatTempErr,2);
	arrReco[(icent+1)*3+icharge] = recoPlusMinusTemp;
	arrRecoStatErrSq[(icent+1)*3+icharge] = pow(recoPlusMinusStatTempErr,2);
	arrCw[(icent+1)*3+icharge] = cWPlusMinusTemp;
	arrCwStatErrSq[(icent+1)*3+icharge] = pow(CwPlusMinusStatTempErr,2);

	icharge = 0;
}

std::vector <float> etaBins;
etaBins.push_back(0.0);
etaBins.push_back(0.25);
etaBins.push_back(0.5);
etaBins.push_back(0.75);
etaBins.push_back(1.0);
etaBins.push_back(1.25);
etaBins.push_back(1.5);
etaBins.push_back(1.75);
etaBins.push_back(2.0);
etaBins.push_back(2.25);
etaBins.push_back(2.5);

const int etaBin = etaBins.size()-1;


TH1F* hPtMinusGenEta[etaBin]  ;
TH1F* hPtPlusGenEta[etaBin] ;
TH1F* hPtPlusMinusGenEta[etaBin] ;

TH1F* hPtMinusGenCutEta[etaBin] ;
TH1F* hPtPlusGenCutEta[etaBin] ;
TH1F* hPtPlusMinusGenCutEta[etaBin] ;

TH1F* hPtMinusRecEta[etaBin] ;
TH1F* hPtPlusRecEta[etaBin] ;
TH1F* hPtPlusMinusRecEta[etaBin] ;


for(int ieta=0; ieta<etaBin; ieta++){

	icharge = 0;

/*	char hPlusGen[200],hPlusGenCut[200],hPlusRec[200];
	char hMinusGen[200],hMinusGenCut[200],hMinusRec[200];
	char hPlusMinusGen[200],hPlusMinusGenCut[200],hPlusMinusRec[200];
*/
	sprintf(hPlusGen,"hPtPlusGenEta%i",ieta);
	sprintf(hMinusGen,"hPtMinusGenEta%i",ieta);
	sprintf(hPlusMinusGen,"hPtPlusMinusGenEta%i",ieta);

	sprintf(hPlusGenCut,"hPtPlusGenCutEta%i",ieta);
	sprintf(hMinusGenCut,"hPtMinusGenCutEta%i",ieta);
	sprintf(hPlusMinusGenCut,"hPtPlusMinusGenCutEta%i",ieta);

	sprintf(hPlusRec,"hPtPlusRecEta%i",ieta);
	sprintf(hMinusRec,"hPtMinusRecEta%i",ieta);
	sprintf(hPlusMinusRec,"hPtPlusMinusRecEta%i",ieta);

	hPtPlusGenEta[ieta] = new TH1F(hPlusGen,hPlusGen,100,0.,100.);
	hPtMinusGenEta[ieta] = new TH1F(hMinusGen,hMinusGen,100,0.,100.);
	hPtPlusMinusGenEta[ieta] = new TH1F(hPlusMinusGen,hPlusMinusGen,100,0.,100.);

	hPtPlusGenCutEta[ieta] = new TH1F(hPlusGenCut,hPlusGenCut,100,0.,100.);
	hPtMinusGenCutEta[ieta] = new TH1F(hMinusGenCut,hMinusGenCut,100,0.,100.);
	hPtPlusMinusGenCutEta[ieta] = new TH1F(hPlusMinusGenCut,hPlusMinusGenCut,100,0.,100.);

	hPtPlusRecEta[ieta] = new TH1F(hPlusRec,hPlusRec,100,0.,100.);
	hPtMinusRecEta[ieta] = new TH1F(hMinusRec,hMinusRec,100,0.,100.);
	hPtPlusMinusRecEta[ieta] = new TH1F(hPlusMinusRec,hPlusMinusRec,100,0.,100.);

	float binLo=etaBins.at(ieta);
	float binUp=etaBins.at(ieta+1);
	//(TString::Format("sBinUp==%f",binUp)).Data();

	TString binCut = "abs(mc_mu_gen_eta)>=";
	binCut+=binLo;
	binCut+= "&&abs(mc_mu_gen_eta)<";
	binCut+= binUp;

	//generator level muons
	tree->Draw((TString("mc_mu_gen_pt>>")+hPlusGen).Data(),binCut+"&&mc_mu_charge==+1");
	tree->Draw((TString("mc_mu_gen_pt>>")+hMinusGen).Data(),binCut+"&&mc_mu_charge==-1");
	tree->Draw((TString("mc_mu_gen_pt>>")+hPlusMinusGen).Data(),binCut+"&&abs(mc_mu_charge)==1");

	//generator level muon tracks in fiducial region
        TString cuts = binCut;
	//cuts+=binCut;
	TString cutsPlus = cuts + "&&mc_mu_charge==1"; 
	TString cutsMinus = cuts + "&&mc_mu_charge==-1"; 

	tree->Draw((TString("mc_mu_gen_pt>>")+hPlusGenCut).Data(),cutsPlus);
	//tree->Draw(hPlusGenCut,cutsPlus);

	tree->Draw((TString("mc_mu_gen_pt>>")+hMinusGenCut).Data(),cutsMinus);

	tree->Draw((TString("mc_mu_gen_pt>>")+hPlusMinusGenCut).Data(),cuts);

	//reconstruced and select W candidates in fiducial region
	if(doMtCorrFact){
		cuts += "&&centrality<=0.8&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&nu_pt>25.&&pt>25.&&mt>40.0&&mt<200.&&ptcone20/pt<0.3";
		cutsPlus += "&&charge==+1&&centrality<=0.8&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&nu_pt>25.&&pt>25.&&mt>40.0&&mt<200.&&ptcone20/pt<0.3";
		cutsMinus += "&&charge==-1&&centrality<=0.8&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&nu_pt>25.&&pt>25.&&mt>40.0&&mt<200.&&ptcone20/pt<0.3";
	}
	else if(doFitCorrFact){
		cuts += "&&centrality<=0.8&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&pt>7.0";
		cutsPlus += "&&charge==+1&&centrality<=0.8&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&pt>7.0";
		cutsMinus += "&&charge==-1&&centrality<=0.8&&truthMatched_muid==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&pt>7.0";
	}

	tree->Draw((TString("mc_mu_gen_pt>>")+hPlusRec).Data(),cutsPlus);

	tree->Draw((TString("mc_mu_gen_pt>>")+hMinusRec).Data(),cutsMinus);

	tree->Draw((TString("mc_mu_gen_pt>>")+hPlusMinusRec).Data(),cuts);

	//calculate Aw correction factor
	float ptPlusGenCutTemp = hPtPlusGenCutEta[ieta]->Integral();
	float ptPlusGenTemp = hPtPlusGenEta[ieta]->Integral();
	float aWPlusTemp = ptPlusGenCutTemp/ptPlusGenTemp;
	float AwPlusStatErrTemp = sqrt(sqrt(ptPlusGenTemp)/ptPlusGenTemp*100 + sqrt(ptPlusGenCutTemp)/ptPlusGenCutTemp*100)*0.01*aWPlusTemp;


	float ptMinusGenCutTemp = hPtMinusGenCutEta[ieta]->Integral();
	float ptMinusGenTemp = hPtMinusGenEta[ieta]->Integral();
	float aWMinusTemp = ptMinusGenCutTemp/ptMinusGenTemp;
	float AwMinusStatErrTemp = sqrt(sqrt(ptMinusGenTemp)/ptMinusGenTemp*100 + sqrt(ptMinusGenCutTemp)/ptMinusGenCutTemp*100)*0.01*aWMinusTemp;


	float ptPlusMinusGenCutTemp = hPtPlusMinusGenCutEta[ieta]->Integral();
	float ptPlusMinusGenTemp = hPtPlusMinusGenEta[ieta]->Integral();
	float aWPlusMinusTemp = ptPlusMinusGenCutTemp/ptPlusMinusGenTemp;
	float AwPlusMinusStatErrTemp = sqrt(sqrt(ptPlusMinusGenTemp)/ptPlusMinusGenTemp*100 + sqrt(ptPlusMinusGenCutTemp)/ptPlusMinusGenCutTemp*100)*0.01*aWPlusMinusTemp;


	//calculate Cw from reco efficiency
	//and trigger efficiency

	float ptPlusRecTemp = hPtPlusRecEta[ieta]->Integral();
//	float ptPlusGenCutTemp = hPtPlusGenCutEta[ieta]->Integral();
	float recoPlusTemp = ptPlusRecTemp/ptPlusGenCutTemp;
	float recoPlusStatTempErr = sqrt(sqrt(ptPlusRecTemp)/ptPlusRecTemp*100 + sqrt(ptPlusGenCutTemp)/ptPlusGenCutTemp*100)*0.01*recoPlusTemp;
	float trigPlusTemp = trigEfficiency(0,ieta+1,0);
	float trigPlusStatTempErr = trigEfficiencyErr(0,ieta+1,0);

	float ptMinusRecTemp = hPtMinusRecEta[ieta]->Integral();
	//float ptMinusGenCutTemp = hPtMinusGenCutEta[ieta]->Integral();
	float recoMinusTemp = ptMinusRecTemp/ptMinusGenCutTemp;
	float recoMinusStatTempErr = sqrt(sqrt(ptMinusRecTemp)/ptMinusRecTemp*100 + sqrt(ptMinusGenCutTemp)/ptMinusGenCutTemp*100)*0.01*recoMinusTemp;
	float trigMinusTemp = trigEfficiency(0,ieta+1,0);
	float trigMinusStatTempErr = trigEfficiencyErr(0,ieta+1,0);

	float ptPlusMinusRecTemp = hPtPlusMinusRecEta[ieta]->Integral();
//	float ptPlusMinusGenCutTemp = hPtPlusMinusGenCutEta[ieta]->Integral();
	float recoPlusMinusTemp = ptPlusMinusRecTemp/ptPlusMinusGenCutTemp;
	float recoPlusMinusStatTempErr = sqrt(sqrt(ptPlusMinusRecTemp)/ptPlusMinusRecTemp*100 + sqrt(ptPlusMinusGenCutTemp)/ptPlusMinusGenCutTemp*100)*0.01*recoPlusMinusTemp;
	float trigPlusMinusTemp = trigEfficiency(0,ieta+1,0);
	float trigPlusMinusStatTempErr = trigEfficiencyErr(0,ieta+1,0);

	float cWPlusTemp = trigPlusTemp*recoPlusTemp;
	float cWMinusTemp = trigMinusTemp*recoMinusTemp;
	float cWPlusMinusTemp = trigPlusMinusTemp*recoPlusMinusTemp;

	cWPlusTemp = trigPlusTemp*recoPlusTemp;
	cWMinusTemp = trigMinusTemp*recoMinusTemp;
	cWPlusMinusTemp = trigPlusMinusTemp*recoPlusMinusTemp;
	float CwPlusStatTempErr = sqrt(recoPlusStatTempErr/recoPlusTemp*100 + trigPlusStatTempErr*100)*0.01*cWPlusTemp;
	float CwMinusStatTempErr = sqrt(recoMinusStatTempErr/recoMinusTemp*100 + trigMinusStatTempErr*100)*0.01*cWMinusTemp;
	float CwPlusMinusStatTempErr = sqrt(recoPlusMinusStatTempErr/recoPlusMinusTemp*100 + trigPlusMinusStatTempErr*100)*0.01*cWPlusMinusTemp;

//don't include since theres no
//acceptance definition when eta binning
/*	AwSum+=aWPlusTemp;
	AwSum+=aWMinusTemp;
	AwSum+=aWPlusMinusTemp;
*/
	TrigSum+=trigPlusTemp;
	TrigSum+=trigMinusTemp;
	TrigSum+=trigPlusMinusTemp;

	RecoSum+=recoPlusTemp;
	RecoSum+=recoMinusTemp;
	RecoSum+=recoPlusMinusTemp;


	CwSum+=cWPlusTemp;
	CwpSum+=cWPlusTemp;
	CwSum+=cWMinusTemp;
	CwmSum+=cWMinusTemp;
	CwSum+=cWPlusMinusTemp;


	std::cout<<setw(82)<<"centrality"<<setw(10)<<"eta"<<setw(10)<<"charge"<<setw(14)<<"WSelection"<<setw(10)<<"trigger"<<setw(10)<<"Cw"<<setw(10)<<"Aw"<<setw(10)<<"\n";
	std::cout<<setw(82)<<"0-0.8"<<setw(10)<<binLo<<"-"<<binUp<<setw(10)<<"+"<<setw(10)<<recoPlusTemp<<setw(10)<<trigPlusTemp<<setw(10)<<cWPlusTemp<<setw(10)<<aWPlusTemp<<"\n";
	std::cout<<setw(82)<<"0-0.8"<<setw(10)<<binLo<<"-"<<binUp<<setw(10)<<"-"<<setw(10)<<recoMinusTemp<<setw(10)<<trigMinusTemp<<setw(10)<<cWMinusTemp<<setw(10)<<aWMinusTemp<<"\n";
	std::cout<<setw(82)<<"0-0.8"<<setw(10)<<binLo<<"-"<<binUp<<setw(10)<<"+/-"<<setw(10)<<recoPlusMinusTemp<<setw(10)<<trigPlusMinusTemp<<setw(10)<<cWPlusMinusTemp<<setw(10)<<aWPlusMinusTemp<<"\n";

	std::cout<<setw(82)<<"centrality"<<setw(10)<<"eta"<<setw(10)<<"charge"<<setw(10)<<"GenEvts"<<setw(10)<<"EvtsInAcc"<<setw(14)<<"EvtsAfterWSel"<<setw(10)<<"\n";
	std::cout<<setw(82)<<"0-0.8"<<setw(10)<<binLo<<"-"<<binUp<<setw(10)<<"+"<<setw(10)<<ptPlusGenTemp<<setw(10)<<ptPlusGenCutTemp<<setw(10)<<ptPlusRecTemp<<"\n";
	std::cout<<setw(82)<<"0-0.8"<<setw(10)<<binLo<<"-"<<binUp<<setw(10)<<"-"<<setw(10)<<ptMinusGenTemp<<setw(10)<<ptMinusGenCutTemp<<setw(10)<<ptMinusRecTemp<<"\n";
	std::cout<<setw(82)<<"0-0.8"<<setw(10)<<binLo<<"-"<<binUp<<setw(10)<<"+/-"<<setw(10)<<ptPlusMinusGenTemp<<setw(10)<<ptPlusMinusGenCutTemp<<setw(10)<<ptPlusMinusRecTemp<<"\n";

	arrAw[(ieta+1)*3+3*centralityBins+icharge] = aWPlusTemp;
	arrAwStatErrSq[(ieta+1)*3+3*centralityBins+icharge] = pow(AwPlusStatErrTemp,2);
	arrTrig[(ieta+1)*3+3*centralityBins+icharge] = trigPlusTemp;
	arrTrigStatErrSq[(ieta+1)*3+3*centralityBins+icharge] = pow(trigPlusStatTempErr,2);
	arrReco[(ieta+1)*3+3*centralityBins+icharge] = recoPlusTemp;
	arrRecoStatErrSq[(ieta+1)*3+3*centralityBins+icharge] = pow(recoPlusStatTempErr,2);
	arrCw[(ieta+1)*3+3*centralityBins+icharge] = cWPlusTemp;
	arrCwStatErrSq[(ieta+1)*3+3*centralityBins+icharge] = pow(CwPlusStatTempErr,2);

	icharge++;
	arrAw[(ieta+1)*3+3*centralityBins+icharge] = aWMinusTemp;
	arrAwStatErrSq[(ieta+1)*3+3*centralityBins+icharge] = pow(AwMinusStatErrTemp,2);
	arrTrig[(ieta+1)*3+3*centralityBins+icharge] = trigMinusTemp;
	arrTrigStatErrSq[(ieta+1)*3+3*centralityBins+icharge] = pow(trigMinusStatTempErr,2);
	arrReco[(ieta+1)*3+3*centralityBins+icharge] = recoMinusTemp;
	arrRecoStatErrSq[(ieta+1)*3+3*centralityBins+icharge] = pow(recoMinusStatTempErr,2);
	arrCw[(ieta+1)*3+3*centralityBins+icharge] = cWMinusTemp;
	arrCwStatErrSq[(ieta+1)*3+3*centralityBins+icharge] = pow(CwMinusStatTempErr,2);

	icharge++;
	arrAw[(ieta+1)*3+3*centralityBins+icharge] = aWPlusMinusTemp;
	arrAwStatErrSq[(ieta+1)*3+3*centralityBins+icharge] = pow(AwPlusMinusStatErrTemp,2);
	arrTrig[(ieta+1)*3+3*centralityBins+icharge] = trigPlusMinusTemp;
	arrTrigStatErrSq[(ieta+1)*3+3*centralityBins+icharge] = pow(trigPlusMinusStatTempErr,2);
	arrReco[(ieta+1)*3+3*centralityBins+icharge] = recoPlusMinusTemp;
	arrRecoStatErrSq[(ieta+1)*3+3*centralityBins+icharge] = pow(recoPlusMinusStatTempErr,2);
	arrCw[(ieta+1)*3+3*centralityBins+icharge] = cWPlusMinusTemp;
	arrCwStatErrSq[(ieta+1)*3+3*centralityBins+icharge] = pow(CwPlusMinusStatTempErr,2);

	icharge = 0;

}

	
	totalBins = globalBin + chargeBins + (globalBin + chargeBins)*centralityBins + (globalBin + chargeBins)*etaBin;
	totalChargeBins = 0.5*(chargeBins + (chargeBins)*centralityBins + (chargeBins)*etaBin);
	totalAwBins = globalBin + chargeBins + (globalBin + chargeBins)*centralityBins ;
	totalAwChargeBins = 0.5*(chargeBins + (chargeBins)*centralityBins );

	AwMean = AwSum/totalAwBins;
	AwpMean = AwpSum/totalAwChargeBins;
	AwmMean = AwmSum/totalAwChargeBins;
	TrigMean = TrigSum/totalBins;
	RecoMean = RecoSum/totalBins;
	CwMean = CwSum/totalBins;
	CwpMean = CwpSum/totalChargeBins;
	CwmMean = CwmSum/totalChargeBins;
	//std::cout << "total bins: " << totalBins << " AwSum: " << AwSum << " Aw Mean: " << AwMean << std::endl;

	for(int i = 0; i < totalBins; i++){
		
		//std::cout << arrAw[i] << std::endl;
		float AwTemp = arrAw[i];
		float diffTempAw = AwTemp-AwMean;
		arrAwSysErrSq[i] = pow(diffTempAw,2);
		AwSysSumSq += arrAwSysErrSq[i];
		AwStatSumSq += arrAwStatErrSq[i];

		float TrigTemp = arrTrig[i];
		float diffTempTrig = TrigTemp-TrigMean;
		arrTrigSysErrSq[i] = pow(diffTempTrig,2);
		TrigSysSumSq += arrTrigSysErrSq[i];
		TrigStatSumSq += arrTrigStatErrSq[i];

		float RecoTemp = arrReco[i];
		float diffTempReco = RecoTemp-RecoMean;
		arrRecoSysErrSq[i] = pow(diffTempReco,2);
		RecoSysSumSq += arrRecoSysErrSq[i];
		RecoStatSumSq += arrRecoStatErrSq[i];

		float CwTemp = arrCw[i];
		float diffTempCw = CwTemp-CwMean;
		arrCwSysErrSq[i] = pow(diffTempCw,2);
		CwSysSumSq += arrCwSysErrSq[i];
		CwStatSumSq += arrCwStatErrSq[i];


	}

	AwSysErr = sqrt(1./(totalBins-1.)*(AwSysSumSq));
	AwStatErr = sqrt(AwStatSumSq);

	TrigSysErr = sqrt(1./(totalBins-1.)*(TrigSysSumSq));
	TrigStatErr = sqrt(TrigStatSumSq);

	RecoSysErr = sqrt(1./(totalBins-1.)*(RecoSysSumSq));
	RecoStatErr = sqrt(RecoStatSumSq);

	CwSysErr = sqrt(1./(totalBins-1.)*(CwSysSumSq));
	CwStatErr = sqrt(CwStatSumSq);

	//std::cout << "Mean Muon Trigger Efficiency: " << TrigMean << " +- " << TrigStatErr << "(stat.) " << TrigSysErr << "(syst.)" << std::endl;
	std::cout << "Mean W Candidate Reconstruction and Selection Efficiency: " << RecoMean << " +- " << RecoStatErr << "(stat.) " << RecoSysErr << "(syst.)" << std::endl;
	std::cout << "Mean Aw+: " << AwpMean << std::endl;
	std::cout << "Mean Cw+: " << CwpMean << std::endl;
	std::cout << "Mean Aw-: " << AwmMean << std::endl;
	std::cout << "Mean Cw-: " << CwmMean << std::endl;
	std::cout << "Mean Aw: " << AwMean << " +- " << AwStatErr << "(stat.) " /*<< AwSysErr << "(syst.)" */<< std::endl;
	std::cout << "Mean Cw: " << CwMean << " +- " << CwStatErr << "(stat.) " /*<< CwSysErr << "(syst.)" */<< std::endl;
	exit(0);
} //eof
