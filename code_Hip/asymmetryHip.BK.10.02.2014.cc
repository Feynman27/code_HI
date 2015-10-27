// Plot the muon distributions for all nucleon
// combinations individually and combined, weighted
// by their cross-sections
// np W-: total (btilde+remnants) cross section in pb   3595.8194838902018      +-   1.6745251619307893
// np W+: total (btilde+remnants) cross section in pb   3588.8202617588995      +-   1.6283306160267506
// pn W-: total (btilde+remnants) cross section in pb   3595.9123198043762      +-   1.6801570220120536
// pn W+: total (btilde+remnants) cross section in pb   3588.4658824657063      +-   1.6237033631234190
// pp W-: total (btilde+remnants) cross section in pb   2900.7899365854364      +-   1.2652663396666060
// pp W+: total (btilde+remnants) cross section in pb   4297.1512199008384      +-   2.1086010346019695
// pp total : 64 mb

#include "TH1.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "AtlasStyle.C"

#include <iostream>
#include <string>

// Forward-Backward Asymmetry
void getFBAsymmetry(TH1F* hFB, TH1F* hf,TH1F* hb){

    std::cout << "hf: " << hf->Integral() << std::endl;
    std::cout << "hb: " << hb->Integral() << std::endl;
    std::cout << "bins: " << hFB->GetNbinsX() << std::endl;
    Int_t xmax = hf->GetNbinsX();
    std::cout << "max bins: " << xmax << std::endl;
    // Currently not including (-0.1--0.1)
    Double_t a,b,num,denom,asymm,error,da,db;
    for(Int_t i = 1, j = xmax, k = hFB->GetNbinsX(); k>0; ++i, --j, --k){
        a = hf->GetBinContent(j);
        b = hb->GetBinContent(i);
        num = a;
        denom = b;
        asymm = num/denom; 
        std::cout << "FB asymm for bin " << k  << " at |eta| = " << hFB->GetBinCenter(k) << " : " << asymm << std::endl;
        hFB->SetBinContent(k,asymm);
        if(denom<1e-20){}
        else{
            da = hf->GetBinError(j)/a;
            db = hb->GetBinError(i)/b;
            error = TMath::Sqrt(da*da + db*db ) * asymm;
            hFB->SetBinError(k,error);
        }
    }//i,j,k

}

TH1F* getChargeAsymmetry(TH1F* h1,TH1F* h2){

    TH1F* hAsymm = (TH1F*)h1->Clone();
    hAsymm->Sumw2();
    TH1F* top = (TH1F*)hAsymm->Clone();
    TH1F* bottom = (TH1F*)hAsymm->Clone();
    top->Add(h1,h2,1.0,-1.0);
    bottom->Add(h1,h2,1.0,1.0);
    hAsymm->Divide(top,bottom);
    Int_t xmax = hAsymm->GetNbinsX();
    Double_t a,b,bot,error,da,db;
    for(Int_t i = 1; i<=xmax; ++i){
        
        a = h1->GetBinContent(i);
        b = h2->GetBinContent(i);
        bot = bottom->GetBinContent(i);
        if(bot<1e-20){}
        else{
            
            da = h1->GetBinError(i);
            db = h2->GetBinError(i);
            error = 2*TMath::Sqrt(a*a*db*db + b*b*da*da )/(bot*bot);
            hAsymm->SetBinError(i,error);
        }
    }

    delete top;
    delete bottom;
    return hAsymm;
    
}

float getEta(float y, float mt){
	float m = 0.105658;
	float coshy = TMath::CosH(y);
	float sinhy = TMath::SinH(y);

	float num = TMath::Sqrt(TMath::Power(mt,2)*TMath::Power(coshy,2)-TMath::Power(m,2))+mt*sinhy;
	float den = TMath::Sqrt(TMath::Power(mt,2)*TMath::Power(coshy,2)-TMath::Power(m,2))-mt*sinhy;
	
	float etaLab = 0.5*TMath::Log(num/den);
	return etaLab;
}

float boostPhi(float pt, float eta, float phi, bool A = false, bool B = false)
{
	// rapidity shift in COM frame
	// in direction of proton
	float yBeta;
	if(A) yBeta = -1.0*0.465;
	else if(B) yBeta = 1.0*0.465;
	else return phi;

	TLorentzVector vmu;
	vmu.SetPtEtaPhiM(pt,eta,phi,0.105658);

	float Beta = TMath::TanH(yBeta);
	// Boost from COM --> ATLAS frame
	vmu.Boost(0.0,0.0,Beta);
	// return eta of muon in ATLAS frame
	float boostedPhi = (float)vmu.Phi();
	return boostedPhi;

}

float boostPt(float pt, float eta, float phi, bool A = false, bool B = false)
{
	// rapidity shift in COM frame
	// in direction of proton
	float yBeta;
	if(A) yBeta = -1.0*0.465;
	else if(B) yBeta = 1.0*0.465;
	else return pt;
	

	TLorentzVector vmu;
	vmu.SetPtEtaPhiM(pt,eta,phi,0.105658);

	float Beta = TMath::TanH(yBeta);
	// Boost from COM --> ATLAS frame
	vmu.Boost(0.0,0.0,Beta);
	// return eta of muon in ATLAS frame
	float boostedPt = (float)vmu.Pt();
	return boostedPt;

}

float boostEta(float pt, float eta, float phi, bool A = false, bool B = false)
{
	// rapidity shift in COM frame
	// in direction of proton
	float yBeta;
	if(A) yBeta = -1.0*0.465;
	else if(B) yBeta = 1.0*0.465;
	else return eta;
	

	TLorentzVector vmu;
	vmu.SetPtEtaPhiM(pt,eta,phi,0.105658);
	//float rapidityCOM = vmu.Rapidity(); // rapidity of muon in COM frame
	//float rapidityLAB = rapidityCOM-yBeta; // rapidity of muon in LAB frame
	//float etaLAB = getEta(rapidityLAB,vmu.Mt());
	//return etaLAB;

	float Beta = TMath::TanH(yBeta);
	//std::cout << "COM eta: " << eta << std::endl;
	// Boost from COM --> ATLAS frame
	vmu.Boost(0.0,0.0,Beta);
	// return eta of muon in ATLAS frame
	float boostedEta = (float)vmu.Eta();
	//std::cout << "boostedEta:" << boostedEta << std::endl;
	return boostedEta;

}

TH1F* fillHisto(TString fileName, double cs, TString hName = "", bool A = false, bool B = false){

	TChain* tree = new TChain("T","T");
	tree->Add(fileName);
	
	std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;

	int muonMother[50], nmu;
	float muonEta[50],muonPt[50],muonPhi[50],neutrinoPt[50],mt[50];
	tree->SetBranchAddress("muonMother",&muonMother);
	tree->SetBranchAddress("muonEta",&muonEta);
	tree->SetBranchAddress("muonPt",&muonPt);
	tree->SetBranchAddress("neutrinoPt",&neutrinoPt);
	tree->SetBranchAddress("muonPhi",&muonPhi);
	tree->SetBranchAddress("mt",&mt);
	tree->SetBranchAddress("nMuons",&nmu);
	
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("nMuons",1);
	tree->SetBranchStatus("muonMother",1);
	tree->SetBranchStatus("muonEta",1);
	tree->SetBranchStatus("muonPt",1);
	tree->SetBranchStatus("neutrinoPt",1);
	tree->SetBranchStatus("muonPhi",1);
	tree->SetBranchStatus("mt",1);

	TH1F* htemp = new TH1F(hName,hName,100,-6.0,6.0);
	htemp->Sumw2();
	for(int iev=0; iev<tree->GetEntries(); ++iev){

		tree->LoadTree(iev);
		tree->GetEntry(iev);

		if(iev%10000==0) std::cout << "Event: " << iev << std::endl;

		for(int imu=0; imu<nmu; ++imu){
			// Boost into lab frame (y_beta = 0.465) 
			float muonBoostedEta = boostEta(muonPt[imu],muonEta[imu],muonPhi[imu], A, B);
			//float muonBoostedPt = boostPt(muonPt[imu],muonEta[imu],muonPhi[imu], A, B);
			//float neutrinoBoostedPt = boostPt(neutrinoPt[imu],neutrinoEta[imu],neutrinoPhi[imu], A, B);
			//float muonBoostedPhi = boostPhi(muonPt[imu],muonEta[imu],muonPhi[imu], A, B);
			//float neutrinoBoostedPhi = boostPhi(neutrinoPt[imu],neutrinoEta[imu],neutrinoPhi[imu], A, B);
			//float dPhi = muonBoostedPhi-neutrinoBoostedPhi;
			//float mtBoosted = TMath::Sqrt(2.0*muonBoostedPt*neutrinoBoostedPt*(1.0-TMath::Cos(dPhi)));
			// Fill histo wt'd by sigma_W/sigma_pp 
			//if(fabs(muonMother[imu])==24) htemp->Fill(muonEta[imu]);
			if(fabs(muonMother[imu])==24 
			   && muonPt[imu] > 25.0
			   && neutrinoPt[imu] > 25.0
			   && mt[imu] > 40.0 	
			  ) {

                htemp->Fill(muonBoostedEta);
            }
		}//imu
	}//iev

	std::cout << "Number of W-->mu events: " << htemp->Integral() << std::endl;
	std::cout << "sigma_W/sigma_pp: " << cs/64.0 << std::endl;
	std::cout << "Cross-section sampled: " << htemp->Integral()/tree->GetEntries()*100.0 << "[%]" << std::endl;
	htemp->Scale(1.0/tree->GetEntries()*cs/64.0);
	return htemp;
}

void run(const std::string& period, TFile* outFile){

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
		std::cout << "Combining periods A and B. Using convention that proton is moving toward positive rapidity." << std::endl;
	}

	TString baseString = "/usatlas/u/tbales/scratch/";
	TString fileName_npWPlus = "PowhegPythia8_MSTW2008NLO_W+_np.08.29.2013.root";
	TString fileName_npWMinus = "PowhegPythia8_MSTW2008NLO_W-_np.08.29.2013.root";
	TString fileName_pnWPlus = "PowhegPythia8_MSTW2008NLO_W+_pn.08.29.2013.root";
	TString fileName_pnWMinus = "PowhegPythia8_MSTW2008NLO_W-_pn.08.29.2013.root";
	TString fileName_ppWPlus = "PowhegPythia8_MSTW2008NLO_W+_pp.08.29.2013.root" ; 	
	TString fileName_ppWMinus = "PowhegPythia8_MSTW2008NLO_W-_pp.08.29.2013.root" ;	

	// Cross-section in mb
	double crossSection_npWPlus = 3588.820*1.0e-9;
	double crossSection_npWMinus = 3595.819*1.0e-9;
	double crossSection_ppWPlus = 4297.151*1.0e-9;
	double crossSection_ppWMinus = 2900.790*1.0e-9;
	double crossSection_pnWPlus = 3588.466*1.0e-9;
	double crossSection_pnWMinus = 3595.912*1.0e-9;

	// Fill histos wt'd by cs
	TH1F* hp_np = new TH1F("hp_np","hp_np",100,-6.0,6.0);
	hp_np->Sumw2();
	TH1F* hn_np = new TH1F("hn_np","hn_np",100,-6.0,6.0);
	hn_np->Sumw2();
	TH1F* hp_ppA = new TH1F("hp_ppA","hp_ppA",100,-6.0,6.0);
	hp_ppA->Sumw2();
	TH1F* hn_ppA = new TH1F("hn_ppA","hn_ppA",100,-6.0,6.0);
	hn_ppA->Sumw2();

	TH1F* hp_pn = new TH1F("hp_pn","hp_pn",100,-6.0,6.0);
	hp_pn->Sumw2();
	TH1F* hn_pn = new TH1F("hn_pn","hn_pn",100,-6.0,6.0);
	hn_pn->Sumw2();
	TH1F* hp_ppB = new TH1F("hp_ppB","hp_ppB",100,-6.0,6.0);
	hp_ppB->Sumw2();
	TH1F* hn_ppB = new TH1F("hn_ppB","hn_ppB",100,-6.0,6.0);
	hn_ppB->Sumw2();

	TString sDataFileIn = "ResultFiles/WAnalysis_fitResultCentChrgEta_";
	if(doPeriodA) sDataFileIn+="PeriodA.09.25.2013.root";
	else if(doPeriodB) sDataFileIn+="PeriodB.09.25.2013.root";
	else {
        // Uncorrected
        //sDataFileIn+="PeriodAB_Raw.07.24.2014.root";
        // Cw correction applied
        //sDataFileIn+="PeriodAB_Corrected.09.12.2014.root";
        // Centrality binning to match Z
        // but periods in MC flipped
        //sDataFileIn+="ZBinning_PeriodAB_Corrected.09.15.2014.root";
        // No fcal re-weighting (negligible)
        sDataFileIn+="PeriodAB_Corrected.09.24.2014.root";
        //sDataFileIn+="PeriodAB_CorrectedFCalReWtd.09.24.2014.root";
	}

	TFile* dataFileIn = new TFile(sDataFileIn,"read");
    if(!dataFileIn->IsOpen()){
        std::cout << "ERROR: Cannot locate file " << sDataFileIn << std::endl;
        exit(0);
    }

	TGraphAsymmErrors* grCentrality = (TGraphAsymmErrors*)dataFileIn->Get("WPtFit_centrality");
	const int nCentralityBins = grCentrality->GetN();
	std::cout << "No. centrality bins: " << nCentralityBins << std::endl;

	TGraphAsymmErrors* grEta = (TGraphAsymmErrors*)dataFileIn->Get("WPtFit_eta");
	const int nEtaBins = grEta->GetN();
	std::cout << "No. eta bins: " << nEtaBins << std::endl;
	
	TGraphAsymmErrors* grPt = (TGraphAsymmErrors*)dataFileIn->Get("WPtFit_pt");
	const int nPtBins = grPt->GetN();
	std::cout << "No. pt bins: " << nPtBins << std::endl;

	float nWPlus[nEtaBins],nWPlusStatErr[nEtaBins],nWMinus[nEtaBins],nWMinusStatErr[nEtaBins];
	float nWPlusFw[nEtaBins],nWPlusFwStatErr[nEtaBins],nWMinusFw[nEtaBins],nWMinusFwStatErr[nEtaBins],nWFw[nEtaBins],nWFwStatErr[nEtaBins];
	float nWPlusBw[nEtaBins],nWPlusBwStatErr[nEtaBins],nWMinusBw[nEtaBins],nWMinusBwStatErr[nEtaBins],nWBw[nEtaBins],nWBwStatErr[nEtaBins];
	float xEtaBins[] = {-2.5,-2.1,-1.55,-1.05,-0.6,-0.1,0.1,0.6,1.05,1.55,2.1,2.5};
	float xAbsEtaBins[] = {0.1,0.6,1.05,1.55,2.1,2.5};
    float nAbsEtaBins = sizeof(xAbsEtaBins)/sizeof(float) - 1;

	TH1F* hDataPlus,*hDataMinus;
	TH1F* hDataFwPlus,*hDataFwMinus, *hDataFw;
	TH1F* hDataBwPlus,*hDataBwMinus, *hDataBw;

	if(doPeriodA) {
		hDataPlus = new TH1F("hDataPerAPlus","hDataPerAPlus",nEtaBins,xEtaBins);
		hDataMinus = new TH1F("hDataPerAMinus","hDataPerAMinus",nEtaBins,xEtaBins);
		hDataFwPlus = new TH1F("hDataFwPerAPlus","hDataFwPerAPlus",nEtaBins,xEtaBins);
		hDataFwMinus = new TH1F("hDataFwPerAMinus","hDataFwPerAMinus",nEtaBins,xEtaBins);
		hDataFw = new TH1F("hDataFwPerA","hDataFwPerA",nEtaBins,xEtaBins);
		hDataBwPlus = new TH1F("hDataBwPerAPlus","hDataBwPerAPlus",nEtaBins,xEtaBins);
		hDataBwMinus = new TH1F("hDataBwPerAMinus","hDataBwPerAMinus",nEtaBins,xEtaBins);
		hDataBw = new TH1F("hDataBwPerA","hDataBwPerA",nEtaBins,xEtaBins);
	}
	else if(doPeriodB) {
		hDataPlus = new TH1F("hDataPerBPlus","hDataPerBPlus",nEtaBins,xEtaBins);
		hDataMinus = new TH1F("hDataPerBMinus","hDataPerBMinus",nEtaBins,xEtaBins);
		hDataFwPlus = new TH1F("hDataFwPerBPlus","hDataFwPerBPlus",nEtaBins,xEtaBins);
		hDataFwMinus = new TH1F("hDataFwPerBMinus","hDataFwPerBMinus",nEtaBins,xEtaBins);
		hDataFw = new TH1F("hDataFwPerB","hDataFwPerB",nEtaBins,xEtaBins);
		hDataBwPlus = new TH1F("hDataBwPerBPlus","hDataBwPerBPlus",nEtaBins,xEtaBins);
		hDataBwMinus = new TH1F("hDataBwPerBMinus","hDataBwPerBMinus",nEtaBins,xEtaBins);
		hDataBw = new TH1F("hDataBwPerB","hDataBwPerB",nEtaBins,xEtaBins);
	}
	else {
		hDataPlus = new TH1F("hDataPlus","hDataPlus",nEtaBins,xEtaBins);
		hDataMinus = new TH1F("hDataMinus","hDataMinus",nEtaBins,xEtaBins);
		hDataFwPlus = new TH1F("hDataFwPlus","hDataFwPlus",nEtaBins,xEtaBins);
		hDataFwMinus = new TH1F("hDataFwMinus","hDataFwMinus",nEtaBins,xEtaBins);
		hDataFw = new TH1F("hDataFw","hDataFw",nEtaBins,xEtaBins);
		hDataBwPlus = new TH1F("hDataBwPlus","hDataBwPlus",nEtaBins,xEtaBins);
		hDataBwMinus = new TH1F("hDataBwMinus","hDataBwMinus",nEtaBins,xEtaBins);
		hDataBw = new TH1F("hDataBw","hDataBw",nEtaBins,xEtaBins);
	}
	hDataPlus->Sumw2();
	hDataMinus->Sumw2();
	hDataFwPlus->Sumw2();
	hDataFwMinus->Sumw2();
	hDataFw->Sumw2();
	hDataBwPlus->Sumw2();
	hDataBwMinus->Sumw2();
	hDataBw->Sumw2();

	// Add up signal candidates from all centrality classes for
	// each eta bin in the data
	for(int ieta=0; ieta<nEtaBins; ++ieta){
	  // Initialize arrays to 0
	  nWPlus[ieta] = 0;
	  nWPlusStatErr[ieta] = 0;
	  nWMinus[ieta] = 0;
	  nWMinusStatErr[ieta] = 0;
	  nWPlusFw[ieta] = 0;
	  nWPlusFwStatErr[ieta] = 0;
	  nWMinusFw[ieta] = 0;
	  nWMinusFwStatErr[ieta] = 0;
	  nWFw[ieta] = 0;
	  nWFwStatErr[ieta] = 0;
	  nWPlusBw[ieta] = 0;
	  nWPlusBwStatErr[ieta] = 0;
	  nWMinusBw[ieta] = 0;
	  nWMinusBwStatErr[ieta] = 0;
	  nWBw[ieta] = 0;
	  nWBwStatErr[ieta] = 0;

	  for(int icent=0; icent<nCentralityBins; ++icent){
	    for(int ipt=0; ipt<nPtBins; ++ipt){
            // Retrieve yields
            TString name= "WPtFit_eta"; name+=ieta; name+="_centrality"; name+=icent; name+="_charge";
            TString namePlus = name+"Plus";
            TString nameMinus = name+"Minus";

            TGraphAsymmErrors* grPlusTemp = (TGraphAsymmErrors*)dataFileIn->Get(namePlus); 
            nWPlus[ieta]+=grPlusTemp->GetY()[ipt];
            nWPlusStatErr[ieta]+=TMath::Power(grPlusTemp->GetEYhigh()[ipt],2);

            TGraphAsymmErrors* grMinusTemp = (TGraphAsymmErrors*)dataFileIn->Get(nameMinus); 
            nWMinus[ieta]+=grMinusTemp->GetY()[ipt];
            nWMinusStatErr[ieta]+=TMath::Power(grMinusTemp->GetEYhigh()[ipt],2);

            if(grEta->GetX()[ieta] < 0.0) {
                // W+
                nWPlusBw[ieta]+=grPlusTemp->GetY()[ipt];
                nWPlusBwStatErr[ieta]+=TMath::Power(grPlusTemp->GetEYhigh()[ipt],2);
                // W-
                nWMinusBw[ieta]+=grMinusTemp->GetY()[ipt];
                nWMinusBwStatErr[ieta]+=TMath::Power(grMinusTemp->GetEYhigh()[ipt],2);
            }
            else if(grEta->GetX()[ieta] > 0.0){
                // W+
                nWPlusFw[ieta]+=grPlusTemp->GetY()[ipt];
                nWPlusFwStatErr[ieta]+=TMath::Power(grPlusTemp->GetEYhigh()[ipt],2);
                // W-
                nWMinusFw[ieta]+=grMinusTemp->GetY()[ipt];
                nWMinusFwStatErr[ieta]+=TMath::Power(grMinusTemp->GetEYhigh()[ipt],2);

            }
	    }//ipt   
		
	  } //icent

	  nWPlusStatErr[ieta]    = TMath::Sqrt(nWPlusStatErr[ieta]);
	  nWMinusStatErr[ieta]   = TMath::Sqrt(nWMinusStatErr[ieta]);
	  nWPlusBwStatErr[ieta]  = TMath::Sqrt(nWPlusBwStatErr[ieta]);
	  nWMinusBwStatErr[ieta] = TMath::Sqrt(nWMinusBwStatErr[ieta]);
	  nWPlusFwStatErr[ieta]  = TMath::Sqrt(nWPlusFwStatErr[ieta]);
	  nWMinusFwStatErr[ieta] = TMath::Sqrt(nWMinusFwStatErr[ieta]);

      // Charge inclusive F/W asymmetry
      if(grEta->GetX()[ieta] < 0.0) {
          // Backward 
          nWBw[ieta] = nWPlusBw[ieta] + nWMinusBw[ieta];
          nWBwStatErr[ieta] = TMath::Sqrt( TMath::Power(nWPlusBwStatErr[ieta],2) + TMath::Power(nWMinusBwStatErr[ieta],2) );
      }
      else if(grEta->GetX()[ieta] > 0.0){
          // Forward 
          nWFw[ieta] = nWPlusFw[ieta] + nWMinusFw[ieta];
          nWFwStatErr[ieta] = TMath::Sqrt( TMath::Power(nWPlusFwStatErr[ieta],2) + TMath::Power(nWMinusFwStatErr[ieta],2) );
      }


	  float binWOver2 = grEta->GetEXhigh()[ieta];
	  float binW = grEta->GetEXhigh()[ieta]*2.0;

	  std::cout << "Total W+ --> mu+ events in eta region " << grEta->GetX()[ieta]-binWOver2 << " - " << grEta->GetX()[ieta]+binWOver2 << 
			": " << nWPlus[ieta] << " +/-" << nWPlusStatErr[ieta] << "(stat.)" << std::endl; 
	  std::cout << "Total W- --> mu- events in eta region " << grEta->GetX()[ieta]-binWOver2 << " - " << grEta->GetX()[ieta]+binWOver2 << 
			": " << nWMinus[ieta] << " +/-" << nWMinusStatErr[ieta] << "(stat.)" << std::endl; 

	  hDataPlus->SetBinContent(ieta+1,nWPlus[ieta]/binW);
	  hDataPlus->SetBinError(ieta+1,nWPlusStatErr[ieta]/binW);

	  hDataMinus->SetBinContent(ieta+1,nWMinus[ieta]/binW);
	  hDataMinus->SetBinError(ieta+1,nWMinusStatErr[ieta]/binW);

      // Backward
      // W+
	  hDataBwPlus->SetBinContent(ieta+1,nWPlusBw[ieta]/binW);
	  hDataBwPlus->SetBinError(ieta+1,nWPlusBwStatErr[ieta]/binW);
      // W-
	  hDataBwMinus->SetBinContent(ieta+1,nWMinusBw[ieta]/binW);
	  hDataBwMinus->SetBinError(ieta+1,nWMinusBwStatErr[ieta]/binW);
      // Charge inclusive
	  hDataBw->SetBinContent(ieta+1,nWBw[ieta]/binW);
	  hDataBw->SetBinError(ieta+1,nWBwStatErr[ieta]/binW);

      // Forward
      // W+
	  hDataFwPlus->SetBinContent(ieta+1,nWPlusFw[ieta]/binW);
	  hDataFwPlus->SetBinError(ieta+1,nWPlusFwStatErr[ieta]/binW);
      // W-
	  hDataFwMinus->SetBinContent(ieta+1,nWMinusFw[ieta]/binW);
	  hDataFwMinus->SetBinError(ieta+1,nWMinusFwStatErr[ieta]/binW);
      // Charge inclusive
	  hDataFw->SetBinContent(ieta+1,nWFw[ieta]/binW);
	  hDataFw->SetBinError(ieta+1,nWFwStatErr[ieta]/binW);

	}//ieta

	//Period A
	hp_np = fillHisto(baseString+fileName_npWPlus,crossSection_npWPlus,hp_np->GetName(),true,false); 	
	hn_np = fillHisto(baseString+fileName_npWMinus,crossSection_npWMinus,hn_np->GetName(),true,false); 	
	hp_ppA = fillHisto(baseString+fileName_ppWPlus,crossSection_ppWPlus,hp_ppA->GetName(),true,false); 	
	hn_ppA = fillHisto(baseString+fileName_ppWMinus,crossSection_ppWMinus,hn_ppA->GetName(),true,false); 	

	// Period B
    // Use this convention when plotting both periods
	hp_ppB = fillHisto(baseString+fileName_ppWPlus,crossSection_ppWPlus,hp_ppB->GetName(),false,true); 	
	hn_ppB = fillHisto(baseString+fileName_ppWMinus,crossSection_ppWMinus,hn_ppB->GetName(),false,true); 	
	hp_pn = fillHisto(baseString+fileName_pnWPlus,crossSection_pnWPlus,hp_pn->GetName(),false,true); 	
	hn_pn = fillHisto(baseString+fileName_pnWMinus,crossSection_pnWMinus,hn_pn->GetName(),false,true); 	

	// Weight pp and np by collisions frequency
	double ppFreq = 82.0/208.0;
	double npFreq = 126.0/208.0;
/*	hp_pp->Scale(ppFreq);
	hn_pp->Scale(ppFreq);
	hp_np->Scale(npFreq);
	hn_np->Scale(npFreq);
*/
	// Add together to simulate A+p (Period A)
	TH1F* hp_Ap  = (TH1F*)hp_ppA->Clone("hp_Ap");
	hp_Ap->Sumw2();
	hp_Ap->Add(hp_ppA,hp_np,ppFreq,npFreq);

	TH1F* hn_Ap  = (TH1F*)hn_ppA->Clone("hn_Ap");
	hn_Ap->Sumw2();
	hn_Ap->Add(hn_ppA,hn_np,ppFreq,npFreq);

	// Add together to simulate p+A (Period B)
	TH1F* hp_pA  = (TH1F*)hp_ppB->Clone("hp_pA");
	hp_pA->Sumw2();
	hp_pA->Add(hp_ppB,hp_pn,ppFreq,npFreq);

	TH1F* hn_pA  = (TH1F*)hn_ppB->Clone("hn_pA");
	hn_pA->Sumw2();
	hn_pA->Add(hn_ppB,hn_pn,ppFreq,npFreq);

	// Charge asymmetry
	TH1F* hp_pAc = (TH1F*)hp_pA->Clone("hp_pAc");
	TH1F* hn_pAc = (TH1F*)hn_pA->Clone("hn_pAc");
	TH1F* hChAsymm_pA = new TH1F("hChAsymm_pA","hChAsymm_pA",100,-6.0,6.0);
	hChAsymm_pA->Sumw2();
	hChAsymm_pA = getChargeAsymmetry(hp_pAc,hn_pAc);	

	TH1F* hp_Apc = (TH1F*)hp_Ap->Clone("hp_Apc");
	TH1F* hn_Apc = (TH1F*)hn_Ap->Clone("hn_Apc");
	TH1F* hChAsymm_Ap = new TH1F("hChAsymm_Ap","hChAsymm_Ap",100,-6.0,6.0);
	hChAsymm_Ap->Sumw2();
	hChAsymm_Ap = getChargeAsymmetry(hp_Apc,hn_Apc);	

	TH1F* hp, *hn ;
	hp = (TH1F*)hDataPlus->Clone("hp");
	hn = (TH1F*)hDataMinus->Clone("hn");
	TH1F* hChAsymm = new TH1F("hChAsymm","Charge Asymmetry",nEtaBins,xEtaBins);
	hChAsymm->Sumw2();
	hChAsymm = getChargeAsymmetry(hp,hn);

    // Forward-Backward asymmetry
    TH1F* hpF, *hpB, *hnF, *hnB, *hF, *hB;
    // W+
    hpF = (TH1F*)hDataFwPlus->Clone("hpF");
    hpB = (TH1F*)hDataBwPlus->Clone("hpB");
	TH1F* hFBAsymmPlus = new TH1F("hFBAsymmPlus","For/Back Asymmetry, mu^{+}",nAbsEtaBins, xAbsEtaBins);
    getFBAsymmetry(hFBAsymmPlus,hpF,hpB);
    // W-
    hnF = (TH1F*)hDataFwMinus->Clone("hnF");
    hnB = (TH1F*)hDataBwMinus->Clone("hnB");
	TH1F* hFBAsymmMinus = new TH1F("hFBAsymmMinus","For/Back Asymmetry, mu^{-}",nAbsEtaBins, xAbsEtaBins);
    getFBAsymmetry(hFBAsymmMinus,hnF,hnB);
    // Charge inclusive
    hF = (TH1F*)hDataFw->Clone("hF");
    hB = (TH1F*)hDataBw->Clone("hB");
	TH1F* hFBAsymm = new TH1F("hFBAsymm","For/Back Asymmetry, mu^{pm}",nAbsEtaBins, xAbsEtaBins);
    getFBAsymmetry(hFBAsymm,hF,hB);

	// Save histos
	outFile->cd();
	hDataPlus->Write();
	hDataMinus->Write();
	if(doPeriodA){
		hp_ppA->Write();
		hn_ppA->Write();
		hp_np->Write();
		hn_np->Write();
		hp_Ap->Write();
		hn_Ap->Write();
		hChAsymm_Ap->Write("hChAsymm_Ap");
		hChAsymm->Write("hChAsymmPerA");
		hFBAsymmPlus->Write("hFBAsymmPlusPerA");
		hFBAsymmMinus->Write("hFBAsymmMinusPerA");
		hFBAsymm->Write("hFBAsymmPerA");
	}
	else if(doPeriodB){
		hp_ppB->Write();
		hn_ppB->Write();
		hp_pn->Write();
		hn_pn->Write();
		hp_pA->Write();
		hn_pA->Write();
		hChAsymm_pA->Write("hChAsymm_pA");
		hChAsymm->Write("hChAsymmPerB");
		hFBAsymmPlus->Write("hFBAsymmPlusPerB");
		hFBAsymmMinus->Write("hFBAsymmMinusPerB");
		hFBAsymm->Write("hFBAsymmPerB");
	}
	else {
		/*hp_ppB->Write();
		hn_ppB->Write();
		hp_pn->Write();
		hn_pn->Write();
		hp_pA->Write();
		hn_pA->Write();
		hChAsymm_pA->Write("hChAsymm_pA");
        */
		hChAsymm->Write("hChAsymmPerAB");
		hFBAsymmPlus->Write("hFBAsymmPlusPerAB");
		hFBAsymmMinus->Write("hFBAsymmMinusPerAB");
		hFBAsymm->Write("hFBAsymmPerAB");
	}

}

int main(){

    SetAtlasStyle();
	TFile* outFile = new TFile("HipChargeAsymmetryHistos.root","recreate");
	//run("A",outFile);	
	//run("B",outFile);	
	run("AB",outFile);	
	outFile->Close();
	return 0;
}
