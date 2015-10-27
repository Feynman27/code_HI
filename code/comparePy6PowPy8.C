#include "TH1.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TF1.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <string>

void fillGenHisto(TString baseString, TString fileNameIn, TH1F* h, bool isMuPlus = false,bool isMuMinus = false){

      h->Sumw2();
      TChain* tree = new TChain("tree","tree");
      tree->Add(baseString+fileNameIn);
      std::cout << "Pathname: " << baseString+fileNameIn << std::endl;
      std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;

      float nuGenPtNt[50], nuPhiGenNt[50];
      float muGenPtNt[50], muPhiGenNt[50];
      float mtGenNt[50];
      int motherNt[50];
      int daughterNt[50];
      int daughterNuNt[50];
      float muEtaGenNt[50];
      float nuEtaGenNt[50];
      float centralityNt;
      float chargeGenNt[50];
      int nmu;



       // --- Set branch adresses ---
  tree->SetBranchAddress("mc_mu_gen_pt", &muGenPtNt);
  tree->SetBranchAddress("mc_nu_gen_pt", &nuGenPtNt);
  tree->SetBranchAddress("mc_mu_gen_mothertype", &motherNt);
  tree->SetBranchAddress("mc_mu_gen_type", &daughterNt);
  tree->SetBranchAddress("mc_nu_gen_type", &daughterNuNt);
  tree->SetBranchAddress("mc_mu_charge", &chargeGenNt);
  tree->SetBranchAddress("mc_mu_gen_eta", &muEtaGenNt);
  tree->SetBranchAddress("mc_nu_gen_eta", &nuEtaGenNt);
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
  tree->SetBranchStatus("mc_nu_gen_type", 1);
  tree->SetBranchStatus("mc_mu_charge", 1);
  tree->SetBranchStatus("mc_mu_gen_eta", 1);
  tree->SetBranchStatus("mc_nu_gen_eta", 1);
  tree->SetBranchStatus("mc_mu_gen_phi", 1);
  tree->SetBranchStatus("mc_nu_gen_phi", 1);


	for(int iev=0; iev<tree->GetEntries(); ++iev){

	    if(iev%10000==0) std::cout << "Event: " << iev << std::endl;
	    tree->LoadTree(iev);
	    tree->GetEntry(iev);

		for(int imu=0; imu<nmu; ++imu){
			if(
			  (fabs(motherNt[imu])==24) 
			  && muGenPtNt[imu]>4.0
			  //&& pt[imu]>25.0
			  //&& mpt>25.0 && mpt<9000.
			  //&& mt[imu]>40.0
			){
			   if(isMuPlus&&chargeGenNt[imu]>0.){
				h->Fill(muEtaGenNt[imu]);
			   }
			   else if(isMuMinus&&chargeGenNt[imu]<0.0){
				h->Fill(muEtaGenNt[imu]);
			   }
		           else continue;
			}
		}//imu
	}//iev

}

void fillRecHisto(TString baseString, TString fileNameIn, TH1F* h, bool isMuPlus = false,bool isMuMinus = false, bool doMc=true){

      h->Sumw2();
      TChain* tree = new TChain("tree","tree");
      tree->Add(baseString+fileNameIn);
      std::cout << "Pathname: " << baseString+fileNameIn << std::endl;
      std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;

      ///reco level variables
      float eLoss[50];
      float scat[50];
      float comp[50];
      float pt[50];
      float mt[50];
      float eta[50];
      float phi[50];
      float mptPhi;
      float charge[50];
      int prompt[50];
      float centrality;
      float mpt;
      float ptcone[50];
      int val[50], truthMatched[50],ZDY[50], matched1[50], matched2[50], matched3[50],matched4[50],matched5[50];
      int nmu,trig1,trig2,trig3,trig4,trig5, mbtrig1, mbtrig2;


      // Set branch address
      tree->SetBranchAddress("pt",&pt);
      if(doMc)tree->SetBranchAddress("prompt",&prompt);
      tree->SetBranchAddress("eLoss", &eLoss);
      tree->SetBranchAddress("ptcone20ID3", &ptcone);
      tree->SetBranchAddress("scat", &scat);
      tree->SetBranchAddress("comp", &comp);
      tree->SetBranchAddress("pt", &pt);
      tree->SetBranchAddress("mt", &mt);
      tree->SetBranchAddress("eta", &eta);
      tree->SetBranchAddress("phi", &phi);
      tree->SetBranchAddress("nu_phi", &mptPhi);
      tree->SetBranchAddress("charge", &charge);
      tree->SetBranchAddress("val", &val); 
      tree->SetBranchAddress("ZDY", &ZDY); 
      tree->SetBranchAddress("centrality", &centrality);
      tree->SetBranchAddress("nu_pt", &mpt);
      tree->SetBranchAddress("mu_muid_n", &nmu);
      if(!doMc){
	  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC",&trig1);
	  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10",&trig2);
	  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20",&trig3);
	  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC_Matched20",&matched1);
	  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10_Matched20",&matched2);
	  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20_Matched20",&matched3);
	  tree->SetBranchAddress("EF_mu4_MSonly_L1TE50",&trig4);
	  tree->SetBranchAddress("EF_mu4_L1VTE50",&trig5);
	  tree->SetBranchAddress("EF_mu4_MSonly_L1TE50_Matched20",&matched4);
	  tree->SetBranchAddress("EF_mu4_L1VTE50_Matched20",&matched5);
  	  tree->SetBranchAddress("EF_mbZdc_a_c_L1VTE50_trk",&mbtrig1);
  	  tree->SetBranchAddress("EF_L1TE50_NoAlg",&mbtrig2);
      }

      // --- Set branch status ---
      tree->SetBranchStatus("*",0) ;
      tree->SetBranchStatus("mu_muid_n", 1);
      tree->SetBranchStatus("eLoss", 1);
      tree->SetBranchStatus("ptcone20ID3", 1);
      tree->SetBranchStatus("scat", 1);
      tree->SetBranchStatus("comp", 1);
      if(doMc)tree->SetBranchStatus("prompt",1);
      tree->SetBranchStatus("pt", 1);
      tree->SetBranchStatus("mt", 1);
      tree->SetBranchStatus("eta", 1);
      tree->SetBranchStatus("phi", 1);
      tree->SetBranchStatus("nu_phi", 1);
      tree->SetBranchStatus("charge", 1);
      tree->SetBranchStatus("val", 1); 
      tree->SetBranchStatus("ZDY", 1); 
      tree->SetBranchStatus("centrality", 1);
      tree->SetBranchStatus("nu_pt", 1);
      if(!doMc) { 
	tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC",1);
  	tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10",1);
  	tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20",1);
  	tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC_Matched20",1);
  	tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10_Matched20",1);
  	tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20_Matched20",1);
  	tree->SetBranchStatus("EF_mu4_MSonly_L1TE50",1);
  	tree->SetBranchStatus("EF_mu4_L1VTE50",1);
  	tree->SetBranchStatus("EF_mu4_MSonly_L1TE50_Matched20",1);
  	tree->SetBranchStatus("EF_mu4_L1VTE50_Matched20",1);
  	tree->SetBranchStatus("EF_mbZdc_a_c_L1VTE50_trk",1);
  	tree->SetBranchStatus("EF_L1TE50_NoAlg",1);
      }
 

	for(int iev=0; iev<tree->GetEntries(); ++iev){

	    if(iev%10000==0) std::cout << "Event: " << iev << std::endl;
	    tree->LoadTree(iev);
	    tree->GetEntry(iev);

		for(int imu=0; imu<nmu; ++imu){
			if(
			  (prompt[imu]==24||!doMc) 
    			  //( ( (trig1&&matched1[imu]) || (trig2&&matched2[imu]) || (trig3&&matched3[imu]) ) || doMc)
    			  //( ( mbtrig1 || mbtrig2) || doMc)
    			  //( ( (trig4&&matched4[imu]) || (trig5&&matched5[imu])) || doMc)
			  &&val[imu]>11
			  //val[imu]>0
			  && fabs(scat[imu])<4.0
			  && fabs(eLoss[imu])<0.5
			  //&& ptcone[imu]/pt[imu] < 0.1
			  //&& ZDY[imu]==0
			  //&& pt[imu]>10.0
			  && pt[imu]>4.0
			  //&& pt[imu]>25.0
			  //&& mpt>25.0 && mpt<9000.
			  //&& mt[imu]>40.0
			){
			   if(isMuPlus&&charge[imu]>0.){
				h->Fill(eta[imu]);
			   }
			   else if(isMuMinus&&charge[imu]<0.0){
				h->Fill(eta[imu]);
			   }
		           else continue;
			}
		}//imu
	}//iev

}

void comparePy6PowPy8(){

	TString baseString = "/usatlas/u/tbales/scratch/";
	TString fInPy6 = "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.04.12.2013.root";
	//PowPy8 samples
	//pp
	TString fPlusPowPy8_pp = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_pp.10.03.2013.root";
	TString fMinusPowPy8_pp = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_pp.10.03.2013.root";
	TString fPlusPowPy8_np = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_np.10.03.2013.root";
	TString fMinusPowPy8_np = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_np.10.03.2013.root";
	TString fPlusPowPy8_pn = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_pn.10.03.2013.root";
	TString fMinusPowPy8_pn = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_pn.10.03.2013.root";
	TString fPlusPowPy8_nn = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_nn.10.03.2013.root";
	TString fMinusPowPy8_nn = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_nn.10.03.2013.root";
	
	TH1F *hPlusPowPy8_pp = new TH1F("hPlusPowPy8_pp","hPlusPowPy8_pp",50,-2.5,2.5); 
	hPlusPowPy8_pp->Sumw2();
	TH1F *hMinusPowPy8_pp = new TH1F("hMinusPowPy8_pp","hMinusPowPy8_pp",50,-2.5,2.5); 
	hMinusPowPy8_pp->Sumw2();
	TH1F *hPlusPowPy8_np = new TH1F("hPlusPowPy8_np","hPlusPowPy8_np",50,-2.5,2.5); 
	hPlusPowPy8_np->Sumw2();
	TH1F *hMinusPowPy8_np = new TH1F("hMinusPowPy8_np","hMinusPowPy8_np",50,-2.5,2.5); 
	hMinusPowPy8_np->Sumw2();
	TH1F *hPlusPowPy8_pn = new TH1F("hPlusPowPy8_pn","hPlusPowPy8_pn",50,-2.5,2.5); 
	hPlusPowPy8_pn->Sumw2();
	TH1F *hMinusPowPy8_pn = new TH1F("hMinusPowPy8_pn","hMinusPowPy8_pn",50,-2.5,2.5); 
	hMinusPowPy8_pn->Sumw2();
	TH1F *hPlusPowPy8_nn = new TH1F("hPlusPowPy8_nn","hPlusPowPy8_nn",50,-2.5,2.5); 
	hPlusPowPy8_nn->Sumw2();
	TH1F *hMinusPowPy8_nn = new TH1F("hMinusPowPy8_nn","hMinusPowPy8_nn",50,-2.5,2.5); 
	hMinusPowPy8_nn->Sumw2();

	TH1F *hPlusPowPy8 = new TH1F("hPlusPowPy8","hPlusPowPy8",50,-2.5,2.5);
	hPlusPowPy8->Sumw2();
	TH1F *hMinusPowPy8 = new TH1F("hMinusPowPy8","hMinusPowPy8",50,-2.5,2.5);
	hPlusPowPy8->Sumw2();

	TH1F *hPlusPy6 = new TH1F("hPlusPy6","hPlusPy6",50,-2.5,2.5);
	hPlusPy6->Sumw2();
	TH1F *hMinusPy6 = new TH1F("hMinusPy6","hMinusPy6",50,-2.5,2.5);
	hPlusPy6->Sumw2();
	
	// Generator level
	TH1F *hGenPlusPowPy8_pp = new TH1F("hGenPlusPowPy8_pp","hGenPlusPowPy8_pp",50,-2.5,2.5); 
	hGenPlusPowPy8_pp->Sumw2();
	TH1F *hGenMinusPowPy8_pp = new TH1F("hGenMinusPowPy8_pp","hGenMinusPowPy8_pp",50,-2.5,2.5); 
	hGenMinusPowPy8_pp->Sumw2();
	TH1F *hGenPlusPowPy8_np = new TH1F("hGenPlusPowPy8_np","hGenPlusPowPy8_np",50,-2.5,2.5); 
	hGenPlusPowPy8_np->Sumw2();
	TH1F *hGenMinusPowPy8_np = new TH1F("hGenMinusPowPy8_np","hGenMinusPowPy8_np",50,-2.5,2.5); 
	hGenMinusPowPy8_np->Sumw2();
	TH1F *hGenPlusPowPy8_pn = new TH1F("hGenPlusPowPy8_pn","hGenPlusPowPy8_pn",50,-2.5,2.5); 
	hGenPlusPowPy8_pn->Sumw2();
	TH1F *hGenMinusPowPy8_pn = new TH1F("hGenMinusPowPy8_pn","hGenMinusPowPy8_pn",50,-2.5,2.5); 
	hGenMinusPowPy8_pn->Sumw2();
	TH1F *hGenPlusPowPy8_nn = new TH1F("hGenPlusPowPy8_nn","hGenPlusPowPy8_nn",50,-2.5,2.5); 
	hGenPlusPowPy8_nn->Sumw2();
	TH1F *hGenMinusPowPy8_nn = new TH1F("hGenMinusPowPy8_nn","hGenMinusPowPy8_nn",50,-2.5,2.5); 
	hGenMinusPowPy8_nn->Sumw2();

	TH1F *hGenPlusPowPy8 = new TH1F("hGenPlusPowPy8","hGenPlusPowPy8",50,-2.5,2.5);
	hGenPlusPowPy8->Sumw2();
	TH1F *hGenMinusPowPy8 = new TH1F("hGenMinusPowPy8","hGenMinusPowPy8",50,-2.5,2.5);
	hGenPlusPowPy8->Sumw2();

	TH1F *hGenPlusPy6 = new TH1F("hGenPlusPy6","hGenPlusPy6",50,-2.5,2.5);
	hGenPlusPy6->Sumw2();
	TH1F *hGenMinusPy6 = new TH1F("hGenMinusPy6","hGenMinusPy6",50,-2.5,2.5);
	hGenPlusPy6->Sumw2();

	float wtpp = 0.155;
	float csPlusPowPy8_pp=2.1120e-9; 
	float csMinusPowPy8_pp=1.2420e-9; 
	float wtnp = 0.478/2.0;
	float csPlusPowPy8_np=1.6690e-9; 
	float csMinusPowPy8_np=1.6540e-9; 
	float wtpn = 0.478/2.0;
	float csPlusPowPy8_pn=1.6700e-9; 
	float csMinusPowPy8_pn=1.6540e-9; 
	float wtnn = 0.367;
	float csPlusPowPy8_nn=1.2530e-9; 
	float csMinusPowPy8_nn=2.0920e-9; 
        // Fill histos	
	fillRecHisto(baseString,fPlusPowPy8_pp,hPlusPowPy8_pp,true);
	fillRecHisto(baseString,fMinusPowPy8_pp,hMinusPowPy8_pp,false,true);
	fillRecHisto(baseString,fPlusPowPy8_np,hPlusPowPy8_np,true);
	fillRecHisto(baseString,fMinusPowPy8_np,hMinusPowPy8_np,false,true);
	fillRecHisto(baseString,fPlusPowPy8_pn,hPlusPowPy8_pn, true);
	fillRecHisto(baseString,fMinusPowPy8_pn,hMinusPowPy8_pn, false,true);
	fillRecHisto(baseString,fPlusPowPy8_nn,hPlusPowPy8_nn, true);
	fillRecHisto(baseString,fMinusPowPy8_nn,hMinusPowPy8_nn, false,true);

	float nTotWPlusPowPy8 = hPlusPowPy8_pp->Integral()+hPlusPowPy8_np->Integral()+hPlusPowPy8_pn->Integral()+hPlusPowPy8_nn->Integral();
	hPlusPowPy8->Add(hPlusPowPy8_pp,hPlusPowPy8_np,wtpp*nTotWPlusPowPy8/hPlusPowPy8_pp->Integral(),wtnp*nTotWPlusPowPy8/hPlusPowPy8_np->Integral());
	hPlusPowPy8->Add(hPlusPowPy8_pn,wtpn*nTotWPlusPowPy8/hPlusPowPy8_pn->Integral());
	hPlusPowPy8->Add(hPlusPowPy8_nn,wtnn*nTotWPlusPowPy8/hPlusPowPy8_nn->Integral());
	hPlusPowPy8->Scale(1.0/hPlusPowPy8->Integral());

	float nTotWMinusPowPy8 = hMinusPowPy8_pp->Integral()+hMinusPowPy8_np->Integral()+hMinusPowPy8_pn->Integral()+hMinusPowPy8_nn->Integral();
	hMinusPowPy8->Add(hMinusPowPy8_pp,hMinusPowPy8_np,wtpp*nTotWMinusPowPy8/hMinusPowPy8_pp->Integral(),wtnp*nTotWMinusPowPy8/hMinusPowPy8_np->Integral());
	hMinusPowPy8->Add(hMinusPowPy8_pn,wtpn*nTotWMinusPowPy8/hMinusPowPy8_pn->Integral());
	hMinusPowPy8->Add(hMinusPowPy8_nn,wtnn*nTotWMinusPowPy8/hMinusPowPy8_nn->Integral());
	hMinusPowPy8->Scale(1.0/hMinusPowPy8->Integral());
	
	fillRecHisto(baseString,fInPy6,hPlusPy6,true);
	hPlusPy6->Scale(1.0/hPlusPy6->Integral());
	
	fillRecHisto(baseString,fInPy6,hMinusPy6,false,true);
	hMinusPy6->Scale(1.0/hMinusPy6->Integral());

        hPlusPowPy8->SetLineColor(7);
        hPlusPy6->SetLineColor(4);
        hMinusPowPy8->SetLineColor(6);
        hMinusPy6->SetLineColor(2);

	// Generator level
	fillGenHisto(baseString,fPlusPowPy8_pp,hGenPlusPowPy8_pp,true);
	fillGenHisto(baseString,fMinusPowPy8_pp,hGenMinusPowPy8_pp,false,true);
	fillGenHisto(baseString,fPlusPowPy8_np,hGenPlusPowPy8_np,true);
	fillGenHisto(baseString,fMinusPowPy8_np,hGenMinusPowPy8_np,false,true);
	fillGenHisto(baseString,fPlusPowPy8_pn,hGenPlusPowPy8_pn, true);
	fillGenHisto(baseString,fMinusPowPy8_pn,hGenMinusPowPy8_pn, false,true);
	fillGenHisto(baseString,fPlusPowPy8_nn,hGenPlusPowPy8_nn, true);
	fillGenHisto(baseString,fMinusPowPy8_nn,hGenMinusPowPy8_nn, false,true);

	nTotWPlusPowPy8 = hGenPlusPowPy8_pp->Integral()+hGenPlusPowPy8_np->Integral()+hGenPlusPowPy8_pn->Integral()+hGenPlusPowPy8_nn->Integral();
	hGenPlusPowPy8->Add(hGenPlusPowPy8_pp,hGenPlusPowPy8_np,wtpp*nTotWPlusPowPy8/hGenPlusPowPy8_pp->Integral(),wtnp*nTotWPlusPowPy8/hGenPlusPowPy8_np->Integral());
	hGenPlusPowPy8->Add(hGenPlusPowPy8_pn,wtpn*nTotWPlusPowPy8/hGenPlusPowPy8_pn->Integral());
	hGenPlusPowPy8->Add(hGenPlusPowPy8_nn,wtnn*nTotWPlusPowPy8/hGenPlusPowPy8_nn->Integral());
	hGenPlusPowPy8->Scale(1.0/hGenPlusPowPy8->Integral());

	nTotWMinusPowPy8 = hGenMinusPowPy8_pp->Integral()+hGenMinusPowPy8_np->Integral()+hGenMinusPowPy8_pn->Integral()+hGenMinusPowPy8_nn->Integral();
	hGenMinusPowPy8->Add(hGenMinusPowPy8_pp,hGenMinusPowPy8_np,wtpp*nTotWMinusPowPy8/hGenMinusPowPy8_pp->Integral(),wtnp*nTotWMinusPowPy8/hGenMinusPowPy8_np->Integral());
	hGenMinusPowPy8->Add(hGenMinusPowPy8_pn,wtpn*nTotWMinusPowPy8/hGenMinusPowPy8_pn->Integral());
	hGenMinusPowPy8->Add(hGenMinusPowPy8_nn,wtnn*nTotWMinusPowPy8/hGenMinusPowPy8_nn->Integral());
	hGenMinusPowPy8->Scale(1.0/hGenMinusPowPy8->Integral());
	
	fillGenHisto(baseString,fInPy6,hGenPlusPy6,true);
	hGenPlusPy6->Scale(1.0/hGenPlusPy6->Integral());
	
	fillGenHisto(baseString,fInPy6,hGenMinusPy6,false,true);
	hGenMinusPy6->Scale(1.0/hGenMinusPy6->Integral());

        hGenPlusPowPy8->SetLineColor(7);
        hGenPlusPy6->SetLineColor(4);
        hGenMinusPowPy8->SetLineColor(6);
        hGenMinusPy6->SetLineColor(2);

	TCanvas *cPlus = new TCanvas("cplus","cplus",600,600);
        hPlusPowPy8->SetMarkerStyle(2);
	hPlusPowPy8->Draw("histe");
        hPlusPy6->SetMarkerStyle(2);
	hPlusPy6->Draw("histesame");

	TCanvas *cMinus = new TCanvas("cminus","cminus",600,600);
	hMinusPowPy8->SetMarkerStyle(2);
	hMinusPowPy8->Draw("histe");
	hMinusPy6->SetMarkerStyle(2);
	hMinusPy6->Draw("histesame");

	TCanvas *c = new TCanvas("c","c",600,600);
	hPlusPowPy8->SetMarkerStyle(2);
	hPlusPy6->SetMarkerStyle(2);
	hMinusPowPy8->SetMarkerStyle(2);
	hMinusPy6->SetMarkerStyle(2);
	hPlusPowPy8->Draw("histe");
	hPlusPy6->Draw("histesame");
	hMinusPowPy8->Draw("histesame");
	hMinusPy6->Draw("histesame");

   TLegend *leg = new TLegend(0.4916107,0.4947552,0.761745,0.791958,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("hPlusPowPy8","W^{+}#rightarrow#mu^{+},PowPy8","l");
   entry->SetLineColor(7);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hMinusPowPy8","W^{-}#rightarrow#mu^{-},PowPy8","l");
   entry->SetLineColor(6);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hPlusPy6","W^{+}#rightarrow#mu^{+},Py6","l");
   entry->SetLineColor(4);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hMinusPy6","W^{-}#rightarrow#mu^{-},Py6","l");
   entry->SetLineColor(2);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();

	TCanvas *cPlusGen = new TCanvas("cplusgen","cplusgen",600,600);
	hGenPlusPowPy8->SetMarkerStyle(2);
	hGenPlusPowPy8->Draw("histe");
	hGenPlusPy6->SetMarkerStyle(2);
	hGenPlusPy6->Draw("histesame");

	TCanvas *cMinusGen = new TCanvas("cminusgen","cminusgen",600,600);
	hGenMinusPowPy8->SetMarkerStyle(2);
	hGenMinusPowPy8->Draw("histe");
	hGenMinusPy6->SetMarkerStyle(2);
	hGenMinusPy6->Draw("histesame");

	TCanvas *cGen = new TCanvas("cGen","cGen",600,600);
	hGenPlusPowPy8->SetMarkerStyle(2);
	hGenPlusPowPy8->Draw("histe");
	hGenPlusPy6->SetMarkerStyle(2);
	hGenPlusPy6->Draw("histesame");
	hGenMinusPowPy8->SetMarkerStyle(2);
	hGenMinusPowPy8->Draw("histesame");
	hGenMinusPy6->SetMarkerStyle(2);
	hGenMinusPy6->Draw("histesame");

	leg->Draw();



}
