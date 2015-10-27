#include "TH1F.h"
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
#include "TF1.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <string>

/*int getCentralityBin(float centrality, bool doAllCentrality = false){
    int index = 0;
    if(doAllCentrality&&(centrality>0.0&&centrality<0.8)) return 0;
    if(centrality>0.0&&centrality<0.05) return index; else ++index;
    if(centrality >0.05&&centrality<0.1) return index; else ++index;
    if(centrality >0.1&&centrality<0.15) return index; else ++index;
    if(centrality >0.15&&centrality<0.2) return index; else ++index;
    if(centrality >0.20&&centrality<0.40) return index; else ++index;
    if(centrality >0.40&&centrality<0.80) return index; else ++index;
    return -1;
}
*/
int getCentralityBin(float centrality, bool doAllCentrality = false){
    int index = 0;
    if(doAllCentrality&&(centrality>0.0&&centrality<0.8)) return 0;
    if(centrality>0.0&&centrality<0.10) return index; else ++index;
    if(centrality >0.1&&centrality<0.80) return index; else ++index;
    return -1;
}


void fillHisto(TString baseString, TString fileNameIn, TH2F* h2DPlus, TH2F* h2DMinus, int index, bool doMc=false){

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

		// Change last argument to false when not binning 
		// in centrality
		int icent = getCentralityBin(centrality, true);
		if(icent<0||icent!=index) continue;
		for(int imu=0; imu<nmu; ++imu){
			if(
			  //(prompt[imu]==24||!doMc) 
    			  //( ( (trig1&&matched1[imu]) || (trig2&&matched2[imu]) || (trig3&&matched3[imu]) ) || doMc)
    			  //( ( mbtrig1 || mbtrig2) || doMc)
    			  ( ( (trig4&&matched4[imu]) || (trig5&&matched5[imu])) || doMc)
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
			   if(charge[imu]>0.) h2DPlus->Fill(eta[imu],phi[imu]); 
	
			   if(charge[imu]<0.) h2DMinus->Fill(eta[imu],phi[imu]);
			}
		}//imu
	}//iev

}

void plotEtaPhiMap(){

      bool doCentrality = false;
      TString baseString = "/usatlas/u/tbales/scratch/";
      TString fileNameInMc = "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.04.12.2013.root";
      TString fileNameInData ="HISingleMuonHardProbesData.04.17.2013.root";
      //TString fileNameInData ="HISingleMuonMinBiasData.07.10.2013.root";
	
      TString fileNameInJ1mu = "MonteCarloFiles/QCD/HISingleMuonJ1_1muPYTHIADataOverlay.04.13.2013.root";
      TString fileNameInJ2mu = "MonteCarloFiles/QCD/HISingleMuonJ2_1muPYTHIADataOverlay.04.13.2013.root";
      TString fileNameInJ3mu = "MonteCarloFiles/QCD/HISingleMuonJ3_1muPYTHIADataOverlay.04.13.2013.root";
      TFile* outFile = new TFile("signalCandidateEtaPhiMap.root","recreate");

      char hNamePlus[50], hNameMinus[50];
      std::vector<double> centralityBins;
      std::vector <float> ncoll;

      centralityBins.push_back(0.00);
      if (doCentrality) {

		centralityBins.push_back(0.05);
		centralityBins.push_back(0.10);
		centralityBins.push_back(0.15);
		centralityBins.push_back(0.20);
		centralityBins.push_back(0.40);
		//centralityBins.push_back(0.60);

		//ncoll
		ncoll.push_back(1683.3); //0-5
		ncoll.push_back(1318.0); //5-10
		ncoll.push_back(1035.4); //10-15
		ncoll.push_back(811.2); //15-20
		ncoll.push_back(440.6); //20-40
		ncoll.push_back(77.8); //40-80

       }

      else  ncoll.push_back(452.0);//0-80
      centralityBins.push_back(0.80);
      const int nCentralityBins = centralityBins.size()-1;

      TH2F* h2DPlus[nCentralityBins], *h2DMinus[nCentralityBins],*h2DMcPlus[nCentralityBins],*h2DMcMinus[nCentralityBins];
      TH2F* h2DJ1muPlus[nCentralityBins],*h2DJ1muMinus[nCentralityBins];
      TH2F* h2DJ2muPlus[nCentralityBins],*h2DJ2muMinus[nCentralityBins];
      TH2F* h2DJ3muPlus[nCentralityBins],*h2DJ3muMinus[nCentralityBins];
      TH2F* h2DJxmuPlus[nCentralityBins],*h2DJxmuMinus[nCentralityBins];

      double etaBins[] = {-2.4,-2.1,-1.85,-1.55,-1.3,-1.05,-0.8,-0.6,-0.35,-0.1,0.1,0.35,0.6,0.8,1.05,1.3,1.55,1.85,2.1,2.4}; 
      int  nEtaBins = sizeof(etaBins)/sizeof(double) - 1;
      const int nPhiBins = 6;
      TH1D* h1DEtaPlus[nCentralityBins][nPhiBins],*h1DEtaMinus[nCentralityBins][nPhiBins],*h1DMcEtaPlus[nCentralityBins][nPhiBins],*h1DMcEtaMinus[nCentralityBins][nPhiBins];
      TH1D* h1DJxmuEtaPlus[nCentralityBins][nPhiBins], *h1DJxmuEtaMinus[nCentralityBins][nPhiBins];
      TH1F* h1DEtaPlusRatio[nCentralityBins][nPhiBins],*h1DEtaMinusRatio[nCentralityBins][nPhiBins];
      TH1F* h1DEtaPlusJxmuRatio[nCentralityBins][nPhiBins],*h1DEtaMinusJxmuRatio[nCentralityBins][nPhiBins];

      for(int icent=0; icent<nCentralityBins; ++icent){

        // data
        sprintf(hNamePlus,"h2DSignalCandidateMuPlusEtaPhiMap_Centrality%i",icent);
        sprintf(hNameMinus,"h2DSignalCandidateMuMinusEtaPhiMap_Centrality%i",icent);
        h2DPlus[icent] = new TH2F(hNamePlus,hNamePlus,nEtaBins,etaBins,nPhiBins,-1*TMath::Pi(),TMath::Pi());
        h2DMinus[icent] = new TH2F(hNameMinus,hNameMinus,nEtaBins,etaBins,nPhiBins,-1*TMath::Pi(),TMath::Pi());
	h2DPlus[icent]->Sumw2();
	h2DMinus[icent]->Sumw2();
	h2DPlus[icent]->GetXaxis()->SetTitle("#eta^{#mu}");
	h2DPlus[icent]->GetYaxis()->SetTitle("#phi^{#mu}");
	h2DMinus[icent]->GetXaxis()->SetTitle("#eta^{#mu}");
	h2DMinus[icent]->GetYaxis()->SetTitle("#phi^{#mu}");

        fillHisto(baseString,fileNameInData,h2DPlus[icent],h2DMinus[icent],icent);

        // MC
        sprintf(hNamePlus,"h2DMcSignalCandidateMuPlusEtaPhiMap_Centrality%i",icent);
        sprintf(hNameMinus,"h2DMcSignalCandidateMuMinusEtaPhiMap_Centrality%i",icent);
        h2DMcPlus[icent] = new TH2F(hNamePlus,hNamePlus,nEtaBins,etaBins,nPhiBins,-1*TMath::Pi(),TMath::Pi());
        h2DMcMinus[icent] = new TH2F(hNameMinus,hNameMinus,nEtaBins,etaBins,nPhiBins,-1*TMath::Pi(),TMath::Pi());
	h2DMcPlus[icent]->Sumw2();
	h2DMcMinus[icent]->Sumw2();
	h2DMcPlus[icent]->GetXaxis()->SetTitle("#eta^{#mu}");
	h2DMcPlus[icent]->GetYaxis()->SetTitle("#phi^{#mu}");
	h2DMcMinus[icent]->GetXaxis()->SetTitle("#eta^{#mu}");
	h2DMcMinus[icent]->GetYaxis()->SetTitle("#phi^{#mu}");

        fillHisto(baseString,fileNameInMc,h2DMcPlus[icent],h2DMcMinus[icent],icent,true);

	// J1mu MC
        sprintf(hNamePlus,"h2DJ1muSignalCandidateMuPlusEtaPhiMap_Centrality%i",icent);
        sprintf(hNameMinus,"h2DJ1SignalCandidateMuMinusEtaPhiMap_Centrality%i",icent);
        h2DJ1muPlus[icent] = new TH2F(hNamePlus,hNamePlus,nEtaBins,etaBins,nPhiBins,-1*TMath::Pi(),TMath::Pi());
        h2DJ1muMinus[icent] = new TH2F(hNameMinus,hNameMinus,nEtaBins,etaBins,nPhiBins,-1*TMath::Pi(),TMath::Pi());
	h2DJ1muPlus[icent]->Sumw2();
	h2DJ1muMinus[icent]->Sumw2();
	fillHisto(baseString,fileNameInJ1mu,h2DJ1muPlus[icent],h2DJ1muMinus[icent],icent,true);

	std::cout << "h2DJ1muPlus: " << h2DJ1muPlus[icent]->Integral() << std::endl;
	std::cout << "h2DJ1muMinus: " << h2DJ1muMinus[icent]->Integral() << std::endl;
	// J2mu MC
        sprintf(hNamePlus,"h2DJ2muSignalCandidateMuPlusEtaPhiMap_Centrality%i",icent);
        sprintf(hNameMinus,"h2DJ2SignalCandidateMuMinusEtaPhiMap_Centrality%i",icent);
        h2DJ2muPlus[icent] = new TH2F(hNamePlus,hNamePlus,nEtaBins,etaBins,nPhiBins,-1*TMath::Pi(),TMath::Pi());
        h2DJ2muMinus[icent] = new TH2F(hNameMinus,hNameMinus,nEtaBins,etaBins,nPhiBins,-1*TMath::Pi(),TMath::Pi());
	h2DJ2muPlus[icent]->Sumw2();
	h2DJ2muMinus[icent]->Sumw2();
	fillHisto(baseString,fileNameInJ2mu,h2DJ2muPlus[icent],h2DJ2muMinus[icent],icent,true);
	std::cout << "h2DJ2muPlus: " << h2DJ2muPlus[icent]->Integral() << std::endl;
	std::cout << "h2DJ2muMinus: " << h2DJ2muMinus[icent]->Integral() << std::endl;
	// J3mu MC
        sprintf(hNamePlus,"h2DJ3muSignalCandidateMuPlusEtaPhiMap_Centrality%i",icent);
        sprintf(hNameMinus,"h2DJ3SignalCandidateMuMinusEtaPhiMap_Centrality%i",icent);
        h2DJ3muPlus[icent] = new TH2F(hNamePlus,hNamePlus,nEtaBins,etaBins,nPhiBins,-1*TMath::Pi(),TMath::Pi());
        h2DJ3muMinus[icent] = new TH2F(hNameMinus,hNameMinus,nEtaBins,etaBins,nPhiBins,-1*TMath::Pi(),TMath::Pi());
	h2DJ3muPlus[icent]->Sumw2();
	h2DJ3muMinus[icent]->Sumw2();
	fillHisto(baseString,fileNameInJ3mu,h2DJ3muPlus[icent],h2DJ3muMinus[icent],icent,true);
	std::cout << "h2DJ3muPlus: " << h2DJ3muPlus[icent]->Integral() << std::endl;
	std::cout << "h2DJ3muMinus: " << h2DJ3muMinus[icent]->Integral() << std::endl;

	//weight the Jx samples according to cross-sections; takes into account prob of find muon+jet (AMI)
	double wtJ1 = 6.6523e-03*1.8770e-4/64.0e-3; double wtJ2 = 1.4941e-02*8.2788e-6/64.0e-3; double wtJ3 = 2.4284e-02*2.9426e-7/64.0e-3;
	double arrCentWidth = centralityBins[icent+1]-centralityBins[icent];
	double evData = 1.03e9; //number of events sampled
	double scaleFactor = arrCentWidth*ncoll[icent]*evData;

        sprintf(hNamePlus,"h2DJxmuSignalCandidateMuPlusEtaPhiMap_Centrality%i",icent);
        sprintf(hNameMinus,"h2DJxSignalCandidateMuMinusEtaPhiMap_Centrality%i",icent);
        h2DJxmuPlus[icent] = new TH2F(hNamePlus,hNamePlus,nEtaBins,etaBins,nPhiBins,-1*TMath::Pi(),TMath::Pi());
        h2DJxmuMinus[icent] = new TH2F(hNameMinus,hNameMinus,nEtaBins,etaBins,nPhiBins,-1*TMath::Pi(),TMath::Pi());

	h2DJxmuPlus[icent]->Add(h2DJ1muPlus[icent],h2DJ2muPlus[icent],wtJ1,wtJ2);
	h2DJxmuPlus[icent]->Add(h2DJ3muPlus[icent],wtJ3);

	h2DJxmuMinus[icent]->Add(h2DJ1muMinus[icent],h2DJ2muMinus[icent],wtJ1,wtJ2);
	h2DJxmuMinus[icent]->Add(h2DJ3muMinus[icent],wtJ3);
	h2DJxmuPlus[icent]->Sumw2();
	h2DJxmuMinus[icent]->Sumw2();
	h2DJxmuPlus[icent]->GetXaxis()->SetTitle("#eta^{#mu}");
	h2DJxmuPlus[icent]->GetYaxis()->SetTitle("#phi^{#mu}");
	h2DJxmuMinus[icent]->GetXaxis()->SetTitle("#eta^{#mu}");
	h2DJxmuMinus[icent]->GetYaxis()->SetTitle("#phi^{#mu}");

	//1D projections of phi slices along eta axis
	//for(int iphi = 0; iphi<nPhiBins; ++iphi){
	//}//iphi
	
      }//icent

      
      // Save histos
      outFile->cd();
	
	for(int icent=0; icent<nCentralityBins; ++icent){
		h2DPlus[icent]->Scale(1.0/h2DPlus[icent]->Integral());
		h2DMinus[icent]->Scale(1.0/h2DMinus[icent]->Integral());
		h2DPlus[icent]->SetMaximum(0.018);
		h2DPlus[icent]->SetMinimum(0.0);
		std::cout << "Writing mu+ histo to file..." << std::endl;
		h2DPlus[icent]->Write();
		std::cout << "Done." << std::endl;
		h2DMinus[icent]->SetMaximum(0.018);
		h2DMinus[icent]->SetMinimum(0.0);
		std::cout << "Writing mu- histo to file..." << std::endl;
		h2DMinus[icent]->Write();
		std::cout << "Done." << std::endl;
			
		h2DMcPlus[icent]->Scale(1.0/h2DMcPlus[icent]->Integral());
		h2DMcMinus[icent]->Scale(1.0/h2DMcMinus[icent]->Integral());
		h2DMcPlus[icent]->SetMaximum(0.018);
		h2DMcPlus[icent]->SetMinimum(0.0);
		std::cout << "Writing MC mu+ histo to file..." << std::endl;
		h2DMcPlus[icent]->Write();
		std::cout << "Done." << std::endl;
		h2DMcMinus[icent]->SetMaximum(0.018);
		h2DMcMinus[icent]->SetMinimum(0.0);
		std::cout << "Writing MC mu- histo to file..." << std::endl;
		h2DMcMinus[icent]->Write();
		std::cout << "Done." << std::endl;
			
		h2DJxmuPlus[icent]->Scale(1.0/h2DJxmuPlus[icent]->Integral());
		h2DJxmuMinus[icent]->Scale(1.0/h2DJxmuMinus[icent]->Integral());
		h2DJxmuPlus[icent]->SetMaximum(0.018);
		h2DJxmuPlus[icent]->SetMinimum(0.0);
		std::cout << "Writing Jx mu+ histo to file..." << std::endl;
		h2DJxmuPlus[icent]->Write();
		std::cout << "Done." << std::endl;
		h2DJxmuMinus[icent]->SetMaximum(0.018);
		h2DJxmuMinus[icent]->SetMinimum(0.0);
		std::cout << "Writing Jx mu- histo to file..." << std::endl;
		h2DJxmuMinus[icent]->Write();
		std::cout << "Done." << std::endl;
	
		// Project eta(x) from phi (y) channels
		for(int iphi = 0; iphi<nPhiBins; ++iphi){

       		   h1DEtaPlus[icent][iphi] = new TH1D();
		   h1DEtaMinus[icent][iphi] = new TH1D();
		   h1DMcEtaPlus[icent][iphi] = new TH1D();
		   h1DMcEtaMinus[icent][iphi] = new TH1D();
		   h1DJxmuEtaPlus[icent][iphi] = new TH1D();
		   h1DJxmuEtaMinus[icent][iphi] = new TH1D();

		   h1DEtaPlus[icent][iphi]->Sumw2();
		   h1DEtaMinus[icent][iphi]->Sumw2();
		
		   h1DEtaPlus[icent][iphi] = h2DPlus[icent]->ProjectionX(Form("muPlusPhiBin%dCentBin%i",iphi+1,icent),iphi+1,iphi+2);
		   h1DEtaPlus[icent][iphi]->GetXaxis()->SetTitle("#eta^{#mu}");
		   h1DEtaPlus[icent][iphi]->Write();
		   h1DEtaMinus[icent][iphi] = h2DMinus[icent]->ProjectionX(Form("muMinusPhiBin%dCentBin%i",iphi+1,icent),iphi+1,iphi+2);
		   h1DEtaMinus[icent][iphi]->GetXaxis()->SetTitle("#eta^{#mu}");
		   h1DEtaMinus[icent][iphi]->Write();

		   h1DMcEtaPlus[icent][iphi]->Sumw2();
		   h1DMcEtaMinus[icent][iphi]->Sumw2();

		   h1DJxmuEtaPlus[icent][iphi]->Sumw2();
		   h1DJxmuEtaMinus[icent][iphi]->Sumw2();
		
		   h1DMcEtaPlus[icent][iphi] = h2DMcPlus[icent]->ProjectionX(Form("muPlusMcPhi%iCentBin%i",iphi,icent),iphi+1,iphi+2);
		   h1DMcEtaPlus[icent][iphi]->GetXaxis()->SetTitle("#eta^{#mu}");
		   h1DMcEtaPlus[icent][iphi]->Write();
		   h1DMcEtaMinus[icent][iphi] = h2DMcMinus[icent]->ProjectionX(Form("muMinusMcPhi%iCentBin%i",iphi,icent),iphi+1,iphi+2);
		   h1DMcEtaMinus[icent][iphi]->GetXaxis()->SetTitle("#eta^{#mu}");
		   h1DMcEtaMinus[icent][iphi]->Write();

		   h1DJxmuEtaPlus[icent][iphi] = h2DJxmuPlus[icent]->ProjectionX(Form("muPlusJxmuPhi%iCentBin%i",iphi,icent),iphi+1,iphi+2);
		   h1DJxmuEtaPlus[icent][iphi]->GetXaxis()->SetTitle("#eta^{#mu}");
		   h1DJxmuEtaPlus[icent][iphi]->Write();
		   h1DJxmuEtaMinus[icent][iphi] = h2DJxmuMinus[icent]->ProjectionX(Form("muMinusJxmuPhi%iCentBin%i",iphi,icent),iphi+1,iphi+2);
		   h1DJxmuEtaMinus[icent][iphi]->GetXaxis()->SetTitle("#eta^{#mu}");
		   h1DJxmuEtaMinus[icent][iphi]->Write();
		   
		   //Data:MC ratio
                   sprintf(hNamePlus,"h1DSignalCandidateMuPlusDataMcRatio_Phi%i_Centrality%i",iphi,icent);
                   sprintf(hNameMinus,"h1DSignalCandidateMuMinusDataMcRatio_Phi%i_Centrality%i",iphi,icent);
		   h1DEtaPlusRatio[icent][iphi] = new TH1F(hNamePlus,hNamePlus,nEtaBins,etaBins);
		   h1DEtaMinusRatio[icent][iphi] = new TH1F(hNameMinus,hNameMinus,nEtaBins,etaBins);
		   h1DEtaPlusRatio[icent][iphi]->Sumw2();
		   h1DEtaMinusRatio[icent][iphi]->Sumw2();

		   h1DEtaPlusRatio[icent][iphi]->Divide(h1DEtaPlus[icent][iphi],h1DMcEtaPlus[icent][iphi]); 
		   h1DEtaMinusRatio[icent][iphi]->Divide(h1DEtaMinus[icent][iphi],h1DMcEtaMinus[icent][iphi]); 

		   h1DEtaPlusRatio[icent][iphi]->Write();
		   h1DEtaMinusRatio[icent][iphi]->Write();

                   sprintf(hNamePlus,"h1DSignalCandidateMuPlusDataJxmuRatio_Phi%i_Centrality%i",iphi,icent);
                   sprintf(hNameMinus,"h1DSignalCandidateMuMinusDataJxmuRatio_Phi%i_Centrality%i",iphi,icent);
		   h1DEtaPlusJxmuRatio[icent][iphi] = new TH1F(hNamePlus,hNamePlus,nEtaBins,etaBins);
		   h1DEtaMinusJxmuRatio[icent][iphi] = new TH1F(hNameMinus,hNameMinus,nEtaBins,etaBins);
		   h1DEtaPlusJxmuRatio[icent][iphi]->Sumw2();
		   h1DEtaMinusJxmuRatio[icent][iphi]->Sumw2();

		   h1DEtaPlusJxmuRatio[icent][iphi]->Divide(h1DEtaPlus[icent][iphi],h1DJxmuEtaPlus[icent][iphi]); 
		   h1DEtaMinusJxmuRatio[icent][iphi]->Divide(h1DEtaMinus[icent][iphi],h1DJxmuEtaMinus[icent][iphi]); 

		   h1DEtaPlusJxmuRatio[icent][iphi]->Write();
		   h1DEtaMinusJxmuRatio[icent][iphi]->Write();

		  TF1 *f1 = new TF1("f1","[0]",-2.5,2.5);
		  f1->FixParameter(0,1.0);
		  h1DEtaPlusRatio[icent][iphi]->Fit("f1");
		  std::cout << "Chi2 for mu+ in centrality bin " << icent << " phi bin " << iphi << " = " << f1->GetChisquare() << "/" << f1->GetNDF() << std::endl;
		  h1DEtaMinusRatio[icent][iphi]->Fit("f1");
		  std::cout << "Chi2 for mu- centrality bin " << icent << " phi bin " << iphi << " = " << f1->GetChisquare() << "/" << f1->GetNDF() << std::endl;

		  delete h1DEtaPlus[icent][iphi];
		  delete h1DEtaMinus[icent][iphi];
		  delete h1DMcEtaPlus[icent][iphi];
		  delete h1DMcEtaMinus[icent][iphi];
		  delete h1DEtaPlusRatio[icent][iphi];
		  delete h1DEtaMinusRatio[icent][iphi];
		  delete h1DEtaPlusJxmuRatio[icent][iphi];
		  delete h1DEtaMinusJxmuRatio[icent][iphi];
		}//iphi

		delete h2DPlus[icent];
		delete h2DMinus[icent];
		delete h2DMcPlus[icent];
		delete h2DMcMinus[icent];
		delete h2DJ1muPlus[icent];
		delete h2DJ1muMinus[icent];
		delete h2DJ2muPlus[icent];
		delete h2DJ2muMinus[icent];
		delete h2DJ3muPlus[icent];
		delete h2DJ3muMinus[icent];
		delete h2DJxmuPlus[icent];
		delete h2DJxmuMinus[icent];
	}//icent	
}
