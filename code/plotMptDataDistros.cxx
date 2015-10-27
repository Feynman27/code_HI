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


int getCentralityBin(float centrality, bool doAllCentrality = false){
    int index = 0;
    if(doAllCentrality&&(centrality>0.0&&centrality<0.8)) return 0;
    if(centrality>0.0&&centrality<0.10) return index; else ++index;
    if(centrality >0.1&&centrality<0.80) return index; else ++index;
    return -1;
}

void fillHisto(TString baseString, TString fileNameIn, TH1F* hPlus, TH1F* hMinus, int index, 
                TString sMpt = "nu_pt", TString sMptPhi = "nu_phi", bool doSystematic = false){

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
      int val[50], ZDY[50], matched1[50], matched2[50], matched3[50],matched4[50],matched5[50];
      int nmu,trig1,trig2,trig3,trig4,trig5, mbtrig1, mbtrig2;


      // Set branch address
      tree->SetBranchAddress("pt",&pt);
      tree->SetBranchAddress("eLoss", &eLoss);
      tree->SetBranchAddress("ptcone20ID3", &ptcone);
      tree->SetBranchAddress("scat", &scat);
      tree->SetBranchAddress("comp", &comp);
      tree->SetBranchAddress("pt", &pt);
      tree->SetBranchAddress("mt", &mt);
      tree->SetBranchAddress("eta", &eta);
      tree->SetBranchAddress("phi", &phi);
      tree->SetBranchAddress(sMptPhi, &mptPhi);
      tree->SetBranchAddress("charge", &charge);
      tree->SetBranchAddress("val", &val); 
      tree->SetBranchAddress("ZDY", &ZDY); 
      tree->SetBranchAddress("centrality", &centrality);
      tree->SetBranchAddress(sMpt, &mpt);
      tree->SetBranchAddress("mu_muid_n", &nmu);
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

      // --- Set branch status ---
      tree->SetBranchStatus("*",0) ;
      tree->SetBranchStatus("mu_muid_n", 1);
      tree->SetBranchStatus("eLoss", 1);
      tree->SetBranchStatus("ptcone20ID3", 1);
      tree->SetBranchStatus("scat", 1);
      tree->SetBranchStatus("comp", 1);
      tree->SetBranchStatus("pt", 1);
      tree->SetBranchStatus("mt", 1);
      tree->SetBranchStatus("eta", 1);
      tree->SetBranchStatus("phi", 1);
      tree->SetBranchStatus(sMptPhi, 1);
      tree->SetBranchStatus("charge", 1);
      tree->SetBranchStatus("val", 1); 
      tree->SetBranchStatus("ZDY", 1); 
      tree->SetBranchStatus("centrality", 1);
      tree->SetBranchStatus(sMpt, 1);
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
 
	for(int iev=0; iev<tree->GetEntries(); ++iev){

	    if(iev%10000==0) std::cout << "Event: " << iev << std::endl;
	    tree->LoadTree(iev);
	    tree->GetEntry(iev);

		// Change last argument to false when not binning 
		// in centrality
		int icent = getCentralityBin(centrality, true);
		if(icent<0||icent!=index) continue;
		for(int imu=0; imu<nmu; ++imu){

            float dPhi = phi[imu]-mptPhi;
            if(dPhi<-1*TMath::Pi()) dPhi += TMath::TwoPi(); if(dPhi>TMath::Pi()) dPhi -= TMath::TwoPi(); //fold btwn [-pi,pi]
            double mtSyst = fabs(mpt)<9000. ? TMath::Sqrt(2.0*mpt*pt[imu]*(1.0-TMath::Cos(dPhi))):-9999.;
            double mtTemp = -9999;
            if(doSystematic) mtTemp = mtSyst;
            else  mtTemp = mt[imu];

			if(
			  //(prompt[imu]==24||!doMc) 
    		   ( ( (trig1&&matched1[imu]) || (trig2&&matched2[imu]) || (trig3&&matched3[imu]) ) )
    			  //( ( mbtrig1 || mbtrig2) || doMc)
    			  //( ( (trig4&&matched4[imu]) || (trig5&&matched5[imu])) || doMc)
			  &&val[imu]>11
			  //val[imu]>0
			  && fabs(scat[imu])<4.0
			  && fabs(eLoss[imu])<0.5
			  && ptcone[imu]/pt[imu] < 0.1
			  && ZDY[imu]==0
			  //&& pt[imu]>10.0
			  && pt[imu]>25.0
			  && mpt<9000.
			  && mtTemp>40.0
			){
			   if(charge[imu]>0.) hPlus->Fill(mpt); 
	
			   if(charge[imu]<0.) hMinus->Fill(mpt);
			}
		}//imu
	}//iev

}

void plotMptDataDistros(){

      bool doCentrality = true;
      TString baseString = "/usatlas/u/tbales/scratch/";
      TString fileNameInData ="HISingleMuonHardProbesData.07.13.2013.root";

      std::vector<double> centralityBins;
      std::vector <float> ncoll;


      TString arrCentLabel[] = {"0-5%","5-10%","10-15%","15-20%","20-40%","40-80%"};
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

      TH1F* h3GeVPlus[nCentralityBins], *h2GeVPlus[nCentralityBins], *h4GeVPlus[nCentralityBins];
      TH1F* h3GeVMinus[nCentralityBins], *h2GeVMinus[nCentralityBins], *h4GeVMinus[nCentralityBins];
      //pp
      TH1F* hpp3GeVPlus[nCentralityBins], *hpp2GeVPlus[nCentralityBins], *hpp4GeVPlus[nCentralityBins];
      TH1F* hpp3GeVMinus[nCentralityBins], *hpp2GeVMinus[nCentralityBins], *hpp4GeVMinus[nCentralityBins];
      //np
      TH1F* hnp3GeVPlus[nCentralityBins], *hnp2GeVPlus[nCentralityBins], *hnp4GeVPlus[nCentralityBins];
      TH1F* hnp3GeVMinus[nCentralityBins], *hnp2GeVMinus[nCentralityBins], *hnp4GeVMinus[nCentralityBins];
      //pn
      TH1F* hpn3GeVPlus[nCentralityBins], *hpn2GeVPlus[nCentralityBins], *hpn4GeVPlus[nCentralityBins];
      TH1F* hpn3GeVMinus[nCentralityBins], *hpn2GeVMinus[nCentralityBins], *hpn4GeVMinus[nCentralityBins];
      //nn
      TH1F* hnn3GeVPlus[nCentralityBins], *hnn2GeVPlus[nCentralityBins], *hnn4GeVPlus[nCentralityBins];
      TH1F* hnn3GeVMinus[nCentralityBins], *hnn2GeVMinus[nCentralityBins], *hnn4GeVMinus[nCentralityBins];

      double wtpp = 0.15;
      double wtnp = 0.48;
      double wtnn = 0.37;


      char hNamePlus[50], hNameMinus[50];
      for(int icent=0; icent<nCentralityBins; ++icent){

        sprintf(hNamePlus,"hMpt3GeVPlus_Centrality%i",icent);
        sprintf(hNameMinus,"hMpt3GeVMinus_Centrality%i",icent);
        h3GeVPlus[icent] = new TH1F(hNamePlus,hNamePlus,25,0.0,100.0);
        h3GeVMinus[icent] = new TH1F(hNameMinus,hNameMinus,25,0.0,100.0);
	    h3GeVPlus[icent]->Sumw2();
	    h3GeVMinus[icent]->Sumw2();
        h3GeVPlus[icent]->GetXaxis()->SetTitle("#slash{p_{T}}");
        h3GeVMinus[icent]->GetXaxis()->SetTitle("#slash{p_{T}}");
        fillHisto(baseString,fileNameInData,h3GeVPlus[icent],h3GeVMinus[icent],icent,"nu_pt","nu_phi");

        sprintf(hNamePlus,"hMpt2GeVPlus_Centrality%i",icent);
        sprintf(hNameMinus,"hMpt2GeVMinus_Centrality%i",icent);
        h2GeVPlus[icent] = new TH1F(hNamePlus,hNamePlus,25,0.0,100.0);
        h2GeVMinus[icent] = new TH1F(hNameMinus,hNameMinus,25,0.0,100.0);
	    h2GeVPlus[icent]->Sumw2();
	    h2GeVMinus[icent]->Sumw2();
        fillHisto(baseString,fileNameInData,h2GeVPlus[icent],h2GeVMinus[icent],icent,"nu_pt2000Nominal","nu_phi2000Nominal",true);

        sprintf(hNamePlus,"hMpt4GeVPlus_Centrality%i",icent);
        sprintf(hNameMinus,"hMpt4GeVMinus_Centrality%i",icent);
        h4GeVPlus[icent] = new TH1F(hNamePlus,hNamePlus,25,0.0,100.0);
        h4GeVMinus[icent] = new TH1F(hNameMinus,hNameMinus,25,0.0,100.0);
	    h4GeVPlus[icent]->Sumw2();
	    h4GeVMinus[icent]->Sumw2();
        fillHisto(baseString,fileNameInData,h4GeVPlus[icent],h4GeVMinus[icent],icent,"nu_pt4000Nominal","nu_phi4000Nominal",true);

        //MC
/*        sprintf(hNamePlus,"hppMpt4GeVPlus_Centrality%i",icent);
        sprintf(hNameMinus,"hppMpt4GeVMinus_Centrality%i",icent);
        hpph4GeVPlus[icent] = new TH1F(hNamePlus,hNamePlus,25,0.0,100.0);
        hpph4GeVMinus[icent] = new TH1F(hNameMinus,hNameMinus,25,0.0,100.0);
	    hpph4GeVPlus[icent]->Sumw2();
	    hpph4GeVMinus[icent]->Sumw2();
        fillHisto(baseString,fileNameInData,hpph4GeVPlus[icent],hpph4GeVMinus[icent],icent,"nu_pt4000Nominal","nu_phi4000Nominal",true);
*/
        // Now that histos are filled for each track threshold, plot on same canvas
        
        TString cNamePlus = "cPlusCent"; cNamePlus+=icent;
        TCanvas* cPlus = new TCanvas(cNamePlus,cNamePlus,600,600);
        h3GeVPlus[icent]->Draw("pe");
        h4GeVPlus[icent]->Draw("pesame");
        h2GeVPlus[icent]->Draw("pesame");

        h4GeVPlus[icent]->SetMarkerColor(kBlue);
        h4GeVPlus[icent]->SetLineColor(kBlue);
        h2GeVPlus[icent]->SetMarkerColor(kRed);
        h2GeVPlus[icent]->SetLineColor(kRed);

        TLegend* leg = new TLegend(0.58,0.22,0.85,0.46);
        leg->SetTextFont(gStyle->GetTextFont());
        leg->SetTextSize(gStyle->GetTextSize());
        leg->SetBorderSize(0);
        leg->SetFillColor(0);
        leg->AddEntry(h3GeVPlus[icent],"p_{T}^{trk}>3GeV","lpe");
        leg->AddEntry(h4GeVPlus[icent],"p_{T}^{trk}>4 GeV","lpe");
        leg->AddEntry(h2GeVPlus[icent],"p_{T}^{trk}>2 GeV","lpe");
        leg->Draw();

        TLatex* tex = new TLatex(0.2,0.2,"#mu^{+}"); 
        tex->SetNDC(); 
        tex->SetTextSize(0.1);
        tex->Draw();
        tex = new TLatex(0.2,0.4,arrCentLabel[icent]); 
        tex->SetNDC(); 
        tex->Draw();

        //cPlus->Print(cNamePlus+"_10_29_2013.pdf");
        
        TString cNameMinus = "cMinusCent"; cNameMinus+=icent;
        TCanvas* cMinus = new TCanvas(cNameMinus,cNameMinus,600,600);
        h3GeVMinus[icent]->Draw("pe");
        h4GeVMinus[icent]->Draw("pesame");
        h2GeVMinus[icent]->Draw("pesame");

        h4GeVMinus[icent]->SetMarkerColor(kBlue);
        h4GeVMinus[icent]->SetLineColor(kBlue);
        h2GeVMinus[icent]->SetMarkerColor(kRed);
        h2GeVMinus[icent]->SetLineColor(kRed);

        leg = new TLegend(0.58,0.22,0.85,0.46);
        leg->SetTextFont(gStyle->GetTextFont());
        leg->SetTextSize(gStyle->GetTextSize());
        leg->SetBorderSize(0);
        leg->SetFillColor(0);
        leg->AddEntry(h3GeVMinus[icent],"p_{T}^{trk}>3GeV","lpe");
        leg->AddEntry(h4GeVMinus[icent],"p_{T}^{trk}>4 GeV","lpe");
        leg->AddEntry(h2GeVMinus[icent],"p_{T}^{trk}>2 GeV","lpe");
        leg->Draw();

        tex = new TLatex(0.2,0.2,"#mu^{-}"); 
        tex->SetNDC(); 
        tex->SetTextSize(0.1);
        tex->Draw();
        tex = new TLatex(0.2,0.4,arrCentLabel[icent]); 
        tex->SetNDC(); 
        tex->Draw();


      }//icent
}
