#include "TH1.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"

#include <iostream>
#include <string>

void weightFcal(TF1* f, TH1F* h){
    
    //TF1* f = wts->GetFunction("f1");
    for(int ibin=1; ibin<h->GetNbinsX(); ++ibin){
       double weight = f->Eval(h->GetBinCenter(ibin)); 
       h->SetBinContent(ibin,h->GetBinContent(ibin)*weight);
    }//bin
}

void fillHisto(TH1F* h, TString s){

	TChain* tree = new TChain("tree","tree");
	tree->Add(s);

	std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;

    float fcalNt;
    float ptNt[50];
    int nmu;
    tree->SetBranchAddress("fcal",&fcalNt);
    tree->SetBranchAddress("pt",&ptNt);
    tree->SetBranchAddress("mu_muid_n", &nmu);
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("fcal",1);
    tree->SetBranchStatus("pt",1);
    tree->SetBranchStatus("mu_muid_n",1);

	for(int iev=0; iev<tree->GetEntries(); ++iev){

		tree->LoadTree(iev);
		tree->GetEntry(iev);
        //if(iev==10000) break; //hack
    
        bool flag = false;
		if(iev%10000==0) std::cout << "Event: " << iev << std::endl;
        for(int imu=0; imu<nmu; ++imu){
            if(ptNt[imu]>10.0){
                flag=true; break;
            }
        }//imu

        if(flag==true) {
            // only fill once per PbPb event
            h->Fill(fcalNt);
            continue;
        }
    }//iev
}

void FCalWeights(){
    
    TFile* outFile = new TFile("fcalHistos.root","recreate");
    TString baseString = "/usatlas/u/tbales/scratch/";

    TString fileNameIn = "HISingleMuonHardProbesData.04.17.2013";

    //PowPy8 samples
    //pp
    TString fileNameIn_ppPlus = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_pp.10.03.2013";
    TString fileNameIn_ppMinus = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_pp.10.03.2013";
    //np
    TString fileNameIn_npPlus = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_np.10.03.2013";
    TString fileNameIn_npMinus = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_np.10.03.2013";
    //pn
    TString fileNameIn_pnPlus = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_pn.10.03.2013";
    TString fileNameIn_pnMinus = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_pn.10.03.2013";
    //nn
    TString fileNameIn_nnPlus = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_nn.10.03.2013";
    TString fileNameIn_nnMinus = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_nn.10.03.2013";

    TH1F *hData, *hpp_pos, *hpn_pos, *hnp_pos, *hnn_pos;
    TH1F *hpp_neg, *hpn_neg, *hnp_neg, *hnn_neg;
    TH1F *hpp,*hpn,*hnn,*hMc; 
    TH1F* hRatio;
    const int nBins = 100;
    const float binLow = 0.0; const float binUpp = 5.0;
    hData=new TH1F("hData","hData",nBins,binLow,binUpp);   hData->Sumw2();
    hpp_pos=new TH1F("hpp_pos","hpp_pos",nBins,binLow,binUpp); hpp_pos->Sumw2();
    hpn_pos=new TH1F("hpn_pos","hpn_pos",nBins,binLow,binUpp); hpn_pos->Sumw2();
    hnp_pos=new TH1F("hnp_pos","hnp_pos",nBins,binLow,binUpp); hnp_pos->Sumw2();
    hnn_pos=new TH1F("hnn_pos","hnn_pos",nBins,binLow,binUpp); hnn_pos->Sumw2();
    hpp_neg=new TH1F("hpp_neg","hpp_neg",nBins,binLow,binUpp); hpp_neg->Sumw2();
    hpn_neg=new TH1F("hpn_neg","hpn_neg",nBins,binLow,binUpp); hpn_neg->Sumw2();
    hnp_neg=new TH1F("hnp_neg","hnp_neg",nBins,binLow,binUpp); hnp_neg->Sumw2();
    hnn_neg=new TH1F("hnn_neg","hnn_neg",nBins,binLow,binUpp); hnn_neg->Sumw2();
    hpp=new TH1F("hpp","hpp",nBins,binLow,binUpp); hpp->Sumw2();
    hpn=new TH1F("hpn","hpn",nBins,binLow,binUpp); hpn->Sumw2();
    hnn=new TH1F("hnn","hnn",nBins,binLow,binUpp); hnn->Sumw2();
    hMc=new TH1F("hMc","hMc",nBins,binLow,binUpp); hMc->Sumw2();
    hRatio=new TH1F("hRatio","hRatio",nBins,binLow,binUpp); hRatio->Sumw2();
                                                              
    fillHisto(hData,baseString+fileNameIn+".root");
    fillHisto(hpp_pos,baseString+fileNameIn_ppPlus+".root");
    fillHisto(hpp_neg,baseString+fileNameIn_ppMinus+".root");
    fillHisto(hpn_pos,baseString+fileNameIn_pnPlus+".root");
    fillHisto(hpn_neg,baseString+fileNameIn_pnMinus+".root");
    fillHisto(hnp_pos,baseString+fileNameIn_npPlus+".root");
    fillHisto(hnp_neg,baseString+fileNameIn_npMinus+".root");
    fillHisto(hnn_pos,baseString+fileNameIn_nnPlus+".root");
    fillHisto(hnn_neg,baseString+fileNameIn_nnMinus+".root");

    // Add pos and neg
    hpp->Add(hpp_pos,hpp_neg);
    hpn->Add(hpn_pos,hpn_neg);
    hpn->Add(hnp_pos); hpn->Add(hnp_neg);
    hnn->Add(hnn_pos,hnn_neg);
    // Weight
    double ppWt = 0.155;
    double pnWt = 0.478;
    double nnWt = 0.367;
    hMc->Add(hpp,hpn,ppWt,pnWt); 
    hMc->Add(hnn,nnWt);
    hRatio->Divide(hData,hMc);
    TF1 *f1 = new TF1("f1", "pol2", 0.0, 5.0);
    f1->SetParameters(0,0.0863814);
    f1->SetParameters(1,8.82331);
    f1->SetParameters(2,-1.41379);

    TH1F* hMCc = new TH1F("hMCc","hMCc",nBins,binLow,binUpp); hMCc->Sumw2();
    hMCc = (TH1F*)hMc->Clone("hMCc");
    hRatio->Fit("f1","R");
    weightFcal(hRatio->GetFunction("f1"),hMCc);

    outFile->cd();
    hMc->Scale(hData->Integral()/hMc->Integral());
    hMc->Write();
    hData->Write();
    hRatio->Write("weights");
    hMCc->Scale(hData->Integral()/hMCc->Integral());
    hMCc->Write("fcalWtd");

}
