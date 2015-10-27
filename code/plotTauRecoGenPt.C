#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooBinning.h"

#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPad.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"

#include "AtlasUtils.C"
#include "AtlasStyle.C"

#include <iostream>
#include <cmath> 

using namespace RooFit;

RooDataSet* fillTauMuonDataSet(const TString& pathName, const TString& fileName, RooArgSet& muonArgSet, TH2F& hpull){
    
    ///Create RooDataSet with N-dimensions
    RooDataSet* set = new RooDataSet(fileName,fileName,muonArgSet);

    ///Load the tree
    TChain* tree = new TChain("tree","tree");
    tree->Add(pathName+fileName);
    int nMu;
    float mc_mt[50],mc_mptPt,centrality,mc_pt[50],mc_ptGen[50],mc_eta[50];
    float pull;
    int nEvents = tree->GetEntries();
    std::cout << "Number of events: " << nEvents << std::endl;
    tree->SetBranchAddress("nTauMu",&nMu);
    tree->SetBranchAddress("mc_mt",&mc_mt);
    tree->SetBranchAddress("mc_mptPt",&mc_mptPt);
    tree->SetBranchAddress("centrality",&centrality);
    tree->SetBranchAddress("mc_pt",&mc_pt);
    tree->SetBranchAddress("mc_ptGen",&mc_ptGen);
    tree->SetBranchAddress("mc_eta",&mc_eta);

    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("nTauMu",1);
    tree->SetBranchStatus("mc_mt",1);
    tree->SetBranchStatus("mc_mptPt",1);
    tree->SetBranchStatus("centrality",1);
    tree->SetBranchStatus("mc_pt",1);
    tree->SetBranchStatus("mc_ptGen",1);
    tree->SetBranchStatus("mc_eta",1);

    for(int iev=0; iev<nEvents; ++iev){

        tree->LoadTree(iev);
        tree->GetEntry(iev);

        muonArgSet.setRealValue("missPt",mc_mptPt);
        muonArgSet.setRealValue("centrality",centrality);
        for(int imu=0; imu<nMu; ++imu){
            
            muonArgSet.setRealValue("muonPt",mc_pt[imu]);
            muonArgSet.setRealValue("muonGenPt",mc_ptGen[imu]);
            muonArgSet.setRealValue("muonEta",mc_eta[imu]);
            muonArgSet.setRealValue("muonMt",mc_mt[imu]);
            pull = (mc_pt[imu]-mc_ptGen[imu])/mc_ptGen[imu];
            muonArgSet.setRealValue("ptPull",pull);
            hpull.Fill(mc_pt[imu],pull);
            ///Add muon to the set
            set->add(muonArgSet);  
        }//imu
    }//iev
   return set;
}

void plotTauRecoGenPt(){
    
    ///W-->tau-->mu Generated sample
	TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
    TString fileNameMCTauIn = "MonteCarloFiles/Wtaumu/HIWtaumuNtuple.07.30.2013";
    //TString fileNameMCTauIn = "MonteCarloFiles/Wtaumu/HIWtaumuNtuple.07.25.2013";
	SetAtlasStyle();
    // --- declare cut variables --- //
	RooRealVar  muonPt("muonPt","p_{T}",0.0,350.0,"GeV");
    RooRealVar  missPt("missPt","p_{T}^{miss}",0.0,250.0,"GeV");
    RooRealVar  muonMt("muonMt","m_{T}",0.0,350.0,"GeV");
	RooRealVar  muonGenPt("muonGenPt","p_{T}^{Gen}",0.0,350.0,"GeV");
    RooRealVar  ptPull("ptPull","ptPull",0.0,+350.0);
  	RooRealVar  centrality("centrality","centrality",0.,1.0);
  	RooRealVar  muonEta("muonEta","muonEta",-3.0,+3.0);
    
    ///Cuts for generator-level tau sample
    TString sCutsTau = "";
    sCutsTau = "abs(muonEta)>0.1&&abs(muonEta)<2.4&&centrality>0.&&centrality<0.8&&muonPt>17.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0";
    std::cout << " Tau cuts: "<< sCutsTau << std::endl;
    RooArgList muonArgList(muonEta,centrality,muonPt,missPt,muonMt);
    RooFormulaVar cuts("cuts","cuts",sCutsTau,muonArgList);
    RooArgSet muonArgSet(muonPt,missPt,muonMt,muonGenPt,centrality,muonEta,ptPull);

    TH2F* hptPull = new TH2F("hptPull","hptPull",200,0.0,200.0,100,-0.1,0.1);
    RooDataSet* dataSet = fillTauMuonDataSet(baseString,fileNameMCTauIn+".root",muonArgSet,*hptPull);dataSet->Print();
    RooBinning b(0.0,200);
    b.addUniform(100,0.0,100.0);
    b.addUniform(25,100.0,150.0);
    b.addUniform(10,150.0,200.0);
    RooBinning b2(100,-0.1,0.1);

    ///create histograms
    TH1F* hSmearedPt = (TH1F*)dataSet->createHistogram("hSmearedPt", muonPt,Binning(b));
//    TH1F* hptPull = (TH1F*)dataSet->createHistogram("hptPull", ptPull,Binning(b2));
//    TH2F* hptPull = (TH2F*)dataSet->createHistogram("hptPull",muonGenPt,Binning(b),YVar(ptPull,Binning(b)));
    TH1F* hGenPt = (TH1F*)dataSet->createHistogram("hGenPt", muonGenPt,Binning(b));
    TProfile* p_Pull = hptPull->ProfileX("p_Pull");
/*    TH1F* hptPull2 = hGenPt;
    hptPull2->Add(hSmearedPt,-1);
    hptPull2->Divide(hGenPt);
*/
    ///Plot histograms on canvas
//    TPad *pad1 = new TPad("pad1","",0,0,1,1);
//    TPad *pad2 = new TPad("pad2","",0,0,1,1);
    ///make pad1 transparent
//    pad1->SetFillStyle(0);
//    TH1F* hdummy = new TH1F("hdummy","hdummy",200,0.0,200.0);
//    ax->SetLabelSize(0.);
//    pad2->Draw();
//    pad2->cd();
//
    TCanvas *c2d = new TCanvas("c2d","c2d",600,400);
    gStyle->SetPalette(1);
    hptPull->Draw("colz");
    p_Pull->Draw("pesame");

    TCanvas *c = new TCanvas("c","c",600,400);
    gStyle->SetOptStat(kFALSE);
    hSmearedPt->SetLineColor(kRed);
    hSmearedPt->Scale(1.0,"width");
    hGenPt->Scale(1.0,"width");
    hSmearedPt->Draw("hist");
    hGenPt->Draw("histsame");
//    Float_t scale = rightmax/gPad->GetUymax();
    
/*    TH1F* p_Pullc = (TH1F*)p_Pull->Clone("p_Pullc");
    Float_t rightmax = 1.1*p_Pullc->GetMaximum();
    Float_t scale = 0.5*hSmearedPt->GetMaximum()/rightmax;
    p_Pullc->SetMarkerColor(kBlue);
//    hptPull->Draw("Y+");
    p_Pullc->Scale(scale);
    p_Pullc->Draw("pesame");
    
    TAxis *ax = p_Pullc->GetXaxis();
    TAxis *ay = p_Pullc->GetYaxis();
//    pad2->Draw();
//    pad2->Update();

//    pad1->Draw();
//    pad1->cd();
    //draw an axis on the right side
    float xmin = ax->GetXmin();
    float ymin = -0.1;
    float xmax = ax->GetXmax();
    float ymax = 0.6;
//    TGaxis *axis = new TGaxis("+L");
//    TGaxis *axis = new TGaxis("axis","axis",xmin,ymin,xmax,ymax,50510,"+L");
//    TGaxis *axis = new TGaxis(8,-0.1,8,0.8,0,0.6,50510,"+L");
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
             gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
    axis->SetLineColor(kBlue);
    axis->SetTextColor(kBlue);
    axis->Draw();
   */ 
}
