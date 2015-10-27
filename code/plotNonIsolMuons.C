#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TChain.h"
#include "TH1F.h"
#include "TMath.h"
#include "TLatex.h"
#include "TTree.h"
#include "TColor.h"
#include "TLegendEntry.h"
#include "TPad.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath> 
#include <sstream>

///Return number of events in this centrality class
double getNEvents(TString fileNameMuonIn){
  TChain* tree = new TChain("tree","tree");
  std::cout << "Getting number of events for " << fileNameMuonIn << std::endl;
  tree->Add(fileNameMuonIn);

  float centralityNt;
  tree->SetBranchAddress("centrality", &centralityNt); 
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("centrality",1);
 
  TH1F* hCent = new TH1F("hCent","hCent",6,0.0,0.8);
  for(int iev=0; iev<tree->GetEntries(); ++iev){

    tree->LoadTree(iev);
    tree->GetEntry(iev);

    hCent->Fill(centralityNt);
  }
  std::cout << "Number of entries: " << hCent->GetEntries() << std::endl; 
  return hCent->GetEntries();
}

//////////////////////////////////////////
//return ratio of hMC over hData
/////////////////////////////////////////
TH1F* getMCDataRatio(TH1F* hRatio, TH1F* hMC, TH1F* hData){
        std::cout << "Calculation MC/Data ratios..." << std::endl;
		hRatio->Divide(hMC,hData);
        std::cout << "Done." << std::endl;
		return hRatio;
}

void fixRatioPlot(TH1F* hRatio){

    for(int i=1; i<=hRatio->GetNbinsX(); ++i){

        float xBin = hRatio->GetBinContent(i);
        if(xBin<=0.0) hRatio->SetBinContent(i,0.0);
        else hRatio->SetBinContent(i,xBin);
    }
}

TH1F* fillHisto(TString fileNameMuonIn , int bins, double xBins[],TString hname,  bool isMc = false, bool doAntiIsolation = false, bool isWMC = false){

  TH1F* hMuTemp = new TH1F(hname,hname,bins,xBins);

  int valNt[30], nmu;
  int trig1,trig2,trig3,trig4,trig5;
  int matched1[50], matched2[50], matched3[50];
  float ptNt[30], centralityNt,etaNt[30],eLossNt[30],scatNt[30];
  float ptconeNt[30], mtNt[30],nu_ptNt;
  int ZDYNt[30], promptNt[30];

  TChain* tree = new TChain("tree","tree");
  std::cout << "Filling the tree for " << fileNameMuonIn << std::endl;
  tree->Add(fileNameMuonIn);
  // --- Set branch adresses ---
  tree->SetBranchAddress("val", &valNt);
  tree->SetBranchAddress("pt", &ptNt);
  tree->SetBranchAddress("eta", &etaNt);
  tree->SetBranchAddress("eLoss",&eLossNt);
  tree->SetBranchAddress("scat",&scatNt);
  tree->SetBranchAddress("mu_muid_n", &nmu);
  tree->SetBranchAddress("centrality", &centralityNt); 
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC",&trig1);
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10",&trig2);
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20",&trig3);
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC_Matched20",&matched1);
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10_Matched20",&matched2);
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20_Matched20",&matched3);
  tree->SetBranchAddress("ptcone20ID3", &ptconeNt);
  tree->SetBranchAddress("nu_pt", &nu_ptNt);
  tree->SetBranchAddress("mt", &mtNt);
  tree->SetBranchAddress("ZDY", &ZDYNt);
  if(isMc) tree->SetBranchAddress("prompt", &promptNt);

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("val", 1);
  if(isMc) tree->SetBranchStatus("prompt", 1);
  tree->SetBranchStatus("pt", 1);
  tree->SetBranchStatus("eta", 1);
  tree->SetBranchStatus("eLoss",1);
  tree->SetBranchStatus("scat",1);
  tree->SetBranchStatus("mu_muid_n", 1);
  tree->SetBranchStatus("centrality", 1); 
  tree->SetBranchStatus("ptcone20ID3", 1);
  tree->SetBranchStatus("nu_pt", 1);
  tree->SetBranchStatus("mt", 1);
  tree->SetBranchStatus("ZDY", 1);
  if(!isMc){
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC", 1); 
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10", 1); 
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20", 1); 
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC_Matched20", 1); 
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10_Matched20", 1); 
      tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20_Matched20", 1); 
  }

 std::cout << "Number of entries: " << tree->GetEntries() << std::endl; 
 
 for(int iev=0; iev<tree->GetEntries(); ++iev){

    tree->LoadTree(iev);
    tree->GetEntry(iev);

  if(trig1||trig2||trig3||isMc){
    for(int imu=0; imu<nmu; ++imu){

          if (
                (valNt[imu]>11) 
                && abs(scatNt[imu])<4.0
                && abs(eLossNt[imu])<0.5
                && abs(etaNt[imu])>0.1 
                && abs(etaNt[imu])<2.4 
                && ptNt[imu]>=4.0
                && centralityNt>=0. 
                && centralityNt<=0.8
                && ( isMc || ( (trig1&&matched1[imu]) ||
                    (trig2&&matched2[imu])
                    ||(trig3&&matched3[imu]) ) )	
                && ZDYNt[imu]==0
                && ( (ptconeNt[imu]/ptNt[imu] >= 0.1/*&&nu_ptNt < 25.&&mtNt<40.*/) || !doAntiIsolation)
                //&& (nu_ptNt < 25. || !doAntiIsolation)
                && ( promptNt[imu]==24|| !isWMC )
		        ){
                   hMuTemp->Fill(ptNt[imu]);
                 }
    }
  }
 }

 return hMuTemp;

}

void plotNonIsolMuons(){

  TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
  TString fileNameMuonIn = baseString+"HardProbesFiles/HISingleMuonHardProbesData.04.17.2013.root";  
  
  //J1 1 muon-filter 
  TString fileNameMCJ1In = baseString+"MonteCarloFiles/QCD/HISingleMuonJ1_1muPYTHIADataOverlay.04.13.2013.root";
  //J2 1 muon-filter 
  TString fileNameMCJ2In = baseString+"MonteCarloFiles/QCD/HISingleMuonJ2_1muPYTHIADataOverlay.04.13.2013.root";
  //J3 1 muon-filter 
  TString fileNameMCJ3In = baseString+"MonteCarloFiles/QCD/HISingleMuonJ3_1muPYTHIADataOverlay.04.13.2013.root";

  //W MC 
  TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.04.12.2013.root";

  ///binning
  double xBins[60] ;

  double xlo = 0.0; double xhi = 50.0;  
  double nbins = 50;
  double binw = (xhi-xlo)/nbins;
  int ib = 0;
  for(double i=xlo ; i<xhi ; i+=binw){
       
            xBins[ib] = i;
            ib++;	
  }
  xlo = 50.0; xhi = 75.0; 
  nbins = 5;
  binw = (xhi-xlo)/nbins;
  for(double i=xlo ; i<xhi ; i+=binw){

        xBins[ib] = i;
        ib++;
   }	
  xlo = 75.0; xhi = 100.0; 
  nbins = 3;
  binw = (xhi-xlo)/nbins;
  for(double i=xlo ; i<xhi ; i+=binw){

        xBins[ib] = i;
        ib++;
  }
  xlo = 100.0; xhi = 250.0; 
  nbins = 2;
  binw = (xhi-xlo)/nbins;
  for(double i=xlo ; i<xhi ; i+=binw){

        xBins[ib] = i;
        ib++;
  }
  int binNum = sizeof(xBins)/sizeof(double) - 1;

  ///Fill anti-isol histo from data muons
  TH1F* hAntiIsolData = fillHisto(fileNameMuonIn,binNum,xBins,"hAntiIsolData",false,true);
  ///Fill histo of preselected(+ Zveto) muons from data 
  TH1F* hDataPS = fillHisto(fileNameMuonIn,binNum,xBins,"hDataPS");
  ///Fill Jxmu preselected(+Zveto) histos from MC
  TH1F* hJ1mu = fillHisto(fileNameMCJ1In,binNum,xBins,"hJ1mu",true,true);
  TH1F* hJ2mu = fillHisto(fileNameMCJ2In,binNum,xBins,"hJ2mu",true,true);
  TH1F* hJ3mu = fillHisto(fileNameMCJ3In,binNum,xBins,"hJ3mu",true,true);
  ///Fill anti-isolated W MC histo with 
  TH1F* hmcWSet = fillHisto(fileNameMCWIn,binNum,xBins,"hmcWSet",true,true,true); 

  ///Weight the Jxmu by cross-sections
  double ncoll = 452.0;
  double wtJ1 = 6.6523e-03*1.8770e-4/64.0e-3; double wtJ2 = 1.4941e-02*8.2788e-6/64.0e-3; double wtJ3 = 2.4284e-02*2.9426e-7/64.0e-3;
  ///Weight by W cross-section
  double wtW = 2.8E-09/64.0e-3*0.95;  
  double arrCentWidth = 0.8; ///0-80%
  double evData = 1.0e9; //number of events sampled in data
  double scaleFactor = arrCentWidth*ncoll*evData;

  ///Get number of events in each mc ds
  double nMCJ1 = getNEvents(fileNameMCJ1In);
  double nMCJ2 = getNEvents(fileNameMCJ2In);
  double nMCJ3 = getNEvents(fileNameMCJ3In);
  double nMCW =  getNEvents(fileNameMCWIn);

  TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",binNum,xBins);
  ///Add the weighted Jxmu samples
  hmcQCDSet->Add(hJ1mu,hJ2mu,(wtJ1)/(nMCJ1)*scaleFactor, (wtJ2)/(nMCJ2)*scaleFactor); 
  hmcQCDSet->Add(hJ3mu,(wtJ3)/(nMCJ3)*scaleFactor);

  ///Scale W anti-isolated sample by cross-section
//  hmcWSet->Scale((wtW)/(nMCW)*scaleFactor);
  ///or use number of estimated non-isolated signal counts from data
  hmcWSet->Scale(249.0/hmcWSet->Integral());


  ///scale factor for scaling JxMu to data in ctrl region 10<pt<20
  std::cout << "Control region : " << hDataPS->GetBinLowEdge(11) << " - " << hDataPS->GetBinLowEdge(21) << std::endl; 
  //double scaleFactorCtrl = hDataPS->Integral(11,21)/hmcQCDSet->Integral(11,21);
  ///Obtained from PS data and PS QCD MC
  double scaleFactorCtrl = 0.4;

  std::cout << "Scale factor in ctrl region of pre-selected data = " << scaleFactorCtrl << std::endl;
  std::cout << "Scale factor in ctrl region of non-isolated data subset = " << hAntiIsolData->Integral(11,21)/hmcQCDSet->Integral(11,21) << std::endl;  
  ///Scale Jxmu MC to data in ctrl region
  hmcQCDSet->Scale(scaleFactorCtrl);
  ///Scale to non-isolated data subset
  //hmcQCDSet->Scale(hAntiIsolData->Integral(11,120)/hmcQCDSet->Integral(11,120));
  ///Scale to non-isolated data subset's ctrl region
  //hmcQCDSet->Scale(hAntiIsolData->Integral(11,21)/hmcQCDSet->Integral(11,21));

  ///Subtract anti-isolated signal from anti-isolated data
  TH1F* hAntiIsolDataSigSub = new TH1F("hAntiIsolDataSigSub","hAntiIsolDataSigSub",binNum,xBins);
//  TH1F* hAntiIsolDataSigSub = (TH1F*)hAntiIsolData->Clone("hAntiIsolDataSigSub");
  hAntiIsolDataSigSub->Add(hAntiIsolData,hmcWSet,+1.0,-1.0);

  TH1F* hAntiIsolDataSigSubc = (TH1F*)hAntiIsolDataSigSub->Clone("hAntiIsolDataSigSub");
  hAntiIsolDataSigSubc->SetMarkerStyle(22);
  hAntiIsolDataSigSubc->SetMarkerColor(kGray+1);

  TH1F* hmcWSetc = (TH1F*)hmcWSet->Clone("hmcWSetc");
  hmcWSetc->SetFillColor(kWhite);

  TH1F* hmcQCDSetc = (TH1F*)hmcQCDSet->Clone("hmcQCDSetc");
  hmcQCDSetc->SetFillColor(kAzure-9);

  TH1F* hDataPSc = (TH1F*)hDataPS->Clone("hDataPSc");
  hDataPSc->SetMarkerStyle(24);

  TCanvas* c = new TCanvas("c","c",600,600);

  c->Divide(1, 2);
  TPad* canvas_up = (TPad*)c->GetListOfPrimitives()->FindObject("c_1");
  TPad* canvas_dw = (TPad*)c->GetListOfPrimitives()->FindObject("c_2");

    double up_height     = 0.8; 
    double dw_correction = 1.60;
    double font_size_dw  = 0.08;
    double dw_height    = (1. - up_height) * dw_correction;

    ///set pad size
    canvas_up->SetPad(0., 1 - up_height, 1., 1.);
    canvas_up->SetFillColor(0);
    canvas_up->SetBorderMode(0);
    canvas_up->SetBorderSize(2);
    canvas_up->SetTickx(1);
    canvas_up->SetTicky(1);
    canvas_up->SetLeftMargin(0.16);
    canvas_up->SetRightMargin(0.05);
    canvas_up->SetTopMargin(0.05);
    canvas_up->SetBottomMargin(0.16);
    canvas_up->SetFrameBorderMode(0);
    canvas_up->SetFrameBorderMode(0);
 
    canvas_dw->SetPad(0., 0., 1., dw_height);
    canvas_dw->Range(-0.3639066,-0.7754386,2.546497,1.31186);
    canvas_dw->SetFillColor(0);
    canvas_dw->SetBorderMode(0);
    canvas_dw->SetBorderSize(2);
    canvas_dw->SetTickx(1);
    canvas_dw->SetTicky(1);
    canvas_dw->SetLeftMargin(0.159396);
    canvas_dw->SetRightMargin(0.05033557);
    canvas_dw->SetTopMargin(0.005681818);
    canvas_dw->SetBottomMargin(0.3715035);
    canvas_dw->SetFrameBorderMode(0);
    canvas_dw->SetFrameBorderMode(0);

    canvas_up->SetFrameFillColor(0);
    canvas_up->SetFillColor(0);
    canvas_dw->SetFillColor(0);
    canvas_dw->SetFrameFillColor(0);


/*  TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
  leg->SetTextFont(gStyle->GetTextFont());
  leg->SetTextSize(gStyle->GetTextSize());
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(hmcQCDSetc, "Jx#mu", "f");
  leg->AddEntry(hmcWSetc, "Anti-Isolated W#rightarrow#mu#nu", "l");
  leg->AddEntry(hDataPSc, "Preselected Data", "pe");
  leg->AddEntry(hAntiIsolData, "Anti-Isolated Data", "pe");
  leg->AddEntry(hAntiIsolDataSigSubc, "Anti-Isolated,Signal Subtracted Data", "pe");
  leg->Draw();
*/
  //TLatex l;
  //l.SetNDC();
  //l.DrawLatex(0.74,0.89,"#frac{#Sigma p_{T}^{trk}(#Delta R<0.2,p_{T}^{trk}>3GeV)}{p_{T}^{#mu}} < 0.1");

  ///Draw spectra
  hmcWSetc->Scale(1.0,"width");
  hmcQCDSetc->Scale(1.0,"width");
  hDataPSc->Scale(1.0,"width");
  hAntiIsolData->Scale(1.0,"width");
  hAntiIsolDataSigSubc->Scale(1.0,"width");

  canvas_up->cd();
 ///MC/Data histo
  TH1F* hRatio = (TH1F*)hAntiIsolDataSigSubc->Clone("hRatio");
  hRatio = getMCDataRatio(hRatio,hmcQCDSetc,hAntiIsolDataSigSubc);
  fixRatioPlot(hRatio);

  hmcQCDSetc->GetXaxis()->SetTitle("p_{T}[GeV]");
  hmcQCDSetc->GetYaxis()->SetTitle("Muons");

  hmcQCDSetc->Draw("hist");
  hDataPSc->Draw("pesame");
  hAntiIsolData->Draw("pesame");
  hmcWSetc->Draw("histsame");
  hAntiIsolDataSigSubc->Draw("pesame");
  //l.Draw(); 
   TLegend *leg = new TLegend(0.340604,0.7692308,0.4748322,0.9423077,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03671329);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("hmcQCDSetc","Anti-Isolated Jx#mu","f");

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#99ccff");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("hmcWSetc","Anti-Isolated W#rightarrow#mu#nu","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("hDataPSc","Preselected Data","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hAntiIsolData","Anti-Isolated Data","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hAntiIsolDataSigSub","Anti-Isolated,Signal Subtracted Data","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#999999");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   leg->Draw();

    ///Draw MC/Data underneath
    canvas_dw->cd();
    /// font size
    hRatio->GetXaxis()->SetRange(11,101);
    hRatio->GetXaxis()->SetTitle("p_{T}[GeV]");
    hRatio->GetXaxis()->SetLabelFont(42);
    hRatio->GetXaxis()->SetLabelSize(0.09);
    hRatio->GetXaxis()->SetTitleSize(0.11);
    hRatio->GetXaxis()->SetTitleOffset(0.8);
    hRatio->GetXaxis()->SetTitleFont(42);
    hRatio->GetYaxis()->SetTitle("MC/Data");
    hRatio->GetYaxis()->SetLabelFont(42);
    hRatio->GetYaxis()->SetLabelSize(0.09);
    hRatio->GetYaxis()->SetTitleSize(0.11);
    hRatio->GetYaxis()->SetTitleOffset(0.7);


    hRatio->Draw();


}

