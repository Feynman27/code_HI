#include "TH1.h"
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"

#include <iostream>
#include <cmath>


TH1F* fillAsymmetry(TString fileNameIn , int bins, float xBins[],TH1F* hGenPlus,TH1F* hGenMinus,TString hname){

  float etaGen[50],phiGen[50],chargeGen[50],ptGen[50],ptGenNu[50];
  int mother[50],daughter[50];
  int ngen;

   TChain* tree = new TChain("tree","tree");
   tree->Add(fileNameIn);

   std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;
        
   tree->SetBranchAddress("mc_mu_gen_mothertype", &mother);
   tree->SetBranchAddress("mc_mu_gen_type", &daughter);
   tree->SetBranchAddress("mc_mu_charge", &chargeGen);
   tree->SetBranchAddress("mc_mu_gen_eta", &etaGen);
   tree->SetBranchAddress("mc_mu_gen_phi", &phiGen);
   tree->SetBranchAddress("mc_mu_gen_pt", &ptGen);
   tree->SetBranchAddress("mc_mu_n", &ngen);
    // --- Set branch status ---
   tree->SetBranchStatus("*",0) ;

   tree->SetBranchStatus("mc_mu_n", 1);
   tree->SetBranchStatus("mc_mu_gen_mothertype", 1);
   tree->SetBranchStatus("mc_mu_gen_type", 1);
   tree->SetBranchStatus("mc_mu_charge", 1);
   tree->SetBranchStatus("mc_mu_gen_eta", 1);
   tree->SetBranchStatus("mc_mu_gen_phi", 1);

   TH1F* hAsymm = new TH1F("hAsymm","hAsymm",bins,xBins);

   for(int iev = 0; iev < tree->GetEntries(); ++iev){
 
    tree->LoadTree(iev);
    tree->GetEntry(iev);

    for(int igen=0; igen<ngen; ++igen){

        ///fill mu+ eta
        if(
            mother[igen]==24
            &&daughter[igen]==-13
            &&chargeGen[igen]>0
            &&mc_mu_gen_pt[igen]>25.0
        ) {
            hGenPlus->Fill(fabs(etaGen[igen]));
        }
        ///fill mu- eta
        else if(
            mother[igen]==-24
            &&daughter[igen]==13
            &&chargeGen[igen]<0
        ) {
            hGenMinus->Fill(fabs(etaGen[igen]));
        }
    } //igen
   }//iev

   hAsymm = (TH1F*)hGenPlus->GetAsymmetry(hGenMinus);
   return hAsymm;
}

void plotChargeAsymmetryMC(){

    TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
    TString fileNamepp = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_pp.06.21.2013.root";
    TString fileNamepn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_pn.06.21.2013.root";
    TString fileNamenp = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_np.06.21.2013.root";
    TString fileNamenn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_nn.06.21.2013.root";

    
    float xBins[] = {0.0,0.1,0.35,0.6,0.8,1.05,1.3,1.55,1.85,2.1,2.4,2.7,3.2,4.0,5.0};
    int bins = sizeof(xBins)/sizeof(float)-1;

    TH1F* hEtaPlus_pp = new TH1F("hEtaPlus_pp","hEtaPlus_pp",bins,xBins);
    TH1F* hEtaMinus_pp = new TH1F("hEtaMinus_pp","hEtaMinus_pp",bins,xBins);
    TH1F* hEtaPlus_pn = new TH1F("hEtaPlus_pn","hEtaPlus_pn",bins,xBins);
    TH1F* hEtaMinus_pn = new TH1F("hEtaMinus_pn","hEtaMinus_pn",bins,xBins);
    TH1F* hEtaPlus_np = new TH1F("hEtaPlus_np","hEtaPlus_np",bins,xBins);
    TH1F* hEtaMinus_np = new TH1F("hEtaMinus_np","hEtaMinus_np",bins,xBins);
    TH1F* hEtaPlus_nn = new TH1F("hEtaPlus_nn","hEtaPlus_nn",bins,xBins);
    TH1F* hEtaMinus_nn = new TH1F("hEtaMinus_nn","hEtaMinus_nn",bins,xBins);

    TH1F* hAsymm_pp = new TH1F("hAsymm_pp","hAsymm_pp",bins,xBins);
    TH1F* hAsymm_pn = new TH1F("hAsymm_pn","hAsymm_pn",bins,xBins);
    TH1F* hAsymm_np = new TH1F("hAsymm_np","hAsymm_np",bins,xBins);
    TH1F* hAsymm_nn = new TH1F("hAsymm_nn","hAsymm_nn",bins,xBins);

    hAsymm_pp = fillAsymmetry(fileNamepp,bins,xBins,hEtaPlus_pp,hEtaMinus_pp,"ppTemp"); 
    hAsymm_np = fillAsymmetry(fileNamenp,bins,xBins,hEtaPlus_np,hEtaMinus_np,"npTemp"); 
    hAsymm_pn = fillAsymmetry(fileNamepn,bins,xBins,hEtaPlus_pn,hEtaMinus_pn,"pnTemp"); 
    hAsymm_nn = fillAsymmetry(fileNamenn,bins,xBins,hEtaPlus_nn,hEtaMinus_nn,"nnTemp"); 

    TCanvas *cpp0 = new TCanvas("cpp0","cpp0",600,600);
    hEtaPlus_pp->Draw("pe");
    hEtaMinus_pp->Draw("pesame");
    TCanvas *cpp1 = new TCanvas("cpp1","cpp1",600,600);
    hAsymm_pp->Draw("pe");
    hAsymm_pp->SetMarkerColor(kBlue);

    TCanvas *cpn0 = new TCanvas("cpn0","cpn0",600,600);
    hEtaPlus_pn->Draw("pe");
    hEtaMinus_pn->Draw("pesame");
    TCanvas *cpn1 = new TCanvas("cpn1","cpn1",600,600);
    hAsymm_pn->Draw("pe");
    hAsymm_pn->SetMarkerColor(kGreen);

    TCanvas *cnp0 = new TCanvas("cnp0","cnp0",600,600);
    hEtaPlus_np->Draw("pe");
    hEtaMinus_np->Draw("pesame");
    TCanvas *cnp1 = new TCanvas("cnp1","cnp1",600,600);
    hAsymm_np->Draw("pe");
    hAsymm_np->SetMarkerColor(kCyan);

    TCanvas *cnn0 = new TCanvas("cnn0","cnn0",600,600);
    hEtaPlus_nn->Draw("pe");
    hEtaMinus_nn->Draw("pesame");
    TCanvas *cnn1 = new TCanvas("cnn1","cnn1",600,600);
    hAsymm_nn->Draw("pe");
    hAsymm_nn->SetMarkerColor(kRed);

    TCanvas *cAsymm = new TCanvas("cAsymm","cAsymm",600,600);
    hAsymm_pp->GetYaxis()->SetRangeUser(-0.9,0.5);
    hAsymm_pp->Draw("pe");
    hAsymm_pn->Draw("pesame");
    hAsymm_np->Draw("pesame");
    hAsymm_nn->Draw("pesame");

    TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->AddEntry(hAsymm_pp,"pp","pe");
    leg->AddEntry(hAsymm_np,"np","pe");
    leg->AddEntry(hAsymm_pn,"pn","pe");
    leg->AddEntry(hAsymm_nn,"nn","pe");
    leg->Draw();

}
