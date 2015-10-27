#include "TString.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TList.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLatex.h"
#include "TClonesArray.h"
#include "TColor.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath> 
#include <sstream>

#include "AtlasUtils.C"
#include "AtlasStyle.C"

TH1F* getAsymmetry(TH1F* h1,TH1F* h2){

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
        if(bot<1e-6){}
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
void fillHistos(TString fileNameIn, int bins, float xBins[], TH1F* hGen, double cs, TString hname){

  float muonPtNt[50],neutrinoPtNt[50],mtNt[50],muonEtaNt[50],muonChargeNt[50];
  float muonPhiNt[50],neutrinoPhiNt[50];
  int muonMotherNt[50],neutrinoMotherNt[50];
  int ngen;
  float centralityNt;

  TChain* tree = new TChain("tree","tree");
  tree->Add(fileNameIn);

  std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;
        
  tree->SetBranchAddress("mc_mu_gen_pt", &muonPtNt);
  tree->SetBranchAddress("mc_nu_gen_pt", &neutrinoPtNt);
  tree->SetBranchAddress("mc_mu_gen_phi", &muonPhiNt);
  tree->SetBranchAddress("mc_nu_gen_phi", &neutrinoPhiNt);
  tree->SetBranchAddress("mc_mu_gen_mothertype", &muonMotherNt);
  tree->SetBranchAddress("mc_nu_gen_mothertype", &neutrinoMotherNt);
  tree->SetBranchAddress("mc_mu_gen_eta", &muonEtaNt);
  tree->SetBranchAddress("mc_mu_charge", &muonChargeNt);
  tree->SetBranchAddress("centrality", &centralityNt);
  tree->SetBranchAddress("mc_mu_n", &ngen);

   // --- Set branch status ---
  tree->SetBranchStatus("*",0) ;
  tree->SetBranchStatus("mc_mu_gen_pt", 1);
  tree->SetBranchStatus("mc_nu_gen_pt", 1);
  tree->SetBranchStatus("mc_mu_gen_phi", 1);
  tree->SetBranchStatus("mc_nu_gen_phi", 1);
  tree->SetBranchStatus("mc_mu_gen_mothertype", 1);
  tree->SetBranchStatus("mc_nu_gen_mothertype", 1);
  tree->SetBranchStatus("mc_mu_gen_eta", 1);
  tree->SetBranchStatus("mc_mu_charge", 1);
  tree->SetBranchStatus("centrality", 1);
  tree->SetBranchStatus("mc_mu_n", 1);

  // Total number of W events generated in sample
  int nW=0;
  for(int iev = 0; iev < tree->GetEntries(); ++iev){
    tree->LoadTree(iev);
    tree->GetEntry(iev);
    if(centralityNt>=0.8) continue; 

    for(int igen=0; igen<ngen; ++igen){

        ///Generator level fiducial cuts
       mtNt[igen] = TMath::Sqrt(2.0*muonPtNt[igen]*neutrinoPtNt[igen]*(1.0-TMath::Cos(muonPhiNt[igen]-neutrinoPhiNt[igen])));
       if(fabs(muonMotherNt[igen])==24) ++nW;
       if(
            fabs(muonMotherNt[igen])==24
            &&fabs(neutrinoMotherNt[igen])==24
            &&fabs(muonEtaNt[igen])>0.0&&fabs(muonEtaNt[igen])<4.0
            &&neutrinoPtNt[igen]>25.0
            &&muonPtNt[igen]>25.0
            &&mtNt[igen]>0.0
        ){ 
          hGen->Fill(fabs(muonEtaNt[igen]));
        }
    } //igen
   }//iev
   
   // Fraction of sampled xsec
   double sf = 1./(nW)*cs;
   hGen->Scale(sf);
   // Scale by eta bw
   hGen->Scale(1.0,"width"); 

}

void plotEtaDistrosAU2CT10(){

    TString fileNameDataOut = "PowhegPythia8_AU2CT10_ChargeEtaDistributions";
    TFile *outFile = new TFile(fileNameDataOut+".root","RECREATE");
    TString base = "/usatlas/u/tbales/scratch/data_files/";
    //PowPy8 samples
    //pp
    TString fileNameIn_ppPlus = base+"MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_pp.10.03.2013.root";
    TString fileNameIn_ppMinus = base+"MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_pp.10.03.2013.root";
    //np
    TString fileNameIn_npPlus = base+"MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_np.10.03.2013.root";
    TString fileNameIn_npMinus = base+"MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_np.10.03.2013.root";
    //pn
    TString fileNameIn_pnPlus = base+"MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_pn.10.03.2013.root";
    TString fileNameIn_pnMinus = base+"MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_pn.10.03.2013.root";
    //nn
    TString fileNameIn_nnPlus = base+"MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_nn.10.03.2013.root";
    TString fileNameIn_nnMinus = base+"MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_nn.10.03.2013.root";

    double WCrossSection_ppPlusnlo = 2.1120e-9/64.0e-3*1.0e9;
    double WCrossSection_ppMinusnlo = 1.2420e-9/64.0e-3*1.0e9;
    double wt_pp = 0.155;
    double WCrossSection_nppnPlusnlo = 1.6690e-9/64.0e-3*1.0e9;
    double WCrossSection_nppnMinusnlo = 1.6540e-9/64.0e-3*1.0e9;
    double wt_nppn = 0.478;
    double WCrossSection_nnPlusnlo = 1.2530e-9/64.0e-3*1.0e9;
    double WCrossSection_nnMinusnlo = 2.0920e-9/64.0e-3*1.0e9;
    double wt_nn = 0.367;

    //float xBinsTheory[] = {0.0,0.1,0.35,0.6,0.8,1.05,1.37,1.52,1.74,2.1,2.50,3.1,3.5,4.0,4.4};
    //D0
    float xBinsTheory[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,3.1,3.5,4.0,4.4};
    int binsTheory = sizeof(xBinsTheory)/sizeof(float)-1;

    TH1F* hEtaPluspp = new TH1F("hEtaPluspp","hEtaPluspp",binsTheory,xBinsTheory);
    hEtaPluspp->Sumw2();
    TH1F* hEtaPlusnn = new TH1F("hEtaPlusnn","hEtaPlusnn",binsTheory,xBinsTheory);
    hEtaPlusnn->Sumw2();
    TH1F* hEtaPlusnp = new TH1F("hEtaPlusnp","hEtaPlusnp",binsTheory,xBinsTheory);
    hEtaPlusnp->Sumw2();
    TH1F* hEtaPluspn = new TH1F("hEtaPluspn","hEtaPluspn",binsTheory,xBinsTheory);
    hEtaPluspn->Sumw2();
    TH1F* hEtaPlusnppn = new TH1F("hEtaPlusnppn","hEtaPlusnppn",binsTheory,xBinsTheory);
    hEtaPlusnppn->Sumw2();
    TH1F* hEtaMinuspp = new TH1F("hEtaMinuspp","hEtaMinuspp",binsTheory,xBinsTheory);
    hEtaMinuspp->Sumw2();
    TH1F* hEtaMinusnn = new TH1F("hEtaMinusnn","hEtaMinusnn",binsTheory,xBinsTheory);
    hEtaMinusnn->Sumw2();
    TH1F* hEtaMinusnp = new TH1F("hEtaMinusnp","hEtaMinusnp",binsTheory,xBinsTheory);
    hEtaMinusnp->Sumw2();
    TH1F* hEtaMinuspn = new TH1F("hEtaMinuspn","hEtaMinuspn",binsTheory,xBinsTheory);
    hEtaMinuspn->Sumw2();
    TH1F* hEtaMinusnppn = new TH1F("hEtaMinusnppn","hEtaMinusnppn",binsTheory,xBinsTheory);
    hEtaMinusnppn->Sumw2();

    fillHistos(fileNameIn_ppPlus,binsTheory,xBinsTheory,hEtaPluspp,WCrossSection_ppPlusnlo,"ppTemp");
    fillHistos(fileNameIn_ppMinus,binsTheory,xBinsTheory,hEtaMinuspp,WCrossSection_ppMinusnlo,"ppTemp");

    fillHistos(fileNameIn_nnPlus,binsTheory,xBinsTheory,hEtaPlusnn,WCrossSection_nnPlusnlo,"nnTemp");
    fillHistos(fileNameIn_nnMinus,binsTheory,xBinsTheory,hEtaMinusnn,WCrossSection_nnMinusnlo,"nnTemp");

    fillHistos(fileNameIn_npPlus,binsTheory,xBinsTheory,hEtaPlusnp,WCrossSection_nppnPlusnlo,"npTemp");
    fillHistos(fileNameIn_npMinus,binsTheory,xBinsTheory,hEtaMinusnp,WCrossSection_nppnMinusnlo,"npTemp");

    fillHistos(fileNameIn_pnPlus,binsTheory,xBinsTheory,hEtaPluspn,WCrossSection_nppnPlusnlo,"pnTemp");
    fillHistos(fileNameIn_pnMinus,binsTheory,xBinsTheory,hEtaMinuspn,WCrossSection_nppnMinusnlo,"pnTemp");

    hEtaPlusnppn->Add(hEtaPlusnp,hEtaPluspn,0.5,0.5); 
    hEtaMinusnppn->Add(hEtaMinusnp,hEtaMinuspn,0.5,0.5);

    TH1F* hAsymmpp = getAsymmetry(hEtaPluspp,hEtaMinuspp);
    TH1F* hAsymmnn = getAsymmetry(hEtaPlusnn,hEtaMinusnn);
    TH1F* hAsymmnppn = getAsymmetry(hEtaPlusnppn,hEtaMinusnppn);

    outFile->cd();
    hEtaPluspp->Write("hEtaPluspp"); hEtaMinuspp->Write("hEtaMinuspp"); hAsymmpp->Write("hAsymmpp");
    hEtaPlusnn->Write("hEtaPlusnn"); hEtaMinusnn->Write("hEtaMinusnn"); hAsymmnn->Write("hAsymmnn");
    hEtaPlusnppn->Write("hEtaPlusnppn"); hEtaMinusnppn->Write("hEtaMinusnppn"); hAsymmnppn->Write("hAsymmnppn");
}
