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
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLatex.h"
#include "TClonesArray.h"
#include "TColor.h"

#include "RooCurve.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooMinuit.h"
#include "RooBinning.h"
#include "RooHistPdf.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath> 
#include <sstream>



#include "AtlasUtils.C"
#include "AtlasStyle.C"

////////////////////////
//Prevent trailing digits
//in figure labeling
////////////////////////
TString format(float value, int precision=3) {
  std::stringstream svalue;
  svalue  << std::setprecision(precision) << value;
  return svalue.str();
}

// Convert TH1 to TGraph
TH1F* convertToHisto(const TGraphAsymmErrors* gr, TString hname, int nBins, float xBins[]){

    TH1F* h = new TH1F(hname,hname,nBins,xBins);
    for(int i = 0; i<gr->GetN(); ++i){

        double x = gr->GetX()[i]; 
        double y = gr->GetY()[i];
        double xerr = gr->GetEXhigh()[i];
        double yerr = gr->GetEYhigh()[i];
        h->SetBinContent(i+1,y);
        h->SetBinError(i+1,yerr);
    }//i

    return h;
}

// Convert TH1 to TGraph
TGraphAsymmErrors* convertToGraph(const TH1F* h){

    TGraphAsymmErrors *gr = new TGraphAsymmErrors(h->GetNbinsX());
    for(int i = 0; i<gr->GetN(); ++i){

        double xgr = h->GetBinCenter(i+1);
        double ygr = h->GetBinContent(i+1);
        double xerr = h->GetBinWidth(i+1)/2.0;
        double yerr = h->GetBinError(i+1);
        gr->SetPoint(i,xgr,ygr);
        gr->SetPointError(i,0.0,0.0,yerr,yerr);
    }//i

    return gr;
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
///////////////////////////////////////////
//writeYieldsToSpreadsheet
//////////////////////////////////////////
void writeYieldsToSpreadsheet(std::ostream& outputFile, TString ieta,TString nObs,TString nBkg, TString sigCounts, 
				TString errStat, TString errSystUncorr, TString errSystCorr){
	std::cout << "Writing yields to spreadsheet..." << std::endl;
	outputFile << ieta << "," << nObs << "," << nBkg << "," <<  sigCounts << "," << errStat << "," << errSystUncorr << "," << errSystCorr << std::endl;
}
void writeAsymmToSpreadsheet(std::ostream& outputFile, TString ieta, TString asymm, TString errStat, TString errSystUnc, TString errSystCor){
	std::cout << "Writing asymmetry to spreadsheet..." << std::endl;
	outputFile << ieta << "," << asymm << "," << errStat << "," << errSystUnc << "," << errSystCor << std::endl;
}

TH1F* scaleToData(TH1F* hEta,double sf, const TGraphAsymmErrors* grAwEta){

    for(int i=1; i<hEta->GetNbinsX(); ++i){

        double ypt = sf*hEta->GetBinContent(i);
        double yErr = sf*hEta->GetBinError(i);
        hEta->SetBinContent(i,ypt);
        hEta->SetBinError(i,yErr);
    }
    return hEta;
}
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

TGraphAsymmErrors* scaleToBinaryColl(TGraphAsymmErrors* gr, double scaleFactor,bool doNcollError = false){

    // Wt'd average relative error on Ncoll in all 6 
    // centrality bins
    const double avgErrNcoll = 0.095;
    for(int igr=0; igr<gr->GetN(); ++igr){
        double xpt = gr->GetX()[igr];
        double xErr = gr->GetEXhigh()[igr];
        double ypt = gr->GetY()[igr]*scaleFactor;
        double yErr = gr->GetEYhigh()[igr]*scaleFactor;
	double yNcollAbs = ypt;
        if(doNcollError) {
	  yNcollAbs*=avgErrNcoll;
	  yErr = TMath::Sqrt(TMath::Power(yNcollAbs,2)+TMath::Power(yErr,2));
	}
	// Add Ncoll and 
        gr->SetPoint(igr,xpt,ypt);
        gr->SetPointError(igr,xErr,xErr,yErr,yErr);
    }
    return gr;
}

TH1F* fillAsymmetry(int bins, float xBins[], TH1F* hGenPlus, TH1F* hGenMinus,double csPlus, double csMinus,TString hname = ""){

  TH1F* hAsymm = new TH1F("hAsymm","hAsymm",bins,xBins);

  hAsymm = getAsymmetry(hGenPlus,hGenMinus);
  return hAsymm;
}
void getEvGen(TString fileNameIn, int nWPos, int nWNeg){

  float muonPtNt[50],neutrinoPtNt[50],mtNt[50],muonEtaNt[50],muonChargeNt[50];
  int muonMotherNt[50],neutrinoMotherNt[50];
  int ngen;

  TChain* tree = new TChain("T","T");
  tree->Add(fileNameIn);

  std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;
        
  tree->SetBranchAddress("muonPt", &muonPtNt);
  tree->SetBranchAddress("neutrinoPt", &neutrinoPtNt);
  tree->SetBranchAddress("muonMother", &muonMotherNt);
  tree->SetBranchAddress("neutrinoMother", &neutrinoMotherNt);
  tree->SetBranchAddress("mt", &mtNt);
  tree->SetBranchAddress("muonEta", &muonEtaNt);
  tree->SetBranchAddress("muonCharge", &muonChargeNt);
  tree->SetBranchAddress("nMuons", &ngen);

   // --- Set branch status ---
  tree->SetBranchStatus("*",0) ;
  tree->SetBranchStatus("muonPt", 1);
  tree->SetBranchStatus("neutrinoPt", 1);
  tree->SetBranchStatus("muonMother", 1);
  tree->SetBranchStatus("neutrinoMother", 1);
  tree->SetBranchStatus("mt", 1);
  tree->SetBranchStatus("muonEta", 1);
  tree->SetBranchStatus("muonCharge", 1);
  tree->SetBranchStatus("nMuons", 1);

  // Total number of W events generated in sample
  nWPos = 0, nWNeg=0;
  for(int iev = 0; iev < tree->GetEntries(); ++iev){
 
    tree->LoadTree(iev);
    tree->GetEntry(iev);

    for(int igen=0; igen<ngen; ++igen){

        ///Generator level fiducial cuts
       if(muonMotherNt[igen]==24) ++nWPos;
       else if(muonMotherNt[igen]==-24) ++nWNeg;
    } //igen
   }//iev
   
   std::cout << "W+ --> mu+: " << nWPos << std::endl;
   std::cout << "W- --> mu-: " << nWNeg << std::endl;
}
void fillMSTWHistos(TString fileNameIn , int bins, float xBins[], TH1F* hGenPlus, TH1F* hGenMinus,double csPlus, double csMinus,TString hname){

  float muonPtNt[50],neutrinoPtNt[50],mtNt[50],muonEtaNt[50],muonChargeNt[50];
  int muonMotherNt[50],neutrinoMotherNt[50];
  int ngen;

  TChain* tree = new TChain("T","T");
  tree->Add(fileNameIn);

  std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;
        
  tree->SetBranchAddress("muonPt", &muonPtNt);
  tree->SetBranchAddress("neutrinoPt", &neutrinoPtNt);
  tree->SetBranchAddress("muonMother", &muonMotherNt);
  tree->SetBranchAddress("neutrinoMother", &neutrinoMotherNt);
  tree->SetBranchAddress("mt", &mtNt);
  tree->SetBranchAddress("muonEta", &muonEtaNt);
  tree->SetBranchAddress("muonCharge", &muonChargeNt);
  tree->SetBranchAddress("nMuons", &ngen);

   // --- Set branch status ---
  tree->SetBranchStatus("*",0) ;
  tree->SetBranchStatus("muonPt", 1);
  tree->SetBranchStatus("neutrinoPt", 1);
  tree->SetBranchStatus("muonMother", 1);
  tree->SetBranchStatus("neutrinoMother", 1);
  tree->SetBranchStatus("mt", 1);
  tree->SetBranchStatus("muonEta", 1);
  tree->SetBranchStatus("muonCharge", 1);
  tree->SetBranchStatus("nMuons", 1);

  // Total number of W events generated in sample
  int nWPos = 0, nWNeg=0;
  for(int iev = 0; iev < tree->GetEntries(); ++iev){
 
    tree->LoadTree(iev);
    tree->GetEntry(iev);

    for(int igen=0; igen<ngen; ++igen){

        ///Generator level fiducial cuts
       if(muonMotherNt[igen]==24) ++nWPos;
       else if(muonMotherNt[igen]==-24) ++nWNeg;
       if(
            fabs(muonMotherNt[igen])==24
            &&fabs(neutrinoMotherNt[igen])==24
            &&fabs(muonEtaNt[igen])>0.1&&fabs(muonEtaNt[igen])<2.4
            &&neutrinoPtNt[igen]>25.0
            &&muonPtNt[igen]>25.0
            &&mtNt[igen]>40.0
        ){ 
        ///fill mu+ eta
        if(muonChargeNt[igen]>0.){
            hGenPlus->Fill(fabs(muonEtaNt[igen]));
        }
        ///fill mu- eta
        else if(muonChargeNt[igen]<0.)
         {
            hGenMinus->Fill(fabs(muonEtaNt[igen]));
         }
    }
    } //igen
   }//iev
   
//   std::cout << hname << " W+:" << hGenPlus->GetEntries() << " cross-section: " << csPlus << " scale factor: " << 1.0/(nEvents/2.)*csPlus << std::endl;
//   std::cout << hname << " W-:" << hGenMinus->GetEntries() << " cross-section: " << csMinus << " scale factor: " << 1.0/(nEvents/2.)*csMinus << std::endl;
   std::cout << "W+ --> mu+: " << nWPos << std::endl;
   std::cout << "W- --> mu-: " << nWNeg << std::endl;
   std::cout << "POW Aw+: " << hGenPlus->GetEntries()/nWPos << std::endl;
   std::cout << "POW Aw-: " << hGenMinus->GetEntries()/nWNeg << std::endl;
   //TH1F* hAsymm = (TH1F*)hGenPlus->Clone("hAsymm");
   TH1F* hGenPlusc = (TH1F*)hGenPlus->Clone("hGenPlusc");
   TH1F* hGenMinusc = (TH1F*)hGenMinus->Clone("hGenMinusc");
   /*for(int i =1; i<hAsymm->GetNbinsX(); ++i){
     float asymm = (hGenPlusc->GetBinContent(i)-hGenMinusc->GetBinContent(i))/(hGenPlusc->GetBinContent(i)+hGenMinusc->GetBinContent(i));
     hAsymm->SetBinContent(i,asymm);
     float asymmErrTop = TMath::Sqrt(hGenPlusc->GetBinContent(i)+hGenMinusc->GetBinContent(i))/(hGenPlusc->GetBinContent(i)-hGenMinusc->GetBinContent(i));
     float asymmErrBottom = TMath::Sqrt(hGenPlusc->GetBinContent(i)+hGenMinusc->GetBinContent(i))/(hGenPlusc->GetBinContent(i)+hGenMinusc->GetBinContent(i));
     float asymmErr = asymmErrTop*TMath::Sqrt(TMath::Power(asymmErrBottom,2)+TMath::Power(asymmErrTop,2));
     hAsymm->SetBinError(i,asymmErr);
     std::cout << "asymm " << asymm << " =?" << hAsymm->GetBinContent(i) << std::endl; 
     
   }*/
//   exit(0);//hack
   // Scale by eta bw
//   double sfPlus = (hGenPlus->GetEntries()+hGenMinus->GetEntries())/(nEvents)*csPlus;
   double sfPlus = 1./(nWPos)*csPlus;
   hGenPlus->Scale(sfPlus);
//   double sfMinus = (hGenMinus->GetEntries()+hGenMinus->GetEntries())/(nEvents)*csMinus;
   std::cout << "AwPOW: " << hGenMinus->GetEntries()/nWNeg << std::endl;
   double sfMinus = 1.0/(nWNeg)*csMinus;
   hGenMinus->Scale(sfMinus);
   hGenPlus->Scale(1.0,"width"); hGenMinus->Scale(1.0,"width");

//   return hAsymm;
}

void fillCT10Histos(TString fileNameIn , int bins, float xBins[], TH1F* hGen, double cs, TString hname){

  float etaGen[50],chargeGen[50];
  float nuGenPt[50], nuPhiGen[50];
  int mother[50],daughter[50];
  int ngen;
  float centrality;
  float muGenPt[50], muPhiGen[50];

   TChain* tree = new TChain("tree","tree");
   tree->Add(fileNameIn);

   std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;
        
  tree->SetBranchAddress("mc_mu_gen_pt", &muGenPt);
  tree->SetBranchAddress("mc_nu_gen_pt", &nuGenPt);
  tree->SetBranchAddress("mc_nu_gen_phi", &nuPhiGen);
  tree->SetBranchAddress("mc_mu_gen_mothertype", &mother);
  tree->SetBranchAddress("mc_mu_gen_type", &daughter);
  tree->SetBranchAddress("mc_mu_charge", &chargeGen);
  tree->SetBranchAddress("mc_mu_gen_eta", &etaGen);
  tree->SetBranchAddress("mc_mu_gen_phi", &muPhiGen);
  tree->SetBranchAddress("mc_mu_n", &ngen);
  tree->SetBranchAddress("centrality", &centrality);

   // --- Set branch status ---
  tree->SetBranchStatus("*",0) ;
  tree->SetBranchStatus("mc_mu_n", 1);
  tree->SetBranchStatus("centrality", 1);
  tree->SetBranchStatus("mc_nu_gen_pt", 1);
  tree->SetBranchStatus("mc_nu_gen_phi", 1);
  tree->SetBranchStatus("mc_mu_gen_pt", 1);
  tree->SetBranchStatus("mc_mu_gen_mothertype", 1);
  tree->SetBranchStatus("mc_mu_gen_type", 1);
  tree->SetBranchStatus("mc_mu_charge", 1);
  tree->SetBranchStatus("mc_mu_gen_eta", 1);
  tree->SetBranchStatus("mc_mu_gen_phi", 1);


   int nW = 0;
   for(int iev = 0; iev < tree->GetEntries(); ++iev){
 
    tree->LoadTree(iev);
    tree->GetEntry(iev);

    for(int igen=0; igen<ngen; ++igen){

        float dPhi = muPhiGen[igen] - nuPhiGen[igen];
        if(dPhi<-1*TMath::Pi()) dPhi += TMath::TwoPi(); if(dPhi>TMath::Pi()) dPhi -= TMath::TwoPi(); //fold btwn [-pi,pi]
        float mtTruth = TMath::Sqrt(2.0*muGenPt[igen]*nuGenPt[igen]*(1.0-TMath::Cos(dPhi)));
        ///Generator level fiducial cuts
        if(fabs(mother[igen])==24&&centrality>0.0&&centrality<0.8) ++nW;
        if(
            fabs(mother[igen])==24
            &&fabs(daughter[igen])==13
            &&centrality>0.0&&centrality<0.8
            &&fabs(etaGen[igen])>0.1&&fabs(etaGen[igen])<2.4
            &&nuGenPt[igen]>25.0
            &&muGenPt[igen]>25.0
            &&mtTruth>40.0
            &&fabs(chargeGen[igen])>0.
        ){ 
            hGen->Fill(fabs(etaGen[igen]));
        }
    } //igen
   }//iev

   std::cout << "W+- --> mu+-: " << nW << std::endl;
   std::cout << "W+- --> mu+- (cut): " << hGen->GetEntries() << std::endl;
   TH1F* hGenc = (TH1F*)hGen->Clone("hGenc");
   // Scale by eta bw
   hGen->Scale(1.0,"width"); 
   //double sf = 1.0/(nW)*cs;
   //hGen->Scale(sf);
}



void fillMcHistos(TString fileNameIn , int bins, float xBins[], TH1F* hGenPlus, TH1F* hGenMinus, double csPlus,double csMinus, TString hname){

  float etaGen[50],chargeGen[50];
  float nuGenPt[50], nuPhiGen[50];
  int mother[50],daughter[50];
  int ngen;
  float centrality;
  float muGenPt[50], muPhiGen[50];

   TChain* tree = new TChain("tree","tree");
   tree->Add(fileNameIn);

   std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;
        
  tree->SetBranchAddress("mc_mu_gen_pt", &muGenPt);
  tree->SetBranchAddress("mc_nu_gen_pt", &nuGenPt);
  tree->SetBranchAddress("mc_nu_gen_phi", &nuPhiGen);
  tree->SetBranchAddress("mc_mu_gen_mothertype", &mother);
  tree->SetBranchAddress("mc_mu_gen_type", &daughter);
  tree->SetBranchAddress("mc_mu_charge", &chargeGen);
  tree->SetBranchAddress("mc_mu_gen_eta", &etaGen);
  tree->SetBranchAddress("mc_mu_gen_phi", &muPhiGen);
  tree->SetBranchAddress("mc_mu_n", &ngen);
  tree->SetBranchAddress("centrality", &centrality);

   // --- Set branch status ---
  tree->SetBranchStatus("*",0) ;
  tree->SetBranchStatus("mc_mu_n", 1);
  tree->SetBranchStatus("centrality", 1);
  tree->SetBranchStatus("mc_nu_gen_pt", 1);
  tree->SetBranchStatus("mc_nu_gen_phi", 1);
  tree->SetBranchStatus("mc_mu_gen_pt", 1);
  tree->SetBranchStatus("mc_mu_gen_mothertype", 1);
  tree->SetBranchStatus("mc_mu_gen_type", 1);
  tree->SetBranchStatus("mc_mu_charge", 1);
  tree->SetBranchStatus("mc_mu_gen_eta", 1);
  tree->SetBranchStatus("mc_mu_gen_phi", 1);


//   float nEvents = tree->GetEntries();
   int nWPos = 0, nWNeg=0;
   for(int iev = 0; iev < tree->GetEntries(); ++iev){
 
    tree->LoadTree(iev);
    tree->GetEntry(iev);

    for(int igen=0; igen<ngen; ++igen){

        float dPhi = muPhiGen[igen] - nuPhiGen[igen];
        if(dPhi<-1*TMath::Pi()) dPhi += TMath::TwoPi(); if(dPhi>TMath::Pi()) dPhi -= TMath::TwoPi(); //fold btwn [-pi,pi]
        float mtTruth = TMath::Sqrt(2.0*muGenPt[igen]*nuGenPt[igen]*(1.0-TMath::Cos(dPhi)));
        ///Generator level fiducial cuts
        if(mother[igen]==24&&centrality>0.0&&centrality<0.8) ++nWPos;
        else if(mother[igen]==-24&&centrality>0.0&&centrality<0.8) ++nWNeg;
        if(
            fabs(mother[igen])==24
            &&fabs(daughter[igen])==13
            &&centrality>0.0&&centrality<0.8
            &&fabs(etaGen[igen])>0.1&&fabs(etaGen[igen])<2.4
            &&nuGenPt[igen]>25.0
            &&muGenPt[igen]>25.0
            &&mtTruth>40.0
            &&fabs(chargeGen[igen])>0.
        ){ 
        ///fill mu+ eta
        if(chargeGen[igen]>0.){
            hGenPlus->Fill(fabs(etaGen[igen]));
        }
        ///fill mu- eta
        else if(chargeGen[igen]<0.)
         {
            hGenMinus->Fill(fabs(etaGen[igen]));
         }
    }
    } //igen
   }//iev

   std::cout << "W+ --> mu+: " << nWPos << std::endl;
   std::cout << "W+ --> mu+ (cut): " << hGenPlus->GetEntries() << std::endl;
   std::cout << "W- --> mu-: " << nWNeg << std::endl;
   std::cout << "W- --> mu- (cut): " << hGenMinus->GetEntries() << std::endl;
   //TH1F* hAsymm = (TH1F*)hGenPlus->Clone("hAsymm");
   TH1F* hGenPlusc = (TH1F*)hGenPlus->Clone("hGenPlusc");
   TH1F* hGenMinusc = (TH1F*)hGenMinus->Clone("hGenMinusc");
   std::cout << "PY6 Aw+: " << hGenPlus->GetEntries()/nWPos << std::endl;
   std::cout << "PY6 Aw-: " << hGenMinus->GetEntries()/nWNeg << std::endl;
   // Scale by eta bw
   hGenPlus->Scale(1.0,"width"); hGenMinus->Scale(1.0,"width");
   double sfPlus = 1.0/(nWPos)*csPlus;
   hGenPlus->Scale(sfPlus);
   double sfMinus = 1.0/(nWNeg)*csMinus;
   hGenMinus->Scale(sfMinus);
   //return hAsymm;
}

///////////////////////////////
//used for mapping negative and positive eta bins
//to absolute eta bin values
//////////////////////////////
int indexNegEta(int imap, int nEtaBins){
     int index = 0.5*(nEtaBins+1.0)-2.0; index-=imap;
     std::cout << "indexNegativeEta:" << index << std::endl;
      if(index>nEtaBins) std::cout << "WARNING: index out of allowed range. EXPECT WRONG RESULTS."<<std::endl;
     return index;
}

int indexPosEta(int imap, int nEtaBins){
      int index = imap+1.0/2.0*(nEtaBins+1.0);
      std::cout << "indexPositiveEta:" << index << std::endl;
      if(index>nEtaBins) std::cout << "WARNING: index out of allowed range. EXPECT WRONG RESULTS."<<std::endl;
      return index;
}


void Write(TFile* outFile, TGraphErrors* gr, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing TGraph to root file..." << std::endl;
	  gr->Write(sFile);
	  std::cout << "Done." << std::endl;
      gDirectory = dir;
}

void fillMcHistos(TTree* treeMC, char arrhMC[], TH1D* h, TString cutsMc){

    treeMC->Draw((TString("(abs(mc_mu_gen_eta))>>")+arrhMC).Data(),cutsMc,"pe");
    //return h;

}

void writeMcHistos(TFile* outFile, TTree* treeMC, char arrhMC[],TH1D* h){

    TDirectory *dir = gDirectory;
	outFile->cd();
    std::cout << "Number of entries in histo = " << h->GetEntries() << std::endl;
	std::cout << "Writing mc histo to root file..." << std::endl;
    h->Write(arrhMC);
	std::cout << "Done." << std::endl;
    gDirectory = dir;

}

///////////////////////////////////////////
//plotEtaDistro
//plots eta distribution separately per centrality class
//////////////////////////////////////////
void plotEtaDistro(TGraphErrors* grEta, double yW, double yWerr, double xEta, double bW, int ieta){

    grEta->SetPoint(ieta,xEta,yW/bW); 
    std::cout << "bW = " << bW << " yWerr = " << yWerr << std::endl;
    grEta->SetPointError(ieta,bW/2.0,yWerr/bW); 

}
///////////////////////////////////////////////////////////////////////////////
// sqr
///////////////////////////////////////////////////////////////////////////////
double sqr(double a) {
  return a*a;
}



int indexIJ(int i, int j, int nj){
  return i*nj+j ;
}

int indexIJK(int i, int j, int k, int nj, int nk){
  return (i*nj+j)*nk+k ;
}

// --- Must use one fit result per centrality class at a time --- //
void plotCT10DataEtaDistros() {

  bool doTheory = false;
  bool doCentrality = false;
  bool doAtlasEta = false;

//    TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_Raw.04.17.2013";
//    TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_CorrectedNoBkgSub.04.17.2013";
//    TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_CorrectedNoTrigMatchNoBkgSub.04.17.2013";
//    TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_v01.04.21.2013";
//    TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_v02.04.21.2013";
//    TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_v03.04.21.2013";
//    TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_v04.04.21.2013";
//    TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_v05.04.21.2013";
//    TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_v06.04.21.2013";
//    TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_CorrectedBkgSub_9etaBins6CentralityBins2ChargeBins.04.21.2013";

    ///main file used for result
    TString fileNameDataIn ;
    if(doAtlasEta) fileNameDataIn =
        //"ResultFiles/WAnalysis_fitResultCentChrgEta_CorrectedBkgSub_19etaBinsNoAbsEta6CentralityBins2ChargeBins.05.08.2013"; 
        //"ResultFiles/WAnalysis_fitResultCentChrgEta_RealEtaBSCorr.05.25.2013"; 
        "systematics/WAnalysis_fitResultCentChrgEta_CwAsymmetry.07.14.2013";
    else fileNameDataIn =
        //"ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr.06.03.2013";
        //"ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorrNuEtaCut.06.08.2013";
        ///Main file
        //"ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr.06.11.2013";
        //"ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr.06.26.2013";
       
        ///Nominal result
        // corrected using Py6
//        "ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr.08.26.2013";
        // corrected using Py8
        //"ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr.10.03.2013";
        // common binning w/ electrons
        "ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr.11.27.2013";
//        "ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr_noMtnoMpt.09.22.2013";
//	"ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr_Qual1.08.30.2013";
//	"ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr_Qual4.08.30.2013";
//	"ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr_Qual5.08.30.2013";
        // Corrected with Aw*Cw
//        "ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr_mu4NoTrigMatch.08.17.2013";
//        "ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr_AwCw.08.14.2013";
//        "ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr_Raw.08.13.2013";
//        "ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr_NoPreSelection.08.13.2013";
//        "ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr_QualityCuts_noDIFVeto.08.13.2013";
//        "ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr_noMPTcut.08.13.2013";
//        "ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr_NoTrigger.08.13.2013";
//          "ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr_noRecGenMatching.08.13.2013";
//        "ResultFiles/WAnalysis_fitResultCentChrgEta_sanityCheck.07.14.2013";
        //////////////////////////////////////////////////
        //Use for correlation in systematics//
        ///////////////////////////////////////////////////
        ///Nominal result used for systematic study (includes tau bkg) 
//        "systematics/WAnalysis_fitResultCentChrgEta_IncIsoConeSizeSyst.07.22.2013";
//        "systematics/WAnalysis_fitResultCentChrgEta_LoosenIsoCut.07.22.2013";
//        "systematics/WAnalysis_fitResultCentChrgEta_2GeVMpt.07.14.2013";
//        "systematics/WAnalysis_fitResultCentChrgEta_4GeVMpt.07.14.2013";
//        "systematics/WAnalysis_fitResultCentChrgEta_MptSmearSystematics.07.20.2013";
//        "systematics/WAnalysis_fitResultCentChrgEta_2GeVMpt.08.04.2013";
//        "systematics/WAnalysis_fitResultCentChrgEta_4GeVMpt.08.04.2013";
//          "systematics/WAnalysis_fitResultCentChrgEta_NominalWithTauBkg.07.30.2013";
//          "systematics/WAnalysis_fitResultCentChrgEta_MCDrivenZBkg.07.30.2013";
//          "systematics/WAnalysis_fitResultCentChrgEta_RaaScaledQCDBkg.07.30.2013";
//          "systematics/WAnalysis_fitResultCentChrgEta_TauBkgNoEmbed.07.30.2013";

        //"WAnalysis_fitResultCentChrgEta";
    ///Closure test with Wmunu MC as dataSet
//        "ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaCorrMcClosure.06.19.2013";
    ///Cw binned over ATLAS eta (i.e. not absolute)
//    TString fileNameDataIn = "ResultFiles/WAnalysis_fitResultCentChrgEta_CorrectedBkgSub_19etaBinsNoAbsEta6CentralityBins2ChargeBins.05.06.2013";
    ///Cw binned in absolute eta
    ///with charge indep cw
    //TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_CorrectedBkgSub_9etaBins6CentralityBins2ChargeBins.04.24.2013";
    ///check run periods
//    TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_periodA.04.26.2013";
//    TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_periodB.04.26.2013";
//    TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_periodC.04.26.2013";
    ///check matching
//    TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_noTriggerMatching.04.28.2013";
    ///check pt cut off
    //TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_pt10GeVUp.04.29.2013";
//    TString fileNameDataIn = "WAnalysis_fitResultCentChrgEta_pt15GeVUp.04.29.2013";

  //SetAtlasStyle();
//  TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
  TString baseString = "/usatlas/u/tbales/scratch/";
  std::cout << "Opening " << fileNameDataIn << std::endl;
  TFile* fDataSet = new TFile(fileNameDataIn+".root", "READ");
  if ( !fDataSet->IsOpen() ) {
    std::cout << fDataSet << " not found!" << std::endl;
    exit(0);
  }
  TTree *tree = (TTree*)fDataSet->Get("tree") ;

  //data overlay
//  TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.03.24.2013";
  TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.04.12.2013.root";
  TString fileNamepp = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_pp.06.21.2013.root";
  TString fileNamepn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_pn.06.21.2013.root";
  TString fileNamenp = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_np.06.21.2013.root";
  TString fileNamenn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_nn.06.21.2013.root";

  // POWHEG+PYTHIA8 MSTW2008NLO
  TString filenameMSTWpp =  baseString+"MonteCarloFiles/Wmunu/PowhegPythia8_MSTW2008NLO_pp.08.10.2013.root";
  TString filenameMSTWnn =  baseString+"MonteCarloFiles/Wmunu/PowhegPythia8_MSTW2008NLO_nn.08.10.2013.root";
  TString filenameMSTWnppn =  baseString+"MonteCarloFiles/Wmunu/PowhegPythia8_MSTW2008NLO_nppn.08.10.2013.root";

  // POWHEG+PYTHIA8 CT10
  //pp
  TString filenameCT10Pluspp =  baseString+"MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_pp.10.03.2013.root";
  TString filenameCT10Minuspp =  baseString+"MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_pp.10.03.2013.root";
  //np
  TString filenameCT10Plusnp =  baseString+"MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_np.10.03.2013.root";
  TString filenameCT10Minusnp =  baseString+"MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_np.10.03.2013.root";
  //pn
  TString filenameCT10Pluspn =  baseString+"MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_pn.10.03.2013.root";
  TString filenameCT10Minuspn =  baseString+"MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_pn.10.03.2013.root";
  //nn
  TString filenameCT10Plusnn =  baseString+"MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_nn.10.03.2013.root";
  TString filenameCT10Minusnn =  baseString+"MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_nn.10.03.2013.root";

  // Aw+- as fcn of eta
  TString filenameAw = "CorrectionFactorFiles/correctionFactorsWEta.08.11.2013.root";
  TFile* fAw = new TFile(filenameAw, "READ");

//  TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonMCWmunu.12.25.2012";
//  TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuMCDataOverlay.2012.11.21";
  TFile* fMcWSet = new TFile(fileNameMCWIn, "READ");

  if ( !fMcWSet->IsOpen() ) {
    std::cout << fMcWSet << " not found!" << std::endl;
    exit(0);
  }

  TTree *treeMC = (TTree*)fMcWSet->Get("tree") ;
  
  // Aw separate for mu+-
  const TGraphAsymmErrors* grAwPlusEta = (TGraphAsymmErrors*)fAw->Get("grAwPlusEtaDistroCent0");
  const TGraphAsymmErrors* grAwMinusEta = (TGraphAsymmErrors*)fAw->Get("grAwMinusEtaDistroCent0");
  ///TGraph holding information on how much mu+,mu- errors cancel
  // expressed as absolute on asymm
  TFile *_fSystematicDistributions = new TFile("systematicDistributions_8_22_2013.root","read");
  TGraphAsymmErrors* grIsoCancellation = (TGraphAsymmErrors*)_fSystematicDistributions->Get("grEtaDiffIso");
  TGraphAsymmErrors* grMptCancellation = (TGraphAsymmErrors*)_fSystematicDistributions->Get("grEtaDiffMpt");

  std::cout << "All files open" << std::endl;

  //yield spreadsheet names	
  TString spreadSheetNameMuPlus = "dataSpreadSheetEtaMuPlus.csv";
  TString spreadSheetNameMuMinus = "dataSpreadSheetEtaMuMinus.csv";
  TString spreadSheetNameAsymm = "dataSpreadSheetEtaAsymm.csv";

  std::ofstream spreadSheetMuPlus;
  spreadSheetMuPlus.open(spreadSheetNameMuPlus);

  std::ofstream spreadSheetMuMinus;
  spreadSheetMuMinus.open(spreadSheetNameMuMinus);

  std::ofstream spreadSheetAsymm;
  spreadSheetAsymm.open(spreadSheetNameAsymm);

	std::vector<double> etaBins;
	/*etaBins.push_back(0.10);
    etaBins.push_back(0.35);
    etaBins.push_back(0.6);
    etaBins.push_back(0.8);
    etaBins.push_back(1.05);
    etaBins.push_back(1.3);
    etaBins.push_back(1.55);
    etaBins.push_back(1.85);
    etaBins.push_back(2.1);
	etaBins.push_back(+2.40);
*/
    etaBins.push_back(0.10);
    etaBins.push_back(0.35);
    etaBins.push_back(0.6);
    etaBins.push_back(0.8);
    etaBins.push_back(1.05);
    etaBins.push_back(1.37);
    etaBins.push_back(1.52);
    etaBins.push_back(1.74);
    etaBins.push_back(2.1);
    etaBins.push_back(+2.40);

	const int nEtaBins = etaBins.size()-1;

    //float arrBinning[] = {0.1,0.35,0.6,0.8,1.05,1.3,1.55,1.85,2.1,2.40};
    float arrBinning[] = {0.1,0.35,0.6,0.8,1.05,1.37,1.52,1.74,2.1,2.40};

	std::vector<double> centralityBins;

	std::vector <float> ncoll;
    std::vector <double> npartBins;

	centralityBins.push_back(0.00);
	if (doCentrality) {

		centralityBins.push_back(0.05);
		centralityBins.push_back(0.10);
		centralityBins.push_back(0.15);
		centralityBins.push_back(0.20);
		centralityBins.push_back(0.40);

		//ncoll
		ncoll.push_back(1683.3); //0-5
		ncoll.push_back(1318.0); //5-10
		ncoll.push_back(1035.4); //10-15
		ncoll.push_back(811.2); //15-20

		ncoll.push_back(440.6); //20-40
		ncoll.push_back(77.8); //40-80
		///npart
		npartBins.push_back(382.16);//0-5
		npartBins.push_back(330.26);//5-10
		npartBins.push_back(281.88);//10-15
		npartBins.push_back(239.52);//15-20
		npartBins.push_back(157.83);//20-40
		npartBins.push_back(45.93);//40-80

	}
	else  {
		ncoll.push_back(452.0);//0-80
		npartBins.push_back(139.5);
	}
	centralityBins.push_back(0.80);

	const int nCentralityBins = centralityBins.size()-1;


//  const int nEtaBins = 9;
//  float arrBinning[] = {0.1,0.35,0.6,0.8,1.05,1.2,1.45,1.7,1.95,2.15,2.40};
  ///arrays for Salgado and eta/charge binned fits
  double y[12] = {0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75} ;
  double noWepPlus[12] = { 0.495, 0.495, 0.493, 0.491, 0.490, 0.479, 0.460, 0.427, 0.380, 0.310, 0.220, 0.120};
  double noWepPlusErr[12] = { 0.061, 0.060, 0.058, 0.051, 0.046, 0.040, 0.038, 0.031, 0.030, 0.020, 0.012, 0.010};
  double noWepMinus[12] = { 0.525, 0.526, 0.528, 0.529, 0.530, 0.527, 0.514, 0.488, 0.440, 0.370, 0.269, 0.155};
  double noWepMinusErr[12] = { 0.068, 0.065, 0.060, 0.055, 0.050, 0.045, 0.040, 0.035, 0.030, 0.025, 0.020, 0.010};
  double WepPlus[12] = {0.510, 0.508, 0.502, 0.494, 0.480, 0.458, 0.430, 0.390, 0.335, 0.260, 0.175, 0.092};
  double WepPlusErr[12] = {0.050, 0.042, 0.037, 0.036, 0.020, 0.025, 0.032, 0.038, 0.030, 0.030, 0.015, 0.010};
  double WepMinus[12] = { 0.540, 0.538, 0.535, 0.530, 0.520, 0.500, 0.478, 0.440, 0.388, 0.310, 0.210, 0.110};
  double WepMinusErr[12] = { 0.051, 0.045, 0.039, 0.029, 0.021, 0.030, 0.037, 0.040, 0.040, 0.032, 0.030, 0.015};

  // --- inclusive asymmetry --- //
  double nWplusFit[nEtaBins] ;
  double nWplusFitErrUp[nEtaBins] ;
  double nWplusFitErrDown[nEtaBins] ;
  double nWminusFit[nEtaBins] ;
  double nWminusFitErrUp[nEtaBins]  ;
  double nWminusFitErrDown[nEtaBins]  ;
 
  //systematic  
  double nWFitSystErrPlus[nEtaBins]  ; 
  double nWFitSystErrMinus[nEtaBins] ; 
  
  int markerStyle[4] = {kFullCircle, kOpenCircle, kFullTriangleUp, kFullTriangleDown };
  int markerColor[3] = {kBlack, kBlue, kRed+1};
  
/*  std::vector<double> centralityBins ;
  centralityBins.push_back(0.0) ;
  centralityBins.push_back(0.8) ;
  std::cout << "centrality " <<centralityBins.at(0)*100.0 << "-" << centralityBins.at(1)*100.0 << std::endl;
*/

  ///Get a TGraph from the fit result
  TGraphAsymmErrors* grWCentrality = (TGraphAsymmErrors*) fDataSet->Get("WPtFit_centrality") ;
  const int nWCentrality =  grWCentrality->GetN() ;
  std::cout << "number of centrality classes: " << nWCentrality << std::endl;

  std::cout << "Plotting the charge asymmetry for centrality class " << centralityBins[0] << "-" <<
    centralityBins[nCentralityBins] << std::endl;

  TString sLow = "" ;
  TString sHigh = "" ;

  sLow += 100*centralityBins[0] ; 
  sHigh += 100*centralityBins[nCentralityBins] ; 

  TString sSel = "" ;
  sSel = sLow; sSel+="-"; sSel+=sHigh; sSel+="%";



  ///pointer to array of centrality values and bin widths
  double* centralityArr = (grWCentrality->GetX()) ;
  double* centralityL   = (grWCentrality->GetEXlow()) ;
  double* centralityH   = (grWCentrality->GetEXhigh()) ;

  ///Get a TGraph that retrieves eta bins fit result
  TGraphAsymmErrors* grWEta = (TGraphAsymmErrors*) fDataSet->Get("WPtFit_eta") ;
  const int nWEta = grWEta->GetN();
  std::cout << "number of eta windows: " << nWEta << std::endl;
  double* etaArr = grWEta->GetX() ;
  double* etaL   = grWEta->GetEXlow();
  double* etaH   = grWEta->GetEXhigh();

  ///total number of TGraphs in fit result
  const int nWGraphs = nWCentrality*nWEta ;

  TObjArray arrWSig = TObjArray(nWGraphs) ;
  TObjArray arrNObs = TObjArray(nWGraphs) ;
  TObjArray arrNBkg = TObjArray(nWGraphs) ;
  TObjArray arrWSigSyst = TObjArray(nWGraphs);
  TObjArray arrWSigSyst1 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst2 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst3 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst4 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst5 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst6 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst7 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst8 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst9 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst10 = TObjArray(nWGraphs);

  TString baseW = "WPtFit" ;
  TString baseWSig = baseW ; baseWSig += "_sig_eta" ; ///stat errors
  TString baseNObs = baseW ; baseNObs += "_nobs_eta" ; ///stat errors
  TString baseNBkg = baseW ; baseNBkg += "_nbkg_eta" ; ///stat errors
  TString baseWSigSyst = baseW;   baseWSigSyst +="_sigSyst_eta"; ///syst. errors
  TString baseWSigSyst1 = baseW;   baseWSigSyst1 +="_sigSyst1_eta"; ///syst. errors
  TString baseWSigSyst2 = baseW;   baseWSigSyst2 +="_sigSyst2_eta"; ///syst. errors
  TString baseWSigSyst3 = baseW;   baseWSigSyst3 +="_sigSyst3_eta"; ///syst. errors
  TString baseWSigSyst4 = baseW;   baseWSigSyst4 +="_sigSyst4_eta"; ///syst. errors
  TString baseWSigSyst5 = baseW;   baseWSigSyst5 +="_sigSyst5_eta"; ///syst. errors
  TString baseWSigSyst6 = baseW;   baseWSigSyst6 +="_sigSyst6_eta"; ///syst. errors
  TString baseWSigSyst7 = baseW;   baseWSigSyst7 +="_sigSyst7_eta"; ///syst. errors
  TString baseWSigSyst8 = baseW;   baseWSigSyst8 +="_sigSyst8_eta"; ///syst. errors
  TString baseWSigSyst9 = baseW;   baseWSigSyst9 +="_sigSyst9_eta"; ///syst. errors
  TString baseWSigSyst10 = baseW;   baseWSigSyst10 +="_sigSyst10_eta"; ///syst. errors

  ///loop over TGraphs in fit result
  for (int i=0; i<nWEta; i++){
    TString searchWI = baseWSig;  searchWI += i;  searchWI += "_centrality";
    TString searchNObsI = baseNObs;  searchNObsI += i;  searchNObsI += "_centrality";
    TString searchNBkgI = baseNBkg;  searchNBkgI += i;  searchNBkgI += "_centrality";
    TString searchSystWI = baseWSigSyst; searchSystWI += i; searchSystWI += "_centrality";
    TString searchSyst1WI = baseWSigSyst1; searchSyst1WI += i; searchSyst1WI += "_centrality";
    TString searchSyst2WI = baseWSigSyst2; searchSyst2WI += i; searchSyst2WI += "_centrality";
    TString searchSyst3WI = baseWSigSyst3; searchSyst3WI += i; searchSyst3WI += "_centrality";
    TString searchSyst4WI = baseWSigSyst4; searchSyst4WI += i; searchSyst4WI += "_centrality";
    TString searchSyst5WI = baseWSigSyst5; searchSyst5WI += i; searchSyst5WI += "_centrality";
    TString searchSyst6WI = baseWSigSyst6; searchSyst6WI += i; searchSyst6WI += "_centrality";
    TString searchSyst7WI = baseWSigSyst7; searchSyst7WI += i; searchSyst7WI += "_centrality";
    TString searchSyst8WI = baseWSigSyst8; searchSyst8WI += i; searchSyst8WI += "_centrality";
    TString searchSyst9WI = baseWSigSyst9; searchSyst9WI += i; searchSyst9WI += "_centrality";
    TString searchSyst10WI = baseWSigSyst10; searchSyst10WI += i; searchSyst10WI += "_centrality";

    for (int j=0; j<nWCentrality; j++){
      TString searchWIJ = searchWI;  searchWIJ += j;
      TString searchNObsIJ = searchNObsI;  searchNObsIJ += j;
      TString searchNBkgIJ = searchNBkgI;  searchNBkgIJ += j;
      TString searchSystWIJ = searchSystWI; searchSystWIJ += j;
      TString searchSyst1WIJ = searchSyst1WI; searchSyst1WIJ += j;
      TString searchSyst2WIJ = searchSyst2WI; searchSyst2WIJ += j;
      TString searchSyst3WIJ = searchSyst3WI; searchSyst3WIJ += j;
      TString searchSyst4WIJ = searchSyst4WI; searchSyst4WIJ += j;
      TString searchSyst5WIJ = searchSyst5WI; searchSyst5WIJ += j;
      TString searchSyst6WIJ = searchSyst6WI; searchSyst6WIJ += j;
      TString searchSyst7WIJ = searchSyst7WI; searchSyst7WIJ += j;
      TString searchSyst8WIJ = searchSyst8WI; searchSyst8WIJ += j;
      TString searchSyst9WIJ = searchSyst9WI; searchSyst9WIJ += j;
      TString searchSyst10WIJ = searchSyst10WI; searchSyst10WIJ += j;

      ///give each TGraph a unique index and fill an array
      const int index = indexIJ(i,j,nWCentrality) ;
      arrWSig[index] = ((TGraphAsymmErrors*) fDataSet->Get(searchWIJ)) ;
      arrNObs[index] = ((TGraphAsymmErrors*) fDataSet->Get(searchNObsIJ)) ;
      arrNBkg[index] = ((TGraphAsymmErrors*) fDataSet->Get(searchNBkgIJ)) ;
      arrWSigSyst[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSystWIJ));
      arrWSigSyst1[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst1WIJ));
      arrWSigSyst2[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst2WIJ));
      arrWSigSyst3[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst3WIJ));
      arrWSigSyst4[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst4WIJ));
      arrWSigSyst5[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst5WIJ));
      arrWSigSyst6[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst6WIJ));
      arrWSigSyst7[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst7WIJ));
      arrWSigSyst8[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst8WIJ));
      arrWSigSyst9[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst9WIJ));
      arrWSigSyst10[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst10WIJ));

      int nWpts = ((TGraphAsymmErrors*) arrWSig[index])->GetN() ;
      for (int k =0; k<nWpts; k++){
        ///pointers to array of y-values in TGraph indexed
        double* ytemp = ((TGraphAsymmErrors*)arrWSig[index])->GetY() ;
        double* ytempNObs = ((TGraphAsymmErrors*)arrNObs[index])->GetY() ;
        double* ytempNBkg = ((TGraphAsymmErrors*)arrNBkg[index])->GetY() ;
        double* ytempL = ((TGraphAsymmErrors*)arrWSig[index])->GetEYlow() ;
        double* ytempH = ((TGraphAsymmErrors*)arrWSig[index])->GetEYhigh() ;
        double* ytempSystL = ((TGraphAsymmErrors*)arrWSigSyst[index])->GetEYlow();
        double* ytempSystH = ((TGraphAsymmErrors*)arrWSigSyst[index])->GetEYhigh();
        double* ytempSyst1 = ((TGraphAsymmErrors*)arrWSigSyst1[index])->GetEYhigh();
        double* ytempSyst2 = ((TGraphAsymmErrors*)arrWSigSyst2[index])->GetEYhigh();
        double* ytempSyst3 = ((TGraphAsymmErrors*)arrWSigSyst3[index])->GetEYhigh();
        double* ytempSyst4 = ((TGraphAsymmErrors*)arrWSigSyst4[index])->GetEYhigh();
        double* ytempSyst5 = ((TGraphAsymmErrors*)arrWSigSyst5[index])->GetEYhigh();
        double* ytempSyst6 = ((TGraphAsymmErrors*)arrWSigSyst6[index])->GetEYhigh();
        double* ytempSyst7 = ((TGraphAsymmErrors*)arrWSigSyst7[index])->GetEYhigh();
        double* ytempSyst8 = ((TGraphAsymmErrors*)arrWSigSyst8[index])->GetEYhigh();
        double* ytempSyst9 = ((TGraphAsymmErrors*)arrWSigSyst9[index])->GetEYhigh();
        double* ytempSyst10 = ((TGraphAsymmErrors*)arrWSigSyst10[index])->GetEYhigh();

        ///Find which Y coordinate corresponds to the nW's in the given eta and cent bin
        ///Only the ytemp matching the charge index in the fit result 
        ///for W^{+},W^{-} will be used (i.e. 102, 103)
      if (ytemp[k]>0&&(abs(ytempL[k])>0||ytempH[k]>0)){
	    
/*        std::cout << " found " << searchWIJ  << "\n" << searchSystWIJ << "\n" 
                    << searchSyst1WIJ << "\n" 
                    << searchSyst2WIJ << "\n" 
                    << searchSyst3WIJ << "\n" 
                    << searchSyst4WIJ << "\n" 
                    << searchSyst5WIJ << "\n" 
                    << searchSyst6WIJ << "\n" 
                    << searchSyst7WIJ << "\n" 
                    << searchSyst8WIJ << "\n" 
                    << searchSyst9WIJ << "\n" 
                    << searchSyst10WIJ << "\n" 
                    << std::endl;
*/
    //	  	    << " " << ytemp[k] << " +"<< ytempH[k]  << " -"<< abs(ytempL[k]) << " (stat.)" << " +" 
   //         << ytempSystH[k]  << " -" << abs(ytempSystL[k]) << " (syst.)" << std::endl;
        }   
      }  //nWpt
    }  //nWCentrality
  }  //nWEta

  double* etaLow = new double[nWEta] ;
  double* etaHigh = new double[nWEta] ;
  double* centralityLow = new double[nWCentrality];
  double* centralityHigh = new double[nWCentrality] ;
  
  ///arrays to hold eta/centrality/charge nW data
  const int nWpt = ((TGraphAsymmErrors*)arrWSig[0])->GetN();
  std::cout << "nWpt = " << nWpt << std::endl;
  //int nWpt=2;
  double* yW = new double[nWEta*nWCentrality*nWpt] ;
  double* eyWL = new double[nWEta*nWCentrality*nWpt] ;
  double* eyWH = new double[nWEta*nWCentrality*nWpt] ;
  double* eyWSystL   = new double[nWEta*nWCentrality*nWpt];
  double* eyWSystH   = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst1 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst2 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst3 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst4 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst5 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst6 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst7 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst8 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst9 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst10 = new double[nWEta*nWCentrality*nWpt];

  double* yWPlus = new double[nWEta*nWCentrality*nWpt] ;
  double* yWMinus = new double[nWEta*nWCentrality*nWpt] ;
  double* eyWPlusL= new double[nWEta]; 
  double* eyWPlusH= new double[nWEta]; 
  double* eyWPlusSystL= new double[nWEta]; 
  double* eyWPlusSystH= new double[nWEta]; 
  double* eyWPlusSyst1= new double[nWEta]; 
  double* eyWPlusSyst2= new double[nWEta]; 
  double* eyWPlusSyst3= new double[nWEta]; 
  double* eyWPlusSyst4= new double[nWEta]; 
  double* eyWPlusSyst5= new double[nWEta]; 
  double* eyWPlusSyst6= new double[nWEta]; 
  double* eyWPlusSyst7= new double[nWEta]; 
  double* eyWPlusSyst8= new double[nWEta]; 
  double* eyWPlusSyst9= new double[nWEta]; 
  double* eyWPlusSyst10= new double[nWEta]; 

  double* eyWMinusL= new double[nWEta]; 
  double* eyWMinusH= new double[nWEta]; 
  double* eyWMinusSystL= new double[nWEta]; 
  double* eyWMinusSystH= new double[nWEta]; 
  double* eyWMinusSyst1= new double[nWEta]; 
  double* eyWMinusSyst2= new double[nWEta]; 
  double* eyWMinusSyst3= new double[nWEta]; 
  double* eyWMinusSyst4= new double[nWEta]; 
  double* eyWMinusSyst5= new double[nWEta]; 
  double* eyWMinusSyst6= new double[nWEta]; 
  double* eyWMinusSyst7= new double[nWEta]; 
  double* eyWMinusSyst8= new double[nWEta]; 
  double* eyWMinusSyst9= new double[nWEta]; 
  double* eyWMinusSyst10= new double[nWEta]; 

  double* eyWPlusSumL= new double[nWEta]; 
  double* eyWPlusSumH= new double[nWEta]; 
  double* eyWPlusSumSystL = new double[nWEta];
  double* eyWPlusSumSystH = new double[nWEta];
  double* eyWPlusSumSyst1 = new double[nWEta];
  double* eyWPlusSumSyst2 = new double[nWEta];
  double* eyWPlusSumSyst3 = new double[nWEta];
  double* eyWPlusSumSyst4 = new double[nWEta];
  double* eyWPlusSumSyst5 = new double[nWEta];
  double* eyWPlusSumSyst6 = new double[nWEta];
  double* eyWPlusSumSyst7 = new double[nWEta];
  double* eyWPlusSumSyst8 = new double[nWEta];
  double* eyWPlusSumSyst9 = new double[nWEta];
  double* eyWPlusSumSyst10 = new double[nWEta];

  double* eyWMinusSumL= new double[nWEta]; 
  double* eyWMinusSumH= new double[nWEta]; 
  double* eyWMinusSumSystL = new double[nWEta];
  double* eyWMinusSumSystH = new double[nWEta];
  double* eyWMinusSumSyst1 = new double[nWEta];
  double* eyWMinusSumSyst2 = new double[nWEta];
  double* eyWMinusSumSyst3 = new double[nWEta];
  double* eyWMinusSumSyst4 = new double[nWEta];
  double* eyWMinusSumSyst5 = new double[nWEta];
  double* eyWMinusSumSyst6 = new double[nWEta];
  double* eyWMinusSumSyst7 = new double[nWEta];
  double* eyWMinusSumSyst8 = new double[nWEta];
  double* eyWMinusSumSyst9 = new double[nWEta];
  double* eyWMinusSumSyst10 = new double[nWEta];

  double* yWPlusTot= new double[nWEta]; 
  double* yNObsPlusTot= new double[nWEta]; 
  double* yNBkgPlusTot= new double[nWEta]; 
  double* yWMinusTot= new double[nWEta];
  double* yNObsMinusTot= new double[nWEta]; 
  double* yNBkgMinusTot= new double[nWEta]; 
  double* yWPlusAbsEtaSum = new double[nEtaBins];
  double* yNObsPlusAbsEtaSum = new double[nEtaBins];
  double* yNBkgPlusAbsEtaSum = new double[nEtaBins];
  double* yWMinusAbsEtaSum = new double[nEtaBins];
  double* yNObsMinusAbsEtaSum = new double[nEtaBins];
  double* yNBkgMinusAbsEtaSum = new double[nEtaBins];
  double* eyWPlusAbsEtaSum = new double[nEtaBins];
  double* eyWMinusAbsEtaSum = new double[nEtaBins];
  double* eyUncorrelatedSystPlusAbsEtaSum = new double[nEtaBins];
  double* eyUncorrelatedSystMinusAbsEtaSum = new double[nEtaBins];
  double* eyCorrelatedSystPlusAbsEtaSum = new double[nEtaBins];
  double* eyCorrelatedSystMinusAbsEtaSum = new double[nEtaBins];
  double* eySyst1WPlusAbsEtaSum = new double[nEtaBins];
  double* eySyst1WMinusAbsEtaSum = new double[nEtaBins];
  double* eySyst2WPlusAbsEtaSum = new double[nEtaBins];
  double* eySyst2WMinusAbsEtaSum = new double[nEtaBins];
  double* eySyst3WPlusAbsEtaSum = new double[nEtaBins];
  double* eySyst3WMinusAbsEtaSum = new double[nEtaBins];
  double* eySyst4WPlusAbsEtaSum = new double[nEtaBins];
  double* eySyst4WMinusAbsEtaSum = new double[nEtaBins];
  double* eySyst5WPlusAbsEtaSum = new double[nEtaBins];
  double* eySyst5WMinusAbsEtaSum = new double[nEtaBins];
  double* eySyst6WPlusAbsEtaSum = new double[nEtaBins];
  double* eySyst6WMinusAbsEtaSum = new double[nEtaBins];
  double* eySyst7WPlusAbsEtaSum = new double[nEtaBins];
  double* eySyst7WMinusAbsEtaSum = new double[nEtaBins];
  double* eySyst8WPlusAbsEtaSum = new double[nEtaBins];
  double* eySyst8WMinusAbsEtaSum = new double[nEtaBins];
  double* eySyst9WPlusAbsEtaSum = new double[nEtaBins];
  double* eySyst9WMinusAbsEtaSum = new double[nEtaBins];
  double* eySyst10WPlusAbsEtaSum = new double[nEtaBins];
  double* eySyst10WMinusAbsEtaSum = new double[nEtaBins];

  //const int characters = 100;
  char* hpName = new char[100];
  char* hnName = new char[100]; 
  char* hpNameMC = new char[100];
  char* hnNameMC = new char[100]; 
  char* hassName = new char[100];
  const int ncent = nWCentrality;

  TH1D* hmAsymm[ncent] ;
  TH1D* hpAsymm[ncent] ;
  TH1D* hmAsymmMC[ncent] ;
  TH1D* hpAsymmMC[ncent] ;
  TH1D* hAsymm[0] ;

  int nassBins = nWEta;
//  double etaPts[] = {0.0, 0.50, 1.0, 1.50, 2.0, 2.5} ;

//  for (int i=0; i<1; i++){
        sprintf(hpName,"hpFit_cent%i",0);
        sprintf(hnName,"hnFit_cent%i",0);
        sprintf(hpNameMC,"hpFit_centMC%i",0);
        sprintf(hnNameMC,"hnFit_centMC%i",0);
        sprintf(hassName,"hassFit_cent%i",0);

        hpAsymm[0] = new TH1D(hpName,hpName,nEtaBins,arrBinning);
        hmAsymm[0] = new TH1D(hnName,hnName,nEtaBins,arrBinning);
        hpAsymmMC[0] = new TH1D(hpNameMC,hpNameMC,nEtaBins,arrBinning);
        hmAsymmMC[0] = new TH1D(hnNameMC,hnNameMC,nEtaBins,arrBinning);
        //hAsymm[0] = new TH1D(hassName,hassName,nEtaBins,arrBinning);
        hAsymm[0] = new TH1D(hassName,hassName,nEtaBins,arrBinning);
//   }      


 TString fileNameDataOut = "signalChargeEtaDistributions";
 TFile *outFile = new TFile(fileNameDataOut+".root","RECREATE");
 TFile *fHisto = new TFile("muonAsymmetry.root","RECREATE");

 TList _grEtaDistroPlusCent;
 TList _grEtaDistroMinusCent;

 TH1D* arrhpD[3];
 TH1D* arrhmD[3];
 TH1D* arrhAsymmD[3];
 TH1D* arrhpMC[nCentralityBins];
 TH1D* arrhmMC[nCentralityBins];

 for (int i=0; i<nWCentrality; i++){
    centralityLow[i] = centralityArr[i] - centralityL[i] ;
    centralityHigh[i] = centralityArr[i] + centralityH[i] ;
    std::cout << "CentralityLow = " << centralityArr[i] << "- " << centralityL[i] << std::endl;
    std::cout << "CentralityHigh = " << centralityArr[i] << "+ " << centralityH[i] << std::endl;

    //_grEtaDistroPlusCent.Add( new TH1D(nWEta));
    //_grEtaDistroMinusCent.Add( new TH1D(nWEta));


  }

double* yWPlusTotCent0_15 = new double[nEtaBins];
double* yWPlusTotCent15_40 = new double[nEtaBins];
double* yWPlusTotCent40_80 = new double[nEtaBins];

double* yWMinusTotCent0_15 = new double[nEtaBins];
double* yWMinusTotCent15_40 = new double[nEtaBins];
double* yWMinusTotCent40_80 = new double[nEtaBins];

if(doCentrality) {
  for(int i=0; i<3; ++i){
    char cHistoPlus[50];
    sprintf(cHistoPlus,"dataEtaDistroPlusCent%i",i);

    char cHistoMinus[50];
    sprintf(cHistoMinus,"dataEtaDistroMinusCent%i",i);

    char cHistoAsymm[50];
    sprintf(cHistoAsymm,"dataChargeAsymmetryCent%i",i);
    arrhpD[i] = new TH1D(cHistoPlus,cHistoPlus,nEtaBins,arrBinning);
    arrhmD[i] = new TH1D(cHistoMinus,cHistoMinus,nEtaBins,arrBinning);
    arrhAsymmD[i] = new TH1D(cHistoAsymm,cHistoAsymm,nEtaBins,arrBinning);
   }
}
 for (int i=0; i<nWEta; i++){
    etaLow[i] = etaArr[i] - etaL[i] ;
    etaHigh[i] = etaArr[i] + etaH[i] ;
    std::cout << "EtaLow = " << etaArr[i] << "- " << etaL[i] << std::endl;
    std::cout << "EtaHigh = " << etaArr[i] << "+ " << etaH[i] << std::endl;
  }


  ///build the data arrays
  ///by running over all event centralities
  ///for a given eta bin i
  std::cout << "Building data arrays" << std::endl;
  for (int i = 0; i<nWEta; i++){
     //reset running total for each eta bin
     yWPlusTot[i] = 0.;
     yNObsPlusTot[i] = 0.;
     yNBkgPlusTot[i] = 0.;
     yWMinusTot[i] = 0.;
     yNObsMinusTot[i] = 0.;
     yNBkgMinusTot[i] = 0.;
     yWPlusTotCent0_15[i] = 0.;
     yWPlusTotCent15_40[i] = 0.;
     yWPlusTotCent40_80[i] = 0.;
     yWMinusTotCent0_15[i] = 0.;
     yWMinusTotCent15_40[i] = 0.;
     yWMinusTotCent40_80[i] = 0.;
     eyWPlusSumL[i] = 0.;
     eyWPlusSumH[i] = 0.;
     eyWMinusSumH[i] = 0.;
     eyWMinusSumL[i] = 0.;
     eyWPlusSumSystL[i] = 0.;
     eyWPlusSumSystH[i] = 0.;
     eyWMinusSumSystH[i] = 0.;
     eyWMinusSumSystL[i] = 0.;
     eyWPlusSumSyst1[i] = 0.;
     eyWPlusSumSyst2[i] = 0.;
     eyWPlusSumSyst3[i] = 0.;
     eyWPlusSumSyst4[i] = 0.;
     eyWPlusSumSyst5[i] = 0.;
     eyWPlusSumSyst6[i] = 0.;
     eyWPlusSumSyst7[i] = 0.;
     eyWPlusSumSyst8[i] = 0.;
     eyWPlusSumSyst9[i] = 0.;
     eyWPlusSumSyst10[i] = 0.;

     eyWMinusSumSyst1[i] = 0.;
     eyWMinusSumSyst2[i] = 0.;
     eyWMinusSumSyst3[i] = 0.;
     eyWMinusSumSyst4[i] = 0.;
     eyWMinusSumSyst5[i] = 0.;
     eyWMinusSumSyst6[i] = 0.;
     eyWMinusSumSyst7[i] = 0.;
     eyWMinusSumSyst8[i] = 0.;
     eyWMinusSumSyst9[i] = 0.;
     eyWMinusSumSyst10[i] = 0.;
 
     
     float xEtaTemp = (etaBins[i+1]-etaBins[i])/2.0 + etaBins[i]; 
     float binWidthTemp = etaBins[i+1]-etaBins[i];
     //loop over each centrality class and add up the
     //number of Ws in this eta bin
     for (int j =0; j<nWCentrality; j++){

       const int index = indexIJ(i,j,nWCentrality) ;
       TString sGraphNamePlus = "sigEtaDistroPlusCent"; sGraphNamePlus+=j;
       TString sGraphNameMinus = "sigEtaDistroMinusCent"; sGraphNameMinus+=j;

       //find the TGrph corresponding to this eta and centrality class
       if (arrWSig[index]) {

        double* yW0 ;
        double* yW0L ;
        double* yW0H ;

	    ///get the number of Ws from this eta/centrality class
        ///this is where we would like to add negative eta bins into absolute 
        ///bins
        yW0 = ((TGraphAsymmErrors*)arrWSig[index])->GetY();
        yW0H = ((TGraphAsymmErrors*)arrWSig[index])->GetEYhigh();
        yW0L = ((TGraphAsymmErrors*)arrWSig[index])->GetEYlow();
	    double* yNObs = ((TGraphAsymmErrors*)arrNObs[index])->GetY();
	    double* yNBkg = ((TGraphAsymmErrors*)arrNBkg[index])->GetY();
	    double* yW0SystL = ((TGraphAsymmErrors*)arrWSigSyst[index])->GetEYlow();
	    double* yW0SystH = ((TGraphAsymmErrors*)arrWSigSyst[index])->GetEYhigh();
	    double* yW0Syst1 = ((TGraphAsymmErrors*)arrWSigSyst1[index])->GetEYhigh();
	    double* yW0Syst2 = ((TGraphAsymmErrors*)arrWSigSyst2[index])->GetEYhigh();
	    double* yW0Syst3 = ((TGraphAsymmErrors*)arrWSigSyst3[index])->GetEYhigh();
	    double* yW0Syst4 = ((TGraphAsymmErrors*)arrWSigSyst4[index])->GetEYhigh();
	    double* yW0Syst5 = ((TGraphAsymmErrors*)arrWSigSyst5[index])->GetEYhigh();
	    double* yW0Syst6 = ((TGraphAsymmErrors*)arrWSigSyst6[index])->GetEYhigh();
	    double* yW0Syst7 = ((TGraphAsymmErrors*)arrWSigSyst7[index])->GetEYhigh();
	    double* yW0Syst8 = ((TGraphAsymmErrors*)arrWSigSyst8[index])->GetEYhigh();
	    double* yW0Syst9 = ((TGraphAsymmErrors*)arrWSigSyst9[index])->GetEYhigh();
	    double* yW0Syst10 = ((TGraphAsymmErrors*)arrWSigSyst10[index])->GetEYhigh();

        int nPts = ((TGraphAsymmErrors*) arrWSig[index])->GetN() ;
        for (int k = 0; k<nPts; k++) {

           const int index2 = indexIJK(i,j,k,nWCentrality,nPts) ;

          if (yW0[k]>0 && (abs(yW0L[k])>0||yW0H[k]>0)){

	       ///mu+
	       if(k == 102){

	         ///number of W+ in this eta/centrality slice
             //yWPlus[index2] = yW0[k];
	         double yErrTempL = yW0L[k];
	         double yErrTempH = yW0H[k];
             double yErrSystTemp = yW0SystL[k];
             double errStatSystTemp = TMath::Sqrt( TMath::Power(yErrSystTemp,2) + TMath::Power(yErrTempL,2)  );

	         ///running total for NW+ in this eta bin	
	         yWPlusTot[i] += yW0[k];
		 yNObsPlusTot[i] +=yNObs[k];
		 yNBkgPlusTot[i] +=yNBkg[k];
             eyWPlusSumL[i] += TMath::Power(yW0L[k],2);
             eyWPlusSumH[i] += TMath::Power(yW0H[k],2);
             eyWPlusSumSystL[i] += TMath::Power(yW0SystL[k],2);
             eyWPlusSumSystH[i] += TMath::Power(yW0SystH[k],2);
             eyWPlusSumSyst1[i] += TMath::Power(yW0Syst1[k],2);
             eyWPlusSumSyst2[i] += TMath::Power(yW0Syst2[k],2);
             eyWPlusSumSyst3[i] += TMath::Power(yW0Syst3[k],2);
             eyWPlusSumSyst4[i] += TMath::Power(yW0Syst4[k],2);
             eyWPlusSumSyst5[i] += TMath::Power(yW0Syst5[k],2);
             eyWPlusSumSyst6[i] += TMath::Power(yW0Syst6[k],2);
             eyWPlusSumSyst7[i] += TMath::Power(yW0Syst7[k],2);
             eyWPlusSumSyst8[i] += TMath::Power(yW0Syst8[k],2);
             eyWPlusSumSyst9[i] += TMath::Power(yW0Syst9[k],2);
             eyWPlusSumSyst10[i] += TMath::Power(yW0Syst10[k],2);

             std::cout << i << ":" << j << ":" << k << ": " << " W^{+} = " << yW0[k] << "+" << yW0H[k] << " " << yW0L[k]
                << " stat. " << std::endl; 
             //   << " +"<<yW0SystH[k]<<"-"<<yW0SystL[k]  << " syst."<< std::endl;
	         std::cout << "Running total W^{+} for eta bin " << i << " = " << yWPlusTot[i] << std::endl;

             ///now that we have the number of W- in this given eta and centrality bin, let's plot the eta distribution
             ///for this particular centrality class
            if(doCentrality){
              if(j>=0&&j<3) {
                  yWPlusTotCent0_15[i] += yW0[k];

              }
              if(j>=3&&j<5) {
                  yWPlusTotCent15_40[i] += yW0[k];
              }
              if(j==5) yWPlusTotCent40_80[i] += yW0[k];

              //plotEtaDistro((TGraphErrors*)_grEtaDistroPlusCent.At(j),yWPlus[index2],yErrTempL,xEtaTemp,binWidthTemp,i);
              if(j==2){
                arrhpD[0]->SetBinContent(i+1,yWPlusTotCent0_15[i]/binWidthTemp); 
                arrhpD[0]->SetBinError(i+1,TMath::Sqrt(yWPlusTotCent0_15[i])/binWidthTemp);

                std::cout << "W+ in 0-15%:" << yWPlusTotCent0_15[i] << std::endl;
                std::cout << "W+ in 0-15 divided by bw:" << yWPlusTotCent0_15[i]/binWidthTemp << std::endl;
              }
              if(j==4){
                arrhpD[1]->SetBinContent(i+1,yWPlusTotCent15_40[i]/binWidthTemp); 
                arrhpD[1]->SetBinError(i+1,TMath::Sqrt(yWPlusTotCent15_40[i])/binWidthTemp);
              }
              if(j==5){
                arrhpD[2]->SetBinContent(i+1,yWPlusTotCent40_80[i]/binWidthTemp); 
                arrhpD[2]->SetBinError(i+1,TMath::Sqrt(yWPlusTotCent40_80[i])/binWidthTemp);
              }

              if(i==nWEta-1)  {

                //Write(outFile, (TH1D*)_grEtaDistroPlusCent.At(j),sGraphNamePlus); 

                double centLow = centralityBins[j]; double centUpp = centralityBins[j+1];
                TString cutsMc = "mc_mu_gen_mothertype==+24";  cutsMc+="&&mc_mu_charge==+1";
                    cutsMc+="&&mc_nu_gen_pt>25.&&mc_mu_gen_pt>25.0"; cutsMc+="&&abs(mc_mu_gen_eta)>";cutsMc+=etaBins[0];
                    cutsMc+="&&abs(mc_mu_gen_eta)<"; cutsMc+=etaBins[nEtaBins]; cutsMc+="&&centrality>";
                    cutsMc+=centLow; cutsMc+="&&centrality<"; cutsMc+=centUpp;
                    cutsMc+="&&TMath::Sqrt(2.0*mc_mu_gen_pt*mc_nu_gen_pt*(1.0-TMath::Cos(mc_mu_gen_phi-mc_nu_gen_phi)))>40.0";
                    std::cout << "cutsMc Plus " << cutsMc << std::endl;

                char cMcHistoPlus[50];
                sprintf(cMcHistoPlus,"mcTruthEtaDistroPlusCent%i",j);
                arrhpMC[j] = new TH1D(cMcHistoPlus,cMcHistoPlus,nEtaBins,arrBinning);
                fillMcHistos(treeMC, cMcHistoPlus,arrhpMC[j], cutsMc);
                //double norm = arrhpD[j]->Integral()/arrhpMC[j]->Integral();
                //arrhpMC[j]->Scale(norm); 
                arrhpMC[j]->Scale(1.0,"width");
                writeMcHistos(outFile, treeMC,cMcHistoPlus, arrhpMC[j]);

                //TH1D* arrhpDc = (TH1D*)arrhpD[j]->Clone(cHistoPlus);
                //arrhpD[j]->Scale(1.0,"width");
                //arrhpDc->Scale(1.0,"width");

                outFile->cd();
                if(j==2){
                    arrhpD[0]->Write();
                }
                if(j==4){
                    arrhpD[1]->Write();
                }
                if(j==5)arrhpD[2]->Write();

            }
          }
 	    } //mu+ index

	    ///mu-
	    if(k == 103){

         //yWMinus[index2] = yW0[k];
	      yWMinusTot[i] += yW0[k];
	      yNObsMinusTot[i] +=yNObs[k];
	      yNBkgMinusTot[i] +=yNBkg[k];
	      double yErrTempL = yW0L[k];
	      double yErrTempH = yW0H[k];
              double yErrSystTemp = yW0SystL[k];
              double errStatSystTemp = TMath::Sqrt( TMath::Power(yErrSystTemp,2) + TMath::Power(yErrTempL,2)  );

		 eyWMinusSumL[i] += TMath::Power(yErrTempL,2);
		 eyWMinusSumH[i] += TMath::Power(yErrTempH,2);
		 eyWMinusSumSystL[i] += TMath::Power(yW0SystL[k],2);
		 eyWMinusSumSystH[i] += TMath::Power(yW0SystH[k],2);
		 eyWMinusSumSystL[i] += TMath::Power(yW0SystL[k],2);
		 eyWMinusSumSystH[i] += TMath::Power(yW0SystH[k],2);
		 eyWMinusSumSyst1[i] += TMath::Power(yW0Syst1[k],2);
		 eyWMinusSumSyst2[i] += TMath::Power(yW0Syst2[k],2);
		 eyWMinusSumSyst3[i] += TMath::Power(yW0Syst3[k],2);
		 eyWMinusSumSyst4[i] += TMath::Power(yW0Syst4[k],2);
		 eyWMinusSumSyst5[i] += TMath::Power(yW0Syst5[k],2);
		 eyWMinusSumSyst6[i] += TMath::Power(yW0Syst6[k],2);
		 eyWMinusSumSyst7[i] += TMath::Power(yW0Syst7[k],2);
		 eyWMinusSumSyst8[i] += TMath::Power(yW0Syst8[k],2);
		 eyWMinusSumSyst9[i] += TMath::Power(yW0Syst9[k],2);
		 eyWMinusSumSyst10[i] += TMath::Power(yW0Syst10[k],2);

		 std::cout << i << ":" << j << ":" << k << ": " << " W^{-} = " << yW0[k] << "+" << yW0H[k] << " " << yW0L[k] <<
			" stat. " <<std::endl; 
		 //       << " +"<<yW0SystH[k]<<"-"<<yW0SystL[k]  << " syst." << std::endl;
		     std::cout << "Running total W^{-} for eta bin " << i << " = " << yWMinusTot[i] << std::endl;


        ///now that we have the number of W- in this given eta and centrality bin, let's plot the eta distribution
        ///for this particular centrality class
        if(doCentrality) {

            //plotEtaDistro((TGraphErrors*)_grEtaDistroMinusCent.At(j),yWMinus[index2],yErrTempH,xEtaTemp,binWidthTemp,i);

              if(j>=0&&j<3) {
                  yWMinusTotCent0_15[i] += yW0[k];
              }
              else if(j>=3&&j<5) yWMinusTotCent15_40[i] += yW0[k];
              else yWMinusTotCent40_80[i] += yW0[k];

              //plotEtaDistro((TGraphErrors*)_grEtaDistroMinusCent.At(j),yWMinus[index2],yErrTempL,xEtaTemp,binWidthTemp,i);
              if(j==2){
                arrhmD[0]->SetBinContent(i+1,yWMinusTotCent0_15[i]/binWidthTemp); 
                arrhmD[0]->SetBinError(i+1,TMath::Sqrt(yWMinusTotCent0_15[i])/binWidthTemp);

                double nSigPlusTemp = arrhpD[0]->GetBinContent(i+1)*binWidthTemp;
                double nSigMinusTemp = arrhmD[0]->GetBinContent(i+1)*binWidthTemp;
                double asymmTemp = ((nSigPlusTemp-nSigMinusTemp)/(nSigPlusTemp+nSigMinusTemp));
                double err = TMath::Sqrt(nSigPlusTemp+nSigMinusTemp);
                double errAsymm = TMath::Sqrt(TMath::Power(err/(nSigPlusTemp-nSigMinusTemp),2)+TMath::Power(err/(nSigPlusTemp+nSigMinusTemp),2))*asymmTemp;
                arrhAsymmD[0]->SetBinContent(i+1,asymmTemp);
                arrhAsymmD[0]->SetBinError(i+1,errAsymm);
              }
              if(j==4){
                
                arrhmD[1]->SetBinContent(i+1,yWMinusTotCent15_40[i]/binWidthTemp); 
                arrhmD[1]->SetBinError(i+1,TMath::Sqrt(yWMinusTotCent15_40[i])/binWidthTemp);
                double nSigPlusTemp = arrhpD[1]->GetBinContent(i+1)*binWidthTemp;
                double nSigMinusTemp = arrhmD[1]->GetBinContent(i+1)*binWidthTemp;
                double asymmTemp = ((nSigPlusTemp-nSigMinusTemp)/(nSigPlusTemp+nSigMinusTemp));
                double err = TMath::Sqrt(nSigPlusTemp+nSigMinusTemp);
                double errAsymm = TMath::Sqrt(TMath::Power(err/(nSigPlusTemp-nSigMinusTemp),2)+TMath::Power(err/(nSigPlusTemp+nSigMinusTemp),2))*asymmTemp;
                arrhAsymmD[1]->SetBinContent(i+1,asymmTemp);
                arrhAsymmD[1]->SetBinError(i+1,errAsymm);

              }
              if(j==5){
                arrhmD[2]->SetBinContent(i+1,yWMinusTotCent40_80[i]/binWidthTemp); 
                arrhmD[2]->SetBinError(i+1,TMath::Sqrt(yWMinusTotCent40_80[i])/binWidthTemp);
                double nSigPlusTemp = arrhpD[2]->GetBinContent(i+1)*binWidthTemp;
                double nSigMinusTemp = arrhmD[2]->GetBinContent(i+1)*binWidthTemp;
                double asymmTemp = ((nSigPlusTemp-nSigMinusTemp)/(nSigPlusTemp+nSigMinusTemp));
                double err = TMath::Sqrt(nSigPlusTemp+nSigMinusTemp);
                double errAsymm = TMath::Sqrt(TMath::Power(err/(nSigPlusTemp-nSigMinusTemp),2)+TMath::Power(err/(nSigPlusTemp+nSigMinusTemp),2))*asymmTemp;
                arrhAsymmD[2]->SetBinContent(i+1,asymmTemp);
                arrhAsymmD[2]->SetBinError(i+1,errAsymm);
              }


            if(i==nWEta-1)  {
                //Write(outFile, (TGraphErrors*)_grEtaDistroMinusCent.At(j),sGraphNameMinus); 
            
                TString cutsMc = "mc_mu_gen_mothertype==-24";  cutsMc+="&&mc_mu_charge==-1";
                    cutsMc+="&&mc_nu_gen_pt>25.&&mc_mu_gen_pt>25.0"; cutsMc+="&&abs(mc_mu_gen_eta)>";cutsMc+=etaBins[0];
                    cutsMc+="&&abs(mc_mu_gen_eta)<"; cutsMc+=etaBins[nEtaBins]; cutsMc+="&&centrality>";
                    cutsMc+=centralityBins[j]; cutsMc+="&&centrality<"; cutsMc+=centralityBins[j+1];
                    cutsMc+="&&TMath::Sqrt(2.0*mc_mu_gen_pt*mc_nu_gen_pt*(1.0-TMath::Cos(mc_mu_gen_phi-mc_nu_gen_phi)))>40.0";
                    std::cout << "cutsMc Minus " << cutsMc << std::endl;

                    char cMcHistoMinus[50];
                    sprintf(cMcHistoMinus,"mcTruthEtaDistroMinusCent%i",j);
                    arrhmMC[j] = new TH1D(cMcHistoMinus,cMcHistoMinus,nEtaBins,arrBinning);
                    fillMcHistos(treeMC, cMcHistoMinus,arrhmMC[j], cutsMc);
                    //double norm = arrhmD[j]->Integral()/arrhmMC[j]->Integral();
                    //arrhmMC[j]->Scale(norm); 
                    arrhmMC[j]->Scale(1.0,"width");
                    writeMcHistos(outFile, treeMC,cMcHistoMinus, arrhmMC[j]);

                    //arrhmD[j]->Scale(1.0,"width");

                    outFile->cd();
                    if(j==2){
                        arrhpD[0]->Write();
                        arrhAsymmD[0]->Write();
                    }
                    if(j==4){
                        arrhpD[1]->Write();
                        arrhAsymmD[1]->Write();
                    }
                    if(j==5){
                        arrhpD[2]->Write();
                        arrhAsymmD[2]->Write();
                    }

            }
        }
	   }///k
      }///yW0 
     } //nPts
    } //arrWSig
   } //centrality

	eyWPlusH[i] = sqrt(eyWPlusSumH[i]);
	eyWPlusL[i] = sqrt(eyWPlusSumL[i]);
    ///correlated systematic errors
	eyWPlusSyst1[i] = sqrt(eyWPlusSumSyst1[i]);
	eyWPlusSyst2[i] = sqrt(eyWPlusSumSyst2[i]);
	eyWPlusSyst3[i] = sqrt(eyWPlusSumSyst3[i]);
	eyWPlusSyst4[i] = sqrt(eyWPlusSumSyst4[i]);
	eyWPlusSyst5[i] = sqrt(eyWPlusSumSyst5[i]);
	eyWPlusSyst6[i] = sqrt(eyWPlusSumSyst6[i]);
	eyWPlusSyst7[i] = sqrt(eyWPlusSumSyst7[i]);
	eyWPlusSyst8[i] = sqrt(eyWPlusSumSyst8[i]);
	eyWPlusSyst9[i] = sqrt(eyWPlusSumSyst9[i]);
    ///Uncorrelated systematic error
	eyWPlusSyst10[i] = sqrt(eyWPlusSumSyst10[i]);

    ///Now calculate the total correlated systematic in this bin
	eyWPlusSystH[i] = sqrt(TMath::Power(eyWPlusSyst1[i],2)+TMath::Power(eyWPlusSyst2[i],2)+TMath::Power(eyWPlusSyst3[i],2)
        +TMath::Power(eyWPlusSyst4[i],2)+TMath::Power(eyWPlusSyst5[i],2)+TMath::Power(eyWPlusSyst6[i],2)+TMath::Power(eyWPlusSyst7[i],2)
        +TMath::Power(eyWPlusSyst8[i],2)+TMath::Power(eyWPlusSyst9[i],2));
	eyWPlusSystL[i] = eyWPlusSystH[i];
    
    std::cout << "total W+ in eta bin" << i << " = " << yWPlusTot[i] << " +- " << eyWPlusH[i] << "(stat.) " << eyWPlusSystH[i] << "(syst.)" << std::endl;

//	writeYieldsToSpreadsheet(spreadSheetMuPlus,i, yWPlusTot[i] , eyWPlusH[i], eyWPlusSystH[i] );

	eyWMinusH[i] = sqrt(eyWMinusSumH[i]);
	eyWMinusL[i] = sqrt(eyWMinusSumL[i]);
    ///correlated systematic error
	eyWMinusSyst1[i] = sqrt(eyWMinusSumSyst1[i]);
	eyWMinusSyst2[i] = sqrt(eyWMinusSumSyst2[i]);
	eyWMinusSyst3[i] = sqrt(eyWMinusSumSyst3[i]);
	eyWMinusSyst4[i] = sqrt(eyWMinusSumSyst4[i]);
	eyWMinusSyst5[i] = sqrt(eyWMinusSumSyst5[i]);
	eyWMinusSyst6[i] = sqrt(eyWMinusSumSyst6[i]);
	eyWMinusSyst7[i] = sqrt(eyWMinusSumSyst7[i]);
	eyWMinusSyst8[i] = sqrt(eyWMinusSumSyst8[i]);
	eyWMinusSyst9[i] = sqrt(eyWMinusSumSyst9[i]);
	eyWMinusSyst10[i] = sqrt(eyWMinusSumSyst10[i]);

    ///Now calculate the total correlated systematic in this bin
	eyWMinusSystH[i] = sqrt(TMath::Power(eyWMinusSyst1[i],2)+TMath::Power(eyWMinusSyst2[i],2)+TMath::Power(eyWMinusSyst3[i],2)
        +TMath::Power(eyWMinusSyst4[i],2)+TMath::Power(eyWMinusSyst5[i],2)+TMath::Power(eyWMinusSyst6[i],2)+TMath::Power(eyWMinusSyst7[i],2)
        +TMath::Power(eyWMinusSyst8[i],2)+TMath::Power(eyWMinusSyst9[i],2));
	eyWMinusSystL[i] = eyWMinusSystH[i];
    
    std::cout << "total W- in eta bin" << i << " = " << yWMinusTot[i] << " +- " << eyWMinusH[i] << "(stat.) " << eyWMinusSystH[i] << "(syst.)" << std::endl;
//	writeYieldsToSpreadsheet(spreadSheetMuMinus,i, yWMinusTot[i] , eyWMinusH[i], eyWMinusSystH[i] );
  } //eta
 


  TGraphAsymmErrors* grassFit = new TGraphAsymmErrors(nEtaBins);
  TGraphAsymmErrors* grassFitSyst = new TGraphAsymmErrors(nEtaBins);
  TGraphAsymmErrors* grass = new TGraphAsymmErrors(nEtaBins); 
  TGraphAsymmErrors* grassSyst[ncent]; 
  int _nPoints = nEtaBins;
  TGraphAsymmErrors* grWp = new TGraphAsymmErrors(_nPoints);
  TGraphAsymmErrors* grWm = new TGraphAsymmErrors(_nPoints);
  TGraphAsymmErrors* grWpUncorrelatedSyst = new TGraphAsymmErrors(_nPoints);
  TGraphAsymmErrors* grWmUncorrelatedSyst = new TGraphAsymmErrors(_nPoints);
  TGraphAsymmErrors* grWpCorrelatedSyst = new TGraphAsymmErrors(_nPoints);
  TGraphAsymmErrors* grWmCorrelatedSyst = new TGraphAsymmErrors(_nPoints);
  TH1F* hAssSyst = new TH1F("hAssSyst","hAssSyst",nEtaBins,arrBinning);

//  for (int i=0; i<1;i++){
     grass = new TGraphAsymmErrors(nEtaBins);
     grassSyst[0] = new TGraphAsymmErrors(nEtaBins);
//  }
  
  const int neta = nWEta;

  double yWPlusTotal=0; 
  double yWMinusTotal=0; 
  ///No. of absolute eta bins
  //const int nAbsEtaBins = (nWEta-1)/2;
  const int nAbsEtaBins = nEtaBins;
  TString sEtaRange[] = {"0.1-0.35","0.35-0.6","0.6-0.8","0.8-1.05","1.05-1.3","1.3-1.55","1.55-1.85","1.85-2.1","2.1-2.4"};
//    for (int j=0; j<nWEta;j++)
    for (int j=0; j<nAbsEtaBins;j++){

      std::cout << "calculating asymmetry for eta bin: " << j << std::endl;

     ///perform the mapping of +-eta bins
     ///to figure out which eta bins share the
     ///same bin number
     if(doAtlasEta){
        int binNeg = indexNegEta(j,nWEta);
        int binPos = indexPosEta(j,nWEta);
        ///get the number of signal candidates in 
        ///each +-eta slice
        double yW0PlusNegEta = yWPlusTot[binNeg];
        double yW0PlusPosEta = yWPlusTot[binPos];
        ///add them together to get the absolute eta yield
        yWPlusAbsEtaSum[j] = yW0PlusNegEta+yW0PlusPosEta;


        ///ditto for W-
        double yW0MinusNegEta = yWMinusTot[binNeg];
        double yW0MinusPosEta = yWMinusTot[binPos];

        yWMinusAbsEtaSum[j] = yW0MinusNegEta+yW0MinusPosEta;

        ///errors
        ///absolute stat error of W+ and W-
        eyWPlusAbsEtaSum[j] = TMath::Sqrt(TMath::Power(TMath::Sqrt(yW0PlusNegEta),2)+TMath::Power(TMath::Sqrt(yW0PlusPosEta),2)); 
        eyWMinusAbsEtaSum[j] = TMath::Sqrt(TMath::Power(TMath::Sqrt(yW0MinusNegEta),2)+TMath::Power(TMath::Sqrt(yW0MinusPosEta),2)); 
        ///systematic error of W+ and W-
        eyUncorrelatedSystPlusAbsEtaSum[j] = TMath::Sqrt(TMath::Power(eyWPlusSystH[binNeg],2)+TMath::Power(eyWPlusSystH[binPos],2)); 
        eyUncorrelatedSystMinusAbsEtaSum[j] = TMath::Sqrt(TMath::Power(eyWMinusSystH[binNeg],2)+TMath::Power(eyWMinusSystH[binPos],2)); 

        std::cout << "Eta bin " << j << " gives W+ yield " << yW0PlusNegEta << "(" << TMath::Sqrt(yW0PlusNegEta)<<")" <<"stat. +" 
        << yW0PlusPosEta << "(" << TMath::Sqrt(yW0PlusPosEta) <<")" <<"stat. = " << yWPlusAbsEtaSum[j] << "(" <<
        eyWPlusAbsEtaSum[j] <<
            ")" <<"stat." << std::endl;


        std::cout << "Eta bin " << j << " gives W- yield " << yW0MinusNegEta << "(" << TMath::Sqrt(yW0MinusNegEta)<<")" <<"stat. +" 
        << yW0MinusPosEta << "(" << TMath::Sqrt(yW0MinusPosEta) <<")" <<"stat. = " << yWMinusAbsEtaSum[j] << "(" <<
        eyWMinusAbsEtaSum[j] <<
            ")" <<"stat." << std::endl;
     }
     else{
         yWPlusAbsEtaSum[j] =yWPlusTot[j];
         yNObsPlusAbsEtaSum[j] =yNObsPlusTot[j];
         yNBkgPlusAbsEtaSum[j] =yNBkgPlusTot[j];
         yWMinusAbsEtaSum[j] = yWMinusTot[j];
         yNObsMinusAbsEtaSum[j] =yNObsMinusTot[j];
         yNBkgMinusAbsEtaSum[j] =yNBkgMinusTot[j];
         eyWPlusAbsEtaSum[j] = eyWPlusH[j];
         eyWMinusAbsEtaSum[j] = eyWMinusH[j];
         eyUncorrelatedSystPlusAbsEtaSum[j] = eyWPlusSyst10[j];
         eyUncorrelatedSystMinusAbsEtaSum[j] = eyWMinusSyst10[j];

         eyCorrelatedSystPlusAbsEtaSum[j] = eyWPlusSystH[j];
         eyCorrelatedSystMinusAbsEtaSum[j] = eyWMinusSystH[j];

        std::cout << "Eta bin " << j << " gives W+ yield " <<  
            yWPlusAbsEtaSum[j] << "(" << eyWPlusAbsEtaSum[j] <<
            ")" <<"stat." << std::endl;
        std::cout << "Eta bin " << j << " gives W+ yield " <<  
            yWMinusAbsEtaSum[j] << "(" << eyWMinusAbsEtaSum[j] <<
            ")" <<"stat." << std::endl;

     }
     double yWPlusTemp = yWPlusAbsEtaSum[j] ;
     double yNObsPlusTemp = yNObsPlusAbsEtaSum[j] ;
     double yNBkgPlusTemp = yNBkgPlusAbsEtaSum[j] ;
     ///add them together to get the absolute eta yield
     yWPlusTotal += yWPlusTemp;

     double yWMinusTemp = yWMinusAbsEtaSum[j] ;
     double yNObsMinusTemp = yNObsMinusAbsEtaSum[j] ;
     double yNBkgMinusTemp = yNBkgMinusAbsEtaSum[j] ;
     yWMinusTotal += yWMinusTemp;


	 ///asymmetry in this eta bin
     double assym = (yWPlusAbsEtaSum[j]-yWMinusAbsEtaSum[j])/(yWPlusAbsEtaSum[j]+yWMinusAbsEtaSum[j]) ;

/*         double enpH = eyWPlusH[j];
         double enpL = eyWPlusL[j] ;
         double enmH = eyWMinusH[j] ;
         double enmL = eyWMinusL[j] ;
*/
         double enpH = eyWPlusAbsEtaSum[j];
         double enpL = eyWPlusAbsEtaSum[j] ;
         double enmH = eyWMinusAbsEtaSum[j] ;
         double enmL = eyWMinusAbsEtaSum[j] ;

	 //relative stat error
         double eStatH0 = sqrt(enpH*enpH  + enmH*enmH ) ;
         double eStatH1 = sqrt( pow(100*eStatH0/(yWPlusAbsEtaSum[j]-yWMinusAbsEtaSum[j]),2) + pow(100*eStatH0/(yWPlusAbsEtaSum[j]+yWMinusAbsEtaSum[j]),2));
         double eStatL0 = sqrt(enpL*enpL  + enmL*enmL) ;
         double eStatL1 = sqrt( pow(100*eStatL0/(yWPlusAbsEtaSum[j]-yWMinusAbsEtaSum[j]),2) + pow(100*eStatL0/(yWPlusAbsEtaSum[j]+yWMinusAbsEtaSum[j]),2));

     ///Add in uncorrelated errors
      double sum = yWPlusAbsEtaSum[j]+yWMinusAbsEtaSum[j];
      double diff = yWPlusAbsEtaSum[j]-yWMinusAbsEtaSum[j];
      double eAsymmUncorrelatedSum = TMath::Sqrt(TMath::Power(enpH,2)+TMath::Power(enmH,2)
                                    +TMath::Power(eyUncorrelatedSystPlusAbsEtaSum[j],2)+TMath::Power(eyUncorrelatedSystMinusAbsEtaSum[j],2))/sum;
      double eAsymmUncorrelatedDiff = eAsymmUncorrelatedSum*sum/diff;
      //double  eAsymmUncorrelated = TMath::Sqrt(pow(eAsymmUncorrelatedDiff,2)+pow(eAsymmUncorrelatedSum,2));
      double eAsymmMCStatErrSum = TMath::Sqrt(TMath::Power(eyUncorrelatedSystPlusAbsEtaSum[j],2)+TMath::Power(eyUncorrelatedSystMinusAbsEtaSum[j],2))/sum;
      double  eAsymmMCStatErrDiff =eAsymmMCStatErrSum*sum/diff;
      double eAsymmMCStatErr = TMath::Sqrt(TMath::Power(eAsymmMCStatErrSum,2)+TMath::Power(eAsymmMCStatErrDiff,2))*100.0;
      double A = yWPlusAbsEtaSum[j];
      double B = yWMinusAbsEtaSum[j];
      double sigmaA = TMath::Sqrt(TMath::Power(eyUncorrelatedSystPlusAbsEtaSum[j],2));
      double sigmaB = TMath::Sqrt(TMath::Power(eyUncorrelatedSystMinusAbsEtaSum[j],2));
      //double eAsymmUncorrelatedDiff = TMath::Sqrt(TMath::Power(sigmaA,2)+TMath::Power(sigmaB,2))/diff;
      std::cout << "sum: " << sum << std::endl;
      std::cout << "diff: " << diff << std::endl;
      std::cout << "Sqrt(sigmaA^2+sigmaB^2) " << TMath::Sqrt(TMath::Power(sigmaA,2)+TMath::Power(sigmaB,2)) << std::endl; 
      ///absolute uncorrelated error on asymmetry
      double eAsymmUncorrelated = 2.0*TMath::Sqrt((B*B*sigmaA*sigmaA+A*A*sigmaB*sigmaB))/TMath::Power((A+B),2);


	 ///additional systematic errors from 
  	 ///inc(dec) Cw fit coefficients by 1 sigma
     /*
	 double systErr1SigmaMuPlus = 0.0023*yWPlusTot[j];
	 double systErr1SigmaMuMinus = 0.0031*yWMinusAbsEtaSum[j];
*/

	 //absolute syst error
     /*    double eSystH0 = sqrt( pow(eyWPlusSystH[j],2)  + pow(eyWMinusSystH[j],2) ) ;
	 //relative 
         double eSystH1 = sqrt( pow(100*eSystH0/(yWPlusAbsEtaSum[j]-yWMinusAbsEtaSum[j]),2) + pow(100*eSystH0/(yWPlusAbsEtaSum[j]+yWMinusAbsEtaSum[j]),2));
         double eSystL0 = sqrt( pow(eyWPlusSystL[j],2) + pow(eyWMinusSystL[j],2)) ;
         double eSystL1 = sqrt( pow(100*eSystL0/(yWPlusAbsEtaSum[j]-yWMinusAbsEtaSum[j]),2) + pow(100*eSystL0/(yWPlusAbsEtaSum[j]+yWMinusAbsEtaSum[j]),2));
	   
	 double eStatPlusHTemp =  enpH;
         double eStatPlusLTemp = enpL;
	 double eTotPlusTemp = sqrt( pow(eStatPlusHTemp,2) + pow(eyWPlusSystH[j],2) );

     double eStatMinusHTemp =  enmH;
     double eStatMinusLTemp = enmL;
	 double eTotMinusTemp = sqrt( pow(eStatMinusHTemp,2) + pow(eyWMinusSystH[j],2) );
     */

     //absolute syst error
     /*double eSystH0 = sqrt( pow(eyUncorrelatedSystPlusAbsEtaSum[j],2)  + pow(eyUncorrelatedSystMinusAbsEtaSum[j],2) ) ;
	 //relative 
     double eSystH1 = sqrt( pow(100*eSystH0/(yWPlusAbsEtaSum[j]-yWMinusAbsEtaSum[j]),2) + pow(100*eSystH0/(yWPlusAbsEtaSum[j]+yWMinusAbsEtaSum[j]),2));
     double eSystL0 = sqrt( pow(eyUncorrelatedSystPlusAbsEtaSum[j],2)  + pow(eyUncorrelatedSystMinusAbsEtaSum[j],2) ) ;
     double eSystL1 = sqrt( pow(100*eSystL0/(yWPlusAbsEtaSum[j]-yWMinusAbsEtaSum[j]),2) + pow(100*eSystL0/(yWPlusAbsEtaSum[j]+yWMinusAbsEtaSum[j]),2));
     */
     ///add in only systematics that don't completely cancel in numerator and denominator
     // Expressed as abs error on asymm
     double eSystL0 = 0.,eSystL1=0., eSystH0 =0., eSystH1=0.;
     eSystH1 = TMath::Sqrt(TMath::Power(grIsoCancellation->GetY()[j],2)+TMath::Power(grMptCancellation->GetY()[j],2));  
     eSystL1 = eSystH1;
	   
	 double eStatPlusHTemp =  enpH;
     double eStatPlusLTemp = enpL;
	 double eTotPlusTemp = sqrt( pow(eStatPlusHTemp,2) + pow(eyUncorrelatedSystPlusAbsEtaSum[j],2) );

     double eStatMinusHTemp =  enmH;
     double eStatMinusLTemp = enmL;
	 double eTotMinusTemp = sqrt( pow(eStatMinusHTemp,2) + pow(eyUncorrelatedSystMinusAbsEtaSum[j],2) );
           
	 //absoluate stat error
     double eStatH  = 0.01*eStatH1*abs(assym);
     double eStatL  = 0.01*eStatL1*abs(assym);

	 ///absoluate stat+uncorrelated+(non-cancelling)correlated syst error
     double eSystH  = TMath::Sqrt(TMath::Power(eSystH1,2)+TMath::Power(eAsymmUncorrelated,2));
     std::cout << "eSystH " << eSystH << std::endl; 
     std::cout << "asymm " << assym << std::endl;
     //eSystH*=abs(assym);
     double eSystL  = TMath::Sqrt(TMath::Power(eSystL1,2)+TMath::Power(eAsymmUncorrelated,2));
     //eSystL*=abs(assym);

         /*double systRelErrFitPlus = nWFitSystErrPlus[i]/(yWPlusTot[j]+yWMinusAbsEtaSum[j])*100; //systematic error  
         double systRelErrFitMinus = nWFitSystErrMinus[i]/(yWPlusTot[j]+yWMinusAbsEtaSum[j])*100; //systematic error 
	 double systErrSum = sqrt(nWFitSystErrPlus[i]*nWFitSystErrPlus[i]+nWFitSystErrMinus[i]*nWFitSystErrMinus[i]) ; 
	 double systErr = 0.01*sqrt( pow(100*systErrSum / (yWPlusTot[j]-yWMinusAbsEtaSum[j]),2) + pow(100*systErrSum / (yWPlusTot[j]+yWMinusAbsEtaSum[j]),2));
	 */
	 //total uncertainty(i.e. stat and syst)
	 double errTotalH = sqrt( pow(eSystH,2) + pow(eStatH,2) );
	 double errTotalL = sqrt(pow(eSystL,2) + pow(eStatL,2) );

     std::cout << "Statistical uncertainty: " << eStatH1/100.0*fabs(assym) << std::endl;  
     std::cout << "MC Statistical uncertainty: " << eAsymmMCStatErr/100.0*fabs(assym) << std::endl; 
     std::cout << "Error from non-cancellation: " << eSystH1*fabs(assym) << std::endl;
     std::cout << "Added in quadrature gives a total uncertainty of : " << errTotalL << std::endl; 

     float xEtaTemp = (etaBins[j+1]-etaBins[j])/2.0 + etaBins[j]; 
     float binWidthTemp = etaBins[j+1]-etaBins[j];
     std::cout << "x point: " << xEtaTemp << std::endl;
     std::cout << "bin width: " << binWidthTemp << std::endl;

     grass->SetPoint(j,xEtaTemp,assym);
     grassSyst[0]->SetPoint(j,xEtaTemp,assym);
     hAssSyst->SetBinContent(j+1,assym);

     grWp->SetPoint(j,xEtaTemp,yWPlusTemp/binWidthTemp);
     grWpUncorrelatedSyst->SetPoint(j,xEtaTemp,yWPlusTemp/binWidthTemp);
     grWpCorrelatedSyst->SetPoint(j,xEtaTemp,yWPlusTemp/binWidthTemp);

     grWm->SetPoint(j,xEtaTemp,yWMinusTemp/binWidthTemp);
     grWmUncorrelatedSyst->SetPoint(j,xEtaTemp,yWMinusTemp/binWidthTemp);
     grWmCorrelatedSyst->SetPoint(j,xEtaTemp,yWMinusTemp/binWidthTemp);

	 grass->SetPointError(j,(etaBins[j+1]-etaBins[j])/2.0, (etaBins[j+1]-etaBins[j])/2.0, eStatL, eStatH) ;
	 grassSyst[0]->SetPointError(j,(etaBins[j+1]-etaBins[j])/2.0, (etaBins[j+1]-etaBins[j])/2.0, errTotalL, errTotalH) ;
     hAssSyst->SetBinError(j+1,errTotalL);
     std::cout << "Stat error on asymmetry: " << eStatL << std::endl;

	 grWp->SetPointError(j,(etaBins[j+1]-etaBins[j])/2.0,
         (etaBins[j+1]-etaBins[j])/2.0,eStatPlusLTemp/binWidthTemp ,eStatPlusHTemp/binWidthTemp ) ;
     std::cout << "W+: " << eStatPlusLTemp << "/" << binWidthTemp << " = " << eStatPlusHTemp << "/" << binWidthTemp 
        << " = " << eStatPlusLTemp/binWidthTemp << std::endl;
	 grWpUncorrelatedSyst->SetPointError(j,(etaBins[j+1]-etaBins[j])/2.0,
         (etaBins[j+1]-etaBins[j])/2.0,eTotPlusTemp/binWidthTemp,eTotPlusTemp/binWidthTemp) ;
     std::cout << "Total error for bin " << j << " = " << eTotPlusTemp/binWidthTemp << std::endl;
	 grWpCorrelatedSyst->SetPointError(j,(etaBins[j+1]-etaBins[j])/2.0,
         (etaBins[j+1]-etaBins[j])/2.0,eyCorrelatedSystPlusAbsEtaSum[j]/binWidthTemp,eyCorrelatedSystPlusAbsEtaSum[j]/binWidthTemp) ;

	 grWm->SetPointError(j,(etaBins[j+1]-etaBins[j])/2.0,
         (etaBins[j+1]-etaBins[j])/2.0,eStatMinusLTemp/binWidthTemp ,eStatMinusHTemp/binWidthTemp ) ;
	 grWmUncorrelatedSyst->SetPointError(j,(etaBins[j+1]-etaBins[j])/2.0,
         (etaBins[j+1]-etaBins[j])/2.0,eTotMinusTemp/binWidthTemp,eTotMinusTemp/binWidthTemp) ;
	 grWmCorrelatedSyst->SetPointError(j,(etaBins[j+1]-etaBins[j])/2.0,
         (etaBins[j+1]-etaBins[j])/2.0,eyCorrelatedSystMinusAbsEtaSum[j]/binWidthTemp,eyCorrelatedSystMinusAbsEtaSum[j]/binWidthTemp) ;

      std::cout << "Asymmetry = " << assym << "+-" << eStatL << "(stat.)" << eSystL << "(syst.)" << std::endl;

	  writeYieldsToSpreadsheet(spreadSheetMuPlus,sEtaRange[j], format(yNObsPlusTemp),format(yNBkgPlusTemp,2),format(yWPlusTemp) , 
					format(eyWPlusH[j],2), format(eyWPlusSyst10[j],2),format(eyWPlusSystL[j],2) );
	  writeYieldsToSpreadsheet(spreadSheetMuMinus,sEtaRange[j],format(yNObsMinusTemp),format(yNBkgMinusTemp,2), format(yWMinusTemp) , 
					format(eyWMinusH[j],2), format(eyWMinusSyst10[j],2),format(eyWMinusSystL[j],2) );

      writeAsymmToSpreadsheet(spreadSheetAsymm,sEtaRange[j],format(assym,2),format(eStatL,2),format(eAsymmUncorrelated,2),format(eSystL1,2));
    } //j 

   spreadSheetMuPlus.close();
   spreadSheetMuMinus.close();

   std::cout << "Total W+ = " << yWPlusTotal << std::endl;
   std::cout << "Total W- = " << yWMinusTotal << std::endl;
 
  //float xBins[] = {0.1,0.35,0.6,0.8,1.05,1.3,1.55,1.85,2.1,2.4};
  float xBins[] = {0.1,0.35,0.6,0.8,1.05,1.37,1.52,1.74,2.1,2.40};
  int bins = sizeof(xBins)/sizeof(float)-1;
  TH1F* hWm, *hWp;
  hWp = convertToHisto(grWp,"hWp",bins,xBins);
  hWm = convertToHisto(grWm,"hWm",bins,xBins);

  ///fill asymmetry for different nucleon components
  //float xBins[] = {0.0,0.1,0.35,0.6,0.8,1.05,1.3,1.55,1.85,2.1,2.4,2.7,3.2,4.0,5.0};
  TH1F* hEtaPlusMcAll = new TH1F("hEtaPlusMcAll","hEtaPlusMcAll",bins,xBins);
  hEtaPlusMcAll->Sumw2();
  TH1F* hEtaMinusMcAll = new TH1F("hEtaMinusMcAll","hEtaMinusMcAll",bins,xBins);
  hEtaMinusMcAll->Sumw2();
  TH1F* hEtaPlus_pp = new TH1F("hEtaPlus_pp","hEtaPlus_pp",bins,xBins);
  hEtaPlus_pp->Sumw2();
  TH1F* hEtaMinus_pp = new TH1F("hEtaMinus_pp","hEtaMinus_pp",bins,xBins);
  hEtaMinus_pp->Sumw2();
  TH1F* hEtaPlus_pn = new TH1F("hEtaPlus_pn","hEtaPlus_pn",bins,xBins);
  hEtaPlus_pn->Sumw2();
  TH1F* hEtaMinus_pn = new TH1F("hEtaMinus_pn","hEtaMinus_pn",bins,xBins);
  hEtaMinus_pn->Sumw2();
  TH1F* hEtaPlus_np = new TH1F("hEtaPlus_np","hEtaPlus_np",bins,xBins);
  hEtaPlus_np->Sumw2();
  TH1F* hEtaMinus_np = new TH1F("hEtaMinus_np","hEtaMinus_np",bins,xBins);
  hEtaMinus_np->Sumw2();
  TH1F* hEtaPlus_nn = new TH1F("hEtaPlus_nn","hEtaPlus_nn",bins,xBins);
  hEtaPlus_nn->Sumw2();
  TH1F* hEtaMinus_nn = new TH1F("hEtaMinus_nn","hEtaMinus_nn",bins,xBins);
  hEtaMinus_nn->Sumw2();
  TH1F* hEtaPlus_crossNucleon = new TH1F("hEtaPlus_crossNucleon","hEtaPlus_crossNucleon",bins,xBins);
  hEtaPlus_crossNucleon->Sumw2();
  TH1F* hEtaMinus_crossNucleon = new TH1F("hEtaMinus_crossNucleon","hEtaMinus_crossNucleon",bins,xBins);
  hEtaMinus_crossNucleon->Sumw2();

  TH1F* hEtaPlusMSTWpp = new TH1F("hEtaPlusMSTWpp","hEtaPlusMSTWpp",bins,xBins);
  hEtaPlusMSTWpp->Sumw2();
  TH1F* hEtaPlusMSTWnn = new TH1F("hEtaPlusMSTWnn","hEtaPlusMSTWnn",bins,xBins);
  hEtaPlusMSTWnn->Sumw2();
  TH1F* hEtaPlusMSTWnppn = new TH1F("hEtaPlusMSTWnppn","hEtaPlusMSTWnppn",bins,xBins);
  hEtaPlusMSTWnppn->Sumw2();
  TH1F* hEtaMinusMSTWpp = new TH1F("hEtaMinusMSTWpp","hEtaMinusMSTWpp",bins,xBins);
  hEtaMinusMSTWpp->Sumw2();
  TH1F* hEtaMinusMSTWnn = new TH1F("hEtaMinusMSTWnn","hEtaMinusMSTWnn",bins,xBins);
  hEtaMinusMSTWnn->Sumw2();
  TH1F* hEtaMinusMSTWnppn = new TH1F("hEtaMinusMSTWnppn","hEtaMinusMSTWnppn",bins,xBins);
  hEtaMinusMSTWnppn->Sumw2();
  TH1F* hEtaPlusMSTW = new TH1F("hEtaPlusMSTW","hEtaPlusMSTW",bins,xBins);
  hEtaPlusMSTW->Sumw2();
  TH1F* hEtaMinusMSTW = new TH1F("hEtaMinusMSTW","hEtaMinusMSTW",bins,xBins);
  hEtaMinusMSTW->Sumw2();

  TH1F* hEtaPlusCT10pp = new TH1F("hEtaPlusCT10pp","hEtaPlusCT10pp",bins,xBins);
  hEtaPlusCT10pp->Sumw2();
  TH1F* hEtaPlusCT10nn = new TH1F("hEtaPlusCT10nn","hEtaPlusCT10nn",bins,xBins);
  hEtaPlusCT10nn->Sumw2();
  TH1F* hEtaPlusCT10np = new TH1F("hEtaPlusCT10np","hEtaPlusCT10np",bins,xBins);
  hEtaPlusCT10np->Sumw2();
  TH1F* hEtaPlusCT10pn = new TH1F("hEtaPlusCT10pn","hEtaPlusCT10pn",bins,xBins);
  hEtaPlusCT10pn->Sumw2();
  TH1F* hEtaMinusCT10pp = new TH1F("hEtaMinusCT10pp","hEtaMinusCT10pp",bins,xBins);
  hEtaMinusCT10pp->Sumw2();
  TH1F* hEtaMinusCT10nn = new TH1F("hEtaMinusCT10nn","hEtaMinusCT10nn",bins,xBins);
  hEtaMinusCT10nn->Sumw2();
  TH1F* hEtaMinusCT10np = new TH1F("hEtaMinusCT10np","hEtaMinusCT10np",bins,xBins);
  hEtaMinusCT10np->Sumw2();
  TH1F* hEtaMinusCT10pn = new TH1F("hEtaMinusCT10pn","hEtaMinusCT10pn",bins,xBins);
  hEtaMinusCT10pn->Sumw2();
  TH1F* hEtaPlusCT10 = new TH1F("hEtaPlusCT10","hEtaPlusCT10",bins,xBins);
  hEtaPlusCT10->Sumw2();
  TH1F* hEtaMinusCT10 = new TH1F("hEtaMinusCT10","hEtaMinusCT10",bins,xBins);
  hEtaMinusCT10->Sumw2();

  TH1F* hAsymmMcAll = new TH1F("hAsymmMcAll","hAsymmMcAll",bins,xBins);
  TH1F* hAsymm_pp = new TH1F("hAsymm_pp","hAsymm_pp",bins,xBins);
  TH1F* hAsymm_pn = new TH1F("hAsymm_pn","hAsymm_pn",bins,xBins);
  TH1F* hAsymm_np = new TH1F("hAsymm_np","hAsymm_np",bins,xBins);
  TH1F* hAsymm_nn = new TH1F("hAsymm_nn","hAsymm_nn",bins,xBins);
//  TH1F* hAsymm_crossNucleon = new TH1F("hAsymm_crossNucleon","hAsymm_crossNucleon",bins,xBins);
  TH1F* hAsymmMSTWnn = new TH1F("hAsymmMSTWnn","hAsymmMSTWnn",bins,xBins);
  TH1F* hAsymmMSTWpp = new TH1F("hAsymmMSTWpp","hAsymmMSTWpp",bins,xBins);
  TH1F* hAsymmMSTWnppn = new TH1F("hAsymmMSTWnppn","hAsymmMSTWnppn",bins,xBins);
//  TH1F* hAsymmMSTW = new TH1F("hAsymmMSTW","hAsymmMSTW",bins,xBins);

  double WCrossSectionPythiaLO = 2.8205e-9/64.0e-3*1.0e9;
  double WCrossSectionMinusPythiaLO = (1.0/2.01)*WCrossSectionPythiaLO;
  double WCrossSectionPlusPythiaLO = (1.0-1.0/2.01)*WCrossSectionPythiaLO;

  ///fill asymmetry for Pythia MRST LO* sample
  fillMcHistos(fileNameMCWIn,bins,xBins,hEtaPlusMcAll,hEtaMinusMcAll,WCrossSectionPlusPythiaLO,WCrossSectionMinusPythiaLO,"AllNucleons"); 
  hAsymmMcAll = getAsymmetry(hEtaPlusMcAll,hEtaMinusMcAll);
  for(int i =1; i<hAsymmMcAll->GetNbinsX(); ++i){
     std::cout << "asymm MRST LO: " << hAsymmMcAll->GetBinContent(i) << std::endl; 
  }

  fillMcHistos(fileNamepp,bins,xBins,hEtaPlus_pp,hEtaMinus_pp,WCrossSectionPlusPythiaLO,WCrossSectionMinusPythiaLO,"ppTemp"); 
  fillMcHistos(fileNamenp,bins,xBins,hEtaPlus_np,hEtaMinus_np,WCrossSectionPlusPythiaLO,WCrossSectionMinusPythiaLO,"npTemp"); 
  fillMcHistos(fileNamepn,bins,xBins,hEtaPlus_pn,hEtaMinus_pn,WCrossSectionPlusPythiaLO,WCrossSectionMinusPythiaLO,"pnTemp"); 
  fillMcHistos(fileNamenn,bins,xBins,hEtaPlus_nn,hEtaMinus_nn,WCrossSectionPlusPythiaLO,WCrossSectionMinusPythiaLO,"nnTemp"); 
  hEtaPlus_crossNucleon->Add(hEtaPlus_np,hEtaPlus_pn);
  hEtaMinus_crossNucleon->Add(hEtaMinus_np,hEtaMinus_pn);
  TH1F* hAsymm_crossNucleon = getAsymmetry(hEtaPlus_crossNucleon,hEtaMinus_crossNucleon);

  // CT10
  double WCrossSection_ppPlusCT10nlo = 2.1120e-9/64.0e-3*1.0e9;
  double WCrossSection_ppMinusCT10nlo = 1.2420e-9/64.0e-3*1.0e9;
  double wt_pp = 0.155;
  double WCrossSection_nppnPlusCT10nlo = (1.6690e-9+1.6700e-9)/2/64.0e-3*1.0e9;
  double WCrossSection_nppnMinusCT10nlo = (1.6540e-9+1.6540e-9)/2/64.0e-3*1.0e9;
  double wt_nppn = 0.478;
  double WCrossSection_npPlusCT10nlo = 1.6690e-9/64.0e-3*1.0e9;
  double WCrossSection_npMinusCT10nlo = 1.6540e-9/64.0e-3*1.0e9;
  double wt_np = 0.478/2.;
  double WCrossSection_pnPlusCT10nlo = 1.6700e-9/64.0e-3*1.0e9;
  double WCrossSection_pnMinusCT10nlo = 1.6540e-9/64.0e-3*1.0e9;
  double wt_pn = 0.478/2.;
  double WCrossSection_nnPlusCT10nlo = 1.2530e-9/64.0e-3*1.0e9;
  double WCrossSection_nnMinusCT10nlo = 2.0920e-9/64.0e-3*1.0e9;
  double wt_nn = 0.367;

  double totalWCrossSection_NNPlusCT10nlo = WCrossSection_ppPlusCT10nlo+WCrossSection_nppnPlusCT10nlo+WCrossSection_nnPlusCT10nlo;
  double totalWCrossSection_NNMinusCT10nlo = WCrossSection_ppMinusCT10nlo+WCrossSection_nppnMinusCT10nlo+WCrossSection_nnMinusCT10nlo;

  double wtdWCrossSectionCT10nloPlus  = (wt_pp*totalWCrossSection_NNPlusCT10nlo
                                            +wt_nppn*totalWCrossSection_NNPlusCT10nlo
                                            +wt_nn*totalWCrossSection_NNPlusCT10nlo)/1.0;
  double wtdWCrossSectionCT10nloMinus  = (wt_pp*totalWCrossSection_NNMinusCT10nlo
                                            +wt_nppn*totalWCrossSection_NNMinusCT10nlo
                                            +wt_nn*totalWCrossSection_NNMinusCT10nlo)/1.0;

  // Fill for all nucleon combos for CT10 NLO
  fillCT10Histos(filenameCT10Plusnn,bins,xBins,hEtaPlusCT10nn,
                    wt_nn*WCrossSection_nnPlusCT10nlo,"nnTemp");
  fillCT10Histos(filenameCT10Pluspp,bins,xBins,hEtaPlusCT10pp,
                    wt_pp*WCrossSection_ppPlusCT10nlo,"ppTemp");
  fillCT10Histos(filenameCT10Pluspn,bins,xBins,hEtaPlusCT10pn,
                    wt_pn*WCrossSection_pnPlusCT10nlo,"pnTemp");
  fillCT10Histos(filenameCT10Plusnp,bins,xBins,hEtaPlusCT10np,
                    wt_np*WCrossSection_npPlusCT10nlo,"npTemp");

  fillCT10Histos(filenameCT10Minusnn,bins,xBins,hEtaMinusCT10nn,
                    wt_nn*WCrossSection_nnMinusCT10nlo,"nnTemp");
  fillCT10Histos(filenameCT10Minuspp,bins,xBins,hEtaMinusCT10pp,
                    wt_pp*WCrossSection_ppMinusCT10nlo,"ppTemp");
  fillCT10Histos(filenameCT10Minuspn,bins,xBins,hEtaMinusCT10pn,
                    wt_pn*WCrossSection_pnMinusCT10nlo,"pnTemp");
  fillCT10Histos(filenameCT10Minusnp,bins,xBins,hEtaMinusCT10np,
                    wt_np*WCrossSection_npMinusCT10nlo,"npTemp");


  // Combine nucleon combinations for CT10
  hEtaPlusCT10->Add(hEtaPlusCT10pp,hEtaPlusCT10nn);
  hEtaPlusCT10->Add(hEtaPlusCT10np);
  hEtaPlusCT10->Add(hEtaPlusCT10pn);
  hEtaMinusCT10->Add(hEtaMinusCT10pp,hEtaMinusCT10nn);
  hEtaMinusCT10->Add(hEtaMinusCT10np);
  hEtaMinusCT10->Add(hEtaMinusCT10pn);

  // normalize to data and take ratio
  hEtaPlusCT10->Scale(hWp->Integral()/hEtaPlusCT10->Integral());
  hEtaMinusCT10->Scale(hWm->Integral()/hEtaMinusCT10->Integral());
  TH1F* hRatioCT10Plus = new TH1F("hRatioCT10Plus","hRatioCT10Plus",bins,xBins);
  hRatioCT10Plus = getMCDataRatio(hRatioCT10Plus,hEtaPlusCT10,hWp);
  TF1* func = new TF1("pol3","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0.0,2.5);
  func->SetParameter(0,0.9769015);
  func->SetParameter(1,-0.413709);
  func->SetParameter(2,0.6514979);
  func->SetParameter(3,-0.2109735);
  hRatioCT10Plus->Fit("pol3");
  TH1F* hRatioCT10Minus = new TH1F("hRatioCT10Minus","hRatioCT10Minus",bins,xBins);
  hRatioCT10Minus = getMCDataRatio(hRatioCT10Minus,hEtaMinusCT10,hWm);
  hRatioCT10Minus->Fit("pol3");

  // Fill asymmetry for POW+PY8 MSTW
  double WCrossSection_ppPlusMSTWnlo = 2.175e-9/64.0e-3*1.0e9;
  double WCrossSection_ppMinusMSTWnlo = 1.330e-9/64.0e-3*1.0e9;
  double WCrossSection_nppnPlusMSTWnlo = 1.743e-9/64.0e-3*1.0e9;
  double WCrossSection_nppnMinusMSTWnlo = 1.744e-9/64.0e-3*1.0e9;
  double WCrossSection_nnPlusMSTWnlo = 1.327e-9/64.0e-3*1.0e9;
  double WCrossSection_nnMinusMSTWnlo = 2.175e-9/64.0e-3*1.0e9;

  double totalWCrossSection_NNPlusMSTWnlo = WCrossSection_ppPlusMSTWnlo+WCrossSection_nppnPlusMSTWnlo+WCrossSection_nnPlusMSTWnlo;
  double totalWCrossSection_NNMinusMSTWnlo = WCrossSection_ppMinusMSTWnlo+WCrossSection_nppnMinusMSTWnlo+WCrossSection_nnMinusMSTWnlo;

  double wtdWCrossSectionMSTWnloPlus  = (wt_pp*totalWCrossSection_NNPlusMSTWnlo
                                            +wt_nppn*totalWCrossSection_NNPlusMSTWnlo
                                            +wt_nn*totalWCrossSection_NNPlusMSTWnlo)/1.0;
  /*double wtdWCrossSectionMSTWnloPlus  = (wt_pp*WCrossSection_ppPlusMSTWnlo
                                            +wt_nppn*WCrossSection_nppnPlusMSTWnlo
                                            +wt_nn*WCrossSection_nnPlusMSTWnlo)/1.0;
  */

  double wtdWCrossSectionMSTWnloMinus  = (wt_pp*totalWCrossSection_NNMinusMSTWnlo
                                            +wt_nppn*totalWCrossSection_NNMinusMSTWnlo
                                            +wt_nn*totalWCrossSection_NNMinusMSTWnlo)/1.0;

  /*double wtdWCrossSectionMSTWnloMinus = (wt_pp*WCrossSection_ppMinusMSTWnlo
                                            +wt_nppn*WCrossSection_nppnMinusMSTWnlo
                                            +wt_nn*WCrossSection_nnMinusMSTWnlo)/1.0;
  */
  // Fill asymmetry for all nucleon combos for MSTW NLO
  fillMSTWHistos(filenameMSTWnn,bins,xBins,hEtaPlusMSTWnn,hEtaMinusMSTWnn,
                    wt_nn*WCrossSection_nnPlusMSTWnlo,wt_nn*WCrossSection_nnMinusMSTWnlo,"nnTemp");
  fillMSTWHistos(filenameMSTWpp,bins,xBins,hEtaPlusMSTWpp,hEtaMinusMSTWpp,
                    wt_pp*WCrossSection_ppPlusMSTWnlo,wt_pp*WCrossSection_ppMinusMSTWnlo,"ppTemp");
  fillMSTWHistos(filenameMSTWnppn,bins,xBins,hEtaPlusMSTWnppn,hEtaMinusMSTWnppn,
                    wt_nppn*WCrossSection_nppnPlusMSTWnlo,wt_nppn*WCrossSection_nppnMinusMSTWnlo,"nppnTemp");


/*  hEtaPlusMSTW->Add(hEtaPlusMSTWpp,hEtaPlusMSTWnn,1./0.501973,1./0.460205);
  hEtaPlusMSTW->Add(hEtaPlusMSTWnppn,1./0.475634);
  hEtaMinusMSTW->Add(hEtaMinusMSTWpp,hEtaMinusMSTWnn,1./0.501973,1./0.460205);
  hEtaMinusMSTW->Add(hEtaMinusMSTWnppn,1./0.475634);
*/
  double nWPlus_pp = hEtaPlusMSTWpp->Integral();
  double nWPlus_nppn = hEtaPlusMSTWnppn->Integral();
  double nWPlus_nn = hEtaPlusMSTWnn->Integral();
  double nWPlus_tot = nWPlus_pp+nWPlus_nppn+nWPlus_nn;

  double nWMinus_pp = hEtaMinusMSTWpp->Integral();
  double nWMinus_nppn = hEtaMinusMSTWnppn->Integral();
  double nWMinus_nn = hEtaMinusMSTWnn->Integral();
  double nWMinus_tot = nWMinus_pp+nWMinus_nppn+nWMinus_nn;
 
  hEtaPlusMSTW->Add(hEtaPlusMSTWpp,hEtaPlusMSTWnn);
  hEtaPlusMSTW->Add(hEtaPlusMSTWnppn);
  hEtaMinusMSTW->Add(hEtaMinusMSTWpp,hEtaMinusMSTWnn);
  hEtaMinusMSTW->Add(hEtaMinusMSTWnppn);

  std::cout << "Fraction of pp events: " << hEtaPlusMSTWpp->Integral()/hEtaPlusMSTW->Integral();
  std::cout << "Fraction of np(pn) events: " << hEtaPlusMSTWnppn->Integral()/hEtaPlusMSTW->Integral();
  std::cout << "Fraction of nn events: " << hEtaPlusMSTWnn->Integral()/hEtaPlusMSTW->Integral();

  TH1F* hAsymmMSTW = getAsymmetry(hEtaPlusMSTW,hEtaMinusMSTW);
  float nMSTWMinusTot=0.;
  for(int i =1; i<=hEtaMinusMSTW->GetNbinsX(); ++i){
     float bw = hEtaMinusMSTW->GetBinWidth(i);
     std::cout << "W- MSTW: " << hEtaMinusMSTW->GetBinContent(i)*bw << std::endl; 
     nMSTWMinusTot +=hEtaMinusMSTW->GetBinContent(i)*bw;
  }
  std::cout << "nMSTWMinusTot: " << nMSTWMinusTot << std::endl;
  
  Int_t ci;   // for color index setting
  gStyle->SetLineScalePS(3.0);
  hAsymm_pp->SetFillColor(kGreen+1); hAsymm_pp->SetFillStyle(3006); 
    hAsymm_pp->SetMarkerColor(kGreen); hAsymm_pp->SetMarkerStyle(2);  hAsymm_pp->SetMarkerSize(1.0);  hAsymm_pp->SetLineWidth(1);
  hAsymm_crossNucleon->SetFillColor(kViolet+1); hAsymm_crossNucleon->SetFillStyle(3011); hAsymm_crossNucleon->SetMarkerSize(1.0); 
    hAsymm_crossNucleon->SetMarkerColor(kViolet); hAsymm_crossNucleon->SetMarkerStyle(2); hAsymm_crossNucleon->SetLineWidth(1);
  hAsymm_nn->SetFillColor(kCyan+2); hAsymm_nn->SetFillStyle(3013); hAsymm_nn->SetLineWidth(1);
    hAsymm_nn->SetMarkerColor(kCyan+1); hAsymm_nn->SetMarkerStyle(2); hAsymm_nn->SetMarkerSize(1.0);
  hAsymmMcAll->SetMarkerColor(kBlue); hAsymmMcAll->SetFillStyle(3005); hAsymmMcAll->SetLineWidth(1);
    hAsymmMcAll->SetMarkerStyle(27); hAsymmMcAll->SetFillColor(kBlue); hAsymmMcAll->SetMarkerSize(1.9);

  hEtaPlusMSTWpp->SetFillColor(kGreen+1); hEtaPlusMSTWpp->SetFillStyle(3006); 
    hEtaPlusMSTWpp->SetMarkerColor(kGreen); hEtaPlusMSTWpp->SetMarkerStyle(2);  hEtaPlusMSTWpp->SetMarkerSize(1.0);  hEtaPlusMSTWpp->SetLineWidth(1);
  hEtaPlusMSTWnppn->SetFillColor(kViolet+1); hEtaPlusMSTWnppn->SetFillStyle(3011); hEtaPlusMSTWnppn->SetMarkerSize(1.0); 
    hEtaPlusMSTWnppn->SetMarkerColor(kViolet); hEtaPlusMSTWnppn->SetMarkerStyle(2); hEtaPlusMSTWnppn->SetLineWidth(1);
  hEtaPlusMSTWnn->SetFillColor(kCyan+2); hEtaPlusMSTWnn->SetFillStyle(3013); hEtaPlusMSTWnn->SetLineWidth(1);
    hEtaPlusMSTWnn->SetMarkerColor(kCyan+1); hEtaPlusMSTWnn->SetMarkerStyle(2); hEtaPlusMSTWnn->SetMarkerSize(1.0);
  hEtaPlusMSTW->SetMarkerColor(kRed); hEtaPlusMSTW->SetFillStyle(3005); hEtaPlusMSTW->SetLineWidth(1);
    hEtaPlusMSTW->SetMarkerStyle(30); hEtaPlusMSTW->SetFillColor(kRed); hEtaPlusMSTW->SetMarkerSize(1.7);

  hEtaMinusMSTWpp->SetFillColor(kGreen+1); hEtaMinusMSTWpp->SetFillStyle(3006); 
    hEtaMinusMSTWpp->SetMarkerColor(kGreen); hEtaMinusMSTWpp->SetMarkerStyle(2);  hEtaMinusMSTWpp->SetMarkerSize(1.0);  hEtaMinusMSTWpp->SetLineWidth(1);
  hEtaMinusMSTWnppn->SetFillColor(kViolet+1); hEtaMinusMSTWnppn->SetFillStyle(3011); hEtaMinusMSTWnppn->SetMarkerSize(1.0); 
    hEtaMinusMSTWnppn->SetMarkerColor(kViolet); hEtaMinusMSTWnppn->SetMarkerStyle(2); hEtaMinusMSTWnppn->SetLineWidth(1);
  hEtaMinusMSTWnn->SetFillColor(kCyan+2); hEtaMinusMSTWnn->SetFillStyle(3013); hEtaMinusMSTWnn->SetLineWidth(1);
    hEtaMinusMSTWnn->SetMarkerColor(kCyan+1); hEtaMinusMSTWnn->SetMarkerStyle(2); hEtaMinusMSTWnn->SetMarkerSize(1.0);
  hEtaMinusMSTW->SetMarkerColor(kBlue); hEtaMinusMSTW->SetFillStyle(3005); hEtaMinusMSTW->SetLineWidth(1);
    hEtaMinusMSTW->SetMarkerStyle(30); hEtaMinusMSTW->SetFillColor(kBlue); hEtaMinusMSTW->SetMarkerSize(1.3);

  hAsymmMSTW->SetMarkerColor(kRed); hAsymmMSTW->SetFillStyle(3006); hAsymmMSTW->SetLineWidth(1);
    hAsymmMSTW->SetMarkerStyle(33); hAsymmMSTW->SetFillColor(kRed); hAsymmMSTW->SetMarkerSize(1.9);

/*  TCanvas *c_ppWsymm = new TCanvas("c_ppWsymm","c_ppWsymm",600,600);
  hAsymmMSTWpp->Draw("pe");
  hAsymm_pp->Draw("pesame");
  c_ppWsymm->Print("ppAsymmMC.root");
*/
  TF1* f0 = new TF1("f0","[0]",0,3);
  TF1* f1 = new TF1("f1","[0]+[1]*x",0,3);
  f0->SetLineStyle(kDashed);
  f1->SetLineStyle(kDashDotted);

  /// Fill theory graphs
  TGraphErrors* grWep   = new TGraphErrors(12);
  TGraphErrors* grNoWep = new TGraphErrors(12);
  if(doTheory){
  grWep->SetLineColor(kAzure-9); grWep->SetFillColor(kAzure-9);
  grNoWep->SetMarkerStyle(kOpenCircle); 
  grNoWep->SetFillColor(kAzure-9);
  for (int i = 0; i < 12; ++i){
    grWep->SetPoint(i,y[i]-0.02,(WepPlus[i]-WepMinus[i])/(WepPlus[i]+WepMinus[i]));
    grWep->SetPointError(i,0,sqrt(WepPlusErr[i]*WepPlusErr[i]+WepMinusErr[i]*WepMinusErr[i])/(WepPlus[i]+WepMinus[i]));
    grNoWep->SetPoint(i,y[i]-0.02,(noWepPlus[i]-noWepMinus[i])/(noWepPlus[i]+noWepMinus[i]));
    grNoWep->SetPointError(i,0,sqrt(noWepPlusErr[i]*noWepPlusErr[i]+noWepMinusErr[i]*noWepMinusErr[i])/(noWepPlus[i]+noWepMinus[i]));
  }
  }

  TCanvas c2 = TCanvas("c2","c2",600,600);
/*  TF1* f2 = new TF1("f2","[0]",0,3);
  TF1* f3 = new TF1("f3","[0]+[1]*x",0,3);
  f2->SetLineStyle(kDashed);
  f3->SetLineStyle(kDashDotted);
*/
     //TCanvas* c3 = new TCanvas("c3","c3",600,600);

//   for (int i=0;i<1;i++){
     TCanvas* c3 = new TCanvas("c3","c3",600,600);
     std::cout << "constructing charge asymmetry... " << std::endl; 
     TGaxis::SetMaxDigits(2);
     hAsymm[0]->SetMarkerStyle(markerStyle[0]); 
     hAsymm[0]->GetYaxis()->SetLabelSize(0.04);; 
     hAsymm[0]->GetXaxis()->SetLabelSize(0.04);; 
     hAsymm[0]->SetMarkerColor(markerColor[0]) ;
     hAsymm[0]->SetLineColor(markerColor[0]) ;
     hAsymm[0]->SetMarkerSize(1.2);

     hAsymm[0]->SetXTitle("|#eta_{#mu}|") ;
     //hAsymm[0]->SetYTitle("W#rightarrow#mu Charge asymmetry");
     hAsymm[0]->SetYTitle("Muon Charge Asymmetry A_{#mu}");
     hAsymm[0]->GetYaxis()->SetNoExponent(kTRUE);

     //hAsymm[0]->SetMaximum(0.9); hAsymm[0]->SetMinimum(-0.9);
     hAsymm[0]->GetXaxis()->SetRangeUser(0.1,2.4);
     hAsymm[0]->GetYaxis()->SetRangeUser(-0.4,0.21);

     hAsymm[0]->Draw("0") ;

     TF1* fd = new TF1("fd","[0]",0,3);
     fd->SetLineStyle(kDashed);
     fd->SetLineColor(kGray+2);
     fd->Draw("lsame");
     grassSyst[0]->SetFillColor(38);
     //gStyle->SetHatchesLineWidth(1);
     grassSyst[0]->SetFillStyle(3004);
     grassSyst[0]->SetMarkerStyle(20);
     ///uncomment to draw systematics
     grassSyst[0]->Draw("2same");
     if (doTheory){
       grNoWep->Draw("3same");
     }

    //hAsymmMcAll->Draw("pe2same");
    //hAsymmMSTW->Draw("pe2same");
   TGraphAsymmErrors *grAsymmMcAll = new TGraphAsymmErrors(hAsymmMcAll->GetNbinsX());
   grAsymmMcAll = convertToGraph(hAsymmMcAll);
   ci = TColor::GetColor("#cc0000");
   grAsymmMcAll->SetLineColor(ci);
   grAsymmMcAll->SetLineStyle(2);
   grAsymmMcAll->SetLineWidth(4);
   grAsymmMcAll->Draw("lsame");

   TGraphAsymmErrors *grAsymmMSTW = new TGraphAsymmErrors(hAsymmMSTW->GetNbinsX());
   grAsymmMSTW = convertToGraph(hAsymmMSTW);
   ci = TColor::GetColor("#0000ff");
   grAsymmMSTW->SetLineColor(ci);
   grAsymmMSTW->SetLineStyle(9);
   grAsymmMSTW->SetLineWidth(4);
   grAsymmMSTW->Draw("lsame");


    // Uncomment to draw individual nucleon asymm from PYTHIA MRST LO*
/*    hAsymm_pp->Draw("pe2same");
    hAsymm_crossNucleon->Draw("pe2same");
    hAsymm_nn->Draw("pe2same");
    */
    TGraphAsymmErrors* grassDummy = (TGraphAsymmErrors*)grass->Clone("grassDummy");
    grassDummy->SetMarkerColor(kWhite); grassDummy->SetMarkerSize(1.9);
//    grassDummy->Draw("pesame");
    grass->SetLineWidth(3);
    grass->Draw("pesame");

    c3->RedrawAxis("g");
    hAsymm[0]->GetXaxis()->SetRangeUser(0.11,2.39);

   TLegend *leg = new TLegend(0.1661074,0.1713287,0.5352349,0.4073427,NULL,"brNDC");
   leg->SetTextFont(gStyle->GetTextFont());
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03846154);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);

    ///save files to outfile 
    outFile->cd();
    grass->Write("grChargeAsymm");
    leg->AddEntry(grassSyst[0],"Data 2011","PEF");
    //leg->AddEntry(grassSyst[0],"MC","PEF");
    if(doTheory) leg->AddEntry(grNoWep,"Salgado","f");
    leg->AddEntry(grAsymmMcAll,"MRST LO*","l");
    leg->AddEntry(grAsymmMSTW,"MSTW NLO","l");
    /*leg->AddEntry(hAsymmMcAll,"PYTHIA (All) (MRST LO*)","pef");
    leg->AddEntry(hAsymm_pp,"pp (MRST LO*)","pef");
    leg->AddEntry(hAsymm_crossNucleon,"np(pn) (MRST LO*)","pef");
    leg->AddEntry(hAsymm_nn,"nn (MRST LO*)","pef");
    */
   TH1F* hAsymmDummyStat = new TH1F("hAsymmDummyStat","hAsymmDummyStat",bins,xBins);
   hAsymmDummyStat->SetMarkerStyle(2);
   leg->Draw() ;
   leg = new TLegend(0.4563758,0.2062937,0.7466443,0.3811189,NULL,"brNDC");
   leg->SetTextFont(gStyle->GetTextFont());
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03846154);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->AddEntry(hAsymmDummyStat,"Stat. uncertainty","pe");
   leg->AddEntry(grassSyst[0],"Total uncertainty","f");
   leg->Draw();

    TLatex *   tex = new TLatex(0.1828859,0.4335664,"#int Ldt #approx 0.14 nb^{-1}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03846154);
   tex->SetLineWidth(2);
   tex->Draw();
   tex = new TLatex(0.4731544,0.4370629,"#sqrt{s_{NN}} = 2.76 TeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04195804);
   tex->SetLineWidth(2);
   tex->Draw();
   tex = new TLatex(0.6073826,0.8566434,"ATLAS Internal");
tex->SetNDC();
   tex->SetTextFont(72);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.4,0.5,"0-80%");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04195804);
   tex->SetLineWidth(2);
   tex->Draw();
    hAsymm[0]->GetXaxis()->SetRangeUser(0.11,2.39);
    c3->Update();
//    c3->Print("chargeAsymmetry"/*+fileNameDataIn+*/".png"); 
//    c3->Print("chargeAsymmetry"/*+fileNameDataIn+*/".eps"); 
//    c3->Print("chargeAsymmetry"/*+fileNameDataIn+*/".pdf"); 
    
    c3->Print("chargeAsymmetry"/*+fileNameDataIn+*/".root"); 

    
//  }

  ///plot the eta distribution
  ///number of W^{+}->mu^{+} and W^{-}->mu^{-}

/*      TH1D* hpMCtemp  = new TH1D("hpMCtemp","hpMCtemp",nEtaBins,arrBinning);
      hpMCtemp->Sumw2();
      TH1D* hmMCtemp = new TH1D("hmMCtemp","hmMCtemp",nEtaBins,arrBinning);
      hpMCtemp->Sumw2();
*/
      c2.cd(); 
      //treeMC->Draw("abs(mc_mu_gen_eta)>>hpMCtemp","mc_mu_charge==1&&mc_mu_gen_mothertype==+24"+cutsMC,"PE");
      //treeMC->Draw("abs(mc_mu_gen_eta)>>hmMCtemp","mc_mu_charge==-1&&mc_mu_gen_mothertype==-24"+cutsMC,"PE");


      //int nMCEvents = hpMCtemp->GetEntries()+hmMCtemp->GetEntries();
 

      TGraphAsymmErrors* grWpUncorrelatedSystc = (TGraphAsymmErrors*)grWpUncorrelatedSyst->Clone("grWpUncorrelatedSystc");
      TGraphAsymmErrors* grWmUncorrelatedSystc = (TGraphAsymmErrors*)grWmUncorrelatedSyst->Clone("grWmUncorrelatedSystc");
      TGraphAsymmErrors* grWpCorrelatedSystc = (TGraphAsymmErrors*)grWpCorrelatedSyst->Clone("grWpCorrelatedSystc");
      TGraphAsymmErrors* grWmCorrelatedSystc = (TGraphAsymmErrors*)grWmCorrelatedSyst->Clone("grWmCorrelatedSystc");
      TGraphAsymmErrors* grWpc = (TGraphAsymmErrors*)grWp->Clone("grWpc");
      TGraphAsymmErrors* grWmc = (TGraphAsymmErrors*)grWm->Clone("grWmc");

      grWpc->GetYaxis()->SetRangeUser(0.0,10.0); /*SetMaximum(0.13); grWpUncorrelatedSystc->SetMinimum(0.0);*/
      grWpc->GetYaxis()->SetTitle("#frac{dN}{d#eta} #frac{1}{N_{events} #LT N_{coll} #GT} 10^{9}"); /*SetMaximum(0.13); grWpUncorrelatedSystc->SetMinimum(0.0);*/
      grWpc->GetXaxis()->SetRangeUser(0.1,2.4);
      grWpc->GetXaxis()->SetTitle("|#eta|"); /*SetMaximum(0.13); grWpUncorrelatedSystc->SetMinimum(0.0);*/

      grWpUncorrelatedSystc->SetLineColor(markerColor[2]);
      grWpUncorrelatedSystc->SetFillColor(kRed-8);
      grWmUncorrelatedSystc->SetLineColor(markerColor[1]);
      grWmUncorrelatedSystc->SetFillColor(kBlue-8);

      grWpCorrelatedSystc->SetLineColor(markerColor[2]);
      grWpCorrelatedSystc->SetFillColor(kRed);
      grWpCorrelatedSystc->SetFillStyle(3004);
      grWmCorrelatedSystc->SetLineColor(markerColor[1]);
      grWmCorrelatedSystc->SetFillColor(kBlue);
      grWmCorrelatedSystc->SetFillStyle(3005);

///uncomment to draw syst. uncert
     TCanvas *c4 = new TCanvas("c4","c4",600,600);

     // Already scaled by eta binwidth
     // Get cross-section from yields in the data
     double nEv = 1.03e9;
     double ncoll_0to80 = 452.0;
     double scaleFactor = 1.0e9/(0.8*nEv*ncoll_0to80);
     grWpc = scaleToBinaryColl(grWpc,scaleFactor);
     grWpCorrelatedSystc = scaleToBinaryColl(grWpCorrelatedSystc,scaleFactor,true);
     grWpUncorrelatedSystc = scaleToBinaryColl(grWpUncorrelatedSystc,scaleFactor);
     grWmc = scaleToBinaryColl(grWmc,scaleFactor);
     TGraphAsymmErrors* grWmDummy = (TGraphAsymmErrors*)grWmc->Clone("grWmDummy");
     grWmDummy->SetMarkerColor(kWhite); grWmDummy->SetMarkerSize(1.8); grWmDummy->SetMarkerStyle(markerStyle[2]);
     TGraphAsymmErrors* grWpDummy = (TGraphAsymmErrors*)grWpc->Clone("grWpDummy");
     grWpDummy->SetMarkerColor(kWhite); grWpDummy->SetMarkerSize(1.8);grWpDummy->SetMarkerStyle(markerStyle[3]);
     grWmCorrelatedSystc = scaleToBinaryColl(grWmCorrelatedSystc,scaleFactor,true);
     grWmUncorrelatedSystc = scaleToBinaryColl(grWmUncorrelatedSystc,scaleFactor);

     grWpc->GetXaxis()->SetRangeUser(0.1,2.4);
     grWpc->GetYaxis()->SetRangeUser(0.,10.);
     grWpc->GetXaxis()->SetLabelSize(0.04);
     grWpc->GetYaxis()->SetLabelSize(0.04);
     grWpc->GetYaxis()->SetTitleSize(0.04);
     grWpc->GetXaxis()->SetTitle("|#eta|");
     grWpc->GetYaxis()->SetTitle("#frac{dN}{d#eta} #frac{1}{N_{events} } #frac{10^{9}}{#LT N_{coll} #GT}");
     grWpc->Draw("ape");
     grWpUncorrelatedSystc->Draw("e2same");
     grWmUncorrelatedSystc->Draw("e2same");
//      grWp->Draw("pesame");
     grWpCorrelatedSystc->Draw("e2same");
     grWmCorrelatedSystc->Draw("e2same");
     grWpc->Draw("pesame");
     grWmc->Draw("pesame");
      //std::cout << "mu^{+}: " << grWp->Integral(1,nEtaBins) << " mu^{-}: " << hmMCtemp->GetEntries() << " MC entries: " << nMCEvents << std::endl;
      //std::cout << yWPlusTotal << std::endl;
      //std::cout << yWMinusTotal << std::endl;
     // Scale by cross-section
     hEtaPlusMcAll->Draw("pesame") ;
     hEtaMinusMcAll->Draw("pesame") ;

     hEtaPlusMSTW->Draw("pesame") ;
     hEtaMinusMSTW->Draw("pesame") ;


     //std::cout << "constructing W^{+}->mu^{+} and W^{-}->mu^{-} eta distribution for centrality bin: " << i << std::endl; 
     hpAsymm[0]->SetMarkerStyle(markerStyle[3]); 
     hpAsymm[0]->SetMarkerColor(markerColor[2]) ;
     hpAsymm[0]->SetLineColor(markerColor[2]) ;
     hpAsymm[0]->SetMarkerSize(2);

     hmAsymm[0]->SetMarkerStyle(markerStyle[2]); 
     hmAsymm[0]->SetMarkerColor(markerColor[1]) ;
     hmAsymm[0]->SetLineColor(markerColor[1]) ;
     hmAsymm[0]->SetMarkerSize(2);

     //std::cout << "constructing Monte Carlo W^{+}->mu^{+} and W^{-}->mu^{-} eta distribution for centrality bin: " << i << std::endl; 
     hEtaPlusMcAll->SetMarkerStyle(30); 
     hEtaPlusMcAll->SetMarkerColor(kRed) ;
     hEtaPlusMcAll->SetMarkerSize(1.9);

     hEtaPlusMSTW->SetMarkerStyle(27); 
     hEtaPlusMSTW->SetMarkerColor(kRed) ;
     hEtaPlusMSTW->SetMarkerSize(1.9);

     hEtaMinusMcAll->SetMarkerStyle(30); 
     hEtaMinusMcAll->SetMarkerColor(kBlue) ;
     hEtaMinusMcAll->SetMarkerSize(1.9);

     hEtaMinusMSTW->SetMarkerStyle(27); 
     hEtaMinusMSTW->SetMarkerColor(kBlue) ;
     hEtaMinusMSTW->SetMarkerSize(1.9);

     grWpc->SetMarkerStyle(markerStyle[3]); 
     grWpUncorrelatedSystc->SetMarkerStyle(markerStyle[3]);
     grWpUncorrelatedSystc->SetMarkerColor(markerColor[2]) ;
     grWpc->SetMarkerColor(markerColor[2]) ;
     grWpc->SetLineColor(markerColor[2]) ;
     grWpc->SetMarkerSize(1.5);

     grWmc->SetMarkerStyle(markerStyle[2]); 
     grWmc->SetMarkerColor(markerColor[1]) ;
     grWmc->SetLineColor(markerColor[1]) ;
     grWmc->SetMarkerSize(1.5);
     grWmUncorrelatedSystc->SetMarkerStyle(markerStyle[2]);
     grWmUncorrelatedSystc->SetMarkerColor(markerColor[1]) ;

     //hpAsymm[0]->Draw() ;

     //grWp[0]->GetXaxis()->SetXTitle("|#eta_{#mu}|") ;
     //grWp[0]->GetYaxis()->SetYTitle("Normalized W#rightarrow#mu Events");

    TGaxis::SetMaxDigits(2);
   leg = new TLegend(0.1795302,0.1905594,0.4161074,0.3234266,NULL,"brNDC");
   leg->SetTextFont(gStyle->GetTextFont());
   leg->SetBorderSize(0);
   leg->SetTextSize(0.02622377);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->AddEntry(grWpUncorrelatedSystc, "W^{+}#rightarrow#mu^{+} Data", "pef");
   leg->AddEntry(grWmUncorrelatedSystc, "W^{-}#rightarrow#mu^{-} Data", "pef");
    leg->Draw();

   leg = new TLegend(0.6577181,0.1730769,0.8104027,0.3391608,NULL,"brNDC");
   leg->SetTextFont(gStyle->GetTextFont());
   leg->SetBorderSize(0);
   leg->SetTextSize(0.02622377);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->AddEntry(hEtaPlusMcAll,"W^{+}#rightarrow#mu^{+} MRST LO*","p");
   leg->AddEntry(hEtaMinusMcAll,"W^{-}#rightarrow#mu^{-} MRST LO*","p");
   leg->Draw();
   leg = new TLegend(0.3959732,0.1713287,0.5486577,0.3374126,NULL,"brNDC");
   leg->SetTextFont(gStyle->GetTextFont());
   leg->SetBorderSize(0);
   leg->SetTextSize(0.02622377);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);

    leg->AddEntry(hEtaPlusMSTW,"W^{+}#rightarrow#mu^{+} MSTW NLO","p");
    leg->AddEntry(hEtaMinusMSTW,"W^{-}#rightarrow#mu^{-} MSTW NLO","p");
    leg->Draw() ;
    tex = new TLatex(0.1862416,0.3653846,"#int Ldt #approx 0.140 nb^{-1}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03846154);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.4899329,0.3671329,"#sqrt{s_{NN}} = 2.76 TeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04195804);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.6073826,0.8566434,"ATLAS");
tex->SetNDC();
   tex->SetTextFont(72);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.7701342,0.8601399,"Internal");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();


    c4->Update();
//    c4->Print("chargeEtaDistr"/*fileNameDataIn+*/".png"); 
//    c4->Print("chargeEtaDistr"/*fileNameDataIn+*/".eps"); 
//    c4->Print("chargeEtaDistr"/*fileNameDataIn+*/".root"); 
//    c4->Print("chargeEtaDistr"/*fileNameDataIn+*/".pdf"); 
    ///write to outfile
    fHisto->cd();
    hAssSyst->Write();
    outFile->cd();
/*    grWpUncorrelatedSystc->Write();
    grWmUncorrelatedSystc->Write();
*/
    grWpc->Write();
    grWmc->Write();
    //hpMCtemp->Write();
    //hmMCtemp->Write();
    /*ATLAS_LABEL(0.6,0.25);
    myText(0.6,0.20,1,"Internal");
    myText(0.2,0.2,1,"#int L #approx 0.140 nb^{-1}");   
*/

    //c4->Print("chargeDistr"+sSel+fileNameDataIn+".root"); 

   TH1F* hEtaPlusData = new TH1F("hEtaPlusData","hEtaPlusData",bins,xBins);
   TH1F* hRatioMSTWPlusUncorrelated = new TH1F("hRatioMSTWPlusUncorrelated","hRatioMSTWPlusUncorrelated",bins,xBins);
   TH1F* hRatioMSTWPlusCorrelated = new TH1F("hRatioMSTWPlusCorrelated","hRatioMSTWPlusCorrelated",bins,xBins);
   TH1F* hRatioCT10PlusUncorrelated = new TH1F("hRatioCT10PlusUncorrelated","hRatioCT10PlusUncorrelated",bins,xBins);

   hRatioMSTWPlusCorrelated->SetFillColor(kRed-8);
   //TGraphErrors* grRatioPlusCorrelated = new TGraphErrors(bins);
   ///Convert TGraph to TH1F
   for(int i=0; i<grWpc->GetN(); ++i){

      float xpt = grWpc->GetX()[i];
      float ypt = grWpc->GetY()[i];
      float errUncorr = grWpUncorrelatedSystc->GetEYhigh()[i];
      float errCorr = grWpCorrelatedSystc->GetEYhigh()[i];
      float errX = grWmCorrelatedSystc->GetEXhigh()[i];
      //MSTW
      float yMod = hEtaPlusMSTW->GetBinContent(i+1);
      float ratio = ypt/yMod;
      //CT10
      float yCT10 = hEtaPlusCT10->GetBinContent(i+1);
      float ratioCT10 = ypt/yCT10;
//      std::cout << "ratio : " << ratio << std::endl;
//      std::cout << "error: " << errCorr/ypt*ratio << std::endl;
      hRatioMSTWPlusCorrelated->SetBinContent(i+1,ratio);
      hRatioMSTWPlusCorrelated->SetBinError(i+1,errCorr/ypt*ratio);
      hEtaPlusData->SetBinContent(i+1,ypt);
      hEtaPlusData->SetBinError(i+1,grWpc->GetEYhigh()[i]);
      hRatioMSTWPlusUncorrelated->SetBinContent(i+1,ratio);
      hRatioMSTWPlusUncorrelated->SetBinError(i+1,errUncorr/ypt*ratio);
      hRatioCT10PlusUncorrelated->SetBinContent(i+1,ratioCT10);
      hRatioCT10PlusUncorrelated->SetBinError(i+1,errUncorr/ypt*ratioCT10);
   }


   TCanvas* cEtaPlus = new TCanvas("cEtaPlus","cEtaPlus",600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cEtaPlus->Range(-0.3658228,-1.043038,2.54557,0.4759494);
   cEtaPlus->SetFillColor(0);
   cEtaPlus->SetBorderMode(0);
   cEtaPlus->SetBorderSize(2);
   cEtaPlus->SetTickx(1);
   cEtaPlus->SetTicky(1);
   cEtaPlus->SetLeftMargin(0.16);
   cEtaPlus->SetRightMargin(0.05);
   cEtaPlus->SetTopMargin(0.05);
   cEtaPlus->SetBottomMargin(0.16);
   cEtaPlus->SetFrameBorderMode(0);
   cEtaPlus->SetFrameBorderMode(0);

   cEtaPlus->Divide(1, 2);
   TPad* canvas_up = (TPad*)cEtaPlus->GetListOfPrimitives()->FindObject("cEtaPlus_1");
   TPad* canvas_dw = (TPad*)cEtaPlus->GetListOfPrimitives()->FindObject("cEtaPlus_2");
    ///Define the size
   double  up_height     = 0.8; 
   double  dw_correction = 1.60;
   double  font_size_dw  = 0.08;
   double  dw_height    = (1. - up_height) * dw_correction;

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


   ///MC/Data histo
   // MSTW
   TH1F* hRatioMSTWPlus = (TH1F*)hEtaPlusData->Clone("hRatioMSTWPlus");
   hRatioMSTWPlus = getMCDataRatio(hRatioMSTWPlus,hEtaPlusData,hEtaPlusMSTW);
   // CT10 
   //TH1F* hRatioCT10Plus = (TH1F*)hEtaPlusData->Clone("hRatioCT10Plus");
   //hRatioCT10Plus = getMCDataRatio(hRatioCT10Plus,hEtaPlusData,hEtaPlusCT10);

   ///Draw eta distros 
   canvas_up->cd();

   // Data
   TH1F* hEtaPlusDummy = new TH1F("hEtaPlusDummy","hEtaPlusDummy",nEtaBins,arrBinning);
   hEtaPlusDummy->GetXaxis()->SetRangeUser(0.11,2.39);
   hEtaPlusDummy->GetXaxis()->SetTitle("|#eta|");
   hEtaPlusDummy->GetYaxis()->SetRangeUser(0.01,10.);
   hEtaPlusDummy->GetYaxis()->SetTitle("#frac{dN}{d#eta} #frac{1}{N_{events} } #frac{10^{9}}{#LT N_{coll} #GT}");
   hEtaPlusDummy->GetXaxis()->SetLabelSize(0.04);
   hEtaPlusDummy->GetYaxis()->SetTitleSize(0.04);
   hEtaPlusDummy->GetYaxis()->SetLabelSize(0.04);
   hEtaPlusDummy->Draw();

   grWpUncorrelatedSystc->Draw("e2same");
   grWpCorrelatedSystc->Draw("e2same");
   grWpc->Draw("pesame");

   // Get ratio of model/data

   // Monte Carlo
   //hEtaPlusMSTWnn->Draw("pesame");
   //hEtaPlusMSTWpp->Draw("pesame");
   //hEtaPlusMSTWnppn->Draw("pesame");
   //hEtaPlusMSTW->Draw("pesame");
   //hEtaPlusMcAll->Draw("pesame") ;

   TGraphAsymmErrors *grEtaPlusMcAll = new TGraphAsymmErrors(hEtaPlusMcAll->GetNbinsX());
   grEtaPlusMcAll = convertToGraph(hEtaPlusMcAll);
   ci = TColor::GetColor("#cc0000");
   grEtaPlusMcAll->SetLineColor(ci);
   grEtaPlusMcAll->SetLineStyle(2);
   grEtaPlusMcAll->SetLineWidth(4);
   grEtaPlusMcAll->Draw("lsame");

   TGraphAsymmErrors *grEtaPlusMSTW = new TGraphAsymmErrors(hEtaPlusMSTW->GetNbinsX());
   grEtaPlusMSTW = convertToGraph(hEtaPlusMSTW);
   ci = TColor::GetColor("#cc0000");
   grEtaPlusMSTW->SetLineColor(ci);
   grEtaPlusMSTW->SetLineStyle(9);
   grEtaPlusMSTW->SetLineWidth(4);
   grEtaPlusMSTW->Draw("lsame");

   TGraphAsymmErrors *grEtaPlusCT10 = new TGraphAsymmErrors(hEtaPlusCT10->GetNbinsX());
   grEtaPlusCT10 = convertToGraph(hEtaPlusCT10);
   ci = TColor::GetColor("#cc0000");
   grEtaPlusCT10->SetLineColor(ci);
   grEtaPlusCT10->SetLineStyle(7);
   grEtaPlusCT10->SetLineWidth(4);
   grEtaPlusCT10->Draw("lsame");

   grWpDummy->Draw("psame");
   grWpc->Draw("pesame");
   hEtaPlusDummy->Draw("sameaxis");
   leg = new TLegend(0.1694631,0.2023601,0.4832215,0.3771853,NULL,"brNDC");
   leg->SetTextFont(gStyle->GetTextFont());
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04545455);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->AddEntry(grWpc,"Data","pe");
   /*leg->AddEntry(hEtaPlusMSTWnn,"nn","pe");
   leg->AddEntry(hEtaPlusMSTWnppn,"np(pn)","pe");
   leg->AddEntry(hEtaPlusMSTWpp,"pp","pe");
   */
   leg->AddEntry(grEtaPlusMSTW,"MSTW NLO","l");
   leg->AddEntry(grEtaPlusCT10,"CT10 NLO","l");
   leg->AddEntry(grEtaPlusMcAll,"MRST LO*","l");
   leg->Draw();
   leg = new TLegend(0.4463087,0.1958042,0.7600671,0.3706294,NULL,"brNDC");
   leg->SetTextFont(gStyle->GetTextFont());
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04545455);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->AddEntry(grWpUncorrelatedSystc,"Uncorr. uncertainty","f");
   leg->AddEntry(grWpCorrelatedSystc,"Corr. uncertainty","f");
   leg->Draw();
   tex = new TLatex(0.6879195,0.7596154,"W^{+}#rightarrow#mu^{+}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.06993007);
   tex->SetLineWidth(2);
   tex->Draw();
//      tex = new TLatex(0.4211409,0.3129371,"POW+PY8 NLO");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
   tex = new TLatex(0.1862416,0.4318182,"#int Ldt #approx 0.14 nb^{-1}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03846154);
   tex->SetLineWidth(2);
   tex->Draw();
   tex = new TLatex(0.4060403,0.4383741,"#sqrt{s_{NN}} = 2.76 TeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03846154);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.6073826,0.8566434,"ATLAS Internal");
tex->SetNDC();
   tex->SetTextFont(72);
   tex->SetLineWidth(2);
   tex->Draw();
   tex = new TLatex(0.62,0.44,"0-80%");
tex->SetNDC();
tex->SetTextSize(0.03846154);
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();


    ///Draw MC/Data underneath
    canvas_dw->cd();
    /// font size
    hRatioMSTWPlusCorrelated->GetXaxis()->SetTitle("|#eta|");
    hRatioMSTWPlusCorrelated->GetXaxis()->SetLabelFont(42);
    hRatioMSTWPlusCorrelated->GetXaxis()->SetLabelSize(0.09);
    hRatioMSTWPlusCorrelated->GetXaxis()->SetTitleSize(0.11);
    hRatioMSTWPlusCorrelated->GetXaxis()->SetTitleOffset(0.8);
    hRatioMSTWPlusCorrelated->GetXaxis()->SetTitleFont(42);
    hRatioMSTWPlusCorrelated->GetYaxis()->SetTitle("Data/MSTW");
    hRatioMSTWPlusCorrelated->GetYaxis()->SetLabelFont(42);
    hRatioMSTWPlusCorrelated->GetYaxis()->SetLabelSize(0.09);
    hRatioMSTWPlusCorrelated->GetYaxis()->SetTitleSize(0.11);
    hRatioMSTWPlusCorrelated->GetYaxis()->SetTitleOffset(0.7);

   TF1* funcRatio = new TF1("funcRatio","[0]",0.0,2.5);
   funcRatio->SetLineStyle(kDashed);
   funcRatio->SetLineColor(kGray+2);
   funcRatio->SetParameter(0,1.0);
   funcRatio->FixParameter(0,1.0);

   hRatioMSTWPlusCorrelated->GetYaxis()->SetRangeUser(0.65,1.45);
   hRatioMSTWPlusCorrelated->GetXaxis()->SetRangeUser(0.11,2.39);
   //hRatioMSTWPlusCorrelated->SetFillColor(kGray+1);
   hRatioMSTWPlusUncorrelated->SetFillColor(kRed-8);

   hRatioMSTWPlusCorrelated->SetFillColor(kRed);
   hRatioMSTWPlusCorrelated->SetFillStyle(3004);

   hRatioMSTWPlusUncorrelated->SetMarkerSize(0);
   hRatioMSTWPlusUncorrelated->Draw("e2");
   hRatioMSTWPlusCorrelated->SetMarkerSize(0);
   hRatioMSTWPlusCorrelated->Draw("e2same");
   funcRatio->Draw("same");
   hRatioMSTWPlus->SetLineColor(kRed);
   hRatioMSTWPlus->SetMarkerStyle(23);
   hRatioMSTWPlus->SetMarkerColor(markerColor[2]);
   hRatioMSTWPlus->SetMarkerSize(1.5);
   TH1F* hRatioPlusDummy = (TH1F*)hRatioMSTWPlus->Clone("hRatioPlusDummy");
   hRatioPlusDummy->SetMarkerStyle(23);
   hRatioPlusDummy->SetMarkerSize(1.8);
   hRatioPlusDummy->SetMarkerColor(kWhite);
   hRatioPlusDummy->Draw("pesame");
   hRatioMSTWPlus->Draw("pesame");

   cEtaPlus->Print("etaPlus.root");
   cEtaPlus->Print("etaPlus.pdf");

   TH1F* hEtaMinusData = new TH1F("hEtaMinusData","hEtaMinusData",bins,xBins);
   TH1F* hRatioMSTWMinusCorrelated = new TH1F("hRatioMSTWMinusCorrelated","hRatioMSTWMinusCorrelated",bins,xBins);
   hRatioMSTWMinusCorrelated->SetFillColor(kBlue-8);
   TH1F* hRatioMSTWMinusUncorrelated = new TH1F("hRatioMSTWMinusUncorrelated","hRatioMSTWMinusUncorrelated",bins,xBins);
   TH1F* hRatioCT10MinusUncorrelated = new TH1F("hRatioCT10MinusUncorrelated","hRatioCT10MinusUncorrelated",bins,xBins);
   //TGraphErrors* grRatioMinusCorrelated = new TGraphErrors(bins);
   ///Convert TGraph to TH1F
   for(int i=0; i<grWmc->GetN(); ++i){

      float ypt = grWmc->GetY()[i];
      float errUncorr = grWmUncorrelatedSystc->GetEYhigh()[i];
      float errCorr = grWmCorrelatedSystc->GetEYhigh()[i];
      float errX = grWmCorrelatedSystc->GetEXhigh()[i];
      float yMod = hEtaMinusMSTW->GetBinContent(i+1);
      float ratio = ypt/yMod;
      //CT10
      float yCT10 = hEtaMinusCT10->GetBinContent(i+1);
      float ratioCT10 = ypt/yCT10;

      hRatioMSTWMinusCorrelated->SetBinContent(i+1,ratio);
      hRatioMSTWMinusCorrelated->SetBinError(i+1,errCorr/ypt*ratio);
      hEtaMinusData->SetBinContent(i+1,ypt);
      hEtaMinusData->SetBinError(i+1,grWmc->GetEYhigh()[i]);
      hRatioMSTWMinusUncorrelated->SetBinContent(i+1,ratio);
      hRatioMSTWMinusUncorrelated->SetBinError(i+1,errUncorr/ypt*ratio);
      hRatioCT10MinusUncorrelated->SetBinContent(i+1,ratioCT10);
      hRatioCT10MinusUncorrelated->SetBinError(i+1,errUncorr/ypt*ratioCT10);
   }

   TCanvas* cEtaMinus = new TCanvas("cEtaMinus","cEtaMinus",600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cEtaMinus->Range(-0.3658228,-1.043038,2.54557,0.4759494);
   cEtaMinus->SetFillColor(0);
   cEtaMinus->SetBorderMode(0);
   cEtaMinus->SetBorderSize(2);
   cEtaMinus->SetTickx(1);
   cEtaMinus->SetTicky(1);
   cEtaMinus->SetLeftMargin(0.16);
   cEtaMinus->SetRightMargin(0.05);
   cEtaMinus->SetTopMargin(0.05);
   cEtaMinus->SetBottomMargin(0.16);
   cEtaMinus->SetFrameBorderMode(0);
   cEtaMinus->SetFrameBorderMode(0);

   cEtaMinus->Divide(1, 2);
   canvas_up = (TPad*)cEtaMinus->GetListOfPrimitives()->FindObject("cEtaMinus_1");
   canvas_dw = (TPad*)cEtaMinus->GetListOfPrimitives()->FindObject("cEtaMinus_2");
    ///Define the size
    up_height     = 0.8; 
    dw_correction = 1.60;
    font_size_dw  = 0.08;
    dw_height    = (1. - up_height) * dw_correction;

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


   ///MC/Data histo
   TH1F* hRatioMSTWMinus = (TH1F*)hEtaMinusData->Clone("hRatioMSTWMinus");
   hRatioMSTWMinus = getMCDataRatio(hRatioMSTWMinus,hEtaMinusData,hEtaMinusMSTW);

   ///Draw eta distros 
   canvas_up->cd();

   // Data
   TH1F* hEtaMinusDummy = new TH1F("hEtaMinusDummy","hEtaMinusDummy",nEtaBins,arrBinning);
   hEtaMinusDummy->GetXaxis()->SetRangeUser(0.11,2.39);
   hEtaMinusDummy->GetXaxis()->SetTitle("|#eta|");
   hEtaMinusDummy->GetYaxis()->SetRangeUser(0.01,10.);
   hEtaMinusDummy->GetYaxis()->SetTitle("#frac{dN}{d#eta} #frac{1}{N_{events} } #frac{10^{9}}{#LT N_{coll} #GT}");
   hEtaMinusDummy->GetXaxis()->SetLabelSize(0.04);
   hEtaMinusDummy->GetYaxis()->SetTitleSize(0.04);
   hEtaMinusDummy->GetYaxis()->SetLabelSize(0.04);
   hEtaMinusDummy->Draw();
/*        float ncounts = 0;
   for (int i =0; i<grWmc->GetN(); ++i){
        float bw = 2.*grWmc->GetEXhigh()[i];
        std::cout << "bw: " << bw << std::endl;
        float ns = grWmc->GetY()[i]*bw;
        std::cout << "number of counts in bin " << i << ": " << ns << std::endl;
        ncounts +=ns;
   }
   std::cout << "total: " << ncounts/0.536 << std::endl; exit(0);
*/   
   grWmUncorrelatedSystc->Draw("e2same");
   grWmCorrelatedSystc->Draw("e2same");
   grWmDummy->Draw("psame");
   grWmc->Draw("pesame");
   // Monte Carlo
   /*hEtaMinusMSTWnn->Draw("pesame");
   hEtaMinusMSTWpp->Draw("pesame");
   hEtaMinusMSTWnppn->Draw("pesame");
   */
   //hEtaMinusMSTW->Draw("pesame");
   //hEtaMinusMcAll->Draw("pesame") ;
   TGraphAsymmErrors *grEtaMinusMcAll = new TGraphAsymmErrors(hEtaMinusMcAll->GetNbinsX());
   grEtaMinusMcAll = convertToGraph(hEtaMinusMcAll);
   ci = TColor::GetColor("#0000ff");
   grEtaMinusMcAll->SetLineColor(ci);
   grEtaMinusMcAll->SetLineStyle(2);
   grEtaMinusMcAll->SetLineWidth(4);
   grEtaMinusMcAll->Draw("lsame");

   TGraphAsymmErrors *grEtaMinusMSTW = new TGraphAsymmErrors(hEtaMinusMSTW->GetNbinsX());
   grEtaMinusMSTW = convertToGraph(hEtaMinusMSTW);
   ci = TColor::GetColor("#0000ff");
   grEtaMinusMSTW->SetLineColor(ci);
   grEtaMinusMSTW->SetLineStyle(9);
   grEtaMinusMSTW->SetLineWidth(4);
   grEtaMinusMSTW->Draw("lsame");

   TGraphAsymmErrors *grEtaMinusCT10 = new TGraphAsymmErrors(hEtaMinusCT10->GetNbinsX());
   grEtaMinusCT10 = convertToGraph(hEtaMinusCT10);
   ci = TColor::GetColor("#0000ff");
   grEtaMinusCT10->SetLineColor(ci);
   grEtaMinusCT10->SetLineStyle(7);
   grEtaMinusCT10->SetLineWidth(4);
   grEtaMinusCT10->Draw("lsame");

   grWmc->Draw("pesame");
   hEtaMinusDummy->Draw("sameaxis");
   leg = new TLegend(0.1694631,0.2023601,0.4832215,0.3771853,NULL,"brNDC");
   leg->SetTextFont(gStyle->GetTextFont());
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04545455);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->AddEntry(grWmc,"Data","pe");
   /*leg->AddEntry(hEtaMinusMSTWnn,"nn","pe");
   leg->AddEntry(hEtaMinusMSTWnppn,"np(pn)","pe");
   leg->AddEntry(hEtaMinusMSTWpp,"pp","pe");
   */
   leg->AddEntry(grEtaMinusMSTW,"MSTW NLO","l");
   leg->AddEntry(grEtaMinusCT10,"CT10 NLO","l");
   leg->AddEntry(grEtaMinusMcAll,"MRST LO*","l");
   leg->Draw();
   leg = new TLegend(0.4463087,0.1958042,0.7600671,0.3706294,NULL,"brNDC");
   leg->SetTextFont(gStyle->GetTextFont());
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04545455);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->AddEntry(grWmUncorrelatedSystc,"Uncorr. uncertainty","f");
   leg->AddEntry(grWmCorrelatedSystc,"Corr. uncertainty","f");
   leg->Draw();
   tex = new TLatex(0.6879195,0.7596154,"W^{-}#rightarrow#mu^{-}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.06993007);
   tex->SetLineWidth(2);
   tex->Draw();
//      tex = new TLatex(0.4211409,0.3129371,"POW+PY8 NLO");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
   tex = new TLatex(0.1862416,0.4318182,"#int Ldt #approx 0.14 nb^{-1}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03846154);
   tex->SetLineWidth(2);
   tex->Draw();
   tex = new TLatex(0.4060403,0.4383741,"#sqrt{s_{NN}} = 2.76 TeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03846154);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.6073826,0.8566434,"ATLAS Internal");
tex->SetNDC();
   tex->SetTextFont(72);
   tex->SetLineWidth(2);
   tex->Draw();
   tex = new TLatex(0.62,0.44,"0-80%");
tex->SetNDC();
tex->SetTextSize(0.03846154);
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();

   ///Draw MC/Data underneath
   canvas_dw->cd();
   /// font size
   hRatioMSTWMinusCorrelated->GetXaxis()->SetTitle("|#eta|");
   hRatioMSTWMinusCorrelated->GetXaxis()->SetLabelFont(42);
   hRatioMSTWMinusCorrelated->GetXaxis()->SetLabelSize(0.09);
   hRatioMSTWMinusCorrelated->GetXaxis()->SetTitleSize(0.11);
   hRatioMSTWMinusCorrelated->GetXaxis()->SetTitleOffset(0.8);
   hRatioMSTWMinusCorrelated->GetXaxis()->SetTitleFont(42);
   hRatioMSTWMinusCorrelated->GetYaxis()->SetTitle("Data/MSTW");
   hRatioMSTWMinusCorrelated->GetYaxis()->SetLabelFont(42);
   hRatioMSTWMinusCorrelated->GetYaxis()->SetLabelSize(0.09);
   hRatioMSTWMinusCorrelated->GetYaxis()->SetTitleSize(0.11);
   hRatioMSTWMinusCorrelated->GetYaxis()->SetTitleOffset(0.7);
   hRatioMSTWMinusCorrelated->GetYaxis()->SetRangeUser(0.65,1.45);
   hRatioMSTWMinusCorrelated->GetXaxis()->SetRangeUser(0.11,2.39);
   //hRatioMSTWMinusCorrelated->SetFillColor(kGray+1);
   hRatioMSTWMinusUncorrelated->SetFillColor(kBlue-8);
   hRatioMSTWMinusCorrelated->SetFillColor(kBlue);
   hRatioMSTWMinusCorrelated->SetFillStyle(3004);
   hRatioMSTWMinusUncorrelated->SetMarkerSize(0);
   hRatioMSTWMinusUncorrelated->Draw("e2");
   hRatioMSTWMinusCorrelated->SetMarkerSize(0);
   hRatioMSTWMinusCorrelated->Draw("e2same");
   funcRatio->Draw("same");
   hRatioMSTWMinus->SetMarkerStyle(22);
   hRatioMSTWMinus->SetMarkerColor(markerColor[1]);
   hRatioMSTWMinus->SetMarkerSize(1.5);
   TH1F* hRatioMinusDummy = (TH1F*)hRatioMSTWMinus->Clone("hRatioMinusDummy");
   hRatioMinusDummy->SetMarkerStyle(22);
   hRatioMinusDummy->SetMarkerSize(1.8);
   hRatioMinusDummy->SetMarkerColor(kWhite);
   hRatioMinusDummy->Draw("pesame");
   hRatioMSTWMinus->Draw("pesame");

   cEtaMinus->Print("etaMinus.root");
   cEtaMinus->Print("etaMinus.pdf");

   outFile->cd();
/*   hRatioCT10PlusUncorrelated->Write();
   hRatioCT10MinusUncorrelated->Write();
*/
   hRatioCT10Plus->Write();
   hRatioCT10Minus->Write();
   hWp->Write();
   hWm->Write();
   hEtaMinusCT10->Write();
   hEtaPlusCT10->Write();
   outFile->Close();
   fHisto->Close();
}
