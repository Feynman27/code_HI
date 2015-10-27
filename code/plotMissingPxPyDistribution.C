///////////////////////////
//This macro plots the missing
//px,py,and pT distributions both with 
//and without centrality binning
//as well as the missing px(y)(t) widths (i.e. std dev) 
//as a function of fcalEt 
//@date: Dec.17, 2012
//////////////////////////

#include "TDirectory.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TList.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>

void plotMissingPxPyDistribution(){

bool doMC = false;

TString fileName;
if(doMC) fileName = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/MonteCarloWmunu/MonteCarlo/*root*";
else fileName =
    "/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.04.12.2013.root";

TChain* tree = new TChain("tree","tree");
int nFiles = tree->Add(fileName);
std::cout <<"Filling the DataSet for "<< fileName << "...Number of files: " << nFiles << std::endl;


//_file0->cd();
#ifdef __CINT__
  gROOT->LoadMacro("AtlasUtils.C");
#endif


char namePx[50], namePy[50];
char namePt[50];
char namePxAll[50], namePyAll[50];

char nameEx[50], nameEy[50];
char nameEt[50];
char nameExAll[50], nameEyAll[50];
//stl container to store lower 
//track thresholds
std::vector<int> missingPThr ;
missingPThr.push_back(500);
//missingPThr.push_back(700);
missingPThr.push_back(1000);
missingPThr.push_back(2000);
missingPThr.push_back(3000); //used in analysis
missingPThr.push_back(5000);
//missingPThr.push_back(7000);
const unsigned int nMissPtBins = missingPThr.size();
const unsigned int nMissEtBins = nMissPtBins;
std::cout << " Number of lower track thresholds: " << nMissPtBins << std::endl;

std::vector<double> centralityBins;
centralityBins.push_back(0.0);
centralityBins.push_back(0.05);
centralityBins.push_back(0.1);
centralityBins.push_back(0.15);
centralityBins.push_back(0.2);
centralityBins.push_back(0.4);
centralityBins.push_back(0.6);
centralityBins.push_back(0.8);
const unsigned int nCentralityBins = centralityBins.size()-1;

double _Infinity = 1.0e30 ;
std::vector<double> fCalBins;
fCalBins.push_back(_Infinity); //limit->Infinity -->0-5
fCalBins.push_back(2.8038); //limit->Infinity -->0-5
fCalBins.push_back(); //1.9073-2.3102 -->10-15%
fCalBins.push_back(0.0236178); //1.5707-1.9073 -->15-20%
fCalBins.push_back(0.0129688); //0.6624-1.5707 -->20-40%
fCalBins.push_back(0.00521983); //0.0438-0.6624 -->40-80% 
const unsigned int nfcalBins = fCalBins.size()-1;
std::cout << " Number of centrality bins: " << nCentralityBins << std::endl;
std::cout << " Number of fcal bins: " << nfcalBins << std::endl;


TH1* hPxAll[nMissPtBins]; TH1* hPyAll[nMissPtBins];
TH1* hPx[nMissPtBins][nCentralityBins]; TH1* hPy[nMissPtBins][nCentralityBins];
TH1* hPt[nMissPtBins][nCentralityBins];

TH1* hExAll[nMissEtBins]; TH1* hEyAll[nMissEtBins];
TH1* hEx[nMissEtBins][nCentralityBins]; TH1* hEy[nMissEtBins][nCentralityBins];
TH1* hEt[nMissEtBins][nCentralityBins];

TCanvas cPt = TCanvas("cPt","cPt",600,600);
TCanvas cPx = TCanvas("cPx","cPx",600,600);
TCanvas cPy = TCanvas("cPy","cPy",600,600);
TCanvas cDummy = TCanvas("cDummy","cDummy",600,600);
TCanvas *cPxCent[nCentralityBins]; TCanvas *cPyCent[nCentralityBins];
TCanvas *cPtCent[nCentralityBins];
TLegend *leg[nCentralityBins];

//collection of Tgraphs
TList _mptMean;
TList _mpxMean;
TList _mpyMean;
TList _mpxSigma;
TList _mpySigma;

TList _metMean;
TList _mexMean;
TList _meyMean;
TList _mexSigma;
TList _meySigma;

TCanvas cMptMean = TCanvas("cMptMean","cMptMean",600,600);
TCanvas cMpxMean = TCanvas("cMpxMean","cMpxMean",600,600);
TCanvas cMpxSigma = TCanvas("cMpxSigma","cMpxSigma",600,600);
TCanvas cMpyMean = TCanvas("cMpyMean","cMpyMean",600,600);
TCanvas cMpySigma = TCanvas("cMpySigma","cMpySigma",600,600);

//TLegend* leg0 = new TLegend(0.189597, 0.725524, 0.444631, 0.90035);
TLegend* leg0 = new TLegend(0.189597, 0.725524, 0.444631, 0.94);
leg0->SetTextFont(gStyle->GetTextFont());
leg0->SetTextSize(gStyle->GetTextSize());
leg0->SetBorderSize(0);
leg0->SetFillColor(0);

TLegend* legMpt = new TLegend(0.23,0.73,0.50,0.91);
legMpt->SetTextFont(gStyle->GetTextFont());
legMpt->SetTextSize(gStyle->GetTextSize());
legMpt->SetBorderSize(0);
legMpt->SetFillColor(0);

char canvasPx[50], canvasPy[50], canvasPt[50];

int lineColor[6] = {kBlue,kRed,kGreen,kBlack,kViolet,49};
int markerStyle[6] = {24,25,26,27,28,30};

//one TGraph per MPT threshold
for(int ipt=0;ipt<nMissPtBins;ipt++){
  _mptMean.Add(new TGraphErrors(nCentralityBins));
  _mpxMean.Add(new TGraphErrors(nCentralityBins));
  _mpyMean.Add(new TGraphErrors(nCentralityBins));
  _mpxSigma.Add(new TGraph(nCentralityBins));
  _mpySigma.Add(new TGraph(nCentralityBins));

  _metMean.Add(new TGraphErrors(nCentralityBins));
  _mexMean.Add(new TGraphErrors(nCentralityBins));
  _meyMean.Add(new TGraphErrors(nCentralityBins));
  _mexSigma.Add(new TGraph(nCentralityBins));
  _meySigma.Add(new TGraph(nCentralityBins));

  for(int ical=0;ical<nfcalBins;ical++){

    sprintf(namePt,"hpt%i_cent%i",ipt,ical);
    sprintf(namePx,"hpx%i_cent%i",ipt,ical);
    sprintf(namePy,"hpy%i_cent%i",ipt,ical);

    sprintf(nameEt,"het%i_cent%i",ipt,ical);
    sprintf(nameEx,"hex%i_cent%i",ipt,ical);
    sprintf(nameEy,"hey%i_cent%i",ipt,ical);

    hPt[ipt][ical] = new TH1F(namePt,namePt,180,0.0,+180.);
    hPx[ipt][ical] = new TH1F(namePx,namePx,180,-90.0,+90.);
    hPy[ipt][ical] = new TH1F(namePy,namePy,180,-90.,+90.);

    hEt[ipt][ical] = new TH1F(nameEt,nameEt,180,0.0,+180.);
    hEx[ipt][ical] = new TH1F(nameEx,nameEx,180,-90.0,+90.);
    hEy[ipt][ical] = new TH1F(nameEy,nameEy,180,-90.,+90.);

    if(ipt==0) {
      sprintf(canvasPt, "cPt_cent%i",ical);
      cPtCent[ical] = new TCanvas(canvasPt,canvasPt,600,600);

      sprintf(canvasPx, "cPx_cent%i",ical);
      cPxCent[ical] = new TCanvas(canvasPx,canvasPx,600,600);

      sprintf(canvasPy, "cPy_cent%i",ical);
      cPyCent[ical] = new TCanvas(canvasPy,canvasPy,600,600);
    }

    double centBinLo = centralityBins.at(ical);
    double centBinUp = centralityBins.at(ical+1);
    //fcal in opposite direction to centrality
    double fcalBinUp = fCalBins.at(ical);
    double fcalBinLo = fCalBins.at(ical+1);
//    double binPt = 0.5*(centBinUp+centBinLo)*100;

    double binPt = fcalBinLo;

    int missPThr = missingPThr.at(ipt);

    //TString centBinCut = "&&Centrality>="; centBinCut+=centBinLo; centBinCut+= "&&Centrality<"; centBinCut+= centBinUp;
    TString fcalBinCut = "fcal>"; fcalBinCut+=fcalBinLo; fcalBinCut+= "&&fcal<="; fcalBinCut+= fcalBinUp;
    TString missPtThrCut = "nu_pt"; missPtThrCut+=missPThr; missPtThrCut+="Nominal";     
    TString missPxThrCut = "nu_px"; missPxThrCut+=missPThr; missPxThrCut+="Nominal";     
    TString missPyThrCut = "nu_py"; missPyThrCut+=missPThr; missPyThrCut+="Nominal"; 
    std::cout << missPtThrCut << std::endl;
    std::cout << missPxThrCut << std::endl;
    std::cout << missPyThrCut << std::endl;

    TString missEtThrCut = "caloMET"; //missEtThrCut+= ">>"; missEtThrCut+= nameEt;    
    TString missExThrCut = "caloMEx"; //missExThrCut+= ">>"; missExThrCut+= nameEx;    
    TString missEyThrCut = "caloMEy"; //missEyThrCut+= ">>"; missEyThrCut+= nameEy;

    TString ptCuts = ""; TString pxCuts = ""; TString pyCuts = "";
    TString etCuts = ""; TString exCuts = ""; TString eyCuts = "";

    if(doMC){
	
	    ptCuts += "abs("; ptCuts += missPtThrCut; ptCuts +=")<300.0"; ptCuts += fcalBinCut;
	    pxCuts += "abs("; pxCuts += missPxThrCut; pxCuts +=")<300.0"; pxCuts += fcalBinCut;
	    pyCuts += "abs("; pyCuts += missPyThrCut; pyCuts +=")<300.0"; pyCuts += fcalBinCut;

    }
    else{
   	/*ptCuts = "(EF_L1TE50_NoAlg||EF_mbZdc_a_c_L1VTE50_trk)";*/ ptCuts += fcalBinCut; ptCuts += "&&abs("; ptCuts += missPtThrCut; ptCuts +=")<300.0";
    	/*pxCuts = "(EF_L1TE50_NoAlg||EF_mbZdc_a_c_L1VTE50_trk)";*/ pxCuts += fcalBinCut; pxCuts += "&&abs("; pxCuts += missPxThrCut; pxCuts +=")<300.0";
    	/*pyCuts = "(EF_L1TE50_NoAlg||EF_mbZdc_a_c_L1VTE50_trk)";*/ pyCuts += fcalBinCut; pyCuts += "&&abs("; pyCuts += missPyThrCut; pyCuts +=")<300.0";
	
	etCuts += fcalBinCut; etCuts += "&&abs("; etCuts += missEtThrCut; etCuts +=")<300.0";
    	exCuts += fcalBinCut; exCuts += "&&abs("; exCuts += missExThrCut; exCuts +=")<300.0";
    	eyCuts += fcalBinCut; eyCuts += "&&abs("; eyCuts += missEyThrCut; eyCuts +=")<300.0";

    }

    std::cout << "cutsPx: " << pxCuts << std::endl; 
    std::cout << "cutsPy: " << pyCuts << std::endl; 
    std::cout << "cutsEx: " << exCuts << std::endl; 
    std::cout << "cutsEy: " << eyCuts << std::endl; 

    cDummy.cd();
    tree->Draw(missPtThrCut+">>"+namePt,pxCuts);
    tree->Draw(missPxThrCut+">>"+namePx,pxCuts);
    tree->Draw(missPyThrCut+">>"+namePy,pyCuts);

    tree->Draw(missEtThrCut+">>"+nameEt,etCuts);
    tree->Draw(missExThrCut+">>"+nameEx,exCuts);
    tree->Draw(missEyThrCut+">>"+nameEy,eyCuts);

    double meanPt = hPt[ipt][ical]->GetMean();
    double errPt = hPt[ipt][ical]->GetMeanError();
    double sigmaPt =  hPt[ipt][ical]->GetRMS();
    ((TGraphErrors*)_mptMean.At(ipt))->SetPoint(ical, binPt, meanPt);
    //((TGraphErrors*)_mptMean.At(ipt))->SetPointError(ical,0.0 , errPt);
    ((TGraphErrors*)_mptMean.At(ipt))->SetPointError(ical,0.0 , 0.0);

    double meanPx = hPx[ipt][ical]->GetMean();
    double errPx = hPx[ipt][ical]->GetMeanError();
    double sigmaPx =  hPx[ipt][ical]->GetRMS();
    ((TGraphErrors*)_mpxMean.At(ipt))->SetPoint(ical, binPt, meanPx);
    //((TGraphErrors*)_mpxMean.At(ipt))->SetPointError(ical,0.0 , errPx);
    ((TGraphErrors*)_mpxMean.At(ipt))->SetPointError(ical,0.0 , 0.0);
    ((TGraph*)_mpxSigma.At(ipt))->SetPoint(ical, binPt, sigmaPx);

    double meanPy = hPy[ipt][ical]->GetMean();
    double errPy = hPy[ipt][ical]->GetMeanError();
    double sigmaPy =  hPy[ipt][ical]->GetRMS();
    ((TGraphErrors*)_mpyMean.At(ipt))->SetPoint(ical, binPt, meanPy);
    //((TGraphErrors*)_mpyMean.At(ipt))->SetPointError(ical,0.0, errPy);
    ((TGraphErrors*)_mpyMean.At(ipt))->SetPointError(ical,0.0, 0.0);
    ((TGraph*)_mpySigma.At(ipt))->SetPoint(ical, binPt, sigmaPy);

    TString strPtOut = "Mean pt with "; strPtOut+=missPThr; strPtOut+="MeV track cut: "; 
    std::cout << strPtOut << "for fcal energy " << fcalBinLo << "-" << fcalBinUp << " " << meanPt << " +- " << sigmaPt << std::endl;

    TString strPxOut = "Mean px with "; strPxOut+=missPThr; strPxOut+="MeV track cut: "; 
    std::cout << strPxOut << meanPx << " +- " << sigmaPx << std::endl;


    TString strPyOut = "Mean py with "; strPyOut+=missPThr; strPyOut+="MeV track cut: "; 
    std::cout << strPyOut << meanPy << " +- " << sigmaPy << std::endl;

    double meanEt = hEt[ipt][ical]->GetMean();
    double errEt = hEt[ipt][ical]->GetMeanError();
    double sigmaEt =  hEt[ipt][ical]->GetRMS();
    ((TGraphErrors*)_metMean.At(ipt))->SetPoint(ical, binPt, meanEt);
    ((TGraphErrors*)_metMean.At(ipt))->SetPointError(ical,0.0 , errEt);

    double meanEx = hEx[ipt][ical]->GetMean();
    double errEx = hEx[ipt][ical]->GetMeanError();
    double sigmaEx =  hEx[ipt][ical]->GetRMS();
    ((TGraphErrors*)_mexMean.At(ipt))->SetPoint(ical, binPt, meanEx);
    ((TGraphErrors*)_mexMean.At(ipt))->SetPointError(ical,0.0 , errEx);
    ((TGraph*)_mexSigma.At(ipt))->SetPoint(ical, binPt, sigmaEx);

    double meanEy = hEy[ipt][ical]->GetMean();
    double errEy = hEy[ipt][ical]->GetMeanError();
    double sigmaEy =  hEy[ipt][ical]->GetRMS();
    ((TGraphErrors*)_meyMean.At(ipt))->SetPoint(ical, binPt, meanEy);
    ((TGraphErrors*)_meyMean.At(ipt))->SetPointError(ical,0.0, errEy);
    ((TGraph*)_meySigma.At(ipt))->SetPoint(ical, binPt, sigmaEy);
    
    std::cout << "Mean Et for fcal energy " << fcalBinLo << "-" << fcalBinUp << " " << meanEt << " +- " << sigmaEt << std::endl;
    std::cout << "Mean Ex for fcal energy " << fcalBinLo << "-" << fcalBinUp << " " << meanEx << " +- " << sigmaEx << std::endl;
    std::cout << "Mean Ey for fcal energy " << fcalBinLo << "-" << fcalBinUp << " " << meanEy << " +- " << sigmaEy << std::endl;


    hPt[ipt][ical]->SetLineColor(lineColor[ipt]);
    hPx[ipt][ical]->SetLineColor(lineColor[ipt]);
    hPy[ipt][ical]->SetLineColor(lineColor[ipt]);

    hPt[ipt][ical]->SetXTitle("#slash{p_{T}}[GeV]");
    hPx[ipt][ical]->SetXTitle("#slash{p_{x}}[GeV]");
    hPy[ipt][ical]->SetXTitle("#slash{p_{y}}[GeV]");


    TString thrLabel = ""; thrLabel+=missPThr; thrLabel+="MeV";
    TString centLabel = ""; centLabel+=centBinLo*100; centLabel+="-"; centLabel+=centBinUp*100; centLabel+="%"; 

    if(ipt==nMissPtBins-1&&ical==0) std::cout<<centLabel<<std::endl;

    cPtCent[ical]->cd();
    if(ipt==0) {
       //TH1F* hdummy = new TH1F("hdummy","hdummy",180,-90.0,+90.);
       TH1F* hdummy = new TH1F("hdummy","hdummy",80,-40.0,+40.);
       hdummy->SetXTitle("#slash{E_{T}},#slash{p_{T}}[GeV]");
       hdummy->GetYaxis()->SetRangeUser(0.,0.2);
       hdummy->Draw();
       myText(0.729866,0.770979,1,centLabel); 
       //leg[ical] = new TLegend(0.189597, 0.725524, 0.444631, 0.90035);
       leg[ical] = new TLegend(0.189597, 0.725524, 0.444631, 0.935);
       leg[ical]->SetTextFont(gStyle->GetTextFont());
       leg[ical]->SetTextSize(gStyle->GetTextSize());
       leg[ical]->SetBorderSize(0);
       leg[ical]->SetFillColor(0);
     }

    hPt[ipt][ical]->DrawNormalized("hist same") ;
    hEt[ipt][ical]->DrawNormalized("hist same") ;

    leg[ical]->AddEntry(hPt[ipt][ical], thrLabel, "l");
    if(ipt==0)leg[ical]->AddEntry(hEt[ipt][ical], "#slash{E_{T}}", "l");
    leg[ical]->Draw();
    
    cPtCent[ical]->Update();

    cPxCent[ical]->cd();

    if(ipt==0) {
       //TH1F* hdummy = new TH1F("hdummy","hdummy",180,-90.0,+90.);
       TH1F* hdummy = new TH1F("hdummy","hdummy",80,-40.0,+40.);
       hdummy->SetXTitle("#slash{E_{x}},#slash{p_{x}}[GeV]");
       hdummy->GetYaxis()->SetRangeUser(0.,0.2);
       hdummy->Draw();
       myText(0.729866,0.770979,1,centLabel); 
     }

    hPx[ipt][ical]->DrawNormalized("hist same") ;
    hEx[ipt][ical]->DrawNormalized("hist same") ;

    leg[ical]->Draw();
    cPxCent[ical]->Update();

    cPyCent[ical]->cd();
    if(ipt==0){ 
       //TH1F* hdummy = new TH1F("hdummy","hdummy",180,-90.0,+90.);
       TH1F* hdummy = new TH1F("hdummy","hdummy",80,-40.0,+40.);
       hdummy->SetXTitle("#slash{E_{y}},#slash{p_{y}}[GeV]");
       hdummy->GetYaxis()->SetRangeUser(0.,0.2);
       hdummy->Draw();
       myText(0.729866,0.770979,1,centLabel); 
    }

    hPy[ipt][ical]->DrawNormalized("hist same") ;
    hEy[ipt][ical]->DrawNormalized("hist same") ;
    leg[ical]->Draw();
    cPyCent[ical]->Update();

    if(ipt==nMissPtBins-1) {
    	cPtCent[ical]->Print("missPtDistr"+centLabel+".pdf");
    	cPxCent[ical]->Print("missPxDistr"+centLabel+".pdf");
    	cPyCent[ical]->Print("missPyDistr"+centLabel+".pdf");
    	cPtCent[ical]->Print("missPtDistr"+centLabel+".root");
    	cPxCent[ical]->Print("missPxDistr"+centLabel+".root");
    	cPyCent[ical]->Print("missPyDistr"+centLabel+".root");
    }
  } //ical
  
  cMptMean.cd();

  (((TGraphErrors*)_mptMean.At(ipt)))->SetMarkerStyle(markerStyle[ipt]);
  (((TGraphErrors*)_mptMean.At(ipt)))->SetMarkerColor(lineColor[ipt]);

  //(((TGraphErrors*)_metMean.At(ipt)))->SetMarkerStyle(markerStyle[ipt]);
  //(((TGraphErrors*)_metMean.At(ipt)))->SetMarkerColor(lineColor[ipt]);

  TString sTitle = ""; sTitle+=missingPThr.at(ipt); sTitle+="MeV";
  legMpt->AddEntry( ((TGraphErrors*)_mptMean.At(ipt)),sTitle,"p");
  legMpt->AddEntry( ((TGraphErrors*)_metMean.At(0)),"#slash{E_{T}}","p");

  if(ipt==0){
    ((TGraphErrors*)_mptMean.At(ipt))->GetXaxis()->SetTitle("#Sigma E_{T}^{fcal}[TeV]");	
    ((TGraphErrors*)_mptMean.At(ipt))->GetYaxis()->SetTitle("#LT #slash{E_{T}} #GT,#LT #slash{p_{T}} #GT");	
    ((TGraphErrors*)_mptMean.At(ipt))->GetXaxis()->SetRangeUser(0.0,3.0);	
    ((TGraphErrors*)_mptMean.At(ipt))->GetYaxis()->SetRangeUser(0.0,70.0);	
    ((TGraphErrors*)_mptMean.At(ipt))->Draw("ap");
    ((TGraphErrors*)_metMean.At(ipt))->Draw("p same");
  }
  else{ ((TGraphErrors*)_mptMean.At(ipt))->Draw("p same");
        //((TGraphErrors*)_metMean.At(ipt))->Draw("p same");
      }

  legMpt->Draw(); cMptMean.Update();

  cMpxMean.cd();

  (((TGraphErrors*)_mpxMean.At(ipt)))->SetMarkerStyle(markerStyle[ipt]);
  (((TGraphErrors*)_mpxMean.At(ipt)))->SetMarkerColor(lineColor[ipt]);

//  (((TGraphErrors*)_mexMean.At(ipt)))->SetMarkerStyle(markerStyle[ipt]);
//  (((TGraphErrors*)_mexMean.At(ipt)))->SetMarkerColor(lineColor[ipt]);

  TString sTitle = ""; sTitle+=missingPThr.at(ipt); sTitle+="MeV";
  if(ipt==0){
    ((TGraphErrors*)_mpxMean.At(ipt))->GetXaxis()->SetTitle("#Sigma E_{T}^{fcal}[TeV]");	
    ((TGraphErrors*)_mpxMean.At(ipt))->GetYaxis()->SetTitle("#LT #slash{p_{x}} #GT");	
    ((TGraphErrors*)_mpxMean.At(ipt))->GetXaxis()->SetRangeUser(0.0,3.0);	
    //((TGraphErrors*)_mpxMean.At(ipt))->GetYaxis()->SetRangeUser(-50.0,70.0);	
    ((TGraphErrors*)_mpxMean.At(ipt))->GetYaxis()->SetRangeUser(-20.0,40.0);	
    ((TGraphErrors*)_mpxMean.At(ipt))->Draw("ap");
    ((TGraphErrors*)_mexMean.At(ipt))->Draw("p same");
  }
  else { ((TGraphErrors*)_mpxMean.At(ipt))->Draw("p same");
         //((TGraphErrors*)_mexMean.At(ipt))->Draw("p same");
  }
  legMpt->Draw(); cMpxMean.Update();

  cMpxSigma.cd();
  
  ((TGraph*)_mpxSigma.At(ipt))->SetMarkerStyle(markerStyle[ipt]);	
  ((TGraph*)_mpxSigma.At(ipt))->SetMarkerColor(lineColor[ipt]);	

//  ((TGraph*)_mexSigma.At(ipt))->SetMarkerStyle(markerStyle[ipt]);	
//  ((TGraph*)_mexSigma.At(ipt))->SetMarkerColor(lineColor[ipt]);	

  if (ipt==0){
    ((TGraph*)_mpxSigma.At(ipt))->GetXaxis()->SetTitle("#Sigma E_{T}^{fcal}[TeV]");	
    ((TGraph*)_mpxSigma.At(ipt))->GetYaxis()->SetTitle("#sigma_{x}");	
    ((TGraph*)_mpxSigma.At(ipt))->GetXaxis()->SetRangeUser(0.0,3.0);	
    ((TGraph*)_mpxSigma.At(ipt))->GetYaxis()->SetRangeUser(0.0,65.0);	
    ((TGraph*)_mpxSigma.At(ipt))->Draw("ap");
    ((TGraph*)_mexSigma.At(ipt))->Draw("p same");
  }
  else { ((TGraph*)_mpxSigma.At(ipt))->Draw("p same");
  	 //((TGraph*)_mexSigma.At(ipt))->Draw("p same");
       }

  legMpt->Draw(); cMpxSigma.Update();

  cMpyMean.cd();

  ((TGraphErrors*)_mpyMean.At(ipt))->SetMarkerStyle(markerStyle[ipt]);	
  ((TGraphErrors*)_mpyMean.At(ipt))->SetMarkerColor(lineColor[ipt]);	

  ((TGraphErrors*)_meyMean.At(ipt))->SetMarkerStyle(markerStyle[ipt]);	
  ((TGraphErrors*)_meyMean.At(ipt))->SetMarkerColor(lineColor[ipt]);	

  if (ipt==0){
    ((TGraphErrors*)_mpyMean.At(ipt))->GetXaxis()->SetTitle("#Sigma E_{T}^{fcal}[TeV]");	
    ((TGraphErrors*)_mpyMean.At(ipt))->GetYaxis()->SetTitle("#LT #slash{p_{y}} #GT");	
    ((TGraphErrors*)_mpyMean.At(ipt))->GetXaxis()->SetRangeUser(0.0,3.0);	
    //((TGraphErrors*)_mpyMean.At(ipt))->GetYaxis()->SetRangeUser(-50.0,70.0);	
    ((TGraphErrors*)_mpyMean.At(ipt))->GetYaxis()->SetRangeUser(-20.0,40.0);	
    ((TGraphErrors*)_mpyMean.At(ipt))->Draw("ap");
    ((TGraphErrors*)_meyMean.At(ipt))->Draw("p same");
  }
  else { ((TGraphErrors*)_mpyMean.At(ipt))->Draw("p same");
         //((TGraphErrors*)_meyMean.At(ipt))->Draw("p same");
  }
  legMpt->Draw(); cMpyMean.Update();

  cMpySigma.cd();
  
  
  ((TGraph*)_mpySigma.At(ipt))->SetMarkerStyle(markerStyle[ipt]);	
  ((TGraph*)_mpySigma.At(ipt))->SetMarkerColor(lineColor[ipt]);	

  ((TGraph*)_meySigma.At(ipt))->SetMarkerStyle(markerStyle[ipt]);	
  ((TGraph*)_meySigma.At(ipt))->SetMarkerColor(lineColor[ipt]);	

  if (ipt==0){
    ((TGraph*)_mpySigma.At(ipt))->GetXaxis()->SetTitle("#Sigma E_{T}^{fcal}[TeV]");	
    ((TGraph*)_mpySigma.At(ipt))->GetYaxis()->SetTitle("#sigma_{y}");	
    ((TGraph*)_mpySigma.At(ipt))->GetXaxis()->SetRangeUser(0.0,3.0);	
    ((TGraph*)_mpySigma.At(ipt))->GetYaxis()->SetRangeUser(0.0,65.0);	
    ((TGraph*)_mpySigma.At(ipt))->Draw("ap");
    ((TGraph*)_meySigma.At(ipt))->Draw("p same");
  }
  else{ ((TGraph*)_mpySigma.At(ipt))->Draw("p same");
        //((TGraph*)_meySigma.At(ipt))->Draw("p same");
  }
  legMpt->Draw(); cMpySigma.Update();


  if(ipt==nMissPtBins-1) {

	cMptMean.Print("MPtMean.pdf");
	cMptMean.Print("MPtMean.root");
	cMpxMean.Print("MPxMean.pdf");
	cMpxMean.Print("MPxMean.root");
	cMpyMean.Print("MPyMean.pdf");
	cMpyMean.Print("MPyMean.root");
	cMpxSigma.Print("MPxSigma.pdf");
	cMpxSigma.Print("MPxSigma.root");
	cMpySigma.Print("MPySigma.pdf");
	cMpySigma.Print("MPySigma.root");
  }
 } //ipt


std::cout << "Done." << std::endl;

}
