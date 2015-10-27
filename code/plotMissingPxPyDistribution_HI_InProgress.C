///////////////////////////
//This macro plots the missing
//px,py,and pT distributions both with 
//and without centrality binning
//as well as the missing px(y)(t) resolution
//as a function of fcalEt 
//@date: Dec.17, 2012
//////////////////////////

#include "TCanvas.h"
#include "TDirectory.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2DErrors.h"
#include "TList.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TChain.h"
#include "THStack.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLatex.h"
#include "TClonesArray.h"

#include "RooDataSet.h"
#include "RooBinning.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHist.h"

#include "WAnalysisHIDep.C"
#include "AtlasUtils.C"
#include "AtlasStyle.C"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath> 
#include <sstream>


RooDataSet* selectFcal(RooDataSet* dataSet, double fcalUp, double fcalLow){
    TString cuts = "fcal>";
    cuts +=fcalLow;
    cuts += "&&fcal<";
    cuts += fcalUp;

    return (RooDataSet*)dataSet->reduce(cuts);
}

RooDataSet* fillDataSet(const TString& pathName, const TString& fileName, RooArgSet& argSet, int mptThresh)
{
  RooDataSet* set = new RooDataSet(fileName,fileName,argSet);
  
  float fcalNt;
  float mpt,mpx,mpy;
  int trig1,trig2;

  // --- Load the MC tree ---
  TChain* tree = new TChain("data","data");
  int nFiles = tree->Add(pathName+fileName);
  
  // --- Set branch adresses ---
  tree->SetBranchAddress("Fcal_Et", &fcalNt);
  tree->SetBranchAddress("EF_L1TE50_NoAlg",&trig1);
  tree->SetBranchAddress("EF_mbZdc_a_c_L1VTE50_trk",&trig2);
  if(mptThresh==500){
    tree->SetBranchAddress("nu_pt500", &mpt);
    tree->SetBranchAddress("nu_px500", &mpx);
    tree->SetBranchAddress("nu_py500", &mpy);
  } else if(mptThresh==1000){
     tree->SetBranchAddress("nu_pt1000", &mpt);
     tree->SetBranchAddress("nu_px1000", &mpx);
     tree->SetBranchAddress("nu_py1000", &mpy);
  }else if(mptThresh==2000){
     tree->SetBranchAddress("nu_pt2000", &mpt);
     tree->SetBranchAddress("nu_px2000", &mpx);
     tree->SetBranchAddress("nu_py2000", &mpy);
  }else if(mptThresh==3000){
     tree->SetBranchAddress("nu_pt3000", &mpt);
     tree->SetBranchAddress("nu_px3000", &mpx);
     tree->SetBranchAddress("nu_py3000", &mpy);
  }else if(mptThresh==5000){
     tree->SetBranchAddress("nu_pt5000", &mpt);
     tree->SetBranchAddress("nu_px5000", &mpx);
     tree->SetBranchAddress("nu_py5000", &mpy);
  }

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("EF_L1TE50_NoAlg",1);
  tree->SetBranchStatus("EF_mbZdc_a_c_L1VTE50_trk",1);
  tree->SetBranchStatus("Fcal_Et", 1);
  tree->SetBranchStatus("nu_pt500", 1);
  tree->SetBranchStatus("nu_px500", 1);
  tree->SetBranchStatus("nu_py500", 1);
  tree->SetBranchStatus("nu_pt1000", 1);
  tree->SetBranchStatus("nu_px1000", 1);
  tree->SetBranchStatus("nu_py1000", 1);
  tree->SetBranchStatus("nu_pt2000", 1);
  tree->SetBranchStatus("nu_px2000", 1);
  tree->SetBranchStatus("nu_py2000", 1);
  tree->SetBranchStatus("nu_pt3000", 1);
  tree->SetBranchStatus("nu_px3000", 1);
  tree->SetBranchStatus("nu_py3000", 1);
  tree->SetBranchStatus("nu_pt5000", 1);
  tree->SetBranchStatus("nu_px5000", 1);
  tree->SetBranchStatus("nu_py5000", 1);
  
  std::cout << "Number of entries: " << tree->GetEntries() << std::endl; 
  for ( int i = 0; i < tree->GetEntries(); i++ ) {
     if (i>100000) return set; // temp hack
    tree->LoadTree(i);
    tree->GetEntry(i);
    if( (trig1)||(trig2) ){
     argSet.setRealValue("fcal",fcalNt);
     argSet.setRealValue("MPt",mpt);
     argSet.setRealValue("MPx",mpx);
     argSet.setRealValue("MPy",mpy);
     set->add(argSet);
    }
  }

  return set;
}


void plotMissingPxPyDistribution_HI_InProgress(){

bool doMC = false;
TString sDate = ".06.10.2013";

TString fileName;
if(doMC) fileName = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/MonteCarloWmunu/MonteCarlo/*root*";
else fileName = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/MinimumBiasFiles/MinimumBias/*root*";
TString pathName = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/MinimumBiasFiles/MinimumBias/";

char namePx[50], namePy[50];
char namePt[50];
char namePxAll[50], namePyAll[50];

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
std::cout << " Number of lower track thresholds: " << nMissPtBins << std::endl;

std::vector<double> centralityBins;
centralityBins.push_back(0.0);
centralityBins.push_back(0.05);
centralityBins.push_back(0.1);
centralityBins.push_back(0.15);
centralityBins.push_back(0.2);
centralityBins.push_back(0.4);
centralityBins.push_back(0.8);
const unsigned int nCentralityBins = centralityBins.size()-1;

double _Infinity = 1.0e30 ;
std::vector<double> fCalBins;
fCalBins.push_back(_Infinity); //limit->Infinity -->0-5
fCalBins.push_back(2.8038); //>2.8038       -->0-5
fCalBins.push_back(2.3102); //2.3102-2.8038 -->5-10
fCalBins.push_back(1.9073); //1.9073-2.3102 -->10-15%
fCalBins.push_back(1.5707); //1.5707-1.9073 -->15-20%
fCalBins.push_back(0.6624); //0.6624-1.5707 -->20-40%
fCalBins.push_back(0.0438); //0.0438-0.6624 -->40-80% 
const unsigned int nFcalBins = fCalBins.size()-1;
std::cout << " Number of centrality bins: " << nCentralityBins << std::endl;
std::cout << " Number of fcal bins: " << nFcalBins << std::endl;


TH1F* hPxAll[nMissPtBins]; TH1F* hPyAll[nMissPtBins];
TH1F* hPx[nMissPtBins][nCentralityBins]; TH1F* hPy[nMissPtBins][nCentralityBins];
TH1F* hPt[nMissPtBins][nCentralityBins];

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

TCanvas cMptMean = TCanvas("cMptMean","cMptMean",600,600);
TCanvas cMpxMean = TCanvas("cMpxMean","cMpxMean",600,600);
TCanvas cMpxSigma = TCanvas("cMpxSigma","cMpxSigma",600,600);
TCanvas cMpyMean = TCanvas("cMpyMean","cMpyMean",600,600);
TCanvas cMpySigma = TCanvas("cMpySigma","cMpySigma",600,600);

TLegend* leg0 = new TLegend(0.189597, 0.725524, 0.444631, 0.935);
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

int lineColor[6] = {kBlue,kRed,kGreen,kCyan,kViolet,kTeal-9};
int markerStyle[6] = {24,25,26,27,28,30};

RooRealVar MPt("MPt","#slash{p_{T}}",0.0,350.0,"GeV");
RooRealVar MPx("MPx","#slash{p_{x}}",0.0,350.0,"GeV");
RooRealVar MPy("MPy","#slash{p_{y}}",0.0,350.0,"GeV");
RooRealVar fcal("fcal","FCal",0.0,4.0,"TeV");

TString ptCuts = ""; TString pxCuts = ""; TString pyCuts = "";

ptCuts = "fabs(MPt)<300.0";
pxCuts = "fabs(MPx)<300.0";
pyCuts = "fabs(MPy)<300.0";

std::cout << "cutsPt: " << ptCuts << std::endl; 
std::cout << "cutsPx: " << pxCuts << std::endl; 
std::cout << "cutsPy: " << pyCuts << std::endl; 

RooArgList ptArgList(MPt,fcal);
RooArgList pxArgList(MPx,fcal);
RooArgList pyArgList(MPy,fcal);

RooArgSet argSet(MPt,MPx,MPy,fcal);

RooFormulaVar cutsPt("cutsPt","cutsPt",ptCuts,MPt);
RooFormulaVar cutsPx("cutsPx","cutsPx",pxCuts,MPx);
RooFormulaVar cutsPy("cutsPy","cutsPy",pyCuts,MPy);

RooDataSet* ds[nMissPtBins];
RooDataSet* dsPt[nMissPtBins][nFcalBins];
RooDataSet* dsPx[nMissPtBins][nFcalBins];
RooDataSet* dsPy[nMissPtBins][nFcalBins];
//one TGraph per MPT threshold
for(int ipt=0;ipt<nMissPtBins;ipt++){
  _mptMean.Add(new TGraphErrors(nCentralityBins));
  _mpxMean.Add(new TGraphErrors(nCentralityBins));
  _mpyMean.Add(new TGraphErrors(nCentralityBins));
  _mpxSigma.Add(new TGraph(nCentralityBins));
  _mpySigma.Add(new TGraph(nCentralityBins));

  int missPThr = missingPThr.at(ipt);
  std::cout << "MPt threshold: " << missPThr << std::endl;
  ds[ipt] = fillDataSet(pathName,"*root*",argSet,missPThr); ds[ipt]->Print();
  for(int ical=0;ical<nFcalBins;ical++){

    std::cout << ipt << ":" << ical << std::endl;
    sprintf(namePt,"hpt%i_cent%i",ipt,ical);
    sprintf(namePx,"hpx%i_cent%i",ipt,ical);
    sprintf(namePy,"hpy%i_cent%i",ipt,ical);
    //hPt[ipt][ical] = new TH1F(namePt,namePt,180,0.0,+180.);
    //hPx[ipt][ical] = new TH1F(namePx,namePx,180,-90.0,+90.);
    //hPy[ipt][ical] = new TH1F(namePy,namePy,180,-90.,+90.);
    RooBinning b = RooBinning(180,-90.,+90.);

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


    RooDataSet *dsSubSet = selectFcal(ds[ipt],fCalBins[ical],fCalBins[ical+1]);dsSubSet->Print();
    dsPt[ipt][ical] = (RooDataSet*)dsSubSet->reduce(cutsPt); dsPt[ipt][ical]->Print();
    dsPx[ipt][ical] = (RooDataSet*)dsSubSet->reduce(cutsPx); dsPx[ipt][ical]->Print();
    dsPy[ipt][ical] = (RooDataSet*)dsSubSet->reduce(cutsPy); dsPy[ipt][ical]->Print();

    /*dsPx[ipt][ical] = fillDataSet(pathName,"*root*",argSet,missPThr); dsPx[ipt][ical]->Print();
    dsPx[ipt][ical] = selectFcal(dsPx[ipt][ical],fCalBins[ical],fCalBins[ical+1]);dsPx[ipt][ical]->Print();
    dsPx[ipt][ical] = (RooDataSet*)dsPx[ipt][ical]->reduce(cutsPx); dsPx[ipt][ical]->Print();

    dsPy[ipt][ical] = fillDataSet(pathName,"*root*",argSet,missPThr); dsPy[ipt][ical]->Print();
    dsPy[ipt][ical] = selectFcal(dsPy[ipt][ical],fCalBins[ical],fCalBins[ical+1]);dsPy[ipt][ical]->Print();
    dsPy[ipt][ical] = (RooDataSet*)dsPy[ipt][ical]->reduce(cutsPy); dsPy[ipt][ical]->Print();
*/
    ///Fill mpt(x)(y) histos for this pT track threshold and centrality bin
    hPt[ipt][ical] = (TH1F*)dsPt[ipt][ical]->createHistogram(namePt,MPt,Binning(b));
    hPx[ipt][ical] = (TH1F*)dsPx[ipt][ical]->createHistogram(namePx,MPx,Binning(b));
    hPy[ipt][ical] = (TH1F*)dsPy[ipt][ical]->createHistogram(namePy,MPy,Binning(b));

    double meanPt = hPt[ipt][ical]->GetMean();
    double sigmaPt =  hPt[ipt][ical]->GetRMS();
    ((TGraphErrors*)_mptMean.At(ipt))->SetPoint(ical, binPt, meanPt);
    ((TGraphErrors*)_mptMean.At(ipt))->SetPointError(ical,0.0 , sigmaPt);

    double meanPx = hPx[ipt][ical]->GetMean();
    double sigmaPx =  hPx[ipt][ical]->GetRMS();
    ((TGraphErrors*)_mpxMean.At(ipt))->SetPoint(ical, binPt, meanPx);
    ((TGraphErrors*)_mpxMean.At(ipt))->SetPointError(ical,0.0 , 0.0);
    ((TGraph*)_mpxSigma.At(ipt))->SetPoint(ical, binPt, sigmaPx);

    TString strPtOut = "Mean pt with "; strPtOut+=missPThr; strPtOut+="MeV track cut: "; 
    std::cout << strPtOut << "for fcal energy " << fcalBinLo << "-" << fcalBinUp << " " << meanPt << " +- " << sigmaPt << std::endl;

    TString strPxOut = "Mean px with "; strPxOut+=missPThr; strPxOut+="MeV track cut: "; 
    std::cout << strPxOut << meanPx << " +- " << sigmaPx << std::endl;

    double meanPy = hPy[ipt][ical]->GetMean();
    double sigmaPy =  hPy[ipt][ical]->GetRMS();
    ((TGraphErrors*)_mpyMean.At(ipt))->SetPoint(ical, binPt, meanPy);
    ((TGraphErrors*)_mpyMean.At(ipt))->SetPointError(ical,0.0, 0.0);
    ((TGraph*)_mpySigma.At(ipt))->SetPoint(ical, binPt, sigmaPy);

    TString strPyOut = "Mean py with "; strPyOut+=missPThr; strPyOut+="MeV track cut: "; 
    std::cout << strPyOut << meanPy << " +- " << sigmaPy << std::endl;


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
       TH1F* hdummy = new TH1F("hdummy","hdummy",180,-90.0,+90.);
       hdummy->SetXTitle("#slash{p_{T}}[GeV]");
       hdummy->GetYaxis()->SetRangeUser(0.,0.095);
       hdummy->Draw();
       myText(0.729866,0.770979,1,centLabel); 
       leg[ical] = new TLegend(0.189597, 0.725524, 0.444631, 0.935);
       leg[ical]->SetTextFont(gStyle->GetTextFont());
       leg[ical]->SetTextSize(gStyle->GetTextSize());
       leg[ical]->SetBorderSize(0);
       leg[ical]->SetFillColor(0);
     }

    hPt[ipt][ical]->DrawNormalized("same") ;

    leg[ical]->AddEntry(hPt[ipt][ical], thrLabel, "l");
    leg[ical]->Draw();
    
    cPtCent[ical]->Update();

    cPxCent[ical]->cd();
     if(ipt==0) {
       TH1F* hdummy = new TH1F("hdummy","hdummy",180,-90.0,+90.);
       hdummy->SetXTitle("#slash{p_{x}}[GeV]");
       hdummy->GetYaxis()->SetRangeUser(0.,0.095);
       hdummy->Draw();
       myText(0.729866,0.770979,1,centLabel); 
         }

    hPx[ipt][ical]->DrawNormalized("same") ;
    //leg[ical]->AddEntry(hPx[ipt][ical], thrLabel, "l");
    leg[ical]->Draw();
    cPxCent[ical]->Update();

    cPyCent[ical]->cd();
    if(ipt==0){ 
       TH1F* hdummy = new TH1F("hdummy","hdummy",180,-90.0,+90.);
       hdummy->SetXTitle("#slash{p_{y}}[GeV]");
       hdummy->GetYaxis()->SetRangeUser(0.,0.095);
       hdummy->Draw();
       myText(0.729866,0.770979,1,centLabel); 
    }

    hPy[ipt][ical]->DrawNormalized("same") ;
    leg[ical]->Draw();
    cPyCent[ical]->Update();

    if(ipt==nMissPtBins-1) {
    	/*cPtCent[ical]->Print("missPtDistr"+ical+".pdf");
    	cPxCent[ical]->Print("missPxDistr"+ical+".pdf");
    	cPyCent[ical]->Print("missPyDistr"+ical+".pdf");

    	cPtCent[ical]->Print("missPtDistrCent"+ical+".pdf");
    	cPxCent[ical]->Print("missPxDistrCent"+ical+".pdf");
    	cPyCent[ical]->Print("missPyDistrCent"+ical+".pdf");
        */
    }
  } //ical
  
  cMptMean.cd();

  (((TGraphErrors*)_mptMean.At(ipt)))->SetMarkerStyle(markerStyle[ipt]);
  (((TGraphErrors*)_mptMean.At(ipt)))->SetMarkerColor(lineColor[ipt]);
  if(ipt==3) {
    ((TGraph*)_mptMean.At(ipt))->SetMarkerSize(1.4);
  }
  TString sTitle = ""; sTitle+=missingPThr.at(ipt); sTitle+="MeV";
  legMpt->AddEntry( ((TGraphErrors*)_mptMean.At(ipt)),sTitle,"p");
  if(ipt==0){
    ((TGraphErrors*)_mptMean.At(ipt))->GetXaxis()->SetTitle("#Sigma E_{T}^{fcal}[TeV]");	
    ((TGraphErrors*)_mptMean.At(ipt))->GetYaxis()->SetTitle("#LT #slash{p_{T}} #GT");	
    ((TGraphErrors*)_mptMean.At(ipt))->GetXaxis()->SetRangeUser(0.0,3.0);	
    ((TGraphErrors*)_mptMean.At(ipt))->GetYaxis()->SetRangeUser(0.0,70.0);	
    ((TGraphErrors*)_mptMean.At(ipt))->Draw("ap");
  }
  else ((TGraphErrors*)_mptMean.At(ipt))->Draw("p same");
  legMpt->Draw(); cMptMean.Update();

  cMpxMean.cd();

  (((TGraphErrors*)_mpxMean.At(ipt)))->SetMarkerStyle(markerStyle[ipt]);
  (((TGraphErrors*)_mpxMean.At(ipt)))->SetMarkerColor(lineColor[ipt]);
  if(ipt==3) {
    ((TGraph*)_mpxMean.At(ipt))->SetMarkerSize(1.4);
  }
  TString sTitle = ""; sTitle+=missingPThr.at(ipt); sTitle+="MeV";
  if(ipt==0){
    ((TGraphErrors*)_mpxMean.At(ipt))->GetXaxis()->SetTitle("#Sigma E_{T}^{fcal}[TeV]");	
    ((TGraphErrors*)_mpxMean.At(ipt))->GetYaxis()->SetTitle("#LT #slash{p_{x}} #GT [GeV]");	
    ((TGraphErrors*)_mpxMean.At(ipt))->GetXaxis()->SetRangeUser(0.0,3.0);	
    ((TGraphErrors*)_mpxMean.At(ipt))->GetYaxis()->SetRangeUser(-1.0,11.0);	
    ((TGraphErrors*)_mpxMean.At(ipt))->Draw("ap");
  }
  else ((TGraphErrors*)_mpxMean.At(ipt))->Draw("p same");
  legMpt->Draw(); cMpxMean.Update();

  cMpxSigma.cd();
  
  ((TGraph*)_mpxSigma.At(ipt))->SetMarkerStyle(markerStyle[ipt]);	
  ((TGraph*)_mpxSigma.At(ipt))->SetMarkerColor(lineColor[ipt]);	
  if(ipt==3) {
    ((TGraph*)_mpxSigma.At(ipt))->SetMarkerSize(1.4);
  }
  if (ipt==0){
    ((TGraph*)_mpxSigma.At(ipt))->GetXaxis()->SetTitle("#Sigma E_{T}^{fcal}[TeV]");	
    ((TGraph*)_mpxSigma.At(ipt))->GetYaxis()->SetTitle("#sigma_{x}[GeV]");	
    ((TGraph*)_mpxSigma.At(ipt))->GetXaxis()->SetRangeUser(0.0,3.0);	
    ((TGraph*)_mpxSigma.At(ipt))->GetYaxis()->SetRangeUser(0.0,45.0);	
    ((TGraph*)_mpxSigma.At(ipt))->Draw("ap");
  }
  else ((TGraph*)_mpxSigma.At(ipt))->Draw("p same");
  legMpt->Draw(); cMpxSigma.Update();

  cMpyMean.cd();

  ((TGraphErrors*)_mpyMean.At(ipt))->SetMarkerStyle(markerStyle[ipt]);	
  ((TGraphErrors*)_mpyMean.At(ipt))->SetMarkerColor(lineColor[ipt]);	
  if(ipt==3) {
    ((TGraph*)_mpuMean.At(ipt))->SetMarkerSize(1.4);
  }
  if (ipt==0){
    ((TGraphErrors*)_mpyMean.At(ipt))->GetXaxis()->SetTitle("#Sigma E_{T}^{fcal}[TeV]");	
    ((TGraphErrors*)_mpyMean.At(ipt))->GetYaxis()->SetTitle("#LT #slash{p_{y}} #GT [GeV]");	
    ((TGraphErrors*)_mpyMean.At(ipt))->GetXaxis()->SetRangeUser(0.0,3.0);	
    ((TGraphErrors*)_mpyMean.At(ipt))->GetYaxis()->SetRangeUser(-1.0,11.0);	
    ((TGraphErrors*)_mpyMean.At(ipt))->Draw("ap");
  }
  else ((TGraphErrors*)_mpyMean.At(ipt))->Draw("p same");
  legMpt->Draw(); cMpyMean.Update();

  cMpySigma.cd();
  
  
  ((TGraph*)_mpySigma.At(ipt))->SetMarkerStyle(markerStyle[ipt]);	
  ((TGraph*)_mpySigma.At(ipt))->SetMarkerColor(lineColor[ipt]);	
  if(ipt==3) {
    ((TGraph*)_mpySigma.At(ipt))->SetMarkerSize(1.4);
  }
  if (ipt==0){
    ((TGraph*)_mpySigma.At(ipt))->GetXaxis()->SetTitle("#Sigma E_{T}^{fcal}[TeV]");	
    ((TGraph*)_mpySigma.At(ipt))->GetYaxis()->SetTitle("#sigma_{y}[GeV]");	
    ((TGraph*)_mpySigma.At(ipt))->GetXaxis()->SetRangeUser(0.0,3.0);	
    ((TGraph*)_mpySigma.At(ipt))->GetYaxis()->SetRangeUser(0.0,45.0);	
    ((TGraph*)_mpySigma.At(ipt))->Draw("ap");
  }
  else ((TGraph*)_mpySigma.At(ipt))->Draw("p same");
  legMpt->Draw(); cMpySigma.Update();


  if(ipt==nMissPtBins-1) {

/*	cMptMean.Print("MPtMean.pdf");
	cMptMean.SaveAs("MPtMean.root");
	cMpxMean.Print("MPxMean.pdf");
	cMpxMean.SaveAs("MPxMean.root");
	cMpyMean.Print("MPyMean.pdf");
	cMpyMean.SaveAs("MPyMean.root");
	cMpxSigma.Print("MPxSigma.pdf");
	cMpxSigma.SaveAs("MPxSigma.root");
	cMpySigma.Print("MPySigma.pdf");
	cMpySigma.SaveAs("MPySigma.root");
    */
  }
 } //ipt


std::cout << "Done." << std::endl;

}
