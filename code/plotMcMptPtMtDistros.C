#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
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
#include "TRandom.h"

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

#include "WPlotterHelper.C"

float smearGeneratorPt(float pt, float smear){

  float mcPtSmeared = pt;
  ///smear the generated pt
  mcPtSmeared +=gRandom->Gaus()*smear;
  return mcPtSmeared;
}

double getEfficiency(double pt,TGraphAsymmErrors* grEff){
    ///array of xBin Centers
    double* xPt = grEff->GetX();
    double* xLo = grEff->GetEXlow();
    double* xUp = grEff->GetEXhigh();
    ///array of efficiencies
    double* yEff = grEff->GetY();
    ///Multiply by efficiency corresponding to bin in hTau
    for(int iEff = 0; iEff<grEff->GetN(); ++iEff){

        ///get bin center
        double binCenter = xPt[iEff];
        ///lower end of bin
        double binLo = xPt[iEff] - xLo[iEff];
        ///upper end of bin
        double binUpp = xPt[iEff] + xUp[iEff];

         if(pt>binLo&&pt<binUpp){

            return yEff[iEff];
         }
    }

    std::cout << "WARNING: No efficiency found for pT " << pt << std::endl;
    return -9999.0;

}
///////////////////////////////////////////////////////////////////////////////
// fillHIEventDataSet(used only for counting number of events in each centrality class)
///////////////////////////////////////////////////////////////////////////////
RooDataSet* fillHIWtauEventDataSet(const TString& pathName, const TString& fileName, RooArgSet& centralityArgSet)
{
  RooDataSet* set = new RooDataSet(fileName,fileName,centralityArgSet);
  
  float centralityNt;

  // --- Load the MC tree ---
  TChain* tree = new TChain("truth","truth");
  int nFiles = tree->Add(pathName+fileName);
  
 
  // --- Set branch adresses ---
  tree->SetBranchAddress("centrality", &centralityNt);
  tree->SetBranchStatus("centrality", 1);
  
  std::cout << "Number of entries: " << tree->GetEntries() << std::endl; 
  for ( int i = 0; i < tree->GetEntries(); i++ ) {
//     if (i>100000) return set; // temp hack
    tree->LoadTree(i);
    tree->GetEntry(i);

    centralityArgSet.setRealValue("centrality",centralityNt);
    set->add(centralityArgSet);
  }

  return set;
}

TH1F* fillWeightedWtauMt(TString pathName, TString fileName, int nBins, float max, TGraphAsymmErrors* grEff){

    TChain* tree = new TChain("truth","truth");
    int nFiles = tree->Add(pathName+fileName);
  
    TH1F* hWtdMt = new TH1F("hWtdMt","hWtdMt",nBins,0.0,max);
    float mc_pt,mc_eta, mc_phi, mc_charge;
    float mc_pdgId,mc_mptPt,mc_mt,centralityNt;
    tree->SetBranchAddress("mc_pdgId",&mc_pdgId);
    tree->SetBranchAddress("mc_pt",&mc_pt);
    tree->SetBranchAddress("mc_eta",&mc_eta);
    tree->SetBranchAddress("mc_mptPt",&mc_mptPt);
    tree->SetBranchAddress("mc_mt",&mc_mt);
    tree->SetBranchAddress("mc_phi",&mc_phi);
    tree->SetBranchAddress("centrality",&centralityNt);
    tree->SetBranchAddress("mc_charge",&mc_charge);

    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("mc_pdgId",1);
    tree->SetBranchStatus("mc_pt",1);
    tree->SetBranchStatus("mc_eta",1);
    tree->SetBranchStatus("mc_mptPt",1);
    tree->SetBranchStatus("mc_mt",1);
    tree->SetBranchStatus("centrality",1);
    tree->SetBranchStatus("mc_charge",1);
    tree->SetBranchStatus("mc_phi",1);

    for(int igen=0; igen<tree->GetEntries(); ++igen){
        tree->GetEntry(igen);

        ///Smear generator pt by 0.1 GeV
        //float ptSmeared = smearGeneratorPt(mc_pt, 0.1);
        double eff = getEfficiency(mc_pt, grEff);
        //std::cout << "efficiency : " << eff << std::endl;
        if(fabs(mc_pdgId)==13
                && fabs(mc_eta)>0.1
                &&fabs(mc_eta)<2.4
                &&mc_mptPt>25.0
                &&mc_pt>25.0
                &&centralityNt<0.8
           ){
        
            hWtdMt->Fill(mc_mt,eff);
         }


    }

    return hWtdMt;
    
}

TH1F* fillWeightedWtauMpt(TString pathName, TString fileName, int nBins, float max, TGraphAsymmErrors* grEff){

    TChain* tree = new TChain("truth","truth");
    int nFiles = tree->Add(pathName+fileName);
  
    TH1F* hWtdMpt = new TH1F("hWtdMpt","hWtdMpt",nBins,0.0,max);
    float mc_pt,mc_eta, mc_phi, mc_charge;
    float mc_pdgId,mc_mptPt,mc_mt,centralityNt;
    tree->SetBranchAddress("mc_pdgId",&mc_pdgId);
    tree->SetBranchAddress("mc_pt",&mc_pt);
    tree->SetBranchAddress("mc_eta",&mc_eta);
    tree->SetBranchAddress("mc_mptPt",&mc_mptPt);
    tree->SetBranchAddress("mc_mt",&mc_mt);
    tree->SetBranchAddress("mc_phi",&mc_phi);
    tree->SetBranchAddress("centrality",&centralityNt);
    tree->SetBranchAddress("mc_charge",&mc_charge);

    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("mc_pdgId",1);
    tree->SetBranchStatus("mc_pt",1);
    tree->SetBranchStatus("mc_eta",1);
    tree->SetBranchStatus("mc_mptPt",1);
    tree->SetBranchStatus("mc_mt",1);
    tree->SetBranchStatus("centrality",1);
    tree->SetBranchStatus("mc_charge",1);
    tree->SetBranchStatus("mc_phi",1);

    for(int igen=0; igen<tree->GetEntries(); ++igen){
        tree->GetEntry(igen);

        ///Smear generator pt by 0.1 GeV
        //float ptSmeared = smearGeneratorPt(mc_pt, 0.1);
        double eff = getEfficiency(mc_pt,grEff);
        if(fabs(mc_pdgId)==13
                && fabs(mc_eta)>0.1
                &&fabs(mc_eta)<2.4
                &&mc_mt>40.0
                &&mc_pt>25.0
                &&centralityNt<0.8
           ){
        
            hWtdMpt->Fill(mc_mptPt,eff);
         }


    }

    return hWtdMpt;
    
}


TH1F* applyCutEfficiency(TH1F* hTau, TGraphAsymmErrors* grEff){

    TH1F* hNew = (TH1F*)hTau->Clone("hNew"); 
    ///array of xBin Centers
    double* xPt = grEff->GetX();
    double* xLo = grEff->GetEXlow();
    double* xUp = grEff->GetEXhigh();
    ///array of efficiencies
    double* yEff = grEff->GetY();
    ///Multiply by efficiency corresponding to bin in hTau
    for(int iEff = 0; iEff<grEff->GetN(); ++iEff){

        ///get bin center
        double binCenter = xPt[iEff];
        ///lower end of bin
        double binLo = xPt[iEff] - xLo[iEff];
        ///upper end of bin
        double binUpp = xPt[iEff] + xUp[iEff];

        for(int itau=1; itau<=hTau->GetNbinsX(); ++itau){
            
            ///Get x bin center
            double tauBin = hTau->GetXaxis()->GetBinCenter(itau);
            ///find corresponding bin in hTau
            if(tauBin>binLo&&tauBin<binUpp){
                double nominalValue = hTau->GetBinContent(itau);
                double eff = yEff[iEff];
                float scaledValue = nominalValue*eff;
                hNew->SetBinContent(itau,scaledValue);
                ///At the moment, assuming eff has no error
                hNew->SetBinError(itau,hTau->GetBinError(itau)*eff);
            }
        } //itau
    } //iEff

    return hNew;
}

RooDataSet* fillGeneratorTauMuonDataSet(const TString& pathName, const TString& fileName, RooArgSet& muonArgSet){

    
  RooDataSet* set = new RooDataSet(fileName,fileName,muonArgSet);
  TChain* tree = new TChain("truth","truth");
  int nFiles = tree->Add(pathName+fileName);
  
  std::cout << "Filling the RooDataSet for " << fileName << "... Number of files: " << nFiles << std::endl;
 
    float mc_pt,mc_eta, mc_phi, mc_charge;
    float mc_pdgId,mc_mptPt,mc_mt,centralityNt;
    tree->SetBranchAddress("mc_pdgId",&mc_pdgId);
    tree->SetBranchAddress("mc_pt",&mc_pt);
    tree->SetBranchAddress("mc_eta",&mc_eta);
    tree->SetBranchAddress("mc_mptPt",&mc_mptPt);
    tree->SetBranchAddress("mc_mt",&mc_mt);
    tree->SetBranchAddress("mc_phi",&mc_phi);
    tree->SetBranchAddress("centrality",&centralityNt);
    tree->SetBranchAddress("mc_charge",&mc_charge);

    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("mc_pdgId",1);
    tree->SetBranchStatus("mc_pt",1);
    tree->SetBranchStatus("mc_eta",1);
    tree->SetBranchStatus("mc_mptPt",1);
    tree->SetBranchStatus("mc_mt",1);
    tree->SetBranchStatus("centrality",1);
    tree->SetBranchStatus("mc_charge",1);
    tree->SetBranchStatus("mc_phi",1);

    for(int igen=0; igen<tree->GetEntries(); ++igen){
        tree->GetEntry(igen);

        muonArgSet.setRealValue("missPt",mc_mptPt);
        muonArgSet.setRealValue("muonMt",mc_mt);
        muonArgSet.setRealValue("muonEta",mc_eta);
        muonArgSet.setRealValue("muonPhi",mc_phi);
        muonArgSet.setRealValue("pdgID",mc_pdgId);
        muonArgSet.setRealValue("centrality",centralityNt);
        ///Smear generator pt by 0.1 GeV
        //float ptSmeared = smearGeneratorPt(mc_pt, 0.1);
        muonArgSet.setRealValue("muonPt",mc_pt);

        if ( mc_charge > 0. ) muonArgSet.setCatLabel("chargeCategory","muPlus");
        else if ( mc_charge < 0.) muonArgSet.setCatLabel("chargeCategory","muMinus");
        set->add(muonArgSet);   
    }

    return set;
}

///////////////////////////////
//plot
//////////////////////////////
void plot(RooDataSet* mcWSet,RooDataSet* mcWtauSet, RooDataSet* mcZSet, RooDataSet* mcJ1Set, RooDataSet* mcJ2Set, RooDataSet* mcJ3Set, 
		RooDataSet* mcWEvents, RooDataSet* mcTauEvents, RooDataSet* mcZEvents, RooDataSet* mcJ1Events,RooDataSet* mcJ2Events, RooDataSet* mcJ3Events, 
        const int iMt, const int iEta, const int iCentrality, int nCentralityBins, 
		double ncoll, double ptLow, double ptUpp, double etaLow, double etaUpp, 
		double centralityLow, double centralityUpp, RooRealVar& muonVar,RooRealVar& centrality, 
        TString sSel, TString sSel2 , TString xAxisLabel, int nBins=100, int max=200.0, double cutLine = 25.0) {


  	RooBinning b = RooBinning(nBins,0.0,max); // 4 GeV per bin

	TH1F* hmcWSet = (TH1F*)mcWSet->createHistogram("hmcWSet",muonVar,Binning(b));

  	// --- Wtau set ---
	TH1F* hmcTauSet = (TH1F*)mcWtauSet->createHistogram("hmcTauSet",muonVar,Binning(b));

    ///apply cut efficiency to Wtau
    //if(iMt==0) hmcTauSet = fillWeightedWtauMpt("/mnt/Lustre/cgrp/atlas_hi/tbalestri/", "MonteCarloFiles/HIWtaumuNtuple.07.09.2013.root",nBins,max,hCutEff); 
    //if(iMt==1) hmcTauSet = applyCutEfficiency(hmcTauSet, hCutEff);
    //if(iMt==2) hmcTauSet = fillWeightedWtauMt("/mnt/Lustre/cgrp/atlas_hi/tbalestri/", "MonteCarloFiles/HIWtaumuNtuple.07.09.2013.root",nBins,max,hCutEff); 

  	// --- Z set ---
	TH1F* hmcZSet = (TH1F*)mcZSet->createHistogram("hmcZSet",muonVar,Binning(b));
  	// --- QCD set ---
	TH1F* hmcJ1Set = (TH1F*)mcJ1Set->createHistogram("hmcJ1Set",muonVar,Binning(b));
	TH1F* hmcJ2Set = (TH1F*)mcJ2Set->createHistogram("hmcJ2Set",muonVar,Binning(b));
	TH1F* hmcJ3Set = (TH1F*)mcJ3Set->createHistogram("hmcJ3Set",muonVar,Binning(b));

	//return weighted histograms
    hmcWSet = getWeightedWHisto(centralityLow, centralityUpp, ncoll,hmcWSet, mcWEvents, nBins,0.0, max); 
    hmcZSet = getWeightedZHisto(centralityLow, centralityUpp, ncoll, hmcZSet, mcZEvents, nBins,0.0, max); 
    hmcTauSet = getWeightedTauHisto(centralityLow, centralityUpp, ncoll, hmcTauSet, mcTauEvents, nBins,0.0, max); 
	TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",nBins,0.0,max);
	hmcQCDSet = getWeightedJxHisto(centralityLow, centralityUpp, ncoll, hmcJ1Set, hmcJ2Set, hmcJ3Set, mcJ1Events, mcJ2Events,mcJ3Events, nBins,0.0, max);

    ///W-->tau--mu events per W-->mu event integrated over all eta and Ncoll
    double bkgFracWtau = 0.0175;
    int sigBinLo = hmcWSet->FindBin(cutLine);
    if(iMt==1) std::cout << "Number of Wtau events in signal phase space after smearing and reweighting: " << hmcTauSet->Integral(sigBinLo,nBins) << std::endl; 
    ///Normalize to expected number of Wtaumu events
    hmcTauSet->Scale(bkgFracWtau*hmcWSet->Integral(sigBinLo,nBins)/hmcTauSet->Integral(sigBinLo,nBins));

    ///Scale factor determined from 
    ///scaling to data in ctrl region (0.1-2.4,0-80%);
    hmcQCDSet->Scale(0.4);

	TH1F* hmcQCDSetc = (TH1F*)hmcQCDSet->Clone("hmcQCDSetc");
	TH1F* hmcZSetc = (TH1F*)hmcZSet->Clone("hmcZSetc");
	TH1F* hmcWSetc = (TH1F*)hmcWSet->Clone("hmcWSetc");
	TH1F* hmcTauSetc = (TH1F*)hmcTauSet->Clone("hmcTauSetc");

	//hmcQCDSetc->SetLineColor(kAzure-9);
	hmcQCDSetc->SetMarkerStyle(23);
	hmcQCDSetc->SetMarkerColor(kAzure-9);
	//hmcZSetc->SetLineColor(kRed);
	hmcZSetc->SetMarkerStyle(29);
	hmcZSetc->SetMarkerColor(kRed);
    //hmcTauSetc->SetLineColor(kYellow-9);
    hmcTauSetc->SetMarkerColor(kGreen);
    hmcTauSetc->SetMarkerStyle(22);

  	TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);

	leg->AddEntry(hmcWSetc, "W#rightarrow#mu#nu", "lpe");
	leg->AddEntry(hmcTauSetc, "W#rightarrow#tau", "lpe");
	leg->AddEntry(hmcZSetc, "Z#rightarrow#mu#mu", "lpe");
	leg->AddEntry(hmcQCDSetc, "QCD", "lpe");

  	TCanvas* cdatamt = new TCanvas("cdatamt","cdatamt",600,600);
  	hmcWSetc->GetXaxis()->SetTitle(xAxisLabel); 
  	if(iMt==0)hmcWSetc->GetYaxis()->SetTitle("dN/d#slash{p_{T}}1/N_{MC} #sigma_{i}/#sigma_{pp} #LT N_{coll} #GT N_{data} #Delta_{bw}");
  	if(iMt==1)hmcWSetc->GetYaxis()->SetTitle("dN/dp_{T} 1/N_{MC} #sigma_{i}/#sigma_{pp} #LT N_{coll} #GT N_{data} #Delta_{bw}");
  	if(iMt==2)hmcWSetc->GetYaxis()->SetTitle("dN/dm_{T} 1/N_{MC} #sigma_{i}/#sigma_{pp} #LT N_{coll} #GT N_{data} #Delta_{bw}");
	hmcWSetc->GetXaxis()->SetRangeUser(cutLine,max); 
	//hmcWSetc->GetXaxis()->SetRangeUser(0.0,max); 
    ///normalize to W MC
    double normQCD = hmcWSetc->Integral()/hmcQCDSetc->Integral();
    std::cout << "QCD norm factor: " << normQCD << std::endl;
    //hmcQCDSetc->Scale(normQCD);
    double normZ = hmcWSetc->Integral()/hmcZSetc->Integral();
    std::cout << "Z norm factor: " << normZ << std::endl;
    //hmcZSetc->Scale(normZ);
	hmcWSetc->Draw("hist pe");
	hmcTauSetc->Draw("hist pesame");
	hmcZSetc->Draw("hist pesame");
	hmcQCDSetc->Draw("hist pesame");
	hmcWSetc->Draw("sameaxis");

 	//TLine *line0 = new TLine(cutLine,0.01,cutLine,881201);
 	TLine *line0 = new TLine(cutLine,0.05,cutLine,1.0799e7);
	line0->SetLineColor(kBlack);
	line0->SetLineStyle(kDashed);
	line0->SetLineWidth(2);
	//line0->Draw();

    TString plotNameLog,plotNameLin;
	TLatex l;
	l.SetNDC();
	l.DrawLatex(0.67,0.66,sSel + "%" );
	l.DrawLatex(0.64,0.59,sSel2);
    TString sSel3;
    if (iMt==0) { sSel3 = "p_{T}^{#mu}>25";plotNameLog = "mpt_"; plotNameLin = plotNameLog;}
    else if (iMt==1) {sSel3 = "";/*"#slash{p_{T}}>25";*/plotNameLog = "pt_"; plotNameLin = plotNameLog;}
    else if(iMt==2) {sSel3 = "p_{T}^{#mu}>25,#slash{p_{T}}>25"; plotNameLog = "mt_"; plotNameLin = plotNameLog;}
	l.DrawLatex(0.70,0.5,sSel3);
	l.SetTextSize(0.034);
//	l.DrawLatex(0.74,0.89,"#sqrt{s_{NN}}=2.76 TeV");
//	l.DrawLatex(0.492,0.89,"#int Ldt #approx 0.140 nb^{-1}"); 

	leg->Draw(); cdatamt->Update();

	plotNameLog += "kinematicsRaw"; plotNameLog+="charge"; plotNameLog+=iMt; plotNameLog+="_eta";
    plotNameLog+=iEta; plotNameLog+="_cent"; plotNameLog+=iCentrality;  
	plotNameLin+="kinematicsRaw"; plotNameLin+="charge"; plotNameLin+=iMt; plotNameLin+="_eta"; plotNameLin+=iEta;
    plotNameLin+="_cent"; plotNameLin+=iCentrality;  

    hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()+50.0); cdatamt->Update();
//	cdatamt->Print(plotNameLin.ReplaceAll("#","")+"_Lin.png");
//	cdatamt->Print(plotNameLin.ReplaceAll("#","")+"_Lin.eps"); 
//	cdatamt->Print(plotNameLin.ReplaceAll("#","")+"_Lin.pdf"); 
	cdatamt->Print(plotNameLin+".root"); 
    cdatamt->SetLogy(true); hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()*1.0e1); cdatamt->Update();
	cdatamt->Print(plotNameLog.ReplaceAll("#","")+"_Log.png");
//	cdatamt->Print(plotNameLog.ReplaceAll("#","")+"_Log.eps"); 
	cdatamt->Print(plotNameLog.ReplaceAll("#","")+"_Log.pdf"); 
	TString plotNameLogRoot = plotNameLog.ReplaceAll("|",",") + ".root";
	cdatamt->Print(plotNameLogRoot); 

	std::cout << "Clean up" << std::endl;

	delete hmcQCDSet;
    delete leg;
	delete cdatamt;	
    delete line0;


}//plot

void plotMcMptPtMtDistros(){

	bool doCharge = false ;
	bool doCentrality = false ;
	bool doEta = false; 
    if(doCharge||doEta||doCentrality){ std::cout << "Cannot bin in eta,centrality,and or charge. " << std::endl; exit(0)}
	bool cutValue = 11; 

	float mtmax = 300.0;
	float ptmax = 300.0;

	SetAtlasStyle();
	TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";

    ///Wmunu
	TString fileNameMCWIn = "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.04.12.2013";

    ///Zmumu
	TString fileNameMCZIn = "MonteCarloFiles/Zmumu/HISingleMuonZmumuPYTHIAHIJINGOverlay.04.13.2013";

	//J1 1 muon-filter 
	TString fileNameMCJ1In = "MonteCarloFiles/QCD/HISingleMuonJ1_1muPYTHIADataOverlay.04.13.2013";
	//J2 1 muon-filter 
	TString fileNameMCJ2In = "MonteCarloFiles/QCD/HISingleMuonJ2_1muPYTHIADataOverlay.04.13.2013";
	//J3 1 muon-filter 
	TString fileNameMCJ3In = "MonteCarloFiles/QCD/HISingleMuonJ3_1muPYTHIADataOverlay.04.13.2013";
    ///Files holding cut efficiencies as a function of pt,mpt,and mt
    //TFile* fEff = new TFile("mcWEffForTauMuonStudy_07_08_2013.root","read");
    //TGraphAsymmErrors* hCutEffPt = (TGraphAsymmErrors*)fEff->Get("pEffWAccCuts");

    ///W-->tau-->mu Generated sample
    TString fileNameMCTauIn = "MonteCarloFiles/Wtaumu/HIWtaumuNtuple.07.25.2013";

    RooRealVar  muonPt("muonPt","p_{T}",0.0,350.0,"GeV");
	RooRealVar  missPt("missPt","p_{T}^{miss}",0.0,250.0,"GeV");
	RooRealVar  muonMt("muonMt","m_{T}",0.0,350.0,"GeV");
	RooRealVar  muonCharge("muonCharge","charge",-2.0,+2.0);
	RooRealVar  isolationMu("isolationMu","isolationMu",0.0,10.0);
  	RooRealVar  centrality("centrality","centrality",0.,1.0);
  	RooRealVar  muonEta("muonEta","muonEta",-3.0,+3.0);
  	RooRealVar  muonPhi("muonPhi","muonPhi",-3.5,+3.5);
    RooRealVar  motherRec("motherRec","motherRec",0.0,50.0,"GeV");
  	RooRealVar  ZDY("ZDY","ZDY",0.0,1.0);
    RooRealVar  muonQuality("muonQuality","muonQuality",0.0,15.0);
    RooRealVar  muonELoss("muonELoss","muonELoss",-10.0,+10.0);
    RooRealVar  muonScat("muonScat","muonScat",-10.0,+10.0);

    RooArgList muonRecArgList(muonEta,muonPt,missPt,muonMt);
    muonRecArgList.add(muonQuality); 
    muonRecArgList.add(muonELoss);
    muonRecArgList.add(muonScat);
    muonRecArgList.add(ZDY);
    muonRecArgList.add(isolationMu);

    RooArgList muonWtauArgList(muonEta,muonPt,missPt,muonMt);


	TString sCutsSigMpt =
          "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>0.0&&missPt<9000.0&&isolationMu<0.1&&muonMt>=40.0&&ZDY==0";

	RooFormulaVar cutsSigMpt("cutsSigMpt", "cutsSigMpt", sCutsSigMpt, muonRecArgList);
	TString sCutsSigPt =
          "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>0.0&&missPt>=25.0&&missPt<9000.0&&isolationMu<0.1&&muonMt>=40.0&&ZDY==0";
	RooFormulaVar cutsSigPt("cutsSigPt", "cutsSigPt", sCutsSigPt, muonRecArgList);
	TString sCutsSigMt =
          "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&isolationMu<0.1&&muonMt>=0.0&&ZDY==0";
	RooFormulaVar cutsSigMt("cutsSigMt", "cutsSigMt", sCutsSigMt, muonRecArgList);


    ///Cuts for generator-level tau sample
	TString sCutsWtauMpt =
          "abs(muonEta)>0.1&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>=0.0&&missPt<9000.0&&muonMt>=40.0";
	RooFormulaVar cutsWtauMpt("cutsWtauMpt", "cutsWtauMpt", sCutsWtauMpt, muonWtauArgList);

	TString sCutsWtauPt =
          "abs(muonEta)>0.1&&abs(muonEta)<2.4&&muonPt>0.0&&missPt>=25.0&&missPt<9000.0&&muonMt>=40.0";
	RooFormulaVar cutsWtauPt("cutsWtauPt", "cutsWtauPt", sCutsWtauPt, muonWtauArgList);

	TString sCutsWtauMt =
          "abs(muonEta)>0.1&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>=25.0&&missPt<9000.0&&muonMt>=0.0";
	RooFormulaVar cutsWtauMt("cutsWtauMt", "cutsWtauMt", sCutsWtauMt, muonWtauArgList);


	RooCategory muonCategory("muonCategory","muonCategory");
	muonCategory.defineType("Unknown",0);
	muonCategory.defineType("HeavyFlavour",1);
	muonCategory.defineType("Tau",2);
	muonCategory.defineType("LightResonance",3); 
	muonCategory.defineType("PionCaloDecay",4);
	muonCategory.defineType("PionEarlyDecay",5);
	muonCategory.defineType("PionLateDecay",6);
	muonCategory.defineType("KaonCaloDecay",7);
	muonCategory.defineType("KaonEarlyDecay",8);
	muonCategory.defineType("KaonLateDecay",9);
	muonCategory.defineType("Fake",10);
	muonCategory.defineType("Z",23);
	muonCategory.defineType("W",24);

  	RooCategory chargeCategory("chargeCategory","chargeCategory") ;
  	chargeCategory.defineType("muMinus",-1) ;
  	chargeCategory.defineType("muPlus",1) ;

  	RooArgSet muonArgSet(muonPt,missPt,isolationMu,muonMt,muonEta,centrality,ZDY,muonCategory,chargeCategory);
	muonArgSet.add(muonPhi);
	muonArgSet.add(muonQuality);
	muonArgSet.add(muonELoss);
	muonArgSet.add(muonScat);
    muonArgSet.add(motherRec);

    RooArgSet centralityArgSet(centrality);  

	// --- Set pt and eta bins ---
	std::vector<double> ptBins;
	ptBins.push_back(0.0);
	//ptBins.push_back(33.5);
	ptBins.push_back(ptmax);
	const int nPtBins = ptBins.size()-1;

	std::vector<double> etaBins;
	etaBins.push_back(0.10);
	if (doEta) {
        etaBins.push_back(0.35);
        etaBins.push_back(0.6);
        etaBins.push_back(0.8);
        etaBins.push_back(1.05);
        etaBins.push_back(1.3);
        etaBins.push_back(1.55);
//        etaBins.push_back(1.73);
        etaBins.push_back(1.85);
        etaBins.push_back(2.1);

	}
	etaBins.push_back(+2.40);

	const int nEtaBins = etaBins.size()-1;
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


  	// --- Fill mc sets ---
    RooDataSet* mcWSet = fillHIMuonDataSet(baseString,fileNameMCWIn+".root",muonArgSet, cutValue, true); mcWSet->Print();
  	mcWSet = (RooDataSet*)mcWSet->reduce(Cut("motherRec==24")); mcWSet->Print();
	RooDataSet* mcWSetMpt = (RooDataSet*)mcWSet->reduce(Cut(cutsSigMpt)); mcWSetMpt->Print();
	RooDataSet* mcWSetPt = (RooDataSet*)mcWSet->reduce(Cut(cutsSigPt)); mcWSetPt->Print();
	RooDataSet* mcWSetMt = (RooDataSet*)mcWSet->reduce(Cut(cutsSigMt)); mcWSetMt->Print();
    RooDataSet* mcWEvents = fillHIEventDataSet(baseString,fileNameMCWIn+".root",centralityArgSet );

    RooDataSet* mcWtauSet = fillHIWTauDataSet(baseString,fileNameMCTauIn+".root",muonArgSet); mcWtauSet->Print();
	RooDataSet* mcWtauSetMpt = (RooDataSet*)mcWtauSet->reduce(Cut(cutsWtauMpt)); mcWtauSetMpt->Print();
	RooDataSet* mcWtauSetPt = (RooDataSet*)mcWtauSet->reduce(Cut(cutsWtauPt)); mcWtauSetPt->Print();
	RooDataSet* mcWtauSetMt = (RooDataSet*)mcWtauSet->reduce(Cut(cutsWtauMt)); mcWtauSetMt->Print();
    RooDataSet* mcWtauEvents = fillHIEventDataSet(baseString,fileNameMCTauIn+".root",centralityArgSet );

    RooDataSet* mcZSet = fillHIMuonDataSet(baseString,fileNameMCZIn+".root",muonArgSet, cutValue, true); mcZSet->Print();
  	mcZSet = (RooDataSet*)mcZSet->reduce(Cut("motherRec==23")); mcZSet->Print();
	RooDataSet* mcZSetMpt = (RooDataSet*)mcZSet->reduce(Cut(cutsSigMpt)); mcZSetMpt->Print();
	RooDataSet* mcZSetPt = (RooDataSet*)mcZSet->reduce(Cut(cutsSigPt)); mcZSetPt->Print();
	RooDataSet* mcZSetMt = (RooDataSet*)mcZSet->reduce(Cut(cutsSigMt)); mcZSetMt->Print();
    RooDataSet* mcZEvents = fillHIEventDataSet(baseString,fileNameMCZIn+".root",centralityArgSet );

    RooDataSet* mcJ1Set = fillHIMuonDataSet(baseString,fileNameMCJ1In+".root",muonArgSet, cutValue, true); mcJ1Set->Print();
	RooDataSet* mcJ1SetMpt = (RooDataSet*)mcJ1Set->reduce(Cut(cutsSigMpt)); mcJ1SetMpt->Print();
	RooDataSet* mcJ1SetPt = (RooDataSet*)mcJ1Set->reduce(Cut(cutsSigPt)); mcJ1SetPt->Print();
	RooDataSet* mcJ1SetMt = (RooDataSet*)mcJ1Set->reduce(Cut(cutsSigMt)); mcJ1SetMt->Print();
    ///number of events in this centrality bin 
    RooDataSet* mcJ1Events = fillHIEventDataSet(baseString,fileNameMCJ1In+".root",centralityArgSet );

    RooDataSet* mcJ2Set = fillHIMuonDataSet(baseString,fileNameMCJ2In+".root",muonArgSet, cutValue, true); mcJ2Set->Print();
	RooDataSet* mcJ2SetMpt = (RooDataSet*)mcJ2Set->reduce(Cut(cutsSigMpt)); mcJ2SetMpt->Print();
	RooDataSet* mcJ2SetPt = (RooDataSet*)mcJ2Set->reduce(Cut(cutsSigPt)); mcJ2SetPt->Print();
	RooDataSet* mcJ2SetMt = (RooDataSet*)mcJ2Set->reduce(Cut(cutsSigMt)); mcJ2SetMt->Print();
    ///number of events in this centrality bin 
    RooDataSet* mcJ2Events = fillHIEventDataSet(baseString,fileNameMCJ2In+".root",centralityArgSet );

    RooDataSet* mcJ3Set = fillHIMuonDataSet(baseString,fileNameMCJ3In+".root",muonArgSet, cutValue, true); mcJ3Set->Print();
	RooDataSet* mcJ3SetMpt = (RooDataSet*)mcJ3Set->reduce(Cut(cutsSigMpt)); mcJ3SetMpt->Print();
	RooDataSet* mcJ3SetPt = (RooDataSet*)mcJ3Set->reduce(Cut(cutsSigPt)); mcJ3SetPt->Print();
	RooDataSet* mcJ3SetMt = (RooDataSet*)mcJ3Set->reduce(Cut(cutsSigMt)); mcJ3SetMt->Print();
    ///number of events in this centrality bin 
    RooDataSet* mcJ3Events = fillHIEventDataSet(baseString,fileNameMCJ3In+".root",centralityArgSet );

	// --- Subdivide in bins ---
	RooDataSet* mcWSubSetMpt[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcWSubSetPt[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcWSubSetMt[nPtBins][nEtaBins][nCentralityBins];

	RooDataSet* mcWtauSubSetMpt[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcWtauSubSetPt[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcWtauSubSetMt[nPtBins][nEtaBins][nCentralityBins];

	RooDataSet* mcZSubSetMpt[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcZSubSetPt[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcZSubSetMt[nPtBins][nEtaBins][nCentralityBins];

	RooDataSet* mcJ1SubSetMpt[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcJ1SubSetPt[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcJ1SubSetMt[nPtBins][nEtaBins][nCentralityBins];

	RooDataSet* mcJ2SubSetMpt[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcJ2SubSetPt[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcJ2SubSetMt[nPtBins][nEtaBins][nCentralityBins];

	RooDataSet* mcJ3SubSetMpt[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcJ3SubSetPt[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcJ3SubSetMt[nPtBins][nEtaBins][nCentralityBins];

    std::cout << "Done initializing dataset. Creating subsets in pt,eta,and centrality..." <<std::endl;
	for ( int i = 0; i < nPtBins; i++ ) {
	  for ( int j = 0; j < nEtaBins; j++ ) {
	    for ( int k = 0; k < nCentralityBins; k++ ){

        	mcWSubSetMpt[i][j][k] = selectPtEtaCentrality(mcWSetMpt , ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 
        	mcWSubSetPt[i][j][k] = selectPtEtaCentrality(mcWSetPt , ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 
        	mcWSubSetMt[i][j][k] = selectPtEtaCentrality(mcWSetMt , ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 

        	mcWtauSubSetMpt[i][j][k] = selectPtEtaCentrality(mcWtauSetMpt , ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 
        	mcWtauSubSetPt[i][j][k] = selectPtEtaCentrality(mcWtauSetPt , ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 
        	mcWtauSubSetMt[i][j][k] = selectPtEtaCentrality(mcWtauSetMt , ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 

        	mcZSubSetMpt[i][j][k] = selectPtEtaCentrality(mcZSetMpt, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 
        	mcZSubSetPt[i][j][k] = selectPtEtaCentrality(mcZSetPt, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 
        	mcZSubSetMt[i][j][k] = selectPtEtaCentrality(mcZSetMt, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 

            std::cout << "Jx shaping.." << std::endl;
		    //use |eta|<2.5 for Jx since mT shape has no eta dep 
            mcJ1SetMpt->Print();
		    mcJ1SubSetMpt[i][j][k] = selectPtEtaCentrality(mcJ1SetMpt, ptBins[i], ptBins[i+1], etaBins[0], etaBins[nEtaBins], centralityBins[k], centralityBins[k+1],true); 
		    mcJ1SubSetPt[i][j][k] = selectPtEtaCentrality(mcJ1SetPt, ptBins[i], ptBins[i+1], etaBins[0], etaBins[nEtaBins], centralityBins[k], centralityBins[k+1],true); 
		    mcJ1SubSetMt[i][j][k] = selectPtEtaCentrality(mcJ1SetMt, ptBins[i], ptBins[i+1], etaBins[0], etaBins[nEtaBins], centralityBins[k], centralityBins[k+1],true); 
            mcJ2SubSetMpt[i][j][k] = selectPtEtaCentrality(mcJ2SetMpt, ptBins[i], ptBins[i+1], etaBins[0], etaBins[nEtaBins], centralityBins[k], centralityBins[k+1],true); 
            mcJ2SubSetPt[i][j][k] = selectPtEtaCentrality(mcJ2SetPt, ptBins[i], ptBins[i+1], etaBins[0], etaBins[nEtaBins], centralityBins[k], centralityBins[k+1],true); 
            mcJ2SubSetMt[i][j][k] = selectPtEtaCentrality(mcJ2SetMt, ptBins[i], ptBins[i+1], etaBins[0], etaBins[nEtaBins], centralityBins[k], centralityBins[k+1],true); 
            mcJ3SubSetMpt[i][j][k] = selectPtEtaCentrality(mcJ3SetMpt, ptBins[i], ptBins[i+1], etaBins[0], etaBins[nEtaBins], centralityBins[k], centralityBins[k+1],true); 
		    mcJ3SubSetPt[i][j][k] = selectPtEtaCentrality(mcJ3SetPt, ptBins[i], ptBins[i+1], etaBins[0], etaBins[nEtaBins], centralityBins[k], centralityBins[k+1],true); 
		    mcJ3SubSetMt[i][j][k] = selectPtEtaCentrality(mcJ3SetMt, ptBins[i], ptBins[i+1], etaBins[0], etaBins[nEtaBins], centralityBins[k], centralityBins[k+1],true); 


	    }
	  }
	}
    std::cout << "Done." <<std::endl;

	for ( int i = 0; i < nPtBins; i++ ) {
	  for ( int j = 0; j < nEtaBins; j++ ) {
	    for ( int k = 0; k < nCentralityBins; k++ ){
        	
            TString sCentLow = "";
            TString sCentUp = "";
            TString sEtaLow = "";
            TString sEtaUp = "";

            sCentLow += 100*centralityBins[k]; //sCentLow.Remove(3);
            sCentUp += 100*centralityBins[k+1]; //sCentUp.Remove(3);

            sEtaLow += etaBins[j];
            sEtaUp += etaBins[j+1];
		
           if(!doCentrality && !doEta && !doCharge){

			plot(mcWSubSetMpt[i][j][k],mcWtauSubSetMpt[i][j][k],mcZSubSetMpt[i][j][k],mcJ1SubSetMpt[i][j][k],mcJ2SubSetMpt[i][j][k],mcJ3SubSetMpt[i][j][k],
				mcWEvents,mcWtauEvents,mcZEvents,mcJ1Events, mcJ2Events, mcJ3Events, 
				0, 0, 0, nCentralityBins, ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],
                missPt,centrality, 
				"#mu^{#pm},0-80" , "0.1 #leq |#eta| < 2.4","#slash{p_{T}}[GeV]",30,120);
    
    	    plot(mcWSubSetPt[i][j][k],mcWtauSubSetPt[i][j][k],mcZSubSetPt[i][j][k],mcJ1SubSetPt[i][j][k],mcJ2SubSetPt[i][j][k],mcJ3SubSetPt[i][j][k],
				mcWEvents,mcWtauEvents,mcZEvents,mcJ1Events, mcJ2Events, mcJ3Events, 
				1, 0, 0, nCentralityBins, ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],
                muonPt,centrality, 
				"#mu^{#pm},0-80" , "0.1 #leq |#eta| < 2.4","p_{T}^{#mu}[GeV]",120,120);
            
	        plot(mcWSubSetMt[i][j][k],mcWtauSubSetMt[i][j][k],mcZSubSetMt[i][j][k],mcJ1SubSetMt[i][j][k],mcJ2SubSetMt[i][j][k],mcJ3SubSetMt[i][j][k],
				mcWEvents,mcWtauEvents,mcZEvents,mcJ1Events, mcJ2Events, mcJ3Events, 
				2, 0, 0, nCentralityBins, ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],
                muonMt,centrality, 
				"#mu^{#pm},0-80" , "0.1 #leq |#eta| < 2.4","m_{T}[GeV]",50,200.0,40.0);

           	}
            
            else std::cout << "Cannot bin for centrality,eta, and/or charge. Must write more code." << std::endl;


        } //i
      } //j
   } //k
}
