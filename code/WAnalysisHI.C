#include "TCanvas.h"
#include "TDirectory.h"
#include "TString.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TList.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLatex.h"
#include "TClonesArray.h"
#include "RooCurve.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooNumber.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooDecay.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooMinuit.h"
#include "RooBinning.h"
#include "RooHistPdf.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooPoisson.h"
#include "RooBifurGauss.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooNumIntConfig.h"
//#include "RooShiftedKeysPdf.h"
// #include "Systematics.C"
#include "EfficiencyCorrection.C"
#include "WAnalysisHIDep.C"
//#include "flavourFitterHIDep.C"
//#include "AtlasUtils.h"
#include "AtlasUtils.C"
#include "AtlasStyle.C"

#include <iostream>
#include <iomanip>
#include <fstream>

////////////////////////////////////////////////////////////////////////////////
// fixplot
////////////////////////////////////////////////////////////////////////////////
void fixplot(RooPlot* fr)
{
 RooHist* hist0 = (RooHist*) fr->getObject(0) ;

 Double_t x,y ;
 for (Int_t i=0 ; i<hist0->GetN() ; i++) {
   hist0->GetPoint(i,x,y) ;
   if (y==0) {
     hist0->SetPointEYhigh(i,0) ;
   }
 }
}


///////////////////////////////////////////////////////////////////////////////
//plot  
///////////////////////////////////////////////////////////////////////////////
void plot( RooDataSet* dataSet, RooDataSet* mcWSet, RooDataSet* mcZSet,RooDataSet* mcQCDSet, RooRealVar& muonMt, 
          FitResult& fitResult, const int iMt, const int iEta, const int iCentrality, const float mtmax, TString fileNameDataIn, TString sSel, TString sSel2 , bool addPercent = true, bool correctSpectra = false)
{

  float mtcutLow = 40.; 

  RooPlot* frameMt=muonMt.frame(Range(0,mtmax)) ;
  frameMt->SetXTitle("#font[52]{m}_{T} [GeV]");

  //RooBinning b = RooBinning(mtmax,0,mtmax); // 1 GeV per bin
  RooBinning b = RooBinning(80,0,mtmax); // 2 GeV per bin
  double binwidth = mtmax/80.0;
  dataSet->plotOn(frameMt, RefreshNorm(), Binning(b), DrawOption("p"), MarkerSize(0.8),DataError(RooAbsData::SumW2));

  frameMt->Print();
  //RooHist* histData = 0;
  double effMt ;

  if(correctSpectra) {
  	correctEfficiencyMt(frameMt, iMt, iCentrality, iEta);
  	effMt = getEfficiencyMt(iMt, iCentrality, iEta);
  }
  else effMt = 1.0;
  //dataSet->plotOn(frameMt, Binning(b), DrawOption("p"), MarkerSize(0.8));
  
  muonMt.setRange("norm",0.0,mtmax) ;  
  mcWSet->plotOn(frameMt,Binning(b), LineColor(kBlack), FillColor(0), DrawOption("fl"), );
  mcQCDSet->plotOn(frameMt,Binning(b), LineColor(kAzure-9), FillColor(kAzure-9), DrawOption("f"),NormRange("norm"));
  mcZSet->plotOn(frameMt,RefreshNorm(), Binning(b),LineColor(kRed), FillColor(kRed), DrawOption("f"));
  dataSet->plotOn(frameMt, RefreshNorm(), Binning(b), DrawOption("p"), MarkerSize(0.8),DataError(RooAbsData::SumW2));


  frameMt->Print();

  muonMt.setRange("signal",mtcutLow,mtmax) ;  

  TString sAll = "(muonMt>0&&muonMt<"; sAll+=mtmax; sAll+=")";
  TString sSig = "(muonMt>="; sSig+= mtcutLow; sSig+="&&muonMt<="; sSig+=mtmax; sSig+=")";

//  double allEvents = dataSet->sumEntries(sAll);
//  double sigEvents = allEvents*int_fracSig;
//  double sigEventsCount = dataSet->sumEntries(sSig);


/*  RooHist* histData = (RooHist*) frameMt->findObject("h_HISingleMuonHP.2012.11.16.root") ;
  double sigEvents = histData->Integral(mtcutLow,mtmax);
*/
  TH1F* h_mtData  = new TH1F("h_mtData", "h_mtData", 80, 0.0, mtmax);
  dataSet->fillHistogram(h_mtData, RooArgList(muonMt));
  double sigEvents = h_mtData->Integral(mtcutLow,mtmax)*binwidth;


  std::cout << " integrated events mT:"<< mtcutLow << "-" << mtmax  << " =  " << sigEvents << " +-" << TMath::Sqrt(sigEvents)/effMt << std::endl;

  fitResult.setSig(iMt, iEta, iCentrality, sigEvents, TMath::Sqrt(sigEvents)/effMt, TMath::Sqrt(sigEvents)/effMt);

  TLegend* leg = new TLegend(0.586, 0.5, 0.857, 0.692);

  TH1F* hdummy0 = new TH1F(); hdummy0->SetMarkerSize(0.8);
  TH1F* hdummy1 = new TH1F();  
  TH1F* hdummy2 = new TH1F(); hdummy2->SetFillColor(kAzure-9); hdummy2->SetLineColor(kAzure-9);
  TH1F* hdummy3 = new TH1F(); hdummy3->SetFillColor(kRed); hdummy3->SetLineColor(kRed);

  leg->SetTextFont(gStyle->GetTextFont());
  leg->SetTextSize(gStyle->GetTextSize());
  leg->SetBorderSize(0);
  leg->SetFillColor(0);

  leg->AddEntry(hdummy0, "Data 2011", "pe");
  leg->AddEntry(hdummy1, "W#rightarrow#mu#nu", "f");
  leg->AddEntry(hdummy2, "QCD", "f");
  leg->AddEntry(hdummy3, "Z#rightarrow#mu#mu", "f");
    
  TCanvas* cdatamt = new TCanvas("cdatamt","datamt",600,600);
  cdatamt->UseCurrentStyle();
  //cdatamt->SetLogy(true);
  frameMt->SetMinimum(0.1);
  //frameMt->SetMinimum(0.5e2);
  frameMt->SetMaximum(1.9e3);
  std::cout << "now printing frameMt" << std::endl;
  frameMt->Print();

  std::cout << "now fixing plot" << std::endl;
  fixplot(frameMt);

  frameMt->Draw();

  leg->Draw();

  TLatex l;
  l.SetNDC();
  l.DrawLatex(0.28,0.75,sSel + ( addPercent ? "%" : "" ));
  l.DrawLatex(0.28,0.65,sSel2);
  l.DrawLatex(0.6,0.75,"#sqrt{s_{NN}}=2.76 TeV");
  l.DrawLatex(0.57,0.85,"#int Ldt #approx 0.140 nb^{-1}"); 

//#ifdef __CINT__
  ATLAS_LABEL(0.17,0.85,1);
  //myText(0.35,0.85, (Color_t)kBlack, (char*)("Preliminary"));
  myText(0.33,0.85, (Color_t)kBlack, (char*)("Internal"));
  //myText(0.35,0.85, (Color_t)kBlack, (char*)("Performance"));
  //myText(0.25,0.70, (Color_t)kRed, (char*)("For Approval"));
  //myText(0.25,0.70, (Color_t)kRed, (char*)("With Approved Parameters"));
//#endif

  std::cout << "now making file names" << std::endl;
  TString plotNameLog = "dataMt_"; plotNameLog+=sSel; plotNameLog+=","; plotNameLog+=sSel2; plotNameLog+=","; plotNameLog+="Log"; //.png";
  TString plotNameLin = "dataMt_"; plotNameLin+=sSel; plotNameLin+="_"; plotNameLin+=sSel2; plotNameLin+="_"; plotNameLin+="Lin"; //.png";
  std::cout << "now saving files" << std::endl;

  cdatamt->Print(plotNameLin.ReplaceAll("#","")+".png");
  cdatamt->Print(plotNameLin.ReplaceAll("#","")+".eps"); 
  cdatamt->Print(plotNameLin.ReplaceAll("#","")+".pdf"); 
  cdatamt->Print(plotNameLin+".root"); 
  
  TCanvas* cdatamtLog = new TCanvas("cdatamtLog","datamtLog",600,600);
  cdatamtLog->UseCurrentStyle();
  cdatamt->SetLogy(true);

  TLegend* leg2 = new TLegend(0.586, 0.5, 0.857, 0.692);
  leg2->SetTextFont(gStyle->GetTextFont());
  leg2->SetTextSize(gStyle->GetTextSize());
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->AddEntry(hdummy0, "Data 2011", "pe");
  leg2->AddEntry(hdummy1, "W#rightarrow#mu#nu", "f");
  leg2->AddEntry(hdummy2, "QCD", "f");
  leg2->AddEntry(hdummy3, "Z#rightarrow#mu#mu", "f");

  frameMt->SetMinimum(1e0) ;
  frameMt->SetMaximum(1e6);
  frameMt->Draw();

  leg2->Draw();

  TLatex lL;
  lL.SetNDC();
  lL.DrawLatex(0.28,0.75,sSel + ( addPercent ? "%" : "" ));
  lL.DrawLatex(0.28,0.65,sSel2);
  lL.DrawLatex(0.6,0.75,"#sqrt{s_{NN}}=2.76 TeV");
  lL.DrawLatex(0.57,0.85,"#int Ldt #approx 0.140 nb^{-1}"); 
  ATLAS_LABEL(0.17,0.85,1);
  myText(0.33,0.85, (Color_t)kBlack, (char*)("Internal"));
 
  cdatamt->Print(plotNameLog.ReplaceAll("#","")+".png");
  cdatamt->Print(plotNameLog.ReplaceAll("#","")+".eps"); 
  cdatamt->Print(plotNameLog.ReplaceAll("#","")+".pdf"); 
  cdatamt->Print(plotNameLog.ReplaceAll("#","")+".root"); 

  std::cout << "Clean up" << std::endl;
  if (leg) delete leg;
  if (leg2) delete leg2;
  if (hdummy0) delete hdummy0;
  if (hdummy1) delete hdummy1;
  if (hdummy2) delete hdummy2;
  if (hdummy3) delete hdummy3;
  if (frameMt) delete frameMt;
  if (cdatamt) delete cdatamt;
  std::cout << "Fit ends here" << std::endl;
}

void WAnalysis(){

  int cutvalue = 11; //2011 DPD 
  bool doCharge =  false ; //turn on when building or using a charge dependent template
  bool doCentrality = false ;
  bool doEta = false ;  //turn on when building eta templates
  bool correctSpectra = false; //corrects spectra using Aw,Cw factors

  int nSigPdf = 1000000;
  TString nSigPdfString;


  if(doCentrality&&doCharge&&doEta){
      nSigPdfString = "CentChrgEta";
    }
  else if(doCentrality&&doCharge&&!doEta){
      nSigPdfString = "CentChrg";
    }
  else if (doEta&&doCentrality&&!doCharge){
      nSigPdfString = "CentEta";
    }

  else if (doEta&&doCharge&&!doCentrality){
      nSigPdfString = "EtaChrg";
   }

  else if (!doEta&&doCharge&&!doCentrality){
      nSigPdfString = "Charge";
   }
  else{ nSigPdfString = "1000000`";}

  float ptmax = 90.0;
  float mtmax = 200.0;
 
  SetAtlasStyle();

  //TString fileNameDataIn = "HISingleMuonHP.2012.10.04"; 
  TString fileNameDataIn = "HISingleMuonHP.2012.11.16"; 

  //data overlay
  TString fileNameMCIn = "HISingleMuonWmunuMCDataOverlay.2012.11.13";
  //HIJING overlay
  TString fileNameZMCIn = "HISingleMuonMC_PYTHIA_HIJING_Zmumu_11.19.2012";
  //J2 1 muon (yujiao)
  TString fileNameQCDMCIn = "HISingleMuonMC_J21Mu_Pythia_2012.11.19";

  TString fileNameFitOut = "WAnalysis_fitResult"; fileNameFitOut+=nSigPdfString; fileNameFitOut+=".root";

  RooRealVar  muonPt("muonPt","p_{T}",0,100,"GeV");
  RooRealVar  missPt("missPt","p_{T}^{miss}",0,180,"GeV");
  RooRealVar  muonMt("muonMt","m_{T}",0,200,"GeV");
  RooRealVar  muonEta("muonEta","muonEta",-3.0,+3.0);
  RooRealVar  muonELoss("muonELoss","muonELoss",-0.5,+0.5);
  RooRealVar  muonScat("muonScat","muonScat",-4.0,+4.0);
  RooRealVar  centrality("centrality","centrality",0.,1.0);
  RooRealVar  muonCharge("muonCharge","charge",-2.0,+2.0);
  RooRealVar  ptCone20("ptCone20","ptCone20",0.0,+2.0);

  RooFormulaVar extraCut("extraCut", "extraCut", "muonPt>25.0&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&missPt>25.0&&ptCone20/muonPt<0.5&&abs(muonEta)<2.5&&centrality<=0.8", RooArgList(muonPt,muonELoss,muonScat,missPt,ptCone20,muonEta,centrality));

  //anti-W cut to estimate QCD bkg
  RooFormulaVar antiWCut("antiWCut", "antiWCut", "muonPt>25.0&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&missPt<25.0&&ptCone20/muonPt<0.5&&abs(muonEta)<2.5&&centrality<=0.8", RooArgList(muonPt,muonELoss,muonScat,missPt,ptCone20,muonEta,centrality));

  //for fitting templates in eta bins
  RooFormulaVar etaCut1("etaCut1", "etaCut1", "abs(muonEta)>=0.0&&abs(muonEta)<0.25", RooArgList(muonEta,muonEta));
  RooFormulaVar etaCut2("etaCut2", "etaCut2", "abs(muonEta)>=0.25&&abs(muonEta)<0.5", RooArgList(muonEta,muonEta));
  RooFormulaVar etaCut3("etaCut3", "etaCut3", "abs(muonEta)>=0.5&&abs(muonEta)<0.75", RooArgList(muonEta,muonEta));
  RooFormulaVar etaCut4("etaCut4", "etaCut4", "abs(muonEta)>=0.75&&abs(muonEta)<1.0", RooArgList(muonEta,muonEta));
  RooFormulaVar etaCut5("etaCut5", "etaCut5", "abs(muonEta)>=1.0&&abs(muonEta)<1.25", RooArgList(muonEta,muonEta));
  RooFormulaVar etaCut6("etaCut6", "etaCut6", "abs(muonEta)>=1.25&&abs(muonEta)<1.5", RooArgList(muonEta,muonEta));
  RooFormulaVar etaCut7("etaCut7", "etaCut7", "abs(muonEta)>=1.5&&abs(muonEta)<1.75", RooArgList(muonEta,muonEta));
  RooFormulaVar etaCut8("etaCut8", "etaCut8", "abs(muonEta)>=1.75&&abs(muonEta)<2.0", RooArgList(muonEta,muonEta));
  RooFormulaVar etaCut9("etaCut9", "etaCut9", "abs(muonEta)>=2.0&&abs(muonEta)<2.25", RooArgList(muonEta,muonEta));
  RooFormulaVar etaCut10("etaCut10", "etaCut10", "abs(muonEta)>=2.25&&abs(muonEta)<2.5", RooArgList(muonEta,muonEta));

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

  RooCategory bkgCutCategory("bkgCutCategory","bkgCutCategory");
  bkgCutCategory.defineType("ZDY",0);
  
  RooArgSet muonArgSet(muonPt,missPt,muonMt,muonEta,muonCategory,chargeCategory,bkgCutCategory,centrality);
  
  // --- Set pt and eta bins ---
  std::vector<double> ptBins;
  ptBins.push_back(0.0);
  ptBins.push_back(ptmax);
  const int nPtBins = ptBins.size()-1;
  std::vector<double> etaBins;
  etaBins.push_back(0.00);
  if (doEta) {
    etaBins.push_back(+0.25);
    etaBins.push_back(+0.50);
    etaBins.push_back(+0.75);
    etaBins.push_back(+1.00);
    etaBins.push_back(+1.25);
    etaBins.push_back(+1.50);
    etaBins.push_back(+1.75);
    etaBins.push_back(+2.00);
    etaBins.push_back(+2.25);
    
  }
  etaBins.push_back(+2.50);

  const int nEtaBins = etaBins.size()-1;
  std::vector<double> centralityBins;

  centralityBins.push_back(0.00);
  if (doCentrality) {
    centralityBins.push_back(0.05);
    centralityBins.push_back(0.10);
    centralityBins.push_back(0.15);
    centralityBins.push_back(0.20);
    centralityBins.push_back(0.40);
    //centralityBins.push_back(0.60);
    
  }
  centralityBins.push_back(0.80);

  const int nCentralityBins = centralityBins.size()-1;
 
  // --- Fill data sets ---
  RooDataSet* dataSet0 = fillHIMuonDataSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/HardProbesFiles/",fileNameDataIn+".root",muonArgSet, cutvalue); dataSet0->Print();
  RooDataSet* dataSet = (RooDataSet*)dataSet0->reduce(extraCut); dataSet->Print(); 
  dataSet = (RooDataSet*)dataSet->reduce(Cut("bkgCutCategory==bkgCutCategory::ZDY"));dataSet->Print();

  // --- Fill mc sets ---
  RooDataSet* mcSet0 = fillHIMuonDataSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameMCIn+".root", muonArgSet, cutvalue, true); mcSet0->Print();
  RooDataSet* mcSet1 = fillHIMuonDataSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Zmumu/",fileNameZMCIn+".root", muonArgSet, cutvalue, true); mcSet1->Print();
  RooDataSet* mcSet2 = fillHIMuonDataSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/QCD/",fileNameQCDMCIn+".root", muonArgSet, cutvalue, true); mcSet1->Print();

  // --- Fill W set ---
  RooDataSet* mcWSet0 = (RooDataSet*)mcSet0->reduce(Cut("muonCategory==muonCategory::W"),Cut("bkgCutCategory==bkgCutCategory::ZDY"), /*Cut(extraCut),*/ EventRange(0, nSigPdf) ); mcWSet0->Print();
  RooDataSet* mcWSet = (RooDataSet*)mcWSet0->reduce(extraCut); mcWSet->Print();

  // --- Fill Z bkg set ---
  RooDataSet* mcZSet0 = (RooDataSet*)mcSet1->reduce(Cut("muonCategory==muonCategory::Z"),EventRange(0, nSigPdf) ); mcZSet0->Print();
  //weight Z mc by per event yield from paper
  RooFormulaVar wFunc("wZ","Z event weight","4.39e-6",muonMt) ;
  // Add column with variable wZ to Z mc set 
  RooRealVar* wZ = (RooRealVar*) mcZSet0->addColumn(wFunc) ; mcZSet0->Print();
  // Instruct dataset d in interpret w as event weight rather than as observable
  RooDataSet* mcZSet = new RooDataSet(mcZSet0->GetName(), mcZSet0->GetTitle(), mcZSet0, *mcZSet0->get(),0, wZ->GetName()); mcZSet->Print();
  //Apply W candidate cuts
  mcZSet = (RooDataSet*)mcZSet->reduce(extraCut);mcZSet->Print();
  mcZSet = (RooDataSet*)mcZSet->reduce(Cut("bkgCutCategory==bkgCutCategory::ZDY"));mcZSet->Print();

  // --- Invert W cuts ---
  RooDataSet* mcQCDSet = (RooDataSet*)mcSet2->reduce(extraCut); mcQCDSet->Print();
  mcQCDSet = (RooDataSet*)mcQCDSet->reduce(Cut("bkgCutCategory==bkgCutCategory::ZDY"));mcQCDSet->Print();

  delete dataSet0;
  delete mcWSet0;
  delete mcZSet0;
  
  RooPlot* frameMt=muonMt.frame(Range(0,mtmax)) ;
  frameMt->SetXTitle("#font[52]{m}_{T} [GeV]");
  frameMt->SetYTitle("Events/GeV");
  dataSet->plotOn(frameMt, DrawOption("p"), MarkerSize(0.8));
  
 
    // --- Subdivide in bins ---
  RooDataSet* dataSubSet[nPtBins][nEtaBins][nCentralityBins];
  RooDataSet* dataQCDSubSet[nPtBins][nEtaBins][nCentralityBins];
  RooDataSet* mcWSubSet[nPtBins][nEtaBins][nCentralityBins];
  RooDataSet* mcZSubSet[nPtBins][nEtaBins][nCentralityBins];

  if(doEta||doCentrality){
  for ( int i = 0; i < nPtBins; i++ ) {
    for ( int j = 0; j < nEtaBins; j++ ) {
      for ( int k = 0; k < nCentralityBins; k++ ) {
        std::cout << "Divide " << i << ":" <<j << ":" <<k << std::endl;
        dataSubSet[i][j][k] = selectPtEtaCentrality( dataSet, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); dataSubSet[i][j][k]->Print();
        mcWSubSet[i][j][k] = selectPtEtaCentrality( dataSet, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); dataSubSet[i][j][k]->Print();
        mcZSubSet[i][j][k] = selectPtEtaCentrality( dataSet, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); dataSubSet[i][j][k]->Print();
      }
    }
  }
 }

    /// --- Open output file ---
    TDirectory *dir = gDirectory;
    TFile outFile(fileNameFitOut,"RECREATE");
    gDirectory = dir;
    TString resultName = "WPtFit";
    FitResult baselineResult(resultName,ptBins,etaBins,centralityBins);

    // --- fill arrays for Aw,Cw ---
    if(correctSpectra) setCorrectionFactors();

    if(!doCentrality && !doEta && !doCharge){
    /// --- fit inclusive spectra --
        plot( dataSet, mcWSet, mcZSet,mcQCDSet,muonMt, baselineResult, 99,0,0, mtmax, fileNameDataIn, "#mu^{#pm},0-80" , "0 #leq |#eta| < 2.5" );

        if(!doCharge) return;
    }

    // --- charge-inclusive spectra over all eta and centrality --- //
    if (doCharge&&!doCentrality&&!doEta) {

      std::cout << "creating muMinus dataset" << std::endl;
      RooDataSet* dataSetMinus = (RooDataSet*) dataSet->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
      RooDataSet* mcQCDSetMinus = (RooDataSet*) mcQCDSet->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
      RooDataSet* mcWSetMinus = (RooDataSet*)mcWSet->reduce(Cut("chargeCategory==chargeCategory::muMinus") ); mcWSetMinus->Print();
      RooDataSet* mcZSetMinus = (RooDataSet*)mcZSet->reduce(Cut("chargeCategory==chargeCategory::muMinus") ); mcZSetMinus->Print();

      std::cout << "creating muPlus dataset" << std::endl;
      RooDataSet* dataSetPlus  = (RooDataSet*) dataSet->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
      RooDataSet* mcQCDSetPlus  = (RooDataSet*) dataSet->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
      RooDataSet* mcWSetPlus = (RooDataSet*)mcWSet->reduce(Cut("chargeCategory==chargeCategory::muPlus") ); mcWSetPlus->Print();
      RooDataSet* mcZSetPlus = (RooDataSet*)mcZSet->reduce(Cut("chargeCategory==chargeCategory::muPlus") ); mcZSetPlus->Print();

      if (dataSet) delete dataSet;
      if (mcWSet)  delete mcWSet;
      if (mcZSet)  delete mcZSet;

      std::cout << "plotting over all ranges for mu-... " << std::endl;
      plot( dataSetMinus, mcWSetMinus, mcZSetMinus, mcQCDSetMinus, muonMt, baselineResult, 100,0,0, mtmax, fileNameDataIn, "#mu^{-},0-80" , "0 #leq |#eta| < 2.5");
      std::cout << "done." << std::endl;
      std::cout << "fitting over all ranges for mu+... " << std::endl;
      plot( dataSetPlus, mcWSetPlus, mcZSetPlus, mcQCDSetPlus, muonMt, baselineResult, 100,0,0, mtmax, fileNameDataIn, "#mu^{-},0-80" , "0 #leq |#eta| < 2.5");
      std::cout << "done." << std::endl;

      return;
    } //charge inclusive spectra

   // --- fitting dataSubSets binned in pT, eta, and/or centrality --- //
    for ( int i = 0; i < nPtBins; i++ ) {
      for ( int j = 0; j < nEtaBins; j++ ) {
        for ( int k = 0; k < nCentralityBins; k++ ) {

          std::cout << " fitting "<<i<<":"<<j<<":"<<k<<std::endl;
          TString sLow = "";
          TString sHigh = "";
          TString sLow2 = "";
          TString sHigh2 = "";
          
	  sLow += 100*centralityBins[k]; sLow.Remove(3);
          sHigh += 100*centralityBins[k+1]; sHigh.Remove(3);
          
	  sLow2 += etaBins[j];
          sHigh2 += etaBins[j+1];

	  if (doCharge) {

            std::cout << "Creating charged datasets." << std::endl;
            RooDataSet* dataSetPlus  = (RooDataSet*) dataSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
            RooDataSet* dataSetMinus = (RooDataSet*) dataSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
            RooDataSet* mcQCDSetPlus  = (RooDataSet*) dataQCDSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
            RooDataSet* mcQCDSetMinus = (RooDataSet*) dataQCDSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
            RooDataSet* mcWSetPlus  = (RooDataSet*) mcWSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
            RooDataSet* mcWSetMinus = (RooDataSet*) mcWSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
            RooDataSet* mcZSetPlus  = (RooDataSet*) mcWSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
            RooDataSet* mcZSetMinus = (RooDataSet*) mcWSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

            TString sSelPlus = "#mu^{+},";
            TString sSelMinus = "#mu^{-},";
            TString sSelPlus2 = sLow2;
            TString sSelMinus2 = sLow2;
            
	    if (doEta&&doCentrality){

              sSelPlus += sLow; sSelPlus+="-"; sSelPlus+= sHigh;
              sSelMinus += sLow; sSelMinus+="-"; sSelMinus+= sHigh;
              sSelMinus2 += " #leq "; sSelMinus2+= "|#eta| "; sSelMinus2+="< "; sSelMinus2 += sHigh2;
              sSelPlus2 += " #leq "; sSelPlus2+= "|#eta| "; sSelPlus2+="< "; sSelPlus2 += sHigh2;

              std::cout << "Binning for eta+charge+centrality." << std::endl;
              std::cout << "plotting for mu+: " << i << ":" << j << ":" << k <<std::endl;

              plot( dataSetPlus, mcWSetPlus, mcZSetPlus,mcQCDSetPlus, muonMt, baselineResult, 102+i,j,k, mtmax, fileNameDataIn, sSelPlus, sSelPlus2);
              
	      std::cout << "Done!" << std::endl;
              std::cout << "plotting for mu-: " << i << ":" << j << ":" << k <<std::endl;

              plot( dataSetMinus, mcWSetMinus, mcZSetMinus,mcQCDSetMinus, muonMt, baselineResult, 103+i,j,k, mtmax, fileNameDataIn, sSelMinus, sSelMinus2);

              std::cout << "Done!" << std::endl;
            }

            if (doEta&&!doCentrality) {

                sSelPlus += sLow; sSelPlus+="-"; sSelPlus+= sHigh;
                sSelMinus += sLow; sSelMinus+="-"; sSelMinus+= sHigh;
                TString sSelMinus2 = sLow2; sSelMinus2 += " #leq "; sSelMinus2+= "|#eta| "; sSelMinus2+="< "; sSelMinus2 += sHigh2;
                TString sSelPlus2 = sLow2; sSelPlus2 += " #leq "; sSelPlus2+= "|#eta| "; sSelPlus2+="< "; sSelPlus2 += sHigh2;

                std::cout << "plottting for mu+: " << i << ":" << j << ":" << k <<std::endl;
                plot( dataSetPlus, mcWSetPlus, mcZSetPlus,mcQCDSetPlus, muonMt, baselineResult, 102+i,j,k, mtmax, fileNameDataIn, sSelPlus, sSelPlus2);
                std::cout << "Done!" << std::endl;

                std::cout << "plotting for mu-: " << i << ":" << j << ":" << k <<std::endl;
                plot( dataSetMinus, mcWSetMinus, mcZSetMinus,mcQCDSetMinus, muonMt, baselineResult, 103+i,j,k, mtmax, fileNameDataIn, sSelMinus, sSelMinus2);
                std::cout << "Done!" << std::endl;

              }  //eta

             if (doCentrality&&!doEta){

                  sSelPlus += sLow; sSelPlus+="-"; sSelPlus+= sHigh;
                  sSelMinus += sLow; sSelMinus+="-"; sSelMinus+= sHigh;
                  sSelMinus2 += " #leq "; sSelMinus2+= "|#eta| "; sSelMinus2+="< "; sSelMinus2 += sHigh2;
                  sSelPlus2 += " #leq "; sSelPlus2+= "|#eta| "; sSelPlus2+="< "; sSelPlus2 += sHigh2;

                  std::cout << "plotting for mu+: " << i << ":" << j << ":" << k <<std::endl;
                  plot( dataSetPlus, mcWSetPlus, mcZSetPlus,mcQCDSetPlus, muonMt, baselineResult, 102+i,j,k, mtmax, fileNameDataIn, sSelPlus, sSelPlus2);
                  std::cout << "Done!" << std::endl;

                  std::cout << "plotting for mu-: " << i << ":" << j << ":" << k <<std::endl;
                  plot( dataSetMinus, mcWSetMinus, mcZSetMinus,mcQCDSetMinus, muonMt, baselineResult, 103+i,j,k, mtmax, fileNameDataIn, sSelMinus, sSelMinus2);
                  std::cout << "Done!" << std::endl;

              }  //centrality
            } //charge

             if (!doCharge&&doCentrality&&doEta) {

                std::cout << "Plotting for centrality and eta classes. Inclusive for charges." << std::endl;
                TString sSel = "#mu^{#pm},"; sSel += sLow; sSel+="-"; sSel+= sHigh;
                TString sSel2 = sLow2; sSel2 += " #leq "; sSel2+= "|#eta| "; sSel2+="< "; sSel2 += sHigh2;

                plot(dataSubSet[i][j][k], mcWSubSet[i][j][k], mcZSubSet[i][j][k],dataQCDSubSet[i][j][k], muonMt, baselineResult, 102+i,j,k, mtmax, fileNameDataIn, sSel, sSel2 );

             }  //eta and centrality
            if ( dataSubSet[i][j][k] ) delete dataSubSet[i][j][k];
        }
      }
    }

    outFile.cd();
    baselineResult.write(outFile);
    outFile.Close();
}

int main() {
  WAnalysis();
}
/// eof













































