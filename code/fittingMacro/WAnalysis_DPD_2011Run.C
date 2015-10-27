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
// #include "Systematics.C"
#include "../EfficiencyCorrection.C"
#include "flavourFitterHIDep_DPD_2011Run.C"
#include "../AtlasUtils.C"
#include "../AtlasStyle.C"

#include <iostream>
#include <iomanip>
#include <fstream>

////////////////////////////////////////////////////////////////////////////////
// fixplot
////////////////////////////////////////////////////////////////////////////////
void fixplot(RooPlot* fr)
{
 RooHist* hist0 = (RooHist*) fr->getObject(0) ;
//  RooHist* hist1 = (RooHist*) fr->getObject(1) ;

 Double_t x,y ;
 for (Int_t i=0 ; i<hist0->GetN() ; i++) {
   hist0->GetPoint(i,x,y) ;
   if (y==0) {
     hist0->SetPointEYhigh(i,0) ;
   }
 }

//  for (Int_t i=0 ; i<hist1->GetN() ; i++) {
//    hist1->GetPoint(i,x,y) ;
//    if (y==0) {
//      hist1->SetPointEYhigh(i,0) ;
//    }
//  }
}

///////////////////////////////////////////////////////////////////////////////
// buildWTemplate
///////////////////////////////////////////////////////////////////////////////
void buildWTemplate( RooDataSet* mcSet, RooRealVar& var, const char* fileName )
{
  std::cout<<"buildWTemplates"<<std::endl;

  RooRealVar nSignal("nSignal","nSignal", mcSet->numEntries() );

  
  // --- Create a workspace ---
  RooWorkspace *w = new RooWorkspace("w","workspace");
  w->import(var);
  
  // --- Build kernel-estimated templates ---
  std::cout << "Building the RooKeysPdf's... This will take some time" << std::endl;
  RooKeysPdf sigPdf("sigPdf","sigPdf",var,*mcSet);

  // --- And write to file ---
  std::cout << "Writing the templates to file " << fileName << std::endl;
  w->import(sigPdf);
  w->import(nSignal);
  w->writeToFile(fileName);
  std::cout << "Done!" << std::endl;
  
  // --- Clean up ---
  delete w;
}

///////////////////////////////////////////////////////////////////////////////
// fit
///////////////////////////////////////////////////////////////////////////////
void fit( RooDataSet* dataSet, RooRealVar& muonPt, 
          FitResult& fitResult, const int iPt, const int iEta, const int iCentrality, const RooAbsPdf& sigPdf, const float ptmax, TString fileNameDataIn, TString sSel, TString sSel2 , bool plotLL = false , bool plotPLL = false , bool addPercent = true)
{
  float ptcut = 30.;
  float ptcutLow = 7.; 
  float ptmaxLow = 20; /// control region
  float chi2 = 0;
  float bkgUncertainty = 0.35172; /// relative uncertainty from 0< cent <1 fit with 30,7,20 GeV.
  bool doCutandCount = false;
  RooPlot* framePt=muonPt.frame(Range(0,ptmax)) ;
  framePt->SetXTitle("#font[52]{p}_{T}^{#mu} [GeV]");
  framePt->SetYTitle("Muons/GeV");
  RooBinning b = RooBinning(ptmax,0,ptmax); // 1 GeV per bin
  dataSet->plotOn(framePt, RefreshNorm(), Binning(b), DrawOption("p"), MarkerSize(0.8));
  
  RooRealVar  exp1K("exp1K","exp1K", -7.27130e-01, -1.0, -0.5); 

  RooExponential exp1Pdf("exp1Pdf","exp1Pdf",muonPt,exp1K);

  RooRealVar  bK("bK","bK",-1.23273e+00, -1.4, -1.0);
  bK.setConstant(kTRUE) ;  //lock

  RooGenericPdf bPdf("bPdf","bPdf","TMath::Exp(bK*sqrt(muonPt))/(muonPt*muonPt*sqrt(muonPt)+0.1)",RooArgSet(bK,muonPt));
  //RooGenericPdf bPdf("bPdf","bPdf","muonPt*TMath::Exp(-1*muonPt*muonPt/(2*sigma*sigma))/sigma*sigma",RooArgSet(muonPt,sigma));
  RooRealVar exp1Ratio("exp1Ratio","exp1Ratio",2.26686e-03, 0.0, 1.0);
  
  ///comment out for power law only fit
  RooAddPdf  bkgPdf("bkgPdf","bkgPdf", RooArgList(exp1Pdf,bPdf), RooArgList(exp1Ratio));
  RooArgList modelComponents;
  modelComponents.add(sigPdf);
  modelComponents.add(bkgPdf);
//  modelComponents.add(bPdf);
  RooRealVar nsig("nsig","Number of W Events",1.0e4,0.0,1e9);
  RooRealVar nbkg("nbkg","nbkg",7.8e7,1.0e2,1.0e11);

  RooAddPdf  mdlPdf("mdlPdf","mdlPdf",modelComponents,RooArgList(nsig,nbkg));
  std::cout << "Creating NLL..." << std::endl;
  //RooAbsReal* nll = mdlPdf.createNLL(*dataSet,Range(ptcutLow,ptmax), NumCPU(4)) ;
  RooAbsReal* nll = mdlPdf.createNLL(*dataSet,Range(ptcutLow,ptmax)) ;
  RooMinuit m1(*nll) ;
  m1.migrad();
  ///release bK
  bK.setConstant(kFALSE) ;
  ///lock exp1K
  exp1K.setConstant(kTRUE) ;
  m1.migrad();
  ///release exp1K
  exp1K.setConstant(kFALSE) ;
  ///lock exp1Ratio
  m1.migrad();
  exp1Ratio.setConstant(kFALSE) ;
  ///re-run
  m1.migrad() ; 
  
  ///lock bK
/*  bK.setConstant(kTRUE) ;
  m1.migrad() ;
  bK.setConstant(kFALSE) ;
  exp1K.setConstant(kTRUE) ;
  m1.migrad();
  exp1K.setConstant(kFALSE) ;
  exp1Ratio.setConstant(kTRUE) ;
  m1.migrad();
  exp1Ratio.setConstant(kFALSE) ;
  bK.setConstant(kTRUE) ;
  m1.migrad() ;
  bK.setConstant(kFALSE) ;
  
*/  
//  m1.migrad() ;
//  exp1Ratio.setConstant(kTRUE) ;
//  m1.migrad();
//  exp1Ratio.setConstant(kFALSE) ;
//  m1.migrad() ;

  cout << "now running minos: " << endl;
  //m1.minos(RooArgSet(bK,exp1K,exp1Ratio)) ;
  m1.minos(RooArgSet(nsig)) ;
  RooFitResult* roofitResult = m1.save();
  //roofitResult->Print();

  //mdlPdf.plotOn(framePt,LineColor(kGreen-8), FillColor(kGreen-8), DrawOption("F"), Range(0,ptmax), Components(sigPdf), Normalization(1.0,RooAbsReal::RelativeExpected));
  mdlPdf.plotOn(framePt,LineColor(kAzure-9), FillColor(kAzure-9), DrawOption("F"), Range(0,ptmax), Components(sigPdf), Normalization(1.0,RooAbsReal::RelativeExpected));
  cout << "now printing framePt: " << endl;
  framePt->Print();
  RooHist* hmdlPdfSig = (RooHist*)framePt->findObject("mdlPdf_Norm[muonPt]_Comp[sigPdf]");
  hmdlPdfSig->SetPoint(0,0,0); 
  hmdlPdfSig->SetPoint(hmdlPdfSig->GetN()-1,ptmax,0.);
  
  mdlPdf.plotOn(framePt,LineColor(kGreen+1),LineStyle(kDashed),LineWidth(2),Range(0,ptmax),Components(bkgPdf),Normalization(1.0,RooAbsReal::RelativeExpected));
//  mdlPdf.plotOn(framePt,LineColor(kGreen+1),LineStyle(kDashed),LineWidth(2),Range(0,ptmax),Components(bPdf),Normalization(1.0,RooAbsReal::RelativeExpected));
  
  std::cout << "now plotting data line" << std::endl;
  mdlPdf.plotOn(framePt,LineColor(kBlack),Range(ptcutLow,ptmax),Normalization(1.0,RooAbsReal::RelativeExpected));
//   TString chi2 = "chi2/ndf = "; chi2+=0.01*floor(100*framePt->chiSquare(5));
  double eventsInFit = dataSet->sumEntries(TString(TString("muonPt>")+=ptcutLow));
  cout << "LL = " << nll->getVal() << " dof = " <<  eventsInFit-2 << " p = " << TMath::Prob(-2*nll->getVal(),eventsInFit-2) <<endl;

  ///calculate the residual for fit wrt data
  /*std::cout << "Plotting residual of likelihood fit wrt data..." << std::endl;
  RooHist* hresid = framePt->residHist("h_"+fileNameDataIn+".root","mdlPdf_Norm[muonPt]") ;
  RooPlot* frameRes = muonPt.frame(Title("Residual Distribution") ) ;
  frameRes->addPlotable(hresid,"P") ;
  frameRes->SetXTitle("#font[52]{p}_{T}^{#mu} [GeV]"); frameRes->SetYTitle("Residual");
  frameRes->SetMarkerSize(0.6);
  std::cout << "Done." << std::endl;*/

  ///calculate pull for fit wrt data
  std::cout << "Plotting pull of likelihood fit wrt data..." << std::endl;
  RooHist* hpull = framePt->pullHist("h_"+fileNameDataIn+".root","mdlPdf_Norm[muonPt]") ;
  RooPlot* framePull = muonPt.frame(Title("Pull Distribution") ) ;
  framePull->addPlotable(hpull,"P") ;
  framePull->SetXTitle("#font[52]{p}_{T}^{#mu} [GeV]"); framePull->SetYTitle("#frac{(y_{i}-#lambda)}{#sigma}");
  framePull->SetMarkerSize(0.6);
  std::cout << "Done." << std::endl;


  ///calculate the residual for signal template wrt data
  std::cout << "Plotting residual of signal template wrt data..." << std::endl;
  RooHist* hresidSig = framePt->residHist("h_"+fileNameDataIn+".root","mdlPdf_Norm[muonPt]_Comp[sigPdf]",true) ;
  RooPlot* frameResSig = muonPt.frame(Title("#frac{(y_{i}-#lambda)}{#sigma}") ) ;
  frameResSig->addPlotable(hresidSig,"P") ;
  frameResSig->SetXTitle("#font[52]{p}_{T}^{#mu} [GeV]"); frameResSig->SetYTitle("#frac{(y_{sig}-#lambda)}{#sigma}");
  frameResSig->SetMarkerSize(0.6);
  std::cout << "Done." << std::endl;

  std::cout << "Plotting residual of background wrt data..." << std::endl;
  RooHist* hresidBkg = framePt->residHist("h_"+fileNameDataIn+".root","mdlPdf_Norm[muonPt]_Comp[bkgPdf]",true) ;
//  RooHist* hresidBkg = framePt->residHist("h_"+fileNameDataIn+".root","mdlPdf_Norm[muonPt]_Comp[bPdf]",true) ;
  RooPlot* frameResBkg = muonPt.frame(Title("#frac{(y_{i}-#lambda)}{#sigma}") ) ;
  frameResBkg->addPlotable(hresidBkg,"P") ;
  frameResBkg->SetXTitle("#font[52]{p}_{T}^{#mu} [GeV]"); frameResBkg->SetYTitle("#frac{(y_{bkg}-#lambda)}{#sigma}");
  frameResBkg->SetMarkerSize(0.6);
  std::cout << "Done." << std::endl;

//  chi2 = framePt->chiSquare("mdlPdf_Norm[muonPt]","h_"+fileNameDataIn+".root",2); 
  //chi2 = framePt->chiSquare(2) ;
  //cout << "chi2/dof = " << framePt->chiSquarePearson(2) << " chi2 = " << framePt->chiSquarePearson(2)*framePt->ndfPearson(2) << " dof = " << framePt->ndfPearson(2) << " p = " << TMath::Prob( framePt->chiSquarePearson(2)*framePt->ndfPearson(2),  framePt->ndfPearson(2)) << endl;

  cout << "nsig = " << nsig.getVal() << " +" << nsig.getAsymErrorHi() << " " << nsig.getAsymErrorLo() << endl;
//  cout << "bK = " << bK.getVal() << " +" << bK.getErrorHi() << " " << bK.getErrorLo() << endl;
  fitResult.setSig(iPt, iEta, iCentrality, nsig.getVal(), nsig.getAsymErrorLo()>0?nsig.getAsymErrorLo():nsig.getErrorLo(), nsig.getAsymErrorHi()>0?nsig.getAsymErrorHi():nsig.getErrorHi());

//  fitResult.setexp1K(iPt, iEta, iCentrality, exp1K.getVal(), exp1K.getError(), exp1K.getError());
//  fitResult.setbK(iPt, iEta, iCentrality, bK.getVal(), bK.getError(), bK.getError());
//  fitResult.setexp1Ratio(iPt, iEta, iCentrality, exp1Ratio.getVal(), exp1Ratio.getError(), exp1Ratio.getError());
//  fitResult.setnbkg(iPt, iEta, iCentrality, nbkg.getVal(), nbkg.getErrorLo()>0?nbkg.getErrorLo():nbkg.getErrorLo(), nbkg.getErrorHi()>0?nbkg.getErrorHi():nbkg.getErrorHi());

  fitResult.setChi2(iPt, iEta, iCentrality, chi2, 0, 0);
  
//   RooAddPdf bkgUpPdf = RooAddPdf(bkgPdf, "bkgUpPdf");
//   RooAddPdf bkgDnPdf = RooAddPdf(bkgPdf, "bkgDnPdf");
//   mdlPdf.plotOn(framePt,LineStyle(kDashed),Range(0,ptmax),Components(bkgPdf),Normalization(1.0+nbkg.getErrorHi()/nbkg.getVal(),RooAbsReal::RelativeExpected));
//   mdlPdf.plotOn(framePt,LineStyle(kDashed),Range(0,ptmax),Components(bkgPdf),Normalization(1.0+nbkg.getErrorLo()/nbkg.getVal(),RooAbsReal::RelativeExpected));
//   bkgPdf.Print();
//   bkgUpPdf.Print();
  
//   cout  << " pt>15: " << dataSet->sumEntries("muonPt>15") 
//         << " pt>20: " << dataSet->sumEntries("muonPt>20") 
//         << " pt>25: " << dataSet->sumEntries("muonPt>25") 
//         << " pt>30: " << dataSet->sumEntries("muonPt>30") << endl;
//   int ptcutBin = hmdlPdfSig->GetXaxis()->FindBin(ptcut);
//   int ptcutBin = (int)floor( ptcut/hmdlPdfSig->GetXaxis()->GetBinWidth(0) ); // 1 GeV bins
//   int ptendBin = hmdlPdfSig->GetN();
//   cout << hmdlPdfSig->GetN() << " " << hmdlPdfSig->GetX()[0] << " to " << hmdlPdfSig->GetX()[hmdlPdfSig->GetN()-1] << endl;
  /// simple root integration since analytic did not work for RooKeysPdf?
//   double sigEffPass = hmdlPdfSig->Integral(ptcutBin, ptendBin)/hmdlPdfSig->Integral(-1,ptendBin);
//   cout  << ptcutBin << ", signal pass " << hmdlPdfSig->Integral(ptcutBin,ptendBin) << " of " << hmdlPdfSig->Integral(-1,ptendBin) <<" gives efficiency " << sigEffPass << endl;
  
//   muonPt.setRange("allRange",ptcutLow,ptmax) ;
//   muonPt.setRange("signalRange",ptcut,ptmax) ;
//   RooAbsReal::defaultIntegratorConfig()->Print("v") ;
//   RooMsgService::instance().addStream(DEBUG,Topic(Integration)) ;
//   RooAbsReal::defaultIntegratorConfig()->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
//   RooAbsReal* int_bkg = bkgPdf.createIntegral(muonPt,Range("signalRange")) ;
//   RooAbsReal* int_bkgAll = bkgPdf.createIntegral(muonPt,Range("allRange")) ;
//   RooAbsReal* int_sig = sigPdf.createIntegral(muonPt,Range("signalRange")) ;
//   RooAbsReal* int_sigAll = sigPdf.createIntegral(muonPt,Range("allRange")) ;
//   double fBkgCut = int_bkg->getVal();
//   double fBkgAll = int_bkgAll->getVal();
//   double fSigCut = int_sig->getVal();
//   double fSigAll = int_sigAll->getVal();
//   double sigEffPass = fSigCut/fSigAll;
//   for (int i = 0; i < 6; ++i){
//     ptcut = 15.+i*5.;
  if(doCutandCount){
  /// estimate by integrating bkgPdf
  TF1* tf1BkgPdf = bkgPdf.asTF(RooArgList(muonPt));
  TF1* tf1SigPdf = sigPdf.asTF(RooArgList(muonPt));
  //TF1* tf1BkgPdf = bPdf.asTF(RooArgList(muonPt));
//   tf1BkgPdf->Print();
//   cout << " par 0 = " << tf1BkgPdf->GetParameter(0) << endl;
//   TF1 tf1BkgPdfUp = TF1(*bkgPdf);
//   TF1 tf1BkgPdfDn = TF1(*bkgPdf);
  
// calculate background intgral
  double bkgInt0 = tf1BkgPdf->Integral(ptcutLow,ptmax) ; //7-90GeV
  double bkgInt1 = tf1BkgPdf->Integral(15.,ptmax) ; 
  double bkgInt2 = tf1BkgPdf->Integral(20.,ptmax) ; 
  double bkgInt3 = tf1BkgPdf->Integral(25.,ptmax) ; 
  double bkgInt4 = tf1BkgPdf->Integral(30.,ptmax) ; 
  double bkgInt5 = tf1BkgPdf->Integral(35.,ptmax) ; 
  double bkgInt6 = tf1BkgPdf->Integral(40.,ptmax) ; 
  cout << " background integral 0: " << bkgInt0 << endl;
  cout << " background integral 1: " << bkgInt1 << endl;
  cout << " background integral 2: " << bkgInt2 << endl;
  cout << " background integral 3: " << bkgInt3 << endl;
  cout << " background integral 4: " << bkgInt4 << endl;
  cout << " background integral 5: " << bkgInt5 << endl;
  cout << " background integral 6: " << bkgInt6 << endl;

// calculate signal integral
  double sigInt0 = tf1SigPdf->Integral(ptcutLow,ptmax) ; //7-90GeV
  double sigInt1 = tf1SigPdf->Integral(15.,ptmax) ; 
  double sigInt2 = tf1SigPdf->Integral(20.,ptmax) ; 
  double sigInt3 = tf1SigPdf->Integral(25.,ptmax) ; 
  double sigInt4 = tf1SigPdf->Integral(30.,ptmax) ; 
  double sigInt5 = tf1SigPdf->Integral(35.,ptmax) ; 
  double sigInt6 = tf1SigPdf->Integral(40.,ptmax) ; 
  cout << " signal integral 0: " << sigInt0 << endl;
  cout << " signal integral 1: " << sigInt1 << endl;
  cout << " signal integral 2: " << sigInt2 << endl;
  cout << " signal integral 3: " << sigInt3 << endl;
  cout << " signal integral 4: " << sigInt4 << endl;
  cout << " signal integral 5: " << sigInt5 << endl;
  cout << " signal integral 6: " << sigInt6 << endl;



  double fBkgCut = tf1BkgPdf->Integral(ptcut,ptmax);
  double fBkgCtrl = tf1BkgPdf->Integral(ptcutLow,ptmaxLow);
  double fBkgAll = tf1BkgPdf->Integral(ptcutLow,ptmax);
  TString sCtrl = "(muonPt>"; sCtrl+=ptcutLow; sCtrl+="&&muonPt<"; sCtrl+=ptmaxLow; sCtrl+=")";
  float nCtrl = dataSet->sumEntries(sCtrl);
  cout  << " integrate pt>" << ptcut << ": " 
        << " background: "      << setprecision(15) << fBkgCut 
        << " control: " << setprecision(15) << fBkgCtrl 
        << " allBkg: " << setprecision(15) << fBkgAll 
        << " nCtrl: " << nCtrl
        << " tf1: " << tf1BkgPdf->Integral(ptcut,ptmax)
        << " all tf1: " << tf1BkgPdf->Integral(ptcutLow,ptmax)
//         << " signal "  << setprecision(15) << fSigCut 
//         << " allSig: " << setprecision(15) << fSigAll 
        << endl;

  /// read in preestimate
  double bkgEffPass = efficiencyWBkg(ptcut,bkgInt0,bkgInt1,bkgInt2,bkgInt3,bkgInt4,bkgInt5,bkgInt6); // pass into signal region
  //double bkgEffPass = efficiencyWBkg(ptcut); // pass into signal region
  double bkgEffCtrl = 1.-efficiencyWBkg(ptmaxLow, bkgInt0,bkgInt1,bkgInt2,bkgInt3,bkgInt4,bkgInt5,bkgInt6); // control region only, efficiency is normalized to [ptcutLow,inf] (make sure to update!)
  double bkgPass = bkgEffPass/bkgEffCtrl*nCtrl;
  //double bkgPass = bkgEffPass*dataSet->sumEntries(TString(TString("muonPt>")+=ptcutLow)); /// assume that at low pt nearly all events are background
  //bkgPass *= (1-591./37909); /// small correction due to W contamination in bkg estimate, derived from data with ptcutLow=7 GeV
  double eventPass = dataSet->sumEntries(TString(TString("muonPt>")+=ptcut));
  double sigPass = eventPass - bkgPass;
  // double estSig  = sigPass/sigEffPass; /// account for W cut out
  double sigEffPass = efficiencyW(ptcut,sigInt0,sigInt1,sigInt2,sigInt3,sigInt4,sigInt5,sigInt6); 
  double estSig = sigPass/sigEffPass;
  cout  << " expected background at pt>" << ptcut << ": " << bkgPass << endl; 
  cout  << " that leaves " << sigPass << " W candidates of total " << eventPass << endl;
  cout  << " accounting for pt cut efficiency nW = " << estSig << endl; 
//   }
  fitResult.setSigN(iPt, iEta, iCentrality, estSig, sqrt(estSig), sqrt(estSig));        /// cut&count with background subtraction
  fitResult.setMcN(iPt, iEta, iCentrality, bkgPass, bkgPass*bkgUncertainty, bkgPass*bkgUncertainty);      /// estimated background with error
  fitResult.setMc(iPt, iEta, iCentrality, eventPass, sqrt(eventPass), sqrt(eventPass)); /// using unused placeholder to store raw counts
  } 

  TLegend* leg = new TLegend(0.58, 0.45, 0.85, 0.73);
  TH1F* hdummy0 = new TH1F(); hdummy0->SetMarkerSize(0.8);
  TH1F* hdummy1 = new TH1F(); hdummy1->SetFillColor(kAzure-9); 
  TH1F* hdummy2 = new TH1F(); hdummy2->SetFillColor(kGray); hdummy2->SetLineColor(kAzure-9); hdummy2->SetLineStyle(kDashed); hdummy2->SetLineWidth(2);
  TH1F* hdummy3 = new TH1F(); hdummy3->SetLineColor(kBlack); hdummy3->SetLineWidth(2);
  leg->SetTextFont(gStyle->GetTextFont());
  leg->SetTextSize(gStyle->GetTextSize());
  leg->SetBorderSize(0);
  leg->SetFillColor(0);

//   leg->AddEntry(hdummy0, "Simulation b#bar{b}", "pe");
  leg->AddEntry(hdummy0, "Data 2011", "pe");
  leg->AddEntry(hdummy1, "W#rightarrow#mu#nu", "f");
  //leg->AddEntry(hdummy2, "Background", "f");
  leg->AddEntry(hdummy2, "Background", "l");
  leg->AddEntry(hdummy3, "Fit", "l");
    
  TCanvas* cdatapt = new TCanvas("cdatapt","datapt",600,600);
  cdatapt->UseCurrentStyle();
  cdatapt->SetLogy(true);
  framePt->SetMinimum(0.1);
//  framePt->SetMaximum(1e5);
  framePt->SetMaximum(1e8);
  std::cout << "now printing framePt" << std::endl;
  framePt->Print();
  std::cout << "now fixing plot" << std::endl;
  fixplot(framePt);
//   framePt->drawAfter("mdlPdf_Norm[muonPt]_Comp[sigPdf]","h_"+fileNameDataIn+".root");
  std::cout << "now creating root file" << std::endl;
  framePt->drawAfter("mdlPdf_Norm[muonPt]","h_"+fileNameDataIn+".root");
//   framePt->Print();
  std::cout << "now drawing framePt" << std::endl;
  framePt->Draw();
  std::cout << "now drawing legend" << std::endl;
  leg->Draw();
  TLatex l;
  l.SetNDC();
  l.DrawLatex(0.28,0.75,sSel + ( addPercent ? "%" : "" ));
  l.DrawLatex(0.28,0.65,sSel2);
  l.DrawLatex(0.6,0.75,"#sqrt{s_{NN}}=2.76 TeV");
  l.DrawLatex(0.57,0.85,"#int Ldt #approx 0.140 nb^{-1}"); 
  //l.DrawLatex(0.6,0.85,"#int L #approx 5 #mub^{-1}"); 
//#ifdef __CINT__
  ATLAS_LABEL(0.17,0.85,1);
  //myText(0.35,0.85, (Color_t)kBlack, (char*)("Preliminary"));
  myText(0.33,0.85, (Color_t)kBlack, (char*)("Internal"));
  //myText(0.35,0.85, (Color_t)kBlack, (char*)("Performance"));
  //myText(0.25,0.70, (Color_t)kRed, (char*)("For Approval"));
  //myText(0.25,0.70, (Color_t)kRed, (char*)("With Approved Parameters"));
//#endif
  //TString plotNameLog = "dataPt_"; plotNameLog+=fileNameDataIn; plotNameLog+=sSel; plotNameLog+=","; plotNameLog+=sSel2; plotNameLog+="Log"; //.png";
  std::cout << "now making file names" << std::endl;
  TString plotNameLog = "dataPt_"; plotNameLog+=sSel; plotNameLog+=","; plotNameLog+=sSel2; plotNameLog+=","; plotNameLog+="Log"; //.png";
  //TString plotNameLog = "dataPt_"; plotNameLog+=sSel; plotNameLog+="_"; plotNameLog+=sSel2; plotNameLog+="_"; plotNameLog+="Log"; //.png";
  std::cout << "now saving files" << std::endl;

/*  cdatapt->Print(plotNameLog.ReplaceAll("#","")+".png");
  cdatapt->Print(plotNameLog.ReplaceAll("#","")+".eps"); 
  cdatapt->Print(plotNameLog.ReplaceAll("#","")+".pdf"); 
  cdatapt->Print(plotNameLog.ReplaceAll("#","")+".root"); 
  //cdatapt->Print("muPt_StdCuts.root"); 
  */
/*  cdatapt->Print(plotNameLog+"_GC.png");
  cdatapt->Print(plotNameLog+"_GC.eps"); 
  cdatapt->Print(plotNameLog+"_GC.pdf"); 
  cdatapt->Print(plotNameLog+"_GC.root"); 
*/
  TCanvas* cPull = new TCanvas("cPull","cPull",600,600) ;
  cPull->cd(); gPad->SetLeftMargin(0.15); framePull->GetYaxis()->SetTitleOffset(1.6); 
  
  //frameRes->SetMinimum(-80.0);
  //frameRes->SetMaximum(80.0);
  framePull->Draw() ;
  TLatex l2;
  l2.SetNDC();
  l2.DrawLatex(0.3,0.75,sSel + ( addPercent ? "%" : "" ));
  l2.DrawLatex(0.3,0.65,sSel2);
  l2.DrawLatex(0.6,0.75,"#sqrt{s_{NN}}=2.76 TeV");
  l2.DrawLatex(0.6,0.85,"#int L dt #approx 0.140 #mub^{-1}"); 
  ATLAS_LABEL(0.17,0.85,1);
  myText(0.33,0.85, (Color_t)kBlack, (char*)("Internal"));
  plotNameLog+="_pull";
/*  cPull->Print(plotNameLog+".root");
  cPull->Print(plotNameLog.ReplaceAll("#","")+".png");
  cPull->Print(plotNameLog.ReplaceAll("#","")+".eps");
  cPull->Print(plotNameLog.ReplaceAll("#","")+".pdf");
*/
  //plot residual of data wrt nll fit
/*  TCanvas* cRes = new TCanvas("cRes","cRes",600,600) ;
  cRes->cd(); gPad->SetLeftMargin(0.15); frameRes->GetYaxis()->SetTitleOffset(1.6); 
  
  frameRes->SetMinimum(-80.0);
  frameRes->SetMaximum(80.0);
  frameRes->Draw() ;
  TLatex l2;
  l2.SetNDC();
  l2.DrawLatex(0.3,0.75,sSel + ( addPercent ? "%" : "" ));
  l2.DrawLatex(0.3,0.65,sSel2);
  l2.DrawLatex(0.6,0.75,"#sqrt{s_{NN}}=2.76 TeV");
  l2.DrawLatex(0.6,0.85,"#int L dt #approx 0.140 #mub^{-1}"); 
  ATLAS_LABEL(0.17,0.85,1);
  myText(0.33,0.85, (Color_t)kBlack, (char*)("Internal"));
  plotNameLog+="_residual";
*/  
  //cRes->Print(plotNameLog+".root");
/*  cRes->Print(plotNameLog,"root");
  cRes->Print(plotNameLog.ReplaceAll("#","")+".png");
  cRes->Print(plotNameLog.ReplaceAll("#","")+".eps");
  cRes->Print(plotNameLog.ReplaceAll("#","")+".pdf");
*/
  //plot residual of data wrt signal template 
  TCanvas* cResSig = new TCanvas("cResSig","cResSig",600,600) ;
  cResSig->cd(); gPad->SetLeftMargin(0.15); frameResSig->GetYaxis()->SetTitleOffset(1.6); 
  frameResSig->SetMinimum(-80.0);
  frameResSig->SetMaximum(80.0);
  frameResSig->Draw() ;
  TLatex l3;
  l3.SetNDC();
  l3.DrawLatex(0.3,0.75,sSel + ( addPercent ? "%" : "" ));
  l3.DrawLatex(0.3,0.65,sSel2);
  l3.DrawLatex(0.6,0.75,"#sqrt{s_{NN}}=2.76 TeV");
  l3.DrawLatex(0.6,0.85,"#int L #approx 0.140 #mub^{-1}"); 
  ATLAS_LABEL(0.17,0.85,1);
  myText(0.33,0.85, (Color_t)kBlack, (char*)("Internal"));
 
  //plotNameLog+="_Signal";
  //cResSig->Print(plotNameLog+"_Signal.root");
  /*cResSig->Print(plotNameLog,"root");
  cResSig->Print(plotNameLog.ReplaceAll("#","")+"_Signal.png");
  cResSig->Print(plotNameLog.ReplaceAll("#","")+"_Signal.eps");
  cResSig->Print(plotNameLog.ReplaceAll("#","")+"_Signal.pdf");
 */ 
  //plot residual of data wrt background function 
  TCanvas* cResBkg = new TCanvas("cResBkg","cResBkg",600,600) ;
  cResBkg->cd(); gPad->SetLeftMargin(0.15); frameResBkg->GetYaxis()->SetTitleOffset(1.6); 
  frameResBkg->SetMinimum(-80.0);
  frameResBkg->SetMaximum(80.0);
  frameResBkg->Draw() ;
  TLatex l4;
  l4.SetNDC();
  l4.DrawLatex(0.3,0.75,sSel + ( addPercent ? "%" : "" ));
  l4.DrawLatex(0.3,0.65,sSel2);
  l4.DrawLatex(0.6,0.75,"#sqrt{s_{NN}}=2.76 TeV");
  l4.DrawLatex(0.6,0.85,"#int L #approx 0.140 #mub^{-1}"); 
  ATLAS_LABEL(0.17,0.85,1);
  myText(0.33,0.85, (Color_t)kBlack, (char*)("Internal"));
 
/*  cResBkg->Print(plotNameLog+"_Bkg.root");
  //cResBkg->Print(plotNameLog,"root");
  cResBkg->Print(plotNameLog.ReplaceAll("#","")+"_Bkg.png");
  cResBkg->Print(plotNameLog.ReplaceAll("#","")+"_Bkg.eps");
  cResBkg->Print(plotNameLog.ReplaceAll("#","")+"_Bkg.pdf");
 */

  //cdatapt->Print(plotNameLog+".pdf"); 
  //cdatapt->Print(plotNameLog+".pdf"); 
//    cdatapt.SetLogy(false);
//  framePt->SetMinimum(0.0);
//  framePt->SetMaximum(500.0);
//  TString plotNameLin = "plots/dataPt_"; plotNameLin+=fileNameDataIn; plotNameLin+=sSel; plotNameLin+="Lin"; //.png";
//  cdatapt.Print(plotNameLin+".png"); 
//  cdatapt.Print(plotNameLin+".eps"); 
//  cdatapt.Print(plotNameLin+".pdf"); 
  
//   framePt->SetAxisRange(30,90);
   //cout << "chi2/dof = " << framePt->chiSquare(2) << " chi2 = " << framePt->chiSquare(2)*(90-ptcutLow-1-2) << " dof = " << 90-ptcutLow-1-2 << " p = "<< TMath::Prob(framePt->chiSquare(2)*(90-ptcutLow-1-2),90-ptcutLow-1-2) << endl;
//   cout << "chi2/dof = " << framePt->chiSquare("mdlPdf_Norm[muonPt]","h_"+fileNameDataIn+".root",5) << " chi2 = " << framePt->chiSquare("mdlPdf_Norm[muonPt]","h_"+fileNameDataIn+".root",5)*(90-ptcutLow-5) << " dof = " << 90-ptcutLow-5 << " p = "<< TMath::Prob(framePt->chiSquare("mdlPdf_Norm[muonPt]","h_"+fileNameDataIn+".root",5)*(90-ptcutLow-5),90-ptcutLow-5) << endl;
  if ( (plotLL) || (plotPLL)) {
    RooPlot* frameLL = nsig.frame(Bins(20),Range(nsig.getVal()-3*nsig.getError(),nsig.getVal()+3*nsig.getError()),Title("LL and profileLL in nsig")) ;
    TCanvas* cdataLL = new TCanvas("cdataLL","dataLL",600,600);
    if (plotLL) {
      std::cout << "Plotting likelihood..." << std::endl;
      double eps = 1e-4;
      nll->plotOn(frameLL,ShiftToZero()) ;
      //nll->plotOn(frameLL,Precision(eps),ShiftToZero()) ;
     }
  // Plot the profile likelihood in nsig
    if (plotPLL) {
      std::cout << "Plotting profile likelihood..." << std::endl;
      RooAbsReal* pll_nsig = nll->createProfile(nsig) ;
      pll_nsig->plotOn(frameLL,LineColor(kRed)) ;
    }
  // Adjust frame maximum for visual clarity
    frameLL->SetMinimum(0) ;
    frameLL->SetMaximum(3) ;
    frameLL->Draw();
    TLatex l2;
    l2.SetNDC();
    l2.DrawLatex(0.35,0.75,sSel + ( addPercent ? "%" : "" ));
    l2.DrawLatex(0.35,0.65,sSel2);
    l2.DrawLatex(0.6,0.75,"#sqrt{s_{NN}}=2.76 TeV");
    l2.DrawLatex(0.6,0.85,"#int L #approx 0.140 #mub^{-1}"); 
    //#ifdef __CINT__
    ATLAS_LABEL(0.17,0.85);
    myText(0.35,0.85, (Color_t)kBlack, (char*)("Internal"));
    //myText(0.35,0.85, (Color_t)kBlack, (char*)("Preliminary"));
    //myText(0.25,0.70, (Color_t)kRed, (char*)("For Approval"));
    //#endif

    TString plotNameLL = "dataLL_"; plotNameLL+=fileNameDataIn; plotNameLL+=sSel; plotNameLL+=sSel2;
    //cdataLL->Print(plotNameLL+".png");
    cdataLL->Print(plotNameLL.ReplaceAll("#","")+".png");
    //cdataLL->Print(plotNameLL+".eps");
    cdataLL->Print(plotNameLL.ReplaceAll("#","")+".eps");
    cdataLL->Print(plotNameLL.ReplaceAll("#","")+".pdf");
    //cdataLL->Print(plotNameLL+".pdf"); 
    //return;
    if (frameLL)  delete frameLL;
    if (cdataLL) delete cdataLL;
   }

  cout << "Clean up" << endl;
  if (roofitResult) delete roofitResult;
  //if (roofitResult2) delete roofitResult2;
  if (leg) delete leg;
  //if (tf1BkgPdf) delete tf1BkgPdf;
  if (hdummy0) delete hdummy0;
  if (hdummy1) delete hdummy1;
  if (hdummy2) delete hdummy2;
  if (hdummy3) delete hdummy3;
  if (framePt) delete framePt;
  if (framePull) delete framePull;
//  if (frameRes)delete frameRes;
//  if (frameResSig)delete frameResSig;
//  if (frameResBkg)delete frameResBkg;
  //if (frameContour) delete frameContour;
  if (cdatapt) delete cdatapt;
  if (cPull)    delete cPull;
//  if (cRes)    delete cRes;
//  if (cResSig)    delete cResSig;
//  if (cResBkg)    delete cResBkg;
  cout << "Fit ends here" << endl;
}

void WAnalysis(){
  int cutvalue = 11; //2011 DPD 
  bool doFit = true; //turn off when building templates
  bool doTemplateSystematics = false;
  bool doCharge =  false ; //turn on when building or using a charge dependent template
  bool doCentrality = false ;
  bool doEta = false ;  //turn on when building eta templates
  bool doBuildTemplates = false;  //turn on when building ANY template 
  bool onlyBuildEtaTemplates = false; //turn on when building ONLY eta dependent templates (charge dependent or independent)

  int nSigPdf = 1000000;
  TString nSigPdfString;

  if(doTemplateSystematics&&doCentrality&&doCharge&&doEta){
      nSigPdfString = "CentChrgEta_Systematics1";
    }
  else if(doTemplateSystematics&&doCentrality&&doCharge&&!doEta){
      nSigPdfString = "CentChrg_Systematics7";
    }
  else if (doTemplateSystematics&&doEta&&doCentrality&&!doCharge){
      nSigPdfString = "EtaCent_Systematics1";
    }

  else if (doTemplateSystematics&&!doEta&&!doCentrality&&!doCharge){
      nSigPdfString = "Inclusive_Systematics";
    }
  else if (doTemplateSystematics&&doEta&&doCharge&&!doCentrality){
      nSigPdfString = "EtaChrg_Systematics10";
   }

  else if (doTemplateSystematics&&!doEta&&doCharge&&!doCentrality){
      nSigPdfString = "Charge_Systematics";
   }
  ///rename when doing simultaneous fits
  else if(doCentrality&&doCharge&&doEta){
      nSigPdfString = "CentChrgEta";
    }
  else if(doCentrality&&doCharge&&!doEta){
      nSigPdfString = "CentChrg";
    }
  ///rename when doing simultaneous fits
  else if (doEta&&doCentrality&&!doCharge){
      nSigPdfString = "CentEta10";
    }

  else if (doEta&&doCharge&&!doCentrality){
      nSigPdfString = "EtaChrg";
   }

  else if (!doEta&&doCharge&&!doCentrality){
      nSigPdfString = "Charge";
   }
  else{ nSigPdfString = "1000000";}

  float ptmax = 90.0;
 
  SetAtlasStyle();

//7.3mub 
  //TString fileNameDataIn = "HISingleMuon_Aug4_v2";//tightened impact paramter (new d0 z0 cuts)

  //TString fileNameDataIn = "HISingleMuonHP.2012.11.22"; 
  TString fileNameDataIn = "HISingleMuonHP.2012.11.26"; //matched to mu4 and mu10 (choose at event level, not muon level)

  TString fileNameMCIn = "HISingleMuonWmunuMCDataOverlay.2012.11.21";

  TString fileNamePdf = "WAnalysis_templates"; fileNamePdf+=nSigPdf; //fileNamePdf+=".root";
  TString fileNameFitOut = "WAnalysis_fitResult"; fileNameFitOut+=nSigPdfString; fileNameFitOut+=".root";

  RooRealVar  muonPt("muonPt","p_{T}",0,100,"GeV");
  RooRealVar  muonEta("muonEta","muonEta",-3.0,+3.0);
  RooRealVar  muonELoss("muonELoss","muonELoss",-0.5,+0.5);
  RooRealVar  muonScat("muonScat","muonScat",-4.0,+4.0);
  RooRealVar  centrality("centrality","centrality",0.,1.0);
  RooRealVar  muonCharge("muonCharge","charge",-2.0,+2.0);
  RooRealVar  ptCone20("ptCone20","ptCone20",0.0,+2.0);

  //RooFormulaVar extraCut("extraCut", "extraCut", "muonPt>7.0&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&ptCone20/muonPt<0.3&&abs(muonEta)<2.5&&centrality<=0.8", RooArgList(muonPt,muonELoss,muonScat,ptCone20,muonEta,centrality));
  RooFormulaVar extraCut("extraCut", "extraCut", "muonPt>7.0&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.5&&centrality<=0.8", RooArgList(muonPt,muonELoss,muonScat,muonEta,centrality));

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

  RooArgSet muonArgSet(muonPt,muonEta,muonCategory,chargeCategory,bkgCutCategory,centrality);
  
    // --- Fill data sets ---
  RooDataSet* dataSet0 = fillHIMuonDataSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/HardProbesFiles/",fileNameDataIn+".root",muonArgSet, cutvalue); dataSet0->Print();
  RooDataSet* dataSet = (RooDataSet*)dataSet0->reduce(extraCut);
  delete dataSet0;
  dataSet = (RooDataSet*)dataSet->reduce(Cut("bkgCutCategory==bkgCutCategory::ZDY"));dataSet->Print();
  
  RooPlot* framePt=muonPt.frame(Range(0,ptmax)) ;
  framePt->SetXTitle("#font[52]{p}_{T}^{#mu} [GeV]");
  framePt->SetYTitle("Muons/GeV");
  if (!doBuildTemplates) dataSet->plotOn(framePt, RefreshNorm(), DrawOption("p"), MarkerSize(0.8));
  
  /// Build templates
  if (doBuildTemplates) {
    std::cout << "Building Templates" << std::endl;
    RooDataSet* mcSet = fillHIMuonDataSet("/tmp/tbalestr/",fileNameMCIn+".root", muonArgSet, cutvalue, true); mcSet->Print();
    //RooDataSet* mcSet = fillHIMuonDataSet("root://eosatlas//eos/atlas/user/t/tbalestr/",fileNameMCIn+".root", muonArgSet, cutvalue, true); mcSet->Print();
    RooDataSet* mcWSet = (RooDataSet*)mcSet->reduce(Cut("muonCategory==muonCategory::W"), /*Cut(extraCut),*/ EventRange(0, nSigPdf) ); mcWSet->Print();
    RooDataSet* mcWSetPlus = (RooDataSet*)mcWSet->reduce(Cut("chargeCategory==chargeCategory::muPlus") ); mcWSetPlus->Print();
    RooDataSet* mcWSetMinus = (RooDataSet*)mcWSet->reduce(Cut("chargeCategory==chargeCategory::muMinus") ); mcWSetMinus->Print();
    RooDataSet* mcWSetEta1 = 0;
    RooDataSet* mcWSetEta2 = 0;
    RooDataSet* mcWSetEta3 = 0;
    RooDataSet* mcWSetEta4 = 0;
    RooDataSet* mcWSetEta5 = 0;
    RooDataSet* mcWSetEta6 = 0;
    RooDataSet* mcWSetEta7 = 0;
    RooDataSet* mcWSetEta8 = 0;
    RooDataSet* mcWSetEta9 = 0;
    RooDataSet* mcWSetEta10 = 0;
    RooDataSet* mcWSetPlusEta1 = 0;
    RooDataSet* mcWSetPlusEta2 = 0;
    RooDataSet* mcWSetPlusEta3 = 0;
    RooDataSet* mcWSetPlusEta4 = 0;
    RooDataSet* mcWSetPlusEta5 = 0;
    RooDataSet* mcWSetPlusEta6 = 0;
    RooDataSet* mcWSetPlusEta7 = 0;
    RooDataSet* mcWSetPlusEta8 = 0;
    RooDataSet* mcWSetPlusEta9 = 0;
    RooDataSet* mcWSetPlusEta10 = 0;
    RooDataSet* mcWSetMinusEta1 = 0;
    RooDataSet* mcWSetMinusEta2 = 0;
    RooDataSet* mcWSetMinusEta3 = 0;
    RooDataSet* mcWSetMinusEta4 = 0;
    RooDataSet* mcWSetMinusEta5 = 0;
    RooDataSet* mcWSetMinusEta6 = 0;
    RooDataSet* mcWSetMinusEta7 = 0;
    RooDataSet* mcWSetMinusEta8 = 0;
    RooDataSet* mcWSetMinusEta9 = 0;
    RooDataSet* mcWSetMinusEta10 = 0; 

    ///commented out for centrality binning of charge asymmetry (17Apr)
    if (doEta) {
       mcWSetEta1 = (RooDataSet*)mcWSet->reduce(etaCut1); mcWSetEta1->Print();
       mcWSetEta2 = (RooDataSet*)mcWSet->reduce(etaCut2); mcWSetEta2->Print();
       mcWSetEta3 = (RooDataSet*)mcWSet->reduce(etaCut3); mcWSetEta3->Print();
       mcWSetEta4 = (RooDataSet*)mcWSet->reduce(etaCut4); mcWSetEta4->Print();
       mcWSetEta5 = (RooDataSet*)mcWSet->reduce(etaCut5); mcWSetEta5->Print();
      ///commented out for centrality binning of charge asymmetry (17Apr)
      mcWSetEta6 = (RooDataSet*)mcWSet->reduce(etaCut6); mcWSetEta6->Print();
      mcWSetEta7 = (RooDataSet*)mcWSet->reduce(etaCut7); mcWSetEta7->Print();
      mcWSetEta8 = (RooDataSet*)mcWSet->reduce(etaCut8); mcWSetEta8->Print();
      mcWSetEta9 = (RooDataSet*)mcWSet->reduce(etaCut9); mcWSetEta9->Print();
      mcWSetEta10 = (RooDataSet*)mcWSet->reduce(etaCut10); mcWSetEta10->Print();

      mcWSetPlusEta1 = (RooDataSet*)mcWSetPlus->reduce(etaCut1); mcWSetPlusEta1->Print();
      mcWSetPlusEta2 = (RooDataSet*)mcWSetPlus->reduce(etaCut2); mcWSetPlusEta2->Print();
      mcWSetPlusEta3 = (RooDataSet*)mcWSetPlus->reduce(etaCut3); mcWSetPlusEta3->Print();
      mcWSetPlusEta4 = (RooDataSet*)mcWSetPlus->reduce(etaCut4); mcWSetPlusEta4->Print();
      mcWSetPlusEta5 = (RooDataSet*)mcWSetPlus->reduce(etaCut5); mcWSetPlusEta5->Print();
      ///commented out for centrality binning of charge asymmetry (17Apr)
       mcWSetPlusEta6 = (RooDataSet*)mcWSetPlus->reduce(etaCut6); mcWSetPlusEta6->Print();
      mcWSetPlusEta7 = (RooDataSet*)mcWSetPlus->reduce(etaCut7); mcWSetPlusEta7->Print();
      mcWSetPlusEta8 = (RooDataSet*)mcWSetPlus->reduce(etaCut8); mcWSetPlusEta8->Print();
      mcWSetPlusEta9 = (RooDataSet*)mcWSetPlus->reduce(etaCut9); mcWSetPlusEta9->Print();
      mcWSetPlusEta10 = (RooDataSet*)mcWSetPlus->reduce(etaCut10); mcWSetPlusEta10->Print();

       mcWSetMinusEta1 = (RooDataSet*)mcWSetMinus->reduce(etaCut1); mcWSetMinusEta1->Print();
       mcWSetMinusEta2 = (RooDataSet*)mcWSetMinus->reduce(etaCut2); mcWSetMinusEta2->Print();
       mcWSetMinusEta3 = (RooDataSet*)mcWSetMinus->reduce(etaCut3); mcWSetMinusEta3->Print();
       mcWSetMinusEta4 = (RooDataSet*)mcWSetMinus->reduce(etaCut4); mcWSetMinusEta4->Print();
       mcWSetMinusEta5 = (RooDataSet*)mcWSetMinus->reduce(etaCut5); mcWSetMinusEta5->Print();

      ///commented out for centrality binning of charge asymmetry (17Apr)
       mcWSetMinusEta6 = (RooDataSet*)mcWSetMinus->reduce(etaCut6); mcWSetMinusEta6->Print();
       mcWSetMinusEta7 = (RooDataSet*)mcWSetMinus->reduce(etaCut7); mcWSetMinusEta7->Print();
       mcWSetMinusEta8 = (RooDataSet*)mcWSetMinus->reduce(etaCut8); mcWSetMinusEta8->Print();
       mcWSetMinusEta9 = (RooDataSet*)mcWSetMinus->reduce(etaCut9); mcWSetMinusEta9->Print();
       mcWSetMinusEta10 = (RooDataSet*)mcWSetMinus->reduce(etaCut10); mcWSetMinusEta10->Print();
      
    }
    mcWSet->plotOn(framePt, DrawOption("p"), MarkerSize(0.8)); 
if (!onlyBuildEtaTemplates) {
    buildWTemplate( mcWSet, muonPt, fileNamePdf+".root" );
    buildWTemplate( mcWSetPlus, muonPt, fileNamePdf+"Plus.root" );
    buildWTemplate( mcWSetMinus, muonPt, fileNamePdf+"Minus.root" );
  }
    if (doEta) {
      std::cout<<"creating Root files for eta templates"<<std::endl;
      buildWTemplate(mcWSetEta1,muonPt,fileNamePdf+"EtaBin1.root");
      buildWTemplate(mcWSetEta2,muonPt,fileNamePdf+"EtaBin2.root");
      buildWTemplate(mcWSetEta3,muonPt,fileNamePdf+"EtaBin3.root");
      buildWTemplate(mcWSetEta4,muonPt,fileNamePdf+"EtaBin4.root");
      buildWTemplate(mcWSetEta5,muonPt,fileNamePdf+"EtaBin5.root");
      buildWTemplate(mcWSetEta6,muonPt,fileNamePdf+"EtaBin6.root");
      buildWTemplate(mcWSetEta7,muonPt,fileNamePdf+"EtaBin7.root");
      buildWTemplate(mcWSetEta8,muonPt,fileNamePdf+"EtaBin8.root");
      buildWTemplate(mcWSetEta9,muonPt,fileNamePdf+"EtaBin9.root");
      buildWTemplate(mcWSetEta10,muonPt,fileNamePdf+"EtaBin10.root");

      buildWTemplate(mcWSetPlusEta1,muonPt,fileNamePdf+"EtaBinPlus1.root");
      buildWTemplate(mcWSetPlusEta2,muonPt,fileNamePdf+"EtaBinPlus2.root");
      buildWTemplate(mcWSetPlusEta3,muonPt,fileNamePdf+"EtaBinPlus3.root");
      buildWTemplate(mcWSetPlusEta4,muonPt,fileNamePdf+"EtaBinPlus4.root");
      buildWTemplate(mcWSetPlusEta5,muonPt,fileNamePdf+"EtaBinPlus5.root");
      buildWTemplate(mcWSetPlusEta6,muonPt,fileNamePdf+"EtaBinPlus6.root");
      buildWTemplate(mcWSetPlusEta7,muonPt,fileNamePdf+"EtaBinPlus7.root");
      buildWTemplate(mcWSetPlusEta8,muonPt,fileNamePdf+"EtaBinPlus8.root");
      buildWTemplate(mcWSetPlusEta9,muonPt,fileNamePdf+"EtaBinPlus9.root");
      buildWTemplate(mcWSetPlusEta10,muonPt,fileNamePdf+"EtaBinPlus10.root");

      buildWTemplate(mcWSetMinusEta1,muonPt,fileNamePdf+"EtaBinMinus1.root");
      buildWTemplate(mcWSetMinusEta2,muonPt,fileNamePdf+"EtaBinMinus2.root");
      buildWTemplate(mcWSetMinusEta3,muonPt,fileNamePdf+"EtaBinMinus3.root");
      buildWTemplate(mcWSetMinusEta4,muonPt,fileNamePdf+"EtaBinMinus4.root");
      buildWTemplate(mcWSetMinusEta5,muonPt,fileNamePdf+"EtaBinMinus5.root");
      buildWTemplate(mcWSetMinusEta6,muonPt,fileNamePdf+"EtaBinMinus6.root");
      buildWTemplate(mcWSetMinusEta7,muonPt,fileNamePdf+"EtaBinMinus7.root");
      buildWTemplate(mcWSetMinusEta8,muonPt,fileNamePdf+"EtaBinMinus8.root");
      buildWTemplate(mcWSetMinusEta9,muonPt,fileNamePdf+"EtaBinMinus9.root");
      buildWTemplate(mcWSetMinusEta10,muonPt,fileNamePdf+"EtaBinMinus10.root");
      
    }
  }
  
  /// Read in templates
  TFile* fPdf      = new TFile(fileNamePdf+".root", "READ");
  TFile* fPdfPlus  = new TFile(fileNamePdf+"Plus.root", "READ");
  TFile* fPdfMinus = new TFile(fileNamePdf+"Minus.root", "READ");
      ///commented out for centrality binning of charge asymmetry (17Apr)
  TFile* fPdfEtaPlus[10];
  TFile* fPdfEtaMinus[10];
  TFile* fPdfEta[10];

  if (doEta || doTemplateSystematics) 
  {
    char fnamePlus2[400];
    char fnameMinus2[400];
    char fname2[400];
    
    for (int l = 0; l < 10; l++) {
       std::sprintf(fnamePlus2,"EtaBinPlus%i",l+1);
       std::sprintf(fnameMinus2,"EtaBinMinus%i",l+1);
       std::sprintf(fname2,"EtaBin%i",l+1);
       fPdfEtaPlus[l]  = new TFile(fileNamePdf+fnamePlus2+".root", "READ");
       fPdfEtaMinus[l]  = new TFile(fileNamePdf+fnameMinus2+".root", "READ");
       fPdfEta[l]  = new TFile(fileNamePdf+fname2+".root", "READ");
       if ( (fPdfEta[l] == 0) || (fPdfEtaPlus[l] == 0) || (fPdfEtaMinus[l] == 0) ) {
         std::cout << "Eta template file not found" << std::endl;
         return;
       }
    }
  }
  if ( (fPdf == 0) || (fPdfPlus == 0) || (fPdfMinus == 0) ) {
    std::cout << "Template file not found" << std::endl;
    return;
  }

  RooWorkspace* w      = (RooWorkspace*) fPdf->Get("w");
  RooWorkspace* wPlus  = (RooWorkspace*) fPdfPlus->Get("w");
  RooWorkspace* wMinus = (RooWorkspace*) fPdfMinus->Get("w");
  RooWorkspace* wEtaPlus[10];
  RooWorkspace* wEtaMinus[10];
  RooWorkspace* wEta[10];

  if (doEta || doTemplateSystematics){ 
    for (int k = 0; k < 10; k++) {
       wEtaPlus[k]  = (RooWorkspace*) fPdfEtaPlus[k]->Get("w");
       wEtaMinus[k]  = (RooWorkspace*) fPdfEtaMinus[k]->Get("w");
       wEta[k]  = (RooWorkspace*) fPdfEta[k]->Get("w");
    }
  }
//retrieve the pdfs
  RooKeysPdf& sigPdf        = *((RooKeysPdf*) w->pdf("sigPdf") );
  RooKeysPdf& sigPdfPlus    = *((RooKeysPdf*) wPlus->pdf("sigPdf") );
  RooKeysPdf& sigPdfMinus   = *((RooKeysPdf*) wMinus->pdf("sigPdf") );
  RooKeysPdf* sigPdfEtaPlus[10];
  RooKeysPdf* sigPdfEtaMinus[10];
  RooKeysPdf* sigPdfEta[10];

  if (doEta||doTemplateSystematics) {
     for (int j = 0; j < 10; j++) {
       sigPdfEtaPlus[j]    = ((RooKeysPdf*) wEtaPlus[j]->pdf("sigPdf") );
       sigPdfEtaMinus[j]    = ((RooKeysPdf*) wEtaMinus[j]->pdf("sigPdf") );
       sigPdfEta[j]    = ((RooKeysPdf*) wEta[j]->pdf("sigPdf") );
    }
  }
  if (doBuildTemplates) {
   if (!onlyBuildEtaTemplates) {
    sigPdf.plotOn(framePt,LineColor(kGreen));
    sigPdfPlus.plotOn(framePt,LineColor(kRed));
    sigPdfMinus.plotOn(framePt,LineColor(kBlue));
   }
    if (doEta||doTemplateSystematics) {
      for (int m = 0; m < 10; m++) {
         sigPdfEtaPlus[m]->plotOn(framePt,LineColor(kRed));
         sigPdfEtaMinus[m]->plotOn(framePt,LineColor(kBlue));
         sigPdfEta[m]->plotOn(framePt,LineColor(kGreen));
        }
     }
    TCanvas* cWpt = new TCanvas("cWpt","Wpt",600,600);
    framePt->Draw();
    cWpt->Print("wPt.png"); cWpt->Print("wPt.eps");
    return;
  }
  
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
  //  centralityBins.push_back(0.60);
    
  }
  centralityBins.push_back(0.80);

  const int nCentralityBins = centralityBins.size()-1;
  
    // --- Subdivide in bins ---
  RooDataSet* dataSubSet[nPtBins][nEtaBins][nCentralityBins];
  if(doEta||doCentrality){
  for ( int i = 0; i < nPtBins; i++ ) {
    for ( int j = 0; j < nEtaBins; j++ ) {
      for ( int k = 0; k < nCentralityBins; k++ ) {
        std::cout << "Divide " << i << ":" <<j << ":" <<k << std::endl;
        dataSubSet[i][j][k] = selectPtEtaCentrality( dataSet, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); dataSubSet[i][j][k]->Print();
      }
    }
  }
 }

    if (doFit) {
    /// --- Open output file ---
    TDirectory *dir = gDirectory;
    TFile outFile(fileNameFitOut,"RECREATE");
    gDirectory = dir;
    TString resultName = "WPtFit";
    FitResult baselineResult(resultName,ptBins,etaBins,centralityBins);


    if(!doCentrality && !doEta && !doCharge){
    /// --- fit inclusive spectra --
        if(!doTemplateSystematics){
        //fit( dataSet, muonPt, baselineResult, 99,0,0, sigPdf, ptmax, fileNameDataIn, "#mu^{#pm}" , "" );
        fit( dataSet, muonPt, baselineResult, 99,0,0, sigPdf, ptmax, fileNameDataIn, "#mu^{#pm},0-80" , "0 #leq |#eta| < 2.5" );
        }

        // --- for inclusive spectra systematics --- //
        if(doTemplateSystematics){
          std::cout<<"fitting systematics for inclusive muon group..."<<std::endl;
          fit( dataSet, muonPt, baselineResult, 99,0,0, sigPdfPlus, ptmax, fileNameDataIn, "#mu^{#pm}" , "systematics" );
          fit( dataSet, muonPt, baselineResult, 99,0,0, sigPdfMinus, ptmax, fileNameDataIn, "#mu^{#pm}" , "systematics" );

          for ( int j = 0; j < 10; j++ ) {
              RooKeysPdf& SigPdfEta = *(sigPdfEta[j]);
              fit( dataSet, muonPt, baselineResult, 99,0,0, SigPdfEta, ptmax, fileNameDataIn, "#mu^{#pm}" , "systematics" );
            }

          //if (dataSet) delete dataSet;
          return;
        } //inclusive template systematics

     if(!doCharge) return;
    }
    // --- fitting for charge-inclusive spectra over all eta and centrality --- //
    if (doCharge&&!doCentrality&&!doEta) {
      std::cout << "creating muMinus dataset" << std::endl;
      RooDataSet* dataSetMinus = (RooDataSet*) dataSet->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
      std::cout << "creating muPlus dataset" << std::endl;
      RooDataSet* dataSetPlus  = (RooDataSet*) dataSet->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
      if (dataSet) delete dataSet;
     

if (!doTemplateSystematics) {
        std::cout << "fitting over all ranges for mu-... " << std::endl;
        fit( dataSetMinus, muonPt, baselineResult, 100,0,0, sigPdfMinus, ptmax, fileNameDataIn, "#mu^{-},0-80" , "0 #leq |#eta| < 2.5");
        //fit( dataSetMinus, muonPt, baselineResult, 100,0,0, sigPdfMinus, ptmax, fileNameDataIn, "#mu^{-}" , "");
        std::cout << "done." << std::endl;
        if (dataSetMinus) delete dataSetMinus;
        std::cout << "fitting over all ranges for mu+... " << std::endl;
        fit( dataSetPlus,  muonPt, baselineResult, 101,0,0, sigPdfPlus, ptmax, fileNameDataIn, "#mu^{+},0-80" , "0 #leq |#eta| < 2.5");
        //fit( dataSetPlus,  muonPt, baselineResult, 101,0,0, sigPdfPlus, ptmax, fileNameDataIn, "#mu^{+}" , "");
        std::cout << "done." << std::endl;
        if (dataSetPlus) delete dataSetPlus;
        }

      if(doTemplateSystematics){
        std::cout<<"fitting systematics..."<<std::endl;
        std::cout<<"fitting systematics for mu-..."<<std::endl;
        fit( dataSetMinus, muonPt, baselineResult, 100,0,0, sigPdfPlus, ptmax, fileNameDataIn, "#mu^{-}" , "systematics");
        fit( dataSetMinus, muonPt, baselineResult, 100,0,0, sigPdf, ptmax, fileNameDataIn, "#mu^{-}" , "systematics");
        for ( int j = 0; j < 10; j++ ) {
            RooKeysPdf& SigPdfEta = *(sigPdfEta[j]);
            fit( dataSetMinus, muonPt, baselineResult, 99,0,0, SigPdfEta, ptmax, fileNameDataIn, "#mu^{-}" , "systematics" );
           }
        //if (dataSetMinus) delete dataSetMinus;
        std::cout<<"Done! Fitting systematics for mu+..."<<std::endl;
        fit( dataSetPlus,  muonPt, baselineResult, 101,0,0, sigPdfMinus, ptmax, fileNameDataIn, "#mu^{+}" , "systematics");
        fit( dataSetPlus,  muonPt, baselineResult, 101,0,0, sigPdf, ptmax, fileNameDataIn, "#mu^{+}" , "systematics");
        for ( int j = 0; j < 10; j++ ) {
            RooKeysPdf& SigPdfEta = *(sigPdfEta[j]);
            fit( dataSetPlus, muonPt, baselineResult, 99,0,0, SigPdfEta, ptmax, fileNameDataIn, "#mu^{+}" , "systematics" );
           }
        //if (dataSetPlus) delete dataSetPlus;
        std::cout<<"Done!"<<std::endl;
        return;
        }
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
         // if (doCentrality){
            sLow += 100*centralityBins[k]; sLow.Remove(3);
            //sLow += centralityBins[k]; sLow.Remove(3);
            sHigh += 100*centralityBins[k+1]; sHigh.Remove(3);
            //sHigh += centralityBins[k+1]; sHigh.Remove(3);
         // }
         // if (doEta){
            sLow2 += etaBins[j];
            sHigh2 += etaBins[j+1];
         // }


	    if (doCharge) {
            std::cout << "Creating charged datasets." << std::endl;
            RooDataSet* dataSetPlus  = (RooDataSet*) dataSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
            RooDataSet* dataSetMinus = (RooDataSet*) dataSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
            //if (dataSubSet[i][j][k]) delete dataSubSet[i][j][k] ;
            std::cout << "datasets created" << std::endl;
            TString sSelPlus = "#mu^{+},";
            TString sSelMinus = "#mu^{-},";
            TString sSelPlus2 = sLow2;
            TString sSelMinus2 = sLow2;
            RooKeysPdf& SigPdfEtaPlus = *(sigPdfEtaPlus[j]);
            RooKeysPdf& SigPdfEtaMinus = *(sigPdfEtaMinus[j]);

            if (doEta&&doCentrality){
              sSelPlus += sLow; sSelPlus+="-"; sSelPlus+= sHigh;
              sSelMinus += sLow; sSelMinus+="-"; sSelMinus+= sHigh;
              sSelMinus2 += " #leq "; sSelMinus2+= "|#eta| "; sSelMinus2+="< "; sSelMinus2 += sHigh2;
              sSelPlus2 += " #leq "; sSelPlus2+= "|#eta| "; sSelPlus2+="< "; sSelPlus2 += sHigh2;
              if (!doTemplateSystematics){
                std::cout << "Binning for eta+charge+centrality." << std::endl;
                std::cout << "Using eta+charge dependent templates." << std::endl;
                std::cout << "fitting for mu+: " << i << ":" << j << ":" << k <<std::endl;
                fit( dataSetPlus, muonPt, baselineResult, 102+i,j,k, SigPdfEtaPlus, ptmax, fileNameDataIn, sSelPlus, sSelPlus2);
                std::cout << "Done!" << std::endl;
                std::cout << "fitting for mu-: " << i << ":" << j << ":" << k <<std::endl;
                fit( dataSetMinus, muonPt, baselineResult, 103+i,j,k, SigPdfEtaMinus, ptmax, fileNameDataIn, sSelMinus, sSelMinus2);
                std::cout << "Done!" << std::endl;
              }
            }
            if (doEta&&!doCentrality) {
                sSelPlus += sLow; sSelPlus+="-"; sSelPlus+= sHigh;
                sSelMinus += sLow; sSelMinus+="-"; sSelMinus+= sHigh;
                TString sSelMinus2 = sLow2; sSelMinus2 += " #leq "; sSelMinus2+= "|#eta| "; sSelMinus2+="< "; sSelMinus2 += sHigh2;
                TString sSelPlus2 = sLow2; sSelPlus2 += " #leq "; sSelPlus2+= "|#eta| "; sSelPlus2+="< "; sSelPlus2 += sHigh2;
                if(!doTemplateSystematics){
                std::cout << "Now fitting with eta+charge dependent templates." << std::endl;
                std::cout << "fitting for mu+: " << i << ":" << j << ":" << k <<std::endl;
                fit( dataSetPlus, muonPt, baselineResult, 102+i,j,k, SigPdfEtaPlus, ptmax, fileNameDataIn, sSelPlus, sSelPlus2);
                //if(dataSetPlus) delete dataSetPlus ;
                std::cout << "Done!" << std::endl;
                std::cout << "fitting for mu-: " << i << ":" << j << ":" << k <<std::endl;
                fit( dataSetMinus, muonPt, baselineResult, 103+i,j,k, SigPdfEtaMinus, ptmax, fileNameDataIn, sSelMinus, sSelMinus2);
                //if (dataSetMinus) delete dataSetMinus;
                std::cout << "Done!" << std::endl;
                }
                if(doTemplateSystematics){
                std::cout << "fitting template systematics for eta and charge binning..." << std::endl;
                sSelPlus+= " systematics";
                std::cout << "Fitting systematics for mu+" << std::endl;

                fit( dataSetPlus, muonPt, baselineResult, 102+i,j,k, sigPdf, ptmax, fileNameDataIn, sSelPlus, sSelPlus2);
                fit( dataSetPlus, muonPt, baselineResult, 102+i,j,k, sigPdfPlus, ptmax, fileNameDataIn, sSelPlus, sSelPlus2);


 		for ( int l = 0; l < 10; l++ ) {
                RooKeysPdf& SigPdfEta = *(sigPdfEta[l]);
                fit( dataSetPlus, muonPt, baselineResult, 102+i,j,k, SigPdfEta, ptmax, fileNameDataIn, sSelPlus, sSelPlus2);
                }
                //if(dataSetPlus) delete dataSetPlus ;
                std::cout << "systematics for mu+ done." << std::endl;
                sSelMinus+= "systematics";
                std::cout << "Fitting systematics for mu-" << std::endl;
                fit( dataSetMinus, muonPt, baselineResult, 103+i,j,k, sigPdf, ptmax, fileNameDataIn, sSelPlus, sSelPlus2);
                fit( dataSetMinus, muonPt, baselineResult, 103+i,j,k, sigPdfPlus, ptmax, fileNameDataIn, sSelPlus, sSelPlus2);
                for ( int l = 0; l < 10; l++ ) {
                RooKeysPdf& SigPdfEta = *(sigPdfEta[l]);
                fit( dataSetMinus, muonPt, baselineResult, 103+i,j,k, SigPdfEta, ptmax, fileNameDataIn, sSelPlus, sSelPlus2);
                }
                //if (dataSetMinus) delete dataSetMinus ;
                std::cout << "systematics for mu- done." << std::endl;
                /*if(k == (nCentralityBins-1)){
                std::cout << "closing templates files..."<< std::endl;
                fPdfEtaPlus[j]->Close();
                fPdfEtaMinus[j]->Close();
                fPdfEta[j]->Close();
                std::cout << "done."<< std::endl;
                }*/
               }
              }  //eta
             if (doCentrality&&!doEta){
                sSelPlus += sLow; sSelPlus+="-"; sSelPlus+= sHigh;
                sSelMinus += sLow; sSelMinus+="-"; sSelMinus+= sHigh;
                sSelMinus2 += " #leq "; sSelMinus2+= "|#eta| "; sSelMinus2+="< "; sSelMinus2 += sHigh2;
                sSelPlus2 += " #leq "; sSelPlus2+= "|#eta| "; sSelPlus2+="< "; sSelPlus2 += sHigh2;
                if (!doTemplateSystematics){
                  std::cout << "Fitting with eta-indep-charge-dep templates" << std::endl;
                  std::cout << "fitting for mu-: " << i << ":" << j << ":" << k <<std::endl;
                  fit( dataSetMinus, muonPt, baselineResult, 103+i,j,k, sigPdfMinus, ptmax, fileNameDataIn, sSelMinus, sSelMinus2 );
                  //if(dataSetMinus) delete dataSetMinus ;
                  std::cout << "Done!" << std::endl;
                  std::cout << "fitting for mu+: " << i << ":" << j << ":" << k <<std::endl;
                  fit( dataSetPlus,  muonPt, baselineResult, 102+i,j,k, sigPdfPlus, ptmax, fileNameDataIn, sSelPlus, sSelPlus2 );
                  //if(dataSetPlus) delete dataSetPlus ;
                  std::cout << "Done!" << std::endl;
                }
                if (doTemplateSystematics){
                   std::cout << "fitting systematics for centrality and charge binning..." << std::endl;
                   sSelMinus+= " systematics";
                   std::cout << "fitting systematics for mu-" << std::endl;
                   fit( dataSetMinus, muonPt, baselineResult, 103+i,j,k, sigPdf, ptmax, fileNameDataIn, sSelMinus, sSelMinus2 );
                   fit( dataSetMinus, muonPt, baselineResult, 103+i,j,k, sigPdfPlus, ptmax, fileNameDataIn, sSelMinus, sSelMinus2 );
                   for ( int l = 0; l < 10; l++ ) {
                     RooKeysPdf& SigPdfEta = *(sigPdfEta[l]);
                     fit( dataSetMinus, muonPt, baselineResult, 103+i,j,k, SigPdfEta, ptmax, fileNameDataIn, sSelMinus, sSelMinus2 );
                   }
                   //if(dataSetMinus) delete dataSetMinus ;
                   std::cout << "fitting systematics for mu- done." << std::endl;
                   sSelPlus+= " systematics";
                   std::cout << "fitting systematics for mu+" << std::endl;
                   fit( dataSetPlus,  muonPt, baselineResult, 102+i,j,k, sigPdf, ptmax, fileNameDataIn, sSelPlus, sSelPlus2 );
   

                   fit( dataSetPlus,  muonPt, baselineResult, 102+i,j,k, sigPdfMinus, ptmax, fileNameDataIn, sSelPlus, sSelPlus2 );
                   for ( int l = 0; l < 10; l++ ) {
                     RooKeysPdf& SigPdfEta = *(sigPdfEta[l]);
                     fit( dataSetPlus, muonPt, baselineResult, 102+i,j,k, SigPdfEta, ptmax, fileNameDataIn, sSelPlus, sSelPlus2 );
                   }
                   //if (dataSetPlus) delete dataSetPlus ;
                   std::cout << "fitting systematics for mu+ done." << std::endl;
                   /*   if(k == (nCentralityBins-1)){
                      std::cout << "closing templates files..." <<std::endl;
                      fPdfEtaPlus[j]->Close();
                      fPdfEtaMinus[j]->Close();
                      fPdfEta[j]->Close();
                      std::cout << "done."<<std::endl;
                  }*/
               }
              }  //centrality
            } //charge

             if (!doCharge&&doCentrality&&doEta) {
                std::cout << "Fitting for centrality and eta classes." << std::endl;
                TString sSel = "#mu^{#pm},"; sSel += sLow; sSel+="-"; sSel+= sHigh;
                TString sSel2 = sLow2; sSel2 += " #leq "; sSel2+= "|#eta| "; sSel2+="< "; sSel2 += sHigh2;
                if (!doTemplateSystematics){
                  ///Change j index to appropriate template
                  ///if fitting for a specific eta window
                  ///(e.g. j+1 is the 0.25-0.5 window)
                  //RooKeysPdf& SigPdfEta = *(sigPdfEta[j+12]);
                  RooKeysPdf& SigPdfEta = *(sigPdfEta[j]);
                  fit(dataSubSet[i][j][k], muonPt, baselineResult, 102+i,j,k, SigPdfEta, ptmax, fileNameDataIn, sSel, sSel2 );
                  std::cout << "Fit completed." << std::endl;
                }
                if (doTemplateSystematics) {
                  std::cout << "fitting systematics for eta and centrality binning..." << std::endl;
                  sSel+=" systematics";
                  fit(dataSubSet[i][j][k], muonPt, baselineResult, 102+i,j,k, sigPdf, ptmax, fileNameDataIn, sSel, sSel2);
                  fit(dataSubSet[i][j][k], muonPt, baselineResult, 102+i,j,k, sigPdfMinus, ptmax, fileNameDataIn, sSel, sSel2);
                  fit(dataSubSet[i][j][k], muonPt, baselineResult, 102+i,j,k, sigPdfPlus, ptmax, fileNameDataIn, sSel, sSel2);
                  for ( int l = 0; l < 10; l++ ) {
                    if (l!=j) {
                      RooKeysPdf& SigPdfEta = *(sigPdfEta[l]);
                      fit(dataSubSet[i][j][k], muonPt, baselineResult, 102+i,j,k, SigPdfEta, ptmax, fileNameDataIn, sSel, sSel2);
                     }
                   }
                  std::cout << "systematics completed." << std::endl;
                }
             }  //eta and centrality
            if ( dataSubSet[i][j][k] ) delete dataSubSet[i][j][k];
        }
      }
    }

    outFile.cd();
    baselineResult.write(outFile);
    outFile.Close();
  }
}

int main() {
  WAnalysis();
}
/// eof













































