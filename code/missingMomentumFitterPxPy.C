#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "TChain.h"
#include "RooBinning.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"
#include "RooGExpModel.h"
#include "RooTruthModel.h"
#include "RooDecay.h"
//#include "RooVoigtian.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooAbsReal.h"
#include "TAxis.h"
//#include "AtlasUtils.C"
//#include "AtlasStyle.C"

using namespace RooFit ;

///////////////////////////////////////////////////////////
//fillDataSet
///////////////////////////////////////////////////////////
RooDataSet* fillDataSet(const TString& fileName, RooArgSet& muonArgSet,float missingPThr,float centBinLo, float centBinUp , bool hasMPt=false) {
       
  RooDataSet* data = new RooDataSet(fileName,fileName,muonArgSet) ;

  float missPtNt;
  float missPxNt;
  float missPyNt;

  float centralityNt;
  int   trigger1 = 0;
  int   trigger2 = 0;
  centralityNt = 0.5;
  float mPx = 0.0;
  float mPy = 0.0;

  TChain* tree = new TChain("data","data");
  int nFiles = tree->Add(fileName);
  std::cout <<"Filling the RooDataSet for "<< fileName << "...Number of files: " << nFiles << std::endl;
  
  
  TString missPxThrCut = "nu_px"; missPxThrCut+=missingPThr;     
  TString missPyThrCut = "nu_py"; missPyThrCut+=missingPThr; 
  TString missPtThrCut = "nu_pt"; missPtThrCut+=missingPThr; 

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("nu*",1);
  tree->SetBranchStatus("Centrality",1);
  tree->SetBranchStatus("EF_L1TE50_NoAlg",1);
  tree->SetBranchStatus("EF_mbZdc_a_c_L1VTE50_trk",1);

  tree->SetBranchAddress("Centrality", &centralityNt);
  tree->SetBranchAddress(missPtThrCut, &missPtNt);
  tree->SetBranchAddress(missPxThrCut, &missPxNt);
  tree->SetBranchAddress(missPyThrCut, &missPyNt);
  if(hasMPt) tree->SetBranchAddress(missPxTruth, &mPx);
  if(hasMPt) tree->SetBranchAddress(missPyTruth, &mPy);
  tree->SetBranchAddress("EF_L1TE50_NoAlg",&trigger1);
  tree->SetBranchAddress("EF_mbZdc_a_c_L1VTE50_trk",&trigger2);
 
  std::cout << "Number of entries: " << tree->GetEntries() << std::endl; 

    for (int i = 0; i < tree->GetEntries(); i++){
     //if(i==30000000) break; //temp hack
     tree->GetEntry(i);
     if (trigger1!=0||trigger2!=0){
      muonArgSet.setRealValue("centrality",centralityNt) ;

      if(centralityNt>=centBinLo&&centralityNt<centBinUp&&missPtNt<200.0) {
        
        muonArgSet.setRealValue("missPt",missPt500Nt);
        muonArgSet.setRealValue("missPx",missPx500Nt);
        muonArgSet.setRealValue("missPy",missPy500Nt);
	muonArgSet.setRealValue("missPxTruth",mPx);
	muonArgSet.setRealValue("missPyTruth",mPy);
        data->add(muonArgSet) ;
      }
     }
  }
    return data ;

}

//////////////////////////
//fit////////////////////
////////////////////////
float fitPx(RooRealVar& missPx,RooRealVar& missPxTruth){

  RooPlot* framePx = missPx.frame(Range(-90.0,+90.0));
  framePx->SetXTitle("#slash{p_{x}} [GeV]");

  //calculate the differene from truth
  RooGenericPdf("resPx","missPx-truth",RooArgList(missPx,truth)) ;

  RooBinning b0(-90.,90.) ;
  b0.addUniform(90,-90.,90.) ; 

//  TString arrMarkerColor[3] = {"kAzure","kBlue","kGreen"};

  //reco
  RooRealVar meanPx("meanPx","meanPx",0.0, -30.0,+30.0) ;
  RooRealVar sigmaPx("sigmaPx","sigmaPx",1.0,0.0,+30.0);

  // Build gaussian p.d.f in terms of missing px and py,mean and sigma
  std::cout << "Fitting for missing Px." << std::endl;
  RooGaussian* gaussPx = new RooGaussian("gaussPx","gaussPx",missPx,meanPx,sigmaPx) ;
  gaussPx->fitTo(*dataset,Range(-40.,40.));

  float resolution = sigmaPx.getVal() ;- 
  
}

void missingMomentumFitterPxPy(){
 #ifdef __CINT__
   gROOT->LoadMacro("AtlasUtils.C");
 #endif


  ///Declare variables
  RooRealVar missPt("missPt","missPt",0.,200.) ;
  RooRealVar  missPx("missPx","missPx",-200.,+200.);
  RooRealVar  missPy("missPy","missPy",-200.,+200.);
  RooRealVar  missPxTruth("missPxTruth","missPxTruth",-200.,+200.);
  RooRealVar  missPyTruth("missPyTruth","missPyTruth",-200.,+200.);
  
  RooRealVar  centrality("centrality","centrality",0.0,1.0);

  // Build gaussian p.d.f in terms of missing px and py,mean and sigma
  RooGaussian* gaussPx500 = new RooGaussian("gaussPx500","gaussPx500",missPx,meanG500,sigmaG500) ;
  RooGaussian* gaussPy500 = new RooGaussian("gaussPy500","gaussPy500",missPy,meanG500,sigmaG500) ;
  RooGaussian* gaussPx1000= new RooGaussian("gaussPx1000","gaussPx1000",missPx,meanG1000,sigmaG1000) ;
  RooGaussian* gaussPy1000= new RooGaussian("gaussPy1000","gaussPy1000",missPy,meanG1000,sigmaG1000) ;
  RooGaussian* gaussPx3000= new RooGaussian("gaussPx3000","gaussPx3000",missPx,meanG3000,sigmaG3000) ;
  RooGaussian* gaussPy3000= new RooGaussian("gaussPy3000","gaussPy3000",missPy,meanG3000,sigmaG3000) ;
 
  //build a truth resolution model (delta function)
  RooTruthModel tmPx("tmPx","truth model",missPx) ;
  RooTruthModel tmPy("tmPy","truth model",missPy) ;

  ///Data2011 from MB HIDPD ntuples 
  TString fileName = "/tmp/tbalestr/MinimumBias/*root*";

  //stl container to store lower 
  //track thresholds
	std::vector<int> missingPThr ;
	missingPThr.push_back(500);
	missingPThr.push_back(1000);
	missingPThr.push_back(3000);
	const unsigned int nMissPtBins = missingPThr.size();
	std::cout << " Number of lower track thresholds: " << nMissPtBins << std::endl;


  std::vector<double> centralityBins;
  centralityBins.push_back(0.00);
  centralityBins.push_back(0.10);
  centralityBins.push_back(0.20);
  centralityBins.push_back(0.40);
  centralityBins.push_back(0.80);
  const unsigned int nCentralityBins = centralityBins.size()-1;

  RooArgSet muonArgSet(missPt, missPx, missPy, centrality,missPxTruth,missPyTruth) ;

  ///create the dataset
  std::cout<<"Constructing datasets..." << std::endl;
  RooDataSet* dataSubSet[nMissPtBins][nCentralityBins];  
  for(int iPt=0;iPt<nMissPtBins;iPt++){
	for(int icent=0;icent<nCentralityBins;icent++){
		
		std::cout << "Divide " << iPt << ":" << icent << std::endl;
		if(dataHasTruMPT){
		  dataSubSet[iPt][icent] = fillDataSet(fileName,muonArgSet,missingPThr[iPt],centralityBins[k], centralityBins[k+1],true); dataSubSet[iPt][icent]->Print();
                }
		else{
		  dataSubSet[iPt][icent] = fillDataSet(fileName,muonArgSet,missingPThr[iPt],centralityBins[k], centralityBins[k+1],false); dataSubSet[iPt][icent]->Print();
		}
	}
  }

 
  ///plot the data 
  RooPlot* framePx=missPx.frame(Range(-90.,90.)) ;
  RooPlot* framePy=missPy.frame(Range(-90.,90.)) ;
  framePx500->SetXTitle("#slash{p_{x}} [GeV]");
  framePy500->SetXTitle("#slash{p_{y}} [GeV]");
  RooBinning b0(-90.,90.) ;
  b0.addUniform(90,-90.,90.) ; 

  TString arrMarkerColor[3] = {"kAzure","kBlue","kGreen"};

  for(int iPt=0;iPt<nMissPtBins;iPt++){
	for(int icent=0;icent<nCentralityBins;icent++){
		
		std::cout << "Fitting " << iPt << ":" << icent << std::endl;
		dataSubSet[iPt][icent]->plotOn(framePx, Binning(b0),DrawOption("p"), MarkerSize(0.8), MarkerColor(arrMarkerColor[iPt])) ;
		dataSubSet[iPt][icent]->plotOn(framePy, Binning(b0),DrawOption("p"), MarkerSize(0.8),MarkerColor(arrMarkerColor[iPt])) ;
		fit(dataSubSet[iPt][icent],missPx,iPt,icent);
		fit(dataSubSet[iPt][icent],missPy,iPt,icent);

	}
  }


  ///fit to gaussian
  std::cout << "Fitting missing Px..." << std::endl;
  gaussPx500->fitTo(*data500,Range(-90.,90.));
  gaussPx500->plotOn(framePx500,NumCPU(4),LineColor(kRed));
  gaussPx1000->fitTo(*data1000,Range(-90.,90.));
  gaussPx1000->plotOn(framePx1000,NumCPU(4),LineColor(kBlue));
  gaussPx3000->fitTo(*data3000,Range(-90.,90.));
  gaussPx3000->plotOn(framePx3000,NumCPU(4),LineColor(kGreen));
  std::cout << "Done." << std::endl;
  gaussPy3000->fitTo(*data3000,Range(-90.,90.));

  std::cout << "Fitting missing Py..." << std::endl;
  gaussPy500->fitTo(*data500,Range(-90.,90.));
  gaussPy500->plotOn(framePy500,NumCPU(4),LineColor(kAzure));
  gaussPy1000->fitTo(*data1000,Range(-90.,90.));
  gaussPy1000->plotOn(framePy1000,NumCPU(4),LineColor(kBlue));
  gaussPy3000->plotOn(framePy3000,NumCPU(4),LineColor(kGreen));
  std::cout << "Done." << std::endl;

  TCanvas* cPx = new TCanvas("cPx","cPx",600,600) ;
  //cPx->SetLogy(true) ;
  gPad->SetLeftMargin(0.15) ; framePx500->GetYaxis()->SetTitleOffset(1.6) ;    //framePx500->SetMaximum(1e7) ;
  framePx->Draw();

  TCanvas* cPy = new TCanvas("cPy","cPy",600,600) ;
  //cPy->SetLogy(true) ;
  gPad->SetLeftMargin(0.15) ; framePy500->GetYaxis()->SetTitleOffset(1.6) ;    //framePy500->SetMaximum(1e7) ;
  framePy->Draw();

}
