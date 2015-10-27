#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TMath.h"

#include <cmath>

using namespace RooFit ;


void zBosonFitter()
{
  const int nPts = 10;
  double arrPoints[nPts] = {2.73954,3.97732,4.9354,5.03371,5.56092,5.78564,4.4642,5.0097,4.31088,2.05925};
  double arrStat[nPts] = {0.430508,0.352002,0.33955,0.322585,0.328123,0.33505,0.302566,0.347129,0.378605,0.356732};
  double arrSyst[nPts] = {0.248308,0.243182,0.277298,0.296689,0.300932,0.317451,0.239293,0.267899,0.299923,0.220807};

  TH1D* h = new TH1D("h","h",10,-2.5,2.5);

  for(int i=1; i<=h->GetNbinsX(); i++){
	h->SetBinContent(i,arrPoints[i-1]);
	double error = TMath::Sqrt(arrStat[i-1]*arrStat[i-1]+arrSyst[i-1]*arrSyst[i-1]);
  	h->SetBinError(i,error);
  }

//  h->SaveAs("rapidityZ.root");

  //TH1* hh = (TH1*) gDirectory->Get("") ;
  RooRealVar yZ("yZ","yZ",-10,10) ;
  RooDataHist data("data","dataset with yZ",yZ,h) ;

  // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
  RooRealVar b("b","b",1.0,0.1,10) ;
  //RooRealVar b("b","b",-2,-10.0,-0.1) ;
  RooRealVar c("c","c",0.1,0.0,10) ;

  //RooGenericPdf mdl("mdl","mdl","TMath::Power(yZ,-1*b)*TMath::Exp(-1*c*TMath::Power(TMath::Log(yZ),2))",RooArgSet(yZ,b,c));
  RooGenericPdf mdl("mdl","mdl","-c*(yZ*yZ)+b",RooArgSet(yZ,b,c));
  
  // Build gaussian p.d.f in terms of x,mean and sigma
  //RooGaussian gauss("gauss","gaussian PDF",yZ,mean,sigma) ;  

  // Construct plot frame in 'x'
  RooPlot* xframe = yZ.frame(Title("p.d.f.")) ;

    // data and the p.d.f in the frame
  data.plotOn(xframe) ;
  mdl.plotOn(xframe) ;
    // F i t   m o d e l   t o   d a t a
  // -----------------------------

  // Fit pdf to data
  mdl.fitTo(data) ;

    // Print values of mean and sigma (that now reflect fitted values and errors)
  b.Print() ;
  c.Print() ;

  xframe->Draw() ;

  
}
