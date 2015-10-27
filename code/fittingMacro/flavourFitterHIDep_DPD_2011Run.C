#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooKeysPdf.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooShiftedKeysPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "TChain.h"
#include "TList.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TStyle.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
// #include "/afs/cern.ch/user/s/sandstro/atlasstyle-00-03-03/AtlasUtils.h"
// #ifndef __CINT__
// #include "/afs/cern.ch/user/s/sandstro/atlasstyle-00-03-03/AtlasStyle.C"
// #endif
//#include "RooShiftedKeysPdf.cxx"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>

using namespace RooFit;
using namespace RooStats;

//ClassImp(RooShiftedKeysPdf);
/*class RooShiftedKeysPdf {
public:
  //ClassDef(RooShiftedKeysPdf,1); // Your description goes here...
    RooShiftedKeysPdf() {} ;
  RooShiftedKeysPdf(const char *name, const char *title, RooAbsReal& x, RooAbsPdf& keysPdf, RooAbsReal& shift, RooAbsReal& stretch) :
  RooAbsPdf(name,title),
  m_xProxy("m_xProxy", "m_xProxy", this, x),
  m_x(static_cast<RooRealVar*>(&x)),
  m_keysPdf(&keysPdf),
  m_shift("m_shift","m_shift", this, shift),
  m_stretch("m_stretch", "m_stretch", this, stretch)
  {
   RooArgSet* deps = keysPdf.getDependents(*m_x);
   m_x = static_cast<RooRealVar*>(deps->find(m_x->GetName()));
   m_minValX = m_x->getMin();
   m_maxValX = m_x->getMax();
   m_mean    = getMean();
  }

RooShiftedKeysPdf(const RooShiftedKeysPdf& other, const char* name) :
RooAbsPdf(other,name),
m_xProxy("xProxy", this, other.m_xProxy),
m_x(other.m_x),
m_keysPdf(other.m_keysPdf),
m_shift("m_shift", this, other.m_shift),
m_stretch("m_stretch", this, other.m_stretch),
m_minValX(other.m_minValX),
m_maxValX(other.m_maxValX),
m_mean(other.m_mean){}
    virtual TObject* clone(const char* newname) const { return new RooShiftedKeysPdf(*this,newname); }
    inline virtual ~RooShiftedKeysPdf() { }
    protected:
    mutable RooRealProxy    m_xProxy;
    mutable RooRealVar*     m_x;
    RooAbsPdf*              m_keysPdf;
    RooRealProxy            m_shift;
    RooRealProxy            m_stretch;
    double                  m_minValX;
    double                  m_maxValX;
    double                  m_mean;
//    Double_t evaluate() const;
		Double_t evaluate() const
{
  double shiftedX = m_xProxy - m_shift;
  double shiftedMean = m_mean - m_shift;
  double newX = (shiftedX-shiftedMean)/m_stretch + shiftedMean;
  
  if (newX<m_minValX) {
  newX = m_minValX;
  }
  if (newX>m_maxValX) {
  newX = m_maxValX;
  }
  m_x->setVal(newX);
  RooArgSet set(*m_x);
  double returnVal = m_keysPdf->getVal( &set );
  return returnVal;
  }
	
//    Double_t getMean() const;
  Double_t getMean() const {
  int nSteps       = 200;
  double step      = (m_maxValX - m_minValX)/nSteps;
  double newX      = m_minValX + step/2.0;
  double integral  = 0;
  double integralX = 0;
  for (int i=0; i<nSteps; i++) {
  m_x->setVal(newX);
  double pdfEval = m_keysPdf->getVal();
  integral  += pdfEval;
  integralX += pdfEval*newX;
  newX += step;
  }
  if (integral==0) {
  std::cerr << "WARNING: RooShiftedKeysPdf::getMean(): integral = 0, this is not causing a crash, but expect wrong results..." << std::endl;
  return 0;
  }
  else {
  std::cout << "Calculated mean of KeysPdf = " << integralX/integral << std::endl;
  return integralX/integral;
  }
  assert(1);
  }


    //private:
    //ClassDef(RooShiftedKeysPdf,1) // Your description goes here...
  };


*/
		       
///////////////////////////////////////////////////////////////////////////////
// FitResult
///////////////////////////////////////////////////////////////////////////////
class FitResult
{
  public:
  
  FitResult( const char* name, std::vector<double>& ptBins, std::vector<double>& etaBins, std::vector<double>& centralityBins ) : _name(name) {
    const int nPtBins = ptBins.size()-1;
    _nEtaBins = etaBins.size()-1;
    _nCentralityBins = centralityBins.size()-1;
    const int nEtaBins = _nEtaBins;
    const int nCentralityBins = _nCentralityBins; 
    double x[nPtBins];
    double xErrLow[nPtBins];
    double xErrUpp[nPtBins];
    double y[nPtBins];
    for ( int i = 0; i < nPtBins; i++ ) {
      x[i] = (ptBins[i]+ptBins[i+1])*0.5;
      xErrLow[i] = x[i]-ptBins[i];
      xErrUpp[i] = ptBins[i+1]-x[i];
      y[i] = 0.0;
    }
    for ( int i = 0; i < nEtaBins; i++ ) {
      double eta = (etaBins[i]+etaBins[i+1])*0.5;
      _eta.SetPoint(i,eta,0.0);
      _eta.SetPointError(i,eta-etaBins[i],etaBins[i+1]-eta,0.0,0.0);
      for ( int j = 0; j < nCentralityBins; j++ ) {
        double centrality = (centralityBins[j]+centralityBins[j+1])*0.5;
        _centrality.SetPoint(j,centrality,0.0);
        _centrality.SetPointError(j,centrality-centralityBins[j],centralityBins[j+1]-centrality,0.0,0.0);

	// -- add objects to TList --//
        _nmc.Add( new TGraphAsymmErrors(nPtBins,x,y,xErrLow,xErrUpp) );
        _nsig.Add( new TGraphAsymmErrors(nPtBins,x,y,xErrLow,xErrUpp) );
        _mc.Add( new TGraphAsymmErrors(nPtBins,x,y,xErrLow,xErrUpp) );
        _toymc.Add( new TGraphAsymmErrors(nPtBins,x,y,xErrLow,xErrUpp) );
        _sig.Add( new TGraphAsymmErrors(nPtBins,x,y,xErrLow,xErrUpp) );
        _bK.Add( new TGraphAsymmErrors(nPtBins,x,y,xErrLow,xErrUpp) );
        _exp1K.Add( new TGraphAsymmErrors(nPtBins,x,y,xErrLow,xErrUpp) );
        _exp1Ratio.Add( new TGraphAsymmErrors(nPtBins,x,y,xErrLow,xErrUpp) );
        _nbkg.Add( new TGraphAsymmErrors(nPtBins,x,y,xErrLow,xErrUpp) );
        _chi2.Add( new TGraphAsymmErrors(nPtBins,x,y,xErrLow,xErrUpp) );
        _shift.Add( new TGraphAsymmErrors(nPtBins,x,y,xErrLow,xErrUpp) );
        _smear.Add( new TGraphAsymmErrors(nPtBins,x,y,xErrLow,xErrUpp) );
        _stretch.Add( new TGraphAsymmErrors(nPtBins,x,y,xErrLow,xErrUpp) );
      }
    }
  }
  ~FitResult() {
    std::cout << "dtor begins" << std::endl;
    for ( int i = 0; i < _sig.GetEntries(); i++ ) {
      if (_nmc.At(i))     delete _nmc.At(i);
      if (_nsig.At(i))    delete _nsig.At(i);
      if (_mc.At(i))      delete _mc.At(i);
      if (_toymc.At(i))   delete _toymc.At(i);
      if (_sig.At(i))     delete _sig.At(i);
      if (_bK.At(i))     delete _bK.At(i);
      if (_exp1K.At(i))     delete _exp1K.At(i);
      if (_exp1Ratio.At(i))     delete _exp1Ratio.At(i);
      if (_nbkg.At(i))     delete _nbkg.At(i);
      if (_chi2.At(i))    delete _chi2.At(i);
      if (_shift.At(i))   delete _shift.At(i);
      if (_smear.At(i))   delete _smear.At(i);
      if (_stretch.At(i)) delete _stretch.At(i);
    };
    std::cout << "graphs deleted" << std::endl;
    // --- Do not delete canvases ---
    for ( int i = 0; i < _canvas.GetEntries(); i++ ) {
      delete _canvas.At(i);
    }
    std::cout << "dtor ends" << std::endl;
  }
  
  void setPt( const int iPt, const int iEta, const int iCentrality, const double value ) {
    const int index = _nCentralityBins*iEta + iCentrality;

    setX((TGraphAsymmErrors*)_nmc.At(index),iPt,value);
    setX((TGraphAsymmErrors*)_nsig.At(index),iPt,value);
    setX((TGraphAsymmErrors*)_mc.At(index),iPt,value);
    setX((TGraphAsymmErrors*)_toymc.At(index),iPt,value);
    setX((TGraphAsymmErrors*)_sig.At(index),iPt,value);
    setX((TGraphAsymmErrors*)_chi2.At(index),iPt,value);
    setX((TGraphAsymmErrors*)_shift.At(index),iPt,value);
    setX((TGraphAsymmErrors*)_smear.At(index),iPt,value);
    setX((TGraphAsymmErrors*)_stretch.At(index),iPt,value);
  }
  
  void setMcN( const int iPt, const int iEta, const int iCentrality, const double value, const double errLow, const double errUpp ) {
    const int index = _nCentralityBins*iEta + iCentrality;
    setY( (TGraphAsymmErrors*)_nmc.At(index), iPt, value, errLow, errUpp );
  }
  void setSigN( const int iPt, const int iEta, const int iCentrality, const double value, const double errLow, const double errUpp ) {
    const int index = _nCentralityBins*iEta + iCentrality;
    setY( (TGraphAsymmErrors*)_nsig.At(index), iPt, value, errLow, errUpp );
  }
  void setMc( const int iPt, const int iEta, const int iCentrality, const double value, const double errLow, const double errUpp ) {
    const int index = _nCentralityBins*iEta + iCentrality;
    setY( (TGraphAsymmErrors*)_mc.At(index), iPt, value, errLow, errUpp );
  }
  void setToyMc( const int iPt, const int iEta, const int iCentrality, const double value, const double errLow, const double errUpp ) {
    const int index = _nCentralityBins*iEta + iCentrality;
    setY( (TGraphAsymmErrors*)_toymc.At(index), iPt, value, errLow, errUpp );
  }
  void setSig( const int iPt, const int iEta, const int iCentrality, const double value, const double errLow, const double errUpp ) {
    const int index = _nCentralityBins*iEta + iCentrality;
    setY( (TGraphAsymmErrors*)_sig.At(index), iPt, value, errLow, errUpp );
  }

  //thomas-functions for plotting nuisance parameters
  void setbK( const int iPt, const int iEta, const int iCentrality, const double value, const double errLow, const double errUpp ) {
    const int index = _nCentralityBins*iEta + iCentrality;
    setY( (TGraphAsymmErrors*)_bK.At(index), iPt, value, errLow, errUpp );
  }
  void setexp1K( const int iPt, const int iEta, const int iCentrality, const double value, const double errLow, const double errUpp ) {
    const int index = _nCentralityBins*iEta + iCentrality;
    setY( (TGraphAsymmErrors*)_exp1K.At(index), iPt, value, errLow, errUpp );
  }
  void setexp1Ratio( const int iPt, const int iEta, const int iCentrality, const double value, const double errLow, const double errUpp ) {
    const int index = _nCentralityBins*iEta + iCentrality;
    setY( (TGraphAsymmErrors*)_exp1Ratio.At(index), iPt, value, errLow, errUpp );
  }
  void setnbkg( const int iPt, const int iEta, const int iCentrality, const double value, const double errLow, const double errUpp ) {
    const int index = _nCentralityBins*iEta + iCentrality;
    setY( (TGraphAsymmErrors*)_nbkg.At(index), iPt, value, errLow, errUpp );
  }
  void setChi2( const int iPt, const int iEta, const int iCentrality, const double value, const double errLow, const double errUpp ) {
    const int index = _nCentralityBins*iEta + iCentrality;
    setY( (TGraphAsymmErrors*)_chi2.At(index), iPt, value, errLow, errUpp );
  }
  void setShift( const int iPt, const int iEta, const int iCentrality, const RooRealVar& var ) {
    const int index = _nCentralityBins*iEta + iCentrality;
    setVar( (TGraphAsymmErrors*)_shift.At(index), iPt, var );
  }
  void setSmear( const int iPt, const int iEta, const int iCentrality, const RooRealVar& var ) {
    const int index = _nCentralityBins*iEta + iCentrality;
    setVar( (TGraphAsymmErrors*)_smear.At(index), iPt, var );
  }
  void setStretch( const int iPt, const int iEta, const int iCentrality, const RooRealVar& var ) {
    const int index = _nCentralityBins*iEta + iCentrality;
    setVar( (TGraphAsymmErrors*)_stretch.At(index), iPt, var );
  }
  
  void addCanvas( TCanvas* c ) {
    _canvas.Add(c);
  }
  
  void write( TFile& f ) {
    TDirectory *dir = gDirectory;
    f.cd();
    _eta.Write( TString(_name)+="_eta" );
    _centrality.Write( TString(_name)+="_centrality" );
    for ( int i = 0; i < _nEtaBins; i++ ) {
      for ( int j = 0; j < _nCentralityBins; j++ ) {
        TString s = "_eta";
        s += i;
        s += "_centrality";
        s += j;
        const int index = i*_nCentralityBins+j;
        ((TGraphAsymmErrors*)_nmc.At(index))->Write( (TString(_name)+="_nmc")+=s );
        ((TGraphAsymmErrors*)_nsig.At(index))->Write( (TString(_name)+="_nsig")+=s );
        ((TGraphAsymmErrors*)_mc.At(index))->Write( (TString(_name)+="_mc")+=s );
        ((TGraphAsymmErrors*)_toymc.At(index))->Write( (TString(_name)+="_toymc")+=s );
        ((TGraphAsymmErrors*)_sig.At(index))->Write( (TString(_name)+="_sig")+=s );
        ((TGraphAsymmErrors*)_bK.At(index))->Write( (TString(_name)+="_bK")+=s );
        ((TGraphAsymmErrors*)_exp1K.At(index))->Write( (TString(_name)+="_exp1K")+=s );
        ((TGraphAsymmErrors*)_exp1Ratio.At(index))->Write( (TString(_name)+="_exp1Ratio")+=s );
        ((TGraphAsymmErrors*)_nbkg.At(index))->Write( (TString(_name)+="_nbkg")+=s );
        ((TGraphAsymmErrors*)_chi2.At(index))->Write( (TString(_name)+="_chi2")+=s );
        ((TGraphAsymmErrors*)_shift.At(index))->Write( (TString(_name)+="_shift")+=s );
        ((TGraphAsymmErrors*)_smear.At(index))->Write( (TString(_name)+="_smear")+=s );
        ((TGraphAsymmErrors*)_stretch.At(index))->Write( (TString(_name)+="_stretch")+=s );
      }
    }
    for ( int i = 0; i < _canvas.GetEntries(); i++ ) {
      TCanvas* c = (TCanvas*)_canvas.At(i);
      c->SetName( (TString(_name)+="_")+=(c->GetName()) );
      c->Write( c->GetName() );
    }
    gDirectory = dir;
  }
  
  private:
  
  void setX( TGraphAsymmErrors* graph, const int iPt, const double value ) {
    double x, y;
    graph->GetPoint(iPt,x,y);
    double xlow = x - graph->GetErrorXlow(iPt);
    double xupp = x + graph->GetErrorXhigh(iPt);
    graph->SetPoint(iPt,value,y);
    graph->SetPointEXlow(iPt,value-xlow);
    graph->SetPointEXhigh(iPt,xupp-value);
  }
  
  void setY( TGraphAsymmErrors* graph, const int iPt, const double value, const double errLow, const double errUpp ) {
    double x, y;
    graph->GetPoint(iPt,x,y);
    graph->SetPoint(iPt,x,value);
    graph->SetPointEYlow(iPt,errLow);
    graph->SetPointEYhigh(iPt,errUpp);
  }
  
  void setVar( TGraphAsymmErrors* graph, const int iPt, const RooRealVar& var ) {
    double value = var.getVal();
    double err   = var.getError();
    setY( graph, iPt, value, err, err );
  }
  
  const char* _name;
  TGraphAsymmErrors _eta;
  TGraphAsymmErrors _centrality;
  TList _nmc;
  TList _nsig;
  TList _mc;
  TList _toymc;
  TList _sig;
  TList _bK;
  TList _exp1K;
  TList _exp1Ratio;
  TList _nbkg;
  TList _chi2;
  TList _shift;
  TList _smear;
  TList _stretch;
  TList _canvas;
  int _nEtaBins;
  int _nCentralityBins;
};


///////////////////////////////////////////////////////////////////////////////
// getPtString
///////////////////////////////////////////////////////////////////////////////
TString getPtString( double ptLow, double ptUpp )
{
  TString str = "pT";
  str += int(ptLow);
  str += "_";
  str += int(ptUpp);
  return str;
}


///////////////////////////////////////////////////////////////////////////////
// getEtaString
///////////////////////////////////////////////////////////////////////////////
TString getEtaString( double etaLow, double etaUpp )
{
  TString str = "eta";
  str += int(100*etaLow);
  str += "_";
  str += int(100*etaUpp);
  return str;
}


///////////////////////////////////////////////////////////////////////////////
// getTemplateFileName
///////////////////////////////////////////////////////////////////////////////
TString getTemplateFileName( double ptLow, double ptUpp, double etaLow, double etaUpp, const char* source, bool doReweightPtSpectrum, bool doMirrorEta, RooRealVar& var )
{
  //TString fileName = "templates/templates";
  TString fileName = "templates";
  fileName += ".";
  fileName += getPtString(ptLow,ptUpp);
  fileName += ".";
  fileName += getEtaString(etaLow,etaUpp);
  fileName += ".";
  fileName += source;
  fileName += ".";
  fileName += "Reweight";
  if ( doReweightPtSpectrum ) fileName += "On";
  else fileName += "Off";
  fileName += "MirrorEta";
  if ( doMirrorEta ) fileName += "On";
  else fileName += "Off";
  fileName += var.getTitle();
  fileName += ".root";
  return fileName;
}

///////////////////////////////////////////////////////////////////////////////
// fillHIMuonDataSet
///////////////////////////////////////////////////////////////////////////////
//RooDataSet* fillHIMuonDataSet(const TString& pathName, const TString& fileName, RooArgSet& muonArgSet, int cutValue = 9,  bool isMc = false )
RooDataSet* fillHIMuonDataSet(const TString& pathName, const TString& fileName, RooArgSet& muonArgSet, int cutValue = 9,  bool isMc = false )
{
  RooDataSet* set = new RooDataSet(fileName,fileName,muonArgSet);
  
  float eLossNt[50];
  float scatNt[50];
  float compNt[50];
  double ptSmNt[50];
  float ptNt[50];
  float etaNt[50];
  float chargeNt[50];
  int promptNt[50];
  float centralityNt;
  float ptcone20Nt[50];
  int valNt[50], ZDYNt[50], efMatched[50];
  int nmu,trig1,trig2,trig3,trig4,trig5;
  

  // --- Load the MC tree ---
  TChain* tree = new TChain("tree","tree");
  int nFiles = tree->Add(pathName+fileName);
  
  std::cout << "Filling the RooDataSet for " << fileName << "... Number of files: " << nFiles << std::endl;
 
  // --- Set branch adresses ---
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC",&trig1);
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10",&trig2);
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20",&trig3);
  tree->SetBranchAddress("EF_mu4_MSonly_L1TE50",&trig4);
  tree->SetBranchAddress("EF_mu4_L1VTE50",&trig5);
  tree->SetBranchAddress("efMatched",&efMatched);
  tree->SetBranchAddress("eLoss", &eLossNt);
  //tree->SetBranchAddress("mujetdR", &dRNt);
  tree->SetBranchAddress("ptcone20", &ptcone20Nt);
  tree->SetBranchAddress("scat", &scatNt);
  tree->SetBranchAddress("comp", &compNt);
  // --- MC_pt_smearing --- //
  if(isMc) tree->SetBranchAddress("pTCB_smeared", &ptSmNt);
  else tree->SetBranchAddress("pt", &ptNt);
  tree->SetBranchAddress("eta", &etaNt);
  //use for MCP ntuples
  tree->SetBranchAddress("charge", &chargeNt);
  if(isMc) tree->SetBranchAddress("prompt", &promptNt);
  tree->SetBranchAddress("val", &valNt); 
  tree->SetBranchAddress("ZDY", &ZDYNt); 
  tree->SetBranchAddress("centrality", &centralityNt);
  tree->SetBranchAddress("mu_muid_n", &nmu);

   // --- Set branch status ---
  tree->SetBranchStatus("*",0) ;
  tree->SetBranchStatus("mu_muid_n", 1);
  tree->SetBranchStatus("eLoss", 1);
  tree->SetBranchStatus("mujetdR", 1);
  tree->SetBranchStatus("ptcone20", 1);
  tree->SetBranchStatus("scat", 1);
  tree->SetBranchStatus("comp", 1);
  // --- MC_pt_smearing --- //
  if(isMc) tree->SetBranchStatus("pTCB_smeared", 1);
  else tree->SetBranchStatus("pt", 1);
  if(isMc) tree->SetBranchStatus("mTCB_smeared", 1);
  tree->SetBranchStatus("eta", 1);
  tree->SetBranchStatus("charge", 1);
  if(isMc) tree->SetBranchStatus("prompt", 1);
  tree->SetBranchStatus("val", 1); 
  tree->SetBranchStatus("ZDY", 1); 
  tree->SetBranchStatus("centrality", 1);
  tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC",1);
  tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10",1);
  tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20",1);
  tree->SetBranchStatus("EF_mu4_MSonly_L1TE50",1);
  tree->SetBranchStatus("EF_mu4_L1VTE50",1);
  tree->SetBranchStatus("efMatched",1);
  
  std::cout << "Number of entries: " << tree->GetEntries() << std::endl; 
  for ( int i = 0; i < tree->GetEntries(); i++ ) {
     if (i>1000000) return set; // temp hack
    tree->LoadTree(i);
    tree->GetEntry(i);
    //std::cout << "aaaa"<<std::endl;

    muonArgSet.setRealValue("centrality",centralityNt);
   // if(trig1||trig2||trig3||trig4||trig5){
    for (int imu = 0; imu<nmu;imu++){
    if(valNt[imu]>cutValue&&fabs(eLossNt[imu])<0.5&&fabs(scatNt[imu])<4.0&&efMatched[imu]==1){
      if (isMc) {
        if ( promptNt[imu] == 0 ) muonArgSet.setCatLabel("muonCategory","Unknown");
        else if ( promptNt[imu] == 1 ) muonArgSet.setCatLabel("muonCategory","HeavyFlavour");
        else if ( promptNt[imu] == 2 ) muonArgSet.setCatLabel("muonCategory","Tau");
        else if ( promptNt[imu] == 3 ) muonArgSet.setCatLabel("muonCategory","LightResonance");
        else if ( promptNt[imu] == 4 ) muonArgSet.setCatLabel("muonCategory","PionCaloDecay");
        else if ( promptNt[imu] == 5 ) muonArgSet.setCatLabel("muonCategory","PionEarlyDecay");
        else if ( promptNt[imu] == 6 ) muonArgSet.setCatLabel("muonCategory","PionLateDecay");
        else if ( promptNt[imu] == 7 ) muonArgSet.setCatLabel("muonCategory","KaonCaloDecay");
        else if ( promptNt[imu] == 8 ) muonArgSet.setCatLabel("muonCategory","KaonEarlyDecay");
        else if ( promptNt[imu] == 9 ) muonArgSet.setCatLabel("muonCategory","KaonLateDecay");
        else if ( promptNt[imu] == 10) muonArgSet.setCatLabel("muonCategory","Fake");
        else if ( promptNt[imu] == 23) muonArgSet.setCatLabel("muonCategory","Z");
        else if ( promptNt[imu] == 24) muonArgSet.setCatLabel("muonCategory","W");
      }

      muonArgSet.setRealValue("muonELoss",eLossNt[imu]);
      muonArgSet.setRealValue("muonScat",scatNt[imu]);
      muonArgSet.setRealValue("composite",compNt[imu]);
      muonArgSet.setRealValue("ptCone20",ptcone20Nt[imu]);
      if(isMc) muonArgSet.setRealValue("muonPt",ptSmNt[imu]);
      else muonArgSet.setRealValue("muonPt",ptNt[imu]);
      muonArgSet.setRealValue("muonEta",etaNt[imu]);
      muonArgSet.setRealValue("muonCharge",chargeNt[imu]);
      muonArgSet.setRealValue("muonCategory",promptNt[imu]);
      if (ZDYNt[imu] == 0) muonArgSet.setCatLabel("bkgCutCategory","ZDY");
      if ( chargeNt[imu] > 0 ) muonArgSet.setCatLabel("chargeCategory","muPlus");
      else if ( chargeNt[imu] < 0) muonArgSet.setCatLabel("chargeCategory","muMinus");
      set->add(muonArgSet);    
    }
   }
  //}//trigger
  }
  return set;
}


///////////////////////////////////////////////////////////////////////////////
// reweightPtSpectrum
///////////////////////////////////////////////////////////////////////////////
RooDataSet* reweightPtSpectrum(RooDataSet* mcSet, RooDataSet* dataSet, const RooArgSet& muonArgSet, double ptLow, double ptHigh, const double binWidth )
{
  const int nSubDivisions = int((ptHigh-ptLow)/binWidth);
  RooDataSet* result = 0;
  
  std::cout << "Reweighting the pt spectrum... Number of sub-divisions: " << nSubDivisions << std::endl;
  
  // --- Fill pT histogram for data and MC --- 
  TH1F* h_ptData  = new TH1F("h_ptData", "h_ptData", nSubDivisions, ptLow, ptHigh);
  TH1F* h_ptMc    = new TH1F("h_ptMc", "h_ptMc", nSubDivisions, ptLow, ptHigh);
  
  dataSet->fillHistogram(h_ptData, *muonArgSet.find("muonPt"));
  mcSet->fillHistogram(h_ptMc, *muonArgSet.find("muonPt"));
  
  for (int iSubDiv=0; iSubDiv<nSubDivisions; iSubDiv++) {
    
    TString cut = "muonPt>=";
    cut += h_ptData->GetXaxis()->GetBinLowEdge(iSubDiv+1);
    cut += "&&muonPt<";
    cut += h_ptData->GetXaxis()->GetBinUpEdge(iSubDiv+1);
    
    if (iSubDiv==0) {
      // --- create new RooDataSet at first pass, don't skim events ---
      result = (RooDataSet*) mcSet->reduce(muonArgSet, cut);
      continue;
    }
    
    // --- Calculate weight (number of events to skim from MC) ---
    RooDataSet* mcSetSubBinned = (RooDataSet*) mcSet->reduce(muonArgSet, cut);
    double ratioData = double(h_ptData->GetBinContent(iSubDiv+1))/double(h_ptData->GetBinContent(1));
    int numEventsKeep = int(ratioData*h_ptMc->GetBinContent(1));
    
    // --- Reduce MC dataset ---
    RooDataSet* mcSetReduced = (RooDataSet*) mcSetSubBinned->reduce(EventRange(0, numEventsKeep));
    result->append(*mcSetReduced);
    
    // --- Debugging info ---
    if (true) {
      std::cout << "cut = " << cut << std::endl;
      std::cout << "mcSetSubBinnend->numEntries() = " << mcSetSubBinned->numEntries() << std::endl;
      std::cout << "Entries to keep = " << numEventsKeep << std::endl;
      std::cout << "Entries kept after reduction = " << mcSetReduced->numEntries() << std::endl;
    }
    
    // --- Checks ---
    if (numEventsKeep>mcSetSubBinned->numEntries()) {
      std::cout << "reweightPtSpectrum(): WARNING: more entries in data than in MC!" << std::endl;
      std::cout << "This reweighting method cannot be applied." << std::endl;
    }
    if (numEventsKeep!=mcSetReduced->numEntries()) {
      std::cout << "ERROR: numEventsKeep!=mcSetReduced->numEntries(): prepare for crash!" << std::endl;
      return 0;
    }
      
    // --- Clean up ---
    delete mcSetSubBinned;
    delete mcSetReduced;
    
  }
  
  delete h_ptData;
  delete h_ptMc;

  return result;
}


///////////////////////////////////////////////////////////////////////////////
// selectPtEta
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectPtEta( RooDataSet* dataSet, double ptLow, double ptUpp, double etaLow, double etaUpp, bool doMirrorEta = false )
{
  TString cut = "muonPt>";
  cut += ptLow;
  cut += "&&muonPt<";
  cut += ptUpp;
  cut += "&&((muonEta>";
  cut += etaLow;
  cut += "&&muonEta<";
  cut += etaUpp;
  cut += ")";
  if ( doMirrorEta ) {
    cut += "||(muonEta>";
    cut += -etaUpp;
    cut += "&&muonEta<";
    cut += -etaLow;
    cut += ")";
  }
  cut += ")";

  //std::cout<<"Selecting PtEta"<<std::endl;
  return (RooDataSet*) dataSet->reduce(cut);
  //RooDataSet* dataSetCut = (RooDataSet*) dataSet->reduce(cut);
  //delete dataSet;
  //return dataSetCut;
}

///////////////////////////////////////////////////////////////////////////////
// selectPtEta
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectPtEtaCentrality( RooDataSet* dataSet, double ptLow, double ptUpp, double etaLow, double etaUpp, double centralityLow, double centralityHigh, bool doMirrorEta = false )
{
  //std::cout<<"Selecting PtEtaCentrality"<<std::endl;
  dataSet = selectPtEta( dataSet, ptLow, ptUpp, etaLow, etaUpp, doMirrorEta);
  //std::cout<<"Selecting centrality"<<std::endl;
  TString cut = "centrality>";
  cut += centralityLow;
  cut += "&&centrality<";
  cut += centralityHigh;
  
  return (RooDataSet*) dataSet->reduce(cut);
  //RooDataSet* dataSetCut = (RooDataSet*) dataSet->reduce(cut);
  //delete dataSet;
  //return dataSetCut;
}


///////////////////////////////////////////////////////////////////////////////
// selectTrueSignal
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectTrueSignal( RooDataSet* mcSet )
{
  ///selectr heavy flavor, taus, light resonances, W, and Z
  return (RooDataSet*) mcSet->reduce("muonCategory>=1&&muonCategory<=3||muonCategory>=23");
}


///////////////////////////////////////////////////////////////////////////////
// selectTrueBackground
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectTrueBackground( RooDataSet* mcSet )
{
  //PionCaloDecay, PionEarlyDecay, PionLateDecay, KaonCaloDecay, KaonEarlyDecay, KaonLateDecay, Fakes
  return (RooDataSet*) mcSet->reduce("muonCategory>=4&&muonCategory<=10");
}


///////////////////////////////////////////////////////////////////////////////
// selectTruePions
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectTruePions( RooDataSet* mcSet )
{
  return (RooDataSet*) mcSet->reduce("muonCategory>=4&&muonCategory<=6");
}


///////////////////////////////////////////////////////////////////////////////
// selectTrueKaons
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectTrueKaons( RooDataSet* mcSet )
{
  return (RooDataSet*) mcSet->reduce("muonCategory>=7&&muonCategory<=9");
}


///////////////////////////////////////////////////////////////////////////////
// selectTrueOthers
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectTrueOthers( RooDataSet* mcSet )
{
  return (RooDataSet*) mcSet->reduce("muonCategory==10");
}


///////////////////////////////////////////////////////////////////////////////
// selectJpsis
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectJpsis( RooDataSet* dataSet, bool isSideBands = false, bool isSameSign = false, double massWindow = 0.10 )
{
  TString cut = "jpsiMass>2.0&&jpsiMass<5.0&&";
  cut += "abs(jpsiMass-3.097)";
  if ( isSideBands ) cut += ">";
  else cut += "<";
  cut += massWindow;
  if ( isSameSign ) cut += "&&jpsiCharge>0.0";
  else cut += "&&jpsiCharge<0.0";
  return (RooDataSet*) dataSet->reduce(cut);
}


///////////////////////////////////////////////////////////////////////////////
// selectKshorts
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectKshorts( RooDataSet* dataSet, bool isSideBands = false, double massWindow = 0.05 )
{
  TString cut = "ksMass>0.0&&ksMass<3.0&&";
  cut += "abs(ksMass-0.498)";
  if ( isSideBands ) cut += ">";
  else cut += "<";
  cut += massWindow;
  return (RooDataSet*) dataSet->reduce(cut);
}


///////////////////////////////////////////////////////////////////////////////
// plotTemplates
///////////////////////////////////////////////////////////////////////////////
TCanvas* plotTemplates( RooRealVar& var, RooDataSet* sigSet, RooDataSet* bkgSet, RooAbsPdf& sigPdf, RooAbsPdf& bkgPdf )
{
  int nBins = var.getBins();
  var.setBins(100);
  TCanvas* canvas = new TCanvas();
  canvas->cd();
  RooPlot* pdfTestFrame = var.frame();
  
  sigSet->plotOn(pdfTestFrame);
  sigPdf.plotOn(pdfTestFrame);
  bkgSet->plotOn(pdfTestFrame,MarkerStyle(kOpenCircle));
  bkgPdf.plotOn(pdfTestFrame,LineStyle(kDashed));
  double binWidth = pdfTestFrame->GetYaxis()->GetBinWidth(1);
  TString yTitle = "Entries / "; yTitle += binWidth;
  pdfTestFrame->GetXaxis()->SetTitle(var.getTitle());
  pdfTestFrame->GetYaxis()->SetTitle(yTitle);
//   pdfTestFrame->GetXaxis()->SetTitle("(#font[52]{p}_{ID}-#font[52]{p}_{MS}-#font[52]{p}_{param})/#font[52]{p}_{ID}");
//   pdfTestFrame->GetYaxis()->SetTitle("Entries / 0.02");
  pdfTestFrame->Draw();
  
  var.setBins(nBins);
  return canvas;
}


///////////////////////////////////////////////////////////////////////////////
// plotDistortions
///////////////////////////////////////////////////////////////////////////////
TCanvas* plotDistortions( RooRealVar& var, RooDataSet* sigSet, RooDataSet* bkgSet, RooAbsPdf& sigPdf, RooAbsPdf& bkgPdf,
                          RooRealVar& parameter, double x0, double x1, double x2, const char* name )
{
  int nBins = var.getBins();
  var.setBins(100);
  TCanvas* canvas = new TCanvas();
  canvas->cd();
  RooPlot* plot = var.frame();
  
  sigSet->plotOn(plot);
  parameter.setVal(x2);
  sigPdf.plotOn(plot, LineColor(kGreen-4));
  parameter.setVal(x1);
  sigPdf.plotOn(plot, LineColor(kCyan-3));
  parameter.setVal(x0);
  sigPdf.plotOn(plot, LineColor(kBlue));
  
  bkgSet->plotOn(plot, MarkerStyle(kOpenCircle));
  parameter.setVal(x2);
  bkgPdf.plotOn(plot, LineColor(kGreen-4), LineStyle(kDashed));
  parameter.setVal(x1);    
  bkgPdf.plotOn(plot, LineColor(kCyan-3), LineStyle(kDashed));
  parameter.setVal(x0);
  bkgPdf.plotOn(plot, LineColor(kBlue), LineStyle(kDashed));
  
  double binWidth = plot->GetYaxis()->GetBinWidth(1);
  TString yTitle = "Entries / "; yTitle += binWidth;
  plot->GetXaxis()->SetTitle(var.getTitle());
  plot->GetYaxis()->SetTitle(yTitle);
//   plot->GetXaxis()->SetTitle("(#font[52]{p}_{ID}-#font[52]{p}_{MS}-#font[52]{p}_{param})/#font[52]{p}_{ID}");
//   plot->GetYaxis()->SetTitle("Entries / 0.02");
  plot->Draw();
  
  TLegend* leg = new TLegend(0.2,0.8,0.5,0.925);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(gStyle->GetTextFont());
  leg->SetTextSize(gStyle->GetTextSize());
  TString str = name;
  str        += " = ";
  str        += x0;
  TLegendEntry* entry = leg->AddEntry("",str,"l");
  entry->SetLineColor(kBlue);
  entry->SetLineWidth(2);
  str  = name;
  str += " = ";
  str += x1;
  entry = leg->AddEntry("",str,"l");
  entry->SetLineColor(kCyan-3);
  entry->SetLineWidth(2);
  str  = name;
  str += " = ";
  str += x2;
  entry = leg->AddEntry("",str,"l");
  entry->SetLineColor(kGreen-4);
  entry->SetLineWidth(2);
  leg->Draw();
  
  var.setBins(nBins);
  return canvas;
}


///////////////////////////////////////////////////////////////////////////////
// plotInvariantMass
///////////////////////////////////////////////////////////////////////////////
TCanvas* plotInvariantMass( RooRealVar& mass, RooDataSet* dataOsSet, RooDataSet* dataSsSet = 0 )
{
  // --- Mass fit ---
  // --- Define peak Gaussian ---
  RooRealVar mean("mean","mean",(mass.getMin()+mass.getMax())/2.0,mass.getMin(),mass.getMax());
  RooRealVar sigma("sigma","sigma",0.05,1.0e-4,0.10);
  RooGaussian gauss("gauss","gauss",mass,mean,sigma);
  // --- Define exponential ---
  RooRealVar alpha("alpha","alpha",0.00);
  RooExponential exp("exp","exp",mass,alpha);
  // --- Sum ---
  RooRealVar nsig("nsig","nsig",100.0,0.0,1.0e6);
  RooRealVar nbkg("nbkg","nbkg",100.0,0.0,1.0e6);
  RooAddPdf sum("sum","sum",RooArgList(gauss,exp),RooArgList(nsig,nbkg));
  // --- Fit ---
  sum.fitTo(*dataOsSet);
  // --- Mass plot ---
  TCanvas* canvas = new TCanvas();
  canvas->cd();
  RooPlot* massPlot = mass.frame();
  dataOsSet->plotOn(massPlot);
  sum.plotOn(massPlot);
  sum.plotOn(massPlot,Components(gauss),LineStyle(kDashed));
  sum.plotOn(massPlot,Components(exp),LineStyle(kDotted));
  if ( dataSsSet != 0 ) {
    dataSsSet->plotOn(massPlot,MarkerStyle(kOpenCircle));
  }
  massPlot->GetXaxis()->SetTitle("#font[52]{M} [GeV]");
  massPlot->Draw();
  
  return canvas;
}


///////////////////////////////////////////////////////////////////////////////
// plotFitResult
///////////////////////////////////////////////////////////////////////////////
TCanvas* plotFitResult( RooRealVar& var, RooDataSet* dataSubSet, RooAbsPdf& mdlPdf, RooAbsPdf& sigConvPdf, RooAbsPdf& bkgConvPdf, double& chi2 )
{
  int nBins = var.getBins();
  var.setBins(100);
  TCanvas* canvas = new TCanvas();
  canvas->cd();
  RooPlot* fitFrame = var.frame();
  dataSubSet->plotOn(fitFrame);
  mdlPdf.plotOn(fitFrame);
  chi2 = fitFrame->chiSquare(4);
  mdlPdf.plotOn(fitFrame,Components(sigConvPdf),LineStyle(kDashed));
  mdlPdf.plotOn(fitFrame,Components(bkgConvPdf),LineStyle(kDotted));
//   double binWidth = fitFrame->GetYaxis()->GetBinWidth(1);
//   TString yTitle = "Entries / "; yTitle += binWidth;
  TString yTitle = "Entries / "; yTitle += 0.02;
  fitFrame->GetXaxis()->SetTitle(var.getTitle());
  fitFrame->GetYaxis()->SetTitle(yTitle);
//   fitFrame->GetXaxis()->SetTitle("(#font[52]{p}_{ID}-#font[52]{p}_{MS}-#font[52]{p}_{param})/#font[52]{p}_{ID}");
//   fitFrame->GetYaxis()->SetTitle("Entries / 0.02");
  fitFrame->Draw();
  TLegend *leg = new TLegend(0.2,0.7,0.5,0.925);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(gStyle->GetTextFont());
  leg->SetTextSize(gStyle->GetTextSize());
  TLegendEntry* entry = leg->AddEntry("","Data 2010 (#sqrt{#font[52]{s}} = 7 TeV)","p");
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1.2);
  entry = leg->AddEntry("","Best-fit template","l");
  entry->SetLineColor(kBlue);
  entry->SetLineWidth(2);
  entry = leg->AddEntry("","Prompt component","l");
  entry->SetLineColor(kBlue);
  entry->SetLineWidth(2);
  entry->SetLineStyle(kDashed);
  entry = leg->AddEntry("","Pion/kaon component","l");
  entry->SetLineColor(kBlue);
  entry->SetLineWidth(2);
  entry->SetLineStyle(kDotted);
  leg->Draw();
  
  var.setBins(nBins);
  return canvas;
}

///////////////////////////////////////////////////////////////////////////////
// plotFitResult
///////////////////////////////////////////////////////////////////////////////
RooPlot* plotFitResultRooPlot( RooRealVar& var, RooDataSet* dataSubSet, RooAbsPdf& mdlPdf, RooAbsPdf& sigConvPdf, RooAbsPdf& bkgConvPdf, double& chi2 )
{
  int nBins = var.getBins();
  var.setBins(100);
//   TCanvas* canvas = new TCanvas();
//   canvas->cd();
  RooPlot* fitFrame = var.frame();
  dataSubSet->plotOn(fitFrame);
  mdlPdf.plotOn(fitFrame);
  chi2 = fitFrame->chiSquare(4);
  mdlPdf.plotOn(fitFrame,Components(sigConvPdf),LineStyle(kDashed));
  mdlPdf.plotOn(fitFrame,Components(bkgConvPdf),LineStyle(kDotted));
//   double binWidth = fitFrame->GetYaxis()->GetBinWidth(1);
//   TString yTitle = "Entries / "; yTitle += binWidth;
  TString yTitle = "Entries / "; yTitle += 0.02;
  fitFrame->GetXaxis()->SetTitle(var.getTitle());
  fitFrame->GetYaxis()->SetTitle(yTitle);
//   fitFrame->GetXaxis()->SetTitle("(#font[52]{p}_{ID}-#font[52]{p}_{MS}-#font[52]{p}_{param})/#font[52]{p}_{ID}");
//   fitFrame->GetYaxis()->SetTitle("Entries / 0.02");
  fitFrame->Draw();
  TLegend *leg = new TLegend(0.2,0.7,0.5,0.925);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(gStyle->GetTextFont());
  leg->SetTextSize(gStyle->GetTextSize());
  TLegendEntry* entry = leg->AddEntry("","Data 2010 (#sqrt{#font[52]{s}} = 7 TeV)","p");
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1.2);
  entry = leg->AddEntry("","Best-fit template","l");
  entry->SetLineColor(kBlue);
  entry->SetLineWidth(2);
  entry = leg->AddEntry("","Prompt component","l");
  entry->SetLineColor(kBlue);
  entry->SetLineWidth(2);
  entry->SetLineStyle(kDashed);
  entry = leg->AddEntry("","Pion/kaon component","l");
  entry->SetLineColor(kBlue);
  entry->SetLineWidth(2);
  entry->SetLineStyle(kDotted);
  leg->Draw();
  
  var.setBins(nBins);
  return fitFrame;
}
///////////////////////////////////////////////////////////////////////////////
// buildTemplates
///////////////////////////////////////////////////////////////////////////////
void buildTemplates( RooDataSet* mcSet, RooRealVar& var, const char* fileName )
{
//   RooRealVar& muonELoss = (RooRealVar&) muonArgSet["muonELoss"];
  
  RooDataSet* sigSet      = selectTrueSignal(mcSet);     sigSet->Print();
  RooDataSet* bkgSet      = selectTrueBackground(mcSet); bkgSet->Print();
  RooDataSet* bkgPionSet  = selectTruePions(mcSet);
  RooDataSet* bkgKaonSet  = selectTrueKaons(mcSet);
  RooDataSet* bkgOtherSet = selectTrueOthers(mcSet);
  
  // --- Number of prompts, pions, kaons, others ---
  RooRealVar nSignal("nSignal","nSignal", sigSet->numEntries() );
  RooRealVar nPions("nPions","nPions",    bkgPionSet->numEntries() );
  RooRealVar nKaons("nKaons","nKaons",    bkgKaonSet->numEntries() );
  RooRealVar nOthers("nOthers","nOthers", bkgOtherSet->numEntries() );
  
  // --- Create a workspace ---
  RooWorkspace *w = new RooWorkspace("w","workspace");
  w->import(var);
  
  // --- Build kernel-estimated templates ---
  std::cout << "Building the RooKeysPdf's... This will take some time" << std::endl;
  RooKeysPdf sigPdf("sigPdf","sigPdf",var,*sigSet);
  RooKeysPdf bkgPionPdf("bkgPionPdf","bkgPionPdf",var,*bkgPionSet);
  RooKeysPdf bkgKaonPdf("bkgKaonPdf","bkgKaonPdf",var,*bkgKaonSet);
  RooKeysPdf bkgOtherPdf("bkgOtherPdf","bkgOtherPdf",var,*bkgOtherSet);

  // --- And write to file ---
  std::cout << "Writing the templates to file " << fileName << std::endl;
  w->import(sigPdf);
  w->import(bkgPionPdf);
  w->import(bkgKaonPdf);
  w->import(bkgOtherPdf);
  w->import(nSignal);
  w->import(nPions);
  w->import(nKaons);
  w->import(nOthers);
  w->writeToFile(fileName);
  std::cout << "Done!" << std::endl;
  
  // --- Clean up ---
  delete sigSet;
  delete bkgSet;
  delete bkgPionSet;
  delete bkgKaonSet;
  delete bkgOtherSet;
  delete w;
}

///comparing function
int compare(const void *a,const void *b){
  return (*(int*)a-*(int*)b);
}

///median fuction after sorting sequence
float getMedian(float arr[],int size)
{
  int middle = size/2;
  float average;
  if (size%2==0) average = static_cast<float>(arr[middle-1]+arr[middle])/2;
  else average = static_cast<float>(arr[middle]);
  
  cout << "Median of array is: " << average << endl;
  return average;
}

///rms fuction for a truncated array
float getRMS(float arr[],int size)
{
  float sum1 = 0.; float sum2 = 0.;
  int nOK = 0;
  for (int i = 1; i<size-1; ++i) {
    if((arr[i]>1e-5) && (1-arr[i]>1e-5)) {
      sum1+=arr[i];
      sum2+=arr[i]*arr[i];
      nOK++;
    }
  }
  float rms = 0.;
  if (sum2/nOK - (sum1*sum1)/(nOK*nOK) > 0.) rms = sqrt(sum2/nOK - (sum1*sum1)/(nOK*nOK));
  return rms;
}

int getMedianIndex(float arr[],int size)
{
  int nOK = 0;
  for (int i = 1; i<size-1; ++i) {
    if((arr[i]>1e-5) && (1-arr[i]>1e-5)) {
      nOK++;
    }
  }
  int aimFor = (int)(nOK/2);
  int index=0;
  int iOK = 0;
  for (int i = 1; i<size-1; ++i) {
    if((arr[i]>1e-5) && (1-arr[i]>1e-5)) {
      iOK++;
      if (iOK == aimFor)  {
	index = i;
	break;
      }
    }
  }
  
  return index;
}

///////////////////////////////////////////////////////////////////////////////
// fit
///////////////////////////////////////////////////////////////////////////////
void fit( RooDataSet* dataSubSet, RooDataSet* refSubSet, RooDataSet* mcSubSet, RooArgSet& muonArgSet, RooRealVar& variable,
          FitResult& fitResult, const int iPt, const int iEta, const int iCentrality,
          TFile* fPdf, double pionWeight = 1.0, bool doParametricModels = false, bool onlyDistortBkg = false )
{
  if (dataSubSet->numEntries()==0) { /// early abort for no data scenario
    fitResult.setSig(iPt, iEta, iCentrality,-1,0,0);
    fitResult.setToyMc(iPt, iEta, iCentrality,-1,0,0);
    fitResult.setMc(iPt, iEta, iCentrality,-1,0,0);
    fitResult.setSigN(iPt,iEta, iCentrality,0,0.0,0.0);
    fitResult.setMcN(iPt, iEta, iCentrality,0,0.0,0.0);
    fitResult.setPt(iPt, iEta, iCentrality,-1);
    return;
  }

  int nToyMC = 8;
  if (pionWeight != 1.0) nToyMC = 0; // no need to to toy mc for the manipulated mc sets
  float* arrSigRatio = new float[nToyMC+1];
  float* arrStretch = new float[nToyMC+1];
  float* arrShift = new float[nToyMC+1];
  float* arrSmear = new float[nToyMC+1];
  float* arrErr   = new float[nToyMC+1];

  // --- Define the RooFit variables ---
//   RooRealVar&  muonELoss = (RooRealVar&) muonArgSet["muonELoss"];
  RooRealVar&  muonPt    = (RooRealVar&) muonArgSet["muonPt"];
  
  RooDataSet* sigRefSet  = selectTrueSignal(refSubSet);
  RooDataSet* bkgRefSet  = selectTrueBackground(refSubSet);
  RooDataSet* sigMCSet  = selectTrueSignal(mcSubSet);
  RooDataSet* bkgMCSet  = selectTrueBackground(mcSubSet);
  int nSig = sigMCSet->numEntries();
  int nBkg = bkgMCSet->numEntries();

  /// if one template dominates, limit the statistics for the toy mc for the other
  if (nSig > 10*nBkg) nSig = 10*nBkg;
  if (nBkg > 10*nSig) nBkg = 10*nSig;

  RooWorkspace* w = (RooWorkspace*) fPdf->Get("w");
  // --- Kernel-estimated models ---
  RooKeysPdf& sigKerPdf   = *((RooKeysPdf*) w->pdf("sigPdf") );
  RooKeysPdf& bkgPionPdf  = *((RooKeysPdf*) w->pdf("bkgPionPdf") );
  RooKeysPdf& bkgKaonPdf  = *((RooKeysPdf*) w->pdf("bkgKaonPdf") );
  RooKeysPdf& bkgOtherPdf = *((RooKeysPdf*) w->pdf("bkgOtherPdf") );
  RooRealVar& nPions  = *((RooRealVar*) w->var("nPions"));
  RooRealVar& nKaons  = *((RooRealVar*) w->var("nKaons"));
  RooRealVar& nOthers = *((RooRealVar*) w->var("nOthers"));
  // --- Parametric models ---
//   RooAbsPdf* parSigPdf = (RooAbsPdf*) w->pdf("parSigPdf_ConvLG");   // --- Signal: Landau (x) Gaussian ---
//   RooAbsPdf* parBkgPdf = (RooAbsPdf*) w->pdf("parBkgPdf_ConvLG");   // --- Signal: Landau (x) Gaussian ---
//   RooAbsPdf* parSigPdf = (RooAbsPdf*) w->pdf("parSigPdf_DG");   // --- Signal: double Gaussian ---
//   RooAbsPdf* parBkgPdf = (RooAbsPdf*) w->pdf("parBkgPdf_DG");
  std::cout << "Done!" << std::endl;
  
  // --- Improve binning for better convolution ---
  variable.setBins(1000);
  
  // --- Shift PDF's ---
  RooRealVar shift("shift","shift",2e-2,-0.3,+0.3);
  RooRealVar stretch("stretch","stretch",1.0,+0.8,+1.3);
  // shift.setConstant(true);
  stretch.setConstant(true);
  
  // --- Background Shifted PDF ---
  // --- Set background contributions as in truth reference sample ---
  nPions.setVal( pionWeight*nPions.getVal() );
  nKaons.setVal( nKaons.getVal() );
  nOthers.setVal( nOthers.getVal() );
  // --- If the MC used to build the templates suffers of low statistics, then bkgOtherPdf will be empty ---
  // --- If this is the case force nOthers to zero to avoid warnings ---
  if ( bkgOtherPdf.getNorm(variable) < 1.0e-5 ) nOthers.setVal(0.0);
  // --- Build the total background pdf ---
  RooAddPdf bkgKerPdf("bkgPdf","bkgPdf",RooArgList(bkgPionPdf,bkgKaonPdf,bkgOtherPdf),RooArgList(nPions,nKaons,nOthers));
  
  RooAbsPdf* sigAbsPdf = 0;
  RooAbsPdf* bkgAbsPdf = 0;
//   RooArgSet* parSigVars = 0;
//   RooArgSet* parBkgVars = 0;
  if ( doParametricModels ) {
    // --- Set the parameters of the parametric PDF's constant ---
//     parSigVars = parSigPdf->getParameters(variable);
//     std::cout << "Printing parSigVars:" << std::endl;
//     parSigVars->Print();
//     parBkgVars = parBkgPdf->getParameters(variable);
//     parSigVars->setAttribAll("Constant",true);
//     parBkgVars->setAttribAll("Constant",true);
//     sigAbsPdf = parSigPdf;
//     bkgAbsPdf = parBkgPdf;
  }
  else {
    sigAbsPdf = &sigKerPdf;
    bkgAbsPdf = &bkgKerPdf;
  }
  RooAbsPdf& sigPdf = *sigAbsPdf;
  RooAbsPdf& bkgPdf = *bkgAbsPdf;
  
  // --- Shifted/stretched PDF's ---
  std::cout<<"RooShiftedKeysPdfBefore"<<std::endl;
  RooShiftedKeysPdf sigShiftPdf("sigShiftPdf","sigShiftPdf",variable,sigPdf,shift,stretch);
  std::cout<<"RooShiftedKeysPdfAfterFirst"<<std::endl;
  RooShiftedKeysPdf bkgShiftPdf("bkgShiftPdf","bkgShiftPdf",variable,bkgPdf,shift,stretch);
  std::cout<<"RooShiftedKeysPdfAfter"<<std::endl;
  
  // --- Smear PDF's ---
  std::cout << " gaussMean = " <<  sigKerPdf.getParameters(variable)->getRealValue("parSig_peak_mean") << std::endl;
//   RooRealVar gaussMean("gaussMean","gaussMean",sigKerPdf.getParameters(variable)->getRealValue("parSig_peak_mean"));
//   RooRealVar gaussMean("gaussMean","gaussMean",0.0,-1.0,1.0);
  RooRealVar gaussMean("gaussMean","gaussMean",0.);
  RooRealVar gaussSigma("gaussSigma","gaussSigma",2.0e-2,1.0e-5,1.0);
  // gaussSigma.setConstant(true);
  RooGaussian gauss("gauss","gauss",variable,gaussMean,gaussSigma);
  RooFFTConvPdf sigConvPdf("sigConvPdf","sigConvPdf",variable,sigShiftPdf,gauss);
  RooFFTConvPdf bkgConvPdf("bkgConvPdf","bkgConvPdf",variable,bkgShiftPdf,gauss);
//   RooFFTConvPdf sigConvPdf("sigConvPdf","sigConvPdf",variable,sigPdf,gauss);
//   RooFFTConvPdf bkgConvPdf("bkgConvPdf","bkgConvPdf",variable,bkgPdf,gauss);
  
  // --- Make a composite model ---
  double margin = 0.;
  RooRealVar sigRatio("sigRatio","sigRatio",0.5,margin,1.0-margin);
  
  RooArgList modelComponents;
  if (onlyDistortBkg) {
    modelComponents.add(sigPdf);
    modelComponents.add(bkgConvPdf);
//   } else if (doParametricModels) {
//     modelComponents.add(sigPdf);
//     modelComponents.add(bkgPdf);
  }  else {
    modelComponents.add(sigConvPdf);
    modelComponents.add(bkgConvPdf);
//     modelComponents.add(sigPdf);
//     modelComponents.add(bkgPdf);
//     modelComponents.add(sigShiftPdf);
//     modelComponents.add(bkgShiftPdf);
  }
  
  RooAddPdf  mdlPdf("mdlPdf","mdlPdf",modelComponents,RooArgList(sigRatio));
  
  // --- Plot examples of shift, stretch, smear ---
//   TCanvas *cTemplate = plotTemplates(variable,sigRefSet,bkgRefSet,sigPdf,bkgPdf);
//   cTemplate->SetName( (((TString("pt")+=iPt)+="_eta")+=iEta)+="Template" );
  // plotDistortions(muonELoss,sigRefSet,bkgRefSet,sigPdf,bkgPdf,nPions,nPions.getVal(),2.0*nPions.getVal(),0.5*nPions.getVal(),"Number of pions");
  // plotDistortions(muonELoss,sigRefSet,bkgRefSet,sigConvPdf,bkgConvPdf,shift,0.0,-0.05,+0.05,"Shift");
  // plotDistortions(muonELoss,sigRefSet,bkgRefSet,sigConvPdf,bkgConvPdf,stretch,1.0,0.95,1.05,"Stretch");
  // plotDistortions(muonELoss,sigRefSet,bkgRefSet,sigConvPdf,bkgConvPdf,gaussSigma,1.0e-4,0.05,0.1,"Smear");

//   RooAbsPdf* prior  = (RooAbsPdf *) w->factory("Uniform::prior(s[0.5,0,1])");
//   RooArgSet argPOI(sigRatio);
//   RooArgSet nui(stretch,shift);  
// //   w->defineSet("poi","sigRatio");
// //   w->defineSet("nuis","shift,stretch");
//   ModelConfig modelConfig = ModelConfig("myConfig");
//   modelConfig.SetWorkspace(*w);
//   modelConfig.SetPdf(*w->pdf("mlpdf"));
//   modelConfig.SetPriorPdf(*w->pdf("prior"));
//   modelConfig.SetParameters(argPOI);
//   modelConfig.SetPdf(mdlPdf);
// //   modelConfig.SetNuisanceParameters(nui);
// //   modelConfig->SetParameters(*w->set("poi"));
// //   modelConfig->SetNuisanceParameters(*w->set("nuis"));
//   w->import(modelConfig);
//   w->import(*dataSubSet);

//   BayesianCalculator bc((RooAbsData&)*dataSubSet, (RooAbsPdf&)mdlPdf, (const RooArgSet&)RooArgSet(sigRatio), *priorPOI , &RooArgSet(shift,stretch));
//   bc.SetTestSize(0.32);
//   bc.SetPriorPdf(priorPOI);
 
  // --- NLL as function of signal ratio ---
  // --- ProfileLikelihoodCalculator usually make intervals: ---
//   ProfileLikelihoodCalculator plc(*dataSubSet,mdlPdf,sigRatio);
//   ProfileLikelihoodCalculator plc(*dataSubSet,modelConfig);
//   plc.SetTestSize(0.32);
  
  // --- Fit and retrieve the results ---
  RooFitResult* roofitOriginalResult = mdlPdf.fitTo(*dataSubSet, Minos(true), SumW2Error(false), Save(true));
  double sigRatioOriginal = sigRatio.getVal();
  double smearOriginal    = gaussSigma.getVal();
  double stretchOriginal  = stretch.getVal();
  double shiftOriginal    = shift.getVal();
  double errOriginal      = sigRatio.getError();
  arrSigRatio[0] = sigRatio.getVal();
  arrSmear[0]    = gaussSigma.getVal();
  arrStretch[0]  = stretch.getVal();
  arrShift[0]    = shift.getVal();
  arrErr[0]      = sigRatio.getError();
  std::cout << " value " << 0 << ": " << sigRatio.getVal() << " status = " << roofitOriginalResult->status() << std::endl;

 // --- Plot result ---
  double chi2 = 0.0;
//   TCanvas *c = new TCanvas("c","c",600,600);
//   c->SetName( ((TString("pt")+=iPt)+="_eta")+=iCentrality );
//   TCanvas *c = plotFitResult(variable, dataSubSet, mdlPdf, sigConvPdf, bkgConvPdf, chi2);
//   RooPlot* fitFrame = plotFitResultRooPlot(variable, dataSubSet, mdlPdf, sigShiftPdf, bkgShiftPdf, chi2);
//   variable.setBins(100);
//   RooPlot* fitFrame = variable.frame();
//   dataSubSet->plotOn(fitFrame);
//   variable.setBins(1000);

//start at i=1 since 0th element is already initialized above
  for (int i = 1; i <= nToyMC; ++i) {
    sigRatio.setVal(0.5);
//     variable.setVal(0.);
    RooDataSet* bkgtemp = bkgConvPdf.generate(variable,nBkg);
    RooDataSet* sigtemp = sigConvPdf.generate(variable,nSig);
//     RooDataSet* bkgtemp = bkgPdf.generate(variable,nBkg);
//     RooDataSet* sigtemp = sigPdf.generate(variable,nSig);
    std::cout << " bg :" <<  bkgtemp->numEntries() << " sig :"<< sigtemp->numEntries() << std::endl;
    RooKeysPdf bkgPdftemp("bkgPdftemp","bkgPdftemp",variable,*bkgtemp);
    RooKeysPdf sigPdftemp("sigPdftemp","sigPdftemp",variable,*sigtemp);
//     RooShiftedKeysPdf bkgShiftPdftemp("bkgShiftPdftemp","bkgShiftPdftemp",variable,bkgPdftemp,shift,stretch);
//     RooShiftedKeysPdf sigShiftPdftemp("sigShiftPdftemp","sigShiftPdftemp",variable,sigPdftemp,shift,stretch);
//     RooFFTConvPdf sigConvPdftemp("sigConvPdftemp","sigConvPdftemp",variable,sigShiftPdftemp,gauss);
//     RooFFTConvPdf bkgConvPdftemp("bkgConvPdftemp","bkgConvPdftemp",variable,bkgShiftPdftemp,gauss);
    RooAddPdf  mdlPdftemp("mdlPdftemp","mdlPdftemp",RooArgList(sigPdftemp,bkgPdftemp),RooArgList(sigRatio));
    RooFitResult* roofitResult = mdlPdftemp.fitTo(*dataSubSet, Minos(true), SumW2Error(false),PrintLevel(-1), Save(true));
    std::cout << " value " << i << ": " << sigRatio.getVal() << " status = " << roofitResult->status() << std::endl;
    arrSigRatio[i] = sigRatio.getVal();
    arrStretch[i]  = stretch.getVal();
    arrShift[i]    = shift.getVal();
    arrSmear[i]    = gaussSigma.getVal();
    arrErr[i]      = sigRatio.getError();
    std::cout << "Clean up." << std::endl;
    if (roofitResult) delete roofitResult;
    if (bkgtemp) delete bkgtemp;
    if (sigtemp) delete sigtemp;
    if(i==nToyMC){std::cout << "Last Entry of nToyMC Loop" << std::endl;}
  }
std::cout << "Outside nToyMC loop"<< std::endl;
//   mdlPdf.plotOn(fitFrame,LineColor(kBlack));
//   mdlPdf.plotOn(fitFrame,Components(sigConvPdf),LineColor(kBlue+1));
//   mdlPdf.plotOn(fitFrame,Components(bkgConvPdf),LineColor(kGreen-8));
//   chi2 = fitFrame->chiSquare(3);
//   fitResult.setChi2(iPt,iCentrality,chi2,0.0,0.0);
//   fitFrame->Draw();
//   TLegend *leg = new TLegend(0.3,0.7,0.6,0.925);
//   leg->SetFillColor(0);
//   leg->SetFillStyle(0);
//   leg->SetBorderSize(0);
//   leg->SetTextFont(gStyle->GetTextFont());
//   leg->SetTextSize(gStyle->GetTextSize());
//   TLegendEntry* entry = leg->AddEntry("","Data 2010 (#sqrt{#font[52]{s}} = 7 TeV)","p");
//   entry->SetMarkerStyle(20);
//   entry->SetMarkerSize(1.2);
//   entry = leg->AddEntry("","Best-fit template","l");
//   entry->SetLineColor(kBlack);
//   entry->SetLineWidth(2);
//   entry = leg->AddEntry("","Prompt component","l");
//   entry->SetLineColor(kBlue+1);
//   entry->SetLineWidth(2);
//   entry = leg->AddEntry("","Pion/kaon component","l");
//   entry->SetLineColor(kGreen-8);
//   entry->SetLineWidth(2);
//   leg->Draw();
//   fitResult.addCanvas(c);
//   if (pionWeight==1.0){
//     c->Print(TString(c->GetName())+=".png"); c->Print(TString(c->GetName())+=".eps");
//   }

  qsort(arrSigRatio,nToyMC+1,sizeof(float),compare); // using qsort for sorting or u can see bubble sort, selection sort...
//   getMedian(arrSigRatio,nToyMC+1);
  std::cout <<"qsort Checkpoint 1" << std::endl;
  int medianIndex = getMedianIndex(arrSigRatio, nToyMC+1);
  std::cout <<"MedianIndex Checkpoint 1" << std::endl;
  float rms = getRMS(arrSigRatio,nToyMC+1);
  std::cout << " median index = " << medianIndex << " " << arrSigRatio[medianIndex]<< "+-" << arrErr[medianIndex] << "+-" << rms << std::endl;

  std::cout <<"MedianIndex Checkpoint 2" << std::endl;

  std::cout  << sigRatio.getVal() << std::endl;
  sigRatio.setVal(sigRatioOriginal);
  std::cout  << sigRatio.getVal() << std::endl;
  stretch.setVal(stretchOriginal);
  std::cout  << "A" << std::endl;
  shift.setVal(shiftOriginal);
  std::cout  << "B" << std::endl;
  gaussSigma.setVal(smearOriginal);
  std::cout  << "C" << std::endl;

  fitResult.setShift(iPt,iEta,iCentrality,shift);
  std::cout  << "D" << std::endl;
  fitResult.setStretch(iPt,iEta,iCentrality,stretch);
  std::cout  << "E" << std::endl;
  fitResult.setSmear(iPt,iEta,iCentrality,gaussSigma);
  std::cout  << "F" << std::endl;
  double fitRatio = sigRatio.getVal();
  std::cout  << fitRatio << std::endl;

  // --- Calculate confidence band ---
  ProfileLikelihoodCalculator plc(*dataSubSet,mdlPdf,sigRatio);
  plc.SetTestSize(0.32);
//   plc.DoGlobalFit();

  double lowerLimit = 0.0;
  double upperLimit = 1.0;
  LikelihoodInterval* interval = plc.GetInterval();
//   double LLlowerLimit = 0.0;
//   double LLupperLimit = 1.0;
  bool limitsOK = interval->FindLimits(sigRatio,lowerLimit,upperLimit);
  if ((limitsOK == false) || (lowerLimit<1e-9)) {
    lowerLimit = sigRatioOriginal - errOriginal > 0. ? sigRatioOriginal - errOriginal : 0.;
  }
  if ((limitsOK == false) || (1.-upperLimit<1e-9)) {
    upperLimit = sigRatioOriginal + errOriginal < 1. ? sigRatioOriginal + errOriginal : 1.;
  }
  std::cout << " got interval " << std::endl;
  std::cout << " limits " << lowerLimit << " " << upperLimit  <<std::endl;
//   std::cout << " limits " << lowerLimit << " " << upperLimit  << " bc: " << bLowerLimit << " " << bUpperLimit  <<std::endl;
  // --- This is a hack (assuming symmetric uncertainties) ---
//   if ( lowerLimit == 0 ) {
//     std::cout << "WARNING: Hacking the lower limit" << std::endl;
//     lowerLimit = fitRatio-(upperLimit-fitRatio);
//     if (lowerLimit<0) lowerLimit = 0.;
//   }
//   if ( upperLimit == 1 ) {
//     std::cout << "WARNING: Hacking the upper limit" << std::endl;
//     upperLimit = fitRatio+(fitRatio-lowerLimit);
//     if (upperLimit>1.) upperLimit = 1.;
//   }

  fitResult.setSig(iPt, iEta, iCentrality,fitRatio,fitRatio-lowerLimit,upperLimit-fitRatio);

  /// set the rms from the toy monte carlo
  fitResult.setToyMc(iPt, iEta, iCentrality,fitRatio,rms,rms);

  /// --- Calculate MC prediction ---
  double mcRatio    = -1.;
  double mcRatioErr = 0.;
  if (refSubSet->numEntries()>0) {
    mcRatio = double(sigRefSet->numEntries())/double(refSubSet->numEntries());
    mcRatioErr = sqrt(double(sigRefSet->numEntries()))/double(refSubSet->numEntries());
  }
  fitResult.setMc(iPt, iEta, iCentrality,mcRatio,mcRatioErr,mcRatioErr);
  
  /// --- Get number of entries ---
  fitResult.setSigN(iPt,iEta, iCentrality,dataSubSet->numEntries(),0.0,0.0);
  fitResult.setMcN(iPt, iEta, iCentrality,refSubSet->numEntries(),0.0,0.0);
  
  // --- Get pt mean ---
  fitResult.setPt(iPt, iEta, iCentrality,dataSubSet->mean(muonPt));
  
  // --- Print summary ---
  std::cout << std::endl;
  std::cout << "-------------" << std::endl;
  std::cout << "S u m m a r y" << std::endl;
  std::cout << "-------------" << std::endl;
  std::cout << "MC events           = " << refSubSet->numEntries() << " ("<< sigRefSet->numEntries() << ")" << std::endl;
  std::cout << "Expected ratio      = " << mcRatio << std::endl;
  std::cout << "pt:eta:centrality   = " << iPt <<":"<< iEta <<":"<<iCentrality << std::endl;
  std::cout << "Events              = " << dataSubSet->numEntries() << std::endl;
  std::cout << "Fitted ratio        = " << fitRatio << std::endl;
  std::cout << "Fit Chi ^2          = " << chi2 << std::endl;
  std::cout << "MC RMS              = " << rms << std::endl;
  std::cout << "1 sigma CL interval = [" << lowerLimit << ", " << upperLimit << "]" << std::endl;

  if (roofitOriginalResult) delete roofitOriginalResult;
//   if (parSigVars) delete parSigVars;
//   if (parBkgVars) delete parBkgVars; 
  if (sigRefSet) delete sigRefSet;
  if (bkgRefSet) delete bkgRefSet;
  if (sigMCSet) delete sigMCSet;
  if (bkgMCSet) delete bkgMCSet;
  if (interval)  delete interval;
  if (w)         delete w;
  if (arrSigRatio)       delete arrSigRatio;
  if (arrStretch)        delete arrStretch;
  if (arrShift)          delete arrShift;
  if (arrSmear)          delete arrSmear;
  if (arrErr)            delete arrErr;
}

TFile* getTemplateFile(double ptLow, double ptUpp, double etaLow, double etaUpp, const char* source, bool doReweightPtSpectrum, bool doMirrorEta, RooRealVar& var){
  TString fileName = getTemplateFileName(ptLow, ptUpp, etaLow, etaUpp, source, doReweightPtSpectrum, doMirrorEta, var );
  // --- Read templates from file ---
  std::cout << "Reading the templates from file " << fileName << "..." << std::endl;
  TFile* fPdf = new TFile(fileName, "READ");
  if ( fPdf == 0 ) {
    std::cout << "Template file not found" << std::endl;
  }
  return fPdf;
}

