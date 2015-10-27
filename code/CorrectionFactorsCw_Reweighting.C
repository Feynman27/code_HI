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
#include "TPad.h"
#include "TProfile.h"
#include "TVector3.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooAbsArg.h"
#include "RooDataSet.h"
#include "RooBinning.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooCategory.h"

#include "AtlasUtils.C"
#include "AtlasStyle.C"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath> 
#include <sstream>

//#include "Systematics.C"
//#include "WPlotterHelper.C"
#include "TriggerEfficiencies.C"
#include "correctionFactorsDep.C"


using namespace RooFit;

///////////////////////////////////////////
//writeToSpreadsheet
//////////////////////////////////////////
void writeToSpreadsheet(std::ostream& outputFile, int ipt, int ieta, int icent, double scaleFactor, double errStatLo,double errStatUp){
	outputFile << ipt << "," << ieta << "," << icent << "," << scaleFactor << "," << errStatLo << "," << errStatUp << std::endl;
}

float getMt(float pt, float eta, float phi, float mptPhi, float mpt){

    TVector3 lvMu;
    lvMu.SetPtEtaPhi(pt,eta,phi) ;
    float phiMuTemp = lvMu.Phi();
    if(phiMuTemp!=phi) {
        std::cout << "ERROR: inconsistency in phi." << std::endl;
        exit(0);
    }
    float phi_munu = phiMuTemp-mptPhi;
    if(phi_munu> TMath::Pi()) {phi_munu -= TMath::TwoPi();}
    if(phi_munu<-1*TMath::Pi()) {phi_munu+=TMath::TwoPi();}
    float mt = (fabs(mpt) < 9999.) ? TMath::Sqrt(2.*fabs(pt)*fabs(mpt)*(1.0-TMath::Cos(phi_munu))) : -9999. ;
    return mt;
}

///////////////////////////////
//used for mapping negative and positive eta bins
//to absolute eta bin values
//////////////////////////////
int indexPositiveEta(int ieta, int nEtaBins){
      int index = ieta-1.0/2.0*(nEtaBins+1.0);
      std::cout << "indexPositiveEta:" << index << std::endl;
      if(index>(nEtaBins-1)/2) std::cout << "WARNING: index out of allowed range. EXPECT WRONG RESULTS."<<std::endl;
      return index;
}
int indexNegativeEta(int ieta, int nEtaBins){
     int index = 0.5*(nEtaBins+1.0)-2.0; index-=ieta;
     std::cout << "indexNegativeEta:" << index << std::endl;
     if(index>(nEtaBins-1)/2) std::cout << "WARNING: index out of allowed range. EXPECT WRONG RESULTS."<<std::endl;
     return index;
}
/////////////////////////////////////////////
//plot residual between fit and Cw(Npart) distribution
//per eta bin
/////////////////////////////////////////////
void plotCwResidual(TGraphAsymmErrors* gr, TGraph* const grResid){

	std::cout << "Plotting residual...." << std::endl;
	//double xPt = etaLo + 0.5*fabs(etaUpp-etaLo);
	TF1* f = gr->GetFunction("f0a");

	double* xPt = gr->GetX();
	double* yPt = gr->GetY();

	for(int igr=0; igr<gr->GetN(); igr++){

		std::cout << "x,y = " << xPt[igr] << "," << yPt[igr] << std::endl;
		double resid = fabs(yPt[igr] - f->Eval(xPt[igr])); 
		grResid->SetPoint(igr,xPt[igr],resid);
	}
}



///////////////////////////////////////////////////////////////////////////////
//plot distro of generator level kinematics 
///////////////////////////////////////////////////////////////////////////////
void plotGeneratorDistros(TH1F* hGen,TH1F* hFid,int ich,TString sSel){

	std::cout << "Now plotting generator level distros..." << std::endl;

	TH1F* hGenc = (TH1F*)hGen->Clone("hGenc");
	TH1F* hFidc = (TH1F*)hFid->Clone("hFidc");

	hGenc->SetFillColor(kYellow);
	hFidc->SetFillColor(kRed);

	//TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
	TLegend* leg = new TLegend(0.619, 0.559, 0.876, 0.939);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);

	leg->AddEntry(hGenc, "gen,all", "f");
	leg->AddEntry(hFidc, "gen,cut", "f");

  	TCanvas* cGen = new TCanvas("cGen","cGen",600,600);
	hGenc->GetXaxis()->SetTitle("#eta_{#mu,gen}[GeV]");
	hGenc->GetYaxis()->SetTitle("Muons/GeV");
	hGenc->GetYaxis()->SetRangeUser(0.0,1.0e7);
	hGenc->GetXaxis()->SetRangeUser(-6.0,6.0);

	hGenc->Scale(1,"width");
	hFidc->Scale(1,"width");

	hGenc->Draw("hist f");
	hFidc->Draw("hist fsame");
	hGenc->Draw("sameaxis");

//	myText(0.196,0.862,kBlack,(char*)sSel);

	leg->Draw();

        //cGen->SetLogy(true); hGenc->GetYaxis()->SetRangeUser(0.1,hGenc->GetMaximum()*7.0e2); cGen->Update();
	cGen->Print(sSel+",GeneratorDistro.pdf");
    TString name = "GeneratorDistro"; name+=ich;
	cGen->Print(name+".root");

}

///////////////////////////////////////////////////////////////////////////////
//plot distro of fiducial kinematic 
///////////////////////////////////////////////////////////////////////////////
void plotFiducialDistros(TH1F* hFid,TH1F* hRec,TString sSel){

	std::cout << "Now plotting fiducial distros..." << std::endl;

	TH1F* hFidc = (TH1F*)hFid->Clone("hFidc");
	TH1F* hRecc = (TH1F*)hRec->Clone("hRecc");

	hFidc->SetFillColor(kRed);
	hRecc->SetFillColor(kBlue);

//	TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
	TLegend* leg = new TLegend(0.619, 0.559, 0.876, 0.939);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);

	leg->AddEntry(hFidc, "gen,cut", "f");
	leg->AddEntry(hRecc, "rec,cut", "f");

  	TCanvas* cFid = new TCanvas("cFid","cFid",600,600);
	hFidc->GetXaxis()->SetTitle("p_{T}^{#mu}[GeV]");
	hFidc->GetYaxis()->SetTitle("Muons/GeV");
	hFidc->GetYaxis()->SetRangeUser(0.1,1.0e5);
	hFidc->GetXaxis()->SetRangeUser(0.0,250.0);

	hFidc->Scale(1,"width");
	hRecc->Scale(1,"width");

	hFidc->Draw("hist f");
	hRecc->Draw("hist fsame");
	hFidc->Draw("sameaxis");

//	myText(0.196,0.862,kBlack,sSel);

	leg->Draw();

        cFid->SetLogy(true); hFidc->GetYaxis()->SetRangeUser(0.1,hFidc->GetMaximum()*7.0e2); cFid->Update();
	cFid->Print(sSel+",FiducialDistro.pdf");

}


///////////////////////////////////////////////////////////////////////////////
//plot Cw as a function of centrality 
///////////////////////////////////////////////////////////////////////////////
void plotCwCentralityDep(TGraphErrors* const grCent,int icent, double npart,float Cw, float CwErr){

	std::cout << "Plotting Cw centrality dependance at centrality bin: " << icent << ": Npart " << npart << std::endl;
	grCent->SetPoint(icent,npart,Cw);
	grCent->SetPointError(icent,0.0,CwErr);

}

///////////////////////////////////////////////////////////////////////////////
//plot Cw as a function of eta 
///////////////////////////////////////////////////////////////////////////////
void plotCwEtaDep(TGraphErrors* const grEta,int ieta, float xEta, float Cw, float CwErr){

	std::cout << "Plotting Cw eta dependance at eta bin: " << ieta << std::endl;
	grEta->SetPoint(ieta,xEta,Cw);
	grEta->SetPointError(ieta,0.0,CwErr);

}
///////////////////////////////////////////////////////////////////////////////
//index ipt,ieta,icent 
///////////////////////////////////////////////////////////////////////////////
int indexIJK(int ich, int i, int j, int k, int nCentrality, int nEta){

	int increment = 0;	
	if(ich==102) increment = 0;
	else if(ich==103) increment = 1;
	else if(ich==104) increment = 2;
	else {
		std::cout << "WARNING: Unable to find charge identifier. Will return -9999." << std::endl;
		return -9999;
	}
	int index =  2*(nCentrality*nEta*i + nCentrality*j + k) + increment ;
	return index;
}

void calcWeightedEfficiency(int i, double x, double binW, int index, int ipt, int ieta, int icent,
			         const TGraphAsymmErrors* grppCw, const TGraphAsymmErrors* grnpCw,
				 const TGraphAsymmErrors* grpnCw, const TGraphAsymmErrors* grnnCw,
                                 TGraphAsymmErrors* grWtdCw){

	double cWpp = grppCw->GetY()[i]; double errcWppUpp = grppCw->GetEYhigh()[i]; double errcWppLow = grppCw->GetEYlow()[i];
	double cWnp = grnpCw->GetY()[i]; double errcWnpUpp = grnpCw->GetEYhigh()[i]; double errcWnpLow = grnpCw->GetEYlow()[i];
	double cWpn = grpnCw->GetY()[i]; double errcWpnUpp = grpnCw->GetEYhigh()[i]; double errcWpnLow = grpnCw->GetEYlow()[i];
	double cWnn = grnnCw->GetY()[i]; double errcWnnUpp = grnnCw->GetEYhigh()[i]; double errcWnnLow = grnnCw->GetEYlow()[i];

	double wtdCw = 0.155*cWpp+0.478*((cWnp+cWpn)/2.0)+0.367*cWnn;
	double wtderrCwUpp = 0.155*errcWppUpp+0.478*((errcWnpUpp+errcWpnUpp)/2.0)+0.367*errcWnnUpp;
	double wtderrCwLow = 0.155*errcWppLow+0.478*((errcWnpLow+errcWpnLow)/2.0)+0.367*errcWnnLow;

	// Fill TGraph of weighted Cw efficiency
	grWtdCw->SetPoint(i,x,wtdCw);
    grWtdCw->SetPointError(i,binW/2.0,binW/2.0,wtderrCwLow,wtderrCwUpp);

}

void calcEfficiency(const TGraphAsymmErrors* const grAwPlus, const TGraphAsymmErrors* const grAwMinus,
                    const TGraphAsymmErrors* const grCwPlus, const TGraphAsymmErrors* const grCwMinus,
                    const TGraphAsymmErrors* const grCwAwPlus, const TGraphAsymmErrors* const grCwAwMinus,
                    double csPlus, double csMinus,
                    int i, double x, double binW, int index, int ipt, int ieta, int icent, 
                    TGraphAsymmErrors* grAw, TGraphAsymmErrors* grCwAw, TGraphAsymmErrors* grCw
                ){

        double aw = (csPlus*grAwPlus->GetY()[i] + csMinus*grAwMinus->GetY()[i])/(csPlus+csMinus); 
        double errAwHi = (csPlus*grAwPlus->GetEYhigh()[i] + csMinus*grAwMinus->GetEYhigh()[i])/(csPlus+csMinus); 
        double errAwLo = (csPlus*grAwPlus->GetEYlow()[i] + csMinus*grAwMinus->GetEYlow()[i])/(csPlus+csMinus) ;
        grAw->SetPoint(i,x,aw);
        grAw->SetPointError(i,binW/2.0,binW/2.0,errAwLo,errAwHi);

        double cw = (csPlus*grCwPlus->GetY()[i] + csMinus*grCwMinus->GetY()[i])/(csPlus+csMinus) ;
        double errCwHi = (csPlus*grCwPlus->GetEYhigh()[i] + csMinus*grCwMinus->GetEYhigh()[i])/(csPlus+csMinus); 
        double errCwLo = (csPlus*grCwPlus->GetEYlow()[i] + csMinus*grCwMinus->GetEYlow()[i])/(csPlus+csMinus) ;
        grCw->SetPoint(i,x,cw);
        grCw->SetPointError(i,binW/2.0,binW/2.0,errCwLo,errCwHi);

        double cwaw = (csPlus*grCwAwPlus->GetY()[i] + csMinus*grCwAwMinus->GetY()[i])/(csPlus+csMinus) ;
        double errCwAwHi = (csPlus*grCwAwPlus->GetEYhigh()[i] + csMinus*grCwAwMinus->GetEYhigh()[i])/(csPlus+csMinus); 
        double errCwAwLo = (csPlus*grCwAwPlus->GetEYlow()[i] + csMinus*grCwAwMinus->GetEYlow()[i])/(csPlus+csMinus) ;
        grCwAw->SetPoint(i,x,cwaw);
        grCwAw->SetPointError(i,binW/2.0,binW/2.0,errCwAwLo,errCwAwHi);


}
///////////////////////////////////////////////////////////////////////////////
//Function for calculating and plotting efficiencies to a TGraph
///////////////////////////////////////////////////////////////////////////////
void calcEfficiency(RooDataSet* mcWGenSet, RooDataSet* mcWFidSet,RooDataSet* mcWRecSet,
                int i, double x, double binW, int index, int ipt, int ieta, int icent, 
                TGraphAsymmErrors* grAw, TGraphAsymmErrors* grCwAw, TGraphAsymmErrors* grCw,
                std::ostream& output, bool doWriteToSpreadsheet=false, bool doSingleTriggerEff=false){

    //number of Wmu in this eta and centrality class
    int tGen = (int)mcWGenSet->sumEntries();
    std::cout << "tGen = " << tGen << std::endl;

    //number of muons at generator level in kinematic fiducial region within this eta and centrality class
    int tGenInFiducial = (int)mcWFidSet->sumEntries();
    std::cout << "tGenInFiducial = " << tGenInFiducial << std::endl;

    double aw, errAwHi, errAwLo;
    Efficiency(tGenInFiducial,tGen,0.683,aw,errAwLo,errAwHi);
    if(aw<0.0) Warning("calcEfficiency","Efficiency  = %d. Expect wrong results", aw);
    if(aw>1.0) Warning("calcEfficiency","Efficiency  = %d. Expect wrong results", aw);
   
    errAwHi = errAwHi - aw;
    errAwLo = aw - errAwLo;
    ///Aw 

    std::cout << "Efficiency Aw: " << aw << " + " << errAwHi << " - " << errAwLo  << std::endl;
    grAw->SetPoint(i,x,aw);
    grAw->SetPointError(i,binW/2.0,binW/2.0,errAwLo,errAwHi);

    ///Retrieve triggere efficiency for this bin
    double trigTemp,trigStatErr;
    ///Use Eta/Centrality averaged trigger
    ///efficiency since not enough statistics
    ///for binning in pt/eta/centrality/charge
    if(doSingleTriggerEff) {
        trigTemp=0.959;
        trigStatErr=0.016;
    }
    else{
        trigTemp = trigEfficiency(index);
        trigStatErr = trigEfficiencyErr(index);
        
        //trigTemp = 1.0;
        //trigStatErr = 0.0;
    }

    ///RecHiQualityWsel
    int pRec = (int)mcWRecSet->sumEntries() ;
    std::cout << "pRec = " << pRec << std::endl;

    ///Calculate efficiency using Bayesian errors 
    ///with a 68.3% C.I.
    double recEff, errRecHi, errRecLo;
    Efficiency(pRec,tGenInFiducial,0.683,recEff,errRecLo,errRecHi);
    if(recEff<0.0) Warning("calcEfficiency","Efficiency %d. Expect wrong results", recEff);
    if(recEff>1.0) Warning("calcEfficiency","Efficiency %d. Expect wrong results", recEff);

    errRecHi = errRecHi - recEff;
    errRecLo = recEff - errRecLo;

    ///Apply trigger efficiency to reconstructed part
    std::cout << "Number of reconstructed muons in fiducial cut space = " << pRec << std::endl;
    std::cout << "Number of generated muons in fiducial cut space = " << tGenInFiducial << std::endl;
    std::cout << "Consistency Check: " << (double)pRec/tGenInFiducial << "=?" << recEff << std::endl;
    std::cout << "Cw before trigger correction for bin = " << ieta << ":" << icent << " = "  << recEff << std::endl;
    std::cout << "Trigger efficiency in bin = " << ieta << ":" << icent << " = "  << trigTemp << std::endl;
    double eff = recEff*trigTemp;
    std::cout << "Cw after trigger correction for bin = " << ieta << ":" << icent << " = "  << eff << std::endl;

    ///Propagate the error
    double effStatErrHi = sqrt( TMath::Power(errRecHi/recEff,2) + TMath::Power(trigStatErr,2) )*eff;
    double effStatErrLo = sqrt( TMath::Power(errRecLo/recEff,2) + TMath::Power(trigStatErr,2) )*eff;

    std::cout << "Cw Efficiency: " << recEff << " + " << effStatErrHi << " - " << effStatErrLo  << std::endl;

    grCw->SetPoint(i,x,eff);
    grCw->SetPointError(i,binW/2.0,binW/2.0,effStatErrLo,effStatErrHi);

    ///Calculate efficiency using Bayesian errors 
    ///with a 68.3% C.I.
    double CwAw, errCwAwHi, errCwAwLo;
    Efficiency(pRec,tGen,0.683,CwAw,errCwAwLo,errCwAwHi);
    if(CwAw<0.0) Warning("calcEfficiency","Efficiency  = %d. Expect wrong results", CwAw);
    if(CwAw>1.0) Warning("calcEfficiency","Efficiency  = %d. Expect wrong results", CwAw);

    ///Return +- value from mean; therefore; must subtract
    ///to get absolute errors
    errCwAwHi = errCwAwHi - CwAw;
    errCwAwLo = CwAw - errCwAwLo;
    ///Apply trigger efficiency to reconstructed part
    CwAw*=trigTemp;

    ///Propagate the error
    double effStatErrHiCwAw = sqrt( TMath::Power(errCwAwHi/CwAw,2) + TMath::Power(trigStatErr,2) )*CwAw;
    double effStatErrLoCwAw = sqrt( TMath::Power(errCwAwLo/CwAw,2) + TMath::Power(trigStatErr,2) )*CwAw;

    std::cout << "EfficiencyCwAw: " << CwAw << " + " << effStatErrHiCwAw << " - " << effStatErrLoCwAw  << std::endl;

    std::cout << "Point " << i << ":" << "x=" << x << " y=" << CwAw << std::endl;
    grCwAw->SetPoint(i,x,CwAw);
    grCwAw->SetPointError(i,binW/2.0,binW/2.0,effStatErrLoCwAw,effStatErrHiCwAw);


    if(doWriteToSpreadsheet) writeToSpreadsheet(output,ipt,ieta,icent,eff,effStatErrLo,effStatErrHi);

}//calcEfficiency


///////////////////////////////////////////////////////////////////////////////
// selectGenPtEta
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectGenPtEta( RooDataSet* dataSet, double ptLow, double ptUpp, double etaLow, double etaUpp, bool doMirrorEta = false )
{
  TString cut = "muonGenPt>";
  cut += ptLow;
  cut += "&&muonGenPt<";
  cut += ptUpp;
  cut += "&&((muEtaGen>=";
  cut += etaLow;
  cut += "&&muEtaGen<";
  cut += etaUpp;
  cut += ")";
  if ( doMirrorEta ) {
    cut += "||(muEtaGen>=";
    cut += -etaUpp;
    cut += "&&muEtaGen<";
    cut += -etaLow;
    cut += ")";
  }
  cut += ")";
  ///cut on neutrino eta
/*  cut += "&&((nuEtaGen>";
  cut += etaLow;
  cut += "&&nuEtaGen<";
  cut += etaUpp;
  cut += ")";
  if ( doMirrorEta ) {
    cut += "||(nuEtaGen>";
    cut += -etaUpp;
    cut += "&&nuEtaGen<";
    cut += -etaLow;
    cut += ")";
  }
  cut += ")";
*/
  dataSet = (RooDataSet*) dataSet->reduce(cut);
  return dataSet;
}
///////////////////////////////////////////////////////////////////////////////
// selectMptMt
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectMptMt( RooDataSet* dataSet,double mptLow, TProfile* _pfx)
{
  std::cout << "Cutting on MPT > " << mptLow << std::endl;
  ///Find the bin corresponding to this mpt cut
  int bin = _pfx->FindBin(mptLow);
  double mtCut = _pfx->GetBinContent(bin);
  std::cout << "This corresponds to and mT cut > " << mtCut << std::endl;
  TString cut = "missPt>";
  cut += mptLow;
  cut += "&&muonMt>";
  cut += mtCut;
  std::cout << "Reducing d.s. with " << cut << std::endl;
  dataSet = (RooDataSet*) dataSet->reduce(cut);
  return dataSet;
}


///////////////////////////////////////////////////////////////////////////////
// selectMpt
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectMpt( RooDataSet* dataSet,double mptLow)
{
  std::cout << "Cutting on MPT." << std::endl;
  TString cut = "missPt>";
  cut += mptLow;
  std::cout << "Reducing d.s. with " << cut << std::endl;
  dataSet = (RooDataSet*) dataSet->reduce(cut);
  return dataSet;
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
  cut += "&&((muonEta>=";
  cut += etaLow;
  cut += "&&muonEta<";
  cut += etaUpp;
  cut += ")";
  if ( doMirrorEta ) {
    cut += "||(muonEta>=";
    cut += -etaUpp;
    cut += "&&muonEta<";
    cut += -etaLow;
    cut += ")";
  }
  cut += ")";

  dataSet = (RooDataSet*) dataSet->reduce(cut);
  return dataSet;
}
///////////////////////////////////////////////////////////////////////////////
//selectCentrality 
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectCentrality( RooDataSet* dataSet,double centralityLow, double centralityUpp){

  //cut on generated eta if this is a generator level cut
  TString cut = "centrality>";
  cut += centralityLow;
  cut += "&&centrality<";
  cut+= centralityUpp;

  dataSet = (RooDataSet*)dataSet->reduce(cut);
  return dataSet;
}
///////////////////////////////////////////////////////////////////////////////
//selectPtEtaCentrality 
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectPtEtaCentrality( RooDataSet* dataSet,double ptLow, double ptUpp, double etaLow, double etaUpp, double centralityLow, double centralityUpp, 
	bool doMirrorEta=false, bool isGen=false){

  //cut on generated eta if this is a generator level cut
  if(isGen) dataSet= selectGenPtEta( dataSet,ptLow,ptUpp, etaLow, etaUpp, doMirrorEta);
  else dataSet = selectPtEta( dataSet,ptLow,ptUpp, etaLow, etaUpp, doMirrorEta);
  TString cut = "centrality>";
  cut += centralityLow;
  cut += "&&centrality<";
  cut+= centralityUpp;

  dataSet = (RooDataSet*)dataSet->reduce(cut);
  return dataSet;
}

double getWeight(const TH1F* h, double etaGen){

    if(etaGen<0.0){
        std::cout << "WARNING: Parametrization only calculated for absolute eta. Use fabs(etaGen) as argument in getWeight()." << std::endl;
        exit(0);
    }
    TF1* f = h->GetFunction("pol3");
    double weight = 1.0/f->Eval(etaGen);
    return weight;
}

///////////////////////////////////////////////////////////////////////////////
// fillHIMuonRecSet
///////////////////////////////////////////////////////////////////////////////
RooDataSet* fillHIMuonRecSet(const TString& pathName, const TString& fileName, RooArgSet& muonArgSet, const TH1F* h, bool doMptSigmaUp = false, bool doMptSigmaDown=false)
{
  RooDataSet* set = new RooDataSet(fileName,fileName,muonArgSet,"weight");
  
  float eLossNt[50];
  float scatNt[50];
  float compNt[50];
  float ptNt[50];
  float mtNt[50];
  float etaNt[50],phiNt[50];
  float etaGen[50];
  float chargeNt[50];
  float centralityNt;
  float nu_ptNt, nu_phiNt;
  float nu_phiNominal;
  float ptconeNt[50];
  int valNt[50], ZDYNt[50], truthMatchedNt[50], promptNt[50];
  int nmu;

  // --- Load the MC tree ---
  TChain* tree = new TChain("tree","tree");
  int nFiles = tree->Add(pathName+fileName);
  
  std::cout << "Filling the reco level RooDataSet for " << fileName << "... Number of files: " << nFiles << std::endl;
 
  // --- Set branch adresses ---
  tree->SetBranchAddress("eLoss", &eLossNt);
  //dR<0.3,pTid>3GeV
  tree->SetBranchAddress("ptcone20ID3", &ptconeNt);
  //vary ptcone variable for systematics
//  tree->SetBranchAddress("ptcone30ID3", &ptconeNt);
  tree->SetBranchAddress("scat", &scatNt);
  tree->SetBranchAddress("comp", &compNt);
  tree->SetBranchAddress("pt", &ptNt);
  tree->SetBranchAddress("eta", &etaNt);
  tree->SetBranchAddress("eta_muid_tru", &etaGen);
  tree->SetBranchAddress("phi", &phiNt);
  tree->SetBranchAddress("charge", &chargeNt);
  tree->SetBranchAddress("prompt", &promptNt);
  tree->SetBranchAddress("val", &valNt); 
  tree->SetBranchAddress("ZDY", &ZDYNt); 
  tree->SetBranchAddress("centrality", &centralityNt);
  if(doMptSigmaDown) {
    std::cout << "WARNING: Using 2GeV track pT threshold relative to the nominal value for mpt." << std::endl;
    tree->SetBranchAddress("nu_pt2000Nominal", &nu_ptNt);
    tree->SetBranchAddress("nu_phi2000Nominal", &nu_phiNt);
  }
  else if(doMptSigmaUp){
    //std::cout << "WARNING: Using smeared mt for systematic study." << std::endl;
    std::cout << "WARNING: Using 4GeV track pT threshold relative to the nominal value for mpt." << std::endl;
    tree->SetBranchAddress("nu_pt4000Nominal", &nu_ptNt);
    tree->SetBranchAddress("nu_phi4000Nominal", &nu_phiNt);
  }
  else {
      std::cout << "Using nominal mpt." << std::endl;
      tree->SetBranchAddress("nu_pt", &nu_ptNt);
  }
  tree->SetBranchAddress("nu_phi", &nu_phiNominal);
  tree->SetBranchAddress("mt", &mtNt);
  tree->SetBranchAddress("truthMatched_muid", &truthMatchedNt);
  tree->SetBranchAddress("mu_muid_n", &nmu);

   // --- Set branch status ---
  tree->SetBranchStatus("*",0) ;
  tree->SetBranchStatus("mu_muid_n", 1);
  tree->SetBranchStatus("truthMatched_muid", 1);
  tree->SetBranchStatus("eLoss", 1);
  tree->SetBranchStatus("ptcone20ID3", 1);
//  tree->SetBranchStatus("ptcone30ID3", 1);
  tree->SetBranchStatus("scat", 1);
  tree->SetBranchStatus("comp", 1);
  tree->SetBranchStatus("pt", 1);
  tree->SetBranchStatus("mt", 1);
  tree->SetBranchStatus("eta", 1);
  tree->SetBranchStatus("eta_muid_tru", 1);
  tree->SetBranchStatus("charge", 1);
  tree->SetBranchStatus("val", 1); 
  tree->SetBranchStatus("ZDY", 1); 
  tree->SetBranchStatus("centrality", 1);
  tree->SetBranchStatus("nu_pt", 1);
  tree->SetBranchStatus("nu_phi", 1);
  tree->SetBranchStatus("phi", 1);
  /*if(doMptSigmaUp) {
      tree->SetBranchStatus("nu_ptMod", 1);
      tree->SetBranchStatus("mtMod", 1);
  }
  */
  if(doMptSigmaDown) {
      tree->SetBranchStatus("nu_pt2000Nominal", 1);
      tree->SetBranchStatus("nu_phi2000Nominal", 1);
  }
  if(doMptSigmaUp) {
      tree->SetBranchStatus("nu_pt4000Nominal", 1);
      tree->SetBranchStatus("nu_phi4000Nominal", 1);
  }
  
  tree->SetBranchStatus("prompt", 1);

  std::cout << "Number of entries: " << tree->GetEntries() << std::endl; 
  for ( int i = 0; i < tree->GetEntries(); i++ ) {
//     if (i>100000) return set; // temp hack
    tree->LoadTree(i);
    tree->GetEntry(i);

    muonArgSet.setRealValue("missPt",nu_ptNt);
    muonArgSet.setRealValue("centrality",centralityNt);

    for (int imu = 0; imu<nmu;imu++){

      if ( promptNt[imu] == 23) muonArgSet.setCatLabel("muonCategory","Z");
      else if ( promptNt[imu] == 24) muonArgSet.setCatLabel("muonCategory","W");

      muonArgSet.setRealValue("muonELoss",eLossNt[imu]);
      muonArgSet.setRealValue("muonScat",scatNt[imu]);
      muonArgSet.setRealValue("muonQuality",valNt[imu]);

      double isolationTemp = -9999.0;
      isolationTemp = ptconeNt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolation",isolationTemp);
      muonArgSet.setRealValue("muonPt",ptNt[imu]);
      if(doMptSigmaUp||doMptSigmaDown){
          double mtTemp = -9999.0;
          mtTemp = getMt(ptNt[imu],etaNt[imu],phiNt[imu],nu_phiNt,nu_ptNt);
          muonArgSet.setRealValue("muonMt",mtTemp);
      }
      ///Nominal mt
      else {
          muonArgSet.setRealValue("muonMt",mtNt[imu]);
      }
      muonArgSet.setRealValue("muonEta",etaNt[imu]);
      muonArgSet.setRealValue("ZDY",ZDYNt[imu]);
      muonArgSet.setRealValue("muonCharge",chargeNt[imu]);
      //muonArgSet.setRealValue("muonCategory",promptNt[imu]);
      muonArgSet.setRealValue("motherRec",promptNt[imu]);
      muonArgSet.setRealValue("muonGenRecMatched",truthMatchedNt[imu]);
      if ( chargeNt[imu] > 0 ) muonArgSet.setCatLabel("chargeCategory","muPlus");
      else if ( chargeNt[imu] < 0) muonArgSet.setCatLabel("chargeCategory","muMinus");
      double w = getWeight(h,fabs(etaGen[imu]));
      set->add(muonArgSet,w,0.0);    
   }
  }

  return set;
}


RooDataSet* fillHIMuonGenSet(const TString& pathName, const TString& fileName, RooArgSet& muonGenSet,const TH1F* h)
{
  RooDataSet* set = new RooDataSet(fileName,fileName,muonGenSet,"weight");
  
  float nuGenPtNt[50], nuPhiGenNt[50];
  float muGenPtNt[50], muPhiGenNt[50];
  float mtGenNt[50];
  int motherNt[50];
  int daughterNt[50];
  int daughterNuNt[50];
  float muEtaGenNt[50];
  float nuEtaGenNt[50];
  float centralityNt;
  float chargeGenNt[50];
  int nmu;

  // --- Load the MC tree ---
  TChain* tree = new TChain("tree","tree");
  int nFiles = tree->Add(pathName+fileName);
  
  std::cout << "Filling the RooDataSet for " << fileName << "... Number of files: " << nFiles << std::endl;
 
  // --- Set branch adresses ---
  tree->SetBranchAddress("mc_mu_gen_pt", &muGenPtNt);
  tree->SetBranchAddress("mc_nu_gen_pt", &nuGenPtNt);
  tree->SetBranchAddress("mc_mu_gen_mothertype", &motherNt);
  tree->SetBranchAddress("mc_mu_gen_type", &daughterNt);
  tree->SetBranchAddress("mc_nu_gen_type", &daughterNuNt);
  tree->SetBranchAddress("mc_mu_charge", &chargeGenNt);
  tree->SetBranchAddress("mc_mu_gen_eta", &muEtaGenNt);
  tree->SetBranchAddress("mc_nu_gen_eta", &nuEtaGenNt);
  tree->SetBranchAddress("mc_mu_gen_phi", &muPhiGenNt);
  tree->SetBranchAddress("mc_nu_gen_phi", &nuPhiGenNt);
  tree->SetBranchAddress("mc_mu_n", &nmu);
  tree->SetBranchAddress("centrality", &centralityNt);

   // --- Set branch status ---
  tree->SetBranchStatus("*",0) ;
  tree->SetBranchStatus("mc_mu_n", 1);
  tree->SetBranchStatus("centrality", 1);
  tree->SetBranchStatus("mc_mu_gen_pt", 1);
  tree->SetBranchStatus("mc_nu_gen_pt", 1);
  tree->SetBranchStatus("mc_mu_gen_mothertype", 1);
  tree->SetBranchStatus("mc_mu_gen_type", 1);
  tree->SetBranchStatus("mc_nu_gen_type", 1);
  tree->SetBranchStatus("mc_mu_charge", 1);
  tree->SetBranchStatus("mc_mu_gen_eta", 1);
  tree->SetBranchStatus("mc_nu_gen_eta", 1);
  tree->SetBranchStatus("mc_mu_gen_phi", 1);
  tree->SetBranchStatus("mc_nu_gen_phi", 1);

  std::cout << "Number of entries: " << tree->GetEntries() << std::endl; 
  for ( int i = 0; i < tree->GetEntries(); i++ ) {
//     if (i>100000) return set; // temp hack
    tree->LoadTree(i);
    tree->GetEntry(i);

    muonGenSet.setRealValue("centrality",centralityNt);

    for (int imu = 0; imu<nmu;imu++){

      float dPhi = nuPhiGenNt[imu]-muPhiGenNt[imu];
      if(dPhi<-1*TMath::Pi()) dPhi += TMath::TwoPi(); if(dPhi>TMath::Pi()) dPhi -= TMath::TwoPi(); //fold btwn [-pi,pi]
      mtGenNt[imu] = TMath::Sqrt(2.0*muGenPtNt[imu]*nuGenPtNt[imu]*(1.0-TMath::Cos(dPhi)));
      muonGenSet.setRealValue("munuGenMt",mtGenNt[imu]);

      muonGenSet.setRealValue("muonGenPt",muGenPtNt[imu]);
      muonGenSet.setRealValue("nuGenPt",nuGenPtNt[imu]);
      muonGenSet.setRealValue("mother",motherNt[imu]);
      muonGenSet.setRealValue("daughter",daughterNt[imu]);
      ///Equivalent number of neutrinos as muons from Ws
      muonGenSet.setRealValue("daughterNu",daughterNuNt[imu]);
      muonGenSet.setRealValue("chargeGen",chargeGenNt[imu]);
      muonGenSet.setRealValue("muEtaGen",muEtaGenNt[imu]);
      muonGenSet.setRealValue("nuEtaGen",nuEtaGenNt[imu]);

      if ( chargeGenNt[imu] > 0. ) muonGenSet.setCatLabel("chargeGenCategory","muPlus");
      else if ( chargeGenNt[imu] < 0.) muonGenSet.setCatLabel("chargeGenCategory","muMinus");
      double w = getWeight(h,fabs(muEtaGenNt[imu]));
      set->add(muonGenSet,w,0.0);    
   }
  }

  return set;
}

void Write(TFile* const outFile, TObject* const gr, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing TGraph to root file..." << std::endl;
	  gr->Write(sFile);
	  std::cout << "Done." << std::endl;
    	  gDirectory = dir;
}


void CorrectionFactorsCw_Reweighting(){

bool doCharge = true;
bool doCentrality = true;
bool doEta = true;
bool doPt = false;
bool doPlotPt = false;
bool doPlotEta = false;
bool doScale = false;
bool doFit = false;
bool doMirrorEta = true; //turn on when using |\eta|
if(!doMirrorEta)std::cout << "Will NOT mirror eta. Please make sure to use ATLAS eta and not |eta|. " << std::endl;
///Turn on for MPT systematics
bool doMptSigmaDown = false;
bool doMptSigmaUp = false;
bool doLoosenIsoCut = false;
if(doLoosenIsoCut) std::cout << "Loosening isolation cut from 0.1 to 0.2 for systematic study." << std::endl;
if(doMptSigmaUp) std::cout << "Doing MPT +1 sigma systematics." << std::endl;
if(doMptSigmaDown) std::cout << "Doing MPT -1 sigma systematics." << std::endl;

/// --- output file ---
TDirectory *dir = gDirectory;
TString fileNameDataOut;

if(doScale) fileNameDataOut = "correctionFactorsWEtaScaled";
else if(doEta&&!doCentrality) fileNameDataOut = "correctionFactorsWEta";
else if(doCentrality&&!doEta) fileNameDataOut = "correctionFactorsWCent";
else if(doCentrality&&doEta) fileNameDataOut = "correctionFactorsWEtaCent";
else fileNameDataOut = "correctionFactorsW";
TFile *outFile = new TFile(fileNameDataOut+".root","RECREATE");

//get scale factors as fcn of eta from root file in directory
//(make graph my turning off cent binning)

TString fileNameEtaGrIn = "correctionFactorsWEta";
TFile* fGr = 0;
if(doScale) fGr = new TFile(fileNameEtaGrIn + ".root", "READ");
TGraphErrors* grToScale = 0;
TGraphErrors* grToScalePlus =0;
TGraphErrors* grToScaleMinus =0;
if(doScale){	

    std::cout << "Getting charge-inclusive eta distribution for scaling..." << std::endl;
    grToScale= (TGraphErrors*) fGr->Get("grCorrectionFactorsWEta");
    std::cout << "Done." << std::endl;

    if(doCharge){
        std::cout << "Getting charge separated eta distribution for scaling..." << std::endl;
        grToScalePlus = (TGraphErrors*) fGr->Get("grCorrectionFactorsWEtaPlus");
        grToScaleMinus = (TGraphErrors*) fGr->Get("grCorrectionFactorsWEtaMinus");
        std::cout << "Done." << std::endl;
    } 
}

gDirectory = dir;

//Correction factors for fitting and mT methods
const int chargeBins = 2;
const int globalBin  = 1;

//gROOT->LoadMacro("AtlasUtils.C");
//read in trigger efficiencies
// Nominal read in
TString fileNameIn = "HISingleMuonWmunuPYTHIADataOverlay.04.12.2013";

//PowPy8 samples
//pp
TString fileNameIn_ppPlus = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_pp.10.03.2013";
TString fileNameIn_ppMinus = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_pp.10.03.2013";
//np
TString fileNameIn_npPlus = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_np.10.03.2013";
TString fileNameIn_npMinus = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_np.10.03.2013";
//pn
TString fileNameIn_pnPlus = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_pn.10.03.2013";
TString fileNameIn_pnMinus = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_pn.10.03.2013";
//nn
TString fileNameIn_nnPlus = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wplusmunu_nn.10.03.2013";
TString fileNameIn_nnMinus = "MonteCarloFiles/HISingleMuonWmunuPYTHIADataOverlay_PowhegPythia8_AU2CT10_Wminmunu_nn.10.03.2013";

// For Diagnostics
//TString fileNameIn = "HISingleMuonWmunuPYTHIADataOverlay_pp.08.29.2013";
//TString fileNameIn = "HISingleMuonWmunuPYTHIADataOverlay.08.29.2013";
//TString fileNameIn = "HISingleMuonWmunuPYTHIADataOverlay_nn.06.21.2013";

TString fileWithWeights = "crossChecks/signalChargeEtaDistributions_CT10DataRatioFit.12.12.2013";
TFile* _fWeights = new TFile(fileWithWeights+".root","read");
//histo with weights
TH1F* hWtPlus = (TH1F*)_fWeights->Get("hRatioCT10Plus");
TH1F* hWtMinus = (TH1F*)_fWeights->Get("hRatioCT10Minus");

TString fileNameMptMtCorrelation = "mptMtCorrelations.07.21.2013";
TFile* _fMptMtCorrelation = NULL;
TProfile* _pfxMptMtCorrelation = NULL;
if(doMptSigmaUp||doMptSigmaDown){
    fileNameIn = "HISingleMuonWmunuPYTHIADataOverlay.07.14.2013";
    _fMptMtCorrelation = new TFile(fileNameMptMtCorrelation+".root","read");
    if(_fMptMtCorrelation!=0) std::cout << "Mpt-Mt correlation file: " << fileNameMptMtCorrelation << " opened." << std::endl;
    else exit(0);
}
std::cout << "Input file: " << fileNameIn << std::endl;


///IMPORTANT:IF CHANGING ETA BINNING, YOU MUST CHANGE
///THE .txt FILE HOLDING THE TRIGGER EFFICIENCIES ACCORDINGLY
std::vector <double> etaBins,etaGenBins;
//etaBins.push_back(0.0);
if(!doMirrorEta){
    etaBins.push_back(-2.4);
    etaGenBins.push_back(-2.5);
}
else{
    etaBins.push_back(0.1);
    etaGenBins.push_back(0.1);
}
if(doEta){
    
/*	etaBins.push_back(-2.35);
	etaBins.push_back(-2.3);
	etaBins.push_back(-2.25);
	etaBins.push_back(-2.2);
	etaBins.push_back(-2.15);
	etaBins.push_back(-2.1);
	etaBins.push_back(-2.05);
	etaBins.push_back(-2.0);
	etaBins.push_back(-1.975);
	etaBins.push_back(-1.950);
	etaBins.push_back(-1.925);
	etaBins.push_back(-1.9);
	etaBins.push_back(-1.8);
	etaBins.push_back(-1.7);
	etaBins.push_back(-1.6);
	etaBins.push_back(-1.5);
	etaBins.push_back(-1.4);
	etaBins.push_back(-1.35);
	etaBins.push_back(-1.3);
	etaBins.push_back(-1.25);
	etaBins.push_back(-1.2);
	etaBins.push_back(-1.15);
	etaBins.push_back(-1.1);
	etaBins.push_back(-1.05);
	etaBins.push_back(-1.0);
	etaBins.push_back(-0.9);
	etaBins.push_back(-0.8);
	etaBins.push_back(-0.7);
	etaBins.push_back(-0.6);
	etaBins.push_back(-0.5);
	etaBins.push_back(-0.4);
	etaBins.push_back(-0.3);
	etaBins.push_back(-0.2);
	etaBins.push_back(-0.1);
	etaBins.push_back(-0.075);
	etaBins.push_back(-0.05);
    etaBins.push_back(-0.025);

	etaBins.push_back(0.0);
*/
if(!doMirrorEta){
    etaBins.push_back(-2.1);
    etaBins.push_back(-1.85);
    etaBins.push_back(-1.55);
    etaBins.push_back(-1.3);
    etaBins.push_back(-1.05);
    etaBins.push_back(-0.8);
    etaBins.push_back(-0.6);
    etaBins.push_back(-0.35);
    etaBins.push_back(-0.1);
    etaBins.push_back(0.1);

    etaGenBins.push_back(-2.1);
    etaGenBins.push_back(-1.85);
    etaGenBins.push_back(-1.55);
    etaGenBins.push_back(-1.3);
    etaGenBins.push_back(-1.05);
    etaGenBins.push_back(-0.8);
    etaGenBins.push_back(-0.6);
    etaGenBins.push_back(-0.35);
    etaGenBins.push_back(-0.1);
    etaGenBins.push_back(0.1);
}
    etaBins.push_back(0.35);
    etaBins.push_back(0.6);
    etaBins.push_back(0.8);
    etaBins.push_back(1.05);
    etaBins.push_back(1.37);
    etaBins.push_back(1.52);
    etaBins.push_back(1.74);
    etaBins.push_back(2.1);

    etaGenBins.push_back(0.35);
    etaGenBins.push_back(0.6);
    etaGenBins.push_back(0.8);
    etaGenBins.push_back(1.05);
    etaGenBins.push_back(1.37);
    etaGenBins.push_back(1.52);
    etaGenBins.push_back(1.74);
    etaGenBins.push_back(2.1);

	/*etaBins.push_back(0.025);
	etaBins.push_back(0.05);
	etaBins.push_back(0.075);
	etaBins.push_back(0.1);
	etaBins.push_back(0.2);
	etaBins.push_back(0.3);
	etaBins.push_back(0.4);
	etaBins.push_back(0.5);
	etaBins.push_back(0.6);
	etaBins.push_back(0.7);
	etaBins.push_back(0.8);
	etaBins.push_back(0.9);
	etaBins.push_back(1.0);
	etaBins.push_back(1.05);
	etaBins.push_back(1.1);
	etaBins.push_back(1.15);
	etaBins.push_back(1.2);
	etaBins.push_back(1.25);
	etaBins.push_back(1.3);
	etaBins.push_back(1.35);
	etaBins.push_back(1.4);
	etaBins.push_back(1.5);
	etaBins.push_back(1.6);
	etaBins.push_back(1.7);
	etaBins.push_back(1.8);
	etaBins.push_back(1.9);
	etaBins.push_back(1.925);
	etaBins.push_back(1.950);
	etaBins.push_back(1.975);
	etaBins.push_back(2.0);
	etaBins.push_back(2.05);
	etaBins.push_back(2.1);
	etaBins.push_back(2.15);
	etaBins.push_back(2.2);
	etaBins.push_back(2.25);
	etaBins.push_back(2.3);
	etaBins.push_back(2.35);
   */ 
	}
etaBins.push_back(2.4);
//extrapolate to 2.5
etaGenBins.push_back(2.5);

const int nEtaBins = etaBins.size()-1;
std::cout << "nEtaBins : " << nEtaBins << std::endl;

//centrality
std::vector <double> centBins;
centBins.push_back(0.0);
if(doCentrality){
	centBins.push_back(0.05);
	centBins.push_back(0.1);
	centBins.push_back(0.15);
	centBins.push_back(0.2);
	centBins.push_back(0.4);
    
}
centBins.push_back(0.8);

const int nCentralityBins = centBins.size()-1;
std::cout << "nCentralityBins : " << nCentralityBins << std::endl;

///MPT lower thresholds as a function of centrality
std::vector <double> mptLow;
double mptLowCut = 25.0;
if(doMptSigmaDown){
    ///-1sigma
    mptLow.push_back(mptLowCut); //0-5%
    mptLow.push_back(mptLowCut); //5-10%
    mptLow.push_back(mptLowCut); //10-15%
    mptLow.push_back(mptLowCut); //15-20%
    mptLow.push_back(mptLowCut); //20-40%
    mptLow.push_back(mptLowCut); //40-80%
//    for(int impt=0; impt<mptLow.size(); ++impt) std::cout << mptLow[impt] << std::endl;
}
else if (doMptSigmaUp){
    ///+1 sigma
    mptLow.push_back(mptLowCut); //0-5%
    mptLow.push_back(mptLowCut); //5-10%
    mptLow.push_back(mptLowCut); //10-15%
    mptLow.push_back(mptLowCut); //15-20%
    mptLow.push_back(mptLowCut); //20-40%
    mptLow.push_back(mptLowCut); //40-80%
//    for(int impt=0; impt<mptLow.size(); ++impt) std::cout << mptLow[impt] << std::endl;
}
//Npart
std::vector <double> npartBins;
if(doCentrality){
	npartBins.push_back(382.16);//0-5
	npartBins.push_back(330.26);//5-10
	npartBins.push_back(281.88);//10-15
	npartBins.push_back(239.52);//15-20
	npartBins.push_back(157.83);//20-40
	npartBins.push_back(45.93);//40-80
   
}
else npartBins.push_back(111.63); //0-80 (weighted by bin width)

const int npartNBins = npartBins.size();

double ptmax = 300.0;
// --- Set pt and eta bins ---
std::vector<double> ptBins;
ptBins.push_back(0.0);
if(doPt){
    ptBins.push_back(30.0);
    ptBins.push_back(35.0);
    ptBins.push_back(40.0);
}
ptBins.push_back(ptmax);
const int nPtBins = ptBins.size()-1;

std::cout << "nPtBins : " << nPtBins << std::endl;

if((nEtaBins==9||(!doMirrorEta))&&nCentralityBins==6) readInputFile(nCentralityBins,nEtaBins,doMirrorEta,"triggerEffZNote_v06.txt");
else if (nEtaBins==38&&nCentralityBins==6) readInputFile(nCentralityBins,nEtaBins,doMirrorEta,"triggerEffZNote_v04.txt");
else std::cout << "WARNING: Trigger bins not initiliazed." << std::endl;
///if not mirroring eta
//readInputFile(nCentralityBins,nEtaBins,"triggerEffWNote_NoAbsEta_v01.txt");

///centrality dep plots
TGraphAsymmErrors* grppAwCentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grppCwAwCentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grppCwEffCentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grppAwCentMinus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grppCwAwCentMinus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grppCwEffCentMinus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grppAwCent = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grppCwAwCent = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grppCwEffCent = new TGraphAsymmErrors(nCentralityBins);

TGraphAsymmErrors* grnpAwCentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grnpCwAwCentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grnpCwEffCentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grnpAwCentMinus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grnpCwAwCentMinus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grnpCwEffCentMinus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grnpAwCent = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grnpCwAwCent = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grnpCwEffCent = new TGraphAsymmErrors(nCentralityBins);

TGraphAsymmErrors* grpnAwCentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grpnCwAwCentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grpnCwEffCentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grpnAwCentMinus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grpnCwAwCentMinus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grpnCwEffCentMinus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grpnAwCent = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grpnCwAwCent = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grpnCwEffCent = new TGraphAsymmErrors(nCentralityBins);

TGraphAsymmErrors* grnnAwCentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grnnCwAwCentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grnnCwEffCentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grnnAwCentMinus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grnnCwAwCentMinus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grnnCwEffCentMinus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grnnAwCent = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grnnCwAwCent = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grnnCwEffCent = new TGraphAsymmErrors(nCentralityBins);

TGraphAsymmErrors* grWtdCwEffCentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grWtdCwEffCentMinus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grWtdCwEffCent = new TGraphAsymmErrors(nCentralityBins);

TGraphAsymmErrors* grWtdCwAwEffCentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grWtdCwAwEffCentMinus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* grWtdCwAwEffCent = new TGraphAsymmErrors(nCentralityBins);


///eta dep plots
TGraphAsymmErrors* grppAwEtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grppCwAwEtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grppCwEffEtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grppAwEtaMinus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grppCwAwEtaMinus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grppCwEffEtaMinus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grppAwEta = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grppCwAwEta = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grppCwEffEta = new TGraphAsymmErrors(nEtaBins);

TGraphAsymmErrors* grnpAwEtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grnpCwAwEtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grnpCwEffEtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grnpAwEtaMinus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grnpCwAwEtaMinus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grnpCwEffEtaMinus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grnpAwEta = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grnpCwAwEta = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grnpCwEffEta = new TGraphAsymmErrors(nEtaBins);

TGraphAsymmErrors* grpnAwEtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grpnCwAwEtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grpnCwEffEtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grpnAwEtaMinus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grpnCwAwEtaMinus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grpnCwEffEtaMinus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grpnAwEta = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grpnCwAwEta = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grpnCwEffEta = new TGraphAsymmErrors(nEtaBins);

TGraphAsymmErrors* grnnAwEtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grnnCwAwEtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grnnCwEffEtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grnnAwEtaMinus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grnnCwAwEtaMinus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grnnCwEffEtaMinus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grnnAwEta = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grnnCwAwEta = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grnnCwEffEta = new TGraphAsymmErrors(nEtaBins);

TGraphAsymmErrors* grWtdCwEffEtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grWtdCwEffEtaMinus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grWtdCwEffEta = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grWtdCwAwEffEtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grWtdCwAwEffEtaMinus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* grWtdCwAwEffEta = new TGraphAsymmErrors(nEtaBins);

//residuals for systematics
//TGraph* grResid = new TGraphErrors(nCentralityBins);
TList _cWResiduals;

// --- declare variables at generator level --- //
  RooRealVar  muonGenPt("muonGenPt","p_{T}",0.0,350.0,"GeV");
  RooRealVar  nuGenPt("nuGenPt","p_{T}^{#nu}",0.0,350.0,"GeV");
  RooRealVar  munuGenMt("munuGenMt","m_{T}",0.0,350.0,"GeV");
  RooRealVar  mother("mother","mother",-30.0,30.0);
  RooRealVar  daughter("daughter","daughter",-20.0,20.0);
  RooRealVar  daughterNu("daughterNu","daughterNu",-20.0,20.0);
  RooRealVar  chargeGen("chargeGen","chargeGen",-2.0,2.0);
  RooRealVar  muEtaGen("muEtaGen","muEtaGen",-10.0,10.0);
  RooRealVar  nuEtaGen("nuEtaGen","nuEtaGen",-10.0,10.0);
  RooRealVar weight = RooRealVar("weight", "weight", 0.0, 5.0);
// --- declare cut variables at reco level --- //
  RooRealVar  muonPt("muonPt","p_{T}",0.0,350.0,"GeV");
  RooRealVar  motherRec("motherRec","motherRec",0.0,250.0);
  RooRealVar  missPt("missPt","p_{T}^{miss}",0.0,350.0,"GeV");
  RooRealVar  muonMt("muonMt","m_{T}",0.0,350.0,"GeV");
  RooRealVar  muonCharge("muonCharge","charge",-2.0,+2.0);
  RooRealVar  isolation("isolation","isolation",0.0,100.0);
  RooRealVar  centrality("centrality","centrality",0.,1.0);
  RooRealVar  muonEta("muonEta","muonEta",-3.0,+3.0);
  RooRealVar  ZDY("ZDY","ZDY",0.0,2.0);
  RooRealVar  muonQuality("muonQuality","muonQuality",0.0,15.0);
  RooRealVar  muonGenRecMatched("muonGenRecMatched","muonGenRecMatched",0.0,2.0);
  RooRealVar  muonELoss("muonELoss","muonELoss",-10.0,+10.0);
  RooRealVar  muonScat("muonScat","muonScat",-10.0,+10.0);

  TString sCutsGen = "abs(mother)==24&&abs(daughter)==13&&centrality>0.0&&centrality<0.8";
  RooArgList muonGenArgList(mother,daughter,centrality);
  RooFormulaVar cutsGen("cutsGen", "cutsGen", sCutsGen, muonGenArgList);
  /// kinematic and geometric fiducial cuts
  ///Note: try to cut on the neutrino's eta in
  ///the pseudorapidity region in which the mpt vector was calculated (i.e. the ID volume for this analysis)
  TString sCutsFid = "";
//  sCutsFid = "muonGenPt>25.0&&abs(muEtaGen)>0.1&&abs(muEtaGen)<2.4&&nuGenPt>25.0&&abs(nuEtaGen)<2.5&&munuGenMt>40.0&&abs(mother)==24&&abs(daughter)==13";
  ///NOMINAL anlysis cuts
  sCutsFid = "muonGenPt>25.0&&abs(muEtaGen)>0.1&&abs(muEtaGen)<2.5&&nuGenPt>25.0&&munuGenMt>40.0&&abs(mother)==24&&abs(daughter)==13&&centrality>0.0&&centrality<0.8";
  //cross-check w/o mpt
  //sCutsFid = "muonGenPt>25.0&&abs(muEtaGen)>0.1&&abs(muEtaGen)<2.5&&nuGenPt>0.0&&munuGenMt>40.0&&abs(mother)==24&&abs(daughter)==13&&centrality>0.0&&centrality<0.8";
  /// Diagnostics
  //sCutsFid = "muonGenPt>25.0&&abs(muEtaGen)>0.1&&abs(muEtaGen)<2.4&&nuGenPt>0.0&&munuGenMt>0.0&&abs(mother)==24&&abs(daughter)==13&&centrality>0.0&&centrality<0.8";
  //sCutsFid = "muonGenPt>25.0&&abs(muEtaGen)>0.1&&abs(muEtaGen)<2.4&&munuGenMt>40.0&&abs(mother)==24&&abs(daughter)==13&&centrality>0.0&&centrality<0.8";

    ///pt check
//  TString sCutsFid = "muonGenPt>40.0&&abs(muEtaGen)<2.4&&nuGenPt>25.0&&munuGenMt>40.0&&abs(mother)==24&&abs(daughter)==13";
//  TString sCutsGeomFid = "abs(mother)==24&&abs(daughter)==13&&abs(muEtaGen)<2.4";
  RooArgList muonGenArgList2(muonGenPt,muEtaGen,nuGenPt,munuGenMt,mother,daughter,centrality);

  RooFormulaVar cutsFid("cutsFid", "cutsFid", sCutsFid, muonGenArgList2);
//  RooFormulaVar cutsGeomFid("cutsGeomFid","cutsGeomFid",sCutsGeomFid,RooArgList(mother,daughter,muEtaGen));
  //reconstruction level cuts

  ///this is the actual correction factor (Cw)
  TString sCutsRecHiQualityWsel = "";
  if(doMptSigmaDown||doMptSigmaUp)
       sCutsRecHiQualityWsel =
            "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&isolation<0.1&&muonMt>40.0&&ZDY==0&&muonGenRecMatched==1";
   else 
       sCutsRecHiQualityWsel =
            // nominal cuts
            "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&isolation<0.1&&muonMt>40.0&&ZDY==0&&muonGenRecMatched==1";
            // cross-check w/o mpt cut
            //"muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>0.0&&missPt<9000.0&&isolation<0.1&&muonMt>40.0&&ZDY==0&&muonGenRecMatched==1";
//           Diagnostics
//            "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&isolation<0.1&&muonMt>40.0&&ZDY==0&&muonGenRecMatched==1";
            //hack for sanity check
//            "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&isolation<0.1&&muonMt>40.0&&ZDY==0"; 
      ///missPt systematic (+-10.85GeV)
      //"muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>14.15&&missPt<9000.0&&isolation<0.1&&muonMt>40.0&&ZDY==0&&muonGenRecMatched==1";
      ///use for isolation systematics
      if(doLoosenIsoCut) sCutsRecHiQualityWsel = "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&isolation<0.2&&muonMt>40.0&&ZDY==0&&muonGenRecMatched==1";

    std::cout << "sCutsRecHiQualityWsel " << sCutsRecHiQualityWsel << std::endl; 
//  TString sCutsRecLoose = "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonGenRecMatched==1";
  ///uncomment for requiring only CB muon
//  TString sCutsRecLoose = "muonGenRecMatched==1&&muonQuality>0";
  ///uncommment for no requirements on reco muons
//  TString sCutsRecLoose = "muonGenRecMatched==1";
  //reco level cuts for systematics
//  TString sCutsRec ="muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&isolation<0.2&&muonMt>40.0&&ZDY==0&&muonGenRecMatched==1";
  
//  TString sCutsRec = "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>33.5&&missPt>25.0&&isolation<0.1&&muonMt>40.0&&ZDY==0";
  RooArgList muonRecArgList5(muonEta,muonPt,missPt,muonGenRecMatched,muonMt);
  muonRecArgList5.add(muonQuality); 
  muonRecArgList5.add(muonELoss);
  muonRecArgList5.add(muonScat);
  muonRecArgList5.add(ZDY);
  muonRecArgList5.add(isolation);

  ///cuts on datasets at RECONSTRUCTION level
  RooFormulaVar cutsRecHiQualityWsel("cutsRecHiQualityWsel", "cutsRecHiQualityWsel",
        sCutsRecHiQualityWsel, muonRecArgList5);

  RooCategory muonCategory("muonCategory","muonCategory");
  muonCategory.defineType("Z",23);
  muonCategory.defineType("W",24);

  RooCategory chargeCategory("chargeCategory","chargeCategory") ;
  chargeCategory.defineType("muMinus",-1) ;
  chargeCategory.defineType("muPlus",1) ;

  RooCategory chargeGenCategory("chargeGenCategory","chargeGenCategory") ;
  chargeGenCategory.defineType("muMinus",-1) ;
  chargeGenCategory.defineType("muPlus",1) ;

  RooArgSet muonGenArgSet(muonGenPt,nuGenPt,munuGenMt,mother,daughter,chargeGenCategory,muEtaGen,nuEtaGen,centrality);
  muonGenArgSet.add(chargeGen);
  muonGenArgSet.add(daughterNu);
  muonGenArgSet.add(weight);

  RooArgSet muonRecArgSet(muonEta,centrality,ZDY,muonCategory,chargeCategory);
  muonRecArgSet.add(muonPt);
  muonRecArgSet.add(motherRec);
  muonRecArgSet.add(missPt);
  muonRecArgSet.add(muonMt);
  muonRecArgSet.add(muonCharge);
  muonRecArgSet.add(isolation);
  muonRecArgSet.add(muonQuality);
  muonRecArgSet.add(muonGenRecMatched);
  muonRecArgSet.add(muonELoss);
  muonRecArgSet.add(muonScat);
  muonRecArgSet.add(muEtaGen);
  muonRecArgSet.add(weight);


  // Create separate d.s. for each nucleon combination
  RooDataSet* mcppWGenSetPlus = fillHIMuonGenSet("/usatlas/u/tbales/scratch/",fileNameIn_ppPlus+".root",muonGenArgSet,hWtPlus); mcppWGenSetPlus->Print();
  double nWplus_pp = mcppWGenSetPlus->sumEntries();
  RooDataSet* mcppWGenSetMinus = fillHIMuonGenSet("/usatlas/u/tbales/scratch/",fileNameIn_ppMinus+".root",muonGenArgSet,hWtMinus); mcppWGenSetMinus->Print();
  double nWminu_pp = mcppWGenSetMinus->sumEntries();
  RooDataSet* mcnpWGenSetPlus = fillHIMuonGenSet("/usatlas/u/tbales/scratch/",fileNameIn_npPlus+".root",muonGenArgSet,hWtPlus); mcnpWGenSetPlus->Print();
  double nWplus_np = mcnpWGenSetPlus->sumEntries();
  RooDataSet* mcnpWGenSetMinus = fillHIMuonGenSet("/usatlas/u/tbales/scratch/",fileNameIn_npMinus+".root",muonGenArgSet,hWtMinus); mcnpWGenSetMinus->Print();
  double nWminu_np = mcnpWGenSetMinus->sumEntries();
  RooDataSet* mcpnWGenSetPlus = fillHIMuonGenSet("/usatlas/u/tbales/scratch/",fileNameIn_pnPlus+".root",muonGenArgSet,hWtPlus); mcpnWGenSetPlus->Print();
  double nWplus_pn = mcpnWGenSetPlus->sumEntries();
  RooDataSet* mcpnWGenSetMinus = fillHIMuonGenSet("/usatlas/u/tbales/scratch/",fileNameIn_pnMinus+".root",muonGenArgSet,hWtMinus); mcpnWGenSetMinus->Print();
  double nWminu_pn = mcpnWGenSetMinus->sumEntries();
  RooDataSet* mcnnWGenSetPlus = fillHIMuonGenSet("/usatlas/u/tbales/scratch/",fileNameIn_nnPlus+".root",muonGenArgSet,hWtPlus); mcnnWGenSetPlus->Print();
  double nWplus_nn = mcnnWGenSetPlus->sumEntries();
  RooDataSet* mcnnWGenSetMinus = fillHIMuonGenSet("/usatlas/u/tbales/scratch/",fileNameIn_nnMinus+".root",muonGenArgSet,hWtMinus); mcnnWGenSetMinus->Print();
  double nWminu_nn = mcnnWGenSetMinus->sumEntries();

  //Total number of events
  double nWplusTot = nWplus_pp+nWplus_np+nWplus_pn+nWplus_nn;
  double wtPlus_pp = 0.155*nWplusTot/nWplus_pp; 
  double wtPlus_np = (0.478/2.0)*nWplusTot/nWplus_np; 
  double wtPlus_pn = (0.478/2.0)*nWplusTot/nWplus_pn; 
  double wtPlus_nn = 0.367*nWplusTot/nWplus_nn; 

  double nWminusTot = nWplus_pp+nWplus_np+nWplus_pn+nWplus_nn;
  double wtMinus_pp = 0.155*nWminusTot/nWplus_pp; 
  double wtMinus_np = (0.478/2.0)*nWminusTot/nWplus_np; 
  double wtMinus_pn = (0.478/2.0)*nWminusTot/nWplus_pn; 
  double wtMinus_nn = 0.367*nWminusTot/nWplus_nn; 

  mcppWGenSetPlus = (RooDataSet*)mcppWGenSetPlus->reduce(Cut(cutsGen)); 
  mcppWGenSetMinus = (RooDataSet*)mcppWGenSetMinus->reduce(Cut(cutsGen)); 
  mcnpWGenSetPlus = (RooDataSet*)mcnpWGenSetPlus->reduce(Cut(cutsGen)); 
  mcnpWGenSetMinus = (RooDataSet*)mcnpWGenSetMinus->reduce(Cut(cutsGen)); 
  mcpnWGenSetPlus = (RooDataSet*)mcpnWGenSetPlus->reduce(Cut(cutsGen)); 
  mcpnWGenSetMinus = (RooDataSet*)mcpnWGenSetMinus->reduce(Cut(cutsGen)); 
  mcnnWGenSetPlus = (RooDataSet*)mcnnWGenSetPlus->reduce(Cut(cutsGen)); 
  mcnnWGenSetMinus = (RooDataSet*)mcnnWGenSetMinus->reduce(Cut(cutsGen)); 

  // Generated muons in fiducial region
  RooDataSet* mcppWFidSetPlus = fillHIMuonGenSet("/usatlas/u/tbales/scratch/",fileNameIn_ppPlus+".root",muonGenArgSet,hWtPlus); 
  mcppWFidSetPlus = (RooDataSet*)mcppWFidSetPlus->reduce(Cut(cutsFid));
  RooDataSet* mcppWFidSetMinus = fillHIMuonGenSet("/usatlas/u/tbales/scratch/",fileNameIn_ppMinus+".root",muonGenArgSet,hWtMinus); 
  mcppWFidSetMinus = (RooDataSet*)mcppWFidSetMinus->reduce(Cut(cutsFid)); 

  RooDataSet* mcnpWFidSetPlus = fillHIMuonGenSet("/usatlas/u/tbales/scratch/",fileNameIn_npPlus+".root",muonGenArgSet,hWtPlus); 
  mcnpWFidSetPlus = (RooDataSet*)mcnpWFidSetPlus->reduce(Cut(cutsFid)); 
  RooDataSet* mcnpWFidSetMinus = fillHIMuonGenSet("/usatlas/u/tbales/scratch/",fileNameIn_npMinus+".root",muonGenArgSet,hWtMinus); 
  mcnpWFidSetMinus = (RooDataSet*)mcnpWFidSetMinus->reduce(Cut(cutsFid)); 

  RooDataSet* mcpnWFidSetPlus = fillHIMuonGenSet("/usatlas/u/tbales/scratch/",fileNameIn_pnPlus+".root",muonGenArgSet,hWtPlus); 
  mcpnWFidSetPlus = (RooDataSet*)mcpnWFidSetPlus->reduce(Cut(cutsFid)); 
  RooDataSet* mcpnWFidSetMinus = fillHIMuonGenSet("/usatlas/u/tbales/scratch/",fileNameIn_pnMinus+".root",muonGenArgSet,hWtMinus); 
  mcpnWFidSetMinus = (RooDataSet*)mcpnWFidSetMinus->reduce(Cut(cutsFid)); 

  RooDataSet* mcnnWFidSetPlus = fillHIMuonGenSet("/usatlas/u/tbales/scratch/",fileNameIn_nnPlus+".root",muonGenArgSet,hWtPlus); 
  mcnnWFidSetPlus = (RooDataSet*)mcnnWFidSetPlus->reduce(Cut(cutsFid)); 
  RooDataSet* mcnnWFidSetMinus = fillHIMuonGenSet("/usatlas/u/tbales/scratch/",fileNameIn_nnMinus+".root",muonGenArgSet,hWtMinus); 
  mcnnWFidSetMinus = (RooDataSet*)mcnnWFidSetMinus->reduce(Cut(cutsFid)); 

  // Aw for each dataset
  // pp
  // Aw+
  double pp_awPlus,pp_errAwPlusLo,pp_errAwPlusHi;
  double pp_fidPlus = mcppWFidSetPlus->sumEntries(); double pp_genPlus = mcppWGenSetPlus->sumEntries();
  Efficiency(pp_fidPlus,pp_genPlus,0.683,pp_awPlus,pp_errAwPlusLo,pp_errAwPlusHi);
  std::cout << " Aw+ = " << mcppWFidSetPlus->sumEntries() << "/" << mcppWGenSetPlus->sumEntries() << " = " << pp_fidPlus/pp_genPlus << std::endl;
  std::cout << "Consistency check: " << pp_fidPlus/pp_genPlus << " =? " << pp_awPlus << std::endl;
  // Aw-
  double pp_awMinus,pp_errAwMinusLo,pp_errAwMinusHi;
  double pp_fidMinus = mcppWFidSetMinus->sumEntries(); double pp_genMinus = mcppWGenSetMinus->sumEntries();
  Efficiency(pp_fidMinus,pp_genMinus,0.683,pp_awMinus,pp_errAwMinusLo,pp_errAwMinusHi);
  std::cout << " Aw- = " << mcppWFidSetMinus->sumEntries() << "/" << mcppWGenSetMinus->sumEntries() << " = " << pp_fidMinus/pp_genMinus << std::endl;
  std::cout << "Consistency check: " << pp_fidMinus/pp_genMinus << " =? " << pp_awMinus << std::endl;
  // np
  // Aw+
  double np_awPlus,np_errAwPlusLo,np_errAwPlusHi;
  double np_fidPlus = mcnpWFidSetPlus->sumEntries(); double np_genPlus = mcnpWGenSetPlus->sumEntries();
  Efficiency(np_fidPlus,np_genPlus,0.683,np_awPlus,np_errAwPlusLo,np_errAwPlusHi);
  std::cout << " Aw+ = " << mcnpWFidSetPlus->sumEntries() << "/" << mcnpWGenSetPlus->sumEntries() << " = " << np_fidPlus/np_genPlus << std::endl;
  std::cout << "Consistency check: " << np_fidPlus/np_genPlus << " =? " << np_awPlus << std::endl;
  // Aw-
  double np_awMinus,np_errAwMinusLo,np_errAwMinusHi;
  double np_fidMinus = mcnpWFidSetMinus->sumEntries(); double np_genMinus = mcnpWGenSetMinus->sumEntries();
  Efficiency(np_fidMinus,np_genMinus,0.683,np_awMinus,np_errAwMinusLo,np_errAwMinusHi);
  std::cout << " Aw- = " << mcnpWFidSetMinus->sumEntries() << "/" << mcnpWGenSetMinus->sumEntries() << " = " << np_fidMinus/np_genMinus << std::endl;
  std::cout << "Consistency check: " << np_fidMinus/np_genMinus << " =? " << np_awMinus << std::endl;
  // pn
  // Aw+
  double pn_awPlus,pn_errAwPlusLo,pn_errAwPlusHi;
  double pn_fidPlus = mcpnWFidSetPlus->sumEntries(); double pn_genPlus = mcpnWGenSetPlus->sumEntries();
  Efficiency(pn_fidPlus,pn_genPlus,0.683,pn_awPlus,pn_errAwPlusLo,pn_errAwPlusHi);
  std::cout << " Aw+ = " << mcpnWFidSetPlus->sumEntries() << "/" << mcpnWGenSetPlus->sumEntries() << " = " << pn_fidPlus/pn_genPlus << std::endl;
  std::cout << "Consistency check: " << pn_fidPlus/pn_genPlus << " =? " << pn_awPlus << std::endl;
  // Aw-
  double pn_awMinus,pn_errAwMinusLo,pn_errAwMinusHi;
  double pn_fidMinus = mcpnWFidSetMinus->sumEntries(); double pn_genMinus = mcpnWGenSetMinus->sumEntries();
  Efficiency(pn_fidMinus,pn_genMinus,0.683,pn_awMinus,pn_errAwMinusLo,pn_errAwMinusHi);
  std::cout << " Aw- = " << mcpnWFidSetMinus->sumEntries() << "/" << mcpnWGenSetMinus->sumEntries() << " = " << pn_fidMinus/pn_genMinus << std::endl;
  std::cout << "Consistency check: " << pn_fidMinus/pn_genMinus << " =? " << pn_awMinus << std::endl;
  // nn
  // Aw+
  double nn_awPlus,nn_errAwPlusLo,nn_errAwPlusHi;
  double nn_fidPlus = mcnnWFidSetPlus->sumEntries(); double nn_genPlus = mcnnWGenSetPlus->sumEntries();
  Efficiency(nn_fidPlus,nn_genPlus,0.683,nn_awPlus,nn_errAwPlusLo,nn_errAwPlusHi);
  std::cout << " Aw+ = " << mcnnWFidSetPlus->sumEntries() << "/" << mcnnWGenSetPlus->sumEntries() << " = " << nn_fidPlus/nn_genPlus << std::endl;
  std::cout << "Consistency check: " << nn_fidPlus/nn_genPlus << " =? " << nn_awPlus << std::endl;
  // Aw-
  double nn_awMinus,nn_errAwMinusLo,nn_errAwMinusHi;
  double nn_fidMinus = mcnnWFidSetMinus->sumEntries(); double nn_genMinus = mcnnWGenSetMinus->sumEntries();
  Efficiency(nn_fidMinus,nn_genMinus,0.683,nn_awMinus,nn_errAwMinusLo,nn_errAwMinusHi);
  std::cout << " Aw- = " << mcnnWFidSetMinus->sumEntries() << "/" << mcnnWGenSetMinus->sumEntries() << " = " << nn_fidMinus/nn_genMinus << std::endl;
  std::cout << "Consistency check: " << nn_fidMinus/nn_genMinus << " =? " << nn_awMinus << std::endl;




  

  RooBinning b = RooBinning(0.0,300.0);
  b.addUniform(90,0.0,90.0);
  b.addUniform(53,90.0,300.0);
  if(doPlotPt){
	// pp
	//entire generator set
	TH1F* hppGenPlus = (TH1F*)mcppWGenSetPlus->createHistogram("hppGenPlus",muonGenPt,Binning(b)); 
	TH1F* hppGenMinus = (TH1F*)mcppWGenSetMinus->createHistogram("hppGenMinus",muonGenPt,Binning(b)); 
	//fiducial generator set
	TH1F* hppFidPlus = (TH1F*)mcppWFidSetPlus->createHistogram("hppFidPlus",muonGenPt,Binning(b)); 
	TH1F* hppFidMinus = (TH1F*)mcppWFidSetMinus->createHistogram("hppFidMinus",muonGenPt,Binning(b)); 

  	plotGeneratorDistros(hppGenPlus,hppFidPlus,102,"W^{+}#rightarrow#mu^{+}#nu,pp");
  	plotGeneratorDistros(hppGenMinus,hppFidMinus,103,"W^{-}#rightarrow#mu^{-}#nu,pp");

	// np
	//entire generator set
	TH1F* hnpGenPlus = (TH1F*)mcnpWGenSetPlus->createHistogram("hnpGenPlus",muonGenPt,Binning(b)); 
	TH1F* hnpGenMinus = (TH1F*)mcnpWGenSetMinus->createHistogram("hnpGenMinus",muonGenPt,Binning(b)); 
	//fiducial generator set
	TH1F* hnpFidPlus = (TH1F*)mcnpWFidSetPlus->createHistogram("hnpFidPlus",muonGenPt,Binning(b)); 
	TH1F* hnpFidMinus = (TH1F*)mcnpWFidSetMinus->createHistogram("hnpFidMinus",muonGenPt,Binning(b)); 

  	plotGeneratorDistros(hnpGenPlus,hnpFidPlus,102,"W^{+}#rightarrow#mu^{+}#nu,np");
  	plotGeneratorDistros(hnpGenMinus,hnpFidMinus,103,"W^{-}#rightarrow#mu^{-}#nu,np");

	// pn
	//entire generator set
	TH1F* hpnGenPlus = (TH1F*)mcpnWGenSetPlus->createHistogram("hpnGenPlus",muonGenPt,Binning(b)); 
	TH1F* hpnGenMinus = (TH1F*)mcpnWGenSetMinus->createHistogram("hpnGenMinus",muonGenPt,Binning(b)); 
	//fiducial generator set
	TH1F* hpnFidPlus = (TH1F*)mcpnWFidSetPlus->createHistogram("hpnFidPlus",muonGenPt,Binning(b)); 
	TH1F* hpnFidMinus = (TH1F*)mcpnWFidSetMinus->createHistogram("hpnFidMinus",muonGenPt,Binning(b)); 

  	plotGeneratorDistros(hpnGenPlus,hpnFidPlus,102,"W^{+}#rightarrow#mu^{+}#nu,pn");
  	plotGeneratorDistros(hpnGenMinus,hpnFidMinus,103,"W^{-}#rightarrow#mu^{-}#nu,pn");

	// nn
	//entire generator set
	TH1F* hnnGenPlus = (TH1F*)mcnnWGenSetPlus->createHistogram("hnnGenPlus",muonGenPt,Binning(b)); 
	TH1F* hnnGenMinus = (TH1F*)mcnnWGenSetMinus->createHistogram("hnnGenMinus",muonGenPt,Binning(b)); 
	//fiducial generator set
	TH1F* hnnFidPlus = (TH1F*)mcnnWFidSetPlus->createHistogram("hnnFidPlus",muonGenPt,Binning(b)); 
	TH1F* hnnFidMinus = (TH1F*)mcnnWFidSetMinus->createHistogram("hnnFidMinus",muonGenPt,Binning(b)); 

  	plotGeneratorDistros(hnnGenPlus,hnnFidPlus,102,"W^{+}#rightarrow#mu^{+}#nu,nn");
  	plotGeneratorDistros(hnnGenMinus,hnnFidMinus,103,"W^{-}#rightarrow#mu^{-}#nu,nn");

  }  

  if(doPlotEta){

    double etaLow, etaUpp;
    etaLow = -6.0; etaUpp=6.0;
    RooBinning b2 = RooBinning(29,etaLow,etaUpp);
    
    b2.addUniform(1,etaLow,-5.0);
    b2.addUniform(1,-5.0,-4.0);
    b2.addUniform(1,-4.0,-3.5);
    b2.addUniform(1,-3.5,-3.0);
    b2.addUniform(1,-3.0,-2.4);
    b2.addUniform(1,-2.4,2.1);
    b2.addUniform(1,-2.1,-1.85);
    b2.addUniform(1,-1.85,-1.55);
    b2.addUniform(1,-1.55,-1.3);
    b2.addUniform(1,-1.3,-1.05);
    b2.addUniform(1,-1.05,-0.8);
    b2.addUniform(1,-0.8,-0.6);
    b2.addUniform(1,-0.6,-0.35);
    b2.addUniform(1,-0.35,-0.1);
    b2.addUniform(1,-0.1,0.1);
    b2.addUniform(1,0.1,0.35);
    b2.addUniform(1,0.35,0.6);
    b2.addUniform(1,0.6,0.8);
    b2.addUniform(1,0.8,1.05);
    b2.addUniform(1,1.05,1.3);
    b2.addUniform(1,1.3,1.55);
    b2.addUniform(1,1.55,1.85);
    b2.addUniform(1,1.85,2.1);
    b2.addUniform(1,2.1,2.4);
    b2.addUniform(1,2.4,3.0);
    b2.addUniform(1,3.0,3.5);
    b2.addUniform(1,3.5,4.0);
    b2.addUniform(1,4.0,5.0);
    b2.addUniform(1,5.0,etaUpp);

    // pp
    // entire generator set
    TH1F* hppGenPlus = (TH1F*)mcppWGenSetPlus->createHistogram("hppGenPlus",muEtaGen,Binning(b2)); 
    TH1F* hppGenMinus = (TH1F*)mcppWGenSetMinus->createHistogram("hppGenMinus",muEtaGen,Binning(b2)); 
    // fiducial generator set
    TH1F* hppFidPlus = (TH1F*)mcppWFidSetPlus->createHistogram("hppFidPlus",muEtaGen,Binning(b2)); 
    TH1F* hppFidMinus = (TH1F*)mcppWFidSetMinus->createHistogram("hppFidMinus",muEtaGen,Binning(b2)); 

    plotGeneratorDistros(hppGenPlus,hppFidPlus,102,"W^{+}#rightarrow#mu^{+}#nu,pp");
    plotGeneratorDistros(hppGenMinus,hppFidMinus,103,"W^{-}#rightarrow#mu^{-}#nu,pp");

    // np
    // entire generator set
    TH1F* hnpGenPlus = (TH1F*)mcnpWGenSetPlus->createHistogram("hnpGenPlus",muEtaGen,Binning(b2)); 
    TH1F* hnpGenMinus = (TH1F*)mcnpWGenSetMinus->createHistogram("hnpGenMinus",muEtaGen,Binning(b2)); 
    // fiducial generator set
    TH1F* hnpFidPlus = (TH1F*)mcnpWFidSetPlus->createHistogram("hnpFidPlus",muEtaGen,Binning(b2)); 
    TH1F* hnpFidMinus = (TH1F*)mcnpWFidSetMinus->createHistogram("hnpFidMinus",muEtaGen,Binning(b2)); 

    plotGeneratorDistros(hnpGenPlus,hnpFidPlus,102,"W^{+}#rightarrow#mu^{+}#nu,np");
    plotGeneratorDistros(hnpGenMinus,hnpFidMinus,103,"W^{-}#rightarrow#mu^{-}#nu,np");

    // pn
    // entire generator set
    TH1F* hpnGenPlus = (TH1F*)mcpnWGenSetPlus->createHistogram("hpnGenPlus",muEtaGen,Binning(b2)); 
    TH1F* hpnGenMinus = (TH1F*)mcpnWGenSetMinus->createHistogram("hpnGenMinus",muEtaGen,Binning(b2)); 
    // fiducial generator set
    TH1F* hpnFidPlus = (TH1F*)mcpnWFidSetPlus->createHistogram("hpnFidPlus",muEtaGen,Binning(b2)); 
    TH1F* hpnFidMinus = (TH1F*)mcpnWFidSetMinus->createHistogram("hpnFidMinus",muEtaGen,Binning(b2)); 

    plotGeneratorDistros(hpnGenPlus,hpnFidPlus,102,"W^{+}#rightarrow#mu^{+}#nu,pn");
    plotGeneratorDistros(hpnGenMinus,hpnFidMinus,103,"W^{-}#rightarrow#mu^{-}#nu,pn");

    // nn
    // entire generator set
    TH1F* hnnGenPlus = (TH1F*)mcnnWGenSetPlus->createHistogram("hnnGenPlus",muEtaGen,Binning(b2)); 
    TH1F* hnnGenMinus = (TH1F*)mcnnWGenSetMinus->createHistogram("hnnGenMinus",muEtaGen,Binning(b2)); 
    // fiducial generator set
    TH1F* hnnFidPlus = (TH1F*)mcnnWFidSetPlus->createHistogram("hnnFidPlus",muEtaGen,Binning(b2)); 
    TH1F* hnnFidMinus = (TH1F*)mcnnWFidSetMinus->createHistogram("hnnFidMinus",muEtaGen,Binning(b2)); 

    plotGeneratorDistros(hnnGenPlus,hnnFidPlus,102,"W^{+}#rightarrow#mu^{+}#nu,nn");
    plotGeneratorDistros(hnnGenMinus,hnnFidMinus,103,"W^{-}#rightarrow#mu^{-}#nu,nn");


  }

  // Datasets at reco level
  // pp
  // mu+
  RooDataSet* mcppSetWRecHiQualityWselPlus =
        fillHIMuonRecSet("/usatlas/u/tbales/scratch/",fileNameIn_ppPlus+".root",muonRecArgSet,hWtPlus,doMptSigmaUp,doMptSigmaDown); 
  mcppSetWRecHiQualityWselPlus =
        (RooDataSet*)mcppSetWRecHiQualityWselPlus->reduce(Cut(cutsRecHiQualityWsel)); 
  mcppSetWRecHiQualityWselPlus = (RooDataSet*)mcppSetWRecHiQualityWselPlus->reduce(Cut("motherRec==24"));
  //mu-
  RooDataSet* mcppSetWRecHiQualityWselMinus =
        fillHIMuonRecSet("/usatlas/u/tbales/scratch/",fileNameIn_ppMinus+".root",muonRecArgSet,hWtMinus,doMptSigmaUp,doMptSigmaDown); 
  mcppSetWRecHiQualityWselMinus =
        (RooDataSet*)mcppSetWRecHiQualityWselMinus->reduce(Cut(cutsRecHiQualityWsel)); 
  mcppSetWRecHiQualityWselMinus = (RooDataSet*)mcppSetWRecHiQualityWselMinus->reduce(Cut("motherRec==24"));

  // np
  // mu+
  RooDataSet* mcnpSetWRecHiQualityWselPlus =
        fillHIMuonRecSet("/usatlas/u/tbales/scratch/",fileNameIn_npPlus+".root",muonRecArgSet,hWtPlus,doMptSigmaUp,doMptSigmaDown); 
  mcnpSetWRecHiQualityWselPlus =
        (RooDataSet*)mcnpSetWRecHiQualityWselPlus->reduce(Cut(cutsRecHiQualityWsel)); 
  mcnpSetWRecHiQualityWselPlus = (RooDataSet*)mcnpSetWRecHiQualityWselPlus->reduce(Cut("motherRec==24"));
  //mu-
  RooDataSet* mcnpSetWRecHiQualityWselMinus =
        fillHIMuonRecSet("/usatlas/u/tbales/scratch/",fileNameIn_npMinus+".root",muonRecArgSet,hWtMinus,doMptSigmaUp,doMptSigmaDown); 
  mcnpSetWRecHiQualityWselMinus =
        (RooDataSet*)mcnpSetWRecHiQualityWselMinus->reduce(Cut(cutsRecHiQualityWsel)); 
  mcnpSetWRecHiQualityWselMinus = (RooDataSet*)mcnpSetWRecHiQualityWselMinus->reduce(Cut("motherRec==24"));

  // pn
  // mu+
  RooDataSet* mcpnSetWRecHiQualityWselPlus =
        fillHIMuonRecSet("/usatlas/u/tbales/scratch/",fileNameIn_pnPlus+".root",muonRecArgSet,hWtPlus,doMptSigmaUp,doMptSigmaDown); 
  mcpnSetWRecHiQualityWselPlus =
        (RooDataSet*)mcpnSetWRecHiQualityWselPlus->reduce(Cut(cutsRecHiQualityWsel)); 
  mcpnSetWRecHiQualityWselPlus = (RooDataSet*)mcpnSetWRecHiQualityWselPlus->reduce(Cut("motherRec==24"));
  //mu-
  RooDataSet* mcpnSetWRecHiQualityWselMinus =
        fillHIMuonRecSet("/usatlas/u/tbales/scratch/",fileNameIn_pnMinus+".root",muonRecArgSet,hWtMinus,doMptSigmaUp,doMptSigmaDown); 
  mcpnSetWRecHiQualityWselMinus =
        (RooDataSet*)mcpnSetWRecHiQualityWselMinus->reduce(Cut(cutsRecHiQualityWsel)); 
  mcpnSetWRecHiQualityWselMinus = (RooDataSet*)mcpnSetWRecHiQualityWselMinus->reduce(Cut("motherRec==24"));

  // nn
  // mu+
  RooDataSet* mcnnSetWRecHiQualityWselPlus =
        fillHIMuonRecSet("/usatlas/u/tbales/scratch/",fileNameIn_nnPlus+".root",muonRecArgSet,hWtPlus,doMptSigmaUp,doMptSigmaDown); 
  mcnnSetWRecHiQualityWselPlus =
        (RooDataSet*)mcnnSetWRecHiQualityWselPlus->reduce(Cut(cutsRecHiQualityWsel)); 
  mcnnSetWRecHiQualityWselPlus = (RooDataSet*)mcnnSetWRecHiQualityWselPlus->reduce(Cut("motherRec==24"));
  //mu-
  RooDataSet* mcnnSetWRecHiQualityWselMinus =
        fillHIMuonRecSet("/usatlas/u/tbales/scratch/",fileNameIn_nnMinus+".root",muonRecArgSet,hWtMinus,doMptSigmaUp,doMptSigmaDown); 
  mcnnSetWRecHiQualityWselMinus =
        (RooDataSet*)mcnnSetWRecHiQualityWselMinus->reduce(Cut(cutsRecHiQualityWsel)); 
  mcnnSetWRecHiQualityWselMinus = (RooDataSet*)mcnnSetWRecHiQualityWselMinus->reduce(Cut("motherRec==24"));

  // --- Subdivide in bins ---
  // pp
  // mu+
  // Generator level subsets (Aw denominator)
  RooDataSet* mcppWGenSubSetPlus[nPtBins][nEtaBins][nCentralityBins];
  // Generator level subsets (Aw numerator, Cw denominator)
  RooDataSet* mcppWFidSubSetPlus[nPtBins][nEtaBins][nCentralityBins];
  // Reconstructed sub-sets (Cw numerator)
  RooDataSet* mcppSubSetWRecHiQualityWselPlus[nPtBins][nEtaBins][nCentralityBins];
  // mu-
  // Generator level subsets (Aw denominator)
  RooDataSet* mcppWGenSubSetMinus[nPtBins][nEtaBins][nCentralityBins];
  // Generator level subsets (Aw numerator, Cw denominator)
  RooDataSet* mcppWFidSubSetMinus[nPtBins][nEtaBins][nCentralityBins];
  // Reconstructed sub-sets (Cw numerator)
  RooDataSet* mcppSubSetWRecHiQualityWselMinus[nPtBins][nEtaBins][nCentralityBins];

  // np
  // mu+
  // Generator level subsets (Aw denominator)
  RooDataSet* mcnpWGenSubSetPlus[nPtBins][nEtaBins][nCentralityBins];
  // Generator level subsets (Aw numerator, Cw denominator)
  RooDataSet* mcnpWFidSubSetPlus[nPtBins][nEtaBins][nCentralityBins];
  // Reconstructed sub-sets (Cw numerator)
  RooDataSet* mcnpSubSetWRecHiQualityWselPlus[nPtBins][nEtaBins][nCentralityBins];
  // mu-
  // Generator level subsets (Aw denominator)
  RooDataSet* mcnpWGenSubSetMinus[nPtBins][nEtaBins][nCentralityBins];
  // Generator level subsets (Aw numerator, Cw denominator)
  RooDataSet* mcnpWFidSubSetMinus[nPtBins][nEtaBins][nCentralityBins];
  // Reconstructed sub-sets (Cw numerator)
  RooDataSet* mcnpSubSetWRecHiQualityWselMinus[nPtBins][nEtaBins][nCentralityBins];

  // pn
  // mu+
  // Generator level subsets (Aw denominator)
  RooDataSet* mcpnWGenSubSetPlus[nPtBins][nEtaBins][nCentralityBins];
  // Generator level subsets (Aw numerator, Cw denominator)
  RooDataSet* mcpnWFidSubSetPlus[nPtBins][nEtaBins][nCentralityBins];
  // Reconstructed sub-sets (Cw numerator)
  RooDataSet* mcpnSubSetWRecHiQualityWselPlus[nPtBins][nEtaBins][nCentralityBins];
  // mu-
  // Generator level subsets (Aw denominator)
  RooDataSet* mcpnWGenSubSetMinus[nPtBins][nEtaBins][nCentralityBins];
  // Generator level subsets (Aw numerator, Cw denominator)
  RooDataSet* mcpnWFidSubSetMinus[nPtBins][nEtaBins][nCentralityBins];
  // Reconstructed sub-sets (Cw numerator)
  RooDataSet* mcpnSubSetWRecHiQualityWselMinus[nPtBins][nEtaBins][nCentralityBins];

  // nn
  // mu+
  // Generator level subsets (Aw denominator)
  RooDataSet* mcnnWGenSubSetPlus[nPtBins][nEtaBins][nCentralityBins];
  // Generator level subsets (Aw numerator, Cw denominator)
  RooDataSet* mcnnWFidSubSetPlus[nPtBins][nEtaBins][nCentralityBins];
  // Reconstructed sub-sets (Cw numerator)
  RooDataSet* mcnnSubSetWRecHiQualityWselPlus[nPtBins][nEtaBins][nCentralityBins];
  // mu-
  // Generator level subsets (Aw denominator)
  RooDataSet* mcnnWGenSubSetMinus[nPtBins][nEtaBins][nCentralityBins];
  // Generator level subsets (Aw numerator, Cw denominator)
  RooDataSet* mcnnWFidSubSetMinus[nPtBins][nEtaBins][nCentralityBins];
  // Reconstructed sub-sets (Cw numerator)
  RooDataSet* mcnnSubSetWRecHiQualityWselMinus[nPtBins][nEtaBins][nCentralityBins];

  // Cut the datasets in pt,eta,cent bins
  for ( int i = 0; i < nPtBins; i++ ) {
    for ( int j = 0; j < nEtaBins; j++ ) {
	//add residual graphs to the list
	_cWResiduals.Add( new TGraph(nCentralityBins) );

      ///set first bool false to not mirror eta
      for ( int k = 0; k < nCentralityBins; k++ ){

        std::cout << "pt :" << ptBins[i] << "-" << ptBins[i+1] << " eta rec:" << etaBins[j] << "-" << etaBins[j+1] << 
            " eta gen:" << etaGenBins[j] << "-" << etaGenBins[j+1] <<
            " cent: " << centBins[k] << "-" << centBins[k+1] << std::endl;
        if(doEta){
	    // pp
            mcppWGenSubSetPlus[i][j][k] = selectPtEtaCentrality(mcppWGenSetPlus ,ptBins[i],ptBins[i+1], etaGenBins[j],
                etaGenBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 
            mcppWGenSubSetMinus[i][j][k] = selectPtEtaCentrality(mcppWGenSetMinus ,ptBins[i],ptBins[i+1], etaGenBins[j],
                etaGenBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 

	    // np
            mcnpWGenSubSetPlus[i][j][k] = selectPtEtaCentrality(mcnpWGenSetPlus ,ptBins[i],ptBins[i+1], etaGenBins[j],
                etaGenBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 
            mcnpWGenSubSetMinus[i][j][k] = selectPtEtaCentrality(mcnpWGenSetMinus ,ptBins[i],ptBins[i+1], etaGenBins[j],
                etaGenBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 

	    // pn
            mcpnWGenSubSetPlus[i][j][k] = selectPtEtaCentrality(mcpnWGenSetPlus ,ptBins[i],ptBins[i+1], etaGenBins[j],
                etaGenBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 
            mcpnWGenSubSetMinus[i][j][k] = selectPtEtaCentrality(mcpnWGenSetMinus ,ptBins[i],ptBins[i+1], etaGenBins[j],
                etaGenBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 

	    // nn
            mcnnWGenSubSetPlus[i][j][k] = selectPtEtaCentrality(mcnnWGenSetPlus ,ptBins[i],ptBins[i+1], etaGenBins[j],
                etaGenBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 
            mcnnWGenSubSetMinus[i][j][k] = selectPtEtaCentrality(mcnnWGenSetMinus ,ptBins[i],ptBins[i+1], etaGenBins[j],
                etaGenBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 
        }
        else{
            mcppWGenSubSetPlus[i][j][k] = selectCentrality(mcppWGenSetPlus,centBins[k], centBins[k+1]); 
            mcppWGenSubSetMinus[i][j][k] = selectCentrality(mcppWGenSetMinus,centBins[k], centBins[k+1]); 

            mcnpWGenSubSetPlus[i][j][k] = selectCentrality(mcnpWGenSetPlus,centBins[k], centBins[k+1]); 
            mcnpWGenSubSetMinus[i][j][k] = selectCentrality(mcnpWGenSetMinus,centBins[k], centBins[k+1]); 

            mcpnWGenSubSetPlus[i][j][k] = selectCentrality(mcpnWGenSetPlus,centBins[k], centBins[k+1]); 
            mcpnWGenSubSetMinus[i][j][k] = selectCentrality(mcpnWGenSetMinus,centBins[k], centBins[k+1]); 

            mcnnWGenSubSetPlus[i][j][k] = selectCentrality(mcnnWGenSetPlus,centBins[k], centBins[k+1]); 
            mcnnWGenSubSetMinus[i][j][k] = selectCentrality(mcnnWGenSetMinus,centBins[k], centBins[k+1]); 
        }

        std::cout << "Entries in pp subset for mu+: " << i << ":" << j << ":" << k << " = " << mcppWGenSubSetPlus[i][j][k]->sumEntries() << std::endl;
        std::cout << "Entries in pp subset for mu-: " << i << ":" << j << ":" << k << " = " << mcppWGenSubSetMinus[i][j][k]->sumEntries() << std::endl;

        std::cout << "Entries in np subset for mu+: " << i << ":" << j << ":" << k << " = " << mcnpWGenSubSetPlus[i][j][k]->sumEntries() << std::endl;
        std::cout << "Entries in np subset for mu-: " << i << ":" << j << ":" << k << " = " << mcnpWGenSubSetMinus[i][j][k]->sumEntries() << std::endl;

        std::cout << "Entries in pn subset for mu+: " << i << ":" << j << ":" << k << " = " << mcpnWGenSubSetPlus[i][j][k]->sumEntries() << std::endl;
        std::cout << "Entries in pn subset for mu-: " << i << ":" << j << ":" << k << " = " << mcpnWGenSubSetMinus[i][j][k]->sumEntries() << std::endl;

        std::cout << "Entries in nn subset for mu+: " << i << ":" << j << ":" << k << " = " << mcnnWGenSubSetPlus[i][j][k]->sumEntries() << std::endl;
        std::cout << "Entries in nn subset for mu-: " << i << ":" << j << ":" << k << " = " << mcnnWGenSubSetMinus[i][j][k]->sumEntries() << std::endl;

        // Generated subsets in fiducial region
        // pp
        // mu+
        mcppWFidSubSetPlus[i][j][k] = selectPtEtaCentrality(mcppWFidSetPlus ,ptBins[i],ptBins[i+1], etaGenBins[j],
            etaGenBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 
        std::cout << " " << i << ":" << j << ":" << k << " = " << mcppWFidSubSetPlus[i][j][k]->sumEntries()<< std::endl;
	// mu-
        mcppWFidSubSetMinus[i][j][k] = selectPtEtaCentrality(mcppWFidSetMinus ,ptBins[i],ptBins[i+1], etaGenBins[j],
            etaGenBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 
        std::cout << " " << i << ":" << j << ":" << k << " = " << mcppWFidSubSetMinus[i][j][k]->sumEntries()<< std::endl;

        // np
        // mu+
        mcnpWFidSubSetPlus[i][j][k] = selectPtEtaCentrality(mcnpWFidSetPlus ,ptBins[i],ptBins[i+1], etaGenBins[j],
            etaGenBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 
        std::cout << " " << i << ":" << j << ":" << k << " = " << mcnpWFidSubSetPlus[i][j][k]->sumEntries()<< std::endl;
	// mu-
        mcnpWFidSubSetMinus[i][j][k] = selectPtEtaCentrality(mcnpWFidSetMinus ,ptBins[i],ptBins[i+1], etaGenBins[j],
            etaGenBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 
        std::cout << " " << i << ":" << j << ":" << k << " = " << mcnpWFidSubSetMinus[i][j][k]->sumEntries()<< std::endl;

        // pn
        // mu+
        mcpnWFidSubSetPlus[i][j][k] = selectPtEtaCentrality(mcpnWFidSetPlus ,ptBins[i],ptBins[i+1], etaGenBins[j],
            etaGenBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 
        std::cout << " " << i << ":" << j << ":" << k << " = " << mcpnWFidSubSetPlus[i][j][k]->sumEntries()<< std::endl;
	// mu-
        mcpnWFidSubSetMinus[i][j][k] = selectPtEtaCentrality(mcpnWFidSetMinus ,ptBins[i],ptBins[i+1], etaGenBins[j],
            etaGenBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 
        std::cout << " " << i << ":" << j << ":" << k << " = " << mcpnWFidSubSetMinus[i][j][k]->sumEntries()<< std::endl;

        // nn
        // mu+
        mcnnWFidSubSetPlus[i][j][k] = selectPtEtaCentrality(mcnnWFidSetPlus ,ptBins[i],ptBins[i+1], etaGenBins[j],
            etaGenBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 
        std::cout << " " << i << ":" << j << ":" << k << " = " << mcnnWFidSubSetPlus[i][j][k]->sumEntries()<< std::endl;
	// mu-
        mcnnWFidSubSetMinus[i][j][k] = selectPtEtaCentrality(mcnnWFidSetMinus ,ptBins[i],ptBins[i+1], etaGenBins[j],
            etaGenBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 
        std::cout << " " << i << ":" << j << ":" << k << " = " << mcnnWFidSubSetMinus[i][j][k]->sumEntries()<< std::endl;

        // Reconstructed subsets in fiducial region
        // pp
        // mu+
        mcppSubSetWRecHiQualityWselPlus[i][j][k] = selectPtEtaCentrality(mcppSetWRecHiQualityWselPlus,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1], doMirrorEta,false);
        // mu-
        mcppSubSetWRecHiQualityWselMinus[i][j][k] = selectPtEtaCentrality(mcppSetWRecHiQualityWselMinus,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1], doMirrorEta,false);

        // np
        // mu+
        mcnpSubSetWRecHiQualityWselPlus[i][j][k] = selectPtEtaCentrality(mcnpSetWRecHiQualityWselPlus,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1], doMirrorEta,false);
        // mu-
        mcnpSubSetWRecHiQualityWselMinus[i][j][k] = selectPtEtaCentrality(mcnpSetWRecHiQualityWselMinus,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1], doMirrorEta,false);

        // pn
        // mu+
        mcpnSubSetWRecHiQualityWselPlus[i][j][k] = selectPtEtaCentrality(mcpnSetWRecHiQualityWselPlus,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1], doMirrorEta,false);
        // mu-
        mcpnSubSetWRecHiQualityWselMinus[i][j][k] = selectPtEtaCentrality(mcpnSetWRecHiQualityWselMinus,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1], doMirrorEta,false);

        // nn
        // mu+
        mcnnSubSetWRecHiQualityWselPlus[i][j][k] = selectPtEtaCentrality(mcnnSetWRecHiQualityWselPlus,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1], doMirrorEta,false);
        // mu-
        mcnnSubSetWRecHiQualityWselMinus[i][j][k] = selectPtEtaCentrality(mcnnSetWRecHiQualityWselMinus,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1], doMirrorEta,false);


          ///Get Mpt/Mt correlation
          if(doMptSigmaDown||doMptSigmaUp) {
                TString sMptMtCorrelation = "h2DMptMtCorrelationCent"; sMptMtCorrelation+=k; sMptMtCorrelation+="_pfx";
                std::cout << "Opening " << sMptMtCorrelation << " TProfilex histo." << std::endl;
                _pfxMptMtCorrelation = (TProfile*)_fMptMtCorrelation->Get(sMptMtCorrelation);
                if(_pfxMptMtCorrelation!=0) std::cout << sMptMtCorrelation << " opened." << std::endl;
                else exit(0);
          }

          if(mcppSubSetWRecHiQualityWselPlus[i][j][k]->sumEntries()>mcppWFidSubSetPlus[i][j][k]->sumEntries()){
            std::cout << "WARNING: Efficiency > 1.0. " << std::endl; exit(0); 
          }
          if(mcppSubSetWRecHiQualityWselMinus[i][j][k]->sumEntries()>mcppWFidSubSetMinus[i][j][k]->sumEntries()){
            std::cout << "WARNING: Efficiency > 1.0. " << std::endl; exit(0); 
	  }

          if(mcnpSubSetWRecHiQualityWselPlus[i][j][k]->sumEntries()>mcnpWFidSubSetPlus[i][j][k]->sumEntries()){
            std::cout << "WARNING: Efficiency > 1.0. " << std::endl; exit(0); 
          }
          if(mcnpSubSetWRecHiQualityWselMinus[i][j][k]->sumEntries()>mcnpWFidSubSetMinus[i][j][k]->sumEntries()){
            std::cout << "WARNING: Efficiency > 1.0. " << std::endl; exit(0); 
	  }

          if(mcpnSubSetWRecHiQualityWselPlus[i][j][k]->sumEntries()>mcpnWFidSubSetPlus[i][j][k]->sumEntries()){
            std::cout << "WARNING: Efficiency > 1.0. " << std::endl; exit(0); 
          }
          if(mcpnSubSetWRecHiQualityWselMinus[i][j][k]->sumEntries()>mcpnWFidSubSetMinus[i][j][k]->sumEntries()){
            std::cout << "WARNING: Efficiency > 1.0. " << std::endl; exit(0); 
	  }

          if(mcnnSubSetWRecHiQualityWselPlus[i][j][k]->sumEntries()>mcnnWFidSubSetPlus[i][j][k]->sumEntries()){
            std::cout << "WARNING: Efficiency > 1.0. " << std::endl; exit(0); 
          }
          if(mcnnSubSetWRecHiQualityWselMinus[i][j][k]->sumEntries()>mcnnWFidSubSetMinus[i][j][k]->sumEntries()){
            std::cout << "WARNING: Efficiency > 1.0. " << std::endl; exit(0); 
	  }
    } //k
  } //j
} //i

  TString spreadSheetName = "correctionFactorSpreadSheet.csv";
  std::ofstream spreadSheet;
  spreadSheet.open(spreadSheetName);

  TString spreadSheetNamePlus = "correctionFactorSpreadSheetPlus.csv";
  std::ofstream spreadSheetPlus;
  spreadSheetPlus.open(spreadSheetNamePlus);

  TString spreadSheetNameMinus = "correctionFactorSpreadSheetMinus.csv";
  std::ofstream spreadSheetMinus;
  spreadSheetMinus.open(spreadSheetNameMinus);

  //////////////////////////////////////////////////////
  ///Loop for plotting efficiencies as fcn of eta in
  ///each centrality bin
  //////////////////////////////////////////////////////
  if(doEta){
   std::cout << "Plotting eta dependence." << std::endl;
   for(int ipt=0; ipt<nPtBins; ipt++){

	for(int icent=0; icent<nCentralityBins; icent++){


            for(int ieta=0; ieta<nEtaBins; ieta++){

              double etaLow = etaBins.at(ieta);
              double etaUpp = etaBins.at(ieta+1); 

              //int index = icent*nEtaBins+ieta; 
	      double etaMed = etaBins[ieta]+(etaBins[ieta+1]-etaBins[ieta])/2.0;
              double binW = etaBins[ieta+1]-etaBins[ieta];
	      int index = -1; 
              
              if(!doMirrorEta){
                  int mappedEta;
                  ///hop over the crack region
                  ///in bin [-0.1,0.1]
                  if(etaMed==0.0) {std::cout << "Skipping crack region." << std::endl; continue; }
                  if(etaMed<0.0){
                      ///number of bins in absolute eta
                      ///used for trigger efficiency
                      mappedEta = indexNegativeEta(ieta,nEtaBins);
                      index = mappedEta*nCentralityBins+icent;
                  }
                   else {
                      mappedEta = indexPositiveEta(ieta,nEtaBins);
                      index = mappedEta*nCentralityBins+icent;
                   }
                }
	          else index = ieta*nCentralityBins+icent; 
              if(index<0.0) {
                std::cout << "WARNING: Invalid index. " << std::endl; break;
              }
                
              std::cout << "Eta median = " << etaMed << std::endl;
	          std::cout << "mu^{#pm}" << std::endl;


                // Calculate efficiencies
                // pp
	    	std::cout << "Calculating efficiency for mu^{+} in pp" << std::endl;
                calcEfficiency(mcppWGenSubSetPlus[ipt][ieta][icent],mcppWFidSubSetPlus[ipt][ieta][icent], mcppSubSetWRecHiQualityWselPlus[ipt][ieta][icent],
                    ieta, etaMed, binW, index,ipt, ieta,icent,grppAwEtaPlus,grppCwAwEtaPlus,
                    grppCwEffEtaPlus,spreadSheetPlus);
	    	std::cout << "Calculating efficiency for mu^{-} in pp" << std::endl;
                calcEfficiency(mcppWGenSubSetMinus[ipt][ieta][icent],mcppWFidSubSetMinus[ipt][ieta][icent], mcppSubSetWRecHiQualityWselMinus[ipt][ieta][icent],
                    ieta, etaMed, binW, index,ipt, ieta,icent,grppAwEtaMinus,grppCwAwEtaMinus,
                    grppCwEffEtaMinus,spreadSheetMinus);

	    	std::cout << "Calculating efficiency for mu^{pm} in pp" << std::endl;
                calcEfficiency(grppAwEtaPlus,grppAwEtaMinus,grppCwEffEtaPlus,grppCwEffEtaMinus,grppCwAwEtaPlus,grppCwAwEtaMinus, 
                    2.1120,1.2420,
                    ieta, etaMed, binW, index,ipt, ieta,icent,
                    grppAwEta,grppCwAwEta,grppCwEffEta);

                // np
	    	std::cout << "Calculating efficiency for mu^{+} in np" << std::endl;
                calcEfficiency(mcnpWGenSubSetPlus[ipt][ieta][icent],mcnpWFidSubSetPlus[ipt][ieta][icent], mcnpSubSetWRecHiQualityWselPlus[ipt][ieta][icent],
                    ieta, etaMed, binW, index,ipt, ieta,icent,grnpAwEtaPlus,grnpCwAwEtaPlus,
                    grnpCwEffEtaPlus,spreadSheetPlus);
	    	std::cout << "Calculating efficiency for mu^{-} in np" << std::endl;
                calcEfficiency(mcnpWGenSubSetMinus[ipt][ieta][icent],mcnpWFidSubSetMinus[ipt][ieta][icent], mcnpSubSetWRecHiQualityWselMinus[ipt][ieta][icent],
                    ieta, etaMed, binW, index,ipt, ieta,icent,grnpAwEtaMinus,grnpCwAwEtaMinus,
                    grnpCwEffEtaMinus,spreadSheetMinus);

	    	std::cout << "Calculating efficiency for mu^{pm} in np" << std::endl;
                calcEfficiency(grnpAwEtaPlus,grnpAwEtaMinus,grnpCwEffEtaPlus,grnpCwEffEtaMinus,grnpCwAwEtaPlus,grnpCwAwEtaMinus, 
                    1.6690,1.6540,
                    ieta, etaMed, binW, index,ipt, ieta,icent,
                    grnpAwEta,grnpCwAwEta,grnpCwEffEta);

                // pn
	    	std::cout << "Calculating efficiency for mu^{+} in pn" << std::endl;
                calcEfficiency(mcpnWGenSubSetPlus[ipt][ieta][icent],mcpnWFidSubSetPlus[ipt][ieta][icent], mcpnSubSetWRecHiQualityWselPlus[ipt][ieta][icent],
                    ieta, etaMed, binW, index,ipt, ieta,icent,grpnAwEtaPlus,grpnCwAwEtaPlus,
                    grpnCwEffEtaPlus,spreadSheetPlus);
	    	std::cout << "Calculating efficiency for mu^{-} in pn" << std::endl;
                calcEfficiency(mcpnWGenSubSetMinus[ipt][ieta][icent],mcpnWFidSubSetMinus[ipt][ieta][icent], mcpnSubSetWRecHiQualityWselMinus[ipt][ieta][icent],
                    ieta, etaMed, binW, index,ipt, ieta,icent,grpnAwEtaMinus,grpnCwAwEtaMinus,
                    grpnCwEffEtaMinus,spreadSheetMinus);

	    	std::cout << "Calculating efficiency for mu^{pm} in pn" << std::endl;
                calcEfficiency(grpnAwEtaPlus,grpnAwEtaMinus,grpnCwEffEtaPlus,grpnCwEffEtaMinus,grpnCwAwEtaPlus,grpnCwAwEtaMinus, 
                    1.6700,1.6540,
                    ieta, etaMed, binW, index,ipt, ieta,icent,
                    grpnAwEta,grpnCwAwEta,grpnCwEffEta);

                // nn
	    	std::cout << "Calculating efficiency for mu^{+} in nn" << std::endl;
                calcEfficiency(mcnnWGenSubSetPlus[ipt][ieta][icent],mcnnWFidSubSetPlus[ipt][ieta][icent], mcnnSubSetWRecHiQualityWselPlus[ipt][ieta][icent],
                    ieta, etaMed, binW, index,ipt, ieta,icent,grnnAwEtaPlus,grnnCwAwEtaPlus,
                    grnnCwEffEtaPlus,spreadSheetPlus);
	    	std::cout << "Calculating efficiency for mu^{-} in nn" << std::endl;
                calcEfficiency(mcnnWGenSubSetMinus[ipt][ieta][icent],mcnnWFidSubSetMinus[ipt][ieta][icent], mcnnSubSetWRecHiQualityWselMinus[ipt][ieta][icent],
                    ieta, etaMed, binW, index,ipt, ieta,icent,grnnAwEtaMinus,grnnCwAwEtaMinus,
                    grnnCwEffEtaMinus,spreadSheetMinus);

	    	std::cout << "Calculating efficiency for mu^{pm} in nn" << std::endl;
                calcEfficiency(grnnAwEtaPlus,grnnAwEtaMinus,grnnCwEffEtaPlus,grnnCwEffEtaMinus,grnnCwAwEtaPlus,grnnCwAwEtaMinus, 
                    1.2530,2.0920,
                    ieta, etaMed, binW, index,ipt, ieta,icent,
                    grnnAwEta,grnnCwAwEta,grnnCwEffEta);
            
                // Calculate Cw weighted by pp,np,pn,nn collision frequencies (15.5%, 47.8%, and 36.7%)
                std::cout << "Calculating weighted efficiency for mu+." << std::endl;
                // Cw
                calcWeightedEfficiency(ieta, etaMed, binW, index,ipt, ieta,icent, 
					grppCwEffEtaPlus, grnpCwEffEtaPlus, grpnCwEffEtaPlus, grnnCwEffEtaPlus, grWtdCwEffEtaPlus);
                //Cw*Aw
                calcWeightedEfficiency(ieta, etaMed, binW, index,ipt, ieta,icent, 
					grppCwAwEtaPlus, grnpCwAwEtaPlus, grpnCwAwEtaPlus, grnnCwAwEtaPlus, grWtdCwAwEffEtaPlus);

                std::cout << "Calculating weighted efficiency for mu-." << std::endl;
                calcWeightedEfficiency(ieta, etaMed, binW, index,ipt, ieta,icent, 
					grppCwEffEtaMinus, grnpCwEffEtaMinus, grpnCwEffEtaMinus, grnnCwEffEtaMinus, grWtdCwEffEtaMinus);
                calcWeightedEfficiency(ieta, etaMed, binW, index,ipt, ieta,icent, 
					grppCwAwEtaMinus, grnpCwAwEtaMinus, grpnCwAwEtaMinus, grnnCwAwEtaMinus, grWtdCwAwEffEtaMinus);

                std::cout << "Calculating weighted efficiency for mu+-." << std::endl;
                calcWeightedEfficiency(ieta, etaMed, binW, index,ipt, ieta,icent, 
					grppCwEffEta, grnpCwEffEta, grpnCwEffEta, grnnCwEffEta, grWtdCwEffEta);
                calcWeightedEfficiency(ieta, etaMed, binW, index,ipt, ieta,icent, 
					grppCwAwEta, grnpCwAwEta, grpnCwAwEta, grnnCwAwEta, grWtdCwAwEffEta);
            
            
            }//ieta

            // Output TGraph names
            // pp
            TString sppEffNameAwPlus = "grpp_AwPlusEtaDistroCent"; sppEffNameAwPlus+=icent;
            TString sppEffNameCwAwPlus = "grpp_CwAwPlusEtaDistroCent"; sppEffNameCwAwPlus+=icent;
            TString sppEffNameCwPlus = "grpp_CwEtaDistroPlusCent"; sppEffNameCwPlus+=icent;
            TString sppEffNameAwMinus = "grpp_AwMinusEtaDistroCent"; sppEffNameAwMinus+=icent;
            TString sppEffNameCwAwMinus = "grpp_CwAwMinusEtaDistroCent"; sppEffNameCwAwMinus+=icent;
            TString sppEffNameCwMinus = "grpp_CwEtaDistroMinusCent"; sppEffNameCwMinus+=icent;

            // np
            TString snpEffNameAwPlus = "grnp_AwPlusEtaDistroCent"; snpEffNameAwPlus+=icent;
            TString snpEffNameCwAwPlus = "grnp_CwAwPlusEtaDistroCent"; snpEffNameCwAwPlus+=icent;
            TString snpEffNameCwPlus = "grnp_CwEtaDistroPlusCent"; snpEffNameCwPlus+=icent;
            TString snpEffNameAwMinus = "grnp_AwMinusEtaDistroCent"; snpEffNameAwMinus+=icent;
            TString snpEffNameCwAwMinus = "grnp_CwAwMinusEtaDistroCent"; snpEffNameCwAwMinus+=icent;
            TString snpEffNameCwMinus = "grnp_CwEtaDistroMinusCent"; snpEffNameCwMinus+=icent;

            // pn
            TString spnEffNameAwPlus = "grpn_AwPlusEtaDistroCent"; spnEffNameAwPlus+=icent;
            TString spnEffNameCwAwPlus = "grpn_CwAwPlusEtaDistroCent"; spnEffNameCwAwPlus+=icent;
            TString spnEffNameCwPlus = "grpn_CwEtaDistroPlusCent"; spnEffNameCwPlus+=icent;
            TString spnEffNameAwMinus = "grpn_AwMinusEtaDistroCent"; spnEffNameAwMinus+=icent;
            TString spnEffNameCwAwMinus = "grpn_CwAwMinusEtaDistroCent"; spnEffNameCwAwMinus+=icent;
            TString spnEffNameCwMinus = "grpn_CwEtaDistroMinusCent"; spnEffNameCwMinus+=icent;

            // nn
            TString snnEffNameAwPlus = "grnn_AwPlusEtaDistroCent"; snnEffNameAwPlus+=icent;
            TString snnEffNameCwAwPlus = "grnn_CwAwPlusEtaDistroCent"; snnEffNameCwAwPlus+=icent;
            TString snnEffNameCwPlus = "grnn_CwEtaDistroPlusCent"; snnEffNameCwPlus+=icent;
            TString snnEffNameAwMinus = "grnn_AwMinusEtaDistroCent"; snnEffNameAwMinus+=icent;
            TString snnEffNameCwAwMinus = "grnn_CwAwMinusEtaDistroCent"; snnEffNameCwAwMinus+=icent;
            TString snnEffNameCwMinus = "grnn_CwEtaDistroMinusCent"; snnEffNameCwMinus+=icent;

	    TString sEffNameCwPlus = "grCwEtaDistroPlusCent"; sEffNameCwPlus+=icent;
	    TString sEffNameCwMinus = "grCwEtaDistroMinusCent"; sEffNameCwMinus+=icent;
	    TString sEffNameCw = "grCwEtaDistroCent"; sEffNameCw+=icent;

	    TString sEffNameCwAwPlus = "grCwAwEtaDistroPlusCent"; sEffNameCwAwPlus+=icent;
	    TString sEffNameCwAwMinus = "grCwAwEtaDistroMinusCent"; sEffNameCwAwMinus+=icent;
	    TString sEffNameCwAw = "grCwAwEtaDistroCent"; sEffNameCwAw+=icent;

            if(doCharge){
                Write(outFile,grWtdCwEffEta,sEffNameCw);
                Write(outFile,grWtdCwAwEffEta,sEffNameCwAw);
		        Write(outFile,grWtdCwEffEtaPlus,sEffNameCwPlus);
		        Write(outFile,grWtdCwEffEtaMinus,sEffNameCwMinus);
                Write(outFile,grWtdCwAwEffEtaPlus,sEffNameCwAwPlus);
		        Write(outFile,grWtdCwAwEffEtaMinus,sEffNameCwAwMinus);

                Write(outFile,grppAwEtaPlus,sppEffNameAwPlus);
                Write(outFile,grppCwAwEtaPlus,sppEffNameCwAwPlus);
                Write(outFile,grppCwEffEtaPlus,sppEffNameCwPlus);
                Write(outFile,grppAwEtaMinus,sppEffNameAwMinus);
                Write(outFile,grppCwAwEtaMinus,sppEffNameCwAwMinus);
                Write(outFile,grppCwEffEtaMinus,sppEffNameCwMinus);

                Write(outFile,grnpAwEtaPlus,snpEffNameAwPlus);
                Write(outFile,grnpCwAwEtaPlus,snpEffNameCwAwPlus);
                Write(outFile,grnpCwEffEtaPlus,snpEffNameCwPlus);
                Write(outFile,grnpAwEtaMinus,snpEffNameAwMinus);
                Write(outFile,grnpCwAwEtaMinus,snpEffNameCwAwMinus);
                Write(outFile,grnpCwEffEtaMinus,snpEffNameCwMinus);

                Write(outFile,grpnAwEtaPlus,spnEffNameAwPlus);
                Write(outFile,grpnCwAwEtaPlus,spnEffNameCwAwPlus);
                Write(outFile,grpnCwEffEtaPlus,spnEffNameCwPlus);
                Write(outFile,grpnAwEtaMinus,spnEffNameAwMinus);
                Write(outFile,grpnCwAwEtaMinus,spnEffNameCwAwMinus);
                Write(outFile,grpnCwEffEtaMinus,spnEffNameCwMinus);

                Write(outFile,grnnAwEtaPlus,snnEffNameAwPlus);
                Write(outFile,grnnCwAwEtaPlus,snnEffNameCwAwPlus);
                Write(outFile,grnnCwEffEtaPlus,snnEffNameCwPlus);
                Write(outFile,grnnAwEtaMinus,snnEffNameAwMinus);
                Write(outFile,grnnCwAwEtaMinus,snnEffNameCwAwMinus);
                Write(outFile,grnnCwEffEtaMinus,snnEffNameCwMinus);
              }
          }//icent
      }//ipt
} //doEta

  //////////////////////////////////////////////////////
  ///Loop for plotting efficiencies as fcn of Npart in
  ///each eta bin
  //////////////////////////////////////////////////////
  if(doCentrality){
      std::cout << "Plotting centrality dependence." << std::endl;
      for(int ipt=0; ipt<nPtBins; ipt++){

          for(int ieta=0; ieta<nEtaBins; ieta++){

            double etabinLo = etaBins.at(ieta);
            double etabinUp = etaBins.at(ieta+1); 
            double xEta = etabinLo+(etabinUp-etabinLo)/2.0; 



	    for(int icent=0; icent<nCentralityBins; icent++){

	       std::cout << "Plotting ipt:ieta:icent = " << ipt << ":" << ieta << ":" << icent << std::endl;
	       double centbinLo = centBins.at(icent);
	       double centbinUp = centBins.at(icent+1); 

	       double npart = npartBins[icent]; 
	       double binW = npartBins[icent+1]-npartBins[icent];
				
	       double etaMed = etaBins[ieta]+(etaBins[ieta+1]-etaBins[ieta])/2.0;
	       int index = -1; 
              
                if(!doMirrorEta){
                  int mappedEta;
                  ///hop over the crack region
                  ///in bin [-0.1,0.1]
                  if(etaMed==0.0) {std::cout << "Skipping crack region." << std::endl; continue; }
                  if(etaMed<0.0){
                      ///number of bins in absolute eta
                      ///used for trigger efficiency
                      mappedEta = indexNegativeEta(ieta,nEtaBins);
                      index = mappedEta*nCentralityBins+icent;
                  }
                  else {
                      mappedEta = indexPositiveEta(ieta,nEtaBins);
                      index = mappedEta*nCentralityBins+icent;
                  }
                }
	        else {
                  index = ieta*nCentralityBins+icent; 
                }
                if(index<0.0) {
                  std::cout << "WARNING: Invalid index. " << std::endl; break;
                }
                

                // Calculate efficiencies
                // pp
	    	std::cout << "Calculating efficiency for mu^{+} in pp" << std::endl;
                calcEfficiency(mcppWGenSubSetPlus[ipt][ieta][icent],mcppWFidSubSetPlus[ipt][ieta][icent], mcppSubSetWRecHiQualityWselPlus[ipt][ieta][icent],
                    icent, npart, binW, index,ipt, ieta,icent,grppAwCentPlus,grppCwAwCentPlus,
                    grppCwEffCentPlus,spreadSheetPlus,true);
	    	std::cout << "Calculating efficiency for mu^{-} in pp" << std::endl;
                calcEfficiency(mcppWGenSubSetMinus[ipt][ieta][icent],mcppWFidSubSetMinus[ipt][ieta][icent], mcppSubSetWRecHiQualityWselMinus[ipt][ieta][icent],
                    icent, npart, binW, index,ipt, ieta,icent,grppAwCentMinus,grppCwAwCentMinus,
                    grppCwEffCentMinus,spreadSheetMinus,true);

	    	std::cout << "Calculating efficiency for mu^{pm} in pp" << std::endl;
                calcEfficiency(grppAwCentPlus,grppAwCentMinus,grppCwEffCentPlus,grppCwEffCentMinus,grppCwAwCentPlus,grppCwAwCentMinus, 
                    2.1120,1.2420,
                    icent, npart, binW, index,ipt, ieta,icent,
                    grppAwCent,grppCwAwCent,grppCwEffCent);

                // np
	    	std::cout << "Calculating efficiency for mu^{+} in np" << std::endl;
                calcEfficiency(mcnpWGenSubSetPlus[ipt][ieta][icent],mcnpWFidSubSetPlus[ipt][ieta][icent], mcnpSubSetWRecHiQualityWselPlus[ipt][ieta][icent],
                    icent, npart, binW, index,ipt, ieta,icent,grnpAwCentPlus,grnpCwAwCentPlus,
                    grnpCwEffCentPlus,spreadSheetPlus,true);
	    	std::cout << "Calculating efficiency for mu^{-} in np" << std::endl;
                calcEfficiency(mcnpWGenSubSetMinus[ipt][ieta][icent],mcnpWFidSubSetMinus[ipt][ieta][icent], mcnpSubSetWRecHiQualityWselMinus[ipt][ieta][icent],
                    icent, npart, binW, index,ipt, ieta,icent,grnpAwCentMinus,grnpCwAwCentMinus,
                    grnpCwEffCentMinus,spreadSheetMinus,true);

	    	std::cout << "Calculating efficiency for mu^{pm} in np" << std::endl;
                calcEfficiency(grnpAwCentPlus,grnpAwCentMinus,grnpCwEffCentPlus,grnpCwEffCentMinus,grnpCwAwCentPlus,grnpCwAwCentMinus, 
                    1.6690,1.6540,
                    icent, npart, binW, index,ipt, ieta,icent,
                    grnpAwCent,grnpCwAwCent,grnpCwEffCent);

                // pn
	    	std::cout << "Calculating efficiency for mu^{+} in pn" << std::endl;
                calcEfficiency(mcpnWGenSubSetPlus[ipt][ieta][icent],mcpnWFidSubSetPlus[ipt][ieta][icent], mcpnSubSetWRecHiQualityWselPlus[ipt][ieta][icent],
                    icent, npart, binW, index,ipt, ieta,icent,grpnAwCentPlus,grpnCwAwCentPlus,
                    grpnCwEffCentPlus,spreadSheetPlus,true);
	    	std::cout << "Calculating efficiency for mu^{-} in pn" << std::endl;
                calcEfficiency(mcpnWGenSubSetMinus[ipt][ieta][icent],mcpnWFidSubSetMinus[ipt][ieta][icent], mcpnSubSetWRecHiQualityWselMinus[ipt][ieta][icent],
                    icent, npart, binW, index,ipt, ieta,icent,grpnAwCentMinus,grpnCwAwCentMinus,
                    grpnCwEffCentMinus,spreadSheetMinus,true);

	    	std::cout << "Calculating efficiency for mu^{pm} in pn" << std::endl;
                calcEfficiency(grpnAwCentPlus,grpnAwCentMinus,grpnCwEffCentPlus,grpnCwEffCentMinus,grpnCwAwCentPlus,grpnCwAwCentMinus, 
                    1.6700,1.6540,
                    icent, npart, binW, index,ipt, ieta,icent,
                    grpnAwCent,grpnCwAwCent,grpnCwEffCent);

                // nn
	    	std::cout << "Calculating efficiency for mu^{+} in nn" << std::endl;
                calcEfficiency(mcnnWGenSubSetPlus[ipt][ieta][icent],mcnnWFidSubSetPlus[ipt][ieta][icent], mcnnSubSetWRecHiQualityWselPlus[ipt][ieta][icent],
                    icent, npart, binW, index,ipt, ieta,icent,grnnAwCentPlus,grnnCwAwCentPlus,
                    grnnCwEffCentPlus,spreadSheetPlus,true);
	    	std::cout << "Calculating efficiency for mu^{-} in nn" << std::endl;
                calcEfficiency(mcnnWGenSubSetMinus[ipt][ieta][icent],mcnnWFidSubSetMinus[ipt][ieta][icent], mcnnSubSetWRecHiQualityWselMinus[ipt][ieta][icent],
                    icent, npart, binW, index,ipt, ieta,icent,grnnAwCentMinus,grnnCwAwCentMinus,
                    grnnCwEffCentMinus,spreadSheetMinus,true);

	    	std::cout << "Calculating efficiency for mu^{pm} in nn" << std::endl;
                calcEfficiency(grnnAwCentPlus,grnnAwCentMinus,grnnCwEffCentPlus,grnnCwEffCentMinus,grnnCwAwCentPlus,grnnCwAwCentMinus, 
                    1.2530,2.0920,
                    icent, npart, binW, index,ipt, ieta,icent,
                    grnnAwCent,grnnCwAwCent,grnnCwEffCent);

                // Calculate Cw weighted by pp,np,pn,nn collision frequencies (15.5%, 47.8%, and 36.7%)
                std::cout << "Calculating weighted efficiency for mu+." << std::endl;
                calcWeightedEfficiency(icent, npart, binW, index,ipt, ieta,icent, 
					grppCwEffCentPlus, grnpCwEffCentPlus, grpnCwEffCentPlus, grnnCwEffCentPlus, grWtdCwEffCentPlus);
                calcWeightedEfficiency(icent, npart, binW, index,ipt, ieta,icent, 
					grppCwAwCentPlus, grnpCwAwCentPlus, grpnCwAwCentPlus, grnnCwAwCentPlus, grWtdCwAwEffCentPlus);

                std::cout << "Calculating weighted efficiency for mu-." << std::endl;
                calcWeightedEfficiency(icent, npart, binW, index,ipt, ieta,icent, 
					grppCwEffCentMinus, grnpCwEffCentMinus, grpnCwEffCentMinus, grnnCwEffCentMinus, grWtdCwEffCentMinus);
                calcWeightedEfficiency(icent, npart, binW, index,ipt, ieta,icent, 
					grppCwAwCentMinus, grnpCwAwCentMinus, grpnCwAwCentMinus, grnnCwAwCentMinus, grWtdCwAwEffCentMinus);

                std::cout << "Calculating weighted efficiency for mu+-." << std::endl;
                calcWeightedEfficiency(icent, npart, binW, index,ipt, ieta,icent, 
					grppCwEffCent, grnpCwEffCent, grpnCwEffCent, grnnCwEffCent, grWtdCwEffCent);
                calcWeightedEfficiency(icent, npart, binW, index,ipt, ieta,icent, 
					grppCwAwCent, grnpCwAwCent, grpnCwAwCent, grnnCwAwCent, grWtdCwAwEffCent);
            
             } //icent

            // Output TGraph names
            // pp
            TString sppEffNameAwPlus = "grpp_AwPlusNpartDistroEta"; sppEffNameAwPlus+=ieta;
            TString sppEffNameCwAwPlus = "grpp_CwAwPlusNpartDistroEta"; sppEffNameCwAwPlus+=ieta;
            TString sppEffNameCwPlus = "grpp_CwNpartDistroPlusEta"; sppEffNameCwPlus+=ieta;
            TString sppEffNameAwMinus = "grpp_AwMinusNpartDistroEta"; sppEffNameAwMinus+=ieta;
            TString sppEffNameCwAwMinus = "grpp_CwAwMinusNpartDistroEta"; sppEffNameCwAwMinus+=ieta;
            TString sppEffNameCwMinus = "grpp_CwNpartDistroMinusEta"; sppEffNameCwMinus+=ieta;

            // np
            TString snpEffNameAwPlus = "grnp_AwPlusNpartDistroEta"; snpEffNameAwPlus+=ieta;
            TString snpEffNameCwAwPlus = "grnp_CwAwPlusNpartDistroEta"; snpEffNameCwAwPlus+=ieta;
            TString snpEffNameCwPlus = "grnp_CwNpartDistroPlusEta"; snpEffNameCwPlus+=ieta;
            TString snpEffNameAwMinus = "grnp_AwMinusNpartDistroEta"; snpEffNameAwMinus+=ieta;
            TString snpEffNameCwAwMinus = "grnp_CwAwMinusNpartDistroEta"; snpEffNameCwAwMinus+=ieta;
            TString snpEffNameCwMinus = "grnp_CwNpartDistroMinusEta"; snpEffNameCwMinus+=ieta;

            // pn
            TString spnEffNameAwPlus = "grpn_AwPlusNpartDistroEta"; spnEffNameAwPlus+=ieta;
            TString spnEffNameCwAwPlus = "grpn_CwAwPlusNpartDistroEta"; spnEffNameCwAwPlus+=ieta;
            TString spnEffNameCwPlus = "grpn_CwNpartDistroPlusEta"; spnEffNameCwPlus+=ieta;
            TString spnEffNameAwMinus = "grpn_AwMinusNpartDistroEta"; spnEffNameAwMinus+=ieta;
            TString spnEffNameCwAwMinus = "grpn_CwAwMinusNpartDistroEta"; spnEffNameCwAwMinus+=ieta;
            TString spnEffNameCwMinus = "grpn_CwNpartDistroMinusEta"; spnEffNameCwMinus+=ieta;

            // nn
            TString snnEffNameAwPlus = "grnn_AwPlusNpartDistroEta"; snnEffNameAwPlus+=ieta;
            TString snnEffNameCwAwPlus = "grnn_CwAwPlusNpartDistroEta"; snnEffNameCwAwPlus+=ieta;
            TString snnEffNameCwPlus = "grnn_CwNpartDistroPlusEta"; snnEffNameCwPlus+=ieta;
            TString snnEffNameAwMinus = "grnn_AwMinusNpartDistroEta"; snnEffNameAwMinus+=ieta;
            TString snnEffNameCwAwMinus = "grnn_CwAwMinusNpartDistroEta"; snnEffNameCwAwMinus+=ieta;
            TString snnEffNameCwMinus = "grnn_CwNpartDistroMinusEta"; snnEffNameCwMinus+=ieta;

	    TString sEffNameCwPlus = "grCwNpartDistroPlusEta"; sEffNameCwPlus+=ieta;
	    TString sEffNameCwMinus = "grCwNpartDistroMinusEta"; sEffNameCwMinus+=ieta;
	    TString sEffNameCw = "grCwNpartDistroEta"; sEffNameCw+=ieta;

	    TString sEffNameCwAwPlus = "grCwAwNpartDistroPlusEta"; sEffNameCwAwPlus+=ieta;
	    TString sEffNameCwAwMinus = "grCwAwNpartDistroMinusEta"; sEffNameCwAwMinus+=ieta;
	    TString sEffNameCwAw = "grCwAwNpartDistroEta"; sEffNameCwAw+=ieta;


		        Write(outFile,grWtdCwAwEffCent,sEffNameCwAw);
                Write(outFile,grWtdCwEffCent,sEffNameCw);
                Write(outFile,grWtdCwEffCentPlus,sEffNameCwPlus);
                Write(outFile,grWtdCwEffCentMinus,sEffNameCwMinus);
                Write(outFile,grWtdCwAwEffCentPlus,sEffNameCwAwPlus);
		        Write(outFile,grWtdCwAwEffCentMinus,sEffNameCwAwMinus);

                Write(outFile,grppAwCentPlus,sppEffNameAwPlus);
                Write(outFile,grppCwAwCentPlus,sppEffNameCwAwPlus);
                Write(outFile,grppCwEffCentPlus,sppEffNameCwPlus);
                Write(outFile,grppAwCentMinus,sppEffNameAwMinus);
                Write(outFile,grppCwAwCentMinus,sppEffNameCwAwMinus);
                Write(outFile,grppCwEffCentMinus,sppEffNameCwMinus);

                Write(outFile,grnpAwCentPlus,snpEffNameAwPlus);
                Write(outFile,grnpCwAwCentPlus,snpEffNameCwAwPlus);
                Write(outFile,grnpCwEffCentPlus,snpEffNameCwPlus);
                Write(outFile,grnpAwCentMinus,snpEffNameAwMinus);
                Write(outFile,grnpCwAwCentMinus,snpEffNameCwAwMinus);
                Write(outFile,grnpCwEffCentMinus,snpEffNameCwMinus);

                Write(outFile,grpnAwCentPlus,spnEffNameAwPlus);
                Write(outFile,grpnCwAwCentPlus,spnEffNameCwAwPlus);
                Write(outFile,grpnCwEffCentPlus,spnEffNameCwPlus);
                Write(outFile,grpnAwCentMinus,spnEffNameAwMinus);
                Write(outFile,grpnCwAwCentMinus,spnEffNameCwAwMinus);
                Write(outFile,grpnCwEffCentMinus,spnEffNameCwMinus);

                Write(outFile,grnnAwCentPlus,snnEffNameAwPlus);
                Write(outFile,grnnCwAwCentPlus,snnEffNameCwAwPlus);
                Write(outFile,grnnCwEffCentPlus,snnEffNameCwPlus);
                Write(outFile,grnnAwCentMinus,snnEffNameAwMinus);
                Write(outFile,grnnCwAwCentMinus,snnEffNameCwAwMinus);
                Write(outFile,grnnCwEffCentMinus,snnEffNameCwMinus);

	   } //ieta

	} //ipt
   } //doCentrality

    std::cout << "Closing spreadsheet..." << std::endl;
    spreadSheet.close();
    spreadSheetPlus.close();
    spreadSheetMinus.close();
    std::cout << "Done." << std::endl;
} 

int main(){
    CorrectionFactorsCw_Reweighting();
    std::cout << "Done running macro." << std::endl;
    return 0;
}
