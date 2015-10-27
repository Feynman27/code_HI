#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooDataSet.h"
#include "RooBinning.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooFormulaVar.h"

#include "TDirectory.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath> 

#include "AtlasUtils.C"
#include "AtlasStyle.C"

#include "TriggerEfficiencies.C"
#include "WPlotterHelper.C"
#include "correctionFactorsDep.C"


using namespace RooFit;

///////////////////////////////////////////
//writeToSpreadsheet
//////////////////////////////////////////
void writeToSpreadsheet(std::ostream& outputFile, int ipt, int ieta, int icent, double scaleFactor, double errStatLo,double errStatUp){
	outputFile << ipt << "," << ieta << "," << icent << "," << scaleFactor << "," << errStatLo << "," << errStatUp << std::endl;
}


/////////////////////////////////////////////
//plot residual between fit and Cw(Npart) distribution
//per eta bin
/////////////////////////////////////////////
void plotCwResidual(TGraphErrors* const gr, TGraph* const grResid){

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
plotGeneratorDistros(TH1F* hGen,TH1F* hFid,int ich,TString sSel){

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

	myText(0.196,0.862,kBlack,sSel);

	leg->Draw();

        //cGen->SetLogy(true); hGenc->GetYaxis()->SetRangeUser(0.1,hGenc->GetMaximum()*7.0e2); cGen->Update();
	cGen->Print(sSel+",GeneratorDistro.pdf");
    TString name = "GeneratorDistro"; name+=ich;
	cGen->Print(name+".root");

}

///////////////////////////////////////////////////////////////////////////////
//plot distro of fiducial kinematic 
///////////////////////////////////////////////////////////////////////////////
plotFiducialDistros(TH1F* hFid,TH1F* hRec,TString sSel){

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

	myText(0.196,0.862,kBlack,sSel);

	leg->Draw();

        cFid->SetLogy(true); hFidc->GetYaxis()->SetRangeUser(0.1,hFidc->GetMaximum()*7.0e2); cFid->Update();
	cFid->Print(sSel+",FiducialDistro.pdf");

}

///////////////////////////////////////////////////////////////////////////////
//plot difference in charge Cw as fcn of centrality 
///////////////////////////////////////////////////////////////////////////////
plotChargeDiffCwNpart(TGraph* const grDiff, TGraphErrors* const grPlus, TGraphErrors* const grMinus, int ieta, int icent, double npart){

	double* yPlus = grPlus->GetY();
	double* yMinus = grMinus->GetY();

	double yTempPlus = yPlus[icent]; 
	double yTempMinus = yMinus[icent]; 
	double yTempPlus = yPlus[icent]; 
	double yTempMinus = yMinus[icent]; 

	double cWdiff = fabs(yTempPlus-yTempMinus);
	std::cout << "Difference in Cw at Eta Bin " << ieta << " and Centrality Bin " << icent << " = " << cWdiff << std::endl;
	grDiff->SetPoint(icent,npart,cWdiff);

}
///////////////////////////////////////////////////////////////////////////////
//plot difference in charge Cw as fcn of eta 
///////////////////////////////////////////////////////////////////////////////
plotChargeDiffCwEta(TGraph* const grDiff, TGraphErrors* const grPlus, TGraphErrors* const grMinus, int ieta, int icent, double xEta) {

	double* yPlus = grPlus->GetY();
	double* yMinus = grMinus->GetY();

	double yTempPlus = yPlus[ieta]; 
	double yTempMinus = yMinus[ieta]; 
	double yTempPlus = yPlus[ieta]; 
	double yTempMinus = yMinus[ieta]; 

	double cWdiff = fabs(yTempPlus-yTempMinus);
	std::cout << "Difference in Cw at Eta Bin " << ieta << " and Centrality Bin " << icent << " = " << cWdiff << std::endl;
	grDiff->SetPoint(ieta,xEta,cWdiff);

}
///////////////////////////////////////////////////////////////////////////////
//plot Cw as a function of centrality 
///////////////////////////////////////////////////////////////////////////////
plotCwCentralityDep(TGraphErrors* const grCent,int icent, double npart,float Cw, float CwErr){

	std::cout << "Plotting Cw centrality dependance at centrality bin: " << icent << ": Npart " << npart << std::endl;
	grCent->SetPoint(icent,npart,Cw);
	grCent->SetPointError(icent,0.0,CwErr);

}

///////////////////////////////////////////////////////////////////////////////
//plot Cw as a function of eta 
///////////////////////////////////////////////////////////////////////////////
plotCwEtaDep(TGraphErrors* const grEta,int ieta, float xEta, float Cw, float CwErr){

	std::cout << "Plotting Cw eta dependance at eta bin: " << ieta << std::endl;
	grEta->SetPoint(ieta,xEta,Cw);
	grEta->SetPointError(ieta,0.0,CwErr);

}
///////////////////////////////////////////////////////////////////////////////
//plot difference in charge A0,A1,and A2
///////////////////////////////////////////////////////////////////////////////
plotChargeDiffAx(TGraph* const grDiff, TGraphErrors* const grPlus, TGraphErrors* const grMinus, int ieta, double xEta){

	double* yPlus = grPlus->GetY();
	double* yMinus = grMinus->GetY();

	double yTempPlus = yPlus[ieta]; 
	double yTempMinus = yMinus[ieta]; 

	double Axdiff = fabs(yTempPlus-yTempMinus);
	std::cout << "Difference in Cw at Eta Bin " << ieta << " = " << Axdiff << std::endl;
	grDiff->SetPoint(ieta,xEta,Axdiff);

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

///////////////////////////////////////////////////////////////////////////////
//Function for calculating and plotting efficiencies to a TGraph
///////////////////////////////////////////////////////////////////////////////
void calcEfficiency(RooDataSet* mcWGenSet, RooDataSet* mcWFidSet,RooDataSet* mcWRecSet0, RooDataSet* mcWRecSet1,
                RooDataSet* mcWRecSet2,RooDataSet* mcWRecSet3,RooDataSet* mcWRecSet4,RooDataSet* mcWRecSet5,
                int i, double x, double binW, int index, int ipt, int ieta, int icent,                 
		        TGraphAsymmErrors* const grWRecEff0, TGraphAsymmErrors* const grWRecEff1,TGraphAsymmErrors* const grWRecEff2,
                TGraphAsymmErrors* const grWRecEff3,TGraphAsymmErrors* const grWRecEff4,TGraphAsymmErrors* const grWRecEff5, 
                std::ostream& output,bool doWriteToSpreadsheet=false, bool doSingleTriggerEff=false){

    //number of Wmu in this eta and centrality class
    int tGen = (int)mcWGenSet->numEntries();
    std::cout << "tGen = " << tGen << std::endl;

    //number of muons at generator level in kinematic fiducial region within this eta and centrality class
    int tGenInFiducial = (int)mcWFidSet->numEntries();
    std::cout << "tGenInFiducial = " << tGenInFiducial << std::endl;

	//number of reconstructed signal muons after pre-selection within the fiducial region in this eta and centrality class
//	double pRecOnly = (double)mcWRecOnlySet->numEntries() ;
//    std::cout << "pRecOnly = " << pRecOnly << std::endl;

    ///NoQualityWselNoIsoNoZveto
    int pRec0 = (int)mcWRecSet0->numEntries() ;
    std::cout << "pRec0 = " << pRec0 << std::endl;

    ///Calculate efficiency using Bayesian errors 
    ///with a 68.3% C.I.
    double recEff0, errRec0Hi, errRec0Lo;
    Efficiency(pRec0,tGenInFiducial,0.683,recEff0,errRec0Lo,errRec0Hi);
    if(recEff0<0.0) Warning("calcEfficiency","Efficiency  = %d. Expect wrong results", recEff0);
    if(recEff0>1.0) Warning("calcEfficiency","Efficiency  = %d. Expect wrong results", recEff0);

    ///Return +- value from mean; therefore; must subtract
    ///to get absolute errors
    errRec0Hi = errRec0Hi - recEff0;
    errRec0Lo = recEff0 - errRec0Lo;
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
    }
    ///Apply trigger efficiency to reconstructed part
    double eff0 = recEff0*trigTemp;

    ///Propagate the error
    double effStatErrHi0 = sqrt( TMath::Power(errRec0Hi/recEff0,2) + TMath::Power(trigStatErr,2) )*eff0;
    double effStatErrLo0 = sqrt( TMath::Power(errRec0Lo/recEff0,2) + TMath::Power(trigStatErr,2) )*eff0;

    std::cout << "Efficiency0: " << recEff0 << " + " << effStatErrHi0 << " - " << effStatErrLo0  << std::endl;

    std::cout << "Point " << i << ":" << "x=" << x << " y=" << recEff0 << std::endl;
    grWRecEff0->SetPoint(i,x,eff0);
    grWRecEff0->SetPointError(i,binW/2.0,binW/2.0,effStatErrLo0,effStatErrHi0);
    //if(doWriteToSpreadsheet) writeToSpreadsheet(output,ipt,ieta,icent,recEff0,effStatErrLo0,effStatErrHi0);

    ///RecNoQualityWselNoIso
    int pRec1 = (int)mcWRecSet1->numEntries() ;
    std::cout << "pRec1 = " << pRec1 << std::endl;

    ///Calculate efficiency using Bayesian errors 
    ///with a 68.3% C.I.
    double recEff1, errRec1Hi, errRec1Lo;
    Efficiency(pRec1,tGenInFiducial,0.683,recEff1,errRec1Lo,errRec1Hi);
    if(recEff1<0.0) Warning("calcEfficiency","Efficiency %d. Expect wrong results", recEff1);
    if(recEff1>1.0) Warning("calcEfficiency","Efficiency %d. Expect wrong results", recEff1);

    errRec1Hi = errRec1Hi - recEff1;
    errRec1Lo = recEff1 - errRec1Lo;

    ///Apply trigger efficiency to reconstructed part
    double eff1 = recEff1*trigTemp;

    ///Propagate the error
    double effStatErrHi1 = sqrt( TMath::Power(errRec1Hi/recEff1,2) + TMath::Power(trigStatErr,2) )*eff1;
    double effStatErrLo1 = sqrt( TMath::Power(errRec1Lo/recEff1,2) + TMath::Power(trigStatErr,2) )*eff1;

    std::cout << "Efficiency1: " << recEff1 << " + " << effStatErrHi1 << " - " << effStatErrLo1  << std::endl;

    grWRecEff1->SetPoint(i,x,eff1);
    grWRecEff1->SetPointError(i,binW/2.0,binW/2.0,effStatErrLo1,effStatErrHi1);
    //if(doWriteToSpreadsheet) writeToSpreadsheet(output,ipt,ieta,icent,recEff1,effStatErrLo1,effStatErrHi1);

    ///RecNoQualityWsel
    int pRec2 = (int)mcWRecSet2->numEntries() ;
    std::cout << "pRec2 = " << pRec2 << std::endl;

    ///Calculate efficiency using Bayesian errors 
    ///with a 68.3% C.I.
    double recEff2, errRec2Hi, errRec2Lo;
    Efficiency(pRec2,tGenInFiducial,0.683,recEff2,errRec2Lo,errRec2Hi);
    if(recEff2<0.0) Warning("calcEfficiency","Efficiency %d. Expect wrong results", recEff2);
    if(recEff2>1.0) Warning("calcEfficiency","Efficiency %d. Expect wrong results", recEff2);

    errRec2Hi = errRec2Hi - recEff2;
    errRec2Lo = recEff2 - errRec2Lo;

    ///Apply trigger efficiency to reconstructed part
    double eff2 = recEff2*trigTemp;

    ///Propagate the error
    double effStatErrHi2 = sqrt( TMath::Power(errRec2Hi/recEff2,2) + TMath::Power(trigStatErr,2) )*eff2;
    double effStatErrLo2 = sqrt( TMath::Power(errRec2Lo/recEff2,2) + TMath::Power(trigStatErr,2) )*eff2;

    std::cout << "Efficiency2: " << recEff2 << " + " << effStatErrHi2 << " - " << effStatErrLo2  << std::endl;

    grWRecEff2->SetPoint(i,x,eff2);
    grWRecEff2->SetPointError(i,binW/2.0,binW/2.0,effStatErrLo2,effStatErrHi2);
    //if(doWriteToSpreadsheet) writeToSpreadsheet(output,ipt,ieta,icent,recEff2,effStatErrLo2,effStatErrHi2);

    ///RecHiQualityWselNoIsoNoZveto
    int pRec3 = (int)mcWRecSet3->numEntries() ;
    std::cout << "pRec3 = " << pRec3 << std::endl;

    ///Calculate efficiency using Bayesian errors 
    ///with a 68.3% C.I.
    double recEff3, errRec3Hi, errRec3Lo;
    Efficiency(pRec3,tGenInFiducial,0.683,recEff3,errRec3Lo,errRec3Hi);
    if(recEff3<0.0) Warning("calcEfficiency","Efficiency %d. Expect wrong results", recEff3);
    if(recEff3>1.0) Warning("calcEfficiency","Efficiency %d. Expect wrong results", recEff3);

    errRec3Hi = errRec3Hi - recEff3;
    errRec3Lo = recEff3 - errRec3Lo;

    ///Apply trigger efficiency to reconstructed part
    double eff3 = recEff3*trigTemp;

    ///Propagate the error
    double effStatErrHi3 = sqrt( TMath::Power(errRec3Hi/recEff3,2) + TMath::Power(trigStatErr,2) )*eff3;
    double effStatErrLo3 = sqrt( TMath::Power(errRec3Lo/recEff3,2) + TMath::Power(trigStatErr,2) )*eff3;

    std::cout << "Efficiency3: " << recEff3 << " + " << effStatErrHi3 << " - " << effStatErrLo3  << std::endl;

    grWRecEff3->SetPoint(i,x,eff3);
    grWRecEff3->SetPointError(i,binW/2.0,binW/2.0,effStatErrLo3,effStatErrHi3);
    //if(doWriteToSpreadsheet) writeToSpreadsheet(output,ipt,ieta,icent,recEff3,effStatErrLo3,effStatErrHi3);

    ///RecHiQualityWselNoIso
    int pRec4 = (int)mcWRecSet4->numEntries() ;
    std::cout << "pRec4 = " << pRec4 << std::endl;

    ///Calculate efficiency using Bayesian errors 
    ///with a 68.3% C.I.
    double recEff4, errRec4Hi, errRec4Lo;
    Efficiency(pRec4,tGenInFiducial,0.683,recEff4,errRec4Lo,errRec4Hi);
    if(recEff4<0.0) Warning("calcEfficiency","Efficiency %d. Expect wrong results", recEff4);
    if(recEff4>1.0) Warning("calcEfficiency","Efficiency %d. Expect wrong results", recEff4);

    errRec4Hi = errRec4Hi - recEff4;
    errRec4Lo = recEff4 - errRec4Lo;

    ///Apply trigger efficiency to reconstructed part
    double eff4 = recEff4*trigTemp;

    ///Propagate the error
    double effStatErrHi4 = sqrt( TMath::Power(errRec4Hi/recEff4,2) + TMath::Power(trigStatErr,2) )*eff4;
    double effStatErrLo4 = sqrt( TMath::Power(errRec4Lo/recEff4,2) + TMath::Power(trigStatErr,2) )*eff4;

    std::cout << "Efficiency4: " << recEff4 << " + " << effStatErrHi4 << " - " << effStatErrLo4  << std::endl;

    grWRecEff4->SetPoint(i,x,eff4);
    grWRecEff4->SetPointError(i,binW/2.0,binW/2.0,effStatErrLo4,effStatErrHi4);
    //if(doWriteToSpreadsheet) writeToSpreadsheet(output,ipt,ieta,icent,recEff4,effStatErrLo4,effStatErrHi4);

    ///RecHiQualityWsel
    int pRec5 = (int)mcWRecSet5->numEntries() ;
    std::cout << "pRec5 = " << pRec5 << std::endl;

    ///Calculate efficiency using Bayesian errors 
    ///with a 68.3% C.I.
    double recEff5, errRec5Hi, errRec5Lo;
    Efficiency(pRec5,tGenInFiducial,0.683,recEff5,errRec5Lo,errRec5Hi);
    if(recEff5<0.0) Warning("calcEfficiency","Efficiency %d. Expect wrong results", recEff5);
    if(recEff5>1.0) Warning("calcEfficiency","Efficiency %d. Expect wrong results", recEff5);

    errRec5Hi = errRec5Hi - recEff5;
    errRec5Lo = recEff5 - errRec5Lo;

    ///Apply trigger efficiency to reconstructed part
    std::cout << "Number of reconstructed muons in fiducial cut space = " << pRec5 << std::endl;
    std::cout << "Number of generated muons in fiducial cut space = " << tGenInFiducial << std::endl;
    std::cout << "Consistency Check: " << (double)pRec5/tGenInFiducial << "=?" << recEff5 << std::endl;
    std::cout << "Cw before trigger correction for bin = " << ieta << ":" << icent << " = "  << recEff5 << std::endl;
    std::cout << "Trigger efficiency in bin = " << ieta << ":" << icent << " = "  << trigTemp << std::endl;
    double eff5 = recEff5*trigTemp;
    std::cout << "Cw after trigger correction for bin = " << ieta << ":" << icent << " = "  << eff5 << std::endl;

    ///Propagate the error
    double effStatErrHi5 = sqrt( TMath::Power(errRec5Hi/recEff5,2) + TMath::Power(trigStatErr,2) )*eff5;
    double effStatErrLo5 = sqrt( TMath::Power(errRec5Lo/recEff5,2) + TMath::Power(trigStatErr,2) )*eff5;

    std::cout << "Efficiency5: " << recEff5 << " + " << effStatErrHi5 << " - " << effStatErrLo5  << std::endl;

    grWRecEff5->SetPoint(i,x,eff5);
    grWRecEff5->SetPointError(i,binW/2.0,binW/2.0,effStatErrLo5,effStatErrHi5);

    if(doWriteToSpreadsheet) writeToSpreadsheet(output,ipt,ieta,icent,eff5,effStatErrLo5,effStatErrHi5);

}//calcEfficiency

///////////////////////////////////////////////////////////////////////////////
// selectNuPt
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectNuPt( RooDataSet* dataSet, double nuPtLow)
{
  std::cout << "Cutting on Neutrino pt." << std::endl;
  TString cut = "nuGenPt>";
  cut += nuPtLow;
  std::cout << "Reducing d.s. with " << cut << std::endl;
  dataSet = (RooDataSet*) dataSet->reduce(cut);
  return dataSet;
}


///////////////////////////////////////////////////////////////////////////////
// selectGenPtEta
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectGenPtEta( RooDataSet* dataSet, double ptLow, double ptUpp, double etaLow, double etaUpp, bool doMirrorEta = false )
{
  TString cut = "muonGenPt>";
  cut += ptLow;
  cut += "&&muonGenPt<";
  cut += ptUpp;
  cut += "&&((muEtaGen>";
  cut += etaLow;
  cut += "&&muEtaGen<";
  cut += etaUpp;
  cut += ")";
  if ( doMirrorEta ) {
    cut += "||(muEtaGen>";
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

  dataSet = (RooDataSet*) dataSet->reduce(cut);
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

///////////////////////////////////////////////////////////////////////////////
// fillHIMuonRecSet
///////////////////////////////////////////////////////////////////////////////
RooDataSet* fillHIMuonRecSet(const TString& pathName, const TString& fileName, RooArgSet& muonArgSet, bool doMptSigmaUp = false, bool doMptSigmaDown=false)
{
  RooDataSet* set = new RooDataSet(fileName,fileName,muonArgSet);
  
  float eLossNt[50];
  float scatNt[50];
  float compNt[50];
  float ptNt[50];
  float mtNt[50];
  float etaNt[50];
  float chargeNt[50];
  float centralityNt;
  float nu_ptNt;
  float ptconeNt[50];
  int valNt[50], ZDYNt[50], truthMatchedNt[50], promptNt[50], barcode_muid_truNt[50],mc_mu_gen_barcodeNt[50];
  float  truthdR_muidNt[50];
  int nmu,ngen;
  
  // --- Load the MC tree ---
  TChain* tree = new TChain("tree","tree");
  int nFiles = tree->Add(pathName+fileName);
  
  std::cout << "Filling the RooDataSet for " << fileName << "... Number of files: " << nFiles << std::endl;
 
  // --- Set branch adresses ---

  tree->SetBranchAddress("mc_mu_n", &ngen);
  tree->SetBranchAddress("truthdR_muid",&truthdR_muidNt);
  tree->SetBranchAddress("barcode_muid_tru",&barcode_muid_truNt);
  tree->SetBranchAddress("mc_mu_gen_barcode",&mc_mu_gen_barcodeNt);
  tree->SetBranchAddress("eLoss", &eLossNt);
  //dR<0.3,pTid>3GeV
  tree->SetBranchAddress("ptcone20ID3", &ptconeNt);
  //vary ptcone variable for systematics
//  tree->SetBranchAddress("ptcone30ID3", &ptconeNt);
  tree->SetBranchAddress("scat", &scatNt);
  tree->SetBranchAddress("comp", &compNt);
  tree->SetBranchAddress("pt", &ptNt);
  tree->SetBranchAddress("mt", &mtNt);
  tree->SetBranchAddress("eta", &etaNt);
  tree->SetBranchAddress("charge", &chargeNt);
  tree->SetBranchAddress("prompt", &promptNt);
  tree->SetBranchAddress("val", &valNt); 
  tree->SetBranchAddress("ZDY", &ZDYNt); 
  tree->SetBranchAddress("centrality", &centralityNt);
  if(doMptSigmaDown) {
    std::cout << "WARNING: Using 2GeV track pT threshold relative to the nominal value for mpt." << std::endl;
    tree->SetBranchAddress("nu_pt2000Nominal", &nu_ptNt);
  }
  if(doMptSigmaUp){
    std::cout << "WARNING: Using 4GeV track pT threshold relative to the nominal value for mpt." << std::endl;
    tree->SetBranchAddress("nu_pt4000Nominal", &nu_ptNt);
  }
  else tree->SetBranchAddress("nu_pt", &nu_ptNt);
  tree->SetBranchAddress("truthMatched_muid", &truthMatchedNt);
  tree->SetBranchAddress("mu_muid_n", &nmu);

   // --- Set branch status ---
  tree->SetBranchStatus("*",0) ;
  tree->SetBranchStatus("mc_mu_n", 1);
  tree->SetBranchStatus("truthdR_muid",1);
  tree->SetBranchStatus("barcode_muid_tru",1);
  tree->SetBranchStatus("mc_mu_gen_barcode",1);

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
  tree->SetBranchStatus("charge", 1);
  tree->SetBranchStatus("val", 1); 
  tree->SetBranchStatus("ZDY", 1); 
  tree->SetBranchStatus("centrality", 1);
  tree->SetBranchStatus("nu_pt", 1);
  if(doMptSigmaDown) tree->SetBranchStatus("nu_pt2000Nominal", 1);
  if(doMptSigmaUp) tree->SetBranchStatus("nu_pt4000Nominal", 1);
  
  tree->SetBranchStatus("prompt", 1);

  std::cout << "Number of entries: " << tree->GetEntries() << std::endl; 
  for ( int i = 0; i < tree->GetEntries(); i++ ) {
//     if (i>100000) return set; // temp hack
    tree->LoadTree(i);
    tree->GetEntry(i);

    muonArgSet.setRealValue("missPt",nu_ptNt);
    muonArgSet.setRealValue("centrality",centralityNt);

  for(int igen=0; igen<ngen; ++igen){
    for (int imu = 0; imu<nmu;imu++){

      if( !(truthMatchedNt[imu]==1&&truthdR_muidNt[imu]<0.2&&barcode_muid_truNt[imu]==mc_mu_gen_barcodeNt[igen]))
         continue; 

      if ( promptNt[imu] == 23) muonArgSet.setCatLabel("muonCategory","Z");
      else if ( promptNt[imu] == 24) muonArgSet.setCatLabel("muonCategory","W");

      muonArgSet.setRealValue("muonELoss",eLossNt[imu]);
      muonArgSet.setRealValue("muonScat",scatNt[imu]);
      muonArgSet.setRealValue("muonQuality",valNt[imu]);

      double isolationTemp = -9999.0;
      isolationTemp = ptconeNt[imu]/ptNt[imu];
      muonArgSet.setRealValue("isolation",isolationTemp);
      muonArgSet.setRealValue("muonPt",ptNt[imu]);
      muonArgSet.setRealValue("muonMt",mtNt[imu]);
      muonArgSet.setRealValue("muonEta",etaNt[imu]);
      muonArgSet.setRealValue("ZDY",ZDYNt[imu]);
      muonArgSet.setRealValue("muonCharge",chargeNt[imu]);
      //muonArgSet.setRealValue("muonCategory",promptNt[imu]);
      muonArgSet.setRealValue("motherRec",promptNt[imu]);
      muonArgSet.setRealValue("muonGenRecMatched",truthMatchedNt[imu]);
      if ( chargeNt[imu] > 0 ) muonArgSet.setCatLabel("chargeCategory","muPlus");
      else if ( chargeNt[imu] < 0) muonArgSet.setCatLabel("chargeCategory","muMinus");
      set->add(muonArgSet);    
   }//imu
  }//igen
  }

  return set;
}


RooDataSet* fillHIMuonGenSet(const TString& pathName, const TString& fileName, RooArgSet& muonGenSet)
{
  RooDataSet* set = new RooDataSet(fileName,fileName,muonGenSet);
  
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

      set->add(muonGenSet);    
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



void fitParam(TGraphErrors* const gr, int ich, TString sPar){

	  std::cout << "Fitting to analytic function for " << sPar << ":" << ich << std::endl;
	  TF1 f0a = TF1("f0a","[0]+[1]*x+[2]*x*x",0.0,2.7);	
	  //TF1 f0a = TF1("f0a","TMath::Sqrt([0]*x)",0.0,500.0);	
	  gr->Fit("f0a");
	  f0a.Draw("same");
	  cout << "chi2 = " << f0a.GetChisquare() << "/" << f0a.GetNDF() << ", p = " << f0a.GetProb() << endl;


}

void SaveGraph(TGraph* const gr, TString sY, TString sBin, TString sLabel, bool doSystematic=false){

	TCanvas* c = new TCanvas("c","c",600,600);
	gr->Draw("ape"); 
	gr->GetYaxis()->SetTitle(sY);
	if(doSystematic) gr->GetYaxis()->SetRangeUser(0.0,0.1);
	else gr->GetYaxis()->SetRangeUser(0.0,1.0);
	gr->GetXaxis()->SetTitle("#LT N_{part} #GT");
	myText(0.33,0.89, (Color_t)kBlack, (char*)(sBin));
	c->Update();
	c->Print(sLabel+".pdf");	

}

void fit(TGraphErrors* const gr, TGraphErrors* const grA0, TGraphErrors* const grA1, TGraphErrors* const grA2, int ich, int ieta, double xEta, bool fixA2, bool fixA1){
	  std::cout << "Fitting to analytic function for " << ich << ":" << ieta << std::endl;
	  //a2 coeffs equal for mu+,mu- within 1sigma
	  //double a2 = (2.84778e-7+4.32554e-7)/2.0;
	  //double a2 = (4.56612e-07+7.19152e-07)/2.0;
	  //mu+-
	  double a2 = 5.72135e-07;
	  TF1* f0a; /*= new TF1("f0a","[0]+[1]*x+[2]*x*x",0.0,500.0);*/	
	  //f0a->FixParameter(2,a2);
 
 	//using average of both charges
/*	 double b0 = (-1.22747e-03+-1.15866e-03)/2.0;
	 double b1 = (-5.47588e-05+-2.42778e-05)/2.0;
	 double b2 = (1.24201e-04+1.09511e-04)/2.0;
	 double b3 = (1.15664e-04+8.42216e-05)/2.0;
	 double a1 = b0+b1*TMath::Cos(xEta)+b2*TMath::Cos(2.0*xEta)+b3*TMath::Cos(3.0*xEta);
*/

	 double b0 = -1.208e-3;
	 double b1 = 3.387e-5;
	 double a1 = b0 + b1*xEta;
	 //temp hack 
	 //double a1 = -0.00116983;
	 //mu+ index
	 if(ich==102) { 
		//double four = -6.934e-4-6.0196e-5*TMath::Cos(xEta)+1.406e-5*TMath::Cos(2.*xEta)-9.778e-6*TMath::Cos(3.*xEta)+5.0098e-5*TMath::Cos(xEta*4.);
		//double a1 = -7.668e-4+1.219e-5*xEta*xEta;
		f0a = new TF1("f0a","[0]+[1]*x+[2]*x*x",0.0,500.0);	
		
	 	if(fixA1) f0a->FixParameter(1,a1);
	 	if(fixA2) f0a->FixParameter(2,a2);
		
 	 }
	  //mu- index
	  else if(ich==103) { 
		//double four = -7.7963e-4-1.0117e-4*TMath::Cos(xEta)-4.2672e-5*TMath::Cos(2.*xEta);
		//double a1 = -8.437e-4+3.879e-5*xEta*xEta;
		f0a = new TF1("f0a","[0]+[1]*x+[2]*x*x",0.0,500.0);	
		
	 	if(fixA1) f0a->FixParameter(1,a1);
	 	if(fixA2) f0a->FixParameter(2,a2);
		
		
	  }

	  //mu+- index
	  else {
	  	f0a = new TF1("f0a","[0]+[1]*x+[2]*x*x",0.0,500.0);	

		
	 	if(fixA1) f0a->FixParameter(1,a1);
	 	if(fixA2) f0a->FixParameter(2,a2);
		
	  }
	  //TF1 f0a = TF1("f0a","TMath::Sqrt([0]*x)",0.0,500.0);	
	  gr->Fit("f0a");
	  f0a->Draw("same");
	  cout << "chi2 = " << f0a->GetChisquare() << "/" << f0a->GetNDF() << ", p = " << f0a->GetProb() << endl;

	  //fill graph pts with parameter value to study eta dependence
	  double a0 = f0a->GetParameter(0);
	  double a0err = f0a->GetParError(0);
	  if(!fixA1) a1 = f0a->GetParameter(1);
	  double a1err = f0a->GetParError(1);
	  if(!fixA2) a2 = f0a->GetParameter(2);
	  double a2err = f0a->GetParError(2);

	  grA0->SetPoint(ieta,xEta,a0);
	  grA0->SetPointError(ieta,0.0,a0err);
	  grA1->SetPoint(ieta,xEta,a1);
	  grA1->SetPointError(ieta,0.0,a1err);
	  grA2->SetPoint(ieta,xEta,a2);
	  grA2->SetPointError(ieta,0.0,a2err);
}

void plotFactor(TGraph2DErrors* const gr, TString sCh){

	std::cout << "Plotting factors..." << std::endl;
	TLatex l;
	TCanvas* c = new TCanvas("c","c",600,600);
	l.SetNDC();
	l.DrawLatex(0.39,0.83,sCh);
	gr->Draw("ape");	
	gr->GetXaxis()->SetTitle("#GT N_{part} #LT");	
	gr->GetYaxis()->SetTitle("|#eta|");
	gr->GetZaxis()->SetTitle(sCh);
	gr->GetXaxis()->SetRangeUser(0.0,100.0);	
	gr->GetYaxis()->SetRangeUser(0.0,2.7);
	gr->GetZaxis()->SetRangeUser(0.0,1.0);

	c->Update();
	c->Print(sCh+".pdf"); c->Print(sCh+".root");
	std::cout << "Done." << std::endl;

}
	
void CorrectionFactorsCw(){

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
TGraphErrors* const grToScale = 0;
TGraphErrors* const grToScalePlus =0;
TGraphErrors* const grToScaleMinus =0;
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

//TString fileNameIn = "HISingleMuonWmunuPYTHIADataOverlay.03.24.2013";
//nominal input mc
//TString fileNameIn = "HISingleMuonWmunuPYTHIADataOverlay.04.12.2013";
TString fileNameIn = "HISingleMuonWmunuPYTHIADataOverlay.07.14.2013";
if(doMptSigmaUp||doMptSigmaDown) fileNameIn = "HISingleMuonWmunuPYTHIADataOverlay.07.14.2013";
std::cout << "Input file: " << fileNameIn << std::endl;
///IMPORTANT:IF CHANGING ETA BINNING, YOU MUST CHANGE
///THE .txt FILE HOLDING THE TRIGGER EFFICIENCIES ACCORDINGLY
std::vector <double> etaBins;
//etaBins.push_back(0.0);
if(!doMirrorEta)etaBins.push_back(-2.4);
else etaBins.push_back(0.1);
//else etaBins.push_back(0.0);
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
}
    etaBins.push_back(0.35);
    etaBins.push_back(0.6);
    etaBins.push_back(0.8);
    etaBins.push_back(1.05);
    etaBins.push_back(1.3);
    etaBins.push_back(1.55);
    etaBins.push_back(1.85);
    etaBins.push_back(2.1);




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
//etaBins.push_back(2.5);

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
/*if(doMptSigmaDown){
    ///-1sigma
    mptLow.push_back(mptLowCut-14.1); //0-5%
    mptLow.push_back(mptLowCut-13.5); //5-10%
    mptLow.push_back(mptLowCut-13.2); //10-15%
    mptLow.push_back(mptLowCut-12.4); //15-20%
    mptLow.push_back(mptLowCut-11.2); //20-40%
    mptLow.push_back(mptLowCut-9.2); //40-80%
    for(int impt=0; impt<mptLow.size(); ++impt) std::cout << mptLow[impt] << std::endl;
}
else if (doMptSigmaUp){
    ///+1 sigma
    mptLow.push_back(mptLowCut+14.1); //0-5%
    mptLow.push_back(mptLowCut+13.5); //5-10%
    mptLow.push_back(mptLowCut+13.2); //10-15%
    mptLow.push_back(mptLowCut+12.4); //15-20%
    mptLow.push_back(mptLowCut+11.2); //20-40%
    mptLow.push_back(mptLowCut+9.2); //40-80%
    for(int impt=0; impt<mptLow.size(); ++impt) std::cout << mptLow[impt] << std::endl;
}*/
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
TGraphAsymmErrors* const grCwNpartPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* const grCwNpartMinus = new TGraphAsymmErrors(nCentralityBins);
TGraph* const grCwNpartDiff = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* const grCwNpart = new TGraphAsymmErrors(nCentralityBins);

TGraphAsymmErrors* const grWRecEff0Cent = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* const grWRecEff0CentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* const grWRecEff0CentMinus = new TGraphAsymmErrors(nCentralityBins);

TGraphAsymmErrors* const grWRecEff1Cent = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* const grWRecEff1CentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* const grWRecEff1CentMinus = new TGraphAsymmErrors(nCentralityBins);

TGraphAsymmErrors* const grWRecEff2Cent = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* const grWRecEff2CentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* const grWRecEff2CentMinus = new TGraphAsymmErrors(nCentralityBins);

TGraphAsymmErrors* const grWRecEff3Cent = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* const grWRecEff3CentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* const grWRecEff3CentMinus = new TGraphAsymmErrors(nCentralityBins);

TGraphAsymmErrors* const grWRecEff4Cent = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* const grWRecEff4CentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* const grWRecEff4CentMinus = new TGraphAsymmErrors(nCentralityBins);

TGraphAsymmErrors* const grWRecEff5Cent = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* const grWRecEff5CentPlus = new TGraphAsymmErrors(nCentralityBins);
TGraphAsymmErrors* const grWRecEff5CentMinus = new TGraphAsymmErrors(nCentralityBins);

///eta dep plots
TGraphAsymmErrors* const grCwEtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* const grCwEtaMinus = new TGraphAsymmErrors(nEtaBins);
TGraph* const grCwEtaDiff = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* const grCwEta = new TGraphAsymmErrors(nEtaBins);

TGraphAsymmErrors* const grWRecEff0Eta = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* const grWRecEff0EtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* const grWRecEff0EtaMinus = new TGraphAsymmErrors(nEtaBins);

TGraphAsymmErrors* const grWRecEff1Eta = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* const grWRecEff1EtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* const grWRecEff1EtaMinus = new TGraphAsymmErrors(nEtaBins);

TGraphAsymmErrors* const grWRecEff2Eta = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* const grWRecEff2EtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* const grWRecEff2EtaMinus = new TGraphAsymmErrors(nEtaBins);

TGraphAsymmErrors* const grWRecEff3Eta = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* const grWRecEff3EtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* const grWRecEff3EtaMinus = new TGraphAsymmErrors(nEtaBins);

TGraphAsymmErrors* const grWRecEff4Eta = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* const grWRecEff4EtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* const grWRecEff4EtaMinus = new TGraphAsymmErrors(nEtaBins);

TGraphAsymmErrors* const grWRecEff5Eta = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* const grWRecEff5EtaPlus = new TGraphAsymmErrors(nEtaBins);
TGraphAsymmErrors* const grWRecEff5EtaMinus = new TGraphAsymmErrors(nEtaBins);

///Pt
TGraphAsymmErrors* const grWRecEff0Pt = new TGraphAsymmErrors(nPtBins);
TGraphAsymmErrors* const grWRecEff0PtPlus = new TGraphAsymmErrors(nPtBins);
TGraphAsymmErrors* const grWRecEff0PtMinus = new TGraphAsymmErrors(nPtBins);

TGraphAsymmErrors* const grWRecEff1Pt = new TGraphAsymmErrors(nPtBins);
TGraphAsymmErrors* const grWRecEff1PtPlus = new TGraphAsymmErrors(nPtBins);
TGraphAsymmErrors* const grWRecEff1PtMinus = new TGraphAsymmErrors(nPtBins);

TGraphAsymmErrors* const grWRecEff2Pt = new TGraphAsymmErrors(nPtBins);
TGraphAsymmErrors* const grWRecEff2PtPlus = new TGraphAsymmErrors(nPtBins);
TGraphAsymmErrors* const grWRecEff2PtMinus = new TGraphAsymmErrors(nPtBins);

TGraphAsymmErrors* const grWRecEff3Pt = new TGraphAsymmErrors(nPtBins);
TGraphAsymmErrors* const grWRecEff3PtPlus = new TGraphAsymmErrors(nPtBins);
TGraphAsymmErrors* const grWRecEff3PtMinus = new TGraphAsymmErrors(nPtBins);

TGraphAsymmErrors* const grWRecEff4Pt = new TGraphAsymmErrors(nPtBins);
TGraphAsymmErrors* const grWRecEff4PtPlus = new TGraphAsymmErrors(nPtBins);
TGraphAsymmErrors* const grWRecEff4PtMinus = new TGraphAsymmErrors(nPtBins);

TGraphAsymmErrors* const grWRecEff5Pt = new TGraphAsymmErrors(nPtBins);
TGraphAsymmErrors* const grWRecEff5PtPlus = new TGraphAsymmErrors(nPtBins);
TGraphAsymmErrors* const grWRecEff5PtMinus = new TGraphAsymmErrors(nPtBins);

//residuals for systematics
//TGraph* grResid = new TGraphErrors(nCentralityBins);
TList _cWResiduals;

//graphs for fit parameter eta dependence
TGraphErrors* const grA0Plus = new TGraphErrors(nEtaBins);
TGraphErrors* const grA1Plus = new TGraphErrors(nEtaBins);
TGraphErrors* const grA2Plus = new TGraphErrors(nEtaBins);

TGraphErrors* const grA0Minus = new TGraphErrors(nEtaBins);
TGraphErrors* const grA1Minus = new TGraphErrors(nEtaBins);
TGraphErrors* const grA2Minus = new TGraphErrors(nEtaBins);

//difference in plus/minus for systematics
TGraph* const grA0Diff = new TGraphErrors(nEtaBins);
TGraph* const grA1Diff = new TGraphErrors(nEtaBins);
TGraph* const grA2Diff = new TGraphErrors(nEtaBins);

TGraphErrors* const grA0 = new TGraphErrors(nEtaBins);
TGraphErrors* const grA1 = new TGraphErrors(nEtaBins);
TGraphErrors* const grA2 = new TGraphErrors(nEtaBins);

// --- declare variables at generator level --- //
  RooRealVar  muonGenPt("muonGenPt","p_{T}",0.0,250.0,"GeV");
  RooRealVar  nuGenPt("nuGenPt","p_{T}^{#nu}",0.0,250.0,"GeV");
  RooRealVar  munuGenMt("munuGenMt","m_{T}",0.0,300.0,"GeV");
  RooRealVar  mother("mother","mother",-30.0,30.0);
  RooRealVar  daughter("daughter","daughter",-20.0,20.0);
  RooRealVar  daughterNu("daughterNu","daughterNu",-20.0,20.0);
  RooRealVar  chargeGen("chargeGen","chargeGen",-2.0,2.0);
  RooRealVar  muEtaGen("muEtaGen","muEtaGen",-10.0,10.0);
  RooRealVar  nuEtaGen("nuEtaGen","nuEtaGen",-10.0,10.0);
// --- declare cut variables at reco level --- //
  RooRealVar  muonPt("muonPt","p_{T}",0.0,250.0,"GeV");
  RooRealVar  motherRec("motherRec","motherRec",0.0,250.0);
  RooRealVar  missPt("missPt","p_{T}^{miss}",0.0,250.0,"GeV");
  RooRealVar  muonMt("muonMt","m_{T}",0.0,250.0,"GeV");
  RooRealVar  muonCharge("muonCharge","charge",-2.0,+2.0);
  RooRealVar  isolation("isolation","isolation",0.0,10.0);
  RooRealVar  centrality("centrality","centrality",0.,1.0);
  RooRealVar  muonEta("muonEta","muonEta",-3.0,+3.0);
  RooRealVar  ZDY("ZDY","ZDY",0.0,1.0);
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
  ///current anlysis cuts
  sCutsFid = "muonGenPt>25.0&&abs(muEtaGen)>0.1&&abs(muEtaGen)<2.4&&nuGenPt>25.0&&munuGenMt>40.0&&abs(mother)==24&&abs(daughter)==13&&centrality>0.0&&centrality<0.8";

    ///pt check
//  TString sCutsFid = "muonGenPt>40.0&&abs(muEtaGen)<2.4&&nuGenPt>25.0&&munuGenMt>40.0&&abs(mother)==24&&abs(daughter)==13";
//  TString sCutsGeomFid = "abs(mother)==24&&abs(daughter)==13&&abs(muEtaGen)<2.4";
  RooFormulaVar cutsFid("cutsFid", "cutsFid", sCutsFid,RooArgList(muonGenPt,muEtaGen,nuGenPt,munuGenMt,mother,daughter,centrality));
//  RooFormulaVar cutsGeomFid("cutsGeomFid","cutsGeomFid",sCutsGeomFid,RooArgList(mother,daughter,muEtaGen));

  //reconstruction level cuts
  ///pt check
  TString sCutsRecNoQualityWselNoIsoNoZveto =  "abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0&&muonGenRecMatched==1"; 
  TString sCutsRecNoQualityWselNoIso =
      "abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0&&muonGenRecMatched==1&&ZDY==0"; 
  TString sCutsRecNoQualityWsel =
      "abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0&&muonGenRecMatched==1&&ZDY==0&&isolation<0.1"; 
  TString sCutsRecHiQualityWselNoIsoNoZveto = 
      "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0&&muonGenRecMatched==1";
  TString sCutsRecHiQualityWselNoIso =
      "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0&&ZDY==0&&muonGenRecMatched==1"; 

  ///this is the actual correction factor (Cw)
  TString sCutsRecHiQualityWsel = "";
//  if(doMptSigmaDown||doMptSigmaUp)
//        sCutsRecHiQualityWsel = "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt<9000.0&&isolation<0.1&&muonMt>40.0&&ZDY==0&&muonGenRecMatched==1";
//   else 
       sCutsRecHiQualityWsel =
            "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&isolation<0.1&&muonMt>40.0&&ZDY==0&&muonGenRecMatched==1";
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
  RooArgList muonRecArgList0(muonEta,muonPt,missPt,muonGenRecMatched,muonMt);

  RooArgList muonRecArgList1(muonEta,muonPt,missPt,muonGenRecMatched,muonMt);
  muonRecArgList1.add(ZDY);

  RooArgList muonRecArgList2(muonEta,muonPt,missPt,muonGenRecMatched,muonMt);
  muonRecArgList2.add(ZDY);
  muonRecArgList2.add(isolation);

  RooArgList muonRecArgList3(muonEta,muonPt,missPt,muonGenRecMatched,muonMt);
  muonRecArgList3.add(muonQuality); 
  muonRecArgList3.add(muonELoss);
  muonRecArgList3.add(muonScat);

  RooArgList muonRecArgList4(muonEta,muonPt,missPt,muonGenRecMatched,muonMt);
  muonRecArgList4.add(muonQuality); 
  muonRecArgList4.add(muonELoss);
  muonRecArgList4.add(muonScat);
  muonRecArgList4.add(ZDY);

  RooArgList muonRecArgList5(muonEta,muonPt,missPt,muonGenRecMatched,muonMt);
  muonRecArgList5.add(muonQuality); 
  muonRecArgList5.add(muonELoss);
  muonRecArgList5.add(muonScat);
  muonRecArgList5.add(ZDY);
  muonRecArgList5.add(isolation);

  ///cuts on datasets at RECONSTRUCTION level
  RooFormulaVar cutsRecNoQualityWselNoIsoNoZveto("cutsRecNoQualityWselNoIsoNoZveto", "cutsRecNoQualityWselNoIsoNoZveto",
        sCutsRecNoQualityWselNoIsoNoZveto, muonRecArgList0);
  RooFormulaVar cutsRecNoQualityWselNoIso("cutsRecNoQualityWselNoIso", "cutsRecNoQualityWselNoIso",
        sCutsRecNoQualityWselNoIso, muonRecArgList1);
  RooFormulaVar cutsRecNoQualityWsel("cutsRecNoQualityWsel", "cutsRecNoQualityWsel",
        sCutsRecNoQualityWsel, muonRecArgList2);
  RooFormulaVar cutsRecHiQualityWselNoIsoNoZveto("cutsRecHiQualityWselNoIsoNoZveto", "cutsRecHiQualityWselNoIsoNoZveto",
        sCutsRecHiQualityWselNoIsoNoZveto, muonRecArgList3);
  RooFormulaVar cutsRecHiQualityWselNoIso("cutsRecHiQualityWselNoIso", "cutsRecHiQualityWselNoIso",
        sCutsRecHiQualityWselNoIso, muonRecArgList4);
  RooFormulaVar cutsRecHiQualityWsel("cutsRecHiQualityWsel", "cutsRecHiQualityWsel",
        sCutsRecHiQualityWsel, muonRecArgList5);

/*  RooArgList muonRecLooseArgList(muonQuality,muonELoss,muonScat,muonEta,muonGenRecMatched);
  RooFormulaVar cutsRecLoose("cutsRecLoose", "cutsRecLoose", sCutsRecLoose, muonRecLooseArgList);
*/

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


  RooDataSet* mcWGenSet = fillHIMuonGenSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonGenArgSet); mcWGenSet->Print();
  mcWGenSet = (RooDataSet*)mcWGenSet->reduce(Cut(cutsGen)); 
  std::cout << "Number of Wmunu evts at generator level : " << mcWGenSet->numEntries() << std::endl;

  //RooDataSet* mcWGenSetPlus = (RooDataSet*)mcWGenSet->reduce("chargeGenCategory==chargeGenCategory::muPlus"); 
  RooDataSet* mcWGenSetPlus = (RooDataSet*)mcWGenSet->reduce("chargeGen==+1.0"); 
  ///cut on neutrino type for consistency check
  ///neutrino
  //mcWGenSetPlus = (RooDataSet*)mcWGenSetPlus->reduce("daughterNu==14");
  //RooDataSet* mcWGenSetMinus = (RooDataSet*)mcWGenSet->reduce("chargeGenCategory==chargeGenCategory::muMinus"); 
  RooDataSet* mcWGenSetMinus = (RooDataSet*)mcWGenSet->reduce("chargeGen==-1.0"); 
  ///anti-neutrino
  //mcWGenSetMinus = (RooDataSet*)mcWGenSetMinus->reduce("daughterNu==-14");

  /*RooDataSet* mcWGenSetPlusc = (RooDataSet*)mcWGenSet->reduce("chargeGen==+1.0"); 
  RooDataSet* mcWGenSetMinusc = (RooDataSet*)mcWGenSet->reduce("chargeGen==-1.0"); 
  std::cout << "Charge consistency check: " << mcWGenSetPlus->numEntries() << " =? " << mcWGenSetPlusc->numEntries() << std::endl;
  std::cout << "Charge consistency check: " << mcWGenSetMinus->numEntries() << " =? " << mcWGenSetMinusc->numEntries() << std::endl;
  */
  /*RooDataSet* mcWGeomFidSet = fillHIMuonGenSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonGenArgSet); 
  mcWGeomFidSet = (RooDataSet*)mcWGeomFidSet->reduce(Cut(cutsGeomFid)); 
  std::cout << "Number of Wmunu evts in purely geometric fiducial region at generator level : " << mcWGeomFidSet->numEntries() << std::endl;
*/

  RooDataSet* mcWFidSet = fillHIMuonGenSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonGenArgSet); 
  mcWFidSet = (RooDataSet*)mcWFidSet->reduce(Cut(cutsFid)); 
  std::cout << "Number of Wmunu evts in fiducial region at generator level : " << mcWFidSet->numEntries() << std::endl;

  //construct charged sets

//  RooDataSet* mcWGeomFidSetPlus = (RooDataSet*)mcWGeomFidSet->reduce("chargeGenCategory==chargeGenCategory::muPlus"); 
//  RooDataSet* mcWGeomFidSetMinus = (RooDataSet*)mcWGeomFidSet->reduce("chargeGenCategory==chargeGenCategory::muMinus"); 

  //RooDataSet* mcWFidSetPlus = (RooDataSet*)mcWFidSet->reduce("chargeGenCategory==chargeGenCategory::muPlus"); 
  RooDataSet* mcWFidSetPlus = (RooDataSet*)mcWFidSet->reduce("chargeGen==+1.0"); 
  ///neutrino
//  mcWFidSetPlus = (RooDataSet*)mcWFidSetPlus->reduce("daughterNu==14");
  //RooDataSet* mcWFidSetMinus = (RooDataSet*)mcWFidSet->reduce("chargeGenCategory==chargeGenCategory::muMinus"); 
  RooDataSet* mcWFidSetMinus = (RooDataSet*)mcWFidSet->reduce("chargeGen==-1.0"); 
  ///anti-neutrino
//  mcWFidSetMinus = (RooDataSet*)mcWFidSetMinus->reduce("daughterNu==-14");

  double fidPlus = mcWFidSetPlus->numEntries(); double genPlus = mcWGenSetPlus->numEntries();
  double fidMinus = mcWFidSetMinus->numEntries(); double genMinus = mcWGenSetMinus->numEntries();

  double awPlus,errAwPlusLo,errAwPlusHi;
  Efficiency(fidPlus,genPlus,0.683,awPlus,errAwPlusLo,errAwPlusHi);
  double awMinus,errAwMinusLo,errAwMinusHi;
  Efficiency(fidMinus,genMinus,0.683,awMinus,errAwMinusLo,errAwMinusHi);
  std::cout << " Aw+ = " << mcWFidSetPlus->numEntries() << "/" << mcWGenSetPlus->numEntries() << " = " << fidPlus/genPlus << std::endl;
  std::cout << "Consistency check: " << fidPlus/genPlus << " =? " << awPlus << std::endl;
  std::cout << " Aw- = " << mcWFidSetMinus->numEntries() << "/" << mcWGenSetMinus->numEntries() << " = " << fidMinus/genMinus << std::endl;
  std::cout << "Consistency check: " << fidMinus/genMinus << " =? " << awMinus << std::endl;
  RooBinning b = RooBinning(0.0,300.0);
  b.addUniform(90,0.0,90.0);
  b.addUniform(53,90.0,300.0);
  if(doPlotPt){
	//Create histos
	//entire generator set
	TH1F* hGenPlus = (TH1F*)mcWGenSetPlus->createHistogram("hGenPlus",muonGenPt,Binning(b)); 
	TH1F* hGenMinus = (TH1F*)mcWGenSetMinus->createHistogram("hGenMinus",muonGenPt,Binning(b)); 
	//fiducial generator set
	TH1F* hFidPlus = (TH1F*)mcWFidSetPlus->createHistogram("hFidPlus",muonGenPt,Binning(b)); 
	TH1F* hFidMinus = (TH1F*)mcWFidSetMinus->createHistogram("hFidMinus",muonGenPt,Binning(b)); 

  	plotGeneratorDistros(hGenPlus,hFidPlus,102,"W^{+}#rightarrow#mu^{+}#nu");
  	plotGeneratorDistros(hGenMinus,hFidMinus,103,"W^{-}#rightarrow#mu^{-}#nu");
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
	//Create histos
	//entire generator set
	TH1F* hGenPlus = (TH1F*)mcWGenSetPlus->createHistogram("hGenPlus",muEtaGen,Binning(b2)); 
	TH1F* hGenMinus = (TH1F*)mcWGenSetMinus->createHistogram("hGenMinus",muEtaGen,Binning(b2)); 
	//fiducial generator set
	TH1F* hFidPlus = (TH1F*)mcWFidSetPlus->createHistogram("hFidPlus",muEtaGen,Binning(b2)); 
	TH1F* hFidMinus = (TH1F*)mcWFidSetMinus->createHistogram("hFidMinus",muEtaGen,Binning(b2)); 

  	plotGeneratorDistros(hGenPlus,hFidPlus,102,"W^{+}#rightarrow#mu^{+}#nu");
  	plotGeneratorDistros(hGenMinus,hFidMinus,103,"W^{-}#rightarrow#mu^{-}#nu");
  }

/*  RooDataSet* mcWRecLooseSet = fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonRecArgSet); 
  mcWRecLooseSet = (RooDataSet*)mcWRecLooseSet->reduce(Cut(cutsRecLoose)); 
  mcWRecLooseSet = (RooDataSet*)mcWRecLooseSet->reduce(Cut("motherRec==24"));
  std::cout << "Number of Wmunu evts reconstructed in fiducial region before W cut selection: " << mcWRecLooseSet->numEntries() << std::endl;
*/

/*  RooDataSet* mcWRecSet = fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonRecArgSet); 
  mcWRecSet = (RooDataSet*)mcWRecSet->reduce(Cut(cutsRec)); 
  mcWRecSet = (RooDataSet*)mcWRecSet->reduce(Cut("motherRec==24"));
  std::cout << "Number of Wmunu evts reconstructed in fiducial region : " << mcWRecSet->numEntries() << std::endl;
*/

  RooDataSet* mcSetWRecNoQualityWselNoIsoNoZveto =
        fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonRecArgSet); 
  mcSetWRecNoQualityWselNoIsoNoZveto =
        (RooDataSet*)mcSetWRecNoQualityWselNoIsoNoZveto->reduce(Cut(cutsRecNoQualityWselNoIsoNoZveto)); 
  std::cout << "Number of Wmunu evts reconstructed in fiducial region for reco set 0: " << mcSetWRecNoQualityWselNoIsoNoZveto->numEntries() << std::endl;
  ///IMPORTANT: The two lines below do not give the same answer (they should). Look into it.
  mcSetWRecNoQualityWselNoIsoNoZveto = (RooDataSet*)mcSetWRecNoQualityWselNoIsoNoZveto->reduce(Cut("motherRec==24"));
  //mcSetWRecNoQualityWselNoIsoNoZveto = (RooDataSet*)mcSetWRecNoQualityWselNoIsoNoZveto->reduce(Cut("motherRec==24"));
  std::cout << "Number of Wmunu evts reconstructed in fiducial region after W mother selection for reco set 0: " << 
            mcSetWRecNoQualityWselNoIsoNoZveto->numEntries() << std::endl;
  RooDataSet* mcSetWRecNoQualityWselNoIso =
        fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonRecArgSet); 
  mcSetWRecNoQualityWselNoIso =
        (RooDataSet*)mcSetWRecNoQualityWselNoIso->reduce(Cut(cutsRecNoQualityWselNoIso)); 
  mcSetWRecNoQualityWselNoIso = (RooDataSet*)mcSetWRecNoQualityWselNoIso->reduce(Cut("motherRec==24"));
  std::cout << "Number of Wmunu evts reconstructed in fiducial region for reco set 1: " << mcSetWRecNoQualityWselNoIso->numEntries() << std::endl;

  RooDataSet* mcSetWRecNoQualityWsel =
        fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonRecArgSet); 
  mcSetWRecNoQualityWsel =
        (RooDataSet*)mcSetWRecNoQualityWsel->reduce(Cut(cutsRecNoQualityWsel)); 
  mcSetWRecNoQualityWsel = (RooDataSet*)mcSetWRecNoQualityWsel->reduce(Cut("motherRec==24"));
  std::cout << "Number of Wmunu evts reconstructed in fiducial region for reco set 2: " << mcSetWRecNoQualityWsel->numEntries() << std::endl;

  RooDataSet* mcSetWRecHiQualityWselNoIsoNoZveto =
        fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonRecArgSet); 
  mcSetWRecHiQualityWselNoIsoNoZveto =
        (RooDataSet*)mcSetWRecHiQualityWselNoIsoNoZveto->reduce(Cut(cutsRecHiQualityWselNoIsoNoZveto)); 
  mcSetWRecHiQualityWselNoIsoNoZveto = (RooDataSet*)mcSetWRecHiQualityWselNoIsoNoZveto->reduce(Cut("motherRec==24"));
  std::cout << "Number of Wmunu evts reconstructed in fiducial region for reco set 3: " << mcSetWRecHiQualityWselNoIsoNoZveto->numEntries() << std::endl;

  RooDataSet* mcSetWRecHiQualityWselNoIso =
        fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonRecArgSet); 
  mcSetWRecHiQualityWselNoIso =
        (RooDataSet*)mcSetWRecHiQualityWselNoIso->reduce(Cut(cutsRecHiQualityWselNoIso)); 
  mcSetWRecHiQualityWselNoIso = (RooDataSet*)mcSetWRecHiQualityWselNoIso->reduce(Cut("motherRec==24"));
  std::cout << "Number of Wmunu evts reconstructed in fiducial region for reco set 4: " << mcSetWRecHiQualityWselNoIso->numEntries() << std::endl;

  RooDataSet* mcSetWRecHiQualityWsel =
        fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonRecArgSet,doMptSigmaUp,doMptSigmaDown); 
  mcSetWRecHiQualityWsel =
        (RooDataSet*)mcSetWRecHiQualityWsel->reduce(Cut(cutsRecHiQualityWsel)); 
  std::cout << "Number of Wmunu evts reconstructed in fiducial region for reco set 5 before W mother selection: " << mcSetWRecHiQualityWsel->numEntries() << std::endl;
  mcSetWRecHiQualityWsel = (RooDataSet*)mcSetWRecHiQualityWsel->reduce(Cut("motherRec==24"));
  //mcSetWRecHiQualityWsel = (RooDataSet*)mcSetWRecHiQualityWsel->reduce(Cut("muonCategory==muonCategory::W"));
  std::cout << "Number of Wmunu evts reconstructed in fiducial region for reco set 5 after W mother selection: " << mcSetWRecHiQualityWsel->numEntries() << std::endl;

  // --- Subdivide in bins ---
///generator level subsets
  RooDataSet* mcWGenSubSet[nPtBins][nEtaBins][nCentralityBins];
  RooDataSet* mcWFidSubSet[nPtBins][nEtaBins][nCentralityBins];
//  RooDataSet* mcWGeomFidSubSet[nPtBins][nEtaBins][nCentralityBins];
//  RooDataSet* mcWGeomFidSubSetPlus[nPtBins][nEtaBins][nCentralityBins];
//  RooDataSet* mcWGeomFidSubSetMinus[nPtBins][nEtaBins][nCentralityBins];
//  RooDataSet* mcWRecLooseSubSet[nPtBins][nEtaBins][nCentralityBins];
//  RooDataSet* mcWRecLooseSubSetPlus[nPtBins][nEtaBins][nCentralityBins];
//  RooDataSet* mcWRecLooseSubSetMinus[nPtBins][nEtaBins][nCentralityBins];
//  RooDataSet* mcWRecSubSet[nPtBins][nEtaBins][nCentralityBins];
///reconstructed sub-sets
  RooDataSet* mcSubSetWRecNoQualityWselNoIsoNoZveto[nPtBins][nEtaBins][nCentralityBins];
  RooDataSet* mcSubSetWRecNoQualityWselNoIso[nPtBins][nEtaBins][nCentralityBins];
  RooDataSet* mcSubSetWRecNoQualityWsel[nPtBins][nEtaBins][nCentralityBins];
  RooDataSet* mcSubSetWRecHiQualityWselNoIsoNoZveto[nPtBins][nEtaBins][nCentralityBins];
  RooDataSet* mcSubSetWRecHiQualityWselNoIso[nPtBins][nEtaBins][nCentralityBins];
  RooDataSet* mcSubSetWRecHiQualityWsel[nPtBins][nEtaBins][nCentralityBins];

///cuts the datasets in pt,eta,cent bins
  for ( int i = 0; i < nPtBins; i++ ) {
    for ( int j = 0; j < nEtaBins; j++ ) {
	//add residual graphs to the list
	_cWResiduals.Add( new TGraph(nCentralityBins) );

      ///set first bool false to not mirror eta
      for ( int k = 0; k < nCentralityBins; k++ ){

        std::cout << "pt :" << ptBins[i] << "-" << ptBins[i+1] << " eta :" << etaBins[j] << "-" << etaBins[j+1] << 
            " cent: " << centBins[k] << "-" << centBins[k+1] << std::endl;

        mcWGenSubSet[i][j][k] = selectPtEtaCentrality(mcWGenSet ,ptBins[i],ptBins[i+1], etaBins[j],
            etaBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 
        std::cout << "mcWGenSubSet in " << i << ":" << j << ":" << k << " = " << mcWGenSubSet[i][j][k]->numEntries() << std::endl;

        mcWFidSubSet[i][j][k] = selectPtEtaCentrality(mcWFidSet ,ptBins[i],ptBins[i+1], etaBins[j],
            etaBins[j+1],centBins[k], centBins[k+1],doMirrorEta,true); 
        std::cout << "mcWFidSubSet in " << i << ":" << j << ":" << k << " = " << mcWFidSubSet[i][j][k]->numEntries()<< std::endl;

        mcSubSetWRecNoQualityWselNoIsoNoZveto[i][j][k] = selectPtEtaCentrality(mcSetWRecNoQualityWselNoIsoNoZveto,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1],doMirrorEta,false);
          std::cout << "mcSubSetWRecNoQualityWselNoIsoNoZveto in "<< i << ":" << j << ":" << k << 
	  	" = " << mcSubSetWRecNoQualityWselNoIsoNoZveto[i][j][k]->numEntries()<< std::endl;

          mcSubSetWRecNoQualityWselNoIso[i][j][k] = selectPtEtaCentrality(mcSetWRecNoQualityWselNoIso,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1],doMirrorEta,false);
          std::cout << "mcSubSetWRecNoQualityWselNoIso in "<< i << ":" << j << ":" << k << 
	  	" = " << mcSubSetWRecNoQualityWselNoIso[i][j][k]->numEntries()<< std::endl;

          mcSubSetWRecNoQualityWsel[i][j][k] = selectPtEtaCentrality(mcSetWRecNoQualityWsel,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1],doMirrorEta,false);
          std::cout << "mcSubSetWRecNoQualityWsel in "<< i << ":" << j << ":" << k << " = " << mcSubSetWRecNoQualityWsel[i][j][k]->numEntries()<< std::endl;

          mcSubSetWRecHiQualityWselNoIsoNoZveto[i][j][k] = selectPtEtaCentrality(mcSetWRecHiQualityWselNoIsoNoZveto,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1],doMirrorEta,false);
          std::cout << "mcSubSetWRecHiQualityWselNoIsoNoZveto in "<< i << ":" << j << ":" << k << 
	  	" = " << mcSubSetWRecHiQualityWselNoIsoNoZveto[i][j][k]->numEntries()<< std::endl;

          mcSubSetWRecHiQualityWselNoIso[i][j][k] = selectPtEtaCentrality(mcSetWRecHiQualityWselNoIso,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1],doMirrorEta,false);
          std::cout << "mcSubSetWRecHiQualityWselNoIso in "<< i << ":" << j << ":" << k << 
	  	" = " << mcSubSetWRecHiQualityWselNoIso[i][j][k]->numEntries()<< std::endl;

          mcSubSetWRecHiQualityWsel[i][j][k] = selectPtEtaCentrality(mcSetWRecHiQualityWsel,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1], doMirrorEta,false);

          ///Use different MPT thresholds for different centrality classes for systematic studies
//          if(doMptSigmaDown||doMptSigmaUp) mcSubSetWRecHiQualityWsel[i][j][k] = selectMpt(mcSubSetWRecHiQualityWsel[i][j][k],mptLow[k]); 
          std::cout << "mcSubSetWRecHiQualityWsel in "<< i << ":" << j << ":" << k << " = " << mcSubSetWRecHiQualityWsel[i][j][k]->numEntries()<< std::endl;
          if(mcSubSetWRecHiQualityWsel[i][j][k]->numEntries()>mcWFidSubSet[i][j][k]->numEntries()){
            std::cout << "WARNING: Efficiency > 1.0. " << std::endl; exit(0); 
          }

    }
  }
}

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

            TString sEffName0 = "grWmunuRecNoQualityWselNoIsoNoZvetoEtaDistroCent"; sEffName0+=icent;
            TString sEffNamePlus0 = "grWmunuRecNoQualityWselNoIsoNoZvetoEtaDistroPlusCent"; sEffNamePlus0+=icent;
            TString sEffNameMinus0 = "grWmunuRecNoQualityWselNoIsoNoZvetoEtaDistroMinusCent"; sEffNameMinus0+=icent;

            TString sEffName1 = "grWmunuRecNoQualityWselNoIsoEtaDistroCent"; sEffName1+=icent;
            TString sEffNamePlus1 = "grWmunuRecNoQualityWselNoIsoEtaDistroPlusCent"; sEffNamePlus1+=icent;
            TString sEffNameMinus1 = "grWmunuRecNoQualityWselNoIsoEtaDistroMinusCent"; sEffNameMinus1+=icent;

            TString sEffName2 = "grWmunuRecNoQualityWselEtaDistroCent"; sEffName2+=icent;
            TString sEffNamePlus2 = "grWmunuRecNoQualityWselEtaDistroPlusCent"; sEffNamePlus2+=icent;
            TString sEffNameMinus2 = "grWmunuRecNoQualityWselEtaDistroMinusCent"; sEffNameMinus2+=icent;

            TString sEffName3 = "grWmunuRecHiQualityWselNoIsoNoZvetoEtaDistroCent"; sEffName3+=icent;
            TString sEffNamePlus3 = "grWmunuRecHiQualityWselNoIsoNoZvetoEtaDistroPlusCent"; sEffNamePlus3+=icent;
            TString sEffNameMinus3 = "grWmunuRecHiQualityWselNoIsoNoZvetoEtaDistroMinusCent"; sEffNameMinus3+=icent;

            TString sEffName4 = "grWmunuRecHiQualityWselNoIsoEtaDistroCent"; sEffName4+=icent;
            TString sEffNamePlus4 = "grWmunuRecHiQualityWselNoIsoEtaDistroPlusCent"; sEffNamePlus4+=icent;
            TString sEffNameMinus4 = "grWmunuRecHiQualityWselNoIsoEtaDistroMinusCent"; sEffNameMinus4+=icent;

            TString sEffName5 = "grWmunuRecHiQualityWselEtaDistroCent"; sEffName5+=icent;
            TString sEffNamePlus5 = "grWmunuRecHiQualityWselEtaDistroPlusCent"; sEffNamePlus5+=icent;
            TString sEffNameMinus5 = "grWmunuRecHiQualityWselEtaDistroMinusCent"; sEffNameMinus5+=icent;

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

	          TH1F* hSet = (TH1F*)mcWGenSubSet[ipt][ieta][icent]->createHistogram("hSet",muEtaGen);
              std::cout << "Consistency check: tgen = " << hSet->GetEntries() << std::endl;
              calcEfficiency(mcWGenSubSet[ipt][ieta][icent],mcWFidSubSet[ipt][ieta][icent],
                    mcSubSetWRecNoQualityWselNoIsoNoZveto[ipt][ieta][icent], 
                    mcSubSetWRecNoQualityWselNoIso[ipt][ieta][icent],
                    mcSubSetWRecNoQualityWsel[ipt][ieta][icent],
                    mcSubSetWRecHiQualityWselNoIsoNoZveto[ipt][ieta][icent],
                    mcSubSetWRecHiQualityWselNoIso[ipt][ieta][icent],
                    mcSubSetWRecHiQualityWsel[ipt][ieta][icent],
                    ieta, etaMed, binW, index,ipt,ieta,icent,
                    grWRecEff0Eta,grWRecEff1Eta,grWRecEff2Eta,grWRecEff3Eta,grWRecEff4Eta,grWRecEff5Eta,spreadSheet);

              if(doCharge){

               RooDataSet* mcWGenSetPlus  = (RooDataSet*) mcWGenSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muPlus"));
               RooDataSet* mcWGenSetMinus  = (RooDataSet*) mcWGenSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muMinus"));

               RooDataSet* mcWFidSetPlus  = (RooDataSet*) mcWFidSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muPlus"));
               RooDataSet* mcWFidSetMinus  = (RooDataSet*) mcWFidSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muMinus"));


                RooDataSet* mcSetWRecNoQualityWselNoIsoNoZvetoPlus =
                        (RooDataSet*)mcSubSetWRecNoQualityWselNoIsoNoZveto[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecNoQualityWselNoIsoNoZvetoMinus =
                        (RooDataSet*)mcSubSetWRecNoQualityWselNoIsoNoZveto[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

                RooDataSet* mcSetWRecNoQualityWselNoIsoPlus =
                        (RooDataSet*)mcSubSetWRecNoQualityWselNoIso[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecNoQualityWselNoIsoMinus =
                        (RooDataSet*)mcSubSetWRecNoQualityWselNoIso[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

                RooDataSet* mcSetWRecNoQualityWselPlus =
                        (RooDataSet*)mcSubSetWRecNoQualityWsel[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecNoQualityWselMinus =
                        (RooDataSet*)mcSubSetWRecNoQualityWsel[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

                RooDataSet* mcSetWRecHiQualityWselNoIsoNoZvetoPlus =
                        (RooDataSet*)mcSubSetWRecHiQualityWselNoIsoNoZveto[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecHiQualityWselNoIsoNoZvetoMinus =
                        (RooDataSet*)mcSubSetWRecHiQualityWselNoIsoNoZveto[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

                RooDataSet* mcSetWRecHiQualityWselNoIsoPlus =
                        (RooDataSet*)mcSubSetWRecHiQualityWselNoIso[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecHiQualityWselNoIsoMinus =
                        (RooDataSet*)mcSubSetWRecHiQualityWselNoIso[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

                std::cout << "Number of muons inclusive for charge = " <<
                       mcSubSetWRecHiQualityWsel[ipt][ieta][icent]->numEntries() << std::endl;
                RooDataSet* mcSetWRecHiQualityWselPlus =
                        (RooDataSet*)mcSubSetWRecHiQualityWsel[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                std::cout << "Number of muons in mu+ d.s. = " << mcSetWRecHiQualityWselPlus->numEntries() << std::endl;
                RooDataSet* mcSetWRecHiQualityWselMinus =
                        (RooDataSet*)mcSubSetWRecHiQualityWsel[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
                std::cout << "Number of muons in mu- d.s. = " << mcSetWRecHiQualityWselMinus->numEntries() << std::endl;

	    	    std::cout << "mu^{+}" << std::endl;
                calcEfficiency(mcWGenSetPlus,mcWFidSetPlus,
                    mcSetWRecNoQualityWselNoIsoNoZvetoPlus,mcSetWRecNoQualityWselNoIsoPlus,mcSetWRecNoQualityWselPlus,
                    mcSetWRecHiQualityWselNoIsoNoZvetoPlus, mcSetWRecHiQualityWselNoIsoPlus, mcSetWRecHiQualityWselPlus,
                    ieta, etaMed, binW, index,ipt, ieta,icent,
                    grWRecEff0EtaPlus,grWRecEff1EtaPlus,grWRecEff2EtaPlus,grWRecEff3EtaPlus,grWRecEff4EtaPlus,grWRecEff5EtaPlus,spreadSheetPlus);

	            std::cout << "mu^{-}" << std::endl;
                calcEfficiency(mcWGenSetMinus,mcWFidSetMinus,
                    mcSetWRecNoQualityWselNoIsoNoZvetoMinus,mcSetWRecNoQualityWselNoIsoMinus,mcSetWRecNoQualityWselMinus,
                    mcSetWRecHiQualityWselNoIsoNoZvetoMinus, mcSetWRecHiQualityWselNoIsoMinus, mcSetWRecHiQualityWselMinus,
                    ieta, etaMed, binW, index,ipt,ieta,icent,
                    grWRecEff0EtaMinus,grWRecEff1EtaMinus,grWRecEff2EtaMinus,grWRecEff3EtaMinus,grWRecEff4EtaMinus,grWRecEff5EtaMinus,spreadSheetMinus);
              } //doCharge

            
            }//ieta

              Write(outFile,grWRecEff0Eta,sEffName0);
              Write(outFile,grWRecEff1Eta,sEffName1);
              Write(outFile,grWRecEff2Eta,sEffName2);
              Write(outFile,grWRecEff3Eta,sEffName3);
              Write(outFile,grWRecEff4Eta,sEffName4);
              Write(outFile,grWRecEff5Eta,sEffName5);
              if(doCharge){
                Write(outFile,grWRecEff0EtaPlus,sEffNamePlus0);
                Write(outFile,grWRecEff1EtaPlus,sEffNamePlus1);
                Write(outFile,grWRecEff2EtaPlus,sEffNamePlus2);
                Write(outFile,grWRecEff3EtaPlus,sEffNamePlus3);
                Write(outFile,grWRecEff4EtaPlus,sEffNamePlus4);
                Write(outFile,grWRecEff5EtaPlus,sEffNamePlus5);

                Write(outFile,grWRecEff0EtaMinus,sEffNameMinus0);
                Write(outFile,grWRecEff1EtaMinus,sEffNameMinus1);
                Write(outFile,grWRecEff2EtaMinus,sEffNameMinus2);
                Write(outFile,grWRecEff3EtaMinus,sEffNameMinus3);
                Write(outFile,grWRecEff4EtaMinus,sEffNameMinus4);
                Write(outFile,grWRecEff5EtaMinus,sEffNameMinus5);

              }
          }//icent
      }//ipt
}

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

/*            TString sEffName = "grWmunuRecEffNpartDistroEta"; sEffName+=ieta;
            TString sEffNamePlus = "grWmunuRecEffNpartDistroPlusEta"; sEffNamePlus+=ieta;
            TString sEffNameMinus = "grWmunuRecEffNpartDistroMinusEta"; sEffNameMinus+=ieta;
*/
            TString sEffName0 = "grWmunuRecNoQualityWselNoIsoNoZvetoNpartDistroEta"; sEffName0+=ieta;
            TString sEffNamePlus0 = "grWmunuRecNoQualityWselNoIsoNoZvetoNpartDistroPlusEta"; sEffNamePlus0+=ieta;
            TString sEffNameMinus0 = "grWmunuRecNoQualityWselNoIsoNoZvetoNpartDistroMinusEta"; sEffNameMinus0+=ieta;

            TString sEffName1 = "grWmunuRecNoQualityWselNoIsoNpartDistroEta"; sEffName1+=ieta;
            TString sEffNamePlus1 = "grWmunuRecNoQualityWselNoIsoNpartDistroPlusEta"; sEffNamePlus1+=ieta;
            TString sEffNameMinus1 = "grWmunuRecNoQualityWselNoIsoNpartDistroMinusEta"; sEffNameMinus1+=ieta;

            TString sEffName2 = "grWmunuRecNoQualityWselNpartDistroEta"; sEffName2+=ieta;
            TString sEffNamePlus2 = "grWmunuRecNoQualityWselNpartDistroPlusEta"; sEffNamePlus2+=ieta;
            TString sEffNameMinus2 = "grWmunuRecNoQualityWselNpartDistroMinusEta"; sEffNameMinus2+=ieta;

            TString sEffName3 = "grWmunuRecHiQualityWselNoIsoNoZvetoNpartDistroEta"; sEffName3+=ieta;
            TString sEffNamePlus3 = "grWmunuRecHiQualityWselNoIsoNoZvetoNpartDistroPlusEta"; sEffNamePlus3+=ieta;
            TString sEffNameMinus3 = "grWmunuRecHiQualityWselNoIsoNoZvetoNpartDistroMinusEta"; sEffNameMinus3+=ieta;

            TString sEffName4 = "grWmunuRecHiQualityWselNoIsoNpartDistroEta"; sEffName4+=ieta;
            TString sEffNamePlus4 = "grWmunuRecHiQualityWselNoIsoNpartDistroPlusEta"; sEffNamePlus4+=ieta;
            TString sEffNameMinus4 = "grWmunuRecHiQualityWselNoIsoNpartDistroMinusEta"; sEffNameMinus4+=ieta;

            TString sEffName5 = "grWmunuRecHiQualityWselNpartDistroEta"; sEffName5+=ieta;
            TString sEffNamePlus5 = "grWmunuRecHiQualityWselNpartDistroPlusEta"; sEffNamePlus5+=ieta;
            TString sEffNameMinus5 = "grWmunuRecHiQualityWselNpartDistroMinusEta"; sEffNameMinus5+=ieta;


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
                

	    std::cout << "mu^{#pm}" << std::endl;

        calcEfficiency(mcWGenSubSet[ipt][ieta][icent],mcWFidSubSet[ipt][ieta][icent],
                    mcSubSetWRecNoQualityWselNoIsoNoZveto[ipt][ieta][icent], 
                    mcSubSetWRecNoQualityWselNoIso[ipt][ieta][icent],
                    mcSubSetWRecNoQualityWsel[ipt][ieta][icent],
                    mcSubSetWRecHiQualityWselNoIsoNoZveto[ipt][ieta][icent],
                    mcSubSetWRecHiQualityWselNoIso[ipt][ieta][icent],
                    mcSubSetWRecHiQualityWsel[ipt][ieta][icent],
                    icent, npart, binW, index,ipt,ieta,icent,
                    grWRecEff0Cent,grWRecEff1Cent,grWRecEff2Cent,grWRecEff3Cent,grWRecEff4Cent,grWRecEff5Cent,spreadSheet,true);



	    if(doCharge) {

		RooDataSet* mcWGenSetPlus  = (RooDataSet*) mcWGenSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muPlus"));
		RooDataSet* mcWGenSetMinus  = (RooDataSet*) mcWGenSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muMinus"));

		RooDataSet* mcWFidSetPlus  = (RooDataSet*) mcWFidSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muPlus"));
		RooDataSet* mcWFidSetMinus  = (RooDataSet*) mcWFidSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muMinus"));

                RooDataSet* mcSetWRecNoQualityWselNoIsoNoZvetoPlus =
                        (RooDataSet*)mcSubSetWRecNoQualityWselNoIsoNoZveto[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecNoQualityWselNoIsoNoZvetoMinus =
                        (RooDataSet*)mcSubSetWRecNoQualityWselNoIsoNoZveto[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

                RooDataSet* mcSetWRecNoQualityWselNoIsoPlus =
                        (RooDataSet*)mcSubSetWRecNoQualityWselNoIso[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecNoQualityWselNoIsoMinus =
                        (RooDataSet*)mcSubSetWRecNoQualityWselNoIso[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

                RooDataSet* mcSetWRecNoQualityWselPlus =
                        (RooDataSet*)mcSubSetWRecNoQualityWsel[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecNoQualityWselMinus =
                        (RooDataSet*)mcSubSetWRecNoQualityWsel[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

                RooDataSet* mcSetWRecHiQualityWselNoIsoNoZvetoPlus =
                        (RooDataSet*)mcSubSetWRecHiQualityWselNoIsoNoZveto[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecHiQualityWselNoIsoNoZvetoMinus =
                        (RooDataSet*)mcSubSetWRecHiQualityWselNoIsoNoZveto[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

                RooDataSet* mcSetWRecHiQualityWselNoIsoPlus =
                        (RooDataSet*)mcSubSetWRecHiQualityWselNoIso[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecHiQualityWselNoIsoMinus =
                        (RooDataSet*)mcSubSetWRecHiQualityWselNoIso[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

                std::cout << "Number of muons inclusive for charge = " <<
                       mcSubSetWRecHiQualityWsel[ipt][ieta][icent]->numEntries() << std::endl;
                RooDataSet* mcSetWRecHiQualityWselPlus =
                        (RooDataSet*)mcSubSetWRecHiQualityWsel[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecHiQualityWselPlusc =
                        (RooDataSet*)mcSubSetWRecHiQualityWsel[ipt][ieta][icent]->reduce(Cut("muonCharge==+1.0"));
                std::cout << "Number of muons in mu+ d.s. = " << mcSetWRecHiQualityWselPlus->numEntries() << std::endl;
                std::cout << "Consistency check. " << mcSetWRecHiQualityWselPlus->numEntries() << " =? " << mcSetWRecHiQualityWselPlusc->numEntries() << std::endl;
                RooDataSet* mcSetWRecHiQualityWselMinus =
                        (RooDataSet*)mcSubSetWRecHiQualityWsel[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
                RooDataSet* mcSetWRecHiQualityWselMinusc =
                        (RooDataSet*)mcSubSetWRecHiQualityWsel[ipt][ieta][icent]->reduce(Cut("muonCharge==-1.0"));
                std::cout << "Number of muons in mu- d.s. = " << mcSetWRecHiQualityWselMinus->numEntries() << std::endl;
                std::cout << "Consistency check. " << mcSetWRecHiQualityWselMinus->numEntries() << " =? " << mcSetWRecHiQualityWselMinusc->numEntries() << std::endl;


			  if(doPlotPt){

				///Create histos
				///fiducial generator set
				TH1F* hFidSubSetPlus = (TH1F*)mcWFidSetPlus->createHistogram("hFidSubSetPlus",muonGenPt,Binning(b)); 
				TH1F* hFidSubSetMinus = (TH1F*)mcWFidSetMinus->createHistogram("hFidSubSetMinus",muonGenPt,Binning(b)); 
				//reco level set in fiducial region
				TH1F* hRecSubSetPlus = (TH1F*)mcWRecSetPlus->createHistogram("hRecSubSetPlus",muonPt,Binning(b)); 
				TH1F* hRecSubSetMinus = (TH1F*)mcWRecSetMinus->createHistogram("hRecSubSetMinus",muonPt,Binning(b)); 

				plotFiducialDistros(hFidSubSetPlus,hRecSubSetPlus,"W^{+}#rightarrow#mu^{+}#nu");
				plotFiducialDistros(hFidSubSetMinus,hRecSubSetMinus,"W^{-}#rightarrow#mu^{-}#nu");
			  }//doPlotPt

		std::cout << "mu^{+}" << std::endl;
                calcEfficiency(mcWGenSetPlus,mcWFidSetPlus,
                        mcSetWRecNoQualityWselNoIsoNoZvetoPlus, 
                        mcSetWRecNoQualityWselNoIsoPlus,
                        mcSetWRecNoQualityWselPlus,
                        mcSetWRecHiQualityWselNoIsoNoZvetoPlus,
                        mcSetWRecHiQualityWselNoIsoPlus,
                        mcSetWRecHiQualityWselPlus,
                  	    icent, npart, binW, index,ipt,ieta,icent,
                        grWRecEff0CentPlus,grWRecEff1CentPlus,grWRecEff2CentPlus,
			            grWRecEff3CentPlus,grWRecEff4CentPlus,grWRecEff5CentPlus,spreadSheetPlus,true);



		std::cout << "mu^{-}" << std::endl;
                calcEfficiency(mcWGenSetMinus,mcWFidSetMinus,
                        mcSetWRecNoQualityWselNoIsoNoZvetoMinus, 
                        mcSetWRecNoQualityWselNoIsoMinus,
                        mcSetWRecNoQualityWselMinus,
                        mcSetWRecHiQualityWselNoIsoNoZvetoMinus,
                        mcSetWRecHiQualityWselNoIsoMinus,
                        mcSetWRecHiQualityWselMinus,
                  	    icent, npart, binW, index,ipt,ieta,icent,
                        grWRecEff0CentMinus,grWRecEff1CentMinus,grWRecEff2CentMinus,
			            grWRecEff3CentMinus,grWRecEff4CentMinus,grWRecEff5CentMinus,spreadSheetMinus,true);

			std::cout << "Plotting difference of Cw+ and Cw- for systematics." << std::endl;

			if(doCentrality&&!doEta) plotChargeDiffCwNpart(grCwNpartDiff,grCwNpartPlus,grCwNpartMinus,ieta, icent,npart);
			if(doEta&&!doCentrality) plotChargeDiffCwEta(grCwEtaDiff,grCwEtaPlus,grCwEtaMinus,ieta,icent,xEta);

			} //doCharge

		  } //icent


          Write(outFile,grWRecEff0Cent,sEffName0);
          Write(outFile,grWRecEff1Cent,sEffName1);
          Write(outFile,grWRecEff2Cent,sEffName2);
          Write(outFile,grWRecEff3Cent,sEffName3);
          Write(outFile,grWRecEff4Cent,sEffName4);
          Write(outFile,grWRecEff5Cent,sEffName5);

          if(doCharge){
            Write(outFile,grWRecEff0CentPlus,sEffNamePlus0);
            Write(outFile,grWRecEff1CentPlus,sEffNamePlus1);
            Write(outFile,grWRecEff2CentPlus,sEffNamePlus2);
            Write(outFile,grWRecEff3CentPlus,sEffNamePlus3);
            Write(outFile,grWRecEff4CentPlus,sEffNamePlus4);
            Write(outFile,grWRecEff5CentPlus,sEffNamePlus5);

            Write(outFile,grWRecEff0CentMinus,sEffNameMinus0);
            Write(outFile,grWRecEff1CentMinus,sEffNameMinus1);
            Write(outFile,grWRecEff2CentMinus,sEffNameMinus2);
            Write(outFile,grWRecEff3CentMinus,sEffNameMinus3);
            Write(outFile,grWRecEff4CentMinus,sEffNameMinus4);
            Write(outFile,grWRecEff5CentMinus,sEffNameMinus5);

          }

		TString sNamePlus = "grCent";sNamePlus+=102;  sNamePlus+="eta";sNamePlus+=ieta;
		TString sNameMinus = "grCent";sNameMinus+=103; sNameMinus+="eta";sNameMinus+=ieta;
		TString sNameDiff = "grCentDiff"; sNameDiff+="eta";sNameDiff+=ieta;
		TString sName = "grCent";sName+=104; sName+="eta";sName+=ieta;

		TString sEtaTemp = ""; sEtaTemp+=etabinLo; sEtaTemp+="<";sEtaTemp+="|#eta|"; sEtaTemp+="<"; sEtaTemp+=etabinUp; 

		if(doFit){
			fit(grCwNpartPlus,grA0Plus,grA1Plus,grA2Plus,102,ieta,xEta,false,false);
			fit(grCwNpartMinus,grA0Minus,grA1Minus,grA2Minus,103,ieta,xEta,false,false);

			std::cout << "Plotting charge difference of a0,a1,and a2 for systematics." << std::endl;
			plotChargeDiffAx(grA0Diff,grA0Plus,grA0Minus,ieta,xEta);
			plotChargeDiffAx(grA1Diff,grA1Plus,grA1Minus,ieta,xEta);
			plotChargeDiffAx(grA2Diff,grA2Plus,grA2Minus,ieta,xEta);

			fit(grCwNpart,grA0,grA1,grA2,104,ieta,xEta,false,false);
			//plot the residuals for systematics
			plotCwResidual(grCwNpart, (TGraph*)_cWResiduals.At(ieta));	

			/*SaveGraph(grCwNpartPlus,"C_{W^{+}}",sEtaTemp,sNamePlus);
			SaveGraph(grCwNpartMinus,"C_{W^{-}}",sEtaTemp,sNameMinus);
			SaveGraph(grCwNpartDiff,"|C_{W^{+}}-C_{W^{-}}|",sEtaTemp,sNameDiff,true);
			SaveGraph(grCwNpart,"C_{W}",sEtaTemp,sName);
			SaveGraph((TGraph*)_cWResiduals.At(ieta),"C_{W}",sEtaTemp,sName+"_Residuals");
            */

			Write(outFile,grCwNpartPlus,sNamePlus);
			Write(outFile,grCwNpartMinus,sNameMinus);
			Write(outFile,grCwNpartDiff,sNameDiff);
			Write(outFile,grCwNpart,sName);
			Write(outFile,(TGraph*)_cWResiduals.At(ieta),sName+"_Residuals");
		}//doFit
	   } //ieta

	} //ipt
   }

  //////////////////////////////////////////////////////
  ///Loop for plotting efficiencies as fcn of pt in
  ///each centrality bin and eta slice
  //////////////////////////////////////////////////////
  if (doPt){
    std::cout << "Plotting pT dependence." << std::endl;
    for(int ieta=0; ieta<nEtaBins; ieta++){
	  for(int icent=0; icent<nCentralityBins; icent++){

            TString sEffName0 = "grWmunuRecNoQualityWselNoIsoNoZvetoPtDistroEta"; sEffName0+=ieta; sEffName0+="Cent";
                sEffName0+=icent;

            TString sEffName0 = "grWmunuRecNoQualityWselNoIsoNoZvetoPtDistroEta"; sEffName0+=ieta; sEffName0+="Cent";sEffName0+=icent;
            TString sEffNamePlus0 = "grWmunuRecNoQualityWselNoIsoNoZvetoPtDistroPlusEta"; sEffNamePlus0+=ieta; sEffNamePlus0+="Cent";sEffNamePlus0+=icent;
            TString sEffNameMinus0 = "grWmunuRecNoQualityWselNoIsoNoZvetoPtDistroMinusEta"; sEffNameMinus0+=ieta; sEffNameMinus0+="Cent";sEffNameMinus0+=icent;

            TString sEffName1 = "grWmunuRecNoQualityWselNoIsoPtDistroEta"; sEffName1+=ieta; sEffName1+="Cent";sEffName1+=icent;
            TString sEffNamePlus1 = "grWmunuRecNoQualityWselNoIsoPtDistroPlusEta"; sEffNamePlus1+=ieta; sEffNamePlus1+="Cent";sEffNamePlus1+=icent;
            TString sEffNameMinus1 = "grWmunuRecNoQualityWselNoIsoPtDistroMinusEta"; sEffNameMinus1+=ieta; sEffNameMinus1+="Cent";sEffNameMinus1+=icent;

            TString sEffName2 = "grWmunuRecNoQualityWselPtDistroEta"; sEffName2+=ieta; sEffName2+="Cent";sEffName2+=icent;
            TString sEffNamePlus2 = "grWmunuRecNoQualityWselPtDistroPlusEta"; sEffNamePlus2+=ieta; sEffNamePlus2+="Cent";sEffNamePlus2+=icent;
            TString sEffNameMinus2 = "grWmunuRecNoQualityWselPtDistroMinusEta"; sEffNameMinus2+=ieta; sEffNameMinus2+="Cent";sEffNameMinus2+=icent;

            TString sEffName3 = "grWmunuRecHiQualityWselNoIsoNoZvetoPtDistroEta"; sEffName3+=ieta; sEffName3+="Cent";sEffName3+=icent;
            TString sEffNamePlus3 = "grWmunuRecHiQualityWselNoIsoNoZvetoPtDistroPlusEta"; sEffNamePlus3+=ieta; sEffNamePlus3+="Cent";sEffNamePlus3+=icent;
            TString sEffNameMinus3 = "grWmunuRecHiQualityWselNoIsoNoZvetoPtDistroMinusEta"; sEffNameMinus3+=ieta; sEffNameMinus3+="Cent";sEffNameMinus3+=icent;

            TString sEffName4 = "grWmunuRecHiQualityWselNoIsoPtDistroEta"; sEffName4+=ieta; sEffName4+="Cent";sEffName4+=icent;
            TString sEffNamePlus4 = "grWmunuRecHiQualityWselNoIsoPtDistroPlusEta"; sEffNamePlus4+=ieta; sEffNamePlus4+="Cent";sEffNamePlus4+=icent;
            TString sEffNameMinus4 = "grWmunuRecHiQualityWselNoIsoPtDistroMinusEta"; sEffNameMinus4+=ieta; sEffNameMinus4+="Cent";sEffNameMinus4+=icent;

            TString sEffName5 = "grWmunuRecHiQualityWselPtDistroEta"; sEffName5+=ieta; sEffName5+="Cent";sEffName5+=icent;
            TString sEffNamePlus5 = "grWmunuRecHiQualityWselPtDistroPlusEta"; sEffNamePlus5+=ieta; sEffNamePlus5+="Cent";sEffNamePlus5+=icent;
            TString sEffNameMinus5 = "grWmunuRecHiQualityWselPtDistroMinusEta"; sEffNameMinus5+=ieta; sEffNameMinus5+="Cent";sEffNameMinus5+=icent;

            for(int ipt=0; ipt<nPtBins; ipt++){

              double etaLow = etaBins.at(ieta);
              double etaUpp = etaBins.at(ieta+1); 

              //int index = icent*nEtaBins+ieta; 
	          double etaMed = etaBins[ieta]+(etaBins[ieta+1]-etaBins[ieta])/2.0;
              double ptMed = ptBins[ipt]+(ptBins[ipt+1]-ptBins[ipt])/2.0;
              double binW = ptBins[ipt+1]-ptBins[ipt];
	          int index = -1; 
              
              ///Assume no pT dependence of trigger efficiency
              ///past turn on region
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
              calcEfficiency(mcWGenSubSet[ipt][ieta][icent],mcWFidSubSet[ipt][ieta][icent],
                    mcSubSetWRecNoQualityWselNoIsoNoZveto[ipt][ieta][icent], 
                    mcSubSetWRecNoQualityWselNoIso[ipt][ieta][icent],
                    mcSubSetWRecNoQualityWsel[ipt][ieta][icent],
                    mcSubSetWRecHiQualityWselNoIsoNoZveto[ipt][ieta][icent],
                    mcSubSetWRecHiQualityWselNoIso[ipt][ieta][icent],
                    mcSubSetWRecHiQualityWsel[ipt][ieta][icent],
                    ipt, ptMed, binW, index,ipt,ieta,icent,
                    grWRecEff0Pt,grWRecEff1Pt,grWRecEff2Pt,grWRecEff3Pt,grWRecEff4Pt,grWRecEff5Pt,spreadSheet,false,true);

              if(doCharge){

               RooDataSet* mcWGenSetPlus  = (RooDataSet*) mcWGenSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muPlus"));
               RooDataSet* mcWGenSetMinus  = (RooDataSet*) mcWGenSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muMinus"));

               RooDataSet* mcWFidSetPlus  = (RooDataSet*) mcWFidSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muPlus"));
               RooDataSet* mcWFidSetMinus  = (RooDataSet*) mcWFidSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muMinus"));


                RooDataSet* mcSetWRecNoQualityWselNoIsoNoZvetoPlus =
                        (RooDataSet*)mcSubSetWRecNoQualityWselNoIsoNoZveto[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecNoQualityWselNoIsoNoZvetoMinus =
                        (RooDataSet*)mcSubSetWRecNoQualityWselNoIsoNoZveto[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

                RooDataSet* mcSetWRecNoQualityWselNoIsoPlus =
                        (RooDataSet*)mcSubSetWRecNoQualityWselNoIso[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecNoQualityWselNoIsoMinus =
                        (RooDataSet*)mcSubSetWRecNoQualityWselNoIso[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

                RooDataSet* mcSetWRecNoQualityWselPlus =
                        (RooDataSet*)mcSubSetWRecNoQualityWsel[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecNoQualityWselMinus =
                        (RooDataSet*)mcSubSetWRecNoQualityWsel[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

                RooDataSet* mcSetWRecHiQualityWselNoIsoNoZvetoPlus =
                        (RooDataSet*)mcSubSetWRecHiQualityWselNoIsoNoZveto[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecHiQualityWselNoIsoNoZvetoMinus =
                        (RooDataSet*)mcSubSetWRecHiQualityWselNoIsoNoZveto[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

                RooDataSet* mcSetWRecHiQualityWselNoIsoPlus =
                        (RooDataSet*)mcSubSetWRecHiQualityWselNoIso[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecHiQualityWselNoIsoMinus =
                        (RooDataSet*)mcSubSetWRecHiQualityWselNoIso[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

                RooDataSet* mcSetWRecHiQualityWselPlus =
                        (RooDataSet*)mcSubSetWRecHiQualityWsel[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecHiQualityWselMinus =
                        (RooDataSet*)mcSubSetWRecHiQualityWsel[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

	    	    std::cout << "mu^{+}" << std::endl;
                calcEfficiency(mcWGenSetPlus,mcWFidSetPlus,
                    mcSetWRecNoQualityWselNoIsoNoZvetoPlus,mcSetWRecNoQualityWselNoIsoPlus,mcSetWRecNoQualityWselPlus,
                    mcSetWRecHiQualityWselNoIsoNoZvetoPlus, mcSetWRecHiQualityWselNoIsoPlus, mcSetWRecHiQualityWselPlus,
                    ipt, ptMed, binW, index,ipt, ieta,icent,
                    grWRecEff0PtPlus,grWRecEff1PtPlus,grWRecEff2PtPlus,grWRecEff3PtPlus,grWRecEff4PtPlus,grWRecEff5PtPlus,spreadSheetPlus,false,true);

	            std::cout << "mu^{-}" << std::endl;
                calcEfficiency(mcWGenSetMinus,mcWFidSetMinus,
                    mcSetWRecNoQualityWselNoIsoNoZvetoMinus,mcSetWRecNoQualityWselNoIsoMinus,mcSetWRecNoQualityWselMinus,
                    mcSetWRecHiQualityWselNoIsoNoZvetoMinus, mcSetWRecHiQualityWselNoIsoMinus, mcSetWRecHiQualityWselMinus,
                    ipt, ptMed, binW, index,ipt,ieta,icent,
                    grWRecEff0PtMinus,grWRecEff1PtMinus,grWRecEff2PtMinus,grWRecEff3PtMinus,grWRecEff4PtMinus,grWRecEff5PtMinus,spreadSheetMinus,false,true);
              } //doCharge

            
            }//ieta

              Write(outFile,grWRecEff0Pt,sEffName0);
              Write(outFile,grWRecEff1Pt,sEffName1);
              Write(outFile,grWRecEff2Pt,sEffName2);
              Write(outFile,grWRecEff3Pt,sEffName3);
              Write(outFile,grWRecEff4Pt,sEffName4);
              Write(outFile,grWRecEff5Pt,sEffName5);
              if(doCharge){
                Write(outFile,grWRecEff0PtPlus,sEffNamePlus0);
                Write(outFile,grWRecEff1PtPlus,sEffNamePlus1);
                Write(outFile,grWRecEff2PtPlus,sEffNamePlus2);
                Write(outFile,grWRecEff3PtPlus,sEffNamePlus3);
                Write(outFile,grWRecEff4PtPlus,sEffNamePlus4);
                Write(outFile,grWRecEff5PtPlus,sEffNamePlus5);

                Write(outFile,grWRecEff0PtMinus,sEffNameMinus0);
                Write(outFile,grWRecEff1PtMinus,sEffNameMinus1);
                Write(outFile,grWRecEff2PtMinus,sEffNameMinus2);
                Write(outFile,grWRecEff3PtMinus,sEffNameMinus3);
                Write(outFile,grWRecEff4PtMinus,sEffNameMinus4);
                Write(outFile,grWRecEff5PtMinus,sEffNameMinus5);

              }
          }//icent
      }//ieta
} //doPt

	if(doFit){
		fitParam(grA0Plus,102,"a0");
		fitParam(grA1Plus,102,"a1");
		fitParam(grA2Plus,102,"a2");

		fitParam(grA0Minus,103,"a0");
		fitParam(grA1Minus,103,"a1");
		fitParam(grA2Minus,103,"a2");

		fitParam(grA0,104,"a0");
		fitParam(grA1,104,"a1");
		fitParam(grA2,104,"a2");
		
		TString sNamePlus = "gr";sNamePlus+=102;
		TString sNamePlusA0 = sNamePlus + "a0";
		TString sNamePlusA1 = sNamePlus + "a1";
		TString sNamePlusA2 = sNamePlus + "a2";

		TString sNameMinus = "gr";sNameMinus+=103;
		TString sNameMinusA0 = sNameMinus + "a0";
		TString sNameMinusA1 = sNameMinus + "a1";
		TString sNameMinusA2 = sNameMinus + "a2";

		TString sNameDiffA0 = "grA0Diff";
		TString sNameDiffA1 = "grA1Diff";
		TString sNameDiffA2 = "grA2Diff";

		TString sName = "gr";sName+=104;
		TString sNameA0 = sName + "a0";
		TString sNameA1 = sName + "a1";
		TString sNameA2 = sName + "a2";

		Write(outFile,grA0Plus,sNamePlusA0);
		Write(outFile,grA1Plus,sNamePlusA1);
		Write(outFile,grA2Plus,sNamePlusA2);

		Write(outFile,grA0Minus,sNameMinusA0);
		Write(outFile,grA1Minus,sNameMinusA1);
		Write(outFile,grA2Minus,sNameMinusA2);

		Write(outFile,grA0Diff,sNameDiffA0);
		Write(outFile,grA1Diff,sNameDiffA1);
		Write(outFile,grA2Diff,sNameDiffA2);

		Write(outFile,grA0,sNameA0);
		Write(outFile,grA1,sNameA1);
		Write(outFile,grA2,sNameA2);

	}
		
	
    std::cout << "Closing spreadsheet..." << std::endl;
	spreadSheet.close();
	spreadSheetPlus.close();
	spreadSheetMinus.close();
    std::cout << "Done." << std::endl;
} 

int main(){
    CorrectionFactorsCw();
    std::cout << "Done running macro." << std::endl;
    return 0;
}
