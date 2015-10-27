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

using namespace RooFit;

///////////////////////////////////////////
//writeToSpreadsheet
//////////////////////////////////////////
void writeToSpreadsheet(std::ostream& outputFile, int ieta, int icent, double scaleFactor, double errStat){
//	outputFile << ieta << " " << icent << " " << scaleFactor << " " << errStat << std::endl;
	outputFile << scaleFactor << " " << errStat << std::endl;
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
plotGeneratorDistros(TH1F* hGen,TH1F* hFid,TString sSel){

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
	hGenc->GetXaxis()->SetTitle("p_{T}^{#mu,gen}[GeV]");
	hGenc->GetYaxis()->SetTitle("Muons/GeV");
	hGenc->GetYaxis()->SetRangeUser(0.1,1.0e7);
	hGenc->GetXaxis()->SetRangeUser(0.0,250.0);

	hGenc->Scale(1,"width");
	hFidc->Scale(1,"width");

	hGenc->Draw("hist f");
	hFidc->Draw("hist fsame");
	hGenc->Draw("sameaxis");

	myText(0.196,0.862,kBlack,sSel);

	leg->Draw();

        cGen->SetLogy(true); hGenc->GetYaxis()->SetRangeUser(0.1,hGenc->GetMaximum()*7.0e2); cGen->Update();
	cGen->Print(sSel+",GeneratorDistro.pdf");

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
//Function for plotting efficiencies as a fcn of eta for rec, rec+ID, rec+ID+trig
///////////////////////////////////////////////////////////////////////////////
void calcEffEta(RooDataSet* mcWGenSet, RooDataSet* mcWFidSet,RooDataSet* mcWRecSet0, RooDataSet* mcWRecSet1,
                RooDataSet* mcWRecSet2,RooDataSet* mcWRecSet3,RooDataSet* mcWRecSet4,RooDataSet* mcWRecSet5,
                int ich, int ipt, int ieta, int icent, int index, int nCentralityBins, int nEtaBins, 
                double etaLow, double etaUpp, 
                TGraphErrors* const grWRecEff0, TGraphErrors* const grWRecEff1,TGraphErrors* const grWRecEff2,
                TGraphErrors* const grWRecEff3,TGraphErrors* const grWRecEff4,TGraphErrors* const grWRecEff5){

    ///rec eff
    double xEta = etaLow+(etaUpp-etaLow)/2.0; 

	//number of Wmu in this eta and centrality class
	double ptGen = (double)mcWGenSet->numEntries();
    std::cout << "ptGen = " << ptGen << std::endl;

	//number of muons at generator level in kinematic fiducial region within this eta and centrality class
	double ptGenInFiducial = (double)mcWFidSet->numEntries();
    std::cout << "ptGenInFiducial = " << ptGenInFiducial << std::endl;

	//number of reconstructed signal muons after pre-selection within the fiducial region in this eta and centrality class
//	double ptRecOnly = (double)mcWRecOnlySet->numEntries() ;
//    std::cout << "ptRecOnly = " << ptRecOnly << std::endl;

    ///NoQualityWselNoIsoNoZveto
    double ptRec0 = (double)mcWRecSet0->numEntries() ;
    std::cout << "ptRec0 = " << ptRec0 << std::endl;

    ///get trig eff for this bin
    double trigTemp = trigEfficiency(index);
    double trigStatErr = trigEfficiencyErr(index);

	double recoEff0 = ptRec0/ptGenInFiducial;
    recoEff0*=trigTemp;

	double recoStatErr0 = sqrt(TMath::Power(sqrt(ptRec0)/ptRec0*100,2) + TMath::Power(trigStatErr*100,2) 
            + TMath::Power(sqrt(ptGenInFiducial)/ptGenInFiducial*100,2))*0.01*recoEff0;

    std::cout << "RecoEff0: " << recoEff0 << " +- " << recoStatErr0 << std::endl;


    grWRecEff0->SetPoint(ieta,xEta,recoEff0);
    grWRecEff0->SetPointError(ieta,0.0,recoStatErr0);

    ///RecNoQualityWselNoIso
    double ptRec1 = (double)mcWRecSet1->numEntries() ;
    std::cout << "ptRec1 = " << ptRec1 << std::endl;

	double recoEff1 = ptRec1/ptGenInFiducial;
    recoEff1*=trigTemp;

	double recoStatErr1 = sqrt(TMath::Power(sqrt(ptRec1)/ptRec1*100,2) + TMath::Power(trigStatErr*100,2) 
            + TMath::Power(sqrt(ptGenInFiducial)/ptGenInFiducial*100,2))*0.01*recoEff1;

    std::cout << "RecoEff1: " << recoEff1 << " +- " << recoStatErr1 << std::endl;

    grWRecEff1->SetPoint(ieta,xEta,recoEff1);
    grWRecEff1->SetPointError(ieta,0.0,recoStatErr1);

    ///RecNoQualityWsel
    double ptRec2 = (double)mcWRecSet2->numEntries() ;
    std::cout << "ptRec2 = " << ptRec2 << std::endl;

	double recoEff2 = ptRec2/ptGenInFiducial;
    recoEff2*=trigTemp;

	double recoStatErr2 = sqrt(TMath::Power(sqrt(ptRec2)/ptRec2*100,2) + TMath::Power(trigStatErr*100,2) 
            + TMath::Power(sqrt(ptGenInFiducial)/ptGenInFiducial*100,2))*0.01*recoEff2;

    std::cout << "RecoEff2: " << recoEff2 << " +- " << recoStatErr2 << std::endl;

    grWRecEff2->SetPoint(ieta,xEta,recoEff2);
    grWRecEff2->SetPointError(ieta,0.0,recoStatErr2);

    ///RecHiQualityWselNoIsoNoZveto
    double ptRec3 = (double)mcWRecSet3->numEntries() ;
    std::cout << "ptRec3 = " << ptRec3 << std::endl;

	double recoEff3 = ptRec3/ptGenInFiducial;
    recoEff3*=trigTemp;

	double recoStatErr3 = sqrt(TMath::Power(sqrt(ptRec3)/ptRec3*100,2) + TMath::Power(trigStatErr*100,2) 
            + TMath::Power(sqrt(ptGenInFiducial)/ptGenInFiducial*100,2))*0.01*recoEff3;

    std::cout << "RecoEff3: " << recoEff3 << " +- " << recoStatErr3 << std::endl;

    grWRecEff3->SetPoint(ieta,xEta,recoEff3);
    grWRecEff3->SetPointError(ieta,0.0,recoStatErr3);

    ///RecHiQualityWselNoIso
    double ptRec4 = (double)mcWRecSet4->numEntries() ;
    std::cout << "ptRec4 = " << ptRec4 << std::endl;

	double recoEff4 = ptRec4/ptGenInFiducial;
    recoEff4*=trigTemp;

	double recoStatErr4 = sqrt(TMath::Power(sqrt(ptRec4)/ptRec4*100,2) + TMath::Power(trigStatErr*100,2) 
            + TMath::Power(sqrt(ptGenInFiducial)/ptGenInFiducial*100,2))*0.01*recoEff4;

    std::cout << "RecoEff4: " << recoEff4 << " +- " << recoStatErr4 << std::endl;

    grWRecEff4->SetPoint(ieta,xEta,recoEff4);
    grWRecEff4->SetPointError(ieta,0.0,recoStatErr4);

    ///RecHiQualityWsel
    double ptRec5 = (double)mcWRecSet5->numEntries() ;
    std::cout << "ptRec5 = " << ptRec5 << std::endl;

	double recoEff5 = ptRec5/ptGenInFiducial;
    recoEff5*=trigTemp;

	double recoStatErr5 = sqrt(TMath::Power(sqrt(ptRec5)/ptRec5*100,2) + TMath::Power(trigStatErr*100,2) 
            + TMath::Power(sqrt(ptGenInFiducial)/ptGenInFiducial*100,2))*0.01*recoEff5;

    std::cout << "RecoEff5: " << recoEff5 << " +- " << recoStatErr5 << std::endl;

    grWRecEff5->SetPoint(ieta,xEta,recoEff5);
    grWRecEff5->SetPointError(ieta,0.0,recoStatErr5);

}//calcEffEta

///////////////////////////////////////////////////////////////////////////////
//calcEffCent 
///////////////////////////////////////////////////////////////////////////////
void calcEffCent(RooDataSet* mcWGenSet, RooDataSet* mcWFidSet,RooDataSet* mcWRecSet0, RooDataSet* mcWRecSet1,
                RooDataSet* mcWRecSet2,RooDataSet* mcWRecSet3,RooDataSet* mcWRecSet4,RooDataSet* mcWRecSet5,
                int ich, int ipt, int ieta, int icent, int index, int nCentralityBins, int nEtaBins, 
                double npart,double xEta, 
                TGraphErrors* const grWRecEff0, TGraphErrors* const grWRecEff1,TGraphErrors* const grWRecEff2,
                TGraphErrors* const grWRecEff3,TGraphErrors* const grWRecEff4,TGraphErrors* const grWRecEff5, 
                TGraphErrors* const grToScale,bool doScale = false){

	//number of Wmu in this eta and centrality class
	double ptGen = (double)mcWGenSet->numEntries();
    std::cout << "ptGen = " << ptGen << std::endl;

	//number of muons at generator level in kinematic fiducial region within this eta and centrality class
	double ptGenInFiducial = (double)mcWFidSet->numEntries();
    std::cout << "ptGenInFiducial = " << ptGenInFiducial << std::endl;

    ///NoQualityWselNoIsoNoZveto
    double ptRec0 = (double)mcWRecSet0->numEntries() ;
    std::cout << "ptRec0 = " << ptRec0 << std::endl;

    ///get trig eff for this bin
    double trigTemp = trigEfficiency(index);
    double trigStatErr = trigEfficiencyErr(index);

	double recoEff0 = ptRec0/ptGenInFiducial;
    recoEff0*=trigTemp;

	double recoStatErr0 = sqrt(TMath::Power(sqrt(ptRec0)/ptRec0*100,2) + TMath::Power(trigStatErr*100,2) 
            + TMath::Power(sqrt(ptGenInFiducial)/ptGenInFiducial*100,2))*0.01*recoEff0;

    std::cout << "RecoEff0: " << recoEff0 << " +- " << recoStatErr0 << std::endl;

    grWRecEff0->SetPoint(icent,npart,recoEff0);
    grWRecEff0->SetPointError(icent,0.0,recoStatErr0);

    ///RecNoQualityWselNoIso
    double ptRec1 = (double)mcWRecSet1->numEntries() ;
    std::cout << "ptRec1 = " << ptRec1 << std::endl;

	double recoEff1 = ptRec1/ptGenInFiducial;
    recoEff1*=trigTemp;

	double recoStatErr1 = sqrt(TMath::Power(sqrt(ptRec1)/ptRec1*100,2) + TMath::Power(trigStatErr*100,2) 
            + TMath::Power(sqrt(ptGenInFiducial)/ptGenInFiducial*100,2))*0.01*recoEff1;

    std::cout << "RecoEff1: " << recoEff1 << " +- " << recoStatErr1 << std::endl;

    grWRecEff1->SetPoint(icent,npart,recoEff1);
    grWRecEff1->SetPointError(icent,0.0,recoStatErr1);

    ///RecNoQualityWsel
    double ptRec2 = (double)mcWRecSet2->numEntries() ;
    std::cout << "ptRec2 = " << ptRec2 << std::endl;

	double recoEff2 = ptRec2/ptGenInFiducial;
    recoEff2*=trigTemp;

	double recoStatErr2 = sqrt(TMath::Power(sqrt(ptRec2)/ptRec2*100,2) + TMath::Power(trigStatErr*100,2) 
            + TMath::Power(sqrt(ptGenInFiducial)/ptGenInFiducial*100,2))*0.01*recoEff2;

    std::cout << "RecoEff2: " << recoEff2 << " +- " << recoStatErr2 << std::endl;

    grWRecEff2->SetPoint(icent,npart,recoEff2);
    grWRecEff2->SetPointError(icent,0.0,recoStatErr2);

    ///RecHiQualityWselNoIsoNoZveto
    double ptRec3 = (double)mcWRecSet3->numEntries() ;
    std::cout << "ptRec3 = " << ptRec3 << std::endl;

	double recoEff3 = ptRec3/ptGenInFiducial;
    recoEff3*=trigTemp;

	double recoStatErr3 = sqrt(TMath::Power(sqrt(ptRec3)/ptRec3*100,2) + TMath::Power(trigStatErr*100,2) 
            + TMath::Power(sqrt(ptGenInFiducial)/ptGenInFiducial*100,2))*0.01*recoEff3;

    std::cout << "RecoEff3: " << recoEff3 << " +- " << recoStatErr3 << std::endl;

    grWRecEff3->SetPoint(icent,npart,recoEff3);
    grWRecEff3->SetPointError(icent,0.0,recoStatErr3);

    ///RecHiQualityWselNoIso
    double ptRec4 = (double)mcWRecSet4->numEntries() ;
    std::cout << "ptRec4 = " << ptRec4 << std::endl;

	double recoEff4 = ptRec4/ptGenInFiducial;
    recoEff4*=trigTemp;

	double recoStatErr4 = sqrt(TMath::Power(sqrt(ptRec4)/ptRec4*100,2) + TMath::Power(trigStatErr*100,2) 
            + TMath::Power(sqrt(ptGenInFiducial)/ptGenInFiducial*100,2))*0.01*recoEff4;

    std::cout << "RecoEff4: " << recoEff4 << " +- " << recoStatErr4 << std::endl;

    grWRecEff4->SetPoint(icent,npart,recoEff4);
    grWRecEff4->SetPointError(icent,0.0,recoStatErr4);

    ///RecHiQualityWsel
    double ptRec5 = (double)mcWRecSet5->numEntries() ;
    std::cout << "ptRec5 = " << ptRec5 << std::endl;

	double recoEff5 = ptRec5/ptGenInFiducial;
    recoEff5*=trigTemp;

	double recoStatErr5 = sqrt(TMath::Power(sqrt(ptRec5)/ptRec5*100,2) + TMath::Power(trigStatErr*100,2) 
            + TMath::Power(sqrt(ptGenInFiducial)/ptGenInFiducial*100,2))*0.01*recoEff5;

    std::cout << "RecoEff5: " << recoEff5 << " +- " << recoStatErr5 << std::endl;

    grWRecEff5->SetPoint(icent,npart,recoEff5);
    grWRecEff5->SetPointError(icent,0.0,recoStatErr5);

    if(doScale){

        double etaBw = abs(etaUpp-etaLow)/2.0;
        TGraphErrors* const grFracNew = new TGraphErrors();

        double* ytemp = ((TGraphErrors*)grToScale)->GetY();
        double* ytempErr = ((TGraphErrors*)grToScale)->GetEY();
        double* xtemp = ((TGraphErrors*)grToScale)->GetX();

        for(int igr=0; igr<grToScale->GetN(); igr++){

            std::cout << "Correction factors pre-scaled for eta " << xtemp[igr]-(etaBw) << "-" << xtemp[igr]+(etaBw) << " = " 
                << ytemp[igr] << std::endl;

            ///centrality binned over centrality integrated
            double scaleFactor = cWTemp/scaleFactorCentInt;
            std::cout << " s.f. = " << cWTemp << "/" << scaleFactorCentInt << " = " << scaleFactor << std::endl;
            ytemp[igr]*=scaleFactor; 
            ytempErr[igr]*=scaleFactor;
            std::cout << "S.F. in eta " << xtemp[igr]-(etaBw) << "-" << xtemp[igr]+(etaBw) 
                << " after scaling to centrality point " << centralityLow << "-" << centralityUpp << " = " 
                << " = " << ytemp[igr] << "+-" << ytempErr[igr] << std::endl;

            grFracNew->SetPoint(igr,xtemp[igr],ytemp[igr]);
            grFracNew->SetPointError(igr,0.0,ytempErr[igr]);

        }

		TString sGraphName = "correctionFactor_charge";sGraphName+=ich; sGraphName+="_eta"; sGraphName+=ieta; sGraphName+="_cent"; sGraphName+=icent;
		Write(outFile, grFracNew,sGraphName);
    }

//    writeToSpreadsheet(output,ieta,icent,cWTemp,CwStatTempErr);

}///calcEffCent

///////////////////////////////////////////////////////////////////////////////
// selectGenPtEta
///////////////////////////////////////////////////////////////////////////////
RooDataSet* selectGenPtEta( RooDataSet* dataSet, double ptLow, double ptUpp, double etaLow, double etaUpp, bool doMirrorEta = false )
{
  TString cut = "muonGenPt>";
  cut += ptLow;
  cut += "&&muonGenPt<";
  cut += ptUpp;
  cut += "&&((etaGen>";
  cut += etaLow;
  cut += "&&etaGen<";
  cut += etaUpp;
  cut += ")";
  if ( doMirrorEta ) {
    cut += "||(etaGen>";
    cut += -etaUpp;
    cut += "&&etaGen<";
    cut += -etaLow;
    cut += ")";
  }
  cut += ")";

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
RooDataSet* fillHIMuonRecSet(const TString& pathName, const TString& fileName, RooArgSet& muonArgSet)
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
  int valNt[50], ZDYNt[50], truthMatchedNt[50], promptNt[50];
  int nmu;

  // --- Load the MC tree ---
  TChain* tree = new TChain("tree","tree");
  int nFiles = tree->Add(pathName+fileName);
  
  std::cout << "Filling the RooDataSet for " << fileName << "... Number of files: " << nFiles << std::endl;
 
  // --- Set branch adresses ---
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
  tree->SetBranchAddress("nu_pt", &nu_ptNt);
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
  tree->SetBranchStatus("charge", 1);
  tree->SetBranchStatus("val", 1); 
  tree->SetBranchStatus("ZDY", 1); 
  tree->SetBranchStatus("centrality", 1);
  tree->SetBranchStatus("nu_pt", 1);
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
      muonArgSet.setRealValue("muonMt",mtNt[imu]);
      muonArgSet.setRealValue("muonEta",etaNt[imu]);
      muonArgSet.setRealValue("ZDY",ZDYNt[imu]);
      muonArgSet.setRealValue("muonCharge",chargeNt[imu]);
      muonArgSet.setRealValue("muonCategory",promptNt[imu]);
      muonArgSet.setRealValue("muonGenRecMatched",truthMatchedNt[imu]);
      if ( chargeNt[imu] > 0 ) muonArgSet.setCatLabel("chargeCategory","muPlus");
      else if ( chargeNt[imu] < 0) muonArgSet.setCatLabel("chargeCategory","muMinus");
      set->add(muonArgSet);    
   }
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
  float etaGenNt[50];
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
  tree->SetBranchAddress("mc_mu_charge", &chargeGenNt);
  tree->SetBranchAddress("mc_mu_gen_eta", &etaGenNt);
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
  tree->SetBranchStatus("mc_mu_charge", 1);
  tree->SetBranchStatus("mc_mu_gen_eta", 1);
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
      muonGenSet.setRealValue("chargeGen",chargeGenNt[imu]);
      muonGenSet.setRealValue("etaGen",etaGenNt[imu]);

      if ( chargeGenNt[imu] > 0 ) muonGenSet.setCatLabel("chargeGenCategory","muPlus");
      else if ( chargeGenNt[imu] < 0) muonGenSet.setCatLabel("chargeGenCategory","muMinus");

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
bool doScale = false;
bool doFit = false;

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
TString fileNameIn = "HISingleMuonWmunuPYTHIADataOverlay.04.12.2013";

///IMPORTANT:IF CHANGING ETA BINNING, YOU MUST CHANGE
///THE .txt FILE HOLDING THE TRIGGER EFFICIENCIES ACCORDINGLY
std::vector <double> etaBins;
//etaBins.push_back(0.0);
etaBins.push_back(0.1);
//etaBins.push_back(-2.4);
if(doEta){

/*	etaBins.push_back(+0.25);
	etaBins.push_back(+0.50);
	etaBins.push_back(+0.75);
	etaBins.push_back(+1.00);
	etaBins.push_back(+1.25);
	etaBins.push_back(+1.50);
	etaBins.push_back(+1.75);
	etaBins.push_back(+2.00);
	etaBins.push_back(+2.25);
*/

    
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

    etaBins.push_back(0.35);
    etaBins.push_back(0.6);
    etaBins.push_back(0.8);
    etaBins.push_back(1.05);
    etaBins.push_back(1.3);
    etaBins.push_back(1.55);
//        etaBins.push_back(1.73);
    etaBins.push_back(1.85);
    etaBins.push_back(2.1);




/*	etaBins.push_back(0.025);
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

if(nEtaBins==9&&nCentralityBins==6) readInputFile(nCentralityBins,nEtaBins,"triggerEffZNote_v06.txt");
else if (nEtaBins==38&&nCentralityBins==6) readInputFile(nCentralityBins,nEtaBins,"triggerEffZNote_v04.txt");
else std::cout << "WARNING: Trigger bins not initiliazed." << std::endl;
///if not mirroring eta
//readInputFile(nCentralityBins,nEtaBins,"triggerEffWNote_NoAbsEta_v01.txt");

///centrality dep plots
TGraphErrors* const grCwNpartPlus = new TGraphErrors(nCentralityBins);
TGraphErrors* const grCwNpartMinus = new TGraphErrors(nCentralityBins);
TGraph* const grCwNpartDiff = new TGraphErrors(nCentralityBins);
TGraphErrors* const grCwNpart = new TGraphErrors(nCentralityBins);

TGraphErrors* const grWRecEff0Cent = new TGraphErrors(nCentralityBins);
TGraphErrors* const grWRecEff0CentPlus = new TGraphErrors(nCentralityBins);
TGraphErrors* const grWRecEff0CentMinus = new TGraphErrors(nCentralityBins);

TGraphErrors* const grWRecEff1Cent = new TGraphErrors(nCentralityBins);
TGraphErrors* const grWRecEff1CentPlus = new TGraphErrors(nCentralityBins);
TGraphErrors* const grWRecEff1CentMinus = new TGraphErrors(nCentralityBins);

TGraphErrors* const grWRecEff2Cent = new TGraphErrors(nCentralityBins);
TGraphErrors* const grWRecEff2CentPlus = new TGraphErrors(nCentralityBins);
TGraphErrors* const grWRecEff2CentMinus = new TGraphErrors(nCentralityBins);

TGraphErrors* const grWRecEff3Cent = new TGraphErrors(nCentralityBins);
TGraphErrors* const grWRecEff3CentPlus = new TGraphErrors(nCentralityBins);
TGraphErrors* const grWRecEff3CentMinus = new TGraphErrors(nCentralityBins);

TGraphErrors* const grWRecEff4Cent = new TGraphErrors(nCentralityBins);
TGraphErrors* const grWRecEff4CentPlus = new TGraphErrors(nCentralityBins);
TGraphErrors* const grWRecEff4CentMinus = new TGraphErrors(nCentralityBins);

TGraphErrors* const grWRecEff5Cent = new TGraphErrors(nCentralityBins);
TGraphErrors* const grWRecEff5CentPlus = new TGraphErrors(nCentralityBins);
TGraphErrors* const grWRecEff5CentMinus = new TGraphErrors(nCentralityBins);

///eta dep plots
TGraphErrors* const grCwEtaPlus = new TGraphErrors(nEtaBins);
TGraphErrors* const grCwEtaMinus = new TGraphErrors(nEtaBins);
TGraph* const grCwEtaDiff = new TGraphErrors(nEtaBins);
TGraphErrors* const grCwEta = new TGraphErrors(nEtaBins);

TGraphErrors* const grWRecEff0Eta = new TGraphErrors(nEtaBins);
TGraphErrors* const grWRecEff0EtaPlus = new TGraphErrors(nEtaBins);
TGraphErrors* const grWRecEff0EtaMinus = new TGraphErrors(nEtaBins);

TGraphErrors* const grWRecEff1Eta = new TGraphErrors(nEtaBins);
TGraphErrors* const grWRecEff1EtaPlus = new TGraphErrors(nEtaBins);
TGraphErrors* const grWRecEff1EtaMinus = new TGraphErrors(nEtaBins);

TGraphErrors* const grWRecEff2Eta = new TGraphErrors(nEtaBins);
TGraphErrors* const grWRecEff2EtaPlus = new TGraphErrors(nEtaBins);
TGraphErrors* const grWRecEff2EtaMinus = new TGraphErrors(nEtaBins);

TGraphErrors* const grWRecEff3Eta = new TGraphErrors(nEtaBins);
TGraphErrors* const grWRecEff3EtaPlus = new TGraphErrors(nEtaBins);
TGraphErrors* const grWRecEff3EtaMinus = new TGraphErrors(nEtaBins);

TGraphErrors* const grWRecEff4Eta = new TGraphErrors(nEtaBins);
TGraphErrors* const grWRecEff4EtaPlus = new TGraphErrors(nEtaBins);
TGraphErrors* const grWRecEff4EtaMinus = new TGraphErrors(nEtaBins);

TGraphErrors* const grWRecEff5Eta = new TGraphErrors(nEtaBins);
TGraphErrors* const grWRecEff5EtaPlus = new TGraphErrors(nEtaBins);
TGraphErrors* const grWRecEff5EtaMinus = new TGraphErrors(nEtaBins);

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
  RooRealVar  chargeGen("chargeGen","chargeGen",-2.0,2.0);
  RooRealVar  etaGen("etaGen","etaGen",-3.0,3.0);
// --- declare cut variables at reco level --- //
  RooRealVar  muonPt("muonPt","p_{T}",0.0,250.0,"GeV");
  RooRealVar  missPt("missPt","p_{T}^{miss}",0.0,250.0,"GeV");
  RooRealVar  muonMt("muonMt","m_{T}",0.0,250.0,"GeV");
  RooRealVar  muonCharge("muonCharge","charge",-2.0,+2.0);
  RooRealVar  isolation("isolation","isolation",0.0,10.0);
  RooRealVar  centrality("centrality","centrality",0.,1.0);
  RooRealVar  muonEta("muonEta","muonEta",-3.0,+3.0);
  RooRealVar  ZDY("ZDY","ZDY",0.0,1.0);
  RooRealVar  muonQuality("muonQuality","muonQuality",0.0,15.0);
  RooRealVar  muonGenRecMatched("muonGenRecMatched","muonGenRecMatched",0.0,1.0);
  RooRealVar  muonELoss("muonELoss","muonELoss",-10.0,+10.0);
  RooRealVar  muonScat("muonScat","muonScat",-10.0,+10.0);

  TString sCutsGen = "abs(mother)==24&&abs(daughter)==13";
  RooArgList muonGenArgList(mother,daughter);
  RooFormulaVar cutsGen("cutsGen", "cutsGen", sCutsGen, muonGenArgList);
  //fiducial cuts
  //TString sCutsFid = "muonGenPt>25.0&&abs(etaGen)<2.4&&nuGenPt>25.0&&munuGenMt>40.0&&abs(mother)==24&&abs(daughter)==13";

  TString sCutsFid = "muonGenPt>25.0&&abs(etaGen)<2.4&&nuGenPt>25.0&&munuGenMt>40.0&&abs(mother)==24&&abs(daughter)==13";
//  TString sCutsGeomFid = "abs(mother)==24&&abs(daughter)==13&&abs(etaGen)<2.4";
  RooFormulaVar cutsFid("cutsFid", "cutsFid", sCutsFid, RooArgList(muonGenPt,etaGen,nuGenPt,munuGenMt,mother,daughter));
//  RooFormulaVar cutsGeomFid("cutsGeomFid","cutsGeomFid",sCutsGeomFid,RooArgList(mother,daughter,etaGen));

  //reconstruction level cuts
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
  TString sCutsRecHiQualityWsel =
      "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&isolation<0.1&&muonMt>40.0&&ZDY==0&&muonGenRecMatched==1";
      ///missPt systematic (+-10.85GeV)
      //"muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>14.15&&missPt<9000.0&&isolation<0.1&&muonMt>40.0&&ZDY==0&&muonGenRecMatched==1";
      ///use for isolation systematics
//      "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&isolation<0.2&&muonMt>40.0&&ZDY==0&&muonGenRecMatched==1";

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

  RooArgSet muonGenArgSet(muonGenPt,nuGenPt,munuGenMt,mother,daughter,chargeGenCategory,etaGen,centrality);

  RooArgSet muonRecArgSet(muonEta,centrality,ZDY,muonCategory,chargeCategory);
  muonRecArgSet.add(muonPt);
  muonRecArgSet.add(missPt);
  muonRecArgSet.add(muonMt);
  muonRecArgSet.add(muonCharge);
  muonRecArgSet.add(isolation);
  muonRecArgSet.add(muonQuality);
  muonRecArgSet.add(muonGenRecMatched);
  muonRecArgSet.add(muonELoss);
  muonRecArgSet.add(muonScat);


  double ptmax = 300.0;
  // --- Set pt and eta bins ---
  std::vector<double> ptBins;
  ptBins.push_back(0.0);
  if(doPt){
  	ptBins.push_back(50.0);
  	//ptBins.push_back(75.0);
  }
  ptBins.push_back(ptmax);
  const int nPtBins = ptBins.size()-1;

  RooDataSet* mcWGenSet = fillHIMuonGenSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonGenArgSet); mcWGenSet->Print();
  mcWGenSet = (RooDataSet*)mcWGenSet->reduce(Cut(cutsGen)); 
  std::cout << "Number of Wmunu evts at generator level : " << mcWGenSet->numEntries() << std::endl;

  RooDataSet* mcWGenSetPlus = (RooDataSet*)mcWGenSet->reduce("chargeGenCategory==chargeGenCategory::muPlus"); 
  RooDataSet* mcWGenSetMinus = (RooDataSet*)mcWGenSet->reduce("chargeGenCategory==chargeGenCategory::muMinus"); 

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

  RooDataSet* mcWFidSetPlus = (RooDataSet*)mcWFidSet->reduce("chargeGenCategory==chargeGenCategory::muPlus"); 
  RooDataSet* mcWFidSetMinus = (RooDataSet*)mcWFidSet->reduce("chargeGenCategory==chargeGenCategory::muMinus"); 

  double fidPlus = mcWFidSetPlus->numEntries(); double genPlus = mcWGenSetPlus->numEntries();
  double fidMinus = mcWFidSetMinus->numEntries(); double genMinus = mcWGenSetMinus->numEntries();

  std::cout << " Aw+ = " << mcWFidSetPlus->numEntries() << "/" << mcWGenSetPlus->numEntries() << " = " << fidPlus/genPlus << std::endl;
  std::cout << " Aw- = " << mcWFidSetMinus->numEntries() << "/" << mcWGenSetMinus->numEntries() << " = " << fidMinus/genMinus << std::endl;

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

  	plotGeneratorDistros(hGenPlus,hFidPlus,"W^{+}#rightarrow#mu^{+}#nu");
  	plotGeneratorDistros(hGenMinus,hFidMinus,"W^{-}#rightarrow#mu^{-}#nu");
  }

/*  RooDataSet* mcWRecLooseSet = fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonRecArgSet); 
  mcWRecLooseSet = (RooDataSet*)mcWRecLooseSet->reduce(Cut(cutsRecLoose)); 
  mcWRecLooseSet = (RooDataSet*)mcWRecLooseSet->reduce(Cut("muonCategory==muonCategory::W"));
  std::cout << "Number of Wmunu evts reconstructed in fiducial region before W cut selection: " << mcWRecLooseSet->numEntries() << std::endl;
*/

/*  RooDataSet* mcWRecSet = fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonRecArgSet); 
  mcWRecSet = (RooDataSet*)mcWRecSet->reduce(Cut(cutsRec)); 
  mcWRecSet = (RooDataSet*)mcWRecSet->reduce(Cut("muonCategory==muonCategory::W"));
  std::cout << "Number of Wmunu evts reconstructed in fiducial region : " << mcWRecSet->numEntries() << std::endl;
*/

  RooDataSet* mcSetWRecNoQualityWselNoIsoNoZveto =
        fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonRecArgSet); 
  mcSetWRecNoQualityWselNoIsoNoZveto =
        (RooDataSet*)mcSetWRecNoQualityWselNoIsoNoZveto->reduce(Cut(cutsRecNoQualityWselNoIsoNoZveto)); 
  mcSetWRecNoQualityWselNoIsoNoZveto = (RooDataSet*)mcSetWRecNoQualityWselNoIsoNoZveto->reduce(Cut("muonCategory==muonCategory::W"));
  std::cout << "Number of Wmunu evts reconstructed in fiducial region for reco set 0: " << mcSetWRecNoQualityWselNoIsoNoZveto->numEntries() << std::endl;

  RooDataSet* mcSetWRecNoQualityWselNoIso =
        fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonRecArgSet); 
  mcSetWRecNoQualityWselNoIso =
        (RooDataSet*)mcSetWRecNoQualityWselNoIso->reduce(Cut(cutsRecNoQualityWselNoIso)); 
  mcSetWRecNoQualityWselNoIso = (RooDataSet*)mcSetWRecNoQualityWselNoIso->reduce(Cut("muonCategory==muonCategory::W"));
  std::cout << "Number of Wmunu evts reconstructed in fiducial region for reco set 1: " << mcSetWRecNoQualityWselNoIso->numEntries() << std::endl;

  RooDataSet* mcSetWRecNoQualityWsel =
        fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonRecArgSet); 
  mcSetWRecNoQualityWsel =
        (RooDataSet*)mcSetWRecNoQualityWsel->reduce(Cut(cutsRecNoQualityWsel)); 
  mcSetWRecNoQualityWsel = (RooDataSet*)mcSetWRecNoQualityWsel->reduce(Cut("muonCategory==muonCategory::W"));
  std::cout << "Number of Wmunu evts reconstructed in fiducial region for reco set 2: " << mcSetWRecNoQualityWsel->numEntries() << std::endl;

  RooDataSet* mcSetWRecHiQualityWselNoIsoNoZveto =
        fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonRecArgSet); 
  mcSetWRecHiQualityWselNoIsoNoZveto =
        (RooDataSet*)mcSetWRecHiQualityWselNoIsoNoZveto->reduce(Cut(cutsRecHiQualityWselNoIsoNoZveto)); 
  mcSetWRecHiQualityWselNoIsoNoZveto = (RooDataSet*)mcSetWRecHiQualityWselNoIsoNoZveto->reduce(Cut("muonCategory==muonCategory::W"));
  std::cout << "Number of Wmunu evts reconstructed in fiducial region for reco set 3: " << mcSetWRecHiQualityWselNoIsoNoZveto->numEntries() << std::endl;

  RooDataSet* mcSetWRecHiQualityWselNoIso =
        fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonRecArgSet); 
  mcSetWRecHiQualityWselNoIso =
        (RooDataSet*)mcSetWRecHiQualityWselNoIso->reduce(Cut(cutsRecHiQualityWselNoIso)); 
  mcSetWRecHiQualityWselNoIso = (RooDataSet*)mcSetWRecHiQualityWselNoIso->reduce(Cut("muonCategory==muonCategory::W"));
  std::cout << "Number of Wmunu evts reconstructed in fiducial region for reco set 4: " << mcSetWRecHiQualityWselNoIso->numEntries() << std::endl;

  RooDataSet* mcSetWRecHiQualityWsel =
        fillHIMuonRecSet("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wmunu/",fileNameIn+".root",muonRecArgSet); 
  mcSetWRecHiQualityWsel =
        (RooDataSet*)mcSetWRecHiQualityWsel->reduce(Cut(cutsRecHiQualityWsel)); 
  mcSetWRecHiQualityWsel = (RooDataSet*)mcSetWRecHiQualityWsel->reduce(Cut("muonCategory==muonCategory::W"));
  std::cout << "Number of Wmunu evts reconstructed in fiducial region for reco set 5: " << mcSetWRecHiQualityWsel->numEntries() << std::endl;


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
        mcWGenSubSet[i][j][k] = selectPtEtaCentrality(mcWGenSet ,ptBins[i],ptBins[i+1], etaBins[j],
            etaBins[j+1],centBins[k], centBins[k+1],true,true); 
        std::cout << "mcWGenSubSet in " << i << ":" << j << ":" << k << " = " << mcWGenSubSet[i][j][k]->numEntries() << std::endl;

        mcWFidSubSet[i][j][k] = selectPtEtaCentrality(mcWFidSet ,ptBins[i],ptBins[i+1], etaBins[j],
            etaBins[j+1],centBins[k], centBins[k+1],true,true); 
        std::cout << "mcWFidSubSet in " << i << ":" << j << ":" << k << " = " << mcWFidSubSet[i][j][k]->numEntries()<< std::endl;
//        mcWGeomFidSubSet[i][j][k] = selectPtEtaCentrality(mcWGeomFidSet,ptBins[i],ptBins[i+1], etaBins[j],
//            etaBins[j+1],centBins[k], centBins[k+1],true,true); 
//        mcWRecLooseSubSet[i][j][k] = selectPtEtaCentrality(mcWRecLooseSet ,ptBins[i],ptBins[i+1],
//            etaBins[j],etaBins[j+1], centBins[k], centBins[k+1],true,false); 
//        mcWRecSubSet[i][j][k] = selectPtEtaCentrality(mcWRecSet ,ptBins[i],ptBins[i+1], etaBins[j],
//            etaBins[j+1],centBins[k], centBins[k+1],true,false); 

          mcSubSetWRecNoQualityWselNoIsoNoZveto[i][j][k] = selectPtEtaCentrality(mcSetWRecNoQualityWselNoIsoNoZveto,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1],true,false);
          std::cout << "mcSubSetWRecNoQualityWselNoIsoNoZveto in "<< i << ":" << j << ":" << k << " = " << mcSubSetWRecNoQualityWselNoIsoNoZveto[i][j][k]->numEntries()<< std::endl;

          mcSubSetWRecNoQualityWselNoIso[i][j][k] = selectPtEtaCentrality(mcSetWRecNoQualityWselNoIso,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1],true,false);
          std::cout << "mcSubSetWRecNoQualityWselNoIso in "<< i << ":" << j << ":" << k << " = " << mcSubSetWRecNoQualityWselNoIso[i][j][k]->numEntries()<< std::endl;

          mcSubSetWRecNoQualityWsel[i][j][k] = selectPtEtaCentrality(mcSetWRecNoQualityWsel,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1],true,false);
          std::cout << "mcSubSetWRecNoQualityWsel in "<< i << ":" << j << ":" << k << " = " << mcSubSetWRecNoQualityWsel[i][j][k]->numEntries()<< std::endl;

          mcSubSetWRecHiQualityWselNoIsoNoZveto[i][j][k] = selectPtEtaCentrality(mcSetWRecHiQualityWselNoIsoNoZveto,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1],true,false);
          std::cout << "mcSubSetWRecHiQualityWselNoIsoNoZveto in "<< i << ":" << j << ":" << k << " = " << mcSubSetWRecHiQualityWselNoIsoNoZveto[i][j][k]->numEntries()<< std::endl;

          mcSubSetWRecHiQualityWselNoIso[i][j][k] = selectPtEtaCentrality(mcSetWRecHiQualityWselNoIso,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1],true,false);
          std::cout << "mcSubSetWRecHiQualityWselNoIso in "<< i << ":" << j << ":" << k << " = " << mcSubSetWRecHiQualityWselNoIso[i][j][k]->numEntries()<< std::endl;

          mcSubSetWRecHiQualityWsel[i][j][k] = selectPtEtaCentrality(mcSetWRecHiQualityWsel,ptBins[i],ptBins[i+1],
                        etaBins[j],etaBins[j+1],centBins[k], centBins[k+1],true,false);
          std::cout << "mcSubSetWRecHiQualityWsel in "<< i << ":" << j << ":" << k << " = " << mcSubSetWRecHiQualityWsel[i][j][k]->numEntries()<< std::endl;

    }
  }
}

	  TString spreadSheetName = "correctionFactorSpreadSheet.csv";
	  std::ofstream spreadSheet;
	  spreadSheet.open(spreadSheetName);

      ///Loop for plotting efficiencies as fcn of eta in
      ///each centrality bin
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
	          int index = ieta*nCentralityBins+icent; 
                
			  std::cout << "mu^{#pm}" << std::endl;
              calcEffEta(mcWGenSubSet[ipt][ieta][icent],mcWFidSubSet[ipt][ieta][icent],
                    mcSubSetWRecNoQualityWselNoIsoNoZveto[ipt][ieta][icent], 
                    mcSubSetWRecNoQualityWselNoIso[ipt][ieta][icent],
                    mcSubSetWRecNoQualityWsel[ipt][ieta][icent],
                    mcSubSetWRecHiQualityWselNoIsoNoZveto[ipt][ieta][icent],
                    mcSubSetWRecHiQualityWselNoIso[ipt][ieta][icent],
                    mcSubSetWRecHiQualityWsel[ipt][ieta][icent],
                    104,ipt,ieta, icent, index, nCentralityBins,nEtaBins,etaLow, etaUpp,
                    grWRecEff0Eta,grWRecEff1Eta,grWRecEff2Eta,grWRecEff3Eta,grWRecEff4Eta,grWRecEff5Eta);

              if(doCharge){

		   		RooDataSet* mcWGenSetPlus  = (RooDataSet*) mcWGenSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muPlus"));
		   		RooDataSet* mcWGenSetMinus  = (RooDataSet*) mcWGenSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muMinus"));

		   		RooDataSet* mcWFidSetPlus  = (RooDataSet*) mcWFidSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muPlus"));
		   		RooDataSet* mcWFidSetMinus  = (RooDataSet*) mcWFidSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muMinus"));

/*		   		RooDataSet* mcWRecSetPlus  = (RooDataSet*) mcWRecSubSet[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   		RooDataSet* mcWRecSetMinus  = (RooDataSet*) mcWRecSubSet[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

                RooDataSet* mcWRecLooseSubSetPlus = (RooDataSet*)mcWRecLooseSubSet[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcWGeomFidSubSetPlus = (RooDataSet*)mcWGeomFidSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muPlus"));

                RooDataSet* mcWRecLooseSubSetMinus = (RooDataSet*)mcWRecLooseSubSet[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
                RooDataSet* mcWGeomFidSubSetMinus =(RooDataSet*)mcWGeomFidSubSet[ipt][ieta][icent]->reduce(Cut("chargeGenCategory==chargeGenCategory::muMinus"));
*/
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
                calcEffEta(mcWGenSetPlus,mcWFidSetPlus,
                    mcSetWRecNoQualityWselNoIsoNoZvetoPlus,mcSetWRecNoQualityWselNoIsoPlus,mcSetWRecNoQualityWselPlus,
                    mcSetWRecHiQualityWselNoIsoNoZvetoPlus, mcSetWRecHiQualityWselNoIsoPlus, mcSetWRecHiQualityWselPlus,
                    102,ipt,ieta, icent, index, nCentralityBins,nEtaBins,etaLow, etaUpp,
                    grWRecEff0EtaPlus,grWRecEff1EtaPlus,grWRecEff2EtaPlus,grWRecEff3EtaPlus,grWRecEff4EtaPlus,grWRecEff5EtaPlus);

				std::cout << "mu^{-}" << std::endl;
                calcEffEta(mcWGenSetMinus,mcWFidSetMinus,
                    mcSetWRecNoQualityWselNoIsoNoZvetoMinus,mcSetWRecNoQualityWselNoIsoMinus,mcSetWRecNoQualityWselMinus,
                    mcSetWRecHiQualityWselNoIsoNoZvetoMinus, mcSetWRecHiQualityWselNoIsoMinus, mcSetWRecHiQualityWselMinus,
                    103,ipt,ieta, icent, index, nCentralityBins,nEtaBins,etaLow, etaUpp,
                    grWRecEff0EtaMinus,grWRecEff1EtaMinus,grWRecEff2EtaMinus,grWRecEff3EtaMinus,grWRecEff4EtaMinus,grWRecEff5EtaMinus);
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

			double centbinLo = centBins.at(icent);
			double centbinUp = centBins.at(icent+1); 

			int index = ieta*nCentralityBins+icent; 
			double npart = npartBins.at(icent); 
			
			std::cout << "mu^{#pm}" << std::endl;
/*			calcEffCent(mcWGenSubSet[ipt][ieta][icent],mcWFidSubSet[ipt][ieta][icent],mcWGeomFidSubSet[ipt][ieta][icent],
                    mcWRecSubSet[ipt][ieta][icent],mcWRecLooseSubSet[ipt][ieta][icent],
                    104,ipt,ieta, icent, index, nCentralityBins,nEtaBins,npart, xEta, 
                    grWRecEffCent, grWRecIDEffCent, grWRecIDTrigEffCent, grToScale, 
                    doScale);
*/

            calcEffCent(mcWGenSubSet[ipt][ieta][icent],mcWFidSubSet[ipt][ieta][icent],
                    mcSubSetWRecNoQualityWselNoIsoNoZveto[ipt][ieta][icent], 
                    mcSubSetWRecNoQualityWselNoIso[ipt][ieta][icent],
                    mcSubSetWRecNoQualityWsel[ipt][ieta][icent],
                    mcSubSetWRecHiQualityWselNoIsoNoZveto[ipt][ieta][icent],
                    mcSubSetWRecHiQualityWselNoIso[ipt][ieta][icent],
                    mcSubSetWRecHiQualityWsel[ipt][ieta][icent],
                    104,ipt,ieta, icent, index, nCentralityBins,nEtaBins,npart, xEta,
                    grWRecEff0Cent,grWRecEff1Cent,grWRecEff2Cent,grWRecEff3Cent,grWRecEff4Cent,grWRecEff5Cent,grToScale,doScale);



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

                RooDataSet* mcSetWRecHiQualityWselPlus =
                        (RooDataSet*)mcSubSetWRecHiQualityWsel[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
                RooDataSet* mcSetWRecHiQualityWselMinus =
                        (RooDataSet*)mcSubSetWRecHiQualityWsel[ipt][ieta][icent]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));


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
                calcEffCent(mcWGenSubSet[ipt][ieta][icent],mcWFidSubSet[ipt][ieta][icent],
                        mcSubSetWRecNoQualityWselNoIsoNoZveto[ipt][ieta][icent], 
                        mcSubSetWRecNoQualityWselNoIso[ipt][ieta][icent],
                        mcSubSetWRecNoQualityWsel[ipt][ieta][icent],
                        mcSubSetWRecHiQualityWselNoIsoNoZveto[ipt][ieta][icent],
                        mcSubSetWRecHiQualityWselNoIso[ipt][ieta][icent],
                        mcSubSetWRecHiQualityWsel[ipt][ieta][icent],
                        102,ipt,ieta, icent, index, nCentralityBins,nEtaBins,npart, xEta,
                        grWRecEff0CentPlus,grWRecEff1CentPlus,grWRecEff2CentPlus,grWRecEff3CentPlus,grWRecEff4CentPlus,grWRecEff5CentPlus,grToScale,doScale);



				std::cout << "mu^{-}" << std::endl;
                calcEffCent(mcWGenSubSet[ipt][ieta][icent],mcWFidSubSet[ipt][ieta][icent],
                        mcSubSetWRecNoQualityWselNoIsoNoZveto[ipt][ieta][icent], 
                        mcSubSetWRecNoQualityWselNoIso[ipt][ieta][icent],
                        mcSubSetWRecNoQualityWsel[ipt][ieta][icent],
                        mcSubSetWRecHiQualityWselNoIsoNoZveto[ipt][ieta][icent],
                        mcSubSetWRecHiQualityWselNoIso[ipt][ieta][icent],
                        mcSubSetWRecHiQualityWsel[ipt][ieta][icent],
                        103,ipt,ieta, icent, index, nCentralityBins,nEtaBins,npart, xEta,
                        grWRecEff0CentMinus,grWRecEff1CentMinus,grWRecEff2CentMinus,grWRecEff3CentMinus,grWRecEff4CentMinus,grWRecEff5CentMinus,grToScale,doScale);

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
    std::cout << "Done." << std::endl;
} 

int main(){
	CorrectionFactorsCw();
    std::cout << "Done running macro." << std::endl;
    return 0;
}
