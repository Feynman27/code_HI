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

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooAbsArg.h"
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

#include "EfficiencyCorrection.C"
#include "Systematics.C"
#include "WPlotterHelper.C"

using namespace RooFit ;
///////////////////////////////
//plot merged histograms
//////////////////////////////
void plotMerged(THStack* hmData,THStack* hmQCD, THStack* hmZ, THStack* hmTau, THStack* hmMcSum, TString sSel, TString sSel2, int chargeIndex, 
                    bool correctSpectra=false, bool doSubtractBkg=false, bool scaleByBinWidth = false, 
                    double sigEventsTot=0.0,double statErrSqTot=0.0,double systErrSqTot=0.0){

    double mtcutLow = 40.0; 
    double mtmax = 300.0;
	//int nBins = 75; double xLo = 0.0;
    RooBinning b = RooBinning(0.0,mtmax); // 4 GeV per bin
    b.addUniform(30,0.0,120.0); //4GeV
    b.addUniform(10,120.0,200.0); //8GeV
    b.addUniform(2,200.0,mtmax); //8GeV
    double* xBins = b.array();
    int nBins = 42;

	TCanvas* cdatamt = new TCanvas("cdatamt","cdatamt",600,600);
    ///Divide canvas for plotting ratios
    std::cout << std::endl;
    std::cout << "TPad..." << std::endl;
    cdatamt->Divide(1, 2);
    TPad* canvas_up = (TPad*)cdatamt->GetListOfPrimitives()->FindObject("cdatamt_1");
    TPad* canvas_dw = (TPad*)cdatamt->GetListOfPrimitives()->FindObject("cdatamt_2");
    std::cout << "Done." << std::endl;

    ///Define the size
    double up_height     = 0.8; 
    double dw_correction = 1.60;
    double font_size_dw  = 0.08;
    double dw_height    = (1. - up_height) * dw_correction;

    ///set pad size
    canvas_up->SetPad(0., 1 - up_height, 1., 1.);
    canvas_up->SetFillColor(0);
    canvas_up->SetBorderMode(0);
    canvas_up->SetBorderSize(2);
    canvas_up->SetTickx(1);
    canvas_up->SetTicky(1);
    canvas_up->SetLeftMargin(0.16);
    canvas_up->SetRightMargin(0.05);
    canvas_up->SetTopMargin(0.05);
    canvas_up->SetBottomMargin(0.16);
    canvas_up->SetFrameBorderMode(0);
    canvas_up->SetFrameBorderMode(0);
 
    canvas_dw->SetPad(0., 0., 1., dw_height);
    canvas_dw->Range(-0.3639066,-0.7754386,2.546497,1.31186);
    canvas_dw->SetFillColor(0);
    canvas_dw->SetBorderMode(0);
    canvas_dw->SetBorderSize(2);
    canvas_dw->SetTickx(1);
    canvas_dw->SetTicky(1);
    canvas_dw->SetLeftMargin(0.159396);
    canvas_dw->SetRightMargin(0.05033557);
    canvas_dw->SetTopMargin(0.005681818);
    canvas_dw->SetBottomMargin(0.3715035);
    canvas_dw->SetFrameBorderMode(0);
    canvas_dw->SetFrameBorderMode(0);

    canvas_up->SetFrameFillColor(0);
    canvas_up->SetFillColor(0);
    canvas_dw->SetFillColor(0);
    canvas_dw->SetFrameFillColor(0);


    ///Merged, uncorrected (both Cw and Aw) histograms
	TH1F* hmDatac =(TH1F*)hmData->GetStack()->Last()->Clone("hmDatac");
	TH1F* hmQCDc =(TH1F*)hmQCD->GetStack()->Last()->Clone("hmQCDc");
	TH1F* hmZc =(TH1F*)hmZ->GetStack()->Last()->Clone("hmZc");
	TH1F* hmTauc =(TH1F*)hmTau->GetStack()->Last()->Clone("hmTauc");
	TH1F* hmMcSumc =(TH1F*)hmMcSum->GetStack()->Last()->Clone("hmMcSumc");

    // Lower mT for signal region
    int binLo = hmDatac->FindBin(40.0);
    int binUpp = hmDatac->FindBin(300.0);
    ///only for uniform binning
    //int binLo = getBinNumber(mtcutLow,nBins,xLo,mtmax);

    ///Integrate signal region of before Aw correction
  	//double sigEventsUncorr = hmDatac->Integral(binLo,nBins); 

	//std::cout << "Number of Events from QCD = "  << hmQCDc->Integral(binLo,nBins) << std::endl;
	//std::cout << "Number of Events from Z bosons = "  << hmZc->Integral(binLo,nBins) - hmQCDc->Integral(binLo,nBins) << std::endl;
    
    ///histo of Z+QCD
    //double bkgEvents = 0.0;
    //if(doSubtractBkg)
	    //bkgEvents = hmZc->Integral(binLo,nBins);

    //double sigEventsSub = sigEventsUncorr-bkgEvents;  
    //



	double Aw ;

    ///WARNING:Update the Aw 
	if(correctSpectra) Aw = getAcceptanceMt(chargeIndex,false,false,false,false);
	else Aw = 1.0;
	std::cout << "Efficiency correction factor = " << Aw << std::endl;

    ///Aw correction
    //std::cout << "Consistency check: " << sigEventsSub << " =? " << sigEventsTot << std::endl;
  	double sigEvents = sigEventsTot/Aw; 
    double totalAbsoluteStatError = TMath::Sqrt(statErrSqTot)/sigEventsTot*100.0;
    double totalPercentSystError = TMath::Sqrt(systErrSqTot)/sigEventsTot*100.0;
    // Events before background subtraction
    double nObs    = hmDatac->Integral(binLo,binUpp);
    // Individual background contributions
    double nQCDTot = hmQCDc->Integral(binLo,binUpp);
    std::cout << "Total number of QCD events: " << nQCDTot << std::endl;
    double nQCDFrac = nQCDTot/nObs*100.;
    double nMuZTot = hmZc->Integral(binLo,binUpp); nMuZTot-=nQCDTot;
    double nMuZFrac = nMuZTot/nObs*100.;
    double nTauTot = hmTauc->Integral(binLo,binUpp); nTauTot-=nMuZTot; nTauTot-=nQCDTot;
    double nTauFrac = nTauTot/nObs*100.;
    double nTotBkgFrac = (nMuZTot+nQCDTot+nTauTot)/nObs*100.;

    std::cout << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "C u m m u l a t i v e  S u m m a r y" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Number of events before background subtraction							= " << nObs << std::endl;
    std::cout << "QCD Background [%]											= " << nQCDFrac << std::endl;
    std::cout << "MuZ Background [%]											= " << nMuZFrac << std::endl;
    std::cout << "Tau Background [%]											= " << nTauFrac << std::endl;
    std::cout << "Total Background [%]											= " << nTotBkgFrac << std::endl;
    std::cout << "Number of reconstructed events after background subtraction and Cw correction, before Aw correction   = " << sigEventsTot << std::endl; 
    std::cout << "Aw Correction factor                                                                                  = " << Aw << std::endl;
    std::cout << "Number of events after Aw correction                                                                  = " << sigEvents << std::endl;
    std::cout << "Statistical uncertainty                                                                               = " << totalAbsoluteStatError << "%" << std::endl;
    std::cout << "Systematic uncertainty                                                                                = " << totalPercentSystError << "%" << std::endl;
    std::cout << std::endl;

    ///MC/Data histo
    TH1F* hRatio = (TH1F*)hmDatac->Clone("hRatio");
    hRatio = getDataMCRatio(hRatio,hmDatac,hmMcSumc);

    ///Draw spectra
    canvas_up->cd();

	hmQCDc->SetFillColor(kAzure-9);
	hmZc->SetFillColor(kRed);
	hmTauc->SetFillColor(kYellow-7);

	hmMcSumc->GetXaxis()->SetTitle("m_{T}[GeV]"); 
	hmMcSumc->GetYaxis()->SetTitle("Muons"); 
	//hmMcSumc->GetYaxis()->SetRangeUser(0.1,hmMcSumc->GetMaximum()+2.5e3); 
	hmMcSumc->GetXaxis()->SetRangeUser(mtcutLow,200.0); 
	//hmData->GetStack()->Last()->Draw("PE");
	hmMcSumc->Draw("hist f");
	hmTauc->Draw("hist fsame");
	hmZc->Draw("hist fsame");
	hmQCDc->Draw("hist fsame");
	hmDatac->Draw("pesame");
	hmMcSumc->Draw("sameaxis");

    if(scaleByBinWidth){
/*        std::cout << "Rebinning..." << std::endl;
        hmDatac->Rebin(nBins,"hmDatacc",xBins);
        hmQCDc->Rebin(nBins,"hmQCDcc",xBins);
        hmZc->Rebin(nBins,"hmZcc",xBins);
        hmTauc->Rebin(nBins,"hmTaucc",xBins);
        hmMcSumc->Rebin(nBins,"hmMcSumcc",xBins);
*/
        std::cout << "Scaling by bin width..." << std::endl;
        hmDatac->Scale(1.0,"width");
        hmQCDc->Scale(1.0,"width");
        hmZc->Scale(1.0,"width");
        hmTauc->Scale(1.0,"width");
        hmMcSumc->Scale(1.0,"width");
    }

	ATLAS_LABEL(0.17,0.89,1);
	myText(0.33,0.89, (Color_t)kBlack, (char*)("Internal"));

	TLegend* leg = new TLegend(0.658, 0.637, 0.928, 0.867);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->AddEntry(hmDatac, "Data 2011", "pe");
	leg->AddEntry(hmMcSumc, "W#rightarrow#mu#nu", "f");
	leg->AddEntry(hmTauc, "W#rightarrow#tau#nu", "f");
	leg->AddEntry(hmZc, "Z#rightarrow#mu#mu", "f");
	leg->AddEntry(hmQCDc, "QCD", "f");
	leg->Draw();

	TLatex l;
	l.SetNDC();
	l.DrawLatex(0.18,0.83,sSel +  "%" );
	l.DrawLatex(0.169,0.767,sSel2);
	l.SetTextSize(0.034);
	l.DrawLatex(0.74,0.89,"#sqrt{s_{NN}}=2.76 TeV");
	l.DrawLatex(0.492,0.89,"#int Ldt #approx 0.14 nb^{-1}"); 
	cdatamt->Update();
	

    hmMcSumc->GetYaxis()->SetRangeUser(0.1,1600.0); cdatamt->Update();

	hmMcSumc->GetXaxis()->SetRangeUser(mtcutLow,200.0); 
    ///Draw MC/Data underneath
    canvas_dw->cd();
    /// font size
    hRatio->GetXaxis()->SetRange(11,51);
    hRatio->GetXaxis()->SetTitle("m_{T}[GeV]");
    hRatio->GetXaxis()->SetLabelFont(42);
    hRatio->GetXaxis()->SetLabelSize(0.09);
    hRatio->GetXaxis()->SetTitleSize(0.11);
    hRatio->GetXaxis()->SetTitleOffset(0.8);
    hRatio->GetXaxis()->SetTitleFont(42);
    hRatio->GetYaxis()->SetTitle("Data/MC");
    hRatio->GetYaxis()->SetLabelFont(42);
    hRatio->GetYaxis()->SetLabelSize(0.09);
    hRatio->GetYaxis()->SetTitleSize(0.11);
    hRatio->GetYaxis()->SetTitleOffset(0.7);


    hRatio->GetYaxis()->SetRangeUser(0.4,1.95);
    hRatio->Draw();

    ///Save the figures
	TString plotNameLog = "dataMt_"; plotNameLog+=sSel2; plotNameLog+="_"; plotNameLog+=sSel; plotNameLog+="_"; plotNameLog+="_Merged"; plotNameLog+="Log"; 
	TString plotNameLin = "dataMt_"; plotNameLin+=sSel2; plotNameLin+="_"; plotNameLin+=sSel; plotNameLin+="_"; plotNameLin+="_Merged"; plotNameLin+="Lin"; 
	cdatamt->Print(plotNameLin.ReplaceAll("#","")+".png");
	cdatamt->Print(plotNameLin.ReplaceAll("#","")+".eps"); 
	cdatamt->Print(plotNameLin.ReplaceAll("#","")+".pdf"); 
	TString plotNameLinRoot = plotNameLin.ReplaceAll("|",",") + ".root";
	cdatamt->Print(plotNameLin); 

    ///Plot on log scale
    canvas_up->cd();
    canvas_up->SetLogy(true); 
    hmMcSumc->GetYaxis()->SetRangeUser(0.11, 3.1e4); cdatamt->Update();
	hmMcSumc->GetXaxis()->SetRangeUser(mtcutLow,200.0); cdatamt->Update(); 
	//hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()*1e3); cdatamt->Update();
	cdatamt->Print(plotNameLog.ReplaceAll("#","")+".png");
	cdatamt->Print(plotNameLog.ReplaceAll("#","")+".eps"); 
	cdatamt->Print(plotNameLog.ReplaceAll("#","")+".pdf"); 
	TString plotNameLogRoot = plotNameLog.ReplaceAll("|",",") + ".root";
	cdatamt->Print(plotNameLogRoot); 
    
    delete leg;
    delete cdatamt;

}
///////////////////////////////
//plot
//////////////////////////////
void plot(RooDataSet* dataSet,RooDataSet* mcWSet, RooDataSet* mcTauSet, RooDataSet* mcZSet, RooDataSet* mcJ1Set, RooDataSet* mcJ2Set, RooDataSet* mcJ3Set, 
		RooDataSet* mcJ1Events,RooDataSet* mcJ2Events, RooDataSet* mcJ3Events, RooDataSet* mcEventsInZSample,
        TFile* fQCDBkg,TGraphAsymmErrors* grBkgTau, TGraph2DErrors* grBkgZ, double totalZinEta, 
        double Cw, double CwStatErr, double qcdSystError, double zBosSystError, 
        THStack* hmData, THStack* hmQCD,THStack* hmZ,THStack* hmTau, THStack* hmMcSum, FitResult& fitResult, 
        const int iMt, const int iEta, const int iCentrality, int nCentralityBins,int nEtaBins, 
		const float mtmax, double ncoll, double ptLow, double ptUpp, double etaLow, double etaUpp, 
		double centralityLow, double centralityUpp, RooRealVar& muonMt, RooRealVar& centrality, 
        TString sSel, TString sSel2 , std::ostream& outputData, 
        double& sigEventsTot,double& statErrSqTot,double& systErrSqTot,bool addPercent = true, 
		bool correctSpectra = true, bool doSubtractBkg = true, bool doPreSelKinematics=false , bool doMirrorEta=true) {


    std::cout << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "B e g i n  W B o s o n  C o u n t i n g" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << std::endl;

	int nBins = 75;
	std::cout << "centrality: " << centralityLow << "-" << centralityUpp << " ncoll: " << ncoll << std::endl;

  	RooBinning b = RooBinning(nBins,0.0,mtmax); // 4 GeV per bin
/*    double binxLo = 0.0;
  	RooBinning b = RooBinning(binxLo,mtmax); // 4 GeV per bin
    b.addUniform(30,binxLo,120.0); //4GeV
    b.addUniform(5,120.0,160.0);
    b.addUniform(4,160.0,200.0); //8GeV
    b.addUniform(2,200.0,mtmax); //8GeV
    double* xBins = b.array();
    nBins = 41;
*/
	if(doPreSelKinematics) { std::cout << "Plotting with preselection cuts" << std::endl; correctSpectra=false; }

	if(correctSpectra) std::cout << "Plotting corrected mT. " << std::endl;
	//initialize histograms	
  	// --- data ---
	TH1F* hdataSet = (TH1F*)dataSet->createHistogram("hdataSet",muonMt,Binning(b));
  	// --- W set ---
	TH1F* hmcWSet = (TH1F*)mcWSet->createHistogram("hmcWSet",muonMt,Binning(b));
  	// --- Z set ---
	TH1F* hmcZSet = (TH1F*)mcZSet->createHistogram("hmcZSet",muonMt,Binning(b));
  	// --- Wtau set ---
	TH1F* hmcTauSet = (TH1F*)mcTauSet->createHistogram("hmcTauSet",muonMt,Binning(b));
  	// --- QCD set ---
	TH1F* hmcJ1Set = (TH1F*)mcJ1Set->createHistogram("hmcJ1Set",muonMt,Binning(b));
	TH1F* hmcJ2Set = (TH1F*)mcJ2Set->createHistogram("hmcJ2Set",muonMt,Binning(b));
	TH1F* hmcJ3Set = (TH1F*)mcJ3Set->createHistogram("hmcJ3Set",muonMt,Binning(b));

	//return correctly weighted QCD histogram
	TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",nBins,0.0,mtmax);
//	TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",nBins,xBins);
	hmcQCDSet = getWeightedJxHisto(centralityLow, centralityUpp, ncoll, hmcJ1Set, hmcJ2Set, hmcJ3Set, mcJ1Events, mcJ2Events,mcJ3Events, nBins,0.0, mtmax);
//	hmcQCDSet = getWeightedJxHisto(centralityLow, centralityUpp, ncoll, hmcJ1Set, hmcJ2Set, hmcJ3Set, mcJ1Events, mcJ2Events,mcJ3Events, nBins,xBins);

    ///return the lower bin number
    ///which you would like to integrate
    ///from in order to obtain signal yield
	double mtcutLow = 40.0;
	int binLo(0); 
	if(doPreSelKinematics) binLo = 1;
	else {
        binLo = getBinNumber(mtcutLow,nBins,0.0,mtmax);
        std::cout << "Consistency check: " << binLo << "=?" << hdataSet->FindBin(mtcutLow) << std::endl;
        binLo = hdataSet->FindBin(mtcutLow);
    }

    std::cout << "Signal bins = " << binLo << "-" << nBins << std::endl;

	TH1F* hdataSetc = (TH1F*)hdataSet->Clone("hdataSetc");

    ///Number of observed events in signal region
    ///before efficiency corrections
  	double sigEventsUncorr = hdataSetc->Integral(binLo,nBins); 

  	std::cout << "Integrated events before correction and background subtraction in mT:"<< mtcutLow << "-" << mtmax  << " =  " 
        << sigEventsUncorr << " +-" << TMath::Sqrt(sigEventsUncorr) << std::endl;

	double xEta = etaLow+(etaUpp-etaLow)/2.0;
	//absolute errors of uncorrected signal events 
	
//////////
//QCD
/////////

    std::cout << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Q C D  B a c k g r o u n d" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << std::endl;

	//graph of QCD bkg fraction as fcn of eta for given charge and centrality class  
	///use inclusive charge set (i.e. iMt=104)
    TString sFracQCD = "fractionQCD"; sFracQCD+="_charge"; sFracQCD+=104; sFracQCD+="_eta"; sFracQCD+=0; sFracQCD+="_cent"; 

	//hack for using 0-10% instead of 0-5%, 5-10%
	/*if(iCentrality==0)*/ sFracQCD+=iCentrality;
	/*else sFracQCD+=iCentrality-1;*/

	TGraphErrors* grBkgQCD = 0;
    if(doSubtractBkg) grBkgQCD = (TGraphErrors*) fQCDBkg->Get(sFracQCD);

	TH1F* hmcQCDSetc = (TH1F*)hmcQCDSet->Clone("hmcQCDSetc");
	
	double survivalFractionQCD ; 

	std::cout << "Getting QCD background fraction..." << std::endl;
	///use inclusive charge set (i.e. iMt=104)
	if(doPreSelKinematics) 

	    ///percent of QCD in W signal region b4 W selection
        survivalFractionQCD = getQCDBkg(104,getQCDBkgBinNumber(iEta,xEta,nEtaBins,doMirrorEta),iCentrality,
                                            grBkgQCD,centralityLow, centralityUpp, doPreSelKinematics) ; 	

	else if(doSubtractBkg) survivalFractionQCD = getQCDBkg(104,getQCDBkgBinNumber(iEta,xEta,nEtaBins,doMirrorEta), 
                                                            iCentrality,grBkgQCD,centralityLow, centralityUpp) ;
    else survivalFractionQCD = 0.0;

    std::cout << "QCD background survivalFractionQCD in centrality bin: " << iCentrality << " eta bin: " << iEta << " = " << survivalFractionQCD << std::endl;
	double dataQCDEvents = survivalFractionQCD*sigEventsUncorr;
    std::cout << "QCD background events in centrality bin: " << iCentrality << " eta bin: " << iEta << " = " << dataQCDEvents << std::endl;

	//relative error of N_corr(sigEvents) 
	/*double errNcorr = effMtSystErr/effMt;*/ 
    ///error in the difference of the QCD bkg frac for mu+ and mu-; relative error shown below
    //double errQCDFrac = qcdSystError/survivalFractionQCD;

	///Error of dataQCDEvents
	double totalUncorrelatedSystError = 0.0; 

    ///Primary statistics for relative stat. error on QCD
    double nJ1CountsAfterSel=0.0, nJ2CountsAfterSel=0.0, nJ3CountsAfterSel=0.0;
    nJ1CountsAfterSel = mcJ1Set->numEntries();
    nJ2CountsAfterSel = mcJ2Set->numEntries();
    nJ3CountsAfterSel = mcJ3Set->numEntries();
    double nJxmuCountsAfterSel = nJ1CountsAfterSel+nJ2CountsAfterSel+nJ3CountsAfterSel;
    /// stat. error on Jx
    double relErrJxmuCountsAfterSel = TMath::Sqrt(nJxmuCountsAfterSel)/nJxmuCountsAfterSel;

    //double relativeErrJxmuCountsAfterSel = relErrJxmuCountsAfterSel*dataQCDEvents;
    ///relative error from QCD Modelling systematic
    ///determined from calculating bkg fraction after
    ///scaling by Raa; value below is difference in 
    ///weighted avg 
    //double relErrJxmuModelling = 0.531;
    
    if(doSubtractBkg) {

        std::cout << "Percent statistical uncertainty from QCD = " << relErrJxmuCountsAfterSel*100. << "%" << std::endl;
    }

    ///Integrate QCD spectrum in signal region
  	double mcQCDEvents = hmcQCDSetc->Integral(binLo,nBins);
	std::cout << "QCD integral in signal region = " <<  mcQCDEvents << std::endl;

    ///Scale factor for scaling to number of
    ///expected QCD events in the data
	double sfQCD = dataQCDEvents/mcQCDEvents;

	if(mcQCDEvents==0) {
		std::cout << "WARNING: 0 QCD MC events in signal region." << std::endl;
		sfQCD = 1.0;
	}
      
    ///Normalize QCD MC shape to
    ///number of expected QCD events in
    ///the data
	hmcQCDSetc->Scale(sfQCD);

    ///Keep a running sum of
    ///QCD bkg histogram for
    ///final figure
	hmQCD->Add(hmcQCDSetc);

	hmcQCDSetc->SetFillColor(kAzure-9);
	std::cout << "Expected number of QCD events in the signal region of the data = " << dataQCDEvents << std::endl;
  	std::cout << "Background percentage from QCD in mT:"<< mtcutLow << "-" << mtmax  << " =  " <<
        dataQCDEvents/sigEventsUncorr*100.0 << "%" <<std::endl;

    totalUncorrelatedSystError = relErrJxmuCountsAfterSel*dataQCDEvents;
///////////
///Z boson
//////////

    std::cout << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Z B o s o n  B a c k g r o u n d" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << std::endl;

	std::cout << "Z integral " << hmcZSet->Integral() << std::endl;
	TH1F* hmcZSetc = (TH1F*)hmcZSet->Clone("hmcZSetc");

	double survivalProb; 
	double errStatSurvivalProbMC; 
	int indexBin = -9999;

    ///return the probability of generated Z 
    ///to survivie W selection cuts
	if(doPreSelKinematics||!doSubtractBkg) {
	    indexBin = iEta*nCentralityBins+iCentrality;
        survivalProb = getZBkg(iMt,iEta,iCentrality,indexBin,grBkgZ,centralityLow, centralityUpp,doPreSelKinematics) ;
        errStatSurvivalProbMC = getZBkgStatErrorMC(iMt,iEta,iCentrality,indexBin,grBkgZ,centralityLow,centralityUpp,doPreSelKinematics)/survivalProb ;
    }
	else if(doSubtractBkg) {

        int indexEta;

        ///bin shift necessary if using
        ///real eta (i.e. not absolute eta)
        if(!doMirrorEta){
          ///hop over the crack region
          ///in bin [-0.1,0.1]
          if(xEta==0.0) std::cout << "WARNING: Bin in crack region. Expect wrong result!" << std::endl;
          if(xEta<0.0){
              ///number of bins in absolute eta
              indexEta = indexNegativeEta(iEta,nEtaBins);
          }
           else {
              indexEta = indexPositiveEta(iEta,nEtaBins);
           }

        }
 
        ///no bin shift if doing |eta|
        else {
            indexEta = iEta;
        }

	    indexBin = indexEta*nCentralityBins+iCentrality;
        survivalProb = getZBkg(iMt,indexEta,iCentrality,indexBin, grBkgZ, centralityLow, centralityUpp) ; //percent of Z in W signal region
        ///absolute error on Z background fraction
        errStatSurvivalProbMC = 
            getZBkgStatErrorMC(iMt,indexEta,iCentrality,indexBin, grBkgZ, centralityLow,centralityUpp)/survivalProb ; //relative of Z in W signal region
    }
    else {
        survivalProb = 0.0;
        errStatSurvivalProbMC = 0.0;
    }

	std::cout <<"Z survival probability = " << survivalProb  << "+/-" << errStatSurvivalProbMC*100. << "%" << std::endl;

    ///Apply this fraction to the number of Zs in this centrality and eta bin 
    ///to afford the number of muons from Zs contaminating sample
    ///This comes from the published Z note
    double nZEventsFromData = getZEventsData(centralityLow, centralityUpp, totalZinEta);
    std::cout << "Number of Z events in this eta/centrality class from the data: " << nZEventsFromData << std::endl;
    ///uncomment for systematic study
/*    std::cout << std::endl;
    std::cout << "WARNING! WARNING! WARNING!" << std::endl;
    std::cout << "Using Z events from Monte Carlo for systematic study." << std::endl;
    nZEventsFromData = getZEventsMC(centralityLow, centralityUpp, ncoll,totalZinEta,mcEventsInZSample);
    std::cout << "Number of Z events in this eta/centrality class from the MC: " << nZEventsFromData << std::endl;
    std::cout << std::endl;
*/
	double bkgZEventFrac = survivalProb*nZEventsFromData/sigEventsUncorr;
    ///relative Stat error in Zs found from data
    double errZEventsFromData = getZCentDataStatErr(centralityLow, centralityUpp);

    std::cout << survivalProb << "*" << nZEventsFromData << "=" <<
            survivalProb*nZEventsFromData << " muons from Zs surviving W selection." << std::endl;
	std::cout << "Fraction of Zs surviving W selection cuts  = " << survivalProb*nZEventsFromData 
            << "/" << sigEventsUncorr << " = " << bkgZEventFrac << std::endl;

    ///Expected number of one-legged muons in the data
	double dataMuZEvents = bkgZEventFrac*sigEventsUncorr ;
    std::cout << "That leaves an estimated " << dataMuZEvents << " one legged-muons in the data signal region." << std::endl;
    ///Stat. error from Z
	double errZCountsStat = 0.0;

    ///Relative percent systematic errors from QCD modeling.
    ///This is used to calculate absolute error in the subtracted yield
    ///for this charge,eta,and centrality class
    double qcdPercentSystErr = getQCDRaaSystErr(iMt, indexBin) ;
    std::cout << "Systematic uncertainty of number of signal events in data due to QCD modeling = " << qcdPercentSystErr << "%" << std::endl;

    ///[%]Relative systematic errors from electroweak bkg.
    ///This is used to calculate absolute error in the subtracted yield
    ///for this charge,eta,and centrality class
    double electroweakPercentSystErr = getElectroweakSystErr(iMt, indexBin) ;
    std::cout << "Systematic uncertainty of number of signal events in data due to electroweak bkg = " << electroweakPercentSystErr << "%" << std::endl;

    if(doSubtractBkg) {

        ///Calculate absolute stat error from statistics in MC and Zs from data
        std::cout << "Error contributions to absolute Z stat: " << errZEventsFromData << "," << errStatSurvivalProbMC << std::endl;
        errZCountsStat = TMath::Sqrt( TMath::Power(errZEventsFromData,2) + TMath::Power(errStatSurvivalProbMC,2))*dataMuZEvents;  

        ///Propagate absolute stat. error from QCD and Z MC 
        totalUncorrelatedSystError = TMath::Sqrt(TMath::Power(totalUncorrelatedSystError,2) + TMath::Power(errZCountsStat,2) );
        std::cout << "Percent statistical uncertainty from Z = " << errZCountsStat*100 << "%" << std::endl;

    }

    ///Number of Z events in the signal region from MC
  	double mcZEvents = hmcZSetc->Integral(binLo,nBins);

    ///Scale factor for scaling MC events
    ///to the expected number of Z background muons
    ///in the data
	double sfZ = dataMuZEvents/mcZEvents;

	if(mcZEvents==0) {
		std::cout << "WARNING: 0 Z MC events in signal region." << std::endl;
		sfZ = 1.0;
	}

    ///Scale the histogram to the expected number of Z background muons
    ///in the data
	hmcZSetc->Scale(sfZ);

    ///Add the QCD portion to afford a "total" background
    ///histogram.
	hmcZSetc->Add(hmcQCDSetc);
	hmcZSetc->SetFillColor(kRed);

	///Keep a running sum of Z+QCD histograms for
    ///the final figure integrated over all eta and centrality
	hmZ->Add(hmcZSetc);

  	std::cout << "Background Z percentage expected in the data from mT:"<< mtcutLow << "-" << mtmax  << " =  " <<
        dataMuZEvents/sigEventsUncorr*100.0 << "%" << std::endl;

//////////////////////////////////////
//Background subtraction 
///////////////////////////////////////

    std::cout << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "S u b t r a c t  B a c k g r o u n d  " << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    std::cout << std::endl;

    ///Subtract the QCD+Z histogram from the data
    ///Expected number of reconstructed signal events in data
    ///after background subtraction.
	double sigEventsSub = 0.0;

    /// relative systematic error from 
    //QCD and EW background 
    double bkgSubPercentSystErr = 0.0;

    ///Histogram with QCD+Z component subtracted off
	TH1F* hdataSetSubtracted = (TH1F*)hdataSetc->Clone("hdataSetSubtracted");
    float dataTauEvents=0.0;
    float errStatTauEvents=0.0;
    if(doSubtractBkg) {
        sigEventsSub = sigEventsUncorr-dataQCDEvents-dataMuZEvents;
        std::cout << "Expected number of reconstructed W events in the data after QCD and Z background subtraction = " 
            << sigEventsSub << std::endl;
        
        ///Fetch %tau contamination per W->mu event in this eta/centralaty bin
        float tauBkgFraction = grBkgTau->GetY()[iEta];
        if(grBkgTau->GetEYhigh()[iEta]>grBkgTau->GetEYlow()[iEta]) errStatTauEvents=grBkgTau->GetEYhigh()[iEta];
        else errStatTauEvents=grBkgTau->GetEYlow()[iEta];
        ///Number of tau events in the data
        dataTauEvents = tauBkgFraction/(1.0+tauBkgFraction)*sigEventsSub;
        ///propagate error in tau bkg fraction from expression above
        errStatTauEvents=sigEventsSub*errStatTauEvents/TMath::Power(1.0+tauBkgFraction,2)/dataTauEvents;
        std::cout << "Number of background events from W-->tau-->mu decays: " << 
            dataTauEvents << " +/- " << errStatTauEvents << "[%]" << std::endl;
        ///absolute uncorrelated systematic error propagation (from MC stat errors)
        totalUncorrelatedSystError = TMath::Sqrt(TMath::Power(totalUncorrelatedSystError,2) + TMath::Power(errStatTauEvents,2) );
        ///convert from absolute to relative
        double totBkgCounts = dataTauEvents+dataQCDEvents+dataMuZEvents;
        totalUncorrelatedSystError /= totBkgCounts;
        ///Subtract off tau background
        sigEventsSub-=dataTauEvents;
        std::cout << "Expected number of reconstructed W events in the data after QCD,Z, and tau background subtraction = " 
            << sigEventsSub << std::endl;

        double sQCD = qcdPercentSystErr/100.0*dataQCDEvents;
        double sEW = electroweakPercentSystErr/100.0*(dataTauEvents+dataMuZEvents);
        ///[%]Relative error in subtracted yield
        bkgSubPercentSystErr = TMath::Sqrt(TMath::Power(sQCD,2)+TMath::Power(sEW,2))/sigEventsSub*100.;
    }
    else sigEventsSub = sigEventsUncorr;

    ///relative statistical errors of the raw number of counts
	double errSigStat = TMath::Sqrt(sigEventsUncorr)/sigEventsUncorr; 
    //std::cout << "Observed signal stat contribution: " << errSigStat << "," << totalAbsoluteStatError <<  std::endl;
    ///absolute stat error on the subtracted yield
    //totalAbsoluteStatError = TMath::Sqrt( TMath::Power(errSigStat,2)+TMath::Power(totalAbsoluteStatError,2));
    double totalRelativeStatError=errSigStat;

    ///Propagate abs correlated syst errror from QCD and Z into subtracted counts for current syst error on the subtracted yield
    double totalPercentSystError = 0.0;
    totalPercentSystError = bkgSubPercentSystErr;
    std::cout << "Current systematic uncertainty on the number of events in the data after QCD, Z, and tau subtraction = +/- " 
            << totalPercentSystError << "%" << std::endl;

    ///Relative systematic errors due to the isolation cut.
    ///This is used to calculate absolute error in the subtracted yield
    ///for this charge,eta,and centrality class
    double isoPercentSystErr = getIsoSystErr(iMt, indexBin);
    std::cout << "Systematic uncertainty of number of signal events in data due to isolation efficiency = " << isoPercentSystErr << "%" << std::endl;

    ///Relative systematic errors due to the MPT cut.
    ///This is used to calculate absolute error in the subtracted yield
    ///for this charge,eta,and centrality class
    double mptPercentSystErr = getMptSystErr(iMt, indexBin) ;
    std::cout << "Systematic uncertainty of number of signal events in data due to MPT resolution = " << mptPercentSystErr << "%" << std::endl;

    ///Relative systematic errors for the QCD modeling in the signal region.
    ///This was obtained by scaling the MC by the PbPb Raa for charged hardrons.
    //double qcdRaaPercentSystErr = getQCDRaaSystErr(iMt, indexBin) ;
    //std::cout << "Systematic uncertainty of number of signal events in data due to QCD modelling in signal region = " << qcdRaaPercentSystErr << "%" << std::endl;

    double cwEtaMirrorPercentSystErr = getCwMirroredEtaSystErr(iMt, indexBin) ;
    cwEtaMirrorPercentSystErr=0.0; ///removed from analysis since more of a cross-check than systematic
    //std::cout << "Systematic uncertainty of number of signal events in data due to mirroring eta in Cw calculation = " << cwEtaMirrorPercentSystErr << "%" << std::endl;

    ///Relative systematic error from tracking efficiency for pT>20GeV (from Petr)
    double trkEffPercentSystErr = 2.0;
    ///Relative systematic errors taken from Z note
    double muonRecoPercentSystErr = 2.5;
    double muonPtResolutionPercentSystErr = 3.0;
    double triggerPercentSystErr = 0.9;
    ///additional systematic for ttbar and Ztautau
    ///based on cross-section calculations
    double additionalBkgSyst = 0.5;

    ///Propagate isolation and MPT cut errors in final relative systematic 
    //error on the number of W events in the data.
    totalPercentSystError = TMath::Sqrt( TMath::Power(totalPercentSystError,2) + TMath::Power(isoPercentSystErr,2) 
        + TMath::Power(mptPercentSystErr,2) + TMath::Power(trkEffPercentSystErr,2)
        + TMath::Power(additionalBkgSyst,2) + TMath::Power(muonRecoPercentSystErr,2)
        + TMath::Power(muonPtResolutionPercentSystErr,2)+ TMath::Power(triggerPercentSystErr,2));
    //std::cout << "Statistical error on the number of events in the data after all background subtraction = +/- " 
    //        << totalRelativeStatError*100. << "%" << std::endl;
    //std::cout << "Current correlated systematic error = +/- " 
    //        << totalPercentSystError << "%" << std::endl;

/////////////////////////
//W-->tau-->mu background
/////////////////////////
    std::cout << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "W --> T a u  B a c k g r o u n d" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << std::endl;

	std::cout << "W-->tau integral " << hmcTauSet->Integral() << std::endl;
	TH1F* hmcTauSetc = (TH1F*)hmcTauSet->Clone("hmcTauSetc");
    ///Number of tau->mu in signal region from MC
    double mcTauEvents = hmcTauSetc->Integral(binLo,nBins);
    ///Scale MC events to expected
    double sfTau = dataTauEvents/mcTauEvents;

	if(mcTauEvents==0) {
		std::cout << "WARNING: 0 Tau MC events in signal region." << std::endl;
		sfTau = 1.0;
	}

    ///Scale the histogram to the expected number of Tau background muons
    ///in the data
	hmcTauSetc->Scale(sfTau);

    ///Add the QCD+Z portion to afford a "total" background
    ///histogram.
	hmcTauSetc->Add(hmcZSetc);
	hmcTauSetc->SetFillColor(kYellow-7);

    if(doSubtractBkg){
      hdataSetSubtracted->Add(hmcTauSetc,-1);
      //Consistency check: 
      std::cout << "Consistency Check: " << std::endl;
      std::cout << sigEventsSub << " =? " << hdataSetSubtracted->Integral(binLo,nBins) << std::endl;
    }
	///Keep a running sum of Z+QCD+Tau histograms for
    ///the final figure integrated over all eta and centrality
	hmTau->Add(hmcTauSetc);

  	std::cout << "Background tau percentage expected in the data from mT:"<< mtcutLow << "-" << mtmax  << " =  " <<
        dataTauEvents/sigEventsUncorr*100.0 << "%" << std::endl;

///////////////
//Cw correction
//////////////

    std::cout << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "E f f i c i e n c y  C o r r e c t i o n" << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
    std::cout << std::endl;
	double effCorrection ; double effCorrectionStatErr; 

	if(correctSpectra) {

        ///Correct background subtracted histogram with Cw
		//hdataSetc = correctEfficiencyMt(hdataSetSubtracted, 104, iCentrality, iEta,Cw);
        hdataSetSubtracted->Scale(1.0/Cw);
		effCorrection = Cw;

        ///scale absolute statistical error accordingly
	    effCorrectionStatErr = CwStatErr/Cw;
        effCorrectionStatErr=0.0;//temp hack
        ///Add in stat error from Cw
        //totalAbsoluteStatError = totalAbsoluteStatError*(1.0/Cw);
        //totalAbsoluteStatError = hdataSetSubtracted->Integral(binLo,nBins)*TMath::Sqrt(TMath::Power(relativeSigSub,2)+TMath::Power(effCorrectionStatErr,2));
        totalUncorrelatedSystError = TMath::Sqrt(TMath::Power(effCorrectionStatErr,2)+TMath::Power(totalUncorrelatedSystError,2));

	    std::cout << "Efficiency correction factor (Cw) for bin " <<  iEta << ":" << iCentrality << " = " << 
            effCorrection << " +/- " << effCorrectionStatErr << "(%) stat." << std::endl;
	}

	else {effCorrection = 1.0;}

    ///Number of background subtracted and corrected events in the data in the signal region
  	double sigEvents = hdataSetSubtracted->Integral(binLo,nBins); 
  	//double sigEvents = sigEventsSub/Cw; 

    //double cWPercentSystErr = 0.0;

    ///Add in systematic uncertainty from Cw statistical error
    //totalPercentSystError = TMath::Sqrt( TMath::Power(totalPercentSystError,2) + TMath::Power(cWPercentSystErr,2) );

    std::cout << "Uncorrelated uncertainty:  " << std::endl;
    std::cout << "QCD:                       " << relErrJxmuCountsAfterSel*100.0 << "[%]" << std::endl; 
    std::cout << "Z:                         " << errZCountsStat*100.0 << "[%]" << std::endl; 
    std::cout << "Tau:                       " << errStatTauEvents*100.0 << "[%]" << std::endl; 
    std::cout << "Cw:                        " << effCorrectionStatErr*100.0 << "[%]" << std::endl; 
    std::cout << "Current uncorrelated statistical error = +/- " 
            << totalUncorrelatedSystError*100. << "%" << std::endl;
    std::cout << "Current statistical error = +/- " 
            << totalRelativeStatError*100. << "%" << std::endl;
    std::cout << "Current systematic error = +/- " 
            << totalPercentSystError << "%" << std::endl;

	///Running sum of data histograms from 
    ///each eta and centrality bin 
	//hmData->Add(hdataSetSubtracted);
    hmData->Add(hdataSetc);
	TH1F* hmDataTemp =(TH1F*)hmData->GetStack()->Last()->Clone("hmDatac");
    std::cout << "Number of total entries in data after bin " << iEta << ":" << iCentrality << " = " << hmDataTemp->Integral() << std::endl;

//////////
//W 
/////////

    std::cout << "-----------------------------------------------------------------" << std::endl;
    std::cout << "S i g n a l  E x p e c t a t i o n  f r o m  M o n t e  C a r l o " << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
    std::cout << std::endl;

    ///We now have a histogram with subtracted data corrected. 
    ///We also have QCD and Z histograms in the same canvas (normalized to the number of 
    ///expected background events). Below, we scale the W MC to the data and see how 
    ///well they agree.

	TH1F* hmcWSetc = (TH1F*)hmcWSet->Clone("hmcWSetc");

    ///Number of W events in MC in signal region
  	double mcWEvents =  hmcWSetc->Integral(binLo,nBins);


    ///Scale factor for scaling W MC events to 
    ///the number of expected signal events in the data.
	double sfW = sigEvents/mcWEvents;

	hmcWSetc->Scale(sfW);

    ///Number of expected W events from the MC
    double mcWExpected = hmcWSetc->Integral(binLo,nBins);
    std::cout << "Number of expected W bosons from the Monte Carlo = " << mcWExpected << std::endl;

    ///Add Z,Tau, and QCD histograms to obtain a 
    ///final MC histogram composed of 
    ///the expected number of Background+Signal events
	hmcWSetc->Add(hmcTauSetc);

	//running sum of combined(QCD+Z+W) MC histos	
    ///Keep a running sum of histograms for final figure
    ///integrated over all eta and centrality
	hmMcSum->Add(hmcWSetc);	

    ///Remove all systematic sources modulo the one of interest for systematic correlation study
    //COMMENT OUT IF NOT DETERMINING CORRELATIONS
//    totalPercentSystError = 0.0;
//    totalPercentSystError = TMath::Sqrt(TMath::Power(isoPercentSystErr,2));
//    totalPercentSystError = TMath::Sqrt(TMath::Power(mptPercentSystErr,2));
//    totalPercentSystError = TMath::Sqrt(TMath::Power(cwEtaMirrorPercentSystErr,2));


    std::cout << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "S u m m a r y" << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "Expected number of events for charge:eta:centrality " << iMt << ":" << iEta << ":" << iCentrality << " = " << mcWExpected << std::endl;
    std::cout << "Number of reconstructed events before background subtraction                                          = " << sigEventsUncorr << std::endl; 
    std::cout << "Number of reconstructed events after background subtraction                                           = " << sigEventsSub << std::endl; 
    std::cout << "Correction factor                                                                                     = " << Cw << std::endl;
    std::cout << "Number of events after efficiency correction                                                          = " << sigEvents << std::endl;
    std::cout << "Statistical uncertainty                                                                               = " << totalRelativeStatError*100. << "%" << std::endl;
    std::cout << "Uncorrelated systematic uncertainty                                                                   = " << totalUncorrelatedSystError*100. << "%" << std::endl;
    std::cout << "Correlated Systematic uncertainty                                                                     = " << totalPercentSystError << "%" << std::endl;
    std::cout << std::endl;


    double totalAbsoluteStatError = totalRelativeStatError*sigEvents;
    double totalAbsoluteSystError = totalPercentSystError/100.0*sigEvents;

    TString sEtaRange[] = {"0.1-0.35","0.35-0.6","0.6-0.8","0.8-1.05","1.05-1.3","1.3-1.55","1.55-1.85","1.85-2.1","2.1-2.4"};
    TString sCentralityRange[] = {"0-5","5-10","10-15","15-20","20-40","40-80"};
    ///Write out summary to spreadsheet
	writeToSpreadsheet(outputData,sEtaRange[iEta], sCentralityRange[iCentrality],format(sigEvents,3),format(sigEventsUncorr,3), format(sigEventsSub,3),
                        format(dataQCDEvents,1),format(dataMuZEvents,1),format(dataTauEvents,1),format(totalRelativeStatError*100.,3), 
                        format(qcdPercentSystErr,2),format(electroweakPercentSystErr,2),format(additionalBkgSyst,2),
                        format(isoPercentSystErr,2),format(mptPercentSystErr,2),
                        format(trkEffPercentSystErr,2),format(muonRecoPercentSystErr,2),format(muonPtResolutionPercentSystErr,2),
                        format(triggerPercentSystErr,2),format(totalPercentSystError,3),format(totalUncorrelatedSystError*100.,3));

    ///Running sum of signal events and 
    ///sum-of-squares for errors for consistency checks
    sigEventsTot+=sigEvents;
    statErrSqTot += TMath::Power(totalAbsoluteStatError,2); 
    systErrSqTot += TMath::Power(totalAbsoluteSystError,2);

	//store sig counts+statistical errors in TGraph
    ///for later use in Rcp and asymmetry macros
    ///Write out summary into TGraph for later access

    std::cout << "Writing yields to TGraph..." << std::endl;
    ///Yield + statistical errors(hi,lo) in given bin
    fitResult.setSig(iMt, iEta, iCentrality, sigEvents, totalAbsoluteStatError, totalAbsoluteStatError);
    fitResult.setNobs(iMt, iEta, iCentrality, sigEventsUncorr, TMath::Sqrt(sigEventsUncorr), TMath::Sqrt(sigEventsUncorr));
    float nBkgTot = dataQCDEvents+dataMuZEvents+dataTauEvents;
    fitResult.setnbkg(iMt, iEta, iCentrality, nBkgTot, TMath::Sqrt(nBkgTot), TMath::Sqrt(nBkgTot));

    ///Yield + systematic errors(hi,lo) in given bin
    fitResult.setSigSyst(iMt, iEta, iCentrality, sigEvents, totalAbsoluteSystError, totalAbsoluteSystError);
    fitResult.setSigSyst1(iMt, iEta, iCentrality, sigEvents, qcdPercentSystErr/100.0*sigEvents, qcdPercentSystErr/100.*sigEvents);
    fitResult.setSigSyst2(iMt, iEta, iCentrality, sigEvents, electroweakPercentSystErr/100.0*sigEvents, electroweakPercentSystErr/100.*sigEvents);
    fitResult.setSigSyst3(iMt, iEta, iCentrality, sigEvents, additionalBkgSyst/100.0*sigEvents, additionalBkgSyst/100.*sigEvents);
    fitResult.setSigSyst4(iMt, iEta, iCentrality, sigEvents, isoPercentSystErr/100.*sigEvents, isoPercentSystErr/100.*sigEvents);
    fitResult.setSigSyst5(iMt, iEta, iCentrality, sigEvents, mptPercentSystErr/100.*sigEvents, mptPercentSystErr/100.*sigEvents);
    fitResult.setSigSyst6(iMt, iEta, iCentrality, sigEvents, trkEffPercentSystErr/100.*sigEvents, trkEffPercentSystErr/100.*sigEvents);
    fitResult.setSigSyst7(iMt, iEta, iCentrality, sigEvents, muonRecoPercentSystErr/100.*sigEvents,muonRecoPercentSystErr/100.*sigEvents);
    fitResult.setSigSyst8(iMt, iEta, iCentrality, sigEvents, muonPtResolutionPercentSystErr/100.*sigEvents, muonPtResolutionPercentSystErr/100.*sigEvents);
    fitResult.setSigSyst9(iMt, iEta, iCentrality, sigEvents, triggerPercentSystErr/100.*sigEvents, triggerPercentSystErr/100.*sigEvents);
    fitResult.setSigSyst10(iMt, iEta, iCentrality, sigEvents, totalUncorrelatedSystError*sigEvents, totalUncorrelatedSystError*sigEvents);

    std::cout << "Done." << std::endl;
////////////////////////////
//Plot intermediate figures
////////////////////////////

	//Now that the histos have been filled and scaled,
	//we now plot them on a canvas
	
//  	TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
	TLegend* leg = new TLegend(0.658, 0.637, 0.928, 0.867);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);

	leg->AddEntry(hdataSetc, "Data 2011", "pe");
	leg->AddEntry(hmcWSetc, "W#rightarrow#mu#nu", "f");
	leg->AddEntry(hmcTauSetc, "W#rightarrow#tau#nu", "f");
	leg->AddEntry(hmcZSetc, "Z#rightarrow#mu#mu", "f");
	leg->AddEntry(hmcQCDSetc, "QCD", "f");

  	TCanvas* cdatamt = new TCanvas("cdatamt","cdatamt",600,600);
  	hmcWSetc->GetXaxis()->SetTitle("m_{T}[GeV]"); 
  	hmcWSetc->GetYaxis()->SetTitle("Muons/4.0 GeV"); 
	//hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()*1.0e2); 
	//hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()+2.0e3); 
	hmcWSetc->GetXaxis()->SetRangeUser(mtcutLow,200.0); 
	hmcWSetc->Draw("hist f");
	hmcTauSetc->Draw("hist fsame");
	hmcZSetc->Draw("hist fsame");
	hmcQCDSetc->Draw("hist fsame");
	hdataSetc->Draw("pesame");
	hmcWSetc->Draw("sameaxis");
        
	TLatex l;
	l.SetNDC();
	l.DrawLatex(0.18,0.83,sSel + ( addPercent ? "%" : "" ));
	l.DrawLatex(0.169,0.767,sSel2);
	l.SetTextSize(0.034);
	l.DrawLatex(0.74,0.89,"#sqrt{s_{NN}}=2.76 TeV");
	l.DrawLatex(0.492,0.89,"#int Ldt #approx 0.14 nb^{-1}"); 

	leg->Draw(); cdatamt->Update();
	ATLAS_LABEL(0.17,0.89,1);
	myText(0.33,0.89, (Color_t)kBlack, (char*)("Internal"));

	//TString plotNameLog = "dataMt_"; plotNameLog+=sSel2; plotNameLog+="_"; plotNameLog+=sSel; plotNameLog+="_"; plotNameLog+="Log"; if(doPreSelKinematics) plotNameLog+="_PreSel";
	TString plotNameLog = "dataMt_"; plotNameLog+="charge"; plotNameLog+=iMt; plotNameLog+="_eta"; plotNameLog+=iEta; plotNameLog+="_cent"; plotNameLog+=iCentrality;  
	if(doPreSelKinematics) plotNameLog+="_PreSel";
	//TString plotNameLin = "dataMt_"; plotNameLin+=sSel2; plotNameLin+="_"; plotNameLin+=sSel; plotNameLin+="_"; plotNameLin+="Lin"; if(doPreSelKinematics) plotNameLin+="_PreSel";
	TString plotNameLin = "dataMt_"; plotNameLin+="charge"; plotNameLin+=iMt; plotNameLin+="_eta"; plotNameLin+=iEta;
    plotNameLin+="_cent"; plotNameLin+=iCentrality;  

    ///Uncomment for saving intermediate plots
  
  /*  hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()+75.0); cdatamt->Update();
	cdatamt->Print(plotNameLin.ReplaceAll("#","")+"_Lin.png");
	cdatamt->Print(plotNameLin.ReplaceAll("#","")+"_Lin.eps"); 
	cdatamt->Print(plotNameLin.ReplaceAll("#","")+"_Lin.pdf"); 
	cdatamt->Print(plotNameLin+".root"); 

  	cdatamt->SetLogy(true); hdataSetc->GetYaxis()->SetRangeUser(0.1,2.1e5); cdatamt->Update();
    cdatamt->SetLogy(true); hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()*3.0e2); cdatamt->Update();
	cdatamt->Print(plotNameLog.ReplaceAll("#","")+"_Log.png");
	cdatamt->Print(plotNameLog.ReplaceAll("#","")+"_Log.eps"); 
	cdatamt->Print(plotNameLog.ReplaceAll("#","")+"_Log.pdf"); 
	TString plotNameLogRoot = plotNameLog.ReplaceAll("|",",") + ".root";
	cdatamt->Print(plotNameLogRoot); 
*/	
	std::cout << "Clean up" << std::endl;

	delete hmcQCDSet;
    delete leg;
	delete cdatamt;	

	std::cout << "Plotting ends here." << std::endl;
}//plot

///////////////////////////////
//WAnalysis
//////////////////////////////

int WAnalysis(){
	
	bool doCharge = true ;
	bool doCentrality = true ;
	bool doEta = true; ///turn off when plotting eta distribution 
	bool correctSpectra = true; //corrects spectra using Aw,Cw factors
	bool doSubtractBkg = true; 
	bool doPreSelKinematics = false; //note: do not correctSpectra when plotting with pre-sel cuts
	bool doMirrorEta = true; ///turn on to bin in |eta|; keep on when producing eta distributions 
    //(current hack in WAnalysisHIDep until better solution is found) 	
    bool doMirrorCw = true; ///turn off to correct using Cw as function of eta (not |eta|)
    bool scaleByBinWidth = false; //turn on to divide by bin width

	//switch for plotting W candidate kinematics
	bool doPlotMt = true; 
	bool doPlotPhi = false; 
	bool doPlotEta = false; ///If turned on, turn off doEta
	bool doPlotAbsoluteEta = false; ///If turned on, turn off doEta; assign muonEta to etaAbs in fillHIMuonDataSet() 
	bool doPlotMPT = false;
	bool doPlotPt =  false;
	bool doPlotIso = false;

    ////////////////////////
    //Systematic switches///
    ////////////////////////
    ///MPT systematics; 
    //Remember to change nu_pt address in WAnalysisHIDep.C
    bool doMptSigmaDown = false;
    bool doMptSigmaUp = false;
    ///Isolation systematic study: remember to change ptconeID variable
    ///in WAnalysisHIDep.C
    bool doIncreaseConeSize = false;
    bool doLoosenIsoCut = false;
    bool doQCDModelSystematics = false;
    bool doWtauSystematics = false;

    if(doMptSigmaUp) std::cout << "Doing MPT +1 sigma systematics." << std::endl;
    if(doMptSigmaDown) std::cout << "Doing MPT -1 sigma systematics." << std::endl;
    if(doIncreaseConeSize) std::cout << "Increasing cone size from 0.2 to 0.3 for systematic study." << std::endl;
    if(doLoosenIsoCut) std::cout << "Loosening isolation cut from 0.1 to 0.2 for systematic study." << std::endl;
    if(doQCDModelSystematics) std::cout << "Doing QCD modelling systematics. " << std::endl;
    if(doWtauSystematics) std::cout << "Doing W-->tau background systematics." << std::endl;
 

	int cutValue = 11;

    double sigEvents=0.0,sigEventsPlus=0.0,sigEventsMinus=0.0;
    double statErr=0.0,statErrPlus=0.0,statErrMinus=0.0;
    double systErr=0.0,systErrPlus=0.0,systErrMinus=0.0;

    if(!doMirrorEta){
        std::cout << "Binning based on ATLAS eta (not |eta|)."<<std::endl;
    }
    if(!doMirrorCw&&correctSpectra){
        std::cout << "Applying Cw based on ATLAS eta (not |eta|)."<<std::endl;
    }
    
    if(doSubtractBkg) std::cout << "Background subtraction turned ON." << std::endl;

    if(!correctSpectra&&!doSubtractBkg) std::cout << "Plotting uncorrected, raw spectra without background subtraction."
        <<std::endl;
    else if(correctSpectra&&!doSubtractBkg) std::cout << "Plotting corrected spectra without background subtraction."
        << std::endl;
    else if(!correctSpectra) std::cout << "Plotting spectra at reconstructed level."
        << std::endl;

    if(!doPlotMt&&!doPlotPhi&&!doPlotEta&&!doPlotAbsoluteEta&&!doPlotMPT&&!doPlotPt&&!doPlotIso){
        std::cout << "Please choose a kinematic variable to plot. " << std::endl; 
        exit(0);
    }

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
	  else{ nSigPdfString = "1000000";}

	float mtmax = 300.0;
	float ptmax = 300.0;
	 
	SetAtlasStyle();
	//TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
	TString baseString = "/usatlas/u/tbales/scratch/";

	//data
	//TString fileNameDataIn = "HardProbesFiles/HISingleMuonHP.12.19.2012";
	TString fileNameDataIn = "HISingleMuonHardProbesData.04.17.2013";
	//TString fileNameDataIn = "HISingleMuonHardProbesData.07.13.2013";
    ///Closure test-->Use Wmunu MC as "data" file
	//TString fileNameDataIn = "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.04.12.2013";
	//data overlay
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


	//HIJING overlay
	TString fileNameMCZIn = "MonteCarloFiles/Zmumu/HISingleMuonZmumuPYTHIAHIJINGOverlay.04.13.2013";
//	TString fileNameMCZIn = "MonteCarloFiles/Zmumu/HISingleMuonMCZmumu.12.30.2012";
	//TString fileNameMCZIn = baseString+"MonteCarloFiles/Zmumu/HISingleMuonMC_PYTHIA_HIJING_Zmumu_11.28.2012";
    
    ///W-->tau-->mu Generated sample
    TString fileNameMCTauIn = "MonteCarloFiles/Wtaumu/HIWtaumuNtuple.07.30.2013";

	//J1 1 muon-filter 
	//TString fileNameMCJ1In = "MonteCarloFiles/QCD/HISingleMuonMcJ1.2012.11.29";
	TString fileNameMCJ1In = "MonteCarloFiles/QCD/HISingleMuonJ1_1muPYTHIADataOverlay.04.13.2013";
	//J2 1 muon-filter 
	//TString fileNameMCJ2In = "MonteCarloFiles/QCD/HISingleMuonMcJ2.2012.11.29";
	TString fileNameMCJ2In = "MonteCarloFiles/QCD/HISingleMuonJ2_1muPYTHIADataOverlay.04.13.2013";
	//J3 1 muon-filter 
	//TString fileNameMCJ3In = "MonteCarloFiles/QCD/HISingleMuonMcJ3.2012.11.29";
	TString fileNameMCJ3In = "MonteCarloFiles/QCD/HISingleMuonJ3_1muPYTHIADataOverlay.04.13.2013";

	//files that hold background values
//	TString fileNameBkgQCDIn = "fractionQCDEtaScaled.12.30.2012";
	///use this file for inclusive charge set
//	TString fileNameBkgQCDIn = "background/fractionQCDEtaScaled.01.21.2013";
    ///11 variable eta bins
	//TString fileNameBkgQCDIn = "background/fractionQCDEtaScaled_variableEtaBinning.04.04.2013";

    TString fileNameMptMtCorrelation = "mptMtCorrelations.07.21.2013";
    TFile* _fMptMtCorrelation = NULL;
    TProfile* _pfxMptMtCorrelation = NULL;
    TString fileNameIn_pp,fileNameIn_pn,fileNameIn_np,fileNameIn_nn;
    if(doMptSigmaDown||doMptSigmaUp){

        _fMptMtCorrelation = new TFile(fileNameMptMtCorrelation+".root","read");
        if(_fMptMtCorrelation!=0) std::cout << "Mpt-Mt correlation file: " << fileNameMptMtCorrelation << " opened." << std::endl;
        else exit(0);

        fileNameDataIn = "HISingleMuonHardProbesData.07.13.2013";
        fileNameIn_pp = "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_pp.07.13.2013";
        fileNameIn_pn = "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_pn.07.13.2013";
        fileNameIn_np = "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_np.07.13.2013";
        fileNameIn_nn = "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_nn.07.13.2013";
        fileNameMCZIn = "MonteCarloFiles/Zmumu/HISingleMuonZmumuPYTHIAHIJINGOverlay.07.13.2013";
        fileNameMCJ1In = "MonteCarloFiles/QCD/HISingleMuonJ1_1muPYTHIADataOverlay.07.13.2013";
        fileNameMCJ2In = "MonteCarloFiles/QCD/HISingleMuonJ2_1muPYTHIADataOverlay.07.13.2013";
        fileNameMCJ3In = "MonteCarloFiles/QCD/HISingleMuonJ3_1muPYTHIADataOverlay.07.13.2013";

/*        fileNameDataIn = "HardProbesFiles/HISingleMuonHardProbesData.07.19.2013";
        fileNameIn_pp = "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_pp.07.19.2013";
        fileNameIn_pn = "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_pn.07.19.2013";
        fileNameIn_np = "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_np.07.19.2013";
        fileNameIn_nn = "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_nn.07.19.2013";
        fileNameMCZIn = "MonteCarloFiles/Zmumu/HISingleMuonZmumuPYTHIAHIJINGOverlay.07.19.2013";
        fileNameMCJ1In = "MonteCarloFiles/QCD/HISingleMuonJ1_1muPYTHIADataOverlay.07.19.2013";
        fileNameMCJ2In = "MonteCarloFiles/QCD/HISingleMuonJ2_1muPYTHIADataOverlay.07.19.2013";
        fileNameMCJ3In = "MonteCarloFiles/QCD/HISingleMuonJ3_1muPYTHIADataOverlay.07.19.2013";
        */
    }
/////////////////////////
///QCD background files
/////////////////////////

    ///use for analysis
	TString fileNameBkgQCDIn = "";
    if(doMptSigmaUp) fileNameBkgQCDIn = "systematics/fractionQCDEtaScaled_4GeVMpt_systematics.08.03.2013";
    else if(doMptSigmaDown) fileNameBkgQCDIn = "systematics/fractionQCDEtaScaled_2GeVMpt_systematics.08.03.2013";
    //else fileNameBkgQCDIn = "background/fractionQCDEtaScaled_9etaBins6centBins.04.21.2013";
    //else fileNameBkgQCDIn = "background/fractionQCDEtaScaled_9etaBins6centBins.08.20.2013";
    ///nominal file used for f_qcd 
    else fileNameBkgQCDIn = "background/fractionQCDEtaScaled_9etaBins6centBins.11.27.2013";
    //cross-check w/o mpt cut
    //else fileNameBkgQCDIn = "crossChecks/fractionQCDEtaScaled_noMptCut.12.10.2013";
    ///isolation efficiency study
    //else fileNameBkgQCDIn = "background/fractionQCDEtaScaled_NoIsoCut_9etaBins6centBins.06.24.2013";
    ///QCD systematics for shape in signal region
    if(doQCDModelSystematics)
        fileNameBkgQCDIn = "systematics/fractionQCDEtaScaled_Raa_06.26.2013";
    ///mpt systematics 
//    TString fileNameBkgQCDIn = "systematics/fractionQCDEtaScaled_MPTup_systematics.04.24.2013";
//    TString fileNameBkgQCDIn = "systematics/fractionQCDEtaScaled_MPTdown_systematics.04.24.2013";
    ///systematics for ptcone
    ///ptcone30ID3<0.1
    if(doIncreaseConeSize)
	    fileNameBkgQCDIn = "systematics/fractionQCDEtaScaled_ptcone30ID3_systematics.06.25.2013";
    ///ptcone20ID3<0.2
    if(doLoosenIsoCut)     
	    fileNameBkgQCDIn = "systematics/fractionQCDEtaScaled_ptcone20ID3_systematics.06.25.2013";

    std::cout << "QCD background file : " << fileNameBkgQCDIn << std::endl;
    TFile *fQCDBkg = new TFile(fileNameBkgQCDIn+".root","READ");
	if ( !fQCDBkg->IsOpen() ) {
	    std::cout << fQCDBkg << " not found!" << std::endl;
	    exit(0);
	}
	//qcd bkg syst as fcn of ncoll
	TString fileNameBkgQCDSystNcoll = "background/fractionQCDCent_6Bins.06.17.2013";
    TFile *fQCDBkgSystNcoll = new TFile(fileNameBkgQCDSystNcoll+".root","READ");
	if ( !fQCDBkgSystNcoll->IsOpen() ) {
	    std::cout << fQCDBkgSystNcoll << " not found!" << std::endl;
	    exit(0);
	}
	//qcd bkg syst as fcn of eta
	TString fileNameBkgQCDSystEta = "background/fractionQCDEta_9etaBins6centBins.04.21.2013";
    ///11 eta bins
//	TString fileNameBkgQCDSystEta = "background/fractionQCDEta_variableEtaBinning.04.04.2013";
    TFile *fQCDBkgSystEta = new TFile(fileNameBkgQCDSystEta+".root","READ");
	if ( !fQCDBkgSystEta->IsOpen() ) {
	    std::cout << fQCDBkgSystEta << " not found!" << std::endl;
	    exit(0);
	}

	//TString fileNameBkgZIn = "fractionZEtaCent.12.30.2012";
	///use this file for inclusive charge set
//	TString fileNameBkgZIn = "background/fractionZEtaCent.02.11.2013";
    //10 eta bins
//	TString fileNameBkgZIn = "background/fractionZEtaCent.03.31.2013";
    ///11 variable eta bins
//	TString fileNameBkgZIn = "background/fractionZEtaCent_variableEtaBinning.04.04.2013";


//////////////////////
///Wtau background files
//////////////////////
	TString fileNameBkgTauIn = "";
    if(doWtauSystematics) fileNameBkgTauIn = "systematics/fractionTauEtaCent_9etaBins6centBins.07.30.2013";
    //else fileNameBkgTauIn = "background/fractionTauEtaCent_9etaBins6centBins.07.30.2013";
    //nominal file
    else fileNameBkgTauIn = "background/fractionTauEtaCent_9etaBins6centBins.11.27.2013";
    //cross-check w/o mpt cut
    //else fileNameBkgTauIn = "crossChecks/fractionTauEtaCent_noMptCut.12.10.2013";
    std::cout << " Wtau background file : " << fileNameBkgTauIn << std::endl;

    TFile *fBkgTau = new TFile(fileNameBkgTauIn+".root","READ");
    if ( !fBkgTau->IsOpen() ) {
        std::cout << fBkgTau << " not found!" << std::endl;
        exit(0);
    }

//////////////////////
///Z background files
//////////////////////

   ///use for analysis
    TString fileNameBkgZIn = "";
    if(doMptSigmaUp) fileNameBkgZIn = "systematics/fractionZEtaChargeCent_4GeVMpt_systematics.08.03.2013";
    else if(doMptSigmaDown) fileNameBkgZIn ="systematics/fractionZEtaChargeCent_2GeVMpt_systematics.08.03.2013";
    //else fileNameBkgZIn = "background/fractionZEtaChargeCent_9etaBins6centBins.04.21.2013";
    //else fileNameBkgZIn = "background/fractionZEtaChargeCent_9etaBins6centBins.06.25.2013";
    //else fileNameBkgZIn = "background/fractionZEtaChargeCent_9etaBins6centBins.07.17.2013";
    //nominal file
    else fileNameBkgZIn = "background/fractionZEtaChargeCent_9etaBins6centBins.11.27.2013";
    // cross-check w/ data-overlay (~20% the stats)
    //else fileNameBkgZIn = "background/fractionZEtaChargeCent_dataOverlay.05.18.2014";
    // cross-check w/o mpt
    //else fileNameBkgZIn = "crossChecks/fractionZEtaChargeCent_noMptCut.12.10.2013";
    ///mpt systematics 
//    TString fileNameBkgZIn ="systematics/fractionZEtaCent_MPTup_systematics.04.24.2013";
//    TString fileNameBkgZIn ="systematics/fractionZEtaCent_MPTdown_systematics.04.24.2013";
    ///systematics with ptcone30ID3<0.1
    if(doIncreaseConeSize)
	 fileNameBkgZIn = "systematics/fractionZEtaCent_ptcone30ID3_systematics.06.25.2013";
    ///systematics with ptcone20ID3<0.2
    if(doLoosenIsoCut)     
	 fileNameBkgZIn = "systematics/fractionZEtaCent_ptcone20ID3_systematics.06.25.2013";

    std::cout << " Z background file : " << fileNameBkgZIn << std::endl;
	TString sFracZ = "fractionZEtaCent";
	TString sFracZPlus = "fractionZEtaCentPlus";
	TString sFracZMinus = "fractionZEtaCentMinus";
	TFile *fZBkg = new TFile(fileNameBkgZIn+".root","READ");
	if ( !fZBkg->IsOpen() ) {
	    std::cout << fZBkg << " not found!" << std::endl;
	    exit(0);
	}

	//files with systmatics for Z
    // Nominal file for analysis
	TString fileNameBkgZSystNcoll = "background/fractionZCent.06.25.2013";
    // cross-check w/ data-overlay
	//TString fileNameBkgZSystNcoll = "background/fractionZCent_dataOverlay.05.18.2014";
    ///10 eta bins
	//TString fileNameBkgZEta = "background/fractionZEta.03.31.2013";
    ///11 eta bins
//	TString fileNameBkgZEta = "background/fractionZEta_variableEtaBinning.04.04.2013";
    ///15 eta bins
	//TString fileNameBkgZEta = "background/fractionZEta.06.25.2013";
    // Nominal file for analysis
	TString fileNameBkgZEta = "background/fractionZEta.11.27.2013";
    // cross-check w/ data-overlay
	//TString fileNameBkgZEta = "background/fractionZEta_dataOverlay.05.18.2014";

	TFile *fBkgZSystNcoll = new TFile(fileNameBkgZSystNcoll+".root","READ");
	TFile *fBkgZEta = new TFile(fileNameBkgZEta+".root","READ");

	///open file that holds TGraph of Cw parametrization values
//	TString fileNameCw = "correctionFactorsW.02.11.2013";
    ///10 eta bins
	TString fileNameCwParam = "correctionFactorsW.02.25.2013";
    ///use for yields
    ///10 eta bins
//	TString fileNameCw = "CorrectionFactorFiles/correctionFactorsW_etaCent_binbybin.04.03.2013";


//////////////////////
///Correction factor (CW,AW) files
//////////////////////
    ///use for analysis
	TString fileNameCw = "";
    //if(doMptSigmaUp)fileNameCw = "systematics/correctionFactorsW_4GeVMpt_systematics.08.03.2013";
    if(doMptSigmaUp)fileNameCw = "systematics/correctionFactorsW_4GeVMpt_systematics.11.20.2013";
    //else if(doMptSigmaDown) fileNameCw = "systematics/correctionFactorsW_2GeVMpt_systematics.08.03.2013";
    else if(doMptSigmaDown) fileNameCw = "systematics/correctionFactorsW_2GeVMpt_systematics.11.20.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges.05.25.2013";
    // Py6
    //else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges.07.12.2013";
    // PowPy8
    //else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges.PowPy8.10.03.2013";
    //else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges.PowPy8.11.27.2013";
    ///NOMINAL input for correction factors
    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_extrapolated.PowPy8.12.02.2013";
    // correct with Cw*Aw
//    else fileNameCw = "CorrectionFactorFiles/correctionFactor_binbybin_9VarEtaBins6CentBins2Charges_AwCw.08.14.2013";
//    cross-check with FCal weighting
    //else fileNameCw = "crossChecks/correctionFactorsW_FcalWeighted.06.12.2014";
    //Diagnostics
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_pp.09.27.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_noMtnoMpt.09.22.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_mpt20nupt15.09.04.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_noTrigger.08.13.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_noPreSelection.08.13.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_QualityCuts_noDIFVeto.08.13.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_noMPTcut.08.13.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_noRecGenMatching.08.13.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_Qual1.08.29.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_Qual2.08.29.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_Qual3.08.29.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_Qual4.08.29.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_Qual5.08.29.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_Qual6.08.29.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_Qual7.08.29.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_Qual8.08.29.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_Qual9.08.29.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_Qual10.08.29.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_Qual11.08.29.2013";
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_Qual12.08.29.2013";

//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_sanityCheck.07.14.2013";
    ///with cut on neutrino eta at generator level
//    else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges.06.08.2013";
    if(!doMirrorCw) fileNameCw = "CorrectionFactorFiles/correctionFactorsWEtaCent_19EtaBinsNoAbsEta6CentBins2Charges.05.26.2013";
    //else fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges.05.03.2013";
    ///new prompt definition
    //else fileNameCw =
    //    "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_IncorrectErrors.05.25.2013";
    //else fileNameCw = "CorrectionFactorFiles/correctionFactorsWEtaCent_9VarEtaBins6CentBins2Charges.05.12.2013";

    //contains only charge-inclusive cw by mistake
//	TString fileNameCw = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges.04.23.2013";
    ///pt inc by 10GeV (check)
	//TString fileNameCw = "CorrectionFactorFiles/correctionFactorsW_pt10up_9VarEtaBins6CentBins2Charges.04.29.2013";
    ///pt inc by 15GeV (check)
	//TString fileNameCw = "CorrectionFactorFiles/correctionFactorsW_pt15up_9VarEtaBins6CentBins2Charges.04.29.2013";
    ///mpt systematics
    //TString fileNameCw = "systematics/correctionFactorsW_MPTup_systematics.04.24.2013";
//    TString fileNameCw = "systematics/correctionFactorsW_MPTdown_systematics.04.24.2013";
    ///use for ptcone30ID3<0.1 systematics
    if(doIncreaseConeSize)
	 fileNameCw = "systematics/correctionFactorsW_ptCone30ID3_systematics.04.24.2013";
    ///use for ptcone20ID3<0.2 systematics
    if(doLoosenIsoCut)     
	 fileNameCw = "systematics/correctionFactorsW_ptCone20ID3_systematics.06.25.2013";
    ///Trigger Efficiencies do not correspond to eta/cent bin; 
    ///these files are only used for systematics so OK for now, but remember to update later!!!
	TString fileNameCwEta = "CorrectionFactorFiles/correctionFactorsWEta_15etaBins.04.07.2013";
	TString fileNameCwCent = "CorrectionFactorFiles/correctionFactorsW_centChargeBinned.04.04.2013";

    //corrections without trigger correction
//	TString fileNameCw = "CorrectionFactorFiles/correctionFactorsW_RecIDNoTrig_binbybin.04.03.2013";
    ///correction with pp binning
//    TString fileNameCw = "correctionFactorsW";
    ///use for isolation systematics
//	TString fileNameCw = "systematics/correctionFactorsW_ptCone30ID3_systematics.04.03.2013";
//	TString fileNameCw = "systematics/correctionFactorsW_ptCone20ID3_systematics.04.03.2013";
	//Cw parametrization for +-1Sigma systematics
//	TString fileNameCw = "systematics/eff/ptSigmaUp/correctionFactorsW_pt1SigmaUp"; //+1pt
//	TString fileNameCw = "systematics/eff/ptSigmaDown/correctionFactorsW_pt1SigmaDown"; //-1pt
//	TString fileNameCw = "systematics/eff/mptSigmaUp/correctionFactorsW_mpt1SigmaUp"; //+1mpt
//	TString fileNameCw = "systematics/eff/mptSigmaDown/correctionFactorsW_mpt1SigmaDown"; //-1mpt

    if(correctSpectra) {
        std::cout << "Correction factors turned ON." << std::endl;
        std::cout << "Filename: " << fileNameCw << std::endl;
    }

	//open files that hold TGraphs of Cw systematic errors
    //use these for Cw parametrization technique
    ///10 eta bins
	TString fileNameCwSyst0 = "systematics/eff/correctionFactorsW.AllFloat_syst.01.19.2013"; 
	TString fileNameCwSyst1 = "systematics/eff/correctionFactorsW.a2_fixed_syst.01.19.2013"; 
	TString fileNameCwSyst2 = "systematics/eff/correctionFactorsW.a1_a2_fixed_syst.01.19.2013"; 
    
    //use these for bin by bin corrections
//	TString fileNameCwSyst0 = "systematics/"; 
//	TString fileNameCwSyst1 = "systematics/"; 

	TFile* fCw = new TFile(fileNameCw+".root","READ");
	TFile* fCwEta = new TFile(fileNameCwEta+".root","READ");
	TFile* fCwCent = new TFile(fileNameCwCent+".root","READ");

    ///files used when using parameterization techinique for Cw
	TFile* fCwParam = new TFile(fileNameCwParam+".root","READ");
	TFile* fCwSyst0 = new TFile(fileNameCwSyst0+".root","READ");
	TFile* fCwSyst1 = new TFile(fileNameCwSyst1+".root","READ");
	TFile* fCwSyst2 = new TFile(fileNameCwSyst2+".root","READ");

	std::cout << "All input files open." << std::endl;

  	//graph of Z fraction for given charge,eta,and centrality class  
	///using only charge inclusive bkg values
  	TGraph2DErrors* grBkgZ = (TGraph2DErrors*) fZBkg->Get(sFracZ);
  	//TGraph2DErrors* grBkgZPlus = (TGraph2DErrors*) fZBkg->Get(sFracZPlus);
  	TGraph2DErrors* grBkgZPlus = (TGraph2DErrors*) fZBkg->Get(sFracZ);
  	//TGraph2DErrors* grBkgZMinus = (TGraph2DErrors*) fZBkg->Get(sFracZMinus);
  	TGraph2DErrors* grBkgZMinus = (TGraph2DErrors*) fZBkg->Get(sFracZ);


	TList _cWsystGraphs; //TList for charge difference in Cw
	TList _cWResiduals; //TList for residuals in Cw
	TList _cWsystGraphAx; //TList for charge difference in Cw parameterization coefficients (a0,a1,a2);

	//TGrphs holding systematic errors 
	//Includes charge diff and deviation of qcd frac from lin
	TGraphErrors* grBkgQCDSystResidNcoll = (TGraphErrors*) fQCDBkgSystNcoll->Get("fitQCDResidual") ;
	TGraphErrors* grBkgQCDSystChDiffNcoll = (TGraphErrors*) fQCDBkgSystNcoll->Get("diffFracQCDNcoll") ;
	TGraphErrors* grBkgQCDSystChDiffEta = (TGraphErrors*) fQCDBkgSystEta->Get("diffFracQCDEta") ;
	
	//charge diff for Z distros
	TGraphErrors* grBkgZSystNcoll = (TGraphErrors*) fBkgZSystNcoll->Get("diffFracZNcoll");
	TGraphErrors* grBkgZNcoll = (TGraphErrors*) fBkgZSystNcoll->Get("fractionZCent");
	TGraphErrors* grBkgZSystEta = (TGraphErrors*) fBkgZEta->Get("diffFracZEta");

    //fraction of generated muons from Zs in particular eta slice
	TGraphErrors* grTotZEta = (TGraphErrors*) fBkgZEta->Get("fractionTotalZEta");
	TGraphErrors* grTotZEtaPlus = (TGraphErrors*) fBkgZEta->Get("fractionTotalZEtaPlus");
	TGraphErrors* grTotZEtaMinus = (TGraphErrors*) fBkgZEta->Get("fractionTotalZEtaMinus");

	TString fileNameFitOut = "WAnalysis_fitResult"; fileNameFitOut+=nSigPdfString; fileNameFitOut+=".root";
	TString spreadSheetName = "dataSpreadSheet.csv";
	TString spreadSheetNameMuPlus = "dataSpreadSheetMuPlus.csv";
	TString spreadSheetNameMuMinus = "dataSpreadSheetMuMinus.csv";

	std::ofstream spreadSheet;
	spreadSheet.open(spreadSheetName);

	std::ofstream spreadSheetMuPlus;
	spreadSheetMuPlus.open(spreadSheetNameMuPlus);

	std::ofstream spreadSheetMuMinus;
	spreadSheetMuMinus.open(spreadSheetNameMuMinus);

        // --- declare cut variables --- //
	 RooRealVar  muonPt("muonPt","p_{T}",0.0,400.0,"GeV");
	 RooRealVar  motherRec("motherRec","motherRec",0.0,400.0);
	 RooRealVar  missPt("missPt","p_{T}^{miss}",0.0,400.0,"GeV");
	 RooRealVar  muonMt("muonMt","m_{T}",0.0,400.0,"GeV");
	 RooRealVar  muonCharge("muonCharge","charge",-2.0,+2.0);
	 RooRealVar  isolationMu("isolationMu","isolationMu",0.0,10.0);
  	 RooRealVar  centrality("centrality","centrality",0.,1.0);
  	 RooRealVar  muonEta("muonEta","muonEta",-3.0,+3.0);
  	 RooRealVar  muonPhi("muonPhi","muonPhi",-3.5,+3.5);
  	 RooRealVar  ZDY("ZDY","ZDY",0.0,2.0);
     RooRealVar  muonQuality("muonQuality","muonQuality",0.0,15.0);
     RooRealVar  muonELoss("muonELoss","muonELoss",-10.0,+10.0);
     RooRealVar  muonScat("muonScat","muonScat",-10.0,+10.0);
     RooRealVar  muonGenRecMatched("muonGenRecMatched","muonGenRecMatched",0.0,2.0);
     ///RooRealVar for weighting W np,pn,pp,nn samples
     RooRealVar  w("w","w",0.0,10.0);

	  ///W selection cuts
      ///match with correction factor that was chosen
//	  TString sCutsSig = "abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0";
//	  TString sCutsSig = "abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0&&ZDY==0";
//	  TString sCutsSig = "abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0&&ZDY==0&&isolationMu<0.1";
//	  TString sCutsSig =
//          "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0";
//	  TString sCutsSig =
//          "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0&&ZDY==0";

      ///Analysis cuts for signal region
	  TString sCutsSig = "";
          if(doMptSigmaDown||doMptSigmaUp)
              sCutsSig=
                "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&isolationMu<0.1&&muonMt>40.0&&ZDY==0";
          else 
              sCutsSig=
                  // NOMINAL
                  "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)>0.1&&abs(muonEta)<2.4&&centrality<0.8&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&isolationMu<0.1&&muonMt>40.0&&ZDY==0";
                  //cross-check w/o mpt cut
                  //"muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)>0.1&&abs(muonEta)<2.4&&centrality<0.8&&muonPt>25.0&&missPt>0.0&&missPt<9000.0&&isolationMu<0.1&&muonMt>40.0&&ZDY==0";
          //                Diagnostics
          //                "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)>0.1&&abs(muonEta)<2.4&&centrality<0.8&&muonPt>25.0&&missPt>20.0&&missPt<9000.0&&isolationMu<0.1&&muonMt>40.0&&ZDY==0";
          //                Sanity check
          //                "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt<9000.0&&isolationMu<0.1&&muonMt>40.0&&ZDY==0";
          ///isolation efficiency study

          //              sCutsSig="muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0&&ZDY==0";
          ///Closure test
          //                sCutsSig="muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&isolationMu<0.1&&muonMt>40.0&&ZDY==0";

     ///isolation systematic with cut of 0.2 instead of 0.1
     if(doLoosenIsoCut)     
          sCutsSig = 
            "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&isolationMu<0.2&&muonMt>40.0&&ZDY==0";

      ///mpt systematic up
      // "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>35.85&&missPt<9000.0&&isolationMu<0.1&&muonMt>40.0&&ZDY==0";

      ///mpt systematic down
//      "muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0&&abs(muonEta)<2.4&&muonPt>25.0&&missPt>14.15&&missPt<9000.0&&isolationMu<0.1&&muonMt>40.0&&ZDY==0";

	  std::cout << " Signals cuts: "<< sCutsSig << std::endl;

      ///Cuts for generator-level tau sample
      TString sCutsTau = "";
      sCutsTau = "abs(muonEta)<2.4&&muonPt>25.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0";
      //cross-check w/o mpt
      //sCutsTau = "abs(muonEta)<2.4&&muonPt>25.0&&missPt>0.0&&missPt<9000.0&&muonMt>40.0";
      std::cout << " Tau cuts: "<< sCutsTau << std::endl;

      ///choose to match with chosen correction factor
    RooArgList muonRecArgList(muonEta,muonPt,missPt,muonMt);
    muonRecArgList.add(muonQuality); 
    muonRecArgList.add(centrality);
    muonRecArgList.add(muonELoss);
    muonRecArgList.add(muonScat);
    muonRecArgList.add(ZDY);
    muonRecArgList.add(isolationMu);
    muonRecArgList.add(muonGenRecMatched);

    RooArgList muonTauArgList(muonEta,muonPt,missPt,muonMt);

	RooFormulaVar cutsSig("cutsSig", "cutsSig", sCutsSig, muonRecArgList);
	RooFormulaVar cutsTau("cutsTau", "cutsTau", sCutsTau, muonTauArgList);

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
    muonArgSet.add(muonGenRecMatched);
    muonArgSet.add(w);
	
	RooArgSet centralityArgSet(centrality);  

	// --- Set pt and eta bins ---
	std::vector<double> ptBins;
	ptBins.push_back(0.0);
	//ptBins.push_back(33.5);
	ptBins.push_back(ptmax);
	const int nPtBins = ptBins.size()-1;

	std::vector<double> etaBins;
    if(!doMirrorEta) etaBins.push_back(-2.4);
	else etaBins.push_back(0.10);

    ///do not bin in eta when plotting 
    ///the distribution
	if (doEta&&!doPlotEta&&!doPlotAbsoluteEta) {

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
        etaBins.push_back(1.37);
        etaBins.push_back(1.52);
//        etaBins.push_back(1.73);
        etaBins.push_back(1.74);
        etaBins.push_back(2.1);



	}
//	etaBins.push_back(+2.50);
	etaBins.push_back(+2.40);

	const int nEtaBins = etaBins.size()-1;
    std::cout << "Plotting for " << nEtaBins << " eta bins." << std::endl;
	std::vector<double> centralityBins;

	std::vector <float> ncoll;
    std::vector <double> npartBins;

	centralityBins.push_back(0.00);
	if (doCentrality) {

        //binning used for optimizing Nsig-Nqcd
        //(0-10,10-40,40-80)
//		centralityBins.push_back(0.10);
//		centralityBins.push_back(0.40);
		//ncoll
//		ncoll.push_back(1500.6); //0-10
//		ncoll.push_back(440.6); //10-40
//		ncoll.push_back(77.8); //40-80
		///npart
//		npartBins.push_back(356.2);//0-10
//		npartBins.push_back(192.1);//10-40
//		npartBins.push_back(45.93);//40-80

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
    std::cout << "Plotting for " << nCentralityBins << " centrality bins." << std::endl;

    ///MPT lower thresholds as a function of centrality
    std::vector <double> mptLow;
    double mptLowCut = 25.0;
    if(doMptSigmaDown&&doCentrality){
        ///-1sigma
        mptLow.push_back(mptLowCut); //0-5%
        mptLow.push_back(mptLowCut); //5-10%
        mptLow.push_back(mptLowCut); //10-15%
        mptLow.push_back(mptLowCut); //15-20%
        mptLow.push_back(mptLowCut); //20-40%
        mptLow.push_back(mptLowCut); //40-80%
        //for(int impt=0; impt<mptLow.size(); ++impt) std::cout << mptLow[impt] << std::endl;
    }
    else if (doMptSigmaUp&&doCentrality){
        ///+1 sigma
        mptLow.push_back(mptLowCut); //0-5%
        mptLow.push_back(mptLowCut); //5-10%
        mptLow.push_back(mptLowCut); //10-15%
        mptLow.push_back(mptLowCut); //15-20%
        mptLow.push_back(mptLowCut); //20-40%
        mptLow.push_back(mptLowCut); //40-80%
    }
	/// --- Open output file ---
	TDirectory *dir = gDirectory;
	TFile outFile(fileNameFitOut,"RECREATE");
	gDirectory = dir;
	TString resultName = "WPtFit";
	FitResult baselineResult(resultName,ptBins,etaBins,centralityBins);

	//stacks for sum of binned histos
	THStack* hmData = new THStack("hmData","hmData");
	THStack* hmDataPlus = new THStack("hmDataPlus","hmDataPlus");
	THStack* hmDataMinus = new THStack("hmDataMinus","hmDataMinus");
	THStack* hmQCD = new THStack("hmQCD","hmQCD");
	THStack* hmQCDPlus = new THStack("hmQCDPlus","hmQCDPlus");
	THStack* hmQCDMinus = new THStack("hmQCDMinus","hmQCDMinus");
	THStack* hmZ = new THStack("hmZ","hmZ");
	THStack* hmZPlus = new THStack("hmZPlus","hmZPlus");
	THStack* hmZMinus = new THStack("hmZMinus","hmZMinus");
	THStack* hmTau = new THStack("hmTau","hmTau");
	THStack* hmTauPlus = new THStack("hmTauPlus","hmTauPlus");
	THStack* hmTauMinus = new THStack("hmTauMinus","hmTauMinus");
	THStack* hmMcSum = new THStack("hmMcSum","hmMcSum");
	THStack* hmMcSumPlus = new THStack("hmMcSumPlus","hmMcSumPlus");
	THStack* hmMcSumMinus = new THStack("hmMcSumMinus","hmMcSumMinus");

	//merged W kinematic histos
	std::vector<std::string> kinBins;
	kinBins.push_back("Phi");
	kinBins.push_back("Eta");
	kinBins.push_back("MPT");
	kinBins.push_back("pT");
	kinBins.push_back("isolation");
	const int nKinBins = kinBins.size();
	
	TList _hstackData;
	TList _hstackDataPlus;
	TList _hstackDataMinus;
	TList _hstackQCD;
	TList _hstackQCDPlus;
	TList _hstackQCDMinus;
	TList _hstackZ;
	TList _hstackZPlus;
	TList _hstackZMinus;
	TList _hstackTau;
	TList _hstackTauPlus;
	TList _hstackTauMinus;
	TList _hstackMcSum;
	TList _hstackMcSumPlus;
	TList _hstackMcSumMinus;

	char shData[50], shDataPlus[50],shDataMinus[50];
	char shQCD[50], shQCDPlus[50],shQCDMinus[50];
	char shZ[50], shZPlus[50],shZMinus[50];
	char shTau[50], shTauPlus[50],shTauMinus[50];
	char shMcSum[50], shMcSumPlus[50],shMcSumMinus[50];
	//add THStacks for each kinematic variable to lists
	for(int ik=0; ik<nKinBins; ik++){
		sprintf(shData,"hstackData%i",ik);
		sprintf(shDataPlus,"hstackDataPlus%i",ik);
		sprintf(shDataMinus,"hstackDataMinus%i",ik);
		sprintf(shQCD,"hstackQCD%i",ik);
		sprintf(shQCDPlus,"hstackQCDPlus%i",ik);
		sprintf(shQCDMinus,"hstackQCDMinus%i",ik);
		sprintf(shZ,"hstackZ%i",ik);
		sprintf(shZPlus,"hstackZPlus%i",ik);
		sprintf(shZMinus,"hstackZMinus%i",ik);
		sprintf(shTau,"hstackTau%i",ik);
		sprintf(shTauPlus,"hstackTauPlus%i",ik);
		sprintf(shTauMinus,"hstackTauMinus%i",ik);
		sprintf(shMcSum,"hstackMcSum%i",ik);
		sprintf(shMcSumPlus,"hstackMcSumPlus%i",ik);
		sprintf(shMcSumMinus,"hstackMcSumMinus%i",ik);

		_hstackData.Add( new THStack((TString(kinBins.at(ik))+shData).Data(), (TString(kinBins.at(ik))+shData).Data() ) );
		_hstackDataPlus.Add(new THStack((TString(kinBins.at(ik))+shDataPlus).Data(), (TString(kinBins.at(ik))+shDataPlus).Data() ));
		_hstackDataMinus.Add(new THStack((TString(kinBins.at(ik))+shDataMinus).Data(), (TString(kinBins.at(ik))+shDataMinus).Data()));

		_hstackQCD.Add(new THStack((TString(kinBins.at(ik))+shQCD).Data(), (TString(kinBins.at(ik))+shQCD).Data() ));
		_hstackQCDPlus.Add(new THStack((TString(kinBins.at(ik))+shQCDPlus).Data(), (TString(kinBins.at(ik))+shQCDPlus).Data() ));
		_hstackQCDMinus.Add(new THStack((TString(kinBins.at(ik))+shQCDMinus).Data(), (TString(kinBins.at(ik))+shQCDMinus).Data() ));

		_hstackZ.Add(new THStack((TString(kinBins.at(ik))+shZ).Data(), (TString(kinBins.at(ik))+shZ).Data() ));
		_hstackZPlus.Add(new THStack((TString(kinBins.at(ik))+shZPlus).Data(), (TString(kinBins.at(ik))+shZPlus).Data() ));
		_hstackZMinus.Add(new THStack((TString(kinBins.at(ik))+shZMinus).Data(), (TString(kinBins.at(ik))+shZMinus).Data() ));

		_hstackTau.Add(new THStack((TString(kinBins.at(ik))+shTau).Data(), (TString(kinBins.at(ik))+shTau).Data() ));
		_hstackTauPlus.Add(new THStack((TString(kinBins.at(ik))+shTauPlus).Data(), (TString(kinBins.at(ik))+shTauPlus).Data() ));
		_hstackTauMinus.Add(new THStack((TString(kinBins.at(ik))+shTauMinus).Data(), (TString(kinBins.at(ik))+shTauMinus).Data() ));

		_hstackMcSum.Add(new THStack((TString(kinBins.at(ik))+shMcSum).Data(), (TString(kinBins.at(ik))+shMcSum).Data() ));
		_hstackMcSumPlus.Add(new THStack((TString(kinBins.at(ik))+shMcSumPlus).Data(), (TString(kinBins.at(ik))+shMcSumPlus).Data() ));
		_hstackMcSumMinus.Add(new THStack((TString(kinBins.at(ik))+shMcSumMinus).Data(), (TString(kinBins.at(ik))+shMcSumMinus).Data() ));
	}

  	// --- Fill data sets ---
        //RooDataSet* dataSet = fillHIMuonDataSet(baseString,fileNameDataIn+".root",muonArgSet, cutValue,false,true,false,false); dataSet->Print(); 	
    ///Main ds
    RooDataSet* dataSet = fillHIMuonDataSet(baseString,fileNameDataIn+".root",muonArgSet,cutValue,false); dataSet->Print(); 	
        ///Closure test
//        RooDataSet* dataSet = fillHIMuonDataSet(baseString,fileNameDataIn+".root",muonArgSet,cutValue,true,false,false,false); dataSet->Print(); 	
	//apply W selection cuts
	dataSet = (RooDataSet*)dataSet->reduce(Cut(cutsSig)); dataSet->Print(); 
        std::cout << "Observed events in data: " << dataSet->numEntries() << std::endl;

  	// --- Fill mc sets ---

    ///Weight W MC ds by N-N collision probability
/*    RooDataSet* mcWSet = fillHIMuonDataSet(baseString,fileNameMCWIn+".root",muonArgSet, cutValue, true); mcWSet->Print();
    ///get total number of signal events from all nucleon combinations for scaling
    double nWtot = mcWSet->numEntries();
    std::cout << "Number of total W MC entries before weighting: " << nWtot << std::endl;
*/
    // pp
    // mu+
    RooDataSet* mcWppSetPlus = fillHIMuonDataSet(baseString,fileNameIn_ppPlus+".root",muonArgSet, cutValue, true); mcWppSetPlus->Print();
    mcWppSetPlus = (RooDataSet*)mcWppSetPlus->reduce(Cut(cutsSig)); mcWppSetPlus->Print();
    mcWppSetPlus = (RooDataSet*)mcWppSetPlus->reduce(Cut("motherRec==24")); mcWppSetPlus->Print();
    double nWppPlus = mcWppSetPlus->numEntries();
    std::cout << "Number of pp-->W+ entries before weighting: " << nWppPlus << std::endl;
    // mu-
    RooDataSet* mcWppSetMinus = fillHIMuonDataSet(baseString,fileNameIn_ppMinus+".root",muonArgSet, cutValue, true); mcWppSetMinus->Print();
    mcWppSetMinus = (RooDataSet*)mcWppSetMinus->reduce(Cut(cutsSig)); mcWppSetMinus->Print();
    mcWppSetMinus = (RooDataSet*)mcWppSetMinus->reduce(Cut("motherRec==24")); mcWppSetMinus->Print();
    double nWppMinus = mcWppSetMinus->numEntries();
    std::cout << "Number of pp-->W- entries before weighting: " << nWppMinus << std::endl;

    // np
    // mu+
    RooDataSet* mcWnpSetPlus = fillHIMuonDataSet(baseString,fileNameIn_npPlus+".root",muonArgSet, cutValue, true); mcWnpSetPlus->Print();
    mcWnpSetPlus = (RooDataSet*)mcWnpSetPlus->reduce(Cut(cutsSig)); mcWnpSetPlus->Print();
    mcWnpSetPlus = (RooDataSet*)mcWnpSetPlus->reduce(Cut("motherRec==24")); mcWnpSetPlus->Print();
    double nWnpPlus = mcWnpSetPlus->numEntries();
    std::cout << "Number of np-->W+ entries before weighting: " << nWnpPlus << std::endl;
    // mu-
    RooDataSet* mcWnpSetMinus = fillHIMuonDataSet(baseString,fileNameIn_npMinus+".root",muonArgSet, cutValue, true); mcWnpSetMinus->Print();
    mcWnpSetMinus = (RooDataSet*)mcWnpSetMinus->reduce(Cut(cutsSig)); mcWnpSetMinus->Print();
    mcWnpSetMinus = (RooDataSet*)mcWnpSetMinus->reduce(Cut("motherRec==24")); mcWnpSetMinus->Print();
    double nWnpMinus = mcWnpSetMinus->numEntries();
    std::cout << "Number of np-->W- entries before weighting: " << nWnpMinus << std::endl;

    // pn
    // mu+
    RooDataSet* mcWpnSetPlus = fillHIMuonDataSet(baseString,fileNameIn_pnPlus+".root",muonArgSet, cutValue, true); mcWpnSetPlus->Print();
    mcWpnSetPlus = (RooDataSet*)mcWpnSetPlus->reduce(Cut(cutsSig)); mcWpnSetPlus->Print();
    mcWpnSetPlus = (RooDataSet*)mcWpnSetPlus->reduce(Cut("motherRec==24")); mcWpnSetPlus->Print();
    double nWpnPlus = mcWpnSetPlus->numEntries();
    std::cout << "Number of pn-->W+ entries before weighting: " << nWpnPlus << std::endl;
    // mu-
    RooDataSet* mcWpnSetMinus = fillHIMuonDataSet(baseString,fileNameIn_pnMinus+".root",muonArgSet, cutValue, true); mcWpnSetMinus->Print();
    mcWpnSetMinus = (RooDataSet*)mcWpnSetMinus->reduce(Cut(cutsSig)); mcWpnSetMinus->Print();
    mcWpnSetMinus = (RooDataSet*)mcWpnSetMinus->reduce(Cut("motherRec==24")); mcWpnSetMinus->Print();
    double nWpnMinus = mcWpnSetMinus->numEntries();
    std::cout << "Number of pn-->W- entries before weighting: " << nWpnMinus << std::endl;
    // nn
    // mu+
    RooDataSet* mcWnnSetPlus = fillHIMuonDataSet(baseString,fileNameIn_nnPlus+".root",muonArgSet, cutValue, true); mcWnnSetPlus->Print();
    mcWnnSetPlus = (RooDataSet*)mcWnnSetPlus->reduce(Cut(cutsSig)); mcWnnSetPlus->Print();
    mcWnnSetPlus = (RooDataSet*)mcWnnSetPlus->reduce(Cut("motherRec==24")); mcWnnSetPlus->Print();
    double nWnnPlus = mcWnnSetPlus->numEntries();
    std::cout << "Number of nn-->W+ entries before weighting: " << nWnnPlus << std::endl;
    // mu-
    RooDataSet* mcWnnSetMinus = fillHIMuonDataSet(baseString,fileNameIn_nnMinus+".root",muonArgSet, cutValue, true); mcWnnSetMinus->Print();
    mcWnnSetMinus = (RooDataSet*)mcWnnSetMinus->reduce(Cut(cutsSig)); mcWnnSetMinus->Print();
    mcWnnSetMinus = (RooDataSet*)mcWnnSetMinus->reduce(Cut("motherRec==24")); mcWnnSetMinus->Print();
    double nWnnMinus = mcWnnSetMinus->numEntries();
    std::cout << "Number of nn-->W- entries before weighting: " << nWnnMinus << std::endl;

    double nWtotPlus = nWnnPlus+nWpnPlus+nWnpPlus+nWppPlus;
    std::cout << "Number of total W+ MC entries before weighting: " << nWtotPlus << std::endl;
    double nWtotMinus = nWnnMinus+nWpnMinus+nWnpMinus+nWppMinus;
    std::cout << "Number of total W+ MC entries before weighting: " << nWtotMinus << std::endl;

    // Weight datasets such that 15.5% are from pp, 47.8% from np+pn, and 36.7% from nn
    // pp
    // mu+
    mcWppSetPlus = fillHIMuonDataSet(baseString,fileNameIn_ppPlus+".root",muonArgSet, cutValue, true,0.155*nWtotPlus/nWppPlus); mcWppSetPlus->Print();
    mcWppSetPlus = (RooDataSet*)mcWppSetPlus->reduce(Cut(cutsSig)); mcWppSetPlus->Print();
    mcWppSetPlus = (RooDataSet*)mcWppSetPlus->reduce(Cut("motherRec==24")); mcWppSetPlus->Print();
    mcWppSetPlus = weightDS(mcWppSetPlus,w); mcWppSetPlus->Print();
    std::cout << "Number of pp-->W+ entries after weighting: " << mcWppSetPlus->numEntries() << std::endl;
    // mu-
    mcWppSetMinus = fillHIMuonDataSet(baseString,fileNameIn_ppMinus+".root",muonArgSet, cutValue, true,0.155*nWtotMinus/nWppMinus); mcWppSetMinus->Print();
    mcWppSetMinus = (RooDataSet*)mcWppSetMinus->reduce(Cut(cutsSig)); mcWppSetMinus->Print();
    mcWppSetMinus = (RooDataSet*)mcWppSetMinus->reduce(Cut("motherRec==24")); mcWppSetMinus->Print();
    mcWppSetMinus = weightDS(mcWppSetMinus,w); mcWppSetMinus->Print();
    std::cout << "Number of pp-->W+ entries after weighting: " << mcWppSetMinus->numEntries() << std::endl;

    // np
    // mu+
    mcWnpSetPlus = fillHIMuonDataSet(baseString,fileNameIn_npPlus+".root",muonArgSet, cutValue, true,0.478/2.*nWtotPlus/nWnpPlus); mcWnpSetPlus->Print();
    mcWnpSetPlus = (RooDataSet*)mcWnpSetPlus->reduce(Cut(cutsSig)); mcWnpSetPlus->Print();
    mcWnpSetPlus = (RooDataSet*)mcWnpSetPlus->reduce(Cut("motherRec==24")); mcWnpSetPlus->Print();
    mcWnpSetPlus = weightDS(mcWnpSetPlus,w); mcWnpSetPlus->Print();
    std::cout << "Number of np-->W+ entries after weighting: " << mcWnpSetPlus->numEntries() << std::endl;
    // mu-
    mcWnpSetMinus = fillHIMuonDataSet(baseString,fileNameIn_npMinus+".root",muonArgSet, cutValue, true,0.478/2.*nWtotMinus/nWnpMinus); mcWnpSetMinus->Print();
    mcWnpSetMinus = (RooDataSet*)mcWnpSetMinus->reduce(Cut(cutsSig)); mcWnpSetMinus->Print();
    mcWnpSetMinus = (RooDataSet*)mcWnpSetMinus->reduce(Cut("motherRec==24")); mcWnpSetMinus->Print();
    mcWnpSetMinus = weightDS(mcWnpSetMinus,w); mcWnpSetMinus->Print();
    std::cout << "Number of np-->W+ entries after weighting: " << mcWnpSetMinus->numEntries() << std::endl;

    // pn
    // mu+
    mcWpnSetPlus = fillHIMuonDataSet(baseString,fileNameIn_pnPlus+".root",muonArgSet, cutValue, true,0.478/2.*nWtotPlus/nWpnPlus); mcWpnSetPlus->Print();
    mcWpnSetPlus = (RooDataSet*)mcWpnSetPlus->reduce(Cut(cutsSig)); mcWpnSetPlus->Print();
    mcWpnSetPlus = (RooDataSet*)mcWpnSetPlus->reduce(Cut("motherRec==24")); mcWpnSetPlus->Print();
    mcWpnSetPlus = weightDS(mcWpnSetPlus,w); mcWpnSetPlus->Print();
    std::cout << "Number of pn-->W+ entries after weighting: " << mcWpnSetPlus->numEntries() << std::endl;
    // mu-
    mcWpnSetMinus = fillHIMuonDataSet(baseString,fileNameIn_pnMinus+".root",muonArgSet, cutValue, true,0.478/2.*nWtotMinus/nWpnMinus); mcWpnSetMinus->Print();
    mcWpnSetMinus = (RooDataSet*)mcWpnSetMinus->reduce(Cut(cutsSig)); mcWpnSetMinus->Print();
    mcWpnSetMinus = (RooDataSet*)mcWpnSetMinus->reduce(Cut("motherRec==24")); mcWpnSetMinus->Print();
    mcWpnSetMinus = weightDS(mcWpnSetMinus,w); mcWpnSetMinus->Print();
    std::cout << "Number of pn-->W+ entries after weighting: " << mcWpnSetMinus->numEntries() << std::endl;

    // nn
    // mu+
    mcWnnSetPlus = fillHIMuonDataSet(baseString,fileNameIn_nnPlus+".root",muonArgSet, cutValue, true,0.367*nWtotPlus/nWnnPlus); mcWnnSetPlus->Print();
    mcWnnSetPlus = (RooDataSet*)mcWnnSetPlus->reduce(Cut(cutsSig)); mcWnnSetPlus->Print();
    mcWnnSetPlus = (RooDataSet*)mcWnnSetPlus->reduce(Cut("motherRec==24")); mcWnnSetPlus->Print();
    mcWnnSetPlus = weightDS(mcWnnSetPlus,w); mcWnnSetPlus->Print();
    std::cout << "Number of nn-->W+ entries after weighting: " << mcWnnSetPlus->numEntries() << std::endl;
    // mu-
    mcWnnSetMinus = fillHIMuonDataSet(baseString,fileNameIn_nnMinus+".root",muonArgSet, cutValue, true,0.367*nWtotMinus/nWnnMinus); mcWnnSetMinus->Print();
    mcWnnSetMinus = (RooDataSet*)mcWnnSetMinus->reduce(Cut(cutsSig)); mcWnnSetMinus->Print();
    mcWnnSetMinus = (RooDataSet*)mcWnnSetMinus->reduce(Cut("motherRec==24")); mcWnnSetMinus->Print();
    mcWnnSetMinus = weightDS(mcWnnSetMinus,w); mcWnnSetMinus->Print();
    std::cout << "Number of nn-->W+ entries after weighting: " << mcWnnSetMinus->numEntries() << std::endl;

    ///Add the signal sets 
    RooDataSet* mcWSetPlus = mcWppSetPlus;
    mcWSetPlus->append(*mcWnpSetPlus); mcWSetPlus->append(*mcWpnSetPlus); mcWSetPlus->append(*mcWnnSetPlus);
    mcWSetPlus->Print();

    RooDataSet* mcWSetMinus = mcWppSetMinus;
    mcWSetMinus->append(*mcWnpSetMinus); mcWSetMinus->append(*mcWpnSetMinus); mcWSetMinus->append(*mcWnnSetMinus);
    mcWSetMinus->Print();

    RooDataSet* mcZSet = fillHIMuonDataSet(baseString,fileNameMCZIn+".root",muonArgSet, cutValue, true); mcZSet->Print();
	mcZSet = (RooDataSet*)mcZSet->reduce(Cut(cutsSig)); mcZSet->Print();
  	mcZSet = (RooDataSet*)mcZSet->reduce(Cut("motherRec==23")); mcZSet->Print();
    RooDataSet* mcZEvents = fillHIEventDataSet(baseString,fileNameMCZIn+".root",centralityArgSet ); mcZEvents->Print();

    RooDataSet* mcTauSet = fillHIWTauDataSet(baseString,fileNameMCTauIn+".root",muonArgSet); mcTauSet->Print();
	mcTauSet = (RooDataSet*)mcTauSet->reduce(Cut(cutsTau)); mcTauSet->Print();

    RooDataSet* mcJ1Set = fillHIMuonDataSet(baseString,fileNameMCJ1In+".root",muonArgSet, cutValue, true); mcJ1Set->Print();
	mcJ1Set = (RooDataSet*)mcJ1Set->reduce(Cut(cutsSig)); mcJ1Set->Print();
    RooDataSet* mcJ1Events = fillHIEventDataSet(baseString,fileNameMCJ1In+".root",centralityArgSet );

    RooDataSet* mcJ2Set = fillHIMuonDataSet(baseString,fileNameMCJ2In+".root",muonArgSet, cutValue, true); mcJ2Set->Print();
	mcJ2Set = (RooDataSet*)mcJ2Set->reduce(Cut(cutsSig)); mcJ2Set->Print();
    RooDataSet* mcJ2Events = fillHIEventDataSet(baseString,fileNameMCJ2In+".root",centralityArgSet );

    RooDataSet* mcJ3Set = fillHIMuonDataSet(baseString,fileNameMCJ3In+".root",muonArgSet, cutValue, true); mcJ3Set->Print();
	mcJ3Set = (RooDataSet*)mcJ3Set->reduce(Cut(cutsSig)); mcJ3Set->Print();
    RooDataSet* mcJ3Events = fillHIEventDataSet(baseString,fileNameMCJ3In+".root",centralityArgSet );

    //apply additional MC selection criteria for reco muons
    ///Closure test; only use on dataSet if it is MC    
  	//dataSet = (RooDataSet*)dataSet->reduce(Cut("motherRec==24")); dataSet->Print();

	// --- Subdivide in bins ---
	RooDataSet* dataSubSet[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcWSubSetPlus[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcWSubSetMinus[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcTauSubSet[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcZSubSet[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcJ1SubSet[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcJ2SubSet[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcJ3SubSet[nPtBins][nEtaBins][nCentralityBins];

	///get charge inclusive A0 distr from Cw parametrization
	TGraph* grCwA0 = (TGraph*) fCwParam->Get("gr104a0");
	TGraph* grCwA1 = (TGraph*) fCwParam->Get("gr104a1");
	TGraph* grCwA2 = (TGraph*) fCwParam->Get("gr104a2");

	_cWsystGraphAx.Add((TGraphErrors*)fCwSyst0->Get("grA2Diff"));
	_cWsystGraphAx.Add((TGraphErrors*)fCwSyst1->Get("grA1Diff"));
	_cWsystGraphAx.Add((TGraphErrors*)fCwSyst2->Get("grA0Diff"));

    ///uncomment when using Cw parametrization technique
/*	for ( int ieta = 0; ieta < nEtaBins; ieta++) {

		//get difference in Cw for each charge for this eta bin
		TString sCwChDiff = "grCentDiffeta"; sCwChDiff+=ieta;
		_cWsystGraphs.Add((TGraphErrors*)fCwSyst2->Get(sCwChDiff));
		//get residuals
		TString sCwResid = "grCent104eta"; sCwResid+=ieta; sCwResid+="_Residuals";
		_cWResiduals.Add((TGraph*)fCwParam->Get(sCwResid));
	}
*/
	for ( int i = 0; i < nPtBins; i++ ) {
	  for ( int j = 0; j < nEtaBins; j++ ) {
	    for ( int k = 0; k < nCentralityBins; k++ ){
        	dataSubSet[i][j][k] = selectPtEtaCentrality(dataSet, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1],centralityBins[k], centralityBins[k+1],doMirrorEta); 
        	mcWSubSetPlus[i][j][k] = selectPtEtaCentrality(mcWSetPlus , ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],doMirrorEta); 
        	mcWSubSetMinus[i][j][k] = selectPtEtaCentrality(mcWSetMinus , ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],doMirrorEta); 
        	mcTauSubSet[i][j][k] = selectPtEtaCentrality(mcTauSet, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],doMirrorEta); 
        	mcZSubSet[i][j][k] = selectPtEtaCentrality(mcZSet, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],doMirrorEta); 
            ///widen bins for Jx1mu so no bins are empty
            mcJ1SubSet[i][j][k] = selectPtEtaCentrality(mcJ1Set, ptBins[i], ptBins[i+1], etaBins[0], etaBins[nEtaBins], 
                centralityBins[k], centralityBins[k+1],doMirrorEta); 
            mcJ2SubSet[i][j][k] = selectPtEtaCentrality(mcJ2Set, ptBins[i], ptBins[i+1], etaBins[0], etaBins[nEtaBins], 
                centralityBins[k], centralityBins[k+1],doMirrorEta); 
            mcJ3SubSet[i][j][k] = selectPtEtaCentrality(mcJ3Set, ptBins[i], ptBins[i+1], etaBins[0], etaBins[nEtaBins], 
                centralityBins[k], centralityBins[k+1],doMirrorEta); 
            /*if(doMptSigmaDown||doMptSigmaUp) {

                TString sMptMtCorrelation = "h2DMptMtCorrelationCent"; sMptMtCorrelation+=k; sMptMtCorrelation+="_pfx";
                std::cout << "Opening " << sMptMtCorrelation << " TProfilex histo." << std::endl;
                _pfxMptMtCorrelation = (TProfile*)_fMptMtCorrelation->Get(sMptMtCorrelation);
                if(_pfxMptMtCorrelation!=0) std::cout << sMptMtCorrelation << " opened." << std::endl;
                else exit(0);
                ///Cut on shifted mpt and mt cuts
                dataSubSet[i][j][k] = selectMptMt(dataSubSet[i][j][k],mptLow[k],_pfxMptMtCorrelation); 
                mcWSubSet[i][j][k]  = selectMptMt(mcWSubSet[i][j][k],mptLow[k],_pfxMptMtCorrelation); 
                mcTauSubSet[i][j][k]  = selectMptMt(mcTauSubSet[i][j][k],mptLow[k],_pfxMptMtCorrelation); 
                mcZSubSet[i][j][k]  = selectMptMt(mcZSubSet[i][j][k],mptLow[k],_pfxMptMtCorrelation); 
                mcJ1SubSet[i][j][k] = selectMptMt(mcJ1SubSet[i][j][k],mptLow[k],_pfxMptMtCorrelation); 
                mcJ2SubSet[i][j][k] = selectMptMt(mcJ2SubSet[i][j][k],mptLow[k],_pfxMptMtCorrelation); 
                mcJ3SubSet[i][j][k] = selectMptMt(mcJ3SubSet[i][j][k],mptLow[k],_pfxMptMtCorrelation); 
            }*/
	    }
	  }
	}


	
	//retrieve the TGraphs
    //use only for Cw parametrization
	/*double *yCwA0Eta = grCwA0->GetY();  
	double *yCwA1Eta = grCwA1->GetY();  
	double *yCwA2Eta = grCwA2->GetY();  
    
	//systematics
    //use only for Cw parametrization
	double *yAx0 = ((TGraphErrors*)_cWsystGraphAx.At(0))->GetY();
	double *yAx1 = ((TGraphErrors*)_cWsystGraphAx.At(1))->GetY();
	double *yAx2 = ((TGraphErrors*)_cWsystGraphAx.At(2))->GetY();
    */
	// --- fill arrays for Aw,Cw ---
    ///change bool to true if wanting charge dependent
    ///Cw scale factor
	//if(correctSpectra) setCorrectionFactors(false);

    /////////////////////////////////////////////
    ///RETRIEVE SYSTEMATIC ERRORS FROM TXT FILES
    ////////////////////////////////////////////

    ///for the isolation cut
    setIsoSystErrors("systematics/IsolationCutSystematics.07.22.2013.txt");
    ///for the mpt cut
    setMptSystErrors("systematics/mptCutSystematics.08.04.2013.txt");
    ///for the qcd shape at high pt
    setQCDRaaSystErrors("systematics/qcdRaaSystematics.07.30.2013.txt");
    ///electroweak bkg (Zmumu, Wtaumu)
    setElectroweakSystErrors("systematics/electroweakSystematics.07.30.2013.txt");
    ///for using |eta| Cw instead of real eta
//    setCwMirroredEtaSystErrors("systematics/CwRealEtaSystematics.05.26.2013.txt");
    setCwMirroredEtaSystErrors("systematics/CwRealEtaSystematics.07.14.2013.txt");

    ///bin-by-bin Cw diff in mu+/- for eta and centrality for systematics
    TGraphErrors* grCwEta = (TGraphErrors*)fCwEta->Get("grCwEtaDiff");
    TGraphErrors* grCwCent = (TGraphErrors*)fCwCent->Get("grCentDiff");

	for ( int i = 0; i < nPtBins; i++ ) {
	  for ( int j = 0; j < nEtaBins; j++ ) {

        ///uncomment efficiencies you'd like to use in yield corrections
/*        TString sCwDistroPlus = "grWmunuRecNoQualityWselNoIsoNoZvetoNpartDistroPlusEta"; sCwDistroPlus+=j;
        TGraphErrors* grCwPlus = (TGraphErrors*)fCw->Get(sCwDistroPlus);
        TString sCwDistroMinus = "grWmunuRecNoQualityWselNoIsoNoZvetoNpartDistroMinusEta"; sCwDistroMinus+=j;
        TGraphErrors* grCwMinus = (TGraphErrors*)fCw->Get(sCwDistroMinus);
        TString sCwDistro = "grWmunuRecNoQualityWselNoIsoNoZvetoNpartDistroEta"; sCwDistro+=j;
        TGraphErrors* grCw = (TGraphErrors*)fCw->Get(sCwDistro);
        
*/
/*        TString sCwDistroPlus = "grWmunuRecNoQualityWselNoIsoNpartDistroPlusEta"; sCwDistroPlus+=j;
        TGraphErrors* grCwPlus = (TGraphErrors*)fCw->Get(sCwDistroPlus);
        TString sCwDistroMinus = "grWmunuRecNoQualityWselNoIsoNpartDistroMinusEta"; sCwDistroMinus+=j;
        TGraphErrors* grCwMinus = (TGraphErrors*)fCw->Get(sCwDistroMinus);
        TString sCwDistro = "grWmunuRecNoQualityWselNoIsoNpartDistroEta"; sCwDistro+=j;
        TGraphErrors* grCw = (TGraphErrors*)fCw->Get(sCwDistro);
*/        
/*        TString sCwDistroPlus = "grWmunuRecNoQualityWselNpartDistroPlusEta"; sCwDistroPlus+=j;
        TGraphErrors* grCwPlus = (TGraphErrors*)fCw->Get(sCwDistroPlus);
        TString sCwDistroMinus = "grWmunuRecNoQualityWselNpartDistroMinusEta"; sCwDistroMinus+=j;
        TGraphErrors* grCwMinus = (TGraphErrors*)fCw->Get(sCwDistroMinus);
        TString sCwDistro = "grWmunuRecNoQualityWselNpartDistroEta"; sCwDistro+=j;
        TGraphErrors* grCw = (TGraphErrors*)fCw->Get(sCwDistro);
*/ 
/*        TString sCwDistroPlus = "grWmunuRecHiQualityWselNoIsoNoZvetoNpartDistroPlusEta"; sCwDistroPlus+=j;
        TGraphErrors* grCwPlus = (TGraphErrors*)fCw->Get(sCwDistroPlus);
        TString sCwDistroMinus = "grWmunuRecHiQualityWselNoIsoNoZvetoNpartDistroMinusEta"; sCwDistroMinus+=j;
        TGraphErrors* grCwMinus = (TGraphErrors*)fCw->Get(sCwDistroMinus);
        TString sCwDistro = "grWmunuRecHiQualityWselNoIsoNoZvetoNpartDistroEta"; sCwDistro+=j;
        TGraphErrors* grCw = (TGraphErrors*)fCw->Get(sCwDistro);
*/ 
/*        TString sCwDistroPlus = "grWmunuRecHiQualityWselNoIsoNpartDistroPlusEta"; sCwDistroPlus+=j;
        TGraphErrors* grCwPlus = (TGraphErrors*)fCw->Get(sCwDistroPlus);
        TString sCwDistroMinus = "grWmunuRecHiQualityWselNoIsoNpartDistroMinusEta"; sCwDistroMinus+=j;
        TGraphErrors* grCwMinus = (TGraphErrors*)fCw->Get(sCwDistroMinus);
        TString sCwDistro = "grWmunuRecHiQualityWselNoIsoNpartDistroEta"; sCwDistro+=j;
        TGraphErrors* grCw = (TGraphErrors*)fCw->Get(sCwDistro);
*/
        ///final Cw
        ///if binning over all eta
        ///use mirrored background and
        ///correction factors
	    double etaMed = etaBins[j]+(etaBins[j+1]-etaBins[j])/2.0;
        TString sCwDistroPlus,sCwDistroMinus,sCwDistro;
        if(!doMirrorCw){

            sCwDistroPlus = "grWmunuRecHiQualityWselNpartDistroPlusEta"; sCwDistroPlus+=j;
            sCwDistroMinus = "grWmunuRecHiQualityWselNpartDistroMinusEta"; sCwDistroMinus+=j;
            sCwDistro = "grWmunuRecHiQualityWselNpartDistroEta"; sCwDistro+=j;
        }
        else if(!doMirrorEta){
          ///hop over the crack region
          ///in bin [-0.1,0.1]
          if(etaMed==0.0) {std::cout << "Skipping crack region." << std::endl; continue; }
          if(etaMed<0.0){
              std::cout << "Negative eta: will use absolute eta values for Cw." << std::endl;
              ///number of bins in absolute eta
              const int index = indexNegativeEta(j,nEtaBins);
              sCwDistroPlus = "grWmunuRecHiQualityWselNpartDistroPlusEta"; sCwDistroPlus+=index;
              sCwDistroMinus = "grWmunuRecHiQualityWselNpartDistroMinusEta"; sCwDistroMinus+=index;
              sCwDistro = "grWmunuRecHiQualityWselNpartDistroEta"; sCwDistro+=index;
          }
           else {
              std::cout << "Positive eta: will use absolute eta values for Cw." << std::endl;
              const int index = indexPositiveEta(j,nEtaBins);
              sCwDistroPlus = "grWmunuRecHiQualityWselNpartDistroPlusEta"; sCwDistroPlus+=index;
              sCwDistroMinus = "grWmunuRecHiQualityWselNpartDistroMinusEta"; sCwDistroMinus+=index;
              sCwDistro = "grWmunuRecHiQualityWselNpartDistroEta"; sCwDistro+=index;
           }
        }
        else{
            // Correct with Cw from Py6
            sCwDistroPlus = "grWmunuRecHiQualityWselNpartDistroPlusEta"; sCwDistroPlus+=j;
            sCwDistroMinus = "grWmunuRecHiQualityWselNpartDistroMinusEta"; sCwDistroMinus+=j;
            sCwDistro = "grWmunuRecHiQualityWselNpartDistroEta"; sCwDistro+=j;
            // Correct with Cw from PowPy8
            sCwDistroPlus = "grCwNpartDistroPlusEta"; sCwDistroPlus+=j;
            sCwDistroMinus = "grCwNpartDistroMinusEta"; sCwDistroMinus+=j;
            
            // Correct with Cw*Aw

            /*sCwDistroPlus = "grCwAwPlusNpartDistroEta"; sCwDistroPlus+=j;
            sCwDistroMinus = "grCwAwMinusNpartDistroEta"; sCwDistroMinus+=j;
            sCwDistro = "grCwAwNpartDistroEta"; sCwDistro+=j;
            */
        }
            

        std::cout << "Fetching " << sCwDistroPlus << " " << sCwDistroMinus << " " << sCwDistro << std::endl;
        TGraphErrors* grCwPlus = (TGraphErrors*)fCw->Get(sCwDistroPlus);
        TGraphErrors* grCwMinus = (TGraphErrors*)fCw->Get(sCwDistroMinus);
        //TGraphErrors* grCw = (TGraphErrors*)fCw->Get(sCwDistro);
 

        //Cw difference for mu+/- as a fcn of eta
        //double syst0 = grCwEta->GetY()[j];
        double syst0 = 0.0;

        ///fraction of Zs generated in this eta bin from MC
        ///(used for estimating number of Zs in this eta bin in the data)

        ///negative eta; use absolute estimation
        double totalZinEta = 0.0;
        double totalZPlusinEta = 0.0;
        double totalZMinusinEta = 0.0;
        if(!doMirrorEta){
          ///hop over the crack region
          ///in bin [-0.1,0.1]
          if(etaMed==0.0) continue;
          if(etaMed<0.0){
              ///number of bins in absolute eta
              const int index = indexNegativeEta(j,nEtaBins);
              totalZinEta = grTotZEta->GetY()[index];
              totalZPlusinEta = grTotZEtaPlus->GetY()[index];
              totalZMinusinEta = grTotZEtaMinus->GetY()[index];
          }
           else {
              
              const int index = indexPositiveEta(j,nEtaBins);
              totalZinEta = grTotZEta->GetY()[index];
              totalZPlusinEta = grTotZEtaPlus->GetY()[index];
              totalZMinusinEta = grTotZEtaMinus->GetY()[index];
           }
        }
      
        else { 
             totalZinEta = grTotZEta->GetY()[j];
             totalZPlusinEta = grTotZEtaPlus->GetY()[j];
             totalZMinusinEta = grTotZEtaMinus->GetY()[j];
        
        }

	    for ( int k = 0; k < nCentralityBins; k++ ){

		std::cout << " plotting "<<i<<":"<<j<<":"<<k<<std::endl;

		TString sCentLow = "";
		TString sCentUp = "";
		TString sEtaLow = "";
		TString sEtaUp = "";

 		sCentLow += format(100*centralityBins[k]); //sCentLow.Remove(3);
        sCentUp += format(100*centralityBins[k+1]); //sCentUp.Remove(3);

		sEtaLow += format(etaBins[j]);
		sEtaUp += format(etaBins[j+1]);

        ///Name of TGraph in file fBkgTau
        TString sFracTau = "tauBkgFractionCent"; sFracTau+=k;
        ///graph of Wtau background fraction as a fcn of eta for 
        ///each centrlity bin
        TGraphAsymmErrors* grBkgTau = (TGraphAsymmErrors*)fBkgTau->Get(sFracTau);

		//Get A0,1,2 in this eta bin
        /*
		double CwA0 = yCwA0Eta[j]; 
		double CwA1 = yCwA1Eta[j]; 
		double CwA2 = yCwA2Eta[j];*/

        ///IMPORTANT: Two techniques are currently being studied
        ///in order to obtain the correction factors. This includes
        ///a parameterization technique and a simple bin-by-bin calculation.
        ///Choose one below when running this macro.

		///calculate parametrized Cw in this eta bin
		//double Cw = CwA0+CwA1*npartBins[k]+CwA2*npartBins[k]*npartBins[k];

        ///get bin-by-bin Cw
        ///index charge by 102 = mu+, 103 = mu-, 104 = mu^pm
        //double Cw = getEfficiencyMt(false,104, j, k);
        //grCw->Print();
        //grCwPlus->Print();
        //grCwMinus->Print();
        //double Cw = grCw->GetY()[k];
        double CwPlus = grCwPlus->GetY()[k];
        double CwMinus = grCwMinus->GetY()[k];
        double CwStatErr=0.0,CwStatErrPlus=0.0,CwStatErrMinus=0.0;
        //if(grCw->GetEYhigh()[k]>grCw->GetEYlow()[k])CwStatErr=grCw->GetEYhigh()[k];
        //else CwStatErr=grCw->GetEYlow()[k];
        if(grCwPlus->GetEYhigh()[k]>grCwPlus->GetEYlow()[k])CwStatErrPlus=grCwPlus->GetEYhigh()[k];
        else CwStatErrPlus=grCwPlus->GetEYlow()[k];
        if(grCwMinus->GetEYhigh()[k]>grCwMinus->GetEYlow()[k])CwStatErrMinus=grCwMinus->GetEYhigh()[k];
        else CwStatErrMinus=grCwMinus->GetEYlow()[k];

        ///uncomment when using Cw parametrization
		//Get TGraph of charge diff of Cw as fcn of Npart for this eta bin	
/*		double *yCwSyst = ((TGraphErrors*)_cWsystGraphs.At(j))->GetY();
		//Get TGraph of Cw fit residuals
		double *yCwResid = ((TGraph*)_cWResiduals.At(j))->GetY();
*/

//		double syst0 = yAx0[j]; double syst1 = yAx1[j]; double syst2 = yAx2[j]; 
        //Cw difference for mu+/- as a fcn of centrality
//        double syst1 = grCwCent->GetY()[k];
        double syst1 = 0.0;
		double syst2 = 0.0; 
        //difference in Cw between mu+,mu-
//		double syst3 = yCwSyst[k]; 
        double syst3 = 0.0;
//        double syst4 = yCwResid[k];
        //read in statistical uncertainty of Cw from txt and use as systematic 
        //(only when NOT using Cw parametrization technique; use stat uncert on Cw instead)
        //double syst4 = getCorrectionFactorError(false,104, j, k) ;
        
        double syst4 = 0.0;/*grCw->GetEY()[k];*/
        double systPlus4 = 0.0;/*grCwPlus->GetEY()[k];*/
        double systMinus4 = 0.0;/*grCwMinus->GetEY()[k];*/
        //std::cout << "Correction factor for Cw in mu+- in eta bin : " << j << " and centrality bin: " << k << " = " << Cw <<
        //    "+/-" << syst4 << std::endl;  
        std::cout << "Correction factor for Cw in mu+ in eta bin : " << j << " and centrality bin: " << k << " = " <<
            CwPlus << "+/-" << syst4 << std::endl;  
        std::cout << "Correction factor for Cw in mu- in eta bin : " << j << " and centrality bin: " << k << " = " <<
        CwMinus << "+/-" << syst4 << std::endl;  

        //Cw systematic errors added in quadrature
	    /*double CwStatErr = TMath::Sqrt( TMath::Power(syst0,2)+TMath::Power(syst1,2)+TMath::Power(syst2,2)+TMath::Power(syst3,2)+TMath::Power(syst4,2) );
	    double CwStatErrPlus = TMath::Sqrt(
           TMath::Power(syst0,2)+TMath::Power(syst1,2)+TMath::Power(syst2,2)+TMath::Power(syst3,2)+TMath::Power(systPlus4,2) );
	    double CwStatErrMinus = TMath::Sqrt(
           TMath::Power(syst0,2)+TMath::Power(syst1,2)+TMath::Power(syst2,2)+TMath::Power(syst3,2)+TMath::Power(systMinus4,2) );
		///debug
		std::cout << "Cw systematic errors for bin " << i<<":"<<j<<":"<<k<<" = " << syst0 << " " << syst1 << " " << syst2 << " " << syst3 << " " << syst4 << std::endl;
		std::cout << "Added in quadrature gives: " << CwStatErr << std::endl;
*/
		//get qcd bkg systematic for this centrality and eta bin
		double *yBkgQCDSystResidNcoll = grBkgQCDSystResidNcoll->GetY();
		double *yBkgQCDSystChDiffNcoll = grBkgQCDSystChDiffNcoll->GetY();
		double *yBkgQCDSystChDiffEta = grBkgQCDSystChDiffEta->GetY();


		//calculate qcd systemetic err for this cent and eta bin
		//hack since used 5 ncoll bins for qcd bkg determination instead of 6
		//if(k==0) syst0 = yBkgQCDSystResidNcoll[k]; else syst0 = yBkgQCDSystResidNcoll[k-1];
        syst0 = 0.0;
		/*if(k==0)*/ syst1 = yBkgQCDSystChDiffNcoll[k]; /*else syst1 = yBkgQCDSystChDiffNcoll[k-1];*/
        //syst1 = 0.0;
        ///get median of eta bin
	    //double xEta = etaBins[j]+(etaBins[j+1]-etaBins[j])/2.0;
        ///uncomment when using >10 eta bins
		//if(doSubtractBkg) syst2 = yBkgQCDSystChDiffEta[getQCDBkgBinNumber(xEta)];
		if(doSubtractBkg){ 

            int indexEta;

            if(!doMirrorEta){
              ///hop over the crack region
              ///in bin [-0.1,0.1]
              if(etaMed==0.0) std::cout << "WARNING: Bin in crack region. Expect wrong result!" << std::endl;
              if(etaMed<0.0){
                  ///number of bins in absolute eta
                  indexEta = indexNegativeEta(j,nEtaBins);
              }
               else {
                  indexEta = indexPositiveEta(j,nEtaBins);
               }

                //syst2 = yBkgQCDSystChDiffEta[indexEta];
                syst2=0.0;
            }

            else {
                //syst2 = yBkgQCDSystChDiffEta[j];
                syst2 = 0.0;
            }
        }
		double qcdSystError = TMath::Sqrt( TMath::Power(syst0,2)+TMath::Power(syst1,2)+TMath::Power(syst2,2) );
		///debug
		//std::cout << "QCD systematic errors for bin " << i<<":"<<j<<":"<<k<<" = " << syst0 << " " << syst1 << " " << syst2 << std::endl;
		//std::cout << "Added in quadrature gives: " << qcdSystError << std::endl;

		//get Z bkg systematic for this centrality and eta bin
		double* yBkgZSystNcoll = grBkgZSystNcoll->GetY();
		double* yBkgZSystEta = grBkgZSystEta->GetY();

		//syst0 = yBkgZSystNcoll[k]; 
        syst0 = 0.0;

        if(!doMirrorEta){
          int indexEta;
          ///hop over the crack region
          ///in bin [-0.1,0.1]
          if(etaMed==0.0) std::cout << "WARNING: Bin in crack region. Expect wrong result!" << std::endl;
          if(etaMed<0.0){
              ///number of bins in absolute eta
              indexEta = indexNegativeEta(j,nEtaBins);
          }
           else {
              indexEta = indexPositiveEta(j,nEtaBins);
           }
            //syst1 = yBkgZSystEta[indexEta];
            syst1 = 0.0;
        }

        else {
            //syst1 = yBkgZSystEta[j];
            syst1 = 0.0;
        }
		double zBosSystError = TMath::Sqrt( TMath::Power(syst0,2)+TMath::Power(syst1,2) );;

		std::cout << "Z systematic errors for bin " << i<<":"<<j<<":"<<k<<" = " << syst0 << " " << syst1  << std::endl;
		std::cout << "Added in quadrature gives: " << zBosSystError << std::endl;

		/// --- plot inclusive spectra --
		if(!doCentrality && !doEta && !doCharge){
			if(doPlotMt){

			plot(dataSubSet[i][j][k],mcWSubSetPlus[i][j][k],mcTauSubSet[i][j][k],mcZSubSet[i][j][k],mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k],
				mcJ1Events, mcJ2Events, mcJ3Events,mcZEvents, 
				fQCDBkg, grBkgTau, grBkgZ,totalZinEta, CwPlus,CwStatErr,qcdSystError,zBosSystError, 
                		hmData, hmQCD, hmZ, hmTau, hmMcSum, baselineResult, 
				99, 0, 0, nCentralityBins,nEtaBins, mtmax, ncoll[k],
				ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],
                		muonMt,centrality, 
				"#mu^{#pm},0-80" , "0 #leq |#eta| < 2.5" , spreadSheet,sigEvents,statErr,systErr,
                true,correctSpectra,doSubtractBkg,doPreSelKinematics,doMirrorEta);
			}

			if(doPlotPhi){ 
					plotWCandidateKinematic(dataSubSet[i][j][k],mcWSubSetPlus[i][j][k],mcTauSubSet[i][j][k],mcZSubSet[i][j][k],
                    mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k],
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgTau, grBkgZ,totalZinEta, CwPlus,CwStatErrMinus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackData.At(0),(THStack*)_hstackQCD.At(0),(THStack*)_hstackZ.At(0),(THStack*)_hstackTau.At(0),(THStack*)_hstackMcSum.At(0),  
					99, 0, 0,nCentralityBins,nEtaBins, -1*TMath::Pi(), -1*TMath::Pi(),TMath::Pi() , ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], 
					centralityBins[k], centralityBins[k+1],muonPhi,"#mu^{#pm},0-80" , "0.1 #leq |#eta| < 2.4", 20.0,"#phi_{#mu}", true, correctSpectra,
					doSubtractBkg,doPreSelKinematics);
			} if(doPlotEta){ 
					plotWCandidateKinematic(dataSubSet[i][j][k],mcWSubSetPlus[i][j][k],mcTauSubSet[i][j][k],mcZSubSet[i][j][k],
                    mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k],
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgTau, grBkgZ,totalZinEta, CwPlus,CwStatErrMinus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackData.At(1),(THStack*)_hstackQCD.At(1),(THStack*)_hstackZ.At(1),(THStack*)_hstackTau.At(1),(THStack*)_hstackMcSum.At(1),  
					99, 0, 0,nCentralityBins,nEtaBins,-2.4, -2.4, 2.4 , ncoll[k], ptBins[i], ptBins[i+1], 
                    etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],
                    muonEta,"#mu^{#pm},0-80" , "0 #leq |#eta| < 2.5", 22.0,"#eta_{#mu}", true, correctSpectra,doSubtractBkg,doPreSelKinematics);
			} 
            if(doPlotAbsoluteEta){ 
					plotWCandidateKinematic(dataSubSet[i][j][k],mcWSubSetPlus[i][j][k],mcTauSubSet[i][j][k],mcZSubSet[i][j][k],
                    mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k],
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgTau, grBkgZ,totalZinEta, CwPlus,CwStatErrMinus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackData.At(1),(THStack*)_hstackQCD.At(1),(THStack*)_hstackZ.At(1),(THStack*)_hstackTau.At(1),(THStack*)_hstackMcSum.At(1),  
					99, 0, 0,nCentralityBins,nEtaBins,0.1, 0.1, 2.4 , ncoll[k], ptBins[i], ptBins[i+1], 
                    etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],
                    muonEta,"#mu^{#pm},0-80" , "0 #leq |#eta| < 2.5", 22.0,"#eta_{#mu}", true, correctSpectra,doSubtractBkg,doPreSelKinematics);
			}
            if(doPlotMPT){ 
					plotWCandidateKinematic(dataSubSet[i][j][k],mcWSubSetPlus[i][j][k],mcTauSubSet[i][j][k],mcZSubSet[i][j][k],
                    mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k],
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgTau, grBkgZ,totalZinEta, CwPlus,CwStatErrMinus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackData.At(2),(THStack*)_hstackQCD.At(2),(THStack*)_hstackZ.At(2),(THStack*)_hstackTau.At(2),(THStack*)_hstackMcSum.At(2),  
					99, 0, 0,nCentralityBins,nEtaBins,0.0, 25.0, 120.0 , ncoll[k], ptBins[i], ptBins[i+1], 
                    etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],
                    missPt, "#mu^{#pm},0-80" , "0 #leq |#eta| < 2.5", 24.0,"#slash{p_{T}}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);
			} if(doPlotPt){ 
					plotWCandidateKinematic(dataSubSet[i][j][k],mcWSubSetPlus[i][j][k],mcTauSubSet[i][j][k],mcZSubSet[i][j][k],
                    mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k],
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgTau, grBkgZ,totalZinEta, CwPlus,CwStatErrMinus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackData.At(3),(THStack*)_hstackQCD.At(3),(THStack*)_hstackZ.At(3),(THStack*)_hstackTau.At(3),(THStack*)_hstackMcSum.At(3),  
					99, 0, 0,nCentralityBins,nEtaBins,0.0, 25.0, 100.0 , ncoll[k], ptBins[i], ptBins[i+1], 
                    etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],
                    muonPt, "#mu^{#pm},0-80" , "0 #leq |#eta| < 2.5", 50.0,"p^{#mu}_{T}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);
			} if(doPlotIso){ 
					plotWCandidateKinematic(dataSubSet[i][j][k],mcWSubSetPlus[i][j][k],mcTauSubSet[i][j][k],mcZSubSet[i][j][k],
                    mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k],
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgTau, grBkgZ,totalZinEta, CwPlus,CwStatErrMinus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackData.At(4),(THStack*)_hstackQCD.At(4),(THStack*)_hstackZ.At(4),(THStack*)_hstackTau.At(4),(THStack*)_hstackMcSum.At(4),  
					99, 0, 0,nCentralityBins,nEtaBins,0.0, 0.0, 0.35 , ncoll[k], ptBins[i], ptBins[i+1], 
                    etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],
                    isolationMu, "#mu^{#pm},0-80" , "0 #leq |#eta| < 2.5", 50.0,"#Sigma p^{ID}_{T}/p^{#mu}_{T}[GeV]", true, 
                    correctSpectra,doSubtractBkg,doPreSelKinematics);
			}

		}

		
		if(doCharge){
		   //bin in charge
	           std::cout << "Creating charged datasets." << std::endl;
		   RooDataSet* dataSetPlus  = (RooDataSet*) dataSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   RooDataSet* dataSetMinus = (RooDataSet*) dataSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
		   RooDataSet* mcWTempSetPlus  = mcWSubSetPlus[i][j][k];
		   RooDataSet* mcWTempSetMinus = mcWSubSetMinus[i][j][k];
		   RooDataSet* mcTauSetPlus  = (RooDataSet*) mcTauSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   RooDataSet* mcTauSetMinus = (RooDataSet*) mcTauSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
		   RooDataSet* mcZSetPlus  = (RooDataSet*) mcZSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   RooDataSet* mcZSetMinus = (RooDataSet*) mcZSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
		   RooDataSet* mcJ1SetPlus  = (RooDataSet*) mcJ1SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   RooDataSet* mcJ1SetMinus = (RooDataSet*) mcJ1SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
		   RooDataSet* mcJ2SetPlus  = (RooDataSet*) mcJ2SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   RooDataSet* mcJ2SetMinus = (RooDataSet*) mcJ2SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
		   RooDataSet* mcJ3SetPlus  = (RooDataSet*) mcJ3SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   RooDataSet* mcJ3SetMinus = (RooDataSet*) mcJ3SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));



	    TString sSelPlus = "#mu^{+},";
	    TString sSelMinus = "#mu^{-},";
	    TString sSelPlusEta = sEtaLow;
	    TString sSelMinusEta = sEtaLow;
	    sSelPlus += sCentLow; sSelPlus+="-"; sSelPlus+= sCentUp;
	    sSelMinus += sCentLow; sSelMinus+="-"; sSelMinus+= sCentUp;
	    sSelMinusEta += "#leq"; sSelMinusEta+= "|#eta|"; sSelMinusEta+="<"; sSelMinusEta += sEtaUp;
	    sSelPlusEta += "#leq"; sSelPlusEta+= "|#eta|"; sSelPlusEta+="<"; sSelPlusEta += sEtaUp;

	    //NOTE: index i is important for correctly matching efficiencies with given bin.
	    //CHANGE WITH CAUTION!
	    //99 = mu^{pm}, 100 = mu-(inclusive), 101 = mu+(inclusive), 102 = mu+ (binned), 103 = mu- (binned), 104 = mu^{pm} (binned)
	    if(!doCentrality&&!doEta){ 
		
                std::cout << "plotting for mu+: " << i << ":" << j << ":" << k <<std::endl;
                plot(dataSetPlus,mcWTempSetPlus,mcTauSetPlus, mcZSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus,mcJ1Events,mcJ2Events,mcJ3Events,mcZEvents, 
                fQCDBkg, grBkgTau, grBkgZ,totalZPlusinEta, CwPlus,CwStatErrPlus,qcdSystError,zBosSystError, 
                hmDataPlus, hmQCDPlus, hmZPlus, hmTauPlus, hmMcSumPlus, baselineResult, 
                101, j, k,nCentralityBins,nEtaBins, mtmax,  ncoll[k],ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],
                muonMt,centrality,sSelPlus , sSelPlusEta,spreadSheetMuPlus, 
                sigEventsPlus,statErrPlus,systErrPlus,
                true,correctSpectra,doSubtractBkg,doPreSelKinematics);

                std::cout << "plotting for mu-: " << i << ":" << j << ":" << k <<std::endl;
                plot(dataSetMinus,mcWTempSetMinus,mcTauSetMinus, mcZSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus,mcJ1Events,mcJ2Events,mcJ3Events,mcZEvents, 
                fQCDBkg, grBkgTau, grBkgZ,totalZMinusinEta, CwMinus,CwStatErrMinus,qcdSystError,zBosSystError, 
                hmDataMinus, hmQCDMinus, hmZMinus, hmTauMinus, hmMcSumMinus, baselineResult, 
                100, j, k,nCentralityBins,nEtaBins, mtmax,  ncoll[k],ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],
                muonMt,centrality,sSelMinus , sSelMinusEta,spreadSheetMuMinus, 
                sigEventsMinus,statErrMinus,systErrMinus,
                true,correctSpectra,doSubtractBkg,doPreSelKinematics);

  	}
	else {
	   if(doPlotMt){

				std::cout << "plotting for mu+: " << i << ":" << j << ":" << k <<std::endl;
				plot(dataSetPlus,mcWTempSetPlus,mcTauSetPlus, mcZSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus,mcJ1Events,mcJ2Events,mcJ3Events,mcZEvents, 
				fQCDBkg, grBkgTau,grBkgZ,totalZPlusinEta, CwPlus,CwStatErrPlus,qcdSystError,zBosSystError, 
                		hmDataPlus,hmQCDPlus,hmZPlus,hmTauPlus, hmMcSumPlus,  baselineResult, 
				102, j, k,nCentralityBins,nEtaBins, mtmax,  ncoll[k],ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], 
                		centralityBins[k], centralityBins[k+1],muonMt,centrality, 
				sSelPlus , sSelPlusEta,spreadSheetMuPlus,
                		sigEventsPlus,statErrPlus,systErrPlus,
                		true,correctSpectra,doSubtractBkg,doPreSelKinematics,doMirrorEta);
	 
				std::cout << "plotting for mu-: " << i << ":" << j << ":" << k <<std::endl;
				plot(dataSetMinus,mcWTempSetMinus,mcTauSetMinus, mcZSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus,
				mcJ1Events,mcJ2Events,mcJ3Events,mcZEvents, 
				fQCDBkg, grBkgTau,grBkgZ,totalZMinusinEta, CwMinus,CwStatErrMinus,qcdSystError,zBosSystError, 
                		hmDataMinus,hmQCDMinus,hmZMinus,hmTauMinus, hmMcSumMinus,  baselineResult, 
				103, j, k,nCentralityBins,nEtaBins, mtmax,  ncoll[k],ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], 
                		centralityBins[k], centralityBins[k+1],muonMt,centrality, 
				sSelMinus , sSelMinusEta,spreadSheetMuMinus,
              			sigEventsMinus,statErrMinus,systErrMinus,
                		true,correctSpectra,doSubtractBkg,doPreSelKinematics,doMirrorEta);


			}
			if(doPlotPhi){ 
					plotWCandidateKinematic(dataSetPlus,mcWTempSetPlus,mcTauSetPlus,mcZSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus,
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgTau, grBkgZ,totalZPlusinEta, CwPlus,CwStatErrPlus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackDataPlus.At(0),(THStack*)_hstackQCDPlus.At(0), 
                    (THStack*)_hstackZPlus.At(0),(THStack*)_hstackTauPlus.At(0),(THStack*)_hstackMcSumPlus.At(0),  
					102+i, j, k,nCentralityBins,nEtaBins, -1*TMath::Pi(), -1*TMath::Pi(),TMath::Pi() , 
                    ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], 
					centralityBins[k+1],muonPhi, "#mu^{+},0-80" , "0 #leq |#eta| < 2.4", 20.0,"#phi_{#mu}", true, correctSpectra,doSubtractBkg,doPreSelKinematics);

					plotWCandidateKinematic(dataSetMinus,mcWTempSetMinus,mcTauSetMinus,mcZSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus,
					mcJ1Events, mcJ2Events, mcJ3Events, fQCDBkg, grBkgTau, grBkgZ,totalZMinusinEta, CwMinus,CwStatErrMinus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackDataMinus.At(0),(THStack*)_hstackQCDMinus.At(0), (THStack*)_hstackZMinus.At(0),(THStack*)_hstackTauMinus.At(0),
                    (THStack*)_hstackMcSumMinus.At(0), 103+i, j, k,nCentralityBins,nEtaBins, -1*TMath::Pi(), -1*TMath::Pi(),TMath::Pi() , 
					ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], 
					centralityBins[k+1],muonPhi,"#mu^{-},0-80" , "0 #leq |#eta| < 2.4", 20.0,"#phi_{#mu}", true, correctSpectra,doSubtractBkg,doPreSelKinematics);

			} if(doPlotEta){ 
                    ///Note: Make sure the nBins argument matches the binning for the Jx histos
                    ///If plotting over both positive and negative eta, change limit to -1.0*etaBins[nEtaBins]
                    ///If plotting over |eta|, change limit to 1.0*etaBins[0*nEtaBins]
					plotWCandidateEta(dataSetPlus,mcWTempSetPlus,mcTauSetPlus,mcZSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus,
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgTau, grBkgZ,totalZPlusinEta, CwPlus,CwStatErrPlus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackDataPlus.At(1),(THStack*)_hstackQCDPlus.At(1), (THStack*)_hstackZPlus.At(1),(THStack*)_hstackTauPlus.At(1),
                    (THStack*)_hstackMcSumPlus.At(1),  
					102, j, k,nCentralityBins,nEtaBins,
                    -1.0*etaBins[nEtaBins], -1.0*etaBins[nEtaBins], etaBins[nEtaBins], 
                    ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],muonEta,
					sSelPlus , "-2.4 #leq #eta < 2.4", nEtaBins,"#eta_{#mu}", true,
                    correctSpectra,doSubtractBkg,doPreSelKinematics,doMirrorEta);

					plotWCandidateEta(dataSetMinus,mcWTempSetMinus,mcTauSetMinus,mcZSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus,
					mcJ1Events, mcJ2Events, mcJ3Events,fQCDBkg, grBkgTau, grBkgZ,totalZMinusinEta, CwMinus,CwStatErrMinus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackDataMinus.At(1),(THStack*)_hstackQCDMinus.At(1), (THStack*)_hstackZMinus.At(1),(THStack*)_hstackTauMinus.At(1),
                    (THStack*)_hstackMcSumMinus.At(1), 103, j, k,nCentralityBins,nEtaBins,
                    -1.0*etaBins[nEtaBins], -1.0*etaBins[nEtaBins], etaBins[nEtaBins], 
					ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],muonEta,
					sSelMinus , "-2.4 #leq #eta < 2.4", nEtaBins,"#eta_{#mu}", true,
                    correctSpectra,doSubtractBkg,doPreSelKinematics,doMirrorEta);

			} 
            if(doPlotAbsoluteEta){ 
                    ///Note: Make sure the nBins argument matches the binning for the Jx histos
                    ///If plotting over both positive and negative eta, change limit to -1.0*etaBins[nEtaBins]
                    ///If plotting over |eta|, change limit to 1.0*etaBins[0*nEtaBins]
					plotWCandidateEta(dataSetPlus,mcWTempSetPlus,mcTauSetPlus,mcZSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus,
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgTau, grBkgZ,totalZPlusinEta, CwPlus,CwStatErrPlus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackDataPlus.At(1),(THStack*)_hstackQCDPlus.At(1), (THStack*)_hstackZPlus.At(1),(THStack*)_hstackTauPlus.At(1),
                    (THStack*)_hstackMcSumPlus.At(1),  
					102, j, k,nCentralityBins,nEtaBins,
                    1.0*etaBins[0*nEtaBins], 1.0*etaBins[0*nEtaBins], etaBins[nEtaBins], 
                    ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],muonEta,
					sSelPlus , "0.1 #leq |#eta| < 2.4", nEtaBins,"#eta_{#mu}", true,
                    correctSpectra,doSubtractBkg,doPreSelKinematics,doMirrorEta);

					plotWCandidateEta(dataSetMinus,mcWTempSetMinus,mcTauSetMinus,mcZSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus,
					mcJ1Events, mcJ2Events, mcJ3Events,
                    fQCDBkg, grBkgTau, grBkgZ,totalZMinusinEta, CwMinus,CwStatErrMinus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackDataMinus.At(1),(THStack*)_hstackQCDMinus.At(1), (THStack*)_hstackZMinus.At(1),(THStack*)_hstackTauMinus.At(1),
                    (THStack*)_hstackMcSumMinus.At(1), 
                    103, j, k,nCentralityBins,nEtaBins,
                    1.0*etaBins[0*nEtaBins], 1.0*etaBins[0*nEtaBins], etaBins[nEtaBins], 
					ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],muonEta,
					sSelMinus , "0.1 #leq |#eta| < 2.4", nEtaBins,"#eta_{#mu}", true,
                    correctSpectra,doSubtractBkg,doPreSelKinematics,doMirrorEta);

			}
            if(doPlotMPT){ 
                //nominal function
					plotWCandidateKinematic(dataSetPlus,mcWTempSetPlus,mcTauSetPlus,mcZSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus,
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgTau, grBkgZ,totalZPlusinEta, CwPlus,CwStatErrPlus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackDataPlus.At(2),(THStack*)_hstackQCDPlus.At(2), (THStack*)_hstackZPlus.At(2),(THStack*)_hstackTauPlus.At(2),(THStack*)_hstackMcSumPlus.At(2),  
					102, j, k,nCentralityBins,nEtaBins,0.0, 0.0, 120.0 , ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], 
                    centralityBins[k], centralityBins[k+1],missPt,
					"#mu^{+},0-80" , "0.1 #leq |#eta| < 2.4", 30,"#slash{p_{T}}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);

					plotWCandidateKinematic(dataSetMinus,mcWTempSetMinus,mcTauSetMinus,mcZSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus,
					mcJ1Events, mcJ2Events, mcJ3Events, fQCDBkg, grBkgTau, grBkgZ,totalZMinusinEta, CwMinus,CwStatErrMinus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackDataMinus.At(2),(THStack*)_hstackQCDMinus.At(2), (THStack*)_hstackZMinus.At(2),(THStack*)_hstackTauMinus.At(2),
                    (THStack*)_hstackMcSumMinus.At(2), 103+i, j,
                    k,nCentralityBins,nEtaBins,0.0, 0.0, 120.0 , 
					ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],missPt,
					"#mu^{-},0-80" , "0.1 #leq |#eta| < 2.4", 30,"#slash{p_{T}}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);
                

                
			} if(doPlotPt){ 
					plotWCandidateKinematic(dataSetPlus,mcWTempSetPlus,mcTauSetPlus,mcZSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus,
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgTau, grBkgZ,totalZPlusinEta, CwPlus,CwStatErrPlus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackDataPlus.At(3),(THStack*)_hstackQCDPlus.At(3), (THStack*)_hstackZPlus.At(3),(THStack*)_hstackTauPlus.At(3),(THStack*)_hstackMcSumPlus.At(3),  
					102, j, k,nCentralityBins,nEtaBins,0.0, 25.0, 100.0 , ncoll[k], ptBins[i], ptBins[i+1], 
                    etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],muonPt,
					"#mu^{+},0-80" , "0.1 #leq |#eta| < 2.4", 50.0,"p^{#mu}_{T}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);

					plotWCandidateKinematic(dataSetMinus,mcWTempSetMinus,mcTauSetMinus,mcZSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus,
					mcJ1Events, mcJ2Events, mcJ3Events, fQCDBkg, grBkgTau, grBkgZ,totalZMinusinEta, CwMinus,CwStatErrMinus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackDataMinus.At(3),(THStack*)_hstackQCDMinus.At(3), (THStack*)_hstackZMinus.At(3),(THStack*)_hstackTauMinus.At(3),
                    (THStack*)_hstackMcSumMinus.At(3), 103, j, k,nCentralityBins,nEtaBins,0.0, 25.0, 100.0 , 
					ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],muonPt,
					"#mu^{-},0-80" , "0.1 #leq |#eta| < 2.4", 50.0,"p^{#mu}_{T}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);
			} if(doPlotIso){ 
					plotWCandidateKinematic(dataSetPlus,mcWTempSetPlus,mcTauSetPlus,mcZSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus,
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgTau, grBkgZ,totalZPlusinEta, CwPlus,CwStatErrPlus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackDataPlus.At(4),(THStack*)_hstackQCDPlus.At(4), (THStack*)_hstackZPlus.At(4),(THStack*)_hstackTauPlus.At(4),(THStack*)_hstackMcSumPlus.At(4),  
					102+i, j, k,nCentralityBins,nEtaBins,0.0, 0.0, 0.35 , ncoll[k], ptBins[i], ptBins[i+1], 
                    etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],isolationMu,
					"#mu^{+},0-80" , "0.1 #leq |#eta| < 2.4", 50.0,"#Sigma p^{ID}_{T}/p^{#mu}_{T}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);

					plotWCandidateKinematic(dataSetMinus,mcWTempSetMinus,mcTauSetMinus,mcZSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus,
					mcJ1Events, mcJ2Events, mcJ3Events, fQCDBkg, grBkgTau, grBkgZ,totalZMinusinEta, CwMinus,CwStatErrMinus,qcdSystError,zBosSystError, 
                    (THStack*)_hstackDataMinus.At(4),(THStack*)_hstackQCDMinus.At(4), (THStack*)_hstackZMinus.At(4),(THStack*)_hstackTauMinus.At(4),
                    (THStack*)_hstackMcSumMinus.At(4), 103+i, j, k,nCentralityBins,nEtaBins,0.0, 0.0, 0.35 , 
					ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],isolationMu,
					"#mu^{-},0-80" , "0.1 #leq |#eta| < 2.4", 50.0,"#Sigma p^{ID}_{T}/p^{#mu}_{T}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);
			}
		   }

		}//charge bin

		if (!doCharge&&doCentrality&&doEta) {

			std::cout << "Plotting for centrality and eta classes for mu^{#pm} " << std::endl;
			TString sSel = "#mu^{#pm},"; sSel += sCentLow; sSel+="-"; sSel+= sCentUp;
			TString sSelEta = sEtaLow; sSelEta += " #leq "; sSelEta+= "|#eta| "; sSelEta+="< "; sSelEta += sEtaUp;

			plot(dataSubSet[i][j][k],mcWSubSetPlus[i][j][k],mcTauSubSet[i][j][k],mcZSubSet[i][j][k],mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k], 
				mcJ1Events, mcJ2Events, mcJ3Events,mcZEvents, 
				fQCDBkg, grBkgTau, grBkgZ,totalZinEta, CwPlus,CwStatErr,qcdSystError,zBosSystError, 
                hmData, hmQCD, hmZ, hmTau, hmMcSum,  baselineResult, 
				104+i, j, k, nCentralityBins,nEtaBins, mtmax, ncoll[k],ptBins[i], ptBins[i+1], etaBins[j],etaBins[j+1], centralityBins[k], centralityBins[k+1],
                muonMt,centrality, sSel , sSelEta, spreadSheet, 
                sigEvents,statErr,systErr, 
                true,correctSpectra,doSubtractBkg,doPreSelKinematics);

		}  //no charge bin
		
	   } //centrality
	  } //eta
	}//pt

    //plot merged histograms
	if(doCharge) {
		
		TString sCentLow = "";
		TString sCentUp = "";
		TString sEtaLow = "";
		TString sEtaUp = "";

 		sCentLow += format(100.0*centralityBins[0]); //sCentLow.Remove(3);
		sCentUp += format(100.0*centralityBins[nCentralityBins]); //sCentUp.Remove(3);

		sEtaLow += format(etaBins[0]);
		sEtaUp += format(etaBins[nEtaBins]);
   		    
		//TString sSelPlus = "W^{+}#rightarrow#mu^{+}#nu,";
		//TString sSelMinus = "W^{-}#rightarrow#mu^{-}#bar{#nu},";
		TString sSelPlus = "#mu^{+}, ";
		TString sSelMinus = "#mu^{-}, ";
		TString sSelPlusEta = sEtaLow;
		TString sSelMinusEta = sEtaLow;
		sSelPlus += sCentLow; sSelPlus+="-"; sSelPlus+= sCentUp;
		sSelMinus += sCentLow; sSelMinus+="-"; sSelMinus+= sCentUp;
		sSelMinusEta += "#leq"; sSelMinusEta+= "|#eta|"; sSelMinusEta+="<"; sSelMinusEta += sEtaUp;
		sSelPlusEta += "#leq"; sSelPlusEta+= "|#eta|"; sSelPlusEta+="<"; sSelPlusEta += sEtaUp;

		std::cout << "Plotting merged histos..." << std::endl;
	        //99 = mu^{pm}, 100 = mu-(inclusive), 101 = mu+(inclusive), 102 = mu+ (binned), 103 = mu- (binned), 104 = mu^{pm} (binned)
		if(doPlotMt){
			plotMerged(hmDataPlus,hmQCDPlus,hmZPlus,hmTauPlus,hmMcSumPlus, sSelPlus , sSelPlusEta, 101,
                        correctSpectra,doSubtractBkg,scaleByBinWidth,sigEventsPlus,statErrPlus,systErrPlus);
			plotMerged(hmDataMinus,hmQCDMinus,hmZMinus,hmTauMinus,hmMcSumMinus, sSelMinus , sSelMinusEta, 100,
                        correctSpectra,doSubtractBkg,scaleByBinWidth,sigEventsMinus,statErrMinus,systErrMinus);
		}
		std::cout << "Plotting merged kinematic histos..." << std::endl;
		if(doPlotPhi){	
			plotMergedKinematic((THStack*)_hstackDataPlus.At(0),(THStack*)_hstackQCDPlus.At(0),
                (THStack*)_hstackZPlus.At(0),(THStack*)_hstackTauPlus.At(0),(THStack*)_hstackMcSumPlus.At(0),
			    -1*TMath::Pi(),TMath::Pi(),20.0,-1*TMath::Pi(),"#phi_{#mu}", sSelPlus , sSelPlusEta, 101, correctSpectra,scaleByBinWidth);

			plotMergedKinematic((THStack*)_hstackDataMinus.At(0),(THStack*)_hstackQCDMinus.At(0),
                (THStack*)_hstackZMinus.At(0),(THStack*)_hstackTauMinus.At(0),(THStack*)_hstackMcSumMinus.At(0),
			    -1*TMath::Pi(),TMath::Pi(),20.0,-1*TMath::Pi(),"#phi_{#mu}", sSelMinus , sSelMinusEta, 100, correctSpectra,scaleByBinWidth);
		}
		if(doPlotEta){

            ///If plotting over both positive and negative eta, change limit to -1.0*etaBins[nEtaBins]
            ///If plotting over |eta|, change limit to 1.0*etaBins[0*nEtaBins]
			plotMergedKinematic((THStack*)_hstackDataPlus.At(1),(THStack*)_hstackQCDPlus.At(1),
            (THStack*)_hstackZPlus.At(1),(THStack*)_hstackTauPlus.At(1),(THStack*)_hstackMcSumPlus.At(1),
			-1.0*etaBins[nEtaBins], etaBins[nEtaBins],nEtaBins,-1.0*etaBins[nEtaBins],"#eta_{#mu}",sSelPlus , sSelPlusEta, 101, correctSpectra,scaleByBinWidth);

			plotMergedKinematic((THStack*)_hstackDataMinus.At(1),(THStack*)_hstackQCDMinus.At(1),
            (THStack*)_hstackZMinus.At(1),(THStack*)_hstackTauMinus.At(1),(THStack*)_hstackMcSumMinus.At(1),
			-1.0*etaBins[nEtaBins], etaBins[nEtaBins],nEtaBins,-1.0*etaBins[nEtaBins],"#eta_{#mu}", sSelMinus , sSelMinusEta, 100, correctSpectra,scaleByBinWidth);
		}
		if(doPlotAbsoluteEta){

            ///If plotting over both positive and negative eta, change limit to -1.0*etaBins[nEtaBins]
            ///If plotting over |eta|, change limit to 1.0*etaBins[0*nEtaBins]
			plotMergedKinematic((THStack*)_hstackDataPlus.At(1),(THStack*)_hstackQCDPlus.At(1),
            (THStack*)_hstackZPlus.At(1),(THStack*)_hstackTauPlus.At(1),(THStack*)_hstackMcSumPlus.At(1),
			1.0*etaBins[0*nEtaBins], etaBins[nEtaBins],nEtaBins,1.0*etaBins[0*nEtaBins],"#eta_{#mu}",sSelPlus , sSelPlusEta, 101, correctSpectra,scaleByBinWidth);

			plotMergedKinematic((THStack*)_hstackDataMinus.At(1),(THStack*)_hstackQCDMinus.At(1),
            (THStack*)_hstackZMinus.At(1),(THStack*)_hstackTauMinus.At(1),(THStack*)_hstackMcSumMinus.At(1),
			1.0*etaBins[0*nEtaBins], etaBins[nEtaBins],nEtaBins,1.0*etaBins[0*nEtaBins],"#eta_{#mu}", sSelMinus , sSelMinusEta, 100, correctSpectra,scaleByBinWidth);
		}
		if(doPlotMPT){
			plotMergedKinematic((THStack*)_hstackDataPlus.At(2),(THStack*)_hstackQCDPlus.At(2),
                (THStack*)_hstackZPlus.At(2),(THStack*)_hstackTauPlus.At(2),(THStack*)_hstackMcSumPlus.At(2),
			    25.0, 100.0,25,0.0,"#slash{p_{T}}[GeV]",sSelPlus , sSelPlusEta, 101, correctSpectra,scaleByBinWidth);
            

			plotMergedKinematic((THStack*)_hstackDataMinus.At(2),(THStack*)_hstackQCDMinus.At(2),
                (THStack*)_hstackZMinus.At(2),(THStack*)_hstackTauMinus.At(2),(THStack*)_hstackMcSumMinus.At(2),
			    25.0, 100.0,25, 0.0, "#slash{p_{T}}[GeV]",sSelMinus , sSelMinusEta, 100, correctSpectra,scaleByBinWidth);
           /* 
            // cross-check w/o mpt cut 
			plotMergedKinematic((THStack*)_hstackDataPlus.At(2),(THStack*)_hstackQCDPlus.At(2),
                (THStack*)_hstackZPlus.At(2),(THStack*)_hstackTauPlus.At(2),(THStack*)_hstackMcSumPlus.At(2),
			    0.0, 100.0,25,0.0,"#slash{p_{T}}[GeV]",sSelPlus , sSelPlusEta, 101, correctSpectra,scaleByBinWidth);
			plotMergedKinematic((THStack*)_hstackDataMinus.At(2),(THStack*)_hstackQCDMinus.At(2),
                (THStack*)_hstackZMinus.At(2),(THStack*)_hstackTauMinus.At(2),(THStack*)_hstackMcSumMinus.At(2),
			    0.0, 100.0,25, 0.0, "#slash{p_{T}}[GeV]",sSelMinus , sSelMinusEta, 100, correctSpectra,scaleByBinWidth);
        */
		}	
		if(doPlotPt){
			plotMergedKinematic((THStack*)_hstackDataPlus.At(3),(THStack*)_hstackQCDPlus.At(3),
                (THStack*)_hstackZPlus.At(3),(THStack*)_hstackTauPlus.At(3),(THStack*)_hstackMcSumPlus.At(3),
			    25.0, 100.0 ,50.0,0.0, "p^{#mu}_{T}[GeV]",sSelPlus , sSelPlusEta, 101, correctSpectra,scaleByBinWidth);

			plotMergedKinematic((THStack*)_hstackDataMinus.At(3),(THStack*)_hstackQCDMinus.At(3),
                (THStack*)_hstackZMinus.At(3),(THStack*)_hstackTauMinus.At(3),(THStack*)_hstackMcSumMinus.At(3),
			    25.0, 100.0 ,50.0,0.0, "p^{#mu}_{T}[GeV]",sSelMinus , sSelMinusEta, 100, correctSpectra,scaleByBinWidth);
		}
		if(doPlotIso){
			plotMergedKinematic((THStack*)_hstackDataPlus.At(4),(THStack*)_hstackQCDPlus.At(4),
                (THStack*)_hstackZPlus.At(4),(THStack*)_hstackTauPlus.At(4),(THStack*)_hstackMcSumPlus.At(4),
			    0.0, 0.35 ,50.0,0.0,"#Sigma p^{ID}_{T}/p^{#mu}_{T}[GeV]",sSelPlus , sSelPlusEta, 101, correctSpectra,scaleByBinWidth);

			plotMergedKinematic((THStack*)_hstackDataMinus.At(4),(THStack*)_hstackQCDMinus.At(4),
                (THStack*)_hstackZMinus.At(4),(THStack*)_hstackTauMinus.At(4),(THStack*)_hstackMcSumMinus.At(4),
			    0.0, 0.35 ,50.0,0.0,"#Sigma p^{ID}_{T}/p^{#mu}_{T}[GeV]",sSelMinus , sSelMinusEta, 100, correctSpectra,scaleByBinWidth);
		}

	}
	else {
		
		TString sCentLow = "";
		TString sCentUp = "";
		TString sEtaLow = "";
		TString sEtaUp = "";

 		sCentLow = 100*centralityBins[0]; //sCentLow.Remove(3);
		sCentUp = 100*centralityBins[nCentralityBins]; //sCentUp.Remove(3);

		sEtaLow += etaBins[0];
		sEtaUp += etaBins[nEtaBins];
 
		TString sSel = "#mu^{#pm},"; sSel += sCentLow; sSel+="-"; sSel+= sCentUp;
		TString sSelEta = sEtaLow; sSelEta += " #leq "; sSelEta+= "|#eta| "; sSelEta+="< "; sSelEta += sEtaUp;
		plotMerged(hmData,hmQCD,hmZ,hmTau,hmMcSum,sSel , sSelEta, 99, correctSpectra,doSubtractBkg);

		if(doPlotPhi){	
			plotMergedKinematic((THStack*)_hstackData.At(0),(THStack*)_hstackQCD.At(0),(THStack*)_hstackZ.At(0),(THStack*)_hstackTau.At(0),
            (THStack*)_hstackMcSum.At(0),
            -1*TMath::Pi(),TMath::Pi(),20.0,
			-1*TMath::Pi(),"#phi_{#mu}",sSel , sSelEta, 99, correctSpectra,scaleByBinWidth);
		}
		if(doPlotEta){
			plotMergedKinematic((THStack*)_hstackData.At(1),(THStack*)_hstackQCD.At(1),(THStack*)_hstackZ.At(1),(THStack*)_hstackTau.At(1),
            (THStack*)_hstackMcSum.At(1),-2.4, 2.4 ,22.0,-2.4,"#eta_{#mu}",
			sSel , sSelEta, 99, correctSpectra,scaleByBinWidth);
		}
		if(doPlotMPT){
			plotMergedKinematic((THStack*)_hstackData.At(2),(THStack*)_hstackQCD.At(2),(THStack*)_hstackZ.At(2),(THStack*)_hstackTau.At(2),
            (THStack*)_hstackMcSum.At(2),25.0, 120.0,24.0,0.0,"#slash{p_{T}}[GeV]",
			sSel , sSelEta, 99, correctSpectra,scaleByBinWidth);
		}	
		if(doPlotPt){
			plotMergedKinematic((THStack*)_hstackData.At(3),(THStack*)_hstackQCD.At(3),(THStack*)_hstackZ.At(3),(THStack*)_hstackTau.At(3),
            (THStack*)_hstackMcSum.At(3),25.0, 100.0 ,50.0,0.0,"p^{#mu}_{T}[GeV]",
			sSel , sSelEta, 99, correctSpectra,scaleByBinWidth);
		}
		if(doPlotIso){
			plotMergedKinematic((THStack*)_hstackData.At(4),(THStack*)_hstackQCD.At(4),(THStack*)_hstackZ.At(4),(THStack*)_hstackTau.At(4),
            (THStack*)_hstackMcSum.At(4),0.0, 0.35 ,50.0,0.0,
			"#Sigma p^{ID}_{T}/p^{#mu}_{T}[GeV]",sSel , sSelEta, 99, correctSpectra,scaleByBinWidth);
		}


	}

	outFile.cd();
	baselineResult.write(outFile);
	outFile.Close();
	spreadSheet.close();
	spreadSheetMuPlus.close();
	spreadSheetMuMinus.close();
	std::cout << "Clean Up." << std::endl;
	//stacks for sum of binned histos
	delete hmData ;
	delete hmDataPlus ;
	delete hmDataMinus ;
	delete hmQCD ;
	delete hmQCDPlus ;
	delete hmQCDMinus ;
	delete hmZ  ;
	delete hmZPlus ;
	delete hmZMinus ;
	delete hmTau  ;
	delete hmTauPlus ;
	delete hmTauMinus ;
	delete hmMcSum ;
	delete hmMcSumPlus ;
	delete hmMcSumMinus ;

	for(int ik=0; ik<nKinBins; ik++){
		delete _hstackData.At(ik);
		delete _hstackDataPlus.At(ik);
		delete _hstackDataMinus.At(ik);
		delete _hstackQCD.At(ik);
		delete _hstackQCDPlus.At(ik);
		delete _hstackQCDMinus.At(ik);
		delete _hstackZ.At(ik);
		delete _hstackZPlus.At(ik);
		delete _hstackZMinus.At(ik);
		delete _hstackTau.At(ik);
		delete _hstackTauPlus.At(ik);
		delete _hstackTauMinus.At(ik);
		delete _hstackMcSum.At(ik);
		delete _hstackMcSumPlus.At(ik);
		delete _hstackMcSumMinus.At(ik);

	}

    delete fQCDBkg;
    delete fQCDBkgSystNcoll;
    delete fQCDBkgSystEta;
    delete fZBkg;
    delete fBkgTau;
    delete fBkgZSystNcoll;
    delete fBkgZEta;
    delete fCw;
    delete fCwEta;
    delete fCwCent;
    delete fCwParam;
    delete fCwSyst0;
    delete fCwSyst1;
    delete fCwSyst2;

    return 0;
}//WAnalysis


int main(){
	int result1 = WAnalysis();
    return 0;
}
