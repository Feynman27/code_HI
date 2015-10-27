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
#include "WPlotterHelper.C"

///////////////////////////////
//plot merged histograms
//////////////////////////////
void plotMerged(THStack* hmData,THStack* hmQCD, THStack* hmZ,THStack* hmMcSum, TString sSel, TString sSel2, int chargeIndex, bool correctSpectra){

	double mtcutLow = 40.0; double mtmax = 300.0;
	double nBins = 75.0;

	TCanvas* cdatamt = new TCanvas("cdatamt","cdatamt",600,600);
	//return summed histos in the THStack
	TH1F* hmDatac =(TH1F*)hmData->GetStack()->Last()->Clone("hmDatac");
	TH1F* hmQCDc =(TH1F*)hmQCD->GetStack()->Last()->Clone("hmQCDc");
	TH1F* hmZc =(TH1F*)hmZ->GetStack()->Last()->Clone("hmZc");
	TH1F* hmMcSumc =(TH1F*)hmMcSum->GetStack()->Last()->Clone("hmMcSumc");

	int binLo = nBins/mtmax*mtcutLow+1;

  	double sigEventsUncorr = hmDatac->Integral(binLo,nBins); //integrate from 40-200 GeV before Aw,Cw correction
  	std::cout << "integrated events before correction in mT:"<< mtcutLow << "-" << mtmax  << " =  " << sigEventsUncorr << /*" +-" << TMath::Sqrt(sigEventsUncorr) <<*/ std::endl;

	double effMt ;
	/*correctSpectra = false; //temp hack
	if(correctSpectra) {
		std::cout << "Correcting for acceptance..." << std::endl;
		hmDatac = correctAcceptanceMt(hmDatac, chargeIndex);
		hmQCDc = correctAcceptanceMt(hmQCDc, chargeIndex,false);
		hmZc = correctAcceptanceMt(hmZc, chargeIndex,false);
		hmMcSumc = correctAcceptanceMt(hmMcSumc, chargeIndex,false);
		effMt = getAcceptanceMt(chargeIndex);
	}*/
	//if(correctSpectra) effMt = getAcceptanceMt(chargeIndex,false,false,false,false);

	///systematics fiducial correction
	///booleans = ptUp,ptDown,MPTup,MPTdown
	if(correctSpectra) effMt = getAcceptanceMt(chargeIndex,false,false,false,false);
	else effMt = 1.0;
	std::cout << "Efficiency correction factor = " << effMt << std::endl;
 
  	double sigEvents = hmDatac->Integral(binLo,nBins); //integrate from 40-200 GeV

  	std::cout << "Corrected number of W candidates in mT BEFORE kinematic and geometrical fiducial correction and background subtraction:"<< mtcutLow << "-" << mtmax  << " =  " << sigEvents 
		<< " +-" << TMath::Sqrt(sigEvents)/effMt << std::endl;

	std::cout << "Number of Events from QCD = "  << hmQCDc->Integral(binLo,nBins) << std::endl;
	std::cout << "Number of Events from Z bosons = "  << hmZc->Integral(binLo,nBins) - hmQCDc->Integral(binLo,nBins) << std::endl;
        //histo of Z+QCD
	double bkgEvents = hmZc->Integral(binLo,nBins);
  	std::cout << "Corrected number of W candidates in mT AFTER background subtraction:"<< mtcutLow << "-" << mtmax  << " =  " << sigEvents-bkgEvents 
		<< " +-" << TMath::Sqrt(sigEvents-bkgEvents) << std::endl;
  	std::cout << "Corrected number of W candidates in mT AFTER fiducial correction:"<< mtcutLow << "-" << mtmax  << " =  " << (sigEvents-bkgEvents)/effMt 
		<< " +-" << TMath::Sqrt(sigEvents-bkgEvents)/effMt << std::endl;


	hmQCDc->SetFillColor(kAzure-9);
	hmZc->SetFillColor(kRed);


	hmMcSumc->GetXaxis()->SetTitle("m_{T}[GeV]"); 
	hmMcSumc->GetYaxis()->SetTitle("Muons/4.0 GeV"); 
	//hmMcSumc->GetYaxis()->SetRangeUser(0.1,hmMcSumc->GetMaximum()+2.5e3); 
	hmMcSumc->GetXaxis()->SetRangeUser(mtcutLow,200.0); 
	//hmData->GetStack()->Last()->Draw("PE");
	hmMcSumc->Draw("hist f");
	hmZc->Draw("hist fsame");
	hmQCDc->Draw("hist fsame");
	hmDatac->Draw("pesame");
	hmMcSumc->Draw("sameaxis");

	ATLAS_LABEL(0.17,0.89,1);
	myText(0.33,0.89, (Color_t)kBlack, (char*)("Internal"));

	TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->AddEntry(hmDatac, "Data 2011", "pe");
	leg->AddEntry(hmMcSumc, "W#rightarrow#mu#nu", "f");
	leg->AddEntry(hmZc, "Z#rightarrow#mu#mu", "f");
	leg->AddEntry(hmQCDc, "QCD", "f");
	leg->Draw();

	TLatex l;
	l.SetNDC();
	l.DrawLatex(0.18,0.83,sSel +  "%" );
	l.DrawLatex(0.169,0.767,sSel2);
	l.SetTextSize(0.034);
	l.DrawLatex(0.74,0.89,"#sqrt{s_{NN}}=2.76 TeV");
	l.DrawLatex(0.492,0.89,"#int Ldt #approx 0.140 nb^{-1}"); 
	cdatamt->Update();
	
	TString plotNameLog = "dataMt_"; plotNameLog+=sSel2; plotNameLog+="_"; plotNameLog+=sSel; plotNameLog+="_"; plotNameLog+="_Merged"; plotNameLog+="Log"; 
	TString plotNameLin = "dataMt_"; plotNameLin+=sSel2; plotNameLin+="_"; plotNameLin+=sSel; plotNameLin+="_"; plotNameLin+="_Merged"; plotNameLin+="Lin"; 

        hmMcSumc->GetYaxis()->SetRangeUser(0.1,1600.0); cdatamt->Update();
	cdatamt->Print(plotNameLin.ReplaceAll("#","")+".png");
	cdatamt->Print(plotNameLin.ReplaceAll("#","")+".eps"); 
	cdatamt->Print(plotNameLin.ReplaceAll("#","")+".pdf"); 
	TString plotNameLinRoot = plotNameLin.ReplaceAll("|",",") + ".root";
	cdatamt->Print(plotNameLin); 

        cdatamt->SetLogy(true); hmMcSumc->GetYaxis()->SetRangeUser(0.1, 3.1e4); cdatamt->Update();
	//hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()*1e3); cdatamt->Update();
	cdatamt->Print(plotNameLog.ReplaceAll("#","")+".png");
	cdatamt->Print(plotNameLog.ReplaceAll("#","")+".eps"); 
	cdatamt->Print(plotNameLog.ReplaceAll("#","")+".pdf"); 
	TString plotNameLogRoot = plotNameLog.ReplaceAll("|",",") + ".root";
	cdatamt->Print(plotNameLogRoot); 
	

}
///////////////////////////////
//plot
//////////////////////////////
void plot(RooDataSet* dataSet,RooDataSet* mcWSet, RooDataSet* mcZSet, RooDataSet* mcJ1Set, RooDataSet* mcJ2Set, RooDataSet* mcJ3Set, 
		RooDataSet* mcJ1Events,RooDataSet* mcJ2Events, RooDataSet* mcJ3Events,TFile* fQCDBkg, TGraph2DErrors* grBkgZ, double Cw, double cWSystErr, double qcdSystError, double zBosSystError, 
		THStack* hmData, THStack* hmQCD,THStack* hmZ, THStack* hmMcSum, FitResult& fitResult, const int iMt, const int iEta, const int iCentrality, int nCentralityBins, 
		const float mtmax, double ncoll, double ptLow, double ptUpp, double etaLow, double etaUpp, 
		double centralityLow, double centralityUpp, RooRealVar& muonMt, RooRealVar& centrality, TString sSel, TString sSel2 , bool addPercent = true, 
		bool correctSpectra = true, bool doSubtractBkg = true, bool doPreSelKinematics=false ) {


	int nBins = 75;
	std::cout << "centrality: " << centralityLow << "-" << centralityUpp << " ncoll: " << ncoll << std::endl;

  	RooBinning b = RooBinning(nBins,0.0,mtmax); // 4 GeV per bin

	if(doPreSelKinematics) { std::cout << "Plotting with preselection cuts" << std::endl; correctSpectra=false; }

	if(correctSpectra) std::cout << "Plotting corrected mT. " << std::endl;
	//initialize histograms	
  	// --- data ---
	TH1F* hdataSet = (TH1F*)dataSet->createHistogram("hdataSet",muonMt,Binning(b));
  	// --- W set ---
	TH1F* hmcWSet = (TH1F*)mcWSet->createHistogram("hmcWSet",muonMt,Binning(b));
  	// --- Z set ---
	TH1F* hmcZSet = (TH1F*)mcZSet->createHistogram("hmcZSet",muonMt,Binning(b));
  	// --- QCD set ---
	TH1F* hmcJ1Set = (TH1F*)mcJ1Set->createHistogram("hmcJ1Set",muonMt,Binning(b));
	TH1F* hmcJ2Set = (TH1F*)mcJ2Set->createHistogram("hmcJ2Set",muonMt,Binning(b));
	TH1F* hmcJ3Set = (TH1F*)mcJ3Set->createHistogram("hmcJ3Set",muonMt,Binning(b));

	//return correctly weighted QCD histogram
	TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",nBins,0.0,mtmax);
	hmcQCDSet = getWeightedJxHisto(centralityLow, centralityUpp, ncoll, hmcJ1Set, hmcJ2Set, hmcJ3Set, mcJ1Events, mcJ2Events,mcJ3Events, nBins,0.0, mtmax);

	double mtcutLow = 40.0;
	int binLo; 
	if(doPreSelKinematics) binLo = 1;
	else binLo = nBins/mtmax*mtcutLow+1;

	TH1F* hdataSetc = (TH1F*)hdataSet->Clone("hdataSetc");

  	double sigEventsUncorr = hdataSetc->Integral(binLo,nBins); //integrate from 40-200 GeV before Aw,Cw correction
  	std::cout << "integrated events before correction in mT:"<< mtcutLow << "-" << mtmax  << " =  " << sigEventsUncorr << " +-" << TMath::Sqrt(sigEventsUncorr) << std::endl;
	double effMt ; double effMtErr; double effMtSystErr;
	//double xEta = etaLow+(etaUpp-etaLow)/2.0;
	if(correctSpectra) {
		//hdataSetc = correctEfficiencyMt(hdataSet, iMt, iCentrality, iEta);
		///use charge inclusive (i.e. iMt = 104)
		hdataSetc = correctEfficiencyMt(hdataSet, 104, iCentrality, iEta,Cw);
		//effMt = getEfficiencyMt(iMt, iCentrality, iEta);
		//effMt = getEfficiencyFitMt(iMt, iCentrality, xEta);
		effMt = Cw;
		//effMtErr = getCorrectionFactorError(iMt, iCentrality, iEta);
		effMtErr = 0.0; //hack
		//Cw syst errs added in quadrature
	        effMtSystErr = cWSystErr;
	}
	else {effMt = 1.0; effMtErr = 0.0; effMtSystErr = 0.0;}

  	double sigEvents = hdataSetc->Integral(binLo,nBins); //integrate from 40-200 GeV
	//relative errors of uncorrected signal events and correction factor
	double err1 = TMath::Sqrt(sigEventsUncorr)/sigEventsUncorr*100.0; double err2 = effMtErr/effMt*100.0;
	//propagated errors
	double errStat = TMath::Sqrt( TMath::Power(err1,2)+TMath::Power(err2,2)  ) *0.01*sigEvents; 
	std::cout << "Efficiency correction factor = " << effMt << "+-" << effMtErr << "(stat.) " << effMtSystErr << "(syst.)" << std::endl;
 
  	//std::cout << "integrated events mT:"<< mtcutLow << "-" << mtmax  << " =  " << sigEvents << " +-" << errStat <<"(stat.) " << effMtSystErr/effMt*sigEvents << "(syst.)" << std::endl;
  	std::cout << "integrated events mT:"<< mtcutLow << "-" << mtmax  << " =  " << sigEvents << " +-" << TMath::Sqrt(sigEventsUncorr)/effMt << "(stat.) " 
		<< effMtSystErr/effMt*sigEvents << "(syst.)" << std::endl;

	//running sum of data histos
	hmData->Add(hdataSetc);
	
	//graph of QCD fraction for given charge,eta,and centrality class  
	///use inclusive charge set (i.e. iMt=104)
        TString sFracQCD = "fractionQCD"; sFracQCD+="_charge"; sFracQCD+=104; sFracQCD+="_eta"; sFracQCD+=0; sFracQCD+="_cent"; 
	//hack for using 0-10% instead of 0-5%, 5-10%
	if(iCentrality==0) sFracQCD+=iCentrality;
	else sFracQCD+=iCentrality-1;
	TGraphErrors* grBkgQCD = (TGraphErrors*) fQCDBkg->Get(sFracQCD);

	TH1F* hmcQCDSetc = (TH1F*)hmcQCDSet->Clone("hmcQCDSetc");
	
	//running sum of qcd histos
	hmQCD->Add(hmcQCDSetc);

	double survivalFractionQCD ; 
	//percent of QCD in W signal region b4 W selection
	///use inclusive charge set (i.e. iMt=104)
	std::cout << "Getting QCD background fraction..." << std::endl;
	if(doPreSelKinematics) survivalFractionQCD = getQCDBkg(104,iEta,iCentrality,grBkgQCD,centralityLow, centralityUpp, doPreSelKinematics) ; 	
	else survivalFractionQCD = getQCDBkg(104, iEta, iCentrality,grBkgQCD,centralityLow, centralityUpp) ;

	double dataQCDEvents = survivalFractionQCD*sigEvents;
	//relative error of N_corr(sigEvents) and QCD bkg fraction(survivalFractionQCD)
	double errNcorr = effMtSystErr/effMt; double errQCDFrac = qcdSystError/survivalFractionQCD;
	//propagated absolute error of dataQCDEvents
	double errQCDCounts = TMath::Sqrt( TMath::Power(errNcorr,2)+ TMath::Power(errQCDFrac,2))*dataQCDEvents;
	///debug
	std::cout << "Relative error in corrected counts: " << errNcorr << " Relative error in qcd fraction: " << errQCDFrac << std::endl;
	std::cout << "That gives an absolute error in the number of QCD bkg counts of " << errQCDCounts << std::endl;
	//error of qcd subtracted counts
	double errCountsQCDSub = TMath::Sqrt(TMath::Power(errNcorr*sigEvents,2) + TMath::Power(errQCDCounts,2) );
	std::cout << "This gives an error in the qcd subtracted yield of " << errCountsQCDSub << std::endl;

  	double mcQCDEvents = hmcQCDSetc->Integral(binLo,nBins);
	std::cout << "QCD integral " <<  mcQCDEvents << std::endl;
	double sfQCD = dataQCDEvents/mcQCDEvents;
	if(mcQCDEvents==0) {
		std::cout << "WARNING: 0 QCD MC events in signal region." << std::endl;
		sfQCD = 1.0;
	}
	hmcQCDSetc->Scale(sfQCD);
	hmcQCDSetc->SetFillColor(kAzure-9);
	std::cout << "QCD events in signal region = " << dataQCDEvents << std::endl;
  	std::cout << "background from QCD mT:"<< mtcutLow << "-" << mtmax  << " =  " << dataQCDEvents/sigEvents*100.0 << "%" <<std::endl;

        //Z boson
	std::cout << "Z integral " << hmcZSet->Integral() << std::endl;
	TH1F* hmcZSetc = (TH1F*)hmcZSet->Clone("hmcZSetc");

	double survivalProb; 
	int index = iEta*nCentralityBins+iCentrality;
	if(doPreSelKinematics) survivalProb = getZBkg(iMt,iEta,iCentrality,index,grBkgZ,centralityLow, centralityUpp,doPreSelKinematics) ;
	else survivalProb = getZBkg(iMt,iEta,iCentrality,index, grBkgZ, centralityLow, centralityUpp) ; //percent of Z in W signal region
	std::cout <<"Z survival probability:"<< survivalProb  << std::endl;
	double bkgZEventFrac = survivalProb*getZEvents(centralityLow, centralityUpp, etaLow, etaUpp)/sigEventsUncorr;
	std::cout << bkgZEventFrac << std::endl;
	double dataZEvents = bkgZEventFrac*sigEvents ;

	//absolute err in number of surviving Z muons
	double errZCounts = TMath::Sqrt( TMath::Power(errNcorr,2)+ TMath::Power(zBosSystError/bkgZEventFrac,2))*dataZEvents;

  	double mcZEvents = hmcZSetc->Integral(binLo,nBins);
	std::cout << dataZEvents << std::endl;
	std::cout << mcZEvents << std::endl;
	//weight Z mc by per event yield from paper
	double sfZ = dataZEvents/mcZEvents;
	if(mcZEvents==0) {
		std::cout << "WARNING: 0 Z MC events in signal region." << std::endl;
		sfZ = 1.0;
	}
	hmcZSetc->Scale(sfZ);
	hmcZSetc->Add(hmcQCDSetc);
	hmcZSetc->SetFillColor(kRed);

	//running sum of Z histos
	hmZ->Add(hmcZSetc);

  	std::cout << "background from Z mT:"<< mtcutLow << "-" << mtmax  << " =  " << dataZEvents/sigEvents*100.0 << "%" << std::endl;

	TH1F* hmcWSetc = (TH1F*)hmcWSet->Clone("hmcWSetc");

        //normalize to data
  	double mcWEvents =  hmcWSetc->Integral(binLo,nBins);
	double sigEventsSub = sigEvents-dataQCDEvents-dataZEvents;
	//absolute error in subtracted W yield
	double sigEventsSubErr = TMath::Sqrt( TMath::Power(errCountsQCDSub,2)+TMath::Power(errZCounts,2)); 

	double err3 = TMath::Sqrt(dataQCDEvents); double err4 = TMath::Sqrt(dataZEvents);
	double errStatSub = TMath::Sqrt( TMath::Power(errStat,2)+TMath::Power(err3,2)+TMath::Power(err4,2)  ) ; 
	double sfW = sigEventsSub/mcWEvents;
	hmcWSetc->Scale(sfW);
	//add bkg contribution
	hmcWSetc->Add(hmcZSetc);

	//running sum of combined MC histos	
	hmMcSum->Add(hmcWSetc);	

	mcWEvents = hmcWSetc->Integral(binLo,nBins);
  	std::cout << "integrated events from W MC mT:"<< mtcutLow << "-" << mtmax  << " =  " << mcWEvents << std::endl;

  	std::cout << "integrated signal after bkg subtraction in mT:"<< mtcutLow << "-" << mtmax  << " =  " << sigEventsSub << "+-" << errStatSub << "(stat.) " 
		<< sigEventsSubErr << "(syst.)" << std::endl;
	
	//store sig counts+statistical errors in TGraph
	if(doSubtractBkg) {
		fitResult.setSig(iMt, iEta, iCentrality, sigEventsSub, errStatSub, errStatSub);
		//store sig counts+systematic errors in TGraph
		fitResult.setSigSyst(iMt, iEta, iCentrality, sigEventsSub, sigEventsSubErr, sigEventsSubErr);
	}
	//else if(doSubtractBkg) fitResult.setSig(iMt, iEta, iCentrality, sigEventsSub, errStatSub, errStatSub);
	//else if(correctSpectra) fitResult.setSig(iMt, iEta, iCentrality, sigEvents, errStat, errStat);
	else fitResult.setSig(iMt, iEta, iCentrality, sigEvents, errStat, errStat);

	//now that the histos have been filled and scaled,
	//we now plot them on a canvas
	
  	TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);

	leg->AddEntry(hdataSetc, "Data 2011", "pe");
	leg->AddEntry(hmcWSetc, "W#rightarrow#mu#nu", "f");
	leg->AddEntry(hmcZSetc, "Z#rightarrow#mu#mu", "f");
	leg->AddEntry(hmcQCDSetc, "QCD", "f");

  	TCanvas* cdatamt = new TCanvas("cdatamt","cdatamt",600,600);
  	hmcWSetc->GetXaxis()->SetTitle("m_{T}[GeV]"); 
  	hmcWSetc->GetYaxis()->SetTitle("Muons/4.0 GeV"); 
	//hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()*1.0e2); 
	//hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()+2.0e3); 
	hmcWSetc->GetXaxis()->SetRangeUser(mtcutLow,200.0); 
	hmcWSetc->Draw("hist f");
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
	l.DrawLatex(0.492,0.89,"#int Ldt #approx 0.140 nb^{-1}"); 

	leg->Draw(); cdatamt->Update();
	ATLAS_LABEL(0.17,0.89,1);
	myText(0.33,0.89, (Color_t)kBlack, (char*)("Internal"));

	//TString plotNameLog = "dataMt_"; plotNameLog+=sSel2; plotNameLog+="_"; plotNameLog+=sSel; plotNameLog+="_"; plotNameLog+="Log"; if(doPreSelKinematics) plotNameLog+="_PreSel";
	TString plotNameLog = "dataMt_"; plotNameLog+="charge"; plotNameLog+=iMt; plotNameLog+="_eta"; plotNameLog+=iEta; plotNameLog+="_cent"; plotNameLog+=iCentrality;  
	if(doPreSelKinematics) plotNameLog+="_PreSel";
	//TString plotNameLin = "dataMt_"; plotNameLin+=sSel2; plotNameLin+="_"; plotNameLin+=sSel; plotNameLin+="_"; plotNameLin+="Lin"; if(doPreSelKinematics) plotNameLin+="_PreSel";
	TString plotNameLin = "dataMt_"; plotNameLin+="charge"; plotNameLin+=iMt; plotNameLin+="_eta"; plotNameLin+=iEta; plotNameLin+="_cent"; plotNameLin+=iCentrality;  

  //uncomment for saving intermediate plots
        hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()+75.0); cdatamt->Update();
	cdatamt->Print(plotNameLin.ReplaceAll("#","")+".png");
	cdatamt->Print(plotNameLin.ReplaceAll("#","")+".eps"); 
	cdatamt->Print(plotNameLin.ReplaceAll("#","")+".pdf"); 
	cdatamt->Print(plotNameLin+".root"); 

//  	cdatamt->SetLogy(true); hdataSetc->GetYaxis()->SetRangeUser(0.1,2.1e5); cdatamt->Update();
        cdatamt->SetLogy(true); hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()*7.0e2); cdatamt->Update();
	cdatamt->Print(plotNameLog.ReplaceAll("#","")+".png");
	cdatamt->Print(plotNameLog.ReplaceAll("#","")+".eps"); 
	cdatamt->Print(plotNameLog.ReplaceAll("#","")+".pdf"); 
	TString plotNameLogRoot = plotNameLog.ReplaceAll("|",",") + ".root";
	cdatamt->Print(plotNameLogRoot); 
	

	std::cout << "Clean up" << std::endl;

	delete gDirectory->FindObject("hdataSet");
	delete gDirectory->FindObject("hmcWSet");
	delete gDirectory->FindObject("hmcZSet");
	delete gDirectory->FindObject("hmcJ1Set");
	delete gDirectory->FindObject("hmcJ2Set");
	delete gDirectory->FindObject("hmcJ3Set");
	delete gDirectory->FindObject("hmcQCDSet");
	delete cdatamt;	

	std::cout << "Plotting ends here." << std::endl;
}//plot

///////////////////////////////
//WAnalysis
//////////////////////////////

void WAnalysis(){
	
	bool doCharge = false ;
	bool doCentrality = true ;
	bool doEta = true; 
	bool correctSpectra = true; //corrects spectra using Aw,Cw factors
	bool doSubtractBkg = true; 
	bool doPreSelKinematics = false; //note: do not correctSpectra when plotting with pre-sel cuts
	int cutValue = 11;

	//switch for plotting W candidate kinematics
	bool doPlotMt = true; 
	bool doPlotPhi = false; 
	bool doPlotEta = false;
	bool doPlotMPT = false;
	bool doPlotPt =  false;
	bool doPlotIso = false;

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
	TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
	//TString baseString = "/tmp/tbalestr/";

	//data
	TString fileNameDataIn = "HardProbesFiles/HISingleMuonHP.12.19.2012";
	//data overlay
	TString fileNameMCWIn = "MonteCarloFiles/Wmunu/HISingleMuonWmunuMCDataOverlay.2012.11.21";
	//HIJING overlay
	TString fileNameMCZIn = "MonteCarloFiles/Zmumu/HISingleMuonMCZmumu.12.30.2012";
	//TString fileNameMCZIn = baseString+"MonteCarloFiles/Zmumu/HISingleMuonMC_PYTHIA_HIJING_Zmumu_11.28.2012";

	//J1 1 muon-filter 
	TString fileNameMCJ1In = "MonteCarloFiles/QCD/HISingleMuonMcJ1.2012.11.29";
	//J2 1 muon-filter 
	TString fileNameMCJ2In = "MonteCarloFiles/QCD/HISingleMuonMcJ2.2012.11.29";
	//J3 1 muon-filter 
	TString fileNameMCJ3In = "MonteCarloFiles/QCD/HISingleMuonMcJ3.2012.11.29";

	//files that hold background values
//	TString fileNameBkgQCDIn = "fractionQCDEtaScaled.12.30.2012";
	///use this file for inclusive charge set
	TString fileNameBkgQCDIn = "background/fractionQCDEtaScaled.01.21.2013";
	///files for pt,mpt cut systematics
	///Change the Z bkg input file when running systematics
//	TString fileNameBkgQCDIn = "systematics/QCD/fractionQCDEtaScaled_1sigmaUpPt.02.08.2013"; //pt +1
//	TString fileNameBkgQCDIn = "systematics/QCD/fractionQCDEtaScaled_1sigmaDownPt.02.08.2013"; //pt -1
//	TString fileNameBkgQCDIn = "systematics/QCD/fractionQCDEtaScaled_1sigmaUpMPt.02.08.2013"; //mpt +1
//	TString fileNameBkgQCDIn = "systematics/QCD/fractionQCDEtaScaled_1sigmaDownMPt.02.08.2013"; //mpt -1

    TFile *fQCDBkg = new TFile(fileNameBkgQCDIn+".root","READ");
	if ( !fQCDBkg->IsOpen() ) {
	    std::cout << fQCDBkg << " not found!" << std::endl;
	    exit(0);
	}
	//qcd bkg syst as fcn of ncoll
	TString fileNameBkgQCDSystNcoll = "background/fractionQCDCent.01.31.2013";
        TFile *fQCDBkgSystNcoll = new TFile(fileNameBkgQCDSystNcoll+".root","READ");
	if ( !fQCDBkgSystNcoll->IsOpen() ) {
	    std::cout << fQCDBkgSystNcoll << " not found!" << std::endl;
	    exit(0);
	}
	//qcd bkg syst as fcn of eta
	TString fileNameBkgQCDSystEta = "background/fractionQCDEta.02.04.2013";
        TFile *fQCDBkgSystEta = new TFile(fileNameBkgQCDSystEta+".root","READ");
	if ( !fQCDBkgSystEta->IsOpen() ) {
	    std::cout << fQCDBkgSystEta << " not found!" << std::endl;
	    exit(0);
	}

	//TString fileNameBkgZIn = "fractionZEtaCent.12.30.2012";
	///use this file for inclusive charge set
	TString fileNameBkgZIn = "background/fractionZEtaCent.02.11.2013";
	///files for pt,mpt cut systematics
//	TString fileNameBkgZIn = "systematics/Zboson/ptSigmaUp/fractionZEtaCent_ptSigmaUp.02.10.2013"; //pt +1
//	TString fileNameBkgZIn = "systematics/Zboson/ptSigmaDown/fractionZEtaCent_ptSigmaDown.02.10.2013"; //pt -1
//	TString fileNameBkgZIn = "systematics/Zboson/mptSigmaUp/fractionZEtaCent_mptSigmaUp.02.10.2013"; //mpt +1
//	TString fileNameBkgZIn = "systematics/Zboson/mptSigmaDown/fractionZEtaCent_mptSigmaDown.02.10.2013"; //mpt -1

	TString sFracZ = "fractionZEtaCent";
	TString sFracZPlus = "fractionZEtaCentPlus";
	TString sFracZMinus = "fractionZEtaCentMinus";
	TFile *fZBkg = new TFile(fileNameBkgZIn+".root","READ");
	if ( !fZBkg->IsOpen() ) {
	    std::cout << fZBkg << " not found!" << std::endl;
	    exit(0);
	}

	//files with systmatics for Z
	TString fileNameBkgZSystNcoll = "background/fractionZCent.02.01.2013";
	TString fileNameBkgZSystEta = "background/fractionZEta.02.01.2013";

	TFile *fBkgZSystNcoll = new TFile(fileNameBkgZSystNcoll+".root","READ");
	TFile *fBkgZSystEta = new TFile(fileNameBkgZSystEta+".root","READ");

	//open file that holds TGraph of Cw parametrization values
//	TString fileNameCw = "correctionFactorsW.02.11.2013";
	TString fileNameCw = "correctionFactorsW.02.25.2013";
	//Cw parametrization for +-1Sigma systematics
//	TString fileNameCw = "systematics/eff/ptSigmaUp/correctionFactorsW_pt1SigmaUp"; //+1pt
//	TString fileNameCw = "systematics/eff/ptSigmaDown/correctionFactorsW_pt1SigmaDown"; //-1pt
//	TString fileNameCw = "systematics/eff/mptSigmaUp/correctionFactorsW_mpt1SigmaUp"; //+1mpt
//	TString fileNameCw = "systematics/eff/mptSigmaDown/correctionFactorsW_mpt1SigmaDown"; //-1mpt

	//open files that hold TGraphs of Cw systematic errors
	TString fileNameCwSyst0 = "systematics/eff/correctionFactorsW.AllFloat_syst.01.19.2013"; 
	TString fileNameCwSyst1 = "systematics/eff/correctionFactorsW.a2_fixed_syst.01.19.2013"; 
	TString fileNameCwSyst2 = "systematics/eff/correctionFactorsW.a1_a2_fixed_syst.01.19.2013"; 

	TFile* fCw = new TFile(fileNameCw+".root","READ");
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
	TGraphErrors* grBkgZSystEta = (TGraphErrors*) fBkgZSystEta->Get("diffFracZEta");

	TString fileNameFitOut = "WAnalysis_fitResult"; fileNameFitOut+=nSigPdfString; fileNameFitOut+=".root";

        // --- declare cut variables --- //
	  RooRealVar  muonPt("muonPt","p_{T}",0.0,350.0,"GeV");
	  RooRealVar  missPt("missPt","p_{T}^{miss}",0.0,250.0,"GeV");
	  RooRealVar  muonMt("muonMt","m_{T}",0.0,350.0,"GeV");
	  RooRealVar  muonCharge("muonCharge","charge",-2.0,+2.0);
	  RooRealVar  isolation20("isolation20","isolation20",0.0,10.0);
  	  RooRealVar  centrality("centrality","centrality",0.,1.0);
  	  RooRealVar  muonEta("muonEta","muonEta",-3.0,+3.0);
  	  RooRealVar  muonPhi("muonPhi","muonPhi",-3.5,+3.5);
  	  RooRealVar  ZDY("ZDY","ZDY",0.0,1.0);

	  //W selection cuts
	  //TString sCutsSig = "muonPt>25.0&&missPt>25.0&&missPt<9000.0&&isolation20<0.3&&muonMt>40.0&&ZDY==0";
	  //for consistency check only; above sCutsSig is used for analysis
	  TString sCutsSig = "muonPt>25.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0&&muonMt<400.0&&ZDY==0";
	  std::cout << " Signals cuts: "<< sCutsSig << std::endl;
	  //TString sCutsSig = "muonPt>25.0&&missPt>(25.0-0.2)&&isolation20<0.4&&muonMt>40.0&&ZDY==0";
	  RooFormulaVar cutsSig("cutsSig", "cutsSig", sCutsSig, RooArgList(muonPt,missPt,isolation20,muonMt,ZDY));

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

  	  RooArgSet muonArgSet(muonPt,missPt,isolation20,muonMt,muonEta,centrality,ZDY,muonCategory,chargeCategory);
	  muonArgSet.add(muonPhi);
	
	  RooArgSet centralityArgSet(centrality);  

	// --- Set pt and eta bins ---
	std::vector<double> ptBins;
	ptBins.push_back(0.0);
	//ptBins.push_back(33.5);
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

	std::vector <float> ncoll;
        std::vector <double> npartBins;

	centralityBins.push_back(0.00);
	if (doCentrality) {

		centralityBins.push_back(0.10);
		centralityBins.push_back(0.40);
		//ncoll
		ncoll.push_back(1500.6); //0-10
		ncoll.push_back(440.6); //10-40
		ncoll.push_back(77.8); //40-80*/
		///npart
		npartBins.push_back(356.2);//0-10
		npartBins.push_back(192.1);//10-40
		npartBins.push_back(45.93);//40-80


	}
	else  {
		ncoll.push_back(452.0);//0-80
		npartBins.push_back(139.5);
	}
	centralityBins.push_back(0.80);

	const int nCentralityBins = centralityBins.size()-1;

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
	TList _hstackMcSum;
	TList _hstackMcSumPlus;
	TList _hstackMcSumMinus;

	char shData[50], shDataPlus[50],shDataMinus[50];
	char shQCD[50], shQCDPlus[50],shQCDMinus[50];
	char shZ[50], shZPlus[50],shZMinus[50];
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

		_hstackMcSum.Add(new THStack((TString(kinBins.at(ik))+shMcSum).Data(), (TString(kinBins.at(ik))+shMcSum).Data() ));
		_hstackMcSumPlus.Add(new THStack((TString(kinBins.at(ik))+shMcSumPlus).Data(), (TString(kinBins.at(ik))+shMcSumPlus).Data() ));
		_hstackMcSumMinus.Add(new THStack((TString(kinBins.at(ik))+shMcSumMinus).Data(), (TString(kinBins.at(ik))+shMcSumMinus).Data() ));
	}

  	// --- Fill data sets ---
        RooDataSet* dataSet = fillHIMuonDataSet(baseString,fileNameDataIn+".root",muonArgSet, cutValue); dataSet->Print(); 	
	//apply W selection cuts
	dataSet = (RooDataSet*)dataSet->reduce(Cut(cutsSig)); dataSet->Print(); std::cout << dataSet->numEntries() << std::endl;

  	// --- Fill mc sets ---
        RooDataSet* mcWSet = fillHIMuonDataSet(baseString,fileNameMCWIn+".root",muonArgSet, cutValue, true); mcWSet->Print();
	mcWSet = (RooDataSet*)mcWSet->reduce(Cut(cutsSig)); mcWSet->Print();

        RooDataSet* mcZSet = fillHIMuonDataSet(baseString,fileNameMCZIn+".root",muonArgSet, cutValue, true); mcZSet->Print();
	mcZSet = (RooDataSet*)mcZSet->reduce(Cut(cutsSig)); mcZSet->Print();

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
  	mcWSet = (RooDataSet*)mcWSet->reduce(Cut("muonCategory==muonCategory::W")); mcWSet->Print();
  	mcZSet = (RooDataSet*)mcZSet->reduce(Cut("muonCategory==muonCategory::Z")); mcZSet->Print();

	// --- Subdivide in bins ---
	RooDataSet* dataSubSet[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcWSubSet[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcZSubSet[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcJ1SubSet[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcJ2SubSet[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* mcJ3SubSet[nPtBins][nEtaBins][nCentralityBins];

	///get charge inclusive A0 distr from Cw parametrization
	TGraph* grCwA0 = (TGraph*) fCw->Get("gr104a0");
	TGraph* grCwA1 = (TGraph*) fCw->Get("gr104a1");
	TGraph* grCwA2 = (TGraph*) fCw->Get("gr104a2");

	_cWsystGraphAx.Add((TGraphErrors*)fCwSyst0->Get("grA2Diff"));
	_cWsystGraphAx.Add((TGraphErrors*)fCwSyst1->Get("grA1Diff"));
	_cWsystGraphAx.Add((TGraphErrors*)fCwSyst2->Get("grA0Diff"));

	for ( int ieta = 0; ieta < nEtaBins; ieta++) {

		//get difference in Cw for each charge for this eta bin
		TString sCwChDiff = "grCentDiffeta"; sCwChDiff+=ieta;
		_cWsystGraphs.Add((TGraphErrors*)fCwSyst2->Get(sCwChDiff));
		//get residuals
		TString sCwResid = "grCent104eta"; sCwResid+=ieta; sCwResid+="_Residuals";
		_cWResiduals.Add((TGraph*)fCw->Get(sCwResid));
	}

	for ( int i = 0; i < nPtBins; i++ ) {
	  for ( int j = 0; j < nEtaBins; j++ ) {
	    for ( int k = 0; k < nCentralityBins; k++ ){
        	dataSubSet[i][j][k] = selectPtEtaCentrality(dataSet, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 
        	mcWSubSet[i][j][k] = selectPtEtaCentrality(mcWSet , ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 
        	mcZSubSet[i][j][k] = selectPtEtaCentrality(mcZSet, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 
		//use |eta|<2.5 for Jx since mT shape has no eta dep 
		mcJ1SubSet[i][j][k] = selectPtEtaCentrality(mcJ1Set, ptBins[i], ptBins[i+1], etaBins[0], etaBins[nEtaBins-1], centralityBins[k], centralityBins[k+1],true); 
		mcJ2SubSet[i][j][k] = selectPtEtaCentrality(mcJ2Set, ptBins[i], ptBins[i+1], etaBins[0], etaBins[nEtaBins-1], centralityBins[k], centralityBins[k+1],true); 
		mcJ3SubSet[i][j][k] = selectPtEtaCentrality(mcJ3Set, ptBins[i], ptBins[i+1], etaBins[0], etaBins[nEtaBins-1], centralityBins[k], centralityBins[k+1],true); 
	    }
	  }
	}


	
	//retrieve the TGraphs
	double *yCwA0Eta = grCwA0->GetY();  
	double *yCwA1Eta = grCwA1->GetY();  
	double *yCwA2Eta = grCwA2->GetY();  
	//systematics
	double *yAx0 = ((TGraphErrors*)_cWsystGraphAx.At(0))->GetY();
	double *yAx1 = ((TGraphErrors*)_cWsystGraphAx.At(1))->GetY();
	double *yAx2 = ((TGraphErrors*)_cWsystGraphAx.At(2))->GetY();
	// --- fill arrays for Aw,Cw ---
	if(correctSpectra) setCorrectionFactors();
	for ( int i = 0; i < nPtBins; i++ ) {
	  for ( int j = 0; j < nEtaBins; j++ ) {
	    for ( int k = 0; k < nCentralityBins; k++ ){

		std::cout << " plotting "<<i<<":"<<j<<":"<<k<<std::endl;
		TString sCentLow = "";
		TString sCentUp = "";
		TString sEtaLow = "";
		TString sEtaUp = "";

		sCentLow += 100*centralityBins[k]; //sCentLow.Remove(3);
		sCentUp += 100*centralityBins[k+1]; //sCentUp.Remove(3);

		sEtaLow += etaBins[j];
		sEtaUp += etaBins[j+1];

		//Get A0,1,2 in this eta bin
		double CwA0 = yCwA0Eta[j]; 
		double CwA1 = yCwA1Eta[j]; 
		double CwA2 = yCwA2Eta[j]; 
		///calculate parametrized Cw in this eta bin
		double Cw = CwA0+CwA1*npartBins[k]+CwA2*npartBins[k]*npartBins[k];
		//Get TGraph of charge diff of Cw as fcn of Npart for this eta bin	
		double *yCwSyst = ((TGraphErrors*)_cWsystGraphs.At(j))->GetY();
		//Get TGraph of Cw fit residuals
		double *yCwResid = ((TGraph*)_cWResiduals.At(j))->GetY();

		double syst0 = yAx0[j]; double syst1 = yAx1[j]; double syst2 = yAx2[j]; 
		double syst3 = yCwSyst[k]; double syst4 = yCwResid[k];


	        double cWSystErr = TMath::Sqrt( TMath::Power(syst0,2)+TMath::Power(syst1,2)+TMath::Power(syst2,2)+TMath::Power(syst3,2)+TMath::Power(syst4,2) );
		///debug
		std::cout << "Cw systematic errors for bin " << i<<":"<<j<<":"<<k<<" = " << syst0 << " " << syst1 << " " << syst2 << " " << syst3 << " " << syst4 << std::endl;
		std::cout << "Added in quadrature gives: " << cWSystErr << std::endl;
		//get qcd bkg systematic for this centrality and eta bin
		double *yBkgQCDSystResidNcoll = grBkgQCDSystResidNcoll->GetY();
		double *yBkgQCDSystChDiffNcoll = grBkgQCDSystChDiffNcoll->GetY();
		double *yBkgQCDSystChDiffEta = grBkgQCDSystChDiffEta->GetY();


		//calculate qcd systemetic err for this cent and eta bin
		//hack since used 5 ncoll bins for qcd bkg determination instead of 6
		if(k==1) syst0 = yBkgQCDSystResidNcoll[0]; else syst0 = yBkgQCDSystResidNcoll[k-1];
		if(k==1) syst1 = yBkgQCDSystChDiffNcoll[0]; else syst1 = yBkgQCDSystChDiffNcoll[k-1];
		syst2 = yBkgQCDSystChDiffEta[j];
		double qcdSystError = TMath::Sqrt( TMath::Power(syst0,2)+TMath::Power(syst1,2)+TMath::Power(syst2,2) );
		///debug
		std::cout << "QCD systematic errors for bin " << i<<":"<<j<<":"<<k<<" = " << syst0 << " " << syst1 << " " << syst2 << std::endl;
		std::cout << "Added in quadrature gives: " << qcdSystError << std::endl;

		//get Z bkg systematic for this centrality and eta bin
		double* yBkgZSystNcoll = grBkgZSystNcoll->GetY();
		double* yBkgZSystEta = grBkgZSystEta->GetY();

		syst0 = yBkgZSystNcoll[k]; syst1 = yBkgZSystEta[j];
		double zBosSystError = TMath::Sqrt( TMath::Power(syst0,2)+TMath::Power(syst1,2) );;

		std::cout << "Z systematic errors for bin " << i<<":"<<j<<":"<<k<<" = " << syst0 << " " << syst1  << std::endl;
		std::cout << "Added in quadrature gives: " << zBosSystError << std::endl;

		/// --- plot inclusive spectra --
		if(!doCentrality && !doEta && !doCharge){
			if(doPlotMt){
			plot(dataSubSet[i][j][k],mcWSubSet[i][j][k],mcZSubSet[i][j][k],mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k],
				mcJ1Events, mcJ2Events, mcJ3Events, 
				fQCDBkg, grBkgZ, Cw,cWSystErr,qcdSystError,zBosSystError, hmData,hmQCD,hmZ,hmMcSum, baselineResult, 
				99, 0, 0, nCentralityBins, mtmax, ncoll[k],ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],muonMt,centrality, 
				"#mu^{#pm},0-80" , "0 #leq |#eta| < 2.5" , true,correctSpectra,doSubtractBkg,doPreSelKinematics);
			}

			if(doPlotPhi){ 
					plotWCandidateKinematic(dataSubSet[i][j][k],mcWSubSet[i][j][k],mcZSubSet[i][j][k],mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k],
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgZ, (THStack*)_hstackData.At(0),(THStack*)_hstackQCD.At(0), (THStack*)_hstackZ.At(0),(THStack*)_hstackMcSum.At(0),  
					99, 0, 0,nCentralityBins, -1*TMath::Pi(), -1*TMath::Pi(),TMath::Pi() , Cw,ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], 
					centralityBins[k], centralityBins[k+1],muonPhi,"#mu^{#pm},0-80" , "0 #leq |#eta| < 2.5", 20.0,"#phi_{#mu}", true, correctSpectra,
					doSubtractBkg,doPreSelKinematics);
			} if(doPlotEta){ 
					plotWCandidateKinematic(dataSubSet[i][j][k],mcWSubSet[i][j][k],mcZSubSet[i][j][k],mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k],
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgZ, (THStack*)_hstackData.At(1),(THStack*)_hstackQCD.At(1), (THStack*)_hstackZ.At(1),(THStack*)_hstackMcSum.At(1),  
					99, 0, 0,nCentralityBins,-2.5, -2.5, 2.5 , Cw,ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],muonEta,
					"#mu^{#pm},0-80" , "0 #leq |#eta| < 2.5", 22.0,"#eta_{#mu}", true, correctSpectra,doSubtractBkg,doPreSelKinematics);
			} if(doPlotMPT){ 
					plotWCandidateKinematic(dataSubSet[i][j][k],mcWSubSet[i][j][k],mcZSubSet[i][j][k],mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k],
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgZ, (THStack*)_hstackData.At(2),(THStack*)_hstackQCD.At(2), (THStack*)_hstackZ.At(2),(THStack*)_hstackMcSum.At(2),  
					99, 0, 0,nCentralityBins,0.0, 25.0, 120.0 , Cw,ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],missPt,
					"#mu^{#pm},0-80" , "0 #leq |#eta| < 2.5", 24.0,"#slash{p_{T}}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);
			} if(doPlotPt){ 
					plotWCandidateKinematic(dataSubSet[i][j][k],mcWSubSet[i][j][k],mcZSubSet[i][j][k],mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k],
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgZ, (THStack*)_hstackData.At(3),(THStack*)_hstackQCD.At(3), (THStack*)_hstackZ.At(3),(THStack*)_hstackMcSum.At(3),  
					99, 0, 0,nCentralityBins,0.0, 25.0, 100.0 , Cw,ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],muonPt,
					"#mu^{#pm},0-80" , "0 #leq |#eta| < 2.5", 50.0,"p^{#mu}_{T}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);
			} if(doPlotIso){ 
					plotWCandidateKinematic(dataSubSet[i][j][k],mcWSubSet[i][j][k],mcZSubSet[i][j][k],mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k],
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgZ, (THStack*)_hstackData.At(4),(THStack*)_hstackQCD.At(4), (THStack*)_hstackZ.At(4),(THStack*)_hstackMcSum.At(4),  
					99, 0, 0,nCentralityBins,0.0, 0.0, 0.35 , Cw,ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],isolation20,
					"#mu^{#pm},0-80" , "0 #leq |#eta| < 2.5", 50.0,"#Sigma p^{ID}_{T}/p^{#mu}_{T}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);
			}

		}

		
		if(doCharge){
		   //bin in charge
	           std::cout << "Creating charged datasets." << std::endl;
		   RooDataSet* dataSetPlus  = (RooDataSet*) dataSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   RooDataSet* dataSetMinus = (RooDataSet*) dataSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
		   RooDataSet* mcWSetPlus  = (RooDataSet*) mcWSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   RooDataSet* mcWSetMinus = (RooDataSet*) mcWSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
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
			plot(dataSetPlus,mcWSetPlus,mcZSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus,mcJ1Events,mcJ2Events,mcJ3Events, 
			fQCDBkg, grBkgZPlus, Cw,cWSystErr,qcdSystError,zBosSystError, hmDataPlus,hmQCDPlus,hmZPlus,hmMcSumPlus, baselineResult, 
			101, j, k,nCentralityBins, mtmax,  ncoll[k],ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],muonMt,centrality, 
			sSelPlus , sSelPlusEta,true,correctSpectra,doSubtractBkg,doPreSelKinematics);

			std::cout << "plotting for mu-: " << i << ":" << j << ":" << k <<std::endl;
			plot(dataSetMinus,mcWSetMinus,mcZSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus,mcJ1Events,mcJ2Events,mcJ3Events, 
			fQCDBkg, grBkgZMinus, Cw,cWSystErr,qcdSystError,zBosSystError, hmDataMinus,hmQCDMinus,hmZMinus,hmMcSumMinus,  baselineResult, 
			100, j, k,nCentralityBins, mtmax,  ncoll[k],ptBins[i], ptBins[i+1], etaBins[j],etaBins[j+1], centralityBins[k], centralityBins[k+1],muonMt,centrality, 
			sSelMinus , sSelMinusEta,true,correctSpectra,doSubtractBkg,doPreSelKinematics);
		    }
		    else {
			if(doPlotMt){
				std::cout << "plotting for mu+: " << i << ":" << j << ":" << k <<std::endl;
				plot(dataSetPlus,mcWSetPlus,mcZSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus,mcJ1Events,mcJ2Events,mcJ3Events, 
				fQCDBkg, grBkgZPlus, Cw,cWSystErr,qcdSystError,zBosSystError, hmDataPlus,hmQCDPlus,hmZPlus,hmMcSumPlus,  baselineResult, 
				102+i, j, k,nCentralityBins, mtmax,  ncoll[k],ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],muonMt,centrality, 
				sSelPlus , sSelPlusEta,true,correctSpectra,doSubtractBkg,doPreSelKinematics);

	 
				std::cout << "plotting for mu-: " << i << ":" << j << ":" << k <<std::endl;
				plot(dataSetMinus,mcWSetMinus,mcZSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus,mcJ1Events,mcJ2Events,mcJ3Events, 
				fQCDBkg, grBkgZMinus, Cw,cWSystErr,qcdSystError,zBosSystError, hmDataMinus,hmQCDMinus,hmZMinus,hmMcSumMinus,  baselineResult, 
				103+i, j, k,nCentralityBins, mtmax,  ncoll[k],ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],muonMt,centrality, 
				sSelMinus , sSelMinusEta, true,correctSpectra,doSubtractBkg,doPreSelKinematics);
			}
			if(doPlotPhi){ 
					plotWCandidateKinematic(dataSetPlus,mcWSetPlus,mcZSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus,
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgZPlus, (THStack*)_hstackDataPlus.At(0),(THStack*)_hstackQCDPlus.At(0), (THStack*)_hstackZPlus.At(0),(THStack*)_hstackMcSumPlus.At(0),  
					102+i, j, k,nCentralityBins, -1*TMath::Pi(), -1*TMath::Pi(),TMath::Pi() , Cw,ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], 
					centralityBins[k+1],muonPhi, "#mu^{+},0-80" , "0 #leq |#eta| < 2.5", 20.0,"#phi_{#mu}", true, correctSpectra,doSubtractBkg,doPreSelKinematics);

					plotWCandidateKinematic(dataSetMinus,mcWSetMinus,mcZSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus,
					mcJ1Events, mcJ2Events, mcJ3Events, fQCDBkg, grBkgZMinus, (THStack*)_hstackDataMinus.At(0),(THStack*)_hstackQCDMinus.At(0), 
					(THStack*)_hstackZMinus.At(0),(THStack*)_hstackMcSumMinus.At(0), 103+i, j, k,nCentralityBins, -1*TMath::Pi(), -1*TMath::Pi(),TMath::Pi() , 
					Cw,ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], 
					centralityBins[k+1],muonPhi,"#mu^{-},0-80" , "0 #leq |#eta| < 2.5", 20.0,"#phi_{#mu}", true, correctSpectra,doSubtractBkg,doPreSelKinematics);

			} if(doPlotEta){ 
					plotWCandidateKinematic(dataSetPlus,mcWSetPlus,mcZSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus,
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgZPlus, (THStack*)_hstackDataPlus.At(1),(THStack*)_hstackQCDPlus.At(1), (THStack*)_hstackZPlus.At(1),(THStack*)_hstackMcSumPlus.At(1),  
					102+i, j, k,nCentralityBins,-2.5, -2.5, 2.5 , Cw,ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],muonEta,
					"#mu^{+},0-80" , "0 #leq |#eta| < 2.5", 22,"#eta_{#mu}", true, correctSpectra,doSubtractBkg,doPreSelKinematics);

					plotWCandidateKinematic(dataSetMinus,mcWSetMinus,mcZSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus,
					mcJ1Events, mcJ2Events, mcJ3Events,fQCDBkg, grBkgZMinus, (THStack*)_hstackDataMinus.At(1),(THStack*)_hstackQCDMinus.At(1), 
					(THStack*)_hstackZMinus.At(1),(THStack*)_hstackMcSumMinus.At(1), 103+i, j, k,nCentralityBins,-2.5, -2.5, 2.5 , 
					Cw,ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],muonEta,
					"#mu^{-},0-80" , "0 #leq |#eta| < 2.5", 22,"#eta_{#mu}", true, correctSpectra,doSubtractBkg,doPreSelKinematics);

			} if(doPlotMPT){ 
					plotWCandidateKinematic(dataSetPlus,mcWSetPlus,mcZSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus,
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgZPlus, (THStack*)_hstackDataPlus.At(2),(THStack*)_hstackQCDPlus.At(2), (THStack*)_hstackZPlus.At(2),(THStack*)_hstackMcSumPlus.At(2),  
					102+i, j, k,nCentralityBins,0.0, 25.0, 120.0 , Cw,ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],missPt,
					"#mu^{+},0-80" , "0 #leq |#eta| < 2.5", 24,"#slash{p_{T}}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);

					plotWCandidateKinematic(dataSetMinus,mcWSetMinus,mcZSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus,
					mcJ1Events, mcJ2Events, mcJ3Events, fQCDBkg, grBkgZMinus, (THStack*)_hstackDataMinus.At(2),(THStack*)_hstackQCDMinus.At(2), 
					(THStack*)_hstackZMinus.At(2),(THStack*)_hstackMcSumMinus.At(2), 103+i, j, k,nCentralityBins,0.0, 25.0, 120.0 , 
					Cw,ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],missPt,
					"#mu^{-},0-80" , "0 #leq |#eta| < 2.5", 24,"#slash{p_{T}}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);
			} if(doPlotPt){ 
					plotWCandidateKinematic(dataSetPlus,mcWSetPlus,mcZSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus,
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgZPlus, (THStack*)_hstackDataPlus.At(3),(THStack*)_hstackQCDPlus.At(3), (THStack*)_hstackZPlus.At(3),(THStack*)_hstackMcSumPlus.At(3),  
					102+i, j, k,nCentralityBins,0.0, 25.0, 100.0 , Cw,ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],muonPt,
					"#mu^{+},0-80" , "0 #leq |#eta| < 2.5", 50.0,"p^{#mu}_{T}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);

					plotWCandidateKinematic(dataSetMinus,mcWSetMinus,mcZSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus,
					mcJ1Events, mcJ2Events, mcJ3Events, fQCDBkg, grBkgZMinus, (THStack*)_hstackDataMinus.At(3),(THStack*)_hstackQCDMinus.At(3), 
					(THStack*)_hstackZMinus.At(3),(THStack*)_hstackMcSumMinus.At(3), 103+i, j, k,nCentralityBins,0.0, 25.0, 100.0 , 
					Cw,ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],muonPt,
					"#mu^{-},0-80" , "0 #leq |#eta| < 2.5", 50.0,"p^{#mu}_{T}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);
			} if(doPlotIso){ 
					plotWCandidateKinematic(dataSetPlus,mcWSetPlus,mcZSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus,
					mcJ1Events, mcJ2Events, mcJ3Events, 
					fQCDBkg, grBkgZPlus, (THStack*)_hstackDataPlus.At(4),(THStack*)_hstackQCDPlus.At(4), (THStack*)_hstackZPlus.At(4),(THStack*)_hstackMcSumPlus.At(4),  
					102+i, j, k,nCentralityBins,0.0, 0.0, 0.35 , Cw,ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],isolation20,
					"#mu^{+},0-80" , "0 #leq |#eta| < 2.5", 50.0,"#Sigma p^{ID}_{T}/p^{#mu}_{T}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);

					plotWCandidateKinematic(dataSetMinus,mcWSetMinus,mcZSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus,
					mcJ1Events, mcJ2Events, mcJ3Events, fQCDBkg, grBkgZMinus, (THStack*)_hstackDataMinus.At(4),(THStack*)_hstackQCDMinus.At(4), 
					(THStack*)_hstackZMinus.At(4),(THStack*)_hstackMcSumMinus.At(4), 103+i, j, k,nCentralityBins,0.0, 0.0, 0.35 , 
					Cw,ncoll[k], ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],isolation20,
					"#mu^{-},0-80" , "0 #leq |#eta| < 2.5", 50.0,"#Sigma p^{ID}_{T}/p^{#mu}_{T}[GeV]", true, correctSpectra,doSubtractBkg,doPreSelKinematics);
			}
		   }

		}//charge bin

		if (!doCharge&&doCentrality&&doEta) {

			std::cout << "Plotting for centrality and eta classes for mu^{#pm} " << std::endl;
			TString sSel = "#mu^{#pm},"; sSel += sCentLow; sSel+="-"; sSel+= sCentUp;
			TString sSelEta = sEtaLow; sSelEta += " #leq "; sSelEta+= "|#eta| "; sSelEta+="< "; sSelEta += sEtaUp;

			plot(dataSubSet[i][j][k],mcWSubSet[i][j][k],mcZSubSet[i][j][k],mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k], 
				mcJ1Events, mcJ2Events, mcJ3Events, 
				fQCDBkg, grBkgZ, Cw,cWSystErr,qcdSystError,zBosSystError, hmData,hmQCD,hmZ,hmMcSum,  baselineResult, 
				104+i, j, k, nCentralityBins, mtmax, ncoll[k],ptBins[i], ptBins[i+1], etaBins[j],etaBins[j+1], centralityBins[k], centralityBins[k+1],muonMt,centrality, 
				sSel , sSelEta, true,correctSpectra,doSubtractBkg,doPreSelKinematics);

		}  //no charge bin
		
	   } //centrality
	  } //eta
	}//pt

	if(doCharge) {
		
		TString sCentLow = "";
		TString sCentUp = "";
		TString sEtaLow = "";
		TString sEtaUp = "";

 		sCentLow += 100*centralityBins[0]; //sCentLow.Remove(3);
		sCentUp += 100*centralityBins[nCentralityBins]; //sCentUp.Remove(3);

		sEtaLow += etaBins[0];
		sEtaUp += etaBins[nEtaBins];
   		    
		TString sSelPlus = "#mu^{+},";
		TString sSelMinus = "#mu^{-},";
		TString sSelPlusEta = sEtaLow;
		TString sSelMinusEta = sEtaLow;
		sSelPlus += sCentLow; sSelPlus+="-"; sSelPlus+= sCentUp;
		sSelMinus += sCentLow; sSelMinus+="-"; sSelMinus+= sCentUp;
		sSelMinusEta += "#leq"; sSelMinusEta+= "|#eta|"; sSelMinusEta+="<"; sSelMinusEta += sEtaUp;
		sSelPlusEta += "#leq"; sSelPlusEta+= "|#eta|"; sSelPlusEta+="<"; sSelPlusEta += sEtaUp;

		std::cout << "Plotting merged histos..." << std::endl;
	        //99 = mu^{pm}, 100 = mu-(inclusive), 101 = mu+(inclusive), 102 = mu+ (binned), 103 = mu- (binned), 104 = mu^{pm} (binned)
		if(doPlotMt){
			plotMerged(hmDataPlus,hmQCDPlus,hmZPlus,hmMcSumPlus, sSelPlus , sSelPlusEta, 101, correctSpectra);
			plotMerged(hmDataMinus,hmQCDMinus,hmZMinus,hmMcSumMinus, sSelMinus , sSelMinusEta, 100, correctSpectra);
		}
		std::cout << "Plotting merged kinematic histos..." << std::endl;
		if(doPlotPhi){	
			plotMergedKinematic((THStack*)_hstackDataPlus.At(0),(THStack*)_hstackQCDPlus.At(0),(THStack*)_hstackZPlus.At(0),(THStack*)_hstackMcSumPlus.At(0),
			-1*TMath::Pi(),TMath::Pi(),20.0,-1*TMath::Pi(),"#phi_{#mu}", sSelPlus , sSelPlusEta, 101, correctSpectra);

			plotMergedKinematic((THStack*)_hstackDataMinus.At(0),(THStack*)_hstackQCDMinus.At(0),(THStack*)_hstackZMinus.At(0),(THStack*)_hstackMcSumMinus.At(0),
			-1*TMath::Pi(),TMath::Pi(),20.0,-1*TMath::Pi(),"#phi_{#mu}", sSelMinus , sSelMinusEta, 100, correctSpectra);
		}
		if(doPlotEta){
			plotMergedKinematic((THStack*)_hstackDataPlus.At(1),(THStack*)_hstackQCDPlus.At(1),(THStack*)_hstackZPlus.At(1),(THStack*)_hstackMcSumPlus.At(1),
			-2.5, 2.5 ,22.0,-2.5,"#eta_{#mu}",sSelPlus , sSelPlusEta, 101, correctSpectra);

			plotMergedKinematic((THStack*)_hstackDataMinus.At(1),(THStack*)_hstackQCDMinus.At(1),(THStack*)_hstackZMinus.At(1),(THStack*)_hstackMcSumMinus.At(1),
			-2.5, 2.5 ,22.0,-2.5,"#eta_{#mu}", sSelMinus , sSelMinusEta, 100, correctSpectra);
		}
		if(doPlotMPT){
			plotMergedKinematic((THStack*)_hstackDataPlus.At(2),(THStack*)_hstackQCDPlus.At(2),(THStack*)_hstackZPlus.At(2),(THStack*)_hstackMcSumPlus.At(2),
			25.0, 200.0,40,0.0,"#slash{p_{T}}[GeV]",sSelPlus , sSelPlusEta, 101, correctSpectra);

			plotMergedKinematic((THStack*)_hstackDataMinus.At(2),(THStack*)_hstackQCDMinus.At(2),(THStack*)_hstackZMinus.At(2),(THStack*)_hstackMcSumMinus.At(2),
			25.0, 200.0,40, 0.0, "#slash{p_{T}}[GeV]",sSelMinus , sSelMinusEta, 100, correctSpectra);
		}	
		if(doPlotPt){
			plotMergedKinematic((THStack*)_hstackDataPlus.At(3),(THStack*)_hstackQCDPlus.At(3),(THStack*)_hstackZPlus.At(3),(THStack*)_hstackMcSumPlus.At(3),
			25.0, 100.0 ,50.0,0.0, "p^{#mu}_{T}[GeV]",sSelPlus , sSelPlusEta, 101, correctSpectra);

			plotMergedKinematic((THStack*)_hstackDataMinus.At(3),(THStack*)_hstackQCDMinus.At(3),(THStack*)_hstackZMinus.At(3),(THStack*)_hstackMcSumMinus.At(3),
			25.0, 100.0 ,50.0,0.0, "p^{#mu}_{T}[GeV]",sSelMinus , sSelMinusEta, 100, correctSpectra);
		}
		if(doPlotIso){
			plotMergedKinematic((THStack*)_hstackDataPlus.At(4),(THStack*)_hstackQCDPlus.At(4),(THStack*)_hstackZPlus.At(4),(THStack*)_hstackMcSumPlus.At(4),
			0.0, 0.35 ,50.0,0.0,"#Sigma p^{ID}_{T}/p^{#mu}_{T}[GeV]",sSelPlus , sSelPlusEta, 101, correctSpectra);

			plotMergedKinematic((THStack*)_hstackDataMinus.At(4),(THStack*)_hstackQCDMinus.At(4),(THStack*)_hstackZMinus.At(4),(THStack*)_hstackMcSumMinus.At(4),
			0.0, 0.35 ,50.0,0.0,"#Sigma p^{ID}_{T}/p^{#mu}_{T}[GeV]",sSelMinus , sSelMinusEta, 100, correctSpectra);
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
		plotMerged(hmData,hmQCD,hmZ,hmMcSum,sSel , sSelEta, 99, correctSpectra);

		if(doPlotPhi){	
			plotMergedKinematic((THStack*)_hstackData.At(0),(THStack*)_hstackQCD.At(0),(THStack*)_hstackZ.At(0),(THStack*)_hstackMcSum.At(0),-1*TMath::Pi(),TMath::Pi(),20.0,
			-1*TMath::Pi(),"#phi_{#mu}",sSel , sSelEta, 99, correctSpectra);
		}
		if(doPlotEta){
			plotMergedKinematic((THStack*)_hstackData.At(1),(THStack*)_hstackQCD.At(1),(THStack*)_hstackZ.At(1),(THStack*)_hstackMcSum.At(1),-2.5, 2.5 ,22.0,-2.5,"#eta_{#mu}",
					sSel , sSelEta, 99, correctSpectra);
		}
		if(doPlotMPT){
			plotMergedKinematic((THStack*)_hstackData.At(2),(THStack*)_hstackQCD.At(2),(THStack*)_hstackZ.At(2),(THStack*)_hstackMcSum.At(2),25.0, 120.0,24.0,0.0,"#slash{p_{T}}[GeV]",
					sSel , sSelEta, 99, correctSpectra);
		}	
		if(doPlotPt){
			plotMergedKinematic((THStack*)_hstackData.At(3),(THStack*)_hstackQCD.At(3),(THStack*)_hstackZ.At(3),(THStack*)_hstackMcSum.At(3),25.0, 100.0 ,50.0,0.0,"p^{#mu}_{T}[GeV]",
					sSel , sSelEta, 99, correctSpectra);
		}
		if(doPlotIso){
			plotMergedKinematic((THStack*)_hstackData.At(4),(THStack*)_hstackQCD.At(4),(THStack*)_hstackZ.At(4),(THStack*)_hstackMcSum.At(4),0.0, 0.35 ,50.0,0.0,
			"#Sigma p^{ID}_{T}/p^{#mu}_{T}[GeV]",sSel , sSelEta, 99, correctSpectra);
		}


	}

	outFile.cd();
	baselineResult.write(outFile);
	outFile.Close();
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
		delete _hstackMcSum.At(ik);
		delete _hstackMcSumPlus.At(ik);
		delete _hstackMcSumMinus.At(ik);

	}

}//WAnalysis


int main(){
	WAnalysis();
}
