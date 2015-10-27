#include "TCanvas.h"
#include "TDirectory.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
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
//#include "WAnalysisHIDep.C"
#include "AtlasUtils.C"
#include "AtlasStyle.C"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath> 
#include <sstream>

#include "EfficiencyCorrection.C"

TString format(float value) {
  std::stringstream svalue;
  svalue  << std::setprecision(2) << value;
  return svalue.str();
}

///////////////////////////////
//Z counts from data
//////////////////////////////
double getZEvents(double centralityLow, double centralityUpp, double etaLow, double etaUpp){

//values calculated from Fig 1&4 of HI Z paper

	if(etaLow>=0.0&&etaUpp<=0.25){
	
		return 567.0/2.0;
		
	}else if(etaLow>=0.25&&etaUpp<=0.5){

		return 567.0/2.0;
		
	}else if(etaLow>=0.5&&etaUpp<=0.75){

		return 475.0/2.0;
		
	}else if(etaLow>=0.75&&etaUpp<=1.0){

		return 475.0/2.0;
		
	}else if(etaLow>=1.0&&etaUpp<=1.25){

		return 497.0/2.0;
		
	}else if(etaLow>=1.25&&etaUpp<=1.5){

		return 497.0/2.0;
		
	}else if(etaLow>=1.5&&etaUpp<=1.75){

		return 414.0/2.0;
		
	}else if(etaLow>=1.75&&etaUpp<=2.0){

		return 414.0/2.0;
		
	}else if(etaLow>=2.0&&etaUpp<=2.25){

		return 240.0/2.0;
		
	}else if(etaLow>=2.25&&etaUpp<=2.5){

		return 240.0/2.0;
	} else if(etaLow>=0.0&&etaUpp<=2.5){
		if(centralityLow>=0.0&&centralityUpp<=0.05) {
			return 338;
		} else if (centralityLow>=0.05&&centralityUpp<=0.1){ 
			return 311.0;
		} else if (centralityLow>=0.1&&centralityUpp<=0.15) {
			return 244.0; 
		} else if (centralityLow>=0.15&&centralityUpp<=0.20) {
			return 195.0;
		} else if (centralityLow>=0.20&&centralityUpp<=0.40) {
			return 400.0;
		} else if (centralityLow>=0.40&&centralityUpp<=0.80) {
			return 136.0;
		} else if (centralityLow>=0.0&&centralityUpp<=0.8) {
			return 1625.0;
		} 
		else {
			std::cout << "WARNING: UNABLE TO RETURN Z COUNT." << std::endl;  	
			return -9999.;
		}
	}
	else {
		std::cout << "WARNING: UNABLE TO RETURN Z COUNT." << std::endl;  	
		return -9999.;
	}
}
///////////////////////////////
//QCD fraction in signal region from data
//////////////////////////////
double getQCDBkg(double centralityLow, double centralityUpp){

	if(centralityLow>=0.0&&centralityUpp<=0.2) {
		return 0.0345947;
	} else if (centralityLow>=0.2&&centralityUpp<=0.4){ 
		return 0.0642742;
	} else if (centralityLow>=0.4&&centralityUpp<=0.8) {
		return 0.118457; 
	} else if (centralityLow>=0.0&&centralityUpp<=0.8) {
		return 0.0826315;
	} else {
		std::cout << "WARNING: UNABLE TO LOCATE CENTRALITY BIN." << std::endl;  	
		return -9999.;
	}
}


/*double getQCDEvents(TTree* tree, TString scutsIsoMPTLo, TString scutsNonIsoMPTUp, TString scutsNonIsoMPTLo, const float mtmax){

	TH1F* hIsoMPTLo = new TH1F("hIsoMPTLo","hIsoMPTLo",50,0.0,mtmax);
	TH1F* hNonIsoMPTUp = new TH1F("hNonIsoMPTUp","hNonIsoMPTUp",50,0.0,mtmax);
	TH1F* hNonIsoMPTLo = new TH1F("hNonIsoMPTLo","hNonIsoMPTLo",50,0.0,mtmax);

	tree->Draw("nu_pt>>hIsoMPTLo",scutsIsoMPTLo,"h");		
	double intIsoMPTLo = hIsoMPTLo->Integral();
	std::cout << "QCD events in region A: " << intIsoMPTLo << std::endl;
	tree->Draw("nu_pt>>hNonIsoMPTUp",scutsNonIsoMPTUp,"hsame");		
	double intNonIsoMPTUp = hNonIsoMPTUp->Integral();
	std::cout << "QCD events in region B: " << intNonIsoMPTUp << std::endl;
	tree->Draw("nu_pt>>hNonIsoMPTLo",scutsNonIsoMPTLo,"hsame");		
	double intNonIsoMPTLo = hNonIsoMPTLo->Integral();
	std::cout << "QCD events in region C: " << intNonIsoMPTLo  << std::endl;


	double IsoMPTUp = intIsoMPTLo*(intNonIsoMPTUp/intNonIsoMPTLo);
	std::cout << "QCD events in region D(signal region): " << IsoMPTUp  << std::endl;

	return IsoMPTUp;
		
}*/
///////////////////////////////
//plot
//////////////////////////////
void plot(TFile* fDataSet, TFile* fMcWSet, TFile* fMcZSet,TFile* fMcJ1Set,TFile* fMcJ2Set, TFile* fMcJ3Set, 
		const int iMt, const int iEta, const int iCentrality, const float mtmax, TString cuts,double ncoll, double ptLow, double ptUpp, double etaLow, double etaUpp, 
		double centralityLow, double centralityUpp, TString sSel, TString sSel2 , bool addPercent = true, 
		bool correctSpectra = true, bool doSubtractBkg = true) {

	float nBins = 50;
	float mptLow = 0.0; float mptUpp = 250.0;
	if(correctSpectra) std::cout << "Plotting corrected mT. " << std::endl;
	//initialize histograms
  	// --- data sets ---
	TH1F* hdataSet = new TH1F("hdataSet","hdataSet",nBins,mptLow,mptUpp);
  	// --- W set ---
	TH1F* hmcWSet = new TH1F("hmcWSet","hmcWSet",nBins,mptLow,mptUpp);
  	// --- Z set ---
	TH1F* hmcZSet = new TH1F("hmcZSet","hmcZSet",nBins,mptLow,mptUpp);
  	// --- QCD set ---
	TH1F* hmcJ1Set = new TH1F("hmcJ1Set","hmcJ1Set",nBins,mptLow,mptUpp);
	TH1F* hmcJ2Set = new TH1F("hmcJ2Set","hmcJ2Set",nBins,mptLow,mptUpp);
	TH1F* hmcJ3Set = new TH1F("hmcJ3Set","hmcJ3Set",nBins,mptLow,mptUpp);
	TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",nBins,mptLow,mptUpp);


	double mtcutLow = 40.0;
	cuts+="&&pt>"; cuts += ptLow; cuts+="&&pt<"; cuts += ptUpp;
	cuts+="&&abs(eta)>"; cuts += etaLow; cuts+="&&abs(eta)<"; cuts += etaUpp;
	cuts+="&&centrality>"; cuts += centralityLow; cuts+="&&centrality<"; cuts += centralityUpp;
//	std::cout << cuts << std::endl;

	double missPtCut = 0.0; 
	double isoCut = 0.3;
	TString sLoPtTrigger = "&&(EF_mu4_MSonly_L1TE50||EF_mu4_L1VTE50)";
	//TString sLoPtTrigger = "&&efMatched&&(EF_mu4_MSonly_L1TE50||EF_mu4_L1VTE50)";
	TString sHiPtTrigger = "&&(EF_mu10_MSonly_EFFS_L1ZDC||EF_mu10_MSonly_EFFS_L1TE10||EF_mu10_MSonly_EFFS_L1TE20)&&(EF_mu10_MSonly_EFFS_L1ZDC_Matched20||EF_mu10_MSonly_EFFS_L1TE10_Matched20||EF_mu10_MSonly_EFFS_L1TE20_Matched20)";
	//TString sHiPtTrigger = "&&efMatched&&(EF_mu10_MSonly_EFFS_L1ZDC||EF_mu10_MSonly_EFFS_L1TE10||EF_mu10_MSonly_EFFS_L1TE20)";
	//TString sHiPtTrigger = "&&(EF_mu10_MSonly_EFFS_L1ZDC||EF_mu10_MSonly_EFFS_L1TE10||EF_mu10_MSonly_EFFS_L1TE20)";
	TString sMissPtLo = "&&nu_pt<"; sMissPtLo+=missPtCut;  TString sMissPtUp = "&&nu_pt>"; sMissPtUp+=missPtCut; 
	TString sIsoCut = "&&ptcone20/pt<"; sIsoCut+=isoCut; TString sNonIsoCut = "&&ptcone20/pt>"; sNonIsoCut+=isoCut;
	TString sIsoCut40 = "&&ptcone40/pt<"; sIsoCut40+=isoCut; TString sNonIsoCut40 = "&&ptcone40/pt>"; sNonIsoCut40+=isoCut;
	TString sMtCutLo = "&&mt>"; sMtCutLo+=mtcutLow; TString sMtCutUp = "&&mt<";sMtCutUp+=mtmax;

	TString scutsSig = cuts + sHiPtTrigger; scutsSig += sMissPtUp; scutsSig += sIsoCut; scutsSig +=sMtCutLo; scutsSig +=sMtCutUp;  
	std::cout << "Signal region: " << scutsSig << std::endl;

	TString scutsMc = cuts; scutsMc += sMissPtUp; scutsMc += sIsoCut; scutsMc +=sMtCutLo; scutsMc +=sMtCutUp; 
//	std::cout << "MC cuts: " << scutsMc << std::endl;

	TString scutsIsoMPTLo = cuts + sHiPtTrigger; scutsIsoMPTLo+=sMissPtLo; scutsIsoMPTLo+=sIsoCut;
//	std::cout << "Region A cuts: " << scutsIsoMPTLo << std::endl;
	TString scutsNonIsoMPTUp = cuts + sHiPtTrigger; scutsNonIsoMPTUp+=sMissPtUp; scutsNonIsoMPTUp+=sNonIsoCut;
//	std::cout << "Region B cuts: " << scutsNonIsoMPTUp << std::endl;
	TString scutsNonIsoMPTLo = cuts + sHiPtTrigger; scutsNonIsoMPTLo +=sMissPtLo; scutsNonIsoMPTLo+=sNonIsoCut;
//	std::cout << "Region C cuts: " << scutsNonIsoMPTLo << std::endl;

	// --- Fill the histograms ---
	TTree *tree = (TTree*)fDataSet->Get("tree");
	tree->Draw("nu_pt>>hdataSet",scutsSig,"pe");		
	//std::cout << "filling data histo with entries: " << hdataSet->GetEntries() << std::endl; exit(0);
	TH1F* hdataSetc = (TH1F*)hdataSet->Clone("hdataSetc");

  	double sigEventsUncorr = hdataSetc->Integral(11,50); //integrate from 40-200 GeV before Aw,Cw correction
	double effMt ;
	if(correctSpectra) {
		hdataSetc = correctEfficiencyMt(hdataSet, iMt, iCentrality, iEta);
		effMt = getEfficiencyMt(iMt, iCentrality, iEta);
	}
	else effMt = 1.0;
 
	double allEvents = hdataSetc->Integral();	
	//double allEvents = hdataSetc->GetEntries();
  	double sigEvents = hdataSetc->Integral(11,50); //integrate from 40-200 GeV



	TTree *treeMcJ1Set = (TTree*)fMcJ1Set->Get("tree");
	TTree *treeMcJ2Set = (TTree*)fMcJ2Set->Get("tree");
	TTree *treeMcJ3Set = (TTree*)fMcJ3Set->Get("tree");

	//event counting histos
	TString sCentrality = "centrality>"; sCentrality+= centralityLow; sCentrality+="&&centrality<"; sCentrality+=centralityUpp;
	treeMcJ1Set->Draw("centrality>>hmcJ1Cent",sCentrality,"hf");
	treeMcJ2Set->Draw("centrality>>hmcJ2Cent",sCentrality,"hf");
	treeMcJ3Set->Draw("centrality>>hmcJ3Cent",sCentrality,"hf");

	TString cutsQCD = scutsMc; //+ "&&prompt==";
	treeMcJ1Set->Draw("nu_pt>>hmcJ1Set",cutsQCD,"hf");		
	treeMcJ2Set->Draw("nu_pt>>hmcJ2Set",cutsQCD,"hf");		
	treeMcJ3Set->Draw("nu_pt>>hmcJ3Set",cutsQCD,"hf");		
	//weight the Jx samples according to cross-sections (AMI)
	double wtJ1 = 6.6523e-03*1.8770e-4/64.0e-3; double wtJ2 = 1.4941e-02*8.2788e-6/64.0e-3; double wtJ3 = 2.4284e-02*2.9426e-7/64.0e-3;
	double arrCentWidth = centralityUpp-centralityLow;
	double evData = 68.7e6; //number of minbias events
	double scaleFactor = arrCentWidth*ncoll*evData;

	hmcQCDSet->Add(hmcJ1Set,hmcJ2Set,(wtJ1)/(hmcJ1Cent->Integral())*scaleFactor, (wtJ2)/(hmcJ2Cent->Integral())*scaleFactor); 
	hmcQCDSet->Add(hmcQCDSet,hmcJ3Set,1.0,(wtJ3)/(hmcJ3Cent->Integral())*scaleFactor);
	TH1F* hmcQCDSetc = (TH1F*)hmcQCDSet->Clone("hmcQCDSetc");

	//double survivalFractionQCD = 0.0975; //percent of QCD in W signal region from fit over ncoll
	double survivalFractionQCD = getQCDBkg(centralityLow, centralityUpp) ; //percent of QCD in W signal region
	double dataQCDEvents = survivalFractionQCD*sigEvents;
  	double mcQCDEvents = hmcQCDSetc->Integral(11,50);
	//number of QCD events from data using "ABCD" method in iso-misspt space 
	//double dataQCDEvents = getQCDEvents(tree, scutsIsoMPTLo, scutsNonIsoMPTUp, scutsNonIsoMPTLo, mtmax ) ;  
        //normalize to data
	hmcQCDSetc->Scale(dataQCDEvents/mcQCDEvents);
	//hmcQCDSetc->Scale(dataQCDEvents/mcQCDEvents/effMt);
	//hmcQCDSetc->Scale(0.0/allEvents);
	hmcQCDSetc->SetFillColor(kAzure-9);

	TTree *treeMcZSet = (TTree*)fMcZSet->Get("tree");
	TString cutsZ = scutsMc + "&&prompt==23";
	treeMcZSet->Draw("nu_pt>>hmcZSet",cutsZ,"hf");		
	TH1F* hmcZSetc = (TH1F*)hmcZSet->Clone("hmcZSetc");
	double survivalProb ; 
	if(iMt==102) survivalProb = 0.0522884; //mu+
	else if(iMt==103) survivalProb = 0.0412539; //mu-
	else survivalProb = 0.09577; //survival probability of Zs averaged over all centrality classes
	double bkgZEventFrac = survivalProb*getZEvents(centralityLow, centralityUpp, etaLow, etaUpp)/sigEventsUncorr;
	//double dataZEvents = survivalProb*getZEvents(centralityLow, centralityUpp) ;
	double dataZEvents = bkgZEventFrac*sigEvents ;
  	double mcZEvents = hmcZSetc->Integral(11,50);
	//weight Z mc by per event yield from paper
	//double mcZEvents = 4.39e-6; 
	hmcZSetc->Scale(dataZEvents/mcZEvents);
	//hmcZSetc->Scale(0.0/mcZEvents);
	//hmcZSetc->Scale(4.39e-6);
	hmcZSetc->Add(hmcZSetc,hmcQCDSetc);
	hmcZSetc->SetFillColor(kRed);

	TTree *treeMcWSet = (TTree*)fMcWSet->Get("tree");
	TString cutsW = scutsMc + "&&prompt==24";
	treeMcWSet->Draw("nu_pt>>hmcWSet",cutsW,"hf");		
	TH1F* hmcWSetc = (TH1F*)hmcWSet->Clone("hmcWSetc");
        //normalize to data
  	double mcWEvents = hmcWSetc->Integral(11,50);
	double sigEventsSub = sigEvents-dataQCDEvents-dataZEvents;
	hmcWSetc->Scale(sigEventsSub/mcWEvents);
	//add bkg contribution
	hmcWSetc->Add(hmcWSetc,hmcZSetc);
	mcWEvents = hmcWSetc->Integral(11,50);



	//now that the histos have been filled and scaled,
	//we now plot them on a canvas
	
	
  	TLegend* leg = new TLegend(0.648, 0.666, 0.918, 0.816);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);

	leg->AddEntry(hdataSetc, "Data 2011", "pe");
	leg->AddEntry(hmcWSetc, "W#rightarrow#mu#nu", "f");
	leg->AddEntry(hmcQCDSetc, "QCD", "f");
	leg->AddEntry(hmcZSetc, "Z#rightarrow#mu#mu", "f");

  	TCanvas* cdatamt = new TCanvas("cdatamt","datamt",600,600);
  	hdataSetc->GetXaxis()->SetTitle("#eta_{#mu}"); 
        float entriesPerBin = (mptUpp-mptLow)/nBins;
	TString sY = "Events/"; TString sentriesPerBin = format(entriesPerBin); sY+=sentriesPerBin;
  	hdataSetc->GetYaxis()->SetTitle(sY); 
  	//if(correctSpectra)hdataSetc->GetYaxis()->SetRangeUser(0.1,2.1e3);
	//hdataSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()+450.0); 
	//hdataSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()+350.0); 
	hdataSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()+50.0); 
	hdataSetc->Draw("pe");
	hmcWSetc->Draw("hfsame");
	hmcZSetc->Draw("hfsame");
	hmcQCDSetc->Draw("hfsame");
	hdataSetc->Draw("pesame");
        
	leg->Draw();

	TLatex l;
	l.SetNDC();
	l.DrawLatex(0.18,0.83,sSel + ( addPercent ? "%" : "" ));
	l.DrawLatex(0.169,0.767,sSel2);
	l.DrawLatex(0.65,0.85,"#sqrt{s_{NN}}=2.76 TeV");
	l.SetTextSize(0.034);
	l.DrawLatex(0.19,0.67,"#int Ldt #approx 0.140 nb^{-1}"); 

	ATLAS_LABEL(0.17,0.89,1);
	myText(0.33,0.89, (Color_t)kBlack, (char*)("Internal"));

	TString plotNameLog = "dataMPT_"; plotNameLog+=sSel; plotNameLog+="_"; plotNameLog+=sSel2; plotNameLog+="_"; plotNameLog+="Log"; //.png";
	TString plotNameLin = "dataMPT_"; plotNameLin+=sSel; plotNameLin+="_"; plotNameLin+=sSel2; plotNameLin+="_"; plotNameLin+="Lin"; //.png";

	cdatamt->Print(plotNameLin.ReplaceAll("#","")+".png");
	cdatamt->Print(plotNameLin.ReplaceAll("#","")+".eps"); 
	cdatamt->Print(plotNameLin.ReplaceAll("#","")+".pdf"); 
	cdatamt->Print(plotNameLin+".root"); 

  	cdatamt->SetLogy(true); hdataSetc->GetYaxis()->SetRangeUser(0.1,2.1e5); cdatamt->Update();
	/*cdatamt->Print(plotNameLog.ReplaceAll("#","")+".png");
	cdatamt->Print(plotNameLog.ReplaceAll("#","")+".eps"); 
	cdatamt->Print(plotNameLog.ReplaceAll("#","")+".pdf"); 
	cdatamt->Print(plotNameLog+".root"); 
*/
	std::cout << "Clean up" << std::endl;
 	//delete cDummy;
	delete hdataSetc;
	delete hmcWSetc;
	delete hmcZSetc;
	delete hmcQCDSetc;
	delete hdataSet;
	delete hmcWSet;
	delete hmcZSet;
	delete hmcQCDSet;
	//delete cdatamt;	
	delete leg;
	delete tree ;
	delete treeMcWSet ;
	delete treeMcZSet ;
	delete treeMcJ1Set ;
	delete treeMcJ2Set ;
	delete treeMcJ3Set ;
	std::cout << "Plotting ends here." << std::endl;
}//plot

///////////////////////////////
//WAnalysis
//////////////////////////////

void WCandidateMptPlotter(){
	
	bool doCharge =  false ;
	bool doCentrality = false ;
	bool doEta = false;
	bool correctSpectra = true; //corrects spectra using Aw,Cw factors
	bool doSubtractBkg = true; 

	int nSigPdf = 1000000;
	TString nSigPdfString;

	  if(doCentrality&&doCharge&&doEta){
	      nSigPdfString = "CentChrgEta";
	    }
	  else if(doCentrality&&doCharge&&!doEta){
	      nSigPdfString = "CentChrgCorrectedBkgSub";
	    }
	  else if (doEta&&doCentrality&&!doCharge){
	      nSigPdfString = "CentEta";
	    }

	  else if (doEta&&doCharge&&!doCentrality){
	      nSigPdfString = "EtaChrgCorrectedBkgSub";
	   }

	  else if (!doEta&&doCharge&&!doCentrality){
	      nSigPdfString = "ChargeCorrectedBkgSub";
	   }
	  else{ nSigPdfString = "1000000CorrectedBkgSub";}

	float mtmax = 200.0;
	float ptmax = 200.0;
	 
	SetAtlasStyle();
	TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
	TString fileNameDataIn = baseString+"HardProbesFiles/HISingleMuonHP.12.09.2012";
	TFile* fDataSet = new TFile(fileNameDataIn+".root", "READ");
	if ( !fDataSet->IsOpen() ) {
	    std::cout << fDataSet << " not found!" << std::endl;
	    exit(0);
	}

	//data overlay
	TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuMCDataOverlay.2012.11.21";
	TFile* fMcWSet = new TFile(fileNameMCWIn+".root", "READ");

	if ( !fMcWSet->IsOpen() ) {
	    std::cout << fMcWSet << " not found!" << std::endl;
	    exit(0);
	}
	//HIJING overlay
	TString fileNameMCZIn = baseString+"MonteCarloFiles/Zmumu/HISingleMuonMC_PYTHIA_HIJING_Zmumu_11.28.2012";
	TFile* fMcZSet = new TFile(fileNameMCZIn+".root", "READ");

	if ( !fMcZSet->IsOpen() ) {
	    std::cout << fMcZSet << " not found!" << std::endl;
	    exit(0);
	}
	//J1 1 muon-filter 
	TString fileNameMCJ1In = baseString+"MonteCarloFiles/QCD/HISingleMuonMcJ1.2012.11.29";
	TFile* fMcJ1Set = new TFile(fileNameMCJ1In+".root", "READ");

	if ( !fMcJ1Set->IsOpen() ) {
	    std::cout << fMcJ1Set << " not found!" << std::endl;
	    exit(0);
	}
	//J2 1 muon-filter 
	TString fileNameMCJ2In = baseString+"MonteCarloFiles/QCD/HISingleMuonMcJ2.2012.11.29";
	TFile* fMcJ2Set = new TFile(fileNameMCJ2In+".root", "READ");

	if ( !fMcJ2Set->IsOpen() ) {
	    std::cout << fMcJ2Set << " not found!" << std::endl;
	    exit(0);
	}
	//J3 1 muon-filter 
	TString fileNameMCJ3In = baseString+"MonteCarloFiles/QCD/HISingleMuonMcJ3.2012.11.29";
	TFile* fMcJ3Set = new TFile(fileNameMCJ3In+".root", "READ");

	if ( !fMcJ3Set->IsOpen() ) {
	    std::cout << fMcJ3Set << " not found!" << std::endl;
	    exit(0);
	}

	// --- Set pt and eta bins ---
	std::vector<double> ptBins;
	ptBins.push_back(25.0);
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
	centralityBins.push_back(0.00);
	if (doCentrality) {
	centralityBins.push_back(0.05);
	centralityBins.push_back(0.10);
	centralityBins.push_back(0.15);
	centralityBins.push_back(0.20);
	centralityBins.push_back(0.40);
	//centralityBins.push_back(0.60);

	//ncoll
	ncoll.push_back(1683.3); //0-5
	ncoll.push_back(1318.0); //5-10
	ncoll.push_back(1035.4); //10-15
	ncoll.push_back(811.2); //15-20

//	ncoll.push_back(1212.0);//0-20
	ncoll.push_back(440.6); //20-40
	ncoll.push_back(77.8); //40-80

	}
	else  ncoll.push_back(361.6);//0-80
	centralityBins.push_back(0.80);
	const int nCentralityBins = centralityBins.size()-1;


	//base muon selection cuts
	TString cuts ="abs(charge)==1";
        cuts+="&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&ZDY==0";
        //cuts+="&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0";
	//std::cout << cuts << std::endl;

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


		/// --- plot inclusive spectra --
		if(!doCentrality && !doEta && !doCharge){
			plot(fDataSet, fMcWSet, fMcZSet, fMcJ1Set, fMcJ2Set, fMcJ3Set,   99, 0, 0, mtmax, cuts,ncoll[i],ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1]			, centralityBins[k], centralityBins[k+1], "#mu^{#pm},0-80" , "0 #leq |#eta| < 2.5", true,correctSpectra,doSubtractBkg);
		}

		
		if(doCharge){

		    TString sSelPlus = "#mu^{+},";
		    TString sSelMinus = "#mu^{-},";
		    TString sSelPlusEta = sEtaLow;
		    TString sSelMinusEta = sEtaLow;
		    sSelPlus += sCentLow; sSelPlus+="-"; sSelPlus+= sCentUp;
		    sSelMinus += sCentLow; sSelMinus+="-"; sSelMinus+= sCentUp;
		    sSelMinusEta += "#leq"; sSelMinusEta+= "|#eta|"; sSelMinusEta+="<"; sSelMinusEta += sEtaUp;
		    sSelPlusEta += "#leq"; sSelPlusEta+= "|#eta|"; sSelPlusEta+="<"; sSelPlusEta += sEtaUp;

		    TString cutsP = cuts + "&&charge==1";
		    std::cout << cutsP << std::endl;
		    TString cutsM = cuts + "&&charge==-1";
		    std::cout << cutsM << std::endl;

		    //NOTE: index i is important for correctly matching efficiencies with given bin.
		    //CHANGE WITH CAUTION!
	            //99 = mu^{pm}, 100 = mu-(inclusive), 101 = mu+(inclusive), 102 = mu+ (binned), 103 = mu- (binned), 104 = mu^{pm} (binned)
		    if(!doCentrality&&!doEta){ 
		    	
			std::cout << "plotting for mu+: " << i << ":" << j << ":" << k <<std::endl;
			plot(fDataSet, fMcWSet, fMcZSet, fMcJ1Set, fMcJ2Set, fMcJ3Set,   101, j, k, mtmax, cutsP, ncoll[i],ptBins[i], ptBins[i+1], etaBins[j],
                        	etaBins[j+1], centralityBins[k], centralityBins[k+1], sSelPlus , sSelPlusEta, true,correctSpectra,doSubtractBkg);
			std::cout << "plotting for mu-: " << i << ":" << j << ":" << k <<std::endl;
			plot(fDataSet, fMcWSet, fMcZSet, fMcJ1Set, fMcJ2Set, fMcJ3Set,   100, j, k, mtmax, cutsM, ncoll[i],ptBins[i], ptBins[i+1], etaBins[j],
                          	etaBins[j+1], centralityBins[k], centralityBins[k+1], sSelMinus , sSelMinusEta, true,correctSpectra,doSubtractBkg);
		    }
		    else {
			std::cout << "plotting for mu+: " << i << ":" << j << ":" << k <<std::endl;
			plot(fDataSet, fMcWSet, fMcZSet, fMcJ1Set, fMcJ2Set, fMcJ3Set,   102+i, j, k, mtmax, cutsP, ncoll[i],ptBins[i], ptBins[i+1], etaBins[j], 
			  	etaBins[j+1], centralityBins[k], centralityBins[k+1], sSelPlus , sSelPlusEta, true,correctSpectra,doSubtractBkg);
		    
		    	std::cout << "plotting for mu-: " << i << ":" << j << ":" << k <<std::endl;
		    	plot(fDataSet, fMcWSet, fMcZSet, fMcJ1Set, fMcJ2Set, fMcJ3Set,   103+i, j, k, mtmax, cutsM, ncoll[i],ptBins[i], ptBins[i+1], etaBins[j], 
				etaBins[j+1], centralityBins[k], centralityBins[k+1], sSelMinus , sSelMinusEta, true,correctSpectra,doSubtractBkg);
		   }

		}//charge bin

		if (!doCharge&&doCentrality&&doEta) {

			std::cout << "Plotting for centrality and eta classes. Inclusive for charges." << std::endl;
			TString sSel = "#mu^{#pm},"; sSel += sCentLow; sSel+="-"; sSel+= sCentUp;
			TString sSelEta = sEtaLow; sSelEta += " #leq "; sSelEta+= "|#eta| "; sSelEta+="< "; sSelEta += sEtaUp;

			plot(fDataSet, fMcWSet, fMcZSet, fMcJ1Set, fMcJ2Set, fMcJ3Set,   104+i, j, k, mtmax, cuts, ncoll[i],ptBins[i], ptBins[i+1], etaBins[j],etaBins[j+1], centralityBins[k], centralityBins[k+1], sSel , sSelEta, true,correctSpectra,doSubtractBkg);

		}  //no charge bin
		
	   } //centrality
	  } //eta
	}//pt

}//WAnalysis


int main(){
	WCandidateEtaPlotter();
}
