#include "TDirectory.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TMath.h"
#include "TF1.h"
#include "TProfile.h"
#include "TLegend.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath> 

#include "AtlasUtils.C"
#include "AtlasStyle.C"



#include "correctionFactorsDep.C"

///////////////////////////////////////////
//This macro calculates the fraction of MUONS
//generated from Zs that survive W selection
//per Z decay;
//This is later in the analysis multiplied
//by the number of Zs produced in a given 
//centrality bin from the data to give
//the number of muons in eta bin i, centrality bin
//j that survive W selection.
//////////////////////////////////////////

///////////////////////////////////////////
//writeToSpreadsheet
//////////////////////////////////////////
void writeToSpreadsheet(std::ostream& outputFile, int ieta, int icent, double fractionZBkg){
//	outputFile << ieta << " " << icent << " " << scaleFactor << " " << errStat << std::endl;
	outputFile << fractionZBkg << std::endl;
}



///////////////////////////////////////////////////////////////////////////////
//plot difference in charge as function of Ncoll 
///////////////////////////////////////////////////////////////////////////////
void plotChargeDiffNcoll(TGraph* grDiff, TGraphAsymmErrors* grPlus, TGraphAsymmErrors* grMinus, int ieta, int icent, double npart){

	double* yPlus = grPlus->GetY();
	double* yMinus = grMinus->GetY();

	double yTempPlus = yPlus[icent]; 
	double yTempMinus = yMinus[icent]; 

	double diff = fabs(yTempPlus-yTempMinus);
	std::cout << "Difference in f_{Z} at eta bin " << ieta << " and centrality bin " << icent << " = " << diff << std::endl;
	grDiff->SetPoint(icent,npart,diff);

}

///////////////////////////////////////////////////////////////////////////////
//plot difference in charge as function of eta 
///////////////////////////////////////////////////////////////////////////////
void plotChargeDiffEta(TGraph* grDiff, TGraphAsymmErrors* grPlus, TGraphAsymmErrors* grMinus, int ieta, int icent, double etaLo, double etaUpp){

	double  xPt = etaLo+fabs(etaUpp-etaLo)/2.0;
	double* yPlus = grPlus->GetY();
	double* yMinus = grMinus->GetY();

	double yTempPlus = yPlus[ieta]; 
	double yTempMinus = yMinus[ieta]; 

	double diff = fabs(yTempPlus-yTempMinus);
	std::cout << "Difference in f_{Z} at eta bin " << ieta << " and centrality bin " << icent << " = " << diff << std::endl;
	grDiff->SetPoint(ieta,xPt,diff);

}


void plotFraction(TGraphErrors* grFracCent, TGraphErrors* grFracEta, TString sCh){

	TLatex l;
	TCanvas* cFracCent = new TCanvas("cFracCent","cFracCent",600,600);
	l.SetNDC();
	l.DrawLatex(0.39,0.83,sCh);
	grFracCent->GetXaxis()->SetTitle("#LT N_{coll} #GT");	
	grFracCent->GetYaxis()->SetTitle("#frac{N_{Z}^{sur}}{N_{Z}^{Gen}}");
	grFracCent->Draw("ape");	

	TF1 f0a = TF1("f0a","[0]",0,1700);
	f0a.SetLineStyle(kDashed);
	grFracCent->Fit("f0a");
	f0a.Draw("same") ; 
	cout << "chi2 = " << f0a.GetChisquare() << "/" << f0a.GetNDF() << ", p = " << f0a.GetProb() << endl;

	cFracCent->Update();
	cFracCent->Print(sCh+"_centrality.pdf"); cFracCent->Print(sCh+"_centrality.root");

	TCanvas* cFracEta = new TCanvas("cFracEta","cFracEta",600,600);
	l.SetNDC();
	l.DrawLatex(0.39,0.83,sCh);
	grFracEta->GetXaxis()->SetTitle("|#eta|");	
	grFracEta->GetYaxis()->SetTitle("#frac{N_{Z}^{sur}}{N_{Z}^{Gen}}");
	grFracEta->Draw("ape");	

	TF1 f0b = TF1("f0b","[0]",0,2.5);
	f0b.SetLineStyle(kDashed);
	grFracEta->Fit("f0b");
	f0b.Draw("same") ; 
	cout << "chi2 = " << f0b.GetChisquare() << "/" << f0b.GetNDF() << ", p = " << f0b.GetProb() << endl;

	cFracEta->Update();
	cFracEta->Print(sCh+"_eta.pdf"); cFracEta->Print(sCh+"_eta.root");
}

void Write(TFile* outFile, TGraph2DErrors* grFrac, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing T2DGraph to root file..." << std::endl;
	  grFrac->Write(sFile);
	  std::cout << "Done." << std::endl;
    	  gDirectory = dir;
}

void Write(TFile* outFile, TObject* grFrac, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing TGraph to root file..." << std::endl;
	  grFrac->Write(sFile);
	  std::cout << "Done." << std::endl;
    	  gDirectory = dir;
}

void plot(TFile* fMcZSet, TGraph2DErrors* grFracEtaCent, TGraphAsymmErrors* grFracCent, TGraphAsymmErrors* grFracEta,TGraphAsymmErrors* grTotZEta,
		const float mtmax, TString cuts, TString sGenChargeCut, int index, double ptLow, double ptUpp, double etaLow, double etaUpp, 
		double centralityLow, double centralityUpp, double ncoll, int ieta, int icent, int iCh, TString sSel, TString
        sSel2, bool doEta, std::ostream& output, 
        bool doMptSigmaDown=false, bool doMptSigmaUp=false, double missPtCut = 25.0, double mtCut = 40.0){

		
		TH1F* hmcZSetFull = new TH1F("hmcZSetFull","hmcZSetFull",80,0.0,mtmax);
		TH1F* hmcZSetCut = new TH1F("hmcZSetCut","hmcZSetCut",80,0.0,mtmax);
		TH1F* hmcZSetPreSelCut = new TH1F("hmcZSetPreSelCut","hmcZSetPreSelCut",80,0.0,mtmax);
		TH1F* hmcZSetPreSelAllFidCut = new TH1F("hmcZSetPreSelAllFidCut","hmcZSetPreSelAllFidCut",80,0.0,mtmax);
		TH1F* hmcZSet = new TH1F("hmcZSet","hmcZSet",80,0.0,mtmax);


		TString sCentCuts = "centrality>"; sCentCuts+=centralityLow; sCentCuts+="&&centrality<";  sCentCuts+=centralityUpp;

		double mtcutLow = mtCut;
		cuts+="&&pt>"; cuts += ptLow; cuts+="&&pt<"; cuts += ptUpp;
		cuts+="&&abs(eta)>"; cuts += etaLow; cuts+="&&abs(eta)<"; cuts += etaUpp;
		cuts+="&&centrality>"; cuts += centralityLow; cuts+="&&centrality<"; cuts += centralityUpp;
		
        ///missPt systematic (+-10.85GeV)
		//double missPtCut = 35.85; 
		//double missPtCut = 14.15; 
		double missPtMax = 9000.0; 
        ///Nominal value used for this analysis
		double isoCut = 0.1;
        //use for ptcone20ID3 systematics
//		double isoCut = 0.2;

		TString sMissPtLo ; 
        TString sMissPtUp ;
		TString sMtCutLo ; 
        TString sMtCutUp;

        if(doMptSigmaDown) {
            std::cout << "WARNING: Using 2GeV lower track pT threshold mpt for systematic study. " << std::endl;
            //sMissPtLo = "&&nu_pt2000Nominal<"; sMissPtLo+=missPtMax;  sMissPtUp = "&&nu_pt2000Nominal>"; sMissPtUp+=missPtCut;
            sMissPtLo = "&&nu_pt2000Nominal<"; sMissPtLo+=missPtMax;  sMissPtUp = "&&nu_pt2000Nominal>"; sMissPtUp+=missPtCut;
            sMtCutLo = "&&TMath::Sqrt(2.0*nu_pt2000Nominal*pt*(1.0-TMath::Cos(phi-nu_phi2000Nominal)))>"; sMtCutLo+=mtcutLow;
        }
        else if(doMptSigmaUp) {
            std::cout << "WARNING: Using smeared mpt for systematic study. " << std::endl;
            sMissPtLo = "&&nu_pt4000Nominal<"; sMissPtLo+=missPtMax;  sMissPtUp = "&&nu_pt4000Nominal>"; sMissPtUp+=missPtCut;
            //sMissPtLo = "&&nu_ptMod<"; sMissPtLo+=missPtMax;  sMissPtUp = "&&nu_ptMod>"; sMissPtUp+=missPtCut;
            sMtCutLo = "&&TMath::Sqrt(2.0*nu_pt4000Nominal*pt*(1.0-TMath::Cos(phi-nu_phi4000Nominal)))>"; sMtCutLo+=mtcutLow;
        }
        
        else {
            std::cout << "Using nominal mpt definition. " << std::endl;
            sMissPtLo = "&&nu_pt<"; sMissPtLo+=missPtMax;  sMissPtUp = "&&nu_pt>"; sMissPtUp+=missPtCut;
            sMtCutLo = "&&mt>"; sMtCutLo+=mtcutLow; sMtCutUp = "&&mt<";sMtCutUp+=mtmax;
        }
        ///Nominal isolation definition used in analysis
		TString sIsoCut = "&&ptcone20ID3/pt<"; sIsoCut+=isoCut; 
        //systematics for ptcone
//		TString sIsoCut = "&&ptcone30ID3/pt<"; sIsoCut+=isoCut; 
		//AMI cross-section and efficiency

		TString scutsMc = cuts; scutsMc += sMissPtUp;scutsMc += sMissPtLo; scutsMc += sIsoCut; scutsMc +=sMtCutLo;  
		double wtZ = 0.99843*2.5743e-10/64.0e-3;

		double mbEvents = 68.7e6;
		double arrCentWidth = centralityUpp-centralityLow;
		double scaleFactor = arrCentWidth*ncoll*mbEvents;

		TTree *treeMcZSet = (TTree*)fMcZSet->Get("tree");

		TString sFull = sCentCuts + "&&prompt==23";
		treeMcZSet->Draw("pt>>hmcZSetFull",sFull,"hf");		

		TH1F* hmcZSetFullc = (TH1F*)hmcZSetFull->Clone("hmcZSetFullc");

		TString sGenCuts = sCentCuts; sGenCuts+="&&abs(mc_mu_gen_mothertype)==23&&abs(mc_mu_gen_type)==13"; 
		TString sFullFid = sGenCuts;
        
        sGenCuts+=sGenChargeCut;

		TString sGenEtaCuts; 

	//	if(doEta) { 
			sGenEtaCuts = "&&abs(mc_mu_gen_eta)>"; sGenEtaCuts+=etaLow; sGenEtaCuts+="&&abs(mc_mu_gen_eta)<"; sGenEtaCuts+=etaUpp;
			sGenCuts+=sGenEtaCuts;
	//	}

		std::cout << "Generated muon cuts : " << sGenCuts << std::endl;

        ///Number of generated muons in eta bin i, centrality bin j
		treeMcZSet->Draw("mc_mu_gen_pt>>hmcZSet",sGenCuts);
		TH1F* hmcZSetc = (TH1F*)hmcZSet->Clone("hmcZSetc");
		hmcZSetc->SetFillColor(kRed);

  		double mcZGenEvents; 
  		double mcMuZGenEvents ;

		if(doEta) { 

			//mcZGenEvents = 0.5*hmcZSet->Integral(); //number of Zs within this eta slice
			//if(mcZGenEvents<1.0) mcZGenEvents=1.0;
			//else mcZGenEvents = ceil(0.5*hmcZSet->Integral());
            

			//std::cout << "Generated Z muons in eta "  << " " << etaLow <<"-"<<etaUpp << " : centrality " << centralityLow << "-" << centralityUpp 
			//	<< " = " << mcZGenEvents <<std::endl;
			//mcMuZGenEvents = mcZGenEvents;
		}
		else {
			//mcZGenEvents = 0.5*hmcZSet->Integral(); //if looking over entire eta space, divide muons froms Zs by 2 to get number of Zs
			//std::cout << "Generated Z muons over all eta : centrality " << centralityLow << "-" << centralityUpp << " = " << mcZGenEvents <<std::endl;
			//mcMuZGenEvents = 2.0*mcZGenEvents;
		}


        double nMuonsFromZGen = hmcZSet->Integral(); ///Number of single muons from Z decays (mu+,mu-, and mu+-)


		hmcZSetFullc->SetFillColor(kRed);

		//TString cutsPreSelZ = cuts + "&&abs(mc_mu_gen_mothertype)==23&&abs(mc_mu_gen_type)==13&&massCB>66.&&massCB<120.";
		//TString cutsPreSelZ = cuts + "&&prompt==23&&massCB>66.&&massCB<120.";
		TString cutsZ = scutsMc; cutsZ+= "&&truthMatched_muid==1"; cutsZ+="&&ZDY==0";  
		std::cout << "W selection cuts: " << cutsZ << std::endl;

		//number of muons from Zs surving muon preselection cuts
		//within this eta slice and centrality class
		//treeMcZSet->Draw("massCB>>hmcZSetPreSelCut",cutsPreSelZ,"hf");
		//treeMcZSet->Draw("massCB>>hmcZSetPreSelAllFidCut",cutsPreSelAllFidZ,"hf");
		std::cout << "Full fiducial cuts = " << sFullFid << std::endl;
		//number of Zs over entire eta space
		treeMcZSet->Draw("mc_mu_gen_pt>>hmcZSetPreSelAllFidCut",sFullFid,"hf");

		//double bkgZpreSel = hmcZSetPreSelCut->Integral();
		//total number of generated Zs
		double bkgZpreSelAllFid = 0.5*hmcZSetPreSelAllFidCut->Integral();
        std::cout << "Number of generated Zs in centrality bin: " << centralityLow << "-" << centralityUpp << " = "
            << bkgZpreSelAllFid << std::endl;
		std::cout << "Percentage of Z bosons in the eta region of interest = " << etaLow << "-" << etaUpp << " " << nMuonsFromZGen/bkgZpreSelAllFid*100.0 << "%" << std::endl; 
		double xpt = (etaUpp-etaLow)/2.0+etaLow;
        ///Fraction of muons found in this eta bin out of all the number of decays 
        double effZMuonsInEta,errZMuonsInEtaUpp,errZMuonsInEtaLow;
        if(doEta) {

            Efficiency(nMuonsFromZGen,bkgZpreSelAllFid,0.683,effZMuonsInEta,errZMuonsInEtaLow, errZMuonsInEtaUpp);
            std::cout << "Sanity check: " << nMuonsFromZGen/bkgZpreSelAllFid << "=?" << effZMuonsInEta << std::endl;
            grTotZEta->SetPoint(ieta,xpt, effZMuonsInEta);
            errZMuonsInEtaUpp = errZMuonsInEtaUpp - effZMuonsInEta;
            errZMuonsInEtaLow = effZMuonsInEta - errZMuonsInEtaLow;
            grTotZEta->SetPointError(ieta,0.0,0.0,errZMuonsInEtaLow,errZMuonsInEtaUpp);
        }
        writeToSpreadsheet(output,ieta,icent,effZMuonsInEta);

/*		std::cout << "estimated percentage of Z muons surviving in centrality " << centralityLow << "-" << centralityUpp << 
			" BEFORE W selection = " << bkgZpreSel/mcZGenEvents*100 << "%" << std::endl; 
*/
        ///Number of generated muons from Zs surviving W selection cuts
		//treeMcZSet->Draw("mc_mu_gen_pt>>hmcZSetCut",cutsZ,"hf");		
		treeMcZSet->Draw("pt>>hmcZSetCut",cutsZ,"hf");		

		TH1F* hmcZSetCutc = (TH1F*)hmcZSetCut->Clone("hmcZSetCutc");
  		//double mcZSurvEvents = hmcZSetCutc->Integral();
  		double mcZSurvEvents = hmcZSetCutc->Integral();
//		std::cout << "Surviving Z muons after W selections in eta "  << " " << etaLow <<"-"<<etaUpp << " : centrality " << centralityLow << "-" << centralityUpp 
//			<< " " << mcZSurvEvents <<std::endl;
        double nMuonsFromZSurv = hmcZSetCutc->Integral();

		hmcZSetCutc->SetFillColor(kYellow);
  	
//		std::cout << "reco Z events in centrality " << centralityLow << "-" << centralityUpp << " = " << mcZGenEvents << std::endl; 
//		std::cout << "Z events (one-legged mu) surviving W selection cuts in centrality " << centralityLow << "-" << centralityUpp << " = " << mcZSurvEvents << std::endl; 
		std::cout << "estimated percentage of generated muons from Z decays surviving in centrality " << centralityLow << "-" << centralityUpp << " = " 
            << nMuonsFromZSurv/nMuonsFromZGen*100 << "%" << std::endl; 

		double e1 = TMath::Sqrt(nMuonsFromZSurv)/nMuonsFromZSurv*100.; double e2 = TMath::Sqrt(nMuonsFromZGen)/nMuonsFromZGen*100.;
		double eFrac = 0.01*TMath::Sqrt(e1*e1 + e2*e2);
		std::cout << nMuonsFromZSurv << " " << TMath::Sqrt(nMuonsFromZSurv) << std::endl;

        double effMuonSurvivors, errLow, errUpp;
        ///Bayesian errors
        Efficiency(nMuonsFromZSurv,nMuonsFromZGen,0.683,effMuonSurvivors,errLow, errUpp);
        ///Get distance from mean to obtain absolute errors
        errUpp = errUpp - effMuonSurvivors;
        errLow = effMuonSurvivors - errLow;

        std::cout << "Sanity check: " << effMuonSurvivors << " =? " << nMuonsFromZSurv/nMuonsFromZGen << std::endl;
		grFracCent->SetPoint(icent,ncoll,effMuonSurvivors);
        grFracCent->SetPointError(icent,0.0,0.0,errLow, errUpp);

		/*grFracCent->SetPoint(icent,ncoll,nMuonsFromZSurv/nMuonsFromZGen);
		grFracCent->SetPointError(icent,0.0,eFrac*nMuonsFromZSurv/nMuonsFromZGen);
        */

        grFracEta->SetPoint(ieta,xpt,effMuonSurvivors);
        grFracEta->SetPointError(ieta,0.0,0.0,errLow, errUpp);
		/*grFracEta->SetPoint(ieta,xpt,nMuonsFromZSurv/nMuonsFromZGen);
		grFracEta->SetPointError(ieta,0.0,eFrac*nMuonsFromZSurv/nMuonsFromZGen);
        */

		grFracEtaCent->SetPoint(index,ncoll,xpt, effMuonSurvivors);
//		grFracEtaCent->SetPoint(index,ncoll,xpt, nMuonsFromZSurv/nMuonsFromZGen);
        double errTemp; 
        if(errUpp>errLow) errTemp = errUpp;
        else errTemp = errLow;
		grFracEtaCent->SetPointError(index,0.0,0.0,errTemp);
        std::cout << "Background fraction in bin index" << index << " = " << grFracEtaCent->GetZ()[index] << " +/- " << grFracEtaCent->GetEZ()[index]<<std::endl;

		TLegend* leg = new TLegend(0.648, 0.722, 0.918, 0.915);
		leg->SetTextFont(gStyle->GetTextFont());
		leg->SetTextSize(gStyle->GetTextSize());
		leg->SetBorderSize(0);
		leg->SetFillColor(0);

		leg->AddEntry(hmcZSetc, "generated", "f");
		leg->AddEntry(hmcZSetCutc, "survivors", "f");

		TCanvas* cdatapt = new TCanvas("cdatapt","datamt",600,600);
		hmcZSetc->GetXaxis()->SetTitle("p_{T}[GeV]"); 
		hmcZSetc->GetYaxis()->SetTitle("Muons/2.5 GeV"); 
		cdatapt->SetLogy(); cdatapt->Update();
		hmcZSetc->Draw("hf");
		hmcZSetCutc->Draw("hfsame");
		hmcZSetc->Draw("sameaxis");
		
		leg->Draw();

		TLatex l;
		l.SetNDC();
		l.DrawLatex(0.53,0.37,sSel + "%" );
		l.DrawLatex(0.57,0.43,sSel2);

		myText(0.53,0.55, (Color_t)kBlack, (char*)("PYTHIA+HIJING"));

		TString plotNameLog = "ptEWBkg_"; plotNameLog+=sSel; plotNameLog+="_"; plotNameLog+=sSel2; plotNameLog+="_"; plotNameLog+="Log"; //.png";

		//cdatapt->Print(plotNameLog.ReplaceAll("|",",")+".pdf"); 
	
	delete hmcZSetFull;
	delete hmcZSetCut;
	//delete cdatapt;	
	delete leg;
	std::cout << "Plotting ends here." << std::endl;
}


void bkgZPlotter()
{
		
	bool doCentrality = true;
	bool doEta = true;
	bool doCharge = true ;
    ///Turn on for MPT systematics
    bool doMptSigmaDown = false;
    bool doMptSigmaUp = false;

    if(doMptSigmaUp) std::cout << "Doing MPT +1 sigma systematics." << std::endl;
    if(doMptSigmaDown) std::cout << "Doing MPT -1 sigma systematics." << std::endl;

	float mtmax = 300.0;
	float ptmax = 300.0;

	//HIJING overlay
	TString baseString = "/usatlas/u/tbales/scratch/";
	//TString baseString = "/tmp/tbalestr/";
	//TString fileNameMCZIn = baseString+"MonteCarloFiles/Zmumu/HISingleMuonMC_PYTHIA_HIJING_Zmumu_11.28.2012";
	//TString fileNameMCZIn = baseString+"MonteCarloFiles/Zmumu/HISingleMuonMCZmumu.12.30.2012";
    ///Nominal input for main analysis (HIJING)
	TString fileNameMCZIn = baseString+"MonteCarloFiles/Zmumu/HISingleMuonZmumuPYTHIAHIJINGOverlay.04.13.2013";
    // New data overlaid sample
	//TString fileNameMCZIn = baseString+"MonteCarloFiles/Zmumu/HISingleMuonZmumuPowPyDataOverlay.05.18.2014";

    TString fileNameMptMtCorrelation = "mptMtCorrelations.07.21.2013";
    TFile* _fMptMtCorrelation = NULL;
    TProfile* _pfxMptMtCorrelation = NULL;
    ///input for systematic study
    if(doMptSigmaDown||doMptSigmaUp) {
        fileNameMCZIn = baseString+"MonteCarloFiles/Zmumu/HISingleMuonZmumuPYTHIAHIJINGOverlay.07.13.2013";
        _fMptMtCorrelation = new TFile(fileNameMptMtCorrelation+".root","read");
        if(_fMptMtCorrelation!=0) std::cout << "Mpt-Mt correlation file: " << fileNameMptMtCorrelation << " opened." << std::endl;
        else exit(0);
    }

	//TString fileNameMCZIn = baseString+"HISingleMuonMCZmumu.12.30.2012";
	TString fileNameOut;

	if(doEta&&doCentrality&&doCharge) fileNameOut = "fractionZEtaChargeCent";
	else if(doEta&&doCentrality) fileNameOut = "fractionZEtaCent";
	else if(doEta) fileNameOut = "fractionZEta";
	else if(doCentrality) fileNameOut = "fractionZCent";
	else fileNameOut = "fractionZ";
	
	/// --- Open output file ---
	TDirectory *dir = gDirectory;
	TFile *outFile = new TFile(fileNameOut+".root","RECREATE");
	gDirectory = dir;

	TString spreadSheetName = "bkgZEtaSpreadSheet.csv";
	std::ofstream spreadSheet;
	spreadSheet.open(spreadSheetName);


	TFile* fMcZSet = new TFile(fileNameMCZIn+".root", "READ");
	std::cout << "Reading input file " << fileNameMCZIn << std::endl;

	if ( !fMcZSet->IsOpen() ) {
	    std::cout << fMcZSet << " not found!" << std::endl;
	    exit(0);
	}

	gROOT->LoadMacro("AtlasUtils.C");


	// --- Set pt and eta bins ---
	std::vector<double> ptBins;
	//systematics
	//double sigmaPt = 8.5;
	//std::cout << "Doing muon pt systematics -1 sigma" << std::endl;
	ptBins.push_back(25.0);
	ptBins.push_back(ptmax);
	const int nPtBins = ptBins.size()-1;

	std::vector<double> etaBins;
	etaBins.push_back(0.10);
	if (doEta) {
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

        etaBins.push_back(0.35);
        etaBins.push_back(0.6);
        etaBins.push_back(0.8);
        etaBins.push_back(1.05);
        etaBins.push_back(1.37);
        etaBins.push_back(1.52);
        etaBins.push_back(1.74);
        etaBins.push_back(2.1);
 
	}
//	etaBins.push_back(+2.50);
	etaBins.push_back(+2.40);

	const int nEtaBins = etaBins.size()-1;
	std::vector<double> centralityBins;
	std::vector <float> ncoll;

	centralityBins.push_back(0.0);
	if(doCentrality){
		centralityBins.push_back(0.05);
		centralityBins.push_back(0.1);
		centralityBins.push_back(0.15);
		centralityBins.push_back(0.2);
		centralityBins.push_back(0.4);

		ncoll.push_back(1683.3); //0-5
		ncoll.push_back(1318.0); //5-10
	//	ncoll.push_back(1500.6); //0-10
		ncoll.push_back(1035.4); //10-15
		ncoll.push_back(811.2); //15-20

	//	ncoll.push_back(1212.0);//0-20
		ncoll.push_back(440.6); //20-40
		ncoll.push_back(77.8); //40-80
	} else  ncoll.push_back(452.0);//0-80

	centralityBins.push_back(0.8);

	const int nCentralityBins = centralityBins.size()-1;

    ///MPT lower thresholds as a function of centrality
    std::vector <double> vMptLowCut;
    double mptLowCut = 25.0;
    //cross-check
    //double mptLowCut = 0.0;
    ///-1sigma
    if(doMptSigmaDown&&doCentrality){
    vMptLowCut.push_back(mptLowCut); //0-5%
    vMptLowCut.push_back(mptLowCut); //5-10%
    vMptLowCut.push_back(mptLowCut); //10-15%
    vMptLowCut.push_back(mptLowCut); //15-20%
    vMptLowCut.push_back(mptLowCut); //20-40%
    vMptLowCut.push_back(mptLowCut); //40-80%
    }
    else if(doMptSigmaDown) vMptLowCut.push_back(mptLowCut); //0-80%
    ///+1 sigma
    else if (doMptSigmaUp&&doCentrality){
    vMptLowCut.push_back(mptLowCut); //0-5%
    vMptLowCut.push_back(mptLowCut); //5-10%
    vMptLowCut.push_back(mptLowCut); //10-15%
    vMptLowCut.push_back(mptLowCut); //15-20%
    vMptLowCut.push_back(mptLowCut); //20-40%
    vMptLowCut.push_back(mptLowCut); //40-80%
    }
    else if (doMptSigmaUp) vMptLowCut.push_back(mptLowCut); //0-80%

	TGraphAsymmErrors* grFracCent = new TGraphAsymmErrors(nCentralityBins);
	TGraphAsymmErrors* grFracEta = new TGraphAsymmErrors(nEtaBins);
	TGraphAsymmErrors* grTotZEta = new TGraphAsymmErrors(nEtaBins);
	TGraphAsymmErrors* grTotZEtaPlus = new TGraphAsymmErrors(nEtaBins);
	TGraphAsymmErrors* grTotZEtaMinus = new TGraphAsymmErrors(nEtaBins);

	TGraphAsymmErrors* grFracCentPlus = new TGraphAsymmErrors(nCentralityBins);
	TGraphAsymmErrors* grFracEtaPlus = new TGraphAsymmErrors(nEtaBins);
	TGraphAsymmErrors* grFracCentMinus = new TGraphAsymmErrors(nCentralityBins);
	TGraphAsymmErrors* grFracEtaMinus = new TGraphAsymmErrors(nEtaBins);

	//plot charge difference for systematics
	TGraph* grDiffEta = new TGraph(nEtaBins);
	TGraph* grDiffNcoll  = new TGraph(nCentralityBins);

	TGraph2DErrors* grFracEtaCent = new TGraph2DErrors(nCentralityBins*nEtaBins);
	TGraph2DErrors* grFracEtaCentPlus = new TGraph2DErrors(nCentralityBins*nEtaBins);
	TGraph2DErrors* grFracEtaCentMinus = new TGraph2DErrors(nCentralityBins*nEtaBins);

	//base muon selection cuts
	TString cuts ="abs(charge)==1";
        cuts+="&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0";
	TString cutsP = cuts+"&&charge==1";
	TString cutsM = cuts+"&&charge==-1";
	//std::cout << cuts << std::endl;
	for ( int i = 0; i < nPtBins; i++ ) {
	  for ( int j = 0; j < nEtaBins; j++ ) {

		TString sNamePlus = "grCent";sNamePlus+=102;  sNamePlus+="eta";sNamePlus+=j;
		TString sNameMinus = "grCent";sNameMinus+=103; sNameMinus+="eta";sNameMinus+=j;
		TString sNameDiff = "grCentDiff"; sNameDiff+="eta";sNameDiff+=j;
		TString sName = "grCent";sName+=104; sName+="eta";sName+=j;

	    for ( int k = 0; k < nCentralityBins; k++ ){

		int index = j*nCentralityBins + k;
		std::cout << " plotting "<<i<<":"<<j<<":"<<k<<std::endl;
		TString sCentLow = "";
		TString sCentUp = "";
		TString sEtaLow = "";
		TString sEtaUp = "";

		sCentLow += 100*centralityBins[k]; //sCentLow.Remove(3);
		sCentUp += 100*centralityBins[k+1]; //sCentUp.Remove(3);

		sEtaLow += etaBins[j];
		sEtaUp += etaBins[j+1];

	    	TString sSel = "Z#rightarrow#mu^{+}#mu^{-} ";
	    	TString sSelEta = sEtaLow;
	    	sSel += sCentLow; sSel+="-"; sSel+= sCentUp;
	    	sSelEta += "#leq"; sSelEta+= "|#eta|"; sSelEta+="<"; sSelEta += sEtaUp;

		std::cout << "plotting for : " << i << ":" << j << ":" << k <<std::endl;
        ///Nominal mt cut
        double mtCut = 40.0;
        if(doMptSigmaDown||doMptSigmaUp) {
            mptLowCut = vMptLowCut[k]; 

            TString sMptMtCorrelation = "h2DMptMtCorrelationCent"; sMptMtCorrelation+=k; sMptMtCorrelation+="_pfx";
            std::cout << "Opening " << sMptMtCorrelation << " TProfilex histo." << std::endl;
            _pfxMptMtCorrelation = (TProfile*)_fMptMtCorrelation->Get(sMptMtCorrelation);
            if(_pfxMptMtCorrelation!=0) std::cout << sMptMtCorrelation << " opened." << std::endl;
            else exit(0);
            //mtCut = _pfxMptMtCorrelation->GetBinContent(_pfxMptMtCorrelation->FindBin(mptLowCut));  
        }
        std::cout << "Lower missing pT cut = " << mptLowCut << " GeV" << std::endl;
        std::cout << "Lower mT cut = " << mtCut << " GeV" << std::endl;

        if(doMptSigmaDown){
            plot(fMcZSet, grFracEtaCent, grFracCent, grFracEta,grTotZEta, mtmax, cuts, "&&abs(mc_mu_charge)==1", index, ptBins[i], ptBins[i+1], etaBins[j], 
			    etaBins[j+1], centralityBins[k], centralityBins[k+1], ncoll[k],j,k,i, sSel , sSelEta,
                doEta,spreadSheet,true,false,mptLowCut,mtCut);
        }
        else if(doMptSigmaUp){
                plot(fMcZSet, grFracEtaCent, grFracCent, grFracEta,grTotZEta, mtmax, cuts, "&&abs(mc_mu_charge)==1", index, ptBins[i], ptBins[i+1], etaBins[j], 
			        etaBins[j+1], centralityBins[k], centralityBins[k+1], ncoll[k],j,k,i, sSel , sSelEta,
                    doEta,spreadSheet,false,true,mptLowCut,mtCut);

        }

		else plot(fMcZSet, grFracEtaCent, grFracCent, grFracEta,grTotZEta, mtmax, cuts, "&&abs(mc_mu_charge)==1", index, ptBins[i], ptBins[i+1], etaBins[j], 
			etaBins[j+1], centralityBins[k], centralityBins[k+1], ncoll[k],j,k,i, sSel , sSelEta,
            doEta,spreadSheet);

		if(doCharge){
            ///Change lower MPT cut as fcn of centrality for systematic study
            if(doMptSigmaDown){
                std::cout << "#mu^{+}" << std::endl;
                plot(fMcZSet,grFracEtaCentPlus, grFracCentPlus, grFracEtaPlus,grTotZEtaPlus, mtmax, cutsP,"&&mc_mu_charge==+1", 
                index, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1], 
                ncoll[k],j,k,102, "#mu^{+},"+sSel , sSelEta, doEta,spreadSheet,true,false,mptLowCut,mtCut);

                std::cout << "#mu^{-}" << std::endl;
                plot(fMcZSet, grFracEtaCentMinus, grFracCentMinus, grFracEtaMinus,grTotZEtaMinus, mtmax, 
                cutsM,"&&mc_mu_charge==-1", index, ptBins[i], ptBins[i+1], etaBins[j], 
                etaBins[j+1], centralityBins[k], centralityBins[k+1], ncoll[k],j,k,103, "#mu^{-},"+sSel , sSelEta,
                doEta,spreadSheet,true,false,mptLowCut,mtCut);
            }
            else if(doMptSigmaUp){
                std::cout << "#mu^{+}" << std::endl;
                plot(fMcZSet,grFracEtaCentPlus, grFracCentPlus, grFracEtaPlus,grTotZEtaPlus, mtmax, cutsP,"&&mc_mu_charge==+1", 
                index, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1], 
                ncoll[k],j,k,102, "#mu^{+},"+sSel , sSelEta, doEta,spreadSheet,false,true,mptLowCut,mtCut);

                std::cout << "#mu^{-}" << std::endl;
                plot(fMcZSet, grFracEtaCentMinus, grFracCentMinus, grFracEtaMinus,grTotZEtaMinus, mtmax, 
                cutsM,"&&mc_mu_charge==-1", index, ptBins[i], ptBins[i+1], etaBins[j], 
                etaBins[j+1], centralityBins[k], centralityBins[k+1], ncoll[k],j,k,103, "#mu^{-},"+sSel , sSelEta,
                doEta,spreadSheet,false,true,mptLowCut,mtCut);
            }
 
           else{ 

                std::cout << "#mu^{+}" << std::endl;
                plot(fMcZSet,grFracEtaCentPlus, grFracCentPlus, grFracEtaPlus,grTotZEtaPlus, mtmax, cutsP,"&&mc_mu_charge==+1", index, ptBins[i], ptBins[i+1], etaBins[j], 
                etaBins[j+1], centralityBins[k], centralityBins[k+1], ncoll[k],j,k,102, "#mu^{+},"+sSel , sSelEta,
                doEta,spreadSheet);

                std::cout << "#mu^{-}" << std::endl;
                plot(fMcZSet, grFracEtaCentMinus, grFracCentMinus, grFracEtaMinus,grTotZEtaMinus, mtmax, cutsM,"&&mc_mu_charge==-1", index, ptBins[i], ptBins[i+1], etaBins[j], 
                etaBins[j+1], centralityBins[k], centralityBins[k+1], ncoll[k],j,k,103, "#mu^{-},"+sSel , sSelEta, doEta,spreadSheet);
           }
			//systematics
	    	        std::cout << "|mu^{-}-mu^{+}|" << std::endl;
			if(doEta){
			  std::cout << "|eta|" << std::endl;
			  plotChargeDiffEta(grDiffEta,grFracEtaPlus,grFracEtaMinus,j,k,etaBins[j], etaBins[j+1]);
			}

			if(doCentrality){
			  std::cout << "<Ncoll>" << std::endl;
			  plotChargeDiffNcoll(grDiffNcoll,grFracCentPlus,grFracCentMinus,j,k,ncoll[k]);
			}

		  }

		} //icent

        ///write ncoll distro per eta bin into TGraph
        if(doCentrality){
            if(doCharge){
			  Write(outFile,grFracCentPlus,sNamePlus);
			  Write(outFile,grFracCentMinus,sNameMinus);
			  Write(outFile,grDiffNcoll,sNameDiff);
            }

			Write(outFile,grFracCent,sName);
        }

	     } //ieta
	}
	
	if(doCharge){
 		TString sChargePlus = "Z,#mu^{+}";
		TString sChargeMinus = "Z,#mu^{-}";
		//plotFraction(grFracCentPlus,grFracEtaPlus,sChargePlus);
		//plotFraction(grFracCentMinus,grFracEtaMinus,sChargeMinus);
        if(doCentrality&&!doEta){
		    Write(outFile, grFracCentPlus, "fractionZCentPlus");
		    Write(outFile, grFracCentMinus, "fractionZCentMinus");
        }
        if(!doCentrality&&doEta){
		    Write(outFile, grFracEtaPlus, "fractionZEtaPlus");
		    Write(outFile, grFracEtaMinus, "fractionZEtaMinus");
		    Write(outFile, grTotZEtaPlus, "fractionTotalZEtaPlus");
		    Write(outFile, grTotZEtaMinus, "fractionTotalZEtaMinus");
        }
        else{
		    Write(outFile, grFracEtaCentPlus, "fractionZEtaCentPlus");
		    Write(outFile, grFracEtaCentMinus, "fractionZEtaCentMinus");
        }
		//systematics
		if(doEta) Write(outFile,grDiffEta,"diffFracZEta"); 
		if(doCentrality) {
			Write(outFile,grDiffNcoll,"diffFracZNcoll");
		}

	}

		TString sCharge = "Z,#mu^{pm}";
		//plotFraction(grFracCent,grFracEta,sCharge);

        if(doCentrality&&!doEta) Write(outFile, grFracCent, "fractionZCent");
		else if(!doCentrality&&doEta){
            Write(outFile, grFracEta, "fractionZEta");
		    Write(outFile, grTotZEta, "fractionTotalZEta");
        }
        else Write(outFile, grFracEtaCent, "fractionZEtaCent");

    std::cout << "Closing spreadsheet..." << std::endl;
	spreadSheet.close();
    std::cout << "Done." << std::endl;
/*	TLegend* leg = new TLegend(0.648, 0.622, 0.918, 0.815);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->AddEntry(grFracMuP,"#mu^{+}","p");
	leg->AddEntry(grFracMuM,"#mu^{-}","p");

	leg->Draw(); 
*/	
}
