#ifndef WAnalysisTool_h
#define WAnalysisTool_h
 
#include "RooKeysPdf.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAbsArg.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooBinning.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TROOT.h"

#include <iostream>
#include <cmath>

using namespace RooFit ;

class WAnalysisTool
{
	private:
		RooDataSet* _dataSet;
		RooBinning _b;
		RooRealVar _var;
		RooWorkspace *w;
		RooRealVar _wt;
		TH1F* _hDataSet;
		THStack* _hsDataSetSum;
		double _cutLow;
		double _cutUpp;
		int _binLo;
		int _binUp;
		double _sigEventsUncorrected;			
		int _nCentralityBins;
		int _nEtaBins;
		int _nPtBins;
		int _charge;
        bool m_doBkgSubtraction;
        float m_NZBkg;
        float m_NTauBkg;
        float m_NBkgA;
        float m_NBkgB;
        bool m_dataSetIsCorrected;
        float m_correctedSum;
        float m_correctedSumStatError;

	public:		

		WAnalysisTool(RooDataSet* dataSet);
		WAnalysisTool(RooDataSet* dataSet,RooRealVar& var);
		// copy ctor
		WAnalysisTool(const WAnalysisTool& CopySource);
	 	WAnalysisTool();	
		~WAnalysisTool();

		void weightDS(double wt); 
		void fit(TFile* fppPdf, TFile* fnpPdf, double cutLow, double cutUpp);
		TH1F* GetHistogram();
		TH1F* GetSummedHistogram();
		RooDataSet* GetDataSet();
		TCanvas* plotOnCanvas(RooRealVar& var);
		float numSignalCandidates();		
        void correctSignalCandidates(float CwA=1.0, float CwB=1.0);
        float getNumCorrectedSignalCandidates();
        float getStatErrorCorrectedSignalCandidates();

		void buildTemplate(RooRealVar& var, const char* fileName);
		void AssignDataSet(RooDataSet* dataSet);
		void AddToStack(TH1F* h);	
		void FillHistogram(RooBinning b, RooRealVar& var, double cutLow, double cutUpp);	
        float getBkgSubtractedSignalCandidates();
        void setZBackground(float NbkgA=0.0, float NbkgB=0.0);
        void setTauBackground(float NbkgA=0.0, float NbkgB=0.0);

		void ApplySelectionCuts(RooFormulaVar cutsSig, int charge = 0, 
				   TString muMinus="chargeCategory==chargeCategory::muMinus",TString muPlus="chargeCategory==chargeCategory::muPlus");

		void SelectPtEtaCentrality(double ptLow, double ptUpp, double etaLow, double etaUpp, double centralityLow, double centralityUpp,
					       TString sPt = "muonPt", TString sEta = "muonEta", TString sCentrality = "centrality");


		void setEtaBins(int etaBins){
			_nEtaBins = etaBins;
		}

		void setCentralityBins(int etaBins){
			_nCentralityBins = etaBins;
		}

		void setPtBins(int etaBins){
			_nPtBins = etaBins;
		}
		
		void setCharge(int charge){
			_charge = charge;
		}
};
#endif
