#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"
#include "TDatime.h"
#include "TGraph.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <sstream>


////////////////////////
//Prevent trailing digits
//in figure labeling
////////////////////////
TString format(float value, int precision) {
  std::stringstream svalue;
  svalue  << std::setprecision(precision) << value;
  return svalue.str();
}

void writeToSpreadSheet(std::ostream& outputFile,TString range, TString syst1, TString syst2, TString syst3, TString syst4,
				TString syst5, TString syst6, TString syst7,TString syst8,TString syst9){

		outputFile << range << "," << syst1 << "," << syst2 << "," << syst3 << "," << syst4 << "," 
			<< syst5 << "," << syst6 << "," << syst7 << "," << syst8 << "," << syst9 << std::endl; 

}

void writeToSpreadSheet(std::ostream& outputFile,std::vector<TString> range, std::vector<TString> syst1, std::vector<TString> syst2, 
			std::vector<TString> syst3, std::vector<TString> syst4, std::vector<TString> syst5, 
			std::vector<TString> syst6, std::vector<TString> syst7, std::vector<TString> syst8,std::vector<TString> syst9){

		        for(int j = 0; j < range.size(); ++j){
                            if(j==range.size()-1){
				outputFile << "\t," << range[j]  << std::endl;
			    }
			    else{
				outputFile << "\t," << range[j] << "," ; 
			    }
			
			}	
		        for(int j = 0; j < syst1.size(); ++j){
                            if(j==syst1.size()-1){
				outputFile << syst1[j]  << std::endl;
			    }
			    else if(j==0) outputFile << "QCD" << "," << syst1[j] << "," ;
			    else{
				outputFile << syst1[j] << "," ; 
			    }
			
			}
		        for(int j = 0; j < syst2.size(); ++j){
                            if(j==syst2.size()-1){
				outputFile << syst2[j]  << std::endl;
			    }
		            else if(j==0) outputFile << "EW" << "," << syst2[j] << "," ;
			    else{
				outputFile << syst2[j] << "," ; 
			    }
			
			}
		        for(int j = 0; j < syst3.size(); ++j){
                            if(j==syst3.size()-1){
				outputFile << syst3[j]  << std::endl;
			    }
		 	    else if(j==0) outputFile << "other bkg" << "," << syst3[j] << "," ;
			    else{
				outputFile << syst3[j] << "," ; 
			    }
			
			}
		        for(int j = 0; j < syst4.size(); ++j){
                            if(j==syst4.size()-1){
				outputFile << syst4[j]  << std::endl;
			    }
		            else if(j==0) outputFile << "isolation" << "," << syst4[j] << "," ;
			    else{
				outputFile << syst4[j] << "," ; 
			    }
			
			}
		        for(int j = 0; j < syst5.size(); ++j){
                            if(j==syst5.size()-1){
				outputFile << syst5[j]  << std::endl;
			    }
		            else if(j==0) outputFile << "mpt reso" << "," << syst5[j] << "," ;
			    else{
				outputFile << syst5[j] << "," ; 
			    }
			
			}
		        for(int j = 0; j < syst6.size(); ++j){
                            if(j==syst6.size()-1){
				outputFile << syst6[j]  << std::endl;
			    }
			    else if(j==0) outputFile << "tracking" << "," << syst6[j] << "," ;
			    else{
				outputFile << syst6[j] << "," ; 
			    }
			
			}
		        for(int j = 0; j < syst7.size(); ++j){
                            if(j==syst7.size()-1){
				outputFile << syst7[j]  << std::endl;
			    }
			    else if(j==0) outputFile << "muon reco" << "," << syst7[j] << "," ;
			    else{
				outputFile << syst7[j] << "," ; 
			    }
			
			}
		        for(int j = 0; j < syst8.size(); ++j){
                            if(j==syst8.size()-1){
				outputFile << syst8[j]  << std::endl;
			    }
			    else if(j==0) outputFile << "pt reso" << "," << syst8[j] << "," ;
			    else{
				outputFile << syst8[j] << "," ; 
			    }
			
			}
		        for(int j = 0; j < syst9.size(); ++j){
                            if(j==syst9.size()-1){
				outputFile << syst9[j]  << std::endl;
			    }
			    else if(j==0) outputFile << "trigger" << "," << syst9[j] << "," ;
			    else{
				outputFile << syst9[j] << "," ; 
			    }
			
			}
}

void fillContainers(TString sFileIn, std::vector <double>& _vSignal,std::vector <double>& _vRaaScalingPercentSyst,
                    std::vector <double>& _vEWPercentSyst,std::vector <double>& _vOtherBkgPercentSyst, std::vector <double>& _vIsolationPercentSyst,
                    std::vector <double>& _vMPTResolPercentSyst,std::vector <double>& _vTrackingPercentSyst,
                    std::vector <double>& _vMuonRecoPercentSyst,std::vector <double>& _vPtResolPercentSyst,std::vector <double>& _vTriggerPercentSyst)
{   
	std::cout << "Opening file " << sFileIn << std::endl;
	std::ifstream s( sFileIn , std::ifstream::in );
	if(!s.is_open()){
		std::cout << "While opening a file an error was encountered" << std::endl;
		exit(1);
    } 
    else{
                std::cout << "File " << sFileIn << " is successfully opened" << std::endl;
    }

	for (unsigned int i = 0; i < _vSignal.size(); ++i){
        s >> _vSignal[i] >> _vRaaScalingPercentSyst[i] >> _vEWPercentSyst[i] >> _vOtherBkgPercentSyst[i]
            >> _vIsolationPercentSyst[i] >> _vMPTResolPercentSyst[i] 
            >> _vTrackingPercentSyst[i] >> _vMuonRecoPercentSyst[i] >> _vPtResolPercentSyst[i] >> _vTriggerPercentSyst[i];
	} //i	

        s.close();
        std::cout << "Done." << std::endl;
	
} //fillContainers()

void plotSystematicUncertainties(){
  

   TFile* outFile = new TFile("systematicErrorHistograms.root","recreate");
   std::vector<double> etaBins;
   etaBins.push_back(0.1);
   etaBins.push_back(0.35);
   etaBins.push_back(0.6);
   etaBins.push_back(0.8);
   etaBins.push_back(1.05);
   etaBins.push_back(1.3);
   etaBins.push_back(1.55);
   etaBins.push_back(1.85);
   etaBins.push_back(2.1);
   etaBins.push_back(2.4);
   const int nEtaBins = etaBins.size()-1;
   const int nCentralityBins = 6;
   const int nBins = nCentralityBins*nEtaBins;

   TString spreadSheetNameCentMuPlus = "systSpreadSheetCentMuPlus.csv";
   TString spreadSheetNameCentMuMinus = "systSpreadSheetCentMuMinus.csv";
   TString spreadSheetNameEtaMuPlus = "systSpreadSheetEtaMuPlus.csv";
   TString spreadSheetNameEtaMuMinus = "systSpreadSheetEtaMuMinus.csv";
	
   std::ofstream spreadSheetCentMuPlus ;
   std::ofstream spreadSheetCentMuMinus;
   std::ofstream spreadSheetEtaMuPlus;
   std::ofstream spreadSheetEtaMuMinus;

   spreadSheetCentMuPlus.open(spreadSheetNameCentMuPlus);
   spreadSheetCentMuMinus.open(spreadSheetNameCentMuMinus);
   spreadSheetEtaMuPlus.open(spreadSheetNameEtaMuPlus);
   spreadSheetEtaMuMinus.open(spreadSheetNameEtaMuMinus);
  
   std::vector<TString> vecSpreadCentRangeMuPlus;
   std::vector<TString> vecSpreadCentSyst1MuPlus;
   std::vector<TString> vecSpreadCentSyst2MuPlus;
   std::vector<TString> vecSpreadCentSyst3MuPlus;
   std::vector<TString> vecSpreadCentSyst4MuPlus;
   std::vector<TString> vecSpreadCentSyst5MuPlus;
   std::vector<TString> vecSpreadCentSyst6MuPlus;
   std::vector<TString> vecSpreadCentSyst7MuPlus;
   std::vector<TString> vecSpreadCentSyst8MuPlus;
   std::vector<TString> vecSpreadCentSyst9MuPlus;

   std::vector<TString> vecSpreadCentRangeMuMinus;
   std::vector<TString> vecSpreadCentSyst1MuMinus;
   std::vector<TString> vecSpreadCentSyst2MuMinus;
   std::vector<TString> vecSpreadCentSyst3MuMinus;
   std::vector<TString> vecSpreadCentSyst4MuMinus;
   std::vector<TString> vecSpreadCentSyst5MuMinus;
   std::vector<TString> vecSpreadCentSyst6MuMinus;
   std::vector<TString> vecSpreadCentSyst7MuMinus;
   std::vector<TString> vecSpreadCentSyst8MuMinus;
   std::vector<TString> vecSpreadCentSyst9MuMinus;

   std::vector<TString> vecSpreadEtaRangeMuPlus;
   std::vector<TString> vecSpreadEtaSyst1MuPlus;
   std::vector<TString> vecSpreadEtaSyst2MuPlus;
   std::vector<TString> vecSpreadEtaSyst3MuPlus;
   std::vector<TString> vecSpreadEtaSyst4MuPlus;
   std::vector<TString> vecSpreadEtaSyst5MuPlus;
   std::vector<TString> vecSpreadEtaSyst6MuPlus;
   std::vector<TString> vecSpreadEtaSyst7MuPlus;
   std::vector<TString> vecSpreadEtaSyst8MuPlus;
   std::vector<TString> vecSpreadEtaSyst9MuPlus;

   std::vector<TString> vecSpreadEtaRangeMuMinus;
   std::vector<TString> vecSpreadEtaSyst1MuMinus;
   std::vector<TString> vecSpreadEtaSyst2MuMinus;
   std::vector<TString> vecSpreadEtaSyst3MuMinus;
   std::vector<TString> vecSpreadEtaSyst4MuMinus;
   std::vector<TString> vecSpreadEtaSyst5MuMinus;
   std::vector<TString> vecSpreadEtaSyst6MuMinus;
   std::vector<TString> vecSpreadEtaSyst7MuMinus;
   std::vector<TString> vecSpreadEtaSyst8MuMinus;
   std::vector<TString> vecSpreadEtaSyst9MuMinus;

   TString sCent[] = {"0-5","5-10","10-15","15-20","20-40","40-80"};
   TString sEta[] = {"0.1-0.35","0.35-0.6","0.6-0.8","0.8-1.05","1.05-1.3","1.3-1.55","1.55-1.85","1.85-2.1","2.1-2.4"};
   TString fInPlus = "muPlusSummary.08.22.2013.txt";
   std::vector <double> _vSignalPlus(nBins),_vRaaScalingPercentSystPlus(nBins),_vEWPercentSystPlus(nBins),
                        _vOtherBkgPercentSystPlus(nBins),_vIsolationPercentSystPlus(nBins),
                       _vMPTResolPercentSystPlus(nBins),_vTrackingPercentSystPlus(nBins),
                       _vMuonRecoPercentSystPlus(nBins),_vPtResolPercentSystPlus(nBins),_vTriggerPercentSystPlus(nBins);


   fillContainers(fInPlus,_vSignalPlus,_vRaaScalingPercentSystPlus,
                    _vEWPercentSystPlus,_vOtherBkgPercentSystPlus,_vIsolationPercentSystPlus,
                    _vMPTResolPercentSystPlus,_vTrackingPercentSystPlus,
                    _vMuonRecoPercentSystPlus,_vPtResolPercentSystPlus,_vTriggerPercentSystPlus);

   TString fInMinus = "muMinusSummary.08.22.2013.txt";
   std::vector <double> _vSignalMinus(nBins),_vRaaScalingPercentSystMinus(nBins),_vEWPercentSystMinus(nBins),
                        _vOtherBkgPercentSystMinus(nBins),_vIsolationPercentSystMinus(nBins),
                       _vMPTResolPercentSystMinus(nBins),_vTrackingPercentSystMinus(nBins),
                       _vMuonRecoPercentSystMinus(nBins),_vPtResolPercentSystMinus(nBins),_vTriggerPercentSystMinus(nBins);

   fillContainers(fInMinus,_vSignalMinus,_vRaaScalingPercentSystMinus,
                    _vEWPercentSystMinus,_vOtherBkgPercentSystMinus,_vIsolationPercentSystMinus,
                    _vMPTResolPercentSystMinus,_vTrackingPercentSystMinus,
                    _vMuonRecoPercentSystMinus,_vPtResolPercentSystMinus,_vTriggerPercentSystMinus);

   double arrCentrality[] = {0.0,5.,10.,15.,20.,40.,80.};
   double arrEta[] = {0.1,0.35,0.6,0.8,1.05,1.3,1.55,1.85,2.1,2.4}; 

   ///mu+
   ///centrality
   TH1F *hRaaScalingAbsoluteSystSumSqCentPlus,*hEWAbsoluteSystSumSqCentPlus,
        *hOtherBkgAbsoluteSystSumSqCentPlus, *hIsolationAbsoluteSystSumSqCentPlus,
        *hMPTResolAbsoluteSystSumSqCentPlus,
        *hTrackingAbsoluteSystSumSqCentPlus, *hMuonRecoAbsoluteSystSumSqCentPlus,
        *hPtResolAbsoluteSystSumSqCentPlus, *hTriggerAbsoluteSystSumSqCentPlus;

   hRaaScalingAbsoluteSystSumSqCentPlus = new
    TH1F("hRaaScalingAbsoluteSystSumSqCentPlus","hRaaScalingAbsoluteSystSumSqCentPlus",nCentralityBins,arrCentrality);
   hEWAbsoluteSystSumSqCentPlus = new
    TH1F("hEWAbsoluteSystSumSqCentPlus","hEWAbsoluteSystSumSqCentPlus",nCentralityBins,arrCentrality);
   hOtherBkgAbsoluteSystSumSqCentPlus = new
    TH1F("hOtherBkgAbsoluteSystSumSqCentPlus","hOtherBkgAbsoluteSystSumSqCentPlus",nCentralityBins,arrCentrality);
   hIsolationAbsoluteSystSumSqCentPlus = new
    TH1F("hIsolationAbsoluteSystSumSqCentPlus","hIsolationAbsoluteSystSumSqCentPlus",nCentralityBins,arrCentrality);
   hMPTResolAbsoluteSystSumSqCentPlus = new
    TH1F("hMPTResolAbsoluteSystSumSqCentPlus","hMPTResolAbsoluteSystSumSqCentPlus",nCentralityBins,arrCentrality);
   hTrackingAbsoluteSystSumSqCentPlus = new
    TH1F("hTrackingAbsoluteSystSumSqCentPlus","hTrackingAbsoluteSystSumSqCentPlus",nCentralityBins,arrCentrality);
   hMuonRecoAbsoluteSystSumSqCentPlus = new
    TH1F("hMuonRecoAbsoluteSystSumSqCentPlus","hMuonRecoAbsoluteSystSumSqCentPlus",nCentralityBins,arrCentrality);
   hPtResolAbsoluteSystSumSqCentPlus = new 
    TH1F("hPtResolAbsoluteSystSumSqCentPlus","hPtResolAbsoluteSystSumSqCentPlus",nCentralityBins,arrCentrality);
   hTriggerAbsoluteSystSumSqCentPlus = new 
    TH1F("hTriggerAbsoluteSystSumSqCentPlus","hTriggerAbsoluteSystSumSqCentPlus",nCentralityBins,arrCentrality);

   ///mu+ Eta
   TH1F *hRaaScalingAbsoluteSystSumSqEtaPlus,*hEWAbsoluteSystSumSqEtaPlus,
        *hOtherBkgAbsoluteSystSumSqEtaPlus, *hIsolationAbsoluteSystSumSqEtaPlus,
        *hMPTResolAbsoluteSystSumSqEtaPlus,
        *hTrackingAbsoluteSystSumSqEtaPlus, *hMuonRecoAbsoluteSystSumSqEtaPlus,
        *hPtResolAbsoluteSystSumSqEtaPlus, *hTriggerAbsoluteSystSumSqEtaPlus;

   hRaaScalingAbsoluteSystSumSqEtaPlus = new
    TH1F("hRaaScalingAbsoluteSystSumSqEtaPlus","hRaaScalingAbsoluteSystSumSqEtaPlus",nEtaBins,arrEta);
   hEWAbsoluteSystSumSqEtaPlus = new
    TH1F("hEWAbsoluteSystSumSqEtaPlus","hEWAbsoluteSystSumSqEtaPlus",nEtaBins,arrEta);
   hOtherBkgAbsoluteSystSumSqEtaPlus = new
    TH1F("hOtherBkgAbsoluteSystSumSqEtaPlus","hOtherBkgAbsoluteSystSumSqEtaPlus",nEtaBins,arrEta);
   hIsolationAbsoluteSystSumSqEtaPlus = new
    TH1F("hIsolationAbsoluteSystSumSqEtaPlus","hIsolationAbsoluteSystSumSqEtaPlus",nEtaBins,arrEta);
   hMPTResolAbsoluteSystSumSqEtaPlus = new
    TH1F("hMPTResolAbsoluteSystSumSqEtaPlus","hMPTResolAbsoluteSystSumSqEtaPlus",nEtaBins,arrEta);
   hTrackingAbsoluteSystSumSqEtaPlus = new
    TH1F("hTrackingAbsoluteSystSumSqEtaPlus","hTrackingAbsoluteSystSumSqEtaPlus",nEtaBins,arrEta);
   hMuonRecoAbsoluteSystSumSqEtaPlus = new
    TH1F("hMuonRecoAbsoluteSystSumSqEtaPlus","hMuonRecoAbsoluteSystSumSqEtaPlus",nEtaBins,arrEta);
   hPtResolAbsoluteSystSumSqEtaPlus = new 
    TH1F("hPtResolAbsoluteSystSumSqEtaPlus","hPtResolAbsoluteSystSumSqEtaPlus",nEtaBins,arrEta);
   hTriggerAbsoluteSystSumSqEtaPlus = new 
    TH1F("hTriggerAbsoluteSystSumSqEtaPlus","hTriggerAbsoluteSystSumSqEtaPlus",nEtaBins,arrEta);

    ///mu- centrality
   TH1F *hRaaScalingAbsoluteSystSumSqCentMinus,*hEWAbsoluteSystSumSqCentMinus,
        *hOtherBkgAbsoluteSystSumSqCentMinus, *hIsolationAbsoluteSystSumSqCentMinus,
        *hMPTResolAbsoluteSystSumSqCentMinus,
        *hTrackingAbsoluteSystSumSqCentMinus, *hMuonRecoAbsoluteSystSumSqCentMinus,
        *hPtResolAbsoluteSystSumSqCentMinus, *hTriggerAbsoluteSystSumSqCentMinus;

   hRaaScalingAbsoluteSystSumSqCentMinus = new
    TH1F("hRaaScalingAbsoluteSystSumSqCentMinus","hRaaScalingAbsoluteSystSumSqCentMinus",nCentralityBins,arrCentrality);
   hEWAbsoluteSystSumSqCentMinus = new
    TH1F("hEWAbsoluteSystSumSqCentMinus","hEWAbsoluteSystSumSqCentMinus",nCentralityBins,arrCentrality);
   hOtherBkgAbsoluteSystSumSqCentMinus = new
    TH1F("hOtherBkgAbsoluteSystSumSqCentMinus","hOtherBkgAbsoluteSystSumSqCentMinus",nCentralityBins,arrCentrality);
   hIsolationAbsoluteSystSumSqCentMinus = new
    TH1F("hIsolationAbsoluteSystSumSqCentMinus","hIsolationAbsoluteSystSumSqCentMinus",nCentralityBins,arrCentrality);
   hMPTResolAbsoluteSystSumSqCentMinus = new
    TH1F("hMPTResolAbsoluteSystSumSqCentMinus","hMPTResolAbsoluteSystSumSqCentMinus",nCentralityBins,arrCentrality);
   hTrackingAbsoluteSystSumSqCentMinus = new
    TH1F("hTrackingAbsoluteSystSumSqCentMinus","hTrackingAbsoluteSystSumSqCentMinus",nCentralityBins,arrCentrality);
   hMuonRecoAbsoluteSystSumSqCentMinus = new
    TH1F("hMuonRecoAbsoluteSystSumSqCentMinus","hMuonRecoAbsoluteSystSumSqCentMinus",nCentralityBins,arrCentrality);
   hPtResolAbsoluteSystSumSqCentMinus = new 
    TH1F("hPtResolAbsoluteSystSumSqCentMinus","hPtResolAbsoluteSystSumSqCentMinus",nCentralityBins,arrCentrality);
   hTriggerAbsoluteSystSumSqCentMinus = new 
    TH1F("hTriggerAbsoluteSystSumSqCentMinus","hTriggerAbsoluteSystSumSqCentMinus",nCentralityBins,arrCentrality);

    ///mu- eta
   TH1F *hRaaScalingAbsoluteSystSumSqEtaMinus,*hEWAbsoluteSystSumSqEtaMinus,
        *hOtherBkgAbsoluteSystSumSqEtaMinus, *hIsolationAbsoluteSystSumSqEtaMinus,
        *hMPTResolAbsoluteSystSumSqEtaMinus,
        *hTrackingAbsoluteSystSumSqEtaMinus, *hMuonRecoAbsoluteSystSumSqEtaMinus,
        *hPtResolAbsoluteSystSumSqEtaMinus, *hTriggerAbsoluteSystSumSqEtaMinus;

   hRaaScalingAbsoluteSystSumSqEtaMinus = new
    TH1F("hRaaScalingAbsoluteSystSumSqEtaMinus","hRaaScalingAbsoluteSystSumSqEtaMinus",nEtaBins,arrEta);
   hEWAbsoluteSystSumSqEtaMinus = new
    TH1F("hEWAbsoluteSystSumSqEtaMinus","hEWAbsoluteSystSumSqEtaMinus",nEtaBins,arrEta);
   hOtherBkgAbsoluteSystSumSqEtaMinus = new
    TH1F("hOtherBkgAbsoluteSystSumSqEtaMinus","hOtherBkgAbsoluteSystSumSqEtaMinus",nEtaBins,arrEta);
   hIsolationAbsoluteSystSumSqEtaMinus = new
    TH1F("hIsolationAbsoluteSystSumSqEtaMinus","hIsolationAbsoluteSystSumSqEtaMinus",nEtaBins,arrEta);
   hMPTResolAbsoluteSystSumSqEtaMinus = new
    TH1F("hMPTResolAbsoluteSystSumSqEtaMinus","hMPTResolAbsoluteSystSumSqEtaMinus",nEtaBins,arrEta);
   hTrackingAbsoluteSystSumSqEtaMinus = new
    TH1F("hTrackingAbsoluteSystSumSqEtaMinus","hTrackingAbsoluteSystSumSqEtaMinus",nEtaBins,arrEta);
   hMuonRecoAbsoluteSystSumSqEtaMinus = new
    TH1F("hMuonRecoAbsoluteSystSumSqEtaMinus","hMuonRecoAbsoluteSystSumSqEtaMinus",nEtaBins,arrEta);
   hPtResolAbsoluteSystSumSqEtaMinus = new 
    TH1F("hPtResolAbsoluteSystSumSqEtaMinus","hPtResolAbsoluteSystSumSqEtaMinus",nEtaBins,arrEta);
   hTriggerAbsoluteSystSumSqEtaMinus = new 
    TH1F("hTriggerAbsoluteSystSumSqEtaMinus","hTriggerAbsoluteSystSumSqEtaMinus",nEtaBins,arrEta);

   
   ///Sum up over all eta, differentially in centrality
   for(int icent=0; icent<nCentralityBins; ++icent){

    double _nSignalCentPlus=0.0,
            nRaaScalingAbsoluteSystSumSqCentPlus=0.0,
            nEWAbsoluteSystSumSqCentPlus=0.0,
            nOtherBkgAbsoluteSystSumSqCentPlus=0.0,
            nIsolationAbsoluteSystSumSqCentPlus=0.0,
            nMPTResolAbsoluteSystSumSqCentPlus=0.,nTrackingAbsoluteSystSumSqCentPlus=0.,
            nMuonRecoAbsoluteSystSumSqCentPlus=0.,nPtResolAbsoluteSystSumSqCentPlus=0.,nTriggerAbsoluteSystSumSqCentPlus=0.;
    double _nSignalCentMinus=0.0,
            nRaaScalingAbsoluteSystSumSqCentMinus=0.0,
            nEWAbsoluteSystSumSqCentMinus=0.0,
            nOtherBkgAbsoluteSystSumSqCentMinus=0.0,
            nIsolationAbsoluteSystSumSqCentMinus=0.0,
            nMPTResolAbsoluteSystSumSqCentMinus=0.,nTrackingAbsoluteSystSumSqCentMinus=0.,
            nMuonRecoAbsoluteSystSumSqCentMinus=0.,nPtResolAbsoluteSystSumSqCentMinus=0.,nTriggerAbsoluteSystSumSqCentMinus=0.;

    for(int ieta=0; ieta<nEtaBins; ++ieta){

        int index = ieta*nCentralityBins+icent;

        ///Sum up over all nCentralityBins for eta bin ieta
        _nSignalCentPlus += _vSignalPlus[index];

        ///Running sum of squares of absolute error on the signal

        ///mu+
        ///systematic on background subtraction
        nRaaScalingAbsoluteSystSumSqCentPlus +=
            TMath::Power(_vRaaScalingPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

        nEWAbsoluteSystSumSqCentPlus +=
            TMath::Power(_vEWPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

        nOtherBkgAbsoluteSystSumSqCentPlus +=
            TMath::Power(_vOtherBkgPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

        ///systematic on isolation cut 
        nIsolationAbsoluteSystSumSqCentPlus +=
            TMath::Power(_vIsolationPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

        ///systematic on isolation cut 
        nMPTResolAbsoluteSystSumSqCentPlus +=
            TMath::Power(_vMPTResolPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

       
        ///systematic on tracking  
        // removed tracking syst
        nTrackingAbsoluteSystSumSqCentPlus +=
            TMath::Power(0.0*_vTrackingPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

        ///systematic on muon reco  
        nMuonRecoAbsoluteSystSumSqCentPlus +=
            TMath::Power(_vMuonRecoPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

        ///systematic on muon pT resolution 
        nPtResolAbsoluteSystSumSqCentPlus +=
            TMath::Power(_vPtResolPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

        ///systematic on trigger efficiency 
        nTriggerAbsoluteSystSumSqCentPlus +=
            TMath::Power(_vTriggerPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

        ///Ditto for mu-
        ///Sum up over all nCentralityBins for eta bin ieta
        _nSignalCentMinus += _vSignalMinus[index];

        ///Running sum of squares of absolute error on the signal

        ///mu+
        ///systematic on background subtraction
        nRaaScalingAbsoluteSystSumSqCentMinus +=
            TMath::Power(_vRaaScalingPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

        nEWAbsoluteSystSumSqCentMinus +=
            TMath::Power(_vEWPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

        nOtherBkgAbsoluteSystSumSqCentMinus +=
            TMath::Power(_vOtherBkgPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

        ///systematic on isolation cut 
        nIsolationAbsoluteSystSumSqCentMinus +=
            TMath::Power(_vIsolationPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

        ///systematic on isolation cut 
        nMPTResolAbsoluteSystSumSqCentMinus +=
            TMath::Power(_vMPTResolPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

       
        ///systematic on tracking  
        // removed tracking syst
        nTrackingAbsoluteSystSumSqCentMinus +=
            TMath::Power(0.0*_vTrackingPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

        ///systematic on muon reco  
        nMuonRecoAbsoluteSystSumSqCentMinus +=
            TMath::Power(_vMuonRecoPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

        ///systematic on muon pT resolution 
        nPtResolAbsoluteSystSumSqCentMinus +=
            TMath::Power(_vPtResolPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

        ///systematic on trigger efficiency 
        nTriggerAbsoluteSystSumSqCentMinus +=
            TMath::Power(_vTriggerPercentSystMinus[index]/100.0*_vSignalMinus[index],2);


    } //ieta

    std::cout << "Number of W^{+} in centrality bin " << icent << " : " << _nSignalCentPlus << std::endl;
    std::cout << "Number of W^{-} in centrality bin " << icent << " : " << _nSignalCentMinus << std::endl;

    double percentSyst1 = 100.0*TMath::Sqrt(nRaaScalingAbsoluteSystSumSqCentPlus)/_nSignalCentPlus;
    double percentSyst2 = 100.0*TMath::Sqrt(nEWAbsoluteSystSumSqCentPlus)/_nSignalCentPlus;
    double percentSyst3 = 100.0*TMath::Sqrt(nOtherBkgAbsoluteSystSumSqCentPlus)/_nSignalCentPlus;
    double percentSyst4 = 100.0*TMath::Sqrt(nIsolationAbsoluteSystSumSqCentPlus)/_nSignalCentPlus;
    double percentSyst5 = 100.0*TMath::Sqrt(nMPTResolAbsoluteSystSumSqCentPlus)/_nSignalCentPlus;
    double percentSyst6 = 100.0*TMath::Sqrt(nTrackingAbsoluteSystSumSqCentPlus)/_nSignalCentPlus;
    double percentSyst7 = 100.0*TMath::Sqrt(nMuonRecoAbsoluteSystSumSqCentPlus)/_nSignalCentPlus;
    double percentSyst8 = 100.0*TMath::Sqrt(nPtResolAbsoluteSystSumSqCentPlus)/_nSignalCentPlus;
    double percentSyst9 = 100.0*TMath::Sqrt(nTriggerAbsoluteSystSumSqCentPlus)/_nSignalCentPlus;

    std::cout << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "W^{+}  S u m m a r y  ( C E N T R A L I T Y )" << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "QCD                                    = " <<
    100.0*TMath::Sqrt(nRaaScalingAbsoluteSystSumSqCentPlus)/_nSignalCentPlus<< std::endl;
    std::cout << "EW                                     = " <<
    100.0*TMath::Sqrt(nEWAbsoluteSystSumSqCentPlus)/_nSignalCentPlus<< std::endl;
    std::cout << "Other background                       = " <<
    100.0*TMath::Sqrt(nOtherBkgAbsoluteSystSumSqCentPlus)/_nSignalCentPlus<< std::endl;
    std::cout << "Isolation                              = " <<
    100.0*TMath::Sqrt(nIsolationAbsoluteSystSumSqCentPlus)/_nSignalCentPlus << std::endl;
    std::cout << "MpT resolution                         = " <<
    100.0*TMath::Sqrt(nMPTResolAbsoluteSystSumSqCentPlus)/_nSignalCentPlus << std::endl; 
    std::cout << "Tracking                               = "
    <<100.0*TMath::Sqrt(nTrackingAbsoluteSystSumSqCentPlus)/_nSignalCentPlus  << std::endl;
    std::cout << "Muon reconstruction                    = " <<
    100.0*TMath::Sqrt(nMuonRecoAbsoluteSystSumSqCentPlus)/_nSignalCentPlus  << std::endl;
    std::cout << "pT resolution                          = "
    <<100.0*TMath::Sqrt(nPtResolAbsoluteSystSumSqCentPlus)/_nSignalCentPlus << std::endl;
    std::cout << "Trigger                                = "   
    <<100.0*TMath::Sqrt(nTriggerAbsoluteSystSumSqCentPlus)/_nSignalCentPlus  << std::endl;
    std::cout << std::endl;


    ///Relative errors (%) in this centrality bin icent
    double errP1 = 100.0*TMath::Sqrt(nRaaScalingAbsoluteSystSumSqCentPlus)/_nSignalCentPlus;
    double errP2 =
        TMath::Sqrt(TMath::Power(errP1,2)+TMath::Power(100.0*TMath::Sqrt(nEWAbsoluteSystSumSqCentPlus)/_nSignalCentPlus,2));
    double errP3 =
        TMath::Sqrt(TMath::Power(errP2,2)+TMath::Power(100.0*TMath::Sqrt(nOtherBkgAbsoluteSystSumSqCentPlus)/_nSignalCentPlus,2));
    double errP4 =
        TMath::Sqrt(TMath::Power(errP3,2)+TMath::Power(100.0*TMath::Sqrt(nIsolationAbsoluteSystSumSqCentPlus)/_nSignalCentPlus,2));
    double errP5 = 
        TMath::Sqrt(TMath::Power(errP4,2)+TMath::Power(100.0*TMath::Sqrt(nMPTResolAbsoluteSystSumSqCentPlus)/_nSignalCentPlus ,2));
    double errP6 = 
        TMath::Sqrt(TMath::Power(errP5,2)+TMath::Power(100.0*TMath::Sqrt(nTrackingAbsoluteSystSumSqCentPlus)/_nSignalCentPlus ,2));
    double errP7 = 
        TMath::Sqrt(TMath::Power(errP6,2)+TMath::Power(100.0*TMath::Sqrt(nMuonRecoAbsoluteSystSumSqCentPlus)/_nSignalCentPlus ,2));
    double errP8 = 
        TMath::Sqrt(TMath::Power(errP7,2)+TMath::Power(100.0*TMath::Sqrt(nPtResolAbsoluteSystSumSqCentPlus)/_nSignalCentPlus ,2));
    double errP9 = 
        TMath::Sqrt(TMath::Power(errP8,2)+TMath::Power(100.0*TMath::Sqrt(nTriggerAbsoluteSystSumSqCentPlus)/_nSignalCentPlus ,2));

    // Save to csv file
    /*writeToSpreadSheet(vecSpreadSheetCentMuPlus,sCent[icent],format(percentSyst1,2),format(percentSyst2,2),
				format(percentSyst3,2),format(percentSyst4,2),format(percentSyst5,2),
				format(percentSyst6,2),format(percentSyst7,2),format(percentSyst8,2),format(percentSyst9,2));
*/

    TString sCentTemp = sCent[icent];
    vecSpreadCentRangeMuPlus.push_back(sCentTemp);
    vecSpreadCentSyst1MuPlus.push_back(format(percentSyst1,1));
    vecSpreadCentSyst2MuPlus.push_back(format(percentSyst2,1));
    vecSpreadCentSyst3MuPlus.push_back(format(percentSyst3,1));
    vecSpreadCentSyst4MuPlus.push_back(format(percentSyst4,1));
    vecSpreadCentSyst5MuPlus.push_back(format(percentSyst5,1));
    vecSpreadCentSyst6MuPlus.push_back(format(percentSyst6,1));
    vecSpreadCentSyst7MuPlus.push_back(format(percentSyst7,1));
    vecSpreadCentSyst8MuPlus.push_back(format(percentSyst8,1));
    vecSpreadCentSyst9MuPlus.push_back(format(percentSyst9,1));

    hRaaScalingAbsoluteSystSumSqCentPlus ->SetBinError(icent+1, errP1);
    hEWAbsoluteSystSumSqCentPlus ->SetBinError(icent+1, errP2);
    hOtherBkgAbsoluteSystSumSqCentPlus ->SetBinError(icent+1, errP3);
    hIsolationAbsoluteSystSumSqCentPlus ->SetBinError(icent+1, errP4);
    hMPTResolAbsoluteSystSumSqCentPlus ->SetBinError(icent+1, errP5);
    hTrackingAbsoluteSystSumSqCentPlus ->SetBinError(icent+1, errP6);
    hMuonRecoAbsoluteSystSumSqCentPlus ->SetBinError(icent+1, errP7);
    hPtResolAbsoluteSystSumSqCentPlus ->SetBinError(icent+1, errP8);
    hTriggerAbsoluteSystSumSqCentPlus ->SetBinError(icent+1, errP9);


    std::cout << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "W^{-}  S u m m a r y  ( C E N T R A L I T Y )" << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "QCD                                    = " <<
    100.0*TMath::Sqrt(nRaaScalingAbsoluteSystSumSqCentMinus)/_nSignalCentMinus<< std::endl;
    std::cout << "EW                                     = " <<
    100.0*TMath::Sqrt(nEWAbsoluteSystSumSqCentMinus)/_nSignalCentMinus<< std::endl;
    std::cout << "Other background                       = " <<
    100.0*TMath::Sqrt(nOtherBkgAbsoluteSystSumSqCentMinus)/_nSignalCentMinus<< std::endl;
    std::cout << "Isolation                              = " <<
    100.0*TMath::Sqrt(nIsolationAbsoluteSystSumSqCentMinus)/_nSignalCentMinus << std::endl;
    std::cout << "MpT resolution                         = " <<
    100.0*TMath::Sqrt(nMPTResolAbsoluteSystSumSqCentMinus)/_nSignalCentMinus << std::endl; 
    std::cout << "Tracking                               = "
    <<100.0*TMath::Sqrt(nTrackingAbsoluteSystSumSqCentMinus)/_nSignalCentMinus  << std::endl;
    std::cout << "Muon reconstruction                    = " <<
    100.0*TMath::Sqrt(nMuonRecoAbsoluteSystSumSqCentMinus)/_nSignalCentMinus  << std::endl;
    std::cout << "pT resolution                          = "
    <<100.0*TMath::Sqrt(nPtResolAbsoluteSystSumSqCentMinus)/_nSignalCentMinus << std::endl;
    std::cout << "Trigger                                = "   
    <<100.0*TMath::Sqrt(nTriggerAbsoluteSystSumSqCentMinus)/_nSignalCentMinus  << std::endl;
    std::cout << std::endl;


    ///Relative errors (%) in this centrality bin icent
    double errM1 = 100.0*TMath::Sqrt(nRaaScalingAbsoluteSystSumSqCentMinus)/_nSignalCentMinus;
    double errM2 =
        TMath::Sqrt(TMath::Power(errM1,2)+TMath::Power(100.0*TMath::Sqrt(nEWAbsoluteSystSumSqCentMinus)/_nSignalCentMinus,2));
    double errM3 =
        TMath::Sqrt(TMath::Power(errM2,2)+TMath::Power(100.0*TMath::Sqrt(nOtherBkgAbsoluteSystSumSqCentMinus)/_nSignalCentMinus,2));
    double errM4 =
        TMath::Sqrt(TMath::Power(errM3,2)+TMath::Power(100.0*TMath::Sqrt(nIsolationAbsoluteSystSumSqCentMinus)/_nSignalCentMinus,2));
    double errM5 = 
        TMath::Sqrt(TMath::Power(errM4,2)+TMath::Power(100.0*TMath::Sqrt(nMPTResolAbsoluteSystSumSqCentMinus)/_nSignalCentMinus ,2));
    double errM6 = 
        TMath::Sqrt(TMath::Power(errM5,2)+TMath::Power(100.0*TMath::Sqrt(nTrackingAbsoluteSystSumSqCentMinus)/_nSignalCentMinus ,2));
    double errM7 = 
        TMath::Sqrt(TMath::Power(errM6,2)+TMath::Power(100.0*TMath::Sqrt(nMuonRecoAbsoluteSystSumSqCentMinus)/_nSignalCentMinus ,2));
    double errM8 = 
        TMath::Sqrt(TMath::Power(errM7,2)+TMath::Power(100.0*TMath::Sqrt(nPtResolAbsoluteSystSumSqCentMinus)/_nSignalCentMinus ,2));
    double errM9 = 
        TMath::Sqrt(TMath::Power(errM8,2)+TMath::Power(100.0*TMath::Sqrt(nTriggerAbsoluteSystSumSqCentMinus)/_nSignalCentMinus ,2));

    percentSyst1 = 100.0*TMath::Sqrt(nRaaScalingAbsoluteSystSumSqCentMinus)/_nSignalCentMinus;
    percentSyst2 = 100.0*TMath::Sqrt(nEWAbsoluteSystSumSqCentMinus)/_nSignalCentMinus;
    percentSyst3 = 100.0*TMath::Sqrt(nOtherBkgAbsoluteSystSumSqCentMinus)/_nSignalCentMinus;
    percentSyst4 = 100.0*TMath::Sqrt(nIsolationAbsoluteSystSumSqCentMinus)/_nSignalCentMinus;
    percentSyst5 = 100.0*TMath::Sqrt(nMPTResolAbsoluteSystSumSqCentMinus)/_nSignalCentMinus;
    percentSyst6 = 100.0*TMath::Sqrt(nTrackingAbsoluteSystSumSqCentMinus)/_nSignalCentMinus;
    percentSyst7 = 100.0*TMath::Sqrt(nMuonRecoAbsoluteSystSumSqCentMinus)/_nSignalCentMinus;
    percentSyst8 = 100.0*TMath::Sqrt(nPtResolAbsoluteSystSumSqCentMinus)/_nSignalCentMinus;
	std::cout << "percentSyst8: " << percentSyst8 << std::endl;
    percentSyst9 = 100.0*TMath::Sqrt(nTriggerAbsoluteSystSumSqCentMinus)/_nSignalCentMinus;
    // Save to csv file

    /*writeToSpreadSheet(vecSpreadSheetCentMuMinus,sCent[icent],format(percentSyst1,2),format(percentSyst2,2),
				format(percentSyst3,2),format(percentSyst4,2),format(percentSyst5,2),
				format(percentSyst6,2),format(percentSyst7,2),format(percentSyst8,2),format(percentSyst9,2));
*/
    vecSpreadCentRangeMuMinus.push_back(sCent[icent]);
    vecSpreadCentSyst1MuMinus.push_back(format(percentSyst1,1));
    vecSpreadCentSyst2MuMinus.push_back(format(percentSyst2,1));
    vecSpreadCentSyst3MuMinus.push_back(format(percentSyst3,1));
    vecSpreadCentSyst4MuMinus.push_back(format(percentSyst4,1));
    vecSpreadCentSyst5MuMinus.push_back(format(percentSyst5,1));
    vecSpreadCentSyst6MuMinus.push_back(format(percentSyst6,1));
    vecSpreadCentSyst7MuMinus.push_back(format(percentSyst7,1));
    vecSpreadCentSyst8MuMinus.push_back(format(percentSyst8,1));
    vecSpreadCentSyst9MuMinus.push_back(format(percentSyst9,1));

    hRaaScalingAbsoluteSystSumSqCentMinus ->SetBinError(icent+1, errM1);
    hEWAbsoluteSystSumSqCentMinus ->SetBinError(icent+1, errM2);
    hOtherBkgAbsoluteSystSumSqCentMinus ->SetBinError(icent+1, errM3);
    hIsolationAbsoluteSystSumSqCentMinus ->SetBinError(icent+1, errM4);
    hMPTResolAbsoluteSystSumSqCentMinus ->SetBinError(icent+1, errM5);
    hTrackingAbsoluteSystSumSqCentMinus ->SetBinError(icent+1, errM6);
    hMuonRecoAbsoluteSystSumSqCentMinus ->SetBinError(icent+1, errM7);
    hPtResolAbsoluteSystSumSqCentMinus ->SetBinError(icent+1, errM8);
    hTriggerAbsoluteSystSumSqCentMinus ->SetBinError(icent+1, errM9);

   } //icent


   ///Sum up over all centrality, differentially in eta
   char cSyst1[50],cSyst2[50],cSyst3[50],cSyst4[50],cSyst5[50],cSyst6[50],cSyst7[50];
/*   TH2D *h2D_RaaScalingPercentSyst[nEtaBins],*h2D_IsolationPercentSyst[nEtaBins],*h2D_MPTResolPercentSyst[nEtaBins];
   TH2D *h2D_TrackingPercentSyst[nEtaBins],*h2D_MuonRecoPercentSyst[nEtaBins],*h2D_PtResolPercentSyst[nEtaBins],*h2D_TriggerPercentSyst[nEtaBins];
   */

/*   TH2D *h2D_RaaScalingPercentSyst,*h2D_IsolationPercentSyst,*h2D_MPTResolPercentSyst;
   TH2D *h2D_TrackingPercentSyst,*h2D_MuonRecoPercentSyst,*h2D_PtResolPercentSyst,*h2D_TriggerPercentSyst;
   sprintf(cSyst1,"chargeCorrelationMatrixRaaScalingEta%i",0);
   h2D_RaaScalingPercentSyst = new TH2D(cSyst1,cSyst1,300,0.0,30.0,300,0.0,30.0); 
   sprintf(cSyst2,"chargeCorrelationMatrixIsolationEta%i",0);
   h2D_IsolationPercentSyst = new TH2D(cSyst2,cSyst2,300,0.0,30.0,300,0.0,30.0); 
   sprintf(cSyst3,"chargeCorrelationMatrixMPTResolEta%i",0);
   h2D_MPTResolPercentSyst = new TH2D(cSyst3,cSyst3,300,0.0,30.0,300,0.0,30.0); 
   sprintf(cSyst4,"chargeCorrelationMatrixTrackingEta%i",0);
   h2D_TrackingPercentSyst = new TH2D(cSyst4,cSyst4,300,0.0,30.0,300,0.0,30.0); 
   sprintf(cSyst5,"chargeCorrelationMatrixMuonRecoEta%i",0);
   h2D_MuonRecoPercentSyst = new TH2D(cSyst5,cSyst5,300,0.0,30.0,300,0.0,30.0); 
   sprintf(cSyst6,"chargeCorrelationMatrixPtResolEta%i",0);
   h2D_PtResolPercentSyst = new TH2D(cSyst6,cSyst6,300,0.0,30.0,300,0.0,30.0); 
   sprintf(cSyst7,"chargeCorrelationMatrixTriggerEta%i",0);
   h2D_TriggerPercentSyst = new TH2D(cSyst7,cSyst7,300,0.0,30.0,300,0.0,30.0); 

   TGraph *_grRaaScalingCorrelationFactor = new TGraph(nEtaBins);
   TGraph *_grIsolationCorrelationFactor = new TGraph(nEtaBins);
   TGraph *_grMPTResolCorrelationFactor = new TGraph(nEtaBins);
   TGraph *_grTrackingCorrelationFactor = new TGraph(nEtaBins);
   TGraph *_grMuonRecoCorrelationFactor = new TGraph(nEtaBins);
   TGraph *_grPtResolCorrelationFactor = new TGraph(nEtaBins);
   TGraph *_grTriggerCorrelationFactor = new TGraph(nEtaBins);
*/
   for(int ieta=0; ieta<nEtaBins; ++ieta){

    double _nSignalEtaPlus=0.0,
            nRaaScalingAbsoluteSystSumSqEtaPlus=0.0,
            nEWAbsoluteSystSumSqEtaPlus=0.0,
            nOtherBkgAbsoluteSystSumSqEtaPlus=0.0,
            nIsolationAbsoluteSystSumSqEtaPlus=0.0,
            nMPTResolAbsoluteSystSumSqEtaPlus=0.,nTrackingAbsoluteSystSumSqEtaPlus=0.,
            nMuonRecoAbsoluteSystSumSqEtaPlus=0.,nPtResolAbsoluteSystSumSqEtaPlus=0.,nTriggerAbsoluteSystSumSqEtaPlus=0.;
    double _nSignalEtaMinus=0.0,
            nRaaScalingAbsoluteSystSumSqEtaMinus=0.0,
            nEWAbsoluteSystSumSqEtaMinus=0.0,
            nOtherBkgAbsoluteSystSumSqEtaMinus=0.0,
            nIsolationAbsoluteSystSumSqEtaMinus=0.0,
            nMPTResolAbsoluteSystSumSqEtaMinus=0.,nTrackingAbsoluteSystSumSqEtaMinus=0.,
            nMuonRecoAbsoluteSystSumSqEtaMinus=0.,nPtResolAbsoluteSystSumSqEtaMinus=0.,nTriggerAbsoluteSystSumSqEtaMinus=0.;

/*    sprintf(cSyst1,"chargeCorrelationMatrixRaaScalingEta%i",ieta);
    h2D_RaaScalingPercentSyst[ieta] = new TH2D(cSyst1,cSyst1,300,0.0,30.0,300,0.0,30.0); 
    sprintf(cSyst2,"chargeCorrelationMatrixIsolationEta%i",ieta);
    h2D_IsolationPercentSyst[ieta] = new TH2D(cSyst2,cSyst2,300,0.0,30.0,300,0.0,30.0); 
    sprintf(cSyst3,"chargeCorrelationMatrixMPTResolEta%i",ieta);
    h2D_MPTResolPercentSyst[ieta] = new TH2D(cSyst3,cSyst3,300,0.0,30.0,300,0.0,30.0); 
    sprintf(cSyst4,"chargeCorrelationMatrixTrackingEta%i",ieta);
    h2D_TrackingPercentSyst[ieta] = new TH2D(cSyst4,cSyst4,300,0.0,30.0,300,0.0,30.0); 
    sprintf(cSyst5,"chargeCorrelationMatrixMuonRecoEta%i",ieta);
    h2D_MuonRecoPercentSyst[ieta] = new TH2D(cSyst5,cSyst5,300,0.0,30.0,300,0.0,30.0); 
    sprintf(cSyst6,"chargeCorrelationMatrixPtResolEta%i",ieta);
    h2D_PtResolPercentSyst[ieta] = new TH2D(cSyst6,cSyst6,300,0.0,30.0,300,0.0,30.0); 
    sprintf(cSyst7,"chargeCorrelationMatrixTriggerEta%i",ieta);
    h2D_TriggerPercentSyst[ieta] = new TH2D(cSyst7,cSyst7,300,0.0,30.0,300,0.0,30.0); 
*/
    for(int icent=0; icent<nCentralityBins; ++icent){

        int index = ieta*nCentralityBins+icent;

        ///Sum up over all nEtaBins for centrality bin icent
        _nSignalEtaPlus += _vSignalPlus[index];

        ///Running sum of squares of absolute error on the signal

        ///mu+
        ///systematic on background subtraction
        nRaaScalingAbsoluteSystSumSqEtaPlus +=
            TMath::Power(_vRaaScalingPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

        nEWAbsoluteSystSumSqEtaPlus +=
            TMath::Power(_vEWPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

        nOtherBkgAbsoluteSystSumSqEtaPlus +=
            TMath::Power(_vOtherBkgPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

        ///systematic on isolation cut 
        nIsolationAbsoluteSystSumSqEtaPlus +=
            TMath::Power(_vIsolationPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

        ///systematic on isolation cut 
        nMPTResolAbsoluteSystSumSqEtaPlus +=
            TMath::Power(_vMPTResolPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

       
        ///systematic on tracking  
        // removed tracking syst
        nTrackingAbsoluteSystSumSqEtaPlus +=
            TMath::Power(0.0*_vTrackingPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

        ///systematic on muon reco  
        nMuonRecoAbsoluteSystSumSqEtaPlus +=
            TMath::Power(_vMuonRecoPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

        ///systematic on muon pT resolution 
        nPtResolAbsoluteSystSumSqEtaPlus +=
            TMath::Power(_vPtResolPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

        ///systematic on trigger efficiency 
        nTriggerAbsoluteSystSumSqEtaPlus +=
            TMath::Power(_vTriggerPercentSystPlus[index]/100.0*_vSignalPlus[index],2);

        ///Ditto for mu-
        ///Sum up over all nEtaralityBins for eta bin ieta
        _nSignalEtaMinus += _vSignalMinus[index];

        ///Running sum of squares of absolute error on the signal

        ///mu+
        ///systematic on background subtraction
        nRaaScalingAbsoluteSystSumSqEtaMinus +=
            TMath::Power(_vRaaScalingPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

        nEWAbsoluteSystSumSqEtaMinus +=
            TMath::Power(_vEWPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

        nOtherBkgAbsoluteSystSumSqEtaMinus +=
            TMath::Power(_vOtherBkgPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

        ///systematic on isolation cut 
        nIsolationAbsoluteSystSumSqEtaMinus +=
            TMath::Power(_vIsolationPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

        ///systematic on isolation cut 
        nMPTResolAbsoluteSystSumSqEtaMinus +=
            TMath::Power(_vMPTResolPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

       
        ///systematic on tracking  
        // removed tracking syst
        nTrackingAbsoluteSystSumSqEtaMinus +=
            TMath::Power(0.0*_vTrackingPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

        ///systematic on muon reco  
        nMuonRecoAbsoluteSystSumSqEtaMinus +=
            TMath::Power(_vMuonRecoPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

        ///systematic on muon pT resolution 
        nPtResolAbsoluteSystSumSqEtaMinus +=
            TMath::Power(_vPtResolPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

        ///systematic on trigger efficiency 
        nTriggerAbsoluteSystSumSqEtaMinus +=
            TMath::Power(_vTriggerPercentSystMinus[index]/100.0*_vSignalMinus[index],2);

    } //icent

    ///TGraph of correlation coefficients 
    double xEta = etaBins[ieta]+(etaBins[ieta+1]-etaBins[ieta])/2.0;

    std::cout << "Number of W^{+} in eta bin " << ieta << " : " << _nSignalEtaPlus << std::endl;
    ///Relative errors (%) in this eta bin ieta
    //

    double errP1 = 100.0*TMath::Sqrt(nRaaScalingAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus;
    double errP2 =
        TMath::Sqrt(TMath::Power(errP1,2)+TMath::Power(100.0*TMath::Sqrt(nEWAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus,2));
    double errP3 =
        TMath::Sqrt(TMath::Power(errP2,2)+TMath::Power(100.0*TMath::Sqrt(nOtherBkgAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus,2));
    double errP4 =
        TMath::Sqrt(TMath::Power(errP3,2)+TMath::Power(100.0*TMath::Sqrt(nIsolationAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus,2));
    double errP5 = 
        TMath::Sqrt(TMath::Power(errP4,2)+TMath::Power(100.0*TMath::Sqrt(nMPTResolAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus ,2));
    double errP6 = 
        TMath::Sqrt(TMath::Power(errP5,2)+TMath::Power(100.0*TMath::Sqrt(nTrackingAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus ,2));
    double errP7 = 
        TMath::Sqrt(TMath::Power(errP6,2)+TMath::Power(100.0*TMath::Sqrt(nMuonRecoAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus ,2));
    double errP8 = 
        TMath::Sqrt(TMath::Power(errP7,2)+TMath::Power(100.0*TMath::Sqrt(nPtResolAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus ,2));
    double errP9 = 
        TMath::Sqrt(TMath::Power(errP8,2)+TMath::Power(100.0*TMath::Sqrt(nTriggerAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus ,2));

    hRaaScalingAbsoluteSystSumSqEtaPlus ->SetBinError(ieta+1, errP1);
    hEWAbsoluteSystSumSqEtaPlus ->SetBinError(ieta+1, errP2);
    hOtherBkgAbsoluteSystSumSqEtaPlus ->SetBinError(ieta+1, errP3);
    hIsolationAbsoluteSystSumSqEtaPlus ->SetBinError(ieta+1, errP4);
    hMPTResolAbsoluteSystSumSqEtaPlus ->SetBinError(ieta+1, errP5);
    hTrackingAbsoluteSystSumSqEtaPlus ->SetBinError(ieta+1, errP6);
    hMuonRecoAbsoluteSystSumSqEtaPlus ->SetBinError(ieta+1, errP7);
    hPtResolAbsoluteSystSumSqEtaPlus ->SetBinError(ieta+1, errP8);
    hTriggerAbsoluteSystSumSqEtaPlus ->SetBinError(ieta+1, errP9);

    std::cout << "Number of W^{-} in eta bin " << ieta << " : " << _nSignalEtaMinus << std::endl;
    ///Relative errors (%) in this eta bin ieta
    double percentSyst1 = 100.0*TMath::Sqrt(nRaaScalingAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus;
    double percentSyst2 = 100.0*TMath::Sqrt(nEWAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus;
    double percentSyst3 = 100.0*TMath::Sqrt(nOtherBkgAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus;
    double percentSyst4 = 100.0*TMath::Sqrt(nIsolationAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus;
    double percentSyst5 = 100.0*TMath::Sqrt(nMPTResolAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus;
    double percentSyst6 = 100.0*TMath::Sqrt(nTrackingAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus;
    double percentSyst7 = 100.0*TMath::Sqrt(nMuonRecoAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus;
    double percentSyst8 = 100.0*TMath::Sqrt(nPtResolAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus;
    double percentSyst9 = 100.0*TMath::Sqrt(nTriggerAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus;
    // Save to csv file

    /*writeToSpreadSheet(vecSpreadSheetEtaMuPlus,sEta[ieta],format(percentSyst1,2),format(percentSyst2,2),
				format(percentSyst3,2),format(percentSyst4,2),format(percentSyst5,2),
				format(percentSyst6,2),format(percentSyst7,2),format(percentSyst8,2),format(percentSyst9,2));

*/
    vecSpreadEtaRangeMuPlus.push_back(sEta[ieta]);
    vecSpreadEtaSyst1MuPlus.push_back(format(percentSyst1,1));
    vecSpreadEtaSyst2MuPlus.push_back(format(percentSyst2,1));
    vecSpreadEtaSyst3MuPlus.push_back(format(percentSyst3,1));
    vecSpreadEtaSyst4MuPlus.push_back(format(percentSyst4,1));
    vecSpreadEtaSyst5MuPlus.push_back(format(percentSyst5,1));
    vecSpreadEtaSyst6MuPlus.push_back(format(percentSyst6,1));
    vecSpreadEtaSyst7MuPlus.push_back(format(percentSyst7,1));
    vecSpreadEtaSyst8MuPlus.push_back(format(percentSyst8,1));
    vecSpreadEtaSyst9MuPlus.push_back(format(percentSyst9,1));

    std::cout << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "W^{+}  S u m m a r y  ( E T A )" << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "QCD                                    = " <<
    100.0*TMath::Sqrt(nRaaScalingAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus<< std::endl;
    std::cout << "EW                                     = " <<
    100.0*TMath::Sqrt(nEWAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus<< std::endl;
    std::cout << "Other background                       = " <<
    100.0*TMath::Sqrt(nOtherBkgAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus<< std::endl;
    std::cout << "Isolation                              = " <<
    100.0*TMath::Sqrt(nIsolationAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus << std::endl;
    std::cout << "MpT resolution                         = " <<
    100.0*TMath::Sqrt(nMPTResolAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus << std::endl; 
    std::cout << "Tracking                               = "
    <<100.0*TMath::Sqrt(nTrackingAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus  << std::endl;
    std::cout << "Muon reconstruction                    = " <<
    100.0*TMath::Sqrt(nMuonRecoAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus  << std::endl;
    std::cout << "pT resolution                          = "
    <<100.0*TMath::Sqrt(nPtResolAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus << std::endl;
    std::cout << "Trigger                                = "   
    <<100.0*TMath::Sqrt(nTriggerAbsoluteSystSumSqEtaPlus)/_nSignalEtaPlus  << std::endl;
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "W^{-}  S u m m a r y  ( E T A )" << std::endl;
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "QCD                                    = " <<
    100.0*TMath::Sqrt(nRaaScalingAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus<< std::endl;
    std::cout << "EW                                     = " <<
    100.0*TMath::Sqrt(nEWAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus<< std::endl;
    std::cout << "Other background                       = " <<
    100.0*TMath::Sqrt(nOtherBkgAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus<< std::endl;
    std::cout << "Isolation                              = " <<
    100.0*TMath::Sqrt(nIsolationAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus << std::endl;
    std::cout << "MpT resolution                         = " <<
    100.0*TMath::Sqrt(nMPTResolAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus << std::endl; 
    std::cout << "Tracking                               = "
    <<100.0*TMath::Sqrt(nTrackingAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus  << std::endl;
    std::cout << "Muon reconstruction                    = " <<
    100.0*TMath::Sqrt(nMuonRecoAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus  << std::endl;
    std::cout << "pT resolution                          = "
    <<100.0*TMath::Sqrt(nPtResolAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus << std::endl;
    std::cout << "Trigger                                = "   
    <<100.0*TMath::Sqrt(nTriggerAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus  << std::endl;
    std::cout << std::endl;

    double errM1 = 100.0*TMath::Sqrt(nRaaScalingAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus;
    double errM2 =
        TMath::Sqrt(TMath::Power(errM1,2)+TMath::Power(100.0*TMath::Sqrt(nEWAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus,2));
    double errM3 =
        TMath::Sqrt(TMath::Power(errM2,2)+TMath::Power(100.0*TMath::Sqrt(nOtherBkgAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus,2));
    double errM4 =
        TMath::Sqrt(TMath::Power(errM3,2)+TMath::Power(100.0*TMath::Sqrt(nIsolationAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus,2));
    double errM5 = 
        TMath::Sqrt(TMath::Power(errM4,2)+TMath::Power(100.0*TMath::Sqrt(nMPTResolAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus ,2));
    double errM6 = 
        TMath::Sqrt(TMath::Power(errM5,2)+TMath::Power(100.0*TMath::Sqrt(nTrackingAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus ,2));
    double errM7 = 
        TMath::Sqrt(TMath::Power(errM6,2)+TMath::Power(100.0*TMath::Sqrt(nMuonRecoAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus ,2));
    double errM8 = 
        TMath::Sqrt(TMath::Power(errM7,2)+TMath::Power(100.0*TMath::Sqrt(nPtResolAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus ,2));
    double errM9 = 
        TMath::Sqrt(TMath::Power(errM8,2)+TMath::Power(100.0*TMath::Sqrt(nTriggerAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus ,2));

    percentSyst1 = 100.0*TMath::Sqrt(nRaaScalingAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus;
    percentSyst2 = 100.0*TMath::Sqrt(nEWAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus;
    percentSyst3 = 100.0*TMath::Sqrt(nOtherBkgAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus;
    percentSyst4 = 100.0*TMath::Sqrt(nIsolationAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus;
    percentSyst5 = 100.0*TMath::Sqrt(nMPTResolAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus;
    // removed tracking syst
    percentSyst6 = 0.0*100.0*TMath::Sqrt(nTrackingAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus;
    percentSyst7 = 100.0*TMath::Sqrt(nMuonRecoAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus;
    percentSyst8 = 100.0*TMath::Sqrt(nPtResolAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus;
    percentSyst9 = 100.0*TMath::Sqrt(nTriggerAbsoluteSystSumSqEtaMinus)/_nSignalEtaMinus;
    // Save to csv file
    /*writeToSpreadSheet(vecSpreadSheetEtaMuMinus,sEta[ieta],format(percentSyst1,2),format(percentSyst2,2),
				format(percentSyst3,2),format(percentSyst4,2),format(percentSyst5,2),
				format(percentSyst6,2),format(percentSyst7,2),format(percentSyst8,2),format(percentSyst9,2));
*/
    vecSpreadEtaRangeMuMinus.push_back(sEta[ieta]);
    vecSpreadEtaSyst1MuMinus.push_back(format(percentSyst1,1));
    vecSpreadEtaSyst2MuMinus.push_back(format(percentSyst2,1));
    vecSpreadEtaSyst3MuMinus.push_back(format(percentSyst3,1));
    vecSpreadEtaSyst4MuMinus.push_back(format(percentSyst4,1));
    vecSpreadEtaSyst5MuMinus.push_back(format(percentSyst5,1));
    vecSpreadEtaSyst6MuMinus.push_back(format(percentSyst6,1));
    vecSpreadEtaSyst7MuMinus.push_back(format(percentSyst7,1));
    vecSpreadEtaSyst8MuMinus.push_back(format(percentSyst8,1));
    vecSpreadEtaSyst9MuMinus.push_back(format(percentSyst9,1));

    hRaaScalingAbsoluteSystSumSqEtaMinus ->SetBinError(ieta+1, errM1);
    hEWAbsoluteSystSumSqEtaMinus ->SetBinError(ieta+1, errM2);
    hOtherBkgAbsoluteSystSumSqEtaMinus ->SetBinError(ieta+1, errM3);
    hIsolationAbsoluteSystSumSqEtaMinus ->SetBinError(ieta+1, errM4);
    hMPTResolAbsoluteSystSumSqEtaMinus ->SetBinError(ieta+1, errM5);
    hTrackingAbsoluteSystSumSqEtaMinus ->SetBinError(ieta+1, errM6);
    hMuonRecoAbsoluteSystSumSqEtaMinus ->SetBinError(ieta+1, errM7);
    hPtResolAbsoluteSystSumSqEtaMinus ->SetBinError(ieta+1, errM8);
    hTriggerAbsoluteSystSumSqEtaMinus ->SetBinError(ieta+1, errM9);

   } //ieta

   // Create csv 
    writeToSpreadSheet(spreadSheetCentMuPlus,vecSpreadCentRangeMuPlus,vecSpreadCentSyst1MuPlus,vecSpreadCentSyst2MuPlus,
			vecSpreadCentSyst3MuPlus,vecSpreadCentSyst4MuPlus,vecSpreadCentSyst5MuPlus,vecSpreadCentSyst6MuPlus,
			vecSpreadCentSyst7MuPlus,vecSpreadCentSyst8MuPlus,vecSpreadCentSyst9MuPlus);
    writeToSpreadSheet(spreadSheetCentMuMinus,vecSpreadCentRangeMuMinus,vecSpreadCentSyst1MuMinus,vecSpreadCentSyst2MuMinus,
			vecSpreadCentSyst3MuMinus,vecSpreadCentSyst4MuMinus,vecSpreadCentSyst5MuMinus,vecSpreadCentSyst6MuMinus,
			vecSpreadCentSyst7MuMinus,vecSpreadCentSyst8MuMinus,vecSpreadCentSyst9MuMinus);
    writeToSpreadSheet(spreadSheetEtaMuPlus,vecSpreadEtaRangeMuPlus,vecSpreadEtaSyst1MuPlus,vecSpreadEtaSyst2MuPlus,
			vecSpreadEtaSyst3MuPlus,vecSpreadEtaSyst4MuPlus,vecSpreadEtaSyst5MuPlus,vecSpreadEtaSyst6MuPlus,
			vecSpreadEtaSyst7MuPlus,vecSpreadEtaSyst8MuPlus,vecSpreadEtaSyst9MuPlus);
    writeToSpreadSheet(spreadSheetEtaMuMinus,vecSpreadEtaRangeMuMinus,vecSpreadEtaSyst1MuMinus,vecSpreadEtaSyst2MuMinus,
			vecSpreadEtaSyst3MuMinus,vecSpreadEtaSyst4MuMinus,vecSpreadEtaSyst5MuMinus,vecSpreadEtaSyst6MuMinus,
			vecSpreadEtaSyst7MuMinus,vecSpreadEtaSyst8MuMinus,vecSpreadEtaSyst9MuMinus);

   ///mu+ eta
   TCanvas *c1 = new TCanvas("c1", "c1",856,64,700,500);
   c1->Range(-0.3658228,-23.10127,2.54557,27.53165);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.16);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.16);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);

   Int_t ci;   // for color index setting

   hTriggerAbsoluteSystSumSqEtaPlus->SetMinimum(-15);
   hTriggerAbsoluteSystSumSqEtaPlus->SetMaximum(25);
   //ci = TColor::GetColor("#00ffff");
   ci = TColor::GetColor("#ffff99");
   hTriggerAbsoluteSystSumSqEtaPlus->SetFillColor(ci);
   hTriggerAbsoluteSystSumSqEtaPlus->GetXaxis()->SetTitle("|#eta|");
   hTriggerAbsoluteSystSumSqEtaPlus->GetYaxis()->SetTitle("Systematic Uncertainties[%]");
   hTriggerAbsoluteSystSumSqEtaPlus->Draw("e2");

   ci = TColor::GetColor("#99ffff");
   hPtResolAbsoluteSystSumSqEtaPlus->SetFillColor(ci);
   hPtResolAbsoluteSystSumSqEtaPlus->Draw("e2same");

   ci = TColor::GetColor("#ff99ff");
   hMuonRecoAbsoluteSystSumSqEtaPlus->SetFillColor(ci);
   hMuonRecoAbsoluteSystSumSqEtaPlus->Draw("e2same");
   
   ci = TColor::GetColor("#ff9999");
   //hTrackingAbsoluteSystSumSqEtaPlus->SetFillColor(ci);
   //hTrackingAbsoluteSystSumSqEtaPlus->Draw("e2same");

   //ci = TColor::GetColor("#cc99ff");
   hMPTResolAbsoluteSystSumSqEtaPlus->SetFillColor(ci);
   hMPTResolAbsoluteSystSumSqEtaPlus->Draw("e2same");

   ci = TColor::GetColor("#9999ff");
   hIsolationAbsoluteSystSumSqEtaPlus->SetFillColor(ci);
   hIsolationAbsoluteSystSumSqEtaPlus->Draw("e2same");

   ci = TColor::GetColor("#ff3333");
   hOtherBkgAbsoluteSystSumSqEtaPlus->SetFillColor(ci);
   hOtherBkgAbsoluteSystSumSqEtaPlus->Draw("e2same");

   ci = TColor::GetColor("#0033ff");
   hEWAbsoluteSystSumSqEtaPlus->SetFillColor(ci);
   hEWAbsoluteSystSumSqEtaPlus->Draw("e2same");

   ci = TColor::GetColor("#99ff99");
   hRaaScalingAbsoluteSystSumSqEtaPlus->SetFillColor(ci);
   hRaaScalingAbsoluteSystSumSqEtaPlus->Draw("e2same");
    
   hTriggerAbsoluteSystSumSqEtaPlus->Draw("sameaxis");

   TLegend *leg = new TLegend(0.1724138,0.690678,0.4425287,0.9279661,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->AddEntry("hRaaScalingAbsoluteSystSumSqEtaPlus","QCD","fe");
   leg->AddEntry("hEWAbsoluteSystSumSqEtaPlus","EW","fe");
   leg->AddEntry("hOtherBkgAbsoluteSystSumSqEtaPlus","Other Background","fe");
   leg->AddEntry("hIsolationAbsoluteSystSumSqEtaPlus","Isolation","fe");
   leg->Draw();
  
   TLegend* leg2 = new TLegend(0.589,0.622,0.859,0.934,NULL,"brNDC");
   leg2->SetBorderSize(0);
   leg2->SetTextSize(0.05);
   leg2->SetLineColor(1);
   leg2->SetLineStyle(1);
   leg2->SetLineWidth(1);
   leg2->SetFillColor(0);
   leg2->SetFillStyle(1001);
   leg2->AddEntry("hMPTResolAbsoluteSystSumSqEtaPlus","#slash{p_{T}} Resolution","fe");
   //leg2->AddEntry("hTrackingAbsoluteSystSumSqEtaPlus","Tracking","fe");
   leg2->AddEntry("hMuonRecoAbsoluteSystSumSqEtaPlus","Muon reconstruction","fe");
   leg2->AddEntry("hPtResolAbsoluteSystSumSqEtaPlus","p_{T}^{#mu} resolution","fe");
   leg2->AddEntry("hTriggerAbsoluteSystSumSqEtaPlus","Trigger","fe");
   leg2->Draw();

   TLatex *   tex = new TLatex(0.3045977,0.2351695,"#mu^{+}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);

   //mu+ cent








   TCanvas *c2 = new TCanvas("c2", "c2",856,64,700,500);
   c2->Range(-0.3658228,-23.10127,2.54557,27.53165);
   c2->SetFillColor(0);
   c2->SetBorderMode(0);
   c2->SetBorderSize(2);
   c2->SetTickx(1);
   c2->SetTicky(1);
   c2->SetLeftMargin(0.16);
   c2->SetRightMargin(0.05);
   c2->SetTopMargin(0.05);
   c2->SetBottomMargin(0.16);
   c2->SetFrameBorderMode(0);
   c2->SetFrameBorderMode(0);

   ci = TColor::GetColor("#ffff99");
   hTriggerAbsoluteSystSumSqCentPlus->SetMinimum(-15);
   hTriggerAbsoluteSystSumSqCentPlus->SetMaximum(25);
   hTriggerAbsoluteSystSumSqCentPlus->SetFillColor(ci);
   hTriggerAbsoluteSystSumSqCentPlus->GetXaxis()->SetTitle("Centrality [%]");
   hTriggerAbsoluteSystSumSqCentPlus->GetYaxis()->SetTitle("Systematic Uncertainties[%]");
   hTriggerAbsoluteSystSumSqCentPlus->Draw("e2");

   ci = TColor::GetColor("#99ffff");
   hPtResolAbsoluteSystSumSqCentPlus->SetFillColor(ci);
   hPtResolAbsoluteSystSumSqCentPlus->Draw("e2same");

   ci = TColor::GetColor("#ff99ff");
   hMuonRecoAbsoluteSystSumSqCentPlus->SetFillColor(ci);
   hMuonRecoAbsoluteSystSumSqCentPlus->Draw("e2same");
   
   ci = TColor::GetColor("#ff9999");
   //hTrackingAbsoluteSystSumSqCentPlus->SetFillColor(ci);
   //hTrackingAbsoluteSystSumSqCentPlus->Draw("e2same");


   //ci = TColor::GetColor("#cc99ff");
   hMPTResolAbsoluteSystSumSqCentPlus->SetFillColor(ci);
   hMPTResolAbsoluteSystSumSqCentPlus->Draw("e2same");

   ci = TColor::GetColor("#9999ff");
   hIsolationAbsoluteSystSumSqCentPlus->SetFillColor(ci);
   hIsolationAbsoluteSystSumSqCentPlus->Draw("e2same");

   ci = TColor::GetColor("#ff3333");
   hOtherBkgAbsoluteSystSumSqCentPlus->SetFillColor(ci);
   hOtherBkgAbsoluteSystSumSqCentPlus->Draw("e2same");

   ci = TColor::GetColor("#0033ff");
   hEWAbsoluteSystSumSqCentPlus->SetFillColor(ci);
   hEWAbsoluteSystSumSqCentPlus->Draw("e2same");

   ci = TColor::GetColor("#99ff99");
   hRaaScalingAbsoluteSystSumSqCentPlus->SetFillColor(ci);
   hRaaScalingAbsoluteSystSumSqCentPlus->Draw("e2same");

   hTriggerAbsoluteSystSumSqCentPlus->Draw("sameaxis");

   leg->Draw();
  
   leg2->Draw();

   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
   c2->Modified();
   c2->cd();
   c2->SetSelected(c2);

   TCanvas *c3 = new TCanvas("c3", "c3",856,64,700,500);
   c3->Range(-0.3658228,-23.10127,2.54557,27.53165);
   c3->SetFillColor(0);
   c3->SetBorderMode(0);
   c3->SetBorderSize(2);
   c3->SetTickx(1);
   c3->SetTicky(1);
   c3->SetLeftMargin(0.16);
   c3->SetRightMargin(0.05);
   c3->SetTopMargin(0.05);
   c3->SetBottomMargin(0.16);
   c3->SetFrameBorderMode(0);
   c3->SetFrameBorderMode(0);

   ci = TColor::GetColor("#ffff99");
   hTriggerAbsoluteSystSumSqCentMinus->SetMinimum(-15);
   hTriggerAbsoluteSystSumSqCentMinus->SetMaximum(25);
   hTriggerAbsoluteSystSumSqCentMinus->SetFillColor(ci);
   hTriggerAbsoluteSystSumSqCentMinus->GetXaxis()->SetTitle("Centrality [%]");
   hTriggerAbsoluteSystSumSqCentMinus->GetYaxis()->SetTitle("Systematic Uncertainties[%]");
   hTriggerAbsoluteSystSumSqCentMinus->Draw("e2");

   ci = TColor::GetColor("#99ffff");
   hPtResolAbsoluteSystSumSqCentMinus->SetFillColor(ci);
   hPtResolAbsoluteSystSumSqCentMinus->Draw("e2same");

   ci = TColor::GetColor("#ff99ff");
   hMuonRecoAbsoluteSystSumSqCentMinus->SetFillColor(ci);
   hMuonRecoAbsoluteSystSumSqCentMinus->Draw("e2same");
   
   ci = TColor::GetColor("#ff9999");
   //hTrackingAbsoluteSystSumSqCentMinus->SetFillColor(ci);
   //hTrackingAbsoluteSystSumSqCentMinus->Draw("e2same");


   //ci = TColor::GetColor("#cc99ff");
   hMPTResolAbsoluteSystSumSqCentMinus->SetFillColor(ci);
   hMPTResolAbsoluteSystSumSqCentMinus->Draw("e2same");

   ci = TColor::GetColor("#9999ff");
   hIsolationAbsoluteSystSumSqCentMinus->SetFillColor(ci);
   hIsolationAbsoluteSystSumSqCentMinus->Draw("e2same");

   ci = TColor::GetColor("#ff3333");
   hOtherBkgAbsoluteSystSumSqCentMinus->SetFillColor(ci);
   hOtherBkgAbsoluteSystSumSqCentMinus->Draw("e2same");

   ci = TColor::GetColor("#0033ff");
   hEWAbsoluteSystSumSqCentMinus->SetFillColor(ci);
   hEWAbsoluteSystSumSqCentMinus->Draw("e2same");

   ci = TColor::GetColor("#99ff99");
   hRaaScalingAbsoluteSystSumSqCentMinus->SetFillColor(ci);
   hRaaScalingAbsoluteSystSumSqCentMinus->Draw("e2same");

   hTriggerAbsoluteSystSumSqCentMinus->Draw("sameaxis");

   leg->Draw();
  
   leg2->Draw();
   TLatex *   tex2 = new TLatex(0.3045977,0.2351695,"#mu^{-}");
   tex2->SetNDC();
   tex2->SetTextFont(42);
   tex2->SetLineWidth(2);
   tex2->Draw();

   c3->Modified();
   c3->cd();
   c3->SetSelected(c2);

   TCanvas *c4 = new TCanvas("c4", "c4",856,64,700,500);
   c4->Range(-0.3658228,-23.10127,2.54557,27.53165);
   c4->SetFillColor(0);
   c4->SetBorderMode(0);
   c4->SetBorderSize(2);
   c4->SetTickx(1);
   c4->SetTicky(1);
   c4->SetLeftMargin(0.16);
   c4->SetRightMargin(0.05);
   c4->SetTopMargin(0.05);
   c4->SetBottomMargin(0.16);
   c4->SetFrameBorderMode(0);
   c4->SetFrameBorderMode(0);

   ci = TColor::GetColor("#ffff99");
   hTriggerAbsoluteSystSumSqEtaMinus->SetMinimum(-15);
   hTriggerAbsoluteSystSumSqEtaMinus->SetMaximum(25);
   hTriggerAbsoluteSystSumSqEtaMinus->SetFillColor(ci);
   hTriggerAbsoluteSystSumSqEtaMinus->GetXaxis()->SetTitle("|#eta|");
   hTriggerAbsoluteSystSumSqEtaMinus->GetYaxis()->SetTitle("Systematic Uncertainties[%]");
   hTriggerAbsoluteSystSumSqEtaMinus->Draw("e2");

   ci = TColor::GetColor("#99ffff");
   hPtResolAbsoluteSystSumSqEtaMinus->SetFillColor(ci);
   hPtResolAbsoluteSystSumSqEtaMinus->Draw("e2same");

   ci = TColor::GetColor("#ff99ff");
   hMuonRecoAbsoluteSystSumSqEtaMinus->SetFillColor(ci);
   hMuonRecoAbsoluteSystSumSqEtaMinus->Draw("e2same");
   
   ci = TColor::GetColor("#ff9999");
   //hTrackingAbsoluteSystSumSqEtaMinus->SetFillColor(ci);
   //hTrackingAbsoluteSystSumSqEtaMinus->Draw("e2same");


   //ci = TColor::GetColor("#cc99ff");
   hMPTResolAbsoluteSystSumSqEtaMinus->SetFillColor(ci);
   hMPTResolAbsoluteSystSumSqEtaMinus->Draw("e2same");

   ci = TColor::GetColor("#9999ff");
   hIsolationAbsoluteSystSumSqEtaMinus->SetFillColor(ci);
   hIsolationAbsoluteSystSumSqEtaMinus->Draw("e2same");

   ci = TColor::GetColor("#ff3333");
   hOtherBkgAbsoluteSystSumSqEtaMinus->SetFillColor(ci);
   hOtherBkgAbsoluteSystSumSqEtaMinus->Draw("e2same");

   ci = TColor::GetColor("#0033ff");
   hEWAbsoluteSystSumSqEtaMinus->SetFillColor(ci);
   hEWAbsoluteSystSumSqEtaMinus->Draw("e2same");

   ci = TColor::GetColor("#99ff99");
   hRaaScalingAbsoluteSystSumSqEtaMinus->SetFillColor(ci);
   hRaaScalingAbsoluteSystSumSqEtaMinus->Draw("e2same");

   hTriggerAbsoluteSystSumSqEtaMinus->Draw("sameaxis");

   leg->Draw();
  
   leg2->Draw();

   tex2->SetNDC();
   tex2->SetTextFont(42);
   tex2->SetLineWidth(2);
   tex2->Draw();

   c4->Modified();
   c4->cd();
   c4->SetSelected(c2);


   TDatime* time = new TDatime();
   TString sDate = "";
   sDate+=time->GetMonth(); sDate+="_"; sDate+=time->GetDay(); sDate+="_"; sDate+=time->GetYear();
   std::cout << "Today is " << sDate << std::endl;


   c1->Print("muPlusSystematicsEta_FullyCorrelated_"+sDate+".pdf");
   c2->Print("muPlusSystematicsCent_FullyCorrelated_"+sDate+".pdf");
   c3->Print("muMinusSystematicsCent_FullyCorrelated_"+sDate+".pdf");
   c4->Print("muMinusSystematicsEta_FullyCorrelated_"+sDate+".pdf");

   spreadSheetCentMuPlus.close();
   spreadSheetCentMuMinus.close();
   spreadSheetEtaMuPlus.close();
   spreadSheetEtaMuMinus.close();

   outFile->Write();
   
//   outFile->Close();
/*   for(int ieta=0; ieta<nEtaBins; ++ieta){
    delete h2D_RaaScalingPercentSyst;
    delete h2D_IsolationPercentSyst;
    delete h2D_MPTResolPercentSyst;
    delete h2D_TrackingPercentSyst;
    delete h2D_MuonRecoPercentSyst;
    delete h2D_PtResolPercentSyst;
    delete h2D_TriggerPercentSyst;

}
*/

}
