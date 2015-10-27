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
#include <vector> 



void Write(TFile* const outFile, TObject* const gr, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing TGraph to root file..." << std::endl;
	  gr->Write(sFile);
	  std::cout << "Done." << std::endl;
    	  gDirectory = dir;
}

void plotZBkgFractionInData()
{
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

    const unsigned int _nEta = etaBins.size()-1;

    std::vector <double> nPartBins;
    nPartBins.push_back(382.16);//0-5
	nPartBins.push_back(330.26);//5-10
	nPartBins.push_back(281.88);//10-15
	nPartBins.push_back(239.52);//15-20
	nPartBins.push_back(157.83);//20-40
	nPartBins.push_back(45.93);//40-80

    const unsigned int _nCentrality = nPartBins.size();
    const unsigned int _nBins = _nCentrality*_nEta;
    ///stl containers to store Z and Raw counts from txt file
    std::vector <double> vZPlus(_nBins); 
    std::vector <double> vZMinus(_nBins); 
    std::vector <double> vRawSigPlus(_nBins); 
    std::vector <double> vRawSigMinus(_nBins); 
    std::vector <double> vTotalEtaSigPlus(_nEta); 
    std::vector <double> vTotalEtaSigMinus(_nEta); 
    std::vector <double> vTotalEtaZPlus(_nEta); 
    std::vector <double> vTotalEtaZMinus(_nEta); 
    std::vector <double> vTotalCentSigPlus(_nCentrality); 
    std::vector <double> vTotalCentSigMinus(_nCentrality); 
    std::vector <double> vTotalCentZPlus(_nCentrality); 
    std::vector <double> vTotalCentZMinus(_nCentrality); 

    TGraphErrors* _grBkgFractionPlusEta = new TGraphErrors(_nEta);
    TGraphErrors* _grBkgFractionMinusEta = new TGraphErrors(_nEta);
    TGraphErrors* _grBkgFractionPlusCent = new TGraphErrors(_nCentrality); 
    TGraphErrors* _grBkgFractionMinusCent = new TGraphErrors(_nCentrality); 


    TString fileNameDataOut = "ZBackgroundFractionsInData";
    TFile *outFile = new TFile(fileNameDataOut+".root","RECREATE");
    TString sFileIn = "rawZBkgFraction.06.27.2013.txt";
    std::ifstream s(sFileIn,std::ifstream::in);

	if(!s.is_open()){
		std::cout << "While opening a file an error was encountered" << std::endl;
		exit(1);
	} else {
		std::cout << "File " << sFileIn << " is successfully opened" << std::endl;
	}

    ///Fill stl containers 
    for(unsigned int i=0; i<vRawSigPlus.size(); ++i){
        s >> vRawSigPlus[i] >> vZPlus[i] >> vRawSigMinus[i] >> vZMinus[i] ;
        std::cout << "Raw counts positive charge:             " << vRawSigPlus[i] << std::endl;
        std::cout << "Z counts :                              " << vZPlus[i] << std::endl;
        std::cout << "Raw counts negative charge:             " << vRawSigMinus[i] << std::endl;
        std::cout << "Z counts :                              " << vZMinus[i] << std::endl;
    }

    ///Sum up in eta
    for(int ieta=0; ieta<_nEta; ++ieta){
        ///initialize to zero
        vTotalEtaSigPlus[ieta] = 0.0;  
        vTotalEtaZPlus[ieta] = 0.0;  
        vTotalEtaSigMinus[ieta] = 0.0;  
        vTotalEtaZMinus[ieta] = 0.0;  
        for(int icent=0; icent<_nCentrality; ++icent){
            int index = ieta*_nCentrality+icent;
            vTotalEtaSigPlus[ieta] += vRawSigPlus[index]; 
            vTotalEtaZPlus[ieta] += vZPlus[index]; 
            vTotalEtaSigMinus[ieta] += vRawSigMinus[index]; 
            vTotalEtaZMinus[ieta] += vZMinus[index]; 
        }//icent
        ///mu+
        double xEta = (etaBins[ieta+1]-etaBins[ieta])/2.0+etaBins[ieta];
        double bkgFractionPlus = vTotalEtaZPlus[ieta]/vTotalEtaSigPlus[ieta];
        double bkgFractionPlusErr = bkgFractionPlus*TMath::Sqrt( TMath::Power(TMath::Sqrt(vTotalEtaZPlus[ieta])/vTotalEtaZPlus[ieta],2) 
            + TMath::Power(TMath::Sqrt(vTotalEtaSigPlus[ieta])/vTotalEtaSigPlus[ieta],2) );
        _grBkgFractionPlusEta->SetPoint(ieta,xEta,bkgFractionPlus);
        _grBkgFractionPlusEta->SetPointError(ieta,0.0,bkgFractionPlusErr);

        ///mu-
        double bkgFractionMinus = vTotalEtaZMinus[ieta]/vTotalEtaSigMinus[ieta];
        double bkgFractionMinusErr = bkgFractionMinus*TMath::Sqrt( TMath::Power(TMath::Sqrt(vTotalEtaZMinus[ieta])/vTotalEtaZMinus[ieta],2) 
            + TMath::Power(TMath::Sqrt(vTotalEtaSigMinus[ieta])/vTotalEtaSigMinus[ieta],2) );
        _grBkgFractionMinusEta->SetPoint(ieta,xEta,bkgFractionMinus);
        _grBkgFractionMinusEta->SetPointError(ieta,0.0,bkgFractionMinusErr);
    }//ieta

    ///Sum up in centrality
    for(int icent=0; icent<_nCentrality; ++icent){
        for(int ieta=0; ieta<_nEta; ++ieta){
            int index = ieta*_nCentrality+icent;
            vTotalCentSigPlus[icent] += vRawSigPlus[index]; 
            vTotalCentZPlus[icent] += vZPlus[index]; 
            vTotalCentSigMinus[icent] += vRawSigMinus[index]; 
            vTotalCentZMinus[icent] += vZMinus[index]; 

        } //ieta
        ///mu+
        double xCent = nPartBins[icent];
        double bkgFractionPlus = vTotalCentZPlus[icent]/vTotalCentSigPlus[icent];
        double bkgFractionPlusErr = bkgFractionPlus*TMath::Sqrt( TMath::Power(TMath::Sqrt(vTotalCentZPlus[icent])/vTotalCentZPlus[icent],2) 
            + TMath::Power(TMath::Sqrt(vTotalCentSigPlus[icent])/vTotalCentSigPlus[icent],2) );
        _grBkgFractionPlusCent->SetPoint(icent,xCent,bkgFractionPlus);
        _grBkgFractionPlusCent->SetPointError(icent,0.0,bkgFractionPlusErr);

        ///mu-
        double bkgFractionMinus = vTotalCentZMinus[icent]/vTotalCentSigMinus[icent];
        double bkgFractionMinusErr = bkgFractionMinus*TMath::Sqrt( TMath::Power(TMath::Sqrt(vTotalCentZMinus[icent])/vTotalCentZMinus[icent],2) 
            + TMath::Power(TMath::Sqrt(vTotalCentSigMinus[icent])/vTotalCentSigMinus[icent],2) );
        _grBkgFractionMinusCent->SetPoint(icent,xCent,bkgFractionMinus);
        _grBkgFractionMinusCent->SetPointError(icent,0.0,bkgFractionMinusErr);

    }//icents

    ///Save TGraphs to file
    Write(outFile,_grBkgFractionPlusEta,"grZBkgFractionPlusEta");
    Write(outFile,_grBkgFractionMinusEta,"grZBkgFractionMinusEta");
    Write(outFile,_grBkgFractionPlusCent,"grZBkgFractionPlusCent");
    Write(outFile,_grBkgFractionMinusCent,"grZBkgFractionMinusCent");

}
