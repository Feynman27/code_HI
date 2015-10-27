//////////////////////////////////////////////////////////////
//This macro plots the relative systematic errors from raising and
//lowering the mpt cut by the respective sigma along with the RMS 
//wrt the nominal yield in that bin
///////////////////////////////////////////////////////////////
#include "TString.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TList.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>
#include <vector>

const unsigned int _nCentrality = 6;
const unsigned int _nEta = 9;
const unsigned int _nBins = _nCentrality*_nEta;
std::vector <double> mptRMSMuPlus (_nBins);
std::vector <double> mptRMSMuMinus (_nBins);
std::vector <double> mptErrHiMuPlus (_nBins);
std::vector <double> mptErrHiMuMinus (_nBins);
std::vector <double> mptErrLoMuPlus (_nBins);
std::vector <double> mptErrLoMuMinus (_nBins);
std::vector <double> WYieldsPlus (_nBins);
std::vector <double> WYieldsMinus (_nBins);

void Write(TFile* const outFile, TObject* const gr, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing TGraph to root file..." << std::endl;
	  gr->Write(sFile);
	  std::cout << "Done." << std::endl;
    	  gDirectory = dir;
}

void setWYields(TString sFileIn = "systematics/mptCutSystematics.04.24.2013.txt"){
    
	unsigned int totalBins =  _nBins;
	std::ifstream s( sFileIn , std::ifstream::in );

	if(!s.is_open()){
		std::cout << "While opening a file an error was encountered" << std::endl;
		exit(1);
	} else {
		std::cout << "File " << sFileIn << " is successfully opened" << std::endl;
	}

    if(s.is_open()){
        for (unsigned int i = 0; i < WYieldsPlus.size(); ++i) {
		    s >> WYieldsPlus[i] >> WYieldsMinus[i];
            std::cout << "MPT systematic : mu+: " << WYieldsPlus[i] << " mu- : " << WYieldsMinus[i] << std::endl;  
        }
        if( (WYieldsPlus.size()!=totalBins) || (WYieldsMinus.size()!=totalBins) ){
            std::cout << "ERROR: mpt container size does not match number of total bins used. Aborting." << std::endl;
            exit(0);
        }
        s.close();
	}
    else std::cout << "unable to open mpt systematics file" << std::endl;
}

float getWYields(int chargeIndex=-1, int binIndex=-1){
    if(chargeIndex==-1||binIndex==-1){
        std::cout << "ERROR READING INPUT. " << std::endl;
        return -1.0;
    }
    else if(chargeIndex==102) {
        std::cout << "Return mpt systematic for mu+ at index: " << binIndex << " = " << WYieldsPlus[binIndex] << std::endl;
        return WYieldsPlus[binIndex];
    }
    else if(chargeIndex==103){
        std::cout << "Return isolation systematic for mu- at index: " << binIndex << " = " << WYieldsMinus[binIndex] << std::endl;
        return WYieldsMinus[binIndex];
    }
    else {

        std::cout << "ERROR READING INPUT. " << std::endl;
        return -1.0;
    }
}


void setMptHiErrors(TString sFileIn = "systematics/mptCutSystematics.04.24.2013.txt"){
    
	unsigned int totalBins =  _nBins;
	std::ifstream s( sFileIn , std::ifstream::in );

	if(!s.is_open()){
		std::cout << "While opening a file an error was encountered" << std::endl;
		exit(1);
	} else {
		std::cout << "File " << sFileIn << " is successfully opened" << std::endl;
	}

    if(s.is_open()){
        for (unsigned int i = 0; i < mptErrHiMuPlus.size(); ++i) {
		    s >> mptErrHiMuPlus[i] >> mptErrHiMuMinus[i];
            std::cout << "MPT systematic : mu+: " << mptErrHiMuPlus[i] << " mu- : " << mptErrHiMuMinus[i] << std::endl;  
        }
        if( (mptErrHiMuPlus.size()!=totalBins) || (mptErrHiMuMinus.size()!=totalBins) ){
            std::cout << "ERROR: mpt container size does not match number of total bins used. Aborting." << std::endl;
            exit(0);
        }
        s.close();
	}
    else std::cout << "unable to open mpt systematics file" << std::endl;
}

float getMPtErrHi(int chargeIndex=-1, int binIndex=-1){
    if(chargeIndex==-1||binIndex==-1){
        std::cout << "ERROR READING INPUT. " << std::endl;
        return -1.0;
    }
    else if(chargeIndex==102) {
        std::cout << "Return mpt systematic for mu+ at index: " << binIndex << " = " << mptErrHiMuPlus[binIndex] << std::endl;
        return mptErrHiMuPlus[binIndex];
    }
    else if(chargeIndex==103){
        std::cout << "Return isolation systematic for mu- at index: " << binIndex << " = " << mptErrHiMuMinus[binIndex] << std::endl;
        return mptErrHiMuMinus[binIndex];
    }
    else {

        std::cout << "ERROR READING INPUT. " << std::endl;
        return -1.0;
    }
}



void setMptLoErrors(TString sFileIn = "systematics/mptCutSystematics.04.24.2013.txt"){
    
	unsigned int totalBins =  _nBins;
	std::ifstream s( sFileIn , std::ifstream::in );

	if(!s.is_open()){
		std::cout << "While opening a file an error was encountered" << std::endl;
		exit(1);
	} else {
		std::cout << "File " << sFileIn << " is successfully opened" << std::endl;
	}

    if(s.is_open()){
        for (unsigned int i = 0; i < mptErrLoMuPlus.size(); ++i) {
		    s >> mptErrLoMuPlus[i] >> mptErrLoMuMinus[i];
            std::cout << "MPT systematic : mu+: " << mptErrLoMuPlus[i] << " mu- : " << mptErrLoMuMinus[i] << std::endl;  
        }
        if( (mptErrLoMuPlus.size()!=totalBins) || (mptErrLoMuMinus.size()!=totalBins) ){
            std::cout << "ERROR: mpt container size does not match number of total bins used. Aborting." << std::endl;
            exit(0);
        }
        s.close();
	}
    else std::cout << "unable to open mpt systematics file" << std::endl;
}

float getMPtErrLo(int chargeIndex=-1, int binIndex=-1){
    if(chargeIndex==-1||binIndex==-1){
        std::cout << "ERROR READING INPUT. " << std::endl;
        return -1.0;
    }
    else if(chargeIndex==102) {
        std::cout << "Return mpt systematic for mu+ at index: " << binIndex << " = " << mptErrLoMuPlus[binIndex] << std::endl;
        return mptErrLoMuPlus[binIndex];
    }
    else if(chargeIndex==103){
        std::cout << "Return isolation systematic for mu- at index: " << binIndex << " = " << mptErrLoMuMinus[binIndex] << std::endl;
        return mptErrLoMuMinus[binIndex];
    }
    else {

        std::cout << "ERROR READING INPUT. " << std::endl;
        return -1.0;
    }
}


void setMptRMSErrors(TString sFileIn = "systematics/mptCutSystematics.04.24.2013.txt"){
    
	unsigned int totalBins =  _nBins;
	std::ifstream s( sFileIn , std::ifstream::in );

	if(!s.is_open()){
		std::cout << "While opening a file an error was encountered" << std::endl;
		exit(1);
	} else {
		std::cout << "File " << sFileIn << " is successfully opened" << std::endl;
	}

    if(s.is_open()){
        for (unsigned int i = 0; i < mptRMSMuPlus.size(); ++i) {
		    s >> mptRMSMuPlus[i] >> mptRMSMuMinus[i];
            std::cout << "MPT systematic : mu+: " << mptRMSMuPlus[i] << " mu- : " << mptRMSMuMinus[i] << std::endl;  
        }
        if( (mptRMSMuPlus.size()!=totalBins) || (mptRMSMuMinus.size()!=totalBins) ){
            std::cout << "ERROR: mpt container size does not match number of total bins used. Aborting." << std::endl;
            exit(0);
        }
        s.close();
	}
    else std::cout << "unable to open mpt systematics file" << std::endl;
}

///return PERCENT systematic error in mpt cuts
///for given charge,centrality, and eta class
float getMptRMSErr(int chargeIndex=-1, int binIndex=-1){
    if(chargeIndex==-1||binIndex==-1){
        std::cout << "ERROR READING INPUT. " << std::endl;
        return -1.0;
    }
    else if(chargeIndex==102) {
        std::cout << "Return mpt systematic for mu+ at index: " << binIndex << " = " << mptRMSMuPlus[binIndex] << std::endl;
        return mptRMSMuPlus[binIndex];
    }
    else if(chargeIndex==103){
        std::cout << "Return isolation systematic for mu- at index: " << binIndex << " = " << mptRMSMuMinus[binIndex] << std::endl;
        return mptRMSMuMinus[binIndex];
    }
    else {

        std::cout << "ERROR READING INPUT. " << std::endl;
        return -1.0;
    }
}

void plotMPtSystematicErrors(){

TString fileNameDataOut = "mptSystematicErrorEtaDistros";
TFile *outFile = new TFile(fileNameDataOut+".root","RECREATE");

setWYields("wYields.06.09.2013.txt");
setMptHiErrors("systematics/mptUpPercentSystError.06.23.2013.txt");
setMptLoErrors("systematics/mptDownPercentSystError.06.23.2013.txt");
//setMptRMSErrors("systematics/mptUpPercentSystError.06.09.2013.txt");

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

std::vector<double> centralityBins;
centralityBins.push_back(0.0);
centralityBins.push_back(0.05);
centralityBins.push_back(0.10);
centralityBins.push_back(0.15);
centralityBins.push_back(0.20);
centralityBins.push_back(0.40);
centralityBins.push_back(0.80);

const int nCentralityBins = centralityBins.size()-1;

TList _grMptRMSMuPlus,_grMptErrHiMuPlus,_grMptErrLoMuPlus,_grWYieldsPlus;
TList _grMptRMSMuMinus,_grMptErrHiMuMinus,_grMptErrLoMuMinus,_grWYieldsMinus;

for(int j=0; j<nCentralityBins; ++j){
     //_grMptRMSMuPlus.Add( new TGraphErrors(nEtaBins) );
     _grMptErrHiMuPlus.Add( new TGraphErrors(nEtaBins) );
     _grMptErrLoMuPlus.Add( new TGraphErrors(nEtaBins) );
     _grWYieldsPlus.Add( new TGraphErrors(nEtaBins) );

     //_grMptRMSMuMinus.Add( new TGraphErrors(nEtaBins) );
     _grMptErrHiMuMinus.Add( new TGraphErrors(nEtaBins) );
     _grMptErrLoMuMinus.Add( new TGraphErrors(nEtaBins) );
     _grWYieldsMinus.Add( new TGraphErrors(nEtaBins) );
  }

for (int ieta=0; ieta<nEtaBins; ++ieta){
    for (int icent=0; icent<nCentralityBins; ++icent){
       int index = ieta*nCentralityBins+icent; 

       double etaX = etaBins[ieta]+(etaBins[ieta+1]-etaBins[ieta])/2.0;

       double wPlus = getWYields(102, index) ;
       ( (TGraphErrors*)_grWYieldsPlus.At(icent) )->SetPoint(ieta,etaX,wPlus);
       ( (TGraphErrors*)_grWYieldsPlus.At(icent) )->SetPointError(ieta,(etaBins[ieta+1]-etaBins[ieta])/2.0,0.0);
       double wMinus = getWYields(103, index) ;
       ( (TGraphErrors*)_grWYieldsMinus.At(icent) )->SetPoint(ieta,etaX,wMinus);
       ( (TGraphErrors*)_grWYieldsMinus.At(icent) )->SetPointError(ieta,(etaBins[ieta+1]-etaBins[ieta])/2.0,0.0);

       double mptErrHiPlus = getMPtErrHi(102, index) ;
       ( (TGraphErrors*)_grMptErrHiMuPlus.At(icent) )->SetPoint(ieta,etaX,0.0);
       ( (TGraphErrors*)_grMptErrHiMuPlus.At(icent))->SetPointError(ieta,0.0,mptErrHiPlus);
       double mptErrHiMinus = getMPtErrHi(103, index) ;
       ( (TGraphErrors*)_grMptErrHiMuMinus.At(icent) )->SetPoint(ieta,etaX,0.0);
       ( (TGraphErrors*)_grMptErrHiMuMinus.At(icent))->SetPointError(ieta,0.0,mptErrHiMinus);

       double mptErrLoPlus = getMPtErrLo(102, index) ;
       ( (TGraphErrors*)_grMptErrLoMuPlus.At(icent) )->SetPoint(ieta,etaX,0.0);
       ( (TGraphErrors*)_grMptErrLoMuPlus.At(icent))->SetPointError(ieta,0.0,mptErrLoPlus);
       double mptErrLoMinus = getMPtErrLo(103, index) ;
       ( (TGraphErrors*)_grMptErrLoMuMinus.At(icent) )->SetPoint(ieta,etaX,0.0);
       ( (TGraphErrors*)_grMptErrLoMuMinus.At(icent))->SetPointError(ieta,0.0,mptErrLoMinus);

       /*double mptErrRMSPlus = getMptRMSErr(102, index) ;
       ( (TGraphErrors*)_grMptRMSMuPlus.At(icent) )->SetPoint(ieta,etaX,0.0);
       ( (TGraphErrors*)_grMptRMSMuPlus.At(icent))->SetPointError(ieta,(etaBins[ieta+1]-etaBins[ieta])/2.0,mptErrRMSPlus);
       double mptErrRMSMinus = getMptRMSErr(103, index) ;
       ( (TGraphErrors*)_grMptRMSMuMinus.At(icent) )->SetPoint(ieta,etaX,0.0);
       ( (TGraphErrors*)_grMptRMSMuMinus.At(icent))->SetPointError(ieta,(etaBins[ieta+1]-etaBins[ieta])/2.0,mptErrRMSMinus);
        */
    }
}
for(int i=0;i<_grMptErrHiMuPlus.GetEntries();++i){
        //TString sGrName0 = "grMptRMSPercentErrorEtaDistroCent"; sGrName0+=i; sGrName0+="Charge";
        //Write(outFile,_grMptRMSMuPlus.At(i),sGrName0+"102");
        //Write(outFile,_grMptRMSMuMinus.At(i),sGrName0+"103");

        TString sGrName1 = "grMptUpPercentErrorEtaDistroCent"; sGrName1+=i; sGrName1+="Charge";
        Write(outFile,_grMptErrHiMuPlus.At(i),sGrName1+"102");
        Write(outFile,_grMptErrHiMuMinus.At(i),sGrName1+"103");

        TString sGrname2 = "grMptDownPercentErrorEtaDistroCent"; sGrname2+=i; sGrname2+="Charge";
        Write(outFile,_grMptErrLoMuPlus.At(i),sGrname2+"102");
        Write(outFile,_grMptErrLoMuMinus.At(i),sGrname2+"103");
  }


for(int j=0; j<_grMptRMSMuPlus.GetEntries(); ++j){
     //delete _grMptRMSMuPlus.At(j);
     delete _grMptErrHiMuPlus.At(j);
     delete _grMptErrLoMuPlus.At(j);
     delete _grWYieldsPlus.At(j);

     //delete _grMptRMSMuMinus.At(j);
     delete _grMptErrHiMuMinus.At(j);
     delete _grMptErrLoMuMinus.At(j);
     delete _grWYieldsMinus.At(j);
  }
}
