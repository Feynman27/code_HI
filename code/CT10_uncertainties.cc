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

void fill(TFile* fout, TString sFileIn, TH1F* hEtaPos, TH1F* hEtaNeg, TH1F* hAsymm){
    std::cout << "Opening file " << sFileIn << std::endl;
    std::ifstream s( sFileIn , std::ifstream::in );

    if(!s.is_open()){
        std::cout << "While opening a file an error was encountered" << std::endl;
        exit(1);
    }
    else{
        std::cout << "File " << sFileIn << " is successfully opened" <<
            std::endl;
    }

    const int nBins = 10;
    std:: cout << "Number of eta bins: " << nBins << std::endl;
    std::vector <double> _vEtaBins(nBins), _vEtaHalfBinWidth(nBins);
    std::vector <double> _vEtaMuPos(nBins), _vEtaMuNeg(nBins), _vAsymmMu(nBins);
    std::vector <double> _vEtaMuPosError(nBins), _vEtaMuNegError(nBins), _vAsymmMuError(nBins);

    double arrBin[nBins+1] = {0.0,0.1,0.35,0.6,0.8,1.05,1.37,1.52,1.85,2.01,2.5};

    hEtaPos = new TH1F("hEtaPosRelErr","mu^{+} CT10", nBins, arrBin);
    hEtaNeg = new TH1F("hEtaNegRelErr","mu^{-} CT10", nBins, arrBin);
    hAsymm = new TH1F("hAsymmAbsErr","asymmetry CT10", nBins, arrBin);

    for (unsigned int i = 0; i < _vEtaBins.size(); ++i){

        s >> _vEtaBins[i] >> _vEtaHalfBinWidth[i] >> _vEtaMuPos[i] >> _vEtaMuPosError[i] >> _vEtaMuNeg[i] 
            >>  _vEtaMuNegError[i] >> _vAsymmMu[i] >> _vAsymmMuError[i] ;

        std::cout << "Eta Central: " << _vEtaBins[i] << " +/- " << _vEtaHalfBinWidth[i] << std::endl;

        float errAbsEtaPos = _vEtaMuPosError[i] ; 
        float errAbsEtaNeg = _vEtaMuNegError[i] ;
        
        hEtaPos->SetBinContent(i+1, errAbsEtaPos); 
        hEtaNeg->SetBinContent(i+1, errAbsEtaNeg); 
        hAsymm->SetBinContent(i+1, _vAsymmMuError[i]); 
    }//i

    fout->cd();
    hEtaPos->Write();
    hEtaNeg->Write();
    hAsymm->Write();

    s.close();
    delete hEtaPos;
    delete hEtaNeg;
    delete hAsymm;

}


void CT10_uncertainties(){
    TString sOutFile = "CT10_uncertainties.root";
    TFile *fOut = new TFile(sOutFile,"RECREATE");
    TString sFileIn = "prediction_only_ct10_error.dat";
    TH1F* hRelErrorsPos, *hRelErrorsNeg, *hAbsErrorsAsymm;
    fill(fOut, sFileIn, hRelErrorsPos,hRelErrorsNeg,hAbsErrorsAsymm);
}
