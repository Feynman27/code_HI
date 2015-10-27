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

double findMaxVariation(double err1, double err2, double err3, double err4){

    string sMax = "none";
    std::cout << " Variations: " << err1 << " " << err2 << " " << err3 << " " << err4 << std::endl;
    double max = -1.0;
    if(err1>max) max = err1;
    if(err2>max) max = err2;
    if(err3>max) max = err3;
    if(err4>max) max = err4;
    if(max==err1) sMax = "ud";
    if(max==err2) sMax = "du";
    if(max==err3) sMax = "uu";
    if(max==err4) sMax = "dd";
    std::cout << "Maximum: " << max << ", " << sMax << std::endl;
    return max;
}

void getAsymmetry(TH1F* hAsymm, TH1F* h1,TH1F* h2){

    TH1F* top = (TH1F*)hAsymm->Clone();
    TH1F* bottom = (TH1F*)hAsymm->Clone();
    top->Add(h1,h2,1.0,-1.0);
    bottom->Add(h1,h2,1.0,1.0);
    hAsymm->Divide(top,bottom);
    Int_t xmax = hAsymm->GetNbinsX();
    Double_t a,b,bot,error,da,db;
    for(Int_t i = 1; i<=xmax; ++i){
        
        a = h1->GetBinContent(i);
        b = h2->GetBinContent(i);
        bot = bottom->GetBinContent(i);
        if(bot<1e-6){}
        else{
            
            da = h1->GetBinError(i);
            db = h2->GetBinError(i);
            error = 2*TMath::Sqrt(a*a*db*db + b*b*da*da )/(bot*bot);
            hAsymm->SetBinError(i,error);
        }
    }

    delete top;
    delete bottom;
    
}

void fillContainers(TFile* fout, TString sFileIn, TH1F* hEtaPos, TH1F* hEtaNeg, TH1F* hAsymm, bool doAsymm=false)
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

    //const int nBins = std::count (std::istreambuf_iterator<char>(s), std::istreambuf_iterator<char>(), '\n' );
    const int nBins = 10;
    std:: cout << "Number of eta bins: " << nBins << std::endl;
    std::vector <double> _vEtaBins(nBins);
    std::vector <double> _vEtaMuPos(nBins), _vEtaMuNeg(nBins), _vAsymmMu(nBins);
    std::vector <double> _vEtaMuPosError(nBins), _vEtaMuNegError(nBins), _vAsymmMuError(nBins);

	for (unsigned int i = 0; i < _vEtaBins.size(); ++i){
        //double input1,input2,input3,input4,input5,input6,input7;

        if(doAsymm) s >> _vEtaBins[i] >> _vEtaMuPos[i] >> _vEtaMuNeg[i] ; 
        else s >> _vEtaBins[i] >> _vEtaMuPos[i] >> _vEtaMuPosError[i] >> _vEtaMuNeg[i] >> _vEtaMuNegError[i] >> _vAsymmMu[i] >> _vAsymmMuError[i] ; 
        //s >> input1 >> input2 >> input3>>input4>>input5>>input6>>input7;
        std::cout << _vEtaBins[i] << std::endl;
    }


    for (unsigned int i = 0; i < nBins; ++i){
        hEtaPos->SetBinContent(i+1,_vEtaMuPos[i]); hEtaPos->SetBinError(i+1,_vEtaMuPosError[i]);
        hEtaNeg->SetBinContent(i+1,_vEtaMuNeg[i]); hEtaNeg->SetBinError(i+1,_vEtaMuNegError[i]);
        if(!doAsymm) hAsymm->SetBinContent(i+1,_vAsymmMu[i]); hAsymm->SetBinError(i+1,_vAsymmMuError[i]);
    }

    if(doAsymm) getAsymmetry(hAsymm,hEtaPos,hEtaNeg); 

/*    fout->cd();
    hEtaPos->Write();
    hEtaNeg->Write();
    hAsymm->Write();
*/
    s.close();
}

void EPS09_plotter(){

    TString sFileOut = "EPS09.root";
    TFile* fOut = new TFile(sFileOut,"recreate");
    TString sFileIn = "prediction.dat";
    TString sFileInUU = "prediction_scaleUU.dat";
    TString sFileInDD = "prediction_scaleDD.dat";
    TH1F *hEtaPos,*hEtaNeg,*hAsymm; 
    TH1F *hEtaPosUU,*hEtaNegUU,*hAsymmUU; 
    TH1F *hEtaPosDD,*hEtaNegDD,*hAsymmDD; 

    const int nBins = 10;
    double arrBin[nBins+1] = {0.0,0.1,0.35,0.6,0.8,1.05,1.37,1.52,1.85,2.01,2.5};
    hEtaPos = new TH1F("hEtaPos","mu^{+} EPS09", nBins, arrBin);
    hEtaNeg = new TH1F("hEtaNeg","mu^{-} EPS09", nBins, arrBin);
    hAsymm = new TH1F("hAsymm","asymmetry EPS09", nBins, arrBin);

    hEtaPosUU = new TH1F("hEtaPosUU","mu^{+} EPS09 scaleUp", nBins, arrBin);
    hEtaNegUU = new TH1F("hEtaNegUU","mu^{-} EPS09 scaleUp", nBins, arrBin);
    hAsymmUU = new TH1F("hAsymmUU","asymmetry EPS09 scaleUp", nBins, arrBin);

    hEtaPosDD = new TH1F("hEtaPosDD","mu^{+} EPS09 scaleDown", nBins, arrBin);
    hEtaNegDD = new TH1F("hEtaNegDD","mu^{-} EPS09 scaleDown", nBins, arrBin);
    hAsymmDD = new TH1F("hAsymmDD","asymmetry EPS09 scaleDown", nBins, arrBin);

    fillContainers(fOut,sFileIn,hEtaPos,hEtaNeg,hAsymm);
    fillContainers(fOut,sFileInUU,hEtaPosUU,hEtaNegUU,hAsymmUU,true);
    fillContainers(fOut,sFileInDD,hEtaPosDD,hEtaNegDD,hAsymmDD,true);

    // set errors 
    for(int ibin=1; ibin<=hAsymm->GetNbinsX(); ++ibin){

        double errorPos = 0.0;
        double errorPosUU = fabs(hEtaPosUU->GetBinContent(ibin)-hEtaPos->GetBinContent(ibin));
        double errorPosDD = fabs(hEtaPosDD->GetBinContent(ibin)-hEtaPos->GetBinContent(ibin));
        errorPos = findMaxVariation(0.0,0.0,errorPosUU,errorPosDD);
        hEtaPos->SetBinError(ibin,errorPos);

        double errorNeg = 0.0;
        double errorNegUU = fabs(hEtaNegUU->GetBinContent(ibin)-hEtaNeg->GetBinContent(ibin));
        double errorNegDD = fabs(hEtaNegDD->GetBinContent(ibin)-hEtaNeg->GetBinContent(ibin));
        errorNeg = findMaxVariation(0.0,0.0,errorNegUU,errorNegDD);
        hEtaNeg->SetBinError(ibin,errorNeg);

        double errorAsymm = 0.0;
        double errorAsymmUU = fabs(hAsymmUU->GetBinContent(ibin)-hAsymm->GetBinContent(ibin));
        double errorAsymmDD = fabs(hAsymmDD->GetBinContent(ibin)-hAsymm->GetBinContent(ibin));
        errorAsymm = findMaxVariation(0.0,0.0,errorAsymmUU,errorAsymmDD);
        hAsymm->SetBinError(ibin,errorAsymm);
    }//ibin
    fOut->cd();
    fOut->Write();
    fOut->Close();
}
