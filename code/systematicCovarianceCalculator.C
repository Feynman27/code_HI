#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"

#include <cmath>
#include <iostream>

///Print Covariance matrix 
void calculateCovMatrix(TGraphAsymmErrors* hNominal, TGraphAsymmErrors* hSyst, TH2F* hCorr){

    if(hNominal->GetN()!=hSyst->GetN()){
        std::cout << "ERROR: Number of bins in systematic and nominal histograms do not agree." << std::endl;
        return;
    }

    ///Loop over bin contents of hSyst
    int nBins = hSyst->GetN();
/*    double binDiff[nBins] = {-9999};
    for(int ibin=0; ibin<nBins; ++ibin){
        
        ///Get the difference between systematic and nominal bin contents
        binDiff[ibin] = hSyst->GetY()[ibin] - hNominal->GetY()[ibin];
    }
*/
    ///Calculate the bin correlations
    for(int ibin = 0; ibin<nBins; ++ibin){

        double x1 = hSyst->GetY()[ibin];
        double mu1 =  hNominal->GetY()[ibin];
        double sigma1 = hSyst->GetEYhigh()[ibin];

        //double residual1 = (x1-mu1)/sigma1;
        double residual1 = (x1-mu1);
        ///Get matrix element V_ibin,ibin2
        for(int ibin2 = 0; ibin2<nBins; ++ibin2){
           double x2 = hSyst->GetY()[ibin2]; 
           double mu2 = hNominal->GetY()[ibin2];
           double sigma2 = hSyst->GetEYhigh()[ibin2];

           //double residual2 = (x2-mu2)/sigma2;
           double residual2 = (x2-mu2);

           double correlation = residual1/residual2;
           std::cout << "Correlation matrix element at position (" << ibin << "," << ibin2 << ") = " << correlation << std::endl;
           hCorr->Fill(residual1,residual2);
        }
    }
}

void systematicCovarianceCalculator(){

    ///Nominal distributions with stat errors only
    TFile* fEtaNominal = new TFile("systematics/signalChargeEtaDistributions_NominalStatErrOnly.07.14.2013.root","read");
    TFile* fCentSyst = new TFile("systematics/binaryScalingDistribution_NominalStatErrOnly.07.14.2013.root","read");

    ///Uncomment systematic of interest
    ///Widen isolation cone
//    TFile* fEtaSyst = new TFile("systematics/signalChargeEtaDistributions_IncIsoConeSizeSyst.07.14.2013.root","read");
//    TFile* fCentNominal = new TFile("systematics/binaryScalingDistribution_IncIsoConeSizeSyst.07.14.2013.root","read");

    ///Loosen isolation cut
//    TFile* fEtaSyst = new TFile("systematics/signalChargeEtaDistributions_LoosenIsoCut.07.14.2013.root","read");
//    TFile* fCentNominal = new TFile("systematics/binaryScalingDistribution_LoosenIsoCut.07.14.2013.root","read");

    ///smear mpt and mt
    TFile* fEtaSyst = new TFile("systematics/signalChargeEtaDistributions_MptSmearSystematics.07.16.2013.root","read");
    TFile* fCentNominal = new TFile("systematics/binaryScalingDistribution_MptSmearSystematics.07.16.2013.root","read");

    if(!(fEtaNominal&&fEtaSyst&&fCentNominal&&fCentSyst)) {
        std::cout << "ERROR OPENING FILE. " << std::endl; exit(0);
    }
    else std::cout << "All files open." << std::endl;

    TGraphAsymmErrors* hEtaNominalPlus = (TGraphAsymmErrors*)fEtaNominal->Get("grWpSystc");
    TGraphAsymmErrors* hEtaSystPlus = (TGraphAsymmErrors*)fEtaSyst->Get("grWpSystc");
    int nEtaBins = hEtaNominalPlus->GetN();
    double etaLo = hEtaNominalPlus->GetX()[0];
    double binWLo = hEtaNominalPlus->GetErrorX(0);
    double etaUpp = hEtaNominalPlus->GetX()[nEtaBins-1];
    double binWUpp = hEtaNominalPlus->GetErrorX(nEtaBins-1);
    std::cout << "Number of eta bins: " << nEtaBins << " ranging from " << etaLo-binWLo << "-" << etaUpp+binWUpp << std::endl;
    TH2F* hEtaCorrPlus = new TH2F("hEtaCorrPlus","hEtaCorrPlus",60,-1.0,1.0,60,-1.0,1.0);

    TGraphAsymmErrors* hEtaNominalMinus = (TGraphAsymmErrors*)fEtaNominal->Get("grWmSystc");
    TGraphAsymmErrors* hEtaSystMinus = (TGraphAsymmErrors*)fEtaSyst->Get("grWmSystc");
    TH2F* hEtaCorrMinus = new TH2F("hEtaCorrMinus","hEtaCorrMinus",nEtaBins,0.0,etaUpp+binWUpp,nEtaBins,0.0,etaUpp+binWUpp);

    TGraphAsymmErrors* hCentNominal = (TGraphAsymmErrors*)fCentNominal->Get("RcpChargeInclusive");
    TGraphAsymmErrors* hCentSyst = (TGraphAsymmErrors*)fCentSyst->Get("RcpChargeInclusive");
    int nCentBins = hCentNominal->GetN();
    double centLo = hCentNominal->GetX()[0];
    binWLo = hCentNominal->GetErrorX(0);
    double centUpp = hCentNominal->GetX()[nCentBins-1];
    binWUpp = hCentNominal->GetErrorX(nCentBins-1);
    std::cout << "Number of centrality bins: " << nCentBins << " ranging from " << centLo-binWLo << "-" << centUpp+binWUpp << std::endl;
    TH2F* hCentCorr = new TH2F("hCentCorr","hCentCorr",nCentBins,0.0,centUpp+binWUpp,nCentBins,0.0,centUpp+binWUpp); 

    std::cout << "Calculate systematic correlations in eta for mu+." << std::endl;
    calculateCovMatrix(hEtaNominalPlus,hEtaSystPlus,hEtaCorrPlus);
    gStyle->SetPalette(1);

    TCanvas* cEtaPlus = new TCanvas("cEtaPlus","cEtaPlus",600,600);
    TCanvas* cEtaMinus = new TCanvas("cEtaMinus","cEtaMinus",600,600);
    TCanvas* cCent = new TCanvas("cCent","cCent",600,600);

    cEtaPlus->cd();
    hEtaCorrPlus->Draw("colz");
    //cEtaPlus.Print("systematicIsoEtaCorrelation.root");

    std::cout << "Calculate systematic correlations in eta for mu-." << std::endl;
    calculateCovMatrix(hEtaNominalMinus,hEtaSystMinus,hEtaCorrMinus);
    cEtaMinus->cd();
    hEtaCorrMinus->Draw("colz");

    std::cout << "Calculate systematic correlations in Npart." << std::endl;
    calculateCovMatrix(hCentNominal,hCentSyst,hCentCorr);
    cCent->cd();
    hCentCorr->Draw("colz");

}
