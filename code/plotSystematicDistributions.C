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
#include "TPad.h"
#include "TLegendEntry.h"

#include "AtlasUtils.C"
#include "AtlasStyle.C"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath> 
#include <sstream>

////////////////////////
//Prevent trailing digits
//in figure labeling
////////////////////////
TString format(float value) {
  std::stringstream svalue;
  svalue  << std::setprecision(3) << value;
  return svalue.str();
}
///////////////////////////////////////////
//writeToSpreadsheet
//////////////////////////////////////////
void writeToSpreadSheet(std::ostream& outputFile, TString range, TString variation1, TString variation2,TString variation3,TString variation4){
	outputFile << range << "," << variation1<< ","  << variation2<< ","  << variation3 << ","  << variation4 << std::endl;
}

void plotSystematicDistributions(){

    TDatime* time = new TDatime();
    TString sDate = "";
    sDate+=time->GetMonth(); sDate+="_"; sDate+=time->GetDay(); sDate+="_"; sDate+=time->GetYear();
    std::cout << "Today is " << sDate << std::endl;

    TString sOutFile = "systematicDistributions_"; sOutFile+=sDate; sOutFile+=".root";
    std::cout << "Writing files to " << sOutFile << std::endl;
    TFile* outFile = new TFile(sOutFile,"recreate");

    TString ssAsymmName = "variationAsymmetry.csv";
    std::ofstream ssAsymm;
    ssAsymm.open(ssAsymmName);

    TString ssCentName = "variationCentrality.csv";
    std::ofstream ssCent;
    ssCent.open(ssCentName);

    TString sEtaRange[] = {"0.1-0.35","0.35-0.6","0.6-0.8","0.8-1.05","1.05-1.3","1.3-1.55","1.55-1.85","1.85-2.1","2.1-2.4"};
    TString sCentralityRange[] = {"0-5","5-10","10-15","15-20","20-40","40-80"};

    TFile* _fEtaNominal = TFile::Open("systematics/signalChargeEtaDistributions_NominalStatErrOnly.07.22.2013.root");
    TFile* _fCentNominal = TFile::Open("systematics/binaryScalingDistribution_NominalStatErrOnly.07.22.2013.root");
    ///includes yields with Wtau bkg subtracted
    TFile* _fEtaNominal2 = TFile::Open("systematics/signalChargeEtaDistributions_NominalWithTauBkg.07.30.2013.root");
    TFile* _fCentNominal2 = TFile::Open("systematics/binaryScalingDistribution_NominalWithTauBkg.07.30.2013.root");

    TFile* _fEtaMptSigUp = TFile::Open("systematics/signalChargeEtaDistributions_MptSigmaUp.08.04.2013.root");
    TFile* _fCentMptSigUp = TFile::Open("systematics/binaryScalingDistribution_MptSigmaUp.08.04.2013.root");
    TFile* _fEtaMptSigDown = TFile::Open("systematics/signalChargeEtaDistributions_MptSigmaDown.08.04.2013.root");
    TFile* _fCentMptSigDown = TFile::Open("systematics/binaryScalingDistribution_MptSigmaDown.08.04.2013.root");
    TFile* _fEtaIncCone = TFile::Open("systematics/signalChargeEtaDistributions_IncIsoConeSizeSyst.07.22.2013.root");
    TFile* _fCentIncCone = TFile::Open("systematics/binaryScalingDistribution_IncIsoConeSizeSyst.07.22.2013.root");
    TFile* _fEtaLoosenIso = TFile::Open("systematics/signalChargeEtaDistributions_LoosenIsoCut.07.22.2013.root");
    TFile* _fCentLoosenIso = TFile::Open("systematics/binaryScalingDistribution_LoosenIsoCut.07.22.2013.root");
    ///includes yields with Wtau bkg subtracted
    ///Compare to Nominal v2
    TFile* _fEtaRaaScaledQCD = TFile::Open("systematics/signalChargeEtaDistributions_RaaScaledQCD.07.30.2013.root");
    TFile* _fEtaZBkg = TFile::Open("systematics/signalChargeEtaDistributions_ZBkg.07.30.2013.root");
    TFile* _fEtaTauBkg = TFile::Open("systematics/signalChargeEtaDistributions_TauBkg.07.30.2013.root");
    TFile* _fCentRaaScaledQCD = TFile::Open("systematics/binaryScalingDistribution_RaaScaledQCD.07.30.2013.root");
    TFile* _fCentZBkg = TFile::Open("systematics/binaryScalingDistribution_ZBkg.07.30.2013.root");
    TFile* _fCentTauBkg = TFile::Open("systematics/binaryScalingDistribution_TauBkg.07.30.2013.root");

    TGraphAsymmErrors *grEtaPlusNominal = (TGraphAsymmErrors*)_fEtaNominal->Get("grWpc");
    const int nEtaBins = grEtaPlusNominal->GetN();
    TGraphAsymmErrors *grEtaPlusMptSigUp = (TGraphAsymmErrors*)_fEtaMptSigUp->Get("grWpc");
    TGraphAsymmErrors *grEtaPlusMptSigDown = (TGraphAsymmErrors*)_fEtaMptSigDown->Get("grWpc");
    TGraphAsymmErrors *grEtaPlusIncCone = (TGraphAsymmErrors*)_fEtaIncCone->Get("grWpc");
    TGraphAsymmErrors *grEtaPlusLoosenIso = (TGraphAsymmErrors*)_fEtaLoosenIso->Get("grWpc");
    TGraphAsymmErrors *grEtaPlusNominal2 = (TGraphAsymmErrors*)_fEtaNominal2->Get("grWpc");
    TGraphAsymmErrors *grEtaPlusRaaScaledQCD = (TGraphAsymmErrors*)_fEtaRaaScaledQCD->Get("grWpc");
    TGraphAsymmErrors *grEtaPlusZBkg = (TGraphAsymmErrors*)_fEtaZBkg->Get("grWpc");
    TGraphAsymmErrors *grEtaPlusTauBkg = (TGraphAsymmErrors*)_fEtaTauBkg->Get("grWpc");

    TGraphAsymmErrors *grEtaMinusNominal = (TGraphAsymmErrors*)_fEtaNominal->Get("grWmc");
    TGraphAsymmErrors *grEtaMinusMptSigUp = (TGraphAsymmErrors*)_fEtaMptSigUp->Get("grWmc");
    TGraphAsymmErrors *grEtaMinusMptSigDown = (TGraphAsymmErrors*)_fEtaMptSigDown->Get("grWmc");
    TGraphAsymmErrors *grEtaMinusIncCone = (TGraphAsymmErrors*)_fEtaIncCone->Get("grWmc");
    TGraphAsymmErrors *grEtaMinusLoosenIso = (TGraphAsymmErrors*)_fEtaLoosenIso->Get("grWmc");
    TGraphAsymmErrors *grEtaMinusNominal2 = (TGraphAsymmErrors*)_fEtaNominal2->Get("grWmc");
    TGraphAsymmErrors *grEtaMinusRaaScaledQCD = (TGraphAsymmErrors*)_fEtaRaaScaledQCD->Get("grWmc");
    TGraphAsymmErrors *grEtaMinusZBkg = (TGraphAsymmErrors*)_fEtaZBkg->Get("grWmc");
    TGraphAsymmErrors *grEtaMinusTauBkg = (TGraphAsymmErrors*)_fEtaTauBkg->Get("grWmc");

    TGraphAsymmErrors *grAsymmNominal = (TGraphAsymmErrors*)_fEtaNominal->Get("grChargeAsymm");
    TGraphAsymmErrors *grAsymmMptSigUp = (TGraphAsymmErrors*)_fEtaMptSigUp->Get("grChargeAsymm");
    TGraphAsymmErrors *grAsymmMptSigDown = (TGraphAsymmErrors*)_fEtaMptSigDown->Get("grChargeAsymm");
    TGraphAsymmErrors *grAsymmIncCone = (TGraphAsymmErrors*)_fEtaIncCone->Get("grChargeAsymm");
    TGraphAsymmErrors *grAsymmLoosenIso = (TGraphAsymmErrors*)_fEtaLoosenIso->Get("grChargeAsymm");
    TGraphAsymmErrors *grAsymmNominal2 = (TGraphAsymmErrors*)_fEtaNominal2->Get("grChargeAsymm");
    TGraphAsymmErrors *grAsymmRaaScaledQCD = (TGraphAsymmErrors*)_fEtaRaaScaledQCD->Get("grChargeAsymm");
    TGraphAsymmErrors *grAsymmZBkg = (TGraphAsymmErrors*)_fEtaZBkg->Get("grChargeAsymm");
    TGraphAsymmErrors *grAsymmTauBkg = (TGraphAsymmErrors*)_fEtaTauBkg->Get("grChargeAsymm");

    double arrEta[] = {0.1,0.35,0.6,0.8,1.05,1.3,1.55,1.85,2.1,2.4};

    TH1F *hEtaPlusNominal = new TH1F("hEtaPlusNominal","hEtaPlusNominal",nEtaBins,arrEta);
    TH1F *hEtaPlusMptSigUp = new TH1F("hEtaPlusMptSigUp","hEtaPlusMptSigUp",nEtaBins,arrEta);
    TH1F *hEtaPlusMptSigDown = new TH1F("hEtaPlusMptSigDown","hEtaPlusMptSigDown",nEtaBins,arrEta);
    TH1F *hEtaPlusIncCone = new TH1F("hEtaPlusIncCone","hEtaPlusIncCone",nEtaBins,arrEta);
    TH1F *hEtaPlusLoosenIso = new TH1F("hEtaPlusLoosenIso","hEtaPlusLoosenIso",nEtaBins,arrEta);
    TH1F *hEtaPlusNominal2 = new TH1F("hEtaPlusNominal2","hEtaPlusNominal2",nEtaBins,arrEta);
    TH1F *hEtaPlusRaaScaledQCD = new TH1F("hEtaPlusRaaScaledQCD","hEtaPlusRaaScaledQCD",nEtaBins,arrEta);
    TH1F *hEtaPlusZBkg = new TH1F("hEtaPlusZBkg","hEtaPlusZBkg",nEtaBins,arrEta);
    TH1F *hEtaPlusTauBkg = new TH1F("hEtaPlusTauBkg","hEtaPlusTauBkg",nEtaBins,arrEta);

    TH1F *hEtaMinusNominal = new TH1F("hEtaMinusNominal","hEtaMinusNominal",nEtaBins,arrEta);
    TH1F *hEtaMinusMptSigUp = new TH1F("hEtaMinusMptSigUp","hEtaMinusMptSigUp",nEtaBins,arrEta);
    TH1F *hEtaMinusMptSigDown = new TH1F("hEtaMinusMptSigDown","hEtaMinusMptSigDown",nEtaBins,arrEta);
    TH1F *hEtaMinusIncCone = new TH1F("hEtaMinusIncCone","hEtaMinusIncCone",nEtaBins,arrEta);
    TH1F *hEtaMinusLoosenIso = new TH1F("hEtaMinusLoosenIso","hEtaMinusLoosenIso",nEtaBins,arrEta);
    TH1F *hEtaMinusNominal2 = new TH1F("hEtaMinusNominal2","hEtaMinusNominal2",nEtaBins,arrEta);
    TH1F *hEtaMinusRaaScaledQCD = new TH1F("hEtaMinusRaaScaledQCD","hEtaMinusRaaScaledQCD",nEtaBins,arrEta);
    TH1F *hEtaMinusZBkg = new TH1F("hEtaMinusZBkg","hEtaMinusZBkg",nEtaBins,arrEta);
    TH1F *hEtaMinusTauBkg = new TH1F("hEtaMinusTauBkg","hEtaMinusTauBkg",nEtaBins,arrEta);

    TH1F *hAsymmNominal = new TH1F("hAsymmNominal","hAsymmNominal",nEtaBins,arrEta);
    TH1F *hAsymmMptSigUp = new TH1F("hAsymmMptSigUp","hAsymmMptSigUp",nEtaBins,arrEta);
    TH1F *hAsymmMptSigDown = new TH1F("hAsymmMptSigDown","hAsymmMptSigDown",nEtaBins,arrEta);
    TH1F *hAsymmIncCone = new TH1F("hAsymmIncCone","hAsymmIncCone",nEtaBins,arrEta);
    TH1F *hAsymmLoosenIso = new TH1F("hAsymmLoosenIso","hAsymmLoosenIso",nEtaBins,arrEta);
    TH1F *hAsymmNominal2 = new TH1F("hAsymmNominal2","hAsymmNominal2",nEtaBins,arrEta);
    TH1F *hAsymmRaaScaledQCD = new TH1F("hAsymmRaaScaledQCD","hAsymmRaaScaledQCD",nEtaBins,arrEta);
    TH1F *hAsymmZBkg = new TH1F("hAsymmZBkg","hAsymmZBkg",nEtaBins,arrEta);
    TH1F *hAsymmTauBkg = new TH1F("hAsymmTauBkg","hAsymmTauBkg",nEtaBins,arrEta);



//    TGraphAsymmErrors *grAsymmMptAvg = new TGraphAsymmErrors(nEtaBins);
//    TGraphAsymmErrors *grAsymmIsoAvg = new TGraphAsymmErrors(nEtaBins);

    TGraphAsymmErrors *grEtaDiffMpt = new TGraphAsymmErrors(nEtaBins);
    TGraphAsymmErrors *grEtaDiffIso = new TGraphAsymmErrors(nEtaBins);
    TGraphAsymmErrors *grEtaDiffRaaScaledQCD = new TGraphAsymmErrors(nEtaBins);
    TGraphAsymmErrors *grEtaDiffZBkg = new TGraphAsymmErrors(nEtaBins);
    TGraphAsymmErrors *grEtaDiffTauBkg = new TGraphAsymmErrors(nEtaBins);

    for(int i=0; i<nEtaBins; ++i){

        hEtaPlusMptSigUp->SetBinContent(i+1,grEtaPlusMptSigUp->GetY()[i]);
        hEtaPlusMptSigDown->SetBinContent(i+1,grEtaPlusMptSigDown->GetY()[i]);
        hEtaPlusIncCone->SetBinContent(i+1,grEtaPlusIncCone->GetY()[i]);
        hEtaPlusLoosenIso->SetBinContent(i+1,grEtaPlusLoosenIso->GetY()[i]);
        hEtaPlusRaaScaledQCD->SetBinContent(i+1,grEtaPlusRaaScaledQCD->GetY()[i]);
        hEtaPlusZBkg->SetBinContent(i+1,grEtaPlusZBkg->GetY()[i]);
        hEtaPlusTauBkg->SetBinContent(i+1,grEtaPlusTauBkg->GetY()[i]);

        hEtaMinusMptSigUp->SetBinContent(i+1,grEtaMinusMptSigUp->GetY()[i]);
        hEtaMinusMptSigDown->SetBinContent(i+1,grEtaMinusMptSigDown->GetY()[i]);
        hEtaMinusIncCone->SetBinContent(i+1,grEtaMinusIncCone->GetY()[i]);
        hEtaMinusLoosenIso->SetBinContent(i+1,grEtaMinusLoosenIso->GetY()[i]);
        hEtaMinusRaaScaledQCD->SetBinContent(i+1,grEtaMinusRaaScaledQCD->GetY()[i]);
        hEtaMinusZBkg->SetBinContent(i+1,grEtaMinusZBkg->GetY()[i]);
        hEtaMinusTauBkg->SetBinContent(i+1,grEtaMinusTauBkg->GetY()[i]);

        hAsymmMptSigUp->SetBinContent(i+1,grAsymmMptSigUp->GetY()[i]);
        hAsymmMptSigDown->SetBinContent(i+1,grAsymmMptSigDown->GetY()[i]);
        hAsymmIncCone->SetBinContent(i+1,grAsymmIncCone->GetY()[i]);
        hAsymmLoosenIso->SetBinContent(i+1,grAsymmLoosenIso->GetY()[i]);
        hAsymmRaaScaledQCD->SetBinContent(i+1,grAsymmRaaScaledQCD->GetY()[i]);
        hAsymmZBkg->SetBinContent(i+1,grAsymmZBkg->GetY()[i]);
        hAsymmTauBkg->SetBinContent(i+1,grAsymmTauBkg->GetY()[i]);
        
        hEtaPlusNominal->SetBinContent(i+1,grEtaPlusNominal->GetY()[i]); 
        hEtaMinusNominal->SetBinContent(i+1,grEtaMinusNominal->GetY()[i]); 
        hAsymmNominal->SetBinContent(i+1,grAsymmNominal->GetY()[i]); 
        hEtaPlusNominal2->SetBinContent(i+1,grEtaPlusNominal2->GetY()[i]); 
        hEtaMinusNominal2->SetBinContent(i+1,grEtaMinusNominal2->GetY()[i]); 
        hAsymmNominal2->SetBinContent(i+1,grAsymmNominal2->GetY()[i]); 


        hEtaPlusMptSigUp->SetBinError(i+1,0.0);
        hEtaPlusMptSigDown->SetBinError(i+1,0.0);
        hEtaPlusIncCone->SetBinError(i+1,0.0);
        hEtaPlusLoosenIso->SetBinError(i+1,0.0);
        hEtaPlusRaaScaledQCD->SetBinError(i+1,0.0);
        hEtaPlusZBkg->SetBinError(i+1,0.0);
        hEtaPlusTauBkg->SetBinError(i+1,0.0);

        hEtaMinusMptSigUp->SetBinError(i+1,0.0);
        hEtaMinusMptSigDown->SetBinError(i+1,0.0);
        hEtaMinusIncCone->SetBinError(i+1,0.0);
        hEtaMinusLoosenIso->SetBinError(i+1,0.0);
        hEtaMinusRaaScaledQCD->SetBinError(i+1,0.0);
        hEtaMinusZBkg->SetBinError(i+1,0.0);
        hEtaMinusTauBkg->SetBinError(i+1,0.0);

        hAsymmMptSigUp->SetBinError(i+1,0.0);
        hAsymmMptSigDown->SetBinError(i+1,0.0);
        hAsymmIncCone->SetBinError(i+1,0.0);
        hAsymmLoosenIso->SetBinError(i+1,0.0);
        hAsymmRaaScaledQCD->SetBinError(i+1,0.0);
        hAsymmZBkg->SetBinError(i+1,0.0);
        hAsymmTauBkg->SetBinError(i+1,0.0);

        
        hEtaPlusNominal->SetBinError(i+1,grEtaPlusNominal->GetEYhigh()[i]); 
        hEtaMinusNominal->SetBinError(i+1,grEtaMinusNominal->GetEYhigh()[i]); 
        hAsymmNominal->SetBinError(i+1,grAsymmNominal->GetEYhigh()[i]); 
        hEtaPlusNominal2->SetBinError(i+1,grEtaPlusNominal2->GetEYhigh()[i]); 
        hEtaMinusNominal2->SetBinError(i+1,grEtaMinusNominal2->GetEYhigh()[i]); 
        hAsymmNominal2->SetBinError(i+1,grAsymmNominal2->GetEYhigh()[i]); 

        ///Fill TGraph with average values
        /*double mptAsymmAvg = (grAsymmMptSigUp->GetY()[i]+grAsymmMptSigDown->GetY()[i])/2.0;
        double isoAsymmAvg = (grAsymmIncCone->GetY()[i]+grAsymmLoosenIso->GetY()[i])/2.0;
        grAsymmMptAvg->SetPoint(i,xEta,mptAsymmAvg); grAsymmMptAvg->SetPointError(i,0.0,0.0,0.0,0.0);
        grAsymmIsoAvg->SetPoint(i,xEta,isoAsymmAvg); grAsymmIsoAvg->SetPointError(i,0.0,0.0,0.0,0.0);
        */

        ///Calculate variations
        double xEta = grAsymmMptSigUp->GetX()[i];
        ///Mpt
        double diffSigUp = fabs(grAsymmMptSigUp->GetY()[i]-grAsymmNominal2->GetY()[i]);
        double diffSigDown = fabs(grAsymmMptSigDown->GetY()[i]-grAsymmNominal2->GetY()[i]);
        double avgDiffPercentage1 = (diffSigUp+diffSigDown)/2.0;  
//        avgDiffPercentage1/=fabs(grAsymmNominal2->GetY()[i]);
        grEtaDiffMpt->SetPoint(i,xEta,avgDiffPercentage1); grEtaDiffMpt->SetPointError(i,0.0,0.0,0.0,0.0);

        ///Isolation efficiency
        diffSigUp= fabs(grAsymmIncCone->GetY()[i]-grAsymmNominal->GetY()[i]);
        diffSigDown= fabs(grAsymmLoosenIso->GetY()[i]-grAsymmNominal->GetY()[i]);
        double avgDiffPercentage2 = (diffSigUp+diffSigDown)/2.0; 
 //       avgDiffPercentage2/=fabs(grAsymmNominal->GetY()[i]);
        grEtaDiffIso->SetPoint(i,xEta,avgDiffPercentage2); grEtaDiffIso->SetPointError(i,0.0,0.0,0.0,0.0);

        ///QCD
        double diffSig1=0.0,diffSig2=0.0,diffSig3=0.0;
        diffSig1= fabs(grAsymmRaaScaledQCD->GetY()[i]-grAsymmNominal2->GetY()[i]); 
//        diffSig1/=fabs(grAsymmNominal2->GetY()[i]); 
        grEtaDiffRaaScaledQCD->SetPoint(i,xEta,diffSig1); grEtaDiffRaaScaledQCD->SetPointError(i,0.0,0.0,0.0,0.0);

        ///Z
        diffSig2= fabs(grAsymmZBkg->GetY()[i]-grAsymmNominal2->GetY()[i]); 
//        diffSig2/=fabs(grAsymmNominal2->GetY()[i]);
        grEtaDiffZBkg->SetPoint(i,xEta,diffSig2); grEtaDiffZBkg->SetPointError(i,0.0,0.0,0.0,0.0);

        ///Tau
        diffSig3= fabs(grAsymmTauBkg->GetY()[i]-grAsymmNominal2->GetY()[i]); 
//        diffSig3/=fabs(grAsymmNominal2->GetY()[i]);
        grEtaDiffTauBkg->SetPoint(i,xEta,diffSig3); grEtaDiffTauBkg->SetPointError(i,0.0,0.0,0.0,0.0);
        
        double avgDiffPercentage3 = (diffSig2+diffSig3)/2.0; 
        avgDiffPercentage3/=fabs(grAsymmNominal2->GetY()[i]); 

        writeToSpreadSheet(ssAsymm,sEtaRange[i],format(100.*avgDiffPercentage1),format(100.*avgDiffPercentage2),format(100.*diffSig1),format(100.*avgDiffPercentage3));

    }


    TGraphAsymmErrors *grCentNominal = (TGraphAsymmErrors*)_fCentNominal->Get("RcpChargeInclusive");
    TGraphAsymmErrors *grCentMptSigUp = (TGraphAsymmErrors*)_fCentMptSigUp->Get("RcpChargeInclusive");
    TGraphAsymmErrors *grCentMptSigDown = (TGraphAsymmErrors*)_fCentMptSigDown->Get("RcpChargeInclusive");
    TGraphAsymmErrors *grCentIncCone = (TGraphAsymmErrors*)_fCentIncCone->Get("RcpChargeInclusive");
    TGraphAsymmErrors *grCentLoosenIso = (TGraphAsymmErrors*)_fCentLoosenIso->Get("RcpChargeInclusive");
    TGraphAsymmErrors *grCentNominal2 = (TGraphAsymmErrors*)_fCentNominal2->Get("RcpChargeInclusive");
    TGraphAsymmErrors *grCentRaaScaledQCD = (TGraphAsymmErrors*)_fCentRaaScaledQCD->Get("RcpChargeInclusive");
    TGraphAsymmErrors *grCentZBkg = (TGraphAsymmErrors*)_fCentZBkg->Get("RcpChargeInclusive");
    TGraphAsymmErrors *grCentTauBkg = (TGraphAsymmErrors*)_fCentTauBkg->Get("RcpChargeInclusive");

    const int nCentBins = grCentNominal->GetN();
    double arrNpart[] = {45.93,157.83,239.52,281.88,330.26,382.16};
    double arrCentBin[6] ;
    arrCentBin[0] = 0.0;
    for(int i=0; i<nCentBins; ++i) {
        arrCentBin[i+1] = arrCentBin[i]+2.*(arrNpart[i] - arrCentBin[i]);
        std::cout << "Bins: " << arrCentBin[i] << "-" << arrCentBin[i+1] << std::endl;
    }


    TH1F* hCentNominal = new TH1F("hCentNominal","hCentNominal",nCentBins,arrCentBin);
    TH1F* hCentMptSigUp = new TH1F("hCentMptSigUp","hCentMptSigUp",nCentBins,arrCentBin);
    TH1F* hCentMptSigDown = new TH1F("hCentMptSigDown","hCentMptSigDown",nCentBins,arrCentBin);
    TH1F* hCentIncCone = new TH1F("hCentIncCone","hCentIncCone",nCentBins,arrCentBin);
    TH1F* hCentLoosenIso = new TH1F("hCentLoosenIso","hCentLoosenIso",nCentBins,arrCentBin);
    TH1F* hCentNominal2 = new TH1F("hCentNominal2","hCentNominal2",nCentBins,arrCentBin);
    TH1F* hCentRaaScaledQCD = new TH1F("hCentRaaScaledQCD","hCentRaaScaledQCD",nCentBins,arrCentBin);
    TH1F* hCentZBkg = new TH1F("hCentZBkg","hCentZBkg",nCentBins,arrCentBin);
    TH1F* hCentTauBkg = new TH1F("hCentTauBkg","hCentTauBkg",nCentBins,arrCentBin);

    TGraphAsymmErrors *grCentDiffMpt = new TGraphAsymmErrors(nEtaBins);
    TGraphAsymmErrors *grCentDiffIso = new TGraphAsymmErrors(nEtaBins);
    TGraphAsymmErrors *grCentDiffRaaScaledQCD = new TGraphAsymmErrors(nEtaBins);
    TGraphAsymmErrors *grCentDiffZBkg = new TGraphAsymmErrors(nEtaBins);
    TGraphAsymmErrors *grCentDiffTauBkg = new TGraphAsymmErrors(nEtaBins);

    for(int i=0; i<nCentBins; ++i){
        hCentNominal->SetBinContent(nCentBins-i,grCentNominal->GetY()[i]);
        hCentMptSigUp->SetBinContent(nCentBins-i,grCentMptSigUp->GetY()[i]);
        hCentMptSigDown->SetBinContent(nCentBins-i,grCentMptSigDown->GetY()[i]);
        hCentIncCone->SetBinContent(nCentBins-i,grCentIncCone->GetY()[i]);
        hCentLoosenIso->SetBinContent(nCentBins-i,grCentLoosenIso->GetY()[i]);
        hCentNominal2->SetBinContent(nCentBins-i,grCentNominal2->GetY()[i]);
        hCentRaaScaledQCD->SetBinContent(nCentBins-i,grCentRaaScaledQCD->GetY()[i]);
        hCentZBkg->SetBinContent(nCentBins-i,grCentZBkg->GetY()[i]);
        hCentTauBkg->SetBinContent(nCentBins-i,grCentTauBkg->GetY()[i]);

        hCentNominal->SetBinError(nCentBins-i,grCentNominal->GetEYhigh()[i]);
        hCentMptSigUp->SetBinError(nCentBins-i,0.0);
        hCentMptSigDown->SetBinError(nCentBins-i,0.0);
        hCentIncCone->SetBinError(nCentBins-i,0.0);
        hCentLoosenIso->SetBinError(nCentBins-i,0.0);
        hCentNominal2->SetBinError(nCentBins-i,grCentNominal2->GetEYhigh()[i]);
        hCentRaaScaledQCD->SetBinError(nCentBins-i,0.0);
        hCentZBkg->SetBinError(nCentBins-i,0.0);
        hCentTauBkg->SetBinError(nCentBins-i,0.0);

        ///Calculate variations
        double xCent = grCentNominal->GetX()[i];
        ///Mpt
        double diffSigUp = fabs(grCentMptSigUp->GetY()[i]-grCentNominal2->GetY()[i]);
        double diffSigDown = fabs(grCentMptSigDown->GetY()[i]-grCentNominal2->GetY()[i]);
        double avgDiffPercentage1 = (diffSigUp+diffSigDown)/2.0;  
//        avgDiffPercentage1/=grCentNominal2->GetY()[i];
        grCentDiffMpt->SetPoint(i,xCent,avgDiffPercentage1); grCentDiffMpt->SetPointError(i,0.0,0.0,0.0,0.0);

        ///Isolation efficiency
        diffSigUp= fabs(grCentIncCone->GetY()[i]-grCentNominal->GetY()[i]);
        diffSigDown= fabs(grCentLoosenIso->GetY()[i]-grCentNominal->GetY()[i]);
        double avgDiffPercentage2 = (diffSigUp+diffSigDown)/2.0; 
 //       avgDiffPercentage2/=grCentNominal->GetY()[i];
        grCentDiffIso->SetPoint(i,xCent,avgDiffPercentage2); grCentDiffIso->SetPointError(i,0.0,0.0,0.0,0.0);

        ///QCD
        double diffSig1=0.0,diffSig2=0.0,diffSig3=0.0;
        diffSig1= fabs(grCentRaaScaledQCD->GetY()[i]-grCentNominal2->GetY()[i]); 
//        diffSig1/=grCentNominal2->GetY()[i]; 
        grCentDiffRaaScaledQCD->SetPoint(i,xCent,diffSig1); grCentDiffRaaScaledQCD->SetPointError(i,0.0,0.0,0.0,0.0);

        ///Z
        diffSig2= fabs(grCentZBkg->GetY()[i]-grCentNominal2->GetY()[i]); 
//        diffSig2/=grCentNominal2->GetY()[i];
        grCentDiffZBkg->SetPoint(i,xCent,diffSig2); grCentDiffZBkg->SetPointError(i,0.0,0.0,0.0,0.0);

        ///Tau
        diffSig3= fabs(grCentTauBkg->GetY()[i]-grCentNominal2->GetY()[i]); 
//        diffSig3/=grCentNominal2->GetY()[i];
        grCentDiffTauBkg->SetPoint(i,xCent,diffSig3); grCentDiffTauBkg->SetPointError(i,0.0,0.0,0.0,0.0);
        
        double avgDiffPercentage3 = (diffSig2+diffSig3)/2.0; 
        avgDiffPercentage3/=grCentNominal2->GetY()[i]; 

        writeToSpreadSheet(ssCent,sCentralityRange[i],format(100.*avgDiffPercentage1),format(100.*avgDiffPercentage2),format(100.*diffSig1),format(100.*avgDiffPercentage3));
    } //i,nCentBins

   ///Eta distributions
   TCanvas *cAsymm1 = new TCanvas("cAsymm1", "cAsymm1",11,142,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cAsymm1->Range(-0.5307818,674.7044,2.797364,3793.449);
   cAsymm1->SetFillColor(0);
   cAsymm1->SetBorderMode(0);
   cAsymm1->SetBorderSize(2);
   cAsymm1->SetTickx(1);
   cAsymm1->SetTicky(1);
   cAsymm1->SetLeftMargin(0.1594828);
   cAsymm1->SetRightMargin(0.05028735);
   cAsymm1->SetTopMargin(0.04872881);
   cAsymm1->SetBottomMargin(0.161017);
   cAsymm1->SetFrameBorderMode(0);
   cAsymm1->SetFrameBorderMode(0);
 
   hAsymmNominal2->GetXaxis()->SetTitle("|#eta|");
   hAsymmNominal2->GetYaxis()->SetTitle("A_{#mu}");
   hAsymmNominal2->GetXaxis()->SetLabelFont(42);
   hAsymmNominal2->GetXaxis()->SetLabelSize(0.05);
   hAsymmNominal2->GetXaxis()->SetTitleSize(0.05);
   hAsymmNominal2->GetXaxis()->SetTitleOffset(1.4);
   hAsymmNominal2->GetXaxis()->SetTitleFont(42);
   hAsymmNominal2->GetYaxis()->SetLabelFont(42);
   hAsymmNominal2->GetYaxis()->SetLabelSize(0.05);
   hAsymmNominal2->GetYaxis()->SetTitleSize(0.05);
   hAsymmNominal2->GetYaxis()->SetTitleOffset(1.4);
   hAsymmNominal2->GetYaxis()->SetTitleFont(42);
   hAsymmNominal2->SetMarkerSize(1.4);
   hAsymmNominal2->GetYaxis()->SetRangeUser(-0.4,+0.3);
   hAsymmNominal2->Draw("hist");
   hAsymmNominal2->Draw("pesame");

   hAsymmMptSigUp->SetMarkerStyle(23);
   hAsymmMptSigUp->SetLineStyle(3);
   hAsymmMptSigUp->SetLineWidth(0);
   hAsymmMptSigUp->SetLineColor(kRed);
   hAsymmMptSigUp->SetMarkerColor(kRed);
   hAsymmMptSigUp->Draw("hist same");
   hAsymmMptSigUp->Draw("pe same");

   hAsymmMptSigDown->SetMarkerStyle(22);
   hAsymmMptSigDown->SetLineStyle(3);
   hAsymmMptSigDown->SetLineWidth(0);
   hAsymmMptSigDown->SetLineColor(kBlue);
   hAsymmMptSigDown->SetMarkerColor(kBlue);
   hAsymmMptSigDown->Draw("hist same");
   hAsymmMptSigDown->Draw("pe same");

/*   hAsymmMptAvg->SetMarkerStyle(24);
   hAsymmMptAvg->SetLineStyle(3);
   hAsymmMptAvg->SetLineWidth(0);
   hAsymmMptAvg->SetMarkerColor(kBlack);
   hAsymmMptAvg->Draw("hist e same");
*/
   TLegend* leg = new TLegend(0.234,0.235,0.5,0.49,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry* entry=leg->AddEntry("hAsymmNominal2","Nominal","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.4);
   entry=leg->AddEntry("hAsymmMptSigUp","#slash{p_{T}} 4GeV","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hAsymmMptSigDown","#slash{p_{T}} 2GeV","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(kBlue);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);

/*   entry=leg->AddEntry("hAsymmMptAvg","Average","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.2);
*/

   leg->Draw();

   ///Eta distros
   TCanvas *cEtaPlus1 = new TCanvas("cEtaPlus1", "cEtaPlus1",11,142,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cEtaPlus1->Range(-0.5307818,674.7044,2.797364,3793.449);
   cEtaPlus1->SetFillColor(0);
   cEtaPlus1->SetBorderMode(0);
   cEtaPlus1->SetBorderSize(2);
   cEtaPlus1->SetTickx(1);
   cEtaPlus1->SetTicky(1);
   cEtaPlus1->SetLeftMargin(0.1594828);
   cEtaPlus1->SetRightMargin(0.05028735);
   cEtaPlus1->SetTopMargin(0.04872881);
   cEtaPlus1->SetBottomMargin(0.161017);
   cEtaPlus1->SetFrameBorderMode(0);
   cEtaPlus1->SetFrameBorderMode(0);
 
   hEtaPlusNominal2->GetXaxis()->SetTitle("|#eta|");
   hEtaPlusNominal2->GetXaxis()->SetLabelFont(42);
   hEtaPlusNominal2->GetXaxis()->SetLabelSize(0.05);
   hEtaPlusNominal2->GetXaxis()->SetTitleSize(0.05);
   hEtaPlusNominal2->GetXaxis()->SetTitleOffset(1.4);
   hEtaPlusNominal2->GetXaxis()->SetTitleFont(42);
   hEtaPlusNominal2->GetYaxis()->SetLabelFont(42);
   hEtaPlusNominal2->GetYaxis()->SetLabelSize(0.05);
   hEtaPlusNominal2->GetYaxis()->SetTitleSize(0.05);
   hEtaPlusNominal2->GetYaxis()->SetTitleOffset(1.4);
   hEtaPlusNominal2->GetYaxis()->SetTitleFont(42);
   hEtaPlusNominal2->SetMarkerSize(1.4);
   hEtaPlusNominal2->Draw("hist");
   hEtaPlusNominal2->Draw("pesame");

   hEtaPlusMptSigUp->SetMarkerStyle(23);
   hEtaPlusMptSigUp->SetLineStyle(3);
   hEtaPlusMptSigUp->SetLineWidth(0);
   hEtaPlusMptSigUp->SetLineColor(kRed);
   hEtaPlusMptSigUp->SetMarkerColor(kRed);
   hEtaPlusMptSigUp->Draw("hist same");
   hEtaPlusMptSigUp->Draw("psame");

   hEtaPlusMptSigDown->SetMarkerStyle(22);
   hEtaPlusMptSigDown->SetLineStyle(3);
   hEtaPlusMptSigDown->SetLineWidth(0);
   hEtaPlusMptSigDown->SetLineColor(kBlue);
   hEtaPlusMptSigDown->SetMarkerColor(kBlue);
   hEtaPlusMptSigDown->Draw("hist same");
   hEtaPlusMptSigDown->Draw("p same");

   leg = new TLegend(0.5732759,0.7351695,0.8433908,0.8855932,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("hEtaPlusNominal2","Nominal","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.4);
   entry=leg->AddEntry("hEtaPlusMptSigUp","#slash{p_{T}} 4GeV","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hEtaPlusMptSigDown","#slash{p_{T}} 2GeV","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(kBlue);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);

   leg->Draw();
   TLatex *   tex = new TLatex(0.2772989,0.3644068,"#mu^{+}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.1101695);
   tex->SetLineWidth(2);
   tex->Draw();

   ///mu-
   TCanvas *cEtaMinus1 = new TCanvas("cEtaMinus1", "cEtaMinus1",11,142,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cEtaMinus1->Range(-0.5307818,674.7044,2.797364,3793.449);
   cEtaMinus1->SetFillColor(0);
   cEtaMinus1->SetBorderMode(0);
   cEtaMinus1->SetBorderSize(2);
   cEtaMinus1->SetTickx(1);
   cEtaMinus1->SetTicky(1);
   cEtaMinus1->SetLeftMargin(0.1594828);
   cEtaMinus1->SetRightMargin(0.05028735);
   cEtaMinus1->SetTopMargin(0.04872881);
   cEtaMinus1->SetBottomMargin(0.161017);
   cEtaMinus1->SetFrameBorderMode(0);
   cEtaMinus1->SetFrameBorderMode(0);
 
   hEtaMinusNominal2->GetXaxis()->SetTitle("|#eta|");
   hEtaMinusNominal2->GetXaxis()->SetLabelFont(42);
   hEtaMinusNominal2->GetXaxis()->SetLabelSize(0.05);
   hEtaMinusNominal2->GetXaxis()->SetTitleSize(0.05);
   hEtaMinusNominal2->GetXaxis()->SetTitleOffset(1.4);
   hEtaMinusNominal2->GetXaxis()->SetTitleFont(42);
   hEtaMinusNominal2->GetYaxis()->SetLabelFont(42);
   hEtaMinusNominal2->GetYaxis()->SetLabelSize(0.05);
   hEtaMinusNominal2->GetYaxis()->SetTitleSize(0.05);
   hEtaMinusNominal2->GetYaxis()->SetTitleOffset(1.4);
   hEtaMinusNominal2->GetYaxis()->SetTitleFont(42);
   hEtaMinusNominal2->SetMarkerSize(1.4);
   hEtaMinusNominal2->Draw("hist");
   hEtaMinusNominal2->Draw("pesame");

   hEtaMinusMptSigUp->SetMarkerStyle(23);
   hEtaMinusMptSigUp->SetLineStyle(3);
   hEtaMinusMptSigUp->SetLineWidth(0);
   hEtaMinusMptSigUp->SetLineColor(kRed);
   hEtaMinusMptSigUp->SetMarkerColor(kRed);
   hEtaMinusMptSigUp->Draw("hist e same");
   hEtaMinusMptSigUp->Draw("p e same");

   hEtaMinusMptSigDown->SetMarkerStyle(22);
   hEtaMinusMptSigDown->SetLineStyle(3);
   hEtaMinusMptSigDown->SetLineWidth(0);
   hEtaMinusMptSigDown->SetLineColor(kBlue);
   hEtaMinusMptSigDown->SetMarkerColor(kBlue);
   hEtaMinusMptSigDown->Draw("hist e same");
   hEtaMinusMptSigDown->Draw("p e same");

   leg = new TLegend(0.3864943,0.7860169,0.6566092,0.9364407,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("hEtaMinusNominal2","Nominal","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.4);
   entry=leg->AddEntry("hEtaMinusMptSigUp","#slash{p_{T}} 4GeV","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hEtaMinusMptSigDown","#slash{p_{T}} 2GeV","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(kBlue);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);

   leg->Draw();
   tex = new TLatex(0.2772989,0.3644068,"#mu^{-}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.1101695);
   tex->SetLineWidth(2);
   tex->Draw();

   ///Isolation systematic in Eta
   ///asymmetry
   TCanvas *cAsymm2 = new TCanvas("cAsymm2", "cAsymm2",11,142,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cAsymm2->Range(-0.5307818,674.7044,2.797364,3793.449);
   cAsymm2->SetFillColor(0);
   cAsymm2->SetBorderMode(0);
   cAsymm2->SetBorderSize(2);
   cAsymm2->SetTickx(1);
   cAsymm2->SetTicky(1);
   cAsymm2->SetLeftMargin(0.1594828);
   cAsymm2->SetRightMargin(0.05028735);
   cAsymm2->SetTopMargin(0.04872881);
   cAsymm2->SetBottomMargin(0.161017);
   cAsymm2->SetFrameBorderMode(0);
   cAsymm2->SetFrameBorderMode(0);

   hAsymmNominal->GetXaxis()->SetTitle("|#eta|");
   hAsymmNominal->GetYaxis()->SetTitle("A_{#mu}");
   hAsymmNominal->GetXaxis()->SetLabelFont(42);
   hAsymmNominal->GetXaxis()->SetLabelSize(0.05);
   hAsymmNominal->GetXaxis()->SetTitleSize(0.05);
   hAsymmNominal->GetXaxis()->SetTitleOffset(1.4);
   hAsymmNominal->GetXaxis()->SetTitleFont(42);
   hAsymmNominal->GetYaxis()->SetLabelFont(42);
   hAsymmNominal->GetYaxis()->SetLabelSize(0.05);
   hAsymmNominal->GetYaxis()->SetTitleSize(0.05);
   hAsymmNominal->GetYaxis()->SetTitleOffset(1.4);
   hAsymmNominal->GetYaxis()->SetTitleFont(42);
   hAsymmNominal->SetMarkerSize(1.4);
   hAsymmNominal->GetYaxis()->SetRangeUser(-0.4,+0.3);
   hAsymmNominal->Draw("hist");
   hAsymmNominal->Draw("pesame");

   hAsymmIncCone->SetMarkerStyle(23);
   hAsymmIncCone->SetLineColor(kRed);
   hAsymmIncCone->SetMarkerColor(kRed);
   hAsymmIncCone->SetLineStyle(3);
   hAsymmIncCone->SetLineWidth(0);
   hAsymmIncCone->Draw("hist e same");
   hAsymmIncCone->Draw("p e same");
   hAsymmLoosenIso->SetMarkerStyle(22);
   hAsymmLoosenIso->SetLineColor(kBlue);
   hAsymmLoosenIso->SetMarkerColor(kBlue);
   hAsymmLoosenIso->SetLineStyle(3);
   hAsymmLoosenIso->SetLineWidth(0);
   hAsymmLoosenIso->Draw("hist e same");
   hAsymmLoosenIso->Draw("p e same");

/*
   hAsymmIsoAvg->SetMarkerStyle(24);
   hAsymmIsoAvg->SetLineStyle(3);
   hAsymmIsoAvg->SetLineWidth(0);
   hAsymmIsoAvg->SetMarkerColor(kBlack);
   hAsymmIsoAvg->Draw("hist e same");
*/

   leg = new TLegend(0.234,0.235,0.5,0.49,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("hAsymmNominal","Nominal","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hAsymmIncCone","Increase isolation cone","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(kRed);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hAsymmLoosenIso","Loosen isolation cut","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(kBlue);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);

/*   entry=leg->AddEntry("hAsymmIsoAvg","Average","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.2);
*/
   leg->Draw();


   TCanvas *cEtaPlus2 = new TCanvas("cEtaPlus2", "cEtaPlus2",11,142,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cEtaPlus2->Range(-0.5307818,674.7044,2.797364,3793.449);
   cEtaPlus2->SetFillColor(0);
   cEtaPlus2->SetBorderMode(0);
   cEtaPlus2->SetBorderSize(2);
   cEtaPlus2->SetTickx(1);
   cEtaPlus2->SetTicky(1);
   cEtaPlus2->SetLeftMargin(0.1594828);
   cEtaPlus2->SetRightMargin(0.05028735);
   cEtaPlus2->SetTopMargin(0.04872881);
   cEtaPlus2->SetBottomMargin(0.161017);
   cEtaPlus2->SetFrameBorderMode(0);
   cEtaPlus2->SetFrameBorderMode(0);

   hEtaPlusNominal->GetXaxis()->SetTitle("|#eta|");
   hEtaPlusNominal->GetXaxis()->SetLabelFont(42);
   hEtaPlusNominal->GetXaxis()->SetLabelSize(0.05);
   hEtaPlusNominal->GetXaxis()->SetTitleSize(0.05);
   hEtaPlusNominal->GetXaxis()->SetTitleOffset(1.4);
   hEtaPlusNominal->GetXaxis()->SetTitleFont(42);
   hEtaPlusNominal->GetYaxis()->SetLabelFont(42);
   hEtaPlusNominal->GetYaxis()->SetLabelSize(0.05);
   hEtaPlusNominal->GetYaxis()->SetTitleSize(0.05);
   hEtaPlusNominal->GetYaxis()->SetTitleOffset(1.4);
   hEtaPlusNominal->GetYaxis()->SetTitleFont(42);
   hEtaPlusNominal->SetMarkerSize(1.4);
   hEtaPlusNominal->Draw("hist");
   hEtaPlusNominal->Draw("pesame");

   hEtaPlusIncCone->SetMarkerStyle(23);
   hEtaPlusIncCone->SetLineColor(kRed);
   hEtaPlusIncCone->SetMarkerColor(kRed);
   hEtaPlusIncCone->SetLineStyle(3);
   hEtaPlusIncCone->SetLineWidth(0);
   hEtaPlusIncCone->Draw("hist e same");
   hEtaPlusIncCone->Draw("p e same");
   hEtaPlusLoosenIso->SetMarkerStyle(22);
   hEtaPlusLoosenIso->SetLineColor(kBlue);
   hEtaPlusLoosenIso->SetMarkerColor(kBlue);
   hEtaPlusLoosenIso->SetLineStyle(3);
   hEtaPlusLoosenIso->SetLineWidth(0);
   hEtaPlusLoosenIso->Draw("hist e same");
   hEtaPlusLoosenIso->Draw("p e same");
   leg = new TLegend(0.4698276,0.75,0.7399425,0.9004237,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("hEtaPlusNominal","Nominal","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hEtaPlusIncCone","Increase isolation cone","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(kRed);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hEtaPlusLoosenIso","Loosen isolation cut","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(kBlue);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);

   leg->Draw();
   tex = new TLatex(0.2772989,0.3644068,"#mu^{+}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.1101695);
   tex->SetLineWidth(2);
   tex->Draw();

   ///mu-
   TCanvas *cEtaMinus2 = new TCanvas("cEtaMinus2", "cEtaMinus2",11,142,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cEtaMinus2->Range(-0.5307818,674.7044,2.797364,3793.449);
   cEtaMinus2->SetFillColor(0);
   cEtaMinus2->SetBorderMode(0);
   cEtaMinus2->SetBorderSize(2);
   cEtaMinus2->SetTickx(1);
   cEtaMinus2->SetTicky(1);
   cEtaMinus2->SetLeftMargin(0.1594828);
   cEtaMinus2->SetRightMargin(0.05028735);
   cEtaMinus2->SetTopMargin(0.04872881);
   cEtaMinus2->SetBottomMargin(0.161017);
   cEtaMinus2->SetFrameBorderMode(0);
   cEtaMinus2->SetFrameBorderMode(0);

   hEtaMinusNominal->GetXaxis()->SetTitle("|#eta|");
   hEtaMinusNominal->GetXaxis()->SetLabelFont(42);
   hEtaMinusNominal->GetXaxis()->SetLabelSize(0.05);
   hEtaMinusNominal->GetXaxis()->SetTitleSize(0.05);
   hEtaMinusNominal->GetXaxis()->SetTitleOffset(1.4);
   hEtaMinusNominal->GetXaxis()->SetTitleFont(42);
   hEtaMinusNominal->GetYaxis()->SetLabelFont(42);
   hEtaMinusNominal->GetYaxis()->SetLabelSize(0.05);
   hEtaMinusNominal->GetYaxis()->SetTitleSize(0.05);
   hEtaMinusNominal->GetYaxis()->SetTitleOffset(1.4);
   hEtaMinusNominal->GetYaxis()->SetTitleFont(42);
   hEtaMinusNominal->SetMarkerSize(1.4);
   hEtaMinusNominal->Draw("hist");
   hEtaMinusNominal->Draw("pesame");

   hEtaMinusIncCone->SetMarkerStyle(23);
   hEtaMinusIncCone->SetLineColor(kRed);
   hEtaMinusIncCone->SetMarkerColor(kRed);
   hEtaMinusIncCone->SetLineStyle(3);
   hEtaMinusIncCone->SetLineWidth(0);
   hEtaMinusIncCone->Draw("hist e same");
   hEtaMinusIncCone->Draw("p e same");
   hEtaMinusLoosenIso->SetMarkerStyle(22);
   hEtaMinusLoosenIso->SetLineColor(kBlue);
   hEtaMinusLoosenIso->SetMarkerColor(kBlue);
   hEtaMinusLoosenIso->SetLineStyle(3);
   hEtaMinusLoosenIso->SetLineWidth(0);
   hEtaMinusLoosenIso->Draw("hist e same");
   hEtaMinusLoosenIso->Draw("p e same");
   leg = new TLegend(0.3864943,0.7860169,0.6566092,0.9364407,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("hEtaMinusNominal","Nominal","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hEtaMinusIncCone","Increase isolation cone","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(kRed);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hEtaMinusLoosenIso","Loosen isolation cut","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(kBlue);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);

   leg->Draw();
   tex = new TLatex(0.2772989,0.3644068,"#mu^{-}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.1101695);
   tex->SetLineWidth(2);
   tex->Draw();

   ///Z and QCD bkg systematics
   ///asymmetry
   TCanvas *cAsymm3 = new TCanvas("cAsymm3", "cAsymm3",11,142,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cAsymm3->Range(-0.5307818,674.7044,2.797364,3793.449);
   cAsymm3->SetFillColor(0);
   cAsymm3->SetBorderMode(0);
   cAsymm3->SetBorderSize(2);
   cAsymm3->SetTickx(1);
   cAsymm3->SetTicky(1);
   cAsymm3->SetLeftMargin(0.1594828);
   cAsymm3->SetRightMargin(0.05028735);
   cAsymm3->SetTopMargin(0.04872881);
   cAsymm3->SetBottomMargin(0.161017);
   cAsymm3->SetFrameBorderMode(0);
   cAsymm3->SetFrameBorderMode(0);

   hAsymmNominal2->GetXaxis()->SetTitle("|#eta|");
   hAsymmNominal2->GetYaxis()->SetTitle("A_{#mu}");
   hAsymmNominal2->GetXaxis()->SetLabelFont(42);
   hAsymmNominal2->GetXaxis()->SetLabelSize(0.05);
   hAsymmNominal2->GetXaxis()->SetTitleSize(0.05);
   hAsymmNominal2->GetXaxis()->SetTitleOffset(1.4);
   hAsymmNominal2->GetXaxis()->SetTitleFont(42);
   hAsymmNominal2->GetYaxis()->SetLabelFont(42);
   hAsymmNominal2->GetYaxis()->SetLabelSize(0.05);
   hAsymmNominal2->GetYaxis()->SetTitleSize(0.05);
   hAsymmNominal2->GetYaxis()->SetTitleOffset(1.4);
   hAsymmNominal2->GetYaxis()->SetTitleFont(42);
   hAsymmNominal2->SetMarkerSize(1.4);
   hAsymmNominal2->GetYaxis()->SetRangeUser(-0.4,+0.3);
   hAsymmNominal2->Draw("hist");
   hAsymmNominal2->Draw("pesame");

   hAsymmZBkg->SetMarkerStyle(23);
   hAsymmZBkg->SetLineColor(kRed);
   hAsymmZBkg->SetMarkerColor(kRed);
   hAsymmZBkg->SetLineStyle(3);
   hAsymmZBkg->SetLineWidth(0);
   hAsymmZBkg->Draw("hist e same");
   hAsymmZBkg->Draw("p e same");

   hAsymmTauBkg->SetMarkerStyle(21);
   hAsymmTauBkg->SetLineColor(kGreen);
   hAsymmTauBkg->SetMarkerColor(kGreen);
   hAsymmTauBkg->SetLineStyle(3);
   hAsymmTauBkg->SetLineWidth(0);
   hAsymmTauBkg->Draw("hist e same");
   hAsymmTauBkg->Draw("p e same");

   hAsymmRaaScaledQCD->SetMarkerStyle(22);
   hAsymmRaaScaledQCD->SetLineColor(kBlue);
   hAsymmRaaScaledQCD->SetMarkerColor(kBlue);
   hAsymmRaaScaledQCD->SetLineStyle(3);
   hAsymmRaaScaledQCD->SetLineWidth(0);
   hAsymmRaaScaledQCD->Draw("hist e same");
   hAsymmRaaScaledQCD->Draw("p e same");

   leg = new TLegend(0.234,0.235,0.5,0.49,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("hAsymmNominal2","Nominal2","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hAsymmZBkg","Z","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(kRed);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hAsymmTauBkg","#tau","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(kGreen);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.2);

   entry=leg->AddEntry("hAsymmRaaScaledQCD","QCD","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(kBlue);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   leg->Draw();



///////////////////////////
//Binary scaling///////////
///////////////////////////
   ///Mpt systematics
   TCanvas *cCent1 = new TCanvas("cCent1", "cCent1",656,152,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cCent1->Range(-83.50784,22.37975,438.4162,32.50633);
   cCent1->SetFillColor(0);
   cCent1->SetBorderMode(0);
   cCent1->SetBorderSize(2);
   cCent1->SetTickx(1);
   cCent1->SetTicky(1);
   cCent1->SetLeftMargin(0.16);
   cCent1->SetRightMargin(0.05);
   cCent1->SetTopMargin(0.05);
   cCent1->SetBottomMargin(0.16);
   cCent1->SetFrameBorderMode(0);
   cCent1->SetFrameBorderMode(0);

   hCentNominal2->SetLineWidth(2);
   hCentNominal2->SetMarkerStyle(20);
   hCentNominal2->SetMarkerSize(1.4);
   hCentNominal2->GetXaxis()->SetTitle("#LT N_{part} #GT");
   hCentNominal2->GetYaxis()->SetTitle("#frac{10^{9}}{#LT N_{coll} #GT}#frac{N^{W#rightarrow#mu,0.1<|#eta|<2.4}}{N_{events}}");
   hCentNominal2->SetMinimum(24);
   hCentNominal2->SetMaximum(32);
   hCentNominal2->GetXaxis()->SetLabelFont(42);
   hCentNominal2->GetXaxis()->SetLabelSize(0.05);
   hCentNominal2->GetXaxis()->SetTitleSize(0.05);
   hCentNominal2->GetXaxis()->SetTitleOffset(1.4);
   hCentNominal2->GetXaxis()->SetTitleFont(42);
   hCentNominal2->GetYaxis()->SetLabelFont(42);
   hCentNominal2->GetYaxis()->SetLabelSize(0.05);
   hCentNominal2->GetYaxis()->SetTitleSize(0.05);
   hCentNominal2->GetYaxis()->SetTitleOffset(1.4);
   hCentNominal2->GetYaxis()->SetTitleFont(42);
   hCentNominal2->GetYaxis()->SetRangeUser(24.0,32.0);
   hCentNominal2->Draw("hist");
   hCentNominal2->Draw("pesame");

   hCentMptSigUp->SetLineWidth(0);
   hCentMptSigUp->SetMarkerStyle(23);
   hCentMptSigUp->SetMarkerSize(1.2);
   hCentMptSigUp->SetLineStyle(3);
   hCentMptSigUp->SetLineColor(kRed);
   hCentMptSigUp->SetMarkerColor(kRed);
   hCentMptSigUp->Draw("histsame");
   hCentMptSigUp->Draw("pesame");

   hCentMptSigDown->SetLineWidth(0);
   hCentMptSigDown->SetMarkerStyle(22);
   hCentMptSigDown->SetMarkerSize(1.2);
   hCentMptSigDown->SetLineStyle(3);
   hCentMptSigDown->SetLineColor(kBlue);
   hCentMptSigDown->SetMarkerColor(kBlue);
   hCentMptSigDown->Draw("histsame");
   hCentMptSigDown->Draw("pesame");

   leg = new TLegend(0.42,0.24,0.69,0.39,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("hCentNominal2","Nominal","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.4);
   entry=leg->AddEntry("hCentMptSigUp","#slash{p_{T}} 4GeV","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hCentMptSigDown","#slash{p_{T}} 2GeV","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(kBlue);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);

   leg->Draw();

   ///Isolation systematics
   TCanvas *cCent2 = new TCanvas("cCent2", "cCent2",1983,107,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cCent2->Range(-73.2111,26.1266,443.9786,33.39859);
   cCent2->SetFillColor(0);
   cCent2->SetBorderMode(0);
   cCent2->SetBorderSize(2);
   cCent2->SetTickx(1);
   cCent2->SetTicky(1);
   cCent2->SetLeftMargin(0.16);
   cCent2->SetRightMargin(0.05);
   cCent2->SetTopMargin(0.05);
   cCent2->SetBottomMargin(0.16);
   cCent2->SetFrameBorderMode(0);
   cCent2->SetFrameBorderMode(0);

   hCentNominal->SetLineWidth(2);
   hCentNominal->SetMarkerStyle(20);
   hCentNominal->SetMarkerSize(1.4);
   hCentNominal->GetXaxis()->SetTitle("#LT N_{part} #GT");
   hCentNominal->GetYaxis()->SetTitle("#frac{10^{9}}{#LT N_{coll} #GT}#frac{N^{W#rightarrow#mu,0.1<|#eta|<2.4}}{N_{events}}");
   hCentNominal->GetXaxis()->SetLabelFont(42);
   hCentNominal->GetXaxis()->SetLabelSize(0.05);
   hCentNominal->GetXaxis()->SetTitleSize(0.05);
   hCentNominal->GetXaxis()->SetTitleOffset(1.4);
   hCentNominal->GetXaxis()->SetTitleFont(42);
   hCentNominal->GetYaxis()->SetLabelFont(42);
   hCentNominal->GetYaxis()->SetLabelSize(0.05);
   hCentNominal->GetYaxis()->SetTitleSize(0.05);
   hCentNominal->GetYaxis()->SetTitleOffset(1.4);
   hCentNominal->GetYaxis()->SetTitleFont(42);
   hCentNominal->GetYaxis()->SetRangeUser(25.5,34.0);
   hCentNominal->Draw("hist");
   hCentNominal->Draw("pesame");

   hCentIncCone->SetLineWidth(0);
   hCentIncCone->SetMarkerStyle(23);
   hCentIncCone->SetMarkerSize(1.2);
   hCentIncCone->SetLineStyle(3);
   hCentIncCone->SetLineColor(kRed);
   hCentIncCone->SetMarkerColor(kRed);
   hCentIncCone->Draw("histsame");
   hCentIncCone->Draw("pesame");
   hCentLoosenIso->SetLineWidth(0);
   hCentLoosenIso->SetMarkerStyle(22);
   hCentLoosenIso->SetMarkerSize(1.2);
   hCentLoosenIso->SetLineStyle(3);
   hCentLoosenIso->SetLineColor(kBlue);
   hCentLoosenIso->SetMarkerColor(kBlue);
   hCentLoosenIso->Draw("histsame");
   hCentLoosenIso->Draw("pesame");

   leg = new TLegend(0.37,0.22,0.64,0.37,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("hCentNominal","Nominal","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.4);
   entry=leg->AddEntry("hCentIncCone","Increase isolation cone","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hCentLoosenIso","Loosen isolation cut","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(kBlue);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   leg->Draw();

   ///Z and QCD background systematics
   TCanvas *cCent3 = new TCanvas("cCent3", "cCent3",1983,107,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cCent3->Range(-73.2111,26.1266,443.9786,33.39859);
   cCent3->SetFillColor(0);
   cCent3->SetBorderMode(0);
   cCent3->SetBorderSize(2);
   cCent3->SetTickx(1);
   cCent3->SetTicky(1);
   cCent3->SetLeftMargin(0.16);
   cCent3->SetRightMargin(0.05);
   cCent3->SetTopMargin(0.05);
   cCent3->SetBottomMargin(0.16);
   cCent3->SetFrameBorderMode(0);
   cCent3->SetFrameBorderMode(0);

   hCentNominal2->SetLineWidth(2);
   hCentNominal2->SetMarkerStyle(20);
   hCentNominal2->SetMarkerSize(1.4);
   hCentNominal2->GetXaxis()->SetTitle("#LT N_{part} #GT");
   hCentNominal2->GetYaxis()->SetTitle("#frac{10^{9}}{#LT N_{coll} #GT}#frac{N^{W#rightarrow#mu,0.1<|#eta|<2.4}}{N_{events}}");
   hCentNominal2->GetXaxis()->SetLabelFont(42);
   hCentNominal2->GetXaxis()->SetLabelSize(0.05);
   hCentNominal2->GetXaxis()->SetTitleSize(0.05);
   hCentNominal2->GetXaxis()->SetTitleOffset(1.4);
   hCentNominal2->GetXaxis()->SetTitleFont(42);
   hCentNominal2->GetYaxis()->SetLabelFont(42);
   hCentNominal2->GetYaxis()->SetLabelSize(0.05);
   hCentNominal2->GetYaxis()->SetTitleSize(0.05);
   hCentNominal2->GetYaxis()->SetTitleOffset(1.4);
   hCentNominal2->GetYaxis()->SetTitleFont(42);
   hCentNominal2->GetYaxis()->SetRangeUser(25.5,34.0);
   hCentNominal2->Draw("hist");
   hCentNominal2->Draw("pesame");

   hCentZBkg->SetLineWidth(0);
   hCentZBkg->SetMarkerStyle(23);
   hCentZBkg->SetMarkerSize(1.2);
   hCentZBkg->SetLineStyle(3);
   hCentZBkg->SetLineColor(kRed);
   hCentZBkg->SetMarkerColor(kRed);
   hCentZBkg->Draw("histsame");
   hCentZBkg->Draw("pesame");
   hCentTauBkg->SetLineWidth(0);
   hCentTauBkg->SetMarkerStyle(21);
   hCentTauBkg->SetMarkerSize(1.2);
   hCentTauBkg->SetLineStyle(3);
   hCentTauBkg->SetLineColor(kGreen);
   hCentTauBkg->SetMarkerColor(kGreen);
   hCentTauBkg->Draw("histsame");
   hCentTauBkg->Draw("pesame");

   hCentRaaScaledQCD->SetLineWidth(0);
   hCentRaaScaledQCD->SetMarkerStyle(22);
   hCentRaaScaledQCD->SetMarkerSize(1.2);
   hCentRaaScaledQCD->SetLineStyle(3);
   hCentRaaScaledQCD->SetLineColor(kBlue);
   hCentRaaScaledQCD->SetMarkerColor(kBlue);
   hCentRaaScaledQCD->Draw("histsame");
   hCentRaaScaledQCD->Draw("pesame");

   leg = new TLegend(0.37,0.22,0.64,0.37,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("hCentNominal2","Nominal","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.4);
   entry=leg->AddEntry("hCentZBkg","Z","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hCentTauBkg","#tau","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(kGreen);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.2);

   entry=leg->AddEntry("hCentRaaScaledQCD","QCD","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(kBlue);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1.2);
   leg->Draw();


   cEtaPlus1->Print("mptSystematicMuPlusEta_"+sDate+".pdf");
   cEtaMinus1->Print("mptSystematicMuMinusEta_"+sDate+".pdf");
   cEtaPlus2->Print("isoSystematicMuPlusEta_"+sDate+".pdf");
   cEtaMinus2->Print("isoSystematicMuMinusEta_"+sDate+".pdf");
   cCent1->Print("mptSystematicCent_"+sDate+".pdf");
   cCent2->Print("isoSystematicCent_"+sDate+".pdf");
   cCent3->Print("bkgSystematicCent_"+sDate+".pdf");

   cAsymm1->Print("mptSystematicMuAsymm_"+sDate+".pdf");
   cAsymm2->Print("isoSystematicMuAsymm_"+sDate+".pdf");
   cAsymm3->Print("bkgSystematicMuAsymm_"+sDate+".pdf");

   ///Save variations into root file
   outFile->cd();
   grEtaDiffIso->Write("grEtaDiffIso");
   grEtaDiffMpt->Write("grEtaDiffMpt");
   //grEtaDiffRaaScaledQCD->Write("grEtaDiffRaaScaledQCD");
   //grEtaDiffZBkg->Write("grEtaDiffZBkg");

   grCentDiffIso->Write("grCentDiffIso");
   grCentDiffMpt->Write("grCentDiffMpt");
   //grCentDiffRaaScaledQCD->Write("grCentDiffRaaScaledQCD");
   //grCentDiffZBkg->Write("grCentDiffZBkg");
   ssAsymm.close();
   ssCent.close();
}
