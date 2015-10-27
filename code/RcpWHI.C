#include "TString.h"
#include "TFile.h"
#include "TColor.h"
#include "TGraphAsymmErrors.h"
#include "TList.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TStyle.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLatex.h"
#include "TClonesArray.h"

#include "RooCurve.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooMinuit.h"
#include "RooBinning.h"
#include "RooHistPdf.h"

#include "RcpDep.C"
#include "Systematics.C"
#include "EfficiencyCorrection.C"
// #include "/home/rikard/atlasstyle-00-03-03/AtlasUtils.h"
// #ifndef __CINT__
// #include "/home/rikard/atlasstyle-00-03-03/AtlasStyle.C"
// #endif
// #include "/afs/cern.ch/user/s/sandstro/atlasstyle-00-03-03/AtlasUtils.h"
// #include "/afs/cern.ch/user/s/sandstro/atlasstyle-00-03-03/AtlasStyle.C"
#include "AtlasUtils.h"
#include "AtlasUtils.C"
#include "AtlasStyle.C"

#include <vector>
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
TString format(float value, int precision=3) {
  std::stringstream svalue;
  svalue  << std::setprecision(precision) << value;
  return svalue.str();
}

///////////////////////////////////////////
//write to outFile
//////////////////////////////////////////
void Write(TFile* outFile, TGraphAsymmErrors* gr, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing TGraph to root file..." << std::endl;
	  gr->Write(sFile);
	  std::cout << "Done." << std::endl;
      gDirectory = dir;
}

void Write(TFile* outFile, TF1* gr, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing TF1 to root file..." << std::endl;
	  gr->Write();
	  std::cout << "Done." << std::endl;
      gDirectory = dir;
}
///////////////////////////////////////////
//writeToSpreadsheet
//////////////////////////////////////////
void writeYieldsToSpreadsheet(std::ostream& outputFile, TString icent,TString nPart, TString nObs,TString nBkg, TString sigCounts, 
				TString errStat, TString errSystUncorr, TString errSystCorr){

	std::cout << "Writing yields to spreadsheet..." << std::endl;
	outputFile << icent << "," << nPart << "," << nObs << "," << nBkg << "," 
		<<  sigCounts << "," << errStat << "," << errSystUncorr << "," << errSystCorr << std::endl;
}

void writeRcpToSpreadsheet(std::ostream& outputFile, int icent, double rcp, double errStat, double errSyst, double errNcoll){
	std::cout << "Writing Rcp to spreadsheet..." << std::endl;
	outputFile << icent << "," << rcp << "," << errStat << "," << errSyst << "," << errNcoll << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// sqr
///////////////////////////////////////////////////////////////////////////////
double sqr(double a) {
  return a*a;
}

///////////////////////////////////////////////////////////////////////////////
// calcSystematic with MC errors
///////////////////////////////////////////////////////////////////////////////
void calcSystematicToyMC(TGraphAsymmErrors& grSyst, TGraphAsymmErrors* grBaseline, TGraphAsymmErrors* grToyMC, TList* systGraphs, double lowBound = -9999, double highBound = 9999) 
{
  for ( int iPoint=0; iPoint < grBaseline->GetN(); iPoint++ ) {
    
    double xBase,yBase;
    grBaseline->GetPoint(iPoint, xBase, yBase);
    grSyst.SetPoint(iPoint, xBase, yBase);
    double xErrorLo = grBaseline->GetErrorXlow(iPoint);
    double xErrorHi = grBaseline->GetErrorXhigh(iPoint);

    double newYErrorLo = (grBaseline->GetErrorYlow(iPoint))*(grBaseline->GetErrorYlow(iPoint));
    double newYErrorHi = (grBaseline->GetErrorYhigh(iPoint))*(grBaseline->GetErrorYhigh(iPoint));
    
    newYErrorLo += (grToyMC->GetErrorYlow(iPoint))*(grToyMC->GetErrorYlow(iPoint));
    newYErrorHi += (grToyMC->GetErrorYhigh(iPoint))*(grToyMC->GetErrorYhigh(iPoint));
    double unbiasedSF = 1./(systGraphs->GetSize()-1); //1/(n-1)sum(xi-xmean)^2
    for (int iGraph=0; iGraph<systGraphs->GetSize(); iGraph++) {

      TGraphAsymmErrors* grTemp = (TGraphAsymmErrors*) systGraphs->At(iGraph);
      
      double x,y;
      grTemp->GetPoint(iPoint, x, y);
            
      newYErrorLo += unbiasedSF*(y-yBase)*(y-yBase);
      newYErrorHi += unbiasedSF*(y-yBase)*(y-yBase);
    }
    newYErrorLo += 0.0001; //1% overall systematic
    newYErrorHi += 0.0001;
    newYErrorLo = sqrt(newYErrorLo);
    newYErrorHi = sqrt(newYErrorHi);
    if ( (yBase+newYErrorHi)>highBound ) {
      newYErrorHi = highBound - yBase;
    }
    if ( (yBase-newYErrorLo)<lowBound ) {
      newYErrorLo = yBase - lowBound;
    }
    
    grSyst.SetPointError(iPoint, xErrorLo, xErrorHi, newYErrorLo, newYErrorHi);
  }
}

int indexIJ(int i, int j, int nj){
  return i*nj+j; 
}

int indexIJK(int i, int j, int k, int nj, int nk){
  return (i*nj+j)*nk+k; 
}

void fillRcp(int icent,double npart,double npartH,double npartL,double yWSum, double eyWSum,double syWSumUncorrelated, double syWSumCorrelated, 
		TGraphAsymmErrors* graphWSum, TGraphAsymmErrors* graphWSumSyst, TGraphAsymmErrors* graphWSumSyst_correlated, double aw=1.0){

      std::cout << "Aw correction factor in centrality bin: " << icent << " = " << aw << std::endl;
      ///stat only
      graphWSum->SetPoint(icent, npart+13, yWSum/aw);
      //double awErr = 0.03;
      double awErr = 0.0;
      double errorTemp = yWSum/aw*TMath::Sqrt(TMath::Power(eyWSum/yWSum,2)+TMath::Power(0.0,2));
      graphWSum->SetPointError(icent, npartL+12.0, npartH+12.0, errorTemp, errorTemp);

      ///stat+uncorrelated systematics
      graphWSumSyst->SetPoint(icent, npart+13, yWSum/aw);
      errorTemp = yWSum/aw*TMath::Sqrt(TMath::Power(syWSumUncorrelated/yWSum,2)+TMath::Power(0.0,2));
      graphWSumSyst->SetPointError(icent, npartL+12.0, npartH+12.0, errorTemp,errorTemp);

      ///Correlated systematics
      graphWSumSyst_correlated->SetPoint(icent, npart, yWSum/aw);
      ///Include aw uncertainty in scaling uncertainty only
      errorTemp = yWSum/aw*TMath::Sqrt(TMath::Power(syWSumCorrelated/yWSum,2)+TMath::Power(awErr,2));
      graphWSumSyst_correlated->SetPointError(icent, npartL+5.0, npartH+5.0, errorTemp,errorTemp);

}//fill Rcp 

///Plot binary scaling at generator level
void fillGeneratorRcp(TString fileNameIn ,double xBins[], int bins, TH1F* hNpart,TH1F* hNpartPlus,TH1F* hNpartMinus){

   float etaGen[50],chargeGen[50];
   float nuGenPt[50], nuPhiGen[50];
   int mother[50],daughter[50];
   int ngen;
   float centrality;
   float muGenPt[50], muPhiGen[50];

   TChain* tree = new TChain("tree","tree");
   tree->Add(fileNameIn);

   std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;
  tree->SetBranchAddress("mc_mu_gen_pt", &muGenPt);
  tree->SetBranchAddress("mc_nu_gen_pt", &nuGenPt);
  tree->SetBranchAddress("mc_nu_gen_phi", &nuPhiGen);
  tree->SetBranchAddress("mc_mu_gen_mothertype", &mother);
  tree->SetBranchAddress("mc_mu_gen_type", &daughter);
  tree->SetBranchAddress("mc_mu_charge", &chargeGen);
  tree->SetBranchAddress("mc_mu_gen_eta", &etaGen);
  tree->SetBranchAddress("mc_mu_gen_phi", &muPhiGen);
  tree->SetBranchAddress("mc_mu_n", &ngen);
  tree->SetBranchAddress("centrality", &centrality);

   // --- Set branch status ---
  tree->SetBranchStatus("*",0) ;
  tree->SetBranchStatus("mc_mu_n", 1);
  tree->SetBranchStatus("centrality", 1);
  tree->SetBranchStatus("mc_nu_gen_pt", 1);
  tree->SetBranchStatus("mc_nu_gen_phi", 1);
  tree->SetBranchStatus("mc_mu_gen_pt", 1);
  tree->SetBranchStatus("mc_mu_gen_mothertype", 1);
  tree->SetBranchStatus("mc_mu_gen_type", 1);
  tree->SetBranchStatus("mc_mu_charge", 1);
  tree->SetBranchStatus("mc_mu_gen_eta", 1);
  tree->SetBranchStatus("mc_mu_gen_phi", 1);

  TH1F* hCent = new TH1F("hCent","hCent",bins,xBins);
  TH1F* hCentPlus = new TH1F("hCentPlus","hCentPlus",bins,xBins);
  TH1F* hCentMinus = new TH1F("hCentMinus","hCentMinus",bins,xBins);

  int signalPosAll = 0;
  int signalPosFid = 0;
  int signalNegAll = 0;
  int signalNegFid = 0;
   for(int iev = 0; iev < tree->GetEntries(); ++iev){
 
      tree->LoadTree(iev);
      tree->GetEntry(iev);

      for(int igen=0; igen<ngen; ++igen){
        if(fabs(mother[igen])==24&&fabs(daughter[igen])==13){
            if((chargeGen[igen])>0.) {++signalPosAll;/*std::cout << "[+]" << std::endl;*/ }
            else if((chargeGen[igen])<0.) {++signalNegAll;/*std::cout << "[-]" << std::endl;*/ }
        }
        float dPhi = nuPhiGen[igen]-muPhiGen[igen];
        if(dPhi<-1*TMath::Pi()) dPhi += TMath::TwoPi(); if(dPhi>TMath::Pi()) dPhi -= TMath::TwoPi(); //fold btwn [-pi,pi]
        float mtTruth = TMath::Sqrt(2.0*muGenPt[igen]*nuGenPt[igen]*(1.0-TMath::Cos(dPhi)));
        if(
            fabs(mother[igen])==24
            &&fabs(daughter[igen])==13
            &&centrality>0.0&&centrality<0.8
            &&fabs(etaGen[igen])>0.1&&fabs(etaGen[igen])<2.4
            &&nuGenPt[igen]>25.0
            &&muGenPt[igen]>25.0
            &&mtTruth>40.0
            &&fabs(chargeGen[igen])>0.
        ){ 
        ///fill centrality distribution for inclusive set
        hCent->Fill(centrality);
        ///fill mu+ eta
        if(chargeGen[igen]>0.){
            hCentPlus->Fill(centrality);
            ++signalPosFid;
        }
        ///fill mu- eta
        else if(chargeGen[igen]<0.)
         {
            hCentMinus->Fill(centrality);
            ++signalNegFid;
         }
      }
      }//igen
    }//iev

   std::cout << std::endl;
   std::cout << "===========>PYTHIA<=========" << std::endl;
   std::cout << "Number of mu^{+} over all space: " << signalPosAll << std::endl;
   std::cout << "Number of mu^{+} in fiducial region: " << signalPosFid << std::endl;
   std::cout << "AW+: " << (float)signalPosFid/signalPosAll << std::endl;
   std::cout << "Number of mu^{-} over all space: " << signalNegAll << std::endl;
   std::cout << "Number of mu^{-} in fiducial region: " << signalNegFid << std::endl;
   std::cout << "AW-: " << (float)signalNegFid/signalNegAll << std::endl;
   std::cout << "===========>PYTHIA<=========" << std::endl;
   std::cout << std::endl;
    ///Now transform from centrlity-->Npart
    TGraphErrors* grNcoll = getNColl();
    double* ncoll = grNcoll->GetY();
    double nEventsSampled = (double)tree->GetEntries(); 
    std::cout << "Number of events in MC file: "<< nEventsSampled << std::endl;
    ///last bin is most peripheral 
    for(int ibin=hCent->GetNbinsX(); ibin>0; --ibin){
       double invBinW = 1.0/hCent->GetBinWidth(ibin);
       std::cout << "binwidth: "<< ibin << " = "  << hCent->GetBinWidth(ibin) << std::endl;
       double scale = 1.0/ncoll[ibin-1];
       std::cout << "Ncoll in bin: "<<ibin << " = "  << ncoll[ibin-1] << std::endl;
       double nWs = hCent->GetBinContent(ibin);
       std::cout << "MC Ws in bin: "<<ibin << " = "  << nWs << std::endl;
       //double scaledWYield = scale; scaledWYield/=nEventsSampled; scaledWYield*=invBinW;  
       double scaledWYield = invBinW;
       scaledWYield*=nWs; 
       hNpart->SetBinContent(ibin,scaledWYield);
       double scaledWYieldErr = TMath::Sqrt(nWs)*invBinW;
       hNpart->SetBinError(ibin,scaledWYieldErr);

       double nWPlus = hCentPlus->GetBinContent(ibin); 
       double scaledWPlusYield = 1.0e9*scale; scaledWPlusYield/=nEventsSampled; scaledWPlusYield*=invBinW;  
       scaledWPlusYield*=nWPlus; 
       hNpartPlus->SetBinContent(ibin,scaledWPlusYield);
       scaledWYieldErr = TMath::Sqrt(nWPlus)*invBinW;
       hNpartPlus->SetBinError(ibin,scaledWYieldErr);

       double nWMinus = hCentMinus->GetBinContent(ibin); 
       double scaledWMinusYield = 1.0e9*scale; scaledWMinusYield/=nEventsSampled; scaledWMinusYield*=invBinW;  
       scaledWMinusYield*=nWMinus; 
       hNpartMinus->SetBinContent(ibin,scaledWMinusYield);
       scaledWYieldErr = TMath::Sqrt(nWMinus)*invBinW;
       hNpartMinus->SetBinError(ibin,scaledWYieldErr);

    }//ibin
}

void getMcCrossSections(double& WCrossSectionPythiaLO,double& WCrossSectionMinusPythiaLO,double& WCrossSectionPlusPythiaLO,
                            double& wtdWCrossSectionMSTWnlo,double& wtdWCrossSectionMSTWnloMinus,double& wtdWCrossSectionMSTWnloPlus){
    ///PYTHIA LO*

    ///fiducial losses from PY6
    ///wtd average of all nucleon combos
/*    double awPlusPy = 0.544; 
    double awMinusPy = 0.536; 
*/
    double awPlusPy = 0.519; 
    double awMinusPy = 0.516; 

    WCrossSectionPythiaLO = 2.8205e-9/64.0e-3*1.0e9;
    WCrossSectionMinusPythiaLO = (1.0/2.01)*WCrossSectionPythiaLO;
    WCrossSectionMinusPythiaLO*=awMinusPy;
    WCrossSectionPlusPythiaLO = (1.0-1.0/2.01)*WCrossSectionPythiaLO;
    WCrossSectionPlusPythiaLO*=awPlusPy;
    WCrossSectionPythiaLO = WCrossSectionMinusPythiaLO+WCrossSectionPlusPythiaLO;


    ///fiducial losses from POWHEG+PY8
    double awPlusPow_nn = 0.548; 
    double awPlusPow_pp = 0.541; 
    double awPlusPow_nppn = 0.543; 
    double awMinusPow_nn = 0.460; 
    double awMinusPow_pp = 0.502; 
    double awMinusPow_nppn = 0.476; 

    ///MSTW2008 NLO 68% C.L. (POWHEG)
    ///pp
    double WCrossSection_ppPlusMSTWnlo = 2.175e-9/64.0e-3*1.0e9;
    double WCrossSection_ppMinusMSTWnlo = 1.330e-9/64.0e-3*1.0e9;
    double WCrossSection_ppMSTWnlo = WCrossSection_ppPlusMSTWnlo+WCrossSection_ppMinusMSTWnlo; 
    double wt_pp = 0.155;
    ///np
    double WCrossSection_nppnPlusMSTWnlo = 1.743e-9/64.0e-3*1.0e9;
    double WCrossSection_nppnMinusMSTWnlo = 1.744e-9/64.0e-3*1.0e9;
    double WCrossSection_nppnMSTWnlo = WCrossSection_nppnPlusMSTWnlo+WCrossSection_nppnMinusMSTWnlo; 
    double wt_nppn = 0.478;
    ///nn
    double WCrossSection_nnPlusMSTWnlo = 1.327e-9/64.0e-3*1.0e9;
    double WCrossSection_nnMinusMSTWnlo = 2.175e-9/64.0e-3*1.0e9;
    double WCrossSection_nnMSTWnlo = WCrossSection_nnPlusMSTWnlo+WCrossSection_nnMinusMSTWnlo; 
    double wt_nn = 0.367;

    double wtdAwPlusPow = wt_pp*awPlusPow_pp+wt_nppn*awPlusPow_nppn+wt_nn*awPlusPow_nn ;
    double wtdAwMinusPow = wt_pp*awMinusPow_pp+wt_nppn*awMinusPow_nppn+wt_nn*awMinusPow_nn ;

    ///weighted average from each nucleon combination
    wtdWCrossSectionMSTWnloPlus  = (wt_pp*WCrossSection_ppPlusMSTWnlo
                                            +wt_nppn*WCrossSection_nppnPlusMSTWnlo
                                            +wt_nn*WCrossSection_nnPlusMSTWnlo)/1.0;
    std::cout << "W^+ --> mu^+ + nu (b4 apply Aw) : " << wtdWCrossSectionMSTWnloPlus << std::endl; 
    wtdWCrossSectionMSTWnloPlus*=wtdAwPlusPow;
    std::cout << "AW^+: " << wtdAwPlusPow << std::endl;
    std::cout << "W^+ --> mu^+ + nu (after apply Aw) : " << wtdWCrossSectionMSTWnloPlus << std::endl; 
    wtdWCrossSectionMSTWnloMinus = (wt_pp*WCrossSection_ppMinusMSTWnlo
                                            +wt_nppn*WCrossSection_nppnMinusMSTWnlo
                                            +wt_nn*WCrossSection_nnMinusMSTWnlo)/1.0;
    std::cout << "W^- --> mu^- + nu~ (b4 apply Aw) : " << wtdWCrossSectionMSTWnloMinus << std::endl; 
    wtdWCrossSectionMSTWnloMinus*=wtdAwMinusPow;
    std::cout << "AW^-: " << wtdAwMinusPow << std::endl;
    std::cout << "W^- --> mu^- + nu~ (b4 apply Aw) : " << wtdWCrossSectionMSTWnloMinus << std::endl; 
    /*double wtdWCrossSectionMSTWnlo      = (wt_pp*WCrossSection_ppMSTWnlo
                                            +wt_nppn*WCrossSection_nppnMSTWnlo
                                            +wt_nn*WCrossSection_nnMSTWnlo)/1.0;
                                            */
    wtdWCrossSectionMSTWnlo = wtdWCrossSectionMSTWnloPlus+wtdWCrossSectionMSTWnloMinus;
    std::cout << "W^+- --> mu^+- + nu: " << wtdWCrossSectionMSTWnlo << std::endl; 

}

void plotChargeRatio(TGraphAsymmErrors* graphWRatio,TGraphAsymmErrors* graphWRatioSyst,TGraphAsymmErrors* graphWRatioSystCorrelated,
                        TGraph* grRatioPythia,TGraph* grRatioPwg,TGraph* grRatioQuark ){

      TCanvas* cWPlusMinus = new TCanvas("cWPlusMinus","cWPlusMinus",600,600);
      TH1F* hDummy = new TH1F("hDummy","hDummy",100,0.0,420.0);

      double WCrossSectionPythiaLO,WCrossSectionMinusPythiaLO,WCrossSectionPlusPythiaLO;
      double wtdWCrossSectionMSTWnlo,wtdWCrossSectionMSTWnloMinus,wtdWCrossSectionMSTWnloPlus;
      getMcCrossSections(WCrossSectionPythiaLO,WCrossSectionMinusPythiaLO,WCrossSectionPlusPythiaLO,
                        wtdWCrossSectionMSTWnlo,wtdWCrossSectionMSTWnloMinus,wtdWCrossSectionMSTWnloPlus);


      // Get overall scaling uncertainty (including cancellations in the ratio)
      double arrAsymm[] = {0.120,0.035,0.130,0.091,0.100,0.004,-0.013,-0.110,-0.300};
      double arrAsymmVariation1[] = {0.0154,0.0098,0.011,0.0087,0.0086,0.0102,0.00976,0.00525,0.00473};
      double arrAsymmVariation2[] = {0.0134365,0.0170867,0.00728431,0.0135928,0.0145589,0.00857352,0.0288549,0.00740872,0.0181329};
      double wt[] = {0.138,0.133,0.103,0.125,0.117,0.107,0.103,0.0812,0.0919};
      int nEtaBins = sizeof(arrAsymm)/sizeof(double);
      double ratioVariation[nEtaBins] ;
      double wtdAvgRatioVariation = 0.0;
      for(int ieta=0; ieta<nEtaBins; ++ieta){
        
        double ratioVariation1 = (1.0/(1.0-arrAsymm[ieta])+1.0/(1+arrAsymm[ieta]))*arrAsymmVariation1[ieta];
        double ratioVariation2 = (1.0/(1.0-arrAsymm[ieta])+1.0/(1+arrAsymm[ieta]))*arrAsymmVariation2[ieta];
        std::cout << "Ratio variation1: " << ratioVariation1 << " Ratio variation2: " << ratioVariation2 << std::endl;
        ratioVariation[ieta] = TMath::Power(ratioVariation1,2);
        ratioVariation[ieta]+=TMath::Power(ratioVariation2,2);
        ratioVariation[ieta] =  TMath::Sqrt(ratioVariation[ieta]);
        double wtdAverage = wt[ieta]*ratioVariation[ieta];
        wtdAvgRatioVariation+=wtdAverage; 
        
      }//ieta

      std::cout << "Total scaling uncertainty in ratio: " << wtdAvgRatioVariation << std::endl;
      TGraphErrors* grRatioScalingUncert = new TGraphErrors(6);
      grRatioScalingUncert->SetFillColor(kGray+1);
      grRatioScalingUncert->SetFillStyle(3244);

      TF1* funcRatioPythia = new TF1("funcRatioPythia","[0]",0.0,420.0);
      funcRatioPythia->SetLineStyle(2);
      funcRatioPythia->SetLineColor(kRed);
      funcRatioPythia->SetParameter(0,WCrossSectionPlusPythiaLO/WCrossSectionMinusPythiaLO);
      funcRatioPythia->FixParameter(0,WCrossSectionPlusPythiaLO/WCrossSectionMinusPythiaLO);

      grRatioPythia->SetPoint(0,10.0,WCrossSectionPlusPythiaLO/WCrossSectionMinusPythiaLO);
      grRatioPythia->SetMarkerStyle(27);
      grRatioPythia->SetMarkerSize(2.5);
      grRatioPythia->SetMarkerColor(kRed);

      TF1* funcRatioPwg = new TF1("funcRatioPwg","[0]",0.0,420.0);
      funcRatioPwg->SetLineStyle(9);
      funcRatioPwg->SetLineColor(kBlue);
      funcRatioPwg->SetParameter(0,wtdWCrossSectionMSTWnloPlus/wtdWCrossSectionMSTWnloMinus);
      funcRatioPwg->FixParameter(0,wtdWCrossSectionMSTWnloPlus/wtdWCrossSectionMSTWnloMinus);

      grRatioPwg->SetPoint(0,27.0,wtdWCrossSectionMSTWnloPlus/wtdWCrossSectionMSTWnloMinus);
      grRatioPwg->SetMarkerStyle(28);
      grRatioPwg->SetMarkerSize(2.3);
      grRatioPwg->SetMarkerColor(kBlue);

      TF1* funcRatioQuark = new TF1("funcRatioQuark","[0]",0.0,420.0);
      funcRatioQuark->SetLineStyle(6);
      funcRatioQuark->SetLineColor(kViolet);
      funcRatioQuark->SetParameter(0,0.87);
      funcRatioQuark->FixParameter(0,0.87);

      grRatioQuark->SetPoint(0,15.0,0.87);

      grRatioQuark->SetMarkerStyle(25);
      grRatioQuark->SetMarkerSize(1.7);
      grRatioQuark->SetMarkerColor(kViolet);

      TGraphAsymmErrors* graphWRatioDummy = (TGraphAsymmErrors*)graphWRatio->Clone("graphWRatioDummy");
      graphWRatioDummy->SetMarkerSize(2.2);
      graphWRatioDummy->SetMarkerStyle(33);
      graphWRatioDummy->SetMarkerColor(kWhite);
      graphWRatio->SetMarkerSize(1.9);
      graphWRatio->SetMarkerStyle(33);

      TGraphAsymmErrors* graphWRatioSystc = (TGraphAsymmErrors*)graphWRatioSyst->Clone("graphWRatioSystc");
      graphWRatioSystc->SetFillColor(kGray+1);
      graphWRatioSystc->SetLineColor(kGray+1);
      //hDummy->GetYaxis()->SetTitle("N^{W^{+}#rightarrow#mu^{+},fid}/N^{W^{-}#rightarrow#mu^{-},fid}");
      hDummy->GetYaxis()->SetTitle("Fiducial Charge Ratio N^{W^{+}}/N^{W^{-}}");
      hDummy->GetYaxis()->SetTitleSize(0.04);
      hDummy->GetYaxis()->SetLabelSize(0.03);
      hDummy->GetXaxis()->SetLabelSize(0.03);
      hDummy->GetXaxis()->SetTitle("#LT N_{part} #GT");

      TGraphAsymmErrors* graphWRatioSystCorrelatedc = (TGraphAsymmErrors*)graphWRatioSystCorrelated->Clone("graphWRatioSystCorrelatedc");
      graphWRatioSystCorrelatedc->SetFillColor(kGray+1);
      graphWRatioSystCorrelatedc->SetFillStyle(3244);
      for(int i = 0; i<graphWRatioSystCorrelatedc->GetN(); ++i){

        double xtemp = graphWRatioSystCorrelatedc->GetX()[i];
        double xtempErr = graphWRatioSystCorrelatedc->GetEXhigh()[i];
        double ytemp = graphWRatioSystCorrelatedc->GetY()[i];
        grRatioScalingUncert->SetPoint(i,xtemp,ytemp);
        grRatioScalingUncert->SetPointError(i,xtempErr,wtdAvgRatioVariation);
      }//i

///uncomment for systematic curve
      hDummy->GetXaxis()->SetRangeUser(0.0,420.0);
      hDummy->GetYaxis()->SetRangeUser(0.0,1.3);
      hDummy->Draw("0");
      cWPlusMinus->Update();

      funcRatioPythia->Draw("lsame");
      funcRatioPwg->Draw("lsame");
      funcRatioQuark->Draw("lsame");
      //graphWRatioSystCorrelatedc->Draw("e2same");
      grRatioScalingUncert->Draw("e2same");
      graphWRatioSystc->Draw("e2same");
      graphWRatioDummy->Draw("pesame");
      graphWRatio->Draw("pesame");
      //grRatioPythia->Draw("psame");
      //grRatioQuark->Draw("psame");
      //grRatioPwg->Draw("psame");
      cWPlusMinus->Update();

   TLegend *leg = new TLegend(0.2315436,0.3216783,0.533557,0.5244755,NULL,"brNDC");
   leg->SetTextFont(gStyle->GetTextFont());
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("NULL","Data 2011","ple");

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#666666");
   entry->SetFillColor(ci);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(1.9);
   entry=leg->AddEntry("NULL","PYTHIA LO*","l");
   entry->SetLineColor(kRed);
   entry->SetLineStyle(2);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(27);
   entry->SetMarkerSize(2);
   entry=leg->AddEntry("NULL","POWHEG NLO","l");
   entry->SetLineColor(kBlue);
   entry->SetLineStyle(9);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(28);
   entry->SetMarkerSize(2.3);
   entry=leg->AddEntry("NULL","u/d quark counts","l");
   entry->SetLineColor(kViolet);
   entry->SetLineStyle(6);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#cc00ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(25);
   entry->SetMarkerSize(1.7);
   leg->Draw();
   TLatex *   tex = new TLatex(0.1828859,0.2255245,"#int Ldt #approx 0.14 nb^{-1}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.5805369,0.2237762,"Pb+Pb #sqrt{s_{NN}} = 2.76 TeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.204698,0.548951,"ATLAS");
tex->SetNDC();
   tex->SetTextFont(72);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.3758389,0.548951,"Preliminary");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
      cWPlusMinus->Print("HI_singleMuonWRatio.pdf"/*+fileNameWFit.ReplaceAll("#","")+".png"*/);
      cWPlusMinus->Print("HI_singleMuonWRatio.root"/*fileNameWFit.ReplaceAll("#","")+".eps"*/);
      //cWPlusMinus->Print("HI_singleMuonWRatio.png"/*fileNameWFit.ReplaceAll("#","")+".pdf"*/);



}

void plotRcp(TCanvas* cRcp, TGraphAsymmErrors* graphWSumc,TGraphAsymmErrors* graphWSumSystc,TGraphAsymmErrors* graphWSumSyst_correlatedc, TString sCharge,
                    TGraphAsymmErrors* graphWPlusSumc,TGraphAsymmErrors* graphWPlusSumSystc,TGraphAsymmErrors* graphWPlusSumSyst_correlatedc,
                    TGraphAsymmErrors* graphWMinusSumc,TGraphAsymmErrors* graphWMinusSumSystc,TGraphAsymmErrors* graphWMinusSumSyst_correlatedc,
                    TH1F* hGenAllNuc, TH1F* hGen_pp, TH1F* hGen_nppn, TH1F* hGen_nn
            ){

//    TCanvas* cRcp = new TCanvas("cRcp","cRcp",600,600);

    TPad *pad1 = new TPad("pad1","",0,0,1,1);
    TPad *pad2 = new TPad("pad2","",0,0,1,1);
    ///make pad1 transparent
    pad1->SetFillStyle(0);

    double arrNCollErr[] = {0.5, 0.9, 1.3, 1.6, 2.6, 5.1, 10.1};
    double arrNcoll[] = {45.93,157.83,239.52,281.88,330.26,382.16};
    int nElements = sizeof(arrNcoll)/sizeof(double) ;

    char* arrCentLabel[] = {"80","40","20","15","10","5","0"};
    double arrCentBin[6];
    arrCentBin[0] = 0.0;

    ///MCFM HERA01 NLO
    double WCrossSection_ppMCFMnlo = 3.351e-9/64.0e-3*1.0e9;
    double WCrossSection_ppPlusMCFMnlo = 2.067e-9/64.0e-3*1.0e9;
    double WCrossSection_ppMinusMCFMnlo = 1.285e-9/64.0e-3*1.0e9;


/*    double AwInclusive = 0.540207;
    double AwPlus = 0.544;
    double AwMinus = 0.536;
    ///accout for fiducial region in cross-section
    WCrossSectionPythiaLO*=AwInclusive;
    WCrossSection_ppMCFMnlo*=AwInclusive;
    WCrossSection_ppPlusMCFMnlo*=AwPlus;
    WCrossSection_ppMinusMCFMnlo*=AwMinus;
    */
    
    double WCrossSectionPythiaLO,WCrossSectionMinusPythiaLO,WCrossSectionPlusPythiaLO;
    double wtdWCrossSectionMSTWnlo,wtdWCrossSectionMSTWnloMinus,wtdWCrossSectionMSTWnloPlus;
    getMcCrossSections(WCrossSectionPythiaLO,WCrossSectionMinusPythiaLO,WCrossSectionPlusPythiaLO,
                        wtdWCrossSectionMSTWnlo,wtdWCrossSectionMSTWnloMinus,wtdWCrossSectionMSTWnloPlus);

    TGraphErrors *grPythia = new TGraphErrors(1);
    TGraph *grPythiaPlus = new TGraphErrors(1);
    TGraph *grPythiaMinus = new TGraphErrors(1);
    TGraphErrors *gr_ppMCFMnlo = new TGraphErrors(1);
    TGraphErrors *gr_ppPlusMCFMnlo = new TGraphErrors(1);
    TGraphErrors *gr_ppMinusMCFMnlo = new TGraphErrors(1);
    TGraph* grMSTWnlo = new TGraphErrors(1);
    TGraph* grMSTWnloPlus = new TGraphErrors(1);
    TGraph* grMSTWnloMinus = new TGraphErrors(1);

    /*for(int i=0; i<6; ++i) {
        arrCentBin[i+1] = arrCentBin[i]+2.*(arrNcoll[i] - arrCentBin[i]);
        //arrCentBin[i] = arrNcoll[i];
        std::cout << "Bins: " << arrCentBin[i] << "-" << arrCentBin[i+1] << std::endl;
        ///Shaded band from PYTHIA@LO* prediction
    }*/
    ///array for plotting bands
    double arrTemp[] = {0.0,45.93,157.83,239.52,281.88,330.26,420.0};
    //for(int i=0;i<7;++i){

    // Pythia LO* prediction
    grPythia->SetPoint(0,15.0,WCrossSectionPythiaLO);
    grPythiaPlus->SetPoint(0,10.0,WCrossSectionPlusPythiaLO);
    grPythiaMinus->SetPoint(0,20.0,WCrossSectionMinusPythiaLO);
    TF1* funcPythia = new TF1("funcPythia","[0]",0.0,420.0);
    TF1* funcPythiaPlus = new TF1("funcPythiaPlus","[0]",0.0,420.0);
    TF1* funcPythiaMinus = new TF1("funcPythiaMinus","[0]",0.0,420.0);
    funcPythia->SetLineStyle(2);
    funcPythia->SetLineColor(kGray+2);
    funcPythia->SetParameter(0,WCrossSectionPythiaLO);
    funcPythia->FixParameter(0,WCrossSectionPythiaLO);
    funcPythiaPlus->SetLineStyle(4);
    funcPythiaPlus->SetLineColor(kRed);
    funcPythiaPlus->SetParameter(0,WCrossSectionPlusPythiaLO);
    funcPythiaPlus->FixParameter(0,WCrossSectionPlusPythiaLO);
    funcPythiaMinus->SetLineStyle(5);
    funcPythiaMinus->SetLineColor(kBlue);
    funcPythiaMinus->SetParameter(0,WCrossSectionMinusPythiaLO);
    funcPythiaMinus->FixParameter(0,WCrossSectionMinusPythiaLO);


    gr_ppMCFMnlo->SetPoint(0,15.0,WCrossSection_ppMCFMnlo);
    gr_ppPlusMCFMnlo->SetPoint(0,15.0,WCrossSection_ppPlusMCFMnlo);
    gr_ppMinusMCFMnlo->SetPoint(0,15.0,WCrossSection_ppMinusMCFMnlo);
    grMSTWnlo->SetPoint(0,15.0,wtdWCrossSectionMSTWnlo);
    grMSTWnloPlus->SetPoint(0,10.0,wtdWCrossSectionMSTWnloPlus);
    grMSTWnloMinus->SetPoint(0,20.0,wtdWCrossSectionMSTWnloMinus);
    TF1* funcMSTWnlo = new TF1("funcMSTWnlo","[0]",0.0,420.0);
    TF1* funcMSTWnloPlus = new TF1("funcMSTWnloPlus","[0]",0.0,420.0);
    TF1* funcMSTWnloMinus = new TF1("funcMSTWnloMinus","[0]",0.0,420.0);
    funcMSTWnlo->SetLineStyle(9);
    funcMSTWnlo->SetLineColor(kGray+2);
    funcMSTWnlo->SetParameter(0,wtdWCrossSectionMSTWnlo);
    funcMSTWnlo->FixParameter(0,wtdWCrossSectionMSTWnlo);
    funcMSTWnloPlus->SetLineStyle(7);
    funcMSTWnloPlus->SetLineColor(kRed);
    funcMSTWnloPlus->SetParameter(0,wtdWCrossSectionMSTWnloPlus);
    funcMSTWnloPlus->FixParameter(0,wtdWCrossSectionMSTWnloPlus);
    funcMSTWnloMinus->SetLineStyle(8);
    funcMSTWnloMinus->SetLineColor(kBlue);
    funcMSTWnloMinus->SetParameter(0,wtdWCrossSectionMSTWnloMinus);
    funcMSTWnloMinus->FixParameter(0,wtdWCrossSectionMSTWnloMinus);



        grPythia->SetPointError(0,0.0,0.0);
        gr_ppMCFMnlo->SetPointError(0,0.0,0.0);
        gr_ppPlusMCFMnlo->SetPointError(0,0.0,0.0);
        gr_ppMinusMCFMnlo->SetPointError(0,0.0,0.0);

    //}

    TH1F* hDummy3 = new TH1F("hDummy3","hDummy3",6,arrCentBin);

    TAxis *ax = hDummy3->GetXaxis();
    TAxis *ay = hDummy3->GetYaxis();
    /*for(int i=0; i<6; ++i) {
        std::cout << "Bin: " << hDummy3->GetBinCenter(i+1) << std::endl;
        std::cout <<  "Centrality " << arrCentLabel[i] << std::endl;
        ax->SetBinLabel(i+1,arrCentLabel[i]);
    }*/

    TH1F* hDummy = new TH1F("hDummy","hDummy",100,0.0,420.0);
    //ax->SetLabelSize(0.);
    ay->SetLabelSize(0.);
    pad2->Draw();
    pad2->cd();
    //hDummy3->Draw("X+");
    pad2->Update();

    pad1->Draw();
    pad1->cd();

    grPythia->SetMarkerColor(kGray+1);
    grPythia->SetMarkerStyle(27);
    grPythia->SetMarkerSize(2.5);
//    grPythia->SetFillStyle(3001);

    gr_ppMCFMnlo->SetMarkerColor(kGray+1);
    gr_ppMCFMnlo->SetMarkerStyle(28);
    gr_ppMCFMnlo->SetMarkerSize(2.3);
    grMSTWnlo->SetMarkerColor(kGray+1);
    grMSTWnlo->SetMarkerStyle(28);
    grMSTWnlo->SetMarkerSize(2.3);
/*    gr_ppMCFMnlo->SetLineWidth(4);
    gr_ppMCFMnlo->SetFillColor(kGray+1);
    gr_ppMCFMnlo->SetFillStyle(3001);
*/
    grPythiaPlus->SetMarkerColor(kRed);
    grPythiaPlus->SetMarkerStyle(27);
    grPythiaPlus->SetMarkerSize(2.5);
    gr_ppPlusMCFMnlo->SetMarkerColor(kRed);
    gr_ppPlusMCFMnlo->SetMarkerStyle(28);
    gr_ppPlusMCFMnlo->SetMarkerSize(2.3);
    grMSTWnloPlus->SetMarkerColor(kRed);
    grMSTWnloPlus->SetMarkerStyle(28);
    grMSTWnloPlus->SetMarkerSize(2.3);
/*    gr_ppPlusMCFMnlo->SetLineWidth(4);
    gr_ppPlusMCFMnlo->SetFillColor(kRed);
    gr_ppPlusMCFMnlo->SetFillStyle(3001);
*/
    grPythiaMinus->SetMarkerColor(kBlue);
    grPythiaMinus->SetMarkerStyle(27);
    grPythiaMinus->SetMarkerSize(2.5);
    gr_ppMinusMCFMnlo->SetMarkerColor(kBlue);
    gr_ppMinusMCFMnlo->SetMarkerStyle(28);
    gr_ppMinusMCFMnlo->SetMarkerSize(2.3);
    grMSTWnloMinus->SetMarkerColor(kBlue);
    grMSTWnloMinus->SetMarkerStyle(28);
    grMSTWnloMinus->SetMarkerSize(2.3);
/*    gr_ppMinusMCFMnlo->SetLineWidth(4);
    gr_ppMinusMCFMnlo->SetFillColor(kBlue);
    gr_ppMinusMCFMnlo->SetFillStyle(3001);
*/
    TGraphAsymmErrors* graphWSumDummy = (TGraphAsymmErrors*)graphWSumc->Clone("graphWSumDummy");
    graphWSumDummy->SetMarkerSize(2.2);
    graphWSumDummy->SetMarkerStyle(33);
    graphWSumDummy->SetMarkerColor(kWhite);
    graphWSumc->SetMarkerSize(1.9);
    graphWSumc->SetMarkerStyle(33);
    graphWSumSystc->SetFillColor(kGray+1);
    graphWSumSyst_correlatedc->SetFillColor(kGray+1);
    graphWSumSyst_correlatedc->SetFillStyle(3244);

    TGraphAsymmErrors* graphWPlusSumDummy = (TGraphAsymmErrors*)graphWPlusSumc->Clone("graphWPlusSumDummy");
    graphWPlusSumDummy->SetMarkerSize(1.9);
    graphWPlusSumDummy->SetMarkerStyle(23);
    graphWPlusSumDummy->SetMarkerColor(kWhite);
    graphWPlusSumc->SetMarkerColor(kRed);
    graphWPlusSumc->SetMarkerStyle(23);
    graphWPlusSumc->SetMarkerSize(1.4);
    graphWPlusSumSystc->SetFillColor(kRed-7);
    graphWPlusSumSyst_correlatedc->SetFillColor(kRed);
    graphWPlusSumSyst_correlatedc->SetFillStyle(3004);

    TGraphAsymmErrors* graphWMinusSumDummy = (TGraphAsymmErrors*)graphWMinusSumc->Clone("graphWMinusSumDummy");
    graphWMinusSumDummy->SetMarkerSize(1.9);
    graphWMinusSumDummy->SetMarkerStyle(22);
    graphWMinusSumDummy->SetMarkerColor(kWhite);
    graphWMinusSumc->SetMarkerColor(kBlue);
    graphWMinusSumc->SetMarkerStyle(22);
    graphWMinusSumc->SetMarkerSize(1.4);
    graphWMinusSumSystc->SetFillColor(kBlue-7);
    graphWMinusSumSyst_correlatedc->SetFillColor(kBlue);
    graphWMinusSumSyst_correlatedc->SetFillStyle(3005);

///uncomment to draw systematic curves

/*    double ymin=0.12,ymax=50.0;
    graphWSumSyst_correlatedc->GetYaxis()->SetRangeUser(ymin,ymax);
    graphWSumSyst_correlatedc->GetXaxis()->SetLabelSize(0.036);
    graphWSumSyst_correlatedc->GetXaxis()->SetTitle("#LT N_{part} #GT");
    //graphWSumSyst_correlatedc->GetYaxis()->SetTitle("R_{CP}");
    graphWSumSyst_correlatedc->GetYaxis()->SetTitle("#frac{10^{9}}{#LT N_{coll} #GT}#frac{N^{W#rightarrow#mu,|#eta|<2.4}}{N_{events}}");
    graphWSumSyst_correlatedc->GetYaxis()->SetLabelSize(0.03);
    graphWSumSyst_correlatedc->GetYaxis()->CenterTitle(true);
    graphWSumSyst_correlatedc->GetYaxis()->SetTitleSize(0.04);
*/    
/*    graphWSumc->GetYaxis()->SetRangeUser(0.6,1.85);
    graphWSumc->GetXaxis()->SetLabelSize(0.036);
    graphWSumc->GetXaxis()->SetTitle("#LT N_{part} #GT");
    graphWSumc->GetYaxis()->SetTitle("R_{CP}");
*/
    double nTot = 10702.0, nTotPlus = 5314.0, nTotMinus = 5388.0;
    double nEv = 1.03e9;
    double ncoll_0to80 = 452.0;
    double norm = 1.0e9/(0.8*nEv*ncoll_0to80);
    nTot*=norm; nTotPlus*=norm; nTotMinus*=norm;
    TF1* funcTot = new TF1("funcTot","[0]",0.0,420.0);
    funcTot->SetLineStyle(kDashed);
    funcTot->SetParameter(0,nTot);
    funcTot->FixParameter(0,nTot);
    TF1* funcTotPlus = new TF1("funcTotPlus","[0]",0.0,420.0);
    funcTotPlus->SetLineStyle(kDashed);
    funcTotPlus->SetParameter(0,nTotPlus);
    funcTotPlus->FixParameter(0,nTotPlus);
    TF1* funcTotMinus = new TF1("funcTotMinus","[0]",0.0,420.0);
    funcTotMinus->SetLineStyle(kDashed);
    funcTotMinus->SetParameter(0,nTotMinus);
    funcTotMinus->FixParameter(0,nTotMinus);
/*    TF1 f0b = TF1("f0b","[0]",0,1);
    f0b.SetParameter(0,1.0);
    f0b.FixParameter(0,1.0);
    f0b.SetLineStyle(kBlack); //f0b.SetLineColor(kYellow);
*/    
//    gStyle->SetHatchesSpacing(0.5);
    gStyle->SetHatchesLineWidth(1);

///uncomment to draw systematic curves
/*
    graphWSumSyst_correlatedc->Draw("ae5");
    graphWPlusSumSyst_correlatedc->Draw("e5same");
    graphWMinusSumSyst_correlatedc->Draw("e5same");
*/
//    graphWSumc->Draw("pesame");

    double ymin=0.0,ymax=37.0;
    hDummy->GetYaxis()->SetRangeUser(ymin,ymax);
    hDummy->GetXaxis()->SetLabelSize(0.036);
    hDummy->GetXaxis()->SetTitle("#LT N_{part} #GT");
    //hDummy->GetYaxis()->SetTitle("R_{CP}");
    hDummy->GetYaxis()->SetTitle("#frac{10^{9}}{#LT N_{coll} #GT}#frac{N^{W#rightarrow#mu,fiducial}}{N_{events}}");
    hDummy->GetYaxis()->SetLabelSize(0.03);
    //hDummy->GetYaxis()->CenterTitle(true);
    hDummy->GetYaxis()->SetTitleSize(0.04);
 
    graphWSumc->SetMarkerSize(1.8);
    hDummy->Draw("0");
    //funcTot->Draw("same");
    funcTotPlus->SetLineColor(kRed);
    //funcTotPlus->Draw("same");
    funcTotMinus->SetLineColor(kBlue);
    //funcTotMinus->Draw("same");
    TFile* outFileTheory = new TFile("theory_npart.root","recreate");
    funcPythia->Write();
    funcPythia->Draw("same");
    //funcPythiaPlus->Draw("same");
    //funcPythiaMinus->Draw("same");
    funcMSTWnlo->Draw("same");
    //funcMSTWnloPlus->Draw("same");
    //funcMSTWnloMinus->Draw("same");
    graphWSumc->Draw("pesame");

    graphWPlusSumc->SetMarkerSize(1.8);

    graphWMinusSumc->SetMarkerSize(1.8);

/*    graphWSumSyst_correlatedc->Draw("[]");
    graphWPlusSumSyst_correlatedc->Draw("[]");
    graphWMinusSumSyst_correlatedc->Draw("[]");
*/
    hGenAllNuc->SetFillColor(kYellow);
    hGen_nppn->SetFillColor(kRed-7);
    hGen_nn->SetFillColor(kAzure-9);
    hGen_pp->SetFillColor(kGray+1);
    hGenAllNuc->Scale(graphWSumc->GetY()[3]/hGenAllNuc->GetMaximum());
    hGen_nppn->Scale(0.478*hGenAllNuc->Integral()/hGen_nppn->Integral());
    hGen_nn->Scale(0.367*hGenAllNuc->Integral()/hGen_nn->Integral());
    hGen_pp->Scale(0.155*hGenAllNuc->Integral()/hGen_pp->Integral());

/*    hGenAllNuc->Draw("histsame");
    hGen_nppn->Draw("histsame");
    hGen_nn->Draw("histsame");
    hGen_pp->Draw("histsame");
    */
    //grPythia->Draw("psame");
    //grPythiaPlus->Draw("psame");
    //grPythiaMinus->Draw("psame");
    //gr_ppMCFMnlo->Draw("psame");
    //gr_ppPlusMCFMnlo->Draw("psame");
    //gr_ppMinusMCFMnlo->Draw("psame");
    //grMSTWnlo->Draw("psame");
    //grMSTWnloPlus->Draw("psame");
    //grMSTWnloMinus->Draw("psame");
    graphWSumSystc->Draw("e2 same");
    graphWPlusSumSystc->Draw("e2 same");
    graphWMinusSumSystc->Draw("e2 same");
    graphWSumSyst_correlatedc->Draw("e2");
    graphWPlusSumSyst_correlatedc->Draw("e2 same");
    graphWMinusSumSyst_correlatedc->Draw("e2 same");
    graphWSumDummy->Draw("pesame");
    graphWPlusSumDummy->Draw("pesame");
    graphWMinusSumDummy->Draw("pesame");
    graphWSumc->Draw("pesame");
    graphWPlusSumc->Draw("pesame");
    graphWMinusSumc->Draw("pesame");
    hDummy->Draw("sameaxis");

/*    std::cout << "Straight line hypothesis." << std::endl;
    graphWSumSystc->Fit("f0a");
    cout << "chi2 = " << f0a.GetChisquare() << "/" << f0a.GetNDF() << ", p = " << f0a.GetProb() << endl;
    f0a.Draw("same") ; 

    graphWPlusSumSystc->Fit("f0a");
    std::cout << "Mu+ " << std::endl;
    cout << "chi2 = " << f0a.GetChisquare() << "/" << f0a.GetNDF() << ", p = " << f0a.GetProb() << endl;
    f0a.Draw("same") ; 

    graphWMinusSumSystc->Fit("f0a");
    std::cout << "Mu- " << std::endl;
    cout << "chi2 = " << f0a.GetChisquare() << "/" << f0a.GetNDF() << ", p = " << f0a.GetProb() << endl;
    f0a.Draw("same") ; 
*/
    TH1F* hdummy0 = new TH1F(); hdummy0->SetLineColor(0);  hdummy0->SetFillColor(kGray+1); hdummy0->SetMarkerStyle(33);
    TH1F* hdummy1 = new TH1F(); hdummy1->SetLineColor(0); hdummy1->SetFillColor(kRed-7);
        hdummy1->SetLineColor(0); hdummy1->SetMarkerColor(kRed); hdummy1->SetMarkerStyle(23);  
    TH1F* hdummy2 = new TH1F(); hdummy2->SetLineColor(0);  hdummy2->SetFillColor(kBlue-7); 
        hdummy2->SetLineColor(0);  hdummy2->SetMarkerColor(kBlue); hdummy2->SetMarkerStyle(22);


   TLegend* leg = new TLegend(0.1694631,0.1730769,0.4848993,0.3741259,NULL,"brNDC");
   leg->SetTextFont(gStyle->GetTextFont());
   leg->SetTextSize(gStyle->GetTextSize());
   leg->SetBorderSize(0);
   leg->SetFillColor(0);
   leg->SetTextSize(0.03321678);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(0);
   TLegendEntry* entry=leg->AddEntry("","W^{#pm}(Data)","p0f");

   int ci;
   ci = TColor::GetColor("#999999");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineStyle(1);
   entry->SetLineWidth(0);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);

   
   entry=leg->AddEntry("funcPythia","W^{#pm}(PYTHIA LO*)","l");

   ci = TColor::GetColor("#666666");
   entry->SetLineColor(ci);
   entry->SetLineStyle(2);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("funcMSTWnlo","W^{#pm}(POWHEG NLO)","l");

   ci = TColor::GetColor("#666666");
   entry->SetLineColor(ci);
   entry->SetLineStyle(9);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();

   leg = new TLegend(0.3909396,0.3094406,0.692953,0.3741259,NULL,"brNDC");
   leg->SetTextFont(gStyle->GetTextFont());
   leg->SetTextSize(gStyle->GetTextSize());
   leg->SetBorderSize(0);
   leg->SetFillColor(0);
   leg->SetTextSize(0.03321678);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(0);

   TLegendEntry* entry2=leg->AddEntry("","W^{+}(Data)","pf");

   ci = TColor::GetColor("#ff6666");
   entry2->SetFillColor(ci);
   entry2->SetFillStyle(1001);
   entry2->SetLineStyle(1);
   entry2->SetLineWidth(0);

   ci = TColor::GetColor("#ff0000");
   entry2->SetMarkerColor(ci);
   entry2->SetMarkerStyle(23);
   entry2->SetMarkerSize(1.2);
   entry2->SetTextFont(42);
   leg->Draw();
   
   leg = new TLegend(0.6157718,0.3111888,0.9278523,0.3758741,NULL,"brNDC");
   leg->SetTextFont(gStyle->GetTextFont());
   leg->SetTextSize(gStyle->GetTextSize());
   leg->SetBorderSize(0);
   leg->SetFillColor(0);
   leg->SetTextSize(0.03321678);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(0);

   entry2=leg->AddEntry("","W^{-}(Data)","pf");

   ci = TColor::GetColor("#6666ff");
   entry2->SetFillColor(ci);
   entry2->SetFillStyle(1001);
   entry2->SetLineStyle(1);
   entry2->SetLineWidth(0);

   ci = TColor::GetColor("#0000ff");
   entry2->SetMarkerColor(ci);
   entry2->SetMarkerStyle(22);
   entry2->SetMarkerSize(1.2);
   entry2->SetTextFont(42);
   leg->Draw();




   TLatex* tex = new TLatex(0.5872483,0.1975524,"#int Ldt #approx 0.14 nb^{-1}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03321678);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.545302,0.2692308,"Pb+Pb #sqrt{s_{NN}} = 2.76 TeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03321678);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(19.30338,33.77392,"ATLAS Preliminary");
   tex->SetTextFont(72);
   tex->SetLineWidth(2);
   tex->Draw();
   

/*    TLegend *legFit = new TLegend(0.1946309,0.1818182,0.4194631,0.3426573,NULL,"brNDC");
    legFit->SetTextFont(gStyle->GetTextFont());
    legFit->SetTextSize(0.02972028);
    legFit->SetBorderSize(0);
    legFit->SetFillColor(0);
    legFit->AddEntry(hdummy0, "W^{#pm}(Data)", "pf");
    legFit->AddEntry(funcPythia,"W^{#pm}(PYTHIA LO*)","l");
    legFit->AddEntry(funcMSTWnlo,"W^{#pm}(POWHEG NLO)","l");
    //legFit->AddEntry(hdummy1, "W^{+}#rightarrow#mu^{+}(Data)", "pf");
    //legFit->AddEntry(hdummy2, "W^{-}#rightarrow#mu^{-}(Data)", "pf");
    legFit->Draw();
    */
    cRcp->Update();
    pad1->cd();

   tex = new TLatex(-25.48387,52.60505,"centrality[%]");
   tex->SetTextAlign(13);
   tex->SetTextFont(42);
   tex->SetTextSize(0.025);
   tex->SetLineWidth(2);
   //tex->Draw();
   tex = new TLatex(0.2902685,0.3618881,"#int Ldt #approx 0.14 nb^{-1}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03146853);
   tex->SetLineWidth(2);
//   tex->Draw();
   tex = new TLatex(0.5184564,0.3636364,"Pb+Pb #sqrt{s_{NN}} = 2.76 TeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03846154);
   tex->SetLineWidth(2);
//   tex->Draw();
    tex = new TLatex(7.707081,32.49048,"ATLAS Preliminary");
   tex->SetTextFont(72);
   tex->SetLineWidth(2);
//   tex->Draw();

/*   legFit = new TLegend(0.465,0.288,0.693,0.342,NULL,"brNDC");
   legFit->SetTextFont(gStyle->GetTextFont());
   legFit->SetBorderSize(0);
   legFit->SetFillColor(0);
   legFit->SetTextSize(0.02622377);
   legFit->AddEntry(hdummy1, "W^{+}(Data)", "pf");
   //legFit->AddEntry(funcPythia,"W^{#pm}(PYTHIA LO*)","l");
   //legFit->AddEntry(funcPythiaPlus,"W^{+}(PYTHIA LO*)","l");
   //legFit->AddEntry(funcPythiaMinus,"W^{-}(PYTHIA LO*)","l");
   legFit->Draw();
   legFit = new TLegend(0.465,0.288,0.693,0.342,NULL,"brNDC");
   legFit->SetTextFont(gStyle->GetTextFont());
   legFit->SetBorderSize(0);
   legFit->SetFillColor(0);
   legFit->SetTextSize(0.02622377);
   legFit->AddEntry(hdummy2, "W^{-}#rightarrow#mu^{-}(Data)", "pf");
   //legFit->AddEntry(funcMSTWnlo,"W^{#pm}(POWHEG NLO)","l");
   //legFit->AddEntry(funcMSTWnloPlus,"W^{+}(POWHEG NLO)","l");
   //legFit->AddEntry(funcMSTWnloMinus,"W^{-}(POWHEG NLO)","l");
//   legFit->AddEntry(hGen_nn,"nn","f");
//   legFit->AddEntry(hGen_pp,"pp","f");
   legFit->Draw();
*/
    cRcp->Print("HI_Rcp,"+sCharge.ReplaceAll("#","")+".png");
    cRcp->Print("HI_Rcp,"+sCharge.ReplaceAll("#","")+".pdf");
    cRcp->Print("HI_Rcp,"+sCharge.ReplaceAll("#","")+".eps");
    cRcp->Print("HI_Rcp,"+sCharge.ReplaceAll("#","")+".root");
}//plot Rcp

void RcpWHI(){
  cout << "Starting Rcp for W Analysis" << endl;
  int ptmin = 2; /// bin
  int ptmax = 5; /// bin
  int etamin = 0; /// bin
  int etamax = 6; /// bin
  int nWSigPdf = 1000000; /// what W template to use and corresponding fit result
  bool doSummaryPlot = true;
  bool doCentralityPlots = false;
  bool doPtPlots = false;
  bool doChargePlots = true;
  bool doEtaCent = false;
  double lowBound  = 0;
  double highBound = 1;

  //yield spreadsheet names	
  TString spreadSheetNameMuPlus = "dataSpreadSheetCentMuPlus.csv";
  TString spreadSheetNameMuMinus = "dataSpreadSheetCentMuMinus.csv";

  TString spreadSheetRcpNameMuInclusive = "dataSpreadSheetRcpInclusive.csv";
  TString spreadSheetRcpNameMuPlus = "dataSpreadSheetRcpPlus.csv";
  TString spreadSheetRcpNameMuMinus = "dataSpreadSheetRcpMinus.csv";

  std::ofstream spreadSheetMuPlus;
  spreadSheetMuPlus.open(spreadSheetNameMuPlus);

  std::ofstream spreadSheetMuMinus;
  spreadSheetMuMinus.open(spreadSheetNameMuMinus);


  std::ofstream spreadSheetRcpMuInclusive;
  spreadSheetRcpMuInclusive.open(spreadSheetRcpNameMuInclusive);

  std::ofstream spreadSheetRcpMuMinus;
  spreadSheetRcpMuMinus.open(spreadSheetRcpNameMuMinus);

  std::ofstream spreadSheetRcpMuPlus;
  spreadSheetRcpMuPlus.open(spreadSheetRcpNameMuPlus);

  ///additional systematic errors from 
  ///inc(dec) Cw fit coefficients by 1 sigma
//  double systErr1SigmaMuPlus = 0.0023;
//  double systErr1SigmaMuMinus = 0.0031;
//  TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
  TString baseString = "/usatlas/u/tbales/scratch/";
  TString fileNameMcAllNuc = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.04.12.2013.root";
  TString fileNameMc_pp = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_pp.06.21.2013.root";
  TString fileNameMc_pn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_pn.06.21.2013.root";
  TString fileNameMc_np = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_np.06.21.2013.root";
  TString fileNameMc_nn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay_nn.06.21.2013.root";

  ///file holding Aw 
  TFile* fAwCent = new TFile("CorrectionFactorFiles/correctionFactorsAwCent.08.01.2013.root","read");

  TString fileNameWFit = "" ;
   
  //fileNameWFit = "WAnalysis_fitResultCentChrgEta.02.06.2013" ;
  //fileNameWFit = "WAnalysis_fitResultCentChrgEta.02.13.2013" ;
//  fileNameWFit = "WAnalysis_fitResultCentChrgEta.04.03.2013";
//  fileNameWFit = "ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr.06.03.2013";
//  fileNameWFit = "ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr.06.11.2013";
  fileNameWFit = 
  // Corrected with Py6
//  "ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr.08.26.2013";
  // Corrected with PowPy8
  //"ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr.10.03.2013";
   // binning to match Iwona's 
  //"ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr.11.27.2013";
  ///Nominal result
  "ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr_extrapolated.12.02.2013";
  // cross-check with FCal wting
  //"crossChecks/WAnalysis_fitResultCentChrgEta_FCalWtd.06.11.2014";
  // cross-check with Z data-overlaid sample 
  //"ResultFiles/WAnalysis_fitResultCentChrgEta_AbsEtaBSCorr_extrapolated_ZDA_05.18.2014";

  //////////////////////////////////////////////////
  //Use for correlation in isolation cut systematic//
  ///////////////////////////////////////////////////
//        "systematics/WAnalysis_fitResultCentChrgEta_IncIsoConeSizeSyst.07.22.2013";
//        "systematics/WAnalysis_fitResultCentChrgEta_LoosenIsoCut.07.22.2013";
//        "systematics/WAnalysis_fitResultCentChrgEta_2GeVMpt.07.14.2013";
//        "systematics/WAnalysis_fitResultCentChrgEta_4GeVMpt.07.14.2013";
//        "systematics/WAnalysis_fitResultCentChrgEta_MptSmearSystematics.07.20.2013";
//        "systematics/WAnalysis_fitResultCentChrgEta_2GeVMpt.08.04.2013";
//        "systematics/WAnalysis_fitResultCentChrgEta_4GeVMpt.08.04.2013";
//          "systematics/WAnalysis_fitResultCentChrgEta_NominalWithTauBkg.07.30.2013";
//          "systematics/WAnalysis_fitResultCentChrgEta_MCDrivenZBkg.07.30.2013";
//          "systematics/WAnalysis_fitResultCentChrgEta_RaaScaledQCDBkg.07.30.2013";
//          "systematics/WAnalysis_fitResultCentChrgEta_TauBkgNoEmbed.07.30.2013";


//  fileNameWFit = "WAnalysis_fitResultCentChrgEta_CorrectedBkgSub_9etaBins6CentralityBins2ChargeBins.04.24.2013";
//  TString fileNameWFit = "WAnalysis_fitResultCentChrgEta_NoCorrection.02.20.2013" ; // w/o correction
//  TString fileNameWFit = "WAnalysis_fitResultCentChrgEta_Raw.02.20.2013" ; // w/o correction or bkg substraction
   
  std::cout << "Reading data from file " << fileNameWFit << std::endl;
  TFile* fDataSet = new TFile(fileNameWFit+".root", "READ");
  if ( !fDataSet->IsOpen() ) {
    std::cout << fDataSet << " not found!" << std::endl;
    exit(0);
  }
  cout << "all files open" << endl;
  
  TFile* outFile = new TFile("binaryScalingDistribution.root","recreate");
  ///TGraphs of Aw as fcn of npart
  TGraphAsymmErrors* grAw = (TGraphAsymmErrors*)fAwCent->Get("grAwNpartDistroEta0");
  TGraphAsymmErrors* grAwPlus = (TGraphAsymmErrors*)fAwCent->Get("grAwNpartDistroPlusEta0");
  TGraphAsymmErrors* grAwMinus = (TGraphAsymmErrors*)fAwCent->Get("grAwNpartDistroMinusEta0");

  ///TGraph of centrality bins used for the fits
  TGraphAsymmErrors* grWCentrality = (TGraphAsymmErrors*) fDataSet->Get("WPtFit_centrality");
  const int nWCentrality = grWCentrality->GetN();  //number of centrality bins
  std::cout << "Centrality bins = " << nWCentrality << std::endl;
  double* centralityArr = (grWCentrality->GetX()); //get centrality
  double* centralityL 	= (grWCentrality->GetEXlow());
  double* centralityH 	= (grWCentrality->GetEXhigh());

  ///TGraph of eta bins used in the fits
  TGraphAsymmErrors* grWEta = (TGraphAsymmErrors*) fDataSet->Get("WPtFit_eta");
  const int nWEta = grWEta->GetN();  //number of eta bins
  cout << "Eta bins = " << nWEta << endl;
  double* etaArr = grWEta->GetX(); //get eta
  double* etaL 	= grWEta->GetEXlow();
  double* etaH 	= grWEta->GetEXhigh();
  const int nWGraphs = nWCentrality*nWEta;  

  TObjArray arrWSig = TObjArray(nWGraphs);
  TObjArray arrNObs = TObjArray(nWGraphs) ;
  TObjArray arrNBkg = TObjArray(nWGraphs) ;
  TObjArray arrWSigSyst = TObjArray(nWGraphs);
  TObjArray arrWSigSyst1 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst2 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst3 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst4 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst5 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst6 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst7 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst8 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst9 = TObjArray(nWGraphs);
  TObjArray arrWSigSyst10 = TObjArray(nWGraphs);

  TString baseW = "WPtFit";
  TString baseWSig = baseW;   baseWSig +="_sig_eta"; //sig+stat. 
  TString baseNObs = baseW ; baseNObs += "_nobs_eta" ; ///stat errors
  TString baseNBkg = baseW ; baseNBkg += "_nbkg_eta" ; ///stat errors
  TString baseWSigSyst = baseW;   baseWSigSyst +="_sigSyst_eta"; //sig+syst.
  TString baseWSigSyst1 = baseW;   baseWSigSyst1 +="_sigSyst1_eta"; ///syst. errors
  TString baseWSigSyst2 = baseW;   baseWSigSyst2 +="_sigSyst2_eta"; ///syst. errors
  TString baseWSigSyst3 = baseW;   baseWSigSyst3 +="_sigSyst3_eta"; ///syst. errors
  TString baseWSigSyst4 = baseW;   baseWSigSyst4 +="_sigSyst4_eta"; ///syst. errors
  TString baseWSigSyst5 = baseW;   baseWSigSyst5 +="_sigSyst5_eta"; ///syst. errors
  TString baseWSigSyst6 = baseW;   baseWSigSyst6 +="_sigSyst6_eta"; ///syst. errors
  TString baseWSigSyst7 = baseW;   baseWSigSyst7 +="_sigSyst7_eta"; ///syst. errors
  TString baseWSigSyst8 = baseW;   baseWSigSyst8 +="_sigSyst8_eta"; ///syst. errors
  TString baseWSigSyst9 = baseW;   baseWSigSyst9 +="_sigSyst9_eta"; ///syst. errors
  TString baseWSigSyst10 = baseW;   baseWSigSyst10 +="_sigSyst10_eta"; ///syst. errors
 
  ///loop over the TGraphs from the WAnalysis output
  ///that store stat. and syst. errors and yields
  for (int icent = 0; icent<nWCentrality; ++icent) {
    for (int ieta = 0; ieta<nWEta; ++ieta) {

       TString searchWI = baseWSig;        searchWI += ieta;      searchWI += "_centrality";
       TString searchNObsI = baseNObs;  searchNObsI += ieta;  searchNObsI += "_centrality";
       TString searchNBkgI = baseNBkg;  searchNBkgI += ieta;  searchNBkgI += "_centrality";
       TString searchSyst1WI = baseWSigSyst1; searchSyst1WI += ieta; searchSyst1WI += "_centrality";
       TString searchSyst2WI = baseWSigSyst2; searchSyst2WI += ieta; searchSyst2WI += "_centrality";
       TString searchSyst3WI = baseWSigSyst3; searchSyst3WI += ieta; searchSyst3WI += "_centrality";
       TString searchSyst4WI = baseWSigSyst4; searchSyst4WI += ieta; searchSyst4WI += "_centrality";
       TString searchSyst5WI = baseWSigSyst5; searchSyst5WI += ieta; searchSyst5WI += "_centrality";
       TString searchSyst6WI = baseWSigSyst6; searchSyst6WI += ieta; searchSyst6WI += "_centrality";
       TString searchSyst7WI = baseWSigSyst7; searchSyst7WI += ieta; searchSyst7WI += "_centrality";
       TString searchSyst8WI = baseWSigSyst8; searchSyst8WI += ieta; searchSyst8WI += "_centrality";
       TString searchSyst9WI = baseWSigSyst9; searchSyst9WI += ieta; searchSyst9WI += "_centrality";
       TString searchSyst10WI = baseWSigSyst10; searchSyst10WI += ieta; searchSyst10WI += "_centrality";
       TString searchSystWI = baseWSigSyst;        searchSystWI += ieta;      searchSystWI += "_centrality";

       TString searchWIJ = searchWI;         searchWIJ += icent;
       TString searchSystWIJ = searchSystWI; searchSystWIJ += icent;
       TString searchNObsIJ = searchNObsI;  searchNObsIJ += icent;
       TString searchNBkgIJ = searchNBkgI;  searchNBkgIJ += icent;
       TString searchSyst1WIJ = searchSyst1WI; searchSyst1WIJ += icent;
       TString searchSyst2WIJ = searchSyst2WI; searchSyst2WIJ += icent;
       TString searchSyst3WIJ = searchSyst3WI; searchSyst3WIJ += icent;
       TString searchSyst4WIJ = searchSyst4WI; searchSyst4WIJ += icent;
       TString searchSyst5WIJ = searchSyst5WI; searchSyst5WIJ += icent;
       TString searchSyst6WIJ = searchSyst6WI; searchSyst6WIJ += icent;
       TString searchSyst7WIJ = searchSyst7WI; searchSyst7WIJ += icent;
       TString searchSyst8WIJ = searchSyst8WI; searchSyst8WIJ += icent;
       TString searchSyst9WIJ = searchSyst9WI; searchSyst9WIJ += icent;
       TString searchSyst10WIJ = searchSyst10WI; searchSyst10WIJ += icent;

      ///create index for the TGraph
      const int index = indexIJ(icent,ieta,nWEta); 

      //get the graphs that contains signal yield with stat. errors
      arrWSig[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchWIJ));
      arrNObs[index] = ((TGraphAsymmErrors*) fDataSet->Get(searchNObsIJ)) ;
      arrNBkg[index] = ((TGraphAsymmErrors*) fDataSet->Get(searchNBkgIJ)) ;
      //get the graphs that contains signal yield with syst. errors
      arrWSigSyst[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSystWIJ));
      arrWSigSyst1[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst1WIJ));
      arrWSigSyst2[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst2WIJ));
      arrWSigSyst3[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst3WIJ));
      arrWSigSyst4[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst4WIJ));
      arrWSigSyst5[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst5WIJ));
      arrWSigSyst6[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst6WIJ));
      arrWSigSyst7[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst7WIJ));
      arrWSigSyst8[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst8WIJ));
      arrWSigSyst9[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst9WIJ));
      arrWSigSyst10[index]  = ((TGraphAsymmErrors*) fDataSet->Get(searchSyst10WIJ));

      int nWpt = ((TGraphAsymmErrors*)arrWSig[index])->GetN();

      for(int k = 0; k<nWpt; k++){

          double* ytemp = ((TGraphAsymmErrors*)arrWSig[index])->GetY();
          double* ytempNObs = ((TGraphAsymmErrors*)arrNObs[index])->GetY() ;
          double* ytempNBkg = ((TGraphAsymmErrors*)arrNBkg[index])->GetY() ;
	  double* ytempL = ((TGraphAsymmErrors*)arrWSig[index])->GetEYlow();
	  double* ytempH = ((TGraphAsymmErrors*)arrWSig[index])->GetEYhigh();
	  double* ytempSystL = ((TGraphAsymmErrors*)arrWSigSyst[index])->GetEYlow();
	  double* ytempSystH = ((TGraphAsymmErrors*)arrWSigSyst[index])->GetEYhigh();
          double* ytempSyst1 = ((TGraphAsymmErrors*)arrWSigSyst1[index])->GetEYhigh();
          double* ytempSyst2 = ((TGraphAsymmErrors*)arrWSigSyst2[index])->GetEYhigh();
          double* ytempSyst3 = ((TGraphAsymmErrors*)arrWSigSyst3[index])->GetEYhigh();
          double* ytempSyst4 = ((TGraphAsymmErrors*)arrWSigSyst4[index])->GetEYhigh();
          double* ytempSyst5 = ((TGraphAsymmErrors*)arrWSigSyst5[index])->GetEYhigh();
          double* ytempSyst6 = ((TGraphAsymmErrors*)arrWSigSyst6[index])->GetEYhigh();
          double* ytempSyst7 = ((TGraphAsymmErrors*)arrWSigSyst7[index])->GetEYhigh();
          double* ytempSyst8 = ((TGraphAsymmErrors*)arrWSigSyst8[index])->GetEYhigh();
          double* ytempSyst9 = ((TGraphAsymmErrors*)arrWSigSyst9[index])->GetEYhigh();
          double* ytempSyst10 = ((TGraphAsymmErrors*)arrWSigSyst10[index])->GetEYhigh();

	      ///for W^{+},W^{-} will be used (i.e. 102, 103)
/*	      if(ytemp[k]>0 && (abs(ytempL[k])>0 || ytempH[k]>0)){
            std::cout << " found " << searchWIJ  << "\n" << searchSystWIJ << "\n" 
                    << searchSyst1WIJ << "\n" 
                    << searchSyst2WIJ << "\n" 
                    << searchSyst3WIJ << "\n" 
                    << searchSyst4WIJ << "\n" 
                    << searchSyst5WIJ << "\n" 
                    << searchSyst6WIJ << "\n" 
                    << searchSyst7WIJ << "\n" 
                    << searchSyst8WIJ << "\n" 
                    << searchSyst9WIJ << "\n" 
                    << searchSyst10WIJ << "\n" 
                    << std::endl;
	      }
*/
      } //nWpt
    } //nWCentrality
  } //nWEta

  //return the number of points in the TGraph of the most central bins 
  int nWpt = ((TGraphAsymmErrors*)arrWSig[0])->GetN();
  cout << "nWpt for arrWSig[0]: " << nWpt << endl;
  double* etaLow  = new double[nWEta]; 
  double* etaHigh = new double[nWEta];
  double* centralityLow  = new double[nWCentrality]; 
  double* centralityHigh = new double[nWCentrality];
  
  /// WAnalysis
  double* yW     = new double[nWEta*nWCentrality*nWpt];
  double* yNObs = new double[nWEta*nWCentrality*nWpt];
  double* yNBkg = new double[nWEta*nWCentrality*nWpt];
  double* eyWL   = new double[nWEta*nWCentrality*nWpt];
  double* eyWH   = new double[nWEta*nWCentrality*nWpt];
  double* eyWSystL   = new double[nWEta*nWCentrality*nWpt];
  double* eyWSystH   = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst1 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst2 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst3 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst4 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst5 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst6 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst7 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst8 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst9 = new double[nWEta*nWCentrality*nWpt];
  double* eyWSyst10 = new double[nWEta*nWCentrality*nWpt];

  double* yW0   =new double[nWpt];
  double* y0NObs =new double[nWpt];
  double* y0NBkg =new double[nWpt];
  double* yW0L  =new double[nWpt];
  double* yW0H  =new double[nWpt];     
  double* yWN0   =new double[nWpt];
  double* yWN0L  =new double[nWpt];
  double* yWN0H  =new double[nWpt];
  double* syWN0  =new double[nWpt];     

  double* yW0SystL =new double[nWpt];
  double* yW0SystH =new double[nWpt];
  double* yW0Syst1 =new double[nWpt];
  double* yW0Syst2 =new double[nWpt];
  double* yW0Syst3 =new double[nWpt];
  double* yW0Syst4 =new double[nWpt];
  double* yW0Syst5 =new double[nWpt];
  double* yW0Syst6 =new double[nWpt];
  double* yW0Syst7 =new double[nWpt];
  double* yW0Syst8 =new double[nWpt];
  double* yW0Syst9 =new double[nWpt];
  double* yW0Syst10=new double[nWpt];

  double* yWN     = new double[nWEta*nWCentrality*nWpt];
  double* eyWNL   = new double[nWEta*nWCentrality*nWpt];
  double* eyWNH   = new double[nWEta*nWCentrality*nWpt];
  double* syWNL   = new double[nWEta*nWCentrality*nWpt];
  double* syWNH   = new double[nWEta*nWCentrality*nWpt];
  
  /// summing up
  double* yWSum = new double[nWCentrality];
  double* yNObsSum = new double[nWCentrality];
  double* yNBkgSum = new double[nWCentrality];
  double* eyWSumL = new double[nWCentrality];
  double* eyWSumH = new double[nWCentrality];
  double* eyWSumSystL = new double[nWCentrality];
  double* eyWSumSystH = new double[nWCentrality];
  double* eyWSumSyst1 = new double[nWCentrality];
  double* eyWSumSyst2 = new double[nWCentrality];
  double* eyWSumSyst3 = new double[nWCentrality];
  double* eyWSumSyst4 = new double[nWCentrality];
  double* eyWSumSyst5 = new double[nWCentrality];
  double* eyWSumSyst6 = new double[nWCentrality];
  double* eyWSumSyst7 = new double[nWCentrality];
  double* eyWSumSyst8= new double[nWCentrality]; 
  double* eyWSumSyst9= new double[nWCentrality]; 
  double* eyWSumSyst10= new double[nWCentrality]; 

  double* yWPlusSum = new double[nWCentrality];
  double* yNObsPlusSum = new double[nWCentrality];
  double* yNBkgPlusSum = new double[nWCentrality];
  double* eyWPlusSumL = new double[nWCentrality];
  double* eyWPlusSumH = new double[nWCentrality];
  double* eyWPlusSumSystL = new double[nWCentrality];
  double* eyWPlusSumSystH = new double[nWCentrality];
  double* eyWPlusSumSyst1 = new double[nWCentrality];
  double* eyWPlusSumSyst2 = new double[nWCentrality];
  double* eyWPlusSumSyst3 = new double[nWCentrality];
  double* eyWPlusSumSyst4 = new double[nWCentrality];
  double* eyWPlusSumSyst5 = new double[nWCentrality];
  double* eyWPlusSumSyst6 = new double[nWCentrality];
  double* eyWPlusSumSyst7 = new double[nWCentrality];
  double* eyWPlusSumSyst8 = new double[nWCentrality];
  double* eyWPlusSumSyst9 = new double[nWCentrality];
  double* eyWPlusSumSyst10 = new double[nWCentrality];
  double* syWPlusSum = new double[nWCentrality];
  double* syRelWPlusSum = new double[nWCentrality];
  double* syWPlusSum2 = new double[nWCentrality]; 
  double* syWPlusSum2_correlated = new double[nWCentrality]; 
  double* syWPlusSumTotal = new double[nWCentrality];
  double* syWPlusSumTotal_correlated = new double[nWCentrality]; 

  double* yWMinusSum = new double[nWCentrality];
  double* yNObsMinusSum = new double[nWCentrality];
  double* yNBkgMinusSum = new double[nWCentrality];
  double* eyWMinusSumL = new double[nWCentrality];
  double* eyWMinusSumH = new double[nWCentrality];
  double* eyWMinusSumSystL = new double[nWCentrality];
  double* eyWMinusSumSystH = new double[nWCentrality];
  double* eyWMinusSumSyst1 = new double[nWCentrality];
  double* eyWMinusSumSyst2 = new double[nWCentrality];
  double* eyWMinusSumSyst3 = new double[nWCentrality];
  double* eyWMinusSumSyst4 = new double[nWCentrality];
  double* eyWMinusSumSyst5 = new double[nWCentrality];
  double* eyWMinusSumSyst6 = new double[nWCentrality];
  double* eyWMinusSumSyst7 = new double[nWCentrality];
  double* eyWMinusSumSyst8 = new double[nWCentrality];
  double* eyWMinusSumSyst9 = new double[nWCentrality];
  double* eyWMinusSumSyst10 = new double[nWCentrality];
  double* syWMinusSum = new double[nWCentrality];
  double* syRelWMinusSum = new double[nWCentrality];
  double* syWMinusSum2 = new double[nWCentrality]; 
  double* syWMinusSum2_correlated = new double[nWCentrality]; 
  double* syWMinusSumTotal = new double[nWCentrality];
  double* syWMinusSumTotal_correlated = new double[nWCentrality]; 

  double* syWSumL = new double[nWCentrality];
  double* syWSumLTemp = new double[nWCentrality];
  double* syWSumL2 = new double[nWCentrality]; 
  double* syWSumL2_correlated = new double[nWCentrality]; 
  double* syWSumLTotal = new double[nWCentrality];
  double* syWSumLTotal_correlated = new double[nWCentrality]; 

  double* syWSumH = new double[nWCentrality];
  double* syWSumHTemp = new double[nWCentrality];
  double* syWSumH2 = new double[nWCentrality];
  double* syWSumH2_correlated = new double[nWCentrality];
  double* syWSumHTotal = new double[nWCentrality];
  double* syWSumHTotal_correlated = new double[nWCentrality];

  double syWSumTemplatesEtaCent[6]  = {0.0	,0.0	,0.0	,0.0	,0.0	,0.0	};
  double syWSumTemplatesCentChrg[6] = {0.0	,0.0	,0.0	,0.0	,0.0	,0.0	};
  double* yWNSum = new double[nWCentrality];
  double* eyWNSumL = new double[nWCentrality];
  double* eyWNSumH = new double[nWCentrality];
  double* syWNSumL = new double[nWCentrality];
  double* syWNSumH = new double[nWCentrality];

  for (int i=0; i < nWCentrality; ++i){
    centralityLow[i] = centralityArr[i]-centralityL[i];
    centralityHigh[i] = centralityArr[i]+centralityH[i];
    std::cout << "CentralityLow = " << centralityArr[i] << " - " << centralityL[i] << std::endl;
    std::cout << "CentralityHigh = " << centralityArr[i] << " + " << centralityH[i] << std::endl;
  }
  for (int i=0; i < nWEta; ++i){
    etaLow[i] = etaArr[i]-etaL[i];
    etaHigh[i] = etaArr[i]+etaH[i];
    std::cout << "EtaLow = " << etaArr[i] << " - " << etaL[i] << std::endl;
    std::cout << "EtaHigh = " << etaArr[i] << " + " << etaH[i] << std::endl;
  }
  
  /// Corrections and systematics from rcol:
  TGraphErrors* grrcolRPC = rcolRPC();  //return TGraph of binary coll ratios
  TGraphErrors* grrcolSysRPC = rcolSysRPC();
  TGraphErrors* grrcolRCP = rcolRCP();
  TGraphErrors* grrcolSysRCP = rcolSysRCP();
  TGraphErrors* grNcoll = getNColl();
  double* rColPC  = grrcolRPC->GetY();  //get the vector of values
  double* erColPC  = grrcolSysRPC->GetY(); // uncertainty on rcolY in percent
  double* rColCP  = grrcolRCP->GetY();  //return RCP values
  double* erColCP  = grrcolSysRCP->GetY(); // uncertainty on rcolY in percent
  double* nColl = grNcoll->GetY(); ///return ncoll for each c.c.
  double* erNcoll = grNcoll->GetEY(); ///return error on ncoll in %
  double* yFactor  = new double[nWCentrality]; 
  double* eyFactor = new double[nWCentrality]; 
  double* yFactorBw  = new double[nWCentrality]; // explicit scaling for bin width
  double* eyFactorBw = new double[nWCentrality]; // explicit scaling for bin width
  double* yFactorBw2  = new double[nWCentrality]; // explicit scaling for bin width
  double* eyFactorBw2 = new double[nWCentrality]; // explicit scaling for bin width

  TGraphAsymmErrors* graphYFactor     = new TGraphAsymmErrors(nWCentrality); 

  for (int j=0; j < nWCentrality; ++j) {

    /*yFactor[j] = rColCP[nWCentrality-1]/rColCP[j] ; //ratio btwn Ncol_40-80% and Ncol_j
    eyFactor[j] = erColCP[j]*0.01*yFactor[j]; //error for each ratio 
    */
    yFactor[j] = 1.0/nColl[j];
    eyFactor[j] = erNcoll[j]*0.01*yFactor[j];
    cout << j << " centrality " << centralityArr[j] << " : " << rColCP[j] << " +- " << rColCP[j]*erColCP[j]*0.01 << "(" << erColCP[j] << "%), " <<
	 centralityL[j]+centralityH[j] <<"/" << (centralityL[nWCentrality-1]+centralityH[nWCentrality-1]) << " -> " << yFactor[j] << " +-" << eyFactor[j]<< endl;

    std::cout << "1 /" << nColl[j] << " = " << yFactor[j] << std::endl;
    //scaling to account for various centrality bin widths
    /*float relBinWidth = (centralityL[nWCentrality-1]+centralityH[nWCentrality-1])/(centralityL[j]+centralityH[j]);
    yFactorBw[j] = yFactor[j]*relBinWidth; //scale the more central bin width to that of the 40-80% bin
    eyFactorBw[j] = eyFactor[j]*relBinWidth;
    */

    float inverseBinWidth = 1.0/(centralityL[j]+centralityH[j]);
    yFactorBw[j] = yFactor[j]*inverseBinWidth; 
    eyFactorBw[j] = eyFactor[j]*inverseBinWidth;

    std::cout << "1.0 / " << centralityL[j]+centralityH[j] << " = " << inverseBinWidth << std::endl;
    std::cout << yFactor[j] << "*" << inverseBinWidth << " = " << yFactorBw[j] << std::endl; 

    const float nEventsSampled = 1.03; ///Number of events (*10^-9) sampled in run
    yFactorBw[j] /=nEventsSampled; //1/(ncoll*nevt*bw)
    std::cout << yFactor[j] << "*" << inverseBinWidth << " = " << yFactorBw[j] << std::endl; 
    eyFactorBw[j] /=nEventsSampled;
    ///Plot the RColl distribution
    graphYFactor->SetPoint(j, 1-centralityArr[j], yFactor[j] );
    graphYFactor->SetPointError(j, centralityL[j], centralityH[j], eyFactor[j], eyFactor[j] ); // systematic

  } //j

  graphYFactor->GetYaxis()->SetRangeUser(0.0,1.05);
  graphYFactor->GetYaxis()->SetTitle("N_{coll}^{40-80}/N_{coll}^{central}");
  graphYFactor->GetXaxis()->SetTitle("1-centrality");
  graphYFactor->Draw("AP");
  
  double arrWPlusTotal = 0.; 
  double arrWMinusTotal = 0.; 


  /// For a given centrality class, add up
  /// all the W candidates in each rapidity window 
  cout << "Building data arrays" << endl;
  for (int i = 0; i<nWCentrality; ++i) {
    for (int j = 0; j<nWEta; ++j) {
      const int index = indexIJ(i,j,nWEta); 

      if (arrWSig[index]) {
        yW0   = ((TGraphAsymmErrors*)arrWSig[index])->GetY();
        yW0L  = ((TGraphAsymmErrors*)arrWSig[index])->GetEYlow();
        yW0H  = ((TGraphAsymmErrors*)arrWSig[index])->GetEYhigh();      

	y0NObs = ((TGraphAsymmErrors*)arrNObs[index])->GetY();
        y0NBkg = ((TGraphAsymmErrors*)arrNBkg[index])->GetY();

	yW0SystL = ((TGraphAsymmErrors*)arrWSigSyst[index])->GetEYlow();
	yW0SystH = ((TGraphAsymmErrors*)arrWSigSyst[index])->GetEYhigh();
	yW0Syst1 = ((TGraphAsymmErrors*)arrWSigSyst1[index])->GetEYhigh();
	yW0Syst2 = ((TGraphAsymmErrors*)arrWSigSyst2[index])->GetEYhigh();
	yW0Syst3 = ((TGraphAsymmErrors*)arrWSigSyst3[index])->GetEYhigh();
	yW0Syst4 = ((TGraphAsymmErrors*)arrWSigSyst4[index])->GetEYhigh();
	yW0Syst5 = ((TGraphAsymmErrors*)arrWSigSyst5[index])->GetEYhigh();
	yW0Syst6 = ((TGraphAsymmErrors*)arrWSigSyst6[index])->GetEYhigh();
	yW0Syst7 = ((TGraphAsymmErrors*)arrWSigSyst7[index])->GetEYhigh();
	yW0Syst8 = ((TGraphAsymmErrors*)arrWSigSyst8[index])->GetEYhigh();
	yW0Syst9 = ((TGraphAsymmErrors*)arrWSigSyst9[index])->GetEYhigh();
	yW0Syst10 = ((TGraphAsymmErrors*)arrWSigSyst10[index])->GetEYhigh();

        for (int k = 0; k<nWpt; ++k) {
          const int index2 = indexIJK(i, j, k, nWEta, nWpt); 
          if (yW0[k]>1 && (abs(yW0L[k])>0 || yW0H[k]>0)) {
      	    cout << "index: " << index << endl;

            std::cout << "k: " << k << std::endl;
            std::cout << "nWpt: " << nWpt << std::endl;
            std::cout << "index2: " << index2 << std::endl;
            std::cout << "yW0Syst1: " << yW0Syst2[k] << std::endl; 
            yW[index2] = yW0[k];
	    yNObs[index2] = y0NObs[k]; 
	    yNBkg[index2] = y0NBkg[k]; 
	        ///stat. error (absolute) from input file
            eyWL[index2] = yW0L[k];
            eyWH[index2] = yW0H[k];
	        ///syst. error(absolute) from input file
            eyWSystL[index2] = yW0SystL[k];
            eyWSystH[index2] = yW0SystH[k];
            eyWSyst1[index2] = yW0Syst1[k];
            eyWSyst2[index2] = yW0Syst2[k];
            eyWSyst3[index2] = yW0Syst3[k];
            eyWSyst4[index2] = yW0Syst4[k];
            eyWSyst5[index2] = yW0Syst5[k];
            // remove tracking eff syst
            eyWSyst6[index2] = 0.0*yW0Syst6[k];
            eyWSyst7[index2] = yW0Syst7[k];
            eyWSyst8[index2] = yW0Syst8[k];
            eyWSyst9[index2] = yW0Syst9[k];
            eyWSyst10[index2] = yW0Syst10[k];

            if (eyWL[index2]==0) eyWL[index2] = eyWH[index2]; //sometimes the minos asymetric error fails to find low bound (?)
            if (eyWH[index2]==0) eyWH[index2] = eyWL[index2]; //sometimes the minos asymetric error fails to find high bound (?)

//            cout << i << ":" << j << ":" << k << ": " << index2 <<" W = " << yW0[k]<< "+"<<yW0H[k]<<"-"<<yW0L[k] << " stat. " << " +"<<yW0SystH[k]<<"-"<<yW0SystL[k]  << " syst." << endl;

	        if(k==102) arrWPlusTotal += yW0[k];
	        if(k==103) arrWMinusTotal += yW0[k];
          } 
          else {
            yW[index2] = 0;
	    yNObs[index2] = 0;
	    yNBkg[index2] = 0;
            eyWL[index2] = 0;
            eyWH[index2] = 0;
            eyWSystL[index2] = 0.;
            eyWSystH[index2] = 0.;
            eyWSyst1[index2] = 0.;
            eyWSyst2[index2] = 0.;
            eyWSyst3[index2] = 0.;
            eyWSyst4[index2] = 0.;
            eyWSyst5[index2] = 0.;
            eyWSyst6[index2] = 0.;
            eyWSyst7[index2] = 0.;
            eyWSyst8[index2] = 0.;
            eyWSyst9[index2] = 0.;
            eyWSyst10[index2] = 0.;
          }
        } 
      }
    } //eta
  } //centrality
       
  std::cout << "Total number of W^{+} in sample = " << arrWPlusTotal << std::endl;
  std::cout << "Total number of W^{-} in sample = " << arrWMinusTotal << std::endl;

  /// sum up per centrality bin
  TString sCentralityRange[] = {"0-5","5-10","10-15","15-20","20-40","40-80"};
  TString sNpart[] = {"45.9","157.8","239.5","281.9","330.3","382.2"};
  for (int icent = 0; icent<nWCentrality; ++icent) {
    yWSum[icent]         = 0.;
    yNObsSum[icent] = 0.;
    yNBkgSum[icent] = 0.;
    eyWSumL[icent]       = 0.;
    eyWSumH[icent]       = 0.;    
    eyWSumSystL[icent]       = 0.;
    eyWSumSystH[icent]       = 0.;    

    yWPlusSum[icent]         = 0.;
    yNObsPlusSum[icent] = 0.;
    yNBkgPlusSum[icent] = 0.;
    eyWPlusSumL[icent]       = 0.;
    eyWPlusSumH[icent]       = 0.;    
    eyWPlusSumSystL[icent]       = 0.;
    eyWPlusSumSystH[icent]       = 0.;    

    eyWSumSyst1[icent] = 0.;
    eyWSumSyst2[icent] = 0.;
    eyWSumSyst3[icent] = 0.;
    eyWSumSyst4[icent] = 0.;
    eyWSumSyst5[icent] = 0.;
    eyWSumSyst6[icent] = 0.;
    eyWSumSyst7[icent] = 0.;
    eyWSumSyst8[icent] = 0.;
    eyWSumSyst9[icent] = 0.;
    eyWSumSyst10[icent] = 0.;

    eyWPlusSumSyst1[icent] = 0.;
    eyWPlusSumSyst2[icent] = 0.;
    eyWPlusSumSyst3[icent] = 0.;
    eyWPlusSumSyst4[icent] = 0.;
    eyWPlusSumSyst5[icent] = 0.;
    eyWPlusSumSyst6[icent] = 0.;
    eyWPlusSumSyst7[icent] = 0.;
    eyWPlusSumSyst8[icent] = 0.;
    eyWPlusSumSyst9[icent] = 0.;
    eyWPlusSumSyst10[icent] = 0.;

    eyWMinusSumSyst1[icent] = 0.;
    eyWMinusSumSyst2[icent] = 0.;
    eyWMinusSumSyst3[icent] = 0.;
    eyWMinusSumSyst4[icent] = 0.;
    eyWMinusSumSyst5[icent] = 0.;
    eyWMinusSumSyst6[icent] = 0.;
    eyWMinusSumSyst7[icent] = 0.;
    eyWMinusSumSyst8[icent] = 0.;
    eyWMinusSumSyst9[icent] = 0.;
    eyWMinusSumSyst10[icent] = 0.;

    yWMinusSum[icent]         = 0.;
    yNObsMinusSum[icent] = 0.;
    yNBkgMinusSum[icent] = 0.;
    eyWMinusSumL[icent]       = 0.;
    eyWMinusSumH[icent]       = 0.;    
    eyWMinusSumSystL[icent]       = 0.;
    eyWMinusSumSystH[icent]       = 0.;    

    syWSumL[icent]       = 0.;
    syWSumLTemp[icent]       = 0.;
    syWSumL2[icent]  = 0.;
    syWSumL2_correlated[icent]  = 0.;
    syWSumLTotal[icent]  = 0.;
    syWSumLTotal_correlated[icent]  = 0.;

    //mu+
    syWPlusSum[icent]       = 0.;
    syRelWPlusSum[icent]       = 0.;
    syWPlusSum2[icent]  = 0.;
    syWPlusSumTotal[icent]  = 0.;
    syWPlusSum2_correlated[icent]  = 0.;    
    syWPlusSumTotal[icent]  = 0.;    
    syWPlusSumTotal_correlated[icent]  = 0.;    
    //mu-
    syWMinusSum[icent]       = 0.;
    syRelWMinusSum[icent]       = 0.;
    syWMinusSum2[icent]  = 0.;
    syWMinusSumTotal[icent]  = 0.;
    syWMinusSum2_correlated[icent]  = 0.;    
    syWMinusSumTotal[icent]  = 0.;    
    syWMinusSumTotal_correlated[icent]  = 0.;    

    syWSumH[icent]       = 0.;    
    syWSumHTemp[icent]       = 0.;    
    syWSumH2[icent]  = 0.;    
    syWSumH2_correlated[icent]  = 0.;    
    syWSumHTotal[icent]  = 0.;    
    syWSumHTotal_correlated[icent]  = 0.;    

    //std::cout << nWEta << std::endl;
    for (int ieta = 0; ieta<nWEta; ++ieta) { 
      for (int k = 102; k<=103; ++k) { 

        const int index2 = indexIJK(icent, ieta, k, nWEta, nWpt); 

	    ///add up W's per eta and centrality bin
        if (yW[index2]>0){
          yWSum[icent] += yW[index2]; //running sum of the W's obtained for the given cent bin
	  yNObsSum[icent] += yNObs[index2];
	  yNBkgSum[icent] += yNBkg[index2];
	  ///running sum of squares
          eyWSumL[icent] += sqr(eyWL[index2]); 
          eyWSumH[icent] += sqr(eyWH[index2]);
          eyWSumSyst1[icent] += TMath::Power(eyWSyst1[index2],2);
          eyWSumSyst2[icent] += TMath::Power(eyWSyst2[index2],2);
          eyWSumSyst3[icent] += TMath::Power(eyWSyst3[index2],2);
          eyWSumSyst4[icent] += TMath::Power(eyWSyst4[index2],2);
          eyWSumSyst5[icent] += TMath::Power(eyWSyst5[index2],2);
          // remove tracking eff syst
          eyWSumSyst6[icent] += TMath::Power(0.0*eyWSyst6[index2],2);
          eyWSumSyst7[icent] += TMath::Power(eyWSyst7[index2],2);
          eyWSumSyst8[icent] += TMath::Power(eyWSyst8[index2],2);
          eyWSumSyst9[icent] += TMath::Power(eyWSyst9[index2],2);
          eyWSumSyst10[icent] += TMath::Power(eyWSyst10[index2],2);

	  if(k==102) {
          	yWPlusSum[icent] += yW[index2]; //running sum of the W's obtained for the given cent bin
	        yNObsPlusSum[icent] += yNObs[index2];
	        yNBkgPlusSum[icent] += yNBkg[index2];
          	eyWPlusSumL[icent] += sqr(eyWL[index2]); 
          	eyWPlusSumH[icent] += sqr(eyWH[index2]);
          	eyWPlusSumSystL[icent] += sqr(eyWSystL[index2]); 
          	eyWPlusSumSystH[icent] += sqr(eyWSystH[index2]);
            eyWPlusSumSyst1[icent] += TMath::Power(eyWSyst1[index2],2);
            eyWPlusSumSyst2[icent] += TMath::Power(eyWSyst2[index2],2);
            eyWPlusSumSyst3[icent] += TMath::Power(eyWSyst3[index2],2);
            eyWPlusSumSyst4[icent] += TMath::Power(eyWSyst4[index2],2);
            eyWPlusSumSyst5[icent] += TMath::Power(eyWSyst5[index2],2);
            // remove tracking eff syst
            eyWPlusSumSyst6[icent] += TMath::Power(0.0*eyWSyst6[index2],2);
            eyWPlusSumSyst7[icent] += TMath::Power(eyWSyst7[index2],2);
            eyWPlusSumSyst8[icent] += TMath::Power(eyWSyst8[index2],2);
            eyWPlusSumSyst9[icent] += TMath::Power(eyWSyst9[index2],2);
            eyWPlusSumSyst10[icent] += TMath::Power(eyWSyst10[index2],2);
            //std::cout << "eyWPlusSumSyst1[icent] " << eyWSyst1[index2] << std::endl; exit(0);
            /*std::cout << std::endl;
            std::cout << "ieta " << ieta << std::endl;
            std::cout << "yWPlus[icent] " << yW[index2] << std::endl; 
            std::cout << "eyWPlusSumSyst10[icent] " << eyWSyst10[index2] << std::endl;  
            std::cout << std::endl;
            */
	  }

	  if(k==103) {
            yWMinusSum[icent] += yW[index2]; //running sum of the W's obtained for the given cent bin
	    yNObsMinusSum[icent] += yNObs[index2];
	    yNBkgMinusSum[icent] += yNBkg[index2];
            eyWMinusSumL[icent] += sqr(eyWL[index2]); 
            eyWMinusSumH[icent] += sqr(eyWH[index2]);
            eyWMinusSumSystL[icent] += sqr(eyWSystL[index2]); 
            eyWMinusSumSystH[icent] += sqr(eyWSystH[index2]);
            eyWMinusSumSyst1[icent] += TMath::Power(eyWSyst1[index2],2);
            eyWMinusSumSyst2[icent] += TMath::Power(eyWSyst2[index2],2);
            eyWMinusSumSyst3[icent] += TMath::Power(eyWSyst3[index2],2);
            eyWMinusSumSyst4[icent] += TMath::Power(eyWSyst4[index2],2);
            eyWMinusSumSyst5[icent] += TMath::Power(eyWSyst5[index2],2);
            // remove tracking eff syst
            eyWMinusSumSyst6[icent] += TMath::Power(0.0*eyWSyst6[index2],2);
            eyWMinusSumSyst7[icent] += TMath::Power(eyWSyst7[index2],2);
            eyWMinusSumSyst8[icent] += TMath::Power(eyWSyst8[index2],2);
            eyWMinusSumSyst9[icent] += TMath::Power(eyWSyst9[index2],2);
            eyWMinusSumSyst10[icent] += TMath::Power(eyWSyst10[index2],2);

            /*std::cout << std::endl;
            std::cout << "ieta " << ieta << std::endl;
            std::cout << "yWMinus[icent] " << yW[index2] << std::endl; 
            std::cout << "eyWMinusSumSyst10[icent] " << eyWSyst10[index2] << std::endl;  
            std::cout << std::endl;
            */
	  }

	  //sanity check
	  std::cout << "Current sum of W candidates in eta bin " << ieta << " and centrality " << icent 
                << " = " << yWPlusSum[icent] << " + " << yWMinusSum[icent] << " = " 
		        << yWSum[icent] << std::endl; 
        }
      }
    } //ieta
//    std::cout << "W+: " << yWPlusSum[icent] << "+" << " W- : " << yWMinusSum[icent] << " = " << yWSum[icent] << std::endl; exit(0);
    //calculate absolute stat and syst errors
    //stat.
    eyWSumL[icent] = sqrt(eyWSumL[icent]);
    eyWSumH[icent] = sqrt(eyWSumH[icent]);
    eyWPlusSumL[icent] = sqrt(eyWPlusSumL[icent]);
    eyWPlusSumH[icent] = sqrt(eyWPlusSumH[icent]);
    eyWMinusSumL[icent] = sqrt(eyWMinusSumL[icent]);
    eyWMinusSumH[icent] = sqrt(eyWMinusSumH[icent]);

    //syst.
	eyWSumSyst1[icent] = sqrt(eyWSumSyst1[icent]);
	eyWSumSyst2[icent] = sqrt(eyWSumSyst2[icent]);
	eyWSumSyst3[icent] = sqrt(eyWSumSyst3[icent]);
	eyWSumSyst4[icent] = sqrt(eyWSumSyst4[icent]);
	eyWSumSyst5[icent] = sqrt(eyWSumSyst5[icent]);
	eyWSumSyst6[icent] = sqrt(eyWSumSyst6[icent]);
	eyWSumSyst7[icent] = sqrt(eyWSumSyst7[icent]);
	eyWSumSyst8[icent] = sqrt(eyWSumSyst8[icent]);
	eyWSumSyst9[icent] = sqrt(eyWSumSyst9[icent]);
	eyWSumSyst10[icent] = sqrt(eyWSumSyst10[icent]);
    ///Now calculate the total correlated systematic in this bin for mu+-
	eyWSumSystH[icent] = sqrt(TMath::Power(eyWSumSyst1[icent],2)+TMath::Power(eyWSumSyst2[icent],2)+TMath::Power(eyWSumSyst3[icent],2)
        +TMath::Power(eyWSumSyst4[icent],2)+TMath::Power(eyWSumSyst5[icent],2)+TMath::Power(eyWSumSyst6[icent],2)
        +TMath::Power(eyWSumSyst7[icent],2)+TMath::Power(eyWSumSyst8[icent],2)+TMath::Power(eyWSumSyst9[icent],2));
	eyWSumSystL[icent] = eyWSumSystH[icent];
    
	eyWPlusSumSyst1[icent] = sqrt(eyWPlusSumSyst1[icent]);
	eyWPlusSumSyst2[icent] = sqrt(eyWPlusSumSyst2[icent]);
	eyWPlusSumSyst3[icent] = sqrt(eyWPlusSumSyst3[icent]);
	eyWPlusSumSyst4[icent] = sqrt(eyWPlusSumSyst4[icent]);
	eyWPlusSumSyst5[icent] = sqrt(eyWPlusSumSyst5[icent]);
	eyWPlusSumSyst6[icent] = sqrt(eyWPlusSumSyst6[icent]);
	eyWPlusSumSyst7[icent] = sqrt(eyWPlusSumSyst7[icent]);
	eyWPlusSumSyst8[icent] = sqrt(eyWPlusSumSyst8[icent]);
	eyWPlusSumSyst9[icent] = sqrt(eyWPlusSumSyst9[icent]);
	eyWPlusSumSyst10[icent] = sqrt(eyWPlusSumSyst10[icent]);

    ///Now calculate the total systematic in this bin for mu+
	eyWPlusSumSystH[icent] = sqrt(TMath::Power(eyWPlusSumSyst1[icent],2)+TMath::Power(eyWPlusSumSyst2[icent],2)+TMath::Power(eyWPlusSumSyst3[icent],2)
        +TMath::Power(eyWPlusSumSyst4[icent],2)+TMath::Power(eyWPlusSumSyst5[icent],2)+TMath::Power(eyWPlusSumSyst6[icent],2)
        +TMath::Power(eyWPlusSumSyst7[icent],2)+TMath::Power(eyWPlusSumSyst8[icent],2)+TMath::Power(eyWPlusSumSyst9[icent],2));
	eyWPlusSumSystL[icent] = eyWPlusSumSystH[icent];
    
	eyWMinusSumSyst1[icent] = sqrt(eyWMinusSumSyst1[icent]);
	eyWMinusSumSyst2[icent] = sqrt(eyWMinusSumSyst2[icent]);
	eyWMinusSumSyst3[icent] = sqrt(eyWMinusSumSyst3[icent]);
	eyWMinusSumSyst4[icent] = sqrt(eyWMinusSumSyst4[icent]);
	eyWMinusSumSyst5[icent] = sqrt(eyWMinusSumSyst5[icent]);
	eyWMinusSumSyst6[icent] = sqrt(eyWMinusSumSyst6[icent]);
	eyWMinusSumSyst7[icent] = sqrt(eyWMinusSumSyst7[icent]);
	eyWMinusSumSyst8[icent] = sqrt(eyWMinusSumSyst8[icent]);
	eyWMinusSumSyst9[icent] = sqrt(eyWMinusSumSyst9[icent]);
	eyWMinusSumSyst10[icent] = sqrt(eyWMinusSumSyst10[icent]);

    ///Now calculate the total correlated systematic in this bin for mu-
	eyWMinusSumSystH[icent] = sqrt(TMath::Power(eyWMinusSumSyst1[icent],2)+TMath::Power(eyWMinusSumSyst2[icent],2)+TMath::Power(eyWMinusSumSyst3[icent],2)
        +TMath::Power(eyWMinusSumSyst4[icent],2)+TMath::Power(eyWMinusSumSyst5[icent],2)+TMath::Power(eyWMinusSumSyst6[icent],2)
        +TMath::Power(eyWMinusSumSyst7[icent],2)+TMath::Power(eyWMinusSumSyst8[icent],2)+TMath::Power(eyWMinusSumSyst9[icent],2));
	eyWMinusSumSystL[icent] = eyWMinusSumSystH[icent];
    
    //add uncorrelated stat and syst in quadrature 
    syWSumL[icent] = sqrt( sqr(eyWSumSyst10[icent]) + sqr(eyWSumL[icent]) ) ;
//    std::cout << "Stat: " << eyWSumL[icent] << " Uncorrelated Syst: " << eyWSumSyst10[icent] << " Quad: " << syWSumL[icent] << std::endl; exit(0);//hack
    syWSumLTemp[icent] = sqrt( sqr(syWSumL[icent]/ yWSum[icent]*100.0)  ) ; 
    syWSumL2[icent] = sqrt( sqr(syWSumLTemp[icent]) ) * 0.01 ; //relative syst error of jth bin 
    ///relative correlated errors pf jth bin
    syWSumL2_correlated[icent] = TMath::Sqrt(TMath::Power(eyWSumSystL[icent]/yWSum[icent],2)+TMath::Power(erNcoll[icent]/100.0,2));

    //add uncorrelated stat and syst in quadrature 
    //syWSumH[icent] = sqrt( sqr(eyWSumSystL[icent]) + sqr(eyWSumH[icent]) ) ;
    syWSumH[icent] = sqrt( sqr(eyWSumSyst10[icent]) + sqr(eyWSumH[icent]) ) ;
    std::cout<<std::endl;
    std::cout << "eyWSumH[icent] " << eyWSumH[icent] << std::endl;
    std::cout << "eyWSumSyst10[icent] " << eyWSumSyst10[icent] << std::endl;
    std::cout << "syWSumH[icent] : " << syWSumH[icent] << std::endl;
    std::cout<<std::endl;
    
    syWSumHTemp[icent] = syWSumH[icent]/ yWSum[icent]; 
    syWSumH2[icent] = syWSumHTemp[icent]; //relative syst error of jth bin 
    std::cout << "Relative error of jth cent bin: " << syWSumH2[icent] << std::endl;
    ///relative correlated errors pf jth bin
    syWSumH2_correlated[icent] = TMath::Sqrt(TMath::Power(eyWSumSystL[icent]/yWSum[icent],2)+TMath::Power(erNcoll[icent]/100.0,2));
    
    ////////
    ///mu+//
    ////////
    std::cout << "mu+ uncorrelated systematic error in bin " << icent << " = " <<  eyWPlusSumSystL[icent] << std::endl;
    ///add uncorrelated absolute systematic error to statistical error
    syWPlusSum[icent] = sqrt( sqr(eyWPlusSumSyst10[icent]) + sqr(eyWPlusSumL[icent])) ;
    ///relative uncorrelated error
    syRelWPlusSum[icent] = sqrt( sqr(syWPlusSum[icent]/ yWPlusSum[icent]*100)) ; 
    syWPlusSum2[icent] = sqrt( sqr(syRelWPlusSum[icent]) ) * 0.01 ;     
    ///correlated relative systematic errors in jth bin
    syWPlusSum2_correlated[icent] = TMath::Sqrt(TMath::Power(eyWPlusSumSystL[icent]/yWPlusSum[icent],2)+TMath::Power(erNcoll[icent]/100.0,2));
    writeYieldsToSpreadsheet(spreadSheetMuPlus,sCentralityRange[icent],sNpart[icent], 
				format(yNObsPlusSum[icent],4),format(yNBkgPlusSum[icent],2),format(yWPlusSum[icent],4) ,
        			format(eyWPlusSumL[icent],2) , format(eyWPlusSumSyst10[icent],2), format(eyWPlusSumSystL[icent],2));

    ////////
    //mu-///
    ////////   
    std::cout << "mu- uncorrelated systematic error in bin " << icent << " = " <<  eyWMinusSumSystL[icent] << std::endl;
    writeYieldsToSpreadsheet(spreadSheetMuMinus,sCentralityRange[icent], sNpart[icent], 
				format(yNObsMinusSum[icent],4),format(yNBkgMinusSum[icent],2),format(yWMinusSum[icent],4) ,
        			format(eyWMinusSumL[icent],2) ,format(eyWMinusSumSyst10[icent],2), format(eyWMinusSumSystL[icent],2) );
    ///add uncorrelated absolute systematic error to statistical error
    syWMinusSum[icent] = sqrt( sqr(eyWMinusSumSyst10[icent]) + sqr(eyWMinusSumL[icent])) ;
    ///relative uncorrelated error
    syRelWMinusSum[icent] = sqrt( sqr(syWMinusSum[icent]/ yWMinusSum[icent]*100)) ; 
    syWMinusSum2[icent] = sqrt( sqr(syRelWMinusSum[icent]) ) * 0.01 ;     
    ///correlated relative systematic errors in jth bin
    syWMinusSum2_correlated[icent] = TMath::Sqrt(TMath::Power(eyWMinusSumSystL[icent]/yWMinusSum[icent],2)+TMath::Power(erNcoll[icent]/100.0,2));

    
    //yWRatio[icent] = yWPlusSum[icent]/yWMinusSum[icent];
  } ///icent
      
  //Here is where the Rcp (Rpc) is calculated
  float yWSumRef = yWSum[nWCentrality-1]; 
  float yWPlusSumRef = yWPlusSum[nWCentrality-1]; 
  float yWMinusSumRef = yWMinusSum[nWCentrality-1]; 

  for (int j = 0; j<nWCentrality; ++j) {

    ///calculate the bin width corrected normalization factor for RCP calculation
/*    float normW = 1./yWSumRef*yFactorBw[j] ; 
    float normWPlus = 1./yWPlusSumRef*yFactorBw[j] ; 
    float normWMinus = 1./yWMinusSumRef*yFactorBw[j] ; 
*/
    ///No division by most central bin
    float normW = yFactorBw[j] ; 
    float normWPlus = yFactorBw[j] ; 
    float normWMinus = yFactorBw[j] ; 

//    cout << "NormW = " << "1" <<"/" << "(" << yWSumRef << ") " << "*" << yFactorBw[j] << " = " << normW << endl;

    yWSum[j] *= normW;
    eyWSumL[j]  *= normW;
    eyWSumH[j]  *= normW;
    eyWSumSystL[j]  *= normW;
    eyWSumSystH[j]  *= normW;

    syWSumL[j]  *= normW; 
    syWSumL2[j] *= yWSum[j]; 
    syWSumLTotal[j] = syWSumL2[j]; 
    syWSumLTotal_correlated[j] = syWSumL2_correlated[j] ; 
    syWSumH[j]  *= normW;
    std::cout << "syWSumH2[j] " << syWSumH2[j] << std::endl;
    syWSumH2[j] *= yWSum[j];
    std::cout << "syWSumH2[j]*= yWSum[j] " << syWSumH2[j] << std::endl;
    syWSumHTotal[j] = syWSumH2[j] ; 
    syWSumHTotal_correlated[j] = syWSumH2_correlated[j] ; 

    ///////
    //mu+//
    ///////
    yWPlusSum[j] *= normWPlus;
    eyWPlusSumL[j]  *= normWPlus;
    eyWPlusSumH[j]  *= normWPlus;
    eyWPlusSumSystL[j]  *= normWPlus;
    eyWPlusSumSystH[j]  *= normWPlus;

    syWPlusSum[j]  *= normWPlus; 
    syWPlusSum2[j] *= yWPlusSum[j]; 
    syWPlusSumTotal[j] = syWPlusSum2[j]; 
    syWPlusSumTotal_correlated[j] = syWPlusSum2_correlated[j]; 

    ///////
    //mu-//
    ///////
    yWMinusSum[j] *= normWMinus;
    eyWMinusSumL[j]  *= normWMinus;
    eyWMinusSumH[j]  *= normWMinus;
    eyWMinusSumSystL[j]  *= normWMinus;
    eyWMinusSumSystH[j]  *= normWMinus;

    syWMinusSum[j]  *= normWMinus; 
    syWMinusSum2[j] *= yWMinusSum[j]; 
    syWMinusSumTotal[j] = syWMinusSum2[j]; 
    syWMinusSumTotal_correlated[j] = syWMinusSum2_correlated[j]; 

    cout << j <<" W average = " << yWSum[j]<< "+"<<eyWSumH[j]<<"-"<<eyWSumL[j] << " stat." << "+" <<eyWSumSystH[j]<<"-" <<eyWSumSystL[j] << endl;
  } //jth centrality
  
  cout << "Filled all arrays" << endl;
 
  if (doCentralityPlots || doSummaryPlot) {
    

    TGraphAsymmErrors* graphWSum      = new TGraphAsymmErrors(nWCentrality);
    TGraphAsymmErrors* graphWSumSyst  = new TGraphAsymmErrors(nWCentrality);
    TGraphAsymmErrors* graphWSumSyst_correlated  = new TGraphAsymmErrors(nWCentrality);

    TGraphAsymmErrors* graphWPlusSum      = new TGraphAsymmErrors(nWCentrality);
    TGraphAsymmErrors* graphWPlusSumSyst  = new TGraphAsymmErrors(nWCentrality);
    TGraphAsymmErrors* graphWPlusSumSyst_correlated  = new TGraphAsymmErrors(nWCentrality);

    TGraphAsymmErrors* graphWMinusSum      = new TGraphAsymmErrors(nWCentrality);
    TGraphAsymmErrors* graphWMinusSumSyst  = new TGraphAsymmErrors(nWCentrality);
    TGraphAsymmErrors* graphWMinusSumSyst_correlated  = new TGraphAsymmErrors(nWCentrality);

    cout << " Centrality \tRCP fit" << endl;
    double Npart[7] = {382.16,	330.26,	281.88,	239.52,	157.83,	45.93} ;
    double nparterr[7] = {0.5, 0.9, 1.3, 1.6, 2.6, 5.1, 10.1};
    double centralityBins[] = {0.0,0.05,0.1,0.15,0.20,0.40,0.80};
    double offset = 10.0;
    double arrNpart[] = {45.93+offset,157.83+offset,239.52+offset,281.88+offset,330.26+offset,382.16+30.0};
    double arrNpartBins[nWCentrality];
    arrNpartBins[0] = 0.0;
    for(int i=0; i<nWCentrality; ++i) {
        arrNpartBins[i+1] = 2.*arrNpart[i] - arrNpartBins[i];
        std::cout << arrNpartBins[i] << "-" << arrNpartBins[i+1] << std::endl;
    }

    ///Generator-level histos
    TH1F* hGenRcpAllNuc = new TH1F("hGenRcpAllNuc","hGenRcpAllNuc",nWCentrality,arrNpartBins);
    TH1F* hGenPlusRcpAllNuc = new TH1F("hGenPlusRcpAllNuc","hGenPlusRcpAllNuc",nWCentrality,arrNpartBins);
    TH1F* hGenMinusRcpAllNuc = new TH1F("hGenMinusRcpAllNuc","hGenMinusRcpAllNuc",nWCentrality,arrNpartBins);
    fillGeneratorRcp(fileNameMcAllNuc,centralityBins,nWCentrality,hGenRcpAllNuc,hGenPlusRcpAllNuc,hGenMinusRcpAllNuc);
    ///pp
    TH1F* hGenRcp_pp = new TH1F("hGenRcp_pp","hGenRcp_pp",nWCentrality,arrNpartBins);
    TH1F* hGenPlusRcp_pp = new TH1F("hGenPlusRcp_pp","hGenPlusRcp_pp",nWCentrality,arrNpartBins);
    TH1F* hGenMinusRcp_pp = new TH1F("hGenMinusRcp_pp","hGenMinusRcp_pp",nWCentrality,arrNpartBins);
    fillGeneratorRcp(fileNameMc_pp,centralityBins,nWCentrality,hGenRcp_pp,hGenPlusRcp_pp,hGenMinusRcp_pp);
    ///np(pn)
    TH1F* hGenRcp_np = new TH1F("hGenRcp_np","hGenRcp_np",nWCentrality,arrNpartBins);
    TH1F* hGenPlusRcp_np = new TH1F("hGenPlusRcp_np","hGenPlusRcp_np",nWCentrality,arrNpartBins);
    TH1F* hGenMinusRcp_np = new TH1F("hGenMinusRcp_np","hGenMinusRcp_np",nWCentrality,arrNpartBins);
    fillGeneratorRcp(fileNameMc_np,centralityBins,nWCentrality,hGenRcp_np,hGenPlusRcp_np,hGenMinusRcp_np);
    TH1F* hGenRcp_pn = new TH1F("hGenRcp_pn","hGenRcp_pn",nWCentrality,arrNpartBins);
    TH1F* hGenPlusRcp_pn = new TH1F("hGenPlusRcp_pn","hGenPlusRcp_pn",nWCentrality,arrNpartBins);
    TH1F* hGenMinusRcp_pn = new TH1F("hGenMinusRcp_pn","hGenMinusRcp_pn",nWCentrality,arrNpartBins);
    fillGeneratorRcp(fileNameMc_pn,centralityBins,nWCentrality,hGenRcp_pn,hGenPlusRcp_pn,hGenMinusRcp_pn);
    TH1F* hGenRcp_nppn = new TH1F("hGenRcp_nppn","hGenRcp_nppn",nWCentrality,arrNpartBins);
    hGenRcp_nppn->Add(hGenRcp_np,hGenRcp_pn);
    TH1F* hGenPlusRcp_nppn = new TH1F("hGenPlusRcp_nppn","hGenPlusRcp_nppn",nWCentrality,arrNpartBins);
    hGenPlusRcp_nppn->Add(hGenPlusRcp_np,hGenPlusRcp_pn);
    TH1F* hGenMinusRcp_nppn = new TH1F("hGenMinusRcp_nppn","hGenMinusRcp_nppn",nWCentrality,arrNpartBins);
    hGenMinusRcp_nppn->Add(hGenMinusRcp_np,hGenMinusRcp_pn);
    ///nn
    TH1F* hGenRcp_nn = new TH1F("hGenRcp_nn","hGenRcp_nn",nWCentrality,arrNpartBins);
    TH1F* hGenPlusRcp_nn = new TH1F("hGenPlusRcp_nn","hGenPlusRcp_nn",nWCentrality,arrNpartBins);
    TH1F* hGenMinusRcp_nn = new TH1F("hGenMinusRcp_nn","hGenMinusRcp_nn",nWCentrality,arrNpartBins);
    fillGeneratorRcp(fileNameMc_nn,centralityBins,nWCentrality,hGenRcp_nn,hGenPlusRcp_nn,hGenMinusRcp_nn);
    for (int j=0; j < nWCentrality; ++j) {

        double npart = Npart[j];
        double npartH = npart*nparterr[j] ; npartH*=0.01;
        double npartL = npart*nparterr[j] ; npartL*=0.01;

        ///Acceptance correction factor Aw in centrality bin j
        /*double aw = grAw->GetY()[j];
        double awPlus = grAwPlus->GetY()[j];
        double awMinus = grAwMinus->GetY()[j];
        */
        double aw = 1.0;
        double awPlus = 1.0;
        double awMinus = 1.0;

    	fillRcp(j,npart,npartH,npartL,yWSum[j],eyWSumL[j],syWSumHTotal[j],syWSumHTotal_correlated[j]*yWSum[j],graphWSum,graphWSumSyst,graphWSumSyst_correlated, aw);

        cout << centralityLow[j] << " - " << centralityHigh[j] << "\t" << yWSum[j]<< " +" 
            << eyWSumL[j] << " -" << eyWSumH[j] << "\t" << "+" 
		    << eyWSumSystL[j] << " (syst.) "<< " -" << eyWSumSystL[j] << " (syst.)" << endl;

        std::cout << "Rcp mu^{#pm}: Point " << j << " = " << yWSum[j] << " +- " << eyWSumL[j] << "(stat)" << 
            eyWSumSystL[j] << "(syst) " << erNcoll[j]*0.01 << "(ncoll) " << std::endl;
            std::cout << std::endl;
        
        writeRcpToSpreadsheet(spreadSheetRcpMuInclusive,j,yWSum[j],eyWSumL[j],eyWSumSystL[j],erNcoll[j]*0.01); 
	    
        if(doChargePlots){
		
            //mu+
		    fillRcp(j,npart,npartH,npartL,yWPlusSum[j],eyWPlusSumL[j],syWPlusSumTotal[j],
                    syWPlusSumTotal_correlated[j]*yWPlusSum[j],graphWPlusSum,graphWPlusSumSyst,graphWPlusSumSyst_correlated,awPlus);
        
            std::cout << "Rcp mu^{+}: Point " << j << " = " << yWPlusSum[j] << " +- " << eyWPlusSumL[j] << "(stat)" << 
                eyWPlusSumSystL[j] << "(syst) " << erNcoll[j]*0.01 << "(ncoll) " << std::endl;
                std::cout << std::endl;

            writeRcpToSpreadsheet(spreadSheetRcpMuPlus,j,yWPlusSum[j],eyWPlusSumL[j],eyWPlusSumSystL[j],erNcoll[j]*0.01); 

		    //mu-
		    fillRcp(j,npart,npartH,npartL,yWMinusSum[j],eyWMinusSumL[j],syWMinusSumTotal[j],
                    syWMinusSumTotal_correlated[j]*yWMinusSum[j],graphWMinusSum,graphWMinusSumSyst,graphWMinusSumSyst_correlated,awMinus);
            
            std::cout << "Rcp mu^{-}: Point " << j << " = " << yWMinusSum[j] << " +- " << eyWMinusSumL[j] << "(stat)" << 
                eyWMinusSumSystL[j] << "(syst) " << erNcoll[j]*0.01 << "(ncoll) " << std::endl;
                std::cout << std::endl;

            writeRcpToSpreadsheet(spreadSheetRcpMuMinus,j,yWMinusSum[j],eyWMinusSumL[j],eyWMinusSumSystL[j],erNcoll[j]*0.01); 
	}
   }///jth centrality

    TGraphAsymmErrors* graphWSumc = (TGraphAsymmErrors*)graphWSum->Clone("graphWSumc");
    TGraphAsymmErrors* graphWSumSystc = (TGraphAsymmErrors*)graphWSumSyst->Clone("graphWSumSystc");
    TGraphAsymmErrors* graphWSumSyst_correlatedc = (TGraphAsymmErrors*)graphWSumSyst_correlated->Clone("graphWSumSyst_correlatedc");

    TGraphAsymmErrors* graphWPlusSumc = (TGraphAsymmErrors*)graphWPlusSum->Clone("graphWPlusSumc");
    TGraphAsymmErrors* graphWPlusSumSystc = (TGraphAsymmErrors*)graphWPlusSumSyst->Clone("graphWPlusSumSystc");
    TGraphAsymmErrors* graphWPlusSumSyst_correlatedc = (TGraphAsymmErrors*)graphWPlusSumSyst_correlated->Clone("graphWPlusSumSyst_correlatedc");

    TGraphAsymmErrors* graphWMinusSumc = (TGraphAsymmErrors*)graphWMinusSum->Clone("graphWMinusSumc");
    TGraphAsymmErrors* graphWMinusSumSystc = (TGraphAsymmErrors*)graphWMinusSumSyst->Clone("graphWMinusSumSystc");
    TGraphAsymmErrors* graphWMinusSumSyst_correlatedc = (TGraphAsymmErrors*)graphWMinusSumSyst_correlated->Clone("graphWMinusSumSyst_correlatedc");

    TString sCh = "W^{#pm}#rightarrow#mu^{#pm}#nu";
    TString sChP = "W^{+}#rightarrow#mu^{+}#nu";
    TString sChM = "W^{-}#rightarrow#mu^{-}#bar{#nu}";

   
    if (doChargePlots) {

      cout << "doing charge plots" << endl;

      TCanvas* cRcp = new TCanvas("cRcp","cRcp",600,600);
      plotRcp(cRcp,graphWSumc,graphWSumSystc,graphWSumSyst_correlatedc,sCh,
                graphWPlusSumc,graphWPlusSumSystc,graphWPlusSumSyst_correlatedc,
                graphWMinusSumc,graphWMinusSumSystc,graphWMinusSumSyst_correlatedc,
                hGenRcpAllNuc,hGenRcp_pp,hGenRcp_nppn,hGenRcp_nn);
      
      Write(outFile,graphWSumc,"RcpChargeInclusive");
      Write(outFile,graphWSumSystc,"RcpChargeInclusive_uncorrelatedErrors");
      Write(outFile,graphWSumSyst_correlatedc,"RcpChargeInclusive_correlatedErrors");
      Write(outFile,graphWPlusSumc,"RcpMuPlus");
      Write(outFile,graphWPlusSumSystc,"RcpMuPlus_uncorrelatedErrors");
      Write(outFile,graphWPlusSumSyst_correlatedc,"RcpMuPlus_correlatedErrors");
      Write(outFile,graphWMinusSumc,"RcpMuMinus");
      Write(outFile,graphWMinusSumSystc,"RcpMuMinus_uncorrelatedErrors");
      Write(outFile,graphWMinusSumSyst_correlatedc,"RcpMuMinus_correlatedErrors");
      //Write(outFileTheory,(TF1*)cRcp->GetPrimitive("funcPythia"),"PYTHIA_LO_NPART");
      //plotRcp(graphWPlusSumc,graphWPlusSumSystc,graphWPlusSumSyst_correlatedc,sChP);
      //plotRcp(graphWMinusSumc,graphWMinusSumSystc,graphWMinusSumSyst_correlatedc,sChM);

      TGraphAsymmErrors* graphWRatio = new TGraphAsymmErrors(nWCentrality);
      TGraphAsymmErrors* graphWRatioSyst = new TGraphAsymmErrors(nWCentrality);
      TGraphAsymmErrors* graphWRatioSystCorrelated = new TGraphAsymmErrors(nWCentrality);
      ///ratio from pythia
      TGraph* grRatioPythia = new TGraph(1);
      //from POWHEG
      TGraph* grRatioPwg = new TGraph(1);
      ///ration from u,d quark counting
      TGraph* grRatioQuark = new TGraph(1);
      
      // Fill data ratios
      for (int j = 0; j<nWCentrality; ++j) {
	    
        double npart = Npart[j];
        double npartH = npart*nparterr[j] ; npartH*=0.01;
        double npartL = npart*nparterr[j] ; npartL*=0.01;

        //float normWPlus = 1./yWPlusSumRef*yFactorBw[j] ; 
        //float normWMinus = 1./yWMinusSumRef*yFactorBw[j] ; 
        float normW = yFactorBw[j] ; 
        float normWPlus = yFactorBw[j] ; 
        float normWMinus = yFactorBw[j] ; 

        double yWPlus = yWPlusSum[j]/normWPlus;
	    std::cout << "W+ in bin " << j << " = " << yWPlus << std::endl;
        double eyWPlusL = eyWPlusSumL[j]/normWPlus;
        double eyWPlusH = eyWPlusSumH[j]/normWPlus;
        double yWMinus = yWMinusSum[j]/normWMinus;
	    std::cout << "W- in bin " << j << " = " << yWMinus << std::endl;
        double eyWMinusL = eyWMinusSumL[j]/normWMinus;
        double eyWMinusH = eyWMinusSumH[j]/normWMinus;
        double yWRatio = yWPlus/yWMinus;
        std::cout << "W+/W- :" << yWRatio << std::endl; 
        double eyWRatioL = yWRatio * 0.01*sqrt( sqr(eyWPlusL/yWPlus*100) + sqr(eyWMinusH/yWMinus*100));
        double eyWRatioH = yWRatio * 0.01*sqrt( sqr(eyWPlusH/yWPlus*100) + sqr(eyWMinusL/yWMinus*100));
        ///stat+uncorrelated systematic error
        double eyTotWRatio = yWRatio * sqrt( sqr(eyWPlusH/yWPlus) + sqr(eyWMinusL/yWMinus) 
                                            + sqr(eyWPlusSumSyst10[j]/yWPlus)+sqr(eyWMinusSumSyst10[j]/yWMinus));
        std::cout << "stat.+uncorrelated syst. on ratio: " << eyTotWRatio << std::endl; 
        ///correlated systematic error

        double syWMinus_correlatedTemp = eyWMinusSumSystL[j]/yWMinusSum[j];
        double syWPlus_correlatedTemp = eyWPlusSumSystL[j]/yWPlusSum[j];
        double eySystCorrelatedRatio = yWRatio * sqrt(sqr(syWMinus_correlatedTemp)+sqr(syWPlus_correlatedTemp));
        //std::cout << "eySystCorrelatedRatio: " << eyWMinusSumSystH[j]/normWMinus << " " << eyWPlusSumSystH[j]/normWPlus << " " 
        //    << sqrt(sqr(eyWMinusSumSystH[j]/normWMinus/yWMinus)+sqr(eyWPlusSumSystH[j]/normWPlus/yWPlus)) << std::endl; exit(0); //hack

        graphWRatio->SetPoint(j, npart+13, yWRatio);
        graphWRatio->SetPointError(j, npartL+12., npartH+12., eyWRatioL, eyWRatioH );
        graphWRatioSyst->SetPoint(j, npart+13, yWRatio);
        graphWRatioSyst->SetPointError(j, npartL+12., npartH+12., eyTotWRatio, eyTotWRatio );
        graphWRatioSystCorrelated->SetPoint(j, npart, yWRatio);
        graphWRatioSystCorrelated->SetPointError(j, npartL+7.0, npartH+7.0, eySystCorrelatedRatio,eySystCorrelatedRatio);
      } //jth centrality

      plotChargeRatio(graphWRatio,graphWRatioSyst,graphWRatioSystCorrelated,grRatioPythia,grRatioPwg,grRatioQuark);
    
      Write(outFile,graphWRatio,"ChargeRatio");
      Write(outFile,graphWRatioSyst,"ChargeRatioUncorrelated");
      Write(outFile,graphWRatioSystCorrelated,"ChargeRatioCorrelated");
  }
}

   outFile->Close();
   spreadSheetMuPlus.close();
   spreadSheetRcpMuPlus.close();
   spreadSheetRcpMuMinus.close();
   spreadSheetMuMinus.close();
   spreadSheetRcpMuInclusive.close();
}
/// eof
