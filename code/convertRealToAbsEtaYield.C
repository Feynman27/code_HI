//////////////////////////////////////////////////////////////
//This macro take the yield in real eta (i.e. not |eta|)
//and combines them into |eta| bins. The output is written
//to a csv file.
///////////////////////////////////////////////////////////////
#include "TString.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TList.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>
#include <vector>


///////////////////////////////////////////
//writeToSpreadsheet
//////////////////////////////////////////
void writeToSpreadsheet(std::ostream& outputFile, int ieta, int icent, double sigPlus, double sigMinus){
	std::cout << "Writing yields to spreadsheet..." << std::endl;
	outputFile << ieta << "," << icent << "," << sigPlus << "," << sigMinus << std::endl;
}

///////////////////////////////////////////
//write TGraphs to root file
//////////////////////////////////////////
void Write(TFile* const outFile, TObject* const gr, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing TGraph to root file..." << std::endl;
	  gr->Write(sFile);
	  std::cout << "Done." << std::endl;
    	  gDirectory = dir;
}

int indexIJ(int i, int j, int nj){
  return i*nj+j ;
}
///////////////////////////////
//used for mapping negative and positive eta bins
//to absolute eta bin values
//////////////////////////////
int indexNegEta(int imap, int nEtaBins){
     int index = 0.5*(nEtaBins+1.0)-2.0; index-=imap;
     std::cout << "indexNegativeEta:" << index << std::endl;
      if(index>nEtaBins) std::cout << "WARNING: index out of allowed range. EXPECT WRONG RESULTS."<<std::endl;
     return index;
}

int indexPosEta(int imap, int nEtaBins){
      int index = imap+1.0/2.0*(nEtaBins+1.0);
      std::cout << "indexPositiveEta:" << index << std::endl;
      if(index>nEtaBins) std::cout << "WARNING: index out of allowed range. EXPECT WRONG RESULTS."<<std::endl;
      return index;
}



void convertRealToAbsEtaYield(){

  ///main file used for result
  TString fileNameDataIn ;
  fileNameDataIn =
        //"ResultFiles/WAnalysis_fitResultCentChrgEta_RealEtaBSCorr.05.15.2013"; 
        //"ResultFiles/WAnalysis_fitResultCentChrgEta_RealEtaBSCorr.05.25.2013"; 
        //"ResultFiles/WAnalysis_fitResultCentChrgEta_RealEtaBSCorr.06.08.2013"; 
        "systematics/WAnalysis_fitResultCentChrgEta_CwAsymmetry.07.14.2013";

  TFile* fDataSet = new TFile(fileNameDataIn+".root", "READ");
  if ( !fDataSet->IsOpen() ) {
    std::cout << fDataSet << " not found!" << std::endl;
    exit(0);
  }

  TString fileNameDataOut = "etaYieldDifferences";
  TFile *outFile = new TFile(fileNameDataOut+".root","RECREATE");
  //yield spreadsheet names	
  TString spreadSheetName = "dataSpreadSheetRealToAbsEta.csv";
  std::ofstream spreadSheet;
  spreadSheet.open(spreadSheetName);



  TGraphAsymmErrors* grWEta = (TGraphAsymmErrors*) fDataSet->Get("WPtFit_eta") ;
  const int nWEta = grWEta->GetN();
  std::cout << "Number of eta windows: " << nWEta << std::endl;

  TGraphAsymmErrors* grWCentrality = (TGraphAsymmErrors*) fDataSet->Get("WPtFit_centrality") ;
  const int nWCentrality =  grWCentrality->GetN() ;
  std::cout << "Number of centrality classes: " << nWCentrality << std::endl;

  ///total number of TGraphs in fit result
  const int nWGraphs = nWCentrality*nWEta ;
  TObjArray arrWSig = TObjArray(nWGraphs) ;
  ///TGraph names from yield file
  TString baseW = "WPtFit" ;
  TString baseWSig = baseW ; baseWSig += "_sig_eta" ;

  double* yWPlus = new double[nWEta*nWCentrality] ;
  double* yWMinus = new double[nWEta*nWCentrality] ;
  double* yWStatPlus = new double[nWEta*nWCentrality] ;
  double* yWStatMinus = new double[nWEta*nWCentrality] ;

  for (int i=0; i<nWEta; i++){
      TString searchWI = baseWSig;  searchWI += i;  searchWI += "_centrality";
      for(int j=0; j<nWCentrality; ++j){
        TString searchWIJ = searchWI;  searchWIJ += j;
        ///Retrieve TGraph of yields per bin
        const int index = indexIJ(i,j,nWCentrality) ;
        arrWSig[index] = ((TGraphAsymmErrors*) fDataSet->Get(searchWIJ)) ;

        ///Pointer to W counts in bin i,j
        double* yW0;
        double* yW0Stat ;

        ///Assign to the W yield in bin i,j
        yW0 = ((TGraphAsymmErrors*)arrWSig[index])->GetY();
        yW0Stat = ((TGraphAsymmErrors*)arrWSig[index])->GetEYlow();
        
        int _nPoints = ((TGraphAsymmErrors*) arrWSig[index])->GetN() ;
        for(int iPoint=0; iPoint<_nPoints; ++iPoint){

            if(yW0[iPoint]>1){
                
	            ///mu+
	            if(iPoint == 102){
                    yWPlus[index] = yW0[iPoint];
                    yWStatPlus[index] = yW0Stat[iPoint];
                    std::cout << "Found NW+ = " << yW0[iPoint] << "+-" << yW0Stat[iPoint] << std::endl;

                }
                ///mu-
	            if(iPoint == 103){
                    yWMinus[index] = yW0[iPoint];
                    yWStatMinus[index] = yW0Stat[iPoint];
                    std::cout << "Found NW- = " << yW0[iPoint] << "+-" << yW0Stat[iPoint] << std::endl;
                }

            }
        }//iPoint
     } //j(cent)
    }//i (eta)

  std::vector<double> etaBins;
  etaBins.push_back(0.10);
  etaBins.push_back(0.35);
  etaBins.push_back(0.6);
  etaBins.push_back(0.8);
  etaBins.push_back(1.05);
  etaBins.push_back(1.3);
  etaBins.push_back(1.55);
  etaBins.push_back(1.85);
  etaBins.push_back(2.1);
  etaBins.push_back(+2.40);

  const int nAbsEtaBins = etaBins.size()-1;

  double* yWPlusAbsEtaSum = new double[nAbsEtaBins*nWCentrality];
  double yWPlusNegEtaSummedOverCent [nAbsEtaBins] = {0.0};
  double yWPlusNegEtaSummedOverCentStat  [nAbsEtaBins] = {0.0};
  double yWPlusPosEtaSummedOverCent [nAbsEtaBins] = {0.0};
  double yWPlusPosEtaSummedOverCentStat [nAbsEtaBins] = {0.0};

  double* yWMinusAbsEtaSum = new double[nAbsEtaBins*nWCentrality];
  double yWMinusNegEtaSummedOverCent [nAbsEtaBins] = {0.0};
  double yWMinusNegEtaSummedOverCentStat  [nAbsEtaBins] = {0.0};
  double yWMinusPosEtaSummedOverCent [nAbsEtaBins] = {0.0};
  double yWMinusPosEtaSummedOverCentStat [nAbsEtaBins] = {0.0};

  TList _grPlus, _grMinus ;
  TGraphErrors* _grPlusAllCent = new TGraphErrors(nAbsEtaBins);
  TGraphErrors* _grMinusAllCent = new TGraphErrors(nAbsEtaBins);

  for(int j=0; j<nWCentrality; ++j){
     _grPlus.Add( new TGraphErrors(nAbsEtaBins) );
     _grMinus.Add( new TGraphErrors(nAbsEtaBins) );
  }

  for (int i=0; i<nAbsEtaBins; i++){

      double xEta = etaBins[i]+(etaBins[i+1]-etaBins[i])/2.0;
      ///Return bin-pair corresponding to 
      ///the ith |eta| bin
      int binNeg = indexNegEta(i,nWEta);
      int binPos = indexPosEta(i,nWEta);
      for(int j=0; j<nWCentrality; ++j){

        int index = indexIJ(i,j,nWCentrality);
        ///get the number of signal candidates in 
        ///each +-eta slice
        int indexNeg = indexIJ(binNeg,j,nWCentrality);
        double yW0PlusNegEta = yWPlus[indexNeg];
        double yW0PlusNegEtaStat = yWStatPlus[indexNeg];
        int indexPos = indexIJ(binPos,j,nWCentrality);
        double yW0PlusPosEta = yWPlus[indexPos];
        double yW0PlusPosEtaStat = yWStatPlus[indexPos];

        ///add them together to get the absolute eta yield
        yWPlusAbsEtaSum[index] = yW0PlusNegEta+yW0PlusPosEta;
        std::cout << yW0PlusNegEta << "+-" << yW0PlusNegEtaStat << " " << yW0PlusPosEta << "+-" << yW0PlusPosEtaStat << std::endl;
        double statErrPlusTemp = TMath::Sqrt(TMath::Power(yW0PlusNegEtaStat,2)+TMath::Power(yW0PlusPosEtaStat,2));
        std::cout << fabs(yW0PlusNegEta-yW0PlusPosEta) << "+-" << statErrPlusTemp << std::endl;
        ///Plot the differences 
        double etaBinDiffPlus = fabs(yW0PlusNegEta-yW0PlusPosEta);

        ///Sum over centrality for each eta bin
        yWPlusNegEtaSummedOverCent[i] += yW0PlusNegEta;
        std::cout << "Running Sum of mu+ negative eta = " << yWPlusNegEtaSummedOverCent[i] << std::endl;
        yWPlusNegEtaSummedOverCentStat[i] += TMath::Power(yW0PlusNegEtaStat,2);
        yWPlusPosEtaSummedOverCent[i] += yW0PlusPosEta;
        std::cout << "Running Sum of mu+ positive eta = " << yWPlusPosEtaSummedOverCent[i] << std::endl;
        yWPlusPosEtaSummedOverCentStat[i] += TMath::Power(yW0PlusPosEtaStat,2);

        std::cout << "Mu+ : " << std::endl;
        std::cout << "Difference at eta " << xEta << " = " << yW0PlusNegEta << "-" << yW0PlusPosEta << " = " << 
            etaBinDiffPlus << " +- " << statErrPlusTemp << std::endl;
        ( (TGraphErrors*)_grPlus.At(j) )->SetPoint(i,xEta,etaBinDiffPlus);
        ( (TGraphErrors*)_grPlus.At(j) )->SetPointError(i,(etaBins[i+1]-etaBins[i])/2.0,statErrPlusTemp);

        ///ditto for W-
        double yW0MinusNegEta = yWMinus[indexNeg];
        double yW0MinusNegEtaStat = yWStatMinus[indexNeg];
        double yW0MinusPosEta = yWMinus[indexPos];
        double yW0MinusPosEtaStat = yWStatMinus[indexPos];

        yWMinusAbsEtaSum[index] = yW0MinusNegEta+yW0MinusPosEta;
        double statErrMinusTemp = TMath::Sqrt(TMath::Power(yW0MinusNegEtaStat,2)+TMath::Power(yW0MinusPosEtaStat,2));

        double etaBinDiffMinus = fabs(yW0MinusNegEta-yW0MinusPosEta);

        ///Sum over centrality for each eta bin
        yWMinusNegEtaSummedOverCent[i] += yW0MinusNegEta;
        std::cout << "Running Sum of mu- negative eta = " << yWMinusNegEtaSummedOverCent[i] << std::endl;
        yWMinusNegEtaSummedOverCentStat[i] += TMath::Power(yW0MinusNegEtaStat,2);
        yWMinusPosEtaSummedOverCent[i] += yW0MinusPosEta;
        std::cout << "Running Sum of mu- positive eta = " << yWMinusPosEtaSummedOverCent[i] << std::endl;
        yWMinusPosEtaSummedOverCentStat[i] += TMath::Power(yW0MinusPosEtaStat,2);

        std::cout << "Mu- : " << std::endl;
        std::cout << "Difference at eta " << xEta << " = " << yW0MinusNegEta << "-" << yW0MinusPosEta << " = " << 
            etaBinDiffMinus << " +- " << statErrMinusTemp << std::endl;
        ( (TGraphErrors*)_grMinus.At(j) )->SetPoint(i,xEta,etaBinDiffMinus);
        ( (TGraphErrors*)_grMinus.At(j) )->SetPointError(i,(etaBins[i+1]-etaBins[i])/2.0,statErrMinusTemp);

        ///Save to a csv file

	    writeToSpreadsheet(spreadSheet,i,j, yWPlusAbsEtaSum[index],yWMinusAbsEtaSum[index]);
 
      } //j(cent)

      double etaBinDiffPlus = fabs(yWPlusNegEtaSummedOverCent[i]-yWPlusPosEtaSummedOverCent[i]); 
      double statErrPlus = TMath::Sqrt(yWPlusNegEtaSummedOverCentStat[i]+yWPlusPosEtaSummedOverCentStat[i]);
      _grPlusAllCent->SetPoint(i,xEta,etaBinDiffPlus);
      _grPlusAllCent->SetPointError(i,(etaBins[i+1]-etaBins[i])/2.0,statErrPlus);
      double etaBinDiffMinus = fabs(yWMinusNegEtaSummedOverCent[i]-yWMinusPosEtaSummedOverCent[i]); 
      double statErrMinus = TMath::Sqrt(yWMinusNegEtaSummedOverCentStat[i]+yWMinusPosEtaSummedOverCentStat[i]);
      _grMinusAllCent->SetPoint(i,xEta,etaBinDiffMinus);
      _grMinusAllCent->SetPointError(i,(etaBins[i+1]-etaBins[i])/2.0,statErrMinus);

  } //i (absEtaBins)

  //Save the TGraphs of eta bin yield differences
  Write(outFile,_grPlusAllCent,"etaYieldDiffAllCentCharge102");
  Write(outFile,_grMinusAllCent,"etaYieldDiffAllCentCharge103");

  for(int i=0;i<_grPlus.GetEntries();++i){
        TString sGrName = "etaYieldDiffCent"; sGrName+=i; sGrName+="Charge";
        Write(outFile,_grPlus.At(i),sGrName+"102");
        Write(outFile,_grMinus.At(i),sGrName+"103");
  }

  std::cout << "Closing spreadsheet." << std::endl;
  spreadSheet.close();
  std::cout << "Done running macro." << std::endl;

  std::cout << "Clean up." << std::endl;

  delete _grPlusAllCent;
  delete _grMinusAllCent;

  for(int i=0;i<_grPlus.GetEntries();++i){
    delete _grPlus.At(i);
    delete _grMinus.At(i);
  }

  /*delete[] yWPlusAbsEtaSum;
  delete[] yWPlusNegEtaSummedOverCent;
  delete[] yWPlusNegEtaSummedOverCentStat;
  delete[] yWPlusPosEtaSummedOverCent;
  delete[] yWPlusPosEtaSummedOverCentStat;

  delete[] yWMinusAbsEtaSum;
  delete[] yWMinusNegEtaSummedOverCent;
  delete[] yWMinusNegEtaSummedOverCentStat;
  delete[] yWMinusPosEtaSummedOverCent;
  delete[] yWMinusPosEtaSummedOverCentStat;
  */
}
