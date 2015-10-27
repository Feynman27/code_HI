//////////////////////////////////////////////////////////////
//This macro take the Cw in real eta (i.e. not |eta|)
//and plots the difference with respect to absolute eta Cw. The output is written
//to a file full of TGraphs.
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
///////////////////////////////
//used for mapping negative and positive eta bins
//to absolute eta bin values
//////////////////////////////
int indexPositiveToAbsEta(int ieta, int nEtaBins){
      int index = ieta-1.0/2.0*(nEtaBins+1.0);
      std::cout << "indexPositiveEta:" << index << std::endl;
      if(index>(nEtaBins-1)/2) std::cout << "WARNING: index out of allowed range. EXPECT WRONG RESULTS."<<std::endl;
      return index;
}
int indexNegativeToAbsEta(int ieta, int nEtaBins){
     int index = 0.5*(nEtaBins+1.0)-2.0; index-=ieta;
     std::cout << "indexNegativeEta:" << index << std::endl;
     if(index>(nEtaBins-1)/2) std::cout << "WARNING: index out of allowed range. EXPECT WRONG RESULTS."<<std::endl;
     return index;
}
void plotCwDifferences(){


   TString fileNameDataOut = "cWDifferencesEtaDistros";
   TFile *outFile = new TFile(fileNameDataOut+".root","RECREATE");

   ///Open files containing Cw values calculated in absolute eta
   TString fileNameCwAbsEta = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges.05.25.2013";
   std::cout << "Filename: " << fileNameCwAbsEta << std::endl;
   TFile* fCwAbsEta = new TFile(fileNameCwAbsEta+".root","READ");

   ///Open file containing Cw values in neg and pos eta
   TString fileNameCwAllEta = "CorrectionFactorFiles/correctionFactorsWEtaCent_19EtaBinsNoAbsEta6CentBins2Charges.05.26.2013";
   std::cout << "Filename: " << fileNameCwAllEta << std::endl;
   TFile* fCwAllEta = new TFile(fileNameCwAllEta+".root","READ");

   std::cout << "All files open." << std::endl;
  
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
   
   const int nWCentrality = 6;
   TList _grPlusNegPosDiff,_grMinusNegPosDiff,_grPlusAbsEtaNegDiff,_grPlusAbsEtaPosDiff,_grMinusAbsEtaNegDiff,_grMinusAbsEtaPosDiff;
   TList _grCwAbsEtaMirrorMinus,_grCwAbsEtaMirrorPlus,_grCwAllEtaMinusc,_grCwAllEtaPlusc;

   for(int j=0; j<nWCentrality; ++j){
     _grPlusNegPosDiff.Add( new TGraphErrors(nAbsEtaBins) );
     _grMinusNegPosDiff.Add( new TGraphErrors(nAbsEtaBins) );
     _grPlusAbsEtaNegDiff.Add( new TGraphErrors(nAbsEtaBins) );
     _grPlusAbsEtaPosDiff.Add( new TGraphErrors(nAbsEtaBins) );
     _grMinusAbsEtaNegDiff.Add( new TGraphErrors(nAbsEtaBins) );
     _grMinusAbsEtaPosDiff.Add( new TGraphErrors(nAbsEtaBins) );
     _grCwAbsEtaMirrorMinus.Add( new TGraphAsymmErrors(19) );
     _grCwAbsEtaMirrorPlus.Add( new TGraphAsymmErrors(19) );
     //_grCwAllEtaMinusc.Add( new TGraphAsymmErrors(19) );
     //_grCwAllEtaPlusc.Add( new TGraphAsymmErrors(19) );
   }

       for (int ieta=0; ieta<nAbsEtaBins; ieta++){
               
	        for (int icent=0; icent<nWCentrality; icent++){

               ///Fetch the TGraphs from their respective TFiles
               TString sCwDistroBase = "grWmunuRecHiQualityWselEtaDistro"; 
               TString sCwDistroPlus = sCwDistroBase+"PlusCent"; sCwDistroPlus+=icent;
               TString sCwDistroMinus = sCwDistroBase+"MinusCent"; sCwDistroMinus+=icent;

               std::cout << "Fetching eta Cw distro for mu+ named " << sCwDistroPlus << std::endl;
               TGraphErrors* grCwAllEtaPlus = (TGraphErrors*)fCwAllEta->Get(sCwDistroPlus);
               TGraphErrors* grCwAbsEtaPlus = (TGraphErrors*)fCwAbsEta->Get(sCwDistroPlus);

               std::cout << "Fetching eta Cw distro for mu- named " << sCwDistroMinus << std::endl;
               TGraphErrors* grCwAllEtaMinus = (TGraphErrors*)fCwAllEta->Get(sCwDistroMinus);
               TGraphErrors* grCwAbsEtaMinus = (TGraphErrors*)fCwAbsEta->Get(sCwDistroMinus);

               int nWEta = grCwAllEtaMinus->GetN();
               std::cout << "Number of eta windows over all eta: " << nWEta << std::endl;

	           double etaX = etaBins[ieta]+(etaBins[ieta+1]-etaBins[ieta])/2.0;
               ///Return bin-pair(neg,pos eta) corresponding to 
               ///the ith |eta| bin (e.g. [-2.4,-2.1] and [2.1,2.4])

               int binNeg = indexNegEta(ieta,nWEta);
               int binPos = indexPosEta(ieta,nWEta);
               
               ///mu+
               double CwPlusEtaNeg = grCwAllEtaPlus->GetY()[binNeg];
               double CwPlusEtaPos = grCwAllEtaPlus->GetY()[binPos];

               ///Place the difference in a TGraph
               double CwPlusNegPosDiff = fabs(CwPlusEtaNeg-CwPlusEtaPos); 
               ( (TGraphErrors*)_grPlusNegPosDiff.At(icent) )->SetPoint(ieta,etaX,CwPlusNegPosDiff);
               ( (TGraphErrors*)_grPlusNegPosDiff.At(icent) )->SetPointError(ieta,(etaBins[ieta+1]-etaBins[ieta])/2.0,0.0);

               ///mu-
               double CwMinusEtaNeg = grCwAllEtaMinus->GetY()[binNeg];
               double CwMinusEtaPos = grCwAllEtaMinus->GetY()[binPos];

               ///Place the difference in a TGraph
               double CwMinusNegPosDiff = fabs(CwMinusEtaNeg-CwMinusEtaPos); 
               ( (TGraphErrors*)_grMinusNegPosDiff.At(icent) )->SetPoint(ieta,etaX,CwMinusNegPosDiff);
               ( (TGraphErrors*)_grMinusNegPosDiff.At(icent) )->SetPointError(ieta,(etaBins[ieta+1]-etaBins[ieta])/2.0,0.0);

               ///Calculate differences wrt to absolute bin
               double CwPlusAbsEta = grCwAbsEtaPlus->GetY()[ieta];

               ///Take the difference of neg and pos eta wrt abs eta bin
               double CwPlusAbsEtaNegDiff = fabs(CwPlusAbsEta-CwPlusEtaNeg);
               double CwPlusAbsEtaPosDiff = fabs(CwPlusAbsEta-CwPlusEtaPos);

               ///Place the difference in a TGraph
               ( (TGraphErrors*)_grPlusAbsEtaNegDiff.At(icent) )->SetPoint(ieta,etaX,CwPlusAbsEtaNegDiff);
               ( (TGraphErrors*)_grPlusAbsEtaNegDiff.At(icent) )->SetPointError(ieta,(etaBins[ieta+1]-etaBins[ieta])/2.0,0.0);

               ( (TGraphErrors*)_grPlusAbsEtaPosDiff.At(icent) )->SetPoint(ieta,etaX,CwPlusAbsEtaPosDiff);
               ( (TGraphErrors*)_grPlusAbsEtaPosDiff.At(icent) )->SetPointError(ieta,(etaBins[ieta+1]-etaBins[ieta])/2.0,0.0);

               ///ditto for mu-
               double CwMinusAbsEta = grCwAbsEtaMinus->GetY()[ieta];
               ///Take the difference of neg and pos eta wrt abs eta bin
               double CwMinusAbsEtaNegDiff = fabs(CwMinusAbsEta-CwMinusEtaNeg);
               double CwMinusAbsEtaPosDiff = fabs(CwMinusAbsEta-CwMinusEtaPos);

               
               std::cout << "I N P U T  P A R A M E T E R S" << std::endl;
               std::cout << "CwPlusAbsEta    = " << CwPlusAbsEta << std::endl;
               std::cout << "CwPlusEtaNeg    = " << CwPlusEtaNeg << std::endl;
               std::cout << "CwPlusEtaPos    = " << CwPlusEtaPos << std::endl;
               std::cout << "CwMinusAbsEta    = " << CwMinusAbsEta << std::endl;
               std::cout << "CwMinusEtaNeg    = " << CwMinusEtaNeg << std::endl;
               std::cout << "CwMinusEtaPos    = " << CwMinusEtaPos << std::endl;
               std::cout << std::endl;


               ///Place the difference in a TGraph
               ( (TGraphErrors*)_grMinusAbsEtaNegDiff.At(icent) )->SetPoint(ieta,etaX,CwMinusAbsEtaNegDiff);
               ( (TGraphErrors*)_grMinusAbsEtaNegDiff.At(icent) )->SetPointError(ieta,(etaBins[ieta+1]-etaBins[ieta])/2.0,0.0);

               ( (TGraphErrors*)_grMinusAbsEtaPosDiff.At(icent) )->SetPoint(ieta,etaX,CwMinusAbsEtaPosDiff);
               ( (TGraphErrors*)_grMinusAbsEtaPosDiff.At(icent) )->SetPointError(ieta,(etaBins[ieta+1]-etaBins[ieta])/2.0,0.0);

            }
       }

       int nWEta = 19;
       std::cout << "Number of eta windows over all eta: " << nWEta << std::endl;
  
       ///Mirror the |eta| distro and plot over all eta
       for (int ieta=0; ieta<nWEta; ieta++){
               
	        for (int icent=0; icent<nWCentrality; icent++){

               ///Fetch the TGraphs from their respective TFiles
               TString sCwDistroBase = "grWmunuRecHiQualityWselEtaDistro"; 
               TString sCwDistroPlus = sCwDistroBase+"PlusCent"; sCwDistroPlus+=icent;
               TString sCwDistroMinus = sCwDistroBase+"MinusCent"; sCwDistroMinus+=icent;

               std::cout << "Fetching eta Cw distro for mu+ named " << sCwDistroPlus << std::endl;
               //TGraphAsymmErrors* grCwAllEtaPlusc = (TGraphAsymmErrors*)fCwAllEta->Get(sCwDistroPlus);
               _grCwAllEtaPlusc.Add((TGraphAsymmErrors*)fCwAllEta->Get(sCwDistroPlus)); 

               TGraphAsymmErrors* grCwAbsEtaPlusc = (TGraphAsymmErrors*)fCwAbsEta->Get(sCwDistroPlus);

               std::cout << "Fetching eta Cw distro for mu- named " << sCwDistroMinus << std::endl;
               //TGraphAsymmErrors* grCwAllEtaMinusc = (TGraphErrors*)fCwAllEta->Get(sCwDistroMinus);
               _grCwAllEtaMinusc.Add((TGraphAsymmErrors*)fCwAllEta->Get(sCwDistroMinus)); 

               TGraphAsymmErrors* grCwAbsEtaMinusc = (TGraphAsymmErrors*)fCwAbsEta->Get(sCwDistroMinus);

	           double etaX = ( (TGraphAsymmErrors*)_grCwAllEtaMinusc.At(icent) )->GetX()[ieta];
               std::cout << "Plotting at eta " << etaX << std::endl;

               int indexAbs;
               if(etaX<0.0) indexAbs = indexNegativeToAbsEta(ieta,nWEta);
               if(etaX>0.0) indexAbs = indexPositiveToAbsEta(ieta,nWEta);
               if(etaX==0.0) continue;

               double yCwAbsPlus = grCwAbsEtaPlusc->GetY()[indexAbs];
               double yCwAbsPlusErrLo = grCwAbsEtaPlusc->GetEYlow()[indexAbs];
               double yCwAbsPlusErrHi = grCwAbsEtaPlusc->GetEYhigh()[indexAbs];
               ( (TGraphAsymmErrors*)_grCwAbsEtaMirrorPlus.At(icent) )->SetPoint(ieta,etaX,yCwAbsPlus);
               ( (TGraphAsymmErrors*)_grCwAbsEtaMirrorPlus.At(icent) )->SetPointError(ieta,0.0,0.0,yCwAbsPlusErrLo,yCwAbsPlusErrHi);
               double yCwAbsMinus = grCwAbsEtaMinusc->GetY()[indexAbs];
               double yCwAbsMinusErrLo = grCwAbsEtaMinusc->GetEYlow()[indexAbs];
               double yCwAbsMinusErrHi = grCwAbsEtaMinusc->GetEYhigh()[indexAbs];
               ( (TGraphAsymmErrors*)_grCwAbsEtaMirrorMinus.At(icent) )->SetPoint(ieta,etaX,yCwAbsMinus);
               ( (TGraphAsymmErrors*)_grCwAbsEtaMirrorMinus.At(icent) )->SetPointError(ieta,0.0,0.0,yCwAbsMinusErrLo,yCwAbsMinusErrHi);

             }
       }

   for(int i=0;i<_grPlusNegPosDiff.GetEntries();++i){
        TString sGrName = "etaCwNegPosDiffCent"; sGrName+=i; sGrName+="Charge";
        Write(outFile,_grPlusNegPosDiff.At(i),sGrName+"102");
        Write(outFile,_grMinusNegPosDiff.At(i),sGrName+"103");

        sGrName = "etaCwAbsEtaNegDiffCent"; sGrName+=i; sGrName+="Charge";
        Write(outFile,_grPlusAbsEtaNegDiff.At(i),sGrName+"102");
        Write(outFile,_grMinusAbsEtaNegDiff.At(i),sGrName+"103");

        sGrName = "etaCwAbsEtaPosDiffCent"; sGrName+=i; sGrName+="Charge";
        Write(outFile,_grPlusAbsEtaPosDiff.At(i),sGrName+"102");
        Write(outFile,_grMinusAbsEtaPosDiff.At(i),sGrName+"103");

        TString sCwDistroBase = "grWmunuRecHiQualityWselAbsEtaDistro"; 
        TString sCwDistroPlus = sCwDistroBase+"PlusCent"; sCwDistroPlus+=i;
        TString sCwDistroMinus = sCwDistroBase+"MinusCent"; sCwDistroMinus+=i;
        Write(outFile,_grCwAbsEtaMirrorMinus.At(i),sCwDistroMinus);
        Write(outFile,_grCwAbsEtaMirrorPlus.At(i),sCwDistroPlus);

        sCwDistroBase = "grWmunuRecHiQualityWselEtaDistro"; 
        sCwDistroPlus = sCwDistroBase+"PlusCent"; sCwDistroPlus+=i;
        sCwDistroMinus = sCwDistroBase+"MinusCent"; sCwDistroMinus+=i;
        Write(outFile,_grCwAllEtaMinusc.At(i),sCwDistroMinus);
        Write(outFile,_grCwAllEtaPlusc.At(i),sCwDistroPlus);
   }

   std::cout << "Clean Up." << std::endl;
   for(int j=0; j<_grPlusNegPosDiff.GetEntries(); ++j){
     _grPlusNegPosDiff.At(j);
     _grMinusNegPosDiff.At(j);
     _grPlusAbsEtaNegDiff.At(j);
     _grPlusAbsEtaPosDiff.At(j);
     _grMinusAbsEtaNegDiff.At(j);
     _grMinusAbsEtaPosDiff.At(j);
     _grCwAbsEtaMirrorMinus.At(j);
     _grCwAbsEtaMirrorPlus.At(j);
     _grCwAllEtaMinusc.At(j);
     _grCwAllEtaPlusc.At(j);
   }
}
