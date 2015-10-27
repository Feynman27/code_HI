#include "TString.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TList.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TLatex.h"
#include "TTree.h"
#include "TClonesArray.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath> 
#include <sstream>
#include <map>

/*#ifdef __CINT__
#pragma link C++ class map<int, int >+;
#pragma link C++ class map<int, double >+;
#endif
*/
using namespace std;

///////////////////////////////////////////
//writeToSpreadsheet
//////////////////////////////////////////
void writeToSpreadsheet(std::ostream& outputFile, int index, int run, double lumi, int muons){
	outputFile << index << "," << run << "," << lumi << "," << muons << std::endl;
}


void countMuonsPerRunD3PD(){

  //gInterpreter->GenerateDictionary("map<int, double >","map.h;map");
  int cutValue = 11;
  float ptmin = 10.0;
  //float ptmin = 4.;
  //TString fileNameMuonIn = "/afs/cern.ch/user/t/tbalestr/public/16.6.2.5/WorkArea/run/HISingleMuon.root";

///7.3 ub^-1 2010 data
  //TString fileNameMuonIn = "/tmp/tbalestr/HISingleMuon_Aug4_v2.root";
  ///2011 PbPb data
  //TString fileNameMuonIn = "/tmp/tbalestr/HISingleMuon_DPD_HP_merged_26May.root";
  //TString fileNameMuonIn = "/tmp/tbalestr/HISingleMuon_DPD_HP_30May.root";

  TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
  TString fileNameMuonIn = baseString+"HardProbesFiles/HISingleMuonHardProbesData.04.17.2013.root";  
  TFile* fMuIn = new TFile(fileNameMuonIn, "READ");
  if ( !fMuIn->IsOpen() ) {
    std::cout << fMuIn << " not found!" << std::endl;
  }

  TFile* outFile = new TFile("muonsPerRun.root","recreate");
  
  ///lumi in ub^-1
  std::vector<double> vLumi;
  vLumi.push_back(0.0137422);
  vLumi.push_back(0.126333);
  vLumi.push_back(1.28605451);
  vLumi.push_back(0.08307752);
  vLumi.push_back(3.948807);
  vLumi.push_back(2.766361);
  vLumi.push_back(3.606919);
  vLumi.push_back(3.809249);
  vLumi.push_back(4.11709);
  vLumi.push_back(4.88581);
  vLumi.push_back(2.50513);
  vLumi.push_back(0.0619227);
  vLumi.push_back(0.9694949);
  vLumi.push_back(0.582203);
  vLumi.push_back(5.42834);
  vLumi.push_back(4.20247);
  vLumi.push_back(0.38549);
  vLumi.push_back(5.50425);
  vLumi.push_back(3.15969);
  vLumi.push_back(4.19951);
  vLumi.push_back(5.24396);
  vLumi.push_back(5.30817);
  vLumi.push_back(0.0433861);
  vLumi.push_back(5.29888);
  vLumi.push_back(5.30346);
  vLumi.push_back(0.757212);
  vLumi.push_back(5.985104);
  vLumi.push_back(5.1515);
  vLumi.push_back(6.57761);
  vLumi.push_back(1.75062);
  vLumi.push_back(5.43771);
  vLumi.push_back(3.13928);
  vLumi.push_back(5.14981);
  vLumi.push_back(4.90751);
  vLumi.push_back(2.521785);
  vLumi.push_back(6.45269);
  vLumi.push_back(3.59951);
  vLumi.push_back(5.33795);
  vLumi.push_back(6.16766);
  vLumi.push_back(1.77379);
  vLumi.push_back(0.923454);
  vLumi.push_back(2.07656);

  std::vector<int> vRun;
  vRun.push_back(193211);
  vRun.push_back(193270);
  vRun.push_back(193291);
  vRun.push_back(193295);
  vRun.push_back(193321);
  vRun.push_back(193403);
  vRun.push_back(193412);
  vRun.push_back(193447);
  vRun.push_back(193463);
  vRun.push_back(193481);
  vRun.push_back(193491);
  vRun.push_back(193492);
  vRun.push_back(193493);
  vRun.push_back(193494);
  vRun.push_back(193546);
  vRun.push_back(193558);
  vRun.push_back(193599);
  vRun.push_back(193604);
  vRun.push_back(193641);
  vRun.push_back(193655);
  vRun.push_back(193662);
  vRun.push_back(193679);
  vRun.push_back(193687);
  vRun.push_back(193718);
  vRun.push_back(193795);
  vRun.push_back(193823);
  vRun.push_back(193825);
  vRun.push_back(193826);
  vRun.push_back(193834);
  vRun.push_back(193890);
  vRun.push_back(194017);
  vRun.push_back(194060);
  vRun.push_back(194061);
  vRun.push_back(194121);
  vRun.push_back(194160);
  vRun.push_back(194163);
  vRun.push_back(194179);
  vRun.push_back(194192);
  vRun.push_back(194193);
  vRun.push_back(194370);
  vRun.push_back(194374);
  vRun.push_back(194382);
  const int nRuns = vRun.size(); 

  TGraphErrors* grLumi = new TGraphErrors(nRuns);
  TGraphErrors* grMuPerLumi = new TGraphErrors(nRuns);

  int valNt[30], EF_mu4_MSonly_L1TE50_Matched20[30],EF_mu4_L1VTE50_Matched20[30],runNt, mu_muid_nNt;
  int EF_mu10_MSonly_EFFS_L1ZDC_Matched20[30], EF_mu10_MSonly_EFFS_L1TE10_Matched20[30],EF_mu10_MSonly_EFFS_L1TE20_Matched20[30];
  int EF_mu10_MSonly_EFFS_L1ZDC, EF_mu10_MSonly_EFFS_L1TE10, EF_mu10_MSonly_EFFS_L1TE20, EF_mu4_MSonly_L1TE50, EF_mu4_L1VTE50;
  float ptNt[30], centralityNt,etaNt[30],eLossNt[30],scatNt[30];

  TChain* tree = new TChain("tree","tree");
  std::cout << "Filling the tree for " << fileNameMuonIn << std::endl;
  tree->Add(fileNameMuonIn);
  // --- Set branch adresses ---
  tree->SetBranchAddress("val", &valNt);
  tree->SetBranchAddress("pt", &ptNt);
  tree->SetBranchAddress("eta", &etaNt);
  tree->SetBranchAddress("eLoss",&eLossNt);
  tree->SetBranchAddress("scat",&scatNt);
  tree->SetBranchAddress("mu_muid_n", &mu_muid_nNt);
  tree->SetBranchAddress("run", &runNt);
  tree->SetBranchAddress("centrality", &centralityNt); 
  tree->SetBranchAddress("EF_mu4_MSonly_L1TE50", &EF_mu4_MSonly_L1TE50); 
  tree->SetBranchAddress("EF_mu4_L1VTE50", &EF_mu4_L1VTE50); 
  tree->SetBranchAddress("EF_mu4_MSonly_L1TE50_Matched20", &EF_mu4_MSonly_L1TE50_Matched20); 
  tree->SetBranchAddress("EF_mu4_L1VTE50_Matched20", &EF_mu4_L1VTE50_Matched20); 
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC", &EF_mu10_MSonly_EFFS_L1ZDC); 
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10", &EF_mu10_MSonly_EFFS_L1TE10); 
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20", &EF_mu10_MSonly_EFFS_L1TE20); 
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC_Matched20", &EF_mu10_MSonly_EFFS_L1ZDC_Matched20); 
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10_Matched20", &EF_mu10_MSonly_EFFS_L1TE10_Matched20); 
  tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20_Matched20", &EF_mu10_MSonly_EFFS_L1TE20_Matched20); 

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("val", 1);
  tree->SetBranchStatus("pt", 1);
  tree->SetBranchStatus("eta", 1);
  tree->SetBranchStatus("eLoss",1);
  tree->SetBranchStatus("scat",1);
  tree->SetBranchStatus("mu_muid_n", 1);
  tree->SetBranchStatus("run", 1);
  tree->SetBranchStatus("centrality", 1); 
  tree->SetBranchStatus("EF_mu4_MSonly_L1TE50", 1); 
  tree->SetBranchStatus("EF_mu4_L1VTE50", 1); 
  tree->SetBranchStatus("EF_mu4_MSonly_L1TE50_Matched20", 1); 
  tree->SetBranchStatus("EF_mu4_L1VTE50_Matched20", 1); 
  tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC", 1); 
  tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10", 1); 
  tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20", 1); 
  tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC_Matched20", 1); 
  tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10_Matched20", 1); 
  tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20_Matched20", 1); 

 std::cout << "Number of entries: " << tree->GetEntries() << std::endl; 
 
 ///map runnumber to #mu
 ///need to used multimap since >1 muon 
 ///value for each run number key
 std::multimap < int, int> runMap;

  ///loop over all events
  for ( int i = 0; i < tree->GetEntries(); ++i ) {
 //   if(i>10000) break;
    tree->LoadTree(i);
    tree->GetEntry(i);

    if(i%10000==0) std::cout << "Event " << i << std::endl;

        for(int imu = 0; imu<mu_muid_nNt; imu++){
            if (
                (valNt[imu]>cutValue) 
                && abs(etaNt[imu])>0.1 
                && abs(etaNt[imu])<2.4 
                && ptNt[imu]>ptmin 
                && centralityNt>=0. 
                && centralityNt<=0.8
                && ( (EF_mu10_MSonly_EFFS_L1ZDC&&EF_mu10_MSonly_EFFS_L1ZDC_Matched20[imu]) ||
                    (EF_mu10_MSonly_EFFS_L1TE10&&EF_mu10_MSonly_EFFS_L1TE10_Matched20[imu])
                    ||(EF_mu10_MSonly_EFFS_L1TE20&&EF_mu10_MSonly_EFFS_L1TE20_Matched20[imu]) )	
		        ){
                ///count number of muons in this run that pass
                runMap.insert(make_pair (runNt,1) );
               }
         } //imu
   } //i
                
   std::cout << std::endl;
   std::cout << "-------------------------  " << std::endl;
   std::cout << " S U M M A R Y " << std::endl;
   std::cout << "-------------------------  " << std::endl;

   TString spreadSheetName = "muonCountingSpreadSheet.csv";
   std::ofstream spreadSheet;
   spreadSheet.open(spreadSheetName);

   double totLumi = 0.0, totalMu = 0;
   for(int iRun = 0; iRun<nRuns; ++iRun){

        int muonsPassed = runMap.count( vRun[iRun] );
        ///Muons per nb^-1
        double muPerLumi = (double)muonsPassed/vLumi[iRun]/1000.0;
        ///Number of muons passed in run iRun
        std::cout << "Index: " << iRun+1 << "\tRun: " << vRun[iRun] << 
            "\tLumi: " << vLumi[iRun]/1000.0 << "\tMuons: " << muonsPassed << "\tMuons/nb^-1: " 
            << muPerLumi << std::endl;

        totLumi += vLumi[iRun]/1000.0;
        totalMu += muonsPassed;

        grMuPerLumi->SetPoint(iRun,vRun[iRun],muPerLumi);
        grMuPerLumi->SetPointError(iRun,0.5, TMath::Sqrt(muonsPassed)/vLumi[iRun]/1000.0);

        grLumi->SetPoint(iRun,totLumi,totalMu);
        grLumi->SetPointError(iRun, 0.0, 0.0);

        writeToSpreadsheet(spreadSheet,iRun+1,vRun[iRun],vLumi[iRun]/1000.,muonsPassed);
   }

   grMuPerLumi->Draw("ape");
   grLumi->Draw("ape");

   outFile->cd();
   grLumi->Write("muonLumiDependence");
   grMuPerLumi->Write("muonPerLumi");

   std::cout << " Total Luminosity: " << totLumi << " nb^-1 " << std::endl;
   std::cout << " Total Muons passed: " << totalMu << std::endl;

   spreadSheet.close();
                
/*
            
  cout << "Number of runs: " << runMap.size() << endl;
  cout << "Run \t" << "Luminosity(ub^-1) \t" << " #muons \t" << "muons/lumi" << endl;
  ///map run lumi to #mu
  //std::map < int, double> lumiMap;
  for (unsigned int j = 0; j<runMap.size(); ++j){
    int runTemp = vRun.at(j);
    double lumiTemp = vLumi.at(j);
    map <int,int>::const_iterator iPairFound = runMap.find(runTemp);
    cout << iPairFound->first << " \t " << lumiTemp << " \t " << " \t " << iPairFound->second << " \t " << " \t " << (double)iPairFound->second/lumiTemp << endl; 
    sum += runMap[runTemp];
  }

  cout << "Total Luminosity = " << totLumi << endl;
  cout << "Total = " << sum << endl;
  */
}
