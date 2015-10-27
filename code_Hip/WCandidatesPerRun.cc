#include "TFile.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TH1.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"

#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

void fillRunLumi(std::string input, std::map<std::string,double>& map, std::map<int,int>& mapRunToMuons, 
                    std::map<int,int>& mapRunToMuonsA,std::map<int,int>& mapRunToMuonsB){

    std::cout << "Reading runs and luminosity from " << input << std::endl;
    std::ifstream s( input.c_str() , std::ifstream::in );

    if(!s.is_open()){
        std::cout << "While opening a file an error was encountered" << std::endl;
        exit(1);
    }
    else{
        std::cout << "File " << input << " is successfully opened" <<
            std::endl;
    }

        std::string line;

        double recLumi = 0.0;
        while(std::getline( s, line )) {
            std::stringstream ss( line );
            std::string run;
            // first column
            std::getline(ss, run, ',');
            std::string sLumi;
            // second column
            std::getline(ss, sLumi,'\n');
            std::stringstream ssLumi(sLumi);
            double lumi;
            ssLumi >> lumi;
            lumi *= 1.0e-6;
            recLumi+=lumi;
            map.insert (std::make_pair <std::string, double>(run, lumi) );

            std::stringstream ssrun(run);
            int runNumber;
            ssrun >> runNumber;
            mapRunToMuons.insert(std::make_pair<int,int>(runNumber, 0) );
            // Period B
            if(runNumber>218589) mapRunToMuonsB.insert (std::make_pair <int, int>(runNumber, 0) ); 
            // Period A
            else mapRunToMuonsA.insert (std::make_pair <int, int>(runNumber, 0) );
        }

        std::cout << "-------------------------------" << std::endl;
        std::cout << "Total recorded luminosity: " << recLumi << "[nb-1]" << std::endl;
        std::cout << "-------------------------------" << std::endl;

}

void convertToHisto(TH1F* h, TGraphErrors* gr){
   for(int i=0; i<gr->GetN(); ++i){
      float ypt = gr->GetY()[i];
      float erry = gr->GetEY()[i];
      h->SetBinContent(i+1,ypt);
      h->SetBinError(i+1,erry);
   }
}

void WCandidatesPerRun(){

    std::map<std::string,double> mapRunToLumi;
    std::map<int,int> mapRunToMuons;
    std::map<int,int> mapRunToMuonsA;
    std::map<int,int> mapRunToMuonsB;
    fillRunLumi("Hip_Lumi_Calc.txt", mapRunToLumi, mapRunToMuons,mapRunToMuonsA,mapRunToMuonsB);

    // output spreadsheet
    TString spreadSheetName = "Hip_WBosons.csv";

    std::ofstream spreadSheet;
    spreadSheet.open(spreadSheetName);

	TChain* tree = new TChain("treeW","treeW");
	tree->Add("$HOME/scratch/Data2013Hip/ntupleW.07.24.2014.root");
    std::cout << "Number of entries in data tree: " << tree->GetEntries() << std::endl;

    int   _run = -9999;
    int   _event = -9999;
    bool  _isPeriodA = false;
    bool  _isPeriodB = false;
	float _ptcone30 = -999.;
	float _pt = -999.;
	tree->SetBranchAddress("run",&_run);
	tree->SetBranchAddress("event",&_event);
	tree->SetBranchAddress("isPeriodA",&_isPeriodA);
	tree->SetBranchAddress("isPeriodB",&_isPeriodB);
	tree->SetBranchAddress("ptcone30",&_ptcone30);
	tree->SetBranchAddress("pt",&_pt);
	// --- Set branch status ---
  	tree->SetBranchStatus("*",0) ;
	tree->SetBranchStatus("run",1);
	tree->SetBranchStatus("event",1);
	tree->SetBranchStatus("isPeriodA",1);
	tree->SetBranchStatus("isPeriodB",1);
	tree->SetBranchStatus("ptcone30",1);
	tree->SetBranchStatus("pt",1);

    int nmu = tree->GetEntries();
    for(int imu =0; imu < nmu; ++imu){

         tree->LoadTree(imu);
         tree->GetEntry(imu);
         // C++11 only
         //thisRun.to_string(_run);
         
	     float isoRatio = _ptcone30/_pt;
         if(isoRatio<0.2){
             if(_isPeriodA) ++mapRunToMuonsA[_run]; 
             if(_isPeriodB) ++mapRunToMuonsB[_run]; 
             ++mapRunToMuons[_run] ; 
             // Write to spreadsheet
             spreadSheet << _run << "," << _event << std::endl;  
         }
         
         
    }//imu

    TGraphErrors* gr = new TGraphErrors();
    std::map<std::string,double>::iterator it;
    std::map<int,int>::iterator it2;
    int ipt;
    int TotalMuon = 0;
    double LintA=0.0, LintB=0.0;
    for( it=mapRunToLumi.begin(), it2=mapRunToMuons.begin(), ipt=0
            ; it!=mapRunToLumi.end()
            ; ++it,++it2, ++ipt){
        
        if(it2->first<=218589) LintA+=it->second;
        else LintB+=it->second;

        std::cout << it->first << " --> " << it->second << "[nb-1]"  << " : " << it2->second << " muons." << std::endl;
        TotalMuon += it2->second;
        double candPerLumi = (double)it2->second/it->second;
        gr->SetPoint(ipt,it2->first,candPerLumi);
        double error = TMath::Sqrt((double)it2->second)/it->second;
        gr->SetPointError(ipt,0,error);
        //h->GetXaxis()->SetBinLabel(ipt+1,it->first.c_str());
        //x = h->GetXaxis()->GetBinCenter(ipt+1);
        //t.DrawText(x,400.0,it->first.c_str());

    }

    TAxis *ax = gr->GetHistogram()->GetXaxis();
    Double_t x1 = ax->GetBinLowEdge(1);
    Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
    //gr->GetHistogram()->GetXaxis()->Set(mapRunToMuons.size(),x1,x2);
    TH1F* h = new TH1F("hDummy","hDummy",mapRunToMuons.size(),x1,x2);
    h->GetXaxis()->SetLabelOffset(99);
    TCanvas* c = new TCanvas("runStability","runStability",800,600);
    c->cd();
    
    std::cout << "==================================" << std::endl;
    std::cout << "Total events : " << TotalMuon << std::endl;
    std::cout << "Lint Period A: " <<  LintA << std::endl;
    std::cout << "Lint Period B: " <<  LintB << std::endl;

    h->GetYaxis()->SetRangeUser(400.0,850.0);
    //h->LabelsDeflate("X");
    //h->LabelsOption("d");
    //h->GetXaxis()->SetLabelSize(0.03);

    h->GetYaxis()->SetTitle("W#rightarrow#mu#nu events/nb^{-1}");
    convertToHisto(h,gr);
    // Plot run-by-run stabiliy
    h->Draw();

    TText t;
    Float_t x;
    t.SetTextAngle(60);
    t.SetTextSize(0.03);
    t.SetTextAlign(33);
    for( it=mapRunToLumi.begin(), ipt=0
            ; it!=mapRunToLumi.end()
            ; ++it,++ipt){
        x = h->GetXaxis()->GetBinCenter(ipt+1);
        t.DrawText(x,399.0,it->first.c_str());

    }
    //gr->Draw("pesame");
    TLegend *leg = new TLegend(0.1934673,0.7867133,0.4635678,0.9370629,NULL,"brNDC");
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->AddEntry(h,"Data #sqrt{s_{NN}} = 5.02 TeV");
    leg->Draw();
    TLatex *   tex = new TLatex(0.3404523,0.7534965,"27.8 nb^{-1}"); tex->SetNDC(); tex->Draw();
    tex = new TLatex(0.2437186,0.756993,"p+Pb"); tex->SetNDC(); tex->Draw();
    tex = new TLatex(0.2575377,0.8986014,"ATLAS Internal"); tex->SetNDC(); tex->Draw();

    c->Draw();
    c->Print("plots/runStability.root"); 
    c->Print("plots/runStability.C"); 

}
