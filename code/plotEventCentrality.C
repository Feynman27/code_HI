#include "TH1.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include <cmath>
#include <iostream>

void plotEventCentrality(){

    TFile* fIn = new TFile("../MinimumBiasFiles/HISingleMuonMinBiasData.07.10.2013.root","READ");
//    TFile* fIn = new TFile("../MinimumBiasFiles/HISingleMuonMB.2012.11.25.root","READ");

    TTree* tIn = (TTree*)fIn->Get("tree");

    float centrality;
    int trig1,trig2;
    tIn->SetBranchAddress("centrality",&centrality);
    tIn->SetBranchAddress("EF_mbZdc_a_c_L1VTE50_trk",&trig1);
    tIn->SetBranchAddress("EF_L1TE50_NoAlg",&trig2);
    tIn->SetBranchStatus("*",0);
    tIn->SetBranchStatus("centrality",1);
    tIn->SetBranchStatus("EF_mbZdc_a_c_L1VTE50_trk",1);
    tIn->SetBranchStatus("EF_L1TE50_NoAlg",1);

    std::vector<double> ncollErr;
    ncollErr.push_back(7.7);
    ncollErr.push_back(7.5);
    ncollErr.push_back(7.4);
    ncollErr.push_back(7.4);
    ncollErr.push_back(7.3);
    ncollErr.push_back(14.2);

    double classes[] = {0.0,0.05,0.1,0.15,0.2,0.4,0.8};
    int nBins = sizeof(classes)/sizeof(double)-1;
    TH1F* hcent = new TH1F("hcent","hcent",nBins,classes);

    std::cout<< "Entries: " << tIn->GetEntries() << std::endl;
    for(int iev=0; iev<tIn->GetEntries(); ++iev){

       tIn->LoadTree(iev); 
       tIn->GetEntry(iev); 

       if(iev%10000==0) std::cout << "Event: " << iev << std::endl;

       ///Did the event pass the MB trigger definition
       if(trig1||trig2) {
//           std::cout << "EF_mbZdc_a_c_L1VTE50_trk:" << trig1 << " " << "EF_L1TE50_NoAlg" << trig2 << std::endl;
           hcent->Fill(centrality);
       }

    }

    std::cout << "Number of centrality classes " << hcent->GetNbinsX() << std::endl;
    /*for(int ibin=1; ibin < hcent->GetNbinsX()+1; ++ibin){
    
        std::cout << "Ncoll error : " << ncollErr[ibin-1] << "%" << std::endl;
        double error = TMath::Sqrt(pow(hcent->GetBinError(ibin),2)+pow(hcent->GetBinContent(ibin)*ncollErr[ibin-1]/100.0,2));
        hcent->SetBinError(ibin,error);
    }
    */

    hcent->Draw("pe");
    hcent->Scale(1.0,"width");
    hcent->Draw("hist e");
}
