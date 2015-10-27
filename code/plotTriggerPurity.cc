#include <iostream>
#include <string>
#include <memory>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TEfficiency.h"

int main(){

    // outfile
    std::unique_ptr<TFile> _fOut(new TFile("triggerPurity.root","recreate"));
    TString pathName = "/usatlas/u/tbales/scratch/";
    TString fileName = "HISingleMuonHardProbesData.04.17.2013.root";

    float binArr[] = {10.,15.,20.,25.,27.,29.,31.,33.,35.,37.,39.,41.,43.,45.,47.,49.,54.,59.,64.,69.,79.,89.,100.,120.,140.,160.,200.};
    const int nbins = sizeof(binArr)/sizeof(float)-1;
    // pt histos for each trig
    TH1F* h_TrigUnion = new TH1F("h_TrigUnion","h_TrigUnion",nbins,binArr); 
    TH1F* h_TrigAll   = new TH1F("h_TrigAll","h_TrigAll",nbins,binArr); 
    TH1F* h_TrigZDC   = new TH1F("h_TrigZDC"  ,"h_TrigZDC",nbins,binArr); 
    TH1F* h_TrigTE10  = new TH1F("h_TrigTE10" ,"h_TrigTE10",nbins,binArr); 
    TH1F* h_TrigTE20  = new TH1F("h_TrigTE20" ,"h_TrigTE20",nbins,binArr); 

    float eLossNt[50];
    float scatNt[50];
    float compNt[50];
    float ptNt[50];
    float mtNt[50];
    float etaNt[50];
    float phiNt[50];
    float chargeNt[50];
    float centralityNt;
    float nu_ptNt;
    float ptconeNt[50];
    int valNt[50], ZDYNt[50], matched1[50], matched2[50], matched3[50];
    int nmu,trig1,trig2,trig3,trig4,trig5;
    int runNt,eventNt;

    std::unique_ptr<TChain> tree(new TChain("tree","tree"));
    int nFiles = tree->Add(pathName+fileName);

    std::cout << "Filling the dataSet for " << fileName << "... Number of files: " << nFiles << std::endl;

    tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC",&trig1);
    tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10",&trig2);
    tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20",&trig3);
    tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1ZDC_Matched20",&matched1);
    tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE10_Matched20",&matched2);
    tree->SetBranchAddress("EF_mu10_MSonly_EFFS_L1TE20_Matched20",&matched3);
    tree->SetBranchAddress("EF_mu4_MSonly_L1TE50",&trig4);
    tree->SetBranchAddress("EF_mu4_L1VTE50",&trig5);
    tree->SetBranchAddress("run", &runNt);
    tree->SetBranchAddress("event", &eventNt);
    tree->SetBranchAddress("eLoss", &eLossNt);
    tree->SetBranchAddress("ptcone20ID3", &ptconeNt);
    tree->SetBranchAddress("scat", &scatNt);
    tree->SetBranchAddress("comp", &compNt);
    tree->SetBranchAddress("pt", &ptNt);
    tree->SetBranchAddress("eta", &etaNt);
    tree->SetBranchAddress("phi", &phiNt);
    tree->SetBranchAddress("charge", &chargeNt);
    tree->SetBranchAddress("val", &valNt); 
    tree->SetBranchAddress("ZDY", &ZDYNt); 
    tree->SetBranchAddress("centrality", &centralityNt);
    tree->SetBranchAddress("nu_pt", &nu_ptNt);
    tree->SetBranchAddress("mt", &mtNt);
    tree->SetBranchAddress("mu_muid_n", &nmu);

    // --- Set branch status ---
    tree->SetBranchStatus("*",0) ;
    tree->SetBranchStatus("run", 1);
    tree->SetBranchStatus("event", 1);
    tree->SetBranchStatus("mu_muid_n", 1);
    tree->SetBranchStatus("eLoss", 1);
    tree->SetBranchStatus("ptcone20ID3", 1);
    tree->SetBranchStatus("scat", 1);
    tree->SetBranchStatus("comp", 1);
    tree->SetBranchStatus("pt", 1);
    tree->SetBranchStatus("mt", 1);
    tree->SetBranchStatus("eta", 1);
    tree->SetBranchStatus("phi", 1);
    tree->SetBranchStatus("charge", 1);
    tree->SetBranchStatus("val", 1); 
    tree->SetBranchStatus("ZDY", 1); 
    tree->SetBranchStatus("centrality", 1);
    tree->SetBranchStatus("nu_pt", 1);
    tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC",1);
    tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10",1);
    tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20",1);
    tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1ZDC_Matched20",1);
    tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE10_Matched20",1);
    tree->SetBranchStatus("EF_mu10_MSonly_EFFS_L1TE20_Matched20",1);
    tree->SetBranchStatus("EF_mu4_MSonly_L1TE50",1);
    tree->SetBranchStatus("EF_mu4_L1VTE50",1);

    std::cout << "Number of entries: " << tree->GetEntries() << std::endl; 
    for ( int i = 0; i < tree->GetEntries(); i++ ) {
        if(i%10000==0) std::cout << "Event: " << i << std::endl;
        tree->LoadTree(i);
        tree->GetEntry(i);
        // 1 of 3 must have fired
        if( (trig1)||(trig2)||(trig3)){
            for (int imu = 0; imu<nmu;imu++){
                // all triggered events
                if(

                   valNt[imu]>11
                   && fabs(etaNt[imu])>0.1 && fabs(etaNt[imu])<2.4
                   && centralityNt >= 0.0 && centralityNt <= 0.8
                   && (matched1[imu]||matched2[imu]||matched3[imu])
                ) h_TrigAll->Fill(ptNt[imu]);
                else continue;
                if(
                   fabs(eLossNt[imu])<0.5
                   && fabs(scatNt[imu])<4.0
                   && ptNt[imu] > 25.0
                   && nu_ptNt > 25.0 && nu_ptNt < 8000.0
                   && (ptconeNt[imu]/ptNt[imu]) <0.1
                   && mtNt[imu] > 40.0 
                   && ZDYNt[imu] == 0 
                  ){ //muon sel
                    h_TrigUnion->Fill(ptNt[imu]);
                    /*if(matched1[imu]) h_TrigZDC->Fill(ptNt[imu]);
                    if(matched2[imu]) h_TrigTE10->Fill(ptNt[imu]);
                    if(matched3[imu]) h_TrigTE20->Fill(ptNt[imu]);
                    */
                 }
            }//imu
        }//trigger
        else continue;
    }//i

    std::unique_ptr<TEfficiency> _eff(new TEfficiency  (*h_TrigUnion,  *h_TrigAll));
    /*std::unique_ptr<TEfficiency> _effZDC(new TEfficiency  (*h_TrigZDC,  *h_TrigUnion));
    std::unique_ptr<TEfficiency> _effTE10(new TEfficiency (*h_TrigTE10, *h_TrigUnion));
    std::unique_ptr<TEfficiency> _effTE20(new TEfficiency (*h_TrigTE20, *h_TrigUnion));
*/
    _fOut->cd();
    _eff->Write("trig_purity");
    /*_effZDC->Write("EF_mu10_MSonly_EFFS_L1ZDC");
    _effTE10->Write("EF_mu10_MSonly_EFFS_L1TE10");
    _effTE20->Write("EF_mu10_MSonly_EFFS_L1TE20");
    */
}
