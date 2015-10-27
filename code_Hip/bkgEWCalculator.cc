#include<memory>
#include<iostream>
#include<vector>
#include<sstream>
#include<map>

#include "TString.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"

#include "BayesianErrors.cc"

int getBinIndex(float var, std::vector<double> vBins){
   for(unsigned int ibin=0; ibin<vBins.size(); ++ibin){
        if(var>vBins[ibin]&&var<vBins[ibin+1]) return ibin;
   }
   return -1;
}

void calculateBkg(std::map<std::string,int>& m_nFid, std::map<std::string,int>& m_nGen, int ipt, int ieta, int icent, std::string beamType, std::string charge, 
        float sigma, float Lint, float& nBkg, float& nBkgErrUpp, float& nBkgErrLow){

    bool isPlus = false;
    bool isMinus = false;
    if(strcmp(charge.c_str(),"+")==0) isPlus=true;
    else if(strcmp(charge.c_str(),"-")==0) isMinus=true;
    std::stringstream binFid,binGen;
    if(isPlus) {
        binGen << "muPlus_" << beamType;
        binFid << "muPlus_" << beamType << "_pt" << ipt << "_eta" << ieta << "_cent" << icent;
    }
    else if(isMinus) {
        binGen << "muMinus_" << beamType ;
        binFid << "muMinus_" << beamType << "_pt" << ipt << "_eta" << ieta << "_cent" << icent;
    }

    double eff=0.0; 
    double eff_err_low=0.0,eff_err_upp=0.0;
    Efficiency(m_nFid[binFid.str()],m_nGen[binGen.str()],0.683,eff,eff_err_low,eff_err_upp);
    //const float Lint = 27.8;
    const float A = 208.0;
    float scale = sigma*Lint*A;
    // fraction of all generated muons passing signal selection
    nBkg = eff*scale;
    nBkgErrLow = nBkg-scale*eff_err_low;
    nBkgErrUpp = scale*eff_err_upp-nBkg;

}

void calculateGen(TString sFile1,TString sFile2, int thisMother, std::string beamType, std::string charge, std::map<std::string,int>& m_nGen){
    bool isPlus = false;
    bool isMinus = false;
    if(strcmp(charge.c_str(),"+")==0) isPlus=true;
    else if(strcmp(charge.c_str(),"-")==0) isMinus=true;

    std::unique_ptr<TChain> _tGen(new TChain("treeGen","treeGen"));
    _tGen->Add(sFile1);
    _tGen->Add(sFile2);
    std::cout << "Number of entries in generator tree: " << _tGen->GetEntries() << std::endl;

    float genPt, genEta, genCentrality, genCharge;
    float fcalGen;
    int genPDG,genMother; 
    bool isGenPeriodA, isGenPeriodB;

    _tGen->SetBranchAddress("muon_gen_pt",&genPt);
    _tGen->SetBranchAddress("muon_gen_eta",&genEta);
    _tGen->SetBranchAddress("centrality",&genCentrality);
    //_tGen->SetBranchAddress("neutrino_gen_pt",&genNeutrinoPt);
    //_tGen->SetBranchAddress("mt_gen_neutrino_muon",&genMt);
    _tGen->SetBranchAddress("muon_gen_charge",&genCharge);
    _tGen->SetBranchAddress("muon_pdg",&genPDG);
    //_tGen->SetBranchAddress("neutrino_pdg",&genNeutrinoPDG);
    _tGen->SetBranchAddress("muon_mother",&genMother);
    _tGen->SetBranchAddress("isPeriodA",&isGenPeriodA);
    _tGen->SetBranchAddress("isPeriodB",&isGenPeriodB);
    _tGen->SetBranchAddress("fcal",&fcalGen);

    // Total number of generated events 
    // over all space
    for(int igen=0; igen<_tGen->GetEntries(); ++igen){
        _tGen->GetEntry(igen);
        if(fabs(genPDG) == 13 && fabs(genMother) == thisMother){

            std::stringstream bin;
            if(genCharge>0.0&&isPlus==true){
             bin << "muPlus_" << beamType ;
             ++m_nGen[bin.str()];
            }
            else if (genCharge<0.0&&isMinus==true){
             bin << "muMinus_" << beamType ;
             ++m_nGen[bin.str()];
            }
        }
    }//generator 
}

void calculateFid(TString sFile1,TString sFile2, int thisMother, std::string beamType, std::string charge, std::map<std::string,int>& m_nFid,   
                    std::vector<double> vPt, std::vector<double> vEta, std::vector<double> vCent){
    bool isPlus = false;
    bool isMinus = false;
    if(strcmp(charge.c_str(),"+")==0) isPlus=true;
    else if(strcmp(charge.c_str(),"-")==0) isMinus=true;

    const float ptLow=vPt[0], ptUpp=vPt[vPt.size()-1];
    const float etaLow=vEta[0], etaUpp=vEta[vEta.size()-1];
    const float centLow=vCent[0], centUpp=vCent[vCent.size()-1];

    std::unique_ptr<TChain> _tRec(new TChain("treeW","treeW"));
    _tRec->Add(sFile1);
    _tRec->Add(sFile2);
    std::cout << "Number of entries in reconstructed tree: " << _tRec->GetEntries() << std::endl;

    float recPt, recEta, recCentrality, recMET, recMt, recCharge,recPtcone30;
    int recPDG, recMother, recMatched; 
    float fcalRec;
    bool isRecPeriodA, isRecPeriodB;
    bool isZCandidate, EF_mu8_matched;

    _tRec->SetBranchAddress("pt",&recPt);
    _tRec->SetBranchAddress("eta",&recEta);
    _tRec->SetBranchAddress("centrality",&recCentrality);
    _tRec->SetBranchAddress("met",&recMET);
    _tRec->SetBranchAddress("mt",&recMt);
    _tRec->SetBranchAddress("charge",&recCharge);
    _tRec->SetBranchAddress("rec_truth_pdg",&recPDG);
    _tRec->SetBranchAddress("rec_truth_mother_pdg",&recMother);
    _tRec->SetBranchAddress("ptcone30",&recPtcone30);
    _tRec->SetBranchAddress("rec_truth_matched",&recMatched);
    _tRec->SetBranchAddress("rec_truth_mother_pdg",&recMother);
    _tRec->SetBranchAddress("rec_truth_pdg",&recPDG);
    _tRec->SetBranchAddress("isPeriodA",&isRecPeriodA);
    _tRec->SetBranchAddress("isPeriodB",&isRecPeriodB);
    _tRec->SetBranchAddress("isZCandidateType1",&isZCandidate);
    _tRec->SetBranchAddress("EF_mu8_matched",&EF_mu8_matched);
    _tRec->SetBranchAddress("fcal",&fcalRec);


    const float metCut = 25.0;
    const float mtCut  = 40.0;
    // Background in signal space
    for(int irec=0; irec<_tRec->GetEntries(); ++irec){
        _tRec->GetEntry(irec);

        float isoRatio = recPtcone30/recPt;
        if(
            recPt>ptLow && recPt<ptUpp
            && recEta>etaLow && recEta<etaUpp
            && recCentrality>centLow  && recCentrality<centUpp 
            && isoRatio <= 0.2
            && recMET > metCut  
            && recMt > mtCut 
            && recMatched==1
            && isZCandidate == false     
            && EF_mu8_matched == 1
            && ( (recCharge>0.0&&isPlus==true) || (recCharge<0.0&&isMinus==true) )

          ) {
            int ipt = getBinIndex(recPt,vPt);
            int ieta = getBinIndex(recEta,vEta);
            int icent = getBinIndex(recCentrality,vCent);
            if(ipt<0 || ieta<0 || icent<0){
                std::cout << "ERROR: Cannot locate bin to which this muon belongs. Expect wrong results." << std::endl;
                return;
            }
            std::stringstream bin;
            if(recCharge>0.0){ 
                bin << "muPlus_" << beamType << "_pt" << ipt << "_eta" << ieta << "_cent" << icent;
                ++m_nFid[bin.str()];
            }
            else if(recCharge<0.0) {
                bin << "muMinus_" << beamType << "_pt" << ipt << "_eta" << ieta << "_cent" << icent;
                ++m_nFid[bin.str()];
            }
        }
    }//reco 

}

int main(){
    std::unique_ptr<TFile> _fOut(new TFile("bkgEW.root","recreate")); 
    std::vector<double> vPt;
    vPt.push_back(25.0);
    vPt.push_back(400.0);
    const int nPtBins = vPt.size()-1;
    std::vector<double> vEta;
    vEta.push_back(-2.5);
    vEta.push_back(-2.1);
    vEta.push_back(-1.55);
    vEta.push_back(-1.05);
    vEta.push_back(-0.6);
    vEta.push_back(-0.1);
    vEta.push_back(0.1);
    vEta.push_back(0.6);
    vEta.push_back(1.05);
    vEta.push_back(1.55);
    vEta.push_back(2.1);
    vEta.push_back(+2.5);
    const int nEtaBins = vEta.size()-1;
    std::vector<double> vCent;
    vCent.push_back(0.00);
    vCent.push_back(0.05);
    vCent.push_back(0.10);
    vCent.push_back(0.20);
    vCent.push_back(0.30); // added to match Z
    vCent.push_back(0.40);
    vCent.push_back(0.60); // added to match Z
    vCent.push_back(0.90);
    const int nCentralityBins = vCent.size()-1;

    TString sBase = "/usatlas/u/tbales/scratch/Data2013Hip/";
    std::map<std::string,int> m_nFid;
    std::map<std::string,int> m_nGen;
    m_nGen.insert(std::make_pair("muPlus_Z_pp",0));
    m_nGen.insert(std::make_pair("muPlus_Z_pn",0));
    m_nGen.insert(std::make_pair("muMinus_Z_pp",0));
    m_nGen.insert(std::make_pair("muMinus_Z_pn",0));

    m_nGen.insert(std::make_pair("muPlus_tau_pp",0));
    m_nGen.insert(std::make_pair("muPlus_tau_pn",0));
    m_nGen.insert(std::make_pair("muMinus_tau_pp",0));
    m_nGen.insert(std::make_pair("muMinus_tau_pn",0));

    std::map<std::string,TGraphAsymmErrors*> m_gr;
    for(int ipt=0; ipt<nPtBins; ++ipt){
        for(int icent=0; icent<nCentralityBins; ++icent){
            std::stringstream grName; 
            grName << "grA_muPlus_Z_pt" << ipt << "_cent" << icent;
            m_gr.insert(std::make_pair(grName.str(),new TGraphAsymmErrors(nEtaBins)));
            grName.str(std::string());
            grName << "grA_muMinus_Z_pt" << ipt << "_cent" << icent;
            m_gr.insert(std::make_pair(grName.str(),new TGraphAsymmErrors(nEtaBins)));
            grName.str(std::string());
            grName << "grA_muPlus_tau_pt" << ipt << "_cent" << icent;
            m_gr.insert(std::make_pair(grName.str(),new TGraphAsymmErrors(nEtaBins)));
            grName.str(std::string());
            grName << "grA_muMinus_tau_pt" << ipt << "_cent" << icent;
            m_gr.insert(std::make_pair(grName.str(),new TGraphAsymmErrors(nEtaBins)));

            grName.str(std::string());
            grName << "grB_muPlus_Z_pt" << ipt << "_cent" << icent;
            m_gr.insert(std::make_pair(grName.str(),new TGraphAsymmErrors(nEtaBins)));
            grName.str(std::string());
            grName << "grB_muMinus_Z_pt" << ipt << "_cent" << icent;
            m_gr.insert(std::make_pair(grName.str(),new TGraphAsymmErrors(nEtaBins)));
            grName.str(std::string());
            grName << "grB_muPlus_tau_pt" << ipt << "_cent" << icent;
            m_gr.insert(std::make_pair(grName.str(),new TGraphAsymmErrors(nEtaBins)));
            grName.str(std::string());
            grName << "grB_muMinus_tau_pt" << ipt << "_cent" << icent;
            m_gr.insert(std::make_pair(grName.str(),new TGraphAsymmErrors(nEtaBins)));
            grName.str(std::string());

            for(int ieta=0; ieta<nEtaBins; ++ieta){
               
                std::stringstream bin;
                // Z-->mumu
                bin << "muPlus_Z_pp_" << "pt" << ipt << "_eta" << ieta << "_cent" << icent;
                m_nFid.insert(std::make_pair(bin.str(),0));
                bin.str(std::string());
                bin << "muPlus_Z_pn_" << "pt" << ipt << "_eta" << ieta << "_cent" << icent;
                m_nFid.insert(std::make_pair(bin.str(),0));

                bin.str(std::string());
                bin << "muMinus_Z_pp_" << "pt" << ipt << "_eta" << ieta << "_cent" << icent;
                m_nFid.insert(std::make_pair(bin.str(),0));
                bin.str(std::string());
                bin << "muMinus_Z_pn_" << "pt" << ipt << "_eta" << ieta << "_cent" << icent;
                m_nFid.insert(std::make_pair(bin.str(),0));

                // W-->tau
                bin.str(std::string());
                bin << "muPlus_tau_pp_" << "pt" << ipt << "_eta" << ieta << "_cent" << icent;
                m_nFid.insert(std::make_pair(bin.str(),0));
                bin.str(std::string());
                bin << "muPlus_tau_pn_" << "pt" << ipt << "_eta" << ieta << "_cent" << icent;
                m_nFid.insert(std::make_pair(bin.str(),0));

                bin.str(std::string());
                bin << "muMinus_tau_pp_" << "pt" << ipt << "_eta" << ieta << "_cent" << icent;
                m_nFid.insert(std::make_pair(bin.str(),0));
                bin.str(std::string());
                bin << "muMinus_tau_pn_" << "pt" << ipt << "_eta" << ieta << "_cent" << icent;
                m_nFid.insert(std::make_pair(bin.str(),0));
            }//ieta
        }//icent
    }//ipt

    // Z-->mumu
    // mu^{+}
    calculateFid(sBase+"ntupleZ_pp.A.09.29.2014.root",sBase+"ntupleZ_pp.B.09.29.2014.root",23,"Z_pp","+", m_nFid, vPt, vEta, vCent);
    calculateGen(sBase+"ntupleZ_pp.A.09.29.2014.root",sBase+"ntupleZ_pp.B.09.29.2014.root",23,"Z_pp","+", m_nGen);
    calculateFid(sBase+"ntupleZ_pn.A.09.29.2014.root",sBase+"ntupleZ_pn.B.09.29.2014.root",23,"Z_pn","+", m_nFid, vPt, vEta, vCent);
    calculateGen(sBase+"ntupleZ_pn.A.09.29.2014.root",sBase+"ntupleZ_pn.B.09.29.2014.root",23,"Z_pn","+", m_nGen);
    // mu^{-}
    calculateFid(sBase+"ntupleZ_pp.A.09.29.2014.root",sBase+"ntupleZ_pp.B.09.29.2014.root",23,"Z_pp","-", m_nFid, vPt, vEta, vCent);
    calculateGen(sBase+"ntupleZ_pp.A.09.29.2014.root",sBase+"ntupleZ_pp.B.09.29.2014.root",23,"Z_pp","-", m_nGen);
    calculateFid(sBase+"ntupleZ_pn.A.09.29.2014.root",sBase+"ntupleZ_pn.B.09.29.2014.root",23,"Z_pn","-", m_nFid, vPt, vEta, vCent);
    calculateGen(sBase+"ntupleZ_pn.A.09.29.2014.root",sBase+"ntupleZ_pn.B.09.29.2014.root",23,"Z_pn","-", m_nGen);

    // W-->tau 
    // mu^{+}
    calculateFid(sBase+"ntupleWtauplusmunu_pp.A.09.29.2014.root",sBase+"ntupleWtauplusmunu_pp.B.09.29.2014.root",15,"tau_pp","+", m_nFid, vPt, vEta, vCent);
    calculateGen(sBase+"ntupleWtauplusmunu_pp.A.09.29.2014.root",sBase+"ntupleWtauplusmunu_pp.B.09.29.2014.root",15,"tau_pp","+", m_nGen);
    calculateFid(sBase+"ntupleWtauplusmunu_pn.A.09.29.2014.root",sBase+"ntupleWtauplusmunu_pn.B.09.29.2014.root",15,"tau_pn","+", m_nFid, vPt, vEta, vCent);
    calculateGen(sBase+"ntupleWtauplusmunu_pn.A.09.29.2014.root",sBase+"ntupleWtauplusmunu_pn.B.09.29.2014.root",15,"tau_pn","+", m_nGen);
    // mu^{-}
    calculateFid(sBase+"ntupleWtauminmunu_pp.A.09.29.2014.root",sBase+"ntupleWtauminmunu_pp.B.09.29.2014.root",15,"tau_pp","-", m_nFid, vPt, vEta, vCent);
    calculateGen(sBase+"ntupleWtauminmunu_pp.A.09.29.2014.root",sBase+"ntupleWtauminmunu_pp.B.09.29.2014.root",15,"tau_pp","-", m_nGen);
    calculateFid(sBase+"ntupleWtauminmunu_pn.A.09.29.2014.root",sBase+"ntupleWtauminmunu_pn.B.09.29.2014.root",15,"tau_pn","-", m_nFid, vPt, vEta, vCent);
    calculateGen(sBase+"ntupleWtauminmunu_pn.A.09.29.2014.root",sBase+"ntupleWtauminmunu_pn.B.09.29.2014.root",15,"tau_pn","-", m_nGen);

    for(auto iter=m_nFid.begin(); iter!=m_nFid.end(); ++iter){
        std::cout << iter->first << "--> " << iter->second << std::endl;
    }
    
    for(auto iter=m_nGen.begin(); iter!=m_nGen.end(); ++iter){
        std::cout << iter->first << "--> " << iter->second << std::endl;
    }

    const double wt_pp = 82.0/208.0;
    const double wt_pn = (208.0-82.0)/208.0;
    // [nb] (AMI)
    const float sigma_Z_pp = 6.4287E-01;
    const float sigma_Z_pn = 6.5651E-01;
    // tau-->mu B.R.
    const float BR = 0.1741;
    const float sigma_tau_plus_pp  = BR*4.5545E-01*4.1918E+00;
    const float sigma_tau_plus_pn  = BR*4.6044E-01*3.4616E+00;
    const float sigma_tau_minus_pp = BR*4.1654E-01*2.7368E+00;
    const float sigma_tau_minus_pn = BR*3.9080E-01*3.4404E+00;
    const float LintA=15.3524, LintB=12.4683;

    for(int ipt=0; ipt<nPtBins; ++ipt){
        for(int icent=0; icent<nCentralityBins; ++icent){
            for(int ieta=0; ieta<nEtaBins; ++ieta){

                float nBkg_pp=0.0,nBkg_pn=0.0,nBkg=0.0;
                float nBkg_pp_err_upp=0.0, nBkg_pp_err_low=0.0, nBkg_pn_err_upp=0.0, nBkg_pn_err_low=0.0, nBkg_err_upp=0.0, nBkg_err_low=0.0;

                // Z-->mumu
                // Period A
                // mu^{+}
                calculateBkg(m_nFid,m_nGen,ipt,ieta,icent,"Z_pp","+",sigma_Z_pp,LintA,nBkg_pp,nBkg_pp_err_upp, nBkg_pp_err_low);
                calculateBkg(m_nFid,m_nGen,ipt,ieta,icent,"Z_pn","+",sigma_Z_pn,LintA,nBkg_pn,nBkg_pn_err_upp, nBkg_pn_err_low);
                nBkg = wt_pp*nBkg_pp + wt_pn*nBkg_pn;
                nBkg_err_upp = wt_pp*nBkg_pp_err_upp + wt_pn*nBkg_pn_err_upp;
                nBkg_err_low = wt_pp*nBkg_pp_err_low + wt_pn*nBkg_pn_err_low;
                
                std::cout << "Number of Period A mu^{+} background events from Z in pt: " << ipt << " eta: " << ieta << " cent: " << icent << " = " 
                    << nBkg << " +/- " << nBkg_err_upp << ", " << nBkg_err_low << std::endl;
                
                std::stringstream grName;
                double x = (vEta[ieta]+vEta[ieta+1])/2.0;
                double xerr = (vEta[ieta+1]-vEta[ieta])/2.0;
                grName << "grA_muPlus_Z_pt" << ipt << "_cent" << icent;
                m_gr[grName.str()]->SetPoint(ieta,x,nBkg);
                m_gr[grName.str()]->SetPointError(ieta,xerr,xerr,nBkg_err_low,nBkg_err_upp);
                grName.str(std::string());

                // B
                calculateBkg(m_nFid,m_nGen,ipt,ieta,icent,"Z_pp","+",sigma_Z_pp,LintB,nBkg_pp,nBkg_pp_err_upp, nBkg_pp_err_low);
                calculateBkg(m_nFid,m_nGen,ipt,ieta,icent,"Z_pn","+",sigma_Z_pn,LintB,nBkg_pn,nBkg_pn_err_upp, nBkg_pn_err_low);
                nBkg = wt_pp*nBkg_pp + wt_pn*nBkg_pn;
                nBkg_err_upp = wt_pp*nBkg_pp_err_upp + wt_pn*nBkg_pn_err_upp;
                nBkg_err_low = wt_pp*nBkg_pp_err_low + wt_pn*nBkg_pn_err_low;
                
                std::cout << "Number of Period B mu^{+} background events from Z in pt: " << ipt << " eta: " << ieta << " cent: " << icent << " = " 
                    << nBkg << " +/- " << nBkg_err_upp << ", " << nBkg_err_low << std::endl;
                grName << "grB_muPlus_Z_pt" << ipt << "_cent" << icent;
                m_gr[grName.str()]->SetPoint(ieta,x,nBkg);
                m_gr[grName.str()]->SetPointError(ieta,xerr,xerr,nBkg_err_low,nBkg_err_upp);
                grName.str(std::string());

                // mu^{-}
                // A
                calculateBkg(m_nFid,m_nGen,ipt,ieta,icent,"Z_pp","-",LintA,sigma_Z_pp,nBkg_pp,nBkg_pp_err_upp, nBkg_pp_err_low);
                calculateBkg(m_nFid,m_nGen,ipt,ieta,icent,"Z_pn","-",LintA,sigma_Z_pn,nBkg_pn,nBkg_pn_err_upp, nBkg_pn_err_low);
                nBkg = wt_pp*nBkg_pp + wt_pn*nBkg_pn;
                nBkg_err_upp = wt_pp*nBkg_pp_err_upp + wt_pn*nBkg_pn_err_upp;
                nBkg_err_low = wt_pp*nBkg_pp_err_low + wt_pn*nBkg_pn_err_low;
                
                std::cout << "Number of Period A mu^{-} background events from Z in pt: " << ipt << " eta: " << ieta << " cent: " << icent << " = " 
                    << nBkg << " +/- " << nBkg_err_upp << ", " << nBkg_err_low << std::endl;
                grName << "grA_muMinus_Z_pt" << ipt << "_cent" << icent;
                m_gr[grName.str()]->SetPoint(ieta,x,nBkg);
                m_gr[grName.str()]->SetPointError(ieta,xerr,xerr,nBkg_err_low,nBkg_err_upp);
                grName.str(std::string());

                // B
                calculateBkg(m_nFid,m_nGen,ipt,ieta,icent,"Z_pp","-",LintB,sigma_Z_pp,nBkg_pp,nBkg_pp_err_upp, nBkg_pp_err_low);
                calculateBkg(m_nFid,m_nGen,ipt,ieta,icent,"Z_pn","-",LintB,sigma_Z_pn,nBkg_pn,nBkg_pn_err_upp, nBkg_pn_err_low);
                nBkg = wt_pp*nBkg_pp + wt_pn*nBkg_pn;
                nBkg_err_upp = wt_pp*nBkg_pp_err_upp + wt_pn*nBkg_pn_err_upp;
                nBkg_err_low = wt_pp*nBkg_pp_err_low + wt_pn*nBkg_pn_err_low;
                
                std::cout << "Number of Period B mu^{-} background events from Z in pt: " << ipt << " eta: " << ieta << " cent: " << icent << " = " 
                    << nBkg << " +/- " << nBkg_err_upp << ", " << nBkg_err_low << std::endl;
                grName << "grB_muMinus_Z_pt" << ipt << "_cent" << icent;
                m_gr[grName.str()]->SetPoint(ieta,x,nBkg);
                m_gr[grName.str()]->SetPointError(ieta,xerr,xerr,nBkg_err_low,nBkg_err_upp);
                grName.str(std::string());

                // W-->tau 
                // mu^{+}
                // A
                calculateBkg(m_nFid,m_nGen,ipt,ieta,icent,"tau_pp","+",sigma_tau_plus_pp,LintA,nBkg_pp,nBkg_pp_err_upp, nBkg_pp_err_low);
                calculateBkg(m_nFid,m_nGen,ipt,ieta,icent,"tau_pn","+",sigma_tau_plus_pn,LintA,nBkg_pn,nBkg_pn_err_upp, nBkg_pn_err_low);
                nBkg = wt_pp*nBkg_pp + wt_pn*nBkg_pn;
                nBkg_err_upp = wt_pp*nBkg_pp_err_upp + wt_pn*nBkg_pn_err_upp;
                nBkg_err_low = wt_pp*nBkg_pp_err_low + wt_pn*nBkg_pn_err_low;
                
                std::cout << "Number of Period A mu^{+} background events from W-->tau in pt: " << ipt << " eta: " << ieta << " cent: " << icent << " = " 
                    << nBkg << " +/- " << nBkg_err_upp << ", " << nBkg_err_low << std::endl;
                grName << "grA_muPlus_tau_pt" << ipt << "_cent" << icent;
                m_gr[grName.str()]->SetPoint(ieta,x,nBkg);
                m_gr[grName.str()]->SetPointError(ieta,xerr,xerr,nBkg_err_low,nBkg_err_upp);
                grName.str(std::string());

                // B
                calculateBkg(m_nFid,m_nGen,ipt,ieta,icent,"tau_pp","+",sigma_tau_plus_pp,LintB,nBkg_pp,nBkg_pp_err_upp, nBkg_pp_err_low);
                calculateBkg(m_nFid,m_nGen,ipt,ieta,icent,"tau_pn","+",sigma_tau_plus_pn,LintB,nBkg_pn,nBkg_pn_err_upp, nBkg_pn_err_low);
                nBkg = wt_pp*nBkg_pp + wt_pn*nBkg_pn;
                nBkg_err_upp = wt_pp*nBkg_pp_err_upp + wt_pn*nBkg_pn_err_upp;
                nBkg_err_low = wt_pp*nBkg_pp_err_low + wt_pn*nBkg_pn_err_low;
                
                std::cout << "Number of Period B mu^{+} background events from W-->tau in pt: " << ipt << " eta: " << ieta << " cent: " << icent << " = " 
                    << nBkg << " +/- " << nBkg_err_upp << ", " << nBkg_err_low << std::endl;
                grName << "grB_muPlus_tau_pt" << ipt << "_cent" << icent;
                m_gr[grName.str()]->SetPoint(ieta,x,nBkg);
                m_gr[grName.str()]->SetPointError(ieta,xerr,xerr,nBkg_err_low,nBkg_err_upp);
                grName.str(std::string());

                // mu^{-}
                // A
                calculateBkg(m_nFid,m_nGen,ipt,ieta,icent,"tau_pp","-",sigma_tau_minus_pp,LintA,nBkg_pp,nBkg_pp_err_upp, nBkg_pp_err_low);
                calculateBkg(m_nFid,m_nGen,ipt,ieta,icent,"tau_pn","-",sigma_tau_minus_pn,LintA,nBkg_pn,nBkg_pn_err_upp, nBkg_pn_err_low);
                nBkg = wt_pp*nBkg_pp + wt_pn*nBkg_pn;
                nBkg_err_upp = wt_pp*nBkg_pp_err_upp + wt_pn*nBkg_pn_err_upp;
                nBkg_err_low = wt_pp*nBkg_pp_err_low + wt_pn*nBkg_pn_err_low;
                
                std::cout << "Number of Period A mu^{-} background events from W-->tau in pt: " << ipt << " eta: " << ieta << " cent: " << icent << " = " 
                    << nBkg << " +/- " << nBkg_err_upp << ", " << nBkg_err_low << std::endl;
                grName << "grA_muMinus_tau_pt" << ipt << "_cent" << icent;
                m_gr[grName.str()]->SetPoint(ieta,x,nBkg);
                m_gr[grName.str()]->SetPointError(ieta,xerr,xerr,nBkg_err_low,nBkg_err_upp);
                grName.str(std::string());

                // B
                calculateBkg(m_nFid,m_nGen,ipt,ieta,icent,"tau_pp","-",sigma_tau_minus_pp,LintB,nBkg_pp,nBkg_pp_err_upp, nBkg_pp_err_low);
                calculateBkg(m_nFid,m_nGen,ipt,ieta,icent,"tau_pn","-",sigma_tau_minus_pn,LintB,nBkg_pn,nBkg_pn_err_upp, nBkg_pn_err_low);
                nBkg = wt_pp*nBkg_pp + wt_pn*nBkg_pn;
                nBkg_err_upp = wt_pp*nBkg_pp_err_upp + wt_pn*nBkg_pn_err_upp;
                nBkg_err_low = wt_pp*nBkg_pp_err_low + wt_pn*nBkg_pn_err_low;
                
                std::cout << "Number of Period B mu^{-} background events from W-->tau in pt: " << ipt << " eta: " << ieta << " cent: " << icent << " = " 
                    << nBkg << " +/- " << nBkg_err_upp << ", " << nBkg_err_low << std::endl;
                grName << "grB_muMinus_tau_pt" << ipt << "_cent" << icent;
                m_gr[grName.str()]->SetPoint(ieta,x,nBkg);
                m_gr[grName.str()]->SetPointError(ieta,xerr,xerr,nBkg_err_low,nBkg_err_upp);
                grName.str(std::string());

            }//ieta
            
            _fOut->cd();
            std::stringstream grName;
            grName << "grA_muPlus_Z_pt" << ipt << "_cent" << icent;
            m_gr[grName.str()]->Write( (grName.str()).c_str());
            delete m_gr[grName.str()];
            grName.str(std::string());

            grName << "grA_muMinus_Z_pt" << ipt << "_cent" << icent;
            m_gr[grName.str()]->Write( (grName.str()).c_str());
            delete m_gr[grName.str()];
            grName.str(std::string());

            grName << "grA_muPlus_tau_pt" << ipt << "_cent" << icent;
            m_gr[grName.str()]->Write( (grName.str()).c_str());
            delete m_gr[grName.str()];
            grName.str(std::string());

            grName << "grA_muMinus_tau_pt" << ipt << "_cent" << icent;
            m_gr[grName.str()]->Write( (grName.str()).c_str());
            delete m_gr[grName.str()];
            grName.str(std::string());

            grName << "grB_muPlus_Z_pt" << ipt << "_cent" << icent;
            m_gr[grName.str()]->Write( (grName.str()).c_str());
            delete m_gr[grName.str()];
            grName.str(std::string());

            grName << "grB_muMinus_Z_pt" << ipt << "_cent" << icent;
            m_gr[grName.str()]->Write( (grName.str()).c_str());
            delete m_gr[grName.str()];
            grName.str(std::string());

            grName << "grB_muPlus_tau_pt" << ipt << "_cent" << icent;
            m_gr[grName.str()]->Write( (grName.str()).c_str());
            delete m_gr[grName.str()];
            grName.str(std::string());

            grName << "grB_muMinus_tau_pt" << ipt << "_cent" << icent;
            m_gr[grName.str()]->Write( (grName.str()).c_str());
            delete m_gr[grName.str()];
            grName.str(std::string());
        }//icent
    }//ipt


}
