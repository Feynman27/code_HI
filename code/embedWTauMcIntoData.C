#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TTimer.h"

#include <vector>
#include <iostream>
#include <cmath>

#include "isGoodEvent.C"

#ifdef __MAKECINT__
#pragma link C++ class vector<vector<int> >+;
#endif

///Corrections to momentum based on independent study of
///MS by MCP group 
void getCorrection(double& delP1,double& delP2,double eta){

    ///Barrel
    if(fabs(eta)>0.&&fabs(eta)<1.05){delP1=1.80/100.;delP2=0.095e-3;}
    ///Transition
    else if(fabs(eta)>=1.05&&fabs(eta)<1.7){delP1=3.17/100.;delP2=0.250e-3;}
    ///End-cap
    else if(fabs(eta)>=1.7&&fabs(eta)<2.0){delP1=1.23/100.;delP2=0.169e-3;}
    ///CSC/No-TRT
    else if(fabs(eta)>=2.0&&fabs(eta)<2.5){delP1=0.52/100.;delP2=0.453e-3;}
    else {delP1=0./100.;delP2=0.e-3;}
}

double transformPt(double pt, double eta){

    double N1 = gRandom->Gaus(0,1);
    double delP1, delP2;
    getCorrection(delP1,delP2,eta);
    pt = pt*(1.0+N1*delP1+N1*delP2*pt);
    return pt;
}

void getParameters(double& p0,double& p1,double& p2,double eta){
     ///Barrel
    if(fabs(eta)>0.&&fabs(eta)<1.05){p0=0.25; p1=3.27/100.;p2=0.168e-3;}
    ///Transition
    else if(fabs(eta)>=1.05&&fabs(eta)<1.7){p0=0.;p1=6.49/100.;p2=0.336e-3;}
    ///End-cap
    else if(fabs(eta)>=1.7&&fabs(eta)<2.0){p0=0.;p1=3.79/100.;p2=0.196e-3;}
    ///CSC/No-TRT
    else if(fabs(eta)>=2.0&&fabs(eta)<2.5){p0=0.15;p1=2.82/100.;p2=0.469e-3;}
    else {p0=0.;p1=0./100.;p2=0.e-3;}
   
}

double smearMuonPt(double eta, double pt){
//    double theta = 2.0*TMath::ATan(TMath::Exp(-1.0*eta));
    //pt = transformPt(pt,eta);
    double p0,p1,p2;
    getParameters(p0,p1,p2,eta);
    double term1 = p0/pt;
    double term2 = p1;
    double term3 = p2*pt;
    double sigmaOverP = TMath::Sqrt(TMath::Power(term1,2)+TMath::Power(term2,2)+TMath::Power(term3,2));
    return sigmaOverP;
}

float AssignCentrality(float fcal_Et){

    float centralityTemp = -1.;
    if(fcal_Et > 2.8038){
        centralityTemp = 0.025; //0-5%
    }
    else if (fcal_Et > 2.3102){
        centralityTemp = 0.075 ; //5-10%
    }
    else if (fcal_Et > 1.9073){
        centralityTemp = 0.125 ; //10-15%
    }
    else if (fcal_Et > 1.5707){
        centralityTemp = 0.175 ; //15-20%
    }
    else if (fcal_Et > 1.0448){
        centralityTemp = 0.25 ; //20-30%
    }
    else if (fcal_Et > 0.6624){
        centralityTemp = 0.35 ; //30-40%
    }
    else if (fcal_Et > 0.3909){
        centralityTemp = 0.45 ;  //40-50%
    }
    else if (fcal_Et > 0.2118){
        centralityTemp = 0.55 ; //50-60%
    }
    else if (fcal_Et > 0.1024){
        centralityTemp = 0.65 ; //60-70%
    }
    else if (fcal_Et > 0.0438){
        centralityTemp = 0.75 ; //70-80%
    }
    else {
        centralityTemp = 0.89;  //80-98%
    }

    return centralityTemp;
}

void embedWTauMcIntoData(){

    TChain* _tTau = new TChain("truth","truth");
    _tTau->Add("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wtaumu/SingleMuonFilter17WtauPYTHIA.NTUP_TRUTH.07.23.2013.root");
//    _tTau->Add("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MonteCarloFiles/Wtaumu/user.tbalestr.mc11_2TeV.209002.Pythia_Wtaunu_1Lepton17.NTUP_TRUTH.07.24.2013.v01.130725000317/*");

/*    std::vector<int>* mc_pdgId = new std::vector<int>();
    std::vector<float>* mc_pt = new std::vector<float>();
    std::vector<float>* mc_eta = new std::vector<float>();
    std::vector<float>* mc_phi = new std::vector<float>();
    std::vector<float>*  mc_charge = new std::vector<float>();
*/
    int mc_n;
    std::vector<int>* mc_pdgId = 0;
    std::vector<float>* mc_pt = 0;
    std::vector<float>* mc_eta = 0;
    std::vector<float>* mc_phi = 0;
    std::vector<float>*  mc_charge = 0;
//    std::vector<float>   *mc_E = 0;
    std::vector<float>   *mc_m = 0;
//    std::vector<float>   *mc_py = 0;
//    std::vector<float>   *mc_pz = 0;
//    std::vector<float>   *mc_px = 0;
    std::vector<std::vector<int> > *mc_child_index = 0;
    std::vector<std::vector<int> > *mc_parent_index = 0;

    _tTau->SetBranchAddress("mc_n",&mc_n);
    _tTau->SetBranchAddress("mc_pt",&mc_pt);
    _tTau->SetBranchAddress("mc_eta",&mc_eta);
    _tTau->SetBranchAddress("mc_phi",&mc_phi);
    _tTau->SetBranchAddress("mc_pdgId",&mc_pdgId);
    _tTau->SetBranchAddress("mc_charge",&mc_charge);
//    _tTau->SetBranchAddress("mc_E",&mc_E);
    _tTau->SetBranchAddress("mc_m",&mc_m);
//    _tTau->SetBranchAddress("mc_py",&mc_py);
//    _tTau->SetBranchAddress("mc_pz",&mc_pz);
//    _tTau->SetBranchAddress("mc_px",&mc_px);
    _tTau->SetBranchAddress("mc_child_index",&mc_child_index);
    _tTau->SetBranchAddress("mc_parent_index",&mc_parent_index);

    _tTau->SetBranchStatus("*",0);
    _tTau->SetBranchStatus("mc_n",1);
    _tTau->SetBranchStatus("mc_pt",1);
    _tTau->SetBranchStatus("mc_eta",1);
    _tTau->SetBranchStatus("mc_phi",1);
    _tTau->SetBranchStatus("mc_pdgId",1);
    _tTau->SetBranchStatus("mc_charge",1);
//    _tTau->SetBranchStatus("mc_E",1);
    _tTau->SetBranchStatus("mc_m",1);
//    _tTau->SetBranchStatus("mc_py",1);
//    _tTau->SetBranchStatus("mc_pz",1);
//    _tTau->SetBranchStatus("mc_px",1);
    _tTau->SetBranchStatus("mc_child_index",1);
    _tTau->SetBranchStatus("mc_parent_index",1);

    TChain* _tMB = new TChain("HeavyIonD3PD","HeavyIonD3PD");
    _tMB->Add("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MinimumBiasFiles/data11_hi.00194193.physics_MinBias.merge.NTUP_HI.f424_m1059/*");
    _tMB->Add("/mnt/Lustre/cgrp/atlas_hi/tbalestri/MinimumBiasFiles/data11_hi.00193558.physics_MinBias.merge.NTUP_HI.f424_m1059/*");
//    _tMB->Add("/mnt/Lustre/cgrp/atlas_hi/zcitron/mb2011/data11_hi.00194193.physics_MinBias.merge.NTUP_HI.f424_m1059/data11_hi*");
//    _tMB->Add("/mnt/Lustre/cgrp/atlas_hi/zcitron/mb2011/data11_hi.00194370.physics_MinBias.merge.NTUP_HI.f424_m1059/data11_hi*");

    int vx_n,RunNumber,EventNumber,lbn,bcid,nTrk;
/*    std::vector<int>* vx_nTracks = new std::vector<int>();
    std::vector<float>* mbtime_timeA  = new std::vector<float>();
    std::vector<float>* mbtime_timeC = new std::vector<float>();
    std::vector<float>* cccEt_Et_Eh_FCal = new std::vector<float>();
    std::vector<float>* trk_theta_wrtPV = new std::vector<float>();
    std::vector<float>* trk_phi_wrtPV = new std::vector<float>();
    std::vector<float>* trk_qoverp_wrtPV = new std::vector<float>();
    std::vector<float>* trk_chi2 = new std::vector<float>();
    std::vector<int>* trk_ndof = new std::vector<int>();
    std::vector<float>* trk_z0_wrtPV = new std::vector<float>();
    std::vector<float>* trk_cov_z0_wrtPV = new std::vector<float>();
    std::vector<float>* trk_cov_theta_wrtPV = new std::vector<float>();
    std::vector<float>* trk_d0_wrtPV = new std::vector<float>();
    std::vector<float>* trk_cov_d0_wrtPV = new std::vector<float>();
    std::vector<float>* trk_eta = new std::vector<float>();

    std::vector<int>* trk_nPixHits = new std::vector<int>();
    std::vector<int>* trk_nSCTHits = new std::vector<int>();
    std::vector<int>* trk_nSCTHoles = new std::vector<int>();
    std::vector<int>* trk_nPixHoles = new std::vector<int>();
    std::vector<int>* trk_expectBLayerHit = new std::vector<int>();
    std::vector<int>* trk_nBLHits = new std::vector<int>();
*/
    std::vector<int>* vx_nTracks = 0;
    float mbtime_timeA , mbtime_timeC ;
    std::vector<float>* cccEt_Et_Eh_FCal = 0;
    std::vector<float>* trk_theta_wrtPV = 0;
    std::vector<float>* trk_phi_wrtPV = 0;
    std::vector<float>* trk_qoverp_wrtPV = 0;
    std::vector<float>* trk_chi2 = 0;
    std::vector<int>* trk_ndof = 0;
    std::vector<float>* trk_z0_wrtPV = 0;
    std::vector<float>* trk_cov_z0_wrtPV = 0;
    std::vector<float>* trk_cov_theta_wrtPV = 0;
    std::vector<float>* trk_d0_wrtPV = 0;
    std::vector<float>* trk_cov_d0_wrtPV = 0;
    std::vector<float>* trk_eta = 0;

    std::vector<int>* trk_nPixHits = 0;
    std::vector<int>* trk_nSCTHits = 0;
    std::vector<int>* trk_nSCTHoles = 0;
    std::vector<int>* trk_nPixHoles = 0;
    std::vector<int>* trk_expectBLayerHit = 0;
    std::vector<int>* trk_nBLHits = 0;


    _tMB->SetBranchAddress("vx_n",&vx_n);
    _tMB->SetBranchAddress("vx_nTracks",&vx_nTracks);
    _tMB->SetBranchAddress("mbtime_timeA",&mbtime_timeA);
    _tMB->SetBranchAddress("mbtime_timeC",&mbtime_timeC);
    _tMB->SetBranchAddress("cccEt_Et_Eh_FCal",&cccEt_Et_Eh_FCal);
    _tMB->SetBranchAddress("trk_theta_wrtPV",&trk_theta_wrtPV);
    _tMB->SetBranchAddress("trk_phi_wrtPV",&trk_phi_wrtPV);
    _tMB->SetBranchAddress("trk_qoverp_wrtPV",&trk_qoverp_wrtPV);
    _tMB->SetBranchAddress("trk_expectBLayerHit",&trk_expectBLayerHit);
    _tMB->SetBranchAddress("trk_nBLHits",&trk_nBLHits);
    _tMB->SetBranchAddress("trk_chi2",&trk_chi2);
    _tMB->SetBranchAddress("trk_ndof",&trk_ndof);
    _tMB->SetBranchAddress("trk_z0_wrtPV",&trk_z0_wrtPV);
    _tMB->SetBranchAddress("trk_cov_z0_wrtPV",&trk_cov_z0_wrtPV);
    _tMB->SetBranchAddress("trk_cov_theta_wrtPV",&trk_cov_theta_wrtPV);
    _tMB->SetBranchAddress("trk_d0_wrtPV",&trk_d0_wrtPV);
    _tMB->SetBranchAddress("trk_cov_d0_wrtPV",&trk_cov_d0_wrtPV);
    _tMB->SetBranchAddress("trk_nPixHits",&trk_nPixHits);
    _tMB->SetBranchAddress("trk_nSCTHits",&trk_nSCTHits);
    _tMB->SetBranchAddress("trk_eta",&trk_eta);
    _tMB->SetBranchAddress("trk_nSCTHoles",&trk_nSCTHoles);
    _tMB->SetBranchAddress("trk_nPixHoles",&trk_nPixHoles);
    _tMB->SetBranchAddress("RunNumber",&RunNumber);
    _tMB->SetBranchAddress("EventNumber",&EventNumber);
    _tMB->SetBranchAddress("lbn",&lbn);
    _tMB->SetBranchAddress("bcid",&bcid);
    _tMB->SetBranchAddress("trk_n",&nTrk);

    _tMB->SetBranchStatus("*",0);
    _tMB->SetBranchStatus("RunNumber",1);
    _tMB->SetBranchStatus("EventNumber",1);
    _tMB->SetBranchStatus("lbn",1);
    _tMB->SetBranchStatus("bcid",1);
    _tMB->SetBranchStatus("vx_n",1);
    _tMB->SetBranchStatus("trk_n",1);
    _tMB->SetBranchStatus("vx_nTracks",1);
    _tMB->SetBranchStatus("mbtime_timeA",1);
    _tMB->SetBranchStatus("mbtime_timeC",1);
    _tMB->SetBranchStatus("cccEt_Et_Eh_FCal",1);
    _tMB->SetBranchStatus("trk_theta_wrtPV",1);
    _tMB->SetBranchStatus("trk_phi_wrtPV",1);
    _tMB->SetBranchStatus("trk_qoverp_wrtPV",1);
    _tMB->SetBranchStatus("trk_expectBLayerHit",1);
    _tMB->SetBranchStatus("trk_nBLHits",1);
    _tMB->SetBranchStatus("trk_chi2",1);
    _tMB->SetBranchStatus("trk_ndof",1);
    _tMB->SetBranchStatus("trk_z0_wrtPV",1);
    _tMB->SetBranchStatus("trk_cov_z0_wrtPV",1);
    _tMB->SetBranchStatus("trk_cov_theta_wrtPV",1);
    _tMB->SetBranchStatus("trk_d0_wrtPV",1);
    _tMB->SetBranchStatus("trk_cov_d0_wrtPV",1);
    _tMB->SetBranchStatus("trk_nPixHits",1);
    _tMB->SetBranchStatus("trk_nSCTHits",1);
    _tMB->SetBranchStatus("trk_eta",1);
    _tMB->SetBranchStatus("trk_nSCTHoles",1);
    _tMB->SetBranchStatus("trk_nPixHoles",1);


    ///Output ntuple variables
    int nTauMu;
    float mc_mt[50], mc_mptPhi,mc_mptPx,mc_mptPy,mc_mptPz,mc_mptPt,centrality;
    float  myMc_pt[50],  myMc_ptNominal[50],myMc_eta[50],myMc_phi[50],myMc_pdgId[50],myMc_charge[50];
    float myMc_E[50],myMc_M[50];
    TFile* _outFile = new TFile("HIWtaumuNtuple.root","recreate");
    TTree* _tOut = new TTree("tree","tree");
    _tOut->Branch("nTauMu",&nTauMu,"nTauMu/I");
    _tOut->Branch("mc_mt",&mc_mt,"mc_mt[nTauMu]/F");
    _tOut->Branch("mc_mptPhi",&mc_mptPhi,"mc_mptPhi/F");
    _tOut->Branch("mc_mptPx",&mc_mptPx,"mc_mptPx/F");
    _tOut->Branch("mc_mptPy",&mc_mptPy,"mc_mptPy/F");
    _tOut->Branch("mc_mptPz",&mc_mptPz,"mc_mptPz/F");
    _tOut->Branch("mc_mptPt",&mc_mptPt,"mc_mptPt/F");
    _tOut->Branch("centrality",&centrality,"centrality/F");
    _tOut->Branch("mc_pt",&myMc_pt,"mc_pt[nTauMu]/F");
    _tOut->Branch("mc_ptGen",&myMc_ptNominal,"mc_ptGen[nTauMu]/F");
    _tOut->Branch("mc_M",&myMc_M,"mc_M[nTauMu]/F");
//    _tOut->Branch("mc_E",&myMc_E,"mc_E[nTauMu]/F");
    _tOut->Branch("mc_eta",&myMc_eta,"mc_eta[nTauMu]/F");
    _tOut->Branch("mc_phi",&myMc_phi,"mc_phi[nTauMu]/F");
    _tOut->Branch("mc_pdgId",&myMc_pdgId,"mc_pdgId[nTauMu]/F");
    _tOut->Branch("mc_charge",&myMc_charge,"mc_charge[nTauMu]/F");

    ///Lower trk threshold for mpt calculation
    double trkThresh = 3.0;

    ///Loop over events in MC
    int nMcEvents = _tTau->GetEntries();
    int nMbEvents = _tMB->GetEntries();
    std::cout << "Number of MC events: " << nMcEvents << std::endl;
    std::cout << "Number of MB events: " << nMbEvents << std::endl;

    int event = -9999, iMB=0; 
    ///Loop over all events in the MC
    for(int iMC=0; iMC<nMcEvents; ++iMC){

        _tTau->GetEntry(iMC);

        ///Number of W->tau->mu events
        nTauMu=0;
        int chargedHadrons=0;

        if(iMC%1000==0) {
            std::cout << "MC Event: " << iMC << std::endl;
//            if(iMC==10000) break;
        }
//        if(iMC==10000) break; //hack

        ///Loop over muons from tau decays
        int iMcMuCandidate=-9999;
        TVector3 vMPT;
        ///Initialize to 0 
        vMPT.SetPtEtaPhi(0.,0.,0.);

        ///Find the muon candidate from W->tau->mu in this MC event
//        std::cout << "Number of mc particles in event:" <<mc_n<< std::endl;
        for(int iTauMu=0; iTauMu<mc_n; ++iTauMu){
            
            ///only take muons within ID acceptance for MPT calculation
            //commented out for testing
            std::vector<int> vecParent = mc_parent_index->at(iTauMu);
            const int nParents = vecParent.size();
            ///make sure we found a parent
            if(nParents==0)continue;
            if(fabs(mc_pdgId->at(iTauMu))==13){
                int parent = mc_parent_index->at(iTauMu)[0];
                if(parent>=0&&parent<mc_n&&fabs(mc_pdgId->at(parent))==15){
                    std::vector<int> vecGrandParent = mc_parent_index->at(parent);
                    const int nGrandParents = vecGrandParent.size();
                    ///make sure there's a grandparent
                    if(nGrandParents==0) continue;
                    int grandParent = mc_parent_index->at(parent)[0];
                    if(grandParent>=0&&grandParent<mc_n&&fabs(mc_pdgId->at(grandParent))==24){
                        iMcMuCandidate = iTauMu;
                    }
                }
             }

             ///Include generated charged hadrons(pions,kaons,protons) in the MPT sum
             ///85% tracking efficiency
             //commented out for testing
             TRandom3 r; 
             double uniform = r.Uniform();
             if((fabs(mc_pdgId->at(iTauMu))==211||fabs(mc_pdgId->at(iTauMu))==321||fabs(mc_pdgId->at(iTauMu))==2212)
                 &&fabs(mc_eta->at(iTauMu))<2.5
                 && uniform<0.85
                 ){
                float ptTemp = mc_pt->at(iTauMu)/1000.;
                TVector3 vHad;
                vHad.SetPtEtaPhi(ptTemp,mc_eta->at(iTauMu),mc_phi->at(iTauMu));                
                vMPT+=vHad;
                ++chargedHadrons;
             }

/*            if( !(fabs(mc_pdgId->at(iTauMu))==13) ) continue;
//            std::cout << "muon index: " << iTauMu << std::endl;
//            std::cout << "muon pdg: " << mc_pdgId->at(iTauMu) << std::endl;

            int parent = mc_parent_index->at(iTauMu)[0];
            if(parent>=mc_n||parent<0) continue;
            if( !(fabs(mc_pdgId->at(parent))==15) ) continue;
//            std::cout << "tau index: " << parent << std::endl;
//            std::cout << "tau pdg: " << mc_pdgId->at(parent) << std::endl;

            int grandParent = mc_parent_index->at(parent)[0];
            if(grandParent>=mc_n||grandParent<0) continue;
            if( !(fabs(mc_pdgId->at(grandParent))==24) ) continue;
//            std::cout << "W index: " << grandParent << std::endl;
//            std::cout << "W pdg: " << mc_pdgId->at(grandParent) << std::endl;

            ///Once we've found it, save the mc index
            iMcMuCandidate = iTauMu; 
            break;
*/
            ///Only use muons
/*            if( ! (fabs(mc_pdgId->at(iTauMu))==13) ) continue;
            std::cout << "daughter index: " << iTauMu << std::endl;
            std::cout << "daughter pdg: " << mc_pdgId->at(iTauMu) << std::endl;
            
            ///index of muon parent
            int parent = mc_parent_index->at(iTauMu)[0];
            ///Sometimes muon parent index is very large???
            if(parent>=mc_n) continue;
            std::cout << "parent index: " << parent << std::endl;
//            std::cout << "parent pdg: " << mc_pdgId->at(parent) << std::endl;
            ///only use muons with a tau parent
            if( !(fabs(mc_pdgId->at(parent))==15) ) continue;
            std::cout << "parent index: " << parent << std::endl;
            ///only use taus from a W parent
            int grandParent = mc_parent_index->at(parent)[0];
            if( !(fabs(mc_pdgId->at(grandParent))==24) )continue;

            ///Once we've found it, save the mc index
            iMcMuCandidate = iTauMu; break;
*/            
        } ///iTauMu
        ///commented out for testing
        if(iMcMuCandidate<0) {
            std::cout <<"Muon candidate index: " << iMcMuCandidate << std::endl;
            std::cout << "WARNING:No W-->tau-->mu candiate found in this event. Going to next event. " << std::endl;
            continue;
        }

        ///Get muon candidate momentum and smear for use in event's MPT calculation
        float muPtTemp = mc_pt->at(iMcMuCandidate)/1000.0; ///GeV

        float nominalPt = muPtTemp;

        ///smear the generated pt
        ///return sigma/pT
        double smear = smearMuonPt(mc_eta->at(iMcMuCandidate),muPtTemp);
        if(muPtTemp<0.0) {
          std::cout << "WARNING: smearing results in pT<0 : " << muPtTemp << std::endl;
          std::cout << "Resetting to nominal pT and resmearing..." << std::endl;
        }

        ///if smearing results in negative value
        ///re-set and re-smear
        while(muPtTemp<0.0){

          ///reset to nominal pT
          muPtTemp = nominalPt;
          ///re-smear
          smear = smearMuonPt(mc_eta->at(iMcMuCandidate),muPtTemp);
          muPtTemp = muPtTemp+muPtTemp*smear*gRandom->Gaus();

          ///If still negative, repeat smearing
          if(muPtTemp>0.0) {
              std::cout << "Smeared pT is now positve:" << muPtTemp << std::endl;
          }
        }

        ///Fill global variables
        ///Fetch a random MB event that passes GRL and good event selection 
        if(nMbEvents<nMcEvents){
            std::cout << "ERROR: MB events must be >= MC events." << std::endl;
            exit(0);
        }

        ///Include smeared muon in mpt vector sum
        TVector3 vMu; 
        vMu.SetPtEtaPhi(muPtTemp,mc_eta->at(iMcMuCandidate),mc_phi->at(iMcMuCandidate));
        if(fabs(mc_eta->at(iMcMuCandidate))<2.5) vMPT+=vMu;

        ///Loop over minimum bias events from data
        bool eventPassed = false;
        TTimer* timer = new TTimer();;
        //commented out for testing
        while(!(eventPassed)){

            ///1 minute timeout
            timer->Start(60.e3);
            ///Start at event 0
            int e1 = iMB; 

            ///Initialize number of good tracks to 0
            int trkCounter = 0;

            ///Make sure we don't go above
            ///the number of events in the MB input file
            ///and in addition, make sure we don't use the
            ///same event twice
            if(e1>=nMbEvents){
                std::cout << "ERROR: Index out of MB event range. " << std::endl;
                exit(0);
            }
            if(e1<nMbEvents&&e1!=event){
                _tMB->GetEntry(e1);
                ///Good run and event
                bool isGoodEventFlag = false;
                isGoodEventFlag = isGoodEvent(RunNumber,lbn,bcid);
                if(
                    isGoodEventFlag==true && 
                    vx_n>1          &&
                    (vx_nTracks->at(0)>2) &&
                    ( (fabs(mbtime_timeA-mbtime_timeC) < 3.) )
                  ) { 
//                    std::cout << "MB event number: " << e1 << std::endl;

                    ///Once we've found a good event, assign the centrality to this MC event
                    float Fcal_Et=(cccEt_Et_Eh_FCal->at(0)+cccEt_Et_Eh_FCal->at(1)+cccEt_Et_Eh_FCal->at(2))/1000000.0;  //Fcal energy in TeV
                    ///Assign centrality to this event
                    centrality = AssignCentrality(Fcal_Et);

//                    std::cout << "Event centrality: " << centrality*100. << "[%]" << std::endl;

                    ///Loop over tracks in event and calculate MPT 
                    for(int itrk=0; itrk<nTrk; ++itrk){
                          
                        ///tight track filter
                        int trkQuality = 0;
                        double trkTheta = trk_theta_wrtPV->at(itrk);
                        double trkEta = -1.*TMath::Log(TMath::Tan(trkTheta/2.)) ;
                        double trkPhi = trk_phi_wrtPV->at(itrk);
                        double trkPt = fabs(1./trk_qoverp_wrtPV->at(itrk))/TMath::CosH(trkEta)/1000. ; ///in GeV

                        if(trk_expectBLayerHit -> at(itrk)==1 && trk_nBLHits -> at(itrk)==0) continue;
                        if(trkPt> 10.0 && TMath::Prob(trk_chi2->at(itrk),trk_ndof->at(itrk)) <= 0.01) continue;

                        float z0     = trk_z0_wrtPV->at(itrk)*sin(trk_theta_wrtPV->at(itrk));
                        float ez0    = sqrt(trk_cov_z0_wrtPV->at(itrk)   *pow(sin(trk_theta_wrtPV->at(itrk)),2) +
                                            trk_cov_theta_wrtPV->at(itrk)*pow(cos(trk_theta_wrtPV->at(itrk))
                                                *trk_z0_wrtPV->at(itrk),2));
                        float d0 = trk_d0_wrtPV->at(itrk);
                        float ed0    = sqrt(trk_cov_d0_wrtPV->at(itrk));

                        if(
                            trk_nPixHits -> at(itrk)     > 1    &&
                            trk_nSCTHits -> at(itrk)       > 5    &&
                            fabs(trk_eta -> at(itrk))<2.5      &&
                            trk_nSCTHoles   -> at(itrk) < 1 &&
                            (trk_nPixHoles-> at(itrk) + trk_nSCTHoles-> at(itrk) )   < 2    &&
                            (trk_nPixHits-> at(itrk)+trk_nSCTHits-> at(itrk) )   > 7    &&
                            fabs(d0/ed0) < 3.0 &&
                            fabs(z0/ez0) < 3.0 &&
                            trkPt>trkThresh
                         ) { 
                            ///Switch on good track locator flag
                            trkQuality = 1; 
                            ///Add 1 to the good-track counter
                            ++trkCounter;
                         }
                         else continue;

                         ///If at least 1 good track was in the event...
                         if(trkQuality==1){

                            ///Construct track 3 vector
                            TVector3 vTrk;
                            vTrk.SetPtEtaPhi(trkPt,trkEta,trkPhi);

                            ///Add it to the running sum of vectors for this event
                            vMPT += vTrk;
                         }


                    } //itrk

                    ///If no good tracks in the event, go 
                    ///to the next MB event
                    if(trkCounter==0){
                      eventPassed = false;
                      ++iMB;
                    }
                    else{
                      eventPassed = true; 
                      event=e1;
                      timer->Stop();
                    }
                    
                 } ///event cuts
                 ///Go to next MB event
                 else ++iMB;
                
            } //event constraints

            ///Go to next MB event
            else ++iMB;
        } ///events passed loop
        
        if(timer->HasTimedOut()==true) {
            timer->Timeout();
            exit(0);
        }

        ///Once we've looped over all the tracks in the event, 
        ///fill the missing momentum kinematics

        ///Flip the vector direction
        //commented out for testing
        mc_mptPhi = vMPT.Phi() + TMath::Pi();
        mc_mptPx = -1.0*vMPT.Px();
        mc_mptPy = -1.0*vMPT.Py();
        mc_mptPz = -1.0*vMPT.Pz();
        mc_mptPt = TMath::Sqrt(TMath::Power(mc_mptPx,2) + TMath::Power(mc_mptPy,2));

        ///Fill muon level variables
        //commented out for testing
        for(int iTauMu=0; iTauMu<mc_n; ++iTauMu){
            
            if( !(fabs(mc_pdgId->at(iTauMu))==13) ) continue;
//            std::cout << "muon index: " << iTauMu << std::endl;
//            std::cout << "muon pdg: " << mc_pdgId->at(iTauMu) << std::endl;
            std::vector<int> vecParent = mc_parent_index->at(iTauMu);
            const int nParents = vecParent.size();
            ///make sure we found a parent
            if(nParents==0)continue;

            int parent = mc_parent_index->at(iTauMu)[0];
            if(parent>=mc_n||parent<0) continue;
            if( !(fabs(mc_pdgId->at(parent))==15) ) continue;
//            std::cout << "tau index: " << parent << std::endl;
//            std::cout << "tau pdg: " << mc_pdgId->at(parent) << std::endl;

            std::vector<int> vecGrandParent = mc_parent_index->at(parent);
            const int nGrandParents = vecGrandParent.size();
            ///make sure there's a grandparent
            if(nGrandParents==0) continue;

            int grandParent = mc_parent_index->at(parent)[0];
            if(grandParent>=mc_n||grandParent<0) continue;
            if( !(fabs(mc_pdgId->at(grandParent))==24) ) continue;
//            std::cout << "W index: " << grandParent << std::endl;
//            std::cout << "W pdg: " << mc_pdgId->at(grandParent) << std::endl;

            
            myMc_pt[nTauMu] = mc_pt->at(iTauMu)/1000.0; ///GeV

            ///smear the generated pt
            nominalPt = myMc_pt[nTauMu];
            myMc_ptNominal[nTauMu] = nominalPt;
            smear = smearMuonPt(mc_eta->at(iTauMu),myMc_pt[nTauMu]);
            myMc_pt[nTauMu] = myMc_pt[nTauMu]+myMc_pt[nTauMu]*smear*gRandom->Gaus();

            if(myMc_pt[nTauMu]<0.0) {
                    std::cout << "WARNING: smearing results in pT<0 : " << myMc_pt[nTauMu] << std::endl;
                    std::cout << "Resetting to nominal pT and resmearing..." << std::endl;
            }

            ///if smearing results in negative value
            ///re-set and re-smear
            while(myMc_pt[nTauMu]<0.0){
                  ///reset to nominal pT
                  myMc_pt[nTauMu] = nominalPt;
                  ///re-smear
                  smear = smearMuonPt(mc_eta->at(iTauMu),myMc_pt[nTauMu]);
                  myMc_pt[nTauMu] = myMc_pt[nTauMu]+smear*gRandom->Gaus()*myMc_pt[nTauMu];

                  ///If still negative, repeat smearing
                  if(myMc_pt[nTauMu]>0.0) {
                      std::cout << "Smeared pT is now positve:" << myMc_pt[nTauMu] << std::endl;
                  }
             }


            myMc_eta[nTauMu] = mc_eta->at(iTauMu);
            myMc_phi[nTauMu] = mc_phi->at(iTauMu);
            myMc_pdgId[nTauMu] = mc_pdgId->at(iTauMu);
            myMc_charge[nTauMu] = mc_charge->at(iTauMu);
//            myMc_E[nTauMu] = mc_E->at(iTauMu);
            myMc_M[nTauMu] = mc_m->at(iTauMu)/1000.;

            double dphi = myMc_phi[nTauMu] - mc_mptPhi; 
            if(dphi<TMath::Pi()) dphi+=TMath::TwoPi();
            if(dphi>TMath::Pi()) dphi-=TMath::TwoPi();

            mc_mt[nTauMu] = (fabs(mc_mptPt) < 9999.) ? TMath::Sqrt( 2.0*myMc_pt[nTauMu]*mc_mptPt*(1.0-TMath::Cos(dphi)) ) :  -9999. ;
            

            ++nTauMu;
        } ///iTauMu (muon level)

        ///Fill tree once per MC event
         _tOut->Fill();
       } ///iMcEvents

       _outFile->cd();
       _outFile->Write();
       _tOut->Print();
       _outFile->Close();

       std::cout << "Clean up. " << std::endl;
//       delete _outFile;
//       delete _tOut;
/*       delete mc_pdgId;
       delete mc_pt;
       delete mc_eta;
       delete mc_phi;
       delete mc_charge;
       delete vx_nTracks;
       delete mbtime_timeA;
       delete mbtime_timeC;
       delete cccEt_Et_Eh_FCal;
       delete trk_theta_wrtPV;
       delete trk_phi_wrtPV;
       delete trk_qoverp_wrtPV;
       delete trk_chi2;
       delete trk_ndof;
       delete trk_z0_wrtPV;
       delete trk_cov_z0_wrtPV;
       delete trk_cov_theta_wrtPV;
       delete trk_d0_wrtPV;
       delete trk_cov_d0_wrtPV;
       delete trk_eta;
       delete trk_nPixHits;
       delete trk_nSCTHits;
       delete trk_nSCTHoles;
       delete trk_nPixHoles;
       delete trk_expectBLayerHit;
       delete trk_nBLHits;
       */
}
