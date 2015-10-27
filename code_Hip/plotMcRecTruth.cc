#include<map>
#include <iostream>
#include<sstream>
#include<vector>
#include<cmath>

#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraphErrors.h"

int getPtIndex(float& pt,const std::vector<double>& vPt){
    for(unsigned int i=0; i<vPt.size()-1; ++i){
        if(pt>vPt[i]&&pt<vPt[i+1]) return i;
    }
//    std::cout << "ERROR: Cannot pT find bin. " << std::endl; 
    return -1;
}
int getEtaIndex(float& eta,const std::vector<double>& vEta){
    for(unsigned int i=0; i<vEta.size()-1; ++i){
        if(eta>vEta[i]&&eta<vEta[i+1]) return i;
    }
    std::cout << "ERROR: Cannot eta find bin. " << std::endl; 
    return -1;

}
int getCentIndex(float& cent,const std::vector<double>& vCent){
    for(unsigned int i=0; i<vCent.size()-1; ++i){
        if(cent>vCent[i]&&cent<vCent[i+1]) return i;
    }
    //std::cout << "ERROR: Cannot cent find bin. " << std::endl; 
    return -1;

}

int main(){

    TFile* _fOut = new TFile("McRecResolutionHistos.root","recreate");
    std::map<std::string,TH1F*> mHistos;
    std::map<std::string,TH1F*> mPt;
    std::map<std::string,TH1F*> mEta;
    std::vector<TGraphErrors*> grPosPtMean;
    std::vector<TGraphErrors*> grPosPtSigma;
    std::vector<TGraphErrors*> grNegPtMean;
    std::vector<TGraphErrors*> grNegPtSigma;

    std::vector<double> vPt;
    vPt.push_back(0.0);
    vPt.push_back(10.0);
    vPt.push_back(20.0);
    vPt.push_back(25.0);
    vPt.push_back(30.0);
    vPt.push_back(35.0);
    vPt.push_back(40.0);
    vPt.push_back(50.0);
    vPt.push_back(60.0);
    vPt.push_back(200.0);
    const int ptBins = vPt.size()-1;
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
    const int etaBins = vEta.size()-1;
    std::vector<double> vCent;
    vCent.push_back(0.00);
    vCent.push_back(0.05);
    vCent.push_back(0.10);
    vCent.push_back(0.20);
    vCent.push_back(0.30); // added to match Z
    vCent.push_back(0.40);
    vCent.push_back(0.60); // added to match Z
    vCent.push_back(0.90);
    const int centBins = vCent.size()-1;

    for(int ipt=0; ipt<ptBins; ++ipt){
        for(int ieta=0; ieta<etaBins; ++ieta){
            for(int icent=0; icent<centBins; ++icent){
                std::stringstream ss;
                ss << "h_pos_pt" << ipt << "_eta" << ieta << "_cent" << icent;  
                TH1F* hPos = new TH1F( (ss.str()).c_str(), (ss.str()).c_str(), 80, -0.2, 0.2);
                mHistos.insert(std::make_pair(ss.str(),hPos));

                ss.str(std::string());
                ss << "h_neg_pt" << ipt << "_eta" << ieta << "_cent" << icent;  
                TH1F* hNeg = new TH1F( (ss.str()).c_str(), (ss.str()).c_str(), 80, -0.2, 0.2);
                mHistos.insert(std::make_pair(ss.str(),hNeg));
            }
        }//ieta
    }//ipt

    float centralityNt,etaNt,ptNt,chargeNt,ptcone30Nt,truth_eta,truth_pt;
    int mother,pdg;
	TChain* _tree = new TChain("treeW","treeW");
    _tree->Add("/usatlas/u/tbales/scratch/Data2013Hip/ntupleWplusmunu_pp.A.PowhegPythia8_AU2CT10.root");

    _tree->SetBranchAddress("centrality",&centralityNt);
    _tree->SetBranchAddress("eta",&etaNt);
    _tree->SetBranchAddress("pt",&ptNt);
    _tree->SetBranchAddress("charge",&chargeNt);
    _tree->SetBranchAddress("ptcone30",&ptcone30Nt);
    _tree->SetBranchAddress("rec_truth_eta",&truth_eta);
    _tree->SetBranchAddress("rec_truth_pt",&truth_pt);
    _tree->SetBranchAddress("rec_truth_mother_pdg",&mother);
    _tree->SetBranchAddress("rec_truth_pdg",&pdg);

    _tree->SetBranchStatus("*",0);
    _tree->SetBranchStatus("centrality",1);
    _tree->SetBranchStatus("eta",1);
    _tree->SetBranchStatus("pt",1);
    _tree->SetBranchStatus("charge",1);
    _tree->SetBranchStatus("ptcone30",1);
    _tree->SetBranchStatus("rec_truth_eta",1);
    _tree->SetBranchStatus("rec_truth_pt",1);
    _tree->SetBranchStatus("rec_truth_mother_pdg",1);
    _tree->SetBranchStatus("rec_truth_pdg",1);

    std::cout << "Number of entries in data tree: " << _tree->GetEntries() << std::endl;
    for(int i=0; i<_tree->GetEntries(); ++i){
       _tree->GetEntry(i);
       truth_pt = truth_pt*1.0e-3;
       float iso = ptcone30Nt/ptNt;
       // isolated muon from a W
       if( (iso>0.2) || (fabs(mother)!=24) || (fabs(pdg)!=13) || (centralityNt < vCent[0]) || (centralityNt > vCent[centBins]) ) continue;
       int ipt = getPtIndex(truth_pt,vPt);
       int ieta = getEtaIndex(truth_eta,vEta);
       int icent = getCentIndex(centralityNt,vCent);
       if(ipt<0||ieta<0||icent<0){
           std::cout << "WARNING: Pt: " << truth_pt << " Eta: " << truth_eta << " Cent: " << centralityNt << std::endl;
           continue;
           //break;
       }

       float relDiff = (ptNt-truth_pt)/truth_pt;
       /*std::cout << "ptNt: " << ptNt << " ptTruth: " << truth_pt << std::endl;
       std::cout << "relDiff: " << relDiff << std::endl;
       */
       std::stringstream hname;
       hname << "h_";
       if(chargeNt>0.) hname << "pos";
       else if(chargeNt<0.) hname << "neg";
       hname << "_pt" << ipt << "_eta" << ieta << "_cent" << icent;
       mHistos[hname.str()]->Fill(relDiff);
       hname.str(std::string());
    }//i

    // Differential in pT and cent, integrated over eta
    for(int icent=0; icent<centBins; ++icent){
        TGraphErrors* grPosMean = new TGraphErrors(ptBins);
        grPosPtMean.push_back(grPosMean);
        TGraphErrors* grNegMean = new TGraphErrors(ptBins);
        grNegPtMean.push_back(grNegMean);

        TGraphErrors* grPosSigma = new TGraphErrors(ptBins);
        grPosPtSigma.push_back(grPosSigma);
        TGraphErrors* grNegSigma = new TGraphErrors(ptBins);
        grNegPtSigma.push_back(grNegSigma);

        for(int ipt=0; ipt<ptBins; ++ipt){
            std::stringstream hnamepos,hnameneg;
            hnamepos << "h_pos_pt" << ipt << "_cent" << icent;  
            TH1F* hPos = new TH1F( (hnamepos.str()).c_str(), (hnamepos.str()).c_str(), 80, -0.2, 0.2);
            mPt.insert(std::make_pair(hnamepos.str(),hPos));

            hnameneg.str(std::string());
            hnameneg << "h_neg_pt" << ipt << "_cent" << icent;  
            TH1F* hNeg = new TH1F( (hnameneg.str()).c_str(), (hnameneg.str()).c_str(), 80, -0.2, 0.2);
            mPt.insert(std::make_pair(hnameneg.str(),hNeg));

            std::stringstream hname2;
            for(int ieta = 0; ieta<etaBins; ++ieta){
                hname2.str(std::string());
                hname2 << "h_pos_pt" << ipt << "_eta" << ieta << "_cent" << icent;
                mPt[hnamepos.str()]->Add(mHistos[hname2.str()]);
                hname2.str(std::string());
                hname2 << "h_neg_pt" << ipt << "_eta" << ieta << "_cent" << icent;
                mPt[hnameneg.str()]->Add(mHistos[hname2.str()]);
            }//ieta
            float ptBin = (vPt[ipt+1]+vPt[ipt])/2.0;
            float ptBinW = vPt[ipt+1]-(vPt[ipt+1]+vPt[ipt])/2.0;
            float mean,sigma;
            mean = mPt[hnamepos.str()]->GetMean();
            sigma = mPt[hnamepos.str()]->GetRMS();
            grPosPtMean[icent]->SetPoint(ipt,ptBin,mean);
            grPosPtMean[icent]->SetPointError(ipt,ptBinW,mPt[hnamepos.str()]->GetMeanError());
            grPosPtSigma[icent]->SetPoint(ipt,ptBin,sigma);
            grPosPtSigma[icent]->SetPointError(ipt,ptBinW,0.0);

            mean = mPt[hnameneg.str()]->GetMean();
            sigma = mPt[hnameneg.str()]->GetRMS();
            grNegPtMean[icent]->SetPoint(ipt,ptBin,mean);
            grNegPtMean[icent]->SetPointError(ipt,ptBinW,mPt[hnameneg.str()]->GetMeanError());
            grNegPtSigma[icent]->SetPoint(ipt,ptBin,sigma);
            grNegPtSigma[icent]->SetPointError(ipt,ptBinW,0.0);

        } //ipt
    }//icent

    std::map<std::string,TH1F*>::iterator it;
    _fOut->cd();
    for(it=mHistos.begin(); it!=mHistos.end(); ++it){
       //it->second->Write(); 
       delete it->second;
    }
    for(it=mPt.begin(); it!=mPt.end(); ++it){
       it->second->Write(); 
       delete it->second;
    }
    for(unsigned int icent=0; icent<grPosPtMean.size(); ++icent){
        std::stringstream ss;
        ss<<"hpos_mean_pt_cent"<<icent;
        grPosPtMean[icent]->Write( (ss.str()).c_str());
        delete grPosPtMean[icent];
        ss.str(std::string());
        ss<<"hpos_sigma_pt_cent"<<icent;
        grPosPtSigma[icent]->Write((ss.str()).c_str());
        delete grPosPtSigma[icent];

        ss.str(std::string());
        ss<<"hneg_mean_pt_cent"<<icent;
        grNegPtMean[icent]->Write( (ss.str()).c_str());
        delete grNegPtMean[icent];
        ss.str(std::string());
        ss<<"hneg_sigma_pt_cent"<<icent;
        grNegPtSigma[icent]->Write((ss.str()).c_str());
        delete grNegPtSigma[icent];
    }

}
