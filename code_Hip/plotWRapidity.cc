#include <map>
#include <vector>
#include <sstream>
#include <iostream>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TH1F.h"

int indexPt(double pt, std::vector<double>& m_vPt){
   int nPtBins = m_vPt.size()-1;
   int index = 0;
   for(int ipt=0; ipt<nPtBins; ++ipt){
        if(pt>=m_vPt.at(ipt)&&pt<m_vPt.at(ipt+1)) return index;
        else ++index;
   }//ipt
   std::cout <<"WARNING: Muon pt: " << pt << " is either outside bin limits or an error has occured." << std::endl;
   return -1;
}

int indexCent(double cent, std::vector<double>& m_vCent){
   int nCentBins = m_vCent.size()-1;
   int index = 0;
   for(int icent=0; icent<nCentBins; ++icent){
        if(cent>=m_vCent.at(icent)&&cent<m_vCent.at(icent+1)) return index;
        else ++index;
   }//icent
   std::cout <<"WARNING: Muon centrality: " << cent << " is either outside bin limits or an error has occured." << std::endl;
   return -1;
}

void plot(TString sFileName, TString charge, TString beamType, std::map<std::string,TH1F*>& m_Wrapidity, std::vector<double>& vPt, std::vector<double>& vEta, std::vector<double>& vCent){
    TChain* _tGen = new TChain("treeGen","treeGen");
    _tGen->Add(sFileName);
    std::cout << "Number of entries in generator tree: " << _tGen->GetEntries() << std::endl;

    int genPDG, genNeutrinoPDG, genMother;
    bool isGenPeriodA, isGenPeriodB;
    float genPt, genEta, genCentrality, genNeutrinoPt, genMt, genCharge,Wrapidity;

    _tGen->SetBranchAddress("muon_gen_pt",&genPt);
    _tGen->SetBranchAddress("muon_gen_eta",&genEta);
    _tGen->SetBranchAddress("centrality",&genCentrality);
    _tGen->SetBranchAddress("neutrino_gen_pt",&genNeutrinoPt);
    _tGen->SetBranchAddress("mt_gen_neutrino_muon",&genMt);
    _tGen->SetBranchAddress("muon_gen_charge",&genCharge);
    _tGen->SetBranchAddress("muon_pdg",&genPDG);
    _tGen->SetBranchAddress("neutrino_pdg",&genNeutrinoPDG);
    _tGen->SetBranchAddress("muon_mother",&genMother);
    _tGen->SetBranchAddress("W_rapidity",&Wrapidity);
    _tGen->SetBranchAddress("isPeriodA",&isGenPeriodA);
    _tGen->SetBranchAddress("isPeriodB",&isGenPeriodB);

    float ptLow=vPt[0],ptUpp=vPt[vPt.size()-1];
    float centLow=vCent[0],centUpp=vCent[vCent.size()-1];
    float etaLow=vEta[0],etaUpp=vEta[vEta.size()-1];
    float metCut = 25.0;
    float mtCut  = 40.0;

    // generator loop
    for(int igen=0; igen<_tGen->GetEntries(); ++igen){
        _tGen->GetEntry(igen);
        if( genPt > ptLow             &&
            genPt < ptUpp &&
            genEta>etaLow   &&
            genEta<etaUpp  &&
            genCentrality>centLow   &&
            genCentrality<centUpp  &&

            genNeutrinoPt > metCut       &&
            genMt > mtCut &&
            //safety
            fabs(genPDG) == 13 &&
            fabs(genMother) == 24
          ){
            int ipt=indexPt(genPt,vPt);
            //int ieta=indexEta(genEta);
            int icent=indexCent(genCentrality,vCent);

            std::stringstream bin;
            bin << "Wrapidity_pt" << ipt << "_eta0" << "_cent" << icent;
            //if(strcmp(charge,"+")==0) bin << "_pos";
            //else if(strcmp(charge,"-")==0) bin << "_neg";
            if(genCharge>0.) bin << "_pos";
            else if(genCharge<0.) bin << "_neg";
            bin << beamType;
            m_Wrapidity[bin.str()]->Fill(Wrapidity);
            
        }
    }//igen
}

void plotWRapidity(){


    TString sFileName = "/usatlas/u/tbales/scratch/Data2013Hip/";
    TFile* _fOut = new TFile("WrapidityHistos.root","recreate");

    std::vector<double> vPt;
    vPt.push_back(25.0);
    vPt.push_back(400.0);
    const int nPtBins = vPt.size()-1;
    std::vector<double> vEta;
    vEta.push_back(-10.0);
    vEta.push_back(+10.0);
    const int nEtaBins = vEta.size()-1;
    std::vector<double> vCent;
    vCent.push_back(0.00);
    vCent.push_back(0.05);
    vCent.push_back(0.10);
    vCent.push_back(0.20);
    vCent.push_back(0.40);
    vCent.push_back(0.90);
    const int nCentBins = vCent.size()-1;

    std::map<std::string,TH1F*> m_Wrapidity;
    TH1F* h=NULL;
    for(int ipt=0; ipt<nPtBins; ++ipt){
        for(int ieta=0; ieta<nEtaBins; ++ieta){
            for(int icent=0; icent<nCentBins; ++icent){

                std::stringstream bin;
                //Period A
                bin << "Wrapidity_pt" << ipt << "_eta0" << "_cent" << icent;
                h = new TH1F((bin.str()+"_pos_ppA").c_str(), "W-->#mu^{+} nu", 100, -10.0,10.0);
                m_Wrapidity.insert(std::make_pair<std::string,TH1F*>(bin.str()+"_pos_ppA",h));
                h = new TH1F((bin.str()+"_neg_ppA").c_str(), "W-->#mu^{-} nu~", 100, -10.0,10.0);
                m_Wrapidity.insert(std::make_pair<std::string,TH1F*>(bin.str()+"_neg_ppA",h));
                h = new TH1F((bin.str()+"_pos_npA").c_str(), "W-->#mu^{+} nu", 100, -10.0,10.0);
                m_Wrapidity.insert(std::make_pair<std::string,TH1F*>(bin.str()+"_pos_npA",h));
                h = new TH1F((bin.str()+"_neg_npA").c_str(), "W-->#mu^{-} nu~", 100, -10.0,10.0);
                m_Wrapidity.insert(std::make_pair<std::string,TH1F*>(bin.str()+"_neg_npA",h));
                //Period B
                h = new TH1F((bin.str()+"_pos_ppB").c_str(), "W-->#mu^{+} nu", 100, -10.0,10.0);
                m_Wrapidity.insert(std::make_pair<std::string,TH1F*>(bin.str()+"_pos_ppB",h));
                h = new TH1F((bin.str()+"_neg_ppB").c_str(), "W-->#mu^{-} nu~", 100, -10.0,10.0);
                m_Wrapidity.insert(std::make_pair<std::string,TH1F*>(bin.str()+"_neg_ppB",h));
                h = new TH1F((bin.str()+"_pos_pnB").c_str(), "W-->#mu^{+} nu", 100, -10.0,10.0);
                m_Wrapidity.insert(std::make_pair<std::string,TH1F*>(bin.str()+"_pos_pnB",h));
                h = new TH1F((bin.str()+"_neg_pnB").c_str(), "W-->#mu^{-} nu~", 100, -10.0,10.0);
                m_Wrapidity.insert(std::make_pair<std::string,TH1F*>(bin.str()+"_neg_pnB",h));
            }//icent
        }//ieta
    }//ipt

    plot(sFileName+"ntupleWminmunu_pp.PerA.PowhegPythia8_AU2CT10.root","-","_ppA",m_Wrapidity,vPt,vEta,vCent);
    plot(sFileName+"ntupleWminmunu_pp.PerB.PowhegPythia8_AU2CT10.root","-","_ppB",m_Wrapidity,vPt,vEta,vCent);
    plot(sFileName+"ntupleWminmunu_pn.PowhegPythia8_AU2CT10.root","-","_pnB",m_Wrapidity,vPt,vEta,vCent);
    plot(sFileName+"ntupleWminmunu_np.PowhegPythia8_AU2CT10.root","-","_npA",m_Wrapidity,vPt,vEta,vCent);

    plot(sFileName+"ntupleWplusmunu_pp.PerA.PowhegPythia8_AU2CT10.root","+","_ppA",m_Wrapidity,vPt,vEta,vCent);
    plot(sFileName+"ntupleWplusmunu_pp.PerB.PowhegPythia8_AU2CT10.root","+","_ppB",m_Wrapidity,vPt,vEta,vCent);
    plot(sFileName+"ntupleWplusmunu_pn.PowhegPythia8_AU2CT10.root","+","_pnB",m_Wrapidity,vPt,vEta,vCent);
    plot(sFileName+"ntupleWplusmunu_np.PowhegPythia8_AU2CT10.root","+","_npA",m_Wrapidity,vPt,vEta,vCent);

    std::map<std::string,TH1F*>::iterator it;
    _fOut->cd();
    for(it=m_Wrapidity.begin(); it!=m_Wrapidity.end(); ++it){
        if(it->second->GetEntries()>0) {
            it->second->Write() ;
            it->second->Print() ;
        }
    }

}
