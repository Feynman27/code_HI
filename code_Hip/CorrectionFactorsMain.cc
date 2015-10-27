#include "EfficiencyCorrection.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include <cmath>

void CalculatedWeightedEfficiencies(TGraphAsymmErrors* gr, double etaLow, double etaUpp, EfficiencyCorrection* ef1, EfficiencyCorrection* ef2, int ipt, int ieta, int icent, std::string charge){

    double weight_pp = 82.0/208.0;
    double weight_pn = (208.0-82.0)/208.0;
    double x = (etaLow+etaUpp)/2.0;
    double ex = fabs(etaLow-etaUpp)/2.0;


    double eff = weight_pp*ef1->GetEfficiency(ipt,ieta,icent,charge) + weight_pn*ef2->GetEfficiency(ipt,ieta,icent,charge);
    double eff_errLow = weight_pp*ef1->GetEfficiencyErrorLow(ipt,ieta,icent,charge) + weight_pn*ef2->GetEfficiencyErrorLow(ipt,ieta,icent,charge);
    double eff_errUp = weight_pp*ef1->GetEfficiencyErrorUp(ipt,ieta,icent,charge) + weight_pn*ef2->GetEfficiencyErrorUp(ipt,ieta,icent,charge);
    gr->SetPoint(ieta,x,eff);
    gr->SetPointError(ieta,ex,ex,eff_errLow,eff_errUp);
}

void CalculateEfficiencies(TString sFileName, TFile* _fOut, EfficiencyCorrection* ef, double ptLow, double ptUpp, double etaLow, double etaUpp, double centLow, double centUpp,bool doFcalReWtd=false){

    TChain* _tGen = new TChain("treeGen","treeGen");
    _tGen->Add(sFileName);
    TChain* _tRec = new TChain("treeW","treeW");
    _tRec->Add(sFileName);
    std::cout << "Number of entries in generator tree: " << _tGen->GetEntries() << std::endl;
    std::cout << "Number of entries in reconstructed tree: " << _tRec->GetEntries() << std::endl;


    float recPt, recEta, recCentrality, recMET, recMt, recCharge,recPtcone30;
    int recPDG, recMother, recMatched; 
    float genPt, genEta, genCentrality, genNeutrinoPt, genMt, genCharge;
    float fcalGen,fcalRec;
    int genPDG, genNeutrinoPDG, genMother;
    bool isGenPeriodA, isGenPeriodB;
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

    _tGen->SetBranchAddress("muon_gen_pt",&genPt);
    _tGen->SetBranchAddress("muon_gen_eta",&genEta);
    _tGen->SetBranchAddress("centrality",&genCentrality);
    _tGen->SetBranchAddress("neutrino_gen_pt",&genNeutrinoPt);
    _tGen->SetBranchAddress("mt_gen_neutrino_muon",&genMt);
    _tGen->SetBranchAddress("muon_gen_charge",&genCharge);
    _tGen->SetBranchAddress("muon_pdg",&genPDG);
    _tGen->SetBranchAddress("neutrino_pdg",&genNeutrinoPDG);
    _tGen->SetBranchAddress("muon_mother",&genMother);
    _tGen->SetBranchAddress("isPeriodA",&isGenPeriodA);
    _tGen->SetBranchAddress("isPeriodB",&isGenPeriodB);
    _tGen->SetBranchAddress("fcal",&fcalGen);


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
            if(doFcalReWtd) ef->fillHistos(genPt,genEta,genCentrality,genCharge,"tot",fcalGen);
            else ef->fillHistos(genPt,genEta,genCentrality,genCharge,"tot");
        }
    }//igen
    // reco loop
    float isoCut = 0.2;
    for(int irec=0; irec<_tRec->GetEntries(); ++irec){

        _tRec->GetEntry(irec);
        float isoRatio = recPtcone30/recPt;
        if( recPt > ptLow             &&
            recPt < ptUpp &&
            recEta>etaLow   &&
            recEta<etaUpp  &&
            recCentrality>centLow   &&
            recCentrality<centUpp  &&

            isoRatio <= isoCut    &&
            recMET > metCut       &&
            recMt > mtCut         &&
            recMatched==1         &&
            //safety
            fabs(recPDG) == 13    &&
            fabs(recMother) == 24 
            &&isZCandidate == 0     
            &&EF_mu8_matched == 1
          ){
            if(doFcalReWtd) ef->fillHistos(recPt,recEta,recCentrality,recCharge,"pass",fcalRec); 
            else ef->fillHistos(recPt,recEta,recCentrality,recCharge,"pass");
        }

    }//irec

    ef->CalculateEfficiency();
    ef->getEtaCentDep();
    ef->writeHistos(_fOut);

}

int main(){
    // Turn on to reweight Cw by fcal
    bool doFCalReWtd=false;

    TString sFileName = "/usatlas/u/tbales/scratch/data_files/Data2013Hip/";
    TFile* _fOut = new TFile("efficiencies.root","recreate");

    std::vector<double> vPt;
    vPt.push_back(25.0);
    vPt.push_back(400.0);
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


    //mu- PeriodA
    EfficiencyCorrection* ef_neg_pp_A = new EfficiencyCorrection(vPt,vEta,vCent,"pp","A","-");
    CalculateEfficiencies(sFileName+"ntupleWminmunu_pp.A.PowhegPythia8_AU2CT10.09.24.2014.root",
            _fOut,ef_neg_pp_A,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],doFCalReWtd);
    EfficiencyCorrection* ef_neg_np_A = new EfficiencyCorrection(vPt,vEta,vCent,"np","A","-");
    CalculateEfficiencies(sFileName+"ntupleWminmunu_pn.A.PowhegPythia8_AU2CT10.09.24.2014.root",
            _fOut,ef_neg_np_A,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],doFCalReWtd);
    
    // mu- PeriodB
    EfficiencyCorrection* ef_neg_pp_B = new EfficiencyCorrection(vPt,vEta,vCent,"pp","B","-");
    CalculateEfficiencies(sFileName+"ntupleWminmunu_pp.B.PowhegPythia8_AU2CT10.09.24.2014.root",
            _fOut,ef_neg_pp_B,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],doFCalReWtd);
    EfficiencyCorrection* ef_neg_pn_B = new EfficiencyCorrection(vPt,vEta,vCent,"pn","B","-");
    CalculateEfficiencies(sFileName+"ntupleWminmunu_pn.B.PowhegPythia8_AU2CT10.09.24.2014.root",
            _fOut,ef_neg_pn_B,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],doFCalReWtd);

    //mu+ PeriodA
    EfficiencyCorrection* ef_pos_pp_A = new EfficiencyCorrection(vPt,vEta,vCent,"pp","A","+");
    CalculateEfficiencies(sFileName+"ntupleWplusmunu_pp.A.PowhegPythia8_AU2CT10.09.24.2014.root",
            _fOut,ef_pos_pp_A,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],doFCalReWtd);
    EfficiencyCorrection* ef_pos_np_A = new EfficiencyCorrection(vPt,vEta,vCent,"np","A","+");
    CalculateEfficiencies(sFileName+"ntupleWplusmunu_pn.A.PowhegPythia8_AU2CT10.09.24.2014.root",
            _fOut,ef_pos_np_A,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],doFCalReWtd);
    
    // mu- PeriodB
    EfficiencyCorrection* ef_pos_pp_B = new EfficiencyCorrection(vPt,vEta,vCent,"pp","B","+");
    CalculateEfficiencies(sFileName+"ntupleWplusmunu_pp.B.PowhegPythia8_AU2CT10.09.24.2014.root",
            _fOut,ef_pos_pp_B,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],doFCalReWtd);
    EfficiencyCorrection* ef_pos_pn_B = new EfficiencyCorrection(vPt,vEta,vCent,"pn","B","+");
    CalculateEfficiencies(sFileName+"ntupleWplusmunu_pn.B.PowhegPythia8_AU2CT10.09.24.2014.root",
            _fOut,ef_pos_pn_B,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],doFCalReWtd);

    for(int ipt=0; ipt<ptBins; ++ipt){
        for(int icent=0; icent<centBins; ++icent){
            TGraphAsymmErrors* gr_NegA = new TGraphAsymmErrors(etaBins);
            TGraphAsymmErrors* gr_PosA = new TGraphAsymmErrors(etaBins);
            TGraphAsymmErrors* gr_NegB = new TGraphAsymmErrors(etaBins);
            TGraphAsymmErrors* gr_PosB = new TGraphAsymmErrors(etaBins);

           for(int ieta=0; ieta<etaBins; ++ieta){
             CalculatedWeightedEfficiencies(gr_NegA,vEta[ieta],vEta[ieta+1],ef_neg_pp_A,ef_neg_np_A,ipt,ieta,icent,"-");
             CalculatedWeightedEfficiencies(gr_NegB,vEta[ieta],vEta[ieta+1],ef_neg_pp_B,ef_neg_pn_B,ipt,ieta,icent,"-");
             CalculatedWeightedEfficiencies(gr_PosA,vEta[ieta],vEta[ieta+1],ef_pos_pp_A,ef_pos_np_A,ipt,ieta,icent,"+");
             CalculatedWeightedEfficiencies(gr_PosB,vEta[ieta],vEta[ieta+1],ef_pos_pp_B,ef_pos_pn_B,ipt,ieta,icent,"+");
         }//ieta
         std::stringstream nameNegA, nameNegB,namePosA,namePosB;
         nameNegA << "effNeg_A_wtd_pt" << ipt << "_cent" << icent;
         nameNegB << "effNeg_B_wtd_pt" << ipt << "_cent" << icent;
         gr_NegA->Write(nameNegA.str().c_str());
         gr_NegB->Write(nameNegB.str().c_str());
         namePosA << "effPos_A_wtd_pt" << ipt << "_cent" << icent;
         namePosB << "effPos_B_wtd_pt" << ipt << "_cent" << icent;
         gr_PosA->Write(namePosA.str().c_str());
         gr_PosB->Write(namePosB.str().c_str());
         delete gr_NegA;
         delete gr_PosA;
         delete gr_NegB;
         delete gr_PosB;
      }//icent
    }//ipt
}
