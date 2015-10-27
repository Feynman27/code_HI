#include "TagAndProbe.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include <cmath>
#include <string>
#include<ctime>

void CalculatedWeightedRecoEfficiencies(TFile* _fOut,std::string ds, std::string period, std::string charge,
                                            TagAndProbe* tp_id_pp, TagAndProbe* tp_id_pn, TagAndProbe* tp_ms_pp, TagAndProbe* tp_ms_pn){

    double weight_pp = 82.0/208.0;
    double weight_pn = (208.0-82.0)/208.0;

    std::string ds_per = (ds+"_"+period);

    std::vector<double> vecPt = tp_id_pp->getPtBins();
    const int nPtBins=vecPt.size()-1;
    std::vector<double> vecEta = tp_id_pp->getEtaBins();
    const int nEtaBins=vecEta.size()-1;
    std::vector<double> vecCent = tp_id_pp->getCentBins();
    const int nCentBins=vecCent.size()-1;
    if(nPtBins==0||nEtaBins==0||nCentBins==0){
        std::cout << "Error: Empty vector." << std::endl;
        return;
    }

    TGraphAsymmErrors* gr_eff_wtd = new TGraphAsymmErrors(nEtaBins);

    for(int ipt=0; ipt<nPtBins; ++ipt){
        for(int icent=0; icent<nCentBins; ++icent){
            std::stringstream effname;
            if(strcmp(charge.c_str(),"+")==0) effname << ds_per << "_Wtd_hPosEff_pt" << ipt << "_cent" << icent; 
            else if(strcmp(charge.c_str(),"-")==0) effname << ds_per << "_Wtd_hNegEff_pt" << ipt << "_cent" << icent; 
            else if(strcmp(charge.c_str(),"+-")==0) effname << ds_per << "_Wtd_hEff_pt" << ipt << "_cent" << icent; 
            else {
                std::cout << "An error has occured." << std::endl;
                return;
            }
            gr_eff_wtd->SetName( (effname.str()).c_str());
            std::cout << "Total weighted efficiency graph name: " << gr_eff_wtd->GetName() << std::endl;

            for(int ieta=0; ieta<nEtaBins; ++ieta){

                // Calculate eff and errors for ID/MS
                // pp
                double eff_id_pp=tp_id_pp->GetEfficiency(ipt,ieta,icent,charge);
                double deff_id_upp_pp= tp_id_pp->GetEfficiencyErrorUp(ipt,ieta,icent,charge);
                double deff_id_low_pp= tp_id_pp->GetEfficiencyErrorLow(ipt,ieta,icent,charge);

                double eff_ms_pp=tp_ms_pp->GetEfficiency(ipt,ieta,icent,charge);
                double deff_ms_upp_pp= tp_ms_pp->GetEfficiencyErrorUp(ipt,ieta,icent,charge);
                double deff_ms_low_pp= tp_ms_pp->GetEfficiencyErrorLow(ipt,ieta,icent,charge);

                // Combined ID/MS
                double eff_pp=0,deff_upp_pp=0,deff_low_pp=0;
                eff_pp = eff_id_pp*eff_ms_pp;
                deff_upp_pp = eff_pp*TMath::Sqrt(TMath::Power(deff_id_upp_pp/eff_id_pp,2) + TMath::Power(deff_ms_upp_pp/eff_ms_pp,2));
                deff_low_pp = eff_pp*TMath::Sqrt(TMath::Power(deff_id_low_pp/eff_id_pp,2) + TMath::Power(deff_ms_low_pp/eff_ms_pp,2));

                // pn
                double eff_id_pn=tp_id_pn->GetEfficiency(ipt,ieta,icent,charge);
                double deff_id_upp_pn= tp_id_pn->GetEfficiencyErrorUp(ipt,ieta,icent,charge);
                double deff_id_low_pn= tp_id_pn->GetEfficiencyErrorLow(ipt,ieta,icent,charge);

                double eff_ms_pn=tp_ms_pn->GetEfficiency(ipt,ieta,icent,charge);
                double deff_ms_upp_pn= tp_ms_pn->GetEfficiencyErrorUp(ipt,ieta,icent,charge);
                double deff_ms_low_pn= tp_ms_pn->GetEfficiencyErrorLow(ipt,ieta,icent,charge);

                // Combined ID/MS
                double eff_pn=0,deff_upp_pn=0,deff_low_pn=0;
                eff_pn = eff_id_pn*eff_ms_pn;
                deff_upp_pn = eff_pn*TMath::Sqrt(TMath::Power(deff_id_upp_pn/eff_id_pn,2) + TMath::Power(deff_ms_upp_pn/eff_ms_pn,2));
                deff_low_pn = eff_pn*TMath::Sqrt(TMath::Power(deff_id_low_pn/eff_id_pn,2) + TMath::Power(deff_ms_low_pn/eff_ms_pn,2));

                // current eta bin
                double eta = (vecEta[ieta+1]+vecEta[ieta])/2.0;
                double deta = fabs(vecEta[ieta+1]-vecEta[ieta])/2.0;

                // weighted eff
                double eff_wtd = (weight_pp*eff_pp+weight_pn*eff_pn)/(weight_pp+weight_pn); 
                double deff_upp_wtd = (weight_pp*deff_upp_pp+weight_pn*deff_upp_pn)/(weight_pp+weight_pn);
                double deff_low_wtd = (weight_pp*deff_low_pp+weight_pn*deff_low_pn)/(weight_pp+weight_pn);

                gr_eff_wtd->SetPoint(ieta,eta,eff_wtd);
                gr_eff_wtd->SetPointEYhigh(ieta,deff_upp_wtd);
                gr_eff_wtd->SetPointEYlow(ieta,deff_low_wtd);
                gr_eff_wtd->SetPointEXhigh(ieta,deta);
                gr_eff_wtd->SetPointEXlow(ieta,deta);
            
            }//ieta

            _fOut->cd();
            gr_eff_wtd->Write();
        }//icent
    }//ipt

    delete gr_eff_wtd;
}

void CalculateRecEffIdMs(TFile* _fOut, std::string ds, std::string period, std::string charge, TagAndProbe* tp_id, TagAndProbe* tp_ms){

    std::string ds_per = (ds+"_"+period);

    std::vector<double> vecPt = tp_id->getPtBins();
    const int nPtBins=vecPt.size()-1;
    std::vector<double> vecEta = tp_id->getEtaBins();
    const int nEtaBins=vecEta.size()-1;
    std::vector<double> vecCent = tp_id->getCentBins();
    const int nCentBins=vecCent.size()-1;
    if(nPtBins==0||nEtaBins==0||nCentBins==0){
        std::cout << "Error: Empty vector." << std::endl;
        return;
    }

    TEfficiency* eff_id=NULL,*eff_ms=NULL;
    TGraphAsymmErrors* gr_eff = new TGraphAsymmErrors(nEtaBins);
    for(int ipt=0; ipt<nPtBins; ++ipt){
        for(int icent=0; icent<nCentBins; ++icent){
            eff_id = tp_id->getCopyEfficiencyHisto(ipt,icent,charge);
            eff_ms = tp_ms->getCopyEfficiencyHisto(ipt,icent,charge);
            std::stringstream effname;
            if(strcmp(charge.c_str(),"+")==0) effname << ds_per << "_Tot_hPosEff_pt" << ipt << "_cent" << icent; 
            else if(strcmp(charge.c_str(),"-")==0) effname << ds_per << "_Tot_hNegEff_pt" << ipt << "_cent" << icent; 
            else if(strcmp(charge.c_str(),"+-")==0) effname << ds_per << "_Tot_hEff_pt" << ipt << "_cent" << icent; 
            else {
                std::cout << "An error has occured." << std::endl;
                return;
            }
            gr_eff->SetName( (effname.str()).c_str());
            std::cout << "Total efficiency graph name: " << gr_eff->GetName() << std::endl;

            for(int ieta=0; ieta<nEtaBins; ++ieta){

                // Calculate eff and errors for ID/MS
                double eff1=eff_id->GetEfficiency(ieta+1);
                double deffupp1 = eff_id->GetEfficiencyErrorUp(ieta+1);
                double defflow1 = eff_id->GetEfficiencyErrorLow(ieta+1);
                double eff2=eff_ms->GetEfficiency(ieta+1);
                double deffupp2 = eff_ms->GetEfficiencyErrorUp(ieta+1);
                double defflow2 = eff_ms->GetEfficiencyErrorLow(ieta+1);

                // Combined ID/MS
                double eff=0,deffupp=0,defflow=0;
                eff = eff1*eff2;
                deffupp = eff*TMath::Sqrt(TMath::Power(deffupp1/eff1,2) + TMath::Power(deffupp2/eff2,2));
                defflow = eff*TMath::Sqrt(TMath::Power(defflow1/eff1,2) + TMath::Power(defflow2/eff2,2));

                double eta = (vecEta[ieta+1]+vecEta[ieta])/2.0;
                double deta = fabs(vecEta[ieta+1]-vecEta[ieta])/2.0;
                gr_eff->SetPoint(ieta,eta,eff);
                gr_eff->SetPointEYhigh(ieta,deffupp);
                gr_eff->SetPointEYlow(ieta,defflow);
                gr_eff->SetPointEXhigh(ieta,deta);
                gr_eff->SetPointEXlow(ieta,deta);
            } //ieta

            // Write PtxCent graphs
            // as fcn of eta
            _fOut->cd();
            gr_eff->Write();
        }//icent
    }//ipt

    delete gr_eff;

}

void CalculateEfficiencies(TString sFileName, TString tName, TFile* _fOut, TagAndProbe* tp, double ptLow, double ptUpp, double etaLow, double etaUpp, double centLow, double centUpp,
                                bool doIdTP=false, bool doMsTP=false, bool doTrigTP=false, bool doMc=false){

    TChain* _tProbe = new TChain(tName,tName);
    _tProbe->Add(sFileName);
    std::cout << "Number of entries in tree: " << _tProbe->GetEntries() << std::endl;

    float centrality,charge,eta,pt;
    int probes,isMatched;
    bool isSS, isOS;
    bool isPeriodA,isPeriodB;

    // Branches
    _tProbe->SetBranchAddress("centrality",&centrality);
    _tProbe->SetBranchAddress("isPeriodA",&isPeriodA);
    _tProbe->SetBranchAddress("isPeriodB",&isPeriodB);
    if(doIdTP){
        _tProbe->SetBranchAddress("probeID_n",&probes);
        _tProbe->SetBranchAddress("rec_IDprobe_charge",&charge);
        _tProbe->SetBranchAddress("rec_IDprobe_eta",&eta);
        _tProbe->SetBranchAddress("rec_IDprobe_pt",&pt);
        _tProbe->SetBranchAddress("rec_IDprobe_isMatched",&isMatched);
        /*_tProbe->SetBranchAddress("rec_IDprobe_isSS",&isSS);
        _tProbe->SetBranchAddress("rec_IDprobe_isOS",&isOS);
        */
    }
    else if(doMsTP){
        _tProbe->SetBranchAddress("probeMS_n",&probes);
        _tProbe->SetBranchAddress("rec_MSprobe_charge",&charge);
        _tProbe->SetBranchAddress("rec_MSprobe_eta",&eta);
        _tProbe->SetBranchAddress("rec_MSprobe_pt",&pt);
        _tProbe->SetBranchAddress("rec_MSprobe_isMatched",&isMatched);
        /*_tProbe->SetBranchAddress("rec_MSprobe_isSS",&isSS);
        _tProbe->SetBranchAddress("rec_MSprobe_isOS",&isOS);
        */
    }
    else if(doTrigTP){
        _tProbe->SetBranchAddress("probeTrig_n",&probes);
        _tProbe->SetBranchAddress("trig_EF_mu8_probe_charge",&charge);
        _tProbe->SetBranchAddress("trig_EF_mu8_probe_eta",&eta);
        _tProbe->SetBranchAddress("trig_EF_mu8_probe_pt",&pt);
        _tProbe->SetBranchAddress("trig_EF_mu8_probe_isMatched",&isMatched);
        /*_tProbe->SetBranchAddress("trig_EF_mu8_probe_isSS",&isSS);
        _tProbe->SetBranchAddress("trig_EF_mu8_probe_isOS",&isOS);
        */
    }
    else {
        std::cout << "Please choose whether doing TP for the ID, MS, or trigger." << std::endl;
        return;
    }
    for(int ipr=0; ipr<_tProbe->GetEntries(); ++ipr){

        _tProbe->GetEntry(ipr);

        if( pt > ptLow             &&
            pt < ptUpp &&
            eta>etaLow   &&
            eta<etaUpp  &&
            centrality>centLow   &&
            centrality<centUpp
            //&& isOS==1
          ){
            tp->fillHistos(pt,eta,centrality,charge,"tot");
            if(isMatched==1) tp->fillHistos(pt,eta,centrality,charge,"pass");
        }
    }//ipr
    tp->CalculateEfficiency();
    tp->getEtaCentDep();
    tp->writeHistos(_fOut);

}

int main(){
    std::clock_t start;
    start=std::clock();

    bool doPeriodA=true;
    bool doPeriodB=true;
    bool doMc = false;
    bool doData = true;
    //TString sFileName = "/usatlas/u/tbales/HipMuonSkimming/WBosonFinder/test.root";
    TString sFileNameA = "/usatlas/u/tbales/scratch/Data2013Hip/ntupleW.periodA.TP.08.25.2014.root";
    TString sFileNameB = "/usatlas/u/tbales/scratch/Data2013Hip/ntupleW.periodB.TP.08.25.2014.root";
    /*TString sFileName_pp_A = "/usatlas/u/tbales/scratch/Data2013Hip/ntupleZ_pp.A.09.08.2014.root";
    TString sFileName_pp_B = "/usatlas/u/tbales/scratch/Data2013Hip/ntupleZ_pp.B.09.08.2014.root";
    TString sFileName_np_A = "/usatlas/u/tbales/scratch/Data2013Hip/ntupleZ_np.09.08.2014.root";
    TString sFileName_pn_B = "/usatlas/u/tbales/scratch/Data2013Hip/ntupleZ_pn.09.08.2014.root";
    */
    // new TP probe definition
    TString sFileName_pp_A = "/usatlas/u/tbales/scratch/Data2013Hip/ntupleZ_pp.A.09.18.2014.root";
    TString sFileName_pp_B = "/usatlas/u/tbales/scratch/Data2013Hip/ntupleZ_pp.B.09.18.2014.root";
    TString sFileName_np_A = "/usatlas/u/tbales/scratch/Data2013Hip/ntupleZ_np.09.18.2014.root";
    TString sFileName_pn_B = "/usatlas/u/tbales/scratch/Data2013Hip/ntupleZ_pn.09.18.2014.root";

    TFile* _fOut = new TFile("tagAndProbe.root","recreate");

    std::vector<double> vPt;
    vPt.push_back(25.0);
    vPt.push_back(100.0);
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
    vCent.push_back(0.40);
    vCent.push_back(0.90);
    const int centBins = vCent.size()-1;

    TString tIDname="tIDprobes",tMSname="tMSprobes",tTrigName="tTrigprobes";
    if(doPeriodA){
        ///////////
        // Period A
        ///////////
        if(doData){
            // T&P ID
            TagAndProbe* tp_recID_A = new TagAndProbe(vPt,vEta,vCent,"data","A","+-","ID");
            CalculateEfficiencies(sFileNameA,tIDname,_fOut,tp_recID_A,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],true,false,false);

            // T&P MS 
            TagAndProbe* tp_recMS_A = new TagAndProbe(vPt,vEta,vCent,"data","A","+-","MS");
            CalculateEfficiencies(sFileNameA,tMSname,_fOut,tp_recMS_A,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],false,true,false);

            // T&P MS/ID
            CalculateRecEffIdMs(_fOut,"data","A","+-",tp_recID_A,tp_recMS_A);

            // T&P trigger 
            TagAndProbe* tp_trig_A = new TagAndProbe(vPt,vEta,vCent,"data","A","+-","EF_mu8");
            CalculateEfficiencies(sFileNameA,tTrigName,_fOut,tp_trig_A,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],false,false,true);

            delete tp_recID_A;
            delete tp_recMS_A;
            delete tp_trig_A;
        }
        if(doMc){
            // pp
            // T&P ID
            TagAndProbe* tp_recID_pp_A = new TagAndProbe(vPt,vEta,vCent,"pp","A","+-","ID");
            CalculateEfficiencies(sFileName_pp_A,tIDname,_fOut,tp_recID_pp_A,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],true,false,false);

            // T&P MS 
            TagAndProbe* tp_recMS_pp_A = new TagAndProbe(vPt,vEta,vCent,"pp","A","+-","MS");
            CalculateEfficiencies(sFileName_pp_A,tMSname,_fOut,tp_recMS_pp_A,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],false,true,false);

            // T&P MS/ID
            CalculateRecEffIdMs(_fOut,"pp","A","+-",tp_recID_pp_A,tp_recMS_pp_A);

            // T&P trigger 
            TagAndProbe* tp_trig_pp_A = new TagAndProbe(vPt,vEta,vCent,"pp","A","+-","EF_mu8");
            CalculateEfficiencies(sFileName_pp_A,tTrigName,_fOut,tp_trig_pp_A,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],false,false,true);

            // np
            // T&P ID
            TagAndProbe* tp_recID_np_A = new TagAndProbe(vPt,vEta,vCent,"np","A","+-","ID");
            CalculateEfficiencies(sFileName_np_A,tIDname,_fOut,tp_recID_np_A,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],true,false,false);

            // T&P MS 
            TagAndProbe* tp_recMS_np_A = new TagAndProbe(vPt,vEta,vCent,"np","A","+-","MS");
            CalculateEfficiencies(sFileName_np_A,tMSname,_fOut,tp_recMS_np_A,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],false,true,false);

            // T&P MS/ID
            CalculateRecEffIdMs(_fOut,"np","A","+-",tp_recID_np_A,tp_recMS_np_A);

            // T&P trigger 
            TagAndProbe* tp_trig_np_A = new TagAndProbe(vPt,vEta,vCent,"np","A","+-","EF_mu8");
            CalculateEfficiencies(sFileName_np_A,tTrigName,_fOut,tp_trig_np_A,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],false,false,true);

            // pPb Reco
            CalculatedWeightedRecoEfficiencies(_fOut,"Pbp","A","+",tp_recID_pp_A,tp_recID_np_A,tp_recMS_pp_A,tp_recMS_np_A);
            CalculatedWeightedRecoEfficiencies(_fOut,"Pbp","A","-",tp_recID_pp_A,tp_recID_np_A,tp_recMS_pp_A,tp_recMS_np_A);
            CalculatedWeightedRecoEfficiencies(_fOut,"Pbp","A","+-",tp_recID_pp_A,tp_recID_np_A,tp_recMS_pp_A,tp_recMS_np_A);

            delete tp_recID_np_A;
            delete tp_recID_pp_A;
            delete tp_recMS_np_A;
            delete tp_recMS_pp_A;
            delete tp_trig_np_A;
            delete tp_trig_pp_A;
            
        }
    }

    if(doPeriodB){
        if(doData){
            ///////////
            // Period B
            ///////////
            // T&P ID
            TagAndProbe* tp_recID_B = new TagAndProbe(vPt,vEta,vCent,"data","B","+-","ID");
            CalculateEfficiencies(sFileNameB,tIDname,_fOut,tp_recID_B,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],true,false,false);

            // T&P MS 
            TagAndProbe* tp_recMS_B = new TagAndProbe(vPt,vEta,vCent,"data","B","+-","MS");
            CalculateEfficiencies(sFileNameB,tMSname,_fOut,tp_recMS_B,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],false,true,false);

            // T&P MS/ID
            CalculateRecEffIdMs(_fOut,"data","B","+-",tp_recID_B,tp_recMS_B);

            // T&P trigger 
            TagAndProbe* tp_trig_B = new TagAndProbe(vPt,vEta,vCent,"data","B","+-","EF_mu8");
            CalculateEfficiencies(sFileNameB,tTrigName,_fOut,tp_trig_B,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],false,false,true);

            delete tp_recID_B;
            delete tp_recMS_B;
            delete tp_trig_B;
        }
        if(doMc){
            // pp
            // T&P ID
            TagAndProbe* tp_recID_pp_B = new TagAndProbe(vPt,vEta,vCent,"pp","B","+-","ID");
            CalculateEfficiencies(sFileName_pp_B,tIDname,_fOut,tp_recID_pp_B,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],true,false,false);

            // T&P MS 
            TagAndProbe* tp_recMS_pp_B = new TagAndProbe(vPt,vEta,vCent,"pp","B","+-","MS");
            CalculateEfficiencies(sFileName_pp_B,tMSname,_fOut,tp_recMS_pp_B,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],false,true,false);

            // T&P MS/ID
            CalculateRecEffIdMs(_fOut,"pp","B","+-",tp_recID_pp_B,tp_recMS_pp_B);

            // T&P trigger 
            TagAndProbe* tp_trig_pp_B = new TagAndProbe(vPt,vEta,vCent,"pp","B","+-","EF_mu8");
            CalculateEfficiencies(sFileName_pp_B,tTrigName,_fOut,tp_trig_pp_B,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],false,false,true);

            // pn
            // T&P ID
            TagAndProbe* tp_recID_pn_B = new TagAndProbe(vPt,vEta,vCent,"pn","B","+-","ID");
            CalculateEfficiencies(sFileName_pn_B,tIDname,_fOut,tp_recID_pn_B,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],true,false,false);

            // T&P MS 
            TagAndProbe* tp_recMS_pn_B = new TagAndProbe(vPt,vEta,vCent,"pn","B","+-","MS");
            CalculateEfficiencies(sFileName_pn_B,tMSname,_fOut,tp_recMS_pn_B,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],false,true,false);

            // T&P MS/ID
            CalculateRecEffIdMs(_fOut,"pn","B","+-",tp_recID_pn_B,tp_recMS_pn_B);

            // T&P trigger 
            TagAndProbe* tp_trig_pn_B = new TagAndProbe(vPt,vEta,vCent,"pn","B","+-","EF_mu8");
            CalculateEfficiencies(sFileName_pn_B,tTrigName,_fOut,tp_trig_pn_B,vPt[0],vPt[vPt.size()-1],vEta[0],vEta[vEta.size()-1],vCent[0],vCent[vCent.size()-1],false,false,true);

            // pPb Reco
            CalculatedWeightedRecoEfficiencies(_fOut,"pPb","B","-",tp_recID_pp_B,tp_recID_pn_B,tp_recMS_pp_B,tp_recMS_pn_B);
            CalculatedWeightedRecoEfficiencies(_fOut,"pPb","B","+",tp_recID_pp_B,tp_recID_pn_B,tp_recMS_pp_B,tp_recMS_pn_B);
            CalculatedWeightedRecoEfficiencies(_fOut,"pPb","B","+-",tp_recID_pp_B,tp_recID_pn_B,tp_recMS_pp_B,tp_recMS_pn_B);

            delete tp_recID_pn_B;
            delete tp_recID_pp_B;
            delete tp_recMS_pn_B;
            delete tp_recMS_pp_B;
            delete tp_trig_pn_B;
            delete tp_trig_pp_B;
        }
    }
    double duration;
    duration = (std::clock()-start) / (double)CLOCKS_PER_SEC;
    std::cout << "Time: " << duration << " [sec]" << std::endl;
}
