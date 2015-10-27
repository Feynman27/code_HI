#include "EfficiencyCorrection.h"

EfficiencyCorrection::EfficiencyCorrection( std::vector<double>& vPt, std::vector<double>& vEta, std::vector<double>& vCent, std::string beam, std::string config, std::string charge,
                                                std::string tpType) : m_vPt(vPt), m_vEta(vEta), m_vCent(vCent),m_doWriteEtaCentDep(false),m_isBothCharges(false),m_isPos(false),m_isNeg(false) {

    std::stringstream beamType, effType; 
    beamType << beam << "_" << config ;
    m_beam=beamType.str();
    /*if(strcmp(beam.c_str(),"data")==0) m_doAll=true;
    else if(strcmp(beam.c_str(),"doAll")==0) m_doAll=true;
    else m_doAll=true;
    */
    m_doAll=true;
    effType << tpType;
    m_eff=effType.str();
    std::cout << "[EfficiencyCorrection::EfficiencyCorrection] Running on " << beam << " collisions." << std::endl;
    std::cout << "[EfficiencyCorrection::EfficiencyCorrection] Period " << config << " configuration." << std::endl;
    std::cout << "[EfficiencyCorrection::EfficiencyCorrection] Initializing histograms for " << m_vPt.size()-1 << " pt bins, " <<  
        m_vEta.size()-1 << " eta bins, " << m_vCent.size()-1 << " centrality bins." << std::endl;
    if(m_doAll) std::cout << "[EfficiencyCorrection::EfficiencyCorrection] Will do mu+, mu-, and mu^{#pm}." << std::endl;

    if(strcmp(charge.c_str(),"+")==0) {
        m_isPos= true;
        std::cout << "[EfficiencyCorrection::EfficiencyCorrection] W+ --> mu+ nu" << std::endl;
    }
    else if(strcmp(charge.c_str(),"-")==0) {
        m_isNeg = true;
        std::cout << "[EfficiencyCorrection::EfficiencyCorrection] W- --> mu- nu~" << std::endl;
    }
    else if(strcmp(charge.c_str(),"+-")==0) {
        m_isBothCharges = true;
        std::cout << "[EfficiencyCorrection::EfficiencyCorrection] W-+ --> mu-+ nu" << std::endl;
    }
    else {

        std::cout << "[EfficiencyCorrection::EfficiencyCorrection] An error has occured. Please specify charge." << std::endl;
        return;
    }
    std::cout << "[EfficiencyCorrection::EfficiencyCorrection] Initializing histograms for T&P analysis for " << tpType << std::endl;

    const int nPt = m_vPt.size();
    double ptBins[nPt];
    const int nEta = m_vEta.size();
    double etaBins[nEta];
    const int nCent = m_vCent.size();
    double centBins[nCent];
    // copy vector-->array for histo binning
    std::copy(m_vEta.begin(), m_vEta.end(), etaBins);
    std::copy(m_vCent.begin(), m_vCent.end(), centBins);
    std::copy(m_vPt.begin(), m_vPt.end(), ptBins);

    for(int ipt=0; ipt<nPt-1; ++ipt){
        for(int ieta=0; ieta<nEta-1; ++ieta){
            for(int icent=0; icent<nCent-1; ++icent){
                std::string hnametot, hnamepass, hnameeff;
                std::string hname_eff_eta, hname_eff_cent;
                if(m_isPos||m_doAll){
                    //mu+
                    getHistoNames(hnametot, hnamepass, hnameeff,ipt,ieta,icent,"+");
                    TH1F* hTotPos = new TH1F(hnametot.c_str(),"mu^{+}, total", nEta-1,etaBins);
                    m_hTotPos.insert(std::make_pair<std::string&, TH1F*& > (hnametot,hTotPos));
                    TH1F* hPassPos = new TH1F(hnamepass.c_str(),"mu^{+}, passed", nEta-1,etaBins);
                    m_hPassPos.insert(std::make_pair<std::string&, TH1F*& > (hnamepass,hPassPos));
                    TEfficiency* hEffPos = new TEfficiency(hnameeff.c_str(),"mu^{+},eff",nEta-1,etaBins);
                    m_hEffPos.insert(std::make_pair<std::string&, TEfficiency*& > (hnameeff,hEffPos));
                    // eta dep
                    getHistoEffNamePtCent(hname_eff_eta, ipt, icent,"+");
                    TEfficiency* hEffPosEta = new TEfficiency(hname_eff_eta.c_str(),"mu^{+},eff,all eta",nEta-1,etaBins);
                    m_hEffPosEta.insert(std::make_pair<std::string&, TEfficiency*& > (hname_eff_eta,hEffPosEta));
                }
                if(m_isNeg||m_doAll){
                    //mu-
                    getHistoNames(hnametot, hnamepass, hnameeff,ipt,ieta,icent,"-");
                    TH1F* hTotNeg = new TH1F(hnametot.c_str(),"mu^{-}, total", nEta-1,etaBins);
                    m_hTotNeg.insert(std::make_pair<std::string&, TH1F*& > (hnametot,hTotNeg));
                    TH1F* hPassNeg = new TH1F(hnamepass.c_str(),"mu^{-}, passed", nEta-1,etaBins);
                    m_hPassNeg.insert(std::make_pair<std::string&, TH1F*& > (hnamepass,hPassNeg));
                    TEfficiency* hEffNeg = new TEfficiency(hnameeff.c_str(),"mu^{-},eff",nEta-1,etaBins);
                    m_hEffNeg.insert(std::make_pair<std::string&, TEfficiency*& > (hnameeff,hEffNeg));
                    // eta dep
                    getHistoEffNamePtCent(hname_eff_eta, ipt, icent,"-");
                    TEfficiency* hEffNegEta = new TEfficiency(hname_eff_eta.c_str(),"mu^{-},eff,all eta",nEta-1,etaBins);
                    m_hEffNegEta.insert(std::make_pair<std::string&, TEfficiency*& > (hname_eff_eta,hEffNegEta));
                }
                if(m_isBothCharges||m_doAll){
                    //mu-
                    getHistoNames(hnametot, hnamepass, hnameeff,ipt,ieta,icent,"+-");
                    TH1F* hTot = new TH1F(hnametot.c_str(),"mu^{#pm}, total", nEta-1,etaBins);
                    m_hTot.insert(std::make_pair<std::string&, TH1F*& > (hnametot,hTot));
                    TH1F* hPass = new TH1F(hnamepass.c_str(),"mu^{#pm}, passed", nEta-1,etaBins);
                    m_hPass.insert(std::make_pair<std::string&, TH1F*& > (hnamepass,hPass));
                    TEfficiency* hEff = new TEfficiency(hnameeff.c_str(),"mu^{#pm},eff",nEta-1,etaBins);
                    m_hEff.insert(std::make_pair<std::string&, TEfficiency*& > (hnameeff,hEff));
                    // eta dep
                    getHistoEffNamePtCent(hname_eff_eta, ipt, icent,"+-");
                    TEfficiency* hEffEta = new TEfficiency(hname_eff_eta.c_str(),"mu^{#pm},eff,all eta",nEta-1,etaBins);
                    m_hEffEta.insert(std::make_pair<std::string&, TEfficiency*& > (hname_eff_eta,hEffEta));
                }
                if(!m_isPos&&!m_isNeg&&!m_isBothCharges&&!m_doAll){ 
                    std::cout << "[EfficiencyCorrection::EfficiencyCorrection] ERROR: No histograms initialized." << std::endl;
                    return;
                }

            }//icent
        }//ieta
    }//ipt
}

EfficiencyCorrection::EfficiencyCorrection( std::vector<double>& vPt, std::vector<double>& vEta, std::vector<double>& vCent, std::string beam, std::string config, std::string charge) : 
    m_vPt(vPt), m_vEta(vEta), m_vCent(vCent),m_doWriteEtaCentDep(false),m_isBothCharges(false),m_isPos(false),m_isNeg(false), m_doAll(false), m_eff("Cw") {

    m_fFcalScF = new TFile("fcalScaleFactors.09.14.2014.root","read") ;

    std::stringstream type; 
    type << beam << "_" << config ;
    m_beam=type.str();
    std::cout << "[EfficiencyCorrection::EfficiencyCorrection] Running on " << beam << " collisions." << std::endl;
    std::cout << "[EfficiencyCorrection::EfficiencyCorrection] Period " << config << " configuration." << std::endl;
    std::cout << "[EfficiencyCorrection::EfficiencyCorrection] Initializing histograms for " << m_vPt.size()-1 << " pt bins, " <<  
        m_vEta.size()-1 << " eta bins, " << m_vCent.size()-1 << " centrality bins." << std::endl;
    std::cout << "[EfficiencyCorrection::EfficiencyCorrection] Run to calculate " << m_eff << " efficiency." << std::endl;

    if(strcmp(charge.c_str(),"+")==0) {
        m_isPos= true;
        std::cout << "[EfficiencyCorrection::EfficiencyCorrection] W+ --> mu+ nu" << std::endl;
    }
    else if(strcmp(charge.c_str(),"-")==0) {
        m_isNeg = true;
        std::cout << "[EfficiencyCorrection::EfficiencyCorrection] W- --> mu- nu~" << std::endl;
    }
    else if(strcmp(charge.c_str(),"+-")==0) {
        m_isBothCharges = true;
        std::cout << "[EfficiencyCorrection::EfficiencyCorrection] W-+ --> mu-+ nu" << std::endl;
    }
    else {

        std::cout << "[EfficiencyCorrection::EfficiencyCorrection] An error has occured. Please specify charge." << std::endl;
        return;
    }

    const int nPt = m_vPt.size();
    double ptBins[nPt];
    const int nEta = m_vEta.size();
    double etaBins[nEta];
    const int nCent = m_vCent.size();
    double centBins[nCent];
    // copy vector-->array for histo binning
    std::copy(m_vEta.begin(), m_vEta.end(), etaBins);
    std::copy(m_vCent.begin(), m_vCent.end(), centBins);
    std::copy(m_vPt.begin(), m_vPt.end(), ptBins);

    for(int ipt=0; ipt<nPt-1; ++ipt){
        for(int ieta=0; ieta<nEta-1; ++ieta){
            for(int icent=0; icent<nCent-1; ++icent){
                std::string hnametot, hnamepass, hnameeff;
                std::string hname_eff_eta, hname_eff_cent;
                if(m_isPos){
                    //mu+
                    getHistoNames(hnametot, hnamepass, hnameeff,ipt,ieta,icent,"+");
                    TH1F* hTotPos = new TH1F(hnametot.c_str(),"mu^{+}, total", nEta-1,etaBins);
                    m_hTotPos.insert(std::make_pair<std::string&, TH1F*& > (hnametot,hTotPos));
                    TH1F* hPassPos = new TH1F(hnamepass.c_str(),"mu^{+}, passed", nEta-1,etaBins);
                    m_hPassPos.insert(std::make_pair<std::string&, TH1F*& > (hnamepass,hPassPos));
                    TEfficiency* hEffPos = new TEfficiency(hnameeff.c_str(),"mu^{+},eff",nEta-1,etaBins);
                    m_hEffPos.insert(std::make_pair<std::string&, TEfficiency*& > (hnameeff,hEffPos));
                    // eta dep
                    getHistoEffNamePtCent(hname_eff_eta, ipt, icent,"+");
                    TEfficiency* hEffPosEta = new TEfficiency(hname_eff_eta.c_str(),"mu^{+},eff,all eta",nEta-1,etaBins);
                    m_hEffPosEta.insert(std::make_pair<std::string&, TEfficiency*& > (hname_eff_eta,hEffPosEta));
                }
                else if(m_isNeg){
                    //mu-
                    getHistoNames(hnametot, hnamepass, hnameeff,ipt,ieta,icent,"-");
                    TH1F* hTotNeg = new TH1F(hnametot.c_str(),"mu^{-}, total", nEta-1,etaBins);
                    m_hTotNeg.insert(std::make_pair<std::string&, TH1F*& > (hnametot,hTotNeg));
                    TH1F* hPassNeg = new TH1F(hnamepass.c_str(),"mu^{-}, passed", nEta-1,etaBins);
                    m_hPassNeg.insert(std::make_pair<std::string&, TH1F*& > (hnamepass,hPassNeg));
                    TEfficiency* hEffNeg = new TEfficiency(hnameeff.c_str(),"mu^{-},eff",nEta-1,etaBins);
                    m_hEffNeg.insert(std::make_pair<std::string&, TEfficiency*& > (hnameeff,hEffNeg));
                    //eta dep
                    getHistoEffNamePtCent(hname_eff_eta, ipt, icent,"-");
                    TEfficiency* hEffNegEta = new TEfficiency(hname_eff_eta.c_str(),"mu^{-},eff,all eta",nEta-1,etaBins);
                    m_hEffNegEta.insert(std::make_pair<std::string&, TEfficiency*& > (hname_eff_eta,hEffNegEta));
                }
                else if(m_isBothCharges){
                    getHistoNames(hnametot, hnamepass, hnameeff,ipt,ieta,icent,"+-");
                    TH1F* hTot = new TH1F(hnametot.c_str(),"mu^{#pm}, total", nEta-1,etaBins);
                    m_hTot.insert(std::make_pair<std::string&, TH1F*& > (hnametot,hTot));
                    TH1F* hPass = new TH1F(hnamepass.c_str(),"mu^{#pm}, passed", nEta-1,etaBins);
                    m_hPass.insert(std::make_pair<std::string&, TH1F*& > (hnamepass,hPass));
                    TEfficiency* hEff = new TEfficiency(hnameeff.c_str(),"mu^{#pm},eff",nEta-1,etaBins);
                    m_hEff.insert(std::make_pair<std::string&, TEfficiency*& > (hnameeff,hEff));
                    // eta dep
                    getHistoEffNamePtCent(hname_eff_eta, ipt, icent,"+-");
                    TEfficiency* hEffEta = new TEfficiency(hname_eff_eta.c_str(),"mu^{#pm},eff,all eta",nEta-1,etaBins);
                    m_hEffEta.insert(std::make_pair<std::string&, TEfficiency*& > (hname_eff_eta,hEffEta));
                }
                else { 
                    std::cout << "[EfficiencyCorrection::EfficiencyCorrection] ERROR: No histograms initialized." << std::endl;
                    return;
                }
            }//icent
        }//ieta

    }//ipt
    std::cout << "[EfficiencyCorrection::EfficiencyCorrection] Histograms initialized and ready for use." << std::endl;

}

EfficiencyCorrection::~EfficiencyCorrection(){

    std::cout << "[EfficiencyCorrection::~EfficiencyCorrection] D'tor" << std::endl;

    delete m_fFcalScF;
    std::map<std::string,TH1F*>::iterator it;
    for(it=m_hPass.begin(); it!=m_hPass.end(); ++it){

        if(it->second->GetEntries()>0){
            
            delete it->second;
            
        }
    }
    for(it=m_hPassPos.begin(); it!=m_hPassPos.end(); ++it){

        if(it->second->GetEntries()>0){
            
            delete it->second;
            
        }
    }
    for(it=m_hPassNeg.begin(); it!=m_hPassNeg.end(); ++it){

        if(it->second->GetEntries()>0){
            
            delete it->second;
            
        }
    }

    std::map<std::string,TEfficiency*>::iterator it2;
    for(it=m_hTot.begin(),it2=m_hEff.begin(); it!=m_hTot.end(); ++it,++it2){

        if(it->second->GetEntries()>0){
            
            delete it->second;
            
            
            delete it2->second;
            
        }
    }

    for(it=m_hTotPos.begin(),it2=m_hEffPos.begin(); it!=m_hTotPos.end(); ++it,++it2){

        if(it->second->GetEntries()>0){
            
            delete it->second;
            
            
            delete it2->second;
            
        }
    }

    for(it=m_hTotNeg.begin(),it2=m_hEffNeg.begin(); it!=m_hTotNeg.end(); ++it,++it2){

        if(it->second->GetEntries()>0){
            
            delete it->second;
            
            
            delete it2->second;
            
        }
    }
    if(m_doWriteEtaCentDep==true) {
        
        for(it2=m_hEffPosEta.begin(); it2!=m_hEffPosEta.end(); ++it2){

            if(it2->second->GetPassedHistogram()->GetEntries()>0){
                
                delete it2->second;
                
            }
        }
        for(it2=m_hEffNegEta.begin(); it2!=m_hEffNegEta.end(); ++it2){

            if(it2->second->GetPassedHistogram()->GetEntries()>0){
                
                delete it2->second;
                
            }
        }
        for(it2=m_hEffEta.begin(); it2!=m_hEffEta.end(); ++it2){

            if(it2->second->GetPassedHistogram()->GetEntries()>0){
                
                delete it2->second;
                
            }
        }
    }
}

std::vector<double> EfficiencyCorrection::getPtBins(){
     std::vector<double> c_vPt (m_vPt) ;
     return c_vPt;
}

std::vector<double> EfficiencyCorrection::getEtaBins(){
    std::vector<double> c_vEta (m_vEta);
    return c_vEta;
}

std::vector<double> EfficiencyCorrection::getCentBins(){
    std::vector<double> c_vCent (m_vCent);
    return c_vCent;
}

TH1* EfficiencyCorrection::getCopyPassedHisto(int ipt, int icent, std::string charge,double weight){

            TH1* copy_passed=0;
            bool doScale = true;
            if(fabs(weight)<1e-9) doScale=false; 
            std::string name;
            getHistoEffNamePtCent(name,ipt,icent,charge);
            if(strcmp(charge.c_str(),"+")==0) {
                copy_passed = m_hEffPosEta[name]->GetCopyPassedHisto();
                if(doScale) copy_passed->Scale(weight/copy_passed->Integral());
            }
            else if(strcmp(charge.c_str(),"-")==0) {
                copy_passed = m_hEffNegEta[name]->GetCopyPassedHisto();
                if(doScale) copy_passed->Scale(weight/copy_passed->Integral());
            }
            else if(strcmp(charge.c_str(),"+-")==0) {
                copy_passed = m_hEffEta[name]->GetCopyPassedHisto();
                if(doScale) copy_passed->Scale(weight/copy_passed->Integral());
            }
            else {
                std::cout << "[EfficiencyCorrection::GetCopyPassedHistogram] ERROR: Cannot find histogram." << std::endl;
                return 0;
            }

            return copy_passed;
}

TH1* EfficiencyCorrection::getCopyTotalHisto(int ipt, int icent, std::string charge,double weight){

            TH1* copy_total=0;
            bool doScale = true;
            if(fabs(weight)<1e-9) doScale=false; 
            std::string name;
            getHistoEffNamePtCent(name,ipt,icent,charge);
            if(strcmp(charge.c_str(),"+")==0) {
                copy_total = m_hEffPosEta[name]->GetCopyTotalHisto();
                if(doScale) copy_total->Scale(weight/copy_total->Integral());
            }
            else if(strcmp(charge.c_str(),"-")==0) {
                copy_total = m_hEffNegEta[name]->GetCopyTotalHisto();
                if(doScale) copy_total->Scale(weight/copy_total->Integral());
            }
            else if(strcmp(charge.c_str(),"+-")==0) {
                copy_total = m_hEffEta[name]->GetCopyTotalHisto();
                if(doScale) copy_total->Scale(weight/copy_total->Integral());
            }
            else {
                std::cout << "[EfficiencyCorrection::GetCopyTotalHistogram] ERROR: Cannot find histogram." << std::endl;
                return 0;
            }

            return copy_total;
}

TEfficiency* EfficiencyCorrection::getCopyEfficiencyHisto(int ipt, int icent, std::string charge){

            TEfficiency* copy_eff=0;
            std::string name;
            getHistoEffNamePtCent(name,ipt,icent,charge);
            if(strcmp(charge.c_str(),"+")==0) {
                copy_eff = m_hEffPosEta[name];
            }
            else if(strcmp(charge.c_str(),"-")==0) {
                copy_eff = m_hEffNegEta[name];
            }
            else if(strcmp(charge.c_str(),"+-")==0) {
                copy_eff = m_hEffEta[name];
            }
            else {
                std::cout << "[EfficiencyCorrection::GetCopyTotalHistogram] ERROR: Cannot find histogram." << std::endl;
                return 0;
            }

            return copy_eff;
}

double EfficiencyCorrection::GetEfficiencyErrorUp(int ipt, int ieta, int icent, std::string charge){

            TEfficiency* peff;
            std::string name;
            getHistoEffName(name,ipt,ieta,icent,charge);
            if(strcmp(charge.c_str(),"+")==0) {
                peff = m_hEffPos[name];
            }
            else if(strcmp(charge.c_str(),"-")==0) {
                peff = m_hEffNeg[name];
            }
            else if(strcmp(charge.c_str(),"+-")==0) {
                peff = m_hEff[name];
            }
            else {
                std::cout << "[EfficiencyCorrection::GetEfficiency] Error. Cannot find charge species." << std::endl;
                return -9999.0;
            }

            for(int ibin=1; ibin<=peff->GetTotalHistogram()->GetNbinsX(); ++ibin){
                double bincontent = peff->GetTotalHistogram()->GetBinContent(ibin);
                if(bincontent>0.) return peff->GetEfficiencyErrorUp(ibin);
                
            }
    
            std::cout << "[EfficiencyCorrection::GetEfficiency] Error. Bins are empty." << std::endl;
            return -9999.0;
}

double EfficiencyCorrection::GetEfficiencyErrorLow(int ipt, int ieta, int icent, std::string charge){

            TEfficiency* peff;
            std::string name;
            getHistoEffName(name,ipt,ieta,icent,charge);
            if(strcmp(charge.c_str(),"+")==0) {
                peff = m_hEffPos[name];
            }
            else if(strcmp(charge.c_str(),"-")==0) {
                peff = m_hEffNeg[name];
            }
            else if(strcmp(charge.c_str(),"+-")==0) {
                peff = m_hEff[name];
            }
            else {
                std::cout << "[EfficiencyCorrection::GetEfficiency] Error. Cannot find charge species." << std::endl;
                return -9999.0;
            }

            for(int ibin=1; ibin<=peff->GetTotalHistogram()->GetNbinsX(); ++ibin){
                double bincontent = peff->GetTotalHistogram()->GetBinContent(ibin);
                if(bincontent>0.) return peff->GetEfficiencyErrorLow(ibin);
            }
    
            std::cout << "[EfficiencyCorrection::GetEfficiency] Error. Bins are empty." << std::endl;
            return -9999.0;
}

double EfficiencyCorrection::GetEfficiency(int ipt, int ieta, int icent, std::string charge){

            TEfficiency* peff;
            std::string name;
            getHistoEffName(name,ipt,ieta,icent,charge);
            if(strcmp(charge.c_str(),"+")==0) {
                peff = m_hEffPos[name];
            }
            else if(strcmp(charge.c_str(),"-")==0) {
                peff = m_hEffNeg[name];
            }
            else if(strcmp(charge.c_str(),"+-")==0) {
                peff = m_hEff[name];
            }
            else {
                std::cout << "[EfficiencyCorrection::GetEfficiency] Error. Cannot find charge species." << std::endl;
                return -9999.0;
            }

            for(int ibin=1; ibin<=peff->GetTotalHistogram()->GetNbinsX(); ++ibin){
                double bincontent = peff->GetTotalHistogram()->GetBinContent(ibin);
                if(bincontent>0.) return peff->GetEfficiency(ibin);
            }
            std::cout << "[EfficiencyCorrection::GetEfficiency] Error. Bins are empty." << std::endl;
            return -9999.0;
    
}


void EfficiencyCorrection::getHistoEffNamePtCent(std::string& hnameeff, int ipt, int icent,std::string charge){

    std::stringstream shnameeff;
    shnameeff.str(std::string());
    if(strcmp(charge.c_str(),"+")==0){
        shnameeff << m_beam<< "_" << m_eff << "_hPosEff" << "_pt" << ipt <<
            "_cent" << icent;
        hnameeff=shnameeff.str();
    }
    else if(strcmp(charge.c_str(),"-")==0){
        shnameeff <<m_beam<< "_" << m_eff << "_hNegEff" << "_pt" << ipt <<
            "_cent" << icent;
        hnameeff=shnameeff.str();
    }
    else if(strcmp(charge.c_str(),"+-")==0){
        shnameeff <<m_beam<< "_" << m_eff << "_hEff" << "_pt" << ipt <<
            "_cent" << icent;
        hnameeff=shnameeff.str();
    }
    
    else std::cout << "[]ERROR: Cannot find key." << std::endl;
}

void EfficiencyCorrection::getHistoEffNamePtEta(std::string& hnameeff, int ipt, int ieta, std::string charge){

    std::stringstream shnameeff;
    shnameeff.str(std::string());
    if(strcmp(charge.c_str(),"+")==0){
        shnameeff << m_beam << "_" << m_eff << "_hPosEff" << "_pt" << ipt <<
            "_eta" << ieta ;
        hnameeff=shnameeff.str();
    }
    else if(strcmp(charge.c_str(),"-")==0){
        shnameeff <<m_beam << "_" << m_eff << "_hNegEff" << "_pt" << ipt <<
            "_eta" << ieta; 
        hnameeff=shnameeff.str();
    }
    else if(strcmp(charge.c_str(),"+-")==0){
        shnameeff <<m_beam << "_" << m_eff << "_hEff" << "_pt" << ipt <<
            "_eta" << ieta; 
        hnameeff=shnameeff.str();
    }
    else std::cout << "[]ERROR: Cannot find key." << std::endl;
}

void EfficiencyCorrection::getHistoEffName(std::string& hnameeff, int ipt, int ieta, int icent,std::string charge){

    std::stringstream shnameeff;
    shnameeff.str(std::string());
    if(strcmp(charge.c_str(),"+")==0){
        shnameeff <<m_beam << "_" << m_eff << "_hPosEff" << "_pt" << ipt <<
            "_eta" << ieta <<
            "_cent" << icent;
        hnameeff=shnameeff.str();
    }
    else if(strcmp(charge.c_str(),"-")==0){
        shnameeff <<m_beam << "_" << m_eff << "_hNegEff" << "_pt" << ipt <<
            "_eta" << ieta <<
            "_cent" << icent;
        hnameeff=shnameeff.str();
    }
    else if(strcmp(charge.c_str(),"+-")==0){
        shnameeff <<m_beam << "_" << m_eff << "_hEff" << "_pt" << ipt <<
            "_eta" << ieta <<
            "_cent" << icent;
        hnameeff=shnameeff.str();
    }
    else std::cout << "[]ERROR: Cannot find key." << std::endl;
}

void EfficiencyCorrection::getHistoNamesPos(std::string& hnametot,  std::string& hnamepass, std::string& hnameeff,
                                                    int ipt, int ieta, int icent){

    std::stringstream shnametot;
    std::stringstream shnamepass;
    std::stringstream shnameeff;
    shnametot.str(std::string());
    shnametot << m_beam << "_hPosTot" << "_pt" << ipt <<
        "_eta" << ieta <<
        "_cent" << icent;
    hnametot=shnametot.str();
    shnamepass.str(std::string());
    shnamepass << m_beam << "_hPosPass" << "_pt" << ipt <<
        "_eta" << ieta <<
        "_cent" << icent;
    hnamepass=shnamepass.str();
    shnameeff.str(std::string());
    shnameeff << m_beam << "_hPosEff" << "_pt" << ipt <<
        "_eta" << ieta <<
        "_cent" << icent;
    hnameeff=shnameeff.str();

}

void EfficiencyCorrection::getHistoNamesNeg(std::string& hnametot,  std::string& hnamepass, std::string& hnameeff,
                                                    int ipt, int ieta, int icent){

    std::stringstream shnametot;
    std::stringstream shnamepass;
    std::stringstream shnameeff;
    shnametot.str(std::string());
    shnametot << m_beam << "_hNegTot" << "_pt" << ipt <<
        "_eta" << ieta <<
        "_cent" << icent;
    hnametot=shnametot.str();
    shnamepass.str(std::string());
    shnamepass << m_beam << "_hNegPass" << "_pt" << ipt <<
        "_eta" << ieta <<
        "_cent" << icent;
    hnamepass=shnamepass.str();
    shnameeff.str(std::string());
    shnameeff << m_beam << "_hNegEff" << "_pt" << ipt <<
        "_eta" << ieta <<
        "_cent" << icent;
    hnameeff=shnameeff.str();

}

void EfficiencyCorrection::getHistoNames(std::string& hnametot,  std::string& hnamepass, std::string& hnameeff,
                                                    int ipt, int ieta, int icent, std::string charge){

    std::stringstream shnametot;
    std::stringstream shnamepass;
    std::stringstream shnameeff;

    shnametot.str(std::string());
    shnamepass.str(std::string());
    shnameeff.str(std::string());
    if(strcmp(charge.c_str(),"+")==0){
        shnametot << m_beam << "_" << m_eff <<  "_hPosTot" << 
            "_pt" << ipt <<
            "_eta" << ieta <<
            "_cent" << icent;
        shnamepass << m_beam <<  "_" << m_eff << "_hPosPass" << 
            "_pt" << ipt <<
            "_eta" << ieta <<
            "_cent" << icent;

        shnameeff << m_beam << "_" << m_eff << "_hPosEff" << 
            "_pt" << ipt <<
            "_eta" << ieta <<
            "_cent" << icent;
    }
    else if(strcmp(charge.c_str(),"-")==0){
        shnametot << m_beam << "_" << m_eff <<  "_hNegTot" << 
            "_pt" << ipt <<
            "_eta" << ieta <<
            "_cent" << icent;
        shnamepass << m_beam <<  "_" << m_eff << "_hNegPass" << 
            "_pt" << ipt <<
            "_eta" << ieta <<
            "_cent" << icent;

        shnameeff << m_beam << "_" << m_eff << "_hNegEff" << 
            "_pt" << ipt <<
            "_eta" << ieta <<
            "_cent" << icent;

    }

    else if(strcmp(charge.c_str(),"+-")==0){

        shnametot << m_beam << "_" << m_eff <<  "_hTot" << 
            "_pt" << ipt <<
            "_eta" << ieta <<
            "_cent" << icent;

        shnamepass << m_beam <<  "_" << m_eff << "_hPass" << 
            "_pt" << ipt <<
            "_eta" << ieta <<
            "_cent" << icent;

        shnameeff << m_beam << "_" << m_eff << "_hEff" << 
            "_pt" << ipt <<
            "_eta" << ieta <<
            "_cent" << icent;

    }
    else {
        std::cout << "[EfficiencyCorrection::getHistoNames] ERROR: Please specify charge with \"+\", \"-\", or \"+-\" " << std::endl;
        return;
    }

    hnametot=shnametot.str();
    hnamepass=shnamepass.str();
    hnameeff=shnameeff.str();
}

int EfficiencyCorrection::indexPt(double pt){
   int nPtBins = m_vPt.size()-1;
   int index = 0;
   for(int ipt=0; ipt<nPtBins; ++ipt){
        if(pt>=m_vPt.at(ipt)&&pt<m_vPt.at(ipt+1)) return index;
        else ++index;
   }//ipt
   std::cout <<"[EfficiencyCorrection::indexPt] WARNING: Muon pt: " << pt << " is either outside bin limits or an error has occured." << std::endl;
   std::cout <<"[EfficiencyCorrection::indexPt] CHECK: " << m_vPt.at(0) << "<=" << pt << "<" << m_vPt.at(nPtBins) << std::endl;
   return -1;
}

int EfficiencyCorrection::indexEta(double eta){
   int nEtaBins = m_vEta.size()-1;
   int index = 0;
   for(int ieta=0; ieta<nEtaBins; ++ieta){
        if(eta>=m_vEta.at(ieta)&&eta<m_vEta.at(ieta+1)) return index;
        else ++index;
   }//ieta
   std::cout <<"[EfficiencyCorrection::indexEta] WARNING: Muon eta: " << eta << " is either outside bin limits or an error has occured." << std::endl;
   std::cout <<"[EfficiencyCorrection::indexEta] CHECK: " << m_vEta.at(0) << "<=" << eta << "<" << m_vEta.at(nEtaBins) << std::endl;
   return -1;
}

int EfficiencyCorrection::indexCentrality(double centrality){
   int nCentralityBins = m_vCent.size()-1;
   int index = 0;
   for(int icent=0; icent<nCentralityBins; ++icent){
        if(centrality>=m_vCent.at(icent)&&centrality<m_vCent.at(icent+1)) return index;
        else ++index;
   }//ieta
   std::cout <<"[EfficiencyCorrection::indexCentrality] WARNING: Muon centrality: " << centrality << " is either outside bin limits or an error has occured." << std::endl;
   std::cout <<"[EfficiencyCorrection::indexCentrality] CHECK: " << m_vCent.at(0) << "<=" << centrality << "<" << m_vCent.at(nCentralityBins) << std::endl;
   return -1;
}

void EfficiencyCorrection::fillHistos(double pt, double eta, double centrality, double charge, std::string type){

    bool isTotLevel = false; 
    bool isPassLevel = false; 
    if(strcmp(type.c_str(),"tot")==0) isTotLevel = true;
    else if(strcmp(type.c_str(),"pass")==0) isPassLevel = true;
    else{
        std::cout << "[EfficiencyCorrection::fillHistos] An error has occured. Please designate \"tot\" or \"pass\". " << std::endl;
        return;
    }

    int ipt   = indexPt(pt);
    int ieta  = indexEta(eta);
    int icent = indexCentrality(centrality);

    std::string hnametot, hnamepass, hnameeff;
    // mu^{+}
    if(charge>0.){  
        getHistoNames(hnametot, hnamepass, hnameeff,ipt,ieta,icent,"+");
        if(isTotLevel){ 
            m_hTotPos[hnametot]->Fill(eta);
        }
        else if(isPassLevel) {
            m_hPassPos[hnamepass]->Fill(eta);
        }
    }

    // mu^{-}
    else if(charge<0.){  
        getHistoNames(hnametot, hnamepass, hnameeff,ipt,ieta,icent,"-");
        if(isTotLevel){ 
            m_hTotNeg[hnametot]->Fill(eta);
        }
        else if(isPassLevel) {
            m_hPassNeg[hnamepass]->Fill(eta);
        }
    }

    if(m_isBothCharges||m_doAll){
        getHistoNames(hnametot, hnamepass, hnameeff,ipt,ieta,icent,"+-");
        if(isTotLevel){ 
            m_hTot[hnametot]->Fill(eta);
        }
        else if(isPassLevel) {
            m_hPass[hnamepass]->Fill(eta);
        }

    }
}
void EfficiencyCorrection::fillHistos(double pt, double eta, double centrality, double charge, std::string type, float fcal){

    bool isTotLevel = false; 
    bool isPassLevel = false; 
    if(strcmp(type.c_str(),"tot")==0) isTotLevel = true;
    else if(strcmp(type.c_str(),"pass")==0) isPassLevel = true;
    else{
        std::cout << "[EfficiencyCorrection::fillHistos] An error has occured. Please designate \"tot\" or \"pass\". " << std::endl;
        return;
    }

    int ipt   = indexPt(pt);
    int ieta  = indexEta(eta);
    int icent = indexCentrality(centrality);

    TF1* _funcScF = (TF1*)((TH1F*)m_fFcalScF->Get("_hRatio"))->GetFunction("_fit");
    double sf = _funcScF->Eval(fcal);

    std::string hnametot, hnamepass, hnameeff;
    // mu^{+}
    if(charge>0.){  
        getHistoNames(hnametot, hnamepass, hnameeff,ipt,ieta,icent,"+");
        if(isTotLevel){ 
            m_hTotPos[hnametot]->Fill(eta,sf);
        }
        else if(isPassLevel) {
            m_hPassPos[hnamepass]->Fill(eta,sf);
        }
    }

    // mu^{-}
    else if(charge<0.){  
        getHistoNames(hnametot, hnamepass, hnameeff,ipt,ieta,icent,"-");
        if(isTotLevel){ 
            m_hTotNeg[hnametot]->Fill(eta,sf);
        }
        else if(isPassLevel) {
            m_hPassNeg[hnamepass]->Fill(eta,sf);
        }
    }

    if(m_isBothCharges||m_doAll){
        getHistoNames(hnametot, hnamepass, hnameeff,ipt,ieta,icent,"+-");
        if(isTotLevel){ 
            m_hTot[hnametot]->Fill(eta,sf);
        }
        else if(isPassLevel) {
            m_hPass[hnamepass]->Fill(eta,sf);
        }

    }
}

void EfficiencyCorrection::writeHistos(TFile* fOut){

    fOut->cd();

    std::map<std::string,TH1F*>::iterator it;
    for(it=m_hPass.begin(); it!=m_hPass.end(); ++it){

        if(it->second->GetEntries()>0){
            std::cout << "Writing " << it->first << std::endl;
            it->second->Write();
            
        }
    }
    for(it=m_hPassPos.begin(); it!=m_hPassPos.end(); ++it){

        if(it->second->GetEntries()>0){
            std::cout << "Writing " << it->first << std::endl;
            it->second->Write();
            
        }
    }
    for(it=m_hPassNeg.begin(); it!=m_hPassNeg.end(); ++it){

        if(it->second->GetEntries()>0){
            std::cout << "Writing " << it->first << std::endl;
            it->second->Write();
            
        }
    }

    std::map<std::string,TEfficiency*>::iterator it2;
    for(it=m_hTot.begin(),it2=m_hEff.begin(); it!=m_hTot.end(); ++it,++it2){

        if(it->second->GetEntries()>0){
            std::cout << "Writing " << it->first << std::endl;
            it->second->Write();
            
            std::cout << "Writing " << it2->first << std::endl;
            it2->second->Write();
            
        }
    }

    for(it=m_hTotPos.begin(),it2=m_hEffPos.begin(); it!=m_hTotPos.end(); ++it,++it2){

        if(it->second->GetEntries()>0){
            std::cout << "Writing " << it->first << std::endl;
            it->second->Write();
            
            std::cout << "Writing " << it2->first << std::endl;
            it2->second->Write();
            
        }
    }

    for(it=m_hTotNeg.begin(),it2=m_hEffNeg.begin(); it!=m_hTotNeg.end(); ++it,++it2){

        if(it->second->GetEntries()>0){
            std::cout << "Writing " << it->first << std::endl;
            it->second->Write();
            
            std::cout << "Writing " << it2->first << std::endl;
            it2->second->Write();
            
        }
    }
    if(m_doWriteEtaCentDep==true) {
        
        for(it2=m_hEffPosEta.begin(); it2!=m_hEffPosEta.end(); ++it2){

            if(it2->second->GetPassedHistogram()->GetEntries()>0){
                std::cout << "Writing " << it2->first << std::endl;
                it2->second->Write();
                
            }
        }
        for(it2=m_hEffNegEta.begin(); it2!=m_hEffNegEta.end(); ++it2){

            if(it2->second->GetPassedHistogram()->GetEntries()>0){
                std::cout << "Writing " << it2->first << std::endl;
                it2->second->Write();
                
            }
        }
        for(it2=m_hEffEta.begin(); it2!=m_hEffEta.end(); ++it2){

            if(it2->second->GetPassedHistogram()->GetEntries()>0){
                std::cout << "Writing " << it2->first << std::endl;
                it2->second->Write();
                
            }
        }
    }
}

void EfficiencyCorrection::CalculateEfficiency(){

    const int nPtBins   = m_vPt.size()-1;
    const int nEtaBins  = m_vEta.size()-1;
    const int nCentBins = m_vCent.size()-1;

    /*TFile* _fScF = new TFile("fcalScaleFactors.09.14.2014.root","read") ;
    TGraph* _grScF = (TGraph*)_fScF->Get("grScaleFactors");
    */
    for(int ipt=0; ipt<nPtBins; ++ipt){
        for(int ieta=0; ieta<nEtaBins; ++ieta){
            for(int icent=0; icent<nCentBins; ++icent){
                // fcal scale factor
                /*if(_grScF->GetN() != nCentBins){
                    std::cout << "[EfficiencyCorrection::CalculateEfficiency] ERROR: Bins in scale factors do not match analysis bins." << std::endl;
                    return;
                }

                double scale_factor = _grScF->GetY()[icent];
                */
                std::string hnametot;
                std::string hnamepass;
                std::string hnameeff;
                //mu+
                getHistoNames(hnametot, hnamepass, hnameeff,ipt,ieta,icent,"+");
                if( (m_isPos||m_doAll) && (m_hTotPos[hnametot]->GetEntries()>0) ){
                    bool cleared1 = m_hEffPos[hnameeff]->SetTotalHistogram(*m_hTotPos[hnametot],"f");
                    bool cleared2 = m_hEffPos[hnameeff]->SetPassedHistogram(*m_hPassPos[hnamepass],"f");
                    if(cleared1==false||cleared2==false){
                        std::cout << "[EfficiencyCorrection::CalculateEfficiency] WARNING: Histogram not replaced." << 
                                "This can lead to inconsistent histograms and useless results or unexpected behaviour." << std::endl;
                    }
                    //m_hEffPos[hnameeff]->SetWeight(scale_factor);
                    
                }
                //mu-
                getHistoNames(hnametot, hnamepass, hnameeff,ipt,ieta,icent,"-");
                if((m_isNeg||m_doAll) && (m_hTotNeg[hnametot]->GetEntries()>0) ){
                    bool cleared1=m_hEffNeg[hnameeff]->SetTotalHistogram(*m_hTotNeg[hnametot],"f");
                    bool cleared2=m_hEffNeg[hnameeff]->SetPassedHistogram(*m_hPassNeg[hnamepass],"f");
                    if(cleared1==false||cleared2==false){
                        std::cout << "[EfficiencyCorrection::CalculateEfficiency] WARNING: Histogram not replaced." << 
                                "This can lead to inconsistent histograms and useless results or unexpected behaviour." << std::endl;
                    }
                    //m_hEffNeg[hnameeff]->SetWeight(scale_factor);
                }
                //mu^pm
                getHistoNames(hnametot, hnamepass, hnameeff,ipt,ieta,icent,"+-");
                if((m_isBothCharges||m_doAll) && (m_hTot[hnametot]->GetEntries()>0) ){
                    bool cleared1=m_hEff[hnameeff]->SetTotalHistogram(*m_hTot[hnametot],"f");
                    bool cleared2=m_hEff[hnameeff]->SetPassedHistogram(*m_hPass[hnamepass],"f");
                    if(cleared1==false||cleared2==false){
                        std::cout << "[EfficiencyCorrection::CalculateEfficiency] WARNING: Histogram not replaced." << 
                                "This can lead to inconsistent histograms and useless results or unexpected behaviour." << std::endl;
                    }
                    //m_hEff[hnameeff]->SetWeight(scale_factor);
                }

            } //icent
        } //ieta
    } //ipt

    //delete _fScF;
}

void EfficiencyCorrection::getEtaCentDep(){

    m_doWriteEtaCentDep = true;
    const int nPtBins   = m_vPt.size()-1;
    const int nEtaBins  = m_vEta.size()-1;
    const int nCentBins = m_vCent.size()-1;

    for(int ipt=0; ipt<nPtBins; ++ipt){
      for(int ieta=0; ieta<nEtaBins; ++ieta){
        for(int icent=0; icent<nCentBins; ++icent){
            std::string currentName, name;
            if(m_isPos||m_doAll){
                //mu+
                getHistoEffName(currentName,ipt,ieta,icent,"+");
                getHistoEffNamePtCent(name,ipt,icent,"+");
                m_hEffPosEta[name]->Add(*m_hEffPos[currentName]); 
                //getHistoEffNamePtEta(name,ipt,ieta,"+");
                //m_hEffPosCent[name]->Add(*m_hEffPos[currentName]); 
            }
            //mu-
            if(m_isNeg||m_doAll){
                getHistoEffName(currentName,ipt,ieta,icent,"-");
                getHistoEffNamePtCent(name,ipt,icent,"-");
                m_hEffNegEta[name]->Add(*m_hEffNeg[currentName]); 
                //getHistoEffNamePtEta(name,ipt,ieta,"-");
                //m_hEffNegCent[name]->Add(*m_hEffNeg[currentName]); 
            }
            if(m_isBothCharges||m_doAll){
                //mu
                getHistoEffName(currentName,ipt,ieta,icent,"+-");
                getHistoEffNamePtCent(name,ipt,icent,"+-");
                m_hEffEta[name]->Add(*m_hEff[currentName]); 
                //getHistoEffNamePtEta(name,ipt,ieta,"+-");
                //m_hEffCent[name]->Add(*m_hEff[currentName]); 
            }
        }
      }
    }
}

