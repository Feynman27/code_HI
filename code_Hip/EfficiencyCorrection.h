#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
//#include <memory>
#include "TF1.h"
#include "TH1.h"
#include "TH1F.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TCollection.h"
#include "TGraph.h"

#ifndef EfficiencyCorrection_h
#define EfficiencyCorrection_h
class EfficiencyCorrection{

    public:
        EfficiencyCorrection(std::vector<double>& vPt, std::vector<double>& vEta, std::vector<double>& vCent,std::string beam,std::string config, std::string charge) ; 
        EfficiencyCorrection(std::vector<double>& vPt, std::vector<double>& vEta, std::vector<double>& vCent,std::string beam,std::string config, std::string charge, std::string tpType) ; 
        virtual ~EfficiencyCorrection();
        std::vector<double> m_vPt;
        std::vector<double> m_vEta;
        std::vector<double> m_vCent;
        bool m_doWriteEtaCentDep;
        std::string m_beam;

        std::map<std::string,TH1F*>  m_hTotPos;
        std::map<std::string,TH1F*>  m_hTotNeg;
        std::map<std::string,TH1F*>  m_hTot;
        std::map<std::string,TH1F*>  m_hPassPos;
        std::map<std::string,TH1F*>  m_hPassNeg;
        std::map<std::string,TH1F*>  m_hPass;
        std::map<std::string,TEfficiency*>  m_hEffPos;
        std::map<std::string,TEfficiency*>  m_hEffNeg;
        std::map<std::string,TEfficiency*>  m_hEff;
        std::map<std::string,TEfficiency*>  m_hEffPosEta;
        std::map<std::string,TEfficiency*>  m_hEffNegEta;
        std::map<std::string,TEfficiency*>  m_hEffEta;
        bool m_isBothCharges ;
        bool m_isPos;
        bool m_isNeg;
        bool m_doAll;
        std::string m_eff;

        TH1* getCopyTotalHisto(int ipt, int icent, std::string charge,double weight=0.0);
        TH1* getCopyPassedHisto(int ipt, int icent, std::string charge,double weight=0.0);
        TEfficiency* getCopyEfficiencyHisto(int ipt, int icent, std::string charge);
        int indexPt(double pt);
        int indexEta(double eta);
        int indexCentrality(double centrality);
        virtual void fillHistos(double pt, double eta, double centrality, double charge, std::string type,float fcal);
        virtual void fillHistos(double pt, double eta, double centrality, double charge, std::string type);
        void writeHistos(TFile* fOut) ;
        void getHistoNamesPos(std::string& hnametot,  std::string& hnamepass, std::string& hnameeff,int ipt, int ieta, int icent);
        void getHistoNamesNeg(std::string& hnametot,  std::string& hnamepass, std::string& hnameeff,int ipt, int ieta, int icent);
        void getHistoNames(std::string& hnametot,  std::string& hnamepass, std::string& hnameeff,int ipt, int ieta, int icent,std::string charge);
        void getHistoEffNamePtCent(std::string& hnameeff, int ipt, int icent,std::string charge);
        void getHistoEffNamePtEta(std::string& hnameeff, int ipt, int ieta, std::string charge);
        void getHistoEffName(std::string& hnameeff, int ipt, int ieta, int icent,std::string charge);
        void CalculateEfficiency();
        void getEtaCentDep();
        std::vector<double> getPtBins();
        std::vector<double> getEtaBins();
        std::vector<double> getCentBins();
        double GetEfficiency(int ipt, int ieta, int icent, std::string charge);
        double GetEfficiencyErrorLow(int ipt, int ieta, int icent, std::string charge);
        double GetEfficiencyErrorUp(int ipt, int ieta, int icent, std::string charge);
    private:
        TFile* m_fFcalScF;
};
#endif

