#include <map>
#include <sstream>
#include <iostream>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TH1.h"
#include "TEfficiency.h"
#include "TStyle.h"
#include "TColor.h"
#include "TPad.h"

void plotVerbose(TEfficiency* teff, TCanvas* c, TString title){

    c->cd();
    std::cout << c->GetName()<<std::endl;
// ------------>Primitives in pad: c_1

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TPad *c_1 = new TPad("c_1", "c_1",0,0.008474576,1,0.9915254);
    c_1->Draw();
    c_1->cd();
    c_1->Range(-3.5,-2446.54,2.827273,3258.188);
    c_1->SetFillColor(0);
    c_1->SetBorderMode(0);
    c_1->SetBorderSize(2);
    c_1->SetTickx(1);
    c_1->SetTicky(1);
    c_1->SetLeftMargin(0.158046);
    c_1->SetRightMargin(0.05172414);
    c_1->SetTopMargin(0.04525862);
    c_1->SetBottomMargin(0.4288793);
    c_1->SetFrameBorderMode(0);
    c_1->SetFrameBorderMode(0);


    ((TH1F*)teff->GetTotalHistogram())->SetFillColor(kRed);
    ((TH1F*)teff->GetTotalHistogram())->Draw("hist");
    ((TH1F*)teff->GetPassedHistogram())->SetFillColor(kYellow);
    ((TH1F*)teff->GetPassedHistogram())->Draw("histsame");
    ((TH1F*)teff->GetTotalHistogram())->Draw("sameaxis");
    double ymax = ((TH1F*)teff->GetTotalHistogram())->GetMaximum();
    ((TH1F*)teff->GetTotalHistogram())->GetYaxis()->SetRangeUser(0.1,ymax+5.0);
    ((TH1F*)teff->GetTotalHistogram())->GetXaxis()->SetTitle("#eta_{lab}^{#mu}");
    ((TH1F*)teff->GetTotalHistogram())->GetYaxis()->SetTitle("Probes");
    
    TLatex *   tex = new TLatex(0.6221264,0.4892241,title);
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04310345);
    tex->SetLineWidth(2);
    tex->Draw();
   
    TPad* c_2 = new TPad("c_2", "c_2",0,0.00862069,1,0.4331897);
    c_2->Draw();
    c_2->cd();
    c_2->Range(-3.5,-0.08581081,2.827273,1.112162);
    c_2->SetFillColor(0);
    c_2->SetBorderMode(0);
    c_2->SetBorderSize(2);
    c_2->SetTickx(1);
    c_2->SetTicky(1);
    c_2->SetLeftMargin(0.158046);
    c_2->SetRightMargin(0.05172414);
    c_2->SetTopMargin(0.01015228);
    c_2->SetBottomMargin(0.2385787);
    c_2->SetFrameBorderMode(0);
    c_2->SetFrameBorderMode(0);

   double bins[] = {-2.5,-2.1,-1.55,-1.05,-0.6,-0.1,0.1,0.6,1.05,1.55,2.1,2.5};
   TH1F* hdummy = new TH1F("hdummy","hdummy",11,bins);
   hdummy->SetMinimum(0.2);
   hdummy->SetMaximum(1.09);
   hdummy->GetXaxis()->SetTitle("#eta^{#mu}");
   hdummy->GetXaxis()->SetLabelSize(0.07);
   hdummy->GetXaxis()->SetTitleSize(0.1);
   hdummy->GetXaxis()->SetTitleOffset(0.87);
   hdummy->GetYaxis()->SetTitle("#epsilon");
   hdummy->GetYaxis()->SetLabelSize(0.07);
   hdummy->GetYaxis()->SetTitleSize(0.09);
   hdummy->GetYaxis()->SetTitleOffset(0.44);
   hdummy->Draw();
   teff->Draw("pesame");

   c_2->Modified();
}

void plotTagAndProbe(){
    bool doData=true;
    bool doMc=false;

    TFile* _file0 = TFile::Open("tagAndProbe_data.08.25.2014.root");
    const int nPtBins = 1;
    const int nEtaBins = 11;
    const int nCentBins = 5;
    std::map<std::string,TGraphAsymmErrors*> graphs_pos_A, graphs_pos_B, graphs_neg_A, graphs_neg_B;
    std::map<std::string,TCanvas*> canvases;
    std::map<std::string,TEfficiency*> m_teff;

    bool isEmptyPosA=true, isEmptyPosB=true, isEmptyNegA=true, isEmptyNegB=true;
    int marker[nCentBins] = {kRed, kBlue, kMagenta, kGreen, kOrange};
    int style[nCentBins] = {20,21,22,23,33};
    TString centrality[nCentBins] = {"0-5%","5-10%","10-20%","20-40%","40-90%"};
    TH1F *h = new TH1F("h","h",100,-3,3);
    h->GetYaxis()->SetRangeUser(0.0,1.1);
    h->GetXaxis()->SetTitle("#eta^{#mu}_{lab}");
    h->GetYaxis()->SetTitle("#epsilon");

    TGraphAsymmErrors *gr=NULL;
    TEfficiency* teff=NULL;
/*    TLegend *leg = new TLegend(0.6623563,0.2457627,0.9324713,0.4872881,NULL,"brNDC");

    TH1F* htotal = new TH1F("htotal","htotal",100,-5.0,5.0);
    TH1F* hpassed = new TH1F("hpassed","hpassed",100,-5.0,5.0);
    TLegend *leg2 = new TLegend(0.6422414,0.7860169,0.9123563,0.9364407,NULL,"brNDC");
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.05);
    leg2->SetLineColor(1);
    leg2->SetLineStyle(1);
    leg2->SetLineWidth(1);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(1001);
    TLegendEntry *entry=leg2->AddEntry("htotal","total","f");

    Int_t ci;      // for color index setting
    ci = TColor::GetColor("#ff0000");
    entry->SetFillColor(ci);
    entry->SetFillStyle(1001);
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(2);
    entry->SetMarkerColor(1);
    entry->SetMarkerStyle(21);
    entry->SetMarkerSize(1);
    entry->SetTextFont(42);
    entry=leg2->AddEntry("hpassed","passed","f");

    ci = TColor::GetColor("#ffff00");
    entry->SetFillColor(ci);
    entry->SetFillStyle(1001);
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(2);
    entry->SetMarkerColor(1);
    entry->SetMarkerStyle(21);
    entry->SetMarkerSize(1);
    entry->SetTextFont(42);
*/
    TCanvas* c = NULL;
    for(int ipt=0; ipt<nPtBins; ++ipt){
        for(int icent=0; icent<nCentBins; ++icent){

            stringstream cNamePosA, cNameNegA;
            stringstream cNamePosB, cNameNegB;
            cNamePosA << "canv_TP_pos_A_cent"<<icent;
            cNameNegA << "canv_TP_neg_A_cent"<<icent;
            cNamePosB << "canv_TP_pos_B_cent"<<icent;
            cNameNegB << "canv_TP_neg_B_cent"<<icent;

            c = new TCanvas((cNamePosA.str()+"_ID").c_str(),(cNamePosA.str()+"_ID").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNamePosA.str()+"_ID",c));
            c = new TCanvas((cNameNegA.str()+"_ID").c_str(),(cNameNegA.str()+"_ID").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNameNegA.str()+"_ID",c));

            c = new TCanvas((cNamePosA.str()+"_MS").c_str(),(cNamePosA.str()+"_MS").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNamePosA.str()+"_MS",c));
            c = new TCanvas((cNameNegA.str()+"_MS").c_str(),(cNameNegA.str()+"_MS").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNameNegA.str()+"_MS",c));

            c = new TCanvas((cNamePosA.str()+"_EF_mu8").c_str(),(cNamePosA.str()+"_EF_mu8").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNamePosA.str()+"_EF_mu8",c));
            c = new TCanvas((cNameNegA.str()+"_EF_mu8").c_str(),(cNameNegA.str()+"_EF_mu8").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNameNegA.str()+"_EF_mu8",c));

            c = new TCanvas((cNamePosB.str()+"_ID").c_str(),(cNamePosB.str()+"_ID").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNamePosB.str()+"_ID",c));
            c = new TCanvas((cNameNegB.str()+"_ID").c_str(),(cNameNegB.str()+"_ID").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNameNegB.str()+"_ID",c));

            c = new TCanvas((cNamePosB.str()+"_MS").c_str(),(cNamePosB.str()+"_MS").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNamePosB.str()+"_MS",c));
            c = new TCanvas((cNameNegB.str()+"_MS").c_str(),(cNameNegB.str()+"_MS").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNameNegB.str()+"_MS",c));

            c = new TCanvas((cNamePosB.str()+"_EF_mu8").c_str(),(cNamePosB.str()+"_EF_mu8").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNamePosB.str()+"_EF_mu8",c));
            c = new TCanvas((cNameNegB.str()+"_EF_mu8").c_str(),(cNameNegB.str()+"_EF_mu8").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNameNegB.str()+"_EF_mu8",c));

            std::stringstream name_pos_A_ID, name_pos_B_ID;
            std::stringstream name_neg_A_ID, name_neg_B_ID;
            std::stringstream name_pos_A_MS, name_pos_B_MS;
            std::stringstream name_neg_A_MS, name_neg_B_MS;
            std::stringstream name_pos_A_EF_mu8, name_pos_B_EF_mu8;
            std::stringstream name_neg_A_EF_mu8, name_neg_B_EF_mu8;

            std::string type;
            if(doData) type="data";
            else if(doMc) type="mc";
            // ID 
            // Period A
            name_pos_A_ID << type << "_A_ID_hPosEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_pos_A_ID.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_pos_A_ID.str(),teff));

            name_neg_A_ID << type << "_A_ID_hNegEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_neg_A_ID.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_neg_A_ID.str(),teff));

            // Period B 
            name_pos_B_ID << type << "_B_ID_hPosEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_pos_B_ID.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_pos_B_ID.str(),teff));

            name_neg_B_ID << type << "_B_ID_hNegEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_neg_B_ID.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_neg_B_ID.str(),teff));

            // MS 
            // Period A
            name_pos_A_MS << type << "_A_MS_hPosEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_pos_A_MS.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_pos_A_MS.str(),teff));

            name_neg_A_MS << type << "_A_MS_hNegEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_neg_A_MS.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_neg_A_MS.str(),teff));

            // Period B 
            name_pos_B_MS << type << "_B_MS_hPosEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_pos_B_MS.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_pos_B_MS.str(),teff));

            name_neg_B_MS << type << "_B_MS_hNegEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_neg_B_MS.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_neg_B_MS.str(),teff));

            // EF_mu8 
            // Period A
            name_pos_A_EF_mu8 << type << "_A_EF_mu8_hPosEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_pos_A_EF_mu8.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_pos_A_EF_mu8.str(),teff));

            name_neg_A_EF_mu8 << type << "_A_EF_mu8_hNegEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_neg_A_EF_mu8.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_neg_A_EF_mu8.str(),teff));

            // Period B 
            name_pos_B_EF_mu8 << type << "_B_EF_mu8_hPosEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_pos_B_EF_mu8.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_pos_B_EF_mu8.str(),teff));

            name_neg_B_EF_mu8 << type << "_B_EF_mu8_hNegEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_neg_B_EF_mu8.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_neg_B_EF_mu8.str(),teff));

            // ID
            plotVerbose(m_teff[name_pos_A_ID.str()], canvases[cNamePosA.str()+"_ID"],"ID(PeriodA), #mu^{+} "+centrality[icent]);
            plotVerbose(m_teff[name_pos_B_ID.str()], canvases[cNamePosB.str()+"_ID"],"ID(PeriodB), #mu^{+} "+centrality[icent]);
            plotVerbose(m_teff[name_neg_A_ID.str()], canvases[cNameNegA.str()+"_ID"],"ID(PeriodA), #mu^{-} "+centrality[icent]);
            plotVerbose(m_teff[name_neg_B_ID.str()], canvases[cNameNegB.str()+"_ID"],"ID(PeriodB), #mu^{-} "+centrality[icent]);

            // MS 
            plotVerbose(m_teff[name_pos_A_MS.str()], canvases[cNamePosA.str()+"_MS"],"MS(PeriodA), #mu^{+} "+centrality[icent]);
            plotVerbose(m_teff[name_pos_B_MS.str()], canvases[cNamePosB.str()+"_MS"],"MS(PeriodB), #mu^{+} "+centrality[icent]);
            plotVerbose(m_teff[name_neg_A_MS.str()], canvases[cNameNegA.str()+"_MS"],"MS(PeriodA), #mu^{-} "+centrality[icent]);
            plotVerbose(m_teff[name_neg_B_MS.str()], canvases[cNameNegB.str()+"_MS"],"MS(PeriodB), #mu^{-} "+centrality[icent]);

            // EF_mu8 
            plotVerbose(m_teff[name_pos_A_EF_mu8.str()], canvases[cNamePosA.str()+"_EF_mu8"],"EF_mu8(PeriodA), #mu^{+} "+centrality[icent]);
            plotVerbose(m_teff[name_pos_B_EF_mu8.str()], canvases[cNamePosB.str()+"_EF_mu8"],"EF_mu8(PeriodB), #mu^{+} "+centrality[icent]);
            plotVerbose(m_teff[name_neg_A_EF_mu8.str()], canvases[cNameNegA.str()+"_EF_mu8"],"EF_mu8(PeriodA), #mu^{-} "+centrality[icent]);
            plotVerbose(m_teff[name_neg_B_EF_mu8.str()], canvases[cNameNegB.str()+"_EF_mu8"],"EF_mu8(PeriodB), #mu^{-} "+centrality[icent]);

        }//icent
    }//ipt

    std::map<std::string,TCanvas*>::iterator it;
    for(it=canvases.begin(); it!=canvases.end(); ++it){

        std::string saveDir = "plots/";
        saveDir+=it->first; 
        it->second->Print((saveDir+".C").c_str());
        it->second->Print((saveDir+".png").c_str());
        it->second->Print((saveDir+".pdf").c_str());
        
    }
}
