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
#include "TROOT.h"
#include "AtlasStyle.C"

void plotVerbose(TEfficiency* teff, TCanvas* c, TString title){

    c->cd();

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
    ((TH1F*)teff->GetTotalHistogram())->GetYaxis()->SetRangeUser(0.1,ymax+250.0);
    ((TH1F*)teff->GetTotalHistogram())->GetXaxis()->SetTitle("#eta_{lab}^{#mu}");
    ((TH1F*)teff->GetTotalHistogram())->GetYaxis()->SetTitle("Events");
    
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

   /*((TH1F*)teff)->GetYaxis()->SetRangeUser(0.0,1.0);
   ((TH1F*)teff)->GetXaxis()->SetTitle("#eta_{lab}^{#mu}");
   ((TH1F*)teff)->GetYaxis()->SetTitle("C_{W}");
   */
   double bins[] = {-2.5,-2.1,-1.55,-1.05,-0.6,-0.1,0.1,0.6,1.05,1.55,2.1,2.5};
   TH1F* hdummy = new TH1F("hdummy","hdummy",11,bins);
   hdummy->SetMinimum(0.2);
   hdummy->SetMaximum(1.1);
   hdummy->GetXaxis()->SetTitle("#eta^{#mu}");
   hdummy->GetXaxis()->SetLabelSize(0.07);
   hdummy->GetXaxis()->SetTitleSize(0.1);
   hdummy->GetXaxis()->SetTitleOffset(0.87);
   hdummy->GetYaxis()->SetTitle("C_{W}");
   hdummy->GetYaxis()->SetLabelSize(0.07);
   hdummy->GetYaxis()->SetTitleSize(0.09);
   hdummy->GetYaxis()->SetTitleOffset(0.44);
   hdummy->Draw();
   teff->Draw("pesame");

   c_2->Modified();
}

int main(){
    //gROOT->LoadMacro("AtlasStyle.C");
    SetAtlasStyle();

    TFile* _file0 = TFile::Open("CorrectionFactors/efficiencies.09.24.2014.root");
    const int nPtBins = 1;
    const int nEtaBins = 11;
    const int nCentBins = 7;
    std::map<std::string,TGraphAsymmErrors*> graphs_pos_A, graphs_pos_B, graphs_neg_A, graphs_neg_B;
    std::map<std::string,TCanvas*> canvases;
    std::map<std::string,TEfficiency*> m_teff;
    TCanvas* c = new TCanvas("canv_pos_A","#mu^{+},PeriodA",700,500);
    canvases.insert(std::make_pair<std::string,TCanvas*>("canv_pos_A",c) );
    c = new TCanvas("canv_pos_B","#mu^{+},PeriodB",700,500);
    canvases.insert(std::make_pair<std::string,TCanvas*>("canv_pos_B",c) );
    c = new TCanvas("canv_neg_A","#mu^{-},PeriodA",700,500);
    canvases.insert(std::make_pair<std::string,TCanvas*>("canv_neg_A",c) );
    c = new TCanvas("canv_neg_B","#mu^{-},PeriodB",700,500);
    canvases.insert(std::make_pair<std::string,TCanvas*>("canv_neg_B",c) );

    bool isEmptyPosA=true, isEmptyPosB=true, isEmptyNegA=true, isEmptyNegB=true;
    int marker[nCentBins] = {kRed, kBlue, kMagenta, kGreen, kOrange,kViolet+7,kCyan};
    int style[nCentBins] = {20,21,22,23,33,29,34};
    TString centrality[nCentBins] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-60%","60-90%"};
    TH1F *h = new TH1F("h","h",100,-3,3);
    h->GetYaxis()->SetRangeUser(0.0,1.1);
    h->GetXaxis()->SetTitle("#eta^{#mu}_{lab}");
    h->GetYaxis()->SetTitle("C_{W}");

    TGraphAsymmErrors *gr=NULL;
    TEfficiency* teff=NULL;
    TLegend *leg = new TLegend(0.6623563,0.2457627,0.9324713,0.4872881,NULL,"brNDC");

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

    for(int ipt=0; ipt<nPtBins; ++ipt){
        for(int icent=0; icent<nCentBins; ++icent){

            std::stringstream cNamePosA, cNameNegA;
            std::stringstream cNamePosB, cNameNegB;
            cNamePosA << "canv_pos_A_cent"<<icent;
            cNameNegA << "canv_neg_A_cent"<<icent;
            cNamePosB << "canv_pos_B_cent"<<icent;
            cNameNegB << "canv_neg_B_cent"<<icent;
            c = new TCanvas((cNamePosA.str()+"_pp").c_str(),(cNamePosA.str()+"_pp").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNamePosA.str()+"_pp",c));
            c = new TCanvas((cNameNegA.str()+"_pp").c_str(),(cNameNegA.str()+"_pp").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNameNegA.str()+"_pp",c));

            c = new TCanvas((cNamePosA.str()+"_np").c_str(),(cNamePosA.str()+"_np").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNamePosA.str()+"_np",c));
            c = new TCanvas((cNameNegA.str()+"_np").c_str(),(cNameNegA.str()+"_np").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNameNegA.str()+"_np",c));

            c = new TCanvas((cNamePosB.str()+"_pp").c_str(),(cNamePosB.str()+"_pp").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNamePosB.str()+"_pp",c));
            c = new TCanvas((cNameNegB.str()+"_pp").c_str(),(cNameNegB.str()+"_pp").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNameNegB.str()+"_pp",c));

            c = new TCanvas((cNamePosB.str()+"_pn").c_str(),(cNamePosB.str()+"_pn").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNamePosB.str()+"_pn",c));
            c = new TCanvas((cNameNegB.str()+"_pn").c_str(),(cNameNegB.str()+"_pn").c_str(),700,500);
            canvases.insert(std::make_pair<std::string,TCanvas*>(cNameNegB.str()+"_pn",c));

            std::stringstream name_pos_A, name_pos_B, name_neg_A, name_neg_B;
            name_pos_A << "effPos_A_wtd_pt" << ipt << "_cent" << icent;
            gr = (TGraphAsymmErrors*)_file0->Get(name_pos_A.str().c_str());
            graphs_pos_A.insert(std::make_pair<std::string,TGraphAsymmErrors*>(name_pos_A.str(),gr) );
            name_neg_A << "effNeg_A_wtd_pt" << ipt << "_cent" << icent;
            gr = (TGraphAsymmErrors*)_file0->Get(name_neg_A.str().c_str());
            graphs_neg_A.insert(std::make_pair<std::string,TGraphAsymmErrors*>(name_neg_A.str(),gr) );
            name_pos_B << "effPos_B_wtd_pt" << ipt << "_cent" << icent;
            gr = (TGraphAsymmErrors*)_file0->Get(name_pos_B.str().c_str());
            graphs_pos_B.insert(std::make_pair<std::string,TGraphAsymmErrors*>(name_pos_B.str(),gr) );
            name_neg_B << "effNeg_B_wtd_pt" << ipt << "_cent" << icent;
            gr = (TGraphAsymmErrors*)_file0->Get(name_neg_B.str().c_str());
            graphs_neg_B.insert(std::make_pair<std::string,TGraphAsymmErrors*>(name_neg_B.str(),gr) );

            std::stringstream name_pos_A_pp, name_pos_B_pp, name_pos_A_pn, name_pos_B_pn;
            std::stringstream name_neg_A_pp, name_neg_B_pp, name_neg_A_pn, name_neg_B_pn;

            name_pos_A_pp << "pp_A_Cw_hPosEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_pos_A_pp.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_pos_A_pp.str(),teff));

            name_pos_A_pn << "np_A_Cw_hPosEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_pos_A_pn.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_pos_A_pn.str(),teff));

            name_neg_A_pp << "pp_A_Cw_hNegEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_neg_A_pp.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_neg_A_pp.str(),teff));

            name_neg_A_pn << "np_A_Cw_hNegEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_neg_A_pn.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_neg_A_pn.str(),teff));

            name_pos_B_pp << "pp_B_Cw_hPosEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_pos_B_pp.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_pos_B_pp.str(),teff));

            name_pos_B_pn << "pn_B_Cw_hPosEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_pos_B_pn.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_pos_B_pn.str(),teff));

            name_neg_B_pp << "pp_B_Cw_hNegEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_neg_B_pp.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_neg_B_pp.str(),teff));

            name_neg_B_pn << "pn_B_Cw_hNegEff_pt"<<ipt<<"_cent"<<icent;
            teff=(TEfficiency*)_file0->Get(name_neg_B_pn.str().c_str());
            m_teff.insert(std::make_pair<std::string,TEfficiency*>(name_neg_B_pn.str(),teff));


            plotVerbose(m_teff[name_pos_A_pp.str()], canvases[cNamePosA.str()+"_pp"],"pp(PeriodA), #mu^{+} "+centrality[icent]);
            plotVerbose(m_teff[name_pos_A_pn.str()], canvases[cNamePosA.str()+"_np"],"np(PeriodA), #mu^{+} "+centrality[icent]);

            plotVerbose(m_teff[name_neg_A_pp.str()], canvases[cNameNegA.str()+"_pp"],"pp(PeriodA), #mu^{-} "+centrality[icent]);
            plotVerbose(m_teff[name_neg_A_pn.str()], canvases[cNameNegA.str()+"_np"],"np(PeriodA), #mu^{-} "+centrality[icent]);

            plotVerbose(m_teff[name_pos_B_pp.str()], canvases[cNamePosB.str()+"_pp"],"pp(PeriodB), #mu^{+} "+centrality[icent]);
            plotVerbose(m_teff[name_pos_B_pn.str()], canvases[cNamePosB.str()+"_pn"],"pn(PeriodB), #mu^{+} "+centrality[icent]);

            plotVerbose(m_teff[name_neg_B_pp.str()], canvases[cNameNegB.str()+"_pp"],"pp(PeriodB), #mu^{-} "+centrality[icent]);
            plotVerbose(m_teff[name_neg_B_pn.str()], canvases[cNameNegB.str()+"_pn"],"pn(PeriodB), #mu^{-} "+centrality[icent]);


            canvases["canv_pos_A"]->cd();
            if(isEmptyPosA) {
                graphs_pos_A[name_pos_A.str()]->SetMarkerColor(marker[icent]);
                graphs_pos_A[name_pos_A.str()]->SetLineColor(marker[icent]);
                graphs_pos_A[name_pos_A.str()]->SetMarkerStyle(style[icent]);
                h->Draw();
                graphs_pos_A[name_pos_A.str()]->Draw("pesame");
                TLatex *   tex = new TLatex(0.3461371,0.8628191,"Period A");
                tex->SetNDC();
                tex->SetTextFont(42);
                tex->SetLineWidth(2);
                tex->Draw();
                tex = new TLatex(0.1933654,0.8628191,"W^{+}#rightarrow#mu^{+}");
                tex->SetNDC();
                tex->SetTextFont(42);
                tex->SetLineWidth(2);
                tex->Draw();

                leg->SetTextFont(gStyle->GetTextFont());
                leg->SetTextSize(gStyle->GetTextSize());
                leg->SetBorderSize(0);
                leg->SetFillColor(0);
                leg->AddEntry(graphs_pos_A[name_pos_A.str()],centrality[icent],"pl");

                isEmptyPosA = false;
            }
            else { 
                graphs_pos_A[name_pos_A.str()]->SetMarkerColor(marker[icent]);
                graphs_pos_A[name_pos_A.str()]->SetLineColor(marker[icent]);
                graphs_pos_A[name_pos_A.str()]->SetMarkerStyle(style[icent]);
                graphs_pos_A[name_pos_A.str()]->Draw("pesame");
                leg->AddEntry(graphs_pos_A[name_pos_A.str()],centrality[icent],"pl");
            }
            leg->Draw();

            canvases["canv_pos_B"]->cd();
            if(isEmptyPosB) { 
                graphs_pos_B[name_pos_B.str()]->SetMarkerColor(marker[icent]);
                graphs_pos_B[name_pos_B.str()]->SetLineColor(marker[icent]);
                graphs_pos_B[name_pos_B.str()]->SetMarkerStyle(style[icent]);
                h->Draw();
                graphs_pos_B[name_pos_B.str()]->Draw("pesame");
                TLatex *   tex = new TLatex(0.3461371,0.8628191,"Period B");
                tex->SetNDC();
                tex->SetTextFont(42);
                tex->SetLineWidth(2);
                tex->Draw();
                tex = new TLatex(0.1933654,0.8628191,"W^{+}#rightarrow#mu^{+}");
                tex->SetNDC();
                tex->SetTextFont(42);
                tex->SetLineWidth(2);
                tex->Draw();
                isEmptyPosB=false;
            }
            else {
                graphs_pos_B[name_pos_B.str()]->SetMarkerColor(marker[icent]);
                graphs_pos_B[name_pos_B.str()]->SetLineColor(marker[icent]);
                graphs_pos_B[name_pos_B.str()]->SetMarkerStyle(style[icent]);
                graphs_pos_B[name_pos_B.str()]->Draw("pesame");
            
            }
            leg->Draw();
            canvases["canv_neg_A"]->cd();
            if(isEmptyNegA){
                graphs_neg_A[name_neg_A.str()]->SetMarkerColor(marker[icent]);
                graphs_neg_A[name_neg_A.str()]->SetLineColor(marker[icent]);
                graphs_neg_A[name_neg_A.str()]->SetMarkerStyle(style[icent]);
                h->Draw();
                graphs_neg_A[name_neg_A.str()]->Draw("pesame");
                TLatex *   tex = new TLatex(0.3461371,0.8628191,"Period A");
                tex->SetNDC();
                tex->SetTextFont(42);
                tex->SetLineWidth(2);
                tex->Draw();
                tex = new TLatex(0.1933654,0.8628191,"W^{-}#rightarrow#mu^{-}");
                tex->SetNDC();
                tex->SetTextFont(42);
                tex->SetLineWidth(2);
                tex->Draw();
                isEmptyNegA=false;
            }
            else {
                graphs_neg_A[name_neg_A.str()]->SetMarkerColor(marker[icent]);
                graphs_neg_A[name_neg_A.str()]->SetLineColor(marker[icent]);
                graphs_neg_A[name_neg_A.str()]->SetMarkerStyle(style[icent]);
                graphs_neg_A[name_neg_A.str()]->Draw("pesame");
            }
            leg->Draw();
            canvases["canv_neg_B"]->cd();
            if(isEmptyNegB){
                graphs_neg_B[name_neg_B.str()]->SetMarkerColor(marker[icent]);
                graphs_neg_B[name_neg_B.str()]->SetLineColor(marker[icent]);
                graphs_neg_B[name_neg_B.str()]->SetMarkerStyle(style[icent]);
                h->Draw();
                graphs_neg_B[name_neg_B.str()]->Draw("pesame");
                TLatex *   tex = new TLatex(0.3461371,0.8628191,"Period B");
                tex->SetNDC();
                tex->SetTextFont(42);
                tex->SetLineWidth(2);
                tex->Draw();
                tex = new TLatex(0.1933654,0.8628191,"W^{-}#rightarrow#mu^{-}");
                tex->SetNDC();
                tex->SetTextFont(42);
                tex->SetLineWidth(2);
                tex->Draw();
                isEmptyNegB=false;
            }
            else {
                graphs_neg_B[name_neg_B.str()]->SetMarkerColor(marker[icent]);
                graphs_neg_B[name_neg_B.str()]->SetLineColor(marker[icent]);
                graphs_neg_B[name_neg_B.str()]->SetMarkerStyle(style[icent]);
                graphs_neg_B[name_neg_B.str()]->Draw("pesame");
            }
            leg->Draw();
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
    
    return 0;
}
