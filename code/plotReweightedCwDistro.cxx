#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TList.h"
#include "TColor.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"

#include <iostream>

void plotOnCanvas(TGraphAsymmErrors* grNomPlus, TGraphAsymmErrors* grReWtPlus, TString cNamePlus, TGraphAsymmErrors* grNomMinus, TGraphAsymmErrors* grReWtMinus, TString cNameMinus, TString sCent){
    
    TCanvas* cPlus = new TCanvas(cNamePlus,cNamePlus,600,600);
    grNomPlus->SetMarkerColor(kRed);  grNomPlus->SetLineColor(kRed); 
    grReWtPlus->SetMarkerColor(kRed); grReWtPlus->SetLineColor(kRed); grReWtPlus->SetMarkerStyle(kOpenCircle);
    grNomPlus->GetYaxis()->SetRangeUser(0.0,1.0);
    grNomPlus->Draw("ape");
    grReWtPlus->Draw("pesame");

    TCanvas* cMinus = new TCanvas(cNameMinus,cNameMinus,600,600);
    grNomMinus->SetMarkerColor(kBlue);  grNomMinus->SetLineColor(kBlue); 
    grReWtMinus->SetMarkerColor(kBlue); grReWtMinus->SetLineColor(kBlue); grReWtMinus->SetMarkerStyle(kOpenCircle);
    grNomMinus->GetYaxis()->SetRangeUser(0.0,1.0);
    grNomMinus->Draw("ape");
    grReWtMinus->Draw("pesame");

    TLegend *leg = new TLegend(0.1845638,0.2045455,0.454698,0.3548951,NULL,"brNDC");
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    leg->AddEntry(grNomPlus,"Nominal","pe");
    leg->AddEntry(grReWtPlus,"ReWeighted","pe");
    cPlus->cd();
    leg->Draw();

    TLatex *   tex = new TLatex(0.5302013,0.2954545,"#mu^{+}");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.1153846);
    tex->SetLineWidth(2);
    tex->Draw();

    tex = new TLatex(0.5285235,0.222028,sCent);
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();

    leg = new TLegend(0.1845638,0.2045455,0.454698,0.3548951,NULL,"brNDC");
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    leg->AddEntry(grNomMinus,"Nominal","pe");
    leg->AddEntry(grReWtMinus,"ReWeighted","pe");
    cMinus->cd();
    leg->Draw();

    tex = new TLatex(0.5302013,0.2954545,"#mu^{-}");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.1153846);
    tex->SetLineWidth(2);
    tex->Draw();

    tex = new TLatex(0.5285235,0.222028,sCent);
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();

    cPlus->Print(cNamePlus+"_CwReWeighted.pdf");
    cMinus->Print(cNameMinus+"_CwReWeighted.pdf");
}

void plotReweightedCwDistro(){

    TString sFileNom = "CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_extrapolated.PowPy8.12.02.2013.root";
    TFile* _fileNom = new TFile(sFileNom,"read");
    if(!_fileNom->IsOpen()){
        std::cout << sFileNom << " not found! " << std::endl; exit(0);
    }
    TString sFileReWt = "crossChecks/correctionFactorsWEtaCent_reweighted.12.12.2013.root";
    TFile* _fileReWt = new TFile(sFileReWt,"read");
    if(!_fileReWt->IsOpen()){
        std::cout << sFileReWt << " not found! " << std::endl; exit(0);
    }

    TList _grNomPlus;
    TList _grReWtPlus;
    TList _grNomMinus;
    TList _grReWtMinus;

    std::vector<double> centralityBins;
    centralityBins.push_back(0.0);
    centralityBins.push_back(0.05);
    centralityBins.push_back(0.10);
    centralityBins.push_back(0.15);
    centralityBins.push_back(0.20);
    centralityBins.push_back(0.40);
    centralityBins.push_back(0.80);
    const int nCentralityBins = centralityBins.size() - 1;
    TString centLabel[] = {"0-5%","5-10%","10-15%","15-20%","20-40%","40-80%"};

    std::vector<double> etaBins;
    etaBins.push_back(0.10);
    etaBins.push_back(0.35);
    etaBins.push_back(0.6);
    etaBins.push_back(0.8);
    etaBins.push_back(1.05);
    etaBins.push_back(1.37);
    etaBins.push_back(1.52);
    etaBins.push_back(1.74);
    etaBins.push_back(2.1);
    etaBins.push_back(+2.40);
    const int nEtaBins = etaBins.size()-1;

    for(int icent=0; icent<nCentralityBins; ++icent){
    
        _fileNom->cd();
        TString graphNamePlus = "grCwEtaDistroPlusCent"; graphNamePlus+=icent;
        _grNomPlus.Add((TGraphAsymmErrors*)_fileNom->Get(graphNamePlus));
        //((TGraphAsymmErrors*)_grNomPlus.At(icent)) = (TGraphAsymmErrors*)_fileNom->Get(graphNamePlus); 

        TString graphNameMinus = "grCwEtaDistroMinusCent"; graphNameMinus+=icent;
        _grNomMinus.Add((TGraphAsymmErrors*)_fileNom->Get(graphNameMinus));
        //(TGraphAsymmErrors*)_grNomMinus.At(icent) = (TGraphAsymmErrors*)_fileNom->Get(graphNameMinus); 

        _fileReWt->cd();
        _grReWtPlus.Add((TGraphAsymmErrors*)_fileReWt->Get(graphNamePlus));
        //(TGraphAsymmErrors*)_grReWtPlus.At(icent) = (TGraphAsymmErrors*)_fileReWt->Get(graphNamePlus);

        _grReWtMinus.Add((TGraphAsymmErrors*)_fileReWt->Get(graphNameMinus));
        //(TGraphAsymmErrors*)_grReWtMinus.At(icent) = (TGraphAsymmErrors*)_fileReWt->Get(graphNameMinus);


        TString cNamePlus = "cPlus"; cNamePlus+=icent;
        TString cNameMinus = "cMinus"; cNameMinus+=icent;
        plotOnCanvas((TGraphAsymmErrors*)_grNomPlus.At(icent),(TGraphAsymmErrors*)_grReWtPlus.At(icent),cNamePlus, 
                (TGraphAsymmErrors*)_grNomMinus.At(icent),(TGraphAsymmErrors*)_grReWtMinus.At(icent),cNameMinus,centLabel[icent]);
    }//icent


    /*for(int i=0; i<_grReWtPlus.GetEntries(); i++ ){
        delete _grNomPlus.At(i);
        delete _grReWtPlus.At(i);
        delete _grNomMinus.At(i);
        delete _grReWtMinus.At(i);
    }
    */
}
