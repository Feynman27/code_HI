#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"

#include <iostream>

void plotDecorCorrectionFactorsMptSystematics(){

    TFile *_fileNom = TFile::Open("CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges.PowPy8.10.03.2013.root");
    TFile *_fileMptUp = TFile::Open("systematics/correctionFactorsW_4GeVMpt_systematics.10.29.2013.root");
    TFile *_fileMptDown = TFile::Open("systematics/correctionFactorsW_2GeVMpt_systematics.10.29.2013.root");
    TFile *_fileNoMpt = TFile::Open("systematics/correctionFactorsW_NoMptCut_systematics.10.29.2013.root");

    const int nCentralityBins = 6;
    TString arrCentLabel[] = {"0-5%","5-10%","10-15%","15-20%","20-40%","40-80%"};
    for(int icent=0; icent<6; ++icent){
        
        TString name = "grCwEtaDistro";
        TString namePlus = name + "PlusCent";
        namePlus+=icent;
        TString nameMinus = name + "MinusCent";
        nameMinus+=icent;

        std::cout << "Opening graph with name " << namePlus << " and "
            << nameMinus << std::endl;

        TGraphAsymmErrors* grPlusNom = (TGraphAsymmErrors*)_fileNom->Get(namePlus);
        TGraphAsymmErrors* grPlusMptUp = (TGraphAsymmErrors*)_fileMptUp->Get(namePlus);
        TGraphAsymmErrors* grPlusMptDown = (TGraphAsymmErrors*)_fileMptDown->Get(namePlus);
        TGraphAsymmErrors* grPlusNoMpt = (TGraphAsymmErrors*)_fileNoMpt->Get(namePlus);

        TGraphAsymmErrors* grMinusNom = (TGraphAsymmErrors*)_fileNom->Get(nameMinus);
        TGraphAsymmErrors* grMinusMptUp = (TGraphAsymmErrors*)_fileMptUp->Get(nameMinus);
        TGraphAsymmErrors* grMinusMptDown = (TGraphAsymmErrors*)_fileMptDown->Get(nameMinus);

        TString cNamePlus = "cPlusCent"; cNamePlus+=icent;
        TCanvas* cPlus = new TCanvas(cNamePlus,cNamePlus,600,600);
        grPlusNom->Draw("ape");
        grPlusMptUp->Draw("pesame");
        grPlusMptDown->Draw("pesame");
//        grPlusNoMpt->Draw("pesame");
        grPlusNom->GetYaxis()->SetRangeUser(0.0,1.0);
        grPlusNom->GetYaxis()->SetTitle("C_{W^{+}}");
        grPlusNom->GetXaxis()->SetTitle("|#eta_{#mu}|");
        grPlusMptUp->SetMarkerColor(kBlue);
        grPlusMptUp->SetLineColor(kBlue);
        grPlusMptDown->SetMarkerColor(kRed);
        grPlusMptDown->SetLineColor(kRed);
        grPlusNoMpt->SetMarkerStyle(kOpenCircle);

        TLegend* leg = new TLegend(0.46,0.22,0.73,0.47);
        leg->SetTextFont(gStyle->GetTextFont());
        leg->SetTextSize(gStyle->GetTextSize());
        leg->SetBorderSize(0);
        leg->SetFillColor(0);
        leg->AddEntry(grPlusNom,"#slash{p_{T}}>25(p_{T}^{trk}>3GeV)","lpe");
        leg->AddEntry(grPlusMptUp,"#slash{p_{T}}>25(p_{T}^{trk}>4 GeV)","lpe");
        leg->AddEntry(grPlusMptDown,"#slash{p_{T}}>25 (p_{T}^{trk}>2 GeV)","lpe");
//        leg->AddEntry(grPlusNoMpt,"No #slash{p_{T}} cut","lpe");
        leg->Draw();

        TLatex* tex = new TLatex(0.2,0.2,"#mu^{+}"); 
        tex->SetNDC(); 
        tex->SetTextSize(0.1);
        tex->Draw();
        tex = new TLatex(0.2,0.4,arrCentLabel[icent]); 
        tex->SetNDC(); 
        tex->Draw();

        cPlus->Print(cNamePlus+"_10_29_2013.pdf");

    }//icent
}
