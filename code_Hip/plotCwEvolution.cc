#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"

#include <iostream>
#include <map>
#include <sstream>

void plot(TCanvas* canvas, std::multimap<std::string,TGraphAsymmErrors*> map, std::string grName, TH1F* hdummy, TString period, TString charge, TString centrality){

    std::multimap<std::string,TGraphAsymmErrors*>::iterator it = map.find(grName);
    if(it!=map.end()){

        canvas->cd();
        hdummy->Draw();
        TLegend *leg = new TLegend(0.1867816,0.2055085,0.4568966,0.4533898,NULL,"brNDC");
        leg->SetBorderSize(0);
        leg->SetTextSize(0.04);
        leg->SetTextFont(42);
        leg->SetLineColor(1);
        leg->SetLineStyle(1);
        leg->SetLineWidth(1);
        leg->SetFillColor(0);
        leg->SetFillStyle(1001);
        TString legLbl[]={"#epsilon_{rec}","#epsilon_{rec}#epsilon_{sel}(isolation)","#epsilon_{rec}#epsilon_{sel}(Z veto)","#epsilon_{rec}#epsilon_{sel}#epsilon_{trig}"};
        int mkSty[]={20,21,23,33};
        int mkCol[]={kGreen,kRed,kBlue,kBlack};
        TString label = period;
        label+=", "; label+=charge; label+=", "; label+=centrality;
        TLatex * tex = new TLatex(0.4891954,0.2944915,label);
        tex->SetNDC();
        tex->SetTextFont(42);
        tex->SetTextSize(0.04);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.4729885,0.2372881,"p_{T}^{#mu}>25GeV, p_{T}^{#nu}>25GeV, m_{T}>40 GeV");
        tex->SetNDC();
        tex->SetTextFont(42);
        tex->SetTextSize(0.03);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.2068966,0.875,"ATLAS Simulation");
        tex->SetNDC();
        tex->SetTextFont(72);
        tex->SetLineWidth(2);
        tex->Draw();
        tex = new TLatex(0.6342529,0.875,"Internal");
        tex->SetNDC();
        tex->SetTextFont(42);
        tex->SetLineWidth(2);
        tex->Draw();

        // Number of graphs named grName
        size_t nGraphs = map.count(grName);
        for(size_t igr = 0; igr<nGraphs; ++igr){
            std::cout << it->first << std::endl;
            it->second->SetMarkerStyle(mkSty[igr]);
            it->second->SetMarkerColor(mkCol[igr]);
            it->second->SetLineColor(mkCol[igr]);
            if(igr==nGraphs-1) it->second->SetMarkerSize(1.6);
            it->second->Draw("pe same");
            leg->AddEntry(it->second,legLbl[igr],"pl");
            ++it;
        }
        leg->Draw();
    }
    else std::cout << "TGraph not found in map." << std::endl;

    std::stringstream ss;
    std::string s;
    ss << canvas->GetName();
    ss >> s;
    TString save = "plots/"+s;
    canvas->Print(save+".C");
    canvas->Print(save+".pdf");
    canvas->Print(save+".png");
}

void plotCwEvolution(){

    TFile* _file0 = TFile::Open("CorrectionFactors/efficiencies.loose.root");
    TFile* _file1 = TFile::Open("CorrectionFactors/efficiencies.isolation.root");
    TFile* _file2 = TFile::Open("CorrectionFactors/efficiencies.isolation.ZVeto.root");
    TFile* _file3 = TFile::Open("CorrectionFactors/efficiencies.isolation.ZVeto.EFmatch.root");

    const int nCentBins=5;
    TString centLabels[nCentBins] = {"0-5%","5-10%","10-20%","20-40%","40-90%"};
    std::map<std::string,TCanvas*> m_canvas;
    std::multimap<std::string,TGraphAsymmErrors*> m_graphs;
    TCanvas* c=NULL;
    TGraphAsymmErrors* gr=NULL;
    //std::vector<TGraphAsymmErrors> v_graphs;
    // Initialize TCanvases
    for(int icent=0; icent<nCentBins; ++icent){
        std::stringstream cName;
        // mu+
        cName << "cPosA_cent" << icent ;
        c = new TCanvas(cName.str().c_str(),cName.str().c_str(),600,600);
        m_canvas.insert(std::make_pair<std::string,TCanvas*>(cName.str(),c));
        cName.str(std::string());

        // mu-
        cName << "cNegA_cent" << icent ;
        c = new TCanvas(cName.str().c_str(),cName.str().c_str(),600,600);
        m_canvas.insert(std::make_pair<std::string,TCanvas*>(cName.str(),c));
        cName.str(std::string());

        // mu+
        cName << "cPosB_cent" << icent ;
        c = new TCanvas(cName.str().c_str(),cName.str().c_str(),600,600);
        m_canvas.insert(std::make_pair<std::string,TCanvas*>(cName.str(),c));
        cName.str(std::string());

        // mu-
        cName << "cNegB_cent" << icent ;
        c = new TCanvas(cName.str().c_str(),cName.str().c_str(),600,600);
        m_canvas.insert(std::make_pair<std::string,TCanvas*>(cName.str(),c));
        cName.str(std::string());

        // Graphs
        
        std::stringstream grName;
        //////
        // mu+
        //////
        // Period A
        // epsilon rec
        grName<<"effPos_A_wtd_pt0_cent"<<icent;
        gr = (TGraphAsymmErrors*)_file0->Get(grName.str().c_str());
        m_graphs.insert(std::make_pair<std::string,TGraphAsymmErrors*>(grName.str(),gr));

        // epsilon sel(isolation) 
        //grName<<"effPos_A_wtd_pt0_cent"<<icent;
        gr = (TGraphAsymmErrors*)_file1->Get(grName.str().c_str());
        m_graphs.insert(std::make_pair<std::string,TGraphAsymmErrors*>(grName.str(),gr));

        // epsilon sel(Zveto+isolation) 
        //grName<<"effPos_A_wtd_pt0_cent"<<icent;
        gr = (TGraphAsymmErrors*)_file2->Get(grName.str().c_str());
        m_graphs.insert(std::make_pair<std::string,TGraphAsymmErrors*>(grName.str(),gr));

        // epsilon trig(EF_mu8) 
        //grName<<"effPos_A_wtd_pt0_cent"<<icent;
        gr = (TGraphAsymmErrors*)_file3->Get(grName.str().c_str());
        m_graphs.insert(std::make_pair<std::string,TGraphAsymmErrors*>(grName.str(),gr));

        //////
        // mu-
        //////
        // Period A
        // epsilon rec
        grName.str(std::string());
        grName<<"effNeg_A_wtd_pt0_cent"<<icent;
        gr = (TGraphAsymmErrors*)_file0->Get(grName.str().c_str());
        m_graphs.insert(std::make_pair<std::string,TGraphAsymmErrors*>(grName.str(),gr));
        grName.str(std::string());

        // epsilon sel(isolation) 
        grName<<"effNeg_A_wtd_pt0_cent"<<icent;
        gr = (TGraphAsymmErrors*)_file1->Get(grName.str().c_str());
        m_graphs.insert(std::make_pair<std::string,TGraphAsymmErrors*>(grName.str(),gr));
        grName.str(std::string());

        // epsilon sel(Zveto+isolation) 
        grName<<"effNeg_A_wtd_pt0_cent"<<icent;
        gr = (TGraphAsymmErrors*)_file2->Get(grName.str().c_str());
        m_graphs.insert(std::make_pair<std::string,TGraphAsymmErrors*>(grName.str(),gr));
        grName.str(std::string());

        // epsilon trig(EF_mu8) 
        grName<<"effNeg_A_wtd_pt0_cent"<<icent;
        gr = (TGraphAsymmErrors*)_file3->Get(grName.str().c_str());
        m_graphs.insert(std::make_pair<std::string,TGraphAsymmErrors*>(grName.str(),gr));
        grName.str(std::string());

        //////
        // mu+
        //////
        // Period B
        // epsilon rec
        grName<<"effPos_B_wtd_pt0_cent"<<icent;
        gr = (TGraphAsymmErrors*)_file0->Get(grName.str().c_str());
        m_graphs.insert(std::make_pair<std::string,TGraphAsymmErrors*>(grName.str(),gr));
        grName.str(std::string());

        // epsilon sel(isolation) 
        grName<<"effPos_B_wtd_pt0_cent"<<icent;
        gr = (TGraphAsymmErrors*)_file1->Get(grName.str().c_str());
        m_graphs.insert(std::make_pair<std::string,TGraphAsymmErrors*>(grName.str(),gr));
        grName.str(std::string());

        // epsilon sel(Zveto+isolation) 
        grName<<"effPos_B_wtd_pt0_cent"<<icent;
        gr = (TGraphAsymmErrors*)_file2->Get(grName.str().c_str());
        m_graphs.insert(std::make_pair<std::string,TGraphAsymmErrors*>(grName.str(),gr));
        grName.str(std::string());

        // epsilon trig(EF_mu8) 
        grName<<"effPos_B_wtd_pt0_cent"<<icent;
        gr = (TGraphAsymmErrors*)_file3->Get(grName.str().c_str());
        m_graphs.insert(std::make_pair<std::string,TGraphAsymmErrors*>(grName.str(),gr));
        grName.str(std::string());

        //////
        // mu-
        //////
        // Period A
        // epsilon rec
        grName<<"effNeg_B_wtd_pt0_cent"<<icent;
        gr = (TGraphAsymmErrors*)_file0->Get(grName.str().c_str());
        m_graphs.insert(std::make_pair<std::string,TGraphAsymmErrors*>(grName.str(),gr));
        grName.str(std::string());

        // epsilon sel(isolation) 
        grName<<"effNeg_B_wtd_pt0_cent"<<icent;
        gr = (TGraphAsymmErrors*)_file1->Get(grName.str().c_str());
        m_graphs.insert(std::make_pair<std::string,TGraphAsymmErrors*>(grName.str(),gr));
        grName.str(std::string());

        // epsilon sel(Zveto+isolation) 
        grName<<"effNeg_B_wtd_pt0_cent"<<icent;
        gr = (TGraphAsymmErrors*)_file2->Get(grName.str().c_str());
        m_graphs.insert(std::make_pair<std::string,TGraphAsymmErrors*>(grName.str(),gr));
        grName.str(std::string());

        // epsilon trig(EF_mu8) 
        grName<<"effNeg_B_wtd_pt0_cent"<<icent;
        gr = (TGraphAsymmErrors*)_file3->Get(grName.str().c_str());
        m_graphs.insert(std::make_pair<std::string,TGraphAsymmErrors*>(grName.str(),gr));
        grName.str(std::string());
    }

    TH1F* hdummy = new TH1F("hdummy","hdummy",50,-2.5,2.5);
    hdummy->GetYaxis()->SetRangeUser(0.0,1.1);
    hdummy->GetYaxis()->SetTitle("#epsilon");
    hdummy->GetXaxis()->SetTitle("#eta_{lab}^{#mu}");
    // plot 
    for(int icent=0; icent<nCentBins; ++icent){
        std::stringstream cName;
        std::stringstream grName;
    
        // mu+ PeriodA
        cName << "cPosA_cent" << icent ;
        grName<<"effPos_A_wtd_pt0_cent"<<icent;
        plot(m_canvas[cName.str()],m_graphs,grName.str(),hdummy,"Period A","#mu^{+}",centLabels[icent]);
        grName.str(std::string()); cName.str(std::string());
        // mu- PeriodA
        cName << "cNegA_cent" << icent ;
        grName<<"effNeg_A_wtd_pt0_cent"<<icent;
        plot(m_canvas[cName.str()],m_graphs,grName.str(),hdummy,"Period A","#mu^{-}",centLabels[icent]);
        grName.str(std::string()); cName.str(std::string());
        // mu+ PeriodB
        cName << "cPosB_cent" << icent ;
        grName<<"effPos_B_wtd_pt0_cent"<<icent;
        plot(m_canvas[cName.str()],m_graphs,grName.str(),hdummy,"Period B","#mu^{+}",centLabels[icent]);
        grName.str(std::string()); cName.str(std::string());
        // mu- PeriodB
        cName << "cNegB_cent" << icent ;
        grName<<"effNeg_B_wtd_pt0_cent"<<icent;
        plot(m_canvas[cName.str()],m_graphs,grName.str(),hdummy,"Period B","#mu^{-}",centLabels[icent]);
        grName.str(std::string()); cName.str(std::string());
    
    }
}
