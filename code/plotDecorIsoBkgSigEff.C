#include "TFile.h"
#include "TList.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TString.h"
#include <iomanip>
#include <vector>

void plotDecorIsoBkgSigEff(){

	TFile* fIn = TFile::Open("isoEfficiencyGraphs.06.26.2013.root");

	gROOT->LoadMacro("AtlasUtils.C");
	const int trkPtCut = 8;
	TString trkPtLo[trkPtCut] = {"05","075","1","2","3","4","5","6"};
	TString titleTrkPtLo[trkPtCut] = {"0.5","0.75","1","2","3","4","5","6"};
	const unsigned int nCones = 6;
	int arrCone[nCones] = {10,15,20,30,40,50};
	TString titleCone[nCones] = {"0.1","0.15","0.2","0.3","0.4","0.5"};
	int nGraphs = trkPtCut*nCones;
	//centrality
	std::vector <float> centBins;
	centBins.push_back(0.0);
	centBins.push_back(0.05);
	centBins.push_back(0.1);
	centBins.push_back(0.15);
	centBins.push_back(0.2);
	centBins.push_back(0.4);
	centBins.push_back(0.8);

	const int nCentralityBins = centBins.size()-1;

	const int isoCutBins = 11;

	int markerColor[nCones] = {kRed,kCyan,kBlue,kGreen,kMagenta,kOrange+5};
    TCanvas *c0[nCentralityBins][trkPtCut];
    TCanvas *c1[nCentralityBins][trkPtCut];

    std::vector <double> vecIsolCut;
    vecIsolCut.push_back(0.05);
    vecIsolCut.push_back(0.1);
    vecIsolCut.push_back(0.2);
    vecIsolCut.push_back(0.3);
    vecIsolCut.push_back(0.4);
    vecIsolCut.push_back(0.5);
    vecIsolCut.push_back(0.7);
    vecIsolCut.push_back(1.0);
    vecIsolCut.push_back(5.0);
    vecIsolCut.push_back(10.0);
    vecIsolCut.push_back(100.0);
	const int nIsolCuts = vecIsolCut.size();

    char cName0[50], cName1[50];
    TList _gr0 ;
    TList _gr1 ;
    ///loop over each centrality class
    for(int icent=0; icent<nCentralityBins; icent++){


	   std::stringstream sCentLow, sCentUpp;
	   sCentLow << std::setprecision(2) << centBins.at(icent)*100;
	   sCentUpp << std::setprecision(2) << centBins.at(icent+1)*100;
	   TString sCent = sCentLow.str(); sCent+="-";sCent+=sCentUpp.str();  sCent+="%";

        ///loop over the track pt
	    for(int itrk=0; itrk<trkPtCut; itrk++){
            sprintf(cName0,"SigBkgEffCent%iTrkPt%i",icent,itrk);
            sprintf(cName1,"EffectiveSignalCent%iTrkPt%i",icent,itrk);
            c0[icent][itrk] = new TCanvas(cName0,cName0,600,600);
            c1[icent][itrk] = new TCanvas(cName1,cName1,600,600);


            TH1F* hDummy1 = new TH1F("hDummy1","hDummy1",50,0,1.1);
            TH1F* hDummy2 = new TH1F("hDummy2","hDummy2",50,0,1.1);;
            
            TLegend *leg = new TLegend(0.168,0.746,0.337,0.946,NULL,"brNDC");
                leg->SetTextFont(gStyle->GetTextFont());
                leg->SetTextSize(gStyle->GetTextSize());
                leg->SetBorderSize(0);
                leg->SetFillColor(0);

            TLegend* leg2 = new TLegend(0.427,0.75,0.584,0.939,NULL,"brNDC");
                leg2->SetTextFont(gStyle->GetTextFont());
                leg2->SetTextSize(gStyle->GetTextSize());
                leg2->SetBorderSize(0);
                leg2->SetFillColor(0);


            TLegend *leg3 = new TLegend(0.5,0.19,0.67,0.39,NULL,"brNDC");
                leg3->SetTextFont(gStyle->GetTextFont());
                leg3->SetTextSize(gStyle->GetTextSize());
                leg3->SetBorderSize(0);
                leg3->SetFillColor(0);

            TLegend* leg4 = new TLegend(0.74,0.19,0.90,0.38,NULL,"brNDC");
                leg4->SetTextFont(gStyle->GetTextFont());
                leg4->SetTextSize(gStyle->GetTextSize());
                leg4->SetBorderSize(0);
                leg4->SetFillColor(0);


            double max = -99999.0, min = 99999.0;
            for(int icone = 0; icone<nCones; ++icone){

                int index = (icent*trkPtCut+itrk)*nCones+icone;
                std::cout << index << std::endl;
                //_gr0.Add( new TGraphErrors(nIsolCuts));
                //_gr1.Add( new TGraph(nIsolCuts));

                TString sGrName = "grIsoEffTrkPt"; sGrName+=trkPtLo[itrk]; sGrName+="cent";
                    sGrName+=icent; sGrName+="ConeRadius"; sGrName+=arrCone[icone];

                ///Graph of effective signal efficiency
                TString sGrNameEffective = sGrName+"SigBkgPurity";

                std::cout << "Graph 1: " << sGrName << std::endl; 
                std::cout << "Graph 2: " << sGrNameEffective << std::endl; 
                _gr0.Add( (TGraphErrors*)fIn->Get(sGrName)  );
                _gr1.Add( (TGraph*)fIn->Get(sGrNameEffective) );

                c0[icent][itrk]->cd();
                if(icone==0) {
                    hDummy1->GetXaxis()->SetTitle("#epsilon_{S}");
                    hDummy1->GetYaxis()->SetTitle("#epsilon_{B}");

                    //((TGraphErrors*)_gr0.At(index))->Draw("ape");

	                hDummy1->GetXaxis()->SetRangeUser(0.645,1.1);
	                hDummy1->GetXaxis()->SetLabelSize(0.04);
	                hDummy1->GetYaxis()->SetRangeUser(0.0,1.1);
	                hDummy1->GetYaxis()->SetLabelSize(0.04);
                    hDummy1->Draw();

                
                    TString sTitle = "p_{T}^{trk}>"; sTitle+=titleTrkPtLo[itrk]; sTitle+="GeV";
                    TLatex* tex = new TLatex(0.678,0.215,sTitle);
                    tex->SetNDC();
                    tex->SetTextFont(42);
                    tex->SetLineWidth(2);
                    tex->Draw();
                    tex = new TLatex(0.76,0.332,sCent);
                    tex->SetNDC();
                    tex->SetTextFont(42);
                    tex->SetLineWidth(2);
                    tex->Draw();

                    //myText(0.67,0.6,kBlack,(char*)sTitle);
                    //myText(0.67,0.7,kBlack,(char*)sCent);

                }
                else ((TGraphErrors*)_gr0.At(index))->Draw("pe same");
                

                ((TGraphErrors*)_gr0.At(index))->SetMarkerColor(markerColor[icone]);
                ((TGraphErrors*)_gr0.At(index))->SetMarkerSize(1.4);

                TString legTitle = "#DeltaR<"; legTitle+=titleCone[icone];
                if(icone<3) leg->AddEntry((TGraphErrors*)_gr0.At(index),legTitle,"p");
                else leg2->AddEntry((TGraphErrors*)_gr0.At(index),legTitle,"p");

                c1[icent][itrk]->cd();
                double maxTemp = ((TGraphErrors*)_gr0.At(index))->GetHistogram()->GetMaximum();
                if(maxTemp>max) max = maxTemp;
                double minTemp = ((TGraphErrors*)_gr0.At(index))->GetHistogram()->GetMinimum();
                if(minTemp<min) min = minTemp;
                std::cout << min << " " << max << std::endl;
                if(icone==0) {
                    //((TGraph*)_gr1.At(index))->Draw("ape");
                    hDummy2->GetXaxis()->SetTitle("i_{#mu} = #frac{#Sigma p_{T}^{trk}(#Delta R, p_{T}^{trk})}{p_{T}^{#mu}}");
                    hDummy2->GetXaxis()->SetTitleSize(0.03);
                    hDummy2->GetXaxis()->SetTitleOffset(1.8);
                    hDummy2->GetXaxis()->SetLabelSize(0.04);
                    hDummy2->GetYaxis()->SetLabelSize(0.04);
                    hDummy2->GetXaxis()->SetRangeUser(0.0,1.0);
                    if(icent==0) hDummy2->GetYaxis()->SetRangeUser(0.0,1200.0);
                    else if(icent==1) hDummy2->GetYaxis()->SetRangeUser(0.0,1000.0);
                    else if(icent==2) hDummy2->GetYaxis()->SetRangeUser(0.0,900.0);
                    else if(icent==3) hDummy2->GetYaxis()->SetRangeUser(0.0,800.0);
                    else if(icent==4) hDummy2->GetYaxis()->SetRangeUser(0.0,2000.0);
                    else if(icent==5) hDummy2->GetYaxis()->SetRangeUser(0.0,700.0);
                    hDummy2->GetYaxis()->SetTitle("N_{eff}");
                    hDummy2->Draw();
//                    c1[icent][itrk]->SetLogy(1);
//                    c1[icent][itrk]->Update();

                
                    TString sTitle = "p_{T}^{trk}>"; sTitle+=titleTrkPtLo[itrk]; sTitle+="GeV";
/*                    myText(0.67,0.3,kBlack,(char*)sTitle);
                    myText(0.67,0.4,kBlack,(char*)sCent);
                    */
                    TLatex *   tex = new TLatex(0.67,0.41,sTitle);
                    tex->SetNDC();
                    tex->SetTextFont(42);
                    tex->SetLineWidth(2);
                    tex->Draw();
                    tex = new TLatex(0.75,0.48,sCent);
                    tex->SetNDC();
                    tex->SetTextFont(42);
                    tex->SetLineWidth(2);
                    tex->Draw();

                }
                else {

                    ((TGraph*)_gr1.At(index))->Draw("pe same");
                }

                ((TGraph*)_gr1.At(index))->SetMarkerColor(markerColor[icone]);
                ((TGraph*)_gr1.At(index))->SetMarkerSize(1.4);
                //TString legTitle = "#DeltaR<"; legTitle+=titleCone[icone];
                if(icone<3) leg3->AddEntry((TGraphErrors*)_gr1.At(index),legTitle,"p");
                else leg4->AddEntry((TGraphErrors*)_gr1.At(index),legTitle,"p");
            }//icone

            c0[icent][itrk]->cd();
            leg->Draw();
            leg2->Draw();

            c1[icent][itrk]->cd();
            leg3->Draw();
            leg4->Draw();

            TString saveName0 = "bkgSigEffIsolationCent"; saveName0+=icent; saveName0+="TrkPt";
            saveName0+=titleTrkPtLo[itrk]; saveName0+="_04_13_2015";  

            TString saveName1 = "effectiveSignalIsolationCent"; saveName1+=icent; saveName1+="TrkPt";
            saveName1+=titleTrkPtLo[itrk]; saveName1+="_04_13_2015";  

            c0[icent][itrk]->Print(saveName0+".pdf");
            c0[icent][itrk]->Print(saveName0+".png");
            c1[icent][itrk]->Print(saveName1+".pdf");
            c1[icent][itrk]->Print(saveName1+".png");
            if(icent==0&&itrk==4) {
                c0[icent][itrk]->Print(saveName0+".C"); c0[icent][itrk]->Print(saveName0+".root");
                c1[icent][itrk]->Print(saveName1+".C"); c1[icent][itrk]->Print(saveName1+".root");
                return 0; //hack
            }

            delete leg;
            delete leg2;
            delete leg3;
            delete leg4;
        } //itrk
    } //icent

    for(int icent=0; icent<nCentralityBins; icent++){
        for(int itrk=0; itrk<trkPtCut; itrk++){
            delete c0[icent][itrk];
            delete c1[icent][itrk];
        }
    }

}        





