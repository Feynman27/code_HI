#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath> 
#include <sstream>
#include <fstream>


TString format(float value) {
        std::stringstream svalue;
        svalue  << std::setprecision(2) << value;
        return svalue.str();
}


void plotDecorCwFcalWtdDistros()
{
   TFile* fIn = TFile::Open("CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges_extrapolated.PowPy8.12.02.2013.root");
   TFile* fIn2 = TFile::Open("crossChecks/correctionFactorsW_FcalWeighted.06.12.2014.root");
   fIn->cd();
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
   etaBins.push_back(+2.50);
   const int nEtaBins = etaBins.size()-1;
   float arrBinning[] = {0.0,0.1,0.35,0.6,0.8,1.05,1.37,1.52,1.74,2.1,2.50};
   TCanvas* cPlus = new TCanvas("cPlus","cPlus",600,600);
   TCanvas* cMinus = new TCanvas("cMinus","cMinus",600,600);

   TH1F* hdummy = new TH1F("hdummy","hdummy",nEtaBins+1,arrBinning);
   hdummy->GetXaxis()->SetRangeUser(0.0,2.5);
   hdummy->GetYaxis()->SetRangeUser(0.0,1.0);
   hdummy->GetYaxis()->SetTitle("C_{W}");
   hdummy->GetXaxis()->SetTitle("|#eta_{#mu}|");
   hdummy->GetXaxis()->SetLabelSize(0.05);

   TLatex *   texP = new TLatex(0.5402685,0.3041958,"#mu^{+}");
   texP->SetTextSize(0.1038462);
   TLatex *   texN = new TLatex(0.5402685,0.3041958,"#mu^{-}");
   texN->SetTextSize(0.1038462);
   TLatex *tex = new TLatex(0.5989933,0.8776224,"ATLAS Internal");
   tex->SetTextFont(72);

   cPlus->cd();
   hdummy->Draw();
   cMinus->cd();
   hdummy->Draw();

   int markerColor[6] = {kRed,kPink+10,kGreen,kViolet,kCyan,kBlue};
   int lineStyle[6] = {kDashed,kDashed,kDashed,kDashed,kDashed,kDashed};
   int markerStyle[6] = {20,21,22,23,29,33};
    
    TGraphErrors* grPlus=0,*grPlus2=0;
    TGraphErrors* grMinus=0,*grMinus2=0;
    TLegend *leg = new TLegend(0.1744966,0.1818182,0.4630872,0.451049,NULL,"brNDC");
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextSize(0.03496503);
    TLegend* leg2 = new TLegend(0.513, 0.6, 0.9295, 0.8636);
    leg2->SetTextFont(gStyle->GetTextFont());
    leg2->SetTextSize(gStyle->GetTextSize());
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);


    TString sCentBin[]={"0-5%","5-10%","10-15%","15-20%","20-40%","40-80%"}

    for(int icent=0; icent<6; ++icent){
        TString sCentPlus = "grCwEtaDistroPlusCent"; sCentPlus+=icent; 
        TString sCentMinus = "grCwEtaDistroMinusCent"; sCentMinus+=icent; 

        fIn->cd();
        std::cout << "Opening " << sCentPlus << std::endl;
        grPlus= (TGraphErrors*) fIn->Get(sCentPlus);
        if(grPlus==0) break;
        std::cout << "Opening " << sCentMinus << std::endl;
        grMinus = (TGraphErrors*)fIn->Get(sCentMinus);
        if(grMinus==0) break;

        grPlus->SetMarkerColor(markerColor[icent]);
        grPlus->SetLineColor(markerColor[icent]);
        grMinus->SetMarkerColor(markerColor[icent]);
        grMinus->SetLineColor(markerColor[icent]);
        grPlus->SetMarkerStyle(markerStyle[icent]);
        grMinus->SetMarkerStyle(markerStyle[icent]);
        
        if(icent==0||icent==4||icent==5){
            cPlus->cd();
            grPlus->Draw("pesame");
            cMinus->cd();
            grMinus->Draw("pesame");

            TString s = sCentBin[icent]; 
            leg->AddEntry(grPlus,s,"pe");
            //leg2->AddEntry(grMinus,s,"pe");
            cPlus->cd(); 
            leg->Draw();
            //leg2->Draw();
            cMinus->cd(); 
            leg->Draw();
            //leg2->Draw();
        }

        fIn2->cd();
        std::cout << "Opening " << sCentPlus << std::endl;
        grPlus2= (TGraphErrors*) fIn2->Get(sCentPlus);
        if(grPlus2==0) break;
        std::cout << "Opening " << sCentMinus << std::endl;
        grMinus2 = (TGraphErrors*)fIn2->Get(sCentMinus);
        if(grMinus2==0) break;

        grPlus2->SetMarkerColor(markerColor[icent]);
        grPlus2->SetLineColor(markerColor[icent]);
        grMinus2->SetMarkerColor(markerColor[icent]);
        grMinus2->SetLineColor(markerColor[icent]);
        grPlus2->SetLineStyle(kDashed);
        grMinus2->SetLineStyle(kDashed);
        grPlus2->SetMarkerStyle(markerStyle[icent]);
        grMinus2->SetMarkerStyle(markerStyle[icent]);

        if(icent==0||icent==4||icent==5){
            cPlus->cd();
            grPlus2->Draw("pesame");
            cMinus->cd();
            grMinus2->Draw("pesame");

            TString s = sCentBin[icent]; 
            //leg->AddEntry(grPlus2,s,"pe");
            //leg2->AddEntry(grMinus2,s,"pe");
            cPlus->cd(); 
            //leg->Draw();
            //leg2->Draw();
            cMinus->cd(); 
            //leg->Draw();
            //leg2->Draw();
        }
   }
    cPlus->cd(); texP->SetNDC(); texP->Draw();
    tex->SetNDC(); tex->Draw();
    cMinus->cd(); texN->SetNDC(); texN->Draw();
    tex->SetNDC(); tex->Draw();

}



