{
TFile *_fQCDEta = TFile::Open("background/fractionQCDEta_9Bins.06.17.2013.root");
TFile *_fQCDCent = TFile::Open("background/fractionQCDCent_6Bins.06.17.2013.root");
TFile *_fZEta = TFile::Open("background/fractionZEta.06.25.2013.root");
TFile *_fZCent = TFile::Open("background/fractionZCent.06.25.2013.root");
TFile *_fTauBkg = TFile::Open("background/fractionTauEtaCent_9etaBins6centBins.07.30.2013.root");

gROOT->LoadMacro("AtlasUtils.C");
TH1F* hDummy0 = new TH1F("hDummy0","hDummy0",5,0.0,2.5);
TH1F* hDummy1 = new TH1F("hDummy1","hDummy1",300,0.0,1800);
TH1F* hDummy2 = new TH1F("hDummy2","hDummy2",5,0.0,2.5);
TH1F* hDummy3 = new TH1F("hDummy3","hDummy3",300,0.0,1800);

_fQCDCent->cd();
TCanvas *c0a = new TCanvas("c0a","c0a",600,600);
hDummy1->GetYaxis()->SetTitle("f_{QCD}");
hDummy1->GetYaxis()->SetRangeUser(0.0,0.07);
hDummy1->GetXaxis()->SetTitle("#LT N_{coll} #GT");
hDummy1->Draw();
TGraph* grMean = (TGraph*)meanFracQCDCent->Clone("grMean");
TGraphErrors* grWtdAvgPlus = (TGraphErrors*)wtdAvgFracQCDPlus->Clone("grWtdAvgPlus");
TGraphErrors* grWtdAvgMinus = (TGraphErrors*)wtdAvgFracQCDMinus->Clone("grWtdAvgMinus");

grMean->SetMarkerColor(kCyan); grMean->SetMarkerStyle(20); grMean->SetMarkerSize(2.5);
myText(0.37,0.52,kBlack,(char*)"#LT f_{QCD} #GT");
grMean->Draw("p same");
fractionQCDCent->Draw("pesame");
TGraphErrors* grCentErr = new TGraphErrors(6);
//(TGraphErrors*)fractionQCDCent->Clone("grCentErr");
double errSyst = fabs( (grWtdAvgPlus->GetY()[0]) - (grWtdAvgMinus->GetY()[0]));
std::cout << errSyst << std::endl;
std::cout << fractionQCDCent->GetN() << std::endl;
///% error on npart
std::vector<double> nCollErr;
nCollErr.push_back(7.7); //0-5
nCollErr.push_back(7.5); //5-10
nCollErr.push_back(7.4); //10-15
nCollErr.push_back(7.4); //15-20
nCollErr.push_back(7.3); //20-40
nCollErr.push_back(14.2); //40-80

for(int i = 0; i<fractionQCDCent->GetN(); ++i){
    double errStat = fractionQCDCent->GetEY()[i];
    double errTotal = TMath::Sqrt(TMath::Power(errStat,2)+TMath::Power(errSyst,2));
    double xPt = fractionQCDCent->GetX()[i];
    double xErr = (fractionQCDCent->GetX()[i])*nCollErr[i]/100.0;
    double yPt = fractionQCDCent->GetY()[i];
    std::cout << errTotal  << std::endl;
    grCentErr->SetPoint(i,xPt,yPt);
    grCentErr->SetPointError(i,xErr,errTotal);
}

grCentErr->SetFillStyle(3001);
//grCentErr->Draw("e2same");
myText(0.29,0.85,kBlack,(char*)"#mu^{#pm}");
myText(0.28,0.75,kBlack,(char*)"0.1<|#eta|<2.4");


TCanvas *c0b = new TCanvas("c0b","c0b",600,600);
hDummy1->Draw();
TGraphErrors* grPlus = (TGraphErrors*)fractionQCDCentPlus->Clone("grPlus");
TGraphErrors* grMinus = (TGraphErrors*)fractionQCDCentMinus->Clone("grMinus");

grPlus->SetMarkerColor(kRed);
grMinus->SetMarkerColor(kBlue);
grPlus->Draw("pesame");
grMinus->Draw("pesame");
TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
leg->SetTextFont(gStyle->GetTextFont());
leg->SetTextSize(gStyle->GetTextSize());
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(grPlus,"#mu^{+}");
leg->AddEntry(grMinus,"#mu^{-}");
leg->Draw();
myText(0.28,0.75,kBlack,(char*)"0.1<|#eta|<2.4");


_fQCDEta->cd();

TCanvas *c1 = new TCanvas("c1","c1",600,600);
hDummy0->GetYaxis()->SetTitle("f_{QCD}");
hDummy0->GetYaxis()->SetRangeUser(0.0,0.07);
hDummy0->GetXaxis()->SetTitle("|#eta|");
hDummy0->Draw();
TGraphErrors* grWtdAvgEtaPlus = (TGraphErrors*)wtdAvgFracQCDPlus->Clone("grWtdAvgEtaPlus");
TGraphErrors* grWtdAvgEtaMinus = (TGraphErrors*)wtdAvgFracQCDMinus->Clone("grWtdAvgEtaMinus");
TGraphErrors* grEtaErr = new TGraphErrors(9);
double errSystEta = fabs( (grWtdAvgEtaPlus->GetY()[0]) - (grWtdAvgEtaMinus->GetY()[0]));
std::cout << errSystEta << std::endl;
std::cout << fractionQCDEta->GetN() << std::endl;
std::vector<double> etaBins;
etaBins.push_back(0.1);
etaBins.push_back(0.35);
etaBins.push_back(0.6);
etaBins.push_back(0.8);
etaBins.push_back(1.05);
etaBins.push_back(1.3);
etaBins.push_back(1.55);
etaBins.push_back(1.85);
etaBins.push_back(2.1);
etaBins.push_back(2.4);

for(int i = 0; i<fractionQCDEta->GetN(); ++i){
    double errStat = fractionQCDEta->GetEY()[i];
    double errTotal = TMath::Sqrt(TMath::Power(errStat,2)+TMath::Power(errSystEta,2));
    double xPt = fractionQCDEta->GetX()[i];
    double xErr = (etaBins[i+1]-etaBins[i])/2.0;
    double yPt = fractionQCDEta->GetY()[i];
    std::cout << errTotal  << std::endl;
    grEtaErr->SetPoint(i,xPt,yPt);
    grEtaErr->SetPointError(i,xErr,errTotal);
}

grEtaErr->SetFillStyle(3001);
//grEtaErr->Draw("e2same");
fractionQCDEta->Draw("pesame");
fractionQCDEta->GetYaxis()->SetTitle("f_{QCD}");
fractionQCDEta->GetXaxis()->SetTitle("|#eta|");
myText(0.5,0.5,kBlack,(char*)"#mu^{#pm}");
myText(0.5,0.4,kBlack,(char*)"0-80%");

TCanvas *c1b = new TCanvas("c1b","c1b",600,600);
hDummy0->Draw();
TGraphErrors* grPlus = (TGraphErrors*)fractionQCDEtaPlus->Clone("grPlus");
TGraphErrors* grMinus = (TGraphErrors*)fractionQCDEtaMinus->Clone("grMinus");
grPlus->SetMarkerColor(kRed);
grMinus->SetMarkerColor(kBlue);
grPlus->Draw("pesame");
grMinus->Draw("pesame");

TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
leg->SetTextFont(gStyle->GetTextFont());
leg->SetTextSize(gStyle->GetTextSize());
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(grPlus,"#mu^{+}");
leg->AddEntry(grMinus,"#mu^{-}");
leg->Draw();
myText(0.28,0.75,kBlack,(char*)"0-80%");



_fZCent->cd();
TCanvas *c2 = new TCanvas("c2","c2",600,600);
hDummy3->Draw();
hDummy3->GetYaxis()->SetTitle("b_{Z}");
hDummy3->GetYaxis()->SetRangeUser(0.0,0.17);
hDummy3->GetXaxis()->SetTitle("#LT N_{coll} #GT");
fractionZCent->Draw("pesame");
myText(0.5,0.5,kBlack,(char*)"#mu^{#pm}");
myText(0.5,0.4,kBlack,(char*)"0.1<|#eta|<2.4");

TCanvas *c2b = new TCanvas("c2b","c2b",600,600);
hDummy3->Draw();
TGraphErrors* grPlus = (TGraphErrors*)fractionZCentPlus->Clone("grPlus");
TGraphErrors* grMinus = (TGraphErrors*)fractionZCentMinus->Clone("grMinus");
grPlus->SetMarkerColor(kRed);
grMinus->SetMarkerColor(kBlue);
grPlus->Draw("pesame");
grMinus->Draw("pesame");
TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
leg->SetTextFont(gStyle->GetTextFont());
leg->SetTextSize(gStyle->GetTextSize());
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(grPlus,"#mu^{+}");
leg->AddEntry(grMinus,"#mu^{-}");
leg->Draw();
myText(0.28,0.75,kBlack,(char*)"0.1<|#eta|<2.4");



_fZEta->cd();
TCanvas *c3 = new TCanvas("c3","c3",600,600);
hDummy2->Draw();
hDummy2->GetYaxis()->SetTitle("b_{Z}");
hDummy2->GetYaxis()->SetRangeUser(0.0,0.17);
hDummy2->GetXaxis()->SetTitle("|#eta|");
fractionZEta->Draw("pesame");
myText(0.5,0.5,kBlack,(char*)"#mu^{#pm}");
myText(0.5,0.4,kBlack,(char*)"0-80%");

TCanvas *c3b = new TCanvas("c3b","c3b",600,600);
hDummy2->Draw();
TGraphErrors* grPlus = (TGraphErrors*)fractionZEtaPlus->Clone("grPlus");
TGraphErrors* grMinus = (TGraphErrors*)fractionZEtaMinus->Clone("grMinus");
grPlus->SetMarkerColor(kRed);
grMinus->SetMarkerColor(kBlue);
grPlus->Draw("pesame");
grMinus->Draw("pesame");
TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
leg->SetTextFont(gStyle->GetTextFont());
leg->SetTextSize(gStyle->GetTextSize());
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(grPlus,"#mu^{+}");
leg->AddEntry(grMinus,"#mu^{-}");
leg->Draw();
myText(0.28,0.75,kBlack,(char*)"0-80%");

////////////////
//W-->tau-->mu
///////////////
   _fTauBkg->cd();
   TCanvas *cTau = new TCanvas("cTau", "cTau",428,133,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cTau->Range(-0.7008091,-0.005093834,2.298318,0.02654155);
   cTau->SetFillColor(0);
   cTau->SetBorderMode(0);
   cTau->SetBorderSize(2);
   cTau->SetTickx(1);
   cTau->SetTicky(1);
   cTau->SetLeftMargin(0.1594828);
   cTau->SetRightMargin(0.05028735);
   cTau->SetTopMargin(0.04872881);
   cTau->SetBottomMargin(0.161017);
   cTau->SetFrameBorderMode(0);
   cTau->SetFrameBorderMode(0);

   tauBkgFractionCent0->SetMarkerSize(1.4);
   tauBkgFractionCent1->SetMarkerSize(1.4);
   tauBkgFractionCent2->SetMarkerSize(1.4);
   tauBkgFractionCent3->SetMarkerSize(1.4);
   tauBkgFractionCent4->SetMarkerSize(1.4);
   tauBkgFractionCent5->SetMarkerSize(1.4);

   tauBkgFractionCent0->SetMarkerColor(2);
   tauBkgFractionCent1->SetMarkerColor(4);
   tauBkgFractionCent2->SetMarkerColor(3);
   tauBkgFractionCent3->SetMarkerColor(6);
   tauBkgFractionCent4->SetMarkerColor(7);
   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#ff6600");
   tauBkgFractionCent5->SetMarkerColor(ci);

   tauBkgFractionCent0->GetXaxis()->SetTitle("|#eta|");
   tauBkgFractionCent0->GetYaxis()->SetTitle("N_{#mu_{#tau}}/N_{#mu_{W}}");
   tauBkgFractionCent0->GetYaxis()->SetRangeUser(0.0,0.025);

   tauBkgFractionCent0->Draw("ape");
   tauBkgFractionCent1->Draw("pesame");
   tauBkgFractionCent2->Draw("pesame");
   tauBkgFractionCent3->Draw("pesame");
   tauBkgFractionCent4->Draw("pesame");
   tauBkgFractionCent5->Draw("pesame");
 

   TLegend *legTau = new TLegend(0.658046,0.2584746,0.9281609,0.5466102,NULL,"brNDC");
   legTau->SetBorderSize(0);
   legTau->SetTextSize(0.05);
   legTau->SetLineColor(1);
   legTau->SetLineStyle(1);
   legTau->SetLineWidth(1);
   legTau->SetFillColor(0);
   legTau->SetFillStyle(1001);
   TLegendEntry *entry=legTau->AddEntry("tauBkgFractionCent0","0-5%","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.4);
   entry=legTau->AddEntry("tauBkgFractionCent1","5-10%","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(4);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.4);
   entry=legTau->AddEntry("tauBkgFractionCent2","10-15%","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(3);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.4);
   entry=legTau->AddEntry("tauBkgFractionCent3","15-20%","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(6);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.4);
   entry=legTau->AddEntry("tauBkgFractionCent4","20-40%","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(7);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.4);
   entry=legTau->AddEntry("tauBkgFractionCent5","40-80%","pe");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff6600");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.4);
   legTau->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);


}
