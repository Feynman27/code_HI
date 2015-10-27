// Trigger efficiencies from Z T&P (provided by Zvi)
// July 29, 2014

#include <iostream>
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLatex.h"

void readInTriggerEff(){


   Double_t bins[22] = {-2.7,-2.4,-2.1,-1.85,-1.55,-1.3,-1.05,-0.8,-0.6,-0.35,-0.1,0.1,0.35,0.6,0.8,1.05,1.3,1.55,1.85,2.1,2.4,2.7};
   const int nEtaBins = sizeof(bins)/sizeof(Double_t)-1;
   std::cout << "Plotting EF_mu8 trigger efficiencies in " << nEtaBins << " pseudorapidity intervals." << std::endl;

   double periodA_MCtrig[21]    = {0,0.919574,0.944285,0.954896,0.951042,0.911578,0.708582,
                                   0.698594,0.757765,0.8339,0.590597,0.817595,0.791514,0.697298,
                                   0.707661,0.894237,0.939925,0.940534,0.946149,0.936891,0};    

   double periodA_trig[21]      = {0,0.969697,0.972973,0.976471,0.971429,0.927711,0.722892,
                                   0.684211,0.763158,0.824176,0.457627,0.834951,0.87931,
                                   0.746667,0.725664,0.898734,0.963855,0.971963,0.924528,0.907407,0};

   double periodB_MCtrig[21]    = {0,0.919166,0.945205,0.954034,0.950955,0.913593,0.709901,0.700557,
                                   0.757269,0.832601,0.594685,0.817325,0.790773,0.696885,0.704813,
                                   0.893053,0.941787,0.940789,0.944381,0.936237,0};

   double periodB_trig[21]      = {0,0.935484,0.93617,0.983333,0.972973,0.850746,0.734177,0.734375,
                                   0.791209,0.785714,0.403846,0.891304,0.846154,0.692308,0.630137,
                                   0.892308,0.951613,0.972973,0.97561,0.932203,0};

   double periodA_trig_down[21] = {0,0.94287,0.948924,0.95542,0.946066,0.895756,0.672494,0.629916,
                                   0.72208,0.782194,0.394262,0.79643,0.846981,0.694829,0.682675,
                                   0.861524,0.939097,0.952561,0.882784,0.863033,0};

   double periodB_trig_down[21] = {0,0.899358,0.893948,0.959305,0.948924,0.804113,0.683045,0.677487,
                                   0.746829,0.734511,0.338207,0.85612,0.806043,0.626979,0.573171,
                                   0.850017,0.918975,0.948924,0.9411,0.89438,0};

   double periodA_trig_up[21]   = {0,0.986616,0.988092,0.989659,0.987397,0.952794,0.769774,0.73534,0.801195,
                                   0.86147,0.521922,0.869185,0.907309,0.794206,0.766024,0.929266,0.980971,
                                   0.985296,0.955615,0.941812,0};

   double periodB_trig_up[21]   = {0,0.962207,0.965953,0.995395,0.988092,0.890512,0.781432,0.786495,0.831385,
                                   0.831582,0.471872,0.920859,0.881253,0.752854,0.684779,0.926637,0.974377,
                                   0.988092,0.993166,0.960239,0};

   TGraphAsymmErrors *grEffA = new TGraphAsymmErrors(nEtaBins);
   TGraphAsymmErrors *grEffA_Mc = new TGraphAsymmErrors(nEtaBins);
   TGraphAsymmErrors *grEffA_DataOverMc = new TGraphAsymmErrors(nEtaBins);
   TGraphAsymmErrors *grEffB = new TGraphAsymmErrors(nEtaBins);
   TGraphAsymmErrors *grEffB_Mc = new TGraphAsymmErrors(nEtaBins);
   TGraphAsymmErrors *grEffB_DataOverMc = new TGraphAsymmErrors(nEtaBins);

   for(int ieta=0; ieta<nEtaBins; ++ieta){

       double xBinWover2 = fabs(bins[ieta+1]-bins[ieta])/2.0;
       double etaMed = bins[ieta]+xBinWover2;
       // Period A
       double effA = periodA_trig[ieta];
       double effA_Mc = periodA_MCtrig[ieta];
       double effA_errUpp = fabs(periodA_trig[ieta]-periodA_trig_up[ieta]);
       double effA_errLow = fabs(periodA_trig[ieta]-periodA_trig_down[ieta]);
       grEffA->SetPoint(ieta,etaMed,effA);
       grEffA->SetPointError(ieta,xBinWover2,xBinWover2,effA_errLow,effA_errUpp);
       grEffA_Mc->SetPoint(ieta,etaMed,effA_Mc);
       grEffA_Mc->SetPointError(ieta,xBinWover2,xBinWover2,0.0,0.0);
       grEffA_DataOverMc->SetPoint(ieta,etaMed,effA/effA_Mc);
       grEffA_DataOverMc->SetPointError(ieta,xBinWover2,xBinWover2,effA_errLow,effA_errUpp);
       // Period B
       double effB = periodB_trig[ieta];
       double effB_Mc = periodB_MCtrig[ieta];
       double effB_errUpp = fabs(periodB_trig[ieta]-periodB_trig_up[ieta]);
       double effB_errLow = fabs(periodB_trig[ieta]-periodB_trig_down[ieta]);
       grEffB->SetPoint(ieta,etaMed,effB);
       grEffB->SetPointError(ieta,xBinWover2,xBinWover2,effB_errLow,effB_errUpp);
       grEffB_Mc->SetPoint(ieta,etaMed,effB_Mc);
       grEffB_Mc->SetPointError(ieta,xBinWover2,xBinWover2,0.0,0.0);
       grEffB_DataOverMc->SetPoint(ieta,etaMed,effB/effB_Mc);
       grEffB_DataOverMc->SetPointError(ieta,xBinWover2,xBinWover2,effB_errLow,effB_errUpp);
   }//ieta

   TH1F* hdummy = new TH1F("hdummy","hdummy",nEtaBins,bins);
   TCanvas *cA = new TCanvas("cA", "",369,106,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cA->Range(0,0,1,1);
   cA->SetFillColor(0);
   cA->SetBorderMode(0);
   cA->SetBorderSize(2);
   cA->SetTickx(1);
   cA->SetTicky(1);
   cA->SetLeftMargin(0.16);
   cA->SetRightMargin(0.05);
   cA->SetTopMargin(0.05);
   cA->SetBottomMargin(0.16);
   cA->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: c_1
   TPad *c_1 = new TPad("c_1", "c_1",0,0.2,1,1);
   c_1->Draw();
   c_1->cd();
   c_1->Range(-3.78,0.09,3.053455,1.150556);
   c_1->SetFillColor(0);
   c_1->SetBorderMode(0);
   c_1->SetBorderSize(2);
   c_1->SetTickx(1);
   c_1->SetTicky(1);
   c_1->SetLeftMargin(0.158046);
   c_1->SetRightMargin(0.05172414);
   c_1->SetTopMargin(0.04766949);
   c_1->SetBottomMargin(0);
   c_1->SetFrameBorderMode(0);
   c_1->SetFrameBorderMode(0);


   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c_1->cd();
   hdummy->SetMinimum(0.09);
   hdummy->SetMaximum(1.1);
   hdummy->SetStats(0);
   hdummy->SetLineWidth(2);
   hdummy->SetMarkerStyle(20);
   hdummy->SetMarkerSize(1.2);
   hdummy->GetXaxis()->SetLabelFont(42);
   hdummy->GetXaxis()->SetLabelSize(0.05);
   hdummy->GetXaxis()->SetTitleSize(0.05);
   hdummy->GetXaxis()->SetTitleOffset(1.4);
   hdummy->GetXaxis()->SetTitleFont(42);
   hdummy->GetYaxis()->SetTitle("#epsilon_{trig}");
   hdummy->GetYaxis()->SetLabelFont(42);
   hdummy->GetYaxis()->SetLabelSize(0.05);
   hdummy->GetYaxis()->SetTitleSize(0.07);
   hdummy->GetYaxis()->SetTitleOffset(1.1);
   hdummy->GetYaxis()->SetTitleFont(42);
   hdummy->GetZaxis()->SetLabelFont(42);
   hdummy->GetZaxis()->SetLabelSize(0.05);
   hdummy->GetZaxis()->SetTitleSize(0.05);
   hdummy->GetZaxis()->SetTitleFont(42);
   hdummy->Draw();
   grEffA->Draw("pesame");
   grEffA_Mc->SetMarkerColor(kBlue);
   grEffA_Mc->SetMarkerStyle(kOpenCircle);
   grEffA_Mc->Draw("pesame");

   TLatex *   tex = new TLatex(0.7567114,0.8720855,"Period A");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();

   TLegend *leg = new TLegend(0.1912752,0.3069948,0.4614094,0.4575777,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);

   TLegendEntry *entry=leg->AddEntry("grEffA","Data","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("grEffA_Mc","MC","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(4);
   entry->SetMarkerSize(1.2);
   leg->Draw();

   c_1->Modified();
   cA->cd();


// ------------>Primitives in pad: c_2
//

   TPad* c_2 = new TPad("c_2", "c_2",0,0,1,0.3241525);
   c_2->Draw();
   c_2->cd();
   c_2->Range(-3.81024,0.5628571,3.077882,1.18);
   c_2->SetFillColor(0);
   c_2->SetBorderMode(0);
   c_2->SetBorderSize(2);
   c_2->SetTickx(1);
   c_2->SetTicky(1);
   c_2->SetLeftMargin(0.158046);
   c_2->SetRightMargin(0.05172414);
   c_2->SetTopMargin(0);
   c_2->SetBottomMargin(0.2222223);
   c_2->SetFrameBorderMode(0);
   c_2->SetFrameBorderMode(0);
   

   TH1F *hdummy2 = new TH1F("hdummy2","Graph",100,-3.24,3.24);
   hdummy2->SetMinimum(0.7);
   hdummy2->SetMaximum(1.18);
   hdummy2->SetDirectory(0);
   hdummy2->SetStats(0);
   hdummy2->SetLineWidth(2);
   hdummy2->SetMarkerStyle(20);
   hdummy2->SetMarkerSize(1.2);
   hdummy2->GetXaxis()->SetTitle("#eta^{#mu}");
   hdummy2->GetXaxis()->SetRange(9,92);
   hdummy2->GetXaxis()->SetLabelFont(42);
   hdummy2->GetXaxis()->SetLabelSize(0.12);
   hdummy2->GetXaxis()->SetTitleSize(0.12);
   hdummy2->GetXaxis()->SetTitleOffset(0.75);
   hdummy2->GetXaxis()->SetTitleFont(42);
   hdummy2->GetYaxis()->SetTitle("Data/MC");
   hdummy2->GetYaxis()->SetNdivisions(409);
   hdummy2->GetYaxis()->SetLabelFont(42);
   hdummy2->GetYaxis()->SetLabelOffset(0.004);
   hdummy2->GetYaxis()->SetLabelSize(0.12);
   hdummy2->GetYaxis()->SetTitleSize(0.11);
   hdummy2->GetYaxis()->SetTickLength(0.02);
   hdummy2->GetYaxis()->SetTitleOffset(0.68);
   hdummy2->GetYaxis()->SetTitleFont(42);
   hdummy2->GetZaxis()->SetLabelFont(42);
   hdummy2->GetZaxis()->SetLabelSize(0.05);
   hdummy2->GetZaxis()->SetTitleSize(0.05);
   hdummy2->GetZaxis()->SetTitleFont(42);
   grEffA_DataOverMc->SetHistogram(hdummy2);

   grEffA_DataOverMc->Draw("ape");
   c_2->Modified();
   cA->cd();
   cA->Modified();
   cA->cd();
   cA->SetSelected(cA);

   TCanvas *cB = new TCanvas("cB", "",369,106,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cB->Range(0,0,1,1);
   cB->SetFillColor(0);
   cB->SetBorderMode(0);
   cB->SetBorderSize(2);
   cB->SetTickx(1);
   cB->SetTicky(1);
   cB->SetLeftMargin(0.16);
   cB->SetRightMargin(0.05);
   cB->SetTopMargin(0.05);
   cB->SetBottomMargin(0.16);
   cB->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: c_1
   TPad *c_1b = new TPad("c_1b", "c_1b",0,0.2,1,1);
   c_1b->Draw();
   c_1b->cd();
   c_1b->Range(-3.78,0.09,3.053455,1.150556);
   c_1b->SetFillColor(0);
   c_1b->SetBorderMode(0);
   c_1b->SetBorderSize(2);
   c_1b->SetTickx(1);
   c_1b->SetTicky(1);
   c_1b->SetLeftMargin(0.158046);
   c_1b->SetRightMargin(0.05172414);
   c_1b->SetTopMargin(0.04766949);
   c_1b->SetBottomMargin(0);
   c_1b->SetFrameBorderMode(0);
   c_1b->SetFrameBorderMode(0);


   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c_1b->cd();
   hdummy->SetMinimum(0.09);
   hdummy->SetMaximum(1.1);
   hdummy->SetStats(0);
   hdummy->SetLineWidth(2);
   hdummy->SetMarkerStyle(20);
   hdummy->SetMarkerSize(1.2);
   hdummy->GetXaxis()->SetLabelFont(42);
   hdummy->GetXaxis()->SetLabelSize(0.05);
   hdummy->GetXaxis()->SetTitleSize(0.05);
   hdummy->GetXaxis()->SetTitleOffset(1.4);
   hdummy->GetXaxis()->SetTitleFont(42);
   hdummy->GetYaxis()->SetTitle("#epsilon_{trig}");
   hdummy->GetYaxis()->SetLabelFont(42);
   hdummy->GetYaxis()->SetLabelSize(0.05);
   hdummy->GetYaxis()->SetTitleSize(0.07);
   hdummy->GetYaxis()->SetTitleOffset(1.1);
   hdummy->GetYaxis()->SetTitleFont(42);
   hdummy->GetZaxis()->SetLabelFont(42);
   hdummy->GetZaxis()->SetLabelSize(0.05);
   hdummy->GetZaxis()->SetTitleSize(0.05);
   hdummy->GetZaxis()->SetTitleFont(42);
   hdummy->Draw();
   grEffB->Draw("pesame");
   grEffB_Mc->SetMarkerColor(kBlue);
   grEffB_Mc->SetMarkerStyle(kOpenCircle);
   grEffB_Mc->Draw("pesame");

   TLatex *   texB = new TLatex(0.7567114,0.8720855,"Period B");
   texB->SetNDC();
   texB->SetTextFont(42);
   texB->SetLineWidth(2);
   texB->Draw();

   TLegend *legB = new TLegend(0.1912752,0.3069948,0.4614094,0.4575777,NULL,"brNDC");
   legB->SetBorderSize(0);
   legB->SetTextSize(0.05);
   legB->SetLineColor(1);
   legB->SetLineStyle(1);
   legB->SetLineWidth(1);
   legB->SetFillColor(0);
   legB->SetFillStyle(1001);

   TLegendEntry *entryB=legB->AddEntry("grEffB","Data","p");
   entryB->SetLineColor(1);
   entryB->SetLineStyle(1);
   entryB->SetLineWidth(1);
   entryB->SetMarkerColor(1);
   entryB->SetMarkerStyle(20);
   entryB->SetMarkerSize(1.2);
   entryB=legB->AddEntry("grEffB_Mc","MC","p");
   entryB->SetLineColor(1);
   entryB->SetLineStyle(1);
   entryB->SetLineWidth(1);
   entryB->SetMarkerColor(1);
   entryB->SetMarkerStyle(4);
   entryB->SetMarkerSize(1.2);
   legB->Draw();

   c_1b->Modified();
   cB->cd();


// ------------>Primitives in pad: c_2b
//

   TPad* c_2b = new TPad("c_2b", "c_2b",0,0,1,0.3241525);
   c_2b->Draw();
   c_2b->cd();
   c_2b->Range(-3.81024,0.5628571,3.077882,1.18);
   c_2b->SetFillColor(0);
   c_2b->SetBorderMode(0);
   c_2b->SetBorderSize(2);
   c_2b->SetTickx(1);
   c_2b->SetTicky(1);
   c_2b->SetLeftMargin(0.158046);
   c_2b->SetRightMargin(0.05172414);
   c_2b->SetTopMargin(0);
   c_2b->SetBottomMargin(0.2222223);
   c_2b->SetFrameBorderMode(0);
   c_2b->SetFrameBorderMode(0);
   
   grEffB_DataOverMc->SetHistogram(hdummy2);

   grEffB_DataOverMc->Draw("ape");
   c_2b->Modified();
   cB->cd();
   cB->Modified();
   cB->cd();
   cB->SetSelected(cB);

}
