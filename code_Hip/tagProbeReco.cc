#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TH1.h"
#include "TH1F.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLegendEntry.h"
#include "TLorentzVector.h"

#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#ifdef __MAKECINT__
#pragma link C++ class vector<vector<bool> >+;
#endif
#ifdef __MAKECINT__
#pragma link C++ class vector<vector<int> >+;
#endif
#ifdef __MAKECINT__
#pragma link C++ class vector<vector<float> >+;
#endif
#ifdef __MAKECINT__
#pragma link C++ class vector<vector<double> >+;
#endif

void plotOnCanvas(TH1F* hpassed, TH1F* htotal,TH1F* hpassed_ss, TH1F* htotal_ss, TGraphAsymmErrors* grEff,
		         TString cLabel,TString yLabelUpp,TString sCent, TString legLabel1,TString legLabel2,TString xLabel,TString yEffLabel,const char* outFileName,
			     bool doEta = false, bool doPt=false, bool doCentrality=false){

   TCanvas *c = new TCanvas(cLabel, cLabel ,2002,92,600,800);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c->Range(0,0,1,1);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetLeftMargin(0.16);
   c->SetRightMargin(0.05);
   c->SetTopMargin(0.05);
   c->SetBottomMargin(0.16);
   c->SetFrameBorderMode(0);
   TPad *c_1 = new TPad("c_1", "c_1",0,0.2,1,1);
   c_1->Draw();
   c_1->cd();
   c_1->Range(-3.508493,-289.2482,2.818471,803.059);
   c_1->SetFillColor(0);
   c_1->SetBorderMode(0);
   c_1->SetBorderSize(2);
   c_1->SetTickx(1);
   c_1->SetTicky(1);
   c_1->SetLeftMargin(0.159396);
   c_1->SetRightMargin(0.05033557);
   c_1->SetTopMargin(0.04857513);
   c_1->SetBottomMargin(0.2648964);
   c_1->SetFrameBorderMode(0);
   c_1->SetFrameBorderMode(0);
   c_1->cd();

   htotal->SetMarkerStyle(20);
   htotal->GetXaxis()->SetLabelFont(42);
   htotal->GetXaxis()->SetLabelSize(0.05);
   htotal->GetXaxis()->SetTitleSize(0.05);
   htotal->GetXaxis()->SetTitleOffset(1.4);
   htotal->GetXaxis()->SetTitleFont(42);
   htotal->GetYaxis()->SetTitle(yLabelUpp);
   htotal->GetYaxis()->SetLabelFont(42);
   htotal->GetYaxis()->SetTitleSize(0.05);
   htotal->GetYaxis()->SetTitleOffset(1.4);
   htotal->GetYaxis()->SetTitleFont(42);
   if(doEta) htotal->GetXaxis()->SetRangeUser(-2.47,2.47);
   else if(doPt) htotal->GetXaxis()->SetRangeUser(20.,100.);
   else if(doCentrality) htotal->GetXaxis()->SetRangeUser(0.0,0.91);
   else {
        std::cout << "ERROR: Must choose x-axis. " << std::endl;
        exit(0);
   }

   htotal->SetMinimum(0.01);
   htotal->SetMaximum(htotal->GetBinContent(htotal->GetMaximumBin())+10.0);
   htotal->Draw("pe");

   htotal_ss->SetMarkerStyle(21);
   htotal_ss->GetXaxis()->SetLabelFont(42);
   htotal_ss->GetXaxis()->SetLabelSize(0.05);
   htotal_ss->GetXaxis()->SetTitleSize(0.05);
   htotal_ss->GetXaxis()->SetTitleOffset(1.4);
   htotal_ss->GetXaxis()->SetTitleFont(42);
   htotal_ss->GetYaxis()->SetTitle(yLabelUpp);
   htotal_ss->GetYaxis()->SetLabelFont(42);
   htotal_ss->GetYaxis()->SetTitleSize(0.05);
   htotal_ss->GetYaxis()->SetTitleOffset(1.4);
   htotal_ss->GetYaxis()->SetTitleFont(42);
   htotal_ss->Draw("pesame");

   hpassed->SetMarkerStyle(4);
   hpassed->SetMarkerSize(1.2);
   hpassed->GetXaxis()->SetLabelFont(42);
   hpassed->GetXaxis()->SetLabelSize(0.05);
   hpassed->GetXaxis()->SetTitleSize(0.05);
   hpassed->GetXaxis()->SetTitleOffset(1.4);
   hpassed->GetXaxis()->SetTitleFont(42);
   hpassed->GetYaxis()->SetLabelFont(42);
   hpassed->GetYaxis()->SetLabelSize(0.05);
   hpassed->GetYaxis()->SetTitleSize(0.05);
   hpassed->GetYaxis()->SetTitleOffset(1.4);
   hpassed->GetYaxis()->SetTitleFont(42);
   hpassed->Draw("pesame");

   hpassed_ss->SetMarkerStyle(25);
   hpassed_ss->SetMarkerSize(1.2);
   hpassed_ss->GetXaxis()->SetLabelFont(42);
   hpassed_ss->GetXaxis()->SetLabelSize(0.05);
   hpassed_ss->GetXaxis()->SetTitleSize(0.05);
   hpassed_ss->GetXaxis()->SetTitleOffset(1.4);
   hpassed_ss->GetXaxis()->SetTitleFont(42);
   hpassed_ss->GetYaxis()->SetLabelFont(42);
   hpassed_ss->GetYaxis()->SetLabelSize(0.05);
   hpassed_ss->GetYaxis()->SetTitleSize(0.05);
   hpassed_ss->GetYaxis()->SetTitleOffset(1.4);
   hpassed_ss->GetYaxis()->SetTitleFont(42);
   hpassed_ss->Draw("pesame");

   TLatex *   tex = new TLatex(0.7567114,0.8720855,sCent);
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
   //leg->SetTextFont(4);
   TLegendEntry *entry=leg->AddEntry("htotal",legLabel1,"p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry=leg->AddEntry("hpassed",legLabel2,"p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(4);
   entry->SetMarkerSize(1.2);
   leg->Draw();

   c_1->Modified();
   c->cd();
// ------------>Primitives in pad: c_2
   TPad* c_2 = new TPad("c_2", "c_2",0,-0,1,0.4119171);
   c_2->Draw();
   c_2->cd();
   c_2->Range(-3.536561,0.2251004,2.841019,1.055221);
   c_2->SetFillColor(0);
   c_2->SetBorderMode(0);
   c_2->SetBorderSize(2);
   c_2->SetTickx(1);
   c_2->SetTicky(1);
   c_2->SetLeftMargin(0.159396);
   c_2->SetRightMargin(0.05033557);
   c_2->SetTopMargin(0.006289308);
   c_2->SetBottomMargin(0.2106918);
   c_2->SetFrameBorderMode(0);
   c_2->SetFrameBorderMode(0);
 
   grEff->GetXaxis()->SetTitle(xLabel);
   grEff->GetXaxis()->SetLabelFont(42);
   grEff->GetXaxis()->SetLabelSize(0.08);
   grEff->GetXaxis()->SetTitleSize(0.08);
   grEff->GetXaxis()->SetTitleOffset(1.07);
   grEff->GetXaxis()->SetTitleFont(42);
   grEff->GetYaxis()->SetTitle(yEffLabel);
   grEff->GetYaxis()->SetLabelFont(42);
   grEff->GetYaxis()->SetLabelSize(0.08);
   grEff->GetYaxis()->SetTitleSize(0.08);
   grEff->GetYaxis()->SetTitleOffset(0.81);
   grEff->GetYaxis()->SetTitleFont(42);
   grEff->GetXaxis()->SetRangeUser(20.,100.);
   if(doEta) grEff->GetXaxis()->SetRangeUser(-2.47,2.47);
   else if(doPt) grEff->GetXaxis()->SetRangeUser(20.,100.);
   else if(doCentrality) grEff->GetXaxis()->SetRangeUser(0.,0.91);
   grEff->SetMinimum(0.7);
   grEff->SetMaximum(1.01);
   grEff->SetLineWidth(2);
   grEff->Draw("ape");
   c_2->Modified();
   c->cd();
   c->Modified();
   c->cd();
   c->SetSelected(c);
   c->Write();
    
} 

void calculateEfficiency(TGraphAsymmErrors* const gr, const TH1* hPassed_OS, const TH1* hPassed_SS, const TH1* hTotal_OS, const TH1* hTotal_SS, int nBins){
    // subtract combinatoric background
    //TH1F* hPassed=NULL,*hTotal=NULL;
    TH1F* hPassed = new TH1F("hPassed","hPassed",hPassed_OS->GetNbinsX(),hPassed_OS->GetXaxis()->GetXmin(),hPassed_OS->GetXaxis()->GetXmax());
    TH1F* hTotal = new TH1F("hTotal","hTotal",hPassed_OS->GetNbinsX(),hPassed_OS->GetXaxis()->GetXmin(),hPassed_OS->GetXaxis()->GetXmax());
    std::cout << "OS probes: " << hTotal_OS->Integral() << std::endl;
    std::cout << "SS probes: " << hTotal_SS->Integral() << std::endl;
    std::cout << "OS probes matched: " << hPassed_OS->Integral() << std::endl;
    std::cout << "SS probes matched: " << hPassed_SS->Integral() << std::endl;
    hPassed->Add(hPassed_OS,hPassed_SS,1.0,-1.0);
    hTotal->Add(hTotal_OS,hTotal_SS,1.0,-1.0);
    std::cout << "Probes after bkg subtraction: " << hTotal->Integral() << std::endl;
    std::cout << "Matched probes after bkg subtraction: " << hPassed->Integral() << std::endl;
	//gr->BayesDivide(hPassed,hTotal);
	gr->BayesDivide(hPassed_OS,hTotal_OS);
}

int getCentralityBin(float centrality, bool doCentrality = false){

	if(!doCentrality) return 0;
	int index = 0;
	if(centrality>=0.0&&centrality<0.05) return index; else ++index;
	if(centrality>=0.05&&centrality<0.10) return index; else ++index;
	if(centrality>=0.10&&centrality<0.20) return index; else ++index;
	if(centrality>=0.20&&centrality<0.40) return index; else ++index;
	if(centrality>=0.40&&centrality<0.90) return index; else ++index;
	return -1;
}

void tagProbeReco(){

	//bool doCentrality = true;
	TFile *outFile = new TFile("tagAndProbeRecEff.root","recreate");
	// 09.22 file uses mu_ms* container; 09.18 file uses mu_me* container
//	TString fileName = "/usatlas/u/tbales/scratch/Data2013Hip/TagAndProbe/HIpSingleMuonHardProbesDataPeriodABMerge_recoTagAndProbe.09.22.2013.root";

//	If using 09.18 file, replace tag_n with tag_isTagMuon->size() since
//	I forgot to add tag_n variable
	//TString fileName = "/usatlas/u/tbales/scratch/Data2013Hip/TagAndProbe/HIpSingleMuonHardProbesDataPeriodABMerge_recoTagAndProbe.09.18.2013.root";
	TString fileName = "$HOME/scratch/Data2013Hip/ntupleW.07.22.2014.root";
	TChain* tree = new TChain("tree","tree");
	tree->Add(fileName);
	std::cout << "Number of entries in data tree: " << tree->GetEntries() << std::endl;
	int tag_n;
    std::vector<bool>* tag_isTagMuon = 0;
    std::vector<float>* tag_pt = 0;
    std::vector<float>* tag_eta = 0;
    std::vector<float>* tag_phi = 0;
    std::vector<float>* tag_charge = 0;
	std::vector<int>* idprobe_n = 0;
	std::vector<std::vector<bool> >* idprobe_isIDProbeMatched = 0;
	std::vector<std::vector<float> >* idprobe_eta = 0;
	std::vector<std::vector<float> >* idprobe_phi = 0;
	std::vector<std::vector<float> >* idprobe_pt = 0;
	std::vector<std::vector<float> >* idprobe_charge = 0;
	std::vector<int>* msprobe_n = 0;
	std::vector<std::vector<bool> >* msprobe_isMSProbeMatched = 0;
	std::vector<std::vector<float> >* msprobe_eta = 0;
	std::vector<std::vector<float> >* msprobe_phi = 0;
	std::vector<std::vector<float> >* msprobe_pt = 0;
	std::vector<std::vector<float> >* msprobe_charge = 0;
	
	float centrality;

	tree->SetBranchAddress("centrality",&centrality);
	tree->SetBranchAddress("tag_isTagMuon",&tag_isTagMuon);
	tree->SetBranchAddress("tag_n",&tag_n);
	tree->SetBranchAddress("tag_pt",&tag_pt);
	tree->SetBranchAddress("tag_eta",&tag_eta);
	tree->SetBranchAddress("tag_phi",&tag_phi);
	tree->SetBranchAddress("tag_charge",&tag_charge);
	tree->SetBranchAddress("idprobe_n",&idprobe_n);
	tree->SetBranchAddress("idprobe_isIDProbeMatched",&idprobe_isIDProbeMatched);
	tree->SetBranchAddress("idprobe_eta",&idprobe_eta);
	tree->SetBranchAddress("idprobe_phi",&idprobe_phi);
	tree->SetBranchAddress("idprobe_pt",&idprobe_pt);
	tree->SetBranchAddress("idprobe_charge",&idprobe_charge);
	tree->SetBranchAddress("msprobe_n",&msprobe_n);
	tree->SetBranchAddress("msprobe_isMSProbeMatched",&msprobe_isMSProbeMatched);
	tree->SetBranchAddress("msprobe_eta",&msprobe_eta);
	tree->SetBranchAddress("msprobe_phi",&msprobe_phi);
	tree->SetBranchAddress("msprobe_pt",&msprobe_pt);
	tree->SetBranchAddress("msprobe_charge",&msprobe_charge);

	// --- Set branch status ---
  	tree->SetBranchStatus("*",0) ;
	tree->SetBranchStatus("centrality",1);
	tree->SetBranchStatus("tag_isTagMuon",1);
	tree->SetBranchStatus("tag_n",1);
	tree->SetBranchStatus("tag_pt",1);
	tree->SetBranchStatus("tag_eta",1);
	tree->SetBranchStatus("tag_phi",1);
	tree->SetBranchStatus("tag_charge",1);
	tree->SetBranchStatus("idprobe_n",1);
	tree->SetBranchStatus("idprobe_isIDProbeMatched",1);
	tree->SetBranchStatus("idprobe_eta",1);
	tree->SetBranchStatus("idprobe_phi",1);
	tree->SetBranchStatus("idprobe_pt",1);
	tree->SetBranchStatus("idprobe_charge",1);
	tree->SetBranchStatus("msprobe_n",1);
	tree->SetBranchStatus("msprobe_isMSProbeMatched",1);
	tree->SetBranchStatus("msprobe_eta",1);
	tree->SetBranchStatus("msprobe_phi",1);
	tree->SetBranchStatus("msprobe_pt",1);
	tree->SetBranchStatus("msprobe_charge",1);
		
    /*std::vector<double> centralityBins;
	centralityBins.push_back(0.0);
	if(doCentrality){
	  centralityBins.push_back(0.05);
	  centralityBins.push_back(0.10);
	  centralityBins.push_back(0.20);
	  centralityBins.push_back(0.40);
	}
	centralityBins.push_back(0.9);
    */
    double xCentralityBins[] = {0.0,0.05,0.10,0.20,0.40,0.9}; 
	const int nCentralityBins = sizeof(xCentralityBins)/sizeof(double) - 1;
	TString centLabels[] = {"0-5%","5-10%","10-20%","20-40%","40-90%"};

	//double xEtaBins[] = {-2.5,-2.1,-1.55,-1.05,-0.6,-0.1,0.1,0.6,1.05,1.55,2.1,2.5};
    // Zvi's binning
    double xEtaBins[] = {-2.47,-2.1,-1.85,-1.55,-1.3,-1.05,-0.8,-0.6,-0.35,-0.1,0.1,0.35,0.6,0.8,1.05,1.3,1.55,1.85,2.1,2.47};
	int nEtaBins = sizeof(xEtaBins)/sizeof(double) - 1;

    // Zvi's binning
    double xPtBins[] = {0,10,15,20,25,30,33,36,38,40,42,45,50,60,70,80,100,150};


/*	double xPtBins[19];
	double ptLowBinEdge = 0.0;
	double ptBinWidth = 3.0; //1GeV binning
	double ptMax = 39.0;
	int ptBinIndex;
	// 4 GeV binning from 0.0,100.0
	for(int ipt = 0; ; ++ipt){

           ptBinIndex = ipt;
	   if(ptLowBinEdge>=ptMax) break;
	   xPtBins[ipt] = ptLowBinEdge;
	   ptLowBinEdge+=ptBinWidth;
	   std::cout << "Pt element " << ipt << ": " << xPtBins[ipt] << std::endl;
		
	}//ipt
	ptBinWidth = 8.0;
	ptLowBinEdge = ptMax;
	ptMax = 71.0;
	
 	for(int ipt = ptBinIndex; ; ++ipt){

           ptBinIndex = ipt;
	   if(ptLowBinEdge>=ptMax) break;
	   xPtBins[ipt] = ptLowBinEdge;
	   ptLowBinEdge+=ptBinWidth;
	   std::cout << "Pt element " << ipt << ": " << xPtBins[ipt] << std::endl;
		
	}//ipt
         
	ptBinWidth = 30.0;
	ptLowBinEdge = ptMax;
	ptMax = 101.0;
	for(int ipt = ptBinIndex; ; ++ipt){

	   if(ptLowBinEdge>ptMax) break;
	   xPtBins[ipt] = ptLowBinEdge;
	   ptLowBinEdge+=ptBinWidth;
	   std::cout << "Pt element " << ipt << ": " << xPtBins[ipt] << std::endl;
		
	}//ipt
*/
	int nPtBins = sizeof(xPtBins)/sizeof(double) - 1;
	std::cout << "Number of pt bins: " << nPtBins << std::endl;
	std::cout << "Number of eta bins: " << nEtaBins << std::endl;
	std::cout << "Number of centrality bins: " << nCentralityBins << std::endl;

	TH1F* hEtaIDProbes,*hEtaMSProbes,*hPtIDProbes,*hPtMSProbes,*hCentIDProbes,*hCentMSProbes;
	TH1F* hEtaIDProbes_SS,*hEtaMSProbes_SS,*hPtIDProbes_SS,*hPtMSProbes_SS,*hCentIDProbes_SS,*hCentMSProbes_SS;
	TH1F* hEtaIDProbesMatched,*hEtaMSProbesMatched,*hPtIDProbesMatched,*hPtMSProbesMatched,*hCentIDProbesMatched,*hCentMSProbesMatched;
	TH1F* hEtaIDProbesMatched_SS,*hEtaMSProbesMatched_SS,*hPtIDProbesMatched_SS,*hPtMSProbesMatched_SS,*hCentIDProbesMatched_SS,*hCentMSProbesMatched_SS;
    TH1F* hTagIdProbeMass, *hTagIdProbeMass_SS, *hTagMsProbeMass, *hTagMsProbeMass_SS;
    TH1F* hTagIdProbeMassMatched, *hTagIdProbeMassMatched_SS, *hTagMsProbeMassMatched, *hTagMsProbeMassMatched_SS;

    hEtaIDProbes 		= new TH1F("hEtaIDProbes","hEtaIDProbes",nEtaBins,xEtaBins);
    hEtaMSProbes 		= new TH1F("hEtaMSProbes","hEtaMSProbes",nEtaBins,xEtaBins);
    hPtIDProbes	 	= new TH1F("hPtIDProbes","hPtIDProbes",nPtBins, xPtBins);
    hPtMSProbes 		= new TH1F("hPtMSProbes","hPtMSProbes",nPtBins,xPtBins);
    hCentIDProbes	 	= new TH1F("hCentIDProbes","hCentIDProbes",nCentralityBins, xCentralityBins);
    hCentMSProbes 		= new TH1F("hCentMSProbes","hCentMSProbes",nCentralityBins,xCentralityBins);
    hEtaIDProbesMatched 	= new TH1F("hEtaIDProbesMatched","hEtaIDProbesMatched",nEtaBins,xEtaBins);
    hEtaMSProbesMatched 	= new TH1F("hEtaMSProbesMatched","hEtaMSProbesMatched",nEtaBins,xEtaBins);
    hPtIDProbesMatched 	= new TH1F("hPtIDProbesMatched","hPtIDProbesMatched",nPtBins, xPtBins);
    hPtMSProbesMatched 	= new TH1F("hPtMSProbesMatched","hPtMSProbesMatched",nPtBins, xPtBins);
    hCentIDProbesMatched 	= new TH1F("hCentIDProbesMatched","hCentIDProbesMatched",nCentralityBins, xCentralityBins);
    hCentMSProbesMatched 	= new TH1F("hCentMSProbesMatched","hCentMSProbesMatched",nCentralityBins, xCentralityBins);

    hEtaIDProbes_SS 		= new TH1F("hEtaIDProbes_SS","hEtaIDProbes_SS",nEtaBins,xEtaBins);
    hEtaMSProbes_SS 		= new TH1F("hEtaMSProbes_SS","hEtaMSProbes_SS",nEtaBins,xEtaBins);
    hPtIDProbes_SS	 	= new TH1F("hPtIDProbes_SS","hPtIDProbes_SS",nPtBins, xPtBins);
    hPtMSProbes_SS 		= new TH1F("hPtMSProbes_SS","hPtMSProbes_SS",nPtBins,xPtBins);
    hCentIDProbes_SS	 	= new TH1F("hCentIDProbes_SS","hCentIDProbes_SS",nCentralityBins, xCentralityBins);
    hCentMSProbes_SS 		= new TH1F("hCentMSProbes_SS","hCentMSProbes_SS",nCentralityBins,xCentralityBins);
    hEtaIDProbesMatched_SS 	= new TH1F("hEtaIDProbesMatched_SS","hEtaIDProbesMatched_SS",nEtaBins,xEtaBins);
    hEtaMSProbesMatched_SS 	= new TH1F("hEtaMSProbesMatched_SS","hEtaMSProbesMatched_SS",nEtaBins,xEtaBins);
    hPtIDProbesMatched_SS 	= new TH1F("hPtIDProbesMatched_SS","hPtIDProbesMatched_SS",nPtBins, xPtBins);
    hPtMSProbesMatched_SS 	= new TH1F("hPtMSProbesMatched_SS","hPtMSProbesMatched_SS",nPtBins, xPtBins);
    hCentIDProbesMatched_SS 	= new TH1F("hCentIDProbesMatched_SS","hCentIDProbesMatched_SS",nCentralityBins, xCentralityBins);
    hCentMSProbesMatched_SS 	= new TH1F("hCentMSProbesMatched_SS","hCentMSProbesMatched_SS",nCentralityBins, xCentralityBins);

    hTagMsProbeMass = new TH1F("hTagMsProbeMass","hTagMsProbeMass",60,60.0,120.0);
    hTagMsProbeMass->Sumw2();
    hTagIdProbeMass = new TH1F("hTagIdProbeMass","hTagIdProbeMass",60,60.0,120.0);
    hTagIdProbeMass->Sumw2();
    hTagMsProbeMass_SS = new TH1F("hTagMsProbeMass_SS","hTagMsProbeMass_SS",60,60.0,120.0);
    hTagMsProbeMass_SS->Sumw2();
    hTagIdProbeMass_SS = new TH1F("hTagIdProbeMass_SS","hTagIdProbeMass_SS",60,60.0,120.0);
    hTagIdProbeMass_SS->Sumw2();

    hTagMsProbeMassMatched = new TH1F("hTagMsProbeMassMatched","hTagMsProbeMassMatched",60,60.0,120.0);
    hTagMsProbeMassMatched->Sumw2();
    hTagIdProbeMassMatched = new TH1F("hTagIdProbeMassMatched","hTagIdProbeMassMatched",60,60.0,120.0);
    hTagIdProbeMassMatched->Sumw2();
    hTagMsProbeMassMatched_SS = new TH1F("hTagMsProbeMassMatched_SS","hTagMsProbeMassMatched_SS",60,60.0,120.0);
    hTagMsProbeMassMatched_SS->Sumw2();
    hTagIdProbeMassMatched_SS = new TH1F("hTagIdProbeMassMatched_SS","hTagIdProbeMassMatched_SS",60,60.0,120.0);
    hTagIdProbeMassMatched_SS->Sumw2();


    TGraphAsymmErrors* _MSEffEta = new TGraphAsymmErrors(nEtaBins);
    TGraphAsymmErrors* _MSEffPt  = new TGraphAsymmErrors(nPtBins) ;
    TGraphAsymmErrors* _MSEffCent  = new TGraphAsymmErrors(nCentralityBins) ;
    TGraphAsymmErrors* _IDEffEta  = new TGraphAsymmErrors(nEtaBins) ;
    TGraphAsymmErrors* _IDEffPt  = new TGraphAsymmErrors(nPtBins) ;
    TGraphAsymmErrors* _IDEffCent  = new TGraphAsymmErrors(nCentralityBins) ;
    TGraphAsymmErrors* _TotEffEta  = new TGraphAsymmErrors(nEtaBins) ;
    TGraphAsymmErrors* _TotEffPt  = new TGraphAsymmErrors(nPtBins) ;
    TGraphAsymmErrors* _TotEffCent  = new TGraphAsymmErrors(nCentralityBins) ;

	//char cEtaIDProbes[30],cEtaMSProbes[30],cPtIDProbes[30],cPtMSProbes[30],cEtaIDProbesMatched[30],cEtaMSProbesMatched[30],cPtIDProbesMatched[30],cPtMSProbesMatched[30];
	for(int iev = 0; iev < tree->GetEntries(); ++iev){
      if(iev%10000==0) std::cout << "Event: " << iev << std::endl;
//	  if(iev==10000) break; //hack
      tree->LoadTree(iev);
      tree->GetEntry(iev);
      
      if(tag_isTagMuon->size()!=tag_n){
        std::cout << "ERROR: Number of tags != tag container size." << std::endl;
        exit(0);
      }

      if(centrality>0.9) continue;
      //if(tag_n>0) std::cout << "No. tags: " << tag_n << std::endl;
	  for(int itag = 0; itag < tag_n; ++itag){

        // safety
        if(tag_isTagMuon->at(itag)==false) continue;

        // construct 4-vector of the tag
        TLorentzVector vTag; vTag.SetPtEtaPhiM(tag_pt->at(itag),tag_eta->at(itag),tag_phi->at(itag),0.10566);
		// loop over probes belonging to each tag
		// ID probes
        //std::cout << "probes " << idprobe_n->at(itag) << std::endl;
		for(int iprobe = 0; iprobe < idprobe_n->at(itag); ++iprobe){
           // construct 4-vector of the id probe
           /*std::cout << "probe pt " << idprobe_pt->at(itag)[iprobe]<< std::endl;
           std::cout << "probe eta " << idprobe_eta->at(itag)[iprobe]<< std::endl;
           std::cout << "probe phi " << idprobe_phi->at(itag)[iprobe]<< std::endl;
           exit(0);
           */

           TLorentzVector vIdProbe; vIdProbe.SetPtEtaPhiM(idprobe_pt->at(itag)[iprobe],idprobe_eta->at(itag)[iprobe],idprobe_phi->at(itag)[iprobe],0.10566);
           // add TP vectors
           TLorentzVector vSum = vIdProbe+vTag;

           // invariant TP mass
           double tpMass = vSum.M();

           if(tag_charge->at(itag)*idprobe_charge->at(itag)[iprobe]>0.0) {

               // fill invariant mass distribution for same-signed pairs
               hTagIdProbeMass_SS->Fill(tpMass);
               if(idprobe_isIDProbeMatched->at(itag)[iprobe]==true) hTagIdProbeMassMatched_SS->Fill(tpMass);
               if(fabs(tpMass-91.1876)>10.0) continue;

               hPtIDProbes_SS->Fill(idprobe_pt->at(itag)[iprobe]);
               hEtaIDProbes_SS->Fill(idprobe_eta->at(itag)[iprobe]);
               hCentIDProbes_SS->Fill(centrality);

                /*std::cout << "tag pt " << vTag.Pt() << std::endl;
                std::cout << "tag eta " << vTag.Eta() << std::endl;
                std::cout << "tag phi " << vTag.Phi() << std::endl;
                */
		       // same-sign probes matched to reconstructed muon
               if(idprobe_isIDProbeMatched->at(itag)[iprobe]==true){
                hPtIDProbesMatched_SS->Fill(idprobe_pt->at(itag)[iprobe]);
                hEtaIDProbesMatched_SS->Fill(idprobe_eta->at(itag)[iprobe]);
                hCentIDProbesMatched_SS->Fill(centrality);
               }

               continue;
           }


           // fill invariant mass distribution for oppositely-signed pairs
           hTagIdProbeMass->Fill(tpMass);
           if(idprobe_isIDProbeMatched->at(itag)[iprobe]==true) hTagIdProbeMassMatched->Fill(tpMass);
           if(fabs(tpMass-91.1876)>10.0) continue;

           //opposite-signed TP pairs
           hPtIDProbes->Fill(idprobe_pt->at(itag)[iprobe]);
           hEtaIDProbes->Fill(idprobe_eta->at(itag)[iprobe]);
           hCentIDProbes->Fill(centrality);
		   // probes matched to reconstructed muon
		   if(idprobe_isIDProbeMatched->at(itag)[iprobe]==true){
		   	hPtIDProbesMatched->Fill(idprobe_pt->at(itag)[iprobe]);
		   	hEtaIDProbesMatched->Fill(idprobe_eta->at(itag)[iprobe]);
		   	hCentIDProbesMatched->Fill(centrality);
		   }
	
		}//iprobe 

	    // MS probes
		for(int iprobe = 0; iprobe < msprobe_n->at(itag); ++iprobe){
           // construct 4-vector of the ms probe

           TLorentzVector vMsProbe; vMsProbe.SetPtEtaPhiM(msprobe_pt->at(itag)[iprobe],msprobe_eta->at(itag)[iprobe],msprobe_phi->at(itag)[iprobe],0.10566);
           // add TP vectors
           TLorentzVector vSum = vMsProbe+vTag;
           // invariant TP mass
           double tpMass = vSum.M();
           if(tag_charge->at(itag)*msprobe_charge->at(itag)[iprobe]>0.0) {
               /*std::cout << "ms same-sign probe found." << std::endl;

               std::cout << "tag pt " << vTag.Pt() << std::endl;
               std::cout << "tag eta " << vTag.Eta() << std::endl;
               std::cout << "tag phi " << vTag.Phi() << std::endl;

               std::cout << "probe pt " << vMsProbe.Pt() << std::endl;
               std::cout << "probe eta " << vMsProbe.Eta() << std::endl;
               std::cout << "probe phi " << vMsProbe.Phi() << std::endl;
               std::cout << "mass: " << tpMass << std::endl;
               */

               // fill invariant mass distribution for same-signed pairs
               hTagMsProbeMass_SS->Fill(tpMass);
               if(msprobe_isMSProbeMatched->at(itag)[iprobe]==true) hTagMsProbeMassMatched_SS->Fill(tpMass);
               if(fabs(tpMass-91.1876)>10.0) continue;

               hPtMSProbes_SS->Fill(msprobe_pt->at(itag)[iprobe]);
               hEtaMSProbes_SS->Fill(msprobe_eta->at(itag)[iprobe]);
               hCentMSProbes_SS->Fill(centrality);

		       // same-sign probes matched to reconstructed muon
               if(msprobe_isMSProbeMatched->at(itag)[iprobe]==true){
                hPtMSProbesMatched_SS->Fill(msprobe_pt->at(itag)[iprobe]);
                hEtaMSProbesMatched_SS->Fill(msprobe_eta->at(itag)[iprobe]);
                hCentMSProbesMatched_SS->Fill(centrality);
               }

               continue;
           }

           // fill invariant mass distribution for opposite-signed pairs
           hTagMsProbeMass->Fill(tpMass);
           if(msprobe_isMSProbeMatched->at(itag)[iprobe]==true) hTagMsProbeMassMatched->Fill(tpMass);
           if(fabs(tpMass-91.1876)>10.0) continue;

           //opposite-signed TP pairs
           hPtMSProbes->Fill(msprobe_pt->at(itag)[iprobe]);
           hEtaMSProbes->Fill(msprobe_eta->at(itag)[iprobe]);
           hCentMSProbes->Fill(centrality);
		   // probes matched to reconstructed muon
		   if(msprobe_isMSProbeMatched->at(itag)[iprobe]==true){
		   	hPtMSProbesMatched->Fill(msprobe_pt->at(itag)[iprobe]);
		   	hEtaMSProbesMatched->Fill(msprobe_eta->at(itag)[iprobe]);
		   	hCentMSProbesMatched->Fill(centrality);
		   }
	
		}//iprobe 

	  }//itag
	
	}//iev

	std::cout << "No. of OS ID Probes" << " = " << hEtaIDProbes->Integral() << std::endl;
	std::cout << "No. of SS ID Probes" << " = " << hEtaIDProbes_SS->Integral() << std::endl;

	std::cout << "No. of OS MS Probes" << " = " << hEtaMSProbes->Integral() << std::endl;
	std::cout << "No. of SS MS Probes" << " = " << hEtaMSProbes_SS->Integral() << std::endl;
    
    nPtBins = hPtIDProbesMatched->GetNbinsX();
/*	  TGraphAsymmErrors* grMSEffEta = new TGraphAsymmErrors(); 
	  TGraphAsymmErrors* grMSEffPt = new TGraphAsymmErrors(); 
	  TGraphAsymmErrors* grIDEffEta = new TGraphAsymmErrors(); 
	  TGraphAsymmErrors* grIDEffPt = new TGraphAsymmErrors(); 
	  TGraphAsymmErrors* grTotEffEta = new TGraphAsymmErrors(); 
	  TGraphAsymmErrors* grTotEffPt = new TGraphAsymmErrors(); 
*/
//	  std::cout << hEtaMSProbesMatched->Integral() << std::endl;
//	  std::cout << hEtaMSProbes->Integral() << std::endl;

    // Efficiency with Bayesian errors
    // MS efficiency
    // Eta dependence
      calculateEfficiency(_MSEffEta, hEtaIDProbesMatched,hEtaIDProbesMatched_SS,hEtaIDProbes,hEtaIDProbes_SS,nEtaBins);
	  _MSEffEta->GetXaxis()->SetTitle("#eta");
	  double xEtaUpp = hEtaIDProbesMatched->GetXaxis()->GetXmax();	
	  double xEtaLow = -1.0*xEtaUpp;	
	  _MSEffEta->GetXaxis()->SetRangeUser(xEtaLow,xEtaUpp);

    // Pt dependence
    calculateEfficiency(_MSEffPt, hPtIDProbesMatched,hPtIDProbesMatched_SS,hPtIDProbes,hPtIDProbes_SS,nPtBins);
	_MSEffPt->GetXaxis()->SetTitle("p_{T}");
	double xPtLow = 20.0;
	double xPtUpp = hPtIDProbesMatched->GetXaxis()->GetXmax();
	_MSEffPt->GetXaxis()->SetRangeUser(xPtLow,xPtUpp);

    // centrality dependence
    calculateEfficiency(_MSEffCent, hCentIDProbesMatched,hCentIDProbesMatched_SS,hCentIDProbes,hCentIDProbes_SS,nCentralityBins);
	_MSEffCent->GetXaxis()->SetTitle("centrality[%] x 10^{-2}");
	double xCentLow = 0.0;
	double xCentUpp = hCentIDProbesMatched->GetXaxis()->GetXmax();
	_MSEffCent->GetXaxis()->SetRangeUser(xCentLow,xCentUpp);

    // ID efficiency
    // Eta dependence
      calculateEfficiency(_IDEffEta, hEtaMSProbesMatched,hEtaMSProbesMatched_SS,hEtaMSProbes,hEtaMSProbes_SS,nEtaBins);
	  _IDEffEta->GetXaxis()->SetTitle("#eta");
	  xEtaUpp = hEtaMSProbesMatched->GetXaxis()->GetXmax();	
	  xEtaLow = -1.0*xEtaUpp;	
	  _IDEffEta->GetXaxis()->SetRangeUser(xEtaLow,xEtaUpp);

    // Pt dependence
    calculateEfficiency(_IDEffPt, hPtMSProbesMatched,hPtMSProbesMatched_SS,hPtMSProbes,hPtMSProbes_SS,nPtBins);
	_IDEffPt->GetXaxis()->SetTitle("p_{T}");
	xPtLow = 20.0;
	xPtUpp = hPtMSProbesMatched->GetXaxis()->GetXmax();
	_IDEffPt->GetXaxis()->SetRangeUser(xPtLow,xPtUpp);

    // centrality dependence
    calculateEfficiency(_IDEffCent, hCentMSProbesMatched,hCentMSProbesMatched_SS,hCentMSProbes,hCentMSProbes_SS,nCentralityBins);
	_IDEffCent->GetXaxis()->SetTitle("p_{T}");
	xCentLow = 0.0;
	xCentUpp = hCentMSProbesMatched->GetXaxis()->GetXmax();
	_IDEffCent->GetXaxis()->SetRangeUser(xCentLow,xCentUpp);

	  //Multiply MS and ID efficiency to get total epsilon_rec
	  // Eta
	  for(int igr = 0; igr<_MSEffEta->GetN(); ++igr){

		double totEffEta = _MSEffEta->GetY()[igr]*_IDEffEta->GetY()[igr]; 

		double dEffEtaUp = TMath::Sqrt(TMath::Power(_MSEffEta->GetEYhigh()[igr]
						/_MSEffEta->GetY()[igr],2)
						+TMath::Power(_IDEffEta->GetEYhigh()[igr]
						/_IDEffEta->GetY()[igr],2))
						*totEffEta;

		double dEffEtaLow = TMath::Sqrt(TMath::Power(_MSEffEta->GetEYlow()[igr]
						/_MSEffEta->GetY()[igr],2)
						+TMath::Power(_IDEffEta->GetEYlow()[igr]
						/_IDEffEta->GetY()[igr],2))
						*totEffEta;

		double etaBinWidthUpp = _MSEffEta->GetEXhigh()[igr];
		double etaBinWidthLow = _MSEffEta->GetEXlow()[igr];

		_TotEffEta->SetPoint(igr,_MSEffEta->GetX()[igr], totEffEta);
		_TotEffEta->SetPointError(igr,etaBinWidthLow,etaBinWidthUpp,dEffEtaLow,dEffEtaUp);

	}//igr

	// Pt	
    for(int igr = 0; igr<_MSEffPt->GetN(); ++igr){

		double totEffPt = _MSEffPt->GetY()[igr]*_IDEffPt->GetY()[igr];

		double dEffPtUp = TMath::Sqrt(TMath::Power(_MSEffPt->GetEYhigh()[igr]
						/_MSEffPt->GetY()[igr],2)
						+TMath::Power(_IDEffPt->GetEYhigh()[igr]
						/_IDEffPt->GetY()[igr],2))*totEffPt;

		double dEffPtLow = TMath::Sqrt(TMath::Power(_MSEffPt->GetEYlow()[igr]
						/_MSEffPt->GetY()[igr],2)
						+TMath::Power(_IDEffPt->GetEYlow()[igr]
						/_IDEffPt->GetY()[igr],2))*totEffPt;

		double ptBinWidthUpp = _MSEffPt->GetEXhigh()[igr];
		double ptBinWidthLow = _MSEffPt->GetEXlow()[igr];

		_TotEffPt->SetPoint(igr,_MSEffPt->GetX()[igr], totEffPt);
		_TotEffPt->SetPointError(igr,ptBinWidthLow,ptBinWidthUpp,dEffPtLow,dEffPtUp);

	}//igr
	 

	// Centrality	
    for(int igr = 0; igr<_MSEffCent->GetN(); ++igr){

		double totEffCent = _MSEffCent->GetY()[igr]*_IDEffCent->GetY()[igr];

		double dEffCentUp = TMath::Sqrt(TMath::Power(_MSEffCent->GetEYhigh()[igr]
						/_MSEffCent->GetY()[igr],2)
						+TMath::Power(_IDEffCent->GetEYhigh()[igr]
						/_IDEffCent->GetY()[igr],2))*totEffCent;

		double dEffCentLow = TMath::Sqrt(TMath::Power(_MSEffCent->GetEYlow()[igr]
						/_MSEffCent->GetY()[igr],2)
						+TMath::Power(_IDEffCent->GetEYlow()[igr]
						/_IDEffCent->GetY()[igr],2))*totEffCent;

		double centBinWidthUpp = _MSEffCent->GetEXhigh()[igr];
		double centBinWidthLow = _MSEffCent->GetEXlow()[igr];

		_TotEffCent->SetPoint(igr,_MSEffCent->GetX()[igr], totEffCent);
		_TotEffCent->SetPointError(igr,centBinWidthLow,centBinWidthUpp,dEffCentLow,dEffCentUp);

	}//igr
	 

	outFile->cd();
	hEtaMSProbesMatched->Scale(1.0,"width");
	hEtaMSProbes->Scale(1.0,"width");
	hEtaIDProbes->Scale(1.0,"width");
	hEtaIDProbesMatched->Scale(1.0,"width");

	hEtaMSProbesMatched_SS->Scale(1.0,"width");
	hEtaMSProbes_SS->Scale(1.0,"width");
	hEtaIDProbes_SS->Scale(1.0,"width");
	hEtaIDProbesMatched_SS->Scale(1.0,"width");

	hPtMSProbesMatched->Scale(1.0,"width");
	hPtMSProbes->Scale(1.0,"width");
	hPtIDProbes->Scale(1.0,"width");
	hPtIDProbesMatched->Scale(1.0,"width");

	hPtMSProbesMatched_SS->Scale(1.0,"width");
	hPtMSProbes_SS->Scale(1.0,"width");
	hPtIDProbes_SS->Scale(1.0,"width");
	hPtIDProbesMatched_SS->Scale(1.0,"width");

	hCentMSProbesMatched->Scale(1.0,"width");
	hCentMSProbes->Scale(1.0,"width");
	hCentIDProbes->Scale(1.0,"width");
	hCentIDProbesMatched->Scale(1.0,"width");

	hCentMSProbesMatched_SS->Scale(1.0,"width");
	hCentMSProbes_SS->Scale(1.0,"width");
	hCentIDProbes_SS->Scale(1.0,"width");
	hCentIDProbesMatched_SS->Scale(1.0,"width");


	plotOnCanvas(hEtaMSProbesMatched,hEtaMSProbes,hEtaMSProbesMatched_SS,hEtaMSProbes_SS,_IDEffEta,
		     "cEtaMSProbes","dN/d#eta", "0-90%","MS Probes","MS Probes Matched",
		     "#eta^{#mu}","#epsilon_{rec}(ID)",outFile->GetName(),true);
	plotOnCanvas(hEtaIDProbesMatched,hEtaIDProbes,hEtaIDProbesMatched_SS,hEtaIDProbes_SS,_MSEffEta,
		     "cEtaIDProbes","dN/d#eta","0-90%","ID Probes","ID Probes Matched",
		     "#eta^{#mu}","#epsilon_{rec}(CB|ID)",outFile->GetName(),true);


	plotOnCanvas(hPtMSProbesMatched,hPtMSProbes,hPtMSProbesMatched_SS,hPtMSProbes_SS,_IDEffPt,
		     "cPtMSProbes","Muons/GeV","0-90%","MS Probes","MS Probes Matched","p_{T}^{#mu}","#epsilon_{rec}(ID)",outFile->GetName(),false,true,false);
	plotOnCanvas(hPtIDProbesMatched,hPtIDProbes,hPtIDProbesMatched_SS,hPtIDProbes_SS,_MSEffPt,
		     "cPtIDProbes","Muons/GeV","0-90%","ID Probes","ID Probes Matched","p_{T}^{#mu}","#epsilon_{rec}(CB|ID)",outFile->GetName(),false,true,false);


	plotOnCanvas(hCentMSProbesMatched,hCentMSProbes,hCentMSProbesMatched_SS,hCentMSProbes_SS,_IDEffCent,
		     "cCentMSProbes","Muons/GeV","0-90%","MS Probes","MS Probes Matched","centrality x10^{-2}","#epsilon_{rec}(ID)",outFile->GetName(),false,false,true);
	plotOnCanvas(hCentIDProbesMatched,hCentIDProbes,hCentIDProbesMatched_SS,hCentIDProbes_SS,_MSEffCent,
		     "cCentIDProbes","Muons/GeV","0-90%","ID Probes","ID Probes Matched","centrality x10^{-2}","#epsilon_{rec}(CB|ID)",outFile->GetName(),false,false,true);

	hEtaMSProbesMatched->Write();
	hEtaMSProbes->Write();
	hEtaIDProbes->Write();
	hEtaIDProbesMatched->Write();

	hEtaMSProbesMatched_SS->Write();
	hEtaMSProbes_SS->Write();
	hEtaIDProbes_SS->Write();
	hEtaIDProbesMatched_SS->Write();

	hPtMSProbesMatched->Write();
	hPtMSProbes->Write();
	hPtIDProbes->Write();
	hPtIDProbesMatched->Write();

	hPtMSProbesMatched_SS->Write();
	hPtMSProbes_SS->Write();
	hPtIDProbes_SS->Write();
	hPtIDProbesMatched_SS->Write();

	hCentMSProbesMatched->Write();
	hCentMSProbes->Write();
	hCentIDProbes->Write();
	hCentIDProbesMatched->Write();

	hCentMSProbesMatched_SS->Write();
	hCentMSProbes_SS->Write();
	hCentIDProbes_SS->Write();
	hCentIDProbesMatched_SS->Write();

    hTagMsProbeMass->Write();
    hTagIdProbeMass->Write();
    hTagMsProbeMass_SS->Write();
    hTagIdProbeMass_SS->Write();

    hTagMsProbeMassMatched->Write();
    hTagIdProbeMassMatched->Write();
    hTagMsProbeMassMatched_SS->Write();
    hTagIdProbeMassMatched_SS->Write();

	_MSEffEta->Write("MSEffEta");
	_IDEffEta->Write("IDEffEta");
	_TotEffEta->Write("TotEffEta");

	_MSEffPt->Write("MSEffPt");
	_IDEffPt->Write("IDEffPt");
	_TotEffPt->Write("TotEffPt");

	_MSEffCent->Write("MSEffCent");
	_IDEffCent->Write("IDEffCent");
	_TotEffCent->Write("TotEffCent");

	// Give back to free store
    delete hEtaMSProbesMatched;
	delete hEtaMSProbes;
	delete hEtaIDProbes;
	delete hEtaIDProbesMatched;
	delete hPtMSProbesMatched;
	delete hPtMSProbes;
	delete hPtIDProbes;
	delete hPtIDProbesMatched;

	delete hEtaMSProbesMatched_SS;
	delete hEtaMSProbes_SS;
	delete hEtaIDProbes_SS;
	delete hEtaIDProbesMatched_SS;
	delete hPtMSProbesMatched_SS;
	delete hPtMSProbes_SS;
	delete hPtIDProbes_SS;
	delete hPtIDProbesMatched_SS;

	delete hCentMSProbesMatched;
	delete hCentMSProbes;
	delete hCentIDProbes;
	delete hCentIDProbesMatched;

	delete hCentMSProbesMatched_SS;
	delete hCentMSProbes_SS;
	delete hCentIDProbes_SS;
	delete hCentIDProbesMatched_SS;

    delete hTagMsProbeMass;
    delete hTagIdProbeMass;

    delete hTagMsProbeMass_SS;
    delete hTagIdProbeMass_SS;

    delete hTagMsProbeMassMatched;
    delete hTagIdProbeMassMatched;

    delete hTagMsProbeMassMatched_SS;
    delete hTagIdProbeMassMatched_SS;

	delete _MSEffEta;
	delete _IDEffEta;
	delete _TotEffEta;
	delete _MSEffPt;
	delete _IDEffPt;
	delete _TotEffPt;

}
