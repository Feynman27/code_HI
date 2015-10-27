#include "TCanvas.h"
#include "TDirectory.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2DErrors.h"
#include "TList.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TChain.h"
#include "THStack.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TMath.h"
#include "TLatex.h"
#include "TClonesArray.h"
#include "TUnixSystem.h"

#include "RooDataSet.h"
#include "RooBinning.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHist.h"

#include "WAnalysisHIDep.C"
#include "AtlasUtils.C"
#include "AtlasStyle.C"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath> 
#include <sstream>

#include "EfficiencyCorrection.C"
#include "Systematics.C"
#include "WPlotterHelper.C"

using namespace RooFit ;
///////////////////////////////////////////////////
//Scale weighted preselected histogram by RAA 
//which corresponds to the current centrality bin
////////////////////////////////////////////////////
TH1F* scaleByRaa(TH1F* hQCD, TH1F* hRAA, double ptSigLo, int icent){
    
    std::cout << "Scaling preselected histogram by the R_AA of centrality bin " << icent << std::endl;
    
        ///Starting bin in signal region from Petr's RAA figure (pt = 21.35 GeV)
        //int binRAALo = 38; //pt = 19.9 - 22.8 GeV
        int binRAALo = 1; //pt = 0.09 GeV
        int binRAAUpp = 49;  //pt bin center = 213.5 GeV

        ///Scale by RAA corresponding to this pt region; 
        ///Start at the bin corresponding to the signal region
        for(int iRAA = binRAALo; iRAA<=binRAAUpp; ++iRAA){

            ///find current pt bin center in RAA histo 
            double ptLo = hRAA->GetXaxis()->GetBinCenter(iRAA);
            double ptUpp = ptLo;
            ///cut off at lower bin width
            ptLo -= 0.5*(hRAA->GetXaxis()->GetBinWidth(iRAA));
            ///cut off at upper bin width
            ptUpp += 0.5*(hRAA->GetXaxis()->GetBinWidth(iRAA));
            std::cout << "Searching RAA in pt region " << ptLo << "-" << ptUpp << std::endl;

            ///Loop over each bin in the preselected,scaled, and weighted QCD pt spectrum
            for(int iQCD=1; iQCD<=hQCD->GetNbinsX(); ++iQCD){
        
                ///Get the pt value from the QCD histo
                double ptQCD = hQCD->GetXaxis()->GetBinCenter(iQCD);

                ///Find where the bin center in the nominal QCD spectrum
                ///falls in the RAA spectrum
                if( (ptQCD>ptLo)  && (ptQCD < ptUpp) ){

                    std::cout << "Bin Found " << "\n" << 
                        "QCD muon pt for bin " << iQCD << " before Raa scaling = " << ptQCD << std::endl;
                    ///Retrieve the RAA in this pt bin
                    double RAA = hRAA->GetBinContent(iRAA);
                    double muonsInPtBin = hQCD->GetBinContent(iQCD);
                    ///Scale this bin by the corresponding RAA
                    double scaledPt = muonsInPtBin*RAA;
                    hQCD->SetBinContent(iQCD,scaledPt);
                    std::cout << "QCD muon pt for bin " << iQCD << " after Raa scaling = " << hQCD->GetBinContent(iQCD) << std::endl;
                }
                else continue;
            }
        }

    ///Return the RAA-scaled preselected histogram
    return hQCD;
}

///////////////////////////////////////////////////////////////////////////////
//plot difference in charge as function of Ncoll 
///////////////////////////////////////////////////////////////////////////////
void plotChargeDiffNcoll(TGraph* grDiff, TGraphErrors* grPlus, TGraphErrors* grMinus, int ieta, int icent, double npart){

	double* yPlus = grPlus->GetY();
	double* yMinus = grMinus->GetY();

	double yTempPlus = yPlus[icent]; 
	double yTempMinus = yMinus[icent]; 

	double diff = fabs(yTempPlus-yTempMinus);
	std::cout << "Difference in f_{QCD} at eta bin " << ieta << " and centrality bin " << icent << " = " << diff << std::endl;
	grDiff->SetPoint(icent,npart,diff);

}

///////////////////////////////////////////////////////////////////////////////
//plot difference in charge as function of eta 
///////////////////////////////////////////////////////////////////////////////
void plotChargeDiffEta(TGraph* grDiff, TGraphErrors* grPlus, TGraphErrors* grMinus, int ieta, int icent, double etaLo, double etaUpp){

	double  xPt = etaLo+fabs(etaUpp-etaLo)/2.0;
	double* yPlus = grPlus->GetY();
	double* yMinus = grMinus->GetY();

	double yTempPlus = yPlus[ieta]; 
	double yTempMinus = yMinus[ieta]; 

	double diff = fabs(yTempPlus-yTempMinus);
	std::cout << "Difference in f_{QCD} at eta bin " << ieta << " and centrality bin " << icent << " = " << diff << std::endl;
	grDiff->SetPoint(ieta,xPt,diff);

}

void PrintScaledBkg(TGraphErrors* grFracNew, TString sCh,TString sEta){
	TCanvas* cFracSc = new TCanvas("cFracSc","cFracSc",600,600);
	grFracNew->GetXaxis()->SetTitle("|#eta|");	
	grFracNew->GetYaxis()->SetTitle("#frac{N_{QCD}^{sig}}{N_{W}}");
	grFracNew->Draw("ape");	
	TLatex l;
	l.SetNDC();
	l.DrawLatex(0.39,0.83,sCh);
	cFracSc->Update();

	TString sLabel = sCh+sEta;
//	cFracSc->Print(sLabel+".pdf");
}

void Write(TFile* outFile, TObject* grFrac, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing TGraph to root file..." << std::endl;
	  grFrac->Write(sFile);
	  std::cout << "Done." << std::endl;
    	  gDirectory = dir;
}

void fitToPol0(TGraphErrors* gr, double xLow, double xUpp){
    TF1* f0 = new TF1("f0","[0]",xLow,xUpp);
    gr->Fit(f0);
}

void plot(TFile* outFile, RooRealVar& muonPt,RooDataSet* dataSet, RooDataSet* mcJ1Set, RooDataSet* mcJ2Set, RooDataSet* mcJ3Set, 
		RooDataSet* psDataSet, RooDataSet* psJ1Set, RooDataSet* psJ2Set, RooDataSet* psJ3Set,
		RooDataSet* mcJ1Events,RooDataSet* mcJ2Events, RooDataSet* mcJ3Events, TH1F* hRAA, 
        TGraphErrors* grFracPt, TGraphErrors* grFrac, TGraphErrors* grFracEta, TGraphErrors* grFracIn,
		double& arrQCDFrac, double& arrSumQCD, double ptLow, double ptUpp, double etaLow, double etaUpp, 
		double centralityLow, double centralityUpp, double ncoll, int ich, int ipt, int icent, int ieta, 
		TString sSel, TString sSel2, double bkgFracInt, bool doScaleByRaa=false, bool doScale=false){

		int nbins = 400;
		double ptmax = 400.0;

		double sigma = 8.5;
		double ptCtrlLo = 10.0; double ptCtrlUpp = 20.0; double ptSigLo = 25.0;  
		if(ptSigLo<25.0) std::cout << "WARNING: Signal lower pT set to " << ptSigLo << std::endl;
		//define bins for control region
//		int ptCtrlBinLo = (nbins/ptmax)*ptCtrlLo+1; int ptCtrlBinUpp = (nbins/ptmax)*ptCtrlUpp+1; 
		//define bins for signal region
//		int ptSigBinLo = (nbins/ptmax)*ptSigLo+1; int ptSigBinUpp = nbins+1;  

		std::cout << "centrality: " << centralityLow << "-" << centralityUpp << " ncoll: " << ncoll << std::endl;

  		RooBinning b = RooBinning(nbins,0.0,ptmax); // 1 GeV per bin

		//hard probes
		double evData = 1.03e9;

		//initialize histograms	
		// --- data ---
		TH1F* hdataSet = (TH1F*)dataSet->createHistogram("hdataSet",muonPt,Binning(b));
		TH1F* hdataSetPS = (TH1F*)psDataSet->createHistogram("hdataSetPS",muonPt,Binning(b));
		// --- QCD set ---
		TH1F* hmcJ1Set = (TH1F*)mcJ1Set->createHistogram("hmcJ1Set",muonPt,Binning(b));
		TH1F* hmcJ2Set = (TH1F*)mcJ2Set->createHistogram("hmcJ2Set",muonPt,Binning(b));
		TH1F* hmcJ3Set = (TH1F*)mcJ3Set->createHistogram("hmcJ3Set",muonPt,Binning(b));
		// --- preselected --- 
		TH1F* hmcJ1SetPS = (TH1F*)psJ1Set->createHistogram("hmcJ1SetPS",muonPt,Binning(b));
		TH1F* hmcJ2SetPS = (TH1F*)psJ2Set->createHistogram("hmcJ2SetPS",muonPt,Binning(b));
		TH1F* hmcJ3SetPS = (TH1F*)psJ3Set->createHistogram("hmcJ3SetPS",muonPt,Binning(b));
                int ptSigBinLo = hdataSet->FindBin(ptSigLo);
                int ptSigBinUpp = hdataSet->FindBin(ptmax);
                int ptCtrlBinLo = hdataSetPS->FindBin(ptCtrlLo);
		int ptCtrlBinUpp = hdataSetPS->FindBin(ptCtrlUpp);
		std::cout << "Consistency checks: " << std::endl;
		std::cout << ptSigBinLo << "=?" << (nbins/ptmax)*ptSigLo+1 << std::endl;
                std::cout << ptSigBinUpp << "=?"<< nbins << std::endl;
                std::cout << ptCtrlBinLo << "=?"<< (nbins/ptmax)*ptCtrlLo+1 << std::endl;
                std::cout << ptCtrlBinUpp << "=?"<< (nbins/ptmax)*ptCtrlUpp+1 << std::endl;

        ///error on QCD from each Jx sample before any weighting or scaling
        double absoluteErrMC =
        TMath::Sqrt(hmcJ1Set->Integral(ptSigBinLo,ptSigBinUpp)+hmcJ2Set->Integral(ptSigBinLo,ptSigBinUpp)
            + hmcJ3Set->Integral(ptSigBinLo,ptSigBinUpp));
        double jxSum = hmcJ1Set->Integral(ptSigBinLo,ptSigBinUpp)+hmcJ2Set->Integral(ptSigBinLo,ptSigBinUpp)
            + hmcJ3Set->Integral(ptSigBinLo,ptSigBinUpp);
        ///relative error on MC after W selection
        double relativeErrMC = absoluteErrMC/jxSum;
        std::cout << "Primary statistic on J1+2+3mu sample: " << jxSum << "+/-" << absoluteErrMC << "(" << relativeErrMC*100.0 << "%)" << std::endl;

        std::cout << "Control Region:" << hdataSet->GetXaxis()->GetBinLowEdge(ptCtrlBinLo) << "-" << hdataSet->GetXaxis()->GetBinLowEdge(ptCtrlBinUpp) << std::endl;
        std::cout << "Signal Region:" << hdataSet->GetXaxis()->GetBinLowEdge(ptSigBinLo) << "-" << hdataSet->GetXaxis()->GetBinLowEdge(ptSigBinUpp)<< std::endl;

		std::cout << "HF muons before W Selection: " << hmcJ1SetPS->Integral() << " " << hmcJ2SetPS->Integral() << " " << hmcJ3SetPS->Integral() << std::endl;
		std::cout << "Muons in data before W selection " << hdataSetPS->Integral() << std::endl;

		//return correctly weighted QCD histogram
		TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",nbins,0.0,ptmax);
		hmcQCDSet = getWeightedJxHisto(centralityLow, centralityUpp, ncoll, hmcJ1Set, hmcJ2Set, hmcJ3Set, mcJ1Events, mcJ2Events,mcJ3Events, nbins,0.0, ptUpp);

		//preselcted
		TH1F* hmcQCDSetPS = new TH1F("hmcQCDSetPS","hmcQCDSetPS",nbins,0.0,ptmax);
		hmcQCDSetPS = getWeightedJxHisto(centralityLow, centralityUpp, ncoll, hmcJ1SetPS, hmcJ2SetPS, hmcJ3SetPS, mcJ1Events, mcJ2Events,mcJ3Events, nbins,0.0, ptUpp);

		int binLo = 1;

		TH1F* hdataSetc = (TH1F*)hdataSet->Clone("hdataSetc");
		TH1F* hdataSetPSc = (TH1F*)hdataSetPS->Clone("hdataSetPSc");
		TH1F* hmcQCDSetc = (TH1F*)hmcQCDSet->Clone("hmcQCDSetc");
		TH1F* hmcQCDSetPSc = (TH1F*)hmcQCDSetPS->Clone("hmcQCDSetPSc");
		
  		double sigEventsPS = hdataSetPSc->Integral(binLo,nbins); 
  		std::cout << "integrated events before W selection:"<< " = " << sigEventsPS << " +-" << TMath::Sqrt(sigEventsPS) << std::endl;
		double sigEvents = hdataSetc->Integral(binLo,nbins); 
  		std::cout << "integrated events after W selection:" << " = " << sigEvents << " +-" << TMath::Sqrt(sigEvents) << std::endl;
	 
		//scale Jx to uncorrected data in control region (10-20GeV);
		//extrapolation to signal region affords QCD bkg estimate
		//before W selection cuts
		std::cout << "Integrating pre-selected samples from bins " << ptCtrlBinLo << "-" << ptCtrlBinUpp << std::endl;
		double scale = hdataSetPS->Integral(ptCtrlBinLo,ptCtrlBinUpp)/hmcQCDSetPSc->Integral(ptCtrlBinLo,ptCtrlBinUpp);

		std::cout << "Total number of pre-selected muons in data = " << hdataSetPS->Integral() << std::endl;
		std::cout << "Total number of pre-selected muons in weighted Jx b4 scaling = " << hmcQCDSetPSc->Integral() << std::endl;
		std::cout << "Scale factor = " << hdataSetPS->Integral(ptCtrlBinLo,ptCtrlBinUpp) << "/" << hmcQCDSetPSc->Integral(ptCtrlBinLo,ptCtrlBinUpp)
			<< " = " << scale << std::endl;
        
        ///Here we choose to scale the pp exected distribution to the data (used in the analysis)
        ///or by the Raa (used for systematics)
        double preSelectedQCD=0.,preSelectedData=0.,bkgQCDPreSel=0.;
        if(doScaleByRaa){
            hmcQCDSetPSc = scaleByRaa(hmcQCDSetPSc,hRAA,ptCtrlUpp,icent);
		    std::cout << "Total number of pre-selected muons in weighted Jx after RAA scaling = " << hmcQCDSetPSc->Integral() << std::endl;
            preSelectedQCD = hmcQCDSetPSc->Integral(ptSigBinLo,ptSigBinUpp);
            bkgQCDPreSel = preSelectedQCD/preSelectedData*100.0;
		    std::cout << "QCD background estimate for centrality class " << centralityLow << "-" << centralityUpp << " and eta class " << etaLow << "-" << 
			    etaUpp << " BEFORE W selection and after RAA scaling = " << bkgQCDPreSel << "%" << std::endl; 

        }
		else {
            hmcQCDSetPSc->Scale(scale);
		    std::cout << "Total number of pre-selected muons in weighted Jx after scaling = " << hmcQCDSetPSc->Integral() << std::endl;
            preSelectedQCD = hmcQCDSetPSc->Integral(ptSigBinLo,ptSigBinUpp);
            bkgQCDPreSel = preSelectedQCD/preSelectedData*100.0;
            std::cout << "QCD background estimate for centrality class " << centralityLow << "-" << centralityUpp << 
            "and eta class " << etaLow << "-" << etaUpp << " BEFORE W selection = " << bkgQCDPreSel << "%" << std::endl;

        }

		double arrCentWidth = centralityUpp-centralityLow;
		double scaleToData = arrCentWidth*ncoll*evData;

		//ratio of Jx+1mu filter cross-section to total pp cross-section; takes into account probability to find muon with a jet (AMI)
		double wtJ1 = 6.6523e-03*1.8770e-4/64.0e-3; double wtJ2 = 1.4941e-02*8.2788e-6/64.0e-3; double wtJ3 = 2.4284e-02*2.9426e-7/64.0e-3;

		std::cout << "Integral of Jx Sets from bins " << ptSigBinLo << "-" << ptSigBinUpp << "= " << 
			hmcJ1SetPS->Integral(ptSigBinLo,ptSigBinUpp) << " " << hmcJ2SetPS->Integral(ptSigBinLo,ptSigBinUpp) << " " << hmcJ1SetPS->Integral(ptSigBinLo,ptSigBinUpp) << std::endl;

		std::cout << "weights = " << wtJ1 << " " << wtJ2 << " " << wtJ3 << std::endl;
		//std::cout << "Jx entries : " << (double)mcJ1Events->numEntries() << " " << (double)mcJ2Events->numEntries() << " " << (double)mcJ3Events->numEntries() << std::endl;

		double nMCJ1 = getNEvents(mcJ1Events, centralityLow, centralityUpp);
		double nMCJ2 = getNEvents(mcJ2Events, centralityLow, centralityUpp);
		double nMCJ3 = getNEvents(mcJ3Events, centralityLow, centralityUpp);

		std::cout << "Jx entries : " << nMCJ1 << " " << nMCJ2 << " " << nMCJ3 << std::endl;
		std::cout << "Scale factor to data = " << scaleToData << std::endl;

		//scale by number of events in data in ctrl region
		std::cout << "Scale factor for QCD MC in ctrl region = " << scale << std::endl;
		//hmcQCDSetc->Scale(scale);
		
        ///temp hack since Petr only provided
        ///0-5% RAA at the moment
        if(doScaleByRaa/*&&icent==0&&centralityUpp<0.11*/){
		    std::cout << "Total number of surviving muons in weighted Jx before RAA scaling = " << hmcQCDSetc->Integral() << std::endl;
            hmcQCDSetc = scaleByRaa(hmcQCDSetc,hRAA,ptCtrlUpp,icent);
		    std::cout << "Total number of surviving muons in weighted Jx after RAA scaling = " << hmcQCDSetc->Integral() << std::endl;
        }
        else hmcQCDSetc->Scale(scale);

		double fractionQCD = hmcQCDSetc->Integral(ptSigBinLo,ptSigBinUpp)/hdataSet->Integral(ptSigBinLo,ptSigBinUpp);
		std::cout << "Number of QCD muons surviving W cuts = " << hmcQCDSetc->Integral(ptSigBinLo,ptSigBinUpp) << std::endl;
		std::cout << "Number of W candidates = " << hdataSet->Integral(ptSigBinLo,ptSigBinUpp) << std::endl;
        
        ///Calculate sqrt(N) errors
        double relativeErrData = TMath::Sqrt(hdataSet->Integral(ptSigBinLo,ptSigBinUpp))/hdataSet->Integral(ptSigBinLo,ptSigBinUpp);
        ///add in quadrature
        double errFractionQCD = fractionQCD*TMath::Sqrt(TMath::Power(relativeErrMC,2)+TMath::Power(relativeErrData,2));

		//scale to fraction of QCD in the current bin
		//and write to outFile
		if(doScale){

			TGraphErrors* grFracNew = new TGraphErrors();

			double* ytemp = ((TGraphErrors*)grFracIn)->GetY();
			double* ytempErr = ((TGraphErrors*)grFracIn)->GetEY();
			double* xtemp = ((TGraphErrors*)grFracIn)->GetX();

			for(int igr=0; igr<grFracIn->GetN(); igr++){
			    
                double etaBw = abs(xtemp[igr+1]-xtemp[igr])/2.0;
	
				std::cout << "QCD background pre-scaled for eta " << xtemp[igr]-(etaBw) << "-" << xtemp[igr]+(etaBw) << " = " 
					<< ytemp[igr] << std::endl;

				double scaleFactor = fractionQCD/bkgFracInt;
				std::cout << "Background s.f. = " << fractionQCD << "/" << bkgFracInt << " = " << scaleFactor << std::endl;
				//std::cout << "N_bkg before applying s.f. = " << ytemp[igr] << std::endl;
				double yScaled = grFracIn->GetY()[igr]*scaleFactor; 
				double yScaledErr = grFracIn->GetEY()[igr]*scaleFactor;
				std::cout << "QCD background for eta " << xtemp[igr]-(etaBw) << "-" << xtemp[igr]+(etaBw) 
					<< " after scaling to centrality point " << centralityLow << "-" << centralityUpp << " = " 
					<< yScaled << "+-" << yScaled << std::endl;

				grFracNew->SetPoint(igr,xtemp[igr],yScaled);
				grFracNew->SetPointError(igr,0.0,yScaledErr);

			}

				PrintScaledBkg(grFracNew,sSel,sSel2);
				TString sGraphName = "fractionQCD_charge";sGraphName+=ich; sGraphName+="_eta"; sGraphName+=ieta; sGraphName+="_cent"; sGraphName+=icent;
				Write(outFile, grFracNew,sGraphName);

		}

		//wt'd fraction of QCD muons in signal region for given eta/cent slice
        double weight = hmcQCDSetc->Integral(ptSigBinLo,ptSigBinUpp)+hdataSetc->Integral(ptSigBinLo,ptSigBinUpp);
		double wtdFrac = weight*fractionQCD;

		//running sum of fraction per centrality class
		arrQCDFrac += wtdFrac;
		std::cout << "Running (wt'd) sum of QCD fraction  = " << arrQCDFrac << std::endl;

		//running sum per centrality classs
		arrSumQCD += weight;
		std::cout << "Running sum of QCD events " << arrSumQCD << std::endl;

		std::cout << "QCD background estimate for centrality class " << centralityLow << "-" << centralityUpp << " and eta class " << etaLow << "-" << 
			etaUpp <<" = " << fractionQCD*100. << "%" << std::endl; 

		double xpt = ncoll ;
		grFrac->SetPoint(icent,xpt,fractionQCD);
		grFrac->SetPointError(icent,0.0,errFractionQCD);
		//grFrac->SetPointError(icent,0.0,eFrac*fractionQCD);

		xpt = (etaUpp-etaLow)/2.0+etaLow;
		grFracEta->SetPoint(ieta,xpt,fractionQCD);
		grFracEta->SetPointError(ieta,0.0,errFractionQCD);
		//grFracEta->SetPointError(ieta,0.0,eFrac*fractionQCD);

		xpt = (ptUpp-ptLow)/2.0+ptLow;
		grFracPt->SetPoint(ipt,xpt,fractionQCD);
		grFracPt->SetPointError(ieta,0.0,errFractionQCD);
		//grFracPt->SetPointError(ieta,0.0,eFrac*fractionQCD);

		TLatex l;
		l.SetNDC();
		TCanvas* cpt = new TCanvas("cdatapt","cdatapt",600,600);
		hdataSetPSc->GetXaxis()->SetTitle("p_{T}[GeV]"); 
		hdataSetPSc->GetYaxis()->SetTitle("Muons/GeV");
        if(doScaleByRaa) hdataSetPSc->GetXaxis()->SetRangeUser(10.0,100.0);
		else hdataSetPSc->GetXaxis()->SetRangeUser(0.0,100.0);
		hdataSetPSc->GetYaxis()->SetRangeUser(0.1,1.0e6);
		hmcQCDSetPSc->SetLineColor(kRed);
		hdataSetPSc->SetMarkerColor(kRed);
		hmcQCDSetc->SetLineColor(kBlue);
		hdataSetc->SetMarkerColor(kBlue);
//		hmcQCDSetPSc->SetMarkerStyle(24);
		//hmcWSetPS->SetLineColor(kBlue);
//		hdataSetPSc->SetLineColor(kBlue);
		hdataSetPSc->Draw("pe ");
		hmcQCDSetPSc->Draw("hist same");
		hdataSetPSc-> SetMinimum(0.01);
		hdataSetPSc-> SetMaximum(7.0e5);

		TLegend* leg = new TLegend(0.56, 0.84, 0.83, 0.67);
		leg->SetTextFont(gStyle->GetTextFont());
		leg->SetTextSize(gStyle->GetTextSize());
		leg->SetBorderSize(0);
		leg->SetFillColor(0);
		leg->AddEntry(hdataSetPSc,"data (PS)","p");
		leg->AddEntry(hmcQCDSetPSc,"J1-3 (PS)","l");
		leg->AddEntry(hdataSetc,"data (W Sel)","p");
		leg->AddEntry(hmcQCDSetc,"J1-3 (W Sel)","l");
		leg->Draw();
//		hmcWSetPS->Draw("hf same");
		double ptCutLine = 25.0;
		cpt->Update();
		//cpt->Range( 0., 10., 0., 1.0e6. );
		TLine *line0 = new TLine(ptCutLine,0.01,ptCutLine,881201);
		line0->SetLineColor(kBlack);
		line0->SetLineStyle(kDashed);
		line0->SetLineWidth(2);
		line0->Draw();
		//l.DrawLatex(0.54,0.87,sSel);
		//l.DrawLatex(0.75,0.87,sSel2);
		leg->Draw();
		cpt->SetLogy(); cpt->Update();
//		cpt->Print("cptQCDBkgEta"+ieta+"Centrality"+icent+".pdf");
//		cpt->Print("cptQCDBkgEta"+ieta+"Centrality"+icent+".root");

		/*TCanvas* cptW = new TCanvas("cdataptW","cdataptW",600,600);
		hdataSet->GetXaxis()->SetTitle("p_{T}[GeV]"); 
		hdataSet->GetYaxis()->SetTitle("Muons/GeV");
		hdataSet->GetXaxis()->SetRangeUser(0.0,100.0);
		*/
//		hmcQCDSetc->SetLineColor(kBlue+2);
//		hmcQCDSetW->SetMarkerStyle(24);
//		hdataSetW->SetLineColor(kBlue);
		hdataSetc->Draw("pe same");
//		hdataSet->SetMinimum(0.0001);
//		hdataSet->SetMaximum(1.0e5);

//		cptW->Update();
		cpt->Update();
		//cptW.Range( 0., 100., 0., 1.0e4. );
		TLine* line = new TLine(ptCutLine,0.0001,ptCutLine,8.0e3);
		line->SetLineColor(kBlack);
		line->SetLineStyle(kDashed);
		line->SetLineWidth(2);
//		line->Draw();
    
		hmcQCDSetc->Draw("hist same");
        gStyle->SetPaintTextFormat("4.1f");
		l.DrawLatex(0.54,0.87,sSel);
		//l.DrawLatex(0.75,0.87,sSel2);
		l.DrawLatex(0.75,0.87,"0.1<|#eta|<2.4");
		leg->Draw();
		cpt->Update();

        TString sSaveName = "cptQCDBkgCharge"; sSaveName+=ich;
        sSaveName+="Eta"; sSaveName+=ieta;
        sSaveName+="Centrality"; sSaveName+=icent;
		cpt->Print(sSaveName+".pdf");
		cpt->Print(sSaveName+".root");
		//cptW->Print("cptW,"+sSel+".pdf");
//		cptW->SetLogy(); cptW->Update();
//		cptW->Print("cptWCuts,"+sSel+","+sSel2+",Log.pdf");

//		cpt->Print("cptWCuts,"+sSel+","+sSel2+",Log.pdf");

} //plot


void plotFraction(TGraphErrors* grFrac, TGraphErrors* grFracEta, TString sCh){
	TLatex l;
	TCanvas* cFrac = new TCanvas("cFrac","cFrac",600,600);
	l.SetNDC();
	l.DrawLatex(0.39,0.83,sCh);
	grFrac->GetXaxis()->SetTitle("#LT N_{coll} #GT");	
	grFrac->GetYaxis()->SetTitle("#frac{N_{QCD}^{sig}}{N_{W}}");
	grFrac->Draw("ape");	
	cFrac->Update();
//	cFrac->Print(sCh+"_centrality.pdf"); cFrac->Print(sCh+ "_centrality" + ".root");

	TCanvas* cFracEta = new TCanvas("cFracEta","cFracEta",600,600);
	l.SetNDC();
	l.DrawLatex(0.39,0.83,sCh);
	grFracEta->GetXaxis()->SetTitle("|#eta|");	
	grFracEta->GetYaxis()->SetTitle("#frac{N_{QCD}^{sig}}{N_{W}}");
	grFracEta->Draw("ape");	
	cFracEta->Update();
//	cFracEta->Print(sCh+"_eta.pdf"); cFracEta->Print(sCh+"_eta" + ".root");
}

void fitFraction(TGraphErrors* gr, TString sName){

	std::cout << "Fitting QCD bkg fraction to analytic function. " << std::endl;
	TF1* f0a = new TF1("f0a","[0]+[1]*x",0.0,2000.0);	
        gr->Fit("f0a");
        f0a->Draw("same");
        std::cout << "chi2 = " << f0a->GetChisquare() << "/" << f0a->GetNDF() << ", p = " << f0a->GetProb() << std::endl;

}

/////////////////////////////////////////////
//plot residual between fit and Ncoll distribution
/////////////////////////////////////////////
void plotResidualNcoll(TGraphErrors* gr, TGraph* grResid){

	//double xPt = etaLo + 0.5*fabs(etaUpp-etaLo);
	TF1* f = gr->GetFunction("f0a");

	double* xPt = gr->GetX();
	double* yPt = gr->GetY();

	for(int igr=0; igr<gr->GetN(); igr++){

		std::cout << "x,y = " << xPt[igr] << "," << yPt[igr] << std::endl;
		double resid = fabs(yPt[igr] - f->Eval(xPt[igr])); 
		grResid->SetPoint(igr,xPt[igr],resid);
	}
}



void bkgQCDPlotter()
{
	double cutValue = 11;
	//switch for scaling eta fraction to centrality fraction
	//produce ete distro first and change fileNameDataGrIn accordingly
	bool doScale = true; 
    bool doScaleByRaa = false;
	bool doCentrality = true;
	bool doPt = false;
	bool doEta = false;
	bool doCharge = true;

    bool doIsolationOptimization = false;
    bool doMptSigmaDown = false;
    bool doMptSigmaUp = false;
    ///Note: remember to change WAnalysisHIDep file accordingly
    bool doIncreaseConeSize = false;
    bool doLoosenIsoCut = false;
    if(doMptSigmaUp) std::cout << "Doing MPT +1 sigma systematics." << std::endl;
    if(doMptSigmaDown) std::cout << "Doing MPT -1 sigma systematics." << std::endl;
    if(doIsolationOptimization) std::cout << "Calculating QCD background with isolation optimization cuts." << std::endl;
    if(doIncreaseConeSize) std::cout << "Increasing cone size from 0.2 to 0.3 for systematic study." << std::endl;
    if(doLoosenIsoCut) std::cout << "Loosening isolation cut from 0.1 to 0.2 for systematic study." << std::endl;
    
    SetAtlasStyle();
	TString date = "Jan20.2012";
//	TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
	TString baseString = "/usatlas/u/tbales/scratch/";
//	TString baseString = "/tmp/tbalestr/";

	//hard probes
//	TString fileNameDataIn = baseString+"HardProbesFiles/HISingleMuonHP.12.19.2012";
    ///Nominal input for main analysis
	TString fileNameDataIn = "HISingleMuonHardProbesData.04.17.2013";
    if(doMptSigmaDown||doMptSigmaUp) fileNameDataIn = "HISingleMuonHardProbesData.07.13.2013";
    std::cout << "Input file for data: " << fileNameDataIn << std::endl;

//	TString fileNameDataIn = "HardProbesFiles/HISingleMuonHP.12.19.2012";
	TString fileNameDataGrIn ;
	//TString fileNameDataGrChargeIn ;
	TString fileNameDataOut;
    TString fileNameRaa = "tb_raa2";

	//input eta distro for scaling cent distro
//	fileNameDataGrIn = "fractionQCDEta.01.20.2013";
//	fileNameDataGrIn = "fractionQCDEta.03.30.2013";

    ///use for analysis eta binning
    //fileNameDataGrIn = "background/fractionQCDEta_9etaBins6centBins.04.21.2013";
    //fileNameDataGrIn = "background/fractionQCDEta_9Bins.06.17.2013";
    fileNameDataGrIn = "background/fractionQCDEta_9Bins.11.27.2013";
    //cross-check w/o mpt cut
    //fileNameDataGrIn = "crossChecks/fractionQCDEta_noMptCut.12.10.2013";
    ///w/o isolation cut applied
//    fileNameDataGrIn = "background/fractionQCDEta_NoIsoCut_9Bins.06.24.2013";
    ///mpt systematics -1sigma
    if(doMptSigmaDown)fileNameDataGrIn = "systematics/fractionQCDEta_2GeVMpt_systematics.08.03.2013";
    ///mpt systematics +1sigma
    if(doMptSigmaUp)fileNameDataGrIn = "systematics/fractionQCDEta_4GeVMpt_systematics.08.03.2013";
    //use for ptcone30ID3<0.1 systematics
	if(doIncreaseConeSize)fileNameDataGrIn = "systematics/fractionQCDEta_ptcone30ID3_systematics.06.25.2013";
    //use for ptcone20ID3<0.2 systematics
	if(doLoosenIsoCut)fileNameDataGrIn = "systematics/fractionQCDEta_ptcone20ID3_systematics.06.25.2013";

	//use this file for systematics
//	fileNameDataGrIn = "fractionQCDEta.02.12.2013";
//	fileNameDataGrIn = "fractionQCDEta_1sigmaUpPt.02.08.2013";
//	fileNameDataGrIn = "fractionQCDEta_1sigmaDownPt.02.08.2013";
//	fileNameDataGrIn = "fractionQCDEta_1sigmaDownMPt.02.08.2013";
//	fileNameDataGrIn = "fractionQCDEta_1sigmaUpMPt.02.08.2013";

	if(doScale) fileNameDataOut = "fractionQCDEtaScaled";
	else if(doEta) fileNameDataOut = "fractionQCDEta";
	else if(doCentrality) fileNameDataOut = "fractionQCDCent";
	else fileNameDataOut = "fractionQCD";
	//minbias
        //TString fileNameDataIn = baseString+"MinimumBiasFiles/HISingleMuonMB.2012.11.25";
	TFile* fGr = 0;

	/// --- Open output file ---
	TDirectory *dir = gDirectory;
	TFile *outFile = new TFile(fileNameDataOut+".root","RECREATE");
	gDirectory = dir;

    fGr = new TFile(fileNameDataGrIn+".root", "READ");
    if ( !fGr->IsOpen() ) {
    std::cout << fGr << " not found!" << std::endl;
    }

    std::cout <<  fileNameDataGrIn+".root" << " opened." << std::endl;

    TFile* fRaa = new TFile(fileNameRaa+".root", "READ");

/*	if(doCharge){
		fGrCharge = new TFile(fileNameDataGrChargeIn+".root", "READ");
		if ( !fGrCharge->IsOpen() ) {
		std::cout << fGrCharge << " not found!" << std::endl;
		return 1;
  		}

  		cout <<  fileNameDataGrChargeIn+".root" << " opened." << endl;
	}
*/
	
//	else {
//		delete fGr;
//	}

	//get qcd frac as fcn of eta graphs from root file in directory
	//(make graph my turning off cent binning)

	TGraphErrors* grFracIn = 0;
	TGraphErrors* grFracInPlus = 0;
	TGraphErrors* grFracInMinus = 0;
    TH1F* hRAA[6] = {0};

    if(doScaleByRaa){
        std::cout << "Filling Raa histograms..." << std::endl;
        hRAA[0] = (TH1F*) fRaa->Get("raa_centr0"); //Use 0-5% for 0-5%
        hRAA[1] = (TH1F*) fRaa->Get("raa_centr1"); //Use 5-10% for 5-10%
        hRAA[2] = (TH1F*) fRaa->Get("raa_centr2"); //Use 10-20% for 10-15%
        hRAA[3] = (TH1F*) fRaa->Get("raa_centr2"); //Use 10-20% for 15-20%
        hRAA[4] = (TH1F*) fRaa->Get("raa_centr3"); //Use 20-30% for 20-40%
        hRAA[5] = (TH1F*) fRaa->Get("raa_centr5"); //Use 40-50% for 40-80%
        std::cout << "Done." << std::endl;
    }

    if(doScale){	

		std::cout << "Getting charge-inclusive eta distribution for scaling..." << std::endl;
        grFracIn= (TGraphErrors*) fGr->Get("fractionQCDEta");
	    std::cout << "Done." << std::endl;

	    if(doCharge){
		    std::cout << "Getting charge separated eta distribution for scaling..." << std::endl;
		    grFracInPlus = (TGraphErrors*) fGr->Get("fractionQCDEtaPlus");
		    grFracInMinus = (TGraphErrors*) fGr->Get("fractionQCDEtaMinus");
		    std::cout << "Done." << std::endl;
	    } 
    }

	//data overlay
//	TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuMCDataOverlay.2012.11.21";
	TString fileNameMCWIn = "HISingleMuonWmunuMCDataOverlay.2012.11.21";

  	// --- Jx set ---

/*	TString fileNameMCJ1In = "MonteCarloFiles/QCD/HISingleMuonMcJ1.2012.11.29";

	TString fileNameMCJ2In = "MonteCarloFiles/QCD/HISingleMuonMcJ2.2012.11.29";

	TString fileNameMCJ3In = "MonteCarloFiles/QCD/HISingleMuonMcJ3.2012.11.29";
*/
    ///Nominal Jxmu MC input files for main analysis
	//J1 1 muon-filter 
	TString fileNameMCJ1In = "MonteCarloFiles/QCD/HISingleMuonJ1_1muPYTHIADataOverlay.04.13.2013";
	//J2 1 muon-filter 
	TString fileNameMCJ2In = "MonteCarloFiles/QCD/HISingleMuonJ2_1muPYTHIADataOverlay.04.13.2013";
	//J3 1 muon-filter 
	TString fileNameMCJ3In = "MonteCarloFiles/QCD/HISingleMuonJ3_1muPYTHIADataOverlay.04.13.2013";
    
    TString fileNameMptMtCorrelation = "mptMtCorrelations.07.21.2013";
    TFile* _fMptMtCorrelation = NULL;
    TProfile* _pfxMptMtCorrelation = NULL;
    if(doMptSigmaDown||doMptSigmaUp){

         fileNameMCJ1In = "MonteCarloFiles/QCD/HISingleMuonJ1_1muPYTHIADataOverlay.07.13.2013"; 
         fileNameMCJ2In = "MonteCarloFiles/QCD/HISingleMuonJ2_1muPYTHIADataOverlay.07.13.2013"; 
         fileNameMCJ3In = "MonteCarloFiles/QCD/HISingleMuonJ3_1muPYTHIADataOverlay.07.13.2013"; 
         _fMptMtCorrelation = new TFile(fileNameMptMtCorrelation+".root","read");
         if(_fMptMtCorrelation!=0) std::cout << "Mpt-Mt correlation file: " << fileNameMptMtCorrelation << " opened." << std::endl;
         else exit(0);
    }

    std::cout << "Jxmu files opened: " << std::endl;
    std::cout << "J1mu: " << fileNameMCJ1In << std::endl;
    std::cout << "J2mu: " << fileNameMCJ2In << std::endl;
    std::cout << "J3mu: " << fileNameMCJ3In << std::endl;

/*	TString fileNameMCJ1In = "MonteCarloFiles/QCD/HISingleMuonJ1_1muPYTHIADataOverlay.03.24.2013";
	//J2 1 muon-filter 
	TString fileNameMCJ2In = "MonteCarloFiles/QCD/HISingleMuonJ2_1muPYTHIADataOverlay.03.24.2013";
	//J3 1 muon-filter 
	TString fileNameMCJ3In = "MonteCarloFiles/QCD/HISingleMuonJ3_1muPYTHIADataOverlay.03.24.2013";
*/
        // --- declare cut variables --- //
	  RooRealVar  muonPt("muonPt","p_{T}",0.0,350.0,"GeV");
	  RooRealVar  missPt("missPt","p_{T}^{miss}",0.0,250.0,"GeV");
	  RooRealVar  muonMt("muonMt","m_{T}",0.0,250.0,"GeV");
	  RooRealVar  muonCharge("muonCharge","charge",-2.0,+2.0);
	  RooRealVar  isolationMu("isolationMu","isolationMu",0.0,10.0);
  	  RooRealVar  centrality("centrality","centrality",0.,1.0);
  	  RooRealVar  muonEta("muonEta","muonEta",-3.0,+3.0);
  	  RooRealVar  muonPhi("muonPhi","muonPhi",-3.5,+3.5);
  	  RooRealVar  ZDY("ZDY","ZDY",0.0,1.0);
      RooRealVar  muonQuality("muonQuality","muonQuality",0.0,15.0);
      RooRealVar  muonELoss("muonELoss","muonELoss",-10.0,+10.0);
      RooRealVar  muonScat("muonScat","muonScat",-10.0,+10.0);


	  //W selection cuts
	  //use upper missPt cut to avoid 9999. values
	  TString sCutsSig = "";
      if(doMptSigmaDown||doMptSigmaUp)
          sCutsSig = "muonPt>4.0&&missPt>25.0&&missPt<9000.0&&isolationMu<0.1&&muonMt>40.0&&ZDY==0&&muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0";
      ///PS + Z veto + isolation
      else if(doIsolationOptimization) 
          sCutsSig =
          "muonPt>4.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0&&ZDY==0&&muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0";
      else 
          ///Nominal cuts used for this analysis
          sCutsSig = "muonPt>4.0&&missPt>25.0&&missPt<9000.0&&isolationMu<0.1&&muonMt>40.0&&ZDY==0&&muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0";
          // cross-check
          //sCutsSig = "muonPt>4.0&&missPt>0.0&&missPt<9000.0&&isolationMu<0.1&&muonMt>40.0&&ZDY==0&&muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0";

      ///systematic for ptcone20ID3
     if(doLoosenIsoCut)     
          sCutsSig = "muonPt>4.0&&missPt>25.0&&missPt<9000.0&&isolationMu<0.2&&muonMt>40.0&&ZDY==0&&muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0";

      std::cout << "Cuts: " << sCutsSig << std::endl;
//	  TString sCutsSig = "muonPt>4.0&&missPt>25.0&&missPt<9000.0&&isolationMu<0.2&&muonMt>40.0&&ZDY==0";
//	  TString sCutsSig = "muonPt>4.0&&missPt>25.0&&missPt<9000.0&&muonMt>40.0&&ZDY==0";
	  //systematics
	  RooFormulaVar cutsSig("cutsSig", "cutsSig", sCutsSig,
          RooArgList(muonPt,missPt,isolationMu,muonMt,ZDY,muonQuality,muonELoss,muonScat));
	  RooFormulaVar cutPS("cutPS","cutPS","ZDY==0&&muonPt>4.0&&muonQuality>11&&abs(muonELoss)<0.5&&abs(muonScat)<4.0", 
            RooArgList(ZDY,muonPt,muonQuality,muonELoss,muonScat));

	  RooCategory muonCategory("muonCategory","muonCategory");
	  muonCategory.defineType("Unknown",0);
	  muonCategory.defineType("HeavyFlavour",1);
	  muonCategory.defineType("Tau",2);
	  muonCategory.defineType("LightResonance",3); 
	  muonCategory.defineType("PionCaloDecay",4);
	  muonCategory.defineType("PionEarlyDecay",5);
	  muonCategory.defineType("PionLateDecay",6);
	  muonCategory.defineType("KaonCaloDecay",7);
	  muonCategory.defineType("KaonEarlyDecay",8);
	  muonCategory.defineType("KaonLateDecay",9);
	  muonCategory.defineType("Fake",10);
	  muonCategory.defineType("Z",23);
	  muonCategory.defineType("W",24);

  	  RooCategory chargeCategory("chargeCategory","chargeCategory") ;
  	  chargeCategory.defineType("muMinus",-1) ;
  	  chargeCategory.defineType("muPlus",1) ;

//  	  RooArgSet muonArgSet(muonPt,missPt,isolationMu,muonMt,muonEta,centrality,ZDY,muonCategory,chargeCategory);
//	  muonArgSet.add(muonPhi);
  	  RooArgSet muonArgSet(muonPt,missPt,isolationMu,muonMt,muonEta,centrality,ZDY,muonCategory,chargeCategory);
	  muonArgSet.add(muonPhi);
	  muonArgSet.add(muonQuality);
	  muonArgSet.add(muonELoss);
	  muonArgSet.add(muonScat);
	
	
	  RooArgSet centralityArgSet(centrality);  

	double ptmax = 400.0;
	// --- Set pt and eta bins ---
	std::vector<double> ptBins;
	ptBins.push_back(0.0);
	if (doPt){
		ptBins.push_back(25.0);
		ptBins.push_back(27.0);
		ptBins.push_back(30.0);
		ptBins.push_back(35.0);
		ptBins.push_back(40.0);
		ptBins.push_back(50.0);
		ptBins.push_back(60.0);
		ptBins.push_back(75.0);
		ptBins.push_back(90.0);
		ptBins.push_back(120.0);
	}
	ptBins.push_back(ptmax);
	const int nPtBins = ptBins.size()-1;



	std::vector<double> etaBins;
	etaBins.push_back(0.10);
	if (doEta) {

/*		etaBins.push_back(+0.25);
		etaBins.push_back(+0.50);
		etaBins.push_back(+0.75);
		etaBins.push_back(+1.00);
		etaBins.push_back(+1.25);
		etaBins.push_back(+1.50);
		etaBins.push_back(+1.75);
		etaBins.push_back(+2.00);
		etaBins.push_back(+2.25);
*/
        etaBins.push_back(0.35);
        etaBins.push_back(0.6);
        etaBins.push_back(0.8);
        etaBins.push_back(1.05);
        etaBins.push_back(1.37);
        etaBins.push_back(1.52);
        etaBins.push_back(1.74);
        etaBins.push_back(2.1);
        
	}
//	etaBins.push_back(+2.50);
	etaBins.push_back(+2.40);
	const int nEtaBins = etaBins.size()-1;

	//centrality
	std::vector <double> centralityBins;
	//ncoll
	std::vector <float> ncoll;

	centralityBins.push_back(0.0);
     	if(doCentrality){

        //binning used for optimizing Nsig-Nqcd
//		centralityBins.push_back(0.10);
//		centralityBins.push_back(0.40);
		//ncoll
//		ncoll.push_back(1500.6); //0-10
//		ncoll.push_back(440.6); //10-40
//		ncoll.push_back(77.8); //40-80*/
	centralityBins.push_back(0.05);
	centralityBins.push_back(0.1);
	centralityBins.push_back(0.15);
	centralityBins.push_back(0.2);
	centralityBins.push_back(0.4);

	ncoll.push_back(1683.3); //0-5
	ncoll.push_back(1318.0); //5-10
//	ncoll.push_back(1500.6); //0-10
	ncoll.push_back(1035.4); //10-15
	ncoll.push_back(811.2); //15-20
	ncoll.push_back(440.6); //20-40
	ncoll.push_back(77.8); //40-80
} else  ncoll.push_back(452.0);//0-80

	centralityBins.push_back(0.8);

	const int nCentralityBins = centralityBins.size()-1;

    ///MPT lower thresholds as a function of centrality
    std::vector <double> mptSyst;
    double mptLowCut = 25.0;
    if(doMptSigmaDown&&doCentrality){
        ///-1sigma
        mptSyst.push_back(mptLowCut); //0-5%
        mptSyst.push_back(mptLowCut); //5-10%
        mptSyst.push_back(mptLowCut); //10-15%
        mptSyst.push_back(mptLowCut); //15-20%
        mptSyst.push_back(mptLowCut); //20-40%
        mptSyst.push_back(mptLowCut); //40-80%
        //for(int impt=0; impt<mptLow.size(); ++impt) std::cout << mptLow[impt] << std::endl;
    }
    else if(doMptSigmaDown) mptSyst.push_back(mptLowCut); //0-80%
    else if (doMptSigmaUp&&doCentrality){
        ///+1 sigma
        mptSyst.push_back(mptLowCut); //0-5%
        mptSyst.push_back(mptLowCut); //5-10%
        mptSyst.push_back(mptLowCut); //10-15%
        mptSyst.push_back(mptLowCut); //15-20%
        mptSyst.push_back(mptLowCut); //20-40%
        mptSyst.push_back(mptLowCut); //40-80%
    }
    else if (doMptSigmaUp) mptSyst.push_back(mptLowCut); //0-80%

	double arrQCDFrac ;
	double arrQCDFracPlus ;
	double arrQCDFracMinus ;

	double arrSumQCD ;
	double arrSumQCDPlus ;
	double arrSumQCDMinus ;

  	// --- Fill data sets ---
        RooDataSet* dataSet = fillHIMuonDataSet(baseString,fileNameDataIn+".root",muonArgSet, cutValue); dataSet->Print(); 	
	//pre-selected dataset
	RooDataSet* psDataSet = (RooDataSet*)dataSet->reduce(Cut(cutPS)); std::cout << "Number of pre-selected muons in data = "<< psDataSet->numEntries() << std::endl;
	//apply W selection cuts
	dataSet = (RooDataSet*)dataSet->reduce(Cut(cutsSig)); dataSet->Print(); std::cout << "After W selection: " << dataSet->numEntries() << std::endl;

        RooDataSet* mcJ1Set = fillHIMuonDataSet(baseString,fileNameMCJ1In+".root",muonArgSet, cutValue, true); mcJ1Set->Print();
	RooDataSet* psJ1Set = (RooDataSet*)mcJ1Set->reduce(Cut(cutPS));
	mcJ1Set = (RooDataSet*)mcJ1Set->reduce(Cut(cutsSig)); mcJ1Set->Print();
        RooDataSet* mcJ1Events = fillHIEventDataSet(baseString,fileNameMCJ1In+".root",centralityArgSet );

        RooDataSet* mcJ2Set = fillHIMuonDataSet(baseString,fileNameMCJ2In+".root",muonArgSet, cutValue, true); mcJ2Set->Print();
	RooDataSet* psJ2Set = (RooDataSet*)mcJ2Set->reduce(Cut(cutPS));
	mcJ2Set = (RooDataSet*)mcJ2Set->reduce(Cut(cutsSig)); mcJ2Set->Print();
        RooDataSet* mcJ2Events = fillHIEventDataSet(baseString,fileNameMCJ2In+".root",centralityArgSet );

        RooDataSet* mcJ3Set = fillHIMuonDataSet(baseString,fileNameMCJ3In+".root",muonArgSet, cutValue, true); mcJ3Set->Print();
	RooDataSet* psJ3Set = (RooDataSet*)mcJ3Set->reduce(Cut(cutPS));
	mcJ3Set = (RooDataSet*)mcJ3Set->reduce(Cut(cutsSig)); mcJ3Set->Print();
        RooDataSet* mcJ3Events = fillHIEventDataSet(baseString,fileNameMCJ3In+".root",centralityArgSet );

	// --- Subdivide in bins ---
	RooDataSet* dataSubSet[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* dataSubSetPS[nPtBins][nEtaBins][nCentralityBins];

	RooDataSet* mcJ1SubSet[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* psJ1SubSet[nPtBins][nEtaBins][nCentralityBins];

	RooDataSet* mcJ2SubSet[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* psJ2SubSet[nPtBins][nEtaBins][nCentralityBins];

	RooDataSet* mcJ3SubSet[nPtBins][nEtaBins][nCentralityBins];
	RooDataSet* psJ3SubSet[nPtBins][nEtaBins][nCentralityBins];

	RooDataSet* mcJ1EventSubSet[nCentralityBins];
	RooDataSet* mcJ2EventSubSet[nCentralityBins];
	RooDataSet* mcJ3EventSubSet[nCentralityBins];

	for ( int i = 0; i < nPtBins; i++ ) {
	  for ( int j = 0; j < nEtaBins; j++ ) {
	    for ( int k = 0; k < nCentralityBins; k++ ){

            ///Get Mpt/Mt correlation
            if(doMptSigmaDown||doMptSigmaUp) {
                TString sMptMtCorrelation = "h2DMptMtCorrelationCent"; sMptMtCorrelation+=k; sMptMtCorrelation+="_pfx";
                std::cout << "Opening " << sMptMtCorrelation << " TProfilex histo." << std::endl;
                _pfxMptMtCorrelation = (TProfile*)_fMptMtCorrelation->Get(sMptMtCorrelation);
                if(_pfxMptMtCorrelation!=0) std::cout << sMptMtCorrelation << " opened." << std::endl;
                else exit(0);
            }

        	dataSubSet[i][j][k] = selectPtEtaCentrality(dataSet, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 
        	dataSubSetPS[i][j][k] = selectPtEtaCentrality(psDataSet, ptBins[0], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 

            mcJ1SubSet[i][j][k] = selectPtEtaCentrality(mcJ1Set, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 
            //if((doMptSigmaDown||doMptSigmaUp))  mcJ1SubSet[i][j][k] = selectMptMt(mcJ1SubSet[i][j][k],mptSyst[k],_pfxMptMtCorrelation); 
            psJ1SubSet[i][j][k] = selectPtEtaCentrality(psJ1Set, ptBins[0], ptBins[nPtBins], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 

            mcJ2SubSet[i][j][k] = selectPtEtaCentrality(mcJ2Set, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 
            //if((doMptSigmaDown||doMptSigmaUp))  mcJ2SubSet[i][j][k] = selectMptMt(mcJ2SubSet[i][j][k],mptSyst[k],_pfxMptMtCorrelation); 
            psJ2SubSet[i][j][k] = selectPtEtaCentrality(psJ2Set, ptBins[0], ptBins[nPtBins], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 

            mcJ3SubSet[i][j][k] = selectPtEtaCentrality(mcJ3Set, ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 
            //if((doMptSigmaDown||doMptSigmaUp))  mcJ3SubSet[i][j][k] = selectMptMt(mcJ3SubSet[i][j][k],mptSyst[k],_pfxMptMtCorrelation); 
            psJ3SubSet[i][j][k] = selectPtEtaCentrality(psJ3Set, ptBins[0], ptBins[nPtBins], etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1],true); 

            //number of Jx events in centrality bin k
            mcJ1EventSubSet[k] = selectCentrality(mcJ1Events,centralityBins[k], centralityBins[k+1]);
            mcJ2EventSubSet[k] = selectCentrality(mcJ2Events,centralityBins[k], centralityBins[k+1]);
            mcJ3EventSubSet[k] = selectCentrality(mcJ3Events,centralityBins[k], centralityBins[k+1]);
	    }
	  }
	}


	TGraphErrors* grFrac = new TGraphErrors(nCentralityBins);
	TGraphErrors* grFracEta = new TGraphErrors(nEtaBins);
	TGraphErrors* grFracPt = new TGraphErrors(nPtBins);
	TGraphErrors* grFracPlus = new TGraphErrors(nCentralityBins);
	TGraphErrors* grFracEtaPlus = new TGraphErrors(nEtaBins);
	TGraphErrors* grFracMinus = new TGraphErrors(nCentralityBins);
	TGraphErrors* grFracEtaMinus = new TGraphErrors(nEtaBins);

	//plot difference in charge for systematics
	TGraph* grDiffEta = new TGraphErrors(nEtaBins);
	TGraph* grDiffNcoll  = new TGraphErrors(nCentralityBins);
	//plot residual of 1st order polynomial with Npart distr
	TGraph* grResid = new TGraphErrors(nCentralityBins);

    ///Analysis values
    double bkgFracInt,bkgFracIntPlus,bkgFracIntMinus;

    double ncollMean = 452.0; //0-80
	TGraph* grFracMean = new TGraphErrors(1);
	TGraph* grFracMeanPlus = new TGraphErrors(1);
	TGraph* grFracMeanMinus = new TGraphErrors(1);

    //use for ptcone systematics
    //ptcone30ID3<0.1
    if(doIncreaseConeSize){
	  bkgFracInt = 0.030;
	  bkgFracIntPlus = 0.028;
	  bkgFracIntMinus = 0.031;
    }
    //ptcone20ID3<0.2
    else if(doLoosenIsoCut){
	  bkgFracInt = 0.046;
	  bkgFracIntPlus = 0.045;
	  bkgFracIntMinus = 0.047;
    }
/*	double bkgFracInt = 0.0456907;
	double bkgFracIntPlus = 0.0448675;
	double bkgFracIntMinus = 0.0465362;
*/

    ///mpt +1 sigma
    else if(doMptSigmaUp){
	 bkgFracInt = 0.0320611;
	 bkgFracIntPlus = 0.0324201;
	 bkgFracIntMinus = 0.0317088;

/*	 bkgFracInt = 0.0354052;
	 bkgFracIntPlus = 0.0342147;
	 bkgFracIntMinus = 0.0366457;
     */
    }
    ///mpt -1 sigma
    else if(doMptSigmaDown){
	 bkgFracInt = 0.037156;
	 bkgFracIntPlus = 0.0359103;
	 bkgFracIntMinus = 0.0384676;
    }
    else if(doIsolationOptimization){
	 bkgFracInt = 0.068708;
	 bkgFracIntPlus = 0.0667931;
	 bkgFracIntMinus = 0.0706622;
    }
    else if(doScaleByRaa){
     ///from pol0 fit to Ncoll
	 bkgFracInt = 0.0507917;
	 bkgFracIntPlus = 0.0507917;
	 bkgFracIntMinus = 0.0507917;

    }
    else{
	 bkgFracInt = 0.0324447;
	 bkgFracIntPlus = 0.0323822;
	 bkgFracIntMinus = 0.0325233;
     
     //cross-check w/o mpt cut
     /*
     bkgFracInt =0.0475482;
     bkgFracIntPlus=0.0467866;
     bkgFracIntMinus =0.0483525;
     */
    }

    std::cout << "Pegged point for inclusive charge = " << bkgFracInt << std::endl;

    grFracMean->SetPoint(0,ncollMean,bkgFracInt);
    grFracMeanPlus->SetPoint(0,ncollMean,bkgFracIntPlus);
    grFracMeanMinus->SetPoint(0,ncollMean,bkgFracIntMinus);

	for ( int i = 0; i < nPtBins; i++ ) {
	  for ( int j = 0; j < nEtaBins; j++ ) {
	    for ( int k = 0; k < nCentralityBins; k++ ){
		
		std::cout << " plotting "<<k<<":"<<j<< std::endl;
        TString sCentLow = "",sCentUp="",sEtaLow="",sEtaUp=""; 
        sCentLow += (100.*centralityBins[k]); //sCentLow.Remove(3);
        sCentUp += (100.*centralityBins[k+1]); //sCentUp.Remove(3);

		sEtaLow += (etaBins[j]);
		sEtaUp += (etaBins[j+1]);

		TString sSel = sCentLow; sSel += "-"; sSel += sCentUp; sSel += "%"; 
        TString sSel2 = sEtaLow; sSel2 += "-"; sSel2 += sEtaUp;

	        //102 = mu+ , 103 = mu- , 104 = mu^{pm}
		std::cout << "mu^{#pm}" << std::endl;
		if(doCentrality) 
			plot(outFile,muonPt,dataSubSet[i][j][k],mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k], 
				dataSubSetPS[i][j][k],psJ1SubSet[i][j][k],psJ2SubSet[i][j][k],psJ3SubSet[i][j][k], 
				mcJ1EventSubSet[k],mcJ2EventSubSet[k],mcJ3EventSubSet[k], hRAA[k],
				grFracPt,grFrac, grFracEta, grFracIn, arrQCDFrac, arrSumQCD,  
				ptBins[i],ptBins[i+1],etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1], ncoll[k],104 ,i
                , k,  j,  sSel, "0.1<|#eta|<2.4", bkgFracInt, 
                doScaleByRaa,doScale);
		else 
			plot(outFile,muonPt,dataSubSet[i][j][k],mcJ1SubSet[i][j][k],mcJ2SubSet[i][j][k],mcJ3SubSet[i][j][k], 
				dataSubSetPS[i][j][k],psJ1SubSet[i][j][k],psJ2SubSet[i][j][k],psJ3SubSet[i][j][k], 
				mcJ1Events,mcJ2Events,mcJ3Events, hRAA[k],
				grFracPt,grFrac, grFracEta, grFracIn, arrQCDFrac, arrSumQCD,  
				ptBins[i],ptBins[i+1],etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1], ncoll[k],104 ,i , k,  j,  sSel, sSel2, bkgFracInt, 
                doScaleByRaa,doScale);

		if(doCharge){ 

	           std::cout << "Creating charged datasets." << std::endl;
		   RooDataSet* dataSetPlus  = (RooDataSet*) dataSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   RooDataSet* dataSetMinus = (RooDataSet*) dataSubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

		   RooDataSet* dataSetPlusPS  = (RooDataSet*) dataSubSetPS[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   RooDataSet* dataSetMinusPS = (RooDataSet*) dataSubSetPS[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

		   RooDataSet* mcJ1SetPlusPS  = (RooDataSet*) psJ1SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   RooDataSet* mcJ1SetMinusPS = (RooDataSet*) psJ1SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
		   RooDataSet* mcJ2SetPlusPS  = (RooDataSet*) psJ2SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   RooDataSet* mcJ2SetMinusPS = (RooDataSet*) psJ2SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
		   RooDataSet* mcJ3SetPlusPS  = (RooDataSet*) psJ3SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   RooDataSet* mcJ3SetMinusPS = (RooDataSet*) psJ3SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));

		   RooDataSet* mcJ1SetPlus  = (RooDataSet*) mcJ1SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   RooDataSet* mcJ1SetMinus = (RooDataSet*) mcJ1SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
		   RooDataSet* mcJ2SetPlus  = (RooDataSet*) mcJ2SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   RooDataSet* mcJ2SetMinus = (RooDataSet*) mcJ2SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));
		   RooDataSet* mcJ3SetPlus  = (RooDataSet*) mcJ3SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muPlus"));
		   RooDataSet* mcJ3SetMinus = (RooDataSet*) mcJ3SubSet[i][j][k]->reduce(Cut("chargeCategory==chargeCategory::muMinus"));


			TString sSelPlus; TString sSelMinus;
			sSelPlus = sSel + ",#mu^{+}";
			sSelMinus = sSel + ",#mu^{-}";

			if(doCentrality){
		    		std::cout << "mu^{+}" << std::endl;
				plot(outFile,muonPt,dataSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus, 
					dataSetPlusPS,mcJ1SetPlusPS,mcJ2SetPlusPS,mcJ3SetPlusPS,
					mcJ1EventSubSet[k],mcJ2EventSubSet[k],mcJ3EventSubSet[k], hRAA[k],
			  		grFracPt,grFracPlus, grFracEtaPlus, grFracInPlus,
			  		arrQCDFracPlus, arrSumQCDPlus, ptBins[i],ptBins[i+1],etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1], ncoll[k], 
					102,i, k,  j,  sSelPlus, sSel2, bkgFracIntPlus ,doScaleByRaa,doScale);
				std::cout << "mu^{-}" << std::endl;
				plot(outFile,muonPt,dataSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus, 
					dataSetMinusPS,mcJ1SetMinusPS,mcJ2SetMinusPS,mcJ3SetMinusPS,
					mcJ1EventSubSet[k],mcJ2EventSubSet[k],mcJ3EventSubSet[k], hRAA[k],
					grFracPt,grFracMinus, grFracEtaMinus, grFracInMinus,
			  		arrQCDFracMinus, arrSumQCDMinus, ptBins[i],ptBins[i+1],etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1], ncoll[k], 
					103,i, k,  j,  sSelMinus, sSel2, bkgFracIntMinus, doScaleByRaa,doScale);

			  std::cout << "<Ncoll>" << std::endl;
			  plotChargeDiffNcoll(grDiffNcoll,grFracPlus,grFracMinus,j,k,ncoll[k]);
			} //centrality

			else /*if(doEta)*/{ 

		    		std::cout << "mu^{+}" << std::endl;
				plot(outFile,muonPt,dataSetPlus,mcJ1SetPlus,mcJ2SetPlus,mcJ3SetPlus, 
					dataSetPlusPS,mcJ1SetPlusPS,mcJ2SetPlusPS,mcJ3SetPlusPS,
					mcJ1Events,mcJ2Events,mcJ3Events, hRAA[k],
			  		grFracPt,grFracPlus, grFracEtaPlus, grFracInPlus,
			  		arrQCDFracPlus, arrSumQCDPlus, ptBins[i],ptBins[i+1],etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1], ncoll[k], 
					102,i, k,  j,  sSelPlus, sSel2, bkgFracIntPlus ,doScaleByRaa,doScale);

		    		std::cout << "mu^{-}" << std::endl;
				plot(outFile,muonPt,dataSetMinus,mcJ1SetMinus,mcJ2SetMinus,mcJ3SetMinus, 
					dataSetMinusPS,mcJ1SetMinusPS,mcJ2SetMinusPS,mcJ3SetMinusPS,
					mcJ1Events,mcJ2Events,mcJ3Events, hRAA[k],
					grFracPt,grFracMinus, grFracEtaMinus, grFracInMinus,
			  		arrQCDFracMinus, arrSumQCDMinus, ptBins[i],ptBins[i+1],etaBins[j], etaBins[j+1], centralityBins[k], centralityBins[k+1], ncoll[k], 
					103,i, k,  j,  sSelMinus, sSel2, bkgFracIntMinus, doScaleByRaa,doScale);

		    	  std::cout << "|mu^{-}-mu^{+}|" << std::endl;
			  std::cout << "|eta|" << std::endl;
			  plotChargeDiffEta(grDiffEta,grFracEtaPlus,grFracEtaMinus,j,k,etaBins[j], etaBins[j+1]);
			} //eta

		} //charge

	  } //k

      }//j
     } //i

//	TF1 f0a = TF1("f0a","[0]",0,1700);
//	f0a.SetLineStyle(kDashed);
//	grFrac->Fit("f0a");


	if((doCentrality||doEta)) {
		TString sCharge = "#mu^{pm}";
		plotFraction(grFrac,grFracEta,sCharge);
		fitFraction(grFrac,sCharge+",fit");
		//systematics
		if(doCentrality) plotResidualNcoll(grFrac,grResid); 
	}
	if((doCentrality||doEta)&&doCharge) {
		TString sChargePlus = "#mu^{+}";
		TString sChargeMinus = "#mu^{-}";
		plotFraction(grFracPlus,grFracEtaPlus,sChargePlus);
		plotFraction(grFracMinus,grFracEtaMinus,sChargeMinus);
		//systematics
		fitFraction(grFracPlus,sChargePlus+",fit");
		fitFraction(grFracMinus,sChargePlus+",fit");
	}

	if(!doScale) {
		Write(outFile,grFrac, "fractionQCDCent");
		Write(outFile,grFracEta, "fractionQCDEta");
		Write(outFile,grFracPt, "fractionQCDPt");
        Write(outFile,grFracMean,"meanFracQCDCent");
		//systematics
		Write(outFile,grResid,"fitQCDResidual");

		if(doCharge) { 

            ///graphs of wtd average of bkg fraction for mu+ and mu-
	        TGraph* grWtdAvgFracPlus = new TGraphErrors(1);
	        TGraph* grWtdAvgFracMinus = new TGraphErrors(1);

            if(doEta){
              grWtdAvgFracPlus->SetPoint(0,(etaBins[0]+etaBins[nEtaBins])/2.0,arrQCDFracPlus/arrSumQCDPlus);
              grWtdAvgFracMinus->SetPoint(0,(etaBins[0]+etaBins[nEtaBins])/2.0,arrQCDFracMinus/arrSumQCDMinus);

              ///Fit to pol0; difference in charges used as systematic
              fitToPol0(grFracEtaPlus,etaBins[0],etaBins[nEtaBins]);
              fitToPol0(grFracEtaMinus,etaBins[0],etaBins[nEtaBins]);

			  Write(outFile,grFracEtaPlus, "fractionQCDEtaPlus");
			  Write(outFile,grFracEtaMinus, "fractionQCDEtaMinus");
			  Write(outFile,grDiffEta,"diffFracQCDEta"); 
            }
            if(doCentrality){
              grWtdAvgFracPlus->SetPoint(0,ncollMean,arrQCDFracPlus/arrSumQCDPlus);
              grWtdAvgFracMinus->SetPoint(0,ncollMean,arrQCDFracMinus/arrSumQCDMinus);

              fitToPol0(grFracPlus,ncoll[0],ncoll[nCentralityBins]);
              fitToPol0(grFracMinus,ncoll[0],ncoll[nCentralityBins]);

			  Write(outFile,grFracPlus, "fractionQCDCentPlus");
			  Write(outFile,grFracMinus, "fractionQCDCentMinus");
		      Write(outFile,grDiffNcoll,"diffFracQCDNcoll");
            }

            Write(outFile,grFracMeanPlus,"meanFracQCDCentPlus");
            Write(outFile,grFracMeanMinus,"meanFracQCDCentMinus");
            Write(outFile,grWtdAvgFracPlus,"wtdAvgFracQCDPlus");
            Write(outFile,grWtdAvgFracMinus,"wtdAvgFracQCDMinus");

		}
        }

    /*if(doScale) {
        char* nameFrom = "fractionQCDEtaScaled.root";
        char* nameTo = "systematics/fractionQCDEtaScaled_MptSmearSystematic.`date +%m.%d.%Y`.root";
        TUnixSystem::Rename(nameFrom,nameTo);
    }*/
}
int main(){
    bkgQCDPlotter();
    std::cout << "Macro completed. " << std::endl;
}
