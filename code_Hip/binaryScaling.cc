// Plot yields as a function of <Npart> for W-->mu (+,-,and +-)
// xsec from AMI [nb]
// W- -->mu- + nu~
//  pn(np): 3.4410E+00
//  pp:     2.7368E+00
// W+ -->mu+ + nu  
//  pn(np): 3.4630E+00
//  pp:     4.1936E+00 
// xsec from POWHEG [nb]
// W- -->mu- + nu~
//  pn(np): 3.595E+00
//  pp:     2.901E+00
// W+ -->mu+ + nu  
//  pn(np): 3.588E+00
//  pp:     4.297E+00 

#include "TH1.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "AtlasStyle.C"

#include <iostream>
#include <string>
#include <memory>

float getModelXSec(TString sFileName, float sigma, double centLow, double centUpp, double etaLow, double etaUpp, std::string charge){
    bool isPlus = false;
    bool isMinus = false;
    if(strcmp(charge.c_str(),"+")==0) isPlus=true;
    else if(strcmp(charge.c_str(),"-")==0) isMinus=true;

    std::unique_ptr<TFile> _fFcalScF(new TFile("fcalScaleFactors.09.14.2014.root","read") );
    TF1* _funcScF = (TF1*)((TH1F*)_fFcalScF->Get("_hRatio"))->GetFunction("_fit");

    std::unique_ptr<TChain> _tGen(new TChain("treeGen","treeGen"));
    _tGen->Add(sFileName);
    std::unique_ptr<TChain> _tRec(new TChain("treeW","treeW"));
    _tRec->Add(sFileName);
    std::cout << "Number of entries in generator tree: " << _tGen->GetEntries() << std::endl;
    std::cout << "Number of entries in reconstructed tree: " << _tRec->GetEntries() << std::endl;

    float recPt, recEta, recCentrality, recMET, recMt, recCharge,recPtcone30;
    int recPDG, recMother, recMatched; 
    float genPt, genEta, genCentrality, genNeutrinoPt, genMt, genCharge;
    float fcalGen,fcalRec;
    int genPDG, genNeutrinoPDG, genMother;
    bool isGenPeriodA, isGenPeriodB;
    bool isRecPeriodA, isRecPeriodB;
    bool isZCandidate, EF_mu8_matched;

    _tRec->SetBranchAddress("pt",&recPt);
    _tRec->SetBranchAddress("eta",&recEta);
    _tRec->SetBranchAddress("centrality",&recCentrality);
    _tRec->SetBranchAddress("met",&recMET);
    _tRec->SetBranchAddress("mt",&recMt);
    _tRec->SetBranchAddress("charge",&recCharge);
    _tRec->SetBranchAddress("rec_truth_pdg",&recPDG);
    _tRec->SetBranchAddress("rec_truth_mother_pdg",&recMother);
    _tRec->SetBranchAddress("ptcone30",&recPtcone30);
    _tRec->SetBranchAddress("rec_truth_matched",&recMatched);
    _tRec->SetBranchAddress("rec_truth_mother_pdg",&recMother);
    _tRec->SetBranchAddress("rec_truth_pdg",&recPDG);
    _tRec->SetBranchAddress("isPeriodA",&isRecPeriodA);
    _tRec->SetBranchAddress("isPeriodB",&isRecPeriodB);
    _tRec->SetBranchAddress("isZCandidateType1",&isZCandidate);
    _tRec->SetBranchAddress("EF_mu8_matched",&EF_mu8_matched);
    _tRec->SetBranchAddress("fcal",&fcalRec);

    _tGen->SetBranchAddress("muon_gen_pt",&genPt);
    _tGen->SetBranchAddress("muon_gen_eta",&genEta);
    _tGen->SetBranchAddress("centrality",&genCentrality);
    _tGen->SetBranchAddress("neutrino_gen_pt",&genNeutrinoPt);
    _tGen->SetBranchAddress("mt_gen_neutrino_muon",&genMt);
    _tGen->SetBranchAddress("muon_gen_charge",&genCharge);
    _tGen->SetBranchAddress("muon_pdg",&genPDG);
    _tGen->SetBranchAddress("neutrino_pdg",&genNeutrinoPDG);
    _tGen->SetBranchAddress("muon_mother",&genMother);
    _tGen->SetBranchAddress("isPeriodA",&isGenPeriodA);
    _tGen->SetBranchAddress("isPeriodB",&isGenPeriodB);
    _tGen->SetBranchAddress("fcal",&fcalGen);


    float metCut = 25.0;
    float mtCut  = 40.0;
    float ptCut  = 25.0;
    // Total number of generated W-->mu 
    // over all space
    int nGen = 0;
    int nFid = 0;
    TH1F* _hGen = new TH1F("_hGen","_hGen",400,0.0,600.0);
    TH1F* _hFid = new TH1F("_hFid","_hFid",400,0.0,600.0);
    for(int igen=0; igen<_tGen->GetEntries(); ++igen){
        _tGen->GetEntry(igen);
        //double sf = _funcScF->Eval(fcalGen);
        double sf = 1.0;
        if(fabs(genPDG) == 13 && fabs(genMother) == 24){ 
            if(genCharge>0.0&&isPlus==true) _hGen->Fill(genPt,sf);
            else if (genCharge<0.0&&isMinus==true) _hGen->Fill(genPt,sf);
            //_hGen->Fill(genPt,sf);
            if(     genPt > ptCut  &&
                    genPt < 600.0  &&
                    genEta>etaLow  && genEta<etaUpp  &&
                    genCentrality>centLow   && genCentrality<centUpp   &&
                    genNeutrinoPt > metCut  &&
                    genMt > mtCut           
                    && ((genCharge<0.0&&isMinus==true)||(genCharge>0.0&&isPlus==true)) 
                    //safety
                    //fabs(genPDG) == 13      &&
                    //fabs(genMother) == 24
          ) _hFid->Fill(genPt,sf);
        }
    }//igen
    // Number of W-->mu events reconstructed
    // in analysis space
    /*for(int irec=0; irec<_tRec->GetEntries(); ++irec){
        _tRec->GetEntry(irec);
        float isoRatio = recPtcone30/recPt;
        if(
            recPt>ptCut && recPt<400.0
            && recEta>etaLow && recEta<etaUpp
            && recCentrality>centLow  && recCentrality<centUpp 
            && isoRatio <= 0.2
            && recMET > metCut && recMET < 400.0 
            && recMt > mtCut && recMt < 400.0
            && recMatched==1
            //safety
            && fabs(recPDG) == 13    
            && fabs(recMother) == 24 
            && isZCandidate == false     
            && EF_mu8_matched == 1
            && ( (recCharge>0.0&&isPlus==true) || (recCharge<0.0&&isMinus==true) )

          ){ 
            //double sf = _funcScF->Eval(fcalRec);
            double sf = 1.0;
            _hFid->Fill(recPt,sf);
        }
        
    }//irec
   */ 
    //const float runLumi = 27.8; // nb^{-1}
    const float sigma_NN = 70e6; // nb
    float Aw = _hFid->Integral()/_hGen->Integral();
    if(Aw>1.0) {std::cout << "ERROR: Aw > 1.0" << std::endl; exit(0);}
    else std::cout << "Aw: " << Aw << std::endl;
    //float Aw = (float)nFid/nGen;
    delete _hGen;
    delete _hFid;
    return Aw*sigma/sigma_NN;
}

void run(const std::string& period, std::string omega, TFile* outFile, bool doCentralityBiasCorr=false){
	
	bool doPeriodA = false;
	bool doPeriodB = false;
    bool doOmega055 = false;
    bool doOmega101 = false;
    TString cName = "c";

    if(strcmp(omega.c_str(),"0.55")==0) { 
        doOmega055=true;
        cName+="_Omega055";
    }
    else if(strcmp(omega.c_str(),"1.01")==0) {
        doOmega101=true;
        cName+="_Omega101";
    }
    else cName+="_Omega0";

	if(period.compare("A") == 0) {
		doPeriodA = true;
		std::cout << "*****************************" << std::endl;
		std::cout << "**Running on Period A data.**" << std::endl;
		std::cout << "*****************************" << std::endl;
	}
	else if(period.compare("B") == 0) {
		doPeriodB = true;
		std::cout << "*****************************" << std::endl;
		std::cout << "**Running on Period B data.**" << std::endl;
		std::cout << "*****************************" << std::endl;
	}
	else {
		std::cout << "*******************************" << std::endl;
		std::cout << "**Running on Period A+B data.**" << std::endl;
		std::cout << "*******************************" << std::endl;
	}

    TString sDataFileIn = "ResultFiles/WAnalysis_fitResultCentChrgEta_";
	if(doPeriodA) sDataFileIn+="PeriodA.09.25.2013.root";
	else if(doPeriodB) sDataFileIn+="PeriodB.09.25.2013.root";
	else {
        // Uncorrected
        //sDataFileIn+="PeriodAB_Raw.07.24.2014.root";
        // Cw correction applied
        //sDataFileIn+="PeriodAB_Corrected.09.12.2014.root";
        // periods in MC reversed
        //sDataFileIn+="ZBinning_PeriodAB_Corrected.09.15.2014.root";
        // no fcal re-weighting (negligible)
        //sDataFileIn+="PeriodAB_Corrected.10.08.2014.root";
        //sDataFileIn+="PeriodAB_CorrectedFCalReWtd.09.24.2014.root";
        sDataFileIn+="PeriodAB_Corrected_EWBkgSub.10.08.2014.root";
    }

	TFile* dataFileIn = new TFile(sDataFileIn,"read");

	TGraphAsymmErrors* grCentrality = (TGraphAsymmErrors*)dataFileIn->Get("WPtFit_centrality");
	const int nCentralityBins = grCentrality->GetN();
	std::cout << "No. centrality bins: " << nCentralityBins << std::endl;

	TGraphAsymmErrors* grEta = (TGraphAsymmErrors*)dataFileIn->Get("WPtFit_eta");
	const int nEtaBins = grEta->GetN();
	std::cout << "No. eta bins: " << nEtaBins << std::endl;
	
	TGraphAsymmErrors* grPt = (TGraphAsymmErrors*)dataFileIn->Get("WPtFit_pt");
	const int nPtBins = grPt->GetN();
	std::cout << "No. pt bins: " << nPtBins << std::endl;

	float nWPlus[nCentralityBins],nWPlusStatErr[nCentralityBins],nWMinus[nCentralityBins],nWMinusStatErr[nCentralityBins];

    TGraphAsymmErrors* grDataPlus = new TGraphAsymmErrors(nCentralityBins);
    TGraphAsymmErrors* grDataMinus = new TGraphAsymmErrors(nCentralityBins);
    TGraphAsymmErrors* grDataRatio = new TGraphAsymmErrors(nCentralityBins);
    TGraphAsymmErrors* grDataInclusive = new TGraphAsymmErrors(nCentralityBins);

    TString sNameDataPlus,sNameDataMinus,sNameDataInclusive, sNameDataRatio;
    if(doPeriodA){
        sNameDataPlus = "muPlusNpartDistroPeriodA";
        sNameDataMinus = "muMinusNpartDistroPeriodA";
        sNameDataRatio = "muChRatioNpartDistroPeriodA";
        sNameDataInclusive = "muInclusiveNpartDistroPeriodA";
    }
    else if(doPeriodB){
        sNameDataPlus = "muPlusNpartDistroPeriodB";
        sNameDataMinus = "muMinusNpartDistroPeriodB";
        sNameDataRatio = "muChRatioNpartDistroPeriodB";
        sNameDataInclusive = "muInclusiveNpartDistroPeriodB";
    }
    else{
        sNameDataPlus = "muPlusNpartDistroPeriodAB";
        sNameDataMinus = "muMinusNpartDistroPeriodAB";
        sNameDataRatio = "muChRatioNpartDistroPeriodAB";
        sNameDataInclusive = "muInclusiveNpartDistroPeriodAB";
    }

    std::vector<double> npart;
    std::vector<double> npart_relerr_upp;
    std::vector<double> npart_relerr_low;
    std::vector<double> ncoll;
    std::vector<double> centBiasCorrection;
    std::vector<double> centBiasCorrection_relerr;
    double meanNcoll;
    if(doOmega055){
        std::cout << "OMEGA = 0.55" << std::endl;
        sNameDataPlus+="Omega055";
        sNameDataMinus+="Omega055";
        sNameDataRatio+="Omega055";
        sNameDataInclusive+="Omega055";

        npart.push_back(20.45); //0-5%
        npart.push_back(16.49); //5-10%
        npart.push_back(13.77); //10-20%
        npart.push_back(11.23); //20-30%
        npart.push_back(9.22); //30-40%
        npart.push_back(6.68); //40-60%
        npart.push_back(3.59); //60-90%
        /*npart.push_back(10.23); //20-40%
        npart.push_back(4.82); //40-90%
        */

        npart_relerr_upp.push_back(6.31); //0-5%
        npart_relerr_upp.push_back(6.06); //5-10%
        npart_relerr_upp.push_back(5.74); //10-20%
        npart_relerr_upp.push_back(5.51); //20-30%
        npart_relerr_upp.push_back(5.40); //30-40%
        npart_relerr_upp.push_back(5.83); //40-60%
        npart_relerr_upp.push_back(5.35); //60-90%
        /*npart_relerr_upp.push_back(5.46); //20-40%
        npart_relerr_upp.push_back(5.61); //40-90%
        */

        npart_relerr_low.push_back(7.09); //0-5%
        npart_relerr_low.push_back(6.06); //5-10%
        npart_relerr_low.push_back(5.91); //10-20%
        npart_relerr_low.push_back(5.92); //20-30%
        npart_relerr_low.push_back(5.81); //30-40%
        npart_relerr_low.push_back(5.66); //40-60%
        npart_relerr_low.push_back(4.85); //60-90%
        /*npart_relerr_low.push_back(5.87); //20-40%
        npart_relerr_low.push_back(5.30); //40-90%
        */

        meanNcoll= 8.535-1.0;  // 0-90%

        centBiasCorrection.push_back(1.11256); //0-5%
        centBiasCorrection.push_back(1.06552); //5-10%
        centBiasCorrection.push_back(1.03368); //10-20%
        centBiasCorrection.push_back(1.00439); //20-30%
        centBiasCorrection.push_back(0.982925); //30-40%
        centBiasCorrection.push_back(0.956645); //40-60%
        centBiasCorrection.push_back(0.869336); //60-90%
        // rel error
        centBiasCorrection_relerr.push_back(2.0); //0-5%
        centBiasCorrection_relerr.push_back(1.0); //5-10%
        centBiasCorrection_relerr.push_back(1.0); //10-20%
        centBiasCorrection_relerr.push_back(0.0); //20-30%
        centBiasCorrection_relerr.push_back(0.0); //30-40%
        centBiasCorrection_relerr.push_back(1.0); //40-60%
        centBiasCorrection_relerr.push_back(3.0); //60-90%
    }
    else if(doOmega101){
        std::cout << "OMEGA = 1.01" << std::endl;

        sNameDataPlus+="Omega101";
        sNameDataMinus+="Omega101";
        sNameDataRatio+="Omega101";
        sNameDataInclusive+="Omega101";

        npart.push_back(22.61); //0-5%
        npart.push_back(17.46); //5-10%
        npart.push_back(14.11); //10-20%
        npart.push_back(11.17); //20-30%
        npart.push_back(8.97); //30-40%
        npart.push_back(6.37); //40-60%
        npart.push_back(3.43); //60-90%

        npart_relerr_upp.push_back(6.55); //0-5%
        npart_relerr_upp.push_back(6.47); //5-10%
        npart_relerr_upp.push_back(6.11); //10-20%
        npart_relerr_upp.push_back(6.09); //20-30%
        npart_relerr_upp.push_back(6.65); //30-40%
        npart_relerr_upp.push_back(4.08); //40-60%
        npart_relerr_upp.push_back(7.57); //60-90%

        npart_relerr_low.push_back(10.88); //0-5%
        npart_relerr_low.push_back(6.19); //5-10%
        npart_relerr_low.push_back(5.61); //10-20%
        npart_relerr_low.push_back(5.60); //20-30%
        npart_relerr_low.push_back(4.46); //30-40%
        npart_relerr_low.push_back(2.45); //40-60%
        npart_relerr_low.push_back(4.54); //60-90%
        meanNcoll= 8.619-1.0; // 0-90%

        centBiasCorrection.push_back(1.08072); //0-5%
        centBiasCorrection.push_back(1.04871); //5-10%
        centBiasCorrection.push_back(1.02705); //10-20%
        centBiasCorrection.push_back(1.00576); //20-30%
        centBiasCorrection.push_back(0.989589); //30-40%
        centBiasCorrection.push_back(0.966759); //40-60%
        centBiasCorrection.push_back(0.880465); //60-90%

        // rel error
        centBiasCorrection_relerr.push_back(2.0); //0-5%
        centBiasCorrection_relerr.push_back(1.0); //5-10%
        centBiasCorrection_relerr.push_back(1.0); //10-20%
        centBiasCorrection_relerr.push_back(0.0); //20-30%
        centBiasCorrection_relerr.push_back(0.0); //30-40%
        centBiasCorrection_relerr.push_back(1.0); //40-60%
        centBiasCorrection_relerr.push_back(2.0); //60-90%
    }
    else {
        std::cout << "Normal Glauber" << std::endl;

        sNameDataPlus+="Glauber";
        sNameDataMinus+="Glauber";
        sNameDataRatio+="Glauber";
        sNameDataInclusive+="Glauber";

        npart.push_back(16.53); //0-5%
        npart.push_back(14.61); //5-10%
        npart.push_back(13.05); //10-20%
        npart.push_back(11.37); //20-30%
        npart.push_back(9.806); //30-40%
        npart.push_back(7.433); //40-60%
        npart.push_back(3.976); //60-90%

        centBiasCorrection.push_back(1.2532); //0-5%
        centBiasCorrection.push_back(1.13585); //5-10%
        centBiasCorrection.push_back(1.06244); //10-20%
        centBiasCorrection.push_back(0.99944); //20-30%
        centBiasCorrection.push_back(0.955975); //30-40%
        centBiasCorrection.push_back(0.911582); //40-60%
        centBiasCorrection.push_back(0.823979); //60-90%

        // rel error
        centBiasCorrection_relerr.push_back(5.0); //0-5%
        centBiasCorrection_relerr.push_back(3.0); //5-10%
        centBiasCorrection_relerr.push_back(1.0); //10-20%
        centBiasCorrection_relerr.push_back(0.0); //20-30%
        centBiasCorrection_relerr.push_back(1.0); //30-40%
        centBiasCorrection_relerr.push_back(2.0); //40-60%
        centBiasCorrection_relerr.push_back(4.0); //60-90%

        npart_relerr_upp.push_back(11.19); //0-5%
        npart_relerr_upp.push_back(8.28); //5-10%
        npart_relerr_upp.push_back(6.29); //10-20%
        npart_relerr_upp.push_back(5.72); //20-30%
        npart_relerr_upp.push_back(5.76); //30-40%
        npart_relerr_upp.push_back(6.00); //40-60%
        npart_relerr_upp.push_back(5.18); //60-90%

        npart_relerr_low.push_back(5.61); //0-5%
        npart_relerr_low.push_back(5.64); //5-10%
        npart_relerr_low.push_back(5.62); //10-20%
        npart_relerr_low.push_back(5.55); //20-30%
        npart_relerr_low.push_back(5.85); //30-40%
        npart_relerr_low.push_back(7.12); //40-60%
        npart_relerr_low.push_back(7.37); //60-90%

        meanNcoll= 8.467-1.0; // 0-90%
    }

    std::vector<double> ncoll_relerr_upp(npart_relerr_upp);
    std::vector<double> ncoll_relerr_low(npart_relerr_low);

    for(unsigned int i=0; i<npart.size(); ++i){
        double ncoll_temp = npart[i]-1.0;
        ncoll.push_back(ncoll_temp);

        // Centrality bias correction
        if(doCentralityBiasCorr){
            std::cout << "Apply centrality bias corrective factor: " << centBiasCorrection[i] << std::endl;
            npart[i] *= centBiasCorrection[i];
            std::cout << "<Ncoll> before correction: " << ncoll[i] << std::endl;
            ncoll[i] = npart[i]-centBiasCorrection[i];
            std::cout << "<Ncoll> after correction: " << ncoll[i] << std::endl;
            npart_relerr_low[i] = TMath::Sqrt(TMath::Power(centBiasCorrection_relerr[i],2)+TMath::Power(npart_relerr_low[i],2));
            npart_relerr_upp[i] = TMath::Sqrt(TMath::Power(centBiasCorrection_relerr[i],2)+TMath::Power(npart_relerr_upp[i],2));
        }
    }//i

    for(int icent=0; icent<nCentralityBins; ++icent){
        // Initialize arrays to 0
        nWPlus[icent] = 0;
        nWPlusStatErr[icent] = 0;
        nWMinus[icent] = 0;
        nWMinusStatErr[icent] = 0;

        // For given centrality, sum over all eta windows
	    for(int ieta=0; ieta<nEtaBins; ++ieta){
        
            for(int ipt=0; ipt<nPtBins; ++ipt){

                // Retrieve yields from input file
                TString name = "WPtFit_eta"; name+=ieta; name+="_centrality"; name+=icent; name+="_charge";
		        TString namePlus = name+"Plus";
		        TString nameMinus = name+"Minus";
                
		        TGraphAsymmErrors* grPlusTemp = (TGraphAsymmErrors*)dataFileIn->Get(namePlus); 
                nWPlus[icent]+=grPlusTemp->GetY()[ipt];
                nWPlusStatErr[icent]+=TMath::Power(grPlusTemp->GetEYhigh()[ipt],2);

		        TGraphAsymmErrors* grMinusTemp = (TGraphAsymmErrors*)dataFileIn->Get(nameMinus); 
                nWMinus[icent]+=grMinusTemp->GetY()[ipt];
                nWMinusStatErr[icent]+=TMath::Power(grMinusTemp->GetEYhigh()[ipt],2);
            }//ipt
            
	    } //ieta 

	  nWPlusStatErr[icent] = TMath::Sqrt(nWPlusStatErr[icent]);
	  nWMinusStatErr[icent] = TMath::Sqrt(nWMinusStatErr[icent]);

      // Place-holder: 10% of npart for this centrality class
	  float centBinWOver2 = grCentrality->GetEXhigh()[icent];
	  float centBinW = grCentrality->GetEXhigh()[icent]*2.0;
	  //float npartBinW = 0.1*npart[icent];
      float npartBinW_low = 1.0e-2*npart_relerr_low[icent]*npart[icent];
      float npartBinW_upp = 1.0e-2*npart_relerr_upp[icent]*npart[icent];

	  std::cout << "Total W+ --> mu+ events in centrality class " << 100.*(grCentrality->GetEXhigh()[icent]-centBinWOver2) 
      << " - " << 100.*(grCentrality->GetX()[icent]+centBinWOver2) << 
		": " << nWPlus[icent] << " +/-" << nWPlusStatErr[icent] << "(stat.)" << std::endl;  
	  std::cout << "Total W- --> mu- events in centrality class " << 100.*(grCentrality->GetX()[icent]-centBinWOver2) 
      << " - " << 100.*(grCentrality->GetX()[icent]+centBinWOver2) << 
		": " << nWMinus[icent] << " +/-" << nWMinusStatErr[icent] << "(stat.)" << std::endl; 

      // Normalize to level of per pA event per binary collision
      float pPbCrossSec = 2.1; //barn
      float runLumi;
      // Place-holders
      if(doPeriodA) runLumi = 0;
      else if(doPeriodB) runLumi = 0;
      //else runLumi = 2.58e10;
      else runLumi = 27.8e9; // barn^{-1}
      float norm = 1.0e9/(centBinW*pPbCrossSec*runLumi*ncoll[icent]);

      grDataPlus->SetPoint(icent,npart[icent],nWPlus[icent]*norm);
      grDataPlus->SetPointError(icent,npartBinW_low,npartBinW_upp,nWPlusStatErr[icent]*norm,nWPlusStatErr[icent]*norm);

      grDataMinus->SetPoint(icent,npart[icent],nWMinus[icent]*norm);
      grDataMinus->SetPointError(icent,npartBinW_low,npartBinW_upp,nWMinusStatErr[icent]*norm,nWMinusStatErr[icent]*norm);

      // W+/W- 
      float chRatio = nWPlus[icent]/nWMinus[icent]; 
      float chRatio_err_upp = chRatio*TMath::Sqrt(TMath::Power(nWPlusStatErr[icent]/nWPlus[icent],2)+TMath::Power(nWMinusStatErr[icent]/nWMinus[icent],2));
      float chRatio_err_low = chRatio_err_upp;
      grDataRatio->SetPoint(icent,npart[icent],chRatio);
      grDataRatio->SetPointError(icent,npartBinW_low,npartBinW_upp,chRatio_err_low,chRatio_err_upp);

      grDataInclusive->SetPoint(icent,npart[icent],(nWPlus[icent]+nWMinus[icent])*norm);
      float errDataInclusive = TMath::Sqrt(TMath::Power(nWMinusStatErr[icent]*norm,2)+TMath::Power(nWPlusStatErr[icent]*norm,2));
      grDataInclusive->SetPointError(icent,npartBinW_low,npartBinW_upp,errDataInclusive,errDataInclusive);


  }//icent

  // For systematic errors
  TGraphAsymmErrors* grDataSystErrorsPlus = (TGraphAsymmErrors*)grDataPlus->Clone();
  TGraphAsymmErrors* grDataSystErrorsMinus = (TGraphAsymmErrors*)grDataMinus->Clone();
  TGraphAsymmErrors* grDataSystErrorsInclusive = (TGraphAsymmErrors*)grDataInclusive->Clone();
  for(int igr=0; igr<grDataSystErrorsInclusive->GetN(); ++igr){
        double yieldInc = grDataSystErrorsInclusive->GetY()[igr];
        double yieldPlus = grDataSystErrorsPlus->GetY()[igr];
        double yieldMinus = grDataSystErrorsMinus->GetY()[igr];
        grDataSystErrorsInclusive->SetPointError(igr, grDataSystErrorsMinus->GetEXlow()[igr], grDataSystErrorsMinus->GetEXhigh()[igr], 
                ncoll_relerr_low[igr]/100.0*yieldInc, ncoll_relerr_upp[igr]/100.0*yieldInc);
        grDataSystErrorsPlus->SetPointError(igr, grDataSystErrorsMinus->GetEXlow()[igr], grDataSystErrorsMinus->GetEXhigh()[igr], 
                ncoll_relerr_low[igr]/100.0*yieldPlus, ncoll_relerr_upp[igr]/100.0*yieldPlus);
        grDataSystErrorsMinus->SetPointError(igr, grDataSystErrorsMinus->GetEXlow()[igr], grDataSystErrorsMinus->GetEXhigh()[igr], 
                ncoll_relerr_low[igr]/100.0*yieldMinus, ncoll_relerr_upp[igr]/100.0*yieldMinus);
  }//igr

  // Save histos
  outFile->cd();
  grDataPlus->Write(sNameDataPlus);
  grDataMinus->Write(sNameDataMinus);
  grDataRatio->Write(sNameDataRatio);
  grDataInclusive->Write(sNameDataInclusive);

  grDataSystErrorsPlus->Write(sNameDataPlus+"_syst");
  grDataSystErrorsMinus->Write(sNameDataMinus+"_syst");
  grDataSystErrorsInclusive->Write(sNameDataInclusive+"_syst");

  // XSec from MC Model (CT10)
  // [nb]
  // AMI
  /*const float sigma_pp_plus = 4.1936; 
  const float sigma_pn_plus = 3.4630; 
  const float sigma_pp_minus = 2.7368; 
  const float sigma_pn_minus = 3.4410; 
  */
  // POWHEG
  const float sigma_pp_plus = 4.297; 
  const float sigma_pn_plus = 3.588; 
  const float sigma_pp_minus = 2.901; 
  const float sigma_pn_minus = 3.595; 
  std::vector<double> centralityBins;
  centralityBins.push_back(0.0);
  centralityBins.push_back(0.9);
  std::vector<double> etaBins;
  etaBins.push_back(-2.5);
  etaBins.push_back(2.5);
  TString sFileName = "/usatlas/u/tbales/scratch/Data2013Hip/";
  std::unique_ptr<TGraphErrors> gr_CT10_plus(new TGraphErrors(1));
  std::unique_ptr<TGraphErrors> gr_CT10_minus(new TGraphErrors(1));
  std::unique_ptr<TGraphErrors> gr_CT10_inclusive(new TGraphErrors(1));
  float nW_CT10_pp_plus = getModelXSec(sFileName+"ntupleWplusmunu_pp.B.PowhegPythia8_AU2CT10.09.24.2014.root",sigma_pp_plus
                                ,centralityBins[0],centralityBins[1],etaBins[0],etaBins[1],"+");
  float nW_CT10_pn_plus = getModelXSec(sFileName+"ntupleWplusmunu_pn.B.PowhegPythia8_AU2CT10.09.24.2014.root",sigma_pn_plus
                                ,centralityBins[0],centralityBins[1],etaBins[0],etaBins[1],"+");
  float nW_CT10_pp_minus = getModelXSec(sFileName+"ntupleWminmunu_pp.B.PowhegPythia8_AU2CT10.09.24.2014.root" ,sigma_pp_minus
                                ,centralityBins[0],centralityBins[1],etaBins[0],etaBins[1],"-");
  float nW_CT10_pn_minus = getModelXSec(sFileName+"ntupleWminmunu_pn.B.PowhegPythia8_AU2CT10.09.24.2014.root" ,sigma_pn_minus
                                ,centralityBins[0],centralityBins[1],etaBins[0],etaBins[1],"-");
  // weight by collision frequency
  double weight_pp = 82.0/208.0;
  double weight_pn = (208.0-82.0)/208.0;
  float nW_CT10_plus = weight_pp*nW_CT10_pp_plus + weight_pn*nW_CT10_pn_plus;
  float nW_CT10_plus_err = 7.14/100.0*nW_CT10_plus;
  float nW_CT10_minus = weight_pp*nW_CT10_pp_minus + weight_pn*nW_CT10_pn_minus;
  float nW_CT10_minus_err = 7.14/100.0*nW_CT10_minus;
  std::cout << "nW_CT10_plus: " << nW_CT10_plus*1e9 << ", nW_CT10_minus: " << nW_CT10_minus*1e9 << std::endl; 
  float nW_CT10_inclusive = nW_CT10_plus+nW_CT10_minus;
  float nW_CT10_inclusive_err = TMath::Sqrt(TMath::Power(nW_CT10_minus_err,2)+TMath::Power(nW_CT10_plus_err,2));

  gr_CT10_plus->SetPoint(0,2.0,nW_CT10_plus*1e9);
  gr_CT10_plus->SetPointError(0,0.0,nW_CT10_plus_err*1e9);
  gr_CT10_minus->SetPoint(0,2.0,nW_CT10_minus*1e9);
  gr_CT10_minus->SetPointError(0,0.0,nW_CT10_minus_err*1e9);
  gr_CT10_inclusive->SetPoint(0,2.0,nW_CT10_inclusive*1e9);
  gr_CT10_inclusive->SetPointError(0,0.0,nW_CT10_inclusive_err*1e9);
  
  outFile->cd();
  // Plot on canvas
  TCanvas* c1 = new TCanvas(cName,cName,600,600);
  c1->cd();
  TPad *c_1 = new TPad("c_1", "c_1",0,0.008474576,1,1);
  c_1->Draw();
  c_1->cd();
  c_1->Range(-1.826128,-24.92296,25.74476,65.96148);
  c_1->SetFillColor(0);
  c_1->SetBorderMode(0);
  c_1->SetBorderSize(2);
  c_1->SetTickx(1);
  c_1->SetTicky(1);
  c_1->SetLeftMargin(0.1577181);
  c_1->SetRightMargin(0.05201342);
  c_1->SetTopMargin(0.01057916);
  c_1->SetBottomMargin(0.2753272);
  c_1->SetFrameBorderMode(0);
  c_1->SetFrameBorderMode(0);

  //TH1F *htemplate = new TH1F("htemplate","htemplate",100,2.5223,24.3107);
  TH1F *htemplate = new TH1F("htemplate","htemplate",100,0.0,24.3107);
  //TH1F *htemplate = new TH1F("htemplate","htemplate",100,0.0,34.3107);
  htemplate->SetMinimum(0.1);
  htemplate->SetMaximum(80);
  htemplate->SetDirectory(0);
  htemplate->SetStats(0);
  htemplate->SetLineWidth(2);
  htemplate->SetMarkerStyle(20);
  htemplate->SetMarkerSize(1.2);
  htemplate->GetXaxis()->SetLabelFont(42);
  htemplate->GetXaxis()->SetLabelSize(0.05);
  htemplate->GetXaxis()->SetTitleSize(0.05);
  htemplate->GetXaxis()->SetTitleOffset(1.4);
  htemplate->GetXaxis()->SetTitleFont(42);
  htemplate->GetYaxis()->SetTitle("#frac{N_{W#rightarrow#mu}^{fiducial}}{#LT N_{coll} #GT} #frac{x10^{9}}{N_{pPb}}");
  htemplate->GetYaxis()->SetLabelFont(42);
  htemplate->GetYaxis()->SetLabelSize(0.03);
  htemplate->GetYaxis()->SetTitleSize(0.03);
  htemplate->GetYaxis()->SetTitleOffset(1.97);
  htemplate->GetYaxis()->SetTitleFont(42);
  htemplate->Draw();

  // mu+-
  grDataSystErrorsInclusive->SetFillColor(kGray+2);
  grDataSystErrorsInclusive->SetFillStyle(3001);
  grDataSystErrorsInclusive->Draw("pe5same");
  grDataInclusive->Draw("pesame");
  // CT10
  gr_CT10_inclusive->SetMarkerColor(kGray+2);
  gr_CT10_inclusive->SetMarkerStyle(21);
  gr_CT10_inclusive->Draw("pesame");

  // mu+
  grDataSystErrorsPlus->SetFillColor(kBlue);
  grDataSystErrorsPlus->SetFillStyle(3001);
  grDataSystErrorsPlus->Draw("pe5same");
  grDataPlus->Draw("pesame");
  // CT10
  gr_CT10_plus->Print();
  gr_CT10_plus->SetMarkerColor(kBlue);
  gr_CT10_plus->SetMarkerStyle(21);
  gr_CT10_plus->Draw("pesame");

  // mu-
  grDataSystErrorsMinus->SetFillColor(kRed);
  grDataSystErrorsMinus->SetFillStyle(3001);
  grDataSystErrorsMinus->Draw("pe5same");
  grDataMinus->Draw("pesame");
  // CT10
  gr_CT10_minus->Print();
  gr_CT10_minus->SetMarkerColor(kRed);
  gr_CT10_minus->SetMarkerStyle(21);
  gr_CT10_minus->Draw("pesame");

  //TLegend *leg = new TLegend(0.614094,0.3088279,0.8708054,0.4992529,NULL,"brNDC");
  TLegend *leg = new TLegend(0.7643098,0.3280423,0.9124579,0.4603175,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  //leg->SetTextSize(0.05);
  leg->SetTextSize(0.03527337);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  TLegendEntry *entry=leg->AddEntry(grDataSystErrorsInclusive,"W^{#pm}#rightarrow#mu^{#pm}","pef");
  entry->SetFillColor(15);
  entry->SetFillStyle(3001);
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(2);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1.2);
  entry->SetTextFont(42);
  entry=leg->AddEntry(grDataSystErrorsPlus,"W^{+}#rightarrow#mu^{+}","pef");
  entry->SetFillColor(4);
  entry->SetFillStyle(3001);
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(2);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1.2);
  entry->SetTextFont(42);
  entry=leg->AddEntry(grDataSystErrorsMinus,"W^{-}#rightarrow#mu^{-}","pef");
  entry->SetFillColor(2);
  entry->SetFillStyle(3001);
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(2);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1.2);
  entry->SetTextFont(42);
  leg->Draw();

  leg = new TLegend(0.533557,0.3352758,0.6040268,0.4498835,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03527337);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  TGraphErrors* grDummy=(TGraphErrors*)grDataSystErrorsInclusive->Clone("grDummy");
  grDummy->SetMarkerColor(kGray+2);
  grDummy->SetMarkerStyle(21);
  entry=leg->AddEntry(grDummy,"W^{#pm}#rightarrow#mu^{#pm} (CT10)","p");
  entry->SetFillColor(2);
  entry->SetFillStyle(3001);
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(2);
  entry->SetMarkerColor(kGray+2);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1.2);
  entry->SetTextFont(42);
  TGraphErrors* grPlusDummy=(TGraphErrors*)grDataSystErrorsPlus->Clone("grPlusDummy");
  grPlusDummy->SetMarkerColor(kBlue);
  grPlusDummy->SetMarkerStyle(21);
  entry=leg->AddEntry(grPlusDummy,"W^{+}#rightarrow#mu^{+} (CT10)","p");
  entry->SetFillColor(2);
  entry->SetFillStyle(3001);
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(2);
  entry->SetMarkerColor(kBlue);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1.2);
  entry->SetTextFont(42);
  TGraphErrors* grMinusDummy=(TGraphErrors*)grDataSystErrorsMinus->Clone("grMinusDummy");
  grMinusDummy->SetMarkerColor(kRed);
  grMinusDummy->SetMarkerStyle(21);
  entry=leg->AddEntry(grMinusDummy,"W^{-}#rightarrow#mu^{-}  (CT10)","p");
  entry->SetFillColor(2);
  entry->SetFillStyle(3001);
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(2);
  entry->SetMarkerColor(kRed);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1.2);
  entry->SetTextFont(42);
  leg->Draw();

  TLatex *   tex = new TLatex(0.6258389,0.9259458,"ATLAS Internal");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.1946309,0.4234355,"p+Pb #sqrt{s_{NN}}=5.02 TeV");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.03526388);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(0.204698,0.3405654,"#int Ldt = 27.8 nb^{-1}");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.03526388);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(0.307047,0.6808619,("Glauber-Gribov (#Omega = "+omega+")").c_str());
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->Draw();

  // W+/W- Ratio
    
  TPad*    c_2 = new TPad("c_2", "c_2",0,0.003527337,1,0.2716049);
  c_2->Draw();
  c_2->cd();
  c_2->Range(-1.826128,0.6060205,25.74476,1.399862);
  c_2->SetFillColor(0);
  c_2->SetBorderMode(0);
  c_2->SetBorderSize(2);
  c_2->SetTickx(1);
  c_2->SetTicky(1);
  c_2->SetLeftMargin(0.1577181);
  c_2->SetRightMargin(0.05201342);
  c_2->SetTopMargin(0.01242356);
  c_2->SetBottomMargin(0.2443553);
  c_2->SetFrameBorderMode(0);
  c_2->SetFrameBorderMode(0);

  //TH1F *hRatioTemplate = new TH1F("hRatioTemplate","hRatioTemplate",100,2.5223,24.3107);
  TH1F *hRatioTemplate = new TH1F("hRatioTemplate","hRatioTemplate",100,0.0,24.3107);
  //TH1F *hRatioTemplate = new TH1F("hRatioTemplate","hRatioTemplate",100,0.0,34.3107);
  hRatioTemplate->SetMinimum(0.8);
  hRatioTemplate->SetMaximum(1.39);
  hRatioTemplate->SetDirectory(0);
  hRatioTemplate->SetStats(0);
  hRatioTemplate->SetLineWidth(2);
  hRatioTemplate->SetMarkerStyle(20);
  hRatioTemplate->SetMarkerSize(1.2);
  hRatioTemplate->GetXaxis()->SetTitle("#LT N_{part} #GT");
  hRatioTemplate->GetXaxis()->SetLabelFont(42);
  hRatioTemplate->GetXaxis()->SetLabelSize(0.09);
  hRatioTemplate->GetXaxis()->SetTitleSize(0.09);
  hRatioTemplate->GetXaxis()->SetTitleOffset(1.22);
  hRatioTemplate->GetXaxis()->SetTitleFont(42);
  hRatioTemplate->GetYaxis()->SetTitle("W+/W-");
  hRatioTemplate->GetYaxis()->SetNdivisions(410);
  hRatioTemplate->GetYaxis()->SetLabelFont(42);
  hRatioTemplate->GetYaxis()->SetLabelSize(0.10);
  hRatioTemplate->GetYaxis()->SetTitleSize(0.14);
  hRatioTemplate->GetYaxis()->SetTitleOffset(0.39);
  hRatioTemplate->GetYaxis()->SetTitleFont(42);
  hRatioTemplate->Draw();

  grDataRatio->Draw("pesame");
  c1->Print("plots/"+cName+".root");
  c1->Write(cName);

  
  delete grDataPlus;
  delete grDataMinus;
  delete grDataInclusive;
  delete grDataSystErrorsPlus;
  delete grDataSystErrorsMinus;
  delete grDataSystErrorsInclusive;
  delete htemplate;
  delete hRatioTemplate;
  delete c1;

}

int main(){
    SetAtlasStyle();

	TFile* outFile = new TFile("HipBinaryScalingHistos.root","recreate");
    std::string omega055 = "0.55";
    std::string omega101 = "1.01";
    // Applying centrality bias correction
	run("AB", omega055, outFile,true);	
	run("AB", omega101, outFile,true);	
	run("AB", "0", outFile,true);	
	outFile->Close();
	return 0;
}
