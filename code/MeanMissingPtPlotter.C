///////////////////////////
//This macro plots the mean missing
//pT as a function of centrality
//for W signal, QCD, and Z/DY bkg
//@date: Oct.22, 2012
//////////////////////////

void MeanMissingPtPlotter() {

//TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
TFile *_fileD = TFile::Open(baseString+"HardProbesFiles/HISingleMuonHardProbesData.04.17.2013.root");

//Slimmed muon ntuple from Wmunu+DataOverlay sample 
TString fileNameMCWIn = "MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.04.12.2013";
TFile *_fileS = TFile::Open(baseString+fileNameMCWIn+".root");

TCanvas cDummy = TCanvas("cDummy","cDummy",600,600);

//stl container to store lower 
//track thresholds
std::vector<int> missingPThr ;
//missingPThr.push_back(500);
//missingPThr.push_back(1000);
missingPThr.push_back(3000);
//missingPThr.push_back(5000);
//missingPThr.push_back(7000);
const unsigned int nMissPtBins = missingPThr.size();
std::cout << " Number of lower track thresholds: " << nMissPtBins << std::endl;

std::vector<double> centralityBins;
centralityBins.push_back(0.0);
centralityBins.push_back(0.05);
centralityBins.push_back(0.1);
centralityBins.push_back(0.15);
centralityBins.push_back(0.2);
centralityBins.push_back(0.4);
centralityBins.push_back(0.8);
const unsigned int nCentralityBins = centralityBins.size()-1;
std::cout << " Number of centrality classes: " << nCentralityBins << std::endl;

std::vector <float> ncoll;
ncoll.push_back(1683.3); //0-5
ncoll.push_back(1318.0); //5-10*/
ncoll.push_back(1035.4); //10-15
ncoll.push_back(811.2); //15-20
ncoll.push_back(440.6); //20-40
ncoll.push_back(77.8); //40-80*/
double pts[] = {77.8,440.6,811.2,1035.4,1318.0,1683.3}  ;
int  binnum = sizeof(pts)/sizeof(double) - 1;

char namePtSig[50],namePtPreSel[50],namePtMC[50] ;
float meanPtSig[nCentralityBins],errPtSig[nCentralityBins], meanPtPreSel[nCentralityBins], errPtPreSel[nCentralityBins], meanPtMC[nCentralityBins],errPtMC[nCentralityBins];

#ifdef __CINT__
  gROOT->LoadMacro("AtlasUtils.C");
#endif

TH1 *hPtSig[nMissPtBins][nCentralityBins],*hPtPreSel[nMissPtBins][nCentralityBins],*hPtMC[nMissPtBins][nCentralityBins] ;

int markerColor[nCentralityBins] = {kBlue,kRed,kGreen,kCyan};
int markerStyle[nCentralityBins] = {kOpenStar,kOpenCircle,kOpenDiamond,kOpenDiamond};

//double pts[nCentralityBins] = {77.8,440.6,811.2,1035.4,1318.0,1683.3} ;
int binFrom[nCentralityBins]  ;
int binTo[nCentralityBins] ;


/*for(int icent=0;icent<nCentralityBins;icent++){
	binFrom[icent] = centralityBins.at(icent)*100.0;
	binTo[icent] = centralityBins.at(icent+1)*100.0;
}*/

TH1F* hMeanSig = new TH1F("hMeanSig","hMeanSig",binnum,pts);
//hMeanSig->Sumw2();

TH1F* hMeanPreSel = new TH1F("hMeanPreSel","hMeanPreSel",binnum,pts);
//hMeanPreSel->Sumw2();
TH1F* hMeanMC = new TH1F("hMeanMC","hMeanMC",binnum,pts);
//hMeanMC->Sumw2();


std::cout << "aaaaa" << std::endl;
TGraphErrors* grMeanPtSig = new TGraphErrors(nCentralityBins); 		
TGraphErrors* grMeanPtPreSel = new TGraphErrors(nCentralityBins); 		
TGraphErrors* grMeanPtMC = new TGraphErrors(nCentralityBins); 		

grMeanPtSig->SetMarkerStyle(kOpenCircle);
grMeanPtSig->SetMarkerColor(kRed);
//grMeanPtMC->SetMarkerStyle(kOpenStar);
//grMeanPtMC->SetMarkerColor(kBlue);
grMeanPtPreSel->SetMarkerStyle(kOpenDiamond);
grMeanPtPreSel->SetMarkerColor(kGreen);

for(int ipt=0; ipt<nMissPtBins;ipt++){
	for(int icent=0;icent<nCentralityBins;icent++){
		
		_fileD->cd();

		sprintf(namePtSig,"hPtSig_%i_cent_%i",ipt,icent);
		sprintf(namePtPreSel,"hPtPreSel_%i_cent_%i",ipt,icent);
		
		hPtSig[ipt][icent] = new TH1F(namePtSig,namePtSig,90,0.0,180.0);
		hPtPreSel[ipt][icent] = new TH1F(namePtPreSel,namePtPreSel,90,0.0,180.0);

		double centBinLo = centralityBins.at(icent);
		double centBinUp = centralityBins.at(icent+1);
//		double centBinCenter = 0.5*(centBinLo+centBinUp);
		double centBinWidth = fabs(centBinUp-centBinLo);
	  	double mtcutLow = 40.0;
		double etaLow = 0.1; double etaUpp = 2.4;

		TString sPtCutPreSel = "&&pt>"; sPtCutPreSel += 0.0; sPtCutPreSel +="&&pt<"; sPtCutPreSel += 300.0;
		TString sPtCutSig = "&&pt>"; sPtCutSig += 25.0; sPtCutSig +="&&pt<"; sPtCutSig += 300.0;
		TString centBinCut = "&&centrality>="; centBinCut+=centBinLo; centBinCut+= "&&centrality<"; centBinCut+= centBinUp;		
		TString sHiPtTrigger =
        "&&((EF_mu10_MSonly_EFFS_L1ZDC&&EF_mu10_MSonly_EFFS_L1ZDC_Matched20)||(EF_mu10_MSonly_EFFS_L1TE10&&EF_mu10_MSonly_EFFS_L1TE10_Matched20)||(EF_mu10_MSonly_EFFS_L1TE20&&EF_mu10_MSonly_EFFS_L1TE20_Matched20))";
		TString sMBTrigger = "&&(EF_L1TE50_NoAlg||EF_mbZdc_a_c_L1VTE50_trk)";

		TString sLoPtTrigger = "&&(EF_mu4_MSonly_L1TE50||EF_mu4_L1VTE50)&&(EF_mu4_MSonly_L1TE50_Matched20||EF_mu4_L1VTE50_Matched20)";
	 	TString preSelCuts = "val>11&&abs(eLoss)<0.5&&abs(scat)<4.0"; preSelCuts+=centBinCut; //preSelCuts+=sLoPtTrigger; 	
		preSelCuts+="&&abs(eta)>"; preSelCuts += etaLow; preSelCuts+="&&abs(eta)<"; preSelCuts += etaUpp; preSelCuts+=sPtCutPreSel;

	 	TString sigCuts = "val>11&&mt>40.0&&mt<300.0&&ptcone20ID3/pt<0.1&&abs(eLoss)<0.5&&abs(scat)<4.0&&ZDY==0"; sigCuts+=centBinCut; sigCuts+=sPtCutSig;
		sigCuts+="&&abs(eta)>"; sigCuts += etaLow; sigCuts+="&&abs(eta)<"; sigCuts += etaUpp; 
		TString sigCutsData = sigCuts + sHiPtTrigger;
		TString sigCutsMC = sigCuts; sigCutsMC+="&&prompt==24";

		std::cout << "Pre-Selection region: " << preSelCuts << "\n" << std::endl;
		std::cout << "Signal region: " << sigCutsData << "\n" << std::endl;

		tree->Draw((TString("nu_pt>>")+namePtSig).Data(),sigCutsData);

	   	meanPtSig[icent] = hPtSig[ipt][icent]->GetMean(); errPtSig[icent] = hPtSig[ipt][icent]->GetMeanError();	
		std::cout << "Mean signal pt: " << ipt << ":" << icent << " "<< ncoll.at(icent) << " " << meanPtSig[icent] << " +- " << errPtSig[icent] << std::endl;

		tree->Draw((TString("nu_pt>>")+namePtPreSel).Data(),preSelCuts);
	   	meanPtPreSel[icent] = hPtPreSel[ipt][icent]->GetMean(); errPtPreSel[icent] = hPtPreSel[ipt][icent]->GetMeanError();	
		std::cout << "Mean PreSel pt: " << ipt << ":" << icent << " "<< ncoll.at(icent) << " " << meanPtPreSel[icent] << " +- " << errPtPreSel[icent] << std::endl;


		//Monte Carlo
		_fileS->cd();
		
		sprintf(namePtMC,"hPtMC_%i_cent_%i",ipt,icent);
		hPtMC[ipt][icent] = new TH1F(namePtMC,namePtMC,90,0.0,180.0);

		tree->Draw((TString("nu_pt>>")+namePtMC).Data(),sigCutsMC);
	   	meanPtMC[icent] = hPtMC[ipt][icent]->GetMean(); errPtMC[icent] = hPtMC[ipt][icent]->GetMeanError();	
		std::cout << "Mean Wmunu MC pt: " << ipt << ":" << icent << " " << ncoll.at(icent) << " " << meanPtMC[icent] << " +- " << errPtMC[icent] << std::endl;

		grMeanPtSig->SetPoint(icent,ncoll.at(icent),meanPtSig[icent]);
		grMeanPtSig->SetPointError(icent,0.0,errPtSig[icent]);

		grMeanPtPreSel->SetPoint(icent,ncoll.at(icent),meanPtPreSel[icent]);
		grMeanPtPreSel->SetPointError(icent,0.0,errPtPreSel[icent]);

		grMeanPtMC->SetPoint(icent,ncoll.at(icent),meanPtMC[icent]);
		grMeanPtMC->SetPointError(icent,0.0,errPtMC[icent]);


	}
    for(int icoll=nCentralityBins; icoll>0; --icoll){
		    hMeanSig->SetBinContent(nCentralityBins-icoll,meanPtSig[icoll-1]);
		    hMeanSig->SetBinError(nCentralityBins-icoll,errPtSig[icoll-1]);

		    hMeanPreSel->SetBinContent(nCentralityBins-icoll,meanPtPreSel[icoll-1]);
		    hMeanPreSel->SetBinError(nCentralityBins-icoll,errPtPreSel[icoll-1]);

		    hMeanMC->SetBinContent(nCentralityBins-icoll,meanPtMC[icoll-1]);
		    hMeanMC->SetBinError(nCentralityBins-icoll,errPtMC[icoll-1]);
		}
}

	TCanvas *c0 = new TCanvas("c0","c0",600,600);


	/*TH1F* hdummy = new TH1F("hdummy","hdummy",4,0,80);
	hdummy->SetXTitle("#LT N_{part} #GT");
	hdummy->SetYTitle("#LT #slash{p_{T}} #GT");
	hdummy->GetYaxis()->SetRangeUser(0.0,55.0);
	hdummy->Draw();


	c0->Update();

*/
	TH1* hMeanMCc = hMeanMC->Clone("hMeanMCc");
	hMeanMCc->SetXTitle("#LT N_{part} #GT");
	hMeanMCc->SetYTitle("#LT #slash{p_{T}} #GT");
	//hMeanMCc->GetYaxis()->SetRangeUser(0.0,1700.0);
	TH1* hMeanSigc = hMeanSig->Clone("hMeanSigc");
	TH1* hMeanPreSelc = hMeanPreSel->Clone("hMeanPreSelc");

	//hMeanMCc->SetOption("hf same");
	hMeanMCc->Draw("histf");

	//grMeanPtSig->Draw("ape");
	//grMeanPtSig->GetYaxis()->SetRangeUser(0.0,55.0);
	//grMeanPtSig->GetXaxis()->SetTitle("#LT N_{part} #GT");
	//grMeanPtSig->GetYaxis()->SetTitle("#LT #slash{p_{T}} #GT");
	c0->Update();
	//hMeanPreSel->SetOption("pe same");
	hMeanPreSelc->Draw("pe same");
	c0->Update();
	//hMeanSig->SetOption("pesame");
	hMeanSigc->Draw("pesame");
	//grMeanPtPreSel->SetMarkerStyle(kFullCircle); grMeanPtPreSel->SetMarkerStyle(22);
	//grMeanPtPreSel->Draw("pe same");
    hMeanMCc->Draw("sameaxis");
	c0->Update();

	hMeanSigc->SetMarkerStyle(kFullCircle);
	hMeanPreSelc->SetMarkerStyle(22);
	//hMeanPreSel->SetFillColor(kRed);
	hMeanMC->SetFillColor(kWhite);
	//grMeanPtMC->Draw("hist same");
	c0->Update();

	TLegend* leg1 = new TLegend(0.18, 0.78, 0.53, 0.92);
	leg1->SetTextFont(gStyle->GetTextFont());
	leg1->SetTextSize(gStyle->GetTextSize());
	leg1->SetBorderSize(0);
	leg1->SetFillColor(0);
	leg1->AddEntry(hMeanSigc, "Signal", "p");
	leg1->AddEntry(hMeanPreSelc, "Pre-Selection", "p");
	leg1->AddEntry(hMeanMC, "W#rightarrow#mu#nu", "f");
	leg1->Draw();
    c0->Print("meanMissingPtCentDep.pdf");
    c0->Print("meanMissingPtCentDep.root");

}
