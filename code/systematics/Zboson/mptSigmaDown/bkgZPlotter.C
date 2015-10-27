///////////////////////////////////////////////////////////////////////////////
//plot difference in charge as function of Ncoll 
///////////////////////////////////////////////////////////////////////////////
void plotChargeDiffNcoll(TGraph* grDiff, TGraphErrors* grPlus, TGraphErrors* grMinus, int ieta, int icent, double npart){

	double* yPlus = grPlus->GetY();
	double* yMinus = grMinus->GetY();

	double yTempPlus = yPlus[icent]; 
	double yTempMinus = yMinus[icent]; 

	double diff = fabs(yTempPlus-yTempMinus);
	std::cout << "Difference in f_{Z} at eta bin " << ieta << " and centrality bin " << icent << " = " << diff << std::endl;
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
	std::cout << "Difference in f_{Z} at eta bin " << ieta << " and centrality bin " << icent << " = " << diff << std::endl;
	grDiff->SetPoint(ieta,xPt,diff);

}


void plotFraction(TGraphErrors* grFracCent, TGraphErrors* grFracEta, TString sCh){

	TLatex l;
	TCanvas* cFracCent = new TCanvas("cFracCent","cFracCent",600,600);
	l.SetNDC();
	l.DrawLatex(0.39,0.83,sCh);
	grFracCent->GetXaxis()->SetTitle("#LT N_{coll} #GT");	
	grFracCent->GetYaxis()->SetTitle("#frac{N_{Z}^{sur}}{N_{Z}^{Gen}}");
	grFracCent->Draw("ape");	

	TF1 f0a = TF1("f0a","[0]",0,1700);
	f0a.SetLineStyle(kDashed);
	grFracCent->Fit("f0a");
	f0a.Draw("same") ; 
	cout << "chi2 = " << f0a.GetChisquare() << "/" << f0a.GetNDF() << ", p = " << f0a.GetProb() << endl;

	cFracCent->Update();
	cFracCent->Print(sCh+"_centrality.pdf"); cFracCent->Print(sCh+"_centrality.root");

	TCanvas* cFracEta = new TCanvas("cFracEta","cFracEta",600,600);
	l.SetNDC();
	l.DrawLatex(0.39,0.83,sCh);
	grFracEta->GetXaxis()->SetTitle("|#eta|");	
	grFracEta->GetYaxis()->SetTitle("#frac{N_{Z}^{sur}}{N_{Z}^{Gen}}");
	grFracEta->Draw("ape");	

	TF1 f0b = TF1("f0b","[0]",0,2.5);
	f0b.SetLineStyle(kDashed);
	grFracEta->Fit("f0b");
	f0b.Draw("same") ; 
	cout << "chi2 = " << f0b.GetChisquare() << "/" << f0b.GetNDF() << ", p = " << f0b.GetProb() << endl;

	cFracEta->Update();
	cFracEta->Print(sCh+"_eta.pdf"); cFracEta->Print(sCh+"_eta.root");
}

void Write(TFile* outFile, TGraph2DErrors* grFrac, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing T2DGraph to root file..." << std::endl;
	  grFrac->Write(sFile);
	  std::cout << "Done." << std::endl;
    	  gDirectory = dir;
}

void Write(TFile* outFile, TObject* grFrac, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing TGraph to root file..." << std::endl;
	  grFrac->Write(sFile);
	  std::cout << "Done." << std::endl;
    	  gDirectory = dir;
}


void plot(TFile* fMcZSet, TGraph2DErrors* grFracEtaCent, TGraphErrors* grFracCent, TGraphErrors* grFracEta,
		const float mtmax, TString cuts, TString sGenChargeCut, int index, double ptLow, double ptUpp, double etaLow, double etaUpp, 
		double centralityLow, double centralityUpp, double ncoll, int ieta, int icent, int iCh, TString sSel, TString sSel2, bool doEta){

		
		TH1F* hmcZSetFull = new TH1F("hmcZSetFull","hmcZSetFull",80,0.0,mtmax);
		TH1F* hmcZSetCut = new TH1F("hmcZSetCut","hmcZSetCut",80,0.0,mtmax);
		TH1F* hmcZSetPreSelCut = new TH1F("hmcZSetPreSelCut","hmcZSetPreSelCut",80,0.0,mtmax);
		TH1F* hmcZSetPreSelAllFidCut = new TH1F("hmcZSetPreSelAllFidCut","hmcZSetPreSelAllFidCut",80,0.0,mtmax);
		TH1F* hmcZSet = new TH1F("hmcZSet","hmcZSet",80,0.0,mtmax);


		TString sCentCuts = "centrality>"; sCentCuts+=centralityLow; sCentCuts+="&&centrality<";  sCentCuts+=centralityUpp;

		double mtcutLow = 40.0;
		cuts+="&&pt>"; cuts += ptLow; cuts+="&&pt<"; cuts += ptUpp;
		TString cutsAllFid = cuts; cutsAllFid+="&&abs(eta)>0.0&&abs(eta)<2.5";
		cuts+="&&abs(eta)>"; cuts += etaLow; cuts+="&&abs(eta)<"; cuts += etaUpp;
		cuts+="&&centrality>"; cuts += centralityLow; cuts+="&&centrality<"; cuts += centralityUpp;
		cutsAllFid+="&&centrality>"; cutsAllFid += centralityLow; cutsAllFid+="&&centrality<"; cutsAllFid += centralityUpp;
		
		//systematics
		double sigmaMPt = 11.05;
		std::cout << "Doing Systematics for MPt -1 sigma " << std::endl;
		double missPtCut = 25.0-sigmaMPt; 
		//double isoCut = 1.0;
		double isoCut = 0.3;
		TString sMissPtLo = "&&nu_pt<"; sMissPtLo+=missPtCut;  TString sMissPtUp = "&&nu_pt>"; sMissPtUp+=missPtCut; 
		TString sIsoCut = "&&ptcone20/pt<"; sIsoCut+=isoCut; 
		TString sMtCutLo = "&&mt>"; sMtCutLo+=mtcutLow; TString sMtCutUp = "&&mt<";sMtCutUp+=mtmax;

		TString scutsMc = cuts; scutsMc += sMissPtUp; scutsMc += sIsoCut; scutsMc +=sMtCutLo; scutsMc +=sMtCutUp; 

		//AMI cross-section and efficiency
		double wtZ = 0.99843*2.5743e-10/64.0e-3;

		double mbEvents = 68.7e6;
		double arrCentWidth = centralityUpp-centralityLow;
		double scaleFactor = arrCentWidth*ncoll*mbEvents;

		TTree *treeMcZSet = (TTree*)fMcZSet->Get("tree");

		TString sFull = sCentCuts + "&&prompt==23";
		treeMcZSet->Draw("pt>>hmcZSetFull",sFull,"hf");		

		TH1F* hmcZSetFullc = (TH1F*)hmcZSetFull->Clone("hmcZSetFullc");

		TString sGenCuts = sCentCuts; sGenCuts+="&&abs(mc_mu_gen_mothertype)==23&&abs(mc_mu_gen_type)==13"; sGenCuts+=sGenChargeCut;
		TString sFullFid = sGenCuts;
		TString sGenEtaCuts; 
		if(doEta) { 
			sGenEtaCuts = "&&abs(mc_mu_gen_eta)>"; sGenEtaCuts+=etaLow; sGenEtaCuts+="&&abs(mc_mu_gen_eta)<"; sGenEtaCuts+=etaUpp;
			sGenCuts+=sGenEtaCuts;
		}

		std::cout << "Generated muon cuts : " << sGenCuts << std::endl;

		treeMcZSet->Draw("mc_mu_gen_pt>>hmcZSet",sGenCuts);
		TH1F* hmcZSetc = (TH1F*)hmcZSet->Clone("hmcZSetc");
		hmcZSetc->SetFillColor(kRed);

  		double mcZGenEvents; 
  		double mcMuZGenEvents ;
		if(doEta) { 
			mcZGenEvents = hmcZSet->Integral(); //number of muons from Zs within this eta slice
			std::cout << "Generated Z muons in eta "  << " " << etaLow <<"-"<<etaUpp << " : centrality " << centralityLow << "-" << centralityUpp << mcZGenEvents <<std::endl;
			mcMuZGenEvents = mcZGenEvents;
		}
		else {
			mcZGenEvents = 0.5*hmcZSet->Integral(); //if looking over entire eta space, divide muons froms Zs by 2 to get number of Zs
			std::cout << "Generated Z muons over all eta : centrality " << centralityLow << "-" << centralityUpp << mcZGenEvents <<std::endl;
			mcMuZGenEvents = 2.0*mcZGenEvents;
		}


		hmcZSetFullc->SetFillColor(kRed);

		TString cutsPreSelZ = cuts + "&&abs(mc_mu_gen_mothertype)==23&&abs(mc_mu_gen_type)==13&&massCB>66.&&massCB<120.";
		//TString cutsPreSelZ = cuts + "&&prompt==23&&massCB>66.&&massCB<120.";
		TString cutsPreSelAllFidZ = cutsAllFid + "&&prompt==23&&massCB>66.&&massCB<120.";
		TString cutsZ = scutsMc; cutsZ+= "&&truthMatched_muid==1"; cutsZ+="&&ZDY==0";  
		std::cout << "W selection cuts: " << cutsZ << std::endl;

		//number of muons from Zs surving muon preselection cuts
		//within this eta slice and centrality class
		treeMcZSet->Draw("massCB>>hmcZSetPreSelCut",cutsPreSelZ,"hf");
		//treeMcZSet->Draw("massCB>>hmcZSetPreSelAllFidCut",cutsPreSelAllFidZ,"hf");
		std::cout << "Full fiducial cuts = " << sFullFid << std::endl;
		//number of muons from Zs over entire eta space
		treeMcZSet->Draw("mc_mu_gen_pt>>hmcZSetPreSelAllFidCut",sFullFid,"hf");

		double bkgZpreSel = hmcZSetPreSelCut->Integral();
		//total number of generated Zs
		double bkgZpreSelAllFid = 0.5*hmcZSetPreSelAllFidCut->Integral();
		std::cout << "Percentage of Z bosons in the eta region of interest = " << etaLow << "-" << etaUpp << " " << bkgZpreSel/bkgZpreSelAllFid*100.0 << "%" << std::endl; 
		std::cout << "estimated percentage of Z muons surviving in centrality " << centralityLow << "-" << centralityUpp << 
			" BEFORE W selection = " << bkgZpreSel/mcZGenEvents*100 << "%" << std::endl; 

		//number of Z muons surviving W cuts (i.e. 1-legged muons)
		//in this eta slice and centrality class
		treeMcZSet->Draw("mc_mu_gen_pt>>hmcZSetCut",cutsZ,"hf");		
		//cross check
		//treeMcZSet->Draw("pt>>hmcZSetCut",cutsZ,"hf");		

		TH1F* hmcZSetCutc = (TH1F*)hmcZSetCut->Clone("hmcZSetCutc");
  		//double mcZSurvEvents = hmcZSetCutc->Integral();
  		double mcZSurvEvents = hmcZSetCutc->Integral();
		std::cout << "Surviving Z muons after W selections in eta "  << " " << etaLow <<"-"<<etaUpp << " : centrality " << centralityLow << "-" << centralityUpp 
			<< mcZSurvEvents <<std::endl;

		hmcZSetCutc->SetFillColor(kYellow);
  	
		std::cout << "reco Z events in centrality " << centralityLow << "-" << centralityUpp << " = " << mcZGenEvents << std::endl; 
		std::cout << "Z events (one-legged mu) surviving W selection cuts in centrality " << centralityLow << "-" << centralityUpp << " = " << mcZSurvEvents << std::endl; 
		std::cout << "estimated percentage of Zs surviving in centrality " << centralityLow << "-" << centralityUpp << " = " << mcZSurvEvents/mcZGenEvents*100 << "%" << std::endl; 

		double e1 = TMath::Sqrt(mcZSurvEvents)/mcZSurvEvents*100.; double e2 = TMath::Sqrt(mcZGenEvents)/mcZGenEvents*100.;
		double eFrac = 0.01*TMath::Sqrt(e1*e1 + e2*e2);
		std::cout << mcZSurvEvents << " " << TMath::Sqrt(mcZSurvEvents) << std::endl;

		grFracCent->SetPoint(icent,ncoll,mcZSurvEvents/mcZGenEvents);
		grFracCent->SetPointError(icent,0.0,eFrac*mcZSurvEvents/mcZGenEvents);

		double xpt = (etaUpp-etaLow)/2.0+etaLow;
		grFracEta->SetPoint(ieta,xpt,mcZSurvEvents/mcZGenEvents);
		grFracEta->SetPointError(ieta,0.0,eFrac*mcZSurvEvents/mcZGenEvents);

		grFracEtaCent->SetPoint(index,ncoll,xpt, mcZSurvEvents/mcZGenEvents);
		grFracEtaCent->SetPointError(index,0.0,0.0, eFrac*mcZSurvEvents/mcZGenEvents);

		TLegend* leg = new TLegend(0.648, 0.722, 0.918, 0.915);
		leg->SetTextFont(gStyle->GetTextFont());
		leg->SetTextSize(gStyle->GetTextSize());
		leg->SetBorderSize(0);
		leg->SetFillColor(0);

		leg->AddEntry(hmcZSetc, "generated", "f");
		leg->AddEntry(hmcZSetCutc, "survivors", "f");

		TCanvas* cdatapt = new TCanvas("cdatapt","datamt",600,600);
		hmcZSetc->GetXaxis()->SetTitle("p_{T}[GeV]"); 
		hmcZSetc->GetYaxis()->SetTitle("Muons/2.5 GeV"); 
		cdatapt->SetLogy(); cdatapt->Update();
		hmcZSetc->Draw("hf");
		hmcZSetCutc->Draw("hfsame");
		hmcZSetc->Draw("sameaxis");
		
		leg->Draw();

		TLatex l;
		l.SetNDC();
		l.DrawLatex(0.53,0.37,sSel + "%" );
		l.DrawLatex(0.57,0.43,sSel2);

		myText(0.53,0.55, (Color_t)kBlack, (char*)("PYTHIA+HIJING"));

		TString plotNameLog = "ptEWBkg_"; plotNameLog+=sSel; plotNameLog+="_"; plotNameLog+=sSel2; plotNameLog+="_"; plotNameLog+="Log"; //.png";

		cdatapt->Print(plotNameLog.ReplaceAll("|",",")+".pdf"); 
	
	delete hmcZSetFullc;
	delete hmcZSetCutc;
	delete hmcZSetFull;
	delete hmcZSetCut;
	//delete cdatapt;	
	delete leg;
	delete treeMcZSet ;
	std::cout << "Plotting ends here." << std::endl;
}


void bkgZPlotter()
{
		
	bool doCentrality = true ;
	bool doEta = true ;
	bool doCharge = false ;

	float mtmax = 300.0;
	float ptmax = 300.0;

	//HIJING overlay
	//TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
	TString baseString = "/tmp/tbalestr/";
	//TString fileNameMCZIn = baseString+"MonteCarloFiles/Zmumu/HISingleMuonMC_PYTHIA_HIJING_Zmumu_11.28.2012";
	//TString fileNameMCZIn = baseString+"MonteCarloFiles/HISingleMuonMCZmumu.12.30.2012";
	TString fileNameMCZIn = baseString+"HISingleMuonMCZmumu.12.30.2012";
	TString fileNameOut;

	if(doEta&&doCentrality&&doCharge) fileNameOut = "fractionZEtaChargeCent";
	else if(doEta) fileNameOut = "fractionZEta";
	else if(doCentrality) fileNameOut = "fractionZCent";
	else fileNameOut = "fractionZ";
	
	/// --- Open output file ---
	TDirectory *dir = gDirectory;
	TFile *outFile = new TFile(fileNameOut+".root","RECREATE");
	TFile* 
	gDirectory = dir;

	TFile* fMcZSet = new TFile(fileNameMCZIn+".root", "READ");
	std::cout << "Reading input file " << fileNameMCZIn << std::endl;

	if ( !fMcZSet->IsOpen() ) {
	    std::cout << fMcZSet << " not found!" << std::endl;
	    exit(0);
	}

	gROOT->LoadMacro("AtlasUtils.C");


	// --- Set pt and eta bins ---
	std::vector<double> ptBins;
	//systematics
	double sigmaPt = 8.8;
	std::cout << "Doing muon pt systematics -1 sigma" << std::endl;
	ptBins.push_back(25.0-sigmaPt);
	ptBins.push_back(ptmax);
	const int nPtBins = ptBins.size()-1;

	std::vector<double> etaBins;
	etaBins.push_back(0.00);
	if (doEta) {
	etaBins.push_back(+0.25);
	etaBins.push_back(+0.50);
	etaBins.push_back(+0.75);
	etaBins.push_back(+1.00);
	etaBins.push_back(+1.25);
	etaBins.push_back(+1.50);
	etaBins.push_back(+1.75);
	etaBins.push_back(+2.00);
	etaBins.push_back(+2.25);

	}
	etaBins.push_back(+2.50);

	const int nEtaBins = etaBins.size()-1;
	std::vector<double> centralityBins;
	std::vector <float> ncoll;

	centralityBins.push_back(0.0);
	if(doCentrality){
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

	//	ncoll.push_back(1212.0);//0-20
		ncoll.push_back(440.6); //20-40
		ncoll.push_back(77.8); //40-80
	} else  ncoll.push_back(361.6);//0-80

	centralityBins.push_back(0.8);

	const int nCentralityBins = centralityBins.size()-1;

	TGraphErrors* grFracCent = new TGraphErrors(nCentralityBins);
	TGraphErrors* grFracEta = new TGraphErrors(nEtaBins);
	TGraphErrors* grFracCentPlus = new TGraphErrors(nCentralityBins);
	TGraphErrors* grFracEtaPlus = new TGraphErrors(nEtaBins);
	TGraphErrors* grFracCentMinus = new TGraphErrors(nCentralityBins);
	TGraphErrors* grFracEtaMinus = new TGraphErrors(nEtaBins);

	//plot charge difference for systematics
	TGraph* grDiffEta = new TGraphErrors(nEtaBins);
	TGraph* grDiffNcoll  = new TGraphErrors(nCentralityBins);

	TGraph2DErrors* grFracEtaCent = new TGraph2DErrors(nCentralityBins*nEtaBins);
	TGraph2DErrors* grFracEtaCentPlus = new TGraph2DErrors(nCentralityBins*nEtaBins);
	TGraph2DErrors* grFracEtaCentMinus = new TGraph2DErrors(nCentralityBins*nEtaBins);

	//base muon selection cuts
	TString cuts ="abs(charge)==1";
        cuts+="&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0";
	TString cutsP = cuts+"&&charge==1";
	TString cutsM = cuts+"&&charge==-1";
	//std::cout << cuts << std::endl;
	for ( int i = 0; i < nPtBins; i++ ) {
	  for ( int j = 0; j < nEtaBins; j++ ) {
	    for ( int k = 0; k < nCentralityBins; k++ ){

		int index = j*nCentralityBins + k;
		std::cout << " plotting "<<i<<":"<<j<<":"<<k<<std::endl;
		TString sCentLow = "";
		TString sCentUp = "";
		TString sEtaLow = "";
		TString sEtaUp = "";

		sCentLow += 100*centralityBins[k]; //sCentLow.Remove(3);
		sCentUp += 100*centralityBins[k+1]; //sCentUp.Remove(3);

		sEtaLow += etaBins[j];
		sEtaUp += etaBins[j+1];

	    	TString sSel = "Z#rightarrow#mu^{+}#mu^{-} ";
	    	TString sSelEta = sEtaLow;
	    	sSel += sCentLow; sSel+="-"; sSel+= sCentUp;
	    	sSelEta += "#leq"; sSelEta+= "|#eta|"; sSelEta+="<"; sSelEta += sEtaUp;

		std::cout << "plotting for : " << i << ":" << j << ":" << k <<std::endl;
		plot(fMcZSet, grFracEtaCent, grFracCent, grFracEta, mtmax, cuts, "&&abs(mc_mu_charge)==1", index, ptBins[i], ptBins[i+1], etaBins[j], 
			etaBins[j+1], centralityBins[k], centralityBins[k+1], ncoll[k],j,k,i, sSel , sSelEta, doEta);

		if(doCharge){
		    std::cout << "#mu^{+}" << std::endl;
		    plot(fMcZSet,grFracEtaCentPlus, grFracCentPlus, grFracEtaPlus, mtmax, cutsP,"&&mc_mu_charge==+1", index, ptBins[i], ptBins[i+1], etaBins[j], 
			etaBins[j+1], centralityBins[k], centralityBins[k+1], ncoll[k],j,k,102, "#mu^{+},"+sSel , sSelEta, doEta);

		    std::cout << "#mu^{-}" << std::endl;
		    plot(fMcZSet, grFracEtaCentMinus, grFracCentMinus, grFracEtaMinus, mtmax, cutsM,"&&mc_mu_charge==-1", index, ptBins[i], ptBins[i+1], etaBins[j], 
			etaBins[j+1], centralityBins[k], centralityBins[k+1], ncoll[k],j,k,103, "#mu^{-},"+sSel , sSelEta, doEta);

			//systematics
	    	        std::cout << "|mu^{-}-mu^{+}|" << std::endl;
			if(doEta){
			  std::cout << "|eta|" << std::endl;
			  plotChargeDiffEta(grDiffEta,grFracEtaPlus,grFracEtaMinus,j,k,etaBins[j], etaBins[j+1]);
			}

			if(doCentrality){
			  std::cout << "<Ncoll>" << std::endl;
			  plotChargeDiffNcoll(grDiffNcoll,grFracCentPlus,grFracCentMinus,j,k,ncoll[k]);
			}

		  }

		} //icent

	     } //ieta
	}
	
	if(doCharge){
 		TString sChargePlus = "Z,#mu^{+}";
		TString sChargeMinus = "Z,#mu^{-}";
		plotFraction(grFracCentPlus,grFracEtaPlus,sChargePlus);
		plotFraction(grFracCentMinus,grFracEtaMinus,sChargeMinus);
		Write(outFile, grFracCentPlus, "fractionZCentPlus");
		Write(outFile, grFracCentMinus, "fractionZCentMinus");
		Write(outFile, grFracEtaPlus, "fractionZEtaPlus");
		Write(outFile, grFracEtaMinus, "fractionZEtaMinus");
		Write(outFile, grFracEtaCentPlus, "fractionZEtaCentPlus");
		Write(outFile, grFracEtaCentMinus, "fractionZEtaCentMinus");
		//systematics
		if(doEta) Write(outFile,grDiffEta,"diffFracZEta"); 
		if(doCentrality) {
			Write(outFile,grDiffNcoll,"diffFracZNcoll");
		}

	}

		TString sCharge = "Z,#mu^{pm}";
		plotFraction(grFracCent,grFracEta,sCharge);
		Write(outFile, grFracCent, "fractionZCent");
		Write(outFile, grFracEta, "fractionZEta");
		Write(outFile, grFracEtaCent, "fractionZEtaCent");

/*	TLegend* leg = new TLegend(0.648, 0.622, 0.918, 0.815);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->AddEntry(grFracMuP,"#mu^{+}","p");
	leg->AddEntry(grFracMuM,"#mu^{-}","p");

	leg->Draw(); 
*/	
}
