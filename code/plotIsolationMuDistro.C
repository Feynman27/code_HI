void Write(TFile* outFile, TH1F* h, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing histo to root file..." << std::endl;
	  h->Write(sFile);
	  std::cout << "Done." << std::endl;
    	  gDirectory = dir;
}

void plotIsolationMuDistro(){

	SetAtlasStyle();
	TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
	
	/// --- Open output file ---
	TDirectory *dir = gDirectory;
	TFile *outFile = new TFile("isoloationMuHistos.root","RECREATE");
	gDirectory = dir;

	//data overlay
	//TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuMCDataOverlay.2012.11.21";
	//TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.2013.03.06";
	TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.03.13.2013";
	TFile* fMcWSet = new TFile(fileNameMCWIn+".root", "READ");

	if ( !fMcWSet->IsOpen() ) {
	    std::cout << fMcWSet << " not found!" << std::endl;
	    exit(0);
	}

  	// --- W set ---

	//TString fileNameMCJ1In = baseString+"MonteCarloFiles/QCD/HISingleMuonMcJ1.2012.11.29";
	///new file with custom isolation cuts
	TString fileNameMCJ1In = baseString+"MonteCarloFiles/QCD/HISingleMuonJ1_1muPYTHIADataOverlay.03.13.2013";
	TFile* fMcJ1Set = new TFile(fileNameMCJ1In+".root", "READ");
	//J2 1 muon-filter 
	//TString fileNameMCJ2In = baseString+"MonteCarloFiles/QCD/HISingleMuonMcJ2.2012.11.29";
	TString fileNameMCJ2In = baseString+"MonteCarloFiles/QCD/HISingleMuonJ2_1muPYTHIADataOverlay.03.13.2013";
	TFile* fMcJ2Set = new TFile(fileNameMCJ2In+".root", "READ");
	//TString fileNameMCJ3In = baseString+"MonteCarloFiles/QCD/HISingleMuonMcJ3.2012.11.29";
	TString fileNameMCJ3In = baseString+"MonteCarloFiles/QCD/HISingleMuonJ3_1muPYTHIADataOverlay.03.13.2013";
	TFile* fMcJ3Set = new TFile(fileNameMCJ3In+".root", "READ");

	TString cuts ="abs(charge)==1";
	cuts+="&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&ZDY==0";

	double ptLow=4.0 ;double ptUpp=100.0 ;double etaLow=0.0; double etaUpp=2.5; double centralityLow=0.0;double centralityUpp=0.8;
	int nbins = 100;
	double ptCtrlLo = 10.0; double ptCtrlUpp = 20.0; double ptSigLo = 25.0; 
	//define bins for control region
	int ptCtrlBinLo = (nbins/ptUpp)*ptCtrlLo+1; int ptCtrlBinUpp = (nbins/ptUpp)*ptCtrlUpp+1; 
	//define bins for signal region
	int ptSigBinLo = (nbins/ptUpp)*ptSigLo+1; int ptSigBinUpp = (nbins/ptUpp)*ptUpp+1;  

	double mtcutLow = 40.0; float mtmax = 400.0;

	cuts+="&&pt>"; cuts += ptLow; cuts+="&&pt<"; cuts += ptUpp;
	cuts+="&&abs(eta)>"; cuts += etaLow; cuts+="&&abs(eta)<"; cuts += etaUpp;

	double missPtCut = 25.0; 
	double missPtMax = 9000.0; 

	TString sMtCutLo = "&&mt>"; sMtCutLo+=mtcutLow; TString sMtCutUp = "&&mt<";sMtCutUp+=mtmax;
	TString sMissPtLo = "&&nu_pt<"; sMissPtLo+=missPtMax; TString sMissPtUp = "&&nu_pt>"; sMissPtUp+=missPtCut;

	TTree *treeMcWSet  = (TTree*)fMcWSet->Get("tree");
	TTree *treeMcJ1Set = (TTree*)fMcJ1Set->Get("tree");
	TTree *treeMcJ2Set = (TTree*)fMcJ2Set->Get("tree");
	TTree *treeMcJ3Set = (TTree*)fMcJ3Set->Get("tree");

	//centrality
	std::vector <float> centBins;
	centBins.push_back(0.0);
//	centBins.push_back(0.05);
	centBins.push_back(0.1);
//	centBins.push_back(0.15);
	centBins.push_back(0.2);
	centBins.push_back(0.4);

	centBins.push_back(0.8);

	const int centralityBins = centBins.size()-1;

	//cone radii around muon
	int cone10 = 10;
	int cone15 = 15;
	int cone20 = 20;
	int cone30 = 30;
	int cone40 = 40;
	const unsigned int nCones = 5;
	int arrCone[nCones] = {10,15,20,30,40};
	//lower track pt used in iso calculation
	const unsigned int nTrkPtCut = 5;
	int trkPtLo[nTrkPtCut] = {2,3,4,5,6};

	//double isoCut[isoCutBins] = {0.2,0.3};
	double arrCentWidth = 0.8;
	double ncoll = 361.6;
	double evData = 1.0e9; //number of hp events
	double scaleFactor = arrCentWidth*ncoll*evData;

	TH1F* hmcJ1 = new TH1F("hmcJ1","hmcJ1",nbins,0.0,ptUpp);
	TH1F* hmcJ2 = new TH1F("hmcJ2","hmcJ2",nbins,0.0,ptUpp);
	TH1F* hmcJ3 = new TH1F("hmcJ3","hmcJ3",nbins,0.0,ptUpp);

	TH1F* hmcJ1Iso = new TH1F("hmcJ1Iso","hmcJ1Iso",nbins,0.0,ptUpp);
	TH1F* hmcJ2Iso = new TH1F("hmcJ2Iso","hmcJ2Iso",nbins,0.0,ptUpp);
	TH1F* hmcJ3Iso = new TH1F("hmcJ3Iso","hmcJ3Iso",nbins,0.0,ptUpp);

	TH1F* hWIsolationMu[nTrkPtCut][nCones] = new TH1F("hWIsolationMu","hWIsolationMu",nbins,0.0,ptUpp);
	TH1F* hJetIsolationMu[nTrkPtCut][nCones] = new TH1F("hJetIsolationMu","hJetIsolationMu",nbins,0.0,ptUpp);

	TH1F* hmcJ1Cent = new TH1F("hmcJ1Cent","hmcJ1Cent",nbins,0.0,ptUpp);
	TH1F* hmcJ2Cent = new TH1F("hmcJ2Cent","hmcJ2Cent",nbins,0.0,ptUpp);
	TH1F* hmcJ3Cent = new TH1F("hmcJ3Cent","hmcJ3Cent",nbins,0.0,ptUpp);

	//graphs of the efficiency as fcn of cone radius round muon
	const int coneRadii = 5;
	int arrConeRad[coneRadii] = {10,15,20,30,40};
	TGraph* _grBkg = new TGraph(coneRadii);
	TGraph* _grSig = new TGraph(coneRadii);

	double wtJ1 = 6.6523e-03*1.8770e-4/64.0e-3; double wtJ2 = 1.4941e-02*8.2788e-6/64.0e-3; double wtJ3 = 2.4284e-02*2.9426e-7/64.0e-3;

	char sJetIsolationMu[200];
	char sWIsolationMu[200];
	for(int icent=0; icent<centralityBins; icent++){
	
		TString sAddCent = "&&centrality>"; sAddCent += centBins.at(icent); sAddCent += "&&centrality<"; sAddCent += centBins.at(icent+1);
		TString selCuts = cuts; selCuts += sMissPtUp; selCuts += sMissPtLo;selCuts +=sMtCutLo; selCuts +=sMtCutUp; selCuts += sAddCent; 
		std::cout << "Selection cuts w/o isolation: " << selCuts << std::endl;
		
		TString sCentrality = "centrality>"; sCentrality+= centBins.at(icent); sCentrality+="&&centrality<"; sCentrality+=centBins.at(icent+1);
		//event counting histos
		treeMcJ1Set->Draw("centrality>>hmcJ1Cent",sCentrality,"hf");
		treeMcJ2Set->Draw("centrality>>hmcJ2Cent",sCentrality,"hf");
		treeMcJ3Set->Draw("centrality>>hmcJ3Cent",sCentrality,"hf");

		double nMcJ1 = hmcJ1Cent->Integral(); double nMcJ2 = hmcJ2Cent->Integral(); double nMcJ3 = hmcJ3Cent->Integral();

	   for(int itrk=0; itrk<nTrkPtCut; itrk++){

			TString sTitle = "p_{T}^{trk}>"; sTitle+=trkPtLo[itrk]; sTitle+="GeV";
		for(int icone=0; icone<nCones; icone++){

			sprintf(sJetIsolationMu,"hJetIsolationMuTrkPt%iCone%i",itrk,icone);
			sprintf(sWIsolationMu,"hWIsolationMuTrkPt%iCone%i",itrk,icone);

			hJetIsolationMu[itrk][icone] = new TH1F(sJetIsolationMu,sJetIsolationMu,nbins,0.0,ptUpp);
			hWIsolationMu[itrk][icone] = new TH1F(sWIsolationMu,sWIsolationMu,nbins,0.0,ptUpp);

			TString sPtCone = "ptcone"; sPtCone+=arrCone[icone]; sPtCone+="ID"; sPtCone+=trkPtLo[itrk];
			std::cout << "cone variable: " << sPtCone << std::endl;
			treeMcJ1Set->Draw(sPtCone+">>hmcJ1Iso",selCuts,"hf");
			treeMcJ2Set->Draw(sPtCone+">>hmcJ2Iso",selCuts,"hf");
			treeMcJ3Set->Draw(sPtCone+">>hmcJ3Iso",selCuts,"hf");
			//weighted muon yield per Jx mc event
			hJetIsolationMu[itrk][icone]->Add(hmcJ1Iso,hmcJ2Iso,(wtJ1)/(nMcJ1)*scaleFactor, (wtJ2)/(nMcJ2)*scaleFactor); 
			hJetIsolationMu[itrk][icone]->Add(hJetIsolationMu[itrk][icone],hmcJ3Iso,1.0,(wtJ3)/(nMcJ3)*scaleFactor);

			treeMcWSet->Draw((TString(sPtCone+">>")+sWIsolationMu).Data(),selCuts+"&&prompt==24","hf");		

			//save histos
			
			Write(outFile,hJetIsolationMu[itrk][icone],sJetIsolationMu );
			Write(outFile,hWIsolationMu[itrk][icone],sWIsolationMu );

		}//icone

       }//lower track pt	
    }//icent

 }
