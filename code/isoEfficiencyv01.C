void Write(TFile* outFile, TObject* grIso, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing TGraph to root file..." << std::endl;
	  grIso->Write(sFile);
	  std::cout << "Done." << std::endl;
    	  gDirectory = dir;
}

void isoEfficiencyv01(){

	SetAtlasStyle();
	TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
	
	/// --- Open output file ---
	TDirectory *dir = gDirectory;
	TFile *outFile = new TFile("isoEfficiencyGraphsv01.root","RECREATE");
	gDirectory = dir;

	//data overlay
	//TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuMCDataOverlay.2012.11.21";
	//TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.2013.03.06";
	TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.03.24.2013";
	TFile* fMcWSet = new TFile(fileNameMCWIn+".root", "READ");

	if ( !fMcWSet->IsOpen() ) {
	    std::cout << fMcWSet << " not found!" << std::endl;
	    exit(0);
	}

  	// --- W set ---

	//TString fileNameMCJ1In = baseString+"MonteCarloFiles/QCD/HISingleMuonMcJ1.2012.11.29";
	///new file with custom isolation cuts
	TString fileNameMCJ1In = baseString+"MonteCarloFiles/QCD/HISingleMuonJ1_1muPYTHIADataOverlay.03.24.2013";
	TFile* fMcJ1Set = new TFile(fileNameMCJ1In+".root", "READ");
	//J2 1 muon-filter 
	//TString fileNameMCJ2In = baseString+"MonteCarloFiles/QCD/HISingleMuonMcJ2.2012.11.29";
	TString fileNameMCJ2In = baseString+"MonteCarloFiles/QCD/HISingleMuonJ2_1muPYTHIADataOverlay.03.24.2013";
	TFile* fMcJ2Set = new TFile(fileNameMCJ2In+".root", "READ");
	//TString fileNameMCJ3In = baseString+"MonteCarloFiles/QCD/HISingleMuonMcJ3.2012.11.29";
	TString fileNameMCJ3In = baseString+"MonteCarloFiles/QCD/HISingleMuonJ3_1muPYTHIADataOverlay.03.24.2013";
	TFile* fMcJ3Set = new TFile(fileNameMCJ3In+".root", "READ");

	TString cuts ="abs(charge)==1";
	cuts+="&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&ZDY==0";

	double ptLow=4.0 ;double ptUpp=400.0 ;double etaLow=0.0; double etaUpp=2.5; double centralityLow=0.0;double centralityUpp=0.8;
	int nbins = 400;
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
//	centBins.push_back(0.2);
	centBins.push_back(0.4);

	centBins.push_back(0.8);

	const int centralityBins = centBins.size()-1;

/*    std::vector<float> cones;
    cones.push_back(0.1);
    cones.push_back(0.15);
    cones.push_back(0.2);
    cones.push_back(0.3);
    cones.push_back(0.4);
    cones.push_back(0.5);
    const int nCones = cones.size();
*/
    //sig and bkg per centrality bin
    double ns[centralityBins] = {2522.0,3749.0,713.0};
    double nb[centralityBins] = {228.0,275.0,72.0};

	//cone radii around muon
	int cone10 = 10;
	int cone15 = 15;
	int cone20 = 20;
	int cone30 = 30;
	int cone40 = 40;
	int cone50 = 50;
	//lower track pt used in iso calculation
	const int trkPtCut = 8;
	TString trkPtLo[trkPtCut] = {"05","075","1","2","3","4","5","6"};
	//int trkPtLo[2] = {2,3};
	//const int trkPtCut = 2;

	//const int isoCutBins = 2;
	const int isoCutBins = 11;
	double isoCut[isoCutBins] = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 5.0, 10.0, 100.0};
	//double isoCut[isoCutBins] = {0.2,0.3};
	double scaleFactor = 1.0e9;

	TH1F* hmcJ1 = new TH1F("hmcJ1","hmcJ1",nbins,0.0,ptUpp);
	TH1F* hmcJ2 = new TH1F("hmcJ2","hmcJ2",nbins,0.0,ptUpp);
	TH1F* hmcJ3 = new TH1F("hmcJ3","hmcJ3",nbins,0.0,ptUpp);

	TH1F* hmcJ1Iso = new TH1F("hmcJ1Iso","hmcJ1Iso",nbins,0.0,ptUpp);
	TH1F* hmcJ2Iso = new TH1F("hmcJ2Iso","hmcJ2Iso",nbins,0.0,ptUpp);
	TH1F* hmcJ3Iso = new TH1F("hmcJ3Iso","hmcJ3Iso",nbins,0.0,ptUpp);

	TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",nbins,0.0,ptUpp);
	TH1F* hmcQCDSetIso10 = new TH1F("hmcQCDSetIso10","hmcQCDSetIso10",nbins,0.0,ptUpp);
	TH1F* hmcQCDSetIso15 = new TH1F("hmcQCDSetIso15","hmcQCDSetIso15",nbins,0.0,ptUpp);
	TH1F* hmcQCDSetIso20 = new TH1F("hmcQCDSetIso20","hmcQCDSetIso20",nbins,0.0,ptUpp);
	TH1F* hmcQCDSetIso30 = new TH1F("hmcQCDSetIso30","hmcQCDSetIso30",nbins,0.0,ptUpp);
	TH1F* hmcQCDSetIso40 = new TH1F("hmcQCDSetIso40","hmcQCDSetIso40",nbins,0.0,ptUpp);
	TH1F* hmcQCDSetIso50 = new TH1F("hmcQCDSetIso50","hmcQCDSetIso50",nbins,0.0,ptUpp);

	TH1F* hmcWSet = new TH1F("hmcWSet","hmcWSet",nbins,0.0,ptUpp);
	TH1F* hmcWSetIso = new TH1F("hmcWSetIso","hmcWSetIso",nbins,0.0,ptUpp);

	TH1F* hmcJ1Cent = new TH1F("hmcJ1Cent","hmcJ1Cent",nbins,0.0,ptUpp);
	TH1F* hmcJ2Cent = new TH1F("hmcJ2Cent","hmcJ2Cent",nbins,0.0,ptUpp);
	TH1F* hmcJ3Cent = new TH1F("hmcJ3Cent","hmcJ3Cent",nbins,0.0,ptUpp);


	
	//graphs of the efficiency as fcn of cone radius round muon
	const int coneRadii = 6;
	int arrConeRad[coneRadii] = {10,15,20,30,40,50};
	TGraph* _grBkg = new TGraph(coneRadii);
	TGraph* _grSig = new TGraph(coneRadii);

	double wtJ1 = 6.6523e-03*1.8770e-4/64.0e-3; double wtJ2 = 1.4941e-02*8.2788e-6/64.0e-3; double wtJ3 = 2.4284e-02*2.9426e-7/64.0e-3;

	for(int icent=0; icent<centralityBins; icent++){
	
        TList _grEffRad10; 
        TList _grEffRad15; 
        TList _grEffRad20; 
        TList _grEffRad30; 
        TList _grEffRad40; 
        TList _grEffRad50; 

        TList _grRatioRad10; 
        TList _grRatioRad15; 
        TList _grRatioRad20; 
        TList _grRatioRad30; 
        TList _grRatioRad40; 
        TList _grRatioRad50; 
        
		TString sAddCent = "&&centrality>"; sAddCent += centBins.at(icent); sAddCent += "&&centrality<"; sAddCent += centBins.at(icent+1);
		TString selCuts = cuts; selCuts += sMissPtUp; selCuts += sMissPtLo;selCuts +=sMtCutLo; selCuts +=sMtCutUp; selCuts += sAddCent; 
		std::cout << "Selection cuts w/o isolation: " << selCuts << std::endl;
		
		TString sCentrality = "centrality>"; sCentrality+= centBins.at(icent); sCentrality+="&&centrality<"; sCentrality+=centBins.at(icent+1);
		//event counting histos
		treeMcJ1Set->Draw("centrality>>hmcJ1Cent",sCentrality,"hf");
		treeMcJ2Set->Draw("centrality>>hmcJ2Cent",sCentrality,"hf");
		treeMcJ3Set->Draw("centrality>>hmcJ3Cent",sCentrality,"hf");

	   for(int itrk=0; itrk<trkPtCut; itrk++){

			TString sTitle = "p_{T}^{trk}>"; sTitle+=trkPtLo[itrk]; sTitle+="GeV";
			TString sGr = "grIsoEffTrkPt"; sGr+=trkPtLo[itrk];  sGr+="cent"; sGr+=icent; 
			TString sGr10 = sGr + "ConeRadius10";  
			TString sGr15 = sGr + "ConeRadius15";  
			TString sGr20 = sGr + "ConeRadius20";  
			TString sGr30 = sGr + "ConeRadius30";  
			TString sGr40 = sGr + "ConeRadius40";  
			TString sGr50 = sGr + "ConeRadius50";  

            TString sGrRatio10 = sGr10 + "SigBkgRatio";  
            TString sGrRatio15 = sGr15 + "SigBkgRatio";  
            TString sGrRatio20 = sGr20 + "SigBkgRatio";  
            TString sGrRatio30 = sGr30 + "SigBkgRatio";  
            TString sGrRatio40 = sGr40 + "SigBkgRatio";  
            TString sGrRatio50 = sGr50 + "SigBkgRatio";  


			_grEffRad10.Add( new TGraph(isoCutBins));
			_grEffRad15.Add( new TGraph(isoCutBins));
			_grEffRad20.Add( new TGraph(isoCutBins));
			_grEffRad30.Add( new TGraph(isoCutBins));
			_grEffRad40.Add( new TGraph(isoCutBins));
			_grEffRad50.Add( new TGraph(isoCutBins));

			_grRatioRad10.Add( new TGraph(isoCutBins));
			_grRatioRad15.Add( new TGraph(isoCutBins));
			_grRatioRad20.Add( new TGraph(isoCutBins));
			_grRatioRad30.Add( new TGraph(isoCutBins));
			_grRatioRad40.Add( new TGraph(isoCutBins));
			_grRatioRad50.Add( new TGraph(isoCutBins));
	
	
		for(int icut=0; icut<isoCutBins; icut++){
			
		
			TString sCone = "&&ptcone"; 
			TString sCone10 = sCone; sCone10+=cone10; sCone10+="ID"; sCone10+=trkPtLo[itrk]; sCone10+="/pt<"; 
			TString sCone15 = sCone; sCone15+=cone15; sCone15+="ID"; sCone15+=trkPtLo[itrk]; sCone15+="/pt<"; 
			TString sCone20 = sCone; sCone20+=cone20; sCone20+="ID"; sCone20+=trkPtLo[itrk]; sCone20+="/pt<"; 
			TString sCone30 = sCone; sCone30+=cone30; sCone30+="ID"; sCone30+=trkPtLo[itrk]; sCone30+="/pt<"; 
			TString sCone40 = sCone; sCone40+=cone40; sCone40+="ID"; sCone40+=trkPtLo[itrk]; sCone40+="/pt<"; 
			TString sCone50 = sCone; sCone50+=cone50; sCone50+="ID"; sCone50+=trkPtLo[itrk]; sCone50+="/pt<"; 


			sCone10+=isoCut[icut];
			sCone15+=isoCut[icut];
			sCone20+=isoCut[icut];
			sCone30+=isoCut[icut];
			sCone40+=isoCut[icut];
			sCone50+=isoCut[icut];

			TString selCuts10 = selCuts + sCone10;
			TString selCuts15 = selCuts + sCone15;
			TString selCuts20 = selCuts + sCone20;
			TString selCuts30 = selCuts + sCone30;
			TString selCuts40 = selCuts + sCone40;
			TString selCuts50 = selCuts + sCone50;
			
			std::cout << "isolation cut for dR 0.1: " << selCuts10 << std::endl;
			std::cout << "isolation cut for dR 0.15: " << selCuts15 << std::endl;
			std::cout << "isolation cut for dR 0.2: " << selCuts20 << std::endl;
			std::cout << "isolation cut for dR 0.3: " << selCuts30 << std::endl;
			std::cout << "isolation cut for dR 0.4: " << selCuts40 << std::endl;
			std::cout << "isolation cut for dR 0.5: " << selCuts50 << std::endl;

			double nMcJ1 = hmcJ1Cent->Integral(); double nMcJ2 = hmcJ2Cent->Integral(); double nMcJ3 = hmcJ3Cent->Integral();

			treeMcJ1Set->Draw("pt>>hmcJ1",selCuts,"hf");
			treeMcJ2Set->Draw("pt>>hmcJ2",selCuts,"hf");
			treeMcJ3Set->Draw("pt>>hmcJ3",selCuts,"hf");
			//weighted muon yield per Jx mc event
			hmcQCDSet->Add(hmcJ1,hmcJ2,(wtJ1)/(nMcJ1)*scaleFactor, (wtJ2)/(nMcJ2)*scaleFactor); 
			hmcQCDSet->Add(hmcQCDSet,hmcJ3,1.0,(wtJ3)/(nMcJ3)*scaleFactor);

			std::cout << "W cuts no iso: " << selCuts << std::endl;
			treeMcWSet->Draw("pt>>hmcWSet",selCuts+"&&prompt==24","hf");		


			treeMcJ1Set->Draw("pt>>hmcJ1Iso",selCuts10,"hf");
			treeMcJ2Set->Draw("pt>>hmcJ2Iso",selCuts10,"hf");
			treeMcJ3Set->Draw("pt>>hmcJ3Iso",selCuts10,"hf");
			//weighted muon yield per Jx mc event
			hmcQCDSetIso10->Add(hmcJ1Iso,hmcJ2Iso,(wtJ1)/(nMcJ1)*scaleFactor, (wtJ2)/(nMcJ2)*scaleFactor); 
			hmcQCDSetIso10->Add(hmcQCDSetIso10,hmcJ3Iso,1.0,(wtJ3)/(nMcJ3)*scaleFactor);

			std::cout << "W cuts after iso: " << selCuts10 << std::endl;
			treeMcWSet->Draw("pt>>hmcWSetIso",selCuts10+"&&prompt==24","hf");		

			double bkgEff10 = hmcQCDSetIso10->Integral(ptSigBinLo,ptSigBinUpp)/hmcQCDSet->Integral(ptSigBinLo,ptSigBinUpp);
			double sigEff10 = hmcWSetIso->Integral(ptSigBinLo,ptSigBinUpp)/hmcWSet->Integral(ptSigBinLo,ptSigBinUpp);

			((TGraph*)_grEffRad10.At(itrk))->SetPoint(icut,sigEff10,bkgEff10);
			std::cout << "dR 0.1 sig fraction at point " << icut << " = " << sigEff10 << "*" << ns[icent] << " = " <<
                sigEff10*ns[icent] << std::endl;
			std::cout << "dR 0.1 bkg fraction at point " << icut << " = " << bkgEff10 << "*" << nb[icent] << " = " <<
                bkgEff10*nb[icent] << std::endl;

            if(sigEff10*ns[icent]>0.0){
                std::cout << "That gives a background:signal ratio of " << (bkgEff10*nb[icent])/(sigEff10*ns[icent])<< std::endl;
			    ((TGraph*)_grRatioRad10.At(itrk))->SetPoint(icut,(bkgEff10*nb[icent])/(sigEff10*ns[icent]),0.1);
            }
            else {
                std::cout <<"WARNING: no signal events in this bin. B/S set to -1.0. " << std::endl;
                ((TGraph*)_grRatioRad10.At(itrk))->SetPoint(icut,-1.0,0.1);
            }

			treeMcJ1Set->Draw("pt>>hmcJ1Iso",selCuts15,"hf");
			treeMcJ2Set->Draw("pt>>hmcJ2Iso",selCuts15,"hf");
			treeMcJ3Set->Draw("pt>>hmcJ3Iso",selCuts15,"hf");
			//weighted muon yield per Jx mc event
			hmcQCDSetIso15->Add(hmcJ1Iso,hmcJ2Iso,(wtJ1)/(nMcJ1)*scaleFactor, (wtJ2)/(nMcJ2)*scaleFactor); 
			hmcQCDSetIso15->Add(hmcQCDSetIso15,hmcJ3Iso,1.0,(wtJ3)/(nMcJ3)*scaleFactor);

			std::cout << "W cuts after iso: " << selCuts15 << std::endl;
			treeMcWSet->Draw("pt>>hmcWSetIso",selCuts15+"&&prompt==24","hf");		

			double bkgEff15 = hmcQCDSetIso15->Integral(ptSigBinLo,ptSigBinUpp)/hmcQCDSet->Integral(ptSigBinLo,ptSigBinUpp);
			double sigEff15 = hmcWSetIso->Integral(ptSigBinLo,ptSigBinUpp)/hmcWSet->Integral(ptSigBinLo,ptSigBinUpp);

			((TGraph*)_grEffRad15.At(itrk))->SetPoint(icut,sigEff15,bkgEff15);
			std::cout << "sig fraction at point " << icut << " = " << sigEff15 << "*" << ns[icent] << " = " <<
                sigEff15*ns[icent] << std::endl;
			std::cout << "bkg fraction at point " << icut << " = " << bkgEff15 << "*" << nb[icent] << " = " <<
                bkgEff15*nb[icent] << std::endl;
            std::cout <<
            if(sigEff15*ns[icent]>0.0){
                std::cout << "That gives a background:signal ratio of " << (bkgEff15*nb[icent])/(sigEff15*ns[icent])<< std::endl;
			    ((TGraph*)_grRatioRad15.At(itrk))->SetPoint(icut,(bkgEff15*nb[icent])/(sigEff15*ns[icent]),0.1);
            }
            else {
                std::cout <<"WARNING: no signal events in this bin. B/S set to -1.0. " << std::endl;
                ((TGraph*)_grRatioRad15.At(itrk))->SetPoint(icut,-1.0,0.1);
            }

			std::cout << "dR 0.15 bkg fraction at point " << icut << " = " << bkgEff15 << std::endl;
			std::cout << "dR 0.15 sig fraction at point " << icut << " = " << sigEff15 << std::endl;


			treeMcJ1Set->Draw("pt>>hmcJ1Iso",selCuts20,"hf");
			treeMcJ2Set->Draw("pt>>hmcJ2Iso",selCuts20,"hf");
			treeMcJ3Set->Draw("pt>>hmcJ3Iso",selCuts20,"hf");
			//weighted muon yield per Jx mc event
			hmcQCDSetIso20->Add(hmcJ1Iso,hmcJ2Iso,(wtJ1)/(nMcJ1)*scaleFactor, (wtJ2)/(nMcJ2)*scaleFactor); 
			hmcQCDSetIso20->Add(hmcQCDSetIso20,hmcJ3Iso,1.0,(wtJ3)/(nMcJ3)*scaleFactor);

			std::cout << "W cuts after iso: " << selCuts20 << std::endl;
			treeMcWSet->Draw("pt>>hmcWSetIso",selCuts20+"&&prompt==24","hf");		

			double bkgEff20 = hmcQCDSetIso20->Integral(ptSigBinLo,ptSigBinUpp)/hmcQCDSet->Integral(ptSigBinLo,ptSigBinUpp);
			double sigEff20 = hmcWSetIso->Integral(ptSigBinLo,ptSigBinUpp)/hmcWSet->Integral(ptSigBinLo,ptSigBinUpp);

			((TGraph*)_grEffRad20.At(itrk))->SetPoint(icut,sigEff20,bkgEff20);
			std::cout << "sig fraction at point " << icut << " = " << sigEff20 << "*" << ns[icent] << " = " <<
                sigEff20*ns[icent] << std::endl;
			std::cout << "bkg fraction at point " << icut << " = " << bkgEff20 << "*" << nb[icent] << " = " <<
                bkgEff20*nb[icent] << std::endl;
            if(sigEff20*ns[icent]>0.0){
                std::cout << "That gives a background:signal ratio of " << (bkgEff20*nb[icent])/(sigEff20*ns[icent])<< std::endl;
			    ((TGraph*)_grRatioRad20.At(itrk))->SetPoint(icut,(bkgEff20*nb[icent])/(sigEff20*ns[icent]),0.1);
            }
            else {
                std::cout <<"WARNING: no signal events in this bin. B/S set to -1.0. " << std::endl;
                ((TGraph*)_grRatioRad20.At(itrk))->SetPoint(icut,-1.0,0.1);
            }

			std::cout << "track Pt " << trkPtLo[itrk]<< " " << hmcQCDSetIso20->Integral(ptSigBinLo,ptSigBinUpp) << "/" << hmcQCDSet->Integral(ptSigBinLo,ptSigBinUpp) 
						<< " = " << bkgEff20 << std::endl;

			std::cout << "track Pt " << trkPtLo[itrk] << " " << hmcWSetIso->Integral(ptSigBinLo,ptSigBinUpp) << "/" << hmcWSet->Integral(ptSigBinLo,ptSigBinUpp) 
						<< " = " << sigEff20 << std::endl;


			std::cout << "dR 0.2 bkg fraction at point " << icut << " = " << bkgEff20 << std::endl;
			std::cout << "dR 0.2 sig fraction at point " << icut << " = " << sigEff20 << std::endl;

			treeMcJ1Set->Draw("pt>>hmcJ1Iso",selCuts30,"hf");
			treeMcJ2Set->Draw("pt>>hmcJ2Iso",selCuts30,"hf");
			treeMcJ3Set->Draw("pt>>hmcJ3Iso",selCuts30,"hf");
			//weighted muon yield per Jx mc event
			hmcQCDSetIso30->Add(hmcJ1Iso,hmcJ2Iso,(wtJ1)/(nMcJ1)*scaleFactor, (wtJ2)/(nMcJ2)*scaleFactor); 
			hmcQCDSetIso30->Add(hmcQCDSetIso30,hmcJ3Iso,1.0,(wtJ3)/(nMcJ3)*scaleFactor);

			treeMcWSet->Draw("pt>>hmcWSetIso",selCuts30+"&&prompt==24","hf");		

			double bkgEff30 = hmcQCDSetIso30->Integral(ptSigBinLo,ptSigBinUpp)/hmcQCDSet->Integral(ptSigBinLo,ptSigBinUpp);
			double sigEff30 = hmcWSetIso->Integral(ptSigBinLo,ptSigBinUpp)/hmcWSet->Integral(ptSigBinLo,ptSigBinUpp);

			((TGraph*)_grEffRad30.At(itrk))->SetPoint(icut,sigEff30,bkgEff30);
			std::cout << "sig fraction at point " << icut << " = " << sigEff30 << "*" << ns[icent] << " = " <<
                sigEff30*ns[icent] << std::endl;
			std::cout << "bkg fraction at point " << icut << " = " << bkgEff30 << "*" << nb[icent] << " = " <<
                bkgEff30*nb[icent] << std::endl;
            if(sigEff30*ns[icent]>0.0){
                std::cout << "That gives a background:signal ratio of " << (bkgEff30*nb[icent])/(sigEff30*ns[icent])<< std::endl;
			    ((TGraph*)_grRatioRad30.At(itrk))->SetPoint(icut,(bkgEff30*nb[icent])/(sigEff30*ns[icent]),0.1);
            }
            else {
                std::cout <<"WARNING: no signal events in this bin. B/S set to -1.0. " << std::endl;
                ((TGraph*)_grRatioRad30.At(itrk))->SetPoint(icut,-1.0,0.1);
            }


			std::cout << "dR 0.3 bkg fraction at point " << icut << " = " << bkgEff30 << std::endl;
			std::cout << "dR 0.3 sig fraction at point " << icut << " = " << sigEff30 << std::endl;

			treeMcJ1Set->Draw("pt>>hmcJ1Iso",selCuts40,"hf");
			treeMcJ2Set->Draw("pt>>hmcJ2Iso",selCuts40,"hf");
			treeMcJ3Set->Draw("pt>>hmcJ3Iso",selCuts40,"hf");
			//weighted muon yield per Jx mc event
			hmcQCDSetIso40->Add(hmcJ1Iso,hmcJ2Iso,(wtJ1)/(nMcJ1)*scaleFactor, (wtJ2)/(nMcJ2)*scaleFactor); 
			hmcQCDSetIso40->Add(hmcQCDSetIso40,hmcJ3Iso,1.0,(wtJ3)/(nMcJ3)*scaleFactor);

			treeMcWSet->Draw("pt>>hmcWSetIso",selCuts40+"&&prompt==24","hf");		

			double bkgEff40 = hmcQCDSetIso40->Integral(ptSigBinLo,ptSigBinUpp)/hmcQCDSet->Integral(ptSigBinLo,ptSigBinUpp);
			double sigEff40 = hmcWSetIso->Integral(ptSigBinLo,ptSigBinUpp)/hmcWSet->Integral(ptSigBinLo,ptSigBinUpp);

			((TGraph*)_grEffRad40.At(itrk))->SetPoint(icut,sigEff40,bkgEff40);
			std::cout << "sig fraction at point " << icut << " = " << sigEff40 << "*" << ns[icent] << " = " <<
                sigEff40*ns[icent] << std::endl;
			std::cout << "bkg fraction at point " << icut << " = " << bkgEff40 << "*" << nb[icent] << " = " <<
                bkgEff40*nb[icent] << std::endl;
            if(sigEff40*ns[icent]>0.0){
                std::cout << "That gives a background:signal ratio of " << (bkgEff40*nb[icent])/(sigEff40*ns[icent])<< std::endl;
			    ((TGraph*)_grRatioRad40.At(itrk))->SetPoint(icut,(bkgEff40*nb[icent])/(sigEff40*ns[icent]),0.1);
            }
            else {
                std::cout <<"WARNING: no signal events in this bin. B/S set to -1.0. " << std::endl;
                ((TGraph*)_grRatioRad40.At(itrk))->SetPoint(icut,-1.0,0.1);
            }



			std::cout << "dR 0.4 bkg fraction at point " << icut << " = " << bkgEff40 << std::endl;
			std::cout << "dR 0.4 sig fraction at point " << icut << " = " << sigEff40 << std::endl;

			treeMcJ1Set->Draw("pt>>hmcJ1Iso",selCuts50,"hf");
			treeMcJ2Set->Draw("pt>>hmcJ2Iso",selCuts50,"hf");
			treeMcJ3Set->Draw("pt>>hmcJ3Iso",selCuts50,"hf");
			//weighted muon yield per Jx mc event
			hmcQCDSetIso50->Add(hmcJ1Iso,hmcJ2Iso,(wtJ1)/(nMcJ1)*scaleFactor, (wtJ2)/(nMcJ2)*scaleFactor); 
			hmcQCDSetIso50->Add(hmcQCDSetIso50,hmcJ3Iso,1.0,(wtJ3)/(nMcJ3)*scaleFactor);

			treeMcWSet->Draw("pt>>hmcWSetIso",selCuts50+"&&prompt==24","hf");		

			double bkgEff50 = hmcQCDSetIso50->Integral(ptSigBinLo,ptSigBinUpp)/hmcQCDSet->Integral(ptSigBinLo,ptSigBinUpp);
			double sigEff50 = hmcWSetIso->Integral(ptSigBinLo,ptSigBinUpp)/hmcWSet->Integral(ptSigBinLo,ptSigBinUpp);

			((TGraph*)_grEffRad50.At(itrk))->SetPoint(icut,sigEff50,bkgEff50);
			std::cout << "sig fraction at point " << icut << " = " << sigEff50 << "*" << ns[icent] << " = " <<
                sigEff50*ns[icent] << std::endl;
			std::cout << "bkg fraction at point " << icut << " = " << bkgEff50 << "*" << nb[icent] << " = " <<
                bkgEff50*nb[icent] << std::endl;
            if(sigEff50*ns[icent]>0.0){
                std::cout << "That gives a background:signal ratio of " << (bkgEff50*nb[icent])/(sigEff50*ns[icent])<< std::endl;
			    ((TGraph*)_grRatioRad50.At(itrk))->SetPoint(icut,(bkgEff50*nb[icent])/(sigEff50*ns[icent]),0.1);
            }
            else {
                std::cout <<"WARNING: no signal events in this bin. B/S set to -1.0. " << std::endl;
                ((TGraph*)_grRatioRad50.At(itrk))->SetPoint(icut,-1.0,0.1);
            }



			std::cout << "dR 0.5 bkg fraction at point " << icut << " = " << bkgEff50 << std::endl;
			std::cout << "dR 0.5 sig fraction at point " << icut << " = " << sigEff50 << std::endl;

		}

	TCanvas* cEff = new TCanvas("cEff","cEff",600,600);

/*	TGraph* grEff10c = ((TGraph*)_grEffRad10.At(itrk))->Clone("grEff10c"); 
	TGraph* grEff15c = ((TGraph*)_grEffRad15.At(itrk))->Clone("grEff15c"); 
	TGraph* grEff20c = ((TGraph*)_grEffRad20.At(itrk))->Clone("grEff20c"); 
	TGraph* grEff30c = ((TGraph*)_grEffRad30.At(itrk))->Clone("grEff30c"); 
	TGraph* grEff40c = ((TGraph*)_grEffRad40.At(itrk))->Clone("grEff40c"); 

	grEff10c->GetYaxis()->SetTitle("background efficiency");
	grEff10c->GetXaxis()->SetTitle("signal efficiency");
	grEff10c->SetMarkerColor(kMagenta); grEff10c->SetMarkerStyle(kFullCircle);
	grEff15c->SetMarkerColor(kOrange+5); grEff15c->SetMarkerStyle(kOpenSquare);
	grEff20c->SetMarkerColor(kRed); grEff20c->SetMarkerStyle(kOpenCircle);
	grEff30c->SetMarkerColor(kGreen); grEff30c->SetMarkerStyle(kOpenTriangleDown);
	grEff40c->SetMarkerColor(kBlue); grEff40c->SetMarkerStyle(kOpenCross);
	
	grEff10c->Draw("ap");
	grEff15c->Draw("psame");
	grEff20c->Draw("psame");
	grEff30c->Draw("psame");
	grEff40c->Draw("psame");

	gROOT->LoadMacro("AtlasUtils.C");
	myText(0.33,0.89, (Color_t)kBlack, sTitle);

	TArrow ar4(1.0,0.062,0.975,0.187,0.05,"|>");
	ar4.SetFillColor(kRed);
	TLegend* leg = new TLegend(0.648, 0.622, 0.918, 0.815);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->AddEntry(grEff10c,"#Delta R 0.1","p");
	leg->AddEntry(grEff15c,"#Delta R 0.15","p");
	leg->AddEntry(grEff20c,"#Delta R 0.2","p");
	leg->AddEntry(grEff30c,"#Delta R 0.3","p");
	leg->AddEntry(grEff40c,"#Delta R 0.4","p");

	leg->Draw(); 
*/	
	Write(outFile,_grEffRad10.At(itrk), sGr10);
	Write(outFile,_grEffRad15.At(itrk), sGr15);
	Write(outFile,_grEffRad20.At(itrk), sGr20);
	Write(outFile,_grEffRad30.At(itrk), sGr30);
	Write(outFile,_grEffRad40.At(itrk), sGr40);
	Write(outFile,_grEffRad50.At(itrk), sGr50);

	Write(outFile,_grRatioRad10.At(itrk), sGrRatio10);
	Write(outFile,_grRatioRad15.At(itrk), sGrRatio15);
	Write(outFile,_grRatioRad20.At(itrk), sGrRatio20);
	Write(outFile,_grRatioRad30.At(itrk), sGrRatio30);
	Write(outFile,_grRatioRad40.At(itrk), sGrRatio40);
	Write(outFile,_grRatioRad50.At(itrk), sGrRatio50);

       }//lower track pt	
    }//icent

 }
