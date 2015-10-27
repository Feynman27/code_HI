TString format(double value) {
  std::stringstream svalue;
  svalue  << std::setprecision(2) << value;
  return svalue.str();
}


void Write(TFile* outFile, TObject* grIso, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Saving TGraph named " << sFile << std::endl;
	  grIso->Write(sFile);
	  std::cout << "Done." << std::endl;
    	  gDirectory = dir;
}

void isoEfficiencyData(){

	bool doWMc = false;
	bool doJx  = true;

	if(doWMc) std::cout << "Running on W mc. " << std::endl;
	if(doJx) std::cout << "Running on Jx+1mu mc. " << std::endl;
	else std::cout << "Running on 2011 PbPb data. " << std::endl;

	SetAtlasStyle();
	TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
	
	/// --- Open output file ---
	TDirectory *dir = gDirectory;
	TFile *outFile = new TFile("isoEfficiencyDeltaRGraphs.root","RECREATE");
	gDirectory = dir;

	//TString fileNameDataIn = baseString+"HardProbesFiles/HISingleMuonHP.03.11.2013";
	TString fileNameDataIn; 
	TString fileNameMCJ1In = baseString+"MonteCarloFiles/QCD/HISingleMuonJ1_1muPYTHIADataOverlay.03.13.2013";
	TString fileNameMCJ2In = baseString+"MonteCarloFiles/QCD/HISingleMuonJ2_1muPYTHIADataOverlay.03.13.2013";
	TString fileNameMCJ3In = baseString+"MonteCarloFiles/QCD/HISingleMuonJ3_1muPYTHIADataOverlay.03.13.2013";
	TFile* fMcJ1Set=0;
	TFile* fMcJ2Set=0;
	TFile* fMcJ3Set=0;
	TFile* fDataSet=0;
	TTree *treeMcJ1Set=0;
	TTree *treeMcJ2Set=0;
	TTree *treeMcJ3Set=0;
	TTree *treeDataSet=0;

	if(doWMc) fileNameDataIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.03.13.2013";
	if(!doJx) {
		fDataSet = new TFile(fileNameDataIn+".root", "READ");
		treeDataSet  = (TTree*)fDataSet->Get("tree");
	}
	else{
		fMcJ1Set = new TFile(fileNameMCJ1In+".root", "READ");
		fMcJ2Set = new TFile(fileNameMCJ2In+".root", "READ");
		fMcJ3Set = new TFile(fileNameMCJ3In+".root", "READ");

		treeMcJ1Set = (TTree*)fMcJ1Set->Get("tree");
		treeMcJ2Set = (TTree*)fMcJ2Set->Get("tree");
		treeMcJ3Set = (TTree*)fMcJ3Set->Get("tree");
	}


	TString cuts ="abs(charge)==1";
	cuts+="&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&ZDY==0";
	if(doWMc) cuts+="&&prompt==24";

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


	//centrality
	std::vector <float> centBins;
	centBins.push_back(0.0);
/*	centBins.push_back(0.05);
	centBins.push_back(0.1);
	centBins.push_back(0.15);
	centBins.push_back(0.2);
	centBins.push_back(0.4);*/

	centBins.push_back(0.8);

	const int centralityBins = centBins.size()-1;

	//cone radii around muon
	int cone10 = 10;
	int cone15 = 15;
	int cone20 = 20;
	int cone30 = 30;
	int cone40 = 40;
	//lower track pt used in iso calculation
	const int trkPtCut = 5;
	int trkPtLo[trkPtCut] = {2,3,4,5,6};
	//int trkPtLo[2] = {2,3};
	//const int trkPtCut = 2;

	//const int isolationMuBins = 2;
	const int isolationMuBins = 4;
	double isoCut[isolationMuBins] = {0.1, 0.2, 0.3, 0.4};
	//double isoCut[isolationMuBins] = {0.2,0.3};
	double scaleFactor = 1.0e9;

	TH1F* hmcDataSet = new TH1F("hmcDataSet","hmcDataSet",nbins,0.0,ptUpp);
	TH1F* hmcDataSetIso = new TH1F("hmcDataSetIso","hmcDataSetIso",nbins,0.0,ptUpp);

	TH1F* hmcJ1 = new TH1F("hmcJ1","hmcJ1",nbins,0.0,ptUpp);
	TH1F* hmcJ2 = new TH1F("hmcJ2","hmcJ2",nbins,0.0,ptUpp);
	TH1F* hmcJ3 = new TH1F("hmcJ3","hmcJ3",nbins,0.0,ptUpp);

	TH1F* hmcJ1Iso = new TH1F("hmcJ1Iso","hmcJ1Iso",nbins,0.0,ptUpp);
	TH1F* hmcJ2Iso = new TH1F("hmcJ2Iso","hmcJ2Iso",nbins,0.0,ptUpp);
	TH1F* hmcJ3Iso = new TH1F("hmcJ3Iso","hmcJ3Iso",nbins,0.0,ptUpp);

	TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",nbins,0.0,ptUpp);
	TH1F* hmcQCDSetIso = new TH1F("hmcQCDSetIso","hmcQCDSetIso",nbins,0.0,ptUpp);

	TH1F* hmcWSet = new TH1F("hmcWSet","hmcWSet",nbins,0.0,ptUpp);
	TH1F* hmcWSetIso = new TH1F("hmcWSetIso","hmcWSetIso",nbins,0.0,ptUpp);

	TH1F* hmcJ1Cent = new TH1F("hmcJ1Cent","hmcJ1Cent",nbins,0.0,ptUpp);
	TH1F* hmcJ2Cent = new TH1F("hmcJ2Cent","hmcJ2Cent",nbins,0.0,ptUpp);
	TH1F* hmcJ3Cent = new TH1F("hmcJ3Cent","hmcJ3Cent",nbins,0.0,ptUpp);

	TList _grEffRad; 
	
	//graphs of the efficiency as fcn of cone radius round muon
	const int nCones = 5;
	int arrConeRad[nCones] = {10,15,20,30,40};
	TString arrConeCut[nCones]; 
	TGraph* _grSig = new TGraph(nCones);

	const unsigned int nGraphs = trkPtCut*isolationMuBins;
	TString sGrName[nGraphs];

	double wtJ1 = 6.6523e-03*1.8770e-4/64.0e-3; double wtJ2 = 1.4941e-02*8.2788e-6/64.0e-3; double wtJ3 = 2.4284e-02*2.9426e-7/64.0e-3;
	for(int icent=0; icent<centralityBins; icent++){
	
		TString sAddCent = "&&centrality>"; sAddCent += centBins.at(icent); sAddCent += "&&centrality<"; sAddCent += centBins.at(icent+1);
		TString selCuts = cuts; selCuts += sMissPtUp; selCuts += sMissPtLo;selCuts +=sMtCutLo; selCuts +=sMtCutUp; selCuts += sAddCent; 
		std::cout << "Selection cuts w/o isolation: " << selCuts << std::endl;
		
		TString sCentrality = "centrality>"; sCentrality+= centBins.at(icent); sCentrality+="&&centrality<"; sCentrality+=centBins.at(icent+1);

	        double nMcJ1; double nMcJ2; double nMcJ3;
		if(doJx) {
		  //event counting histos
		  treeMcJ1Set->Draw("centrality>>hmcJ1Cent",sCentrality,"hf");
		  treeMcJ2Set->Draw("centrality>>hmcJ2Cent",sCentrality,"hf");
		  treeMcJ3Set->Draw("centrality>>hmcJ3Cent",sCentrality,"hf");
	          nMcJ1 = hmcJ1Cent->Integral();  nMcJ2 = hmcJ2Cent->Integral();  nMcJ3 = hmcJ3Cent->Integral();
		}
	   //lower pT threshold on the tracks
	   for(int itrk=0; itrk<trkPtCut; itrk++){

			TString sTitle = "p_{T}^{trk}>"; sTitle+=trkPtLo[itrk]; sTitle+="GeV";
	
		//isolation_mu upper limits
		for(int icut=0; icut<isolationMuBins; icut++){
			
			//unique index for given lower track pT
			//threshold and isolation_mu cut
			int index = itrk*isolationMuBins+icut;

			TString sGr = "grIsoEffTrkPt"; sGr+=trkPtLo[itrk]; 
			sGr += "IsolationMu"; 
			sGr+=icut; 
			sGrName[index] = sGr;

			//1 TGraph per trk pT threshold and
			//isolation_mu cut
			_grEffRad.Add( new TGraph(nCones));

			std::cout << "cuts w/ no isolation requirement: " << selCuts << std::endl;
			if(doJx){

				treeMcJ1Set->Draw("pt>>hmcJ1",selCuts,"hf");
				treeMcJ2Set->Draw("pt>>hmcJ2",selCuts,"hf");
				treeMcJ3Set->Draw("pt>>hmcJ3",selCuts,"hf");
				//weighted muon yield per Jx mc event
				hmcDataSet->Add(hmcJ1,hmcJ2,(wtJ1)/(nMcJ1)*scaleFactor, (wtJ2)/(nMcJ2)*scaleFactor); 
				hmcDataSet->Add(hmcDataSet,hmcJ3,1.0,(wtJ3)/(nMcJ3)*scaleFactor);

			}
			else{ 
			   treeDataSet->Draw("pt>>hmcDataSet",selCuts,"hf");		
			}

			//different cone radii around the muon
			for (int iCone = 0; iCone < nCones; iCone++){
			  
			  //construct isolation cut for given cone size 
			  TString sCone = "&&ptcone"; 
			  sCone+=arrConeRad[iCone]; sCone+="ID"; sCone+=trkPtLo[itrk]; sCone+="/pt<"; sCone+=isoCut[icut]; 
			  std::cout << "Isolation cut: " << sCone << std::endl;

			  TString selCutsCone = selCuts + sCone;
			  std::cout << "cuts after isolation requirement: " << selCutsCone << std::endl;
		          if(doJx){	
			        treeMcJ1Set->Draw("pt>>hmcJ1Iso",selCutsCone,"hf");
		         	treeMcJ2Set->Draw("pt>>hmcJ2Iso",selCutsCone,"hf");
		        	treeMcJ3Set->Draw("pt>>hmcJ3Iso",selCutsCone,"hf");
				//weighted muon yield per Jx mc event
				hmcDataSetIso->Add(hmcJ1Iso,hmcJ2Iso,(wtJ1)/(nMcJ1)*scaleFactor, (wtJ2)/(nMcJ2)*scaleFactor); 
				hmcDataSetIso->Add(hmcDataSetIso,hmcJ3Iso,1.0,(wtJ3)/(nMcJ3)*scaleFactor);
			  }
			  else treeDataSet->Draw("pt>>hmcDataSetIso",selCutsCone,"hf");		

			  //calculate the efficiency 
			  double sigEff = hmcDataSetIso->Integral()/hmcDataSet->Integral();

			  //plot the efficiency as a function of cone size
			  ((TGraph*)_grEffRad.At(index))->SetPoint(iCone,arrConeRad[iCone],sigEff);
			  std::cout << "dR sig fraction at point for cone radius " << arrConeRad[iCone] << " = " << sigEff << std::endl;
			
			}//cone size 

			//save the TGraph
			
			Write(outFile,_grEffRad.At(index), sGrName[index]);

		}//isolation_mu upper limit

       }//lower track pt	
    }//icent

 }
