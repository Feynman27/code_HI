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

void isoEfficiencyDatav02(){

	SetAtlasStyle();
	TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
	
	/// --- Open output file ---
	TDirectory *dir = gDirectory;
	TFile *outFile = new TFile("isoEfficiencyDataGraphsv02.root","RECREATE");
	gDirectory = dir;

	TString fileNameDataIn = baseString+"HardProbesFiles/HISingleMuonHP.03.11.2013";
	TFile* fDataSet = new TFile(fileNameDataIn+".root", "READ");

	if ( !fDataSet->IsOpen() ) {
	    std::cout << fDataSet << " not found!" << std::endl;
	    exit(0);
	}

  	// --- W set ---


	double ptLow=4.0 ;double ptUpp=100.0 ;double etaLow=0.0; double etaUpp=2.5; double centralityLow=0.0;double centralityUpp=0.8;
	int nbins = 100;
	double ptCtrlLo = 10.0; double ptCtrlUpp = 20.0; double ptSigLo = 25.0; 
	//define bins for control region
	int ptCtrlBinLo = (nbins/ptUpp)*ptCtrlLo+1; int ptCtrlBinUpp = (nbins/ptUpp)*ptCtrlUpp+1; 
	//define bins for signal region
	int ptSigBinLo = (nbins/ptUpp)*ptSigLo+1; int ptSigBinUpp = (nbins/ptUpp)*ptUpp+1;  

	double mtcutLow = 40.0; float mtmax = 400.0;

	double missPtCut = 25.0; 
	double missPtMax = 9000.0; 

	float eLoss[50];
	float scat[50];
	float pt[50];
	float mt[50];
	float eta[50];
	float phi[50];
	float charge[50];
	float centrality;
	float nu_pt;
	float ptcone10ID2[50];
	float ptcone10ID3[50];
	float ptcone10ID4[50];
	float ptcone10ID5[50];
	float ptcone10ID6[50];

	float ptcone15ID2[50];
	float ptcone15ID3[50];
	float ptcone15ID4[50];
	float ptcone15ID5[50];
	float ptcone15ID6[50];

	float ptcone20ID2[50];
	float ptcone20ID3[50];
	float ptcone20ID4[50];
	float ptcone20ID5[50];
	float ptcone20ID6[50];

	float ptcone30ID2[50];
	float ptcone30ID3[50];
	float ptcone30ID4[50];
	float ptcone30ID5[50];
	float ptcone30ID6[50];

	float ptcone40ID2[50];
	float ptcone40ID3[50];
	float ptcone40ID4[50];
	float ptcone40ID5[50];
	float ptcone40ID6[50];

	int val[50], ZDY[50];
	int nmu;

	TTree *tree  = (TTree*)fDataSet->Get("tree");
	tree->SetBranchAddress("charge",&charge);
	tree->SetBranchAddress("val",&val);
	tree->SetBranchAddress("eLoss",&eLoss);
	tree->SetBranchAddress("scat",&scat);
	tree->SetBranchAddress("ZDY",&ZDY);
	tree->SetBranchAddress("pt",&pt);
	tree->SetBranchAddress("eta",&eta);
	tree->SetBranchAddress("nu_pt",&nu_pt);
	tree->SetBranchAddress("mt",&mt);
	tree->SetBranchAddress("centrality",&centrality);
	tree->SetBranchAddress("ptcone10ID2",&ptcone10ID2);
	tree->SetBranchAddress("ptcone15ID2",&ptcone15ID2);
	tree->SetBranchAddress("ptcone20ID2",&ptcone10ID2);
	tree->SetBranchAddress("ptcone30ID2",&ptcone10ID2);
	tree->SetBranchAddress("ptcone40ID2",&ptcone10ID2);
	tree->SetBranchAddress("ptcone10ID3",&ptcone10ID3);
	tree->SetBranchAddress("ptcone15ID3",&ptcone15ID3);
	tree->SetBranchAddress("ptcone20ID3",&ptcone10ID3);
	tree->SetBranchAddress("ptcone30ID3",&ptcone10ID3);
	tree->SetBranchAddress("ptcone40ID3",&ptcone10ID3);
	tree->SetBranchAddress("ptcone10ID4",&ptcone10ID4);
	tree->SetBranchAddress("ptcone15ID4",&ptcone15ID4);
	tree->SetBranchAddress("ptcone20ID4",&ptcone10ID4);
	tree->SetBranchAddress("ptcone30ID4",&ptcone10ID4);
	tree->SetBranchAddress("ptcone40ID4",&ptcone10ID4);
	tree->SetBranchAddress("ptcone10ID5",&ptcone10ID5);
	tree->SetBranchAddress("ptcone15ID5",&ptcone15ID5);
	tree->SetBranchAddress("ptcone20ID5",&ptcone10ID5);
	tree->SetBranchAddress("ptcone30ID5",&ptcone10ID5);
	tree->SetBranchAddress("ptcone40ID5",&ptcone10ID5);
	tree->SetBranchAddress("ptcone10ID6",&ptcone10ID6);
	tree->SetBranchAddress("ptcone15ID6",&ptcone15ID6);
	tree->SetBranchAddress("ptcone20ID6",&ptcone10ID6);
	tree->SetBranchAddress("ptcone30ID6",&ptcone10ID6);
	tree->SetBranchAddress("ptcone40ID6",&ptcone10ID6);
        tree->SetBranchAddress("mu_muid_n", &nmu);

        tree->SetBranchStatus("*",0);
        tree->SetBranchStatus("mu_muid_n", 1);
	tree->SetBranchStatus("charge",1);
	tree->SetBranchStatus("val",1);
	tree->SetBranchStatus("eLoss",1);
	tree->SetBranchStatus("scat",1);
	tree->SetBranchStatus("ZDY",1);
	tree->SetBranchStatus("pt",1);
	tree->SetBranchStatus("eta",1);
	tree->SetBranchStatus("nu_pt",1);
	tree->SetBranchStatus("mt",1);
	tree->SetBranchStatus("centrality",1);
	tree->SetBranchStatus("ptcone10ID2",1);
	tree->SetBranchStatus("ptcone15ID2",1);
	tree->SetBranchStatus("ptcone20ID2",1);
	tree->SetBranchStatus("ptcone30ID2",1);
	tree->SetBranchStatus("ptcone40ID2",1);
	tree->SetBranchStatus("ptcone10ID3",1);
	tree->SetBranchStatus("ptcone15ID3",1);
	tree->SetBranchStatus("ptcone20ID3",1);
	tree->SetBranchStatus("ptcone30ID3",1);
	tree->SetBranchStatus("ptcone40ID3",1);
	tree->SetBranchStatus("ptcone10ID4",1);
	tree->SetBranchStatus("ptcone15ID4",1);
	tree->SetBranchStatus("ptcone20ID4",1);
	tree->SetBranchStatus("ptcone30ID4",1);
	tree->SetBranchStatus("ptcone40ID4",1);
	tree->SetBranchStatus("ptcone10ID5",1);
	tree->SetBranchStatus("ptcone15ID5",1);
	tree->SetBranchStatus("ptcone20ID5",1);
	tree->SetBranchStatus("ptcone30ID5",1);
	tree->SetBranchStatus("ptcone40ID5",1);
	tree->SetBranchStatus("ptcone10ID6",1);
	tree->SetBranchStatus("ptcone15ID6",1);
	tree->SetBranchStatus("ptcone20ID6",1);
	tree->SetBranchStatus("ptcone30ID6",1);
	tree->SetBranchStatus("ptcone40ID6",1);


        
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

	//lower track pt used in iso calculation
	const int trkPtCut = 5;
	int trkPtLo[trkPtCut] = {2,3,4,5,6};

	const int isolationMuBins = 4;
	double isoCut[isolationMuBins] = {0.1, 0.2, 0.3, 0.4};

	const int nCones = 5;
	const int nConeVars = trkPtCut*nCones;
	TH1F* hmcDataSet = new TH1F("hmcDataSet","hmcDataSet",nbins,0.0,ptUpp);
	TH1F* hptcone10ID2[isolationMuBins] ; 
	TH1F* hptcone15ID2[isolationMuBins] ; 	
        TH1F* hptcone20ID2[isolationMuBins] ; 
	TH1F* hptcone30ID2[isolationMuBins] ; 
	TH1F* hptcone40ID2[isolationMuBins] ; 
	TH1F* hptcone10ID3[isolationMuBins] ; 
	TH1F* hptcone15ID3[isolationMuBins] ; 
	TH1F* hptcone20ID3[isolationMuBins] ; 
	TH1F* hptcone30ID3[isolationMuBins] ; 
	TH1F* hptcone40ID3[isolationMuBins] ; 
	TH1F* hptcone10ID4[isolationMuBins] ; 
	TH1F* hptcone15ID4[isolationMuBins] ; 
	TH1F* hptcone20ID4[isolationMuBins] ; 
	TH1F* hptcone30ID4[isolationMuBins] ; 
	TH1F* hptcone40ID4[isolationMuBins] ; 
	TH1F* hptcone10ID5[isolationMuBins] ; 
	TH1F* hptcone15ID5[isolationMuBins] ;
	TH1F* hptcone20ID5[isolationMuBins] ; 
	TH1F* hptcone30ID5[isolationMuBins] ; 
	TH1F* hptcone40ID5[isolationMuBins] ; 
	TH1F* hptcone10ID6[isolationMuBins] ; 
	TH1F* hptcone15ID6[isolationMuBins] ; 
	TH1F* hptcone20ID6[isolationMuBins] ; 
	TH1F* hptcone30ID6[isolationMuBins] ; 
	TH1F* hptcone40ID6[isolationMuBins] ; 

	TList _grEffID2; 
	TList _grEffID3; 
	TList _grEffID4; 
	TList _grEffID5; 
	TList _grEffID6; 
	
	//graphs of the efficiency as fcn of cone radius round muon
	const int nCones = 5;
	int arrConeRad[nCones] = {10,15,20,30,40};
	int arrEff[nCones] ;
	TString arrConeCut[nCones]; 
	TGraph* _grSig = new TGraph(nCones);

	const unsigned int nGraphs = trkPtCut*isolationMuBins;
	TString sGrName[nGraphs];
			     for(int icut=0; icut<isolationMuBins; icut++){

				hptcone10ID2[icut] = new TH1F("hptcone10ID2","hptcone10ID2",nbins,0.0,ptUpp);
				hptcone10ID3[icut] = new TH1F("hptcone10ID3","hptcone10ID3",nbins,0.0,ptUpp);
				hptcone10ID4[icut] = new TH1F("hptcone10ID4","hptcone10ID4",nbins,0.0,ptUpp);
				hptcone10ID5[icut] = new TH1F("hptcone10ID5","hptcone10ID5",nbins,0.0,ptUpp);
				hptcone10ID6[icut] = new TH1F("hptcone10ID6","hptcone10ID6",nbins,0.0,ptUpp);

				hptcone15ID2[icut] = new TH1F("hptcone15ID2","hptcone15ID2",nbins,0.0,ptUpp);
				hptcone15ID3[icut] = new TH1F("hptcone15ID3","hptcone15ID3",nbins,0.0,ptUpp);
				hptcone15ID4[icut] = new TH1F("hptcone15ID4","hptcone15ID4",nbins,0.0,ptUpp);
				hptcone15ID5[icut] = new TH1F("hptcone15ID5","hptcone15ID5",nbins,0.0,ptUpp);
				hptcone15ID6[icut] = new TH1F("hptcone15ID6","hptcone15ID6",nbins,0.0,ptUpp);

				hptcone20ID2[icut] = new TH1F("hptcone20ID2","hptcone20ID2",nbins,0.0,ptUpp);
				hptcone20ID3[icut] = new TH1F("hptcone20ID3","hptcone20ID3",nbins,0.0,ptUpp);
				hptcone20ID4[icut] = new TH1F("hptcone20ID4","hptcone20ID4",nbins,0.0,ptUpp);
				hptcone20ID5[icut] = new TH1F("hptcone20ID5","hptcone20ID5",nbins,0.0,ptUpp);
				hptcone20ID6[icut] = new TH1F("hptcone20ID6","hptcone20ID6",nbins,0.0,ptUpp);

				hptcone30ID2[icut] = new TH1F("hptcone30ID2","hptcone30ID2",nbins,0.0,ptUpp);
				hptcone30ID3[icut] = new TH1F("hptcone30ID3","hptcone30ID3",nbins,0.0,ptUpp);
				hptcone30ID4[icut] = new TH1F("hptcone30ID4","hptcone30ID4",nbins,0.0,ptUpp);
				hptcone30ID5[icut] = new TH1F("hptcone30ID5","hptcone30ID5",nbins,0.0,ptUpp);
				hptcone30ID6[icut] = new TH1F("hptcone30ID6","hptcone30ID6",nbins,0.0,ptUpp);

				hptcone40ID2[icut] = new TH1F("hptcone40ID2","hptcone40ID2",nbins,0.0,ptUpp);
				hptcone40ID3[icut] = new TH1F("hptcone40ID3","hptcone40ID3",nbins,0.0,ptUpp);
				hptcone40ID4[icut] = new TH1F("hptcone40ID4","hptcone40ID4",nbins,0.0,ptUpp);
				hptcone40ID5[icut] = new TH1F("hptcone40ID5","hptcone40ID5",nbins,0.0,ptUpp);
				hptcone40ID6[icut] = new TH1F("hptcone40ID6","hptcone40ID6",nbins,0.0,ptUpp);

}

	std::cout << "Entries : " << tree->GetEntries() << std::endl;
	for(int iev=0; iev<tree->GetEntries(); iev++){
		
		tree->LoadTree(iev);
		tree->GetEntry(iev);

		if(iev%10000==0) std::cout << "Event: " << iev << std::endl; 

		for(int imu=0;imu<nmu;imu++){
			
			for(int icent=0; icent<centralityBins; icent++){

			   if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400) 
					hmcDataSet->Fill(pt[imu]);	

			     for(int icut=0; icut<isolationMuBins; icut++){

			        //ptcone10
				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone10ID2[imu]/pt[imu]<isoCut[icut]) {hptcone10ID2[icut]->Fill(pt[imu]);	}

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone10ID3[imu]/pt[imu]<isoCut[icut]) hptcone10ID3[icut]->Fill(pt[imu]);

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone10ID4[imu]/pt[imu]<isoCut[icut]) hptcone10ID4[icut]->Fill(pt[imu]);

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone10ID5[imu]/pt[imu]<isoCut[icut]) hptcone10ID5[icut]->Fill(pt[imu]);

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone10ID6[imu]/pt[imu]<isoCut[icut]) hptcone10ID6[icut]->Fill(pt[imu]);

				//ptcone15
				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone15ID2[imu]/pt[imu]<isoCut[icut]) hptcone15ID2[icut]->Fill(pt[imu]);	

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone15ID3[imu]/pt[imu]<isoCut[icut]) hptcone15ID3[icut]->Fill(pt[imu]);

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone15ID4[imu]/pt[imu]<isoCut[icut]) hptcone15ID4[icut]->Fill(pt[imu]);

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone15ID5[imu]/pt[imu]<isoCut[icut]) hptcone15ID5[icut]->Fill(pt[imu]);

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone15ID6[imu]/pt[imu]<isoCut[icut]) hptcone15ID6[icut]->Fill(pt[imu]);

				//ptcone20
				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone20ID2[imu]/pt[imu]<isoCut[icut]) hptcone20ID2[icut]->Fill(pt[imu]);	

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone20ID3[imu]/pt[imu]<isoCut[icut]) hptcone20ID3[icut]->Fill(pt[imu]);

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone20ID4[imu]/pt[imu]<isoCut[icut]) hptcone20ID4[icut]->Fill(pt[imu]);

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone20ID5[imu]/pt[imu]<isoCut[icut]) hptcone20ID5[icut]->Fill(pt[imu]);

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone20ID6[imu]/pt[imu]<isoCut[icut]) hptcone20ID6[icut]->Fill(pt[imu]);

				//ptcone30
				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone30ID2[imu]/pt[imu]<isoCut[icut]) hptcone30ID2[icut]->Fill(pt[imu]);	

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone30ID3[imu]/pt[imu]<isoCut[icut]) hptcone30ID3[icut]->Fill(pt[imu]);

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone30ID4[imu]/pt[imu]<isoCut[icut]) hptcone30ID4[icut]->Fill(pt[imu]);

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone30ID5[imu]/pt[imu]<isoCut[icut]) hptcone30ID5[icut]->Fill(pt[imu]);

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone30ID6[imu]/pt[imu]<isoCut[icut]) hptcone30ID6[icut]->Fill(pt[imu]);
				//ptcone40	
				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone40ID2[imu]/pt[imu]<isoCut[icut]) hptcone40ID2[icut]->Fill(pt[imu]);	

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone40ID3[imu]/pt[imu]<isoCut[icut]) hptcone40ID3[icut]->Fill(pt[imu]);

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone40ID4[imu]/pt[imu]<isoCut[icut]) hptcone40ID4[icut]->Fill(pt[imu]);

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone40ID5[imu]/pt[imu]<isoCut[icut]) hptcone40ID5[icut]->Fill(pt[imu]);

				if(abs(charge[imu])==1&&val[imu]>11&&abs(eLoss[imu])<0.5&&abs(scat[imu])<4.0&&ZDY[imu]==0&&pt[imu]>ptLow
				   &&pt[imu]<ptUpp&&abs(eta[imu])>etaLow&&abs(eta[imu])<etaUpp&&centrality>centBins.at(icent)&&centrality<centBins.at(icent+1)
				   &&nu_pt>25&&nu_pt<9000&&mt[imu]>40&&mt[imu]<400
				   &&ptcone40ID6[imu]/pt[imu]<isoCut[icut]) hptcone40ID6[icut]->Fill(pt[imu]);

			}//icut
	   } //icent
	} //imu
      } //iev

	std::cout << "Done fill histograms. Now taking ratios..." << std::endl;

	for(int icent=0; icent<centralityBins; icent++){

		//isolation_mu upper limits
		for(int icut=0; icut<isolationMuBins; icut++){
			
			//1 TGraph per trk pT threshold and
			//isolation_mu cut
			_grEffID2.Add( new TGraph(nCones));

			//calculate the efficiency  
			double sigEffCone10ID2 = hptcone10ID2[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone15ID2 = hptcone15ID2[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone20ID2 = hptcone20ID2[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone30ID2 = hptcone30ID2[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone40ID2 = hptcone40ID2[icut]->Integral()/hmcDataSet->Integral();

			arrEff[] = {sigEffCone10ID2,sigEffCone15ID2,sigEffCone20ID2,sigEffCone30ID2,sigEffCone40ID2};

			//plot the efficiency as a function of cone size
			_grEffID2.Add(new TGraph(nCones,arrConeRad,arrEff));

			//save the TGraph
			
			TString sGr = "grIsoEffTrkPt"; 
			TString sGrID2=sGr+trkPtLo[0]; 
			sGrID2 += "IsolationMu"; 
			sGrID2+=icut; 

			Write(outFile,_grEffID2.At(icut), sGrID2);

		      //3GeV trk thr
			_grEffID3.Add( new TGraph(nCones));

			//calculate the efficiency  
			double sigEffCone10ID3 = hptcone10ID3[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone15ID3 = hptcone15ID3[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone20ID3 = hptcone20ID3[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone30ID3 = hptcone30ID3[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone40ID3 = hptcone40ID3[icut]->Integral()/hmcDataSet->Integral();

			arrEff[] = {sigEffCone10ID3,sigEffCone15ID3,sigEffCone20ID3,sigEffCone30ID3,sigEffCone40ID3};

			//plot the efficiency as a function of cone size
			_grEffID3.Add(new TGraph(nCones,arrConeRad,arrEff));

			//save the TGraph
			
			TString sGr = "grIsoEffTrkPt"; 
			TString sGrID3=sGr+trkPtLo[1]; 
			sGrID3 += "IsolationMu"; 
			sGrID3+=icut; 

			Write(outFile,_grEffID3.At(icut), sGrID3);

		        //4GeV trk thr
			_grEffID4.Add( new TGraph(nCones));

			//calculate the efficiency  
			double sigEffCone10ID4 = hptcone10ID4[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone15ID4 = hptcone15ID4[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone20ID4 = hptcone20ID4[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone30ID4 = hptcone30ID4[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone40ID4 = hptcone40ID4[icut]->Integral()/hmcDataSet->Integral();

			arrEff[] = {sigEffCone10ID4,sigEffCone15ID4,sigEffCone20ID4,sigEffCone30ID4,sigEffCone40ID4};

			//plot the efficiency as a function of cone size
			_grEffID4.Add(new TGraph(nCones,arrConeRad,arrEff));

			//save the TGraph
			
			TString sGr = "grIsoEffTrkPt"; 
			TString sGrID4=sGr+trkPtLo[2]; 
			sGrID4 += "IsolationMu"; 
			sGrID4+=icut; 

			Write(outFile,_grEffID4.At(icut), sGrID4);

		        //5GeV trk thr
			_grEffID5.Add( new TGraph(nCones));

			//calculate the efficiency  
			double sigEffCone10ID5 = hptcone10ID5[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone15ID5 = hptcone15ID5[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone20ID5 = hptcone20ID5[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone30ID5 = hptcone30ID5[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone40ID5 = hptcone40ID5[icut]->Integral()/hmcDataSet->Integral();

			arrEff[] = {sigEffCone10ID5,sigEffCone15ID5,sigEffCone20ID5,sigEffCone30ID5,sigEffCone40ID5};

			//plot the efficiency as a function of cone size
			_grEffID5.Add(new TGraph(nCones,arrConeRad,arrEff));

			//save the TGraph
			
			TString sGr = "grIsoEffTrkPt"; 
			TString sGrID5=sGr+trkPtLo[3]; 
			sGrID5 += "IsolationMu"; 
			sGrID5+=icut; 

			Write(outFile,_grEffID5.At(icut), sGrID5);

		        //6GeV trk thr
			_grEffID6.Add( new TGraph(nCones));

			//calculate the efficiency  
			double sigEffCone10ID6 = hptcone10ID6[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone15ID6 = hptcone15ID6[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone20ID6 = hptcone20ID6[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone30ID6 = hptcone30ID6[icut]->Integral()/hmcDataSet->Integral();
			double sigEffCone40ID6 = hptcone40ID6[icut]->Integral()/hmcDataSet->Integral();

			arrEff[] = {sigEffCone10ID6,sigEffCone15ID6,sigEffCone20ID6,sigEffCone30ID6,sigEffCone40ID6};

			//plot the efficiency as a function of cone size
			_grEffID6.Add(new TGraph(nCones,arrConeRad,arrEff));

			//save the TGraph
			
			TString sGrID6=sGr+trkPtLo[4]; 
			sGrID6 += "IsolationMu"; 
			sGrID6+=icut; 

			Write(outFile,_grEffID6.At(icut), sGrID6);

		}//isolation_mu upper limit

    }//icent

 }
