

void plotZdEta(){
        TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";

	TFile *fIn = new TFile(baseString+"MonteCarloFiles/Zmumu/HISingleMuonMCZmumu.12.30.2012.root","READ");
	TFile *fOut = new TFile("dEtaZmumu.root","recreate");
	TTree* tree = (TTree*)fIn->Get("tree");

	int nmu;
	int mother[50],daughter[50];
	float eta[50];
	TH1* heta = new TH1F("heta","heta",150,0.0,3.5);

	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("mc_mu_gen_eta",1);
	tree->SetBranchStatus("mc_mu_gen_mothertype",1);
	tree->SetBranchStatus("mc_mu_gen_type",1);
	tree->SetBranchStatus("mc_mu_n",1);


	tree->SetBranchAddress("mc_mu_gen_eta",&eta);
	tree->SetBranchAddress("mc_mu_gen_mothertype",&mother);
	tree->SetBranchAddress("mc_mu_gen_type",&daughter);
	tree->SetBranchAddress("mc_mu_n",&nmu);
	
	int nEvents = tree->GetEntries();

	for(int iev=0;iev<nEvents;iev++){
	  tree->LoadTree(iev);
	  tree->GetEntry(iev);
	  if(iev%10000==0) std::cout << "Event: " << iev << " Muons : " << nmu << std::endl;
	  for(int imu=0; imu<nmu; imu++){
		double eta1 = eta[imu];
		for(int jmu=imu+1; jmu<nmu; jmu++){
			double eta2 = eta[jmu];
			if(abs(mother[imu])==23&&abs(daughter[imu])==13&&abs(mother[jmu])==23&&abs(daughter[jmu])==13) {
				double deta = abs(eta1-eta2); heta->Fill(deta);
				//if(deta==0.0)std::cout << deta << std::endl;
			}
		}
	  }
	}

	fOut->cd();
	fOut->Write();
	fOut->Close();
}
