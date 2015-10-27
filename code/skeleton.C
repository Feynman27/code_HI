void plotMptMtCorrelation(){
      
      TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";
	  TString fileNameIn = baseString+ ".root";

      ///reco level variables
      float eLossNt[50];
      float scatNt[50];
      float compNt[50];
      float ptNt[50];
      float mtNt[50];
      float etaNt[50];
      float phiNt[50];
      float mptPhi;
      float chargeNt[50];
      int promptNt[50];
      float centralityNt;
      float mpt;
      float ptconeNt[50];
      int valNt[50], truthMatchedNt[50],ZDYNt[50], matched1[50], matched2[50], matched3[50];
      int nmu,trig1,trig2,trig3,trig4,trig5;

      TChain* tree = new TChain("tree","tree");
      tree->Add(fileNameIn);
      std::cout << "Number of entries in tree: " << tree->GetEntries() << std::endl;
      tree->SetBranchAddress("pt",&ptNt);

      tree->SetBranchAddress("prompt",&promptNt);
      tree->SetBranchAddress("eLoss", &eLossNt);
      tree->SetBranchAddress("ptcone20ID3", &ptconeNt);
      tree->SetBranchAddress("scat", &scatNt);
      tree->SetBranchAddress("comp", &compNt);
      tree->SetBranchAddress("pt", &ptNt);
      tree->SetBranchAddress("mt", &mtNt);
      tree->SetBranchAddress("eta", &etaNt);
      tree->SetBranchAddress("phi", &phiNt);
      tree->SetBranchAddress("nu_phi", &mptPhi);
      tree->SetBranchAddress("charge", &chargeNt);
      tree->SetBranchAddress("val", &valNt); 
      tree->SetBranchAddress("ZDY", &ZDYNt); 
      tree->SetBranchAddress("centrality", &centralityNt);
      tree->SetBranchAddress("nu_pt", &mpt);
      tree->SetBranchAddress("mu_muid_n", &nmu);


      // --- Set branch status ---
      tree->SetBranchStatus("*",0) ;
      tree->SetBranchStatus("mu_muid_n", 1);
      tree->SetBranchStatus("eLoss", 1);
      tree->SetBranchStatus("ptcone20ID3", 1);
      tree->SetBranchStatus("scat", 1);
      tree->SetBranchStatus("comp", 1);
      tree->SetBranchStatus("prompt",1);
      tree->SetBranchStatus("pt", 1);
      tree->SetBranchStatus("mt", 1);
      tree->SetBranchStatus("eta", 1);
      tree->SetBranchStatus("phi", 1);
      tree->SetBranchStatus("nu_phi", 1);
      tree->SetBranchStatus("charge", 1);
      tree->SetBranchStatus("val", 1); 
      tree->SetBranchStatus("ZDY", 1); 
      tree->SetBranchStatus("centrality", 1);
      tree->SetBranchStatus("nu_pt", 1);

      for(int iev = 0; iev < tree->GetEntries(); ++iev){
 
        tree->LoadTree(iev);
        tree->GetEntry(iev);

        for(int imu=0; imu<nmu; ++imu){

        }//imu
       }//iev
}
