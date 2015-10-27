

void pythia8(Int_t nev = 100, Int_t ndeb = 1) {
	
   	char *p8dataenv = gSystem->Getenv("PYTHIA8DATA");
   	if (!p8dataenv) {
      	char *p8env = gSystem->Getenv("PYTHIA8");
      	if (!p8env) {
         	Error("pythia8.C",
               		"Environment variable PYTHIA8 must contain path to pythia directory!");
         	return;
      	}
      		TString p8d = p8env;
      		p8d += "/xmldoc";
      		gSystem->Setenv("PYTHIA8DATA", p8d);
   	}

	char* path = gSystem->ExpandPathName("$PYTHIA8DATA") ;
	if (gSystem->AccessPathName(path)) {
		Warning("pythia8.C","Environment variable PYTHIA8DATA must contain path to pythia8100/xmldoc directory !");
		return;
	}

	//load libraries
	gSystem->Load("$PYTHIA8/lib/libpythia8") ;
	gSystem->Load("libEG") ;
	gSystem->Load("libEGPythia8");

	//histos
	TH1F* ptH = new TH1F("ptH","pt",100,0.0,100.0);

	//array of particles
	TCloneArrary* particles = new TClonesArray("TParticle",1000);

	//create pythia8 object
	TPythia8* = new TPythia8();

	//configure
	pythia8->ReadString("HardQCD:qqbar2bbbar = on");	//turn on hf production from quark-anti-quark 
	pythia8->ReadString("HardQCD:qqbar2ccbar = on");
	pythia8->ReadString("5:onMode = off"); //turn off bottom decay channels
	pythia8->ReadString("5:onIfMatch = 13") ; //turn on muon decay channel
	pythia8->ReadString("4:onMode = off"); //turn off charm decay channels
	pythia8->ReadString("4:onIfMatch = -13") ; //turn on muon decay channel

	//initialize

	pythia8->Initialize(2212/*p*/,2212/*p*/,7000./*TeV*/);

	//event loop
	for (Int_t iev = 0; iev < nev; iev++){
		pythia8->GenerateEvent();
		if (iev < ndeb) pythia8->EventListing();
		pythia8->ImportParticles(particles,"All") ;
		Int_t np = particles->GetEntriesFast();
		//particle loop
		for (Int_t ip = 0; ip < np; ip++) {
			TParticle* part = (TParticle*) particles->At(ip);
			Int_t ist = part->GetStatusCode();
			Int_t pdg = part->GetPdgCode();
			if (ist != 1) continue;
			Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
			if(charge==0) continue;
			Float_t pt = part->Pt();

			if (pt>0.) ptH->Fill(pt,1./(2.*pt));
		}
	}

	pythia8->PrintStatistics();

	TCanvas *c1 = new TCanvas("c1","c1",600,600);
	gPad->SetLogy();
	ptH->Scale(5./Float_t(nev));
	ptH->Draw();
	ptH->SetXTitle("p_{T} [GeV/c]");
	ptH->SetYTitle("dN/dp_{t}^{2} [GeV/c]^{-2}");


}
