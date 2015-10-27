void plotWtdJxSples(){
	
	bool doEta = false;
	bool doCentrality = false;
	SetAtlasStyle();
	TString baseString = "/usatlas/u/tbales/scratch/MonteCarloFiles/QCD/";

    ///data
	TString fileNameDataIn = "/usatlas/u/tbales/scratch/HISingleMuonHardProbesData.04.17.2013";
	TFile* fDataSet = new TFile(fileNameDataIn+".root", "READ");

	//J1 1 muon-filter 
	TString fileNameMCJ1In = baseString+"HISingleMuonJ1_1muPYTHIADataOverlay.04.13.2013";
	//TString fileNameMCJ1In = baseString+"HISingleMuonMCQCDJ1_1mu.03.04.2013";
	TFile* fMcJ1Set = new TFile(fileNameMCJ1In+".root", "READ");

	if ( !fMcJ1Set->IsOpen() ) {
	    std::cout << fMcJ1Set << " not found!" << std::endl;
	    exit(0);
	}
	//J2 1 muon-filter 
	TString fileNameMCJ2In = baseString+"HISingleMuonJ2_1muPYTHIADataOverlay.04.13.2013";
	//TString fileNameMCJ2In = baseString+"HISingleMuonMCQCDJ2_1mu.03.04.2013";
	TFile* fMcJ2Set = new TFile(fileNameMCJ2In+".root", "READ");

	if ( !fMcJ2Set->IsOpen() ) {
	    std::cout << fMcJ2Set << " not found!" << std::endl;
	    exit(0);
	}
	//J3 1 muon-filter 
	TString fileNameMCJ3In = baseString+"HISingleMuonJ3_1muPYTHIADataOverlay.04.13.2013";
	//TString fileNameMCJ3In = baseString+"HISingleMuonMCQCDJ3_1mu.03.04.2013";
	TFile* fMcJ3Set = new TFile(fileNameMCJ3In+".root", "READ");

	if ( !fMcJ3Set->IsOpen() ) {
	    std::cout << fMcJ3Set << " not found!" << std::endl;
	    exit(0);
	}

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


	std::vector<double> centralityBins;

	centralityBins.push_back(0.00);
	if (doCentrality) {

		centralityBins.push_back(0.05);
		centralityBins.push_back(0.10);
		centralityBins.push_back(0.15);
		centralityBins.push_back(0.20);
		centralityBins.push_back(0.40);
		//centralityBins.push_back(0.60);
	}
	centralityBins.push_back(0.80);

	const int nCentralityBins = centralityBins.size()-1;

	TString cuts ="abs(charge)==1";
        cuts+="&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0";

	double mtmax = 300.0; double ncoll = 361.6; double ptLow = 0.0; double ptUpp=300.0; double etaLow = 0.0; double etaUpp = 2.5;
	double centralityLow=0.0; double centralityUpp = 0.8;

	int nBins = 300;
  	// --- QCD set ---
	TH1F* hmcJ1Set = new TH1F("hmcJ1Set","hmcJ1Set",nBins,0.0,mtmax);
	TH1F* hmcJ2Set = new TH1F("hmcJ2Set","hmcJ2Set",nBins,0.0,mtmax);
	TH1F* hmcJ3Set = new TH1F("hmcJ3Set","hmcJ3Set",nBins,0.0,mtmax);
	TH1F* hmcJ1Cent = new TH1F("hmcJ1Cent","hmcJ1Cent",nBins,0.0,mtmax);
	TH1F* hmcJ2Cent = new TH1F("hmcJ2Cent","hmcJ2Cent",nBins,0.0,mtmax);
	TH1F* hmcJ3Cent = new TH1F("hmcJ3Cent","hmcJ3Cent",nBins,0.0,mtmax);
	TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",nBins,0.0,mtmax);
	TH1F* hDataSet = new TH1F("hDataSet","hDataSet",nBins,0.0,mtmax);

	double mtcutLow = 40.0;
	double missPtCut = 25.0; 
	double isoCut = 0.3;

	TString sMissPtLo = "&&nu_pt<"; sMissPtLo+=missPtCut;  TString sMissPtUp = "&&nu_pt>"; sMissPtUp+=missPtCut; 
	TString sIsoCut = "&&ptcone20/pt<"; sIsoCut+=isoCut; TString sNonIsoCut = "&&ptcone20/pt>"; sNonIsoCut+=isoCut;
	TString sIsoCut40 = "&&ptcone40/pt<"; sIsoCut40+=isoCut; TString sNonIsoCut40 = "&&ptcone40/pt>"; sNonIsoCut40+=isoCut;
	TString sMtCutLo = "&&mt>"; sMtCutLo+=mtcutLow; TString sMtCutUp = "&&mt<";sMtCutUp+=mtmax;
	TString sPtCut = "&&pt>"; sPtCut += ptLow; sPtCut +="&&pt<"; sPtCut += ptUpp;

	cuts+="&&centrality>"; cuts += centralityLow; cuts+="&&centrality<"; cuts += centralityUpp;
	//assume QCD bkg fraction is constant over eta to increase statistics in MC
	TString cutsQCD = cuts+"&&abs(eta)>0.0&&abs(eta)<2.5";

    TString sHiPtTrigger =
    "&&((EF_mu10_MSonly_EFFS_L1ZDC&&EF_mu10_MSonly_EFFS_L1ZDC_Matched20)||(EF_mu10_MSonly_EFFS_L1TE10&&EF_mu10_MSonly_EFFS_L1TE10_Matched20)||(EF_mu10_MSonly_EFFS_L1TE20&&EF_mu10_MSonly_EFFS_L1TE20_Matched20))"; 
 
    TString cutsData = cutsQCD+sHiPtTrigger;
	//cutsQCD += sMissPtUp; scutsQCD += sIsoCut; scutsQCD +=sMtCutLo; scutsQCD +=sMtCutUp; scutsQCD +="&&ZDY==0"; scutsQCD += sPtCut;
	//std::cout << "QCD cuts " << scutsQCD << std::endl;

	TTree *treeMcJ1Set = (TTree*)fMcJ1Set->Get("tree");
	TTree *treeMcJ2Set = (TTree*)fMcJ2Set->Get("tree");
	TTree *treeMcJ3Set = (TTree*)fMcJ3Set->Get("tree");
	TTree *treeDataSet = (TTree*)fDataSet->Get("tree");

	treeMcJ1Set->Draw("pt>>hmcJ1Set",cutsQCD,"hf");		
	treeMcJ2Set->Draw("pt>>hmcJ2Set",cutsQCD,"hf");		
	treeMcJ3Set->Draw("pt>>hmcJ3Set",cutsQCD,"hf");		
	treeDataSet->Draw("pt>>hDataSet",cutsData,"pe");		

	//event counting histos
	TString sCentrality = "centrality>"; sCentrality+= centralityLow; sCentrality+="&&centrality<"; sCentrality+=centralityUpp;
	treeMcJ1Set->Draw("centrality>>hmcJ1Cent",sCentrality,"hf");
	treeMcJ2Set->Draw("centrality>>hmcJ2Cent",sCentrality,"hf");
	treeMcJ3Set->Draw("centrality>>hmcJ3Cent",sCentrality,"hf");
	double nMCJ1 = hmcJ1Cent->Integral(); double nMCJ2 = hmcJ2Cent->Integral(); double nMCJ3 = hmcJ3Cent->Integral();

	//weight the Jx samples according to cross-sections; takes into account prob of find muon+jet (AMI)
	double wtJ1 = 6.6523e-03*1.8770e-4/64.0e-3; double wtJ2 = 1.4941e-02*8.2788e-6/64.0e-3; double wtJ3 = 2.4284e-02*2.9426e-7/64.0e-3;
	double arrCentWidth = centralityUpp-centralityLow;
	double evData = 1.03e9; //number of hp events
	double scaleFactor = arrCentWidth*ncoll*evData;

    TCanvas* c0 = new TCanvas("c0","c0",600,600);
    TH1F* hmcJ1Setc = (TH1F*)hmcJ1Set->Clone("hmcJ1Setc");
    hmcJ1Setc->Scale((wtJ1)/(nMCJ1)*scaleFactor);
    hmcJ1Setc->SetFillColor(kYellow);
    TH1F* hmcJ2Setc = (TH1F*)hmcJ2Set->Clone("hmcJ2Setc");
    hmcJ2Setc->Scale((wtJ2)/(nMCJ2)*scaleFactor);
    hmcJ2Setc->SetFillColor(kRed);
    TH1F* hmcJ3Setc = (TH1F*)hmcJ3Set->Clone("hmcJ3Setc");
    hmcJ3Setc->Scale((wtJ3)/(nMCJ3)*scaleFactor);
    hmcJ3Setc->SetFillColor(kGreen);

    hmcJ1Setc->Draw("hist");
    hmcJ2Setc->Draw("histsame");
    hmcJ3Setc->Draw("histsame");

    TCanvas* c1 = new TCanvas("c1","c1",600,600);
	//hmcQCDSet->Add(hmcJ1Set,hmcJ2Set,(wtJ1)/(nMCJ1)*scaleFactor, (wtJ2)/(nMCJ2)*scaleFactor); 
	//hmcQCDSet->Add(hmcJ3Set,(wtJ3)/(nMCJ3)*scaleFactor);
    hmcQCDSet->Add(hmcJ1Setc,hmcJ2Setc); hmcQCDSet->Add(hmcJ3Setc);
	TH1F* hmcQCDSetc = (TH1F*)hmcQCDSet->Clone("hmcQCDSetc");

	TH1F* hDataSetc = (TH1F*)hDataSet->Clone("hDataSetc");
	
	Int_t ci;   // for color index setting
        hmcQCDSet->GetXaxis()->SetRangeUser(10.0,120.0);
	hmcQCDSet->SetFillColor(15);
	ci = TColor::GetColor("#99ccff");
	hmcQCDSetc->SetFillColor(ci);
	hmcQCDSet->Draw("hist");
        int binCtrlDataLo = hDataSetc->FindBin(10.0); int binCtrlDataUp = hDataSetc->FindBin(20.0);
        int binCtrlMCLo = hmcQCDSet->FindBin(10.0); int binCtrlMCUp = hmcQCDSetc->FindBin(20.0);
	double rescaleFactor = hDataSetc->Integral(binCtrlDataLo,binCtrlDataUp)/hmcQCDSet->Integral(binCtrlMCLo,binCtrlMCUp);
        hmcQCDSetc->Scale(rescaleFactor);
	hmcQCDSetc->Draw("histsame");
	hDataSetc->Draw("pesame");
	hmcQCDSet->Draw("sameaxis");
   TLegend* leg = new TLegend(0.192953,0.3216783,0.5620805,0.5576923,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04545455);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->AddEntry(hDataSetc,"Data 2011","pe");
   leg->AddEntry(hmcQCDSet,"QCD MC (Scaled to #LT N_{coll} #GT","f");
   leg->AddEntry(hmcQCDSetc,"QCD MC (Re-scaled to ctrl region","f");
   leg->Draw();

   TLatex *   tex = new TLatex(0.1862416,0.2342657,"#int Ldt #approx 0.140 nb^{-1}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03846154);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.4848993,0.236014,"#sqrt{s_{NN}} = 2.76 TeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04195804);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.6073826,0.8566434,"ATLAS");
tex->SetNDC();
   tex->SetTextFont(72);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.7701342,0.8601399,"Internal");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();



   c1->SetLogy();
}
