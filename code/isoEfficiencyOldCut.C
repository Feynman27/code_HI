void Write(TFile* outFile, TObject* grIso, TString sFile){
	  TDirectory *dir = gDirectory;
	  outFile->cd();
	  std::cout << "Writing TGraph to root file..." << std::endl;
	  grIso->Write(sFile);
	  std::cout << "Done." << std::endl;
    	  gDirectory = dir;
}



void isoEfficiencyOldCut(){

	SetAtlasStyle();
	TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";


	TDirectory *dir = gDirectory;
	TFile *outFile = new TFile("isoEfficiencyBuiltInVarGraphs.root","RECREATE");
	gDirectory = dir;
	//data overlay
	TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuMCDataOverlay.2012.11.21";
	//TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuPYTHIADataOverlay.03.17.2013";
	TFile* fMcWSet = new TFile(fileNameMCWIn+".root", "READ");

	if ( !fMcWSet->IsOpen() ) {
	    std::cout << fMcWSet << " not found!" << std::endl;
	    exit(0);
	}

  	// --- W set ---

	TString fileNameMCJ1In = baseString+"MonteCarloFiles/QCD/HISingleMuonMcJ1.2012.11.29";
	///new file with custom isolation cuts
	//TString fileNameMCJ1In = baseString+"MonteCarloFiles/QCD/HISingleMuonJ1_1muPYTHIADataOverlay.03.17.2013";
	TFile* fMcJ1Set = new TFile(fileNameMCJ1In+".root", "READ");
	//J2 1 muon-filter 
	TString fileNameMCJ2In = baseString+"MonteCarloFiles/QCD/HISingleMuonMcJ2.2012.11.29";
	//TString fileNameMCJ2In = baseString+"MonteCarloFiles/QCD/HISingleMuonJ2_1muPYTHIADataOverlay.03.17.2013";
	TFile* fMcJ2Set = new TFile(fileNameMCJ2In+".root", "READ");
	TString fileNameMCJ3In = baseString+"MonteCarloFiles/QCD/HISingleMuonMcJ3.2012.11.29";
	//TString fileNameMCJ3In = baseString+"MonteCarloFiles/QCD/HISingleMuonJ3_1muPYTHIADataOverlay.03.17.2013";
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
	std::cout << "Signal region from bin " << ptSigBinLo << " - " << ptSigBinUpp << std::endl;

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
//	centBins.push_back(0.1);
//	centBins.push_back(0.15);
//	centBins.push_back(0.2);
//	centBins.push_back(0.4);

	centBins.push_back(0.8);

	//ncoll
	std::vector <float> ncoll;

//	ncoll.push_back(1683.3); //0-5
//	ncoll.push_back(1318.0); //5-10
//	ncoll.push_back(1500.6); //0-10
//	ncoll.push_back(1035.4); //10-15
//	ncoll.push_back(811.2); //15-20
//	ncoll.push_back(923.3); //10-20

//	ncoll.push_back(1212.0);//0-20
//	ncoll.push_back(440.6); //20-40
//	ncoll.push_back(77.8); //40-80
        ncoll.push_back(361.6);//0-80

	const int centralityBins = centBins.size()-1;

	int cone20 = 20;
	int cone30 = 30;
	int cone40 = 40;
	const int cutBins = 11;
	double isoCut[cutBins] = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 5.0, 10.0, 100.0};
	//double isoCut[cutBins] = {0.2,0.3};
	double evData = 1.0e9;

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

	TGraph* grEff20 = new TGraph(cutBins);
	TGraph* grEff30 = new TGraph(cutBins);
	TGraph* grEff40 = new TGraph(cutBins);


	double wtJ1 = 6.6523e-03*1.8770e-4/64.0e-3; double wtJ2 = 1.4941e-02*8.2788e-6/64.0e-3; double wtJ3 = 2.4284e-02*2.9426e-7/64.0e-3;

	for(int icent=0; icent<centralityBins; icent++){

	        double arrCentWidth = centBins.at(icent+1)-centBins.at(icent);
	        double scaleFactor = arrCentWidth*ncoll[icent]*evData;
		TString sGr = "grIsoEffBuiltInVar"; sGr+="cent"; sGr+=icent; 
	 	TString sGr20 = sGr + "ConeRadius20";
		TString sGr30 = sGr + "ConeRadius30";  
		TString sGr40 = sGr + "ConeRadius40";  

		for(int icut=0; icut<cutBins; icut++){
		
			TString sAddCent = "&&centrality>"; sAddCent += centBins.at(icent); sAddCent += "&&centrality<"; sAddCent += centBins.at(icent+1);
			TString selCuts = cuts; selCuts += sMissPtUp; selCuts += sMissPtLo;selCuts +=sMtCutLo; selCuts +=sMtCutUp; selCuts += sAddCent; 
			std::cout << "Selection cuts w/o isolation: " << selCuts << std::endl;

			TString sCone = "&&ptcone"; 
			TString sCone20 = sCone; sCone20+=cone20; sCone20+="/pt<"; 
			TString sCone30 = sCone; sCone30+=cone30; sCone30+="/pt<"; 
			TString sCone40 = sCone; sCone40+=cone40; sCone40+="/pt<"; 
			sCone20+=isoCut[icut];
			sCone30+=isoCut[icut];
			sCone40+=isoCut[icut];

			TString selCuts20 = selCuts + sCone20;
			TString selCuts30 = selCuts + sCone30;
			TString selCuts40 = selCuts + sCone40;
			
			std::cout << "isolation cut for dR 0.2: " << selCuts20 << std::endl;
			std::cout << "isolation cut for dR 0.3: " << selCuts30 << std::endl;
			std::cout << "isolation cut for dR 0.4: " << selCuts40 << std::endl;

			TString sCentrality = "centrality>"; sCentrality+= centBins.at(icent); sCentrality+="&&centrality<"; sCentrality+=centBins.at(icent+1);
			//event counting histos
			treeMcJ1Set->Draw("centrality>>hmcJ1Cent",sCentrality,"hf");
			treeMcJ2Set->Draw("centrality>>hmcJ2Cent",sCentrality,"hf");
			treeMcJ3Set->Draw("centrality>>hmcJ3Cent",sCentrality,"hf");

			
			double nMcJ1 = hmcJ1Cent->Integral(); double nMcJ2 = hmcJ2Cent->Integral(); double nMcJ3 = hmcJ3Cent->Integral();

			treeMcJ1Set->Draw("pt>>hmcJ1",selCuts,"hf");
			treeMcJ2Set->Draw("pt>>hmcJ2",selCuts,"hf");
			treeMcJ3Set->Draw("pt>>hmcJ3",selCuts,"hf");
			//weighted muon yield per Jx mc event
			hmcQCDSet->Add(hmcJ1,hmcJ2,(wtJ1)/(nMcJ1)*scaleFactor, (wtJ2)/(nMcJ2)*scaleFactor); 
			hmcQCDSet->Add(hmcQCDSet,hmcJ3,1.0,(wtJ3)/(nMcJ3)*scaleFactor);

			std::cout << "W cuts no iso: " << selCuts << std::endl;
			treeMcWSet->Draw("pt>>hmcWSet",selCuts+"&&prompt==24","hf");		

			treeMcJ1Set->Draw("pt>>hmcJ1Iso",selCuts20,"hf");
			treeMcJ2Set->Draw("pt>>hmcJ2Iso",selCuts20,"hf");
			treeMcJ3Set->Draw("pt>>hmcJ3Iso",selCuts20,"hf");
			//weighted muon yield per Jx mc event
			hmcQCDSetIso->Add(hmcJ1Iso,hmcJ2Iso,(wtJ1)/(nMcJ1)*scaleFactor, (wtJ2)/(nMcJ2)*scaleFactor); 
			hmcQCDSetIso->Add(hmcQCDSetIso,hmcJ3Iso,1.0,(wtJ3)/(nMcJ3)*scaleFactor);

			std::cout << "W cuts after iso: " << selCuts20 << std::endl;
			treeMcWSet->Draw("pt>>hmcWSetIso",selCuts20+"&&prompt==24","hf");		

			double bkgEff20 = hmcQCDSetIso->Integral(ptSigBinLo,ptSigBinUpp)/hmcQCDSet->Integral(ptSigBinLo,ptSigBinUpp);
			double sigEff20 = hmcWSetIso->Integral(ptSigBinLo,ptSigBinUpp)/hmcWSet->Integral(ptSigBinLo,ptSigBinUpp);

			std::cout <<  hmcQCDSetIso->Integral(ptSigBinLo,ptSigBinUpp) << "/" << hmcQCDSet->Integral(ptSigBinLo,ptSigBinUpp) 
						<< " = " << bkgEff20 << std::endl;

			std::cout <<  hmcWSetIso->Integral(ptSigBinLo,ptSigBinUpp) << "/" << hmcWSet->Integral(ptSigBinLo,ptSigBinUpp) 
						<< " = " << sigEff20 << std::endl;



			grEff20->SetPoint(icut,sigEff20,bkgEff20);
			std::cout << "dR 0.2 bkg fraction at point " << icut << " = " << bkgEff20 << std::endl;
			std::cout << "dR 0.2 sig fraction at point " << icut << " = " << sigEff20 << std::endl;

			treeMcJ1Set->Draw("pt>>hmcJ1Iso",selCuts30,"hf");
			treeMcJ2Set->Draw("pt>>hmcJ2Iso",selCuts30,"hf");
			treeMcJ3Set->Draw("pt>>hmcJ3Iso",selCuts30,"hf");
			//weighted muon yield per Jx mc event
			hmcQCDSetIso->Add(hmcJ1Iso,hmcJ2Iso,(wtJ1)/(nMcJ1)*scaleFactor, (wtJ2)/(nMcJ2)*scaleFactor); 
			hmcQCDSetIso->Add(hmcQCDSetIso,hmcJ3Iso,1.0,(wtJ3)/(nMcJ3)*scaleFactor);

			treeMcWSet->Draw("pt>>hmcWSetIso",selCuts30+"&&prompt==24","hf");		

			double bkgEff30 = hmcQCDSetIso->Integral(ptSigBinLo,ptSigBinUpp)/hmcQCDSet->Integral(ptSigBinLo,ptSigBinUpp);
			double sigEff30 = hmcWSetIso->Integral(ptSigBinLo,ptSigBinUpp)/hmcWSet->Integral(ptSigBinLo,ptSigBinUpp);

			grEff30->SetPoint(icut,sigEff30,bkgEff30);
			std::cout << "dR 0.3 bkg fraction at point " << icut << " = " << bkgEff30 << std::endl;
			std::cout << "dR 0.3 sig fraction at point " << icut << " = " << sigEff30 << std::endl;

			treeMcJ1Set->Draw("pt>>hmcJ1Iso",selCuts40,"hf");
			treeMcJ2Set->Draw("pt>>hmcJ2Iso",selCuts40,"hf");
			treeMcJ3Set->Draw("pt>>hmcJ3Iso",selCuts40,"hf");
			//weighted muon yield per Jx mc event
			hmcQCDSetIso->Add(hmcJ1Iso,hmcJ2Iso,(wtJ1)/(nMcJ1)*scaleFactor, (wtJ2)/(nMcJ2)*scaleFactor); 
			hmcQCDSetIso->Add(hmcQCDSetIso,hmcJ3Iso,1.0,(wtJ3)/(nMcJ3)*scaleFactor);

			treeMcWSet->Draw("pt>>hmcWSetIso",selCuts40+"&&prompt==24","hf");		

			double bkgEff40 = hmcQCDSetIso->Integral(ptSigBinLo,ptSigBinUpp)/hmcQCDSet->Integral(ptSigBinLo,ptSigBinUpp);
			double sigEff40 = hmcWSetIso->Integral(ptSigBinLo,ptSigBinUpp)/hmcWSet->Integral(ptSigBinLo,ptSigBinUpp);

			grEff40->SetPoint(icut,sigEff40,bkgEff40);
			std::cout << "dR 0.4 bkg fraction at point " << icut << " = " << bkgEff40 << std::endl;
			std::cout << "dR 0.4 sig fraction at point " << icut << " = " << sigEff40 << std::endl;

		}
	

	TCanvas* cEff = new TCanvas("cEff","cEff",600,600);

	TGraph* grEff20c = grEff20->Clone("grEff20c"); 
	TGraph* grEff30c = grEff30->Clone("grEff30c"); 
	TGraph* grEff40c = grEff40->Clone("grEff40c"); 

	grEff20c->GetYaxis()->SetTitle("background efficiency");
	grEff20c->GetXaxis()->SetTitle("signal efficiency");
	grEff20c->SetMarkerColor(kRed); grEff20c->SetMarkerStyle(kOpenCircle);
	grEff30c->SetMarkerColor(kGreen); grEff30c->SetMarkerStyle(kOpenTriangleDown);
	grEff40c->SetMarkerColor(kBlue); grEff40c->SetMarkerStyle(kOpenCross);
	
	grEff20c->Draw("ap");
	grEff30c->Draw("psame");
	grEff40c->Draw("psame");

	TArrow ar4(1.0,0.062,0.975,0.187,0.05,"|>");
	ar4.SetFillColor(kRed);
	TLegend* leg = new TLegend(0.648, 0.622, 0.918, 0.815);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->AddEntry(grEff20c,"#Delta R 0.2","p");
	leg->AddEntry(grEff30c,"#Delta R 0.3","p");
	leg->AddEntry(grEff40c,"#Delta R 0.4","p");

	leg->Draw(); 

	Write(outFile,grEff20c, sGr20);
	Write(outFile,grEff30c, sGr30);
	Write(outFile,grEff40c, sGr40);
	
	}//icent

 }
