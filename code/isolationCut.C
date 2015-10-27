{
	SetAtlasStyle();
	TString date = "Dec3.2012";
	TString baseString = "/mnt/Lustre/cgrp/atlas_hi/tbalestri/";

	//data overlay
//	TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonWmunuMCDataOverlay.2012.11.21";
	TString fileNameMCWIn = baseString+"MonteCarloFiles/Wmunu/HISingleMuonMCWmunu.12.25.2012";
	TFile* fMcWSet = new TFile(fileNameMCWIn+".root", "READ");

	if ( !fMcWSet->IsOpen() ) {
	    std::cout << fMcWSet << " not found!" << std::endl;
	    exit(0);
	}

  	// --- W set ---

	TString fileNameMCJ1In = baseString+"MonteCarloFiles/QCD/HISingleMuonMcJ1.2012.11.29";
	TFile* fMcJ1Set = new TFile(fileNameMCJ1In+".root", "READ");
	//J2 1 muon-filter 
	TString fileNameMCJ2In = baseString+"MonteCarloFiles/QCD/HISingleMuonMcJ2.2012.11.29";
	TFile* fMcJ2Set = new TFile(fileNameMCJ2In+".root", "READ");
	TString fileNameMCJ3In = baseString+"MonteCarloFiles/QCD/HISingleMuonMcJ3.2012.11.29";
	TFile* fMcJ3Set = new TFile(fileNameMCJ3In+".root", "READ");
	TString cuts ="abs(charge)==1";
	cuts+="&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&ZDY==0";
	double ptLow=4.0 ;double ptUpp=100.0 ;double etaLow=0.0; double etaUpp=2.5; double centralityLow=0.0;double centralityUpp=0.8;
	TH1F* hmcJ1Set = new TH1F("hmcJ1Set","hmcJ1Set",100,0.0,100.0);
	TH1F* hmcJ2Set = new TH1F("hmcJ2Set","hmcJ2Set",100,0.0,100.0);
	TH1F* hmcJ3Set = new TH1F("hmcJ3Set","hmcJ3Set",100,0.0,100.0);
	TH1F* hmcJ1Cent = new TH1F("hmcJ1Cent","hmcJ1Cent",100,0.0,100.0);
	TH1F* hmcJ2Cent = new TH1F("hmcJ2Cent","hmcJ2Cent",100,0.0,100.0);
	TH1F* hmcJ3Cent = new TH1F("hmcJ3Cent","hmcJ3Cent",100,0.0,100.0);
	TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",100,0.0,100.0);

	double mtcutLow = 40.0; float mtmax = 400.0;

	cuts+="&&pt>"; cuts += ptLow; cuts+="&&pt<"; cuts += ptUpp;
	cuts+="&&abs(eta)>"; cuts += etaLow; cuts+="&&abs(eta)<"; cuts += etaUpp;
	 //cuts+="&&centrality>"; cuts += centralityLow; cuts+="&&centrality<"; cuts += centralityUpp;

	double missPtCut = 25.0; 
	double isoCut10 = 0.1;
	double isoCut20 = 0.2;
	double isoCut30 = 0.3;
	int cone = 20;

	TString sCone = "&&ptcone"; sCone+=cone; sCone+="/pt<";
	TString sIsoCut10 = sCone; sIsoCut10+=isoCut10;
	TString sIsoCut20 = sCone; sIsoCut20+=isoCut20;
	TString sIsoCut30 = sCone; sIsoCut30+=isoCut30;

	TString sMtCutLo = "&&mt>"; sMtCutLo+=mtcutLow; TString sMtCutUp = "&&mt<";sMtCutUp+=mtmax;
	TString sMissPtLo = "&&nu_pt<"; sMissPtLo+=missPtCut; TString sMissPtUp = "&&nu_pt>"; sMissPtUp+=missPtCut;

	TString scutsMc = cuts; scutsMc += sMissPtUp;scutsMc +=sMtCutLo; /*scutsMc +=sMtCutUp;*/
	TString scutsMc10 = cuts; scutsMc10 += sMissPtUp; scutsMc10 +=sMtCutLo; scutsMc10 +=sMtCutUp; scutsMc10 += sIsoCut10;
	TString scutsMc20 = cuts; scutsMc20 += sMissPtUp; scutsMc20 +=sMtCutLo; scutsMc20 +=sMtCutUp; scutsMc20 += sIsoCut20;
	TString scutsMc30 = cuts; scutsMc30 += sMissPtUp; scutsMc30 +=sMtCutLo; scutsMc30 +=sMtCutUp; scutsMc30 += sIsoCut30;

	//centrality
	std::vector <float> centBins;
	centBins.push_back(0.0);
	centBins.push_back(0.05);
	centBins.push_back(0.1);
	centBins.push_back(0.15);
	centBins.push_back(0.2);
	centBins.push_back(0.4);
	centBins.push_back(0.8);

	const int centralityBins = centBins.size()-1;

	//ncoll
	std::vector <float> ncoll;
	ncoll.push_back(1683.3); //0-5
	ncoll.push_back(1318.0); //5-10
	ncoll.push_back(1035.4); //10-15
	ncoll.push_back(811.2); //15-20
	ncoll.push_back(440.6); //20-40
	ncoll.push_back(77.8); //40-80

	double mbEvents = 68.7e6;


	TH1F* hmcWSet = new TH1F("hmcWSet","hmcWSet",100,0.0,100.0);
	TH1F* hmcWCent = new TH1F("hmcWCent","hmcWCent",100,0.0,100.0);

	TTree *treeMcWSet = (TTree*)fMcWSet->Get("tree");

	TTree *treeMcJ1Set = (TTree*)fMcJ1Set->Get("tree");
	TTree *treeMcJ2Set = (TTree*)fMcJ2Set->Get("tree");
	TTree *treeMcJ3Set = (TTree*)fMcJ3Set->Get("tree");

//	double eventWt = (623803*1.8770e-4+624874*8.2788e-6+639291*2.9426e-7)/(1.00e6*64.0e-3); //events taken from grid over events generated times relative cross-section

	//ratio of Jx cross-section to total pp cross-section
	double wtJsum = 1.8770e-4+8.2788e-6+2.9426e-7;
	//double wtJ1 = 1.8770e-4/wtJsum; double wtJ2 = 8.2788e-6/wtJsum; double wtJ3 = 2.9426e-7/wtJsum;
	double wtJ1 = 1.8770e-4/64.0e-3; double wtJ2 = 8.2788e-6/64.0e-3; double wtJ3 = 2.9426e-7/64.0e-3;
	double wtW = 2.8205e-9/64.0e-3;
	double minBiasLoPt = 5.02730000000000000e+04; //integral from 7-20 in mb stream

	for(int icent=0; icent<centralityBins; icent++){
		
		TString sAddCent = "&&centrality>"; sAddCent += centBins.at(icent); sAddCent += "&&centrality<"; sAddCent += centBins.at(icent+1);
		scutsMc += sAddCent; scutsMc10+= sAddCent; scutsMc20+= sAddCent; scutsMc30+= sAddCent;
		TString cutsW = scutsMc + "&&prompt==24";
		double arrCentWidth = centBins.at(icent+1)-centBins.at(icent);
		double scaleFactor = arrCentWidth*ncoll.at(icent)*mbEvents;

		TString sCentrality = "centrality>"; sCentrality+= centBins.at(icent); sCentrality+="&&centrality<"; sCentrality+=centBins.at(icent+1);

		std::cout << sCentrality << std::endl;
		std::cout << scutsMc  << std::endl;
		std::cout << cutsW  << std::endl; 

		//J1-J3

		//no cuts
		treeMcJ1Set->Draw("pt>>hmcJ1Set",sCentrality,"hf");
		treeMcJ2Set->Draw("pt>>hmcJ2Set",sCentrality,"hf");
		treeMcJ3Set->Draw("pt>>hmcJ3Set",sCentrality,"hf");

		treeMcJ1Set->Draw("centrality>>hmcJ1Cent",sCentrality,"hf");
		treeMcJ2Set->Draw("centrality>>hmcJ2Cent",sCentrality,"hf");
		treeMcJ3Set->Draw("centrality>>hmcJ3Cent",sCentrality,"hf");

		//muon yield per Jx mc event
		hmcQCDSet->Add(hmcJ1Set,hmcJ2Set,(wtJ1)/(hmcJ1Cent->Integral()), (wtJ2)/(hmcJ2Cent->Integral())); 
		hmcQCDSet->Add(hmcQCDSet,hmcJ3Set,1.0,(wtJ3)/(hmcJ3Cent->Integral()));
		TH1F* hmcQCDSetNoCuts = (TH1F*)hmcQCDSet->Clone("hmcQCDSetNoCuts");

		//muon yield per HI collision in N_minbias events
		hmcQCDSetNoCuts->Scale(scaleFactor);

		//isolation cut < 0.1
/*		treeMcJ1Set->Draw("pt>>hmcJ1Set",scutsMc10,"hf");
		treeMcJ2Set->Draw("pt>>hmcJ2Set",scutsMc10,"hf");
		treeMcJ3Set->Draw("pt>>hmcJ3Set",scutsMc10,"hf");
		hmcQCDSet->Add(hmcJ1Set,hmcJ2Set,(wtJ1)/(hmcJ1Cent->Integral()), (wtJ2)/(hmcJ2Cent->Integral())); 
		hmcQCDSet->Add(hmcQCDSet,hmcJ3Set,1.0,(wtJ3)/(hmcJ3Cent->Integral()));
		TH1F* hmcQCDSet10 = (TH1F*)hmcQCDSet->Clone("hmcQCDSet10");
		hmcQCDSet10->Scale(scaleFactor);

		//isolation cut < 0.2
		treeMcJ1Set->Draw("pt>>hmcJ1Set",scutsMc20,"hf");
		treeMcJ2Set->Draw("pt>>hmcJ2Set",scutsMc20,"hf");
		treeMcJ3Set->Draw("pt>>hmcJ3Set",scutsMc20,"hf");
		hmcQCDSet->Add(hmcJ1Set,hmcJ2Set,(wtJ1)/(hmcJ1Cent->Integral()), (wtJ2)/(hmcJ2Cent->Integral())); 
		hmcQCDSet->Add(hmcQCDSet,hmcJ3Set,1.0,(wtJ3)/(hmcJ3Cent->Integral()));
		TH1F* hmcQCDSet20 = (TH1F*)hmcQCDSet->Clone("hmcQCDSet20");
		hmcQCDSet20->Scale(scaleFactor);
*/
		//isolation cut < 0.3
		treeMcJ1Set->Draw("pt>>hmcJ1Set",scutsMc30,"hf");
		treeMcJ2Set->Draw("pt>>hmcJ2Set",scutsMc30,"hf");
		treeMcJ3Set->Draw("pt>>hmcJ3Set",scutsMc30,"hf");

		hmcQCDSet->Add(hmcJ1Set,hmcJ2Set,(wtJ1)/(hmcJ1Cent->Integral()), (wtJ2)/(hmcJ2Cent->Integral())); 
		std::cout << scutsMc30 << std::endl;
		hmcQCDSet->Add(hmcQCDSet,hmcJ3Set,1.0,(wtJ3)/(hmcJ3Cent->Integral()));
		TH1F* hmcQCDSet30 = (TH1F*)hmcQCDSet->Clone("hmcQCDSet30");
		hmcQCDSet30->Scale(scaleFactor);

		//Wmunu

		//no cuts
		treeMcWSet->Draw("pt>>hmcWCent",sCentrality+"&&prompt==24","hf");		
		TH1F* hmcWSetNoCuts = (TH1F*)hmcWCent->Clone("hmcWSetNoCuts");
		hmcWSetNoCuts->Scale(scaleFactor/(hmcWCent->Integral()));

		cutsW = scutsMc10 + "&&prompt==24";
		treeMcWSet->Draw("pt>>hmcWSet",cutsW,"hf");		
		TH1F* hmcWSet10 = (TH1F*)hmcWSet->Clone("hmcWSet10");
		hmcWSet10->Scale(scaleFactor/(hmcWCent->Integral()));

		cutsW = scutsMc20 + "&&prompt==24";
		treeMcWSet->Draw("pt>>hmcWSet",cutsW,"hf");		
		TH1F* hmcWSet20 = (TH1F*)hmcWSet->Clone("hmcWSet20");
		hmcWSet20->Scale(scaleFactor/(hmcWCent->Integral()));

		cutsW = scutsMc30 + "&&prompt==24";
		treeMcWSet->Draw("pt>>hmcWSet",cutsW,"hf");		
		TH1F* hmcWSet30 = (TH1F*)hmcWSet->Clone("hmcWSet30");
		hmcWSet30->Scale(scaleFactor/(hmcWCent->Integral()));

		TLegend* leg = new TLegend(0.648, 0.622, 0.918, 0.815);
		leg->SetTextFont(gStyle->GetTextFont());
		leg->SetTextSize(gStyle->GetTextSize());
		leg->SetBorderSize(0);
		leg->SetFillColor(0);
		leg->AddEntry(hmcQCDSetNoCuts,"J1-3","l");
		leg->AddEntry(hmcWSetNoCuts,"Wmunu","l");

		TString sCone = "#Delta R <"; sCone+=cone;
		TString sSel = ""; sSel += 100*centBins.at(icent); sSel.Remove(3); sSel += "-"; sSel +=100*centBins.at(icent+1); sSel.Remove(3);sSel +="%";  

		TLatex l;
		l.SetNDC();
		TCanvas* cpt = new TCanvas("cdatapt","cdatapt",600,600);
		hmcQCDSetNoCuts->GetXaxis()->SetTitle("p_{T}[GeV]"); 
		hmcQCDSetNoCuts->GetYaxis()->SetTitle("Events/GeV");
		hmcQCDSetNoCuts->SetLineColor(kRed);
		hmcWSetNoCuts->SetLineColor(kBlue);
		hmcQCDSetNoCuts->Draw("hf");
		hmcWSetNoCuts->Draw("hf same");
		l.DrawLatex(0.18,0.83,"no selection cuts");
		l.DrawLatex(0.713,0.5,sSel);
		leg->Draw();
		cpt->Print("cpt"+sSel+date+".pdf");

		TCanvas* cpt30 = new TCanvas("cdatapt30","cdatapt30",600,600);
		hmcQCDSet30->GetXaxis()->SetTitle("p_{T}[GeV]"); 
		hmcQCDSet30->GetYaxis()->SetTitle("Events/GeV");
		hmcQCDSet30->SetLineColor(kRed);
		hmcWSet30->SetLineColor(kBlue);
		hmcWSet30->Draw("hf ");
		hmcQCDSet30->Draw("hfsame");
		l.DrawLatex(0.18,0.83,sCone+",isolation cut < 0.3");
		l.DrawLatex(0.15,0.73,sSel);
		leg->Draw();
		std::cout << "QCD background with isolation cut < 0.3 = " << hmcQCDSet30->Integral(26,100)/hmcWSet30->Integral(26,100)*100 << "%" << std::endl;
		cpt30->Print("cpt30"+sCone+sSel+date+".pdf");

		//break;
		/*TCanvas* cpt20 = new TCanvas("cdatapt20","cdatapt20",600,600);
		hmcQCDSet20->GetXaxis()->SetTitle("p_{T}[GeV]"); 
		hmcQCDSet20->GetYaxis()->SetTitle("Events/GeV");
		hmcQCDSet20->SetLineColor(kRed);
		hmcWSet20->SetLineColor(kBlue);
		hmcQCDSet20->Draw("hf");
		hmcWSet20->Draw("hf same");
		l.DrawLatex(0.18,0.83,sCone+",isolation cut < 0.2");
		l.DrawLatex(0.15,0.73,sSel);
		leg->Draw();
		std::cout << "QCD background with isolation cut < 0.2 = " << hmcQCDSet20->Integral(26,100)/hmcWSet20->Integral(26,100)*100 << "%" << std::endl;
		cpt20->Print("cpt20"+sCone+sSel+date+".pdf");

		TCanvas* cpt10 = new TCanvas("cdatapt10","cdatapt10",600,600);
		hmcQCDSet10->GetXaxis()->SetTitle("p_{T}[GeV]"); 
		hmcQCDSet10->GetYaxis()->SetTitle("Events/GeV");
		hmcQCDSet10->SetLineColor(kRed);
		hmcWSet10->SetLineColor(kBlue);
		hmcQCDSet10->Draw("hf");
		hmcWSet10->Draw("hf same");
		l.DrawLatex(0.18,0.83,sCone+",isolation cut < 0.1");
		l.DrawLatex(0.15,0.73,sSel);
		leg->Draw();
		std::cout << "QCD background with isolation cut < 0.1 = " << hmcQCDSet10->Integral(26,100)/hmcWSet10->Integral(26,100)*100 << "%" << std::endl;
		cpt10->Print("cpt10"+sCone+sSel+date+".pdf");
*/

	}	
	
}
