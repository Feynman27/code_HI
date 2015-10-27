{
	TCanvas *cWpt = new TCanvas("cWpt","cWpt",600,600);
	TFile *fData = TFile::Open("HISingleMuon_hp_aug28_v2.root");
	TFile *fMC = TFile::Open("HISingleMuon_mcWmunuDataOverlay_aug28_v2.root");

	TH1F *hD = new TH1F("hD","hD",100,0,100);
	TH1F *hMC = new TH1F("hMC","hMC",100,0,100);

	fMC->cd();
	tree->Draw("ptW>>hMC","val>10&&pt>25.0&&nu_pt>20.0&&abs(scat)<4.0&&abs(eLoss)<0.5&&ZDY!=1&&centrality<=0.8&&mt>40.0&&mt<120.0&&prompt==24");
	fData->cd();
	tree->Draw("ptW>>hD","val>11&&pt>25.0&&nu_pt>20.0&&abs(scat)<4.0&&abs(eLoss)<0.5&&ZDY!=1&&centrality<=0.8&&mt>40.0&&mt<120.0");
	hMC->DrawNormalized("same",7498);

	TLegend* leg0 = new TLegend(0.75, 0.65, 0.9, 0.8);
	leg0->SetTextFont(gStyle->GetTextFont());
	leg0->SetTextSize(gStyle->GetTextSize());
	leg0->SetBorderSize(0);
	leg0->SetFillColor(0);
	leg0->AddEntry("hD","Data","l");
	leg0->AddEntry("hMC","MC","p");
	leg0->Draw();
	///eta
	TCanvas *cWeta = new TCanvas("cWeta","cWeta",600,600);
	TH1F *hDeta = new TH1F("hDeta","hDeta",60,-3.0,3.0);
	TH1F *hMCeta = new TH1F("hMCeta","hMCeta",60,-3.0,3.0);
 	fMC->cd();
	tree->Draw("etaW>>hMCeta","prompt==24&&val>10&&abs(eta)<2.5&&pt>25.0&&nu_pt>20.0&&abs(scat)<4.0&&abs(eLoss)<0.5&&ZDY!=1&&centrality<=0.8&&mt>40.0&&mt<120.0");
	fData->cd();
	tree->Draw("etaW>>hDeta","val>11&&pt>25.0&&nu_pt>20.0&&abs(scat)<4.0&&abs(eLoss)<0.5&&ZDY!=1&&centrality<=0.8&&mt>40.0&&mt<120.0");
	hMCeta->DrawNormalized("same",7498);

	TLegend* leg1 = new TLegend(0.75, 0.65, 0.9, 0.8);
	leg1->SetTextFont(gStyle->GetTextFont());
	leg1->SetTextSize(gStyle->GetTextSize());
	leg1->SetBorderSize(0);
	leg1->SetFillColor(0);
	leg1->AddEntry("hMCeta","MC","l");
	leg1->AddEntry("hDeta","Data","l");
	leg1->Draw();

}
