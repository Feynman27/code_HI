{
	gROOT->LoadMacro("AtlasUtils.C");
	TFile* fMC = new TFile("correctionFactorsW_allFloat.01.13.2013.root","READ");
	char cGrPlus[50],cGrMinus[50],cTitle[50], cTitle2[50],cGr[50];
	std::vector<double> etaBins;
	etaBins.push_back(0.00);
	etaBins.push_back(+0.25);
	etaBins.push_back(+0.50);
	etaBins.push_back(+0.75);
	etaBins.push_back(+1.00);
	etaBins.push_back(+1.25);
	etaBins.push_back(+1.50);
	etaBins.push_back(+1.75);
	etaBins.push_back(+2.00);
	etaBins.push_back(+2.25);
	etaBins.push_back(+2.50);

	const int nEtaBins = etaBins.size()-1;

	
  	TCanvas* c = new TCanvas("c","c",600,600);
  	TCanvas* c2 = new TCanvas("c2","c2",600,600);

	for(int igr = 0; igr<nEtaBins; igr++){	

		c->cd();
		sprintf(cGrPlus,"grCent102eta%i",igr);
		sprintf(cGrMinus,"grCent103eta%i",igr);
		sprintf(cGr,"grCent104eta%i",igr);
		sprintf(cTitle,"grCentPlusMinuseta%i",igr);
		sprintf(cTitle2,"grCentInclusiveEta%i",igr);

		TString sEta = "";  sEta+= etaBins.at(igr); sEta+="<|#eta|<"; sEta+=etaBins.at(igr+1);
		std::cout << sEta << std::endl;

		TGraphErrors* grPlus = (TGraphErrors*)fMC->Get(cGrPlus);
		TGraphErrors* grMinus = (TGraphErrors*)fMC->Get(cGrMinus);
		TGraphErrors* gr = (TGraphErrors*)fMC->Get(cGr);

		grPlus->SetMarkerColor(kRed);
		grMinus->SetMarkerColor(kBlue);

		grPlus->Draw("ape");
		grMinus->Draw("pe same");

		grPlus->GetYaxis()->SetTitle("C_{W^{#pm}}");
		myText(0.33,0.89, (Color_t)kBlack, (char*)sEta);
		TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
		leg->SetTextFont(gStyle->GetTextFont());
		leg->SetTextSize(gStyle->GetTextSize());
		leg->SetBorderSize(0);
		leg->SetFillColor(0);
		leg->AddEntry(grPlus,"#mu^{+}","pe");
		leg->AddEntry(grMinus,"#mu^{-}","pe");
		leg->Draw();

		TString sTitle = cTitle;
		c->Print(sTitle+".pdf");

		c2->cd();
		gr->Draw("ape");
		myText(0.33,0.89, (Color_t)kBlack, (char*)sEta);
		TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
		leg->SetTextFont(gStyle->GetTextFont());
		leg->SetTextSize(gStyle->GetTextSize());
		leg->SetBorderSize(0);
		leg->SetFillColor(0);
		leg->AddEntry(gr,"#mu^{#pm}","pe");
		leg->Draw();
		TString sTitle2 = cTitle2;
		c2->Print(sTitle2+".pdf");

	}
}
