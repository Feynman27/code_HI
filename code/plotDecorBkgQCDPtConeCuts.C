{
    TFile* _fEta0 = TFile::Open("background/fractionQCDEta.03.30.2013.root");
    TFile* _fEta1 = TFile::Open("systematics/fractionQCDEta_ptcone20ID3_systematics.04.02.2013.root");
    TFile* _fEta2 = TFile::Open("systematics/fractionQCDEta_ptcone30ID3_systematics.04.02.2013.root");

    TGraphErrors* grEta0 = (TGraphErrors*)_fEta0->Get("fractionQCDEta");
    TGraphErrors* grEta1 = (TGraphErrors*)_fEta1->Get("fractionQCDEta");
    TGraphErrors* grEta2 = (TGraphErrors*)_fEta2->Get("fractionQCDEta");

    TCanvas* cEta = new TCanvas("cEta","cEta",700,600);
	TH1F* hdummy = new TH1F("hdummy","hdummy",10,0.0,2.6);
	hdummy->GetXaxis()->SetRangeUser(0.0,2.6);
	hdummy->GetYaxis()->SetRangeUser(0.0,0.1);
	hdummy->GetYaxis()->SetTitle("f_{QCD}");
	hdummy->GetXaxis()->SetTitle("|#eta|");

    hdummy->Draw();
    grEta0->Draw("pesame");
    grEta1->SetMarkerColor(kBlue);
    grEta1->Draw("pesame");
    grEta2->SetMarkerColor(kRed);
    grEta2->Draw("pesame");

    TLegend* leg = new TLegend(0.2,0.6,0.5,0.9);
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    leg->AddEntry(grEta0,"Optimal","pe");
    leg->AddEntry(grEta1,"#uparrow isolation_{#mu}","pe");
    leg->AddEntry(grEta2,"#uparrow #Delta R","pe");
    leg->Draw();

    TFile* _fCent0 = TFile::Open("background/fractionQCDCent.03.30.2013.root");
    TFile* _fCent1 = TFile::Open("systematics/fractionQCDCent_ptcone20ID3_systematics.04.02.2013.root");
    TFile* _fCent2 = TFile::Open("systematics/fractionQCDCent_ptcone30ID3_systematics.04.02.2013.root");

    TGraphErrors* grCent0 = (TGraphErrors*)_fCent0->Get("fractionQCDCent");
    TGraphErrors* grCent1 = (TGraphErrors*)_fCent1->Get("fractionQCDCent");
    TGraphErrors* grCent2 = (TGraphErrors*)_fCent2->Get("fractionQCDCent");

    TCanvas* cCent = new TCanvas("cCent","cCent",700,600);
	TH1F* hdummy2 = new TH1F("hdummy2","hdummy2",10,0.0,1600.0);
	hdummy2->GetXaxis()->SetRangeUser(0.0,1600.0);
	hdummy2->GetYaxis()->SetRangeUser(0.0,0.1);
	hdummy2->GetYaxis()->SetTitle("f_{QCD}");
	hdummy2->GetXaxis()->SetTitle("#LT N_{coll} #GT");

    hdummy2->Draw();
    grCent0->Draw("pesame");
    grCent1->SetMarkerColor(kBlue);
    grCent1->Draw("pesame");
    grCent2->SetMarkerColor(kRed);
    grCent2->Draw("pesame");

    TLegend* leg = new TLegend(0.2,0.6,0.5,0.9);
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    leg->AddEntry(grCent0,"Optimal","pe");
    leg->AddEntry(grCent1,"#uparrow isolation_{#mu}","pe");
    leg->AddEntry(grCent2,"#uparrow #Delta R","pe");
    leg->Draw();

    /////////////////////////////
    //Z
    /////////////////////////////

    TFile* _fEtaZ0 = TFile::Open("background/fractionZEta.03.31.2013.root");
    TFile* _fEtaZ1 = TFile::Open("systematics/fractionZEta_ptcone20ID3_systematics.04.02.2013.root");
    TFile* _fEtaZ2 = TFile::Open("systematics/fractionZEta_ptcone30ID3_systematics.04.02.2013.root");

    TGraphErrors* grEta0 = (TGraphErrors*)_fEtaZ0->Get("fractionZEta");
    TGraphErrors* grEta1 = (TGraphErrors*)_fEtaZ1->Get("fractionZEta");
    TGraphErrors* grEta2 = (TGraphErrors*)_fEtaZ2->Get("fractionZEta");

    TCanvas* cEtaZ = new TCanvas("cEtaZ","cEtaZ",700,600);
	TH1F* hdummy = new TH1F("hdummy","hdummy",10,0.0,2.6);
	hdummy->GetXaxis()->SetRangeUser(0.0,2.6);
	hdummy->GetYaxis()->SetRangeUser(0.0,0.19);
	hdummy->GetYaxis()->SetTitle("f_{Z}");
	hdummy->GetXaxis()->SetTitle("|#eta|");

    hdummy->Draw();
    grEta0->Draw("pesame");
    grEta1->SetMarkerColor(kBlue);
    grEta1->Draw("pesame");
    grEta2->SetMarkerColor(kRed);
    grEta2->Draw("pesame");

    TLegend* leg = new TLegend(0.2,0.6,0.5,0.9);
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    leg->AddEntry(grEta0,"Optimal","pe");
    leg->AddEntry(grEta1,"#uparrow isolation_{#mu}","pe");
    leg->AddEntry(grEta2,"#uparrow #Delta R","pe");
    leg->Draw();

    TFile* _fCentZ0 = TFile::Open("background/fractionZCent.03.31.2013.root");
    TFile* _fCentZ1 = TFile::Open("systematics/fractionZCent_ptcone20ID3_systematics.04.02.2013.root");
    TFile* _fCentZ2 = TFile::Open("systematics/fractionZCent_ptcone30ID3_systematics.04.02.2013.root");

    TGraphErrors* grCent0 = (TGraphErrors*)_fCentZ0->Get("fractionZCent");
    TGraphErrors* grCent1 = (TGraphErrors*)_fCentZ1->Get("fractionZCent");
    TGraphErrors* grCent2 = (TGraphErrors*)_fCentZ2->Get("fractionZCent");

    TCanvas* cCentZ = new TCanvas("cCentZ","cCentZ",700,600);
	TH1F* hdummy2 = new TH1F("hdummy2","hdummy2",10,0.0,1600.0);
	hdummy2->GetXaxis()->SetRangeUser(0.0,1600.0);
	hdummy2->GetYaxis()->SetRangeUser(0.0,0.19);
	hdummy2->GetYaxis()->SetTitle("f_{Z}");
	hdummy2->GetXaxis()->SetTitle("#LT N_{coll} #GT");

    hdummy2->Draw();
    grCent0->Draw("pesame");
    grCent1->SetMarkerColor(kBlue);
    grCent1->Draw("pesame");
    grCent2->SetMarkerColor(kRed);
    grCent2->Draw("pesame");

    TLegend* leg = new TLegend(0.2,0.6,0.5,0.9);
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetTextSize(gStyle->GetTextSize());
    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    leg->AddEntry(grCent0,"Optimal","pe");
    leg->AddEntry(grCent1,"#uparrow isolation_{#mu}","pe");
    leg->AddEntry(grCent2,"#uparrow #Delta R","pe");
    leg->Draw();




    //delete cEta;
}
