#include "TEfficiency.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TLine.h"
#include "TLegend.h"

TEfficiency* writeEfficiency(TFile* fEff, const TH1& hPassed, const TH1& hTotal, TString sEff){
	if(TEfficiency::CheckConsistency(hPassed,hTotal)) {

		TEfficiency* pEff = 0;
		fEff->cd();
		pEff = new TEfficiency(hPassed,hTotal);
		pEff->Write(sEff);
		pEff->Draw("AP");
		return pEff;
	}
}
void dataCutFlow(){

	//TFile::Open("../HardProbesFiles/HISingleMuonHP.12.19.2012.root");
	TString fileNameIn = "/usatlas/u/tbales/scratch/HISingleMuonHardProbesData.04.17.2013";
	TFile* pFile = new TFile("dataCutRatios.root","recreate");
	TFile* fIn = new TFile(fileNameIn+".root", "READ");

	float mtmax = 400.0;float ptmax = 400.0; 
//	TString cuts = "abs(charge)==1&&val>11&&abs(eLoss)<0.5&&abs(scat)<4.0&&abs(eta)>0.1&&abs(eta)<2.4&&centrality>0&&centrality<0.8&&((EF_mu10_MSonly_EFFS_L1ZDC&&EF_mu10_MSonly_EFFS_L1ZDC_Matched20)||(EF_mu10_MSonly_EFFS_L1TE10&&EF_mu10_MSonly_EFFS_L1TE10_Matched20)||(EF_mu10_MSonly_EFFS_L1TE20&&EF_mu10_MSonly_EFFS_L1TE20_Matched20))" ;
        TString trig = "((EF_mu10_MSonly_EFFS_L1ZDC&&EF_mu10_MSonly_EFFS_L1ZDC_Matched20)||(EF_mu10_MSonly_EFFS_L1TE10&&EF_mu10_MSonly_EFFS_L1TE10_Matched20)||(EF_mu10_MSonly_EFFS_L1TE20&&EF_mu10_MSonly_EFFS_L1TE20_Matched20))";
        TString ps = "&&abs(eta)>0.1&&abs(eta)<2.4&&centrality>0.&&centrality<0.8&&abs(scat)<4.&&abs(eLoss)<0.5&&val>11&&abs(charge)==1.0";
        TString cuts = trig+ps;

        std::cout << "P.S.: " << cuts << std::endl;
        std::cout << std::endl;
	double ptLow=25.0;double mtcutLow = 40.0;
	double missPtCut = 25.0; double isoCut = 0.1;
	TString sMissPtUp = "&&nu_pt>"; sMissPtUp+=missPtCut;
	TString sIsoCut = "&&ptcone20ID3/pt<"; sIsoCut+=isoCut;
	//TString sIsoCut = "&&ptcone20/pt<"; sIsoCut+=0.3;
	TString sMtCutLo = "&&mt>"; sMtCutLo+=mtcutLow;

//	TString cutsZ = cuts+"&&ZDY==0";
	TString cutsIso = cuts + sIsoCut;
        std::cout << "PS+Iso: " << cutsIso << std::endl;
        std::cout << std::endl;
	TString cutsZ = cutsIso+"&&ZDY==0";
        std::cout << "PS+Iso+Zveto: " << cutsZ << std::endl;
        std::cout << std::endl;
	TString cutsMPT = cutsZ+sMissPtUp;
        std::cout << "PS+Iso+Zveto+mpt: " << cutsMPT << std::endl;
        std::cout << std::endl;
	TString cutsMt = cutsMPT + sMtCutLo;
        std::cout << "PS+Iso+Zveto+mpt+mt: " << cutsMt << std::endl;
        std::cout << std::endl;

	double xBins[61] ;
	
	double xlo = 0.0; double xhi = 50.0;  
	double nbins = 50;
	double binw = (xhi-xlo)/nbins;
	int ib = 0;
	for(double i=xlo ; i<xhi ; i+=binw){
		
			xBins[ib] = i;
			ib++;	
	}
	xlo = 50.0; xhi = 75.0; 
	nbins = 5;
	binw = (xhi-xlo)/nbins;
	for(double i=xlo ; i<xhi ; i+=binw){

		xBins[ib] = i;
		ib++;
	}	
	xlo = 75.0; xhi = 100.0; 
	nbins = 3;
	binw = (xhi-xlo)/nbins;
	for(double i=xlo ; i<xhi ; i+=binw){

		xBins[ib] = i;
		ib++;
	}
	xlo = 100.0; xhi = 400.0; 
	nbins = 3;
	binw = (xhi-xlo)/nbins;
	for(double i=xlo ; i<xhi ; i+=binw){

		xBins[ib] = i;
		ib++;
	}
	int binNum = sizeof(xBins)/sizeof(double) - 1;
        std::cout << "No. of bins: " << binNum << std::endl;

	TH1F* h0 = new TH1F("h0","h0",binNum,xBins);
	TH1F* h1 = new TH1F("h1","h1",binNum,xBins);
	TH1F* h2 = new TH1F("h2","h2",binNum,xBins);
	TH1F* h3 = new TH1F("h3","h3",binNum,xBins);
	TH1F* h4 = new TH1F("h4","h4",binNum,xBins);

	tree->Draw("pt>>h0",cuts,"pe");
	TH1F* h0c= (TH1F*)h0->Clone("h0c");
	h0c->Scale(1,"width");
	h0c->SetMarkerColor(kGreen);

	tree->Draw("pt>>h1",cutsIso,"pesame");
	TH1F* h1c= (TH1F*)h1->Clone("h1c");
	h1c->Scale(1,"width");
	h1c->SetMarkerColor(kYellow);

	tree->Draw("pt>>h4",cutsZ,"pesame");
	TH1F* h4c= (TH1F*)h4->Clone("h4c");
	h4c->Scale(1,"width");
	h4c->SetMarkerColor(kBlue);

	tree->Draw("pt>>h2",cutsMPT,"pesame");
	TH1F* h2c= (TH1F*)h2->Clone("h2c");
	h2c->Scale(1,"width");
	h2c->SetMarkerColor(kViolet+5);

	tree->Draw("pt>>h3",cutsMt,"pesame");
	TH1F* h3c= (TH1F*)h3->Clone("h3c");
	h3c->Scale(1,"width");
	h3c->SetMarkerColor(kCyan);

	TEfficiency* pEff0 =writeEfficiency(pFile,*h1c,*h0c,"pEff0");
	fIn->cd();
        TEfficiency* pEff1 = writeEfficiency(pFile,*h2c,*h0c,"pEff1");
	fIn->cd();
	TEfficiency* pEff2 =writeEfficiency(pFile,*h3c,*h0c,"pEff2");
	fIn->cd();
	TEfficiency* pEff3 =writeEfficiency(pFile,*h4c,*h0c,"pEff3");
	fIn->cd();
	
	TCanvas *cptData = new TCanvas("cptData","cptData",600,700);
    //frame1->SetFillColor(0);
    //frame1->Draw();
	//TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
	TPad *pad1 = new TPad("pad1","pad1",0.0,0.219,1.0,0.778);
	pad1->SetBottomMargin(0);
	pad1->Draw();
	pad1->cd();
    TFrame* frame1 = new TFrame(10.0,-0.958607,100.0,5.176091);
	h0c->GetXaxis()->SetTitle("p_{T}^{#mu}[GeV]"); 
	h0c->GetYaxis()->SetRangeUser(0.15,1.0e7); 
	h0c->GetXaxis()->SetRangeUser(10.0,100.0); 
	h0c->Draw("pe");
	h1c->Draw("pesame");
	h4c->Draw("pesame");
	h2c->Draw("pesame");
	h3c->Draw("pesame");
    
	cptData->SetLogy(); cptData->Update();

	TLegend* leg = new TLegend(0.659, 0.73, 0.9295, 0.8636);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->AddEntry("h0c","Pre-Selection","p");
	leg->AddEntry("h1c","Isolation_{#mu}","p");
	leg->AddEntry("h4c","Z veto","p");
	leg->AddEntry("h2c","#slash{p_{T}}","p");
	leg->AddEntry("h3c","m_{T}","p");
	leg->Draw();

	double ptCutLine = 25.0;
    TLine *line0 = new TLine(ptCutLine,0.01,ptCutLine,881201);
    line0->SetLineColor(kBlack); line0->SetLineStyle(kDashed);
    line0->SetLineWidth(2);
    line0->Draw();



	cptData->cd();
    //frame2->SetFillColor(0);
    //frame2->Draw();
	//TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
	TPad *pad2 = new TPad("pad2","pad2",0.0,0.0476,1,0.301);
	pad2->SetTopMargin(0);
	pad2->Draw();
	pad2->cd();
    TFrame* frame2 = new TFrame(9.873000,0.0,100.92,1.05);

	/*pEff0->GetPaintedGraph()->SetMinimum(0.0);
	pEff0->GetPaintedGraph()->SetMaximum(1.1);
	pEff0->GetPaintedGraph()->GetXaxis()->SetTitle("p_{T}^{#mu}[GeV]");	
	*/
	pEff0->SetMarkerColor(kYellow);
	pEff1->SetMarkerColor(kViolet+5);
	pEff2->SetMarkerColor(kCyan);
	pEff3->SetMarkerColor(kBlue);
	pEff0->Draw("ap");
	pEff3->Draw("psame");
	pEff1->Draw("psame");
	pEff2->Draw("psame");
	pEff0->Draw("sameaxis");

    cptData->Print("dataCutFlow.root");

    int binCtrlLo = h0->FindBin(10.0);
    int binSigLo = h4->FindBin(25.0);
    std::cout << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "S u m m a r y" << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "Number of preselected events                                   = " << h0->Integral(binCtrlLo,binNum) << std::endl; 
    std::cout << "Number of events after isolation                               = " << h1->Integral(binCtrlLo,binNum) << std::endl; 
    std::cout << "Number of events after Z veto                                  = " << h4->Integral(binCtrlLo,binNum) << std::endl; 
    std::cout << "Number of events at pt>25GeV                                   = " << h4->Integral(binSigLo,binNum) << std::endl; 
    std::cout << "Number of events at mpt>25GeV                                  = " << h2->Integral(binSigLo,binNum) << std::endl; 
    std::cout << "Number of events at mt>40GeV                                   = " << h3->Integral(binSigLo,binNum) << std::endl; 




}
