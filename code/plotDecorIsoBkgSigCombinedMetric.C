#include <ctime>
#include <iomanip>

void plotDecorIsoBkgSigCombinedMetric(){
	//TFile* fIn = TFile::Open("isoBkgSigEfficiencyCustomCut.03.18.2013.root");
	//TFile* fIn = TFile::Open("isoBkgSigEfficiencyCustomCut.03.20.2013.root");
	TFile* fIn = TFile::Open("isoBkgSigEfficiencyCustomCut.03.25.2013.root");
	//TFile* fInW = TFile::Open("isoEfficiencyDeltaRWmunu.03.13.2013.root");
	//TFile* fInW = TFile::Open("isoEfficiencyDeltaRQCD.03.13.2013.root");

        // current date/time based on current system
//        std::time_t now = time(0);
//	std::tm *ltm = localtime(&now);
   
//	std::cout << "month: " << ltm->tm_mon << " day: " << ltm->tm_mday << " year: " << 1900+ltm->tm_year << std::endl; exit(0); 
//        std::string sTime = "."; sTime += (char*)ltm->tm_mon; sTime+=".";  sTime+= (char*)ltm->tm_mday; sTime+="."; //sTime+=1900+ltm->tm_year;sTime+="."; 
//	std::cout << sTime << std::endl; exit(0);
        // convert now to string form
        //char* dt = ctime(&now);

	gROOT->LoadMacro("AtlasUtils.C");
    TGaxis::SetMaxDigits(2);
	//const int trkPtCut = 8;
	const int trkPtCut = 5;
	//TString trkPtLo[trkPtCut] = {"05","075","1","2","3","4","5","6"};
	TString trkPtLo[trkPtCut] = {"05","075","1","2","3"};
	TString titleTrkPtLo[trkPtCut] = {"0.5","0.75","1","2","3","4","5","6"};
	const unsigned int nCones = 6;
	int arrCone[nCones] = {10,15,20,30,40,50};
	int nGraphs = trkPtCut*nCones;
	//centrality
	std::vector <float> centBins;
	centBins.push_back(0.0);
	centBins.push_back(0.1);
//	centBins.push_back(0.2);
	centBins.push_back(0.4);
	centBins.push_back(0.8);

	const int nCentralityBins = centBins.size()-1;


	const int isoCutBins = 11;

	int markerColor[8] = {kRed,kCyan,kBlue,kGreen,kMagenta,kOrange+5,kAzure-9, kYellow};

      for(int icent=0; icent<nCentralityBins; icent++){

	TCanvas* cRatio = new TCanvas("cRatio","cRatio",700,700);
        //TLegend* leg2 = new TLegend(0.2,0.4,0.5,0.9);
	TLegend* leg2 = new TLegend(0.2,0.65,0.5,0.92);
        leg2->SetTextFont(gStyle->GetTextFont());
        leg2->SetTextSize(gStyle->GetTextSize());
	leg2->SetTextSize(0.03);
        leg2->SetBorderSize(0);
        leg2->SetFillColor(0);


	TGraph* grRatioOpt = new TGraph(1);
	TGraph* grRatioOptSyst = new TGraph(2);
    TList _grRatio10 ;
    TList _grRatio15 ;
    TList _grRatio20 ;
    TList _grRatio30 ;
    TList _grRatio40 ;
    TList _grRatio50 ;


	TH1F* hdummy2 = new TH1F("hdummy2","hdummy2",6,0.0,0.75);
	//hdummy2->GetYaxis()->SetRangeUser(0.1,1.5e7);
	if(icent==0)hdummy2->GetYaxis()->SetRangeUser(0.0,0.15);
	else if(icent==1)hdummy2->GetYaxis()->SetRangeUser(0.0,0.15);
	else if(icent==2)hdummy2->GetYaxis()->SetRangeUser(0.0,0.15);
	hdummy2->GetXaxis()->SetRangeUser(0.0,0.6);
	hdummy2->GetYaxis()->SetNoExponent();
	hdummy2->GetXaxis()->SetNoExponent();
    TGaxis::SetMaxDigits(2);
	//hdummy2->GetYaxis()->SetTitle("#epsilon_{S} N_{S}^{D}-#epsilon_{B} N_{B}^{D}");
	hdummy2->GetYaxis()->SetTitle("N_{B}^{D}/N_{S}^{D}");
	hdummy2->GetXaxis()->SetTitle("#Delta R");

	std::stringstream sCentLow, sCentUpp;
	sCentLow << std::setprecision(2) << centBins.at(icent)*100;
	sCentUpp << std::setprecision(2) << centBins.at(icent+1)*100;
	TString sCent = sCentLow.str(); sCent+="-";sCent+=sCentUpp.str();  sCent+="%";

	for(int itrk=0; itrk<trkPtCut; itrk++){
	  
	     TString sTitle = "p_{T}^{trk}>"; sTitle+=titleTrkPtLo[itrk]; sTitle+="GeV";
	     TString cTitle = "grIsoEffTrkPt"; cTitle+=trkPtLo[itrk]; cTitle+="cent"; cTitle+=icent;
	     TString cTitle2 = "grIsoEff";  cTitle2+="cent"; cTitle2+=icent;

	     TString sGr = "grIsoEffTrkPt"; sGr+=trkPtLo[itrk]; sGr+="cent"; sGr+=icent; sGr+="ConeRadius"; 
	     TString sGr10 = sGr; sGr10+=10; 
	     TString sGr15 = sGr; sGr15+=15; 
	     TString sGr20 = sGr; sGr20+=20; 
	     TString sGr30 = sGr; sGr30+=30; 
	     TString sGr40 = sGr; sGr40+=40; 
	     TString sGr50 = sGr; sGr50+=50; 

	     TString sGrRatio10 = sGr10; sGrRatio10+="SigBkgRatio"; 
	     std::cout << "Opening " << sGrRatio10 << std::endl;
	     TString sGrRatio15 = sGr15; sGrRatio15+="SigBkgRatio"; 
	     std::cout << "Opening " << sGrRatio15 << std::endl;
	     TString sGrRatio20 = sGr20; sGrRatio20+="SigBkgRatio"; 
	     std::cout << "Opening " << sGrRatio20 << std::endl;
	     TString sGrRatio30 = sGr30; sGrRatio30+="SigBkgRatio"; 
	     std::cout << "Opening " << sGrRatio30 << std::endl;
	     TString sGrRatio40 = sGr40; sGrRatio40+="SigBkgRatio"; 
	     std::cout << "Opening " << sGrRatio40 << std::endl;
	     TString sGrRatio50 = sGr50; sGrRatio50+="SigBkgRatio"; 
	     std::cout << "Opening " << sGrRatio50 << std::endl;

        double* ytemp10 = ((TGraph*)fIn->Get(sGrRatio10))->GetY(); 
        double* xtemp10 = ((TGraph*)fIn->Get(sGrRatio10))->GetX(); 
        double* ytemp15 = ((TGraph*)fIn->Get(sGrRatio15))->GetY(); 
        double* xtemp15 = ((TGraph*)fIn->Get(sGrRatio15))->GetX(); 
        double* ytemp20 = ((TGraph*)fIn->Get(sGrRatio20))->GetY(); 
        double* xtemp20 = ((TGraph*)fIn->Get(sGrRatio20))->GetX(); 
        double* ytemp30 = ((TGraph*)fIn->Get(sGrRatio30))->GetY(); 
        double* xtemp30 = ((TGraph*)fIn->Get(sGrRatio30))->GetX(); 
        double* ytemp40 = ((TGraph*)fIn->Get(sGrRatio40))->GetY(); 
        double* xtemp40 = ((TGraph*)fIn->Get(sGrRatio40))->GetX(); 
        double* ytemp50 = ((TGraph*)fIn->Get(sGrRatio50))->GetY(); 
        double* xtemp50 = ((TGraph*)fIn->Get(sGrRatio50))->GetX(); 
        //take only points up to isolation_mu<0.5 (i.e. 50% of muon pt)
        const int nPoints = 5;
        _grRatio10.Add(new TGraph(nPoints));
        _grRatio15.Add(new TGraph(nPoints));
        _grRatio20.Add(new TGraph(nPoints));
        _grRatio30.Add(new TGraph(nPoints));
        _grRatio40.Add(new TGraph(nPoints));
        _grRatio50.Add(new TGraph(nPoints));

        for(int igr=0; igr<nPoints; igr++){
         ( (TGraph*)_grRatio10.At(itrk))->SetPoint(igr,ytemp10[igr],xtemp10[igr]);
         ( (TGraph*)_grRatio15.At(itrk))->SetPoint(igr,ytemp15[igr],xtemp15[igr]);
         ( (TGraph*)_grRatio20.At(itrk))->SetPoint(igr,ytemp20[igr],xtemp20[igr]);
         ( (TGraph*)_grRatio30.At(itrk))->SetPoint(igr,ytemp30[igr],xtemp30[igr]);
         ( (TGraph*)_grRatio40.At(itrk))->SetPoint(igr,ytemp40[igr],xtemp40[igr]);
         ( (TGraph*)_grRatio50.At(itrk))->SetPoint(igr,ytemp50[igr],xtemp50[igr]);
        }

	    //plot all graphs on same plot for a given centrality bin
	     if(itrk==4){
             double xOpt = ((TGraph*)_grRatio20.At(itrk))->GetX()[1]; double yOpt = ((TGraph*)_grRatio20.At(itrk))->GetY()[1];
             std::cout << "Optimal bkg eff = " << yOpt << " Optimal sig eff = " << xOpt << std::endl;
             grRatioOpt->SetMarkerColor(kBlack);
             grRatioOpt->SetMarkerStyle(29);
             grRatioOpt->SetMarkerSize(2.8);
             grRatioOptSyst->SetMarkerStyle(34);
             grRatioOptSyst->SetMarkerSize(2.3);
             grRatioOptSyst->SetMarkerColor(14);
             grRatioOpt->SetPoint(0,xOpt,yOpt);
             xOpt = ((TGraph*)_grRatio30.At(itrk))->GetX()[1]; yOpt = ((TGraph*)_grRatio30.At(itrk))->GetY()[1];
             grRatioOptSyst->SetPoint(0,xOpt,yOpt);
             xOpt = ((TGraph*)_grRatio20.At(itrk))->GetX()[2]; yOpt = ((TGraph*)_grRatio20.At(itrk))->GetY()[2];
             grRatioOptSyst->SetPoint(1,xOpt,yOpt);
	     }	     

         cRatio->cd();
         if(itrk==0){
                hdummy2->Draw();
	            hdummy2->GetXaxis()->SetDecimals(kTRUE);
         }
    

         cRatio->Update();
	     ((TGraph*)_grRatio10.At(itrk))->SetMarkerColor(markerColor[itrk]);
	     ((TGraph*)_grRatio10.At(itrk))->SetMarkerStyle(8);
	     ((TGraph*)_grRatio15.At(itrk))->SetMarkerColor(markerColor[itrk]);
	     ((TGraph*)_grRatio15.At(itrk))->SetMarkerStyle(8);
	     ((TGraph*)_grRatio20.At(itrk))->SetMarkerColor(markerColor[itrk]);
	     ((TGraph*)_grRatio20.At(itrk))->SetMarkerStyle(8);
	     ((TGraph*)_grRatio30.At(itrk))->SetMarkerColor(markerColor[itrk]);
	     ((TGraph*)_grRatio30.At(itrk))->SetMarkerStyle(8);
	     ((TGraph*)_grRatio40.At(itrk))->SetMarkerColor(markerColor[itrk]);
	     ((TGraph*)_grRatio40.At(itrk))->SetMarkerStyle(8);
	     ((TGraph*)_grRatio50.At(itrk))->SetMarkerColor(markerColor[itrk]);
	     ((TGraph*)_grRatio50.At(itrk))->SetMarkerStyle(8);

	     ((TGraph*)_grRatio10.At(itrk))->Draw("psame");
	     ((TGraph*)_grRatio15.At(itrk))->Draw("psame");
	     ((TGraph*)_grRatio20.At(itrk))->Draw("psame");
	     ((TGraph*)_grRatio30.At(itrk))->Draw("psame");
	     ((TGraph*)_grRatio40.At(itrk))->Draw("psame");
	     ((TGraph*)_grRatio50.At(itrk))->Draw("psame");
         if(itrk>=4){
               grRatioOptSyst->Draw("psame");
               grRatioOpt->Draw("psame");
         }
	     cRatio->SetGrid();
	     cRatio->Update();
	     leg2->AddEntry((TGraph*)_grRatio10.At(itrk),sTitle,"p");

	     //delete c0;
	} //itrk

	cRatio->cd();
	//cRatio->SetLogy();
	leg2->Draw();
	myText(0.6,0.80,kBlack,(char*)sCent);
	cRatio->Print(cTitle2+"Ratio.03.25.2013.pdf");
	//delete cAll;
	//delete hdummy;
	//delete c0;
     }//icent
}
