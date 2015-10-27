#include <iomanip>

void plotDecorEffPlotter()
{
    bool drawOnlyCw = true;
    TString date = "NuEtaCutAbsEta9bins_06_08_2013";
//    TString date = "RealEta19Bins_05_26_2013";
    //TString date = "RealEta9Bins_05_14_2013";
    ///Cw evolution 
    TFile* _file0 =
//    TFile::Open("CorrectionFactorFiles/correctionFactorsWEtaCent_9VarEtaBins6CentBins2Charges.05.12.2013.root");
    //TFile::Open("CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges.05.03.2013.root");
//    TFile::Open("CorrectionFactorFiles/correctionFactorsWEtaCent_19EtaBinsAbsEta6CentBins2Charges.05.12.2013.root");
//    TFile::Open("CorrectionFactorFiles/correctionFactorsWEtaCent_19EtaBinsNoAbsEta6CentBins2Charges.05.12.2013.root");
//    TFile::Open("CorrectionFactorFiles/correctionFactorsW_binbybin_9EtaBinsNoAbsEta6CentBins2Charges.05.14.2013.root");
//    TFile::Open("CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges.05.25.2013.root");
//    TFile::Open("CorrectionFactorFiles/correctionFactorsWEtaCent_19EtaBinsNoAbsEta6CentBins2Charges.05.26.2013.root");
    ///with neutrino eta gen cut
    TFile::Open("CorrectionFactorFiles/correctionFactorsW_binbybin_9VarEtaBins6CentBins2Charges.06.08.2013.root");

    std::vector <double> etaBins;
    etaBins.push_back(0.0);

	etaBins.push_back(0.025);
	etaBins.push_back(0.05);
	etaBins.push_back(0.075);
	etaBins.push_back(0.1);
	etaBins.push_back(0.2);
	etaBins.push_back(0.3);
	etaBins.push_back(0.4);
	etaBins.push_back(0.5);
	etaBins.push_back(0.6);
	etaBins.push_back(0.7);
	etaBins.push_back(0.8);
	etaBins.push_back(0.9);
	etaBins.push_back(1.0);
	etaBins.push_back(1.05);
	etaBins.push_back(1.1);
	etaBins.push_back(1.15);
	etaBins.push_back(1.2);
	etaBins.push_back(1.25);
	etaBins.push_back(1.3);
	etaBins.push_back(1.35);
	etaBins.push_back(1.4);
	etaBins.push_back(1.5);
	etaBins.push_back(1.6);
	etaBins.push_back(1.7);
	etaBins.push_back(1.8);
	etaBins.push_back(1.9);
	etaBins.push_back(1.925);
	etaBins.push_back(1.950);
	etaBins.push_back(1.975);
	etaBins.push_back(2.0);
	etaBins.push_back(2.05);
	etaBins.push_back(2.1);
	etaBins.push_back(2.15);
	etaBins.push_back(2.2);
	etaBins.push_back(2.25);
	etaBins.push_back(2.3);
	etaBins.push_back(2.35);
    
    etaBins.push_back(2.4);
    const int nEtaBins = etaBins.size()-1;
    std::cout << "nEtaBins : " << nEtaBins << std::endl;

    //centrality
    std::vector <double> centBins;
    centBins.push_back(0.0);
    centBins.push_back(0.05);
    centBins.push_back(0.1);
    centBins.push_back(0.15);
    centBins.push_back(0.2);
    centBins.push_back(0.4);
    centBins.push_back(0.8);

    const int nCentralityBins = centBins.size()-1;
    std::cout << "nCentralityBins : " << nCentralityBins << std::endl;

    //Npart
    std::vector <double> npartBins;
    npartBins.push_back(382.16);//0-5
    npartBins.push_back(330.26);//5-10
    npartBins.push_back(281.88);//10-15
    npartBins.push_back(239.52);//15-20
    npartBins.push_back(157.83);//20-40
    npartBins.push_back(45.93);//40-80

    const int npartNBins = npartBins.size();

    TList _grPlus0 ;
    TList _grPlus1 ;
    TList _grPlus2 ;
    TList _grPlus3 ;
    TList _grPlus4 ;
    TList _grPlus5 ;

    TList _grMinus0 ;
    TList _grMinus1 ;
    TList _grMinus2 ;
    TList _grMinus3 ;
    TList _grMinus4 ;
    TList _grMinus5 ;

    TList _canvasPlus,_canvasMinus;
    
    gROOT->LoadMacro("AtlasUtils.C");

    TH1F* hdummy = new TH1F("hdummy","hdummy",20,0.0,2.5);
    //TH1F* hdummy = new TH1F("hdummy","hdummy",40,-2.5,2.5);
    hdummy->GetXaxis()->SetRangeUser(0.0,2.45);
    //hdummy->GetXaxis()->SetRangeUser(-2.45,2.45);
    hdummy->GetXaxis()->SetTitle("|#eta|");
    //hdummy->GetXaxis()->SetTitle("#eta");
    hdummy->GetYaxis()->SetRangeUser(0.0,1.45);
    //hdummy->GetYaxis()->SetTitle("N_{k}/N_{gen}");
    hdummy->GetYaxis()->SetTitle("#epsilon_{W#rightarrow#mu}");

    for(int icent=0; icent<nCentralityBins; ++icent){
        
        TString sRecEffEta0 = "grWmunuRecNoQualityWselNoIsoNoZvetoEtaDistro";
        //TString sRecEffEta0 = "grWmunuRecGenMatchedEtaDistro";
        TString sRecEffEtaPlus0 = sRecEffEta0+"PlusCent"; sRecEffEtaPlus0+=icent;
        TString sRecEffEtaMinus0 = sRecEffEta0+"MinusCent"; sRecEffEtaMinus0+=icent;

        _grPlus0.Add( (TGraphErrors*)_file0->Get(sRecEffEtaPlus0) );
        _grMinus0.Add( (TGraphErrors*)_file0->Get(sRecEffEtaMinus0) );

        TString sRecEffEta1 = "grWmunuRecNoQualityWselNoIsoEtaDistro";
        //TString sRecEffEta1 = "grWmunuRecPreSelectedEtaDistro";
        TString sRecEffEtaPlus1 = sRecEffEta1+"PlusCent"; sRecEffEtaPlus1+=icent;
        TString sRecEffEtaMinus1 = sRecEffEta1+"MinusCent"; sRecEffEtaMinus1+=icent;

        _grPlus1.Add( (TGraphErrors*)_file0->Get(sRecEffEtaPlus1) );
        _grMinus1.Add( (TGraphErrors*)_file0->Get(sRecEffEtaMinus1) );

        TString sRecEffEta2 = "grWmunuRecNoQualityWselEtaDistro";
        //TString sRecEffEta2 = "grWmunuRecPreSelectedPtCutEtaDistro";
        TString sRecEffEtaPlus2 = sRecEffEta2+"PlusCent"; sRecEffEtaPlus2+=icent;
        TString sRecEffEtaMinus2 = sRecEffEta2+"MinusCent"; sRecEffEtaMinus2+=icent;

        _grPlus2.Add( (TGraphErrors*)_file0->Get(sRecEffEtaPlus2) );
        _grMinus2.Add( (TGraphErrors*)_file0->Get(sRecEffEtaMinus2) );

        TString sRecEffEta3 = "grWmunuRecHiQualityWselNoIsoNoZvetoEtaDistro";
        //TString sRecEffEta3 = "grWmunuRecPreSelWselNoIsoNoZvetoEtaDistro";
        TString sRecEffEtaPlus3 = sRecEffEta3+"PlusCent"; sRecEffEtaPlus3+=icent;
        TString sRecEffEtaMinus3 = sRecEffEta3+"MinusCent"; sRecEffEtaMinus3+=icent;

        _grPlus3.Add( (TGraphErrors*)_file0->Get(sRecEffEtaPlus3) );
        _grMinus3.Add( (TGraphErrors*)_file0->Get(sRecEffEtaMinus3) );

        TString sRecEffEta4 = "grWmunuRecHiQualityWselNoIsoNoZvetoEtaDistro";
        //TString sRecEffEta4 = "grWmunuRecPreSelWselNoIsoEtaDistro";
        TString sRecEffEtaPlus4 = sRecEffEta4+"PlusCent"; sRecEffEtaPlus4+=icent;
        TString sRecEffEtaMinus4 = sRecEffEta4+"MinusCent"; sRecEffEtaMinus4+=icent;

        _grPlus4.Add( (TGraphErrors*)_file0->Get(sRecEffEtaPlus4) );
        _grMinus4.Add( (TGraphErrors*)_file0->Get(sRecEffEtaMinus4) );

        TString sRecEffEta5 = "grWmunuRecHiQualityWselEtaDistro";
        //TString sRecEffEta5 = "grWmunuRecPreSelWselEtaDistro";
        TString sRecEffEtaPlus5 = sRecEffEta5+"PlusCent"; sRecEffEtaPlus5+=icent;
        TString sRecEffEtaMinus5 = sRecEffEta5+"MinusCent"; sRecEffEtaMinus5+=icent;

        _grPlus5.Add( (TGraphErrors*)_file0->Get(sRecEffEtaPlus5) );
        _grMinus5.Add( (TGraphErrors*)_file0->Get(sRecEffEtaMinus5) );

         ((TGraphErrors*)_grPlus0.At(icent))->SetMarkerColor(kOrange+2);
         ((TGraphErrors*)_grPlus1.At(icent))->SetMarkerColor(kBlue);
         ((TGraphErrors*)_grPlus2.At(icent))->SetMarkerColor(kRed);
         ((TGraphErrors*)_grPlus3.At(icent))->SetMarkerColor(kMagenta);
         ((TGraphErrors*)_grPlus3.At(icent))->SetMarkerSize(1.5);
         ((TGraphErrors*)_grPlus4.At(icent))->SetMarkerColor(kAzure+10);
         ((TGraphErrors*)_grPlus5.At(icent))->SetMarkerColor(kSpring);

         ((TGraphErrors*)_grMinus0.At(icent))->SetMarkerColor(kOrange+2);
         ((TGraphErrors*)_grMinus1.At(icent))->SetMarkerColor(kBlue);
         ((TGraphErrors*)_grMinus2.At(icent))->SetMarkerColor(kRed);
         ((TGraphErrors*)_grMinus3.At(icent))->SetMarkerColor(kMagenta);
         ((TGraphErrors*)_grMinus3.At(icent))->SetMarkerSize(1.5);
         ((TGraphErrors*)_grMinus4.At(icent))->SetMarkerColor(kAzure+10);
         ((TGraphErrors*)_grMinus5.At(icent))->SetMarkerColor(kSpring);

         //TLegend* leg = new TLegend(0.19,0.18,0.46,0.28);
         //TLegend* leg = new TLegend(0.21,0.8,0.49,0.9);
         TLegend* leg = new TLegend(0.2,0.76,0.49,0.9);
         leg->SetTextFont(gStyle->GetTextFont());
         //leg->SetTextSize(0.042);
         leg->SetTextSize(0.032);
         leg->SetBorderSize(0);
         leg->SetFillColor(0);
         
         leg->AddEntry((TGraphErrors*)_grPlus0.At(icent),"k_{0}:|#eta|<2.4,p_{T}^{#mu,rec}>25,#slash{p_{T}}>25,m_{T}>40","pe");
//         leg->AddEntry((TGraphErrors*)_grPlus1.At(icent),"k_{1}:k_{0}+Z veto","pe");
//         leg->AddEntry((TGraphErrors*)_grPlus2.At(icent),"k_{2}:k_{1}+isolation","pe");

/*         leg->AddEntry((TGraphErrors*)_grPlus0.At(icent),"k_{0}:rec-gen matched","pe");
         leg->AddEntry((TGraphErrors*)_grPlus1.At(icent),"k_{1}:k_{0}+PS","pe");
         leg->AddEntry((TGraphErrors*)_grPlus2.At(icent),"k_{2}:k_{1}+p_{T}^{#mu,rec}>25","pe");
*/         
         //TLegend* leg2 = new TLegend(0.57,0.18,0.83,0.28);
         TLegend* leg2 = new TLegend(0.5,0.8,0.77,0.9);
         leg2->SetTextFont(gStyle->GetTextFont());
         leg2->SetTextSize(0.032);
         leg2->SetBorderSize(0);
         leg2->SetFillColor(0);


         leg->AddEntry((TGraphErrors*)_grPlus3.At(icent),"k_{0}:k_{0}+PS","pe");
         leg->AddEntry((TGraphErrors*)_grPlus4.At(icent),"k_{1}:k_{0}+Z veto","pe");
         leg->AddEntry((TGraphErrors*)_grPlus5.At(icent),"k_{2}:k_{1}+isolation","pe");

/*         leg->AddEntry((TGraphErrors*)_grPlus3.At(icent),"k_{1}:k_{0}+#slash{p_{T}}>25,m_{T}>40","pe");
         leg->AddEntry((TGraphErrors*)_grPlus4.At(icent),"k_{2}:k_{1}+Z veto","pe");
         leg->AddEntry((TGraphErrors*)_grPlus5.At(icent),"k_{3}:k_{2}+isolation","pe");
*/

         std::stringstream sCentLow, sCentUpp;
         sCentLow << std::setprecision(2) << centBins.at(icent)*100;
         sCentUpp << std::setprecision(2) << centBins.at(icent+1)*100;
         TString sCent = sCentLow.str(); sCent+="-";sCent+=sCentUpp.str();  sCent+="%";

         TLatex l;
         l.SetNDC();
         _canvasPlus.Add(new TCanvas());
         ((TCanvas*)_canvasPlus.At(icent))->cd();
         hdummy->Draw();
         if(!drawOnlyCw){
         ((TGraphErrors*)_grPlus0.At(icent))->Draw("pesame");
           //((TGraphErrors*)_grPlus1.At(icent))->Draw("pesame");
           //((TGraphErrors*)_grPlus2.At(icent))->Draw("pesame");
           ((TGraphErrors*)_grPlus3.At(icent))->Draw("pesame");
           ((TGraphErrors*)_grPlus4.At(icent))->Draw("pesame");
         }
         ((TGraphErrors*)_grPlus5.At(icent))->Draw("pesame");
         l.SetTextSize(0.038);
         l.DrawLatex(0.42,0.73,"#epsilon_{W#rightarrow#mu} =
          #frac{N_{k}}{N_{gen}#left{p_{T}^{#mu,gen}>25,p_{T}^{#nu,gen}>25,m_{T}^{#mu#nu}>40#right}}");
         l.Draw();
         leg->Draw();
         //leg2->Draw();
         //myText(0.8,0.7,kBlack,(char*)"#mu^{+}");
         //myText(0.77,0.67,kBlack,(char*)sCent);
         myText(0.85,0.81,kBlack,(char*)"#mu^{+}");
         myText(0.83,0.86,kBlack,(char*)sCent);
         l.SetTextSize(0.038);
         l.DrawLatex(0.42,0.73,"#epsilon_{W#rightarrow#mu} =
          #frac{N_{k}}{N_{gen}#left{p_{T}^{#mu,gen}>25,p_{T}^{#nu,gen}>25,m_{T}^{#mu#nu}>40#right}}");
         l.Draw();

         ((TCanvas*)_canvasPlus.At(icent))->Draw();
         //TString sNamePlus = "efficienciesPlusCent";sNamePlus+=icent;sNamePlus+=date;
         TString sNamePlus = "efficienciesTotalPlusCent";sNamePlus+=icent;sNamePlus+=date;
         ((TCanvas*)_canvasPlus.At(icent))->Print(sNamePlus+".pdf");
         ((TCanvas*)_canvasPlus.At(icent))->Print(sNamePlus+".root");

         _canvasMinus.Add(new TCanvas());
         ((TCanvas*)_canvasMinus.At(icent))->cd();
         hdummy->Draw();
         if(!drawOnlyCw){
         ((TGraphErrors*)_grMinus0.At(icent))->Draw("pesame");
           //((TGraphErrors*)_grMinus1.At(icent))->Draw("pesame");
           //((TGraphErrors*)_grMinus2.At(icent))->Draw("pesame");
           ((TGraphErrors*)_grMinus3.At(icent))->Draw("pesame");
           ((TGraphErrors*)_grMinus4.At(icent))->Draw("pesame");
         }
         ((TGraphErrors*)_grMinus5.At(icent))->Draw("pesame");
         l.SetTextSize(0.038);
         l.DrawLatex(0.42,0.73,"#epsilon_{W#rightarrow#mu} =
          #frac{N_{k}}{N_{gen}#left{p_{T}^{#mu,gen}>25,p_{T}^{#nu,gen}>25,m_{T}^{#mu#nu}>40#right}}");
         l.Draw();
         leg->Draw();
         //leg2->Draw();
         ((TCanvas*)_canvasMinus.At(icent))->Draw();
         //myText(0.8,0.7,kBlack,(char*)"#mu^{-}");
         //myText(0.77,0.67,kBlack,(char*)sCent);
         myText(0.85,0.81,kBlack,(char*)"#mu^{-}");
         myText(0.83,0.86,kBlack,(char*)sCent);
         //myText(0.254,0.854,kBlack,(char*)"#epsilon_{W#rightarrow#mu}=#frac{N_{k}}{N_{gen}#left{p_{T}^{#mu,gen}>25,p_{T}^{#nu,gen}>25,m_{T}^{#mu#nu}>40#right}}");

         //TString sNameMinus = "efficienciesMinusCent";sNameMinus+=icent;sNameMinus+=date;
         TString sNameMinus = "efficienciesTotalMinusCent";sNameMinus+=icent;sNameMinus+=date;
         ((TCanvas*)_canvasMinus.At(icent))->Print(sNameMinus+".pdf");
         ((TCanvas*)_canvasMinus.At(icent))->Print(sNameMinus+".root");

         //delete leg;
         //delete leg2;
         //delete _canvasPlus.At(icent);
         //delete _canvasMinus.At(icent);
         //break;
    }//icent

    //delete hdummy;
}
