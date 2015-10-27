#include<iostream>
#include<memory>
#include "TFile.h"
#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"

void fillHistos(const TString& pathName, const TString& fileName, TH1F* h){
    std::unique_ptr<TChain> tree(new TChain("treeW","treeW"));
    int nFiles=tree->Add(pathName+fileName);
    std::cout << "Filling  " << fileName << "... Number of files: " << nFiles << std::endl;
    std::cout << "Number of W-->mu candidates: " << tree->GetEntries() << std::endl;

    float _fcal;
    tree->SetBranchAddress("fcal",&_fcal);
    for(int i=0; i<tree->GetEntries(); ++i){
        tree->GetEntry(i);
        if(i%10000==0) std::cout << "Event: " << i << std::endl;
        h->Fill(_fcal);
    }//i
}

int main(){

/*    std::unique_ptr<TFile> _fData(new TFile); 

    std::unique_ptr<TFile> _fWmin_pp_A(new TFile); 
    std::unique_ptr<TFile> _fWmin_pn_A(new TFile); 
    std::unique_ptr<TFile> _fWmin_pp_B(new TFile); 
    std::unique_ptr<TFile> _fWmin_pn_B(new TFile); 

    std::unique_ptr<TFile> _fWplus_pp_A(new TFile); 
    std::unique_ptr<TFile> _fWplus_pn_A(new TFile); 
    std::unique_ptr<TFile> _fWplus_pp_B(new TFile); 
    std::unique_ptr<TFile> _fWplus_pn_B(new TFile); 
*/
    const float binUpp = 150.0;
    const float binLow = 2.555; //90%

    std::unique_ptr<TFile> _fOut(new TFile("fcalHistos.root","recreate") );
    TH1F* _hData=new TH1F("_hData","_hData",binUpp,binLow,binUpp) ;
    _hData->Sumw2();
    TH1F* _hWmin_pp_A=new TH1F("_hWmin_pp_A","_hWmin_pp_A",binUpp,binLow,binUpp) ;
    _hWmin_pp_A->Sumw2();
    TH1F* _hWmin_pp_B=new TH1F("_hWmin_pp_B","_hWmin_pp_B",binUpp,binLow,binUpp) ;
    _hWmin_pp_B->Sumw2();
    TH1F* _hWmin_pn_A=new TH1F("_hWmin_pn_A","_hWmin_pn_A",binUpp,binLow,binUpp) ;
    _hWmin_pn_A->Sumw2();
    TH1F* _hWmin_pn_B=new TH1F("_hWmin_pn_B","_hWmin_pn_B",binUpp,binLow,binUpp) ;
    _hWmin_pn_B->Sumw2();
    TH1F* _hWplus_pp_A=new TH1F("_hWplus_pp_A","_hWplus_pp_A",binUpp,binLow,binUpp) ;
    _hWplus_pp_A->Sumw2();
    TH1F* _hWplus_pp_B=new TH1F("_hWplus_pp_B","_hWplus_pp_B",binUpp,binLow,binUpp) ;
    _hWplus_pp_B->Sumw2();
    TH1F* _hWplus_pn_A=new TH1F("_hWplus_pn_A","_hWplus_pn_A",binUpp,binLow,binUpp) ;
    _hWplus_pn_A->Sumw2();
    TH1F* _hWplus_pn_B=new TH1F("_hWplus_pn_B","_hWplus_pn_B",binUpp,binLow,binUpp) ;
    _hWplus_pn_B->Sumw2();
    // sum
    TH1F* _hMc_sum=new TH1F("_hMc_sum","_hMc_sum",binUpp,binLow,binUpp) ;
    _hMc_sum->Sumw2();
    TH1F* _hMc_sum_scaled=new TH1F("_hMc_sum_scaled","_hMc_sum_scaled",binUpp,binLow,binUpp) ;
    _hMc_sum_scaled->Sumw2();
    TH1F* _hRatio=new TH1F("_hRatio","_hRatio",binUpp,binLow,binUpp) ;
    _hRatio->Sumw2();

    TString baseString = "/usatlas/u/tbales/scratch/Data2013Hip/";
    fillHistos(baseString,"ntupleW.07.24.2014.root",_hData);
    fillHistos(baseString,"ntupleWminmunu_pn.A.PowhegPythia8_AU2CT10.09.24.2014.root",_hWmin_pn_A);
    fillHistos(baseString,"ntupleWminmunu_pn.B.PowhegPythia8_AU2CT10.09.24.2014.root",_hWmin_pn_B);
    fillHistos(baseString,"ntupleWminmunu_pp.A.PowhegPythia8_AU2CT10.09.24.2014.root",_hWmin_pp_A);
    fillHistos(baseString,"ntupleWminmunu_pp.B.PowhegPythia8_AU2CT10.09.24.2014.root",_hWmin_pp_B);
    fillHistos(baseString,"ntupleWplusmunu_pn.A.PowhegPythia8_AU2CT10.09.24.2014.root",_hWplus_pn_A);
    fillHistos(baseString,"ntupleWplusmunu_pn.B.PowhegPythia8_AU2CT10.09.24.2014.root",_hWplus_pn_B);
    fillHistos(baseString,"ntupleWplusmunu_pp.A.PowhegPythia8_AU2CT10.09.24.2014.root",_hWplus_pp_A);
    fillHistos(baseString,"ntupleWplusmunu_pp.B.PowhegPythia8_AU2CT10.09.24.2014.root",_hWplus_pp_B);


	// Weight pp and np by collisions frequency
	double ppFreq = 82.0/208.0;
	double npFreq = 126.0/208.0;
    _hWmin_pn_A->Add(_hWmin_pn_B);
    _hWmin_pp_A->Add(_hWmin_pp_B);
    _hWmin_pn_A->Add(_hWmin_pn_A,_hWmin_pp_A,npFreq,ppFreq);

    _hWplus_pn_A->Add(_hWplus_pn_B);
    _hWplus_pp_A->Add(_hWplus_pp_B);
    _hWplus_pn_A->Add(_hWplus_pn_A,_hWplus_pp_A,npFreq,ppFreq);

    // W+/W- ~ 1.3
    _hMc_sum->Add(_hWmin_pn_A,_hWplus_pn_A,1.0,1.3);
    
    // norm to data
    _hMc_sum->Scale(_hData->Integral()/_hMc_sum->Integral());

    _fOut->cd();
    _hData->Write();
    _hMc_sum->Write();


/*    for(int ifcal=0; ifcal<nBins; ++ifcal){
        double fcal_low = fcal_bins[ifcal];
        double fcal_upp = fcal_bins[ifcal+1];

        int ibin_low = _hData->FindBin(fcal_low); 
        int ibin_upp = _hData->FindBin(fcal_upp); 

        double scale_factor = _hData->Integral(ibin_low,ibin_upp)/_hMc_sum->Integral(ibin_low,ibin_upp);

        std::cout << "Scale factors in fcal bin: " << ifcal << " = " << scale_factor << std::endl;
        std::cout << "Will be applied to mc histo bins: " << ibin_low+1 << " - " << ibin_upp << std::endl;
        for(int ibin=ibin_low+1; ibin<=ibin_upp; ++ibin){
            double fcal_scaled= scale_factor*_hMc_sum->GetBinContent(ibin);
            _hMc_sum_scaled->SetBinContent(ibin,fcal_scaled);
        }//ibin
    }//ifcal
*/
    _hRatio->Divide(_hData,_hMc_sum);

    TF1* _fit = new TF1("_fit","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
    _fit->SetParameter(0,0.739155);
    _fit->SetParameter(1,-0.0380782);
    _fit->SetParameter(2,0.00254976);
    _fit->SetParameter(3,-3.85359e-05);
    _fit->SetParameter(4,1.78379e-07);

    _hRatio->Fit(_fit);
    for(int ibin=1; ibin<=_hMc_sum->GetNbinsX(); ++ibin){
        double x = _hMc_sum->GetBinCenter(ibin);
        double scale_factor = _fit->Eval(x);
        double fcal_scaled = scale_factor*_hMc_sum->GetBinContent(ibin);
        _hMc_sum_scaled->SetBinContent(ibin,fcal_scaled);
    }//ibin
    
    std::vector<double> fcal_bins;
    //fcal_bins.push_back(_hData->GetXaxis()->GetXmax()); //0-5%
    fcal_bins.push_back(100.0); //0-5%
    fcal_bins.push_back(66.03); //5-10%
    fcal_bins.push_back(53.74); //10-20%
    fcal_bins.push_back(40.04); //20-30%
    fcal_bins.push_back(31.07); //30-40%
    fcal_bins.push_back(24.10); //40-60%
    fcal_bins.push_back(13.41); //60-90%
    fcal_bins.push_back(2.555);
    const int nBins = fcal_bins.size()-1;

    std::vector<double> vCent;
    vCent.push_back(0.00);
    vCent.push_back(0.05);
    vCent.push_back(0.10);
    vCent.push_back(0.20);
    vCent.push_back(0.30); // added to match Z
    vCent.push_back(0.40);
    vCent.push_back(0.60); // added to match Z
    vCent.push_back(0.90);

    TGraph* grRatioCent = new TGraph(nBins);
    for(int ifc=0; ifc<nBins; ++ifc){
        double binCenter = (fcal_bins[ifc]+fcal_bins[ifc+1])/2.0;
        double scale_factor  = _fit->Eval(binCenter);
        double centrality = (vCent[ifc]+vCent[ifc+1])/2.0;
        grRatioCent->SetPoint(ifc,centrality,scale_factor);
    }//ifc

    _hMc_sum_scaled->Write();
    _hRatio->Write();
    grRatioCent->Write("grScaleFactors");

    delete grRatioCent;
    delete _hRatio;
    delete _hMc_sum_scaled;
    delete _hWmin_pn_A;
    delete _hWmin_pn_B;
    delete _hWmin_pp_A;
    delete _hWmin_pp_B;
    delete _hWplus_pn_A;
    delete _hWplus_pn_B;
    delete _hWplus_pp_A;
    delete _hWplus_pp_B;
    delete _hData;
    delete _hMc_sum;

    return 0;
}
