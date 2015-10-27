#include "RooBinning.h"
using namespace RooFit;

void toyBinning(){

    double mtmax=300.0;
    RooBinning b = RooBinning(0.0,mtmax);
    b.addUniform(30,0.0,120.0);
    b.addUniform(10,120.0,200.0);
    b.addUniform(2,200.0,mtmax);
    double* xBins = b.array();
    int nBins=42;
    for(int i=0; i<nBins+1; ++i){
        std::cout << xBins[i] << std::endl;
    }
    //TH1F* h = new TH1F("h","h",nBins,xBins);
    TH1F* h = new TH1F("h","h",50,0.0,200.0);
    ///fill
    for (i=0;i<1600;i++) h->Fill(gRandom->Landau(38.77,14.06));
    std::cout << "Sampling point: " << h->GetRandom() << std::endl;
    TCanvas *c = new TCanvas("c","c",600,600);
    h->Draw();

}
