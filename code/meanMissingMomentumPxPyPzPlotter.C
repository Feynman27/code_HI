///////////////////////////
//This macro plots the mean and sigma 
//for the event missing
//px,py as a function of FCal_et
//@date: Dec 17, 2012
//////////////////////////

void meanMissingPxPyPlotter(){

TFile *fOut = new TFile("meanMissPxy_MB_.root","recreate") ;
TFile *fileName = new TFile("/tmp/tbalestr/"); 
TChain* data = new TChain("tree","tree");
int nFiles = data->Add(fileName);
std::cout <<"Filling the for "<< fileName << "...Number of files: " << nFiles << std::endl;

TCanvas cDummy = TCanvas("cDummy","cDummy",600,600);

#ifdef __CINT__
  gROOT->LoadMacro("AtlasUtils.C");
#endif

//stl container to store lower 
//track thresholds used in mpt calc
std::vector<int> missingPThr ;
missingPThr.push_back(500);
missingPThr.push_back(700);
missingPThr.push_back(2000);
missingPThr.push_back(3000);
missingPThr.push_back(4000);
//missingPThr.push_back(7000);
const unsigned int nMissPtBins = missingPThr.size();
std::cout << " Number of mpt trk thresholds: " << nMissPtBins << std::endl;

std::vector<double> centralityBins;
centralityBins.push_back(0.00);
centralityBins.push_back(0.05);
centralityBins.push_back(0.10);
centralityBins.push_back(0.15);
centralityBins.push_back(0.20);
centralityBins.push_back(0.40);
centralityBins.push_back(0.80);

const int nCentralityBins = centralityBins.size()-1;

TH1 *hPtSig[nMissPtBins][nCentralityBins];

TH1 *hpx[10],*hpy[10],*hpz[10] ;
char namePx[200], namePy[200], namePz[200];
int nbins = 7;

double pts[]  = {0.0,0.05,0.1,0.15,0.2,0.4,0.8} ;
std::cout << "initializing histograms..." << std::endl;
for(int ih=0;ih<10;ih++){
  for(int icent;icent<4;icent++){
    sprintf(namePx,"hpx%i_cent%i",ih,icent);
    sprintf(namePy,"hpy%i_cent%i",ih,icent);
    sprintf(namePz,"hpz%i_cent%i",ih,icent);
    hpx[ih] = new TH1F(namePx,namePx,nbins,pts);
    hpy[ih] = new TH1F(namePy,namePy,nbins,pts);
    hpz[ih] = new TH1F(namePz,namePz,nbins,pts);

    /*data->Draw("nu_px500>>"+namePx,"Centrality>=0.0&&Centrality<0.1&&nu_pt500<180.");
    meanPx[ih][icent] = hpx0_cent0->GetMean();
    errPx[ih][icent]  = hpx0_cent0->GetMeanError();
    std::cout << "Mean px with 500MeV track cut: " << meanPx0_cent0 << " +- " << errPx0_cent0 << std::endl;
    */
  }
 }
std::cout << "done." << std::endl;
///////////////////////
//mean pX vs centrality
//////////////////////
const int ntrkCuts = 10;
const int ncent = 4;
///0-10%
data->Draw("nu_px>>hpx0_cent0","Centrality>=0.0&&Centrality<0.1&&nu_pt500<300.");
float meanPx0_cent0 = hpx0_cent0->GetMean();
float errPx0_cent0  = hpx0_cent0->GetRMS();

data->Draw("caloMEx>>hpx9_cent0","Centrality>=0.0&&Centrality<0.1");
float meanEx0_cent0 = hex0_cent0->GetMean();
float errEx0_cent0  = hex0_cent0->GetMeanError();

///10-20%
data->Draw("nu_px>>hpx0_cent1","Centrality>=0.1&&Centrality<0.2&&nu_pt500<300.");
float meanPx0_cent1 = hpx0_cent1->GetMean();
float errEx0_cent1  = hex0_cent1->GetRMS();

data->Draw("caloMEx>>hpx9_cent1","Centrality>=0.1&&Centrality<0.2");
float meanEx0_cent1 = hex0_cent1->GetMean();
float errEx0_cent1  = hex0_cent1->GetMeanError();

//20-40%
data->Draw("nu_px>>hpx0_cent0","Centrality>=0.0&&Centrality<0.1&&nu_pt500<300.");
float meanPx0_cent2 = hpx0_cent2->GetMean();
float errEx0_cent2  = hex0_cent2->GetRMS();

data->Draw("caloMEx>>hpx9_cent2","Centrality>=0.0&&Centrality<0.1");
float meanEx0_cent2 = hex0_cent2->GetMean();
float errEx0_cent2  = hex0_cent2->GetMeanError();

//40-80%
data->Draw("nu_px500>>hpx0_cent3","Centrality>=0.1&&Centrality<0.2&&nu_pt500<180.");
float meanPx0_cent3 = hpx0_cent3->GetMean();
float errEx0_cent3  = hex0_cent3->GetMeanError();

data->Draw("nu_px600>>hpx1_cent3","Centrality>=0.1&&Centrality<0.2&&nu_pt600<180.");
float meanPx0_cent3 = hpx0_cent3->GetMean();
float errEx0_cent3  = hex0_cent3->GetMeanError();


float meanPx3000[4] = {meanPx9_cent0,meanPx9_cent1,meanPx9_cent2,meanPx9_cent3};
float errPx3000[4]  = {errPx9_cent0,errPx9_cent1,errPx9_cent2,errPx9_cent3} ;

float meanEx3000[4] = {meanEx9_cent0,meanEx9_cent1,meanEx9_cent2,meanEx9_cent3};
float errEx3000[4]  = {errEx9_cent0,errEx9_cent1,errEx9_cent2,errEx9_cent3} ;

float centrality[4] = {5.0,15.0,30.0,60.0};
float centralityBw[4] = {5.0,5.0,10.0,20.0};

TGraphErrors* grPx3000 = new TGraphErrors(4);grPx3000->SetNameTitle("grPx3000","grPx3000") ;
TGraphErrors* grEx3000 = new TGraphErrors(4);grEx3000->SetNameTitle("grEx3000","grEx3000") ;

grEx3000->SetMarkerColor(kMagenta+5);
grPx3000->SetMarkerColor(kOrange+5);


///0-10%
std::cout << "Plotting mean Pz as function of centrality." << std::endl;
std::cout << "Centrality 0-10%: " << std::endl;
data->Draw("nu_pz500>>hpz0_cent0","Centrality>=0.0&&Centrality<0.1&&nu_pt500<180.");
float meanPz0_cent0 = hpz0_cent0->GetMean();
float errPz0_cent0  = hpz0_cent0->GetMeanError();
std::cout << "Mean pz with 500MeV track cut: " << meanPz0_cent0 << " +- " << errPz0_cent0 << std::endl;

data->Draw("nu_pz600>>hpz1_cent0","Centrality>=0.0&&Centrality<0.1&&nu_pt600<180.");
float meanPz1_cent0 = hpz1_cent0->GetMean();
float errPz1_cent0  = hpz1_cent0->GetMeanError();
std::cout << "Mean pz with 600MeV track cut: " << meanPz1_cent0 << " +- " << errPz1_cent0 << std::endl;

data->Draw("nu_pz700>>hpz2_cent0","Centrality>=0.0&&Centrality<0.1&&nu_pt700<180.");
float meanPz2_cent0 = hpz2_cent0->GetMean();
float errPz2_cent0  = hpz2_cent0->GetMeanError();
std::cout << "Mean pz with 700MeV track cut: " << meanPz2_cent0 << " +- " << errPz2_cent0 << std::endl;

data->Draw("nu_pz800>>hpz3_cent0","Centrality>=0.0&&Centrality<0.1&&nu_pt800<180.");
float meanPz3_cent0 = hpz3_cent0->GetMean();
float errPz3_cent0  = hpz3_cent0->GetMeanError();
std::cout << "Mean pz with 800MeV track cut: " << meanPz3_cent0 << " +- " << errPz3_cent0 << std::endl;

data->Draw("nu_pz900>>hpz4_cent0","Centrality>=0.0&&Centrality<0.1&&nu_pt900<180.");
float meanPz4_cent0 = hpz4_cent0->GetMean();
float errPz4_cent0  = hpz4_cent0->GetMeanError();
std::cout << "Mean pz with 900MeV track cut: " << meanPz4_cent0 << " +- " << errPz4_cent0 << std::endl;

data->Draw("nu_pz1000>>hpz5_cent0","Centrality>=0.0&&Centrality<0.1&&nu_pt1000<180.");
float meanPz5_cent0 = hpz5_cent0->GetMean();
float errPz5_cent0  = hpz5_cent0->GetMeanError();
std::cout << "Mean pz with 1000MeV track cut: " << meanPz5_cent0 << " +- " << errPz5_cent0 << std::endl;

data->Draw("nu_pz1500>>hpz6_cent0","Centrality>=0.0&&Centrality<0.1&&nu_pt1500<180.");
float meanPz6_cent0 = hpz6_cent0->GetMean();
float errPz6_cent0  = hpz6_cent0->GetMeanError();
std::cout << "Mean pz with 1500MeV track cut: " << meanPz6_cent0 << " +- " << errPz6_cent0 << std::endl;

data->Draw("nu_pz2000>>hpz7_cent0","Centrality>=0.0&&Centrality<0.1&&nu_pt2000<180.");
float meanPz7_cent0 = hpz7_cent0->GetMean();
float errPz7_cent0  = hpz7_cent0->GetMeanError();
std::cout << "Mean pz with 2000MeV track cut: " << meanPz7_cent0 << " +- " << errPz7_cent0 << std::endl;

data->Draw("nu_pz2500>>hpz8_cent0","Centrality>=0.0&&Centrality<0.1&&nu_pt2500<180.");
float meanPz8_cent0 = hpz8_cent0->GetMean();
float errPz8_cent0  = hpz8_cent0->GetMeanError();
std::cout << "Mean pz with 2500MeV track cut: " << meanPz8_cent0 << " +- " << errPz8_cent0 << std::endl;

data->Draw("nu_pz3000>>hpz9_cent0","Centrality>=0.0&&Centrality<0.1&&nu_pt3000<180.");
float meanPz9_cent0 = hpz9_cent0->GetMean();
float errPz9_cent0  = hpz9_cent0->GetMeanError();
std::cout << "Mean pz with 3000MeV track cut: " << meanPz9_cent0 << " +- " << errPz9_cent0 << std::endl;

///10-20%
std::cout << "Centrality 10-20%: " << std::endl;
data->Draw("nu_pz500>>hpz0_cent1","Centrality>=0.1&&Centrality<0.2&&nu_pt500<180.");
float meanPz0_cent1 = hpz0_cent1->GetMean();
float errPz0_cent1  = hpz0_cent1->GetMeanError();
std::cout << "Mean pz with 500MeV track cut: " << meanPz0_cent1 << " +- " << errPz0_cent1 << std::endl;

data->Draw("nu_pz600>>hpz1_cent1","Centrality>=0.1&&Centrality<0.2&&nu_pt600<180.");
float meanPz1_cent1 = hpz1_cent1->GetMean();
float errPz1_cent1  = hpz1_cent1->GetMeanError();
std::cout << "Mean pz with 600MeV track cut: " << meanPz1_cent1 << " +- " << errPz1_cent1 << std::endl;

data->Draw("nu_pz700>>hpz2_cent1","Centrality>=0.1&&Centrality<0.2&&nu_pt700<180.");
float meanPz2_cent1 = hpz2_cent1->GetMean();
float errPz2_cent1  = hpz2_cent1->GetMeanError();
std::cout << "Mean pz with 700MeV track cut: " << meanPz2_cent1 << " +- " << errPz2_cent1 << std::endl;

data->Draw("nu_pz800>>hpz3_cent1","Centrality>=0.1&&Centrality<0.2&&nu_pt800<180.");
float meanPz3_cent1 = hpz3_cent1->GetMean();
float errPz3_cent1  = hpz3_cent1->GetMeanError();
std::cout << "Mean pz with 800MeV track cut: " << meanPz3_cent1 << " +- " << errPz3_cent1 << std::endl;

data->Draw("nu_pz900>>hpz4_cent1","Centrality>=0.1&&Centrality<0.2&&nu_pt900<180.");
float meanPz4_cent1 = hpz4_cent1->GetMean();
float errPz4_cent1  = hpz4_cent1->GetMeanError();
std::cout << "Mean pz with 900MeV track cut: " << meanPz4_cent1 << " +- " << errPz4_cent1 << std::endl;

data->Draw("nu_pz1000>>hpz5_cent1","Centrality>=0.1&&Centrality<0.2&&nu_pt1000<180.");
float meanPz5_cent1 = hpz5_cent1->GetMean();
float errPz5_cent1  = hpz5_cent1->GetMeanError();
std::cout << "Mean pz with 1000MeV track cut: " << meanPz5_cent1 << " +- " << errPz5_cent1 << std::endl;

data->Draw("nu_pz1500>>hpz6_cent1","Centrality>=0.1&&Centrality<0.2&&nu_pt1500<180.");
float meanPz6_cent1 = hpz6_cent1->GetMean();
float errPz6_cent1  = hpz6_cent1->GetMeanError();
std::cout << "Mean pz with 1500MeV track cut: " << meanPz6_cent1 << " +- " << errPz6_cent1 << std::endl;

data->Draw("nu_pz2000>>hpz7_cent1","Centrality>=0.1&&Centrality<0.2&&nu_pt2000<180.");
float meanPz7_cent1 = hpz7_cent1->GetMean();
float errPz7_cent1  = hpz7_cent1->GetMeanError();
std::cout << "Mean pz with 2000MeV track cut: " << meanPz7_cent1 << " +- " << errPz7_cent1 << std::endl;

data->Draw("nu_pz2500>>hpz8_cent1","Centrality>=0.1&&Centrality<0.2&&nu_pt2500<180.");
float meanPz8_cent1 = hpz8_cent1->GetMean();
float errPz8_cent1  = hpz8_cent1->GetMeanError();
std::cout << "Mean pz with 2500MeV track cut: " << meanPz8_cent1 << " +- " << errPz8_cent1 << std::endl;

data->Draw("nu_pz3000>>hpz9_cent1","Centrality>=0.1&&Centrality<0.2&&nu_pt3000<180.");
float meanPz9_cent1 = hpz9_cent1->GetMean();
float errPz9_cent1  = hpz9_cent1->GetMeanError();
std::cout << "Mean pz with 3000MeV track cut: " << meanPz9_cent1 << " +- " << errPz9_cent1 << std::endl;

///20-40%
std::cout << "Centrality 20-40%: " << std::endl;
data->Draw("nu_pz500>>hpz0_cent2","Centrality>=0.2&&Centrality<0.4&&nu_pt500<180.");
float meanPz0_cent2 = hpz0_cent2->GetMean();
float errPz0_cent2  = hpz0_cent2->GetMeanError();
std::cout << "Mean pz with 500MeV track cut: " << meanPz0_cent2 << " +- " << errPz0_cent2 << std::endl;

data->Draw("nu_pz600>>hpz1_cent2","Centrality>=0.2&&Centrality<0.4&&nu_pt600<180.");
float meanPz1_cent2 = hpz1_cent2->GetMean();
float errPz1_cent2  = hpz1_cent2->GetMeanError();
std::cout << "Mean pz with 600MeV track cut: " << meanPz1_cent2 << " +- " << errPz1_cent2 << std::endl;

data->Draw("nu_pz700>>hpz2_cent2","Centrality>=0.2&&Centrality<0.4&&nu_pt700<180.");
float meanPz2_cent2 = hpz2_cent2->GetMean();
float errPz2_cent2  = hpz2_cent2->GetMeanError();
std::cout << "Mean pz with 700MeV track cut: " << meanPz2_cent2 << " +- " << errPz2_cent2 << std::endl;

data->Draw("nu_pz800>>hpz3_cent2","Centrality>=0.2&&Centrality<0.4&&nu_pt800<180.");
float meanPz3_cent2 = hpz3_cent2->GetMean();
float errPz3_cent2  = hpz3_cent2->GetMeanError();
std::cout << "Mean pz with 800MeV track cut: " << meanPz3_cent2 << " +- " << errPz3_cent2 << std::endl;

data->Draw("nu_pz900>>hpz4_cent2","Centrality>=0.2&&Centrality<0.4&&nu_pt900<180.");
float meanPz4_cent2 = hpz4_cent2->GetMean();
float errPz4_cent2  = hpz4_cent2->GetMeanError();
std::cout << "Mean pz with 900MeV track cut: " << meanPz4_cent2 << " +- " << errPz4_cent2 << std::endl;

data->Draw("nu_pz1000>>hpz5_cent2","Centrality>=0.2&&Centrality<0.4&&nu_pt1000<180.");
float meanPz5_cent2 = hpz5_cent2->GetMean();
float errPz5_cent2  = hpz5_cent2->GetMeanError();
std::cout << "Mean pz with 1000MeV track cut: " << meanPz5_cent2 << " +- " << errPz5_cent2 << std::endl;

data->Draw("nu_pz1500>>hpz6_cent2","Centrality>=0.2&&Centrality<0.4&&nu_pt1500<180.");
float meanPz6_cent2 = hpz6_cent2->GetMean();
float errPz6_cent2  = hpz6_cent2->GetMeanError();
std::cout << "Mean pz with 1500MeV track cut: " << meanPz6_cent2 << " +- " << errPz6_cent2 << std::endl;

data->Draw("nu_pz2000>>hpz7_cent2","Centrality>=0.2&&Centrality<0.4&&nu_pt2000<180.");
float meanPz7_cent2 = hpz7_cent2->GetMean();
float errPz7_cent2  = hpz7_cent2->GetMeanError();
std::cout << "Mean pz with 2000MeV track cut: " << meanPz7_cent2 << " +- " << errPz7_cent2 << std::endl;

data->Draw("nu_pz2500>>hpz8_cent2","Centrality>=0.2&&Centrality<0.4&&nu_pt2500<180.");
float meanPz8_cent2 = hpz8_cent2->GetMean();
float errPz8_cent2  = hpz8_cent2->GetMeanError();
std::cout << "Mean pz with 2500MeV track cut: " << meanPz8_cent2 << " +- " << errPz8_cent2 << std::endl;

data->Draw("nu_pz3000>>hpz9_cent2","Centrality>=0.2&&Centrality<0.4&&nu_pt3000<180.");
float meanPz9_cent2 = hpz9_cent2->GetMean();
float errPz9_cent2  = hpz9_cent2->GetMeanError();
std::cout << "Mean pz with 3000MeV track cut: " << meanPz9_cent2 << " +- " << errPz9_cent2 << std::endl;

///40-80%
std::cout << "Centrality 40-80%: " << std::endl;
data->Draw("nu_pz500>>hpz0_cent3","Centrality>=0.4&&Centrality<0.8&&nu_pt500<180.");
float meanPz0_cent3 = hpz0_cent3->GetMean();
float errPz0_cent3  = hpz0_cent3->GetMeanError();
std::cout << "Mean pz with 500MeV track cut: " << meanPz0_cent3 << " +- " << errPz0_cent3 << std::endl;

data->Draw("nu_pz600>>hpz1_cent3","Centrality>=0.4&&Centrality<0.8&&nu_pt600<180.");
float meanPz1_cent3 = hpz1_cent3->GetMean();
float errPz1_cent3  = hpz1_cent3->GetMeanError();
std::cout << "Mean pz with 600MeV track cut: " << meanPz1_cent3 << " +- " << errPz1_cent3 << std::endl;

data->Draw("nu_pz700>>hpz2_cent3","Centrality>=0.4&&Centrality<0.8&&nu_pt700<180.");
float meanPz2_cent3 = hpz2_cent3->GetMean();
float errPz2_cent3  = hpz2_cent3->GetMeanError();
std::cout << "Mean pz with 700MeV track cut: " << meanPz2_cent3 << " +- " << errPz2_cent3 << std::endl;

data->Draw("nu_pz800>>hpz3_cent3","Centrality>=0.4&&Centrality<0.8&&nu_pt800<180.");
float meanPz3_cent3 = hpz3_cent3->GetMean();
float errPz3_cent3  = hpz3_cent3->GetMeanError();
std::cout << "Mean pz with 800MeV track cut: " << meanPz3_cent3 << " +- " << errPz3_cent3 << std::endl;

data->Draw("nu_pz900>>hpz4_cent3","Centrality>=0.4&&Centrality<0.8&&nu_pt900<180.");
float meanPz4_cent3 = hpz4_cent3->GetMean();
float errPz4_cent3  = hpz4_cent3->GetMeanError();
std::cout << "Mean pz with 900MeV track cut: " << meanPz4_cent3 << " +- " << errPz4_cent3 << std::endl;

data->Draw("nu_pz1000>>hpz5_cent3","Centrality>=0.4&&Centrality<0.8&&nu_pt1000<180.");
float meanPz5_cent3 = hpz5_cent3->GetMean();
float errPz5_cent3  = hpz5_cent3->GetMeanError();
std::cout << "Mean pz with 1000MeV track cut: " << meanPz5_cent3 << " +- " << errPz5_cent3 << std::endl;

data->Draw("nu_pz1500>>hpz6_cent3","Centrality>=0.4&&Centrality<0.8&&nu_pt1500<180.");
float meanPz6_cent3 = hpz6_cent3->GetMean();
float errPz6_cent3  = hpz6_cent3->GetMeanError();
std::cout << "Mean pz with 1500MeV track cut: " << meanPz6_cent3 << " +- " << errPz6_cent3 << std::endl;

data->Draw("nu_pz2000>>hpz7_cent3","Centrality>=0.4&&Centrality<0.8&&nu_pt2000<180.");
float meanPz7_cent3 = hpz7_cent3->GetMean();
float errPz7_cent3  = hpz7_cent3->GetMeanError();
std::cout << "Mean pz with 2000MeV track cut: " << meanPz7_cent3 << " +- " << errPz7_cent3 << std::endl;

data->Draw("nu_pz2500>>hpz8_cent3","Centrality>=0.4&&Centrality<0.8&&nu_pt2500<180.");
float meanPz8_cent3 = hpz8_cent3->GetMean();
float errPz8_cent3  = hpz8_cent3->GetMeanError();
std::cout << "Mean pz with 2500MeV track cut: " << meanPz8_cent3 << " +- " << errPz8_cent3 << std::endl;

data->Draw("nu_pz3000>>hpz9_cent3","Centrality>=0.4&&Centrality<0.8&&nu_pt3000<180.");
float meanPz9_cent3 = hpz9_cent3->GetMean();
float errPz9_cent3  = hpz9_cent3->GetMeanError();
std::cout << "Mean pz with 3000MeV track cut: " << meanPz9_cent3 << " +- " << errPz9_cent3 << std::endl;


float meanPz500[4] = {meanPz0_cent0,meanPz0_cent1,meanPz0_cent2,meanPz0_cent3};
float errPz500[4]  = {errPz0_cent0,errPz0_cent1,errPz0_cent2,errPz0_cent3} ;

float meanPz600[4] = {meanPz1_cent0,meanPz1_cent1,meanPz1_cent2,meanPz1_cent3};
float errPz600[4]  = {errPz1_cent0,errPz1_cent1,errPz1_cent2,errPz1_cent3} ;

float meanPz700[4] = {meanPz2_cent0,meanPz2_cent1,meanPz2_cent2,meanPz2_cent3};
float errPz700[4]  = {errPz2_cent0,errPz2_cent1,errPz2_cent2,errPz2_cent3} ;

float meanPz800[4] = {meanPz3_cent0,meanPz3_cent1,meanPz3_cent2,meanPz3_cent3};
float errPz800[4]  = {errPz3_cent0,errPz3_cent1,errPz3_cent2,errPz3_cent3} ;

float meanPz900[4] = {meanPz4_cent0,meanPz4_cent1,meanPz4_cent2,meanPz4_cent3};
float errPz900[4]  = {errPz4_cent0,errPz4_cent1,errPz4_cent2,errPz4_cent3} ;

float meanPz1000[4] = {meanPz5_cent0,meanPz5_cent1,meanPz5_cent2,meanPz5_cent3};
float errPz1000[4]  = {errPz5_cent0,errPz5_cent1,errPz5_cent2,errPz5_cent3} ;

float meanPz1500[4] = {meanPz6_cent0,meanPz6_cent1,meanPz6_cent2,meanPz6_cent3};
float errPz1500[4]  = {errPz6_cent0,errPz6_cent1,errPz6_cent2,errPz6_cent3} ;

float meanPz2000[4] = {meanPz7_cent0,meanPz7_cent1,meanPz7_cent2,meanPz7_cent3};
float errPz2000[4]  = {errPz7_cent0,errPz7_cent1,errPz7_cent2,errPz7_cent3} ;

float meanPz2500[4] = {meanPz8_cent0,meanPz8_cent1,meanPz8_cent2,meanPz8_cent3};
float errPz2500[4]  = {errPz8_cent0,errPz8_cent1,errPz8_cent2,errPz8_cent3} ;

float meanPz3000[4] = {meanPz9_cent0,meanPz9_cent1,meanPz9_cent2,meanPz9_cent3};
float errPz3000[4]  = {errPz9_cent0,errPz9_cent1,errPz9_cent2,errPz9_cent3} ;

float centrality[4] = {5.0,15.0,30.0,60.0};
float centralityBw[4] = {5.0,5.0,10.0,20.0};
TGraphErrors* grPz500 = new TGraphErrors(4); grPz500->SetNameTitle("grPz500","grPz500") ;
TGraphErrors* grPz600 = new TGraphErrors(4);grPz600->SetNameTitle("grPz600","grPz600") ;
TGraphErrors* grPz700 = new TGraphErrors(4);grPz700->SetNameTitle("grPz700","grPz700") ;
TGraphErrors* grPz800 = new TGraphErrors(4);grPz800->SetNameTitle("grPz800","grPz800") ;
TGraphErrors* grPz900 = new TGraphErrors(4);grPz900->SetNameTitle("grPz900","grPz900") ;
TGraphErrors* grPz1000 = new TGraphErrors(4);grPz1000->SetNameTitle("grPz1000","grPz1000") ;
TGraphErrors* grPz1500 = new TGraphErrors(4);grPz1500->SetNameTitle("grPz1500","grPz1500") ;
TGraphErrors* grPz2000 = new TGraphErrors(4);grPz2000->SetNameTitle("grPz2000","grPz2000") ;
TGraphErrors* grPz2500 = new TGraphErrors(4);grPz2500->SetNameTitle("grPz2500","grPz2500") ;
TGraphErrors* grPz3000 = new TGraphErrors(4);grPz3000->SetNameTitle("grPz3000","grPz3000") ;

//grPz500->SetMarkerStyle(kOpenCircle);
grPz500->SetMarkerColor(kRed);
grPz600->SetMarkerColor(kBlue);
grPz700->SetMarkerColor(kGreen);
grPz800->SetMarkerColor(kMagenta);
grPz900->SetMarkerColor(kOrange);
grPz1000->SetMarkerColor(kRed+5);
grPz1500->SetMarkerColor(kBlue+5);
grPz2000->SetMarkerColor(kGreen+5);
grPz2500->SetMarkerColor(kMagenta+5);
grPz3000->SetMarkerColor(kOrange+5);



for(int ipx=0; ipx<4;ipx++){
  grPx500->SetPoint(ipx,centrality[ipx],meanPx500[ipx]);
  grPx500->SetPointError(ipx,centralityBw[ipx],errPx500[ipx]);

  grPx600->SetPoint(ipx,centrality[ipx],meanPx600[ipx]);
  grPx600->SetPointError(ipx,centralityBw[ipx],errPx600[ipx]);

  grPx700->SetPoint(ipx,centrality[ipx],meanPx700[ipx]);
  grPx700->SetPointError(ipx,centralityBw[ipx],errPx700[ipx]);

  grPx800->SetPoint(ipx,centrality[ipx],meanPx800[ipx]);
  grPx800->SetPointError(ipx,centralityBw[ipx],errPx800[ipx]);

  grPx900->SetPoint(ipx,centrality[ipx],meanPx900[ipx]);
  grPx900->SetPointError(ipx,centralityBw[ipx],errPx900[ipx]);

  grPx1000->SetPoint(ipx,centrality[ipx],meanPx1000[ipx]);
  grPx1000->SetPointError(ipx,centralityBw[ipx],errPx1000[ipx]);

  grPx1500->SetPoint(ipx,centrality[ipx],meanPx1500[ipx]);
  grPx1500->SetPointError(ipx,centralityBw[ipx],errPx1500[ipx]);

  grPx2000->SetPoint(ipx,centrality[ipx],meanPx2000[ipx]);
  grPx2000->SetPointError(ipx,centralityBw[ipx],errPx2000[ipx]);

  grPx2500->SetPoint(ipx,centrality[ipx],meanPx2500[ipx]);
  grPx2500->SetPointError(ipx,centralityBw[ipx],errPx2500[ipx]);

  grPx3000->SetPoint(ipx,centrality[ipx],meanPx3000[ipx]);
  grPx3000->SetPointError(ipx,centralityBw[ipx],errPx3000[ipx]);
}
for(int ipy=0; ipy<4;ipy++){
  grPy500->SetPoint(ipy,centrality[ipy],meanPy500[ipy]);
  grPy500->SetPointError(ipy,centralityBw[ipy],errPy500[ipy]);

  grPy600->SetPoint(ipy,centrality[ipy],meanPy600[ipy]);
  grPy600->SetPointError(ipy,centralityBw[ipy],errPy600[ipy]);

  grPy700->SetPoint(ipy,centrality[ipy],meanPy700[ipy]);
  grPy700->SetPointError(ipy,centralityBw[ipy],errPy700[ipy]);

  grPy800->SetPoint(ipy,centrality[ipy],meanPy800[ipy]);
  grPy800->SetPointError(ipy,centralityBw[ipy],errPy800[ipy]);

  grPy900->SetPoint(ipy,centrality[ipy],meanPy900[ipy]);
  grPy900->SetPointError(ipy,centralityBw[ipy],errPy900[ipy]);

  grPy1000->SetPoint(ipy,centrality[ipy],meanPy1000[ipy]);
  grPy1000->SetPointError(ipy,centralityBw[ipy],errPy1000[ipy]);

  grPy1500->SetPoint(ipy,centrality[ipy],meanPy1500[ipy]);
  grPy1500->SetPointError(ipy,centralityBw[ipy],errPy1500[ipy]);

  grPy2000->SetPoint(ipy,centrality[ipy],meanPy2000[ipy]);
  grPy2000->SetPointError(ipy,centralityBw[ipy],errPy2000[ipy]);

  grPy2500->SetPoint(ipy,centrality[ipy],meanPy2500[ipy]);
  grPy2500->SetPointError(ipy,centralityBw[ipy],errPy2500[ipy]);

  grPy3000->SetPoint(ipy,centrality[ipy],meanPy3000[ipy]);
  grPy3000->SetPointError(ipy,centralityBw[ipy],errPy3000[ipy]);
}

for(int ipz=0; ipz<4;ipz++){
  grPz500->SetPoint(ipz,centrality[ipz],meanPz500[ipz]);
  grPz500->SetPointError(ipz,centralityBw[ipz],errPz500[ipz]);

  grPz600->SetPoint(ipz,centrality[ipz],meanPz600[ipz]);
  grPz600->SetPointError(ipz,centralityBw[ipz],errPz600[ipz]);

  grPz700->SetPoint(ipz,centrality[ipz],meanPz700[ipz]);
  grPz700->SetPointError(ipz,centralityBw[ipz],errPz700[ipz]);

  grPz800->SetPoint(ipz,centrality[ipz],meanPz800[ipz]);
  grPz800->SetPointError(ipz,centralityBw[ipz],errPz800[ipz]);

  grPz900->SetPoint(ipz,centrality[ipz],meanPz900[ipz]);
  grPz900->SetPointError(ipz,centralityBw[ipz],errPz900[ipz]);

  grPz1000->SetPoint(ipz,centrality[ipz],meanPz1000[ipz]);
  grPz1000->SetPointError(ipz,centralityBw[ipz],errPz1000[ipz]);

  grPz1500->SetPoint(ipz,centrality[ipz],meanPz1500[ipz]);
  grPz1500->SetPointError(ipz,centralityBw[ipz],errPz1500[ipz]);

  grPz2000->SetPoint(ipz,centrality[ipz],meanPz2000[ipz]);
  grPz2000->SetPointError(ipz,centralityBw[ipz],errPz2000[ipz]);

  grPz2500->SetPoint(ipz,centrality[ipz],meanPz2500[ipz]);
  grPz2500->SetPointError(ipz,centralityBw[ipz],errPz2500[ipz]);

  grPz3000->SetPoint(ipz,centrality[ipz],meanPz3000[ipz]);
  grPz3000->SetPointError(ipz,centralityBw[ipz],errPz3000[ipz]);
}

std::cout << "Saving graphs..." << std::endl;
fOut->cd();
grPx500->Write();
grPx600->Write();
grPx700->Write();
grPx800->Write();
grPx900->Write();
grPx1000->Write();
grPx1500->Write();
grPx2000->Write();
grPx2500->Write();
grPx3000->Write();
grPy500->Write();
grPy600->Write();
grPy700->Write();
grPy800->Write();
grPy900->Write();
grPy1000->Write();
grPy1500->Write();
grPy2000->Write();
grPy2500->Write();
grPy3000->Write();
grPz500->Write();
grPz600->Write();
grPz700->Write();
grPz800->Write();
grPz900->Write();
grPz1000->Write();
grPz1500->Write();
grPz2000->Write();
grPz2500->Write();
grPz3000->Write();
fOut->Write();
fOut->Close();
std::cout << "Done." << std::endl;

TCanvas c0 = TCanvas("c0","c0",600,600);
grPx500->GetYaxis()->SetRangeUser(0.0,11.0) ;
grPx500->GetYaxis()->SetTitle("#LT#slash{p_{x}}#RT") ;
grPx500->GetXaxis()->SetTitle("centrality (%)") ;
grPx500->Draw("ape")
grPx600->Draw("pesame")
grPx700->Draw("pesame")
grPx800->Draw("pesame")
grPx900->Draw("pesame")
grPx1000->Draw("pesame")
grPx1500->Draw("pesame")
grPx2000->Draw("pesame")
grPx2500->Draw("pesame")
grPx3000->Draw("pesame")

TLegend* leg0 = new TLegend(0.55, 0.5, 0.9, 0.9);
leg0->SetTextFont(gStyle->GetTextFont());
leg0->SetTextSize(gStyle->GetTextSize());
leg0->SetBorderSize(0);
leg0->SetFillColor(0);
  leg0->AddEntry(grPx500, "500MeV", "p");
  leg0->AddEntry(grPx600, "600MeV", "p");
  leg0->AddEntry(grPx700, "700MeV", "p");
  leg0->AddEntry(grPx800, "800MeV", "p");
  leg0->AddEntry(grPx900, "900MeV", "p");
  leg0->AddEntry(grPx1000, "1000MeV", "p");
  leg0->AddEntry(grPx1500, "1500MeV", "p");
  leg0->AddEntry(grPx2000, "2000MeV", "p");
  leg0->AddEntry(grPx2500, "2500MeV", "p");
  leg0->AddEntry(grPx3000, "3000MeV", "p");
  leg0->Draw();
c0->Print("meanMissingPx_trkCuts_28June.pdf") ;

TCanvas c1 = TCanvas("c1","c1",600,600);
grPy500->GetYaxis()->SetRangeUser(0.0,11.0) ;
grPy500->GetYaxis()->SetTitle("#LT#slash{p_{x}}#RT") ;
grPy500->GetXaxis()->SetTitle("centrality (%)") ;
grPy500->Draw("ape")
grPy600->Draw("pesame")
grPy700->Draw("pesame")
grPy800->Draw("pesame")
grPy900->Draw("pesame")
grPy1000->Draw("pesame")
grPy1500->Draw("pesame")
grPy2000->Draw("pesame")
grPy2500->Draw("pesame")
grPy3000->Draw("pesame")

TLegend* leg0 = new TLegend(0.55, 0.5, 0.9, 0.9);
leg0->SetTextFont(gStyle->GetTextFont());
leg0->SetTextSize(gStyle->GetTextSize());
leg0->SetBorderSize(0);
leg0->SetFillColor(0);
  leg0->AddEntry(grPy500, "500MeV", "p");
  leg0->AddEntry(grPy600, "600MeV", "p");
  leg0->AddEntry(grPy700, "700MeV", "p");
  leg0->AddEntry(grPy800, "800MeV", "p");
  leg0->AddEntry(grPy900, "900MeV", "p");
  leg0->AddEntry(grPy1000, "1000MeV", "p");
  leg0->AddEntry(grPy1500, "1500MeV", "p");
  leg0->AddEntry(grPy2000, "2000MeV", "p");
  leg0->AddEntry(grPy2500, "2500MeV", "p");
  leg0->AddEntry(grPy3000, "3000MeV", "p");
  leg0->Draw();
c1->Print("meanMissingPy_trkCuts_28June.pdf") ;

TCanvas c2 = TCanvas("c2","c2",600,600);
grPz500->GetYaxis()->SetRangeUser(0.0,11.0) ;
grPz500->GetYaxis()->SetTitle("#LT#slash{p_{x}}#RT") ;
grPz500->GetXaxis()->SetTitle("centrality (%)") ;
grPz500->Draw("ape")
grPz600->Draw("pesame")
grPz700->Draw("pesame")
grPz800->Draw("pesame")
grPz900->Draw("pesame")
grPz1000->Draw("pesame")
grPz1500->Draw("pesame")
grPz2000->Draw("pesame")
grPz2500->Draw("pesame")
grPz3000->Draw("pesame")

TLegend* leg0 = new TLegend(0.55, 0.5, 0.9, 0.9);
leg0->SetTextFont(gStyle->GetTextFont());
leg0->SetTextSize(gStyle->GetTextSize());
leg0->SetBorderSize(0);
leg0->SetFillColor(0);
  leg0->AddEntry(grPz500, "500MeV", "p");
  leg0->AddEntry(grPz600, "600MeV", "p");
  leg0->AddEntry(grPz700, "700MeV", "p");
  leg0->AddEntry(grPz800, "800MeV", "p");
  leg0->AddEntry(grPz900, "900MeV", "p");
  leg0->AddEntry(grPz1000, "1000MeV", "p");
  leg0->AddEntry(grPz1500, "1500MeV", "p");
  leg0->AddEntry(grPz2000, "2000MeV", "p");
  leg0->AddEntry(grPz2500, "2500MeV", "p");
  leg0->AddEntry(grPz3000, "3000MeV", "p");
  leg0->Draw();
c2->Print("meanMissingPz_trkCuts_28June.pdf") ;
//myText(0.4,0.75,1,"#int L #approx 5 #mub^{-1}");

}
