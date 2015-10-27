#include "TVectorF.h"
#include "TGraphErrors.h"

int __npoints = 6; //number of centrality bins
TVectorF xp = TVectorF(__npoints);
TVectorF xpe = TVectorF(__npoints);

const unsigned int _nCentrality = 6;
const unsigned int _nEta = 9;
const unsigned int _nBins = _nCentrality*_nEta;
double isoSystMuPlus[_nBins] ;
double isoSystMuMinus[_nBins];
std::vector <double> mptSystMuPlus (_nBins);
std::vector <double> mptSystMuMinus (_nBins);
std::vector <double> qcdRaaSystMuPlus (_nBins);
std::vector <double> qcdRaaSystMuMinus (_nBins);
std::vector <double> electroweakSystMuPlus (_nBins);
std::vector <double> electroweakSystMuMinus (_nBins);
std::vector <double> cwEtaSystMuPlus (_nBins);
std::vector <double> cwEtaSystMuMinus (_nBins);

//////////////////////////////////////
//Uncertainty in using eta-mirrored Cw
//////////////////////////////////////
void setCwMirroredEtaSystErrors(TString sFileIn = "systematics/cwMirroredEtaSystematics.05.15.2013.txt"){
    
	unsigned int totalBins =  _nBins;
	std::ifstream s( sFileIn , std::ifstream::in );

	if(!s.is_open()){
		std::cout << "While opening a file an error was encountered" << std::endl;
		exit(1);
	} else {
		std::cout << "File " << sFileIn << " is successfully opened" << std::endl;
	}

    if(s.is_open()){
        for (unsigned int i = 0; i < cwEtaSystMuPlus.size(); ++i) {
		    s >> cwEtaSystMuPlus[i] >> cwEtaSystMuMinus[i];
            std::cout << "Cw Eta systematic : mu+: " << cwEtaSystMuPlus[i] << " mu- : " << cwEtaSystMuMinus[i] << std::endl;  
        }
        if( (cwEtaSystMuPlus.size()!=totalBins) || (cwEtaSystMuMinus.size()!=totalBins) ){
            std::cout << "ERROR: Cw Eta container size does not match number of total bins used. Aborting." << std::endl;
            exit(0);
        }
        s.close();
	}
    else std::cout << "unable to open cwEta systematics file" << std::endl;
}

///return PERCENT systematic error for Cw Eta scaling
///for given charge,centrality, and eta class
float getCwMirroredEtaSystErr(int chargeIndex=-1, int binIndex=-1){
    if(chargeIndex==-1||binIndex==-1){
        std::cout << "ERROR READING Cw Eta SYSTEMATICS. " << std::endl;
        return -1.0;
    }
    else if(chargeIndex==102) {
        std::cout << "Return cwEta systematic for mu+ at index: " << binIndex << " = " << cwEtaSystMuPlus[binIndex] << std::endl;
        return cwEtaSystMuPlus[binIndex];
    }
    else if(chargeIndex==103){
        std::cout << "Return isolation systematic for mu- at index: " << binIndex << " = " << cwEtaSystMuMinus[binIndex] << std::endl;
        return cwEtaSystMuMinus[binIndex];
    }
    else {

        std::cout << "ERROR READING Cw Eta SYSTEMATICS. " << std::endl;
        return -1.0;
    }
}


/////////////////////////////////////////////////////
//Uncertainty in electroweak bkg(Z->mumu,W->tau->mu) 
/////////////////////////////////////////////////////
void setElectroweakSystErrors(TString sFileIn = "systematics/electroweakSystematics.07.28.2013.txt"){
    
	unsigned int totalBins =  _nBins;
	std::ifstream s( sFileIn , std::ifstream::in );

	if(!s.is_open()){
		std::cout << "While opening a file an error was encountered" << std::endl;
		exit(1);
	} else {
		std::cout << "File " << sFileIn << " is successfully opened" << std::endl;
	}

    if(s.is_open()){
        for (unsigned int i = 0; i < electroweakSystMuPlus.size(); ++i) {
		    s >> electroweakSystMuPlus[i] >> electroweakSystMuMinus[i];
            std::cout << "Electroweak systematic : mu+: " << electroweakSystMuPlus[i] << " mu- : " << electroweakSystMuMinus[i] << std::endl;  
        }
        if( (electroweakSystMuPlus.size()!=totalBins) || (electroweakSystMuMinus.size()!=totalBins) ){
            std::cout << "ERROR: QCD Raa container size does not match number of total bins used. Aborting." << std::endl;
            exit(0);
        }
        s.close();
	}
    else std::cout << "unable to open electroweak systematics file" << std::endl;
}
//////////////////////////////////////
//Uncertainty in QCD shape at high pT 
//////////////////////////////////////
void setQCDRaaSystErrors(TString sFileIn = "systematics/qcdRaaSystematics.05.13.2013.txt"){
    
	unsigned int totalBins =  _nBins;
	std::ifstream s( sFileIn , std::ifstream::in );

	if(!s.is_open()){
		std::cout << "While opening a file an error was encountered" << std::endl;
		exit(1);
	} else {
		std::cout << "File " << sFileIn << " is successfully opened" << std::endl;
	}

    if(s.is_open()){
        for (unsigned int i = 0; i < qcdRaaSystMuPlus.size(); ++i) {
		    s >> qcdRaaSystMuPlus[i] >> qcdRaaSystMuMinus[i];
            std::cout << "QCD Raa systematic : mu+: " << qcdRaaSystMuPlus[i] << " mu- : " << qcdRaaSystMuMinus[i] << std::endl;  
        }
        if( (qcdRaaSystMuPlus.size()!=totalBins) || (qcdRaaSystMuMinus.size()!=totalBins) ){
            std::cout << "ERROR: QCD Raa container size does not match number of total bins used. Aborting." << std::endl;
            exit(0);
        }
        s.close();
	}
    else std::cout << "unable to open qcdRaa systematics file" << std::endl;
}

//////////////////////////////////////////////////////
///return PERCENT systematic error for electroweak bkg 
///for given charge,centrality, and eta class
//////////////////////////////////////////////////////
float getElectroweakSystErr(int chargeIndex=-1, int binIndex=-1){
    if(chargeIndex==-1||binIndex==-1){
        std::cout << "ERROR READING ELECTROWEAK SYSTEMATICS. " << std::endl;
        return -1.0;
    }
    else if(chargeIndex==102) {
        std::cout << "Return electroweak systematic for mu+ at index: " << binIndex << " = " << electroweakSystMuPlus[binIndex] << std::endl;
        return electroweakSystMuPlus[binIndex];
    }
    else if(chargeIndex==103){
        std::cout << "Return electroweak systematic for mu- at index: " << binIndex << " = " << electroweakSystMuMinus[binIndex] << std::endl;
        return electroweakSystMuMinus[binIndex];
    }
    else {

        std::cout << "ERROR READING ELECTROWEAK SYSTEMATICS. " << std::endl;
        return -1.0;
    }
}
//////////////////////////////////////////////////////
///return PERCENT systematic error for QCD Raa scaling
///for given charge,centrality, and eta class
/////////////////////////////////////////////////////////
float getQCDRaaSystErr(int chargeIndex=-1, int binIndex=-1){
    if(chargeIndex==-1||binIndex==-1){
        std::cout << "ERROR READING QCD Raa SYSTEMATICS. " << std::endl;
        return -1.0;
    }
    else if(chargeIndex==102) {
        std::cout << "Return qcdRaa systematic for mu+ at index: " << binIndex << " = " << qcdRaaSystMuPlus[binIndex] << std::endl;
        return qcdRaaSystMuPlus[binIndex];
    }
    else if(chargeIndex==103){
        std::cout << "Return qcdRaa systematic for mu- at index: " << binIndex << " = " << qcdRaaSystMuMinus[binIndex] << std::endl;
        return qcdRaaSystMuMinus[binIndex];
    }
    else {

        std::cout << "ERROR READING QCD Raa SYSTEMATICS. " << std::endl;
        return -1.0;
    }
}

//////////////////////////////////////
//Uncertainty in MPT due to resolution 
//////////////////////////////////////
void setMptSystErrors(TString sFileIn = "systematics/mptCutSystematics.04.24.2013.txt"){
    
	unsigned int totalBins =  _nBins;
	std::ifstream s( sFileIn , std::ifstream::in );

	if(!s.is_open()){
		std::cout << "While opening a file an error was encountered" << std::endl;
		exit(1);
	} else {
		std::cout << "File " << sFileIn << " is successfully opened" << std::endl;
	}

    if(s.is_open()){
        for (unsigned int i = 0; i < mptSystMuPlus.size(); ++i) {
		    s >> mptSystMuPlus[i] >> mptSystMuMinus[i];
            std::cout << "MPT systematic : mu+: " << mptSystMuPlus[i] << " mu- : " << mptSystMuMinus[i] << std::endl;  
        }
        if( (mptSystMuPlus.size()!=totalBins) || (mptSystMuMinus.size()!=totalBins) ){
            std::cout << "ERROR: mpt container size does not match number of total bins used. Aborting." << std::endl;
            exit(0);
        }
        s.close();
	}
    else std::cout << "unable to open mpt systematics file" << std::endl;
}

///return PERCENT systematic error in mpt cuts
///for given charge,centrality, and eta class
float getMptSystErr(int chargeIndex=-1, int binIndex=-1){
    if(chargeIndex==-1||binIndex==-1){
        std::cout << "ERROR READING MPT SYSTEMATICS. " << std::endl;
        return -1.0;
    }
    else if(chargeIndex==102) {
        std::cout << "Return mpt systematic for mu+ at index: " << binIndex << " = " << mptSystMuPlus[binIndex] << std::endl;
        return mptSystMuPlus[binIndex];
    }
    else if(chargeIndex==103){
        std::cout << "Return isolation systematic for mu- at index: " << binIndex << " = " << mptSystMuMinus[binIndex] << std::endl;
        return mptSystMuMinus[binIndex];
    }
    else {

        std::cout << "ERROR READING MPT SYSTEMATICS. " << std::endl;
        return -1.0;
    }
}

//////////////////////////////////////
//Uncertainty in isolation efficiency 
//////////////////////////////////////
void setIsoSystErrors(TString sFileIn = "systematics/IsolationCutSystematics.04.03.2013.txt") {
    
	std::cout << "Opening file " << sFileIn << std::endl;
	std::ifstream s( sFileIn , std::ifstream::in );
	if(!s.is_open()){
		std::cout << "While opening a file an error was encountered" << std::endl;
		exit(1);
        } else {
                std::cout << "File " << sFileIn << " is successfully opened" << std::endl;
        }

    unsigned int i = 0;
	
    std::cout << "Filling iso systematic arrays for mu+/-.... " << std::endl;
	while(!s.eof()) {
        //read in systematic errors; ordered in txt
        //by looping over all centrality classes per eta bin
        s >> isoSystMuPlus[i] >> isoSystMuMinus[i];
        std::cout << "Isolation systematic : mu+: " << isoSystMuPlus[i] << " mu- : " << isoSystMuMinus[i] << std::endl;  
		++i;
	}	

        std::cout << "Done." << std::endl;

}

///return PERCENT systematic error in isolation cuts
///for given charge,centrlity, and eta class
float getIsoSystErr(int chargeIndex=-1, int binIndex=-1){
    if(chargeIndex==-1||binIndex==-1){
        std::cout << "ERROR READING ISO SYSTEMATICS. " << std::endl;
        return -1.0;
    }
    else if(chargeIndex==102) {
        std::cout << "Return isolation systematic for mu+ at index: " << binIndex << " = " << isoSystMuPlus[binIndex] << std::endl;
        return isoSystMuPlus[binIndex];
    }
    else if(chargeIndex==103){
        std::cout << "Return isolation systematic for mu- at index: " << binIndex << " = " << isoSystMuMinus[binIndex] << std::endl;
        return isoSystMuMinus[binIndex];
    }
    else {

        std::cout << "ERROR READING ISO SYSTEMATICS. " << std::endl;
        return -1.0;
    }
}

void setpoints() {
  //xp[0] = 95;
  //xpe[0] = 5.;
  //xp[1] = 85;
  //xpe[1] = 5.;
  //xp[2] = 70;
  //xpe[2] = 10.;
  //xp[3] = 40;
  //xpe[3] = 20.;
 /* xp[0] = 70;
  xpe[0] = 10.;
  xp[1] = 50;
  xpe[1] = 10.;
  xp[2] = 30;
  xpe[2] = 10.;
  xp[3] = 17.5;
  xpe[3] = 2.5;
  xp[4] = 12.5;
  xpe[4] = 2.5;
  xp[5] = 7.5;
  xpe[5] = 2.5;
  xp[6] = 2.5;
  xpe[6] = 2.5;
  */
  xp[0] = 2.5;
  xpe[0] = 2.5;
  xp[1] = 7.5;
  xpe[1] = 2.5;
  xp[2] = 12.5;
  xpe[2] = 2.5;
  xp[3] = 17.5;
  xpe[3] = 2.5;
  xp[4] = 30.0;
  xpe[4] = 10.0;
  xp[5] = 60.0;
  xpe[5] = 20.0;
}


TGraphErrors* yield() {
  TVectorF yp(__npoints);
  TVectorF ype(__npoints);
  yp[0] = 188;
  yp[1] = 160;
  yp[2] = 180;
  yp[3] = 93;
  ype[0] = 22;
  ype[1] = 19;
  ype[2] = 17;
  ype[3] = 12;
  TGraphErrors* func = new TGraphErrors(xp,yp,xpe,ype);
  return func;
}


TGraphAsymmErrors* efficiency() {
  setpoints();
  TVectorF yp(__npoints);
  TVectorF ypeL(__npoints);
  TVectorF ypeH(__npoints);
  yp[0] = 0.380615;
  yp[1] = 0.383965;
  yp[2] = 0.384777;
  yp[3] = 0.396718;
  ypeL[0] = 0.0029118;
  ypeL[1] = 0.00355478;
  ypeL[2] = 0.00259136;
  ypeL[3] = 0.00141997;
  ypeH[0] = 0.00291754;
  ypeH[1] = 0.00356306;
  ypeH[2] = 0.00259573;
  ypeH[3] = 0.00142113;
  
  TGraphAsymmErrors* func = new TGraphAsymmErrors(xp,yp,xpe,xpe,ypeL,ypeH);
  return func;
}


TGraphErrors* rcolRCP() {
  setpoints();
  TVectorF yp(__npoints);
  TVectorF ype(__npoints);
//   yp[0] = 19.5;
//   yp[1] = 11.9;
//   yp[2] = 5.7;
//   yp[3] = 1;
//   yp[0] = 1518.;
//   yp[1] = 933.;
//   yp[2] = 450.;
//   yp[3] = 79.;
//   yp[0] = 1466.; // new not yet approved numbers - March 22 2011
//   yp[1] = 929.;
//   yp[2] = 443.;
//   yp[3] = 78.;
//   yp[0] = 1500.6; // new numbers - April 29 2011
//   yp[1] = 923.3;
//   yp[2] = 440.6;
//   yp[3] = 77.8;
  //May 11 2011
  //2010 Data
/*  yp[0] = 19.3;  
  yp[1] = 11.9; 
  yp[2] = 5.7;  
  yp[3] = 1.0; 
*/
  //2011 Data i/40-80
//  yp[0] = 21.6 ; //0-5 
  /*yp[0] = 19.3 ; //0-10
  yp[1] = 13.3 ;//10-15
  yp[2] = 10.4 ;//15-20
  yp[3] = 5.7 ;//20-40
  yp[4] = 1.0 ; //40-80
 */ 
   yp[0] = 21.6 ; //0-5 
   yp[1] = 16.9 ; //5-10
   yp[2] = 13.3 ;//10-15
   yp[3] = 10.4 ;//15-20
   yp[4] = 5.7 ;//20-40
   yp[5] = 1.0 ; //40-80
  
  //2011 Data i/60-80
  //27 March, 2012
  /*yp[0] = 63.5; //0-5 
  yp[1] = 49.7; //5-10
  yp[2] = 39.1;  //10-15
  yp[3] = 30.6; //15-20
  yp[4] = 16.6; //20-40
  yp[5] = 4.9; //40-60
  yp[6] = 1.0; //60-80
 */
  double err = 1.;
  ype[0] = err;
  ype[1] = err;
  ype[2] = err;
  ype[3] = err;
  ype[4] = err;
  ype[5] = err;
  //ype[6] = err;
  TGraphErrors* func = new TGraphErrors(xp,yp,xpe,ype);
  return func;

}

TGraphErrors* rcolRPC() {
  setpoints();
  TVectorF yp(__npoints);
  TVectorF ype(__npoints);

// -- Run2010 -- //
/*  yp[0] = 1.0; // rpc
  yp[1] = 1.6;
  yp[2] = 3.4;
  yp[3] = 19.3;  
  double err = 1.;
  ype[0] = err;
  ype[1] = err;
  ype[2] = err;
  ype[3] = err;
*/
// --- Run2011 --- //
  yp[0] = 1.0; // rpc
  yp[1] = 1.3; //5-10
  yp[2] = 1.6; //10-15
  yp[3] = 2.1; //15-20
  yp[4] = 3.82; //20-40
  yp[5] = 13.0; //40-60
  //yp[6] = 63.5; //60-80
  double err = 1.;
  ype[0] = err;
  ype[1] = err;
  ype[2] = err;
  ype[3] = err;
  ype[4] = err;
  ype[5] = err;
  //ype[6] = err;


  TGraphErrors* func = new TGraphErrors(xp,yp,xpe,ype);
  return func;

}

TGraphErrors* signExSys() {
  TVectorF yp(__npoints);
  TVectorF ype(__npoints);
  yp[0] = 5.8;
  yp[1] = 1.3;
  yp[2] = 6.6;
  yp[3] = 3.2;
  double err = .1;
  ype[0] = err;
  ype[1] = err;
  ype[2] = err;
  ype[3] = err;
  TGraphErrors* func = new TGraphErrors(xp,yp,xpe,ype);
  return func;

}

TGraphErrors* effSys() {
  TVectorF yp(__npoints);
  TVectorF ype(__npoints);
  yp[0] = 6.2;
  yp[1] = 3.8;
  yp[2] = 2.2;
  yp[3] = 1.6;
  double err = .1;
  ype[0] = err;
  ype[1] = err;
  ype[2] = err;
  ype[3] = err;
  TGraphErrors* func = new TGraphErrors(xp,yp,xpe,ype);
  return func;

}

TGraphErrors* rcolSysRCP() {
  setpoints();
  TVectorF yp(__npoints);
  TVectorF ype(__npoints);
  // --- 2010 --- //
/*  yp[0] = 7.5; // rcp
  yp[1] = 6.3;
  yp[2] = 4.7;
  yp[3] = 0.01;
*/
  // --- 2011 --- //40-80
//  yp[0] = 7.5; // rcp
  yp[0] = 7.5;
  yp[1] = 7.5;
  yp[2] = 6.3;
  yp[3] = 6.3;
  yp[4] = 4.7;
  yp[5] = 0.0;
 
  // --- 2011 --- //60-80
/*  yp[0] = 11.7; //0-5
  yp[1] = 11.6; //5-10
  yp[2] = 10.5; //10-15
  yp[3] = 10.5; //15-20
  yp[4] = 9.6; //20-40
  yp[5] = 5.0; //40-60
  yp[6] = 0.0; //60-80; 
  */
  //yp[6] = 13.2; //60-80; ATLAS-CONF-2011-079
 /* 
  yp[0] = 0; //0-5
  yp[1] = 0; //5-10
  yp[2] = 0; //10-15
  yp[3] = 0; //15-20
  yp[4] = 0; //20-40
  yp[5] =0; //40-60
  yp[6] = 0; //60-80
*/

  double err = .01;
  //ype[0] = err;
  //ype[1] = err;
  //ype[2] = err;
  //ype[3] = err;
  ype[0] = err;
  ype[1] = err;
  ype[2] = err;
  ype[3] = err;
  ype[4] = err;
  ype[5] = err;
  //ype[6] = err;
  
  TGraphErrors* func = new TGraphErrors(xp,yp,xpe,ype);
  return func;

}



TGraphErrors* rcolSysRPC() {
  setpoints();
  TVectorF yp(__npoints);
  TVectorF ype(__npoints);
/*  yp[0] = 0.01;// rpc
  yp[1] = 1.1;
  yp[2] = 2.9;
  yp[3] = 7.5;
  double err = .01;
  ype[0] = err;
  ype[1] = err;
  ype[2] = err;
  ype[3] = err;
*/
  //yp[0] = 7.7;// rpc
  yp[0] = 0.0;// rpc
  yp[1] = 7.8; //5-10
  yp[2] = 7.4; //10-15
  yp[3] = 7.4; //15-20
  yp[4] = 7.4; //20-40
  yp[5] = 7.5; //40-80
 // yp[5] = 9.8; //40-60
  //yp[6] = 11.7; //60-80
  double err = .01;
  ype[0] = err;
  ype[1] = err;
  ype[2] = err;
  ype[3] = err;
  ype[4] = err;
  ype[5] = err;
  //ype[6] = err;

  TGraphErrors* func = new TGraphErrors(xp,yp,xpe,ype);
  return func;

}

TGraphErrors* totSys() {
  TVectorF yp(__npoints);
  TVectorF ype(__npoints);
  yp[0] = 10;
  yp[1] = 6.2;
  yp[2] = 7.7;
  yp[3] = 3.6;
  double err = .1;
  ype[0] = err;
  ype[1] = err;
  ype[2] = err;
  ype[3] = err;
  TGraphErrors* func = new TGraphErrors(xp,yp,xpe,ype);
  std::cout << " got total systematics " << std::endl;
  return func;
}


void draw( int func ) {
 
  // book base histogram
  TH1* h = new TH1F("func","func",100,0,100);
  TGraphErrors* graph = 0;
  if( func == 0 ) graph = yield();
  if( func == 1 ) graph = rcolRCP();
  if( func == 2 ) graph = signExSys();
  if( func == 3 ) graph = effSys();
  if( func == 4 ) graph = rcolSysRCP();
  if( func == 5 ) graph = totSys();

  // get hist range
  double max = 0;
  for( int i=0;i<__npoints;++i ) {
    if( graph->GetY()[i] > max ) max = graph->GetY()[i];
  }
  h->SetMaximum(1.3*max);
  h->SetXTitle("1-Centrality %");  
  TString cname = "func";
  cname += func;
  TCanvas* c;
  c = new TCanvas(cname,cname,500,500);
  h->Draw();
  graph->Draw("PE");
  
}

void sysCheck() {
  TGraphErrors* graph1 = signExSys();
  TGraphErrors* graph2 = effSys();
  TGraphErrors* graph3 = rcolSysRCP();
  TGraphErrors* graph4 = totSys();

  for( int i=0;i<__npoints;++i ) {
    double val = graph1->GetY()[i]*graph1->GetY()[i];  
    val += graph2->GetY()[i]*graph2->GetY()[i];  
    val += graph3->GetY()[i]*graph3->GetY()[i];  
    val = sqrt(val);
    std::cout << " tot sys Fabio: " << graph4->GetY()[i]
	      << " val " << val << std::endl;
  }

}

void drawSystematics() {
  setpoints();
  TH1* h = new TH1F("func","func",100,0,100);
  h->SetMaximum(13);
  h->SetXTitle("1-Centrality %");  
  h->SetYTitle("Systematical Error %");  

  TGraphErrors* graph1 = signExSys();
  TGraphErrors* graph2 = effSys();
  TGraphErrors* graph3 = rcolSysRCP();
  TGraphErrors* graph4 = totSys();

  TString cname = "systematics";
  TCanvas* c = new TCanvas(cname,cname,500,500);
  h->Draw();

  graph1->SetMarkerStyle(20);
  graph2->SetMarkerStyle(20);
  graph3->SetMarkerStyle(20);
  graph4->SetMarkerStyle(20);

  graph2->SetMarkerColor(2);
  graph3->SetMarkerColor(3);
  graph4->SetMarkerColor(4);

  graph2->SetLineColor(2);
  graph3->SetLineColor(3);
  graph4->SetLineColor(4);

  graph1->Draw("P SAME");
  graph2->Draw("P SAME");
  graph3->Draw("P SAME");
  graph4->Draw("P SAME");
  TLegend* leg = new TLegend(0.2,0.65,0.53,0.93);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  
  leg->AddEntry(graph1,"Sign.Ex.Sys.","P");
  leg->AddEntry(graph2,"Eff.Sys.","P");
  leg->AddEntry(graph3,"R_{col} Sys.","P");
  leg->AddEntry(graph4,"Tot. Sys.","P");
  leg->Draw("SAME");

   c->Print(cname+".gif","gif");
   c->Print(cname+".eps","eps");
}
