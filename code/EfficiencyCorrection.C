double range[7] = {2,3,4,5,6,7,9999};
double point0[7] = { 0.0180369, 0.209694, 0.546498, 0.804184, 0.882564, 0.910473, 0.926754};
double point1[7] = { 0.0181116, 0.199667, 0.534592, 0.795927, 0.853862, 0.92936,  0.917184};
double point2[7] = { 0.0181071, 0.202095, 0.53292 , 0.798361, 0.881928, 0.872611, 0.915318};
double point3[7] = { 0.0195342, 0.209801, 0.550089, 0.813898, 0.896122, 0.917946,0.9339};
double ptW[7] = {0., 15., 20., 25., 30., 35., 40.}; 
// double nWMC[7] = {2158693., 1960099., 1783494., 1545781., 1245191., 877533., 469581.}; // mc10 7TeV W MC
//double bkgFit[7] = {0.000976709682906985, 4.93857326511932e-05, 1.47227683658063e-05, 5.36735001850178e-06, 2.22838742066616e-06, 1.0154715998808e-06, 4.9634802745296e-07}; // from data fit, first point is pt>7 (not 0) 
//double bkgFit[7] = {0.000943129, 4.69479e-05, 1.38545e-05, 5.01177e-06, 2.06662e-06, 9.35929e-07, 4.54865e-07}; // from data2011 fit, first point is pt>7 (not 0) 
//double bkgFit[7] = {0.000934801, 5.07704e-05, 1.51085e-05, 5.50527e-06, 2.28537e-06, 1.04136e-06, 5.08968e-07}; // from data2011 fit, first point is pt>7 (not 0) 
//double nWMC[7] = {385053., 351434., 320936., 279377., 225197., 155941., 74137.}; // 2TeV W MC 

//double nWMC[7] = {2.73397, 2.47046, 2.30576, 2.14106, 1.97637, 1.81167, 1.64697}; // 2TeV W MC 
//double nWMC[7] = {8.36152000000000000e+05, 7.85244000000000000e+05, 7.25483000000000000e+05, 6.41633000000000000e+05, 5.29964000000000000e+05, 3.84168000000000000e+05,2.05801000000000000e+05 }; // 2TeV W MC11 
// double bkgFit[7] = {0.00107169636084781, 5.60326128613431e-05, 1.70146394646687e-05, 6.30404397409549e-06, 2.65556117365143e-06, 1.22625320034725e-06, 6.06712871044546e-07}; // fit with April28 data, 2TeV W MC 
// double bkgFit[7] = {0.000979803784013707, 4.9540105512599e-05, 1.47758797428469e-05, 5.39112110411487e-06, 2.23995286719929e-06, 1.02144812670629e-06, 4.99587007391913e-07}; // fit with May2 (repro) data, 2TeV W MC , should first entry be 0.000994579663756554 ?
double etaVal[7] = {-2.5,-2.0,-1.7,-1.05,1.05,2.0,2.5};
double pt_min[999];
double pt_max[999];
double eta_Min[999];
double eta_Max[999];
double centrality_bin[999]; 
double sct_eff[999]; 
double sct_err[999];

void readInputFile(TString sFileIn = "eff_may10.txt"){
//   TFile* fIn = new TFile(sFileIn, "READ");
//   if ( !fIn->IsOpen() ) {
//     std::cout << fIn << " not found!" << std::endl;
//     return;
//   }
  //std::ifstream s(sFileIn);
  std::ifstream s( sFileIn , std::ifstream::in );
  if(!s.is_open()){
    std::cout << "While opening a file an error is encountered" << std::endl;
    exit(1);
  } else {
    std::cout << "File " << sFileIn << " is successfully opened" << std::endl;
  }
//   pt_min[GeV] pt_max[GeV] eta_Min eta_Max centrality  SCT-eff. SCT-eff error
//   4 5 -2.5  -2  1 0.905146  0.00316707
  unsigned int i = 0;
  while(!s.eof())
  {
    s >> pt_min[i] >> pt_max[i] >> eta_Min[i] >> eta_Max[i] >> centrality_bin[i] >> sct_eff[i] >> sct_err[i];
//     cout << i << " " << pt_min[i] << " " << pt_max[i] << " " << eta_Min[i] << " " << eta_Max[i] << " " << centrality_bin[i] << " " << sct_eff[i] << " " << sct_err[i] <<  endl;
    ++i;
  }
}

double efficiencyW( double pt,double sig0,double sig1,double sig2,double sig3,double sig4,double sig5,double sig6 ) {
 int index = -1;
 double nWMC[7] = {sig0,sig1,sig2,sig3,sig4,sig5,sig6};
 for( int i=0;i<7;++i ){
   if( pt <= ptW[i] )  {
     index = i;
     break;
   }
 }
 std::cout << " index " << index << " pt " << pt << std::endl;

 if( index == -1 ) {
   std::cout << " error finding bin " << pt << std::endl;
   return 0;
 }

 return nWMC[index]/nWMC[0];
}

double efficiencyWBkg( double pt ,double bkg0,double bkg1,double bkg2,double bkg3,double bkg4,double bkg5,double bkg6) {
 int index = -1;
 double bkgFit[7] = {bkg0,bkg1,bkg2,bkg3,bkg4,bkg5,bkg6};
 for( int i=0;i<7;++i ){
   if( pt <= ptW[i] )  {
     index = i;
     break;
   }
 }
 std::cout << " index " << index << " pt " << pt << std::endl;

 if( index == -1 ) {
   std::cout << " error finding bin " << pt << std::endl;
   return 0;
 }

 return bkgFit[index]/bkgFit[0];
}


double efficiency( double pt, int centralityBin ) {
 int index = -1;
 for( int i=0;i<7;++i ){
   if( pt < range[i] )  {
     index = i;
     break;
   }
 }
 //std::cout << " index " << i << " pt " << pt << std::endl;

 if( index == -1 ) {
   std::cout << " error finding bin1 " << pt << std::endl;
   return 0;
 }

 if( centralityBin == 0 ) return point0[index];
 if( centralityBin == 1 ) return point1[index];
 if( centralityBin == 2 ) return point2[index];
 if( centralityBin == 3 ) return point3[index];
 std::cout << " error finding bin "<< std::endl;

 return 0;
}

double efficiencyLowPt( double pt, double eta, int centralityBin ) {
  int index = -1;
  for( int i=0;i<999;++i ){
    if ( ( pt > pt_min[i] ) &&  ( pt < pt_max[i] ) 
     && ( eta > eta_Min[i]) && ( eta < eta_Max[i])
     && (centralityBin == centrality_bin[i]-1) )
    {
      index = i;
      break;
    }
  }
//   std::cout << " index " << index << " pt " << pt << std::endl;

 if( index == -1 ) {
//    std::cout << " error finding bin1 " << pt << std::endl;
   return 1;
 } else {
  return sct_eff[index];
 }
}

double efficiencyErrLowPt( double pt, double eta, int centralityBin ) {
  int index = -1;
  for( int i=0;i<999;++i ){
    if ( ( pt > pt_min[i] ) &&  ( pt < pt_max[i] ) 
     && ( eta > eta_Min[i]) && ( eta < eta_Max[i])
     && (centralityBin == centrality_bin[i]-1) )
    {
      index = i;
      break;
    }
  }
//   std::cout << " index " << index << " pt " << pt << std::endl;

 if( index == -1 ) {
//    std::cout << " error finding bin1 " << pt << std::endl;
   return 0;
 } else {
  return sct_err[index];
 }
}

//double efficiencyHighPt( int etaBin, int centralityBin ) {
//
// if( centralityBin == 0 ) return highPtPoint0[index];
// if( centralityBin == 1 ) return highPtPoint1[index];
// if( centralityBin == 2 ) return highPtPoint2[index];
// if( centralityBin == 3 ) return highPtPoint3[index];
// std::cout << " error finding bin "<< std::endl;
//
// return 0;
//}

double integratedEfficiency( TH1* h, int centralityBin ) {

 double effSum = 0;
 double entries = 0;
 for( int i=1;i<=h->GetNbinsX(); ++i ){
   double eff = efficiency( h->GetBinCenter(i), centralityBin );
   //std::cout << " pt " << h->GetBinCenter(i) << " eff " << eff << " entries " << h->GetBinContent(i) << std::endl;
   effSum += h->GetBinContent(i)*eff;
   entries += h->GetBinContent(i);
 }

 return effSum/entries;
}

TH1* correctEfficiency( TH1* h, int centralityBin ) {

 TH1* hnew = (TH1*)h->Clone();
 hnew->SetName(TString("new") + h->GetName());
 for( int i=1;i<=h->GetNbinsX(); ++i ){
   if( h->GetBinContent(i) <= 0 ) hnew->SetBinContent(i,h->GetBinContent(i) );
   else{
     double eff = efficiency( h->GetBinCenter(i), centralityBin );
     hnew->SetBinContent(i,h->GetBinContent(i)/eff );
   }
 }
 return hnew;
}

////////////////////////////////////////////////////
//Aw,Cw corrections
///////////////////////////////////////////////////

const int nCentrality = 6;
const int nEta = 10;

//Cw,Aw taken from reco, trigger, and event selection efficiencies
//Use macro CorrectionFactors.C to see these

const int nBinsCharge = nCentrality*nEta*2;
const int nBins = nCentrality*nEta;
float arrCwChargeSep[nBinsCharge];
float arrCwChargeSepErr[nBinsCharge];
float arrCw[nBins];
float arrCwErr[nBins];
float arrAw[nBinsCharge];
float arrAWxCW[nBinsCharge];

//void setCorrectionFactors(TString sFileIn = "CorrectionFactors.10.09.2012.txt") {
//void setCorrectionFactors(TString sFileIn = "CorrectionFactors.11.22.2012.txt") {
//void setCorrectionFactors(TString sFileIn = "CorrectionFactors.12.09.2012.txt") {
//void setCorrectionFactors(TString sFileIn = "CorrectionFactors.12.24.2012.txt") {
//void setCorrectionFactors(TString sFileIn = "CorrectionFactors.12.25.2012.txt") {
//void setCorrectionFactors( bool doCharge = false, TString sFileIn = "CorrectionFactors.01.12.2013.txt") {

void setCorrectionFactors( bool doCharge = false, TString sFileIn =
    "CorrectionFactorFiles/correctionFactorSpreadSheetRecIDTrig_trkPt3.03.29.2013.txt") {
    
    ////use for systematics ptcone30ID3<0.1
/*void setCorrectionFactors( bool doCharge = false, TString sFileIn =
    "correctionFactorSpreadSheet_ptCone30ID3_systematics.03.31.2013.txt"){
*/
    ////use for systematics ptcone20ID3<0.2
/*void setCorrectionFactors( bool doCharge = false, TString sFileIn =
    "correctionFactorSpreadSheet_ptCone20ID3_systematics.03.31.2013.txt"){
*/    

//void setCorrectionFactors( bool doCharge = false, TString sFileIn = "correctionFactorSpreadSheetRecID.03.27.2013.txt") {
	std::cout << "Opening file " << sFileIn << std::endl;

	std::ifstream s( sFileIn , std::ifstream::in );
	if(!s.is_open()){
		std::cout << "While opening a file an error was encountered" << std::endl;
		exit(1);
        } else {
                std::cout << "File " << sFileIn << " is successfully opened" << std::endl;
        }

	unsigned int i = 0;
	
    std::cout << "Filling arrays.... " << std::endl;
	while(!s.eof()) {
		//s >> arrCw[i] >> arrCwErr[i] >>arrAw[i]; 
        if(doCharge) s >> arrCwChargeSep[i] >> arrCwChargeSepErr[i];
		else s >> arrCw[i] >> arrCwErr[i]; 
		std::cout << "Cw for index " << i << " = " << arrCw[i] << " +- " << arrCwErr[i] << std::endl;
		//arrAWxCW[i] = arrCw[i]*arrAw[i];
		i++;
	}	

        std::cout << "Done." << std::endl;
        //hack
        //exit(0);
}

double getCorrectionFactorError(bool doCharge, int uniqueID, int centralityBin, int etaBin ) {

//NOTE: Efficiencies include unbinned set and binned sets 
//for charge, charge+eta, and charge+centrality
//The index corresponds to the output of CorrectionFactors.C

 int index = -1;

 if(uniqueID==101) return 1.0; //mu+
 //if(uniqueID==101) return arrCW[0]; //mu+
 if(uniqueID==100) return 1.0; //mu-
 //if(uniqueID==100) return arrCW[1]; //mu-
 //if(uniqueID==99)  return arrCW[2];  //mu+- 
 if(uniqueID==99)  return 1.0;  //mu+- 

int nCh = -1;
if(!doCharge) {
    index = centralityBin+etaBin*(centralityBin);
    return arrCwErr[index];
}
else if(doCharge) {
    
    nCh = 2;

  for(int icharge = 0; icharge<nCh; icharge++){
	for (int ieta = 0; ieta<nEta; ieta++){
		for(int icent = 0; icent<nCentrality; icent++){

			//uniqueID = 102 is mu+, 103 = mu-, 104 = mu^{pm}
                index = 2*(icent+ieta*(centralityBin))+icharge;
			    if((uniqueID==102+icharge)&&(centralityBin==icent)&&(etaBin==ieta)) return arrCwChargeSepErr[index];
       }
	}
  }
 }
 std::cout << " WARNING: error finding bin "<< std::endl;

 return 0;
}

double getEfficiencyMt(bool doCharge, int uniqueID, int centralityBin, int etaBin ) {

//NOTE: Efficiencies include unbinned set and binned sets 
//for charge, charge+eta, and charge+centrality
//The index corresponds to the output of CorrectionFactors.C

 int index = -1;

 if(uniqueID==101) return 1.0; //mu+
 //if(uniqueID==101) return arrCW[0]; //mu+
 if(uniqueID==100) return 1.0; //mu-
 //if(uniqueID==100) return arrCW[1]; //mu-
 //if(uniqueID==99)  return arrCW[2];  //mu+- 
 if(uniqueID==99)  return 1.0;  //mu+- 

int nCh = -1;
if(!doCharge) {
    index = centralityBin+etaBin*(centralityBin);
    std::cout << "Cw for index " << index << " = " << arrCw[index] << std::endl;
    return arrCw[index];
}
else if(doCharge) {
    
    nCh = 2;

  for(int icharge = 0; icharge<nCh; icharge++){
	for (int ieta = 0; ieta<nEta; ieta++){
		for(int icent = 0; icent<nCentrality; icent++){

			//uniqueID = 102 is mu+, 103 = mu-, 104 = mu^{pm}
                index = 2*(icent+ieta*(centralityBin))+icharge;
			    if((uniqueID==102+icharge)&&(centralityBin==icent)&&(etaBin==ieta)) return arrCwChargeSep[index];
       }
	}
  }
 }
 std::cout << " WARNING: error finding bin "<< std::endl;

 return 0;
}

double getAcceptanceMt( int chargeID, bool doPtSigmaUp=false, bool doPtSigmaDown=false, bool doMPTSigmaUp=false, bool doMPTSigmaDown=false) {

     //99 = mu^{pm}, 100 = mu-(inclusive), 101 = mu+(inclusive)
     if(doPtSigmaUp){
        std::cout << "Returning systematic fiducial acceptance factor for +1Sigma muonPt." << std::endl;
     	if(chargeID==101) return 0.387397;
     	else if(chargeID==100) return 0.392361;
	else return 0.906;
     }else if(doPtSigmaDown){
        std::cout << "Returning systematic fiducial acceptance factor for -1Sigma muonPt." << std::endl;
     	if(chargeID==101) return 0.64816;
     	else if(chargeID==100) return 0.638914;
	else return 0.906;

     }else if(doMPTSigmaUp){
        std::cout << "Returning systematic fiducial acceptance factor for +1Sigma MPT." << std::endl;
     	if(chargeID==101) return 0.298481;
     	else if(chargeID==100) return 0.326034;
	else return 0.906;

     }else if(doMPTSigmaDown){
        std::cout << "Returning systematic fiducial acceptance factor for -1Sigma MPT." << std::endl;
     	if(chargeID==101) return 0.653929;
     	else if(chargeID==100) return 0.635758;
	else return 0.906;

     }
     else if(chargeID==101) return 0.544; //mu+
     else if(chargeID==100) return 0.536; //mu-
     else return 0.906;
     
}

/*void correctEfficiencyMt( RooPlot* fr, int chargeID, int centralityBin , int etaBin) {

 RooHist* h = (RooHist*) fr->getObject(0) ;

 Double_t x,y;
 //histo bins start at 1; loop over mT bins
 for (Int_t i=0 ; i<h->GetN() ; i++) {
   h->GetPoint(i,x,y) ;
   if( y <= 0 ) h->SetPoint(i,x,y);
   else{
     double eff = applyCorrectionFactors( chargeID, centralityBin , etaBin);
     eff = 0.476834*0.788243; //temporary average of Aw,Cw until eta/centrality bins are one-to-one
     h->SetPoint(i,x,y/eff );
     h->SetPointEYhigh(i,TMath::Sqrt(h->GetN())/eff);
     h->SetPointEYlow(i,TMath::Sqrt(h->GetN())/eff);
   }
 }
 //return h;
}
*/
TH1F* correctAcceptanceMt( TH1F* h, int chargeID, bool bError = true) {

 TH1F* hnew = (TH1F*)h->Clone("hnew");
 hnew->SetName(TString("new") + h->GetName());
 //histo bins start at 1; loop over mT bins
 for (Int_t i=1 ; i<=h->GetNbinsX() ; i++) {
   if( h->GetBinContent(i) <= 0 ) hnew->SetBinContent(i,h->GetBinContent(i) );
   else{

      //99 = mu^{pm}, 100 = mu-(inclusive), 101 = mu+(inclusive)
     //double eff = applyAcceptanceFactors( chargeID, centralityBin , etaBin);
     double eff;
     //if(chargeID==101) eff = 0.957;
     //else if(chargeID==100) eff = 0.854;
     if(chargeID==101) eff = 0.5816;
     else if(chargeID==100) eff = 0.5787;
     else eff = 0.906;

     hnew->SetBinContent(i,h->GetBinContent(i)/eff );
     if(bError) hnew->SetBinError(i,h->GetBinError(i)/eff);
   }
 } 

 return hnew;

 //return h;
}

double getA0(int chargeID, double xEta, bool doCharge){
  if(!doCharge){

	//temp hack with a1 = constant
/*	if(xEta>=0.0&&xEta<=0.25){
		return 0.713877;		
	}
	else if(xEta>=0.25&&xEta<=0.5){
		return 0.848649; 
	}
	else if(xEta>=0.5&&xEta<=0.75){
		return 0.868517;
	}
	else if(xEta>=0.75&&xEta<=1.0){
		return 0.893854;
	}
	else if(xEta>=1.0&&xEta<=1.25){
		return 0.757972;
	
	}
	else if(xEta>=1.25&&xEta<=1.5){
		return 0.781857;
	}
	else if(xEta>=1.5&&xEta<=1.75){
		return 0.766436;
	}
	else if(xEta>=1.75&&xEta<=2.0){
		return 0.765991;
	}
	else if(xEta>=2.0&&xEta<=2.25){
		return 0.660769;
	}
	else if(xEta>=2.25&&xEta<=2.5){
		return 0.594477;
	}
	else {
		std::cout << "WARNING: CANNOT FIND A0. A0=-9999" << std::endl;
		return -9999.;
	}
*/
	if(xEta>=0.0&&xEta<=0.25){
		//return 0.69325;		
		return 0.720296;		
	}
	else if(xEta>=0.25&&xEta<=0.5){
		//return 0.843794; 
		return 0.853614; 
	}
	else if(xEta>=0.5&&xEta<=0.75){
		//return 0.883841;
		return 0.871791;
	}
	else if(xEta>=0.75&&xEta<=1.0){
		//return 0.919016;
		return 0.89553;
	}
	else if(xEta>=1.0&&xEta<=1.25){
		//return 0.775817;
		return 0.757985;
	
	}
	else if(xEta>=1.25&&xEta<=1.5){
		//return 0.783077;
		return 0.78026;
	}
	else if(xEta>=1.5&&xEta<=1.75){
		//return 0.75507;
		return 0.763268;
	}
	else if(xEta>=1.75&&xEta<=2.0){
		//return 0.7536;
		return 0.761311;
	}
	else if(xEta>=2.0&&xEta<=2.25){
		//return 0.654838;
		return 0.652516;
	}
	else if(xEta>=2.25&&xEta<=2.5){
		//return 0.59988;
		return 0.584525;
	}
	else {
		std::cout << "WARNING: CANNOT FIND A0. A0=-9999" << std::endl;
		return -9999.;
	}
  }

  else if(chargeID==102){
	if(xEta>=0.0&&xEta<=0.25){
		return 6.94341e-01;		
	}
	else if(xEta>=0.25&&xEta<=0.5){
		return 8.47020e-01; 
	}
	else if(xEta>=0.5&&xEta<=0.75){
		return 8.86220e-01;
	}
	else if(xEta>=0.75&&xEta<=1.0){
		return 9.27995e-01;
	}
	else if(xEta>=1.0&&xEta<=1.25){
		return 8.07948e-01;
	
	}
	else if(xEta>=1.25&&xEta<=1.5){
		return 8.24931e-01;
	}
	else if(xEta>=1.5&&xEta<=1.75){
		return 7.85357e-01;
	}
	else if(xEta>=1.75&&xEta<=2.0){
		return 7.69695e-01;
	}
	else if(xEta>=2.0&&xEta<=2.25){
		return 6.49063e-01;
	}
	else if(xEta>=2.25&&xEta<=2.5){
		return 5.81055e-01;
	}
	else {
		std::cout << "WARNING: CANNOT FIND A0. A0=-9999" << std::endl;
		return -9999.;
	}
  }
  else if(chargeID==103){
	if(xEta>=0.0&&xEta<=0.25){
		return 6.76635e-01 ;		
	}
	else if(xEta>=0.25&&xEta<=0.5){
		return 8.22981e-01 ;		
	}
	else if(xEta>=0.5&&xEta<=0.75){

		return 8.66772e-01;		
	}
	else if(xEta>=0.75&&xEta<=1.0){

		return 9.17658e-01;		
	}
	else if(xEta>=1.0&&xEta<=1.25){

		return 7.83378e-01;		
	
	}
	else if(xEta>=1.25&&xEta<=1.5){

		return 8.07017e-01;		
	}
	else if(xEta>=1.5&&xEta<=1.75){

		return 7.89193e-01;		
	}
	else if(xEta>=1.75&&xEta<=2.0){

		return 7.70116e-01;		
	}
	else if(xEta>=2.0&&xEta<=2.25){

		return 6.48344e-01;		
	}
	else if(xEta>=2.25&&xEta<=2.5){

		return 5.71289e-01;		

	}
	else {
		std::cout << "WARNING: CANNOT FIND A0. A0=-9999" << std::endl;
		return -9999.;
	}
  }
	else {
		std::cout << "WARNING: CANNOT FIND A0. A0=-9999" << std::endl;
		return -9999.;
	}
}

double getA1(int chargeID, double xEta){

/*	 double b0 = (-1.22747e-03+-1.15866e-03)/2.0;
	 double b1 = (-5.47588e-05+-2.42778e-05)/2.0;
	 double b2 = (1.24201e-04+1.09511e-04)/2.0;
	 double b3 = (1.15664e-04+8.42216e-05)/2.0;
	 double a1 = b0+b1*TMath::Cos(xEta)+b2*TMath::Cos(2.0*xEta)+b3*TMath::Cos(3.0*xEta);
*/

/*	 double b0 = -1.17977e-03;
	 double b1 = 7.08965e-05;
	 double b2 = 6.04041e-05;
	 double a1 = b0+b1*TMath::Cos(3.0*xEta)+b2*TMath::Cos(4.0*xEta);
*/
	 double b0 = -1.208e-3;
	 double b1 = 3.387e-5;
	 //systematics +-1sigma
	 //double b0 = -1.208e-3+3.252e-05;
	 //double b0 = -1.208e-3-3.252e-05;
	 //double b1 = 3.387e-5+2.448e-05;
	 //double b1 = 3.387e-5-2.448e-05;
	 double a1 = b0 + b1*xEta;
	 //temp hack
	 //a1 = -0.00116983;
	 return a1;

	/*if(chargeID==102) return -7.668e-4+1.219e-5*xEta*xEta;
	else if(chargeID==103) return -8.437e-4+3.879e-5*xEta*xEta;
	else {
		std::cout << "WARNING: CANNOT FIND A1. A1=-9999" << std::endl;
		return -9999.;
	}*/
}

//const coeff as fcn of eta
double getA2(){
	
	double a2 = 5.72135e-07;
	//systematics +-1sigma
	//double a2 = 5.72135e-07 +1.729e-07;
	//double a2 = 5.72135e-07 -1.729e-07;
	//return (4.56612e-07+7.19152e-07)/2.0;
	return a2;
}

//Cw parametrized as fcn of eta and centrality
double getParamCorrectionFactor(double a0, double a1, double a2, double npart){

		return a0 + a1*npart + a2*npart*npart; 

}


//used for fitting Cw
double getEfficiencyFitMt(int chargeID, int centralityBin , double xEta) {

 std::vector <double> npartBins;
 npartBins.push_back(382.16);//0-5
 npartBins.push_back(330.26);//5-10
 npartBins.push_back(281.88);//10-15
 npartBins.push_back(239.52);//15-20
 npartBins.push_back(157.83);//20-40
 npartBins.push_back(45.93);//40-80

 double a0 = getA0(chargeID,xEta,false);
 double a1 = getA1(chargeID,xEta);
 double a2 = getA2();
     //double eff = applyCorrectionFactors( chargeID, centralityBin , etaBin);
 double eff = getParamCorrectionFactor(a0,a1,a2,npartBins.at(centralityBin));
 std::cout << "Cw correction for bin " << chargeID << ":" << xEta << ":" << centralityBin << 
 		" = " << a0 << "+" << a1 << "(" << npartBins.at(centralityBin) << ")" << "+" << a2 << "(" << npartBins.at(centralityBin)*npartBins.at(centralityBin) << ")" 
		<< " = " << eff << std::endl;

 return eff;
     
}

TH1F* correctEfficiencyMt( TH1F* h, int chargeID, int centralityBin, int etaBin , double Cw) {

 TH1F* hnew = (TH1F*)h->Clone("hnew");
 hnew->SetName(TString("new") + h->GetName());
 //histo bins start at 1; loop over mT bins
 std::cout << "Number of bins in input histogram: " << h->GetNbinsX() << std::endl;
 for (Int_t i=1 ; i<=h->GetNbinsX() ; ++i) {
   /*if( h->GetBinContent(i) <= 0 ) { 
       //std::cout << "WARNING: Input histogram bin content negative or zero." << std::endl;
       hnew->SetBinContent(i,h->GetBinContent(i) );
   }*/

     //return efficiency and stat. error of efficiency factor
     //double eff = applyCorrectionFactors( chargeID, npartBins.at(centralityBin) , etaBin);
     //double eff = getEfficiencyFitMt(chargeID, centralityBin, xEta);
     double eff = Cw;
     //double errEff = getCorrectionFactorError( chargeID, centralityBin , etaBin);
     //double eff = getEfficiencyMt( chargeID, centralityBin , etaBin);
     if(i==1) std::cout << "Cw for bin " << ":" << etaBin << ":" << centralityBin << " = " << eff << std::endl;
     //double errEff = 0.0;

     //double err1 = TMath::Sqrt(h->GetBinContent(i))/h->GetBinContent(i)*100.0;
     //double err2 = errEff/eff*100.0; 
     //scale the histogram
     std::cout << "Bin content BEFORE correction = " << h->GetBinContent(i) << std::endl;
     hnew->SetBinContent(i,h->GetBinContent(i)/eff );
     double error = TMath::Sqrt(h->GetBinContent(i))/eff;
     hnew->SetBinError(i,error);
     std::cout << "Bin content AFTER correction = " << hnew->GetBinContent(i) << "+/-" << hnew->GetBinError(i) << std::endl;
 } 

 return hnew;

 //return h;
}
