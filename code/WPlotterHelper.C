#include "TriggerEfficiencies.C"

//////////////////////////////////////////
//Helper functions for WAnalysis
/////////////////////////////////////////

//////////////////////////////////////////
///Weight W MC ds by N-N collision probability
/////////////////////////////////////////
/*RooDataSet* weightDS( RooDataSet* ds, RooRealVar& weight, double wt){

    ds->addColumn(weight);
    weight.setVal(wt);
    //weight = wt;
    ds->setWeightVar(weight.GetName());
    return ds;
}
*/
//////////////////////////////////////////
//return ratio of hMC over hData
/////////////////////////////////////////
TH1F* getDataMCRatio(TH1F* hRatio, TH1F* hnum, TH1F* hden){
        std::cout << "Calculation Data/MC ratios..." << std::endl;
		hRatio->Divide(hnum,hden);
        std::cout << "Done." << std::endl;
		return hRatio;
}
//////////////////////////////////////////
//Return bin number for x in UNIFORMLY
//binned histogram
/////////////////////////////////////////
int getBinNumber(const float x, int nBins, float xLo, const float xUp){

    ///returns bin number of x for
    ///UNIFORM binning
    double binW = (double)(xUp-xLo)/nBins; 
    int binLo(0);
    for(int i=0; i<nBins; ++i){
        if(x<(xLo+i*binW)) return binLo; 
        ++binLo;
    }

    if(binLo>nBins) {
        std::cout << "ERROR: Cannot find bin number. " << std::endl;
        return -9999;
    }
    return binLo;


}///get bin number
///////////////////////////////////////////
//writeToSpreadsheet
//////////////////////////////////////////
void writeToSpreadsheet(std::ostream& outputFile, int ieta, int icent, double sigCounts, double errStat, double errSyst){
	outputFile << ieta << "," << icent << "," << sigCounts << "," << errStat << "," << errSyst << std::endl;
}
///////////////////////////////////////////
//writeToSpreadsheet
//////////////////////////////////////////
void writeToSpreadsheet(std::ostream& outputFile, TString ieta, TString icent, TString sigCounts, TString sigCountsRaw, TString sigCountsBS,                
                            TString bkgCount1, TString bkgCounts2, TString bkgCounts3,
                            TString errStat, TString syst1, TString syst2, TString syst3, TString syst4, 
                            TString syst5, TString syst6, TString syst7,TString syst8,TString syst9, TString errSystQuad,TString errSystUncorrelated){

	outputFile << ieta << "," << icent << "," << sigCounts << "," << sigCountsRaw << "," <<  sigCountsBS << "," 
            << bkgCount1 << "," << bkgCounts2 << "," << bkgCounts3 << ","
            << errStat << ","
            << syst1 << "," << syst2 << "," << syst3 << "," << syst4 << "," << 
            syst5 << "," << syst6 << "," << syst7 << "," << syst8 << "," << syst9 << "," 
            << errSystQuad << "," << errSystUncorrelated << std::endl;
}
///////////////////////////////
//used for mapping negative and positive eta bins
//to absolute eta bin values
//////////////////////////////
int indexPositiveEta(int ieta, int nEtaBins){
      int index = ieta-1.0/2.0*(nEtaBins+1.0);
      std::cout << "indexPositiveEta:" << index << std::endl;
      if(index>(nEtaBins-1)/2) std::cout << "WARNING: index out of allowed range. EXPECT WRONG RESULTS."<<std::endl;
      return index;
}
int indexNegativeEta(int ieta, int nEtaBins){
     int index = 0.5*(nEtaBins+1.0)-2.0; index-=ieta;
     std::cout << "indexNegativeEta:" << index << std::endl;
     if(index>(nEtaBins-1)/2) std::cout << "WARNING: index out of allowed range. EXPECT WRONG RESULTS."<<std::endl;
     return index;
}

///////////////////////////////
//return absolute eta bin number
//////////////////////////////
int getTauBkgBinNumber(int ieta, float eta, int nEtaBins, bool doMirrorEta=true){

    ///cover binning in more than 10 bins
    if(nEtaBins<10) return ieta;
    int index = 0;
    ///for MIRROR eta, performing mapping
    ///of negative to positive
    if(!doMirrorEta||nEtaBins>9){
      
      std::cout << "No. of eta bins >9. Will shift bin indices." << std::endl;

      ///hop over the crack region
      ///in bin [-0.1,0.1]
      if(eta==0.0) std::cout << "WARNING: Bin in crack region. Expect wrong result!" << std::endl;
      if(eta<0.0){
          ///number of bins in absolute eta
          return indexNegativeEta(ieta,nEtaBins);
      }
       else {
          return indexPositiveEta(ieta,nEtaBins);
       }
    }
 
    if(fabs(eta)>0.0&&fabs(eta)<0.25) return index; else ++index;
    if(fabs(eta)>0.25&&fabs(eta)<0.5) return index; else ++index;
    if(fabs(eta)>0.5&&fabs(eta)<0.75) return index; else ++index;
    if(fabs(eta)>0.75&&fabs(eta)<1.0) return index; else ++index;
    if(fabs(eta)>1.0&&fabs(eta)<1.25) return index; else ++index;
    if(fabs(eta)>1.25&&fabs(eta)<1.5) return index; else ++index;
    if(fabs(eta)>1.5&&fabs(eta)<1.75) return index; else ++index;
    if(fabs(eta)>1.75&&fabs(eta)<2.0) return index; else ++index;
    if(fabs(eta)>2.0&&fabs(eta)<2.25) return index; else ++index;
    if(fabs(eta)>2.25&&fabs(eta)<2.5) return index; 
    else {std::cout << "ERROR INDEXING Z BKG. EXPECT WRONG RESULTS"
                << std::endl;
                return -1;
         }
} //getTauBkgBinNumber


///////////////////////////////
//return absolute eta bin number
//////////////////////////////
int getZBkgBinNumber(int ieta, float eta, int nEtaBins, bool doMirrorEta=true){

    ///cover binning in more than 10 bins
    if(nEtaBins<10) return ieta;
    int index = 0;
    ///for MIRROR eta, performing mapping
    ///of negative to positive
    if(!doMirrorEta||nEtaBins>9){
      
      std::cout << "No. of eta bins >9. Will shift bin indices." << std::endl;

      ///hop over the crack region
      ///in bin [-0.1,0.1]
      if(eta==0.0) std::cout << "WARNING: Bin in crack region. Expect wrong result!" << std::endl;
      if(eta<0.0){
          ///number of bins in absolute eta
          return indexNegativeEta(ieta,nEtaBins);
      }
       else {
          return indexPositiveEta(ieta,nEtaBins);
       }
    }
 
    if(fabs(eta)>0.0&&fabs(eta)<0.25) return index; else ++index;
    if(fabs(eta)>0.25&&fabs(eta)<0.5) return index; else ++index;
    if(fabs(eta)>0.5&&fabs(eta)<0.75) return index; else ++index;
    if(fabs(eta)>0.75&&fabs(eta)<1.0) return index; else ++index;
    if(fabs(eta)>1.0&&fabs(eta)<1.25) return index; else ++index;
    if(fabs(eta)>1.25&&fabs(eta)<1.5) return index; else ++index;
    if(fabs(eta)>1.5&&fabs(eta)<1.75) return index; else ++index;
    if(fabs(eta)>1.75&&fabs(eta)<2.0) return index; else ++index;
    if(fabs(eta)>2.0&&fabs(eta)<2.25) return index; else ++index;
    if(fabs(eta)>2.25&&fabs(eta)<2.5) return index; 
    else {std::cout << "ERROR INDEXING Z BKG. EXPECT WRONG RESULTS"
                << std::endl;
                return -1;
         }
} //getZBkgBinNumber


///////////////////////////////
//return QCD background for bins 
//with good statistics
//////////////////////////////
int getQCDBkgBinNumber(int ieta, float eta, int nEtaBins, bool doMirrorEta=true){

    ///cover binning in more than 10 bins
    if(nEtaBins<10) return ieta;
    int index = 0;
    ///for MIRROR eta, performing mapping
    ///of negative to positive
    if(!doMirrorEta||nEtaBins>9){
      
      std::cout << "No. of eta bins >9. Will shift bin indices." << std::endl;

      ///hop over the crack region
      ///in bin [-0.1,0.1]
      if(eta==0.0) std::cout << "WARNING: Bin in crack region. Expect wrong result!" << std::endl;
      if(eta<0.0){
          ///number of bins in absolute eta
          return indexNegativeEta(ieta,nEtaBins);
      }
       else {
          return indexPositiveEta(ieta,nEtaBins);
       }
    }
 
    if(fabs(eta)>0.0&&fabs(eta)<0.25) return index; else ++index;
    if(fabs(eta)>0.25&&fabs(eta)<0.5) return index; else ++index;
    if(fabs(eta)>0.5&&fabs(eta)<0.75) return index; else ++index;
    if(fabs(eta)>0.75&&fabs(eta)<1.0) return index; else ++index;
    if(fabs(eta)>1.0&&fabs(eta)<1.25) return index; else ++index;
    if(fabs(eta)>1.25&&fabs(eta)<1.5) return index; else ++index;
    if(fabs(eta)>1.5&&fabs(eta)<1.75) return index; else ++index;
    if(fabs(eta)>1.75&&fabs(eta)<2.0) return index; else ++index;
    if(fabs(eta)>2.0&&fabs(eta)<2.25) return index; else ++index;
    if(fabs(eta)>2.25&&fabs(eta)<2.5) return index; 
    else {std::cout << "ERROR INDEXING QCD BKG. EXPECT WRONG RESULTS"
                << std::endl;
                return -1;
         }
} //getQCDBkgBinNumber



///////////////////////////////////////////
//return number of HI collisions in sample
//////////////////////////////////////////
double getNEvents(RooDataSet* ds, double centralityLow, double centralityUpp){
	//event counting dataset
    RooDataSet* mcEventsTemp = selectCentrality(ds,centralityLow,centralityUpp); mcEventsTemp->Print();
	return mcEventsTemp->numEntries();
}


///////////////////////////////////////////
//return Ncoll weighted W histogram
//////////////////////////////////////////
TH1F* getWeightedWHisto(double centralityLow,double centralityUpp, double ncoll, TH1F* hmc1,
			RooDataSet* mcNEvents, int nBins, double xLo, double varMax){
	std::cout << "Weighting W sample..." << std::endl;
	std::cout << "centrality: " << centralityLow << "-" << centralityUpp << " ncoll: " << ncoll << " nBins: " << nBins << " varMax " << varMax << std::endl;
	//return the number of MC collisions in this centrality bin
	double nMC = getNEvents(mcNEvents, centralityLow, centralityUpp);

	//TH1F* hWtd = new TH1F("hWtd","hWtd",nBins,xLo,varMax);
	//weight the sample according to pp cross-sections (AMI)
	double crossSec = 2.8E-09/64.0e-3*0.95;  
	double arrCentWidth = centralityUpp-centralityLow;
	//double evData = 68.7e6; //number of minbias events
	double evData = 1.03e9; //number of hp events
	double scaleFactor = arrCentWidth*ncoll*evData;
    ///Scale by cross section
    //hmc1->Scale(crossSec);
    ///Scale to number of collisions in all events of this centrality class
    hmc1->Scale(crossSec*scaleFactor/nMC);
    return hmc1;

}

///////////////////////////////////////////
//return Ncoll weighted Wtau histogram
//////////////////////////////////////////
TH1F* getWeightedTauHisto(double centralityLow,double centralityUpp, double ncoll, TH1F* hmc1,
			RooDataSet* mcNEvents, int nBins, double xLo, double varMax){
	std::cout << "Weighting W-->tau sample..." << std::endl;
	std::cout << "centrality: " << centralityLow << "-" << centralityUpp << " ncoll: " << ncoll << " nBins: " << nBins << " varMax " << varMax << std::endl;
	//return the number of MC collisions in this centrality bin
	double nMC = getNEvents(mcNEvents, centralityLow, centralityUpp);

	//weight the sample according to pp cross-sections (AMI)
	double crossSec = 2.817E-09/64.0e-3*0.13754;  
	double arrCentWidth = centralityUpp-centralityLow;
	//double evData = 68.7e6; //number of minbias events
	double evData = 1.03e9; //number of hp events
	double scaleFactor = arrCentWidth*ncoll*evData;
    ///Scale by cross section
    //hmc1->Scale(crossSec);
    ///Scale to number of collisions in all events of this centrality class
    hmc1->Scale(crossSec*scaleFactor/nMC);
    return hmc1;

}


///////////////////////////////////////////
//return Ncoll weighted Z histogram
//////////////////////////////////////////

TH1F* getWeightedZHisto(double centralityLow,double centralityUpp, double ncoll, TH1F* hmc1,
			RooDataSet* mcNEvents, int nBins, double xLo, double varMax){
	std::cout << "Weighting Z sample..." << std::endl;
	std::cout << "centrality: " << centralityLow << "-" << centralityUpp << " ncoll: " << ncoll << " nBins: " << nBins << " varMax " << varMax << std::endl;
	//return the number of MC collisions in this centrality bin
	double nMC = getNEvents(mcNEvents, centralityLow, centralityUpp);

	//TH1F* hWtd = new TH1F("hWtd","hWtd",nBins,xLo,varMax);
	//weight the sample according to pp cross-section (AMI)
	double crossSec = 2.6E-10/64.0e-3*9.99E-01 ; 
	double arrCentWidth = centralityUpp-centralityLow;
	//double evData = 68.7e6; //number of minbias events
	double evData = 1.03e9; //number of hp events
	double scaleFactor = arrCentWidth*ncoll*evData;
    ///Scale by cross section
    //hmc1->Scale(crossSec);
    ///Scale to number of collisions in all events of this centrality class
    hmc1->Scale(scaleFactor/nMC*crossSec);
    return hmc1;

}

///////////////////////////////////////////
//return weighted Jx histogram(currently only J1-J3) 
//with variable binning
//////////////////////////////////////////
TH1F* getWeightedJxHisto(double centralityLow,double centralityUpp, double ncoll, TH1F* hmcJ1Set, TH1F* hmcJ2Set, TH1F* hmcJ3Set, 
			RooDataSet* mcJ1Events,RooDataSet* mcJ2Events, RooDataSet* mcJ3Events, int nBins, const double* xBins=0){

	std::cout << "Weighting Jx samples..." << std::endl;
	//return the number of MC collisions in this centrality bin
	double nMCJ1 = getNEvents(mcJ1Events, centralityLow, centralityUpp);
	double nMCJ2 = getNEvents(mcJ2Events, centralityLow, centralityUpp);
	double nMCJ3 = getNEvents(mcJ3Events, centralityLow, centralityUpp);

	TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",nBins,xBins);
	//weight the Jx samples according to cross-sections; takes into account prob of find muon+jet (AMI)
	double wtJ1 = 6.6523e-03*1.8770e-4/64.0e-3; double wtJ2 = 1.4941e-02*8.2788e-6/64.0e-3; double wtJ3 = 2.4284e-02*2.9426e-7/64.0e-3;
	double arrCentWidth = centralityUpp-centralityLow;
	//double evData = 68.7e6; //number of minbias events
	double evData = 1.03e9; //number of hp events
	double scaleFactor = arrCentWidth*ncoll*evData;

	hmcQCDSet->Add(hmcJ1Set,hmcJ2Set,(wtJ1)/(nMCJ1)*scaleFactor, (wtJ2)/(nMCJ2)*scaleFactor); 
	hmcQCDSet->Add(hmcJ3Set,(wtJ3)/(nMCJ3)*scaleFactor);

	std::cout << "J1 coeff : " << (wtJ1)/(nMCJ1)*scaleFactor << " J2 coeff : " << (wtJ2)/(nMCJ2)*scaleFactor << " J2 coeff : " << (wtJ3)/(nMCJ3)*scaleFactor << std::endl;
	std::cout << "MCevents: " << nMCJ1 << " " << nMCJ2 << " " << nMCJ3 << std::endl;

	return hmcQCDSet;
}


///////////////////////////////////////////
//return weighted Jx histogram(currently only J1-J3) 
//////////////////////////////////////////
TH1F* getWeightedJxHisto(double centralityLow,double centralityUpp, double ncoll, TH1F* hmcJ1Set, TH1F* hmcJ2Set, TH1F* hmcJ3Set, 
			RooDataSet* mcJ1Events,RooDataSet* mcJ2Events, RooDataSet* mcJ3Events, int nBins, double xLo, double varMax){

	std::cout << "Weighting Jx samples..." << std::endl;
	std::cout << "centrality: " << centralityLow << "-" << centralityUpp << " ncoll: " << ncoll << " nBins: " << nBins << " varMax " << varMax << std::endl;
	//return the number of MC collisions in this centrality bin
	double nMCJ1 = getNEvents(mcJ1Events, centralityLow, centralityUpp);
	double nMCJ2 = getNEvents(mcJ2Events, centralityLow, centralityUpp);
	double nMCJ3 = getNEvents(mcJ3Events, centralityLow, centralityUpp);

	TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",nBins,xLo,varMax);
	//weight the Jx samples according to cross-sections; takes into account prob of find muon+jet (AMI)
	double wtJ1 = 6.6523e-03*1.8770e-4/64.0e-3; double wtJ2 = 1.4941e-02*8.2788e-6/64.0e-3; double wtJ3 = 2.4284e-02*2.9426e-7/64.0e-3;
	double arrCentWidth = centralityUpp-centralityLow;
	//double evData = 68.7e6; //number of minbias events
	double evData = 1.03e9; //number of hp events
	double scaleFactor = arrCentWidth*ncoll*evData;

	hmcQCDSet->Add(hmcJ1Set,hmcJ2Set,(wtJ1)/(nMCJ1)*scaleFactor, (wtJ2)/(nMCJ2)*scaleFactor); 
	hmcQCDSet->Add(hmcJ3Set,(wtJ3)/(nMCJ3)*scaleFactor);

	std::cout << "J1 coeff : " << (wtJ1)/(nMCJ1)*scaleFactor << " J2 coeff : " << (wtJ2)/(nMCJ2)*scaleFactor << " J2 coeff : " << (wtJ3)/(nMCJ3)*scaleFactor << std::endl;
	std::cout << "MCevents: " << nMCJ1 << " " << nMCJ2 << " " << nMCJ3 << std::endl;

	return hmcQCDSet;
}

///////////////////////////////////////////
//return weighted Jx histogram for variable-binned input(currently only J1-J3) 
//////////////////////////////////////////
TH1F* getWeightedJxVariableBinnedHisto(double xBins[], int arrLength, double centralityLow,double centralityUpp, double ncoll, TH1F* hmcJ1Set, TH1F* hmcJ2Set, TH1F* hmcJ3Set, 
			RooDataSet* mcJ1Events,RooDataSet* mcJ2Events, RooDataSet* mcJ3Events, int nBins, double xLo, double varMax){

	std::cout << "Weighting Jx samples..." << std::endl;
	std::cout << "centrality: " << centralityLow << "-" << centralityUpp << " ncoll: " << ncoll << " nBins: " << nBins << " varMax " << varMax << std::endl;
	//return the number of MC collisions in this centrality bin
	double nMCJ1 = getNEvents(mcJ1Events, centralityLow, centralityUpp);
	double nMCJ2 = getNEvents(mcJ2Events, centralityLow, centralityUpp);
	double nMCJ3 = getNEvents(mcJ3Events, centralityLow, centralityUpp);

    //double xBins2[20] =
    //{-2.4,-2.1,-1.85,-1.55,-1.3,-1.05,-0.8,-0.6,-0.35,-0.1,0.1,0.35,0.6,0.8,1.05,1.3,1.55,1.85,2.1,2.4}; 
	TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",arrLength,xBins);

	//weight the Jx samples according to cross-sections; takes into account prob of find muon+jet (AMI)
	double wtJ1 = 6.6523e-03*1.8770e-4/64.0e-3; double wtJ2 = 1.4941e-02*8.2788e-6/64.0e-3; double wtJ3 = 2.4284e-02*2.9426e-7/64.0e-3;
	double arrCentWidth = centralityUpp-centralityLow;
	//double evData = 68.7e6; //number of minbias events
	double evData = 1.03e9; //number of hp events
	double scaleFactor = arrCentWidth*ncoll*evData;

    std::cout << "Bins must be equal. Let's check: " << hmcQCDSet->GetNbinsX() << " =? " << hmcJ1Set->GetNbinsX() << " =? "
        << hmcJ2Set->GetNbinsX() << " =? " << hmcJ3Set->GetNbinsX() << " =? " << std::endl;
	hmcQCDSet->Add(hmcJ1Set,hmcJ2Set,(wtJ1)/(nMCJ1)*scaleFactor, (wtJ2)/(nMCJ2)*scaleFactor); 
	hmcQCDSet->Add(hmcJ3Set,(wtJ3)/(nMCJ3)*scaleFactor);

	std::cout << "J1 coeff : " << (wtJ1)/(nMCJ1)*scaleFactor << " J2 coeff : " << (wtJ2)/(nMCJ2)*scaleFactor << " J2 coeff : " << (wtJ3)/(nMCJ3)*scaleFactor << std::endl;
	std::cout << "MCevents: " << nMCJ1 << " " << nMCJ2 << " " << nMCJ3 << std::endl;

	return hmcQCDSet;
}
///////////////////////////////////////////
//correctly weight the Jx samples
//////////////////////////////////////////
double weightedJxEvents(TTree* treeMcJ1Set, TTree* treeMcJ2Set, TTree* treeMcJ3Set, TString cutsQCD, double ncoll, double centralityLow, double centralityUpp, double mtcutLow, double mtmax ){
	
	int nBins = 50;
	// --- QCD set ---
	TH1F* hmcQCDTemp = new TH1F("hmcQCDTemp","hmcQCDTemp",nBins,0.0,mtmax);
	TH1F* hmcJ1Temp = new TH1F("hmcJ1Temp","hmcJ1Temp",nBins,0.0,mtmax);
	TH1F* hmcJ2Temp = new TH1F("hmcJ2Temp","hmcJ2Temp",nBins,0.0,mtmax);
	TH1F* hmcJ3Temp = new TH1F("hmcJ3Temp","hmcJ3Temp",nBins,0.0,mtmax);
	TH1F* hmcJ1CentTemp = new TH1F("hmcJ1CentTemp","hmcJ1CentTemp",nBins,0.0,mtmax);
	TH1F* hmcJ2CentTemp = new TH1F("hmcJ2CentTemp","hmcJ2CentTemp",nBins,0.0,mtmax);
	TH1F* hmcJ3CentTemp = new TH1F("hmcJ3CentTemp","hmcJ3CentTemp",nBins,0.0,mtmax);

	treeMcJ1Set->Draw("mt>>hmcJ1Temp",cutsQCD,"hf");		
	treeMcJ2Set->Draw("mt>>hmcJ2Temp",cutsQCD,"hf");		
	treeMcJ3Set->Draw("mt>>hmcJ3Temp",cutsQCD,"hf");		

	//event counting histos
	TString sCentrality = "centrality>"; sCentrality+= centralityLow; sCentrality+="&&centrality<"; sCentrality+=centralityUpp;
	treeMcJ1Set->Draw("centrality>>hmcJ1CentTemp",sCentrality,"hf");
	treeMcJ2Set->Draw("centrality>>hmcJ2CentTemp",sCentrality,"hf");
	treeMcJ3Set->Draw("centrality>>hmcJ3CentTemp",sCentrality,"hf");
	double nMCJ1 = hmcJ1CentTemp->Integral(); double nMCJ2 = hmcJ2CentTemp->Integral(); double nMCJ3 = hmcJ3CentTemp->Integral();

	//weight the Jx samples according to cross-sections; takes into account prob of find muon+jet (AMI)
	double wtJ1 = 6.6523e-03*1.8770e-4/64.0e-3; double wtJ2 = 1.4941e-02*8.2788e-6/64.0e-3; double wtJ3 = 2.4284e-02*2.9426e-7/64.0e-3;
	double arrCentWidth = centralityUpp-centralityLow;
	double evData = 1.03e9; //number of hp events
	double scaleFactor = arrCentWidth*ncoll*evData;

	hmcQCDTemp->Add(hmcJ1Temp,hmcJ2Temp,(wtJ1)/(nMCJ1)*scaleFactor, (wtJ2)/(nMCJ2)*scaleFactor); 
	hmcQCDTemp->Add(hmcJ3Temp,(wtJ3)/(nMCJ3)*scaleFactor);

//	int binLo = nBins/mtmax*mtcutLow+1;
    int binLo = getBinNumber(mtcutLow,nBins,0.0,mtmax);
  	double events = hmcQCDTemp->Integral(binLo,nBins); //integrate from 40-200 GeV before Aw,Cw correction
	return events;

}

///////////////////////////////////////////
//returns fraction of Z background for given
//Z percentage in signal region determined from MC
//studies 
//////////////////////////////////////////
double getCorrectedSignalYield(double yieldUC, int iMt, int iCentrality, int iEta){

	double eff ;
	//eff = getEfficiencyMt(iMt, iCentrality, iEta);
	eff = getEfficiencyMt(false,104, iCentrality, iEta);
	return yieldUC/eff;
}

///////////////////////////////////////////
//counts in signal region given by scutsSig
//////////////////////////////////////////
double getSigEvents(TTree *tree, double mtcutLow, double mtmax, TString scutsSig ){

  	std::cout << "Calculating signal events for region: "<< scutsSig << std::endl;
	int nBins = 50;
	TH1F* hTemp = new TH1F("hTemp","hTemp",nBins,0.0,mtmax);
	// --- Fill the histograms ---
	tree->Draw("mt>>hTemp",scutsSig,"pe");		

//	int binLo = nBins/mtmax*mtcutLow+1;
    int binLo = getBinNumber(mtcutLow,nBins,0.0,mtmax);
  	double sigEvents = hTemp->Integral(binLo,nBins); //integrate from 40-200 GeV before Aw,Cw correction

	return sigEvents;

}

double getZEtaCentSliceMC(double centralityLow, double centralityUpp, double ncoll,double frac, RooDataSet* mcNEvents){

	    double nMC = getNEvents(mcNEvents, centralityLow, centralityUpp);
	    double crossSec = 2.6E-10/64.0e-3*9.99E-01 ; 
	    double arrCentWidth = centralityUpp-centralityLow;
	    double evData = 1.03e9; //number of hp events
	    double scaleFactor = arrCentWidth*ncoll*evData;

        ///Yields in each centrality bin taken from Z MC
		if(centralityLow>=0.0&&centralityUpp<=0.05) {
			return frac*(scaleFactor/nMC*crossSec)*17470.;
            std::cout << frac << "*" << "338.0 Zs in centrality " << centralityLow << centralityUpp << std::endl;
		} else if (centralityLow>=0.05&&centralityUpp<=0.1){ 
			return frac*(scaleFactor/nMC*crossSec)*24974.;
		} else if (centralityLow>=0.1&&centralityUpp<=0.15) {
			return frac*(scaleFactor/nMC*crossSec)*25135.;
		} else if (centralityLow>=0.15&&centralityUpp<=0.20) {
			return frac*(scaleFactor/nMC*crossSec)*24918.;
		} else if (centralityLow>=0.20&&centralityUpp<=0.40) {
			return frac*(scaleFactor/nMC*crossSec)*104028.;
		} else if (centralityLow>=0.40&&centralityUpp<=0.80) {
			return frac*(scaleFactor/nMC*crossSec)*234548.;
		} else {
			return -1.0;
		} 
	
}

double getZEtaCentSliceData(double centralityLow, double centralityUpp, double frac){

		if(centralityLow>=0.0&&centralityUpp<=0.05) {
			return frac*349.0;
            std::cout << frac << "*" << "338.0 Zs in centrality " << centralityLow << centralityUpp << std::endl;
		} else if (centralityLow>=0.05&&centralityUpp<=0.1){ 
			return frac*320.0;
		} else if (centralityLow>=0.1&&centralityUpp<=0.15) {
			return frac*251.0; 
		} else if (centralityLow>=0.15&&centralityUpp<=0.20) {
			return frac*201.0;
		} else if (centralityLow>=0.20&&centralityUpp<=0.40) {
			return frac*411.0;
		} else if (centralityLow>=0.40&&centralityUpp<=0.80) {
			return frac*140.0;
		} else {
			return -1.0;
		} 
	
}

//////////////////////////////////////////////
//Return [%] stat error from
//Zs found in the data in this centrality bin
///////////////////////////////////////////////
double getZCentDataStatErr(double centralityLow, double centralityUpp){

		if(centralityLow>=0.0&&centralityUpp<=0.05) {
			return TMath::Sqrt(349.0)/349.0;
		} else if (centralityLow>=0.05&&centralityUpp<=0.1){ 
			return TMath::Sqrt(320.0)/320.0;
		} else if (centralityLow>=0.1&&centralityUpp<=0.15) {
			return TMath::Sqrt(251.0)/251.0; 
		} else if (centralityLow>=0.15&&centralityUpp<=0.20) {
			return TMath::Sqrt(201.0)/201.0;
		} else if (centralityLow>=0.20&&centralityUpp<=0.40) {
			return TMath::Sqrt(411.0)/411.0;
		} else if (centralityLow>=0.40&&centralityUpp<=0.80) {
			return TMath::Sqrt(140.0)/140.0;
		} else {
			return -1.0;
		} 
	
}

///////////////////////////////
//Z counts from MC 
//////////////////////////////
double getZEventsMC(double centralityLow, double centralityUpp, double ncoll, double frac, RooDataSet* mcNEvents){

    // ---Eta/Centrality Slices --- //

    ///frac = alpha_{Z} in note
    std::cout << "Fraction of generated muons from Z decays in this eta slice = " << frac << std::endl;
	double ZEtaCentSlice = getZEtaCentSliceMC(centralityLow,centralityUpp,ncoll,frac,mcNEvents);
    //std::cout << "Fraction of Zs generated in this eta slice :" << frac << std::endl;
    std::cout << " Number of Z events in this eta slice from data: " << ZEtaCentSlice << std::endl;
	if(ZEtaCentSlice>0) return ZEtaCentSlice;
	else {
            std::cout << "WARNING: error in returning number of Z events in data. Expect wrong results." << std::endl;
            return -1.0;
        }

}///getZEvents



///////////////////////////////
//Z counts from data
//////////////////////////////
double getZEventsData(double centralityLow, double centralityUpp, double frac){

    // ---Eta/Centrality Slices --- //
    ///values calculated from Fig 1&4 of HI Z paper
    //frac = alpha_{Z} in analysis note

    //double frac = grZEta->GetY()[ieta];

    std::cout << "Fraction of generated muons from Z decays in this eta slice = " << frac << std::endl;
	double ZEtaCentSlice = getZEtaCentSliceData(centralityLow,centralityUpp,frac);
    //std::cout << "Fraction of Zs generated in this eta slice :" << frac << std::endl;
    std::cout << " Number of Z events in this eta slice from data: " << ZEtaCentSlice << std::endl;
	if(ZEtaCentSlice>0) return ZEtaCentSlice;
	else {
            std::cout << "WARNING: error in returning number of Z events in data. Expect wrong results." << std::endl;
            return -1.0;
        }

}///getZEvents

///////////////////////////////
//Zmumu fraction in signal region from MC study
//(called b_{Z} in analysis note
//////////////////////////////
double getZBkg(int iMt, int ieta, int icent, int index, TGraph2DErrors* grZIn, double centralityLow, double centralityUpp, bool isPreSel = false){

    //99 = mu^{pm}, 100 = mu-(inclusive), 101 = mu+(inclusive), 102 = mu+ (binned), 103 = mu- (binned), 104 = mu^{pm} (binned)
    if(ieta>9) std::cout << "WARNING: Eta bin has no affiliated background fraction. EXPECT WRONG RESULTS." <<
        std::endl;

	//return pre-selected efficiency
	//0-80%, mu^{pm}
	if(isPreSel){
		return 0.986346;
	}


	//one-legged mu+
	if(iMt==102){

		std::cout << "Z background fraction in bin for mu^{+}:" << iMt << ":" << ieta << ":" << icent << " = " << grZIn->GetZ()[index] << std::endl;
  		return grZIn->GetZ()[index];

		if(centralityLow>=0.0&&centralityUpp<=0.05) {
			return 0.0685714 ;
		} else if (centralityLow>=0.05&&centralityUpp<=0.1){ 
			return 0.059448;
		} else if (centralityLow>=0.1&&centralityUpp<=0.15){ 
			return 0.0408602;
		} else if (centralityLow>=0.15&&centralityUpp<=0.2){ 
			return 0.044843;
		} else if (centralityLow>=0.2&&centralityUpp<=0.4){ 
			return 0.0622969;
		} else if (centralityLow>=0.4&&centralityUpp<=0.8) {
			return 0.0479873;
		} else if (centralityLow>=0.0&&centralityUpp<=0.8) {
			return 0.0522884;
		} else {
			std::cout << "WARNING: UNABLE TO LOCATE CENTRALITY BIN." << std::endl;  	
			return -9999.;
		}
	} else if (iMt==103){

		std::cout << "Z background fraction in bin for mu^{-}:" << iMt << ":" << ieta << ":" << icent << " = " << grZIn->GetZ()[index] << std::endl;
		return grZIn->GetZ()[index];

		if(centralityLow>=0.0&&centralityUpp<=0.05) {
			return 0.0342857;
		} else if (centralityLow>=0.05&&centralityUpp<=0.1){ 
			return 0.0467091;
		} else if (centralityLow>=0.1&&centralityUpp<=0.15){ 
			return 0.0387097;
		} else if (centralityLow>=0.15&&centralityUpp<=0.2){ 
			return 0.058296;
		} else if (centralityLow>=0.2&&centralityUpp<=0.4){ 
			return 0.0390033;
		} else if (centralityLow>=0.4&&centralityUpp<=0.8) {
			return 0.0407096;
		} else if (centralityLow>=0.0&&centralityUpp<=0.8) {
			return 0.0412539;
		} else {
			std::cout << "WARNING: UNABLE TO LOCATE CENTRALITY BIN." << std::endl;  	
			return -9999.;
		}
	} else if (iMt==104){

		std::cout << "Z background fraction in bin for mu^{pm}:" << iMt << ":" << ieta << ":" << icent << " = " << grZIn->GetZ()[index] << std::endl;
		return grZIn->GetZ()[index];

        } else if (centralityLow>=0.0&&centralityUpp<=0.8) {
			return 0.0935423;
	} else {
			std::cout << "WARNING: UNABLE TO LOCATE CENTRALITY BIN." << std::endl;  	
			return -9999.;
	  }
}


///////////////////////////////
//statistical error in Zmumu fraction from MC
//////////////////////////////
double getZBkgStatErrorMC(int iMt, int ieta, int icent, int index, TGraph2DErrors* grZIn, double centralityLow, double centralityUpp, bool isPreSel = false){

	//one-legged mu+
	if(iMt==102){

		std::cout << "Z background fraction MC stat error in bin for mu^{+}:" << iMt << ":" << ieta << ":" << icent <<
        "= " << grZIn->GetEZ()[index] << std::endl;
  		return grZIn->GetEZ()[index];

	} else if (iMt==103){

		std::cout << "Z background fraction MC stat error in bin for mu^{-}:" << iMt << ":" << ieta << ":" << icent << " = " << grZIn->GetEZ()[index] << std::endl;
		return grZIn->GetEZ()[index];

	} else if (iMt==104){

		std::cout << "Z background fraction MC stat error in bin for mu^{#pm}:" << iMt << ":" << ieta << ":" << icent << " = " << grZIn->GetEZ()[index] << std::endl;
		return grZIn->GetEZ()[index];
    } 
 
    else {
			std::cout << "WARNING: UNABLE TO LOCATE CENTRALITY BIN." << std::endl;  	
			return -9999.;
	  }
}


///////////////////////////////
//Zmumu fraction in signal region 
//binned over centrality and
//averaged over eta
//////////////////////////////
double getZBkgCent(int iMt, int ieta, int icent, int index, TGraphErrors* grZIn, double centralityLow, double centralityUpp, bool isPreSel = false){

  		return grZIn->GetY()[icent];
}

///////////////////////////////
//QCD fraction in signal region from data
//////////////////////////////
double getQCDBkg(int iMt, int ieta, int icent, TGraphErrors* grQCDIn, double centralityLow, double centralityUpp, bool isPreSel = false){

  //99 = mu^{pm}, 100 = mu-(inclusive), 101 = mu+(inclusive), 102 = mu+ (binned), 103 = mu- (binned), 104 = mu^{pm} (binned)

  if(isPreSel){
		if(centralityLow>=0.0&&centralityUpp<=0.10) {
			return 0.176637;
		} else if(centralityLow>=0.10&&centralityUpp<=0.20) {
			return 0.204476;
		} else if (centralityLow>=0.2&&centralityUpp<=0.4){ 
			return 0.243826;
		} else if (centralityLow>=0.4&&centralityUpp<=0.8) {
			return 0.315636; 
		} else if (centralityLow>=0.0&&centralityUpp<=0.8) {
			return 0.211384;
			//return 0.0852526; //bkg w/o mt cut
		} else {
			std::cout << "WARNING: UNABLE TO LOCATE CENTRALITY BIN." << std::endl;  	
			return -9999.;
		}
  }  


  if(iMt==99) {
	if (centralityLow>=0.0&&centralityUpp<=0.8) {
		return 0.0328355; 
	} else {
		std::cout << "WARNING: UNABLE TO LOCATE CENTRALITY BIN." << std::endl;  	
		return -9999.;
	}

  } else if(iMt == 100){
	if (centralityLow>=0.0&&centralityUpp<=0.8) {
		return 0.0357009; 
	} else {
		std::cout << "WARNING: UNABLE TO LOCATE CENTRALITY BIN." << std::endl;  	
		return -9999.;
	}
  } else if(iMt == 101){
	if (centralityLow>=0.0&&centralityUpp<=0.8) {
		return 0.0300548; 
	} else {
		std::cout << "WARNING: UNABLE TO LOCATE CENTRALITY BIN." << std::endl;  	
		return -9999.;
	}
  }  else if(iMt==102){

	std::cout << "QCD background fraction in bin for mu^{+}:" << iMt << ":" << ieta << ":" << icent << " = " << grQCDIn->GetY()[ieta] << std::endl;
  	return grQCDIn->GetY()[ieta];

	if(centralityLow>=0.0&&centralityUpp<=0.10) {
		return 0.0198221;
	} else if(centralityLow>=0.10&&centralityUpp<=0.20) {
		return 0.0124922;
	} else if (centralityLow>=0.2&&centralityUpp<=0.4){ 
		return 0.036892;
	} else if (centralityLow>=0.4&&centralityUpp<=0.8) {
		return 0.0423467; 
	} else if (centralityLow>=0.0&&centralityUpp<=0.8) {
		return 0.0300548; 
	} else {
		std::cout << "WARNING: UNABLE TO LOCATE CENTRALITY BIN." << std::endl;  	
		return -9999.;
	}
  } else if (iMt==103){ //mu^{-}

	std::cout << "QCD background fraction in bin for mu^{-}:" << iMt << ":" << ieta << ":" << icent << " = " << grQCDIn->GetY()[ieta] << std::endl;
  	return grQCDIn->GetY()[ieta];

	if(centralityLow>=0.0&&centralityUpp<=0.10) {
		return 0.0171457;
	} else if(centralityLow>=0.10&&centralityUpp<=0.20) {
		return 0.0203586;
	} else if (centralityLow>=0.2&&centralityUpp<=0.4){ 
		return 0.0371438;
	} else if (centralityLow>=0.4&&centralityUpp<=0.8) {
		return 0.0574238; 
	} else if (centralityLow>=0.0&&centralityUpp<=0.8) {
		return 0.0357009; 
	}  else {
		std::cout << "WARNING: UNABLE TO LOCATE CENTRALITY BIN." << std::endl;  	
		return -9999.;
	}
  } else if(iMt==104){

	std::cout << "QCD background fraction in bin for mu^{pm}:" << iMt << ":" << ieta << ":" << icent << " = " << grQCDIn->GetY()[ieta] << std::endl;
    //if(xEta<0.0) ieta = nEtaBins-ieta;
  	return grQCDIn->GetY()[ieta];
  } else {
		std::cout << "WARNING: UNABLE TO LOCATE CENTRALITY BIN." << std::endl;  	
		return -9999.;
	}

}

///////////////////////////////
//QCD fraction in signal region averaged over eta
//as a function of centrality
//////////////////////////////
double getQCDBkgCent(int iMt, int icent, TGraphErrors* grQCDIn, double centralityLow, double centralityUpp, bool isPreSel = false){

	//hack for using 0-10% instead of 0-5%, 5-10%
    //if(icent!=0) --icent; 
    return grQCDIn->GetY()[icent];
}

////////////////////////
//Prevent trailing digits
//in figure labeling
////////////////////////
TString format(float value) {
  std::stringstream svalue;
  svalue  << std::setprecision(2) << value;
  return svalue.str();
}

////////////////////////
//Prevent trailing digits
//in figure labeling
////////////////////////
TString format(float value, int precision) {
  std::stringstream svalue;
  svalue  << std::setprecision(precision) << value;
  return svalue.str();
}
//////////////////////////////////////////
//Plot the last element in the THStack
//produced from the individual histograms
//from each bin for MC and Data
///////////////////////////////////////////

void plotMergedKinematic(THStack* hmData,THStack* hmQCD, THStack* hmZ,THStack* hmTau, THStack* hmMcSum,
                const float kinVarLow, const float kinVarMax,int nBins, float xLo, TString sXtitle, 
				TString sSel, TString sSel2, int chargeIndex, bool correctSpectra, bool scaleBW=false){

	TCanvas* cKin = new TCanvas("cKin","cKin",600,600);
	//return summed histos in the THStack

    ///Divide canvas for plotting ratios
    std::cout << std::endl;
    std::cout << "TPad..." << std::endl;
    cKin->Divide(1, 2);
    TPad* canvas_up = (TPad*)cKin->GetListOfPrimitives()->FindObject("cKin_1");
    TPad* canvas_dw = (TPad*)cKin->GetListOfPrimitives()->FindObject("cKin_2");
    std::cout << "Done." << std::endl;

    ///Define the size
    double up_height     = 0.8; 
    double dw_correction = 1.60;
    double font_size_dw  = 0.08;
    double dw_height    = (1. - up_height) * dw_correction;

    ///set pad size
    canvas_up->SetPad(0., 1 - up_height, 1., 1.);
    canvas_up->SetFillColor(0);
    canvas_up->SetBorderMode(0);
    canvas_up->SetBorderSize(2);
    canvas_up->SetTickx(1);
    canvas_up->SetTicky(1);
    canvas_up->SetLeftMargin(0.16);
    canvas_up->SetRightMargin(0.05);
    canvas_up->SetTopMargin(0.05);
    canvas_up->SetBottomMargin(0.16);
    canvas_up->SetFrameBorderMode(0);
    canvas_up->SetFrameBorderMode(0);
 
    canvas_dw->SetPad(0., 0., 1., dw_height);
    canvas_dw->Range(-0.3639066,-0.7754386,2.546497,1.31186);
    canvas_dw->SetFillColor(0);
    canvas_dw->SetBorderMode(0);
    canvas_dw->SetBorderSize(2);
    canvas_dw->SetTickx(1);
    canvas_dw->SetTicky(1);
    canvas_dw->SetLeftMargin(0.159396);
    canvas_dw->SetRightMargin(0.05033557);
    canvas_dw->SetTopMargin(0.005681818);
    canvas_dw->SetBottomMargin(0.3715035);
    canvas_dw->SetFrameBorderMode(0);
    canvas_dw->SetFrameBorderMode(0);

    canvas_up->SetFrameFillColor(0);
    canvas_up->SetFillColor(0);
    canvas_dw->SetFillColor(0);
    canvas_dw->SetFrameFillColor(0);

	TH1F* hmDatac =(TH1F*)hmData->GetStack()->Last()->Clone("hmDatac");
	TH1F* hmQCDc =(TH1F*)hmQCD->GetStack()->Last()->Clone("hmQCDc");
	TH1F* hmZc =(TH1F*)hmZ->GetStack()->Last()->Clone("hmZc");
	TH1F* hmTauc =(TH1F*)hmTau->GetStack()->Last()->Clone("hmTauc");
	TH1F* hmMcSumc =(TH1F*)hmMcSum->GetStack()->Last()->Clone("hmMcSumc");

	///int binLo = nBins/kinVarMax*kinVarLow+1;
    ///get the bin corresponding to kinVarLow

    ///note:for uniform binning only
    int binLo = getBinNumber(kinVarLow,nBins,xLo,kinVarMax);


	std::cout << "integrating from bin " << binLo << " to " << nBins << std::endl;
  	double sigEventsUncorr = hmDatac->Integral(binLo,nBins); //integrate from 40-200 GeV before Aw,Cw correction
  	std::cout << "integrated events before correction in range:"<< kinVarLow << "-" << kinVarMax  << " =  " << sigEventsUncorr << " +-" << TMath::Sqrt(sigEventsUncorr) << std::endl;

	double effMt ;
	correctSpectra = false; //temp hack
	if(correctSpectra) {
		std::cout << "Correcting for acceptance..." << std::endl;
		hmDatac = correctAcceptanceMt(hmDatac, chargeIndex);
		hmQCDc = correctAcceptanceMt(hmQCDc, chargeIndex,false);
		hmZc = correctAcceptanceMt(hmZc, chargeIndex,false);
		hmMcSumc = correctAcceptanceMt(hmMcSumc, chargeIndex,false);
		effMt = getAcceptanceMt(chargeIndex);
	}
	else effMt = 1.0;
	std::cout << "Efficiency correction factor = " << effMt << std::endl;
 

	hmQCDc->SetFillColor(kAzure-9);
	hmZc->SetFillColor(kRed);
	hmTauc->SetFillColor(kYellow-7);


	std::cout << "Parameter mean = " << hmDatac->GetMean() << "+- " << hmDatac->GetMeanError() << " RMS = " << hmDatac->GetRMS() << std::endl;
	std::cout << "|Data-MC|/Data*100 = " << fabs(hmDatac->GetMean()-hmMcSumc->GetMean())/hmDatac->GetMean()*100 << "%" << std::endl;

    ///MC/Data histo
    TH1F* hRatio = (TH1F*)hmDatac->Clone("hRatio");
    hRatio = getDataMCRatio(hRatio,hmDatac,hmMcSumc);

    ///Draw spectra
    canvas_up->cd();



	hmMcSumc->GetXaxis()->SetTitle(sXtitle); 
    float entriesPerBin = (float)(kinVarMax-xLo)/nBins;
	TString sY = "Muons/"; TString sentriesPerBin = format(entriesPerBin); sY+=sentriesPerBin;

    //divide by binwidth
    if(scaleBW){
        hmDatac->Scale(1.0,"width");
        hmQCDc->Scale(1.0,"width");
        hmZc->Scale(1.0,"width");
        hmMcSumc->Scale(1.0,"width");
    }


	//hmMcSumc->GetYaxis()->SetTitle(sY); 
	hmMcSumc->GetYaxis()->SetTitle("Muons"); 
	//hmMcSumc->GetYaxis()->SetTitle("dN^{#mu}/d #eta"); 
	//hmMcSumc->GetYaxis()->SetRangeUser(0.1,1600.0); 
	hmMcSumc->GetXaxis()->SetRangeUser(kinVarLow,kinVarMax); 
	hmMcSumc->Draw("hist f");
	hmTauc->Draw("hist fsame");
	hmZc->Draw("hist fsame");
	hmQCDc->Draw("hist fsame");
	hmDatac->Draw("pesame");
	hmMcSumc->Draw("sameaxis");

    //canvas_up->SetLogy(true);

	ATLAS_LABEL(0.17,0.89,1);
	myText(0.33,0.89, (Color_t)kBlack, (char*)("Internal"));

	TLegend* leg = new TLegend(0.658, 0.637, 0.928, 0.867);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->AddEntry(hmDatac, "Data 2011", "pe");
	leg->AddEntry(hmMcSumc, "W#rightarrow#mu#nu", "f");
	leg->AddEntry(hmTauc, "W#rightarrow#tau#nu", "f");
	leg->AddEntry(hmZc, "Z#rightarrow#mu#mu", "f");
	leg->AddEntry(hmQCDc, "QCD", "f");
	leg->Draw();

	TLatex l;
	l.SetNDC();
	l.DrawLatex(0.18,0.83,sSel +  "%" );
	//l.DrawLatex(0.169,0.767,"-2.4<#eta<2.4");
	l.DrawLatex(0.169,0.767,"0.1<|#eta|<2.4");
	l.SetTextSize(0.034);
	l.DrawLatex(0.74,0.89,"#sqrt{s_{NN}}=2.76 TeV");
	l.DrawLatex(0.492,0.89,"#int Ldt #approx 0.14 nb^{-1}"); 
	cKin->Update();
	
	TString plotNameLog = "data,"; plotNameLog+=sSel; plotNameLog+=","; /*plotNameLog+=sSel2; plotNameLog+="_";*/  plotNameLog+=sXtitle; plotNameLog+=",Log";
	TString plotNameLin = "data,"; plotNameLin+=sSel; plotNameLin+=","; /*plotNameLin+=sSel2; plotNameLin+="_";*/ plotNameLin+=sXtitle; plotNameLin+=",Lin"; 
    
    hmMcSumc->GetYaxis()->SetRangeUser(0.1,hmMcSumc->GetMaximum()+7.7e2); cKin->Update();
//    hmMcSumc->GetYaxis()->SetRangeUser(0.009,375.0); cKin->Update(); ///use for phi
//    hmMcSumc->GetYaxis()->SetRangeUser(0.009,500.0); cKin->Update(); ///use for eta
//    hmMcSumc->GetYaxis()->SetRangeUser(0.009,750.0); cKin->Update(); //use for |eta|
    //hmMcSumc->GetYaxis()->SetRangeUser(0.009,2.0e4); cKin->Update();
    //hmMcSumc->GetYaxis()->SetRangeUser(0.1,hmMcSumc->GetMaximum()+4e2); cKin->Update();

    ///Draw MC/Data underneath
    canvas_dw->cd();
    /// font size
    hRatio->GetXaxis()->SetTitle(sXtitle);
    hRatio->GetXaxis()->SetLabelFont(42);
    hRatio->GetXaxis()->SetLabelSize(0.09);
    hRatio->GetXaxis()->SetTitleSize(0.11);
    hRatio->GetXaxis()->SetTitleOffset(0.8);
    hRatio->GetXaxis()->SetTitleFont(42);
    hRatio->GetYaxis()->SetTitle("Data/MC");
    hRatio->GetYaxis()->SetLabelFont(42);
    hRatio->GetYaxis()->SetLabelSize(0.09);
    hRatio->GetYaxis()->SetTitleSize(0.11);
    hRatio->GetYaxis()->SetTitleOffset(0.7);

    hRatio->GetYaxis()->SetRangeUser(0.4,2.0);

    /*hRatio->GetXaxis()->SetLabelSize(font_size_dw);
    hRatio->GetXaxis()->SetTitleSize(font_size_dw);
    hRatio->GetYaxis()->SetLabelSize(font_size_dw);
    hRatio->GetYaxis()->SetTitleSize(font_size_dw);
    hRatio->GetYaxis()->SetTitleOffset(0.7);
    hRatio->GetXaxis()->SetTitleOffset(1.4);
    */
    
    hRatio->Draw();


	
	plotNameLin.ReplaceAll("#","");
	plotNameLin.ReplaceAll("[","");
	plotNameLin.ReplaceAll("]","");
	cKin->Print(plotNameLin.ReplaceAll("|",",") +".png");
	cKin->Print(plotNameLin.ReplaceAll("|",",") +".eps"); 
	cKin->Print(plotNameLin.ReplaceAll("|",",")+ ".pdf"); 
	TString plotNameLinRoot = plotNameLin.ReplaceAll("|",",") + ".root";
	cKin->Print(plotNameLin); 


	plotNameLog.ReplaceAll("#","");
	plotNameLog.ReplaceAll("[","");
	plotNameLog.ReplaceAll("]","");

    canvas_up->cd();
    
    hmMcSumc->GetYaxis()->SetRangeUser(0.11,1.6e4); cKin->Update();
    //canvas_up->SetLogy(true); /*hmMcSumc->GetYaxis()->SetRangeUser(0.1, hmMcSumc->GetMaximum()*1e3);*/ 
    canvas_up->Update();

    canvas_dw->cd();
    hRatio->GetXaxis()->SetRangeUser(kinVarLow,kinVarMax);
    hRatio->Draw();
	//hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()*1e3); cKin->Update();
	cKin->Print(plotNameLog.ReplaceAll("|",",") +".png");
	cKin->Print(plotNameLog.ReplaceAll("|",",") +".eps"); 
	cKin->Print(plotNameLog.ReplaceAll("|",",") +".pdf"); 
	TString plotNameLogRoot = plotNameLog.ReplaceAll("|",",") + ".root";
	cKin->Print(plotNameLogRoot); 
	

} ///merged histograms

/////////////////////////////////////////////
//Generic function for plotting
//kinematic variables for the W. This 
//is obsolete, and it is recommended to use the
//more specialized functions below.
//////////////////////////////////////////////
void plotWCandidateKinematic( RooDataSet* dataSet, RooDataSet* mcWSet, RooDataSet* mcTauSet, RooDataSet* mcZSet, 
                RooDataSet* mcJ1Set, RooDataSet* mcJ2Set, RooDataSet* mcJ3Set, 
				RooDataSet* mcJ1Events,RooDataSet* mcJ2Events, RooDataSet* mcJ3Events, 
                TFile* fQCDBkg, TGraphAsymmErrors* grBkgTau,TGraph2DErrors* grBkgZ, double totalZinEta, 
                double Cw, double cWSystErr, double qcdSystError, double zBosSystError, 
				THStack* hmData, THStack* hmQCD,THStack* hmZ, THStack* hmTau, THStack* hmMcSum, 
                const int iMt, const int iEta, const int
                iCentrality, int nCentralityBins, int nEtaBins, float xLo,
				float kinVarLow, float kinVarMax ,double ncoll, double ptLow , double ptUpp, double etaLow, double etaUpp, double centralityLow, 
                double centralityUpp, RooRealVar& kinVar, TString sSel, TString sSel2 , int nBins, TString sXtitle, 
				bool addPercent = true, bool correctSpectra = true, bool doSubtractBkg = true, bool
                doPreSelKinematics=false, bool doMirrorEta=true) {


	std::cout << "Plotting for kinematic variable " << sXtitle << " in range " << kinVarLow << " to " << kinVarMax <<std::endl;

  	RooBinning b = RooBinning(nBins,xLo,kinVarMax); 

	//initialize histograms	
  	// --- data ---
	TH1F* hdataSet = (TH1F*)dataSet->createHistogram("hdataSet",kinVar ,Binning(b));
  	// --- W set ---
	TH1F* hmcWSet = (TH1F*)mcWSet->createHistogram("hmcWSet",kinVar,Binning(b));
  	// --- Z set ---
	TH1F* hmcZSet = (TH1F*)mcZSet->createHistogram("hmcZSet",kinVar,Binning(b));
    // --- Wtau set ---
	TH1F* hmcTauSet = (TH1F*)mcTauSet->createHistogram("hmcTauSet",kinVar,Binning(b));
  	// --- QCD set ---
	TH1F* hmcJ1Set = (TH1F*)mcJ1Set->createHistogram("hmcJ1Set",kinVar,Binning(b));
	TH1F* hmcJ2Set = (TH1F*)mcJ2Set->createHistogram("hmcJ2Set",kinVar,Binning(b));
	TH1F* hmcJ3Set = (TH1F*)mcJ3Set->createHistogram("hmcJ3Set",kinVar,Binning(b));

	///Return correctly weighted QCD histogram
    std::cout << "Weighting QCD eta distributions with Jxmu samples..." << std::endl;

	TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",nBins, xLo, kinVarMax);
	hmcQCDSet = getWeightedJxHisto(centralityLow, centralityUpp, ncoll, hmcJ1Set, hmcJ2Set, hmcJ3Set, mcJ1Events, mcJ2Events,mcJ3Events, nBins,xLo, kinVarMax);
    std::cout << "Bins must be equal. Let's check: " << hmcQCDSet->GetNbinsX() << " =? " << hmcJ1Set->GetNbinsX() << " =? " << hmcJ2Set->GetNbinsX() 
        << " =? " << hmcJ3Set->GetNbinsX() << std::endl;
    std::cout << "Done." << std::endl;

    ///return the lower bin number
    ///which you would like to integrate
    ///from in order to obtain signal yield
    ///Note:For uniform binning only
    int binLo(0); 
	if(doPreSelKinematics) binLo = 1;
	else binLo = getBinNumber(kinVarLow,nBins,xLo,kinVarMax);
    ///For variable binning, declare a RooBinning object
    ///and call RooBinning::binNumber(double x) method
    //int binLo = b.binNumber(kinVarLow)+1;

    std::cout << "Signal bins = " << binLo << "-" << nBins << std::endl;
	std::cout << "Corresponding to eta range " << kinVarLow << " to " << kinVarMax << std::endl;

	TH1F* hdataSetc = (TH1F*)hdataSet->Clone("hdataSetc");

    ///Number of observed events in signal region
    ///before efficiency corrections
  	double sigEventsUncorr = hdataSetc->Integral(binLo,nBins); 
  	std::cout << "Integrated events before correction in eta:"<< kinVarLow << "-" << kinVarMax  << " =  " 
        << sigEventsUncorr << " +-" << TMath::Sqrt(sigEventsUncorr) << std::endl;
	
	double xEta = etaLow+(etaUpp-etaLow)/2.0;
	//relative errors of uncorrected signal events and correction factor
	double err1 = TMath::Sqrt(sigEventsUncorr); 
    ///statistical errors of the raw number of counts
    double errSigStat = err1;
	
//////////
//QCD
/////////

    TString sFracQCD = "fractionQCD"; sFracQCD+="_charge"; sFracQCD+=104; sFracQCD+="_eta"; sFracQCD+=0; sFracQCD+="_cent"; 
	//hack for using 0-10% instead of 0-5%, 5-10%
	/*if(iCentrality==0)*/ sFracQCD+=iCentrality;
	//else sFracQCD+=iCentrality-1;
	TGraphErrors* grBkgQCD = 0;
    if(doSubtractBkg) grBkgQCD = (TGraphErrors*) fQCDBkg->Get(sFracQCD);

	TH1F* hmcQCDSetc = (TH1F*)hmcQCDSet->Clone("hmcQCDSetc");
	
	double survivalFractionQCD ; 

	std::cout << "Getting QCD background fraction..." << std::endl;
	///use inclusive charge set (i.e. iMt=104)
	if(doPreSelKinematics) 

	    ///percent of QCD in W signal region b4 W selection
        survivalFractionQCD = getQCDBkg(104,getQCDBkgBinNumber(iEta,xEta,nEtaBins,doMirrorEta),iCentrality,
                                            grBkgQCD,centralityLow, centralityUpp, doPreSelKinematics) ; 	

	else if(doSubtractBkg) survivalFractionQCD = getQCDBkg(104,getQCDBkgBinNumber(iEta,xEta,nEtaBins,doMirrorEta), 
                                                            iCentrality,grBkgQCD,centralityLow, centralityUpp) ;
    else survivalFractionQCD = 0.0;

	double dataQCDEvents = survivalFractionQCD*sigEventsUncorr;

	//relative error of N_corr(sigEvents) 
    ///error in the difference of the QCD bkg frac for mu+ and mu-; relative error shown below
    //double errQCDFrac = qcdSystError/survivalFractionQCD;

	///Error of dataQCDEvents
	double errQCDCountsSyst = 0.0;
	double errQCDCountsStat = 0.0;
	double totalPercentStatError = 0.0; 
    double totalPercentSystError = 0.0;

/*    if(doSubtractBkg) {
        errQCDCountsSyst = errQCDFrac*dataQCDEvents;
        errQCDCountsStat = TMath::Sqrt(dataQCDEvents);

        std::cout << "Percent systematic uncertainty of the data from QCD = " << errQCDCountsSyst/sigEventsUncorr*100 << "%" << std::endl;
        totalPercentSystError = errQCDCountsSyst/sigEventsUncorr*100.0;
    }
*/
    std::cout << "Current systematic error on the number of events in the data after QCD subtraction = +/- " 
            << totalPercentSystError << "%" << std::endl;

    ///Integrate QCD spectrum in signal region
  	double mcQCDEvents = hmcQCDSetc->Integral(binLo,nBins);
	std::cout << "QCD integral in signal region = " <<  mcQCDEvents << std::endl;

    ///Scale factor for scaling to number of
    ///expected QCD events in the data
	double sfQCD = dataQCDEvents/mcQCDEvents;

	if(mcQCDEvents==0) {
		std::cout << "WARNING: 0 QCD MC events in signal region." << std::endl;
		sfQCD = 1.0;
	}
      
    ///Normalize QCD MC shape to
    ///number of expected QCD events in
    ///the data
	hmcQCDSetc->Scale(sfQCD);

    ///Keep a running sum of
    ///QCD bkg histogram for
    ///final figure
	hmQCD->Add(hmcQCDSetc);

	hmcQCDSetc->SetFillColor(kAzure-9);
	std::cout << "Expected number of QCD events in the signal region of the data = " << dataQCDEvents << std::endl;
  	std::cout << "Background percentage from QCD in kinematic region:"<< kinVarLow << "-" << kinVarMax  << " =  " <<
        dataQCDEvents/sigEventsUncorr*100.0 << "%" <<std::endl;

///////////
///Z boson
//////////

	std::cout << "Z integral " << hmcZSet->Integral() << std::endl;
	TH1F* hmcZSetc = (TH1F*)hmcZSet->Clone("hmcZSetc");

	double survivalProb; 
	int indexBin = -9999;

    ///return the probability of generated Z 
    ///to survivie W selection cuts
	if(doPreSelKinematics) {
	    indexBin = iEta*nCentralityBins+iCentrality;
        survivalProb = getZBkg(iMt,iEta,iCentrality,indexBin,grBkgZ,centralityLow, centralityUpp,doPreSelKinematics) ;
    }
	else if(doSubtractBkg) {

        int indexEta;

        ///bin shift necessary if using
        ///real eta (i.e. not absolute eta)
        if(!doMirrorEta){
          ///hop over the crack region
          ///in bin [-0.1,0.1]
          if(xEta==0.0) std::cout << "WARNING: Bin in crack region. Expect wrong result!" << std::endl;
          if(xEta<0.0){
              ///number of bins in absolute eta
              indexEta = indexNegativeEta(iEta,nEtaBins);
          }
           else {
              indexEta = indexPositiveEta(iEta,nEtaBins);
           }

        }
 
        ///no bin shift if doing |eta|
        else {
            indexEta = iEta;
        }

	    indexBin = indexEta*nCentralityBins+iCentrality;
        survivalProb = getZBkg(iMt,indexEta,iCentrality,indexBin, grBkgZ, centralityLow, centralityUpp) ; //percent of Z in W signal region
    }
    else survivalProb = 0.0;

	std::cout <<"Z survival probability = " << survivalProb  << std::endl;

    ///Apply this fraction to the number of Zs in this centrality and eta bin 
    ///This comes from the published Z note
	double bkgZEventFrac = survivalProb*getZEventsData(centralityLow, centralityUpp, totalZinEta)/sigEventsUncorr;
	std::cout << "Fraction of Zs surviving W selection cuts  = " << bkgZEventFrac << std::endl;

    ///Expected number of one-legged muons in the data
	double dataZEvents = bkgZEventFrac*sigEventsUncorr ;
    std::cout << "That leaves an estimated " << dataZEvents << " one legged-muons in the data signal region." << std::endl;

	///Absolute systematic in the number of Z events calculated from 
    ///the difference in the bkg frac for mu+ and mu- one-legged muons.
    ///Propagate into data stat. and QCD systematic errors
	double errZCountsStat = 0.0;
	double errZCountsSyst = 0.0;

    if(doSubtractBkg) {

        ///stat. error
        errZCountsStat = TMath::Sqrt(dataZEvents);

        ///Propagate syst. error 
        errZCountsSyst = zBosSystError/bkgZEventFrac*dataZEvents;  
        std::cout << "Percent systematic uncertainty from Z = " << errZCountsSyst/sigEventsUncorr*100 << "%" << std::endl;

        totalPercentSystError = TMath::Sqrt( TMath::Power(totalPercentSystError,2) +
            TMath::Power(errZCountsSyst/sigEventsUncorr*100,2));
    }

    std::cout << "Current systematic error on the number of events in the data after QCD and Z subtraction = +/- " 
            << totalPercentSystError << "%" << std::endl;

    ///Number of Z events in the signal region from MC
  	double mcZEvents = hmcZSetc->Integral(binLo,nBins);

    ///Scale factor for scaling MC events
    ///to the expected number of Z background muons
    ///in the data
	double sfZ = dataZEvents/mcZEvents;

	if(mcZEvents==0) {
		std::cout << "WARNING: 0 Z MC events in signal region." << std::endl;
		sfZ = 1.0;
	}

    ///Scale the histogram to the expected number of Z background muons
    ///in the data
	hmcZSetc->Scale(sfZ);

    ///Add the QCD portion to afford a "total" background
    ///histogram.
	hmcZSetc->Add(hmcQCDSetc);
	hmcZSetc->SetFillColor(kRed);

	///Keep a running sum of Z+QCD histograms for
    ///the final figure integrated over all eta and centrality
	hmZ->Add(hmcZSetc);

  	std::cout << "Background Z percentage expected in the data from :"<< kinVarLow << "-" << kinVarMax  << " =  " <<
        dataZEvents/sigEventsUncorr*100.0 << "%" << std::endl;

////////////////////////
//Background subtraction
////////////////////////

    ///Histogram with QCD+Z component subtracted off
	TH1F* hdataSetSubtracted = (TH1F*)hdataSetc->Clone("hdataSetSubtracted");
    float dataTauEvents=0.0;

    ///Subtract the QCD+Z histogram from the data
    ///Expected number of reconstructed signal events in data
    ///after background subtraction.
	double sigEventsSub = sigEventsUncorr-dataQCDEvents-dataZEvents;
    std::cout << "Expected number of reconstructed W events in the data after background subtraction = " 
        << sigEventsSub << std::endl;

    ///Fetch %tau contamination per W->mu event in this eta/centralaty bin
    float tauBkgFraction = grBkgTau->GetY()[iEta];
    ///Number of tau events in the data
    dataTauEvents = tauBkgFraction/(1.0+tauBkgFraction)*sigEventsSub;
    std::cout << "Number of background events from W-->tau-->mu decays: " << dataTauEvents << std::endl;
    ///Subtract off tau background
    sigEventsSub-=dataTauEvents;
    std::cout << "Expected number of reconstructed W events in the data after QCD,Z, and tau background subtraction = "
        << sigEventsSub << std::endl;

    ///Propagate stat errror from QCD and Z into data counts for current stat error on the subtracted yield
    totalPercentStatError = TMath::Sqrt(TMath::Power(errSigStat,2) +
        TMath::Power(errQCDCountsStat,2) + TMath::Power(errZCountsStat,2) )/sigEventsSub*100.0;

    ///Absolute systematic errors due to the isolation cut.
    ///This is used to calculate absolute error in the subtracted yield
    ///for this charge,eta,and centrality class
    double isoPercentSystErr = getIsoSystErr(iMt, indexBin);
    std::cout << "Systematic uncertainty of number of signal events in data due to isolation efficiency = " << isoPercentSystErr << "%" << std::endl;

    ///Relative systematic errors due to the MPT cut.
    ///This is used to calculate absolute error in the subtracted yield
    ///for this charge,eta,and centrality class
    double mptPercentSystErr = getMptSystErr(iMt, indexBin) ;
    std::cout << "Systematic uncertainty of number of signal events in data due to MPT resolution = " << mptPercentSystErr << "%" << std::endl;

    ///Propagate isolation and MPT cut errors in final relative systematic 
    //error on the number of W events in the data.
    totalPercentSystError = TMath::Sqrt( TMath::Power(totalPercentSystError,2) + TMath::Power(isoPercentSystErr,2) + 
        TMath::Power(mptPercentSystErr,2) );

    std::cout << "Current statistical error on the number of events in the data after all background subtraction = +/- " 
            << totalPercentStatError << "%" << std::endl;
    std::cout << "Current systematic error on the number of events in the data after QCD,Z,isolation, and MPT errors propagated = +/- " 
            << totalPercentSystError << "%" << std::endl;

/////////////////////////
//W-->tau-->mu background
/////////////////////////
    std::cout << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "W --> T a u  B a c k g r o u n d" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << std::endl;

	std::cout << "W-->tau integral " << hmcTauSet->Integral() << std::endl;
	TH1F* hmcTauSetc = (TH1F*)hmcTauSet->Clone("hmcTauSetc");
    ///Number of tau->mu in signal region from MC
    double mcTauEvents = hmcTauSetc->Integral(binLo,nBins);
    ///Scale MC events to expected
    double sfTau = dataTauEvents/mcTauEvents;

	if(mcTauEvents==0) {
		std::cout << "WARNING: 0 Tau MC events in signal region." << std::endl;
		sfTau = 1.0;
	}

    ///Scale the histogram to the expected number of Tau background muons
    ///in the data
	hmcTauSetc->Scale(sfTau);

    ///Add the QCD+Z portion to afford a "total" background
    ///histogram.
	hmcTauSetc->Add(hmcZSetc);
	hmcTauSetc->SetFillColor(kYellow-7);

    if(doSubtractBkg){
      hdataSetSubtracted->Add(hmcTauSetc,-1);
      //Consistency check: 
      std::cout << "Consistency Check: " << std::endl;
      std::cout << sigEventsSub << " =? " << hdataSetSubtracted->Integral(binLo,nBins) << std::endl;
    }
	///Keep a running sum of Z+QCD+Tau histograms for
    ///the final figure integrated over all eta and centrality
	hmTau->Add(hmcTauSetc);

  	std::cout << "Background tau percentage expected in the data :" <<
        dataTauEvents/sigEventsUncorr*100.0 << "%" << std::endl;


///////////////
//Cw correction
//////////////

	double effCorrection ; double effCorrectionSystErr; 

	if(correctSpectra) {

        ///Correct background subtracted histogram with Cw
		//hdataSetc = correctEfficiencyMt(hdataSetSubtracted, 104, iCentrality, iEta,Cw);
        hdataSetSubtracted->Scale(1.0/Cw);
		effCorrection = Cw;

        ///Take systematic error on Cw
        ///as statistical error from W MC
        ///in eta/centrality bin
	    effCorrectionSystErr = cWSystErr;

	    std::cout << "Efficiency correction factor = " << effCorrection << " +/- " << effCorrectionSystErr << std::endl;
	}

	else {effCorrection = 1.0; effCorrectionSystErr = 0.0;}

    ///Number of background subtracted and corrected events in the data in the signal region
  	double sigEvents = hdataSetSubtracted->Integral(binLo,nBins); 
  	//double sigEvents = sigEventsSub/Cw; 

    double cWPercentSystErr = effCorrectionSystErr/effCorrection*100.0;

    ///Add in systematic uncertainty from Cw statistical error
    totalPercentSystError = TMath::Sqrt( TMath::Power(totalPercentSystError,2) + TMath::Power(cWPercentSystErr,2) );

    std::cout << "Current statistical error on the number of events in the data after all background subtraction = +/- " 
            << totalPercentStatError << "%" << std::endl;
    std::cout << "Current systematic error on the number of events in the data after QCD,Z,isolation, MPT, and Cw errors propagated = +/- " 
            << totalPercentSystError << "%" << std::endl;

	///Running sum of data histograms from 
    ///each eta and centrality bin 
	//hmData->Add(hdataSetSubtracted);
	hmData->Add(hdataSetc);
	TH1F* hmDataTemp =(TH1F*)hmData->GetStack()->Last()->Clone("hmDatac");
    std::cout << "Number of total entries in data after bin " << iEta << ":" << iCentrality << " = " << hmDataTemp->Integral() << std::endl;

//////////
//W 
/////////

    ///We now have a histogram with subtracted data corrected. 
    ///We also have QCD and Z histograms in the same canvas (normalized to the number of 
    ///expected background events). Below, we scale the W MC to the data and see how 
    ///well they agree.

	TH1F* hmcWSetc = (TH1F*)hmcWSet->Clone("hmcWSetc");

    ///Number of W events in MC in signal region
  	double mcWEvents =  hmcWSetc->Integral(binLo,nBins);


    ///Scale factor for scaling W MC events to 
    ///the number of expected signal events in the data.
	double sfW = sigEvents/mcWEvents;

	hmcWSetc->Scale(sfW);

    ///Number of expected W events from the MC
    double mcWExpected = hmcWSetc->Integral(binLo,nBins);
    std::cout << "Number of expected W bosons from the Monte Carlo = " << mcWExpected << std::endl;

    ///Add Z,Tau, and QCD histograms to obtain a 
    ///final MC histogram composed of 
    ///the expected number of Background+Signal events
	hmcWSetc->Add(hmcTauSetc);

	//running sum of combined(QCD+Z+W) MC histos	
    ///Keep a running sum of histograms for final figure
    ///integrated over all eta and centrality
	hmMcSum->Add(hmcWSetc);	

    std::cout << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "S u m m a r y" << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "Expected number of events for charge:eta:centrality" << iMt << ":" << iEta << ":" << iCentrality << " = " << mcWExpected << std::endl;
    std::cout << "Number of reconstructed events before background subtraction                                          = " << sigEventsUncorr << std::endl; 
    std::cout << "Number of reconstructed events after background subtraction                                           = " << sigEventsSub << std::endl; 
    std::cout << "Correction factor                                                                                     = " << Cw << std::endl;
    std::cout << "Number of events after efficiency correction                                                          = " << sigEvents << std::endl;
    std::cout << "Statistical uncertainty                                                                               = " << totalPercentStatError << "%" << std::endl;
    std::cout << "Systematic uncertainty                                                                                = " << totalPercentSystError << "%" << std::endl;
    std::cout << std::endl;

////////////////////////////
//Plot intermediate figures
////////////////////////////

	//Now that the histos have been filled and scaled,
	//we now plot them on a canvas
	
  	/*TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);

	leg->AddEntry(hdataSetc, "Data 2011", "pe");
	leg->AddEntry(hmcWSetc, "W#rightarrow#mu#nu", "f");
	leg->AddEntry(hmcZSetc, "Z#rightarrow#mu#mu", "f");
	leg->AddEntry(hmcQCDSetc, "QCD", "f");

  	hmcWSetc->GetXaxis()->SetTitle(sXtitle); 
        float entriesPerBin = (kinVarMax)/nBins;
	TString sY = "Events/"; TString sentriesPerBin = format(entriesPerBin); sY+=sentriesPerBin;
  	hmcWSetc->GetYaxis()->SetTitle(sY); 
	hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()+1.0e3); 
	hmcWSetc->Draw("hist f");
	hmcZSetc->Draw("hist fsame");
	hmcQCDSetc->Draw("hist fsame");
	hdataSetc->Draw("pesame");
	hmcWSetc->Draw("sameaxis");
        
	leg->Draw();

	TLatex l;
	l.SetNDC();
	l.DrawLatex(0.18,0.83,sSel + ( addPercent ? "%" : "" ));
	l.DrawLatex(0.169,0.767,sSel2);
	l.DrawLatex(0.65,0.89,"#sqrt{s_{NN}}=2.76 TeV");
	l.SetTextSize(0.034);
	l.DrawLatex(0.421,0.820,"#int Ldt #approx 0.14 nb^{-1}"); 

	ATLAS_LABEL(0.17,0.89,1);
	myText(0.33,0.89, (Color_t)kBlack, (char*)("Internal"));

	TString plotNameLog = "data"; plotNameLog+=sXtitle; plotNameLog+=sSel; plotNameLog+="_"; plotNameLog+=sSel2; plotNameLog+="_"; plotNameLog+="Log";
	if(doPreSelKinematics) plotNameLog+="_PreSel"; 
	TString plotNameLin = "data"; plotNameLin+=sXtitle; plotNameLin+=sSel; plotNameLin+="_"; plotNameLin+=sSel2; plotNameLin+="_"; plotNameLin+="Lin"; //.png";
	if(doPreSelKinematics) plotNameLin+="_PreSel"; 

	plotNameLin.ReplaceAll("/",",");plotNameLin.ReplaceAll("#","");
	cdata->Print(plotNameLin+".png");
	cdata->Print(plotNameLin+".eps"); 
	cdata->Print(plotNameLin+".pdf"); 
	cdata->Print(plotNameLin+".root"); 

	plotNameLog.ReplaceAll("/",",");plotNameLog.ReplaceAll("#","");
  	cdata->SetLogy(true); hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()*1e3); cdata->Update();
	cdata->Print(plotNameLog.ReplaceAll("#","")+".png");
	cdata->Print(plotNameLog.ReplaceAll("#","")+".eps"); 
	cdata->Print(plotNameLog.ReplaceAll("#","")+".pdf"); 
	cdata->Print(plotNameLog+".root"); 
*/
	std::cout << "Clean up" << std::endl;
    delete hmcQCDSet;

}//plotKinematics

// Calculate bin index corresponding
// to this eta and centrality bin
int getBinIndex(int ieta, int icent, int nCentralityBins, int nEtaBins, double etaBins[], bool doMirrorEta = true){

	       double etaMed = etaBins[ieta]+(etaBins[ieta+1]-etaBins[ieta])/2.0;
	       int index = -1; 
              
                if(!doMirrorEta){
                  int mappedEta;
                  ///hop over the crack region
                  ///in bin [-0.1,0.1]
                  if(etaMed==0.0) {std::cout << "Skipping crack region." << std::endl; return index; }
                  if(etaMed<0.0){
                      ///number of bins in absolute eta
                      ///used for trigger efficiency
                      mappedEta = indexNegativeEta(ieta,nEtaBins);
                      index = mappedEta*nCentralityBins+icent;
                  }
                  else {
                      mappedEta = indexPositiveEta(ieta,nEtaBins);
                      index = mappedEta*nCentralityBins+icent;
                  }
                }
	        else {
                  index = ieta*nCentralityBins+icent; 
                }
                if(index<0.0) {
                  std::cout << "ERROR: Invalid index. " << std::endl; exit(0);
                }
                
                return index;

}


////////////////////////////////////////////////////
//Plotter function for the eta distribution
//of muon signal candidates. Note: When plotting
//the eta distribution , only one eta bin is used
//since the if not the MC
//distribution is normalized to the data in every eta bin,
//thereby producing an MC distribution that exactly matches
//the data (which is not what we want)
//////////////////////////////////////////////////////
void plotWCandidateEta( RooDataSet* dataSet, RooDataSet* mcWSet, RooDataSet* mcTauSet, RooDataSet* mcZSet, 
                RooDataSet* mcJ1Set, RooDataSet* mcJ2Set, RooDataSet* mcJ3Set, 
				RooDataSet* mcJ1Events,RooDataSet* mcJ2Events, RooDataSet* mcJ3Events, 
                TFile* fQCDBkg, TGraphAsymmErrors* grBkgTau,TGraph2DErrors* grBkgZ, double totalZinEta, 
                double Cw, double cWSystErr, double qcdSystError, double zBosSystError,
				THStack* hmData, THStack* hmQCD,THStack* hmZ,THStack* hmTau, THStack* hmMcSum, const int iMt, const int iEta, 
                const int iCentrality, int nCentralityBins,int nEtaBins, float xLo,
				float kinVarLow, float kinVarMax , double ncoll, double ptLow , double ptUpp, 
                double etaLow, double etaUpp, double centralityLow, double centralityUpp, 
				RooRealVar& kinVar, TString sSel, TString sSel2 , int nBins, TString sXtitle, 
				bool addPercent = true, bool correctSpectra = true, bool doSubtractBkg = true, bool
                doPreSelKinematics=false, bool doMirrorEta=true) {


	std::cout << "Plotting for kinematic variable " << sXtitle << " in range " << kinVarLow << " to " << kinVarMax <<std::endl;
    // If plotting over all eta, change xBins and binning accordingly
    ///match bins to analysis bins
  	RooBinning b = RooBinning(kinVarLow,kinVarMax); 
    // Uncomment to plot over all eta
//    double xBins[] = {kinVarLow,-2.1,-1.85,-1.55,-1.3,-1.05,-0.8,-0.6,-0.35,-0.1,0.1,0.35,0.6,0.8,1.05,1.3,1.55,1.85,2.1,kinVarMax}; 
    ///use for plotting |eta|
    //double xBins[] = {kinVarLow,0.35,0.6,0.8,1.05,1.3,1.55,1.85,2.1,kinVarMax}; 
    double xBins[] = {kinVarLow,0.35,0.6,0.8,1.05,1.37,1.52,1.74,2.1,kinVarMax}; 
    //double xBins[] = {xLo,0.35,0.6,0.8,1.05,1.3,1.55,1.85,2.1,kinVarMax}; 
    int  binnum = sizeof(xBins)/sizeof(double) - 1;
    nBins = binnum;

    
    // Initialize trigger efficiencies
    if((nBins==9||nBins==19||(!doMirrorEta))&&nCentralityBins==6) readInputFile(nCentralityBins,9,doMirrorEta,"triggerEffZNote_v06.txt");
    else {
        std::cout << "There was error in reading in the trigger efficiencies. Check that you have the correct number of eta and centrality bins. " << std::endl;
        exit(0);
    }
    // Arrary of trigger efficiencies for each eta/centrality bin
    float trigTemp[nBins];
    for(int ieta = 0; ieta < nBins; ++ieta){

        // Turn on last argument for absolute eta
        int index = getBinIndex(ieta,iCentrality,nCentralityBins,nBins,xBins,true);
        if(index<0) {
            std::cout << "WARNING: Negative bin index. You are either in the crack region or something is wrong." << std::endl; 
        }
        else {
            trigTemp[ieta] = trigEfficiency(index);
            std::cout << "Trigger efficiency at eta bin: " << ieta << " centrality bin: " << iCentrality << " = " << trigTemp[ieta] << std::endl;
        }
    }//ieta

    ///IMPORTANT: Be sure to use exact limits or error will occur when 
    ///adding histograms (i.e. use kinVarLow and kinVarMax directly and NOT a number; e.g. 0.1,2.4)
/*    b.addUniform(1,kinVarLow,-2.1);
    b.addUniform(1,-2.1,-1.85);
    b.addUniform(1,-1.85,-1.55);
    b.addUniform(1,-1.55,-1.3);
    b.addUniform(1,-1.3,-1.05);
    b.addUniform(1,-1.05,-0.8);
    b.addUniform(1,-0.8,-0.6);
    b.addUniform(1,-0.6,-0.35);
    b.addUniform(1,-0.35,-0.1);
    b.addUniform(1,-0.1,0.1);
    b.addUniform(1,0.1,0.35);
*/
   b.addUniform(1,kinVarLow,0.35);
    b.addUniform(1,0.35,0.6);
    b.addUniform(1,0.6,0.8);
    b.addUniform(1,0.8,1.05);
    b.addUniform(1,1.05,1.37);
    b.addUniform(1,1.37,1.52);
    b.addUniform(1,1.52,1.74);
    b.addUniform(1,1.74,2.1);
    b.addUniform(1,2.1,kinVarMax);

	//initialize histograms	
  	// --- data ---
	TH1F* hdataSet = (TH1F*)dataSet->createHistogram("hdataSet",kinVar,Binning(b));
    hdataSet->Scale(1.0,"width");
    std::cout << "entries in data:" << dataSet->numEntries() << " : " << hdataSet->Integral() << std::endl;
  	// --- W set ---
	TH1F* hmcWSet = (TH1F*)mcWSet->createHistogram("hmcWSet",kinVar,Binning(b));
    hmcWSet->Scale(1.0,"width");
    std::cout << "entries in W MC:" << mcWSet->numEntries() << " : " << hmcWSet->Integral() << std::endl;
  	// --- Z set ---
	TH1F* hmcZSet = (TH1F*)mcZSet->createHistogram("hmcZSet",kinVar,Binning(b));
    //hmcZSet->Scale(1.0,"width");
    std::cout << "entries in Z MC:" << mcZSet->numEntries() << " : " << hmcZSet->Integral() << std::endl;
    // --- Wtau set ---
	TH1F* hmcTauSet = (TH1F*)mcTauSet->createHistogram("hmcTauSet",kinVar,Binning(b));
    //hmcTauSet->Scale(1.0,"width");
  	// --- QCD set ---
	TH1F* hmcJ1Set = (TH1F*)mcJ1Set->createHistogram("hmcJ1Set",kinVar,Binning(b));
    std::cout << "entries in J1mu MC:" << mcJ1Set->numEntries() << " : " << hmcJ1Set->Integral() << std::endl;
	TH1F* hmcJ2Set = (TH1F*)mcJ2Set->createHistogram("hmcJ2Set",kinVar,Binning(b));
    std::cout << "entries in J2mu MC:" << mcJ2Set->numEntries() << " : " << hmcJ2Set->Integral() << std::endl;
	TH1F* hmcJ3Set = (TH1F*)mcJ3Set->createHistogram("hmcJ3Set",kinVar,Binning(b));
    std::cout << "entries in J3mu MC:" << mcJ3Set->numEntries() << " : " << hmcJ3Set->Integral() << std::endl;

	//return correctly weighted QCD histogram
    std::cout << "Weighting QCD eta distributions with Jxmu samples..." << std::endl;
	//TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",nBins, xBins);
	TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",binnum, xBins);
    std::cout << "Bins must be equal. Let's check: " << hmcQCDSet->GetNbinsX() << " =? " << hmcJ1Set->GetNbinsX() << " =? " << hmcJ2Set->GetNbinsX() 
        << " =? " << hmcJ3Set->GetNbinsX() << std::endl;
	hmcQCDSet = getWeightedJxVariableBinnedHisto(xBins,hmcQCDSet->GetNbinsX(),centralityLow, centralityUpp, ncoll, hmcJ1Set, hmcJ2Set,
        hmcJ3Set, mcJ1Events, mcJ2Events,mcJ3Events, nBins,kinVarLow, kinVarMax);
    //hmcQCDSet->Scale(1.0,"width");
    std::cout << "Done." << std::endl;



    ///return the lower bin number
    ///which you would like to integrate
    ///from in order to obtain signal yield
    int binLo = b.binNumber(kinVarLow)+1;

    std::cout << "Signal bins = " << binLo << "-" << nBins << std::endl;
//	std::cout << "nBins: " << nBins << "\n" << "Binning from " << xLo << " to " << kinVarMax << std::endl;
	std::cout << "Corresponding to eta range " << kinVarLow << " to " << kinVarMax << std::endl;
//	std::cout << "nBins: " << nBins << "\n" << "Integrating from " << kinVarLow << " to " << kinVarMax << std::endl;
//	std::cout << "lower cut bin corresponding to eta bin median " << xEta << " = " << binLo << std::endl;

	TH1F* hdataSetc = (TH1F*)hdataSet->Clone("hdataSetc");

    ///Number of observed events in signal region
    ///before efficiency corrections
  	double sigEventsUncorr = hdataSetc->Integral(binLo,nBins); 
  	std::cout << "Integrated events before correction in eta:"<< kinVarLow << "-" << kinVarMax  << " =  " 
        << sigEventsUncorr << " +-" << TMath::Sqrt(sigEventsUncorr) << std::endl;
	
	double xEta = etaLow+(etaUpp-etaLow)/2.0;
	//relative errors of uncorrected signal events and correction factor
	double err1 = TMath::Sqrt(sigEventsUncorr); 
    ///statistical errors of the raw number of counts
    double errSigStat = err1;
	
//////////
//QCD
/////////

	//graph of QCD bkg fraction as fcn of eta for given charge and centrality class  
	///use inclusive charge set (i.e. iMt=104)
    TString sFracQCD = "fractionQCD"; sFracQCD+="_charge"; sFracQCD+=104; sFracQCD+="_eta"; sFracQCD+=0; sFracQCD+="_cent"; 

	sFracQCD+=iCentrality;

	TGraphErrors* grBkgQCD = 0;
    if(doSubtractBkg) grBkgQCD = (TGraphErrors*) fQCDBkg->Get(sFracQCD);

	TH1F* hmcQCDSetc = (TH1F*)hmcQCDSet->Clone("hmcQCDSetc");
	// Primary statistics for relative stat. error on QCD
    double nJ1CountsAfterSel=0.0, nJ2CountsAfterSel=0.0, nJ3CountsAfterSel=0.0;
    nJ1CountsAfterSel = mcJ1Set->numEntries();
    nJ2CountsAfterSel = mcJ2Set->numEntries();
    nJ3CountsAfterSel = mcJ3Set->numEntries();
    double nJxmuCountsAfterSel = nJ1CountsAfterSel+nJ2CountsAfterSel+nJ3CountsAfterSel;
    /// stat. error on Jx
    double relErrJxmuCountsAfterSel = TMath::Sqrt(nJxmuCountsAfterSel)/nJxmuCountsAfterSel;


	std::cout << "Getting QCD background fraction..." << std::endl;
    // Get f_qcd in this eta and centrality bin
	///use inclusive charge set (i.e. iMt=104)
    double survivalFractionQCD = 0.0;
	if(doPreSelKinematics) 
        survivalFractionQCD = getQCDBkg(104,getQCDBkgBinNumber(iEta,xEta,nEtaBins,doMirrorEta),iCentrality,grBkgQCD,centralityLow, centralityUpp, doPreSelKinematics) ; 	
	double dataQCDEvents = survivalFractionQCD*sigEventsUncorr;

    // Loop over each eta bin from the data distro 
    // and normalize QCD bin content to relative expected 
    double dataQCDEventsTot = 0;
  	for (int ieta = 1;  ieta <= hdataSetc->GetNbinsX(); ++ieta){

        if(!doSubtractBkg) break;
        // Retrieve QCD bkg fraction in this eta/centrlity bin
        survivalFractionQCD = getQCDBkg(104,getQCDBkgBinNumber(ieta-1,hdataSetc->GetBinCenter(ieta),hdataSetc->GetNbinsX(),doMirrorEta),
                                iCentrality,grBkgQCD,centralityLow, centralityUpp) ;
        // jump over crack region [-0.1,0.1]
        if(hdataSetc->GetBinCenter(ieta)==0.0) survivalFractionQCD=0.0;

        std::cout << "survivalFractionQCD in centrality bin: " << iCentrality << " eta bin: " << ieta-1 << " = " << survivalFractionQCD << std::endl;
        // Number of expected QCD events in this eta/centrality bin
        dataQCDEvents = survivalFractionQCD*hdataSetc->GetBinContent(ieta);
        std::cout << "Number of qcd events in centrality: " << iCentrality << " eta: " << ieta-1 << " = " << dataQCDEvents << std::endl;
        // Set eta bin content for qcd bkg histo
        hmcQCDSetc->SetBinContent(ieta,dataQCDEvents);
        hmcQCDSetc->SetBinError(ieta,relErrJxmuCountsAfterSel*dataQCDEvents);
        dataQCDEventsTot+=dataQCDEvents;
    }   
    ///error in the difference of the QCD bkg frac for mu+ and mu-; relative error shown below
//    double errQCDFrac = qcdSystError/survivalFractionQCD;

	///Error of dataQCDEvents
	double errQCDCountsSyst = 0.0;
	double errQCDCountsStat = 0.0;
	double totalPercentStatError = 0.0; 
    double totalPercentSystError = 0.0;

/*    if(doSubtractBkg) {
        errQCDCountsSyst = errQCDFrac*dataQCDEvents;
        errQCDCountsStat = TMath::Sqrt(dataQCDEvents);

        std::cout << "Percent systematic uncertainty of the data from QCD = " << errQCDCountsSyst/sigEventsUncorr*100 << "%" << std::endl;
        totalPercentSystError = errQCDCountsSyst/sigEventsUncorr*100.0;
    }
*/
    std::cout << "Current systematic error on the number of events in the data after QCD subtraction = +/- " 
            << totalPercentSystError << "%" << std::endl;

    ///Integrate QCD spectrum in signal region
  	double mcQCDEvents = hmcQCDSetc->Integral(binLo,nBins);
	std::cout << "QCD integral in signal region = " <<  mcQCDEvents << std::endl;

    ///Scale factor for scaling to number of
    ///expected QCD events in the data
	double sfQCD = dataQCDEvents/mcQCDEvents;

	if(mcQCDEvents==0) {
		std::cout << "WARNING: 0 QCD MC events in signal region." << std::endl;
		sfQCD = 1.0;
	}
 
    ///Normalize QCD MC shape to
    ///number of expected QCD events in
    ///the data
	//hmcQCDSetc->Scale(sfQCD);

    ///Keep a running sum of
    ///QCD bkg histogram for
    ///final figure
	hmQCD->Add(hmcQCDSetc);

	hmcQCDSetc->SetFillColor(kAzure-9);
	std::cout << "Expected number of QCD events in the signal region of the data = " << dataQCDEvents << std::endl;
  	std::cout << "Background percentage from QCD in eta :"<< kinVarLow << "-" << kinVarMax  << " =  " <<
        dataQCDEvents/sigEventsUncorr*100.0 << "%" <<std::endl;

///////////
///Z boson
//////////

	TH1F* hmcZSetc = (TH1F*)hmcZSet->Clone("hmcZSetc");

	double survivalProb; 
	int indexBin = -9999;

    ///return the probability of generated Z 
    ///to survivie W selection cuts
    double dataMuZEventsTot = 0;
	if(doPreSelKinematics) survivalProb = getZBkg(iMt,iEta,iCentrality,indexBin,grBkgZ,centralityLow, centralityUpp,doPreSelKinematics) ;
	else if(doSubtractBkg) {
        int indexEta;
        if(!doMirrorEta){
          ///hop over the crack region
          ///in bin [-0.1,0.1]
          if(xEta==0.0) std::cout << "WARNING: Bin in crack region. Expect wrong result!" << std::endl;
          if(xEta<0.0){
              ///number of bins in absolute eta
              indexEta = indexNegativeEta(iEta,nEtaBins);
	          indexBin = indexEta*nCentralityBins+iCentrality;
          }
           else {
              indexEta = indexPositiveEta(iEta,nEtaBins);
	          indexBin = indexEta*nCentralityBins+iCentrality;
           }
        }
        else {
          indexEta = iEta;
          indexBin = iEta*nCentralityBins+iCentrality;
        }
        // Loop over each eta bin from the data distro 
        // and normalize Z bin content to relative expected 

        double errZEventsFromData = getZCentDataStatErr(centralityLow, centralityUpp);
  	    for (int ieta = 1;  ieta <= hdataSetc->GetNbinsX(); ++ieta){

            indexEta=getZBkgBinNumber(ieta-1,hdataSetc->GetBinCenter(ieta),hdataSetc->GetNbinsX());
            indexBin = (indexEta)*nCentralityBins+iCentrality;

            // Retrieve Z bkg fraction in this eta/centrlity bin
            survivalProb = getZBkg(iMt,indexEta,iCentrality,indexBin, grBkgZ, centralityLow, centralityUpp) ; //percent of Z in W signal region
            // jump over crack region [-0.1,0.1]
            if(hdataSetc->GetBinCenter(ieta)==0.0) survivalProb=0.0;
            std::cout << "survivalFractionZ in centrality bin: " << iCentrality << " eta bin: " << ieta-1 << " = " << survivalProb << std::endl;

            double errStatSurvivalProbMC = 
                getZBkgStatErrorMC(iMt,indexEta,iCentrality,indexBin, grBkgZ, centralityLow,centralityUpp)/survivalProb ; //relative of Z in W signal region

            ///Apply this fraction to the number of Zs in this centrality and eta bin 
            ///to afford the number of muons from Zs contaminating sample
            ///This comes from the published Z note
            double nZEventsFromData = getZEventsData(centralityLow, centralityUpp, totalZinEta);

            ///Expected number of one-legged muons in the data
            double dataMuZEvents = survivalProb*nZEventsFromData;
            //Split in half if plotting full eta
            if(hdataSetc->GetNbinsX()>9) dataMuZEvents/=2.0;
            double errZCountsStat = TMath::Sqrt( TMath::Power(errZEventsFromData,2) + TMath::Power(errStatSurvivalProbMC,2))*dataMuZEvents;  
            // divide by bin width
            dataMuZEvents/=hdataSetc->GetBinWidth(ieta);
            errZCountsStat/=hdataSetc->GetBinWidth(ieta);
            std::cout << "dataMuZEvents in centrality bin: " << iCentrality << " eta bin: " << ieta-1 << " = " << dataMuZEvents << std::endl;

            // Set eta bin content for qcd bkg histo
            hmcZSetc->SetBinContent(ieta,dataMuZEvents);
            hmcZSetc->SetBinError(ieta,errZCountsStat);
            dataMuZEventsTot+=dataMuZEvents;
        } //ieta 
    }
    else survivalProb = 0.0;

	double errZCountsStat = 0.0;
	double errZCountsSyst = 0.0;

/*    if(doSubtractBkg) {

        ///stat. error
        errZCountsStat = TMath::Sqrt(dataZEvents);

        ///Propagate syst. error 
        errZCountsSyst = zBosSystError/bkgZEventFrac*dataZEvents;  
        std::cout << "Percent systematic uncertainty from Z = " << errZCountsSyst/sigEventsUncorr*100 << "%" << std::endl;

        totalPercentSystError = TMath::Sqrt( TMath::Power(totalPercentSystError,2) +
            TMath::Power(errZCountsSyst/sigEventsUncorr*100,2));
    }
*/
    std::cout << "Current systematic error on the number of events in the data after QCD and Z subtraction = +/- " 
            << totalPercentSystError << "%" << std::endl;

    ///Number of Z events in the signal region from MC
  	double mcZEvents = hmcZSetc->Integral(binLo,nBins);

    ///Scale factor for scaling MC events
    ///to the expected number of Z background muons
    ///in the data
	double sfZ = dataMuZEventsTot/mcZEvents;

	if(mcZEvents==0) {
		std::cout << "WARNING: 0 Z MC events in signal region." << std::endl;
		sfZ = 1.0;
	}

    ///Scale the histogram to the expected number of Z background muons
    ///in the data
	//hmcZSetc->Scale(sfZ);

    ///Add the QCD portion to afford a "total" background
    ///histogram.
	hmcZSetc->Add(hmcQCDSetc);
	hmcZSetc->SetFillColor(kRed);

	///Keep a running sum of Z+QCD histograms for
    ///the final figure integrated over all eta and centrality
	hmZ->Add(hmcZSetc);

  	std::cout << "Background Z percentage expected in the data from :"<< kinVarLow << "-" << kinVarMax  << " =  " <<
        dataMuZEventsTot/sigEventsUncorr*100.0 << "%" << std::endl;

////////////////////////
//Background subtraction
////////////////////////

    ///Histogram with QCD+Z component subtracted off
	TH1F* hdataSetSubtracted = (TH1F*)hdataSetc->Clone("hdataSetSubtracted");

    ///Subtract the QCD+Z histogram from the data
    ///Expected number of reconstructed signal events in data
    ///after background subtraction.
	double sigEventsSubTot = sigEventsUncorr-dataQCDEventsTot-dataMuZEventsTot;
    std::cout << "Expected number of reconstructed W events in the data after background subtraction = " 
        << sigEventsSubTot << std::endl;
    double dataTauEventsTot=0;
  	for (int ieta = 1;  ieta <= hdataSetc->GetNbinsX(); ++ieta){
        
        float dataTauEvents=0.0;

	    // Number of qcd,Z bkg-subtracted sig events in this eta bin
        double sigEventsSub = hdataSetc->GetBinContent(ieta)-hmcZSetc->GetBinContent(ieta);

        ///Fetch %tau contamination per W->mu event in this eta/centralaty bin
        int index = getTauBkgBinNumber(ieta-1,hdataSetc->GetBinCenter(ieta),hdataSetc->GetNbinsX());
        float tauBkgFraction = grBkgTau->GetY()[index];
        // jump over crack region [-0.1,0.1]
        if(hdataSetc->GetBinCenter(ieta)==0.0) tauBkgFraction=0.0;
        std::cout << "Background fraction from W-->tau-->mu decays: " << tauBkgFraction << std::endl;

        ///Number of tau events in the data in this eta bin
        dataTauEvents = tauBkgFraction/(1.0+tauBkgFraction)*sigEventsSub;
        double errStatTauEvents=0.;
        if(grBkgTau->GetEYhigh()[index]>grBkgTau->GetEYlow()[index]) errStatTauEvents=grBkgTau->GetEYhigh()[index];
        else errStatTauEvents=grBkgTau->GetEYlow()[index];
        errStatTauEvents=sigEventsSub*errStatTauEvents/TMath::Power(1.0+tauBkgFraction,2)/dataTauEvents;
        std::cout << "Number of tau events in centrality: " << iCentrality << " eta: " << ieta-1 << " = " << dataTauEvents << std::endl;
        hmcTauSet->SetBinContent(ieta,dataTauEvents);
        hmcTauSet->SetBinError(ieta,errStatTauEvents);
        dataTauEventsTot+=dataTauEvents;
    }//ieta
    ///Subtract off tau background
    sigEventsSubTot-=dataTauEventsTot;
    std::cout << "Expected number of reconstructed W events in the data after QCD,Z, and tau background subtraction = "
        << sigEventsSubTot << std::endl;

    ///Absolute systematic errors due to the isolation cut.
    ///This is used to calculate absolute error in the subtracted yield
    ///for this charge,eta,and centrality class
    double isoPercentSystErr = getIsoSystErr(iMt, indexBin);
    std::cout << "Systematic uncertainty of number of signal events in data due to isolation efficiency = " << isoPercentSystErr << "%" << std::endl;

    ///Relative systematic errors due to the MPT cut.
    ///This is used to calculate absolute error in the subtracted yield
    ///for this charge,eta,and centrality class
    double mptPercentSystErr = getMptSystErr(iMt, indexBin) ;
    std::cout << "Systematic uncertainty of number of signal events in data due to MPT resolution = " << mptPercentSystErr << "%" << std::endl;

/////////////////////////
//W-->tau-->mu background
/////////////////////////
    std::cout << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "W --> T a u  B a c k g r o u n d" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << std::endl;

	std::cout << "W-->tau integral " << hmcTauSet->Integral() << std::endl;
	TH1F* hmcTauSetc = (TH1F*)hmcTauSet->Clone("hmcTauSetc");
    ///Number of tau->mu in signal region from MC
    double mcTauEvents = hmcTauSetc->Integral(binLo,nBins);
    ///Scale MC events to expected
    double sfTau = dataTauEventsTot/mcTauEvents;

	if(mcTauEvents==0) {
		std::cout << "WARNING: 0 Tau MC events in signal region." << std::endl;
		sfTau = 1.0;
	}

    ///Scale the histogram to the expected number of Tau background muons
    ///in the data
	//hmcTauSetc->Scale(sfTau);

    ///Add the QCD+Z portion to afford a "total" background
    ///histogram.
	hmcTauSetc->Add(hmcZSetc);
	hmcTauSetc->SetFillColor(kYellow-7);

    if(doSubtractBkg){
      hdataSetSubtracted->Add(hmcTauSetc,-1);
      //Consistency check: 
      std::cout << "Consistency Check: " << std::endl;
      std::cout << sigEventsSubTot << " =? " << hdataSetSubtracted->Integral(binLo,nBins) << std::endl;
    }
	///Keep a running sum of Z+QCD+Tau histograms for
    ///the final figure integrated over all eta and centrality
	hmTau->Add(hmcTauSetc);

  	std::cout << "Background tau percentage expected in the data from mT:" <<
        dataTauEventsTot/sigEventsUncorr*100.0 << "%" << std::endl;

///////////////
//Cw correction
//////////////

	double effCorrection ; double effCorrectionSystErr; 

	if(correctSpectra) {

        ///Correct background subtracted histogram with Cw
		//hdataSetc = correctEfficiencyMt(hdataSetSubtracted, 104, iCentrality, iEta,Cw);
        hdataSetSubtracted->Scale(1.0/Cw);
		effCorrection = Cw;

        ///Take systematic error on Cw
        ///as statistical error from W MC
        ///in eta/centrality bin
	    effCorrectionSystErr = cWSystErr;

	    std::cout << "Efficiency correction factor = " << effCorrection << " +/- " << effCorrectionSystErr << std::endl;
	}

	else {effCorrection = 1.0; effCorrectionSystErr = 0.0;}

    ///Number of background subtracted and corrected events in the data in the signal region
  	double sigEvents = hdataSetSubtracted->Integral(binLo,nBins); 
  	//double sigEvents = sigEventsSub/Cw; 

    double cWPercentSystErr = effCorrectionSystErr/effCorrection*100.0;

    ///Add in systematic uncertainty from Cw statistical error
    totalPercentSystError = TMath::Sqrt( TMath::Power(totalPercentSystError,2) + TMath::Power(cWPercentSystErr,2) );

    std::cout << "Current statistical error on the number of events in the data after all background subtraction = +/- " 
            << totalPercentStatError << "%" << std::endl;
    std::cout << "Current systematic error on the number of events in the data after QCD,Z,isolation, MPT, and Cw errors propagated = +/- " 
            << totalPercentSystError << "%" << std::endl;

	///Running sum of data histograms from 
    ///each eta and centrality bin 
	//hmData->Add(hdataSetSubtracted);
    hmData->Add(hdataSetc);
	TH1F* hmDataTemp =(TH1F*)hmData->GetStack()->Last()->Clone("hmDatac");
    std::cout << "Number of total entries in data after bin " << iEta << ":" << iCentrality << " = " << hmDataTemp->Integral() << std::endl;

//////////
//W 
/////////

    ///We now have a histogram with subtracted data corrected. 
    ///We also have QCD,tau, and Z histograms in the same canvas (normalized to the number of 
    ///expected background events). Below, we scale the W MC to the data and see how 
    ///well they agree.

	TH1F* hmcWSetc = (TH1F*)hmcWSet->Clone("hmcWSetc");
 
    // Jump over crack region [-0.1,0.1]
    for(int ieta = 1; ieta <= hdataSetc->GetNbinsX(); ++ieta){
       if(hmcWSetc->GetBinCenter(ieta)==0.0) hmcWSetc->SetBinContent(ieta,0.0);
       // Shape changes due to trigger-induced losses 
       std::cout << "W MC bin content : " << hmcWSetc->GetBinContent(ieta) << std::endl;
       std::cout << "Trigger Eff of eta bin : " << ieta << "  = " << trigTemp[ieta-1]<< std::endl;
       hmcWSetc->SetBinContent(ieta,hmcWSetc->GetBinContent(ieta)*trigTemp[ieta-1]);
       std::cout << "Bin content of MC after trigger losses applied: " << hmcWSetc->GetBinContent(ieta) << std::endl; 
    } //ieta

    ///Number of W events in MC in signal region
  	double mcWEvents =  hmcWSetc->Integral(binLo,nBins);


    ///Scale factor for scaling W MC events to 
    ///the number of expected signal events in the data.
	double sfW = sigEvents/mcWEvents;

	hmcWSetc->Scale(sfW);

    ///Number of expected W events from the MC
    double mcWExpected = hmcWSetc->Integral(binLo,nBins);
    std::cout << "Number of expected W bosons from the Monte Carlo = " << mcWExpected << std::endl;
    ///Add Z,Tau, and QCD histograms to obtain a 
    ///final MC histogram composed of 
    ///the expected number of Background+Signal events
	hmcWSetc->Add(hmcTauSetc);

	//running sum of combined(QCD+Z+W) MC histos	
    ///Keep a running sum of histograms for final figure
    ///integrated over all eta and centrality
	hmMcSum->Add(hmcWSetc);	


    std::cout << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "S u m m a r y" << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << "Expected number of events for charge:eta:centrality" << iMt << ":" << iEta << ":" << iCentrality << " = " << mcWExpected << std::endl;
    std::cout << "Number of reconstructed events before background subtraction                                          = " << sigEventsUncorr << std::endl; 
    std::cout << "Number of reconstructed events after background subtraction                                           = " << sigEventsSubTot << std::endl; 
    std::cout << "Correction factor                                                                                     = " << Cw << std::endl;
    std::cout << "Number of events after efficiency correction                                                          = " << sigEvents << std::endl;
    std::cout << "Statistical uncertainty                                                                               = " << totalPercentStatError << "%" << std::endl;
    std::cout << "Systematic uncertainty                                                                                = " << totalPercentSystError << "%" << std::endl;
    std::cout << std::endl;

/////////////////////////////
//Plot intermediate figures
/////////////////////////////

	//Now that the histos have been filled and scaled,
	//we now plot them on a canvas
	TLegend* leg = new TLegend(0.658, 0.637, 0.928, 0.867);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);

	leg->AddEntry(hdataSetc, "Data 2011", "pe");
	leg->AddEntry(hmcWSetc, "W#rightarrow#mu#nu", "f");
	leg->AddEntry(hmcTauSetc, "W#rightarrow#tau#nu", "f");
	leg->AddEntry(hmcZSetc, "Z#rightarrow#mu#mu", "f");
	leg->AddEntry(hmcQCDSetc, "QCD", "f");

  	hmcWSetc->GetXaxis()->SetTitle(sXtitle); 
	TString sY = "Muons";
  	hmcWSetc->GetYaxis()->SetTitle(sY); 
	hmcWSetc->GetYaxis()->SetRangeUser(0.1,110.0); 

    std::cout << "aaaaaaaaaaa"<< std::endl;
  	TCanvas* cdata = new TCanvas("cdata","cdata",600,600);
	hmcWSetc->Draw("hist f");
	hmcTauSetc->Draw("hist fsame");
	hmcZSetc->Draw("hist fsame");
	hmcQCDSetc->Draw("hist fsame");
	hdataSetc->Draw("pesame");
	hmcWSetc->Draw("sameaxis");
    std::cout << "bbbbbbbbbb"<< std::endl;
    
        
	leg->Draw();
	TLatex l;
	l.SetNDC();
	l.DrawLatex(0.18,0.83,sSel + ( addPercent ? "%" : "" ));
	l.DrawLatex(0.169,0.767,sSel2);
	l.SetTextSize(0.034);
	l.DrawLatex(0.74,0.89,"#sqrt{s_{NN}}=2.76 TeV");
	l.DrawLatex(0.492,0.89,"#int Ldt #approx 0.14 nb^{-1}"); 

	leg->Draw(); cdata->Update();
	ATLAS_LABEL(0.17,0.89,1);
	myText(0.33,0.89, (Color_t)kBlack, (char*)("Internal"));

	//TString plotNameLog = "dataMt_"; plotNameLog+=sSel2; plotNameLog+="_"; plotNameLog+=sSel; plotNameLog+="_"; plotNameLog+="Log"; if(doPreSelKinematics) plotNameLog+="_PreSel";
	TString plotNameLog = "dataEta_"; plotNameLog+="charge"; plotNameLog+=iMt; plotNameLog+="_eta"; plotNameLog+=iEta; plotNameLog+="_cent"; plotNameLog+=iCentrality;  
	if(doPreSelKinematics) plotNameLog+="_PreSel";
	//TString plotNameLin = "dataMt_"; plotNameLin+=sSel2; plotNameLin+="_"; plotNameLin+=sSel; plotNameLin+="_"; plotNameLin+="Lin"; if(doPreSelKinematics) plotNameLin+="_PreSel";
	TString plotNameLin = "dataEta_"; plotNameLin+="charge"; plotNameLin+=iMt; plotNameLin+="_eta"; plotNameLin+=iEta;
    plotNameLin+="_cent"; plotNameLin+=iCentrality;  

    ///uncomment for saving intermediate plots
  
    //hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()+75.0); cdata->Update();
/*    hmcWSetc->GetYaxis()->SetRangeUser(0.1, 90.0); cdata->Update();
	cdata->Print(plotNameLin.ReplaceAll("#","")+"_Lin.png");
	cdata->Print(plotNameLin.ReplaceAll("#","")+"_Lin.eps"); 
	cdata->Print(plotNameLin.ReplaceAll("#","")+"_Lin.pdf"); 
	cdata->Print(plotNameLin+".root"); 

//  	cdata->SetLogy(true); hdataSetc->GetYaxis()->SetRangeUser(0.1,2.1e5); cdata->Update();
    cdata->SetLogy(true); hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()*7.0e2); cdata->Update();
	cdata->Print(plotNameLog.ReplaceAll("#","")+"_Log.png");
	cdata->Print(plotNameLog.ReplaceAll("#","")+"_Log.eps"); 
	cdata->Print(plotNameLog.ReplaceAll("#","")+"_Log.pdf"); 
	TString plotNameLogRoot = plotNameLog.ReplaceAll("|",",") + ".root";
	cdata->Print(plotNameLogRoot); 
*/	

	std::cout << "Clean up" << std::endl;
    delete hmcQCDSet;
    delete cdata;

}///plotEta

////////////////////////////////////////////////////
//Plotter function for the missing pt (MPT) distribution
//for events with a muon signal candidate. 
//WORK IN PROGRESS (I.E. NOT READY FOR USE) 
//////////////////////////////////////////////////////
void plotWCandidateMPT( RooDataSet* dataSet, RooDataSet* mcWSet, RooDataSet* mcZSet, RooDataSet* mcJ1Set, RooDataSet* mcJ2Set, RooDataSet* mcJ3Set, 
				RooDataSet* mcJ1Events,RooDataSet* mcJ2Events, RooDataSet* mcJ3Events, TFile* fQCDBkg, TGraph2DErrors* grBkgZ, double totalZinEta, 
				THStack* hmData, THStack* hmQCD,THStack* hmZ, THStack* hmMcSum, const int iMt, const int iEta, const int
                iCentrality, int nCentralityBins,int nEtaBins, float xLo,
				float kinVarLow, float kinVarMax ,double Cw,double ncoll, double ptLow , double ptUpp, double etaLow, double etaUpp, double centralityLow, double centralityUpp, 
				RooRealVar& kinVar, TString sSel, TString sSel2 , float nBins, TString sXtitle, 
				bool addPercent = true, bool correctSpectra = true, bool doSubtractBkg = true, bool
                doPreSelKinematics=false,bool doMirrorEta=true) {


	std::cout << "Plotting for kinematic variable " << sXtitle << " in range " << kinVarLow << " to " << kinVarMax <<std::endl;
//	std::cout << "Z s.f. = " << sfZ << " QCD s.f. = " << sfQCD << " W s.f. = " << sfW << std::endl;

//	int binLo = nBins/kinVarMax*kinVarLow+1;
    int binLo = getBinNumber(kinVarLow,nBins,xLo,kinVarMax);


	std::cout << "nBins: " << nBins << " Binning from " << xLo << " to " << kinVarMax << std::endl;
	std::cout << "lower cut bin : " << binLo << std::endl;
  	RooBinning b = RooBinning(nBins,xLo,kinVarMax); 
  	//RooBinning b = RooBinning(xLo,kinVarMax); 

	//initialize histograms	
  	// --- data ---
	TH1F* hdataSet = (TH1F*)dataSet->createHistogram("hdataSet",kinVar ,Binning(b));
  	// --- W set ---
	TH1F* hmcWSet = (TH1F*)mcWSet->createHistogram("hmcWSet",kinVar,Binning(b));
  	// --- Z set ---
	TH1F* hmcZSet = (TH1F*)mcZSet->createHistogram("hmcZSet",kinVar,Binning(b));
  	// --- QCD set ---
	TH1F* hmcJ1Set = (TH1F*)mcJ1Set->createHistogram("hmcJ1Set",kinVar,Binning(b));
	TH1F* hmcJ2Set = (TH1F*)mcJ2Set->createHistogram("hmcJ2Set",kinVar,Binning(b));
	TH1F* hmcJ3Set = (TH1F*)mcJ3Set->createHistogram("hmcJ3Set",kinVar,Binning(b));

	//return correctly weighted QCD histogram
	TH1F* hmcQCDSet = new TH1F("hmcQCDSet","hmcQCDSet",nBins,xLo,kinVarMax);
	hmcQCDSet = getWeightedJxHisto(centralityLow, centralityUpp, ncoll, hmcJ1Set, hmcJ2Set, hmcJ3Set, mcJ1Events, mcJ2Events,mcJ3Events, nBins,xLo, kinVarMax);

	TH1F* hdataSetc = (TH1F*)hdataSet->Clone("hdataSetc");

  	double sigEventsUncorr = hdataSetc->Integral(binLo,nBins); //integrate from 40-200 GeV before Aw,Cw correction
  	std::cout << "integrated events before correction in mpt signal space = " << sigEventsUncorr << " +-" << TMath::Sqrt(sigEventsUncorr) << std::endl;

	
	double effMt ; double effMtErr; double effMtSystErr;
	double xEta = etaLow+(etaUpp-etaLow)/2.0;
	//return efficiency correction factor (Cw)
	if(correctSpectra) {
		//hdataSetc = correctEfficiencyMt(hdataSet, iMt, iCentrality, iEta);
		///use charge inclusive (i.e. iMt = 104)
		hdataSetc = correctEfficiencyMt(hdataSet, 104, iCentrality, iEta,Cw);
		//effMt = getEfficiencyMt(iMt, iCentrality, iEta);
		//effMt = getEfficiencyFitMt(iMt, iCentrality, xEta);
		effMt = Cw;
		//effMtErr = getCorrectionFactorError(iMt, iCentrality, iEta);
		effMtErr = 0.0; //hack
		//Cw syst errs added in quadrature
	    effMtSystErr = 0.0;
	}
	else {effMt = 1.0; effMtErr = 0.0; effMtSystErr = 0.0;}

  	double sigEvents = hdataSetc->Integral(binLo,nBins); 
	//relative errors of uncorrected signal events and correction factor
	double err1 = TMath::Sqrt(sigEventsUncorr)/sigEventsUncorr*100.0; double err2 = effMtErr/effMt*100.0;
	//propagated stat errors from signal and Cw
	double errStat = TMath::Sqrt( TMath::Power(err1,2)+TMath::Power(err2,2)  ) *0.01*sigEvents; 
	std::cout << "Efficiency correction factor = " << effMt << "+-" << effMtErr << "(stat.) " << effMtSystErr << "(syst.)" << std::endl;
    //double errStatFromCw = effMtSystErr/effMt*sigEvents;
 	//running sum of data histos
	hmData->Add(hdataSetc);
	TH1F* hmDataTemp =(TH1F*)hmData->GetStack()->Last()->Clone("hmDatac");
    std::cout << "Number of total entries in data after bin " << iEta << ":" << iCentrality << " = " << hmDataTemp->Integral() << std::endl;
	
    TString sFracQCD = "fractionQCD"; sFracQCD+="_charge"; sFracQCD+=104; sFracQCD+="_eta"; sFracQCD+=0; sFracQCD+="_cent"; 
	//hack for using 0-10% instead of 0-5%, 5-10%
	/*if(iCentrality==0)*/ sFracQCD+=iCentrality;
	//else sFracQCD+=iCentrality-1;
	TGraphErrors* grBkgQCD = 0;
    if(doSubtractBkg) grBkgQCD = (TGraphErrors*) fQCDBkg->Get(sFracQCD);

	TH1F* hmcQCDSetc = (TH1F*)hmcQCDSet->Clone("hmcQCDSetc");
	
	double survivalFractionQCD ; 
	//percent of QCD in W signal region b4 W selection
	///use inclusive charge set (i.e. iMt=104)
	std::cout << "Getting QCD background fraction..." << std::endl;
	if(doPreSelKinematics) survivalFractionQCD = getQCDBkg(104,getQCDBkgBinNumber(iEta,xEta,nEtaBins),iCentrality,grBkgQCD,centralityLow, centralityUpp, doPreSelKinematics) ; 	
	else if(doSubtractBkg) survivalFractionQCD = getQCDBkg(104, getQCDBkgBinNumber(iEta,xEta,nEtaBins), iCentrality,grBkgQCD,centralityLow, centralityUpp) ;
    else survivalFractionQCD = 0.0;

	double dataQCDEvents = survivalFractionQCD*sigEventsUncorr;
	//relative error of N_corr(sigEvents) and QCD bkg fraction(survivalFractionQCD)
	//double errNcorr = effMtSystErr/effMt; double errQCDFrac = qcdSystError/survivalFractionQCD;
	//propagated absolute error of dataQCDEvents
	//double errQCDCounts = 0.0;
    //if(doSubtractBkg) errQCDCounts = TMath::Sqrt( TMath::Power(errNcorr,2)+ TMath::Power(errQCDFrac,2))*dataQCDEvents;
	///debug
	//std::cout << "Relative error in corrected counts: " << errNcorr << " Relative error in qcd fraction: " << errQCDFrac << std::endl;
	//std::cout << "That gives an absolute error in the number of QCD bkg counts of " << errQCDCounts << std::endl;

	//error of qcd subtracted counts
	double errCountsQCDSub = 0.0; 
    //if(doSubtractBkg) errCountsQCDSub = TMath::Sqrt(TMath::Power(errNcorr*sigEvents,2) + TMath::Power(errQCDCounts,2) );
	//std::cout << "This gives an error in the qcd subtracted yield of " << errCountsQCDSub << std::endl;

  	double mcQCDEvents = hmcQCDSetc->Integral(binLo,nBins);
	std::cout << "QCD integral " <<  mcQCDEvents << std::endl;
	double sfQCD = dataQCDEvents/mcQCDEvents;
	if(mcQCDEvents==0) {
		std::cout << "WARNING: 0 QCD MC events in signal region." << std::endl;
		sfQCD = 1.0;
	}

	hmcQCDSetc->Scale(sfQCD);
    std::cout << "Number of QCD events in MC = " <<  hmcQCDSetc->Integral(binLo,nBins) << std::endl;
	hmQCD->Add(hmcQCDSetc);
	hmcQCDSetc->SetFillColor(kAzure-9);
	std::cout << "QCD events in signal region = " << dataQCDEvents << std::endl;
  	std::cout << "background from QCD :" << " =  " << dataQCDEvents/sigEvents*100.0 << "%" <<std::endl;

    ///Z boson
	std::cout << "Z integral " << hmcZSet->Integral() << std::endl;
	TH1F* hmcZSetc = (TH1F*)hmcZSet->Clone("hmcZSetc");

	double survivalProb; 
	int index = iEta*nCentralityBins+iCentrality;
	if(doPreSelKinematics) survivalProb = getZBkg(iMt,iEta,iCentrality,index,grBkgZ,centralityLow, centralityUpp,doPreSelKinematics) ;
	else if(doSubtractBkg) {
        int indexEta;
        if(!doMirrorEta){
          ///hop over the crack region
          ///in bin [-0.1,0.1]
          if(xEta==0.0) std::cout << "WARNING: Bin in crack region. Expect wrong result!" << std::endl;
          if(xEta<0.0){
              ///number of bins in absolute eta
              indexEta = indexNegativeEta(iEta,nEtaBins);
	          index = indexEta*nCentralityBins+iCentrality;
          }
           else {
              indexEta = indexPositiveEta(iEta,nEtaBins);
	          index = indexEta*nCentralityBins+iCentrality;
           }
        }

        else {
          indexEta = iEta;
          index = iEta*nCentralityBins+iCentrality;
        }

        survivalProb = getZBkg(iMt,indexEta,iCentrality,index, grBkgZ, centralityLow, centralityUpp) ; //percent of Z in W signal region
    }
    else survivalProb = 0.0;

	std::cout <<"Z survival probability:"<< survivalProb  << std::endl;
	double bkgZEventFrac = survivalProb*getZEventsData(centralityLow, centralityUpp, totalZinEta)/sigEventsUncorr;
	std::cout << "Fraction of Zs surviving cuts :" << bkgZEventFrac << std::endl;
	double dataZEvents = bkgZEventFrac*sigEventsUncorr ;

	//absolute err in number of surviving Z muons
	double errZCounts = 0.0;
    //if(doSubtractBkg) errZCounts = TMath::Sqrt( TMath::Power(errNcorr,2)+ TMath::Power(zBosSystError/bkgZEventFrac,2))*dataZEvents;

  	double mcZEvents = hmcZSetc->Integral(binLo,nBins);
	std::cout << dataZEvents << std::endl;
	std::cout << mcZEvents << std::endl;
	//weight Z mc by per event yield from paper
	double sfZ = dataZEvents/mcZEvents;
	if(mcZEvents==0) {
		std::cout << "WARNING: 0 Z MC events in signal region." << std::endl;
		sfZ = 1.0;
	}
	hmcZSetc->Scale(sfZ);
	hmcZSetc->Add(hmcQCDSetc);
	hmcZSetc->SetFillColor(kRed);

	//running sum of Z histos
	hmZ->Add(hmcZSetc);

  	std::cout << "background from Z MPT:"<< kinVarLow << "-" << kinVarMax  << " =  " << dataZEvents/sigEvents*100.0 << "%" << std::endl;

	TH1F* hmcWSetc = (TH1F*)hmcWSet->Clone("hmcWSetc");

    ///normalize to data
  	double mcWEvents =  hmcWSetc->Integral(binLo,nBins);
	double sigEventsSub = sigEvents-dataQCDEvents-dataZEvents;

    ///for this charge,eta,and centrality class
    double errIsoSub = getIsoSystErr(iMt, index)/100.0*sigEventsSub;
    std::cout << "Absolute systematic error in W yield due to isolation cuts: +- " << errIsoSub << std::endl;

	//absolute error in subtracted W yield
    //from each syst src added in quadrature
	double sigEventsSubErr = TMath::Sqrt(TMath::Power(errCountsQCDSub,2)+TMath::Power(errZCounts,2)+TMath::Power(errIsoSub,2)); 
	//double sigEventsSubErr = TMath::Sqrt(TMath::Power(errCountsQCDSub,2)+TMath::Power(errZCounts,2)+TMath::Power(0.0,2)); 

	double err3 = TMath::Sqrt(dataQCDEvents); double err4 = TMath::Sqrt(dataZEvents);
	double errStatSub = TMath::Sqrt( TMath::Power(errStat,2)+TMath::Power(err3,2)+TMath::Power(err4,2)  ) ; 
	double sfW = sigEventsSub/mcWEvents;
	hmcWSetc->Scale(sfW);
	//add bkg contribution
	hmcWSetc->Add(hmcZSetc);

	//running sum of combined(QCD+Z+W) MC histos	
	hmMcSum->Add(hmcWSetc);	

  	std::cout << "integrated signal after bkg subtraction in MPT :"<< kinVarLow << "-" << kinVarMax  << " =  " << sigEventsSub << "+-" << errStatSub << "(stat.) " 
		<< sigEventsSubErr << "(syst.)" << std::endl;

	mcWEvents = hmcWSetc->Integral(binLo,nBins);

//  	std::cout << "integrated events from W MC MPT:"<< kinVarLow << "-" << kinVarMax  << " =  " << mcWEvents << std::endl;

//  	std::cout << "integrated signal after bkg subtraction in MPT:"<< kinVarLow << "-" << kinVarMax  << " =  " << sigEventsSub << "+-" << errStatSub << std::endl;

	//now that the histos have been filled and scaled,
	//we now plot them on a canvas
	
  /*	TLegend* leg = new TLegend(0.659, 0.713, 0.9295, 0.8636);
	leg->SetTextFont(gStyle->GetTextFont());
	leg->SetTextSize(gStyle->GetTextSize());
	leg->SetBorderSize(0);
	leg->SetFillColor(0);

	leg->AddEntry(hdataSetc, "Data 2011", "pe");
	leg->AddEntry(hmcWSetc, "W#rightarrow#mu#nu", "f");
	leg->AddEntry(hmcZSetc, "Z#rightarrow#mu#mu", "f");
	leg->AddEntry(hmcQCDSetc, "QCD", "f");

  	hmcWSetc->GetXaxis()->SetTitle(sXtitle); 
        float entriesPerBin = (kinVarMax)/nBins;
	TString sY = "Events/"; TString sentriesPerBin = format(entriesPerBin); sY+=sentriesPerBin;
  	hmcWSetc->GetYaxis()->SetTitle(sY); 
	hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()+1.0e3); 
	hmcWSetc->Draw("hist f");
	hmcZSetc->Draw("hist fsame");
	hmcQCDSetc->Draw("hist fsame");
	hdataSetc->Draw("pesame");
	hmcWSetc->Draw("sameaxis");
        
	leg->Draw();

	TLatex l;
	l.SetNDC();
	l.DrawLatex(0.18,0.83,sSel + ( addPercent ? "%" : "" ));
	l.DrawLatex(0.169,0.767,sSel2);
	l.DrawLatex(0.65,0.89,"#sqrt{s_{NN}}=2.76 TeV");
	l.SetTextSize(0.034);
	l.DrawLatex(0.421,0.820,"#int Ldt #approx 0.14 nb^{-1}"); 

	ATLAS_LABEL(0.17,0.89,1);
	myText(0.33,0.89, (Color_t)kBlack, (char*)("Internal"));

	TString plotNameLog = "data"; plotNameLog+=sXtitle; plotNameLog+=sSel; plotNameLog+="_"; plotNameLog+=sSel2; plotNameLog+="_"; plotNameLog+="Log";
	if(doPreSelKinematics) plotNameLog+="_PreSel"; 
	TString plotNameLin = "data"; plotNameLin+=sXtitle; plotNameLin+=sSel; plotNameLin+="_"; plotNameLin+=sSel2; plotNameLin+="_"; plotNameLin+="Lin"; //.png";
	if(doPreSelKinematics) plotNameLin+="_PreSel"; 

	plotNameLin.ReplaceAll("/",",");plotNameLin.ReplaceAll("#","");
	cdata->Print(plotNameLin+".png");
	cdata->Print(plotNameLin+".eps"); 
	cdata->Print(plotNameLin+".pdf"); 
	cdata->Print(plotNameLin+".root"); 

	plotNameLog.ReplaceAll("/",",");plotNameLog.ReplaceAll("#","");
  	cdata->SetLogy(true); hmcWSetc->GetYaxis()->SetRangeUser(0.1,hmcWSetc->GetMaximum()*1e3); cdata->Update();
	cdata->Print(plotNameLog.ReplaceAll("#","")+".png");
	cdata->Print(plotNameLog.ReplaceAll("#","")+".eps"); 
	cdata->Print(plotNameLog.ReplaceAll("#","")+".pdf"); 
	cdata->Print(plotNameLog+".root"); 
*/
	std::cout << "Clean up" << std::endl;
    delete hmcQCDSet;
	/*delete gDirectory->FindObject("hdataSet");
	delete gDirectory->FindObject("hmcWSet");
	delete gDirectory->FindObject("hmcZSet");
	delete gDirectory->FindObject("hmcJ1Set");
	delete gDirectory->FindObject("hmcJ2Set");
	delete gDirectory->FindObject("hmcJ3Set");
	delete gDirectory->FindObject("hmcQCDSet");
*/
}///plotKinematics


// --- EOF --- //
