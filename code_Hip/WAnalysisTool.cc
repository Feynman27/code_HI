#define WAnalysisTool_cc
#include "WAnalysisTool.h"

WAnalysisTool::WAnalysisTool(RooDataSet* dataSet):
    _dataSet(dataSet),m_doBkgSubtraction(false),m_NZBkg(0.0),m_NTauBkg(0.0),m_NBkgA(0.0),m_NBkgB(0.0), m_dataSetIsCorrected(false), m_correctedSum(0)
{
	
//	_hDataSet = new TH1F();
	w = new RooWorkspace("w","workspace");
	_hsDataSetSum = new THStack();
	std::cout << "[WAnalysisTool::WAnalysisTool] Dataset initialized." << std::endl;
	_dataSet->Print();
	
}

WAnalysisTool::WAnalysisTool(RooDataSet* dataSet,RooRealVar& var):_dataSet(dataSet),_var(var)
{
	
	_hsDataSetSum = new THStack();
	std::cout << "[WAnalysisTool::WAnalysisTool] Dataset initialized." << std::endl;
	_dataSet->Print();
	
}
// Copy Ctor
WAnalysisTool::WAnalysisTool(const WAnalysisTool& CopySource){

/*	if(CopySource._hDataSet != NULL){
		// ensure deep copy by first allocating own _hDataSet
		_hDataSet = new TH1F();
		// copy from source into local THStack
		_hDataSet = (TH1F*)_hDataSet->Clone();
	}
*/
	if(CopySource._hsDataSetSum != NULL){
		// ensure deep copy by first allocating own _hsDataSetSum
		_hsDataSetSum = new THStack();
		// copy from source into local THStack
		_hsDataSetSum = (THStack*)CopySource._hsDataSetSum->Clone();
	}
}

WAnalysisTool::WAnalysisTool(){

}



void WAnalysisTool::fit(TFile* fppPdf, TFile* fnpPdf, double cutLow, double cutUpp){
	
	// Retrieve the workspace from files for pp,np(pn) 
	RooWorkspace* wpp = (RooWorkspace*) fppPdf->Get("w");
	RooWorkspace* wnp = (RooWorkspace*) fnpPdf->Get("w");

	// Retrieve the signal pdf for pp and np(pn) from the workspace
	RooKeysPdf sigppPdf = *((RooKeysPdf*)wpp->pdf("sigPdf"));
	RooKeysPdf signpPdf = *((RooKeysPdf*)wnp->pdf("sigPdf"));

	RooRealVar ppSig("ppSig","Number of W events from pp",350,0,1e9);
	RooRealVar npSig("npSig","Number of W events from np",350,0,1e9);

	// Build model
	RooArgList modelComponents;
	modelComponents.add(sigppPdf);
	modelComponents.add(signpPdf);
	// Add pdfs for merged pp,np(pn) model
	RooAddPdf mdlPdf("mdlPdf","mdlPdf",modelComponents,RooArgList(ppSig,npSig));
	
	// Create nll object
	std::cout << "[WAnalysisTool::fit] Creating NLL." << std::endl;
	RooAbsReal* nll = mdlPdf.createNLL(*_dataSet,Range(cutLow,cutUpp), NumCPU(4));

	// Interface to minuit
	RooMinuit m1(*nll) ;

	// Run migrad minimization
	m1.migrad();

	// Calculate minos errors
	std::cout << "[WAnalysisTool::fit] Calculating errors with minos." << std::endl;
	m1.minos(RooArgSet(ppSig,npSig));
	RooFitResult* roofitResult = m1.save();

	std::cout << "[WAnalysisTool::fit] Number of W events from pp: " << ppSig.getVal() << " +" << ppSig.getAsymErrorHi() << " -" << ppSig.getAsymErrorLo() << std::endl;
	std::cout << "[WAnalysisTool::fit] Number of W events from np: " << npSig.getVal() << " +" << npSig.getAsymErrorHi() << " -" << npSig.getAsymErrorLo() << std::endl;
	std::cout << std::endl;
	
}

// Kernel estimation of pdf
void WAnalysisTool::buildTemplate(RooRealVar& var, const char* fileName){
	
	// Initialize with number of entries in mcSet
	RooRealVar nSignal("nSignal","nSignal",_dataSet->numEntries());

	// Create workspace that may be retrieved for later fitting
	//RooWorkspace *w = new RooWorkspace("w","workspace");
	w->import(var);
	
	// Build template pdf using kernel estimation
	std::cout << "Building the RooKeysPdf's... This will take some time" << std::endl;
	RooKeysPdf sigPdf("sigPdf","sigPdf",var,*_dataSet);

	// Write to file
	std::cout << "Writing the templates to file " << fileName << std::endl;
	w->import(sigPdf);
	w->import(nSignal);
	w->writeToFile(fileName);
	std::cout << "Done." << std::endl;
	
}

void WAnalysisTool::weightDS(double wt)
{
	RooRealVar dsWt("dsWt","dsWt",0.0,10.0);
	dsWt.setVal(wt);
	std::cout << "DS weight: " << dsWt.getVal() << std::endl;
	_dataSet->addColumn(dsWt);
	RooDataSet* wSet = new RooDataSet(_dataSet->GetName(),_dataSet->GetTitle(),_dataSet,*_dataSet->get(),0,dsWt.GetName()) ;
	std::cout << "[WAnalysisTool::WAnalysisTool] Weighting dataset with weight " << dsWt.getVal() << std::endl;
	_dataSet = wSet;
	_dataSet->Print();
}

void WAnalysisTool::AssignDataSet(RooDataSet* dataSet){
	std::cout << "[WAnalysisTool::AssignDataSet] Reassigning dataset." << std::endl;
	_dataSet = dataSet;
}

void WAnalysisTool::AddToStack(TH1F* h){
	_hsDataSetSum->Add(h);
	// Assign _hDataSet to last added histogram
	//_hDataSet = (TH1F*)_hsDataSetSum->GetStack()->Last();
	//std::cout << "[WAnalysisTool::AddToStack] Current integral of cumulative histogram:" 
	//		<< (TH1F*)_hsDataSetSum->GetStack()->Last()->Integral() << std::endl;
}


void WAnalysisTool::FillHistogram(RooBinning b, RooRealVar& var, double cutLow, double cutUpp)
{

	//_var = var;
	_cutLow = cutLow;
	_cutUpp = cutUpp;
	_b = b;

	std::cout << "[WAnalysisTool::FillHistogram] Initializing histograms..." << std::endl;
	
	TString hname = "h"+(TString)_dataSet->GetTitle();
    // Check if hist with same name exists
    // and if so, delete it 
    TString checkName = hname; checkName +="__"; checkName+=(TString)var.GetName();
    if ( gROOT->FindObject(checkName) ) { 
        std::cout << "found object: " << checkName << std::endl;
        delete gROOT->FindObject(checkName);
    } 
	_hDataSet = (TH1F*)_dataSet->createHistogram(hname, var , Binning(b));
	std::cout << "binning in range: [" << _cutLow << "," << _cutUpp << "]" << std::endl;
	std::cout << "bins in " << checkName << ": " << _hDataSet->GetNbinsX() << std::endl;

	// Find lower bin boundary in signal region
	_binLo = _hDataSet->FindBin(_cutLow);
	_binUp = _hDataSet->FindBin(_cutUpp);
	std::cout << "Bin low: " << _binLo << " Bin upp: " << _binUp << std::endl;
	
	std::cout << "[WAnalysisTool::FillHistogram] Number of observed signal candidates: " << _hDataSet->Integral(_binLo,_binUp) << std::endl;
		
}

TCanvas* WAnalysisTool::plotOnCanvas(RooRealVar& var){
	std::cout << "[WAnalysisTool::plotOnCanvas]" << std::endl;

	TCanvas* canvas = new TCanvas();
	canvas->cd();
	_hDataSet->Draw("pe");
	_hDataSet->GetXaxis()->SetTitle(var.getTitle());
	return canvas;
}

void WAnalysisTool::ApplySelectionCuts(RooFormulaVar cutsSig, int charge, TString muMinus,TString muPlus)
{
	std::cout << "[WAnalysisTool::ApplySelectionCuts] Applying selection cuts." << std::endl; 
	_dataSet = (RooDataSet*)_dataSet->reduce(Cut(cutsSig)); _dataSet->Print();
	if(charge<0){  
		_dataSet = (RooDataSet*)_dataSet->reduce(Cut(muMinus));
		std::cout << "[WAnalysisTool::ApplySelectionCuts] W- --> mu-." << std::endl;
		_dataSet->Print();
	}
	else if(charge>0) {
		_dataSet = (RooDataSet*)_dataSet->reduce(Cut(muPlus));
		std::cout << "[WAnalysisTool::ApplySelectionCuts] W+ --> mu+." << std::endl;
		_dataSet->Print();
	}
	else {
		std::cout << "[WAnalysisTool::ApplySelectionCuts] No cut on charge applied." << std::endl;
		_dataSet->Print();
	}
}

void WAnalysisTool::SelectPtEtaCentrality(double ptLow, double ptUpp, double etaLow, double etaUpp, double centralityLow, double centralityUpp,
					TString sPt , TString sEta , TString sCentrality ){
	std::cout << "[WAnalysisTool::SelectPtEtaCentrality]" << std::endl;
	TString cut = sPt;           cut+=">";  cut+=ptLow;
	cut+="&&"; cut+=sPt;         cut+="<";  cut+=ptUpp;
	cut+="&&", cut+=sEta;        cut+=">";  cut+=etaLow;
	cut+="&&"; cut+=sEta; 	     cut+="<";  cut+=etaUpp;
	cut+="&&"; cut+=sCentrality; cut+=">";  cut+=centralityLow;
	cut+="&&"; cut+=sCentrality; cut+="<";  cut+=centralityUpp;	
	std::cout << "[WAnalysisTool::SelectPtEtaCentrality] Pt:Eta:Centrality cuts: " << cut << std::endl;
	_dataSet = (RooDataSet*)_dataSet->reduce(Cut(cut));
	std::cout << "pt: " << ptLow << " - " << ptUpp << "[GeV] " << std::endl;
	std::cout << "eta: " << etaLow << " - " << etaUpp << std::endl;
	std::cout << "centrality: " << centralityLow*100.0 << " - " << centralityUpp*100 << "[%]" << std::endl;
	_dataSet->Print();
}

float WAnalysisTool::numSignalCandidates(){

	return (float)_dataSet->numEntries();
}

float WAnalysisTool::getNumCorrectedSignalCandidates(){
    if(m_dataSetIsCorrected == false){
        std::cout << "[WAnalysisTool::getNumCorrectedSignalCandidates()] ERROR: Dataset is not corrected. Must first call WAnalysisTool::correctSignalCandidates()." << std::endl;
        return -9999.0;
    }
    return m_correctedSum;
}

float WAnalysisTool::getStatErrorCorrectedSignalCandidates(){
    if(m_dataSetIsCorrected == false){
        std::cout << "[WAnalysisTool::getStatErrorCorrectedSignalCandidates()] ERROR: Dataset is not corrected. Must first call WAnalysisTool::correctSignalCandidates()." << std::endl;
        return -9999.0;
    }
        return m_correctedSumStatError;
}   

void WAnalysisTool::setZBackground(float NbkgA, float NbkgB){
    if(!m_doBkgSubtraction) m_doBkgSubtraction = true;
    std::cout << "[WAnalysisTool::setZBackground] Period A: " << NbkgA << ", Period B: " << NbkgB << std::endl;
    m_NZBkg = NbkgA+NbkgB;
    m_NBkgA+=NbkgA;
    m_NBkgB+=NbkgB;
}
void WAnalysisTool::setTauBackground(float NbkgA, float NbkgB){
    if(!m_doBkgSubtraction) m_doBkgSubtraction = true;
    std::cout << "[WAnalysisTool::setTauBackground] Period A: " << NbkgA << ", Period B: " << NbkgB << std::endl;
    m_NTauBkg = NbkgA+NbkgB;
    m_NBkgA+=NbkgA;
    m_NBkgB+=NbkgB;
}

float WAnalysisTool::getBkgSubtractedSignalCandidates(){
    if(m_doBkgSubtraction==false) {
        std::cout << "[WAnalysisTool::getBkgSubtractedSignalCandidates] WARNING: Background = 0.0. Please set it by calling setEWBackground(float NbkgA, float NbkgB)." << std::endl;
    }
    return numSignalCandidates()-m_NBkgA-m_NBkgB;
}

void WAnalysisTool::correctSignalCandidates(float CwA, float CwB){
    std::cout << "[WAnalysisTool::correctSignalCandidates] Applying Cw factors to signal yields." << std::endl;

    m_dataSetIsCorrected = true;
    float correctedSum=0.0;
    float isPeriodA = 0.0, isPeriodB = 0.0;
    float numEntriesA=0.0, numEntriesB=0.0;
    for(int i=0; i<_dataSet->numEntries(); ++i){
        isPeriodA=_dataSet->get(i)->getRealValue("isPeriodA");
        isPeriodB=_dataSet->get(i)->getRealValue("isPeriodB");
        if(isPeriodA==1.0) {
            correctedSum+=1.0/CwA;
            numEntriesA+=1.0;
        }
        else if(isPeriodB==1.0){
            correctedSum+=1.0/CwB;
            numEntriesB+=1.0;
        }
    }//i
    m_correctedSum = correctedSum;
    //std::cout << "[WAnalysisTool::correctSignalCandidates] DEBUG: " << correctedSum << " =? " << numEntriesA/CwA + numEntriesB/CwB << std::endl;
    m_correctedSumStatError = TMath::Sqrt(numEntriesA)/CwA + TMath::Sqrt(numEntriesB)/CwB;
    if(m_doBkgSubtraction) {
        numEntriesA -= m_NBkgA;
        numEntriesB -= m_NBkgB;
        m_correctedSum = numEntriesA/CwA + numEntriesB/CwB;
    }
    std::cout << "[WAnalysisTool::correctSignalCandidates] Period A: B/4 correction, NW = " << numEntriesA 
        << ", " << " CW = " << CwA << ", " << "After correction, NW/CW = " << numEntriesA/CwA << std::endl;
    std::cout << "[WAnalysisTool::correctSignalCandidates] Period B: B/4 correction, NW = " << numEntriesB 
        << ", " << " CW = " << CwB << ", " << "After correction, NW/CW = " << numEntriesB/CwB << std::endl;
}


TH1F* WAnalysisTool::GetSummedHistogram(){
	std::cout << "[WAnalysisTool::GetSummedHistogram] Retrieving summed histogram." << std::endl;
	return (TH1F*)_hsDataSetSum->GetStack()->Last();
	
}

TH1F* WAnalysisTool::GetHistogram(){
	std::cout << "[WAnalysisTool::GetHistogram] Retrieving histogram." << std::endl;
	return _hDataSet;
	
}

RooDataSet* WAnalysisTool::GetDataSet(){
	std::cout << "[WAnalysisTool::GetDataSet] Retrieving dataset." << std::endl;
	return _dataSet;
}

WAnalysisTool::~WAnalysisTool(){
	std::cout << "[WAnalysisTool::~WAnalysisTool] Dtor. " << std::endl;	
	delete _hsDataSetSum;
	delete w;
	//delete _hDataSet;
}
