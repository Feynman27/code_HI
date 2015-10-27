#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TList.h"
#include "TAxis.h"

#include <iostream>
#include <vector>
#include <cmath>

class WBookKeeper

{
	public: 
		
		WBookKeeper( const char* name, std::vector<double>& ptBins, std::vector<double>& etaBins, std::vector<double>& centralityBins) 
		 		: _name(name), _nPtBins(ptBins.size()-1), _nEtaBins(etaBins.size()-1),_nCentralityBins(centralityBins.size()-1) {

			std::cout << "[WBookKeeper::WBookKeeper] Constructing WBookKeeper object." << std::endl;
			// Containers for yields as function of muon pt
			float x[_nPtBins] ;
			float xErrUpp[_nPtBins] ;
			float xErrLow[_nPtBins] ;
			float y[_nPtBins] ;

			// Store TGraph of pt, eta, and centrality bins
			for(int ipt=0; ipt<_nPtBins; ++ipt){
				float binHalfW = (ptBins.at(ipt+1)-ptBins.at(ipt))/2.0;
				xErrUpp[ipt] = binHalfW;
				xErrLow[ipt] = binHalfW;
				float centralPoint = ptBins.at(ipt)+binHalfW;
				x[ipt] = centralPoint;
				y[ipt] = 0.0;
				_pt.SetPoint(ipt, centralPoint,0.0);
				_pt.SetPointError(ipt,binHalfW,binHalfW,0.0,0.0);
			}//ipt
			for(int ieta=0; ieta<_nEtaBins; ++ieta){
				float binHalfW = (etaBins.at(ieta+1)-etaBins.at(ieta))/2.0;
				float centralPoint = etaBins.at(ieta)+binHalfW;
				_eta.SetPoint(ieta, centralPoint,0.0);
				_eta.SetPointError(ieta,binHalfW,binHalfW,0.0,0.0);

			
				for(int icent=0; icent<_nCentralityBins; ++icent){

					float binHalfW = (centralityBins.at(icent+1)-centralityBins.at(icent))/2.0;
					float centralPoint = centralityBins.at(icent)+binHalfW;
					_centrality.SetPoint(icent, centralPoint,0.0);
					_centrality.SetPointError(icent,binHalfW,binHalfW,0.0,0.0);

					// Add TGraph objects to TList
					// Yields as function of pt per eta/centrality bin
					_sigPlus.Add( new TGraphAsymmErrors(_nPtBins,x,y,xErrLow,xErrUpp) );
					_sigMinus.Add( new TGraphAsymmErrors(_nPtBins,x,y,xErrLow,xErrUpp) );
					_sig.Add( new TGraphAsymmErrors(_nPtBins,x,y,xErrLow,xErrUpp) );

				}//icent

			}//ieta

		} //ctor 
		
		~WBookKeeper(){

			std::cout << "[WBookKeeper::~WBookKeeper] Dtor." << std::endl;
			for(int i = 0; i<_sig.GetEntries(); ++i){
	
				if(_sig.At(i)) delete _sig.At(i);
				if(_sigPlus.At(i)) delete _sigPlus.At(i);
				if(_sigMinus.At(i)) delete _sigMinus.At(i);
			}//i
		}//dtor

		// Method to set TGraph point to yield
		// in this eta/centrality bin
		void setSig(const int charge, const int ipt, const int ieta, const int icent, const double yield, const double errLow, const double errUpp){
			// Retrieve index for this eta/centrality bin
			const int index = _nCentralityBins*ieta+icent;
			// Set pt point ipt to current yield
			if(charge>0) {
				setY( (TGraphAsymmErrors*)_sigPlus.At(index), ipt, yield, errLow, errUpp);	
				((TGraphAsymmErrors*)_sigPlus.At(index))->GetXaxis()->SetTitle("pt");
			}
			else if(charge<0) {
				setY( (TGraphAsymmErrors*)_sigMinus.At(index), ipt, yield, errLow, errUpp);	
				((TGraphAsymmErrors*)_sigMinus.At(index))->GetXaxis()->SetTitle("pt");
			}
			else {
				setY( (TGraphAsymmErrors*)_sig.At(index), ipt, yield, errLow, errUpp);	
				((TGraphAsymmErrors*)_sig.At(index))->GetXaxis()->SetTitle("pt");
			}
		}		

		// Method for writing TGraphs to outfile after 
		// obtaining yields
		void write(TFile& f, const int charge){
			TDirectory *dir = gDirectory;
			f.cd();
			_pt.Write( TString(_name)+="_pt");
			_eta.Write( TString(_name)+="_eta");
			_centrality.Write( TString(_name)+="_centrality");
			for(int ieta = 0; ieta<_nEtaBins; ++ieta){
			  for(int icent = 0; icent<_nCentralityBins; ++icent){
	                    int index = ieta*_nCentralityBins+icent;
			    TString baseName = TString(_name);
		            baseName+="_eta"; baseName+=ieta;
			    baseName+="_centrality"; baseName+=icent;	
				
			    if(charge>0) {
				baseName+="_chargePlus";
				((TGraphAsymmErrors*)_sigPlus.At(index))->Write(baseName);
			    }	
			    else if(charge<0){ 
				baseName+="_chargeMinus";
				((TGraphAsymmErrors*)_sigMinus.At(index))->Write(baseName);
			    }
			    else {
				((TGraphAsymmErrors*)_sig.At(index))->Write(baseName);
			    }
			  }//icent
			}//eta	

			gDirectory = dir;
		}

	private:
		
		void setY(TGraphAsymmErrors* gr, const int ipt, const double yield, const double errLow, const double errUpp){
			
			double x,y;
			//Fill buffer of x and y with values from TGraph	
			gr->GetPoint(ipt,x,y);
			gr->SetPoint(ipt,x,yield);
			gr->SetPointEYlow(ipt,errLow);
			gr->SetPointEYhigh(ipt,errUpp);
		}		

		const char* _name;
		TGraphAsymmErrors _pt;
		TGraphAsymmErrors _eta;
		TGraphAsymmErrors _centrality;
		TList _sig;
		TList _sigPlus;
		TList _sigMinus;
		int _charge;
		int _nPtBins;
		int _nEtaBins ;
		int _nCentralityBins ;
};
