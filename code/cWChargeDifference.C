//This macro calculates the difference in the Cw charge-separated distribution
//as a function of Npart over 10 eta bins 

{
	TFile* fCw = new TFile("correctionFactorsW.a1_a2_fixed.01.17.2013.root","UPDATE");	
	//TGraphErrors *grPlus = fCw->Get("gr102");
	int nEtaBins = 10; int nCentralityBins = 6;
	double sumOfSqrs = 0;

	for(int ieta = 0; ieta<nEtaBins; ieta++){

		TString sGrPlus = "grCent";
		TString sGrMinus = "grCent";
		sGrPlus += 102; sGrPlus+="eta"; sGrPlus+=ieta;
		sGrMinus += 103; sGrMinus+="eta"; sGrMinus+=ieta;
		std::cout << sGrPlus << std::endl;
		std::cout << sGrMinus << std::endl;
		TGraphErrors* grPlus = fCw->Get(sGrPlus);
		TGraphErrors* grMinus = fCw->Get(sGrMinus);
		double* yPlus = grPlus->GetY();
		double* yMinus = grMinus->GetY();
		for(int icent = 0; icent<nCentralityBins; icent++){
			
			double yTempPlus = yPlus[icent]; 
			double yTempMinus = yMinus[icent]; 
			std::cout << "Difference in Cw at Eta Bin " << ieta << " and Centrality Bin " << icent << " = " << fabs(yTempPlus-yTempMinus) << std::endl;
		}
	}
}
