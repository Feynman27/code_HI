double eta_min[200] = {0};
double eta_max[200] = {0};
double cent_min[200] = {0};
int    eta_bin[200] = {0};
double cent_max[200] = {0};
int    cent_bin[200] = {0};
//index 0 if not binned for charge
int    charge[200] = {0};
//double trig_eff[222] = {0};
std::vector <double> trig_eff (54);
//std::vector <double> trig_eff (228);
//std::vector <double> trig_eff (456);
//double trig_err[222] = {0};
std::vector <double> trig_err (54);
//std::vector <double> trig_err (228);
//std::vector <double> trig_err (456);

//void readInputFile(TString sFileIn = "TrigEff_Sep22.txt"){
//void readInputFile(TString sFileIn = "TrigEff.Znote.txt"){
//void readInputFile(TString sFileIn = "triggerEffZNote_v03.txt"){
void readInputFile(int nCentralityBins,int nEtaBins, bool doMirrorEta=true, TString sFileIn = "triggerEffZNote_v04.txt"){
    
	int totalBins =  nCentralityBins*nEtaBins;
	//ifstream s(sFileIn);
	std::ifstream s( sFileIn , std::ifstream::in );

	if(!s.is_open()){
		std::cout << "While opening a file an error was encountered" << std::endl;
		exit(1);
	} else {
		std::cout << "File " << sFileIn << " is successfully opened" << std::endl;
	}

//	unsigned int i = 0;
//	while(!s.eof()) 
    if(s.is_open()){
        for (int i = 0; i < trig_eff.size(); ++i) {
		    s >> trig_eff[i] >> trig_err[i];
            std::cout << "index: " << i << " trigger eff: " << trig_eff[i] << " +- " << trig_err[i] << std::endl;
        }
        if( (trig_eff.size()==totalBins) || (trig_err.size()==totalBins) || !doMirrorEta){
            std:: cout << "Trigger container filled successfully. " << std::endl;
        }

        else{
            std::cout << "ERROR: trigger container size does not match number of total bins used. Aborting." << std::endl;
            exit(0);
        }
        s.close();
        //exit(0);
		//s >> pt_min[i] >> pt_max[i] >> eta_min[i] >> eta_max[i] >> cent_min[i] >> cent_max[i] >> trig_eff[i] >> trig_err[i];
//		s >> charge[i] >> eta_bin[i] >> cent_bin[i] >> trig_eff[i] >> trig_err[i];
        ///reads in trig eff indexed by ieta*nCentrality + icent
//		s >> trig_eff[i] >> trig_err[i];
//		std::cout << i << " " << charge[i] << " " << eta_bin[i] << " " << cent_bin[i] << " " << trig_eff[i] << " " << trig_err[i] << std::endl;
		//s >> charge[i] >> eta_min[i] >> eta_max[i] >> cent_bin[i] >> trig_eff[i] >> trig_err[i];
//		++i;
	}
    else std::cout << "unable to open trigger eff file" << std::endl;
}

//bin 0 is over whole range
//double trigEfficiency(int chargeBin, int etaBin, int centralityBin){
double trigEfficiency(int index=-1){

    std::cout << "Index: " << index << std::endl;
    if(index>trig_eff.size()){
        std::cout << "ERROR: Attempting to access element outside of trigger container size. Aborting." << std::endl;
        exit(0);
    }


//	std::cout << chargeBin << ":" << etaBin << ":" << centralityBin << std::endl;

/*	int index = -1;

	for(int i=0;i<200;++i){

		if( chargeBin == charge[i] && etaBin == eta_bin[i] && centralityBin == cent_bin[i] )
		{
//			std::cout << i << " " << charge[i] << " " << eta_bin[i] << " " << cent_bin[i] << " " << trig_eff[i] << " " << trig_err[i] << std::endl;
			index = i;
			break;
		}
	}
*/
	if(index==-1){
		return -999.0;
	} else {
		std::cout << "trigger eff at index : " << index << " = " << trig_eff[index] << std::endl;
		return trig_eff[index];
	}

}

//double trigEfficiencyErr(int chargeBin, int etaBin, int centralityBin){
double trigEfficiencyErr(int index = -1){

    std::cout << "Index: " << index << std::endl;
    if(index>trig_err.size()){
        std::cout << "ERROR: Attempting to access element outside of trigger container size. Aborting." << std::endl;
        exit(0);
    }

/*	int index = -1;



	for(int i=0;i<200;i++){

		if( chargeBin == charge[i] && etaBin == eta_bin[i] && centralityBin == cent_bin[i])
		{
			//std::cout << i << " " << charge[i] << " " << eta_bin[i] << " " << cent_bin[i] << " " << trig_eff[i] << " " << trig_err[i] << std::endl;
			index = i;
			break;
		}

	}
*/
	if(index==-1){
		return -999.0;
	} else {
		return trig_err[index];
	}

}


