double eta_min[200] = {0};
double eta_max[200] = {0};
double cent_min[200] = {0};
int    eta_bin[200] = {0};
double cent_max[200] = {0};
int    cent_bin[200] = {0};
//index 0 if not binned for charge
int    charge[200] = {0};
double trig_eff[200] = {0};
double trig_err[200] = {0};

//void readInputFile(TString sFileIn = "TrigEff_Sep22.txt"){
void readInputFile(TString sFileIn = "/afs/cern.ch/user/t/tbalestr/public/17.2.2.2/code/TrigEff.Znote.txt"){
//void readInputFile(TString sFileIn = "triggerEffZNote_v02.txt"){
	
	//ifstream s(sFileIn);
	std::ifstream s( sFileIn , std::ifstream::in );

	if(!s.is_open()){
		std::cout << "While opening a file an error was encountered" << std::endl;
		exit(1);
	} else {
		std::cout << "File " << sFileIn << " is successfully opened" << std::endl;
	}

	unsigned int i = 0;
	while(!s.eof()) {
		//s >> pt_min[i] >> pt_max[i] >> eta_min[i] >> eta_max[i] >> cent_min[i] >> cent_max[i] >> trig_eff[i] >> trig_err[i];
		s >> charge[i] >> eta_bin[i] >> cent_bin[i] >> trig_eff[i] >> trig_err[i];
//		std::cout << i << " " << charge[i] << " " << eta_bin[i] << " " << cent_bin[i] << " " << trig_eff[i] << " " << trig_err[i] << std::endl;
		//s >> charge[i] >> eta_min[i] >> eta_max[i] >> cent_bin[i] >> trig_eff[i] >> trig_err[i];
		++i;
	}
}

//bin 0 is over whole range
double trigEfficiency(int chargeBin, int etaBin, int centralityBin){

//	std::cout << chargeBin << ":" << etaBin << ":" << centralityBin << std::endl;

	int index = -1;

	for(int i=0;i<200;++i){

		if( chargeBin == charge[i] && etaBin == eta_bin[i] && centralityBin == cent_bin[i] )
		{
//			std::cout << i << " " << charge[i] << " " << eta_bin[i] << " " << cent_bin[i] << " " << trig_eff[i] << " " << trig_err[i] << std::endl;
			index = i;
			break;
		}
	}

	if(index==-1){
		return -999.0;
	} else {
//		std::cout << "trigger eff: " << trig_eff[index] << std::endl;
		return trig_eff[index];
	}

}

double trigEfficiencyErr(int chargeBin, int etaBin, int centralityBin){

	int index = -1;



	for(int i=0;i<200;i++){

		if( chargeBin == charge[i] && etaBin == eta_bin[i] && centralityBin == cent_bin[i])
		{
			//std::cout << i << " " << charge[i] << " " << eta_bin[i] << " " << cent_bin[i] << " " << trig_eff[i] << " " << trig_err[i] << std::endl;
			index = i;
			break;
		}

	}

	if(index==-1){
		return 1;
	} else {
		return trig_err[index];
	}

}


