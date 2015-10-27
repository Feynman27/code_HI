#include "TagAndProbe.h"

TagAndProbe::TagAndProbe( std::vector<double>& vPt, std::vector<double>& vEta, std::vector<double>& vCent, std::string beam, std::string config, std::string charge,
                            std::string tpType) : EfficiencyCorrection(vPt,vEta,vCent,beam,config,charge,tpType){

    if(strcmp(beam.c_str(),"data")==0) std::cout << "[TagAndProbe::TagAndProbe] Running T&P analysis with data." << std::endl;
    else if(strcmp(beam.c_str(),"pp")==0) std::cout << "[TagAndProbe::TagAndProbe] Running T&P analysis with p+p MC." << std::endl;
    else if(strcmp(beam.c_str(),"pn")==0) std::cout << "[TagAndProbe::TagAndProbe] Running T&P analysis with p+n MC." << std::endl;
    else if(strcmp(beam.c_str(),"np")==0) std::cout << "[TagAndProbe::TagAndProbe] Running T&P analysis with n+p MC." << std::endl;
    else {
        std::cout << "[TagAndProbe::TagAndProbe] ERROR: Cannot identify beam type. Please select data, pp, np, or pn. " << std::endl;
        return;
    }

    std::cout << "[TagAndProbe::TagAndProbe] Performing T&P analysis for " << tpType << std::endl;

}

void TagAndProbe::fillTPHistos(double pt, double eta, double centrality, double charge, std::string type){
    fillHistos(pt,eta,centrality,charge,type);
}
