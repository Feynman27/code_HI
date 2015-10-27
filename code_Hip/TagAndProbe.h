#include "EfficiencyCorrection.h"

#ifndef TagAndProbe_h
#define TagAndProbe_h
class TagAndProbe:public EfficiencyCorrection{
    public:
        TagAndProbe(std::vector<double>& vPt, std::vector<double>& vEta, std::vector<double>& vCent, std::string beam, std::string config, std::string charge, std::string tpType) ;
        void fillTPHistos(double pt, double eta, double centrality, double charge, std::string type);
        ~TagAndProbe(){};
};
#endif
