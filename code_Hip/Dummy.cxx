#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
#include "TH1F.h"

int main(){

    std::string hnametot;
    std::stringstream shnametot;
    shnametot.str(std::string());
    shnametot << "_hPosTot" << "_pt" << 0 <<
        "_eta" << 0 <<
        "_cent" << 0;
    hnametot=shnametot.str();
    
    std::map<std::string, TH1F* >  m_hPass;
    TH1F* hTotPos = new TH1F("h","mu^{+}, total", 10,0.0,10.0);
    m_hPass.insert(std::make_pair<std::string&, TH1F*& > (hnametot,hTotPos));
}
