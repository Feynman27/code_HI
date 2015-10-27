#include "TVectorF.h"
#include "TGraphErrors.h"

int __Npoints = 6; //number of centrality bins
TVectorF etap = TVectorF(__Npoints);
TVectorF etape = TVectorF(__Npoints);

void setEtaPoints() {
  etap[0] = 2.5;
  etape[0] = 2.5;
  etap[1] = 7.5;
  etape[1] = 2.5;
  etap[2] = 12.5;
  etape[2] = 2.5;
  etap[3] = 17.5;
  etape[3] = 2.5;
  etap[4] = 30.0;
  etape[4] = 10.0;
  etap[5] = 60.0;
  etape[5] = 20.0;

}

TGraphErrors* getNColl() {
  setEtaPoints();
  TVectorF yp(__Npoints);
  TVectorF ype(__Npoints);
  yp[0] = 1683.3 ; //0-5 
  yp[1] = 1318.0 ; //5-10
  yp[2] = 1035.4 ;//10-15
  yp[3] = 811.2 ;//15-20
  yp[4] = 440.6 ;//20-40
  yp[5] = 77.8 ; //40-80
 
  ype[0] = 7.7;
  ype[1] = 7.5;
  ype[2] = 7.4;
  ype[3] = 7.4;
  ype[4] = 7.3;
  ype[5] = 11.6;
  TGraphErrors* func = new TGraphErrors(etap,yp,etape,ype);
  return func;

}
