#include "TH2F"

void toy1()
  
{
  float nu_range = 150; // GeV
  float smear = 50; // GeV
  int stat = 100000;
  
  double px_n, py_n, pt_n, phi;
  double px_r, py_r, pt_r;
  double px_s, py_s, pt_s;
  
  
  TH2F *cor = new TH2F("cor",";pt_nu;pt_smear",200,0,200,200,0,200);
  
  for (int i=0; i<stat; i++)
    {
      pt_n = gRandom->Rndm()*nu_range;
      phi  = gRandom->Rndm()*2.*3.14159;
      px_n = pt_n*cos(phi);
      py_n = pt_n*sin(phi);

      px_r = gRandom->Gaus()*smear;
      py_r = gRandom->Gaus()*smear;


      px_s = px_n+px_r;
      py_s = py_n+py_r;

      pt_s = pow(pow(px_s,2)+pow(py_s,2),0.5); 

      cor->Fill(pt_n,pt_s);
    }


  gStyle->SetPalette(1);
  cor->Draw("colz");
  TProfile *cor_p = (TProfile*) cor->ProfileX("cor_p");
  cor_p->SetMarkerSize(0.5);
  cor_p->SetMarkerStyle(8);
  //  cor_p->Rebin(5);
  cor_p->Draw("same");

}
