{
/*TFile *_file0 = TFile::Open("user.thomas.MuonJetsHardProbes.D3PD.193825.500MeV.v2.18Jun_2012.193825.120618210344/HISingleMuon.root");
TFile *_file1 = TFile::Open("user.thomas.MuonJetsHardProbes.D3PD.193825.600MeV.18Jun_2012.193825.120618174613/HISingleMuon.root");
TFile *_file2 = TFile::Open("user.thomas.MuonJetsHardProbes.D3PD.193825.700MeV.18Jun_2012.193825.120618174725/HISingleMuon.root");
TFile *_file3 = TFile::Open("user.thomas.MuonJetsHardProbes.D3PD.193825.800MeV.18Jun_2012.193825.120618174851/HISingleMuon.root");
TFile *_file4 = TFile::Open("user.thomas.MuonJetsHardProbes.D3PD.193825.900MeV.18Jun_2012.193825.120618175123/HISingleMuon.root");
TFile *_file5 = TFile::Open("user.thomas.MuonJetsHardProbes.D3PD.193825.1000MeV.18Jun_2012.193825.120618175234/HISingleMuon.root");
TFile *_file6 = TFile::Open("user.thomas.MuonJetsHardProbes.D3PD.193825.1500MeV.18Jun_2012.193825.120618175536/HISingleMuon.root");
TFile *_file7 = TFile::Open("user.thomas.MuonJetsHardProbes.D3PD.193825.2000MeV.18Jun_2012.193825.120618175708/HISingleMuon.root");
TFile *_file8 = TFile::Open("user.thomas.MuonJetsHardProbes.D3PD.193825.2500MeV.18Jun_2012.193825.120618175831/HISingleMuon.root");
TFile *_file9 = TFile::Open("user.thomas.MuonJetsHardProbes.D3PD.193825.3000MeV.18Jun_2012.193825.120618175942/HISingleMuon.root");
*/
TFile *_file0 = TFile::Open("/tmp/tbalestr/500MeV/HISingleMuon.root");
TFile *_file1 = TFile::Open("/tmp/tbalestr/600MeV/HISingleMuon.root");
TFile *_file2 = TFile::Open("/tmp/tbalestr/700MeV/HISingleMuon.root");
TFile *_file3 = TFile::Open("/tmp/tbalestr/800MeV/HISingleMuon.root");
TFile *_file4 = TFile::Open("/tmp/tbalestr/900MeV/HISingleMuon.root");
TFile *_file5 = TFile::Open("/tmp/tbalestr/1000MeV/HISingleMuon.root");
TFile *_file6 = TFile::Open("/tmp/tbalestr/1500MeV/HISingleMuon.root");
TFile *_file7 = TFile::Open("/tmp/tbalestr/2000MeV/HISingleMuon.root");
TFile *_file8 = TFile::Open("/tmp/tbalestr/2500MeV/HISingleMuon.root");
TFile *_file9 = TFile::Open("/tmp/tbalestr/3000MeV/HISingleMuon.root");

TCanvas c0 = TCanvas("c0","c0",600,600);
#ifdef __CINT__
  gROOT->LoadMacro("AtlasUtils.C");
#endif


 _file0->cd();
TH1F* h0 = new TH1F("h0","h0",90,0.,180.);
tree->Draw("nu_pt>>h0","val>11&&centrality<=0.8&&abs(scat)<4.0&&abs(eLoss)<0.5&&hasHiPtMuon==1&&pt<90.");
std::cout << "Mean: " << h0->GetMean() << std::endl;
//h0->GetYaxis()->SetRangeUser(0.0,100);
h0->SetXTitle("#slash{p_{T}}");
myText(0.55,0.75,1,"PYTHIA+HIJING");
myText(0.55,0.55,1,"p_{T}^{trk}>500MeV");
c0.Print("MCmissPt500_21Jun.pdf");

//h0->SetLineColor(kRed+5);
//h0->Draw();

 _file1->cd();
TH1F* h1 = new TH1F("h1","h1",90,0.,180.);
tree->Draw("nu_pt>>h1","val>11&&centrality<=0.8&&abs(scat)<4.0&&abs(eLoss)<0.5&&hasHiPtMuon==1&&pt<90.","same");
h1->SetXTitle("#slash{p_{T}}");
myText(0.55,0.75,1,"PYTHIA+HIJING");
myText(0.55,0.55,1,"p_{T}^{trk}>600MeV");
c0.Print("MCmissPt600_21Jun.pdf") ;
//h1->SetLineColor(kRed);
//h1->Draw("same");

_file2->cd();
TH1F* h2 = new TH1F("h2","h2",90,0.,180.);
tree->Draw("nu_pt>>h2","val>11&&centrality<=0.8&&abs(scat)<4.0&&abs(eLoss)<0.5&&hasHiPtMuon==1&&pt<90.","same");
h2->SetXTitle("#slash{p_{T}}");
myText(0.55,0.75,1,"PYTHIA+HIJING");
myText(0.55,0.55,1,"p_{T}^{trk}>700MeV");
c0.Print("MCmissPt700_21Jun.pdf");
//h2->SetLineColor(kBlue);
//h2->Draw("same");

_file3->cd();
TH1F* h3 = new TH1F("h3","h3",90,0.,180.);
tree->Draw("nu_pt>>h3","val>11&&centrality<=0.8&&abs(scat)<4.0&&abs(eLoss)<0.5&&hasHiPtMuon==1&&pt<90.","same");
h3->SetXTitle("#slash{p_{T}}");
myText(0.55,0.75,1,"PYTHIA+HIJING");
myText(0.55,0.55,1,"p_{T}^{trk}>800MeV");
c0.Print("MCmissPt800_21Jun.pdf");
//h3->SetLineColor(kGreen);
//h3->Draw("same");

_file4->cd();
TH1F* h4 = new TH1F("h4","h4",90,0.,180.);
tree->Draw("nu_pt>>h4","val>11&&centrality<=0.8&&abs(scat)<4.0&&abs(eLoss)<0.5&&hasHiPtMuon==1&&pt<90.","same");
h4->SetXTitle("#slash{p_{T}}");
myText(0.55,0.75,1,"PYTHIA+HIJING");
myText(0.55,0.55,1,"p_{T}^{trk}>900MeV");
c0.Print("MCmissPt900_21Jun.pdf") ;
//h4->SetLineColor(kOrange);
//h4->Draw("same");

 _file5->cd();
TH1F* h5 = new TH1F("h5","h5",90,0.,180.);
tree->Draw("nu_pt>>h5","val>11&&centrality<=0.8&&abs(scat)<4.0&&abs(eLoss)<0.5&&hasHiPtMuon==1&&pt<90.","same");
h5->SetXTitle("#slash{p_{T}}");
myText(0.55,0.75,1,"PYTHIA+HIJING");
myText(0.55,0.55,1,"p_{T}^{trk}>1000MeV");
c0.Print("MCmissPt1000_21Jun.pdf") ;
//h5->SetLineColor(kMagenta);
//h5->Draw("same");

_file6->cd();
TH1F* h6 = new TH1F("h6","h6",90,0.,180.);
tree->Draw("nu_pt>>h6","val>11&&centrality<=0.8&&abs(scat)<4.0&&abs(eLoss)<0.5&&hasHiPtMuon==1&&pt<90.","same");
h6->SetXTitle("#slash{p_{T}}");
myText(0.55,0.75,1,"PYTHIA+HIJING");
myText(0.55,0.55,1,"p_{T}^{trk}>1500MeV");
c0.Print("MCmissPt1500_21Jun.pdf") ;
//h6->SetLineColor(kBlue-9);
//h6->Draw("same");

_file7->cd();
TH1F* h7 = new TH1F("h7","h7",90,0.,180.);
tree->Draw("nu_pt>>h7","val>11&&centrality<=0.8&&abs(scat)<4.0&&abs(eLoss)<0.5&&hasHiPtMuon==1&&pt<90.","same");
h7->SetXTitle("#slash{p_{T}}");
myText(0.55,0.75,1,"PYTHIA+HIJING");
myText(0.55,0.55,1,"p_{T}^{trk}>2000MeV");
c0.Print("MCmissPt2000_21Jun.pdf") ;
//h7->SetLineColor(kCyan);
//h7->Draw("same");

_file8->cd();
TH1F* h8 = new TH1F("h8","h8",90,0.,180.);
tree->Draw("nu_pt>>h8","val>11&&centrality<=0.8&&abs(scat)<4.0&&abs(eLoss)<0.5&&hasHiPtMuon==1&&pt<90.","same");
h8->SetXTitle("#slash{p_{T}}");
myText(0.55,0.75,1,"PYTHIA+HIJING");
myText(0.55,0.55,1,"p_{T}^{trk}>2500MeV");
c0.Print("MCmissPt2500_21Jun.pdf") ;
//h8->SetLineColor(kYellow+3);
//h8->Draw("same");

_file9->cd();
TH1F* h9 = new TH1F("h9","h9",90,0.,180.);
tree->Draw("nu_pt>>h9","val>11&&centrality<=0.8&&abs(scat)<4.0&&abs(eLoss)<0.5&&hasHiPtMuon==1&&pt<90.","same");
h9->SetXTitle("#slash{p_{T}}");
myText(0.55,0.75,1,"PYTHIA+HIJING");
myText(0.55,0.55,1,"p_{T}^{trk}>3000MeV");
c0.Print("MCmissPt3000_21Jun.pdf") ;
//h9->SetLineColor(kBlack);
//h9->Draw("same");

/*TLegend* leg0 = new TLegend(0.55, 0.5, 0.9, 0.9);
leg0->SetTextFont(gStyle->GetTextFont());
leg0->SetTextSize(gStyle->GetTextSize());
leg0->SetBorderSize(0);
leg0->SetFillColor(0);
  leg0->AddEntry(h0, "500MeV", "l");
  leg0->AddEntry(h1, "600MeV", "l");
  leg0->AddEntry(h0, "700MeV", "l");
  leg0->AddEntry(h1, "800MeV", "l");
  leg0->AddEntry(h0, "900MeV", "l");
  leg0->AddEntry(h1, "1000MeV", "l");
  leg0->AddEntry(h0, "1500MeV", "l");
  leg0->AddEntry(h1, "2000MeV", "l");
  leg0->AddEntry(h0, "2500MeV", "l");
  leg0->AddEntry(h1, "3000MeV", "l");
  //leg0->Draw();
//myText(0.4,0.75,1,"#int L #approx 5 #mub^{-1}");
*/
}
