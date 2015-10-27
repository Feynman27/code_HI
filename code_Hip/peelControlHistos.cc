
void peelControlHistos(){

    TFile* _fIn = TFile::Open("$HOME/scratch/Data2013Hip/skimmed.07.22.2014.root");
    //TFile* _fIn = TFile::Open("$HOME/scratch/Data2013Hip/skimmed.mergedPeriods.06.28.2014.root");
    TFile* _fOut = new TFile("skimmedCtrlHistos.root","recreate");
    
    _fIn->cd();
    TH1I* h0 = (TH1I*)gDirectory->Get("hmu_medium");
    TH1I* h1 = (TH1I*)gDirectory->Get("hmu_author");
    TH1I* h2 = (TH1I*)gDirectory->Get("hmu_nPixHits");
    TH1I* h3 = (TH1I*)gDirectory->Get("hmu_nBLHits");
    TH1I* h4 = (TH1I*)gDirectory->Get("hmu_nSCTHits");
    TH1I* h5 = (TH1I*)gDirectory->Get("hmu_nPixSCTHoles");
    TH1I* h6 = (TH1I*)gDirectory->Get("hmu_isCombinedMuon");
    TH1F* h7 = (TH1F*)gDirectory->Get("hmu_EF_mu8_dR");

    TH1F* h8 = (TH1F*)gDirectory->Get("hmu_momentumBalanceSignificance");
    TH1F* h9 = (TH1F*)gDirectory->Get("hmu_id_d0_exPV");
    TH1F* h10 = (TH1F*)gDirectory->Get("hmu_id_z0_exPVsinTheta");
    TH1F* h11 = (TH1F*)gDirectory->Get("hmu_matchchi2Overndof");
    TH1F* h12 = (TH1F*)gDirectory->Get("hmu_scatteringCurvatureSignificance");
    TH1F* h13 = (TH1F*)gDirectory->Get("hmu_z0_vx_z_diff");

    TCanvas* c0 = new TCanvas("c_mu_medium","c_mu_medium",600,600);
    TCanvas* c1 = new TCanvas("c_mu_author","c_mu_author",600,600);
    TCanvas* c2 = new TCanvas("c_mu_nPixHits","c_mu_nPixHits",600,600);
    TCanvas* c3 = new TCanvas("c_mu_nBLHits","c_mu_nBLHits",600,600);
    TCanvas* c4 = new TCanvas("c_mu_nSCTHits","c_mu_nSCTHits",600,600);
    TCanvas* c5 = new TCanvas("c_mu_nPixSCTHoles","c_mu_nPixSCTHoles",600,600);
    TCanvas* c6 = new TCanvas("c_mu_isCombinedMuon","c_mu_isCombinedMuon",600,600);
    TCanvas* c7 = new TCanvas("c_mu_EF_mu8_dR","c_mu_EF_mu8_dR",600,600);
    c7->SetLogy(1);
    TCanvas* c8 = new TCanvas("c_mu_momentumBalanceSignificance","c_mu_momentumBalanceSignificance",600,600);
    TCanvas* c9 = new TCanvas("c_mu_id_d0_exPV","c_mu_id_d0_exPV",600,600);
    c9->SetLogy(1);
    TCanvas* c10 = new TCanvas("c_mu_id_z0_exPVsinTheta","c_mu_id_z0_exPVsinTheta",600,600);
    c10->SetLogy(1);
    TCanvas* c11 = new TCanvas("c_mu_matchchi2Overndof","c_mu_matchchi2Overndof",600,600);
    c11->SetLogy(1);
    TCanvas* c12 = new TCanvas("c_mu_scatteringCurvatureSignificance","c_mu_scatteringCurvatureSignificance",600,600);


    h0->GetXaxis()->SetRangeUser(0,2);
    h1->GetXaxis()->SetRangeUser(0,30);
    h2->GetXaxis()->SetRangeUser(0,10);
    h3->GetXaxis()->SetRangeUser(0,10);
    h4->GetXaxis()->SetRangeUser(0,20);
    h5->GetXaxis()->SetRangeUser(0,5);
    h6->GetXaxis()->SetRangeUser(0,2);
    h7->GetXaxis()->SetRangeUser(0,1.0);
    h8->GetXaxis()->SetRangeUser(0,10);
    h9->GetXaxis()->SetRangeUser(0,10);
    h10->GetXaxis()->SetRangeUser(0,10);
    h11->GetXaxis()->SetRangeUser(0,25);
    h12->GetXaxis()->SetRangeUser(0,10);

    h7->GetYaxis()->SetRangeUser(0.1,1e6);
    h9->GetYaxis()->SetRangeUser(0.1,1e6);
    h10->GetYaxis()->SetRangeUser(0.1,1e6);
    h11->GetYaxis()->SetRangeUser(0.1,1e6);

    h0->SetXTitle("medium");
    h1->SetXTitle("author");
    h2->SetXTitle("nPixHits");
    h3->SetXTitle("nBLHits");
    h4->SetXTitle("nSCTHits");
    h5->SetXTitle("nPixSCTHoles");
    h6->SetXTitle("isCombinedMuon");
    h7->SetXTitle("EF_mu8_dR");  
    h8->SetXTitle("momentumBalanceSignificance");  
    h9->SetXTitle("id_d0_exPV");  
    h10->SetXTitle("id_z0_exPVsinTheta");
    h11->SetXTitle("matchchi2Overndof");
    h12->SetXTitle("scatteringCurvatureSignificance");
    h13->SetXTitle("z0_vx_z_diff");

    h0->SetFillColor(kYellow);
    h1->SetFillColor(kYellow);
    h2->SetFillColor(kYellow);
    h3->SetFillColor(kYellow);
    h4->SetFillColor(kYellow);
    h5->SetFillColor(kYellow);
    h6->SetFillColor(kYellow);
    h7->SetFillColor(kYellow);  
    h8->SetFillColor(kYellow);  
    h9->SetFillColor(kYellow);  
    h10->SetFillColor(kYellow);
    h11->SetFillColor(kYellow);
    h12->SetFillColor(kYellow);
    h13->SetFillColor(kYellow);

    h0->SetYTitle("Muons");
    h1->SetYTitle("Muons");
    h2->SetYTitle("Muons");
    h3->SetYTitle("Muons");
    h4->SetYTitle("Muons");
    h5->SetYTitle("Muons");
    h6->SetYTitle("Muons");
    h7->SetYTitle("Muons");  
    h8->SetYTitle("Muons");  
    h9->SetYTitle("Muons");  
    h10->SetYTitle("Muons");
    h11->SetYTitle("Muons");
    h12->SetYTitle("Muons");
    h13->SetYTitle("Muons");

    c0->cd(); h0->Draw();  c0->Print("plots/hmu_medium.pdf");
    c1->cd(); h1->Draw();c1->Print("plots/hmu_author.pdf");
    c2->cd(); h2->Draw();c2->Print("plots/hmu_nPixHits.pdf");
    c3->cd(); h3->Draw();c3->Print("plots/hmu_nBLHits.pdf");
    c4->cd(); h4->Draw();c4->Print("plots/hmu_nSCTHits.pdf");
    c5->cd(); h5->Draw();c5->Print("plots/hmu_nPixSCTHoles.pdf");
    c6->cd(); h6->Draw();c6->Print("plots/hmu_isCombinedMuon.pdf");
    c7->cd(); h7->Draw();c7->Print("plots/hmu_EF_mu8_dR.pdf");
    c8->cd(); h8->Draw();c8->Print("plots/hmu_momentumBalanceSignificance.pdf");
    c9->cd(); h9->Draw();c9->Print("plots/hmu_id_d0_exPV.pdf");
    c10->cd(); h10->Draw();c10->Print("plots/hmu_id_z0_exPVsinTheta.pdf");
    c11->cd(); h11->Draw();c11->Print("plots/hmu_matchchi2Overndof.pdf");
    c12->cd(); h12->Draw();c12->Print("plots/hmu_scatteringCurvatureSignificance.pdf");

    _fOut->cd();

    c0->Write();
    c1->Write();
    c2->Write();
    c3->Write();
    c4->Write();
    c5->Write();
    c6->Write();
    c7->Write();  
    c8->Write();  
    c9->Write();  
    c10->Write();
    c11->Write();
    c12->Write();

    h0->Write();
    h1->Write();
    h2->Write();
    h3->Write();
    h4->Write();
    h5->Write();
    h6->Write();
    h7->Write();  
    h8->Write();  
    h9->Write();  
    h10->Write();
    h11->Write();
    h12->Write();
    h13->Write();
    _fOut->Close();
}
