void plots() {
    gSystem->Load("libTPORec.so");
//    TFile f("Batch-TPORecevent_200035_nutauCC.root");
    TFile f("Batch-TPORecevent_200026_nueCC.root");
 
 #if 0
    TCanvas *c1 = new TCanvas("c1", "dE/dx", 800, 800);
    c1->Divide(2, 1);
    c1->cd(1);
    TH1D* h_mu_dedx;
    f.GetObject("mu_dedx",h_mu_dedx);
    h_mu_dedx->Draw();
    c1->cd(2);
    TH1D* h_tau_dedx;
    f.GetObject("tau_dedx",h_tau_dedx);
    h_tau_dedx->Draw();
    c1->Update();
#endif

    TTree *tt;f.GetObject("RecoEvent",tt);

    tt->Draw("fPOFullEvent->TotalEvis()");
    tt->Draw("fPOFullEvent->TotalET()");
    tt->Draw("fPOFullEvent->TotalET():fPOFullEvent->TotalEvis()","","COLZ");
//    tt->Draw("fPOFullEvent->TotalET():fPOFullEvent->TotalEvis()","fPOFullEvent->TotalEvis()<1500","COLZ");
//    tt->Draw("primary_n_charged");
//    tt->Draw("nhits_tracker_first");
//    tt->Draw("nhits_tracker_first:primary_n_charged", "", "COLZ");
    tt->Draw("nhits_tau");
    tt->Draw("fPOFullEvent->TotalEvis()", "nhits_tau>4");
//    tt->Draw("fPOFullEvent->TotalET()", "nhits_tau>4");

    f.Close();
} 
