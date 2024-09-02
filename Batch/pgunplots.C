TTree* pgtree;

void plotenergy() {
    TH1F *h1 = new TH1F("h1", "Energy Distribution", 100, 0, 3000); // 100 bins between 0 and 3000 GeV
    TH1F *h2 = new TH1F("h2", "Energy Distribution", 100, 0, 3000); // 100 bins between 0 and 3000 GeV

    // Fill the histograms with the desired selections
    pgtree->Draw("m_energy>>h1", "abs(m_pdg_id)==11");
    pgtree->Draw("m_energy>>h2", "abs(m_pdg_id)==111");

    // Set the colors for the histograms
    h1->SetLineColor(kRed); // Set h1 color to red
    h2->SetLineColor(kBlue); // Set h2 color to blue

    // Create a canvas and draw the histograms
    h2->Draw(); 
    h1->Draw("same");

    // Save the canvas as an image if desired
 //   c1->SaveAs("pgun_h1.png");
}

void plotperpdgid(int ityp,  std::string what, float min, float max, std::string cut, std::string f="", std::string opt="") {

    int nbin = 25;
    TH1F *h1 = new TH1F("h1", (std::string(what)+" for electrons").c_str(), nbin, min, max);
    TH1F *h2 = new TH1F("h2", (std::string(what)+" for pi0").c_str(), nbin, min, max);
    TH1F *h3 = new TH1F("h3", (std::string(what)+" for pi+-").c_str(), nbin, min, max);
    
    std::string drawCommand1 = what + ">>h1";
    pgtree->Draw(drawCommand1.c_str(), ("abs(m_pdg_id)==11" + cut).c_str());
    std::string drawCommand2 = what + ">>h2";
    pgtree->Draw(drawCommand2.c_str(), ("abs(m_pdg_id)==111" + cut).c_str());
    std::string drawCommand3 = what + ">>h3";
    pgtree->Draw(drawCommand3.c_str(), ("abs(m_pdg_id)==211" + cut).c_str());

#if 0
    int iwhat_value;
    float what_value;
    int pdg_id;

    // Set the branch addresses
    if(ityp == 0) {
    pgtree->SetBranchAddress(what, &iwhat_value);
    } else {
    pgtree->SetBranchAddress(what, &what_value);
    }
    pgtree->SetBranchAddress("m_pdg_id", &pdg_id);

     // Loop over the entries and fill the histogram
    Long64_t nentries = pgtree->GetEntries();

    for (Long64_t i = 0; i < nentries; i++) {
        pgtree->GetEntry(i);

        if(ityp==0) what_value = iwhat_value;

        // Apply the condition and fill the histogram
        if (abs(pdg_id) == 11) {
            h1->Fill(what_value);
        }
        if (abs(pdg_id) == 111) {
            h2->Fill(what_value);
        }
        if (abs(pdg_id) == 211) {
            h3->Fill(what_value);
        }
    }
#endif

    TCanvas *c1 = new TCanvas("c1", "Longitudinal energy profile", 1024 , 800);
    gStyle->SetStatFontSize(0.05);  // Set the font size for the statistics box

    if(opt == "same") {
        h1->SetFillColor(kBlue);
        h1->SetFillStyle(1001);

        h2->SetFillColor(kGreen);
        h2->SetFillStyle(1001);

        h3->SetFillColor(kRed);
        h3->SetFillStyle(1001);
        h1->Draw("hist");
        h2->Draw("hist,same");
        h3->Draw("hist,same");
    }
    else
    {
        c1->Divide(1, 3);
        c1->cd(1);
        h1->Draw();
        h1->GetXaxis()->SetLabelSize(0.075); // Set X-axis label size
        h1->GetYaxis()->SetLabelSize(0.075); // Set Y-axis label size
        c1->cd(2);
        h2->Draw();
        h2->GetXaxis()->SetLabelSize(0.075); // Set X-axis label size
        h2->GetYaxis()->SetLabelSize(0.075); // Set Y-axis label size
        c1->cd(3);
        h3->Draw();
        h3->GetXaxis()->SetLabelSize(0.075); // Set X-axis label size
        h3->GetYaxis()->SetLabelSize(0.075); // Set Y-axis label size
    }
    c1->Modified();
    c1->Update();
    if(f=="")
    c1->SaveAs(("pgunplots/" + std::string(what)+".png").c_str());
    else 
    c1->SaveAs(("pgunplots/" + f).c_str());
}


void plotperpdgid2D(std::string what, float min, float max, std::string cut, std::string f="") {

    TH2F *h1 = new TH2F("h1", (std::string(what)+" for electrons").c_str(), 20, min, max, 20, 8,16);
    TH2F *h2 = new TH2F("h2", (std::string(what)+" for pi0").c_str(), 20, min, max, 20, 8, 16);
    TH2F *h3 = new TH2F("h3", (std::string(what)+" for pi+-").c_str(), 20, min, max, 20, 8, 16);
    
    std::string drawCommand1 = what + ">>h1";
    pgtree->Draw(drawCommand1.c_str(), ("abs(m_pdg_id)==11" + cut).c_str());
    std::string drawCommand2 = what + ">>h2";
    pgtree->Draw(drawCommand2.c_str(), ("abs(m_pdg_id)==111" + cut).c_str());
    std::string drawCommand3 = what + ">>h3";
    pgtree->Draw(drawCommand3.c_str(), ("abs(m_pdg_id)==211" + cut).c_str());

       TCanvas *c1 = new TCanvas("c1", "Longitudinal energy profile", 1024 , 800);
    gStyle->SetStatFontSize(0.05);  // Set the font size for the statistics box

    c1->Divide(2,1);
    c1->cd(1);
    h1->Draw();
    h1->GetXaxis()->SetLabelSize(0.075);  // Set X-axis label size
    h1->GetYaxis()->SetLabelSize(0.075);  // Set Y-axis label size   
    c1->cd(2);
    h2->Draw();
    h2->GetXaxis()->SetLabelSize(0.075);  // Set X-axis label size
    h2->GetYaxis()->SetLabelSize(0.075);  // Set Y-axis label size   
     c1->Modified();
    c1->Update();
    if(f=="")
    c1->SaveAs(("pgunplots/" + std::string("scat")+".png").c_str());
    else 
    c1->SaveAs(("pgunplots/" + f).c_str());
}


void pgunplots() {
    gSystem->Load("libTPORec.so");
    TFile f("Batch-TPORecevent_1000000001_v2.root");

    pgtree = (TTree*)f.Get("ParticleGun");

//    pgtree->Draw("m_pdg_id");

    plotperpdgid(0,"m_pdg_id", -300,300,"");
    plotperpdgid(1,"m_energy", 0,3000,"");
    plotperpdgid(1,"ep_chi2_per_ndf",0,1000," && m_energy>100");
    std::string finalcut = " && m_energy>100  && ep_b > 0";
    plotperpdgid(1,"ep_a", -10, 20.,finalcut);
    plotperpdgid(1,"ep_b", 0, 1.,finalcut);
    plotperpdgid(1,"ep_c", -2, 5,finalcut);
    plotperpdgid(1,"ep_a/ep_b-1/ep_b-log(m_energy)-4", -5, 5,finalcut, "c.png", "same");

    plotperpdgid2D("log(ep_E0):(ep_a-1)/(ep_b)", 0., 20., finalcut);
}
