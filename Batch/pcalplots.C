TTree* pgtree;

void pcalplots() {
    gSystem->Load("LibTPORec.so");

    TFile f("Batch-TPORecevent_88800100_0_99999999.root");
//    TFile f("Batch-TPORecevent_9990000_0_99999999.root");

    pgtree = (TTree*)f.Get("ParticleGun");

    TH1F *h1 = new TH1F("h1", "Visible Energy Distribution", 100, 0, 200);
    pgtree->Draw("ef_evis>>h1");
    h1->SetLineColor(kRed); // Set h1 color to red
    h1->GetXaxis()->SetTitle("Visible Energy [GeV]");
    h1->GetYaxis()->SetTitle("Entries");
    TCanvas *c1 = new TCanvas("c1", "Visible Energy Distribution", 800, 600);
    gStyle->SetStatFontSize(0.05);  // Set the font size for the statistics box
    h1->Draw();
    c1->Modified();
    c1->Update();
    c1->SaveAs("visible_energy.png");

    TH1F *h2 = new TH1F("h2"," Transverse Energy Distribution", 100, 0, 10);
    pgtree->Draw("ef_et>>h2");
    h2->SetLineColor(kBlue); // Set h2 color to blue
    h2->GetXaxis()->SetTitle("Transverse Energy [GeV]");
    h2->GetYaxis()->SetTitle("Entries");
    TCanvas *c2 = new TCanvas("c2", "Transverse Energy Distribution", 800, 600);
    gStyle->SetStatFontSize(0.05);  // Set the font size for the statistics box
    h2->Draw();
    c2->Modified();
    c2->Update();
    c2->SaveAs("transverse_energy.png");

    TH1F *h3 = new TH1F("h3"," FASER Calorimeter Energy Distribution", 100, -3, 3);
    pgtree->Draw("ef_fasercal_x>>h3");
    h3->SetLineColor(kGreen+2); // Set h3 color to green
    h3->GetXaxis()->SetTitle("FASER Calorimeter Energy [GeV]");
    h3->GetYaxis()->SetTitle("Entries");
    TCanvas *c3 = new TCanvas("c3", "FASER Calorimeter Energy Distribution", 800, 600);
    h3->Draw();
    c3->Modified();
    c3->Update();
    c3->SaveAs("fasercal_energy_x.png");

    TH1F *h4 = new TH1F("h4"," FASER Calorimeter Energy Distribution", 100, -3, 3);
    pgtree->Draw("ef_fasercal_y>>h4");
    h4->SetLineColor(kMagenta); // Set h4 color to magenta
    h4->GetXaxis()->SetTitle("FASER Calorimeter Energy [GeV]");
    h4->GetYaxis()->SetTitle("Entries");
    TCanvas *c4 = new TCanvas("c4", "FASER Calorimeter Energy Distribution", 800, 600);
    h4->Draw();
    c4->Modified();
    c4->Update();
    c4->SaveAs("fasercal_energy_y.png");

    TH1F *h5 = new TH1F("h5"," FASER Calorimeter Energy Distribution", 100, 0, 100);
    pgtree->Draw("ef_fasercal_z>>h5");
    h5->SetLineColor(kCyan); // Set h5 color to cyan
    h5->GetXaxis()->SetTitle("FASER Calorimeter Energy [GeV]");
    h5->GetYaxis()->SetTitle("Entries");
    TCanvas *c5 = new TCanvas("c5", "FASER Calorimeter Energy Distribution", 800, 600);
    h5->Draw();
    c5->Modified();
    c5->Update();
    c5->SaveAs("fasercal_energy_z.png");

    float b = 5;
    TH1F *h6 = new TH1F("h6"," Calorimetric Energy Distribution", 100, -b, b);
    pgtree->Draw("ef_fasercal_x+4.16*ef_ecal_x+38.76*ef_hcal_x>>h6");
    h6->SetLineColor(kOrange); // Set h6 color to orange
    h6->GetXaxis()->SetTitle("Calorimetric Energy [GeV]");
    h6->GetYaxis()->SetTitle("Entries");
    TCanvas *c6 = new TCanvas("c6", "Calorimetric Energy Distribution", 800, 600);
    h6->Draw();
    c6->Modified();
    c6->Update();
    c6->SaveAs("comp_energy_x.png");

    TH1F *h7 = new TH1F("h7"," Calorimetric Energy Distribution", 100, -b, b);
    pgtree->Draw("ef_fasercal_y+4.16*ef_ecal_y+38.76*ef_hcal_y>>h7");
    h7->SetLineColor(kOrange); // Set h7 color to orange
    h7->GetXaxis()->SetTitle("Calorimetric Energy [GeV]");
    h7->GetYaxis()->SetTitle("Entries");
    TCanvas *c7 = new TCanvas("c7", "Calorimetric Energy Distribution", 800, 600);
    h7->Draw();
    c7->Modified();
    c7->Update();
    c7->SaveAs("comp_energy_y.png");

    TH1F *h8 = new TH1F("h8"," Calorimetric Energy Distribution", 100, 0, 200);
    pgtree->Draw("ef_fasercal_z+4.16*ef_ecal_z+38.76*ef_hcal_z>>h8");
    h8->SetLineColor(kOrange); // Set h8 color to orange
    h8->GetXaxis()->SetTitle("Calorimetric Energy [GeV]");
    h8->GetYaxis()->SetTitle("Entries");
    TCanvas *c8 = new TCanvas("c8", "Calorimetric Energy Distribution", 800, 600);
    h8->Draw();
    c8->Modified();
    c8->Update();
    c8->SaveAs("comp_energy_z.png");

    TH1F *h9 = new TH1F("h9"," normalized energy", 100, 0, 2.0);
    pgtree->Draw("(1.2*ef_fasercal_z+4.16*ef_ecal_z+38.76*ef_hcal_z)/m_jet_energy>>h9");
    h9->SetLineColor(kBlack); // Set h9 color to black
    h9->GetXaxis()->SetTitle("normalized energy");
    h9->GetYaxis()->SetTitle("Entries");
    TCanvas *c9 = new TCanvas("c9", "normalized energy", 800, 600);
    h9->Draw();
    c9->Modified();
    c9->Update();
    c9->SaveAs("comp_energy_normalized.png");

    TH1F *h10 = new TH1F("h10"," normalized px",100,-0.3,0.3);
    pgtree->Draw("((ef_fasercal_x+4.16*ef_ecal_x+38.76*ef_hcal_x)/m_jet_energy)>>h10");
    h10->SetLineColor(kBlack); // Set h10 color to black
    h10->GetXaxis()->SetTitle("normalized px");
    h10->GetYaxis()->SetTitle("Entries");
    TCanvas *c10 = new TCanvas("c10", "normalized px", 800, 600);
    h10->Draw();
    c10->SetLogy();
    c10->Modified();
    c10->Update();
    c10->SaveAs("comp_px_normalized.png");

    TH1F *h11 = new TH1F("h11"," normalized py",100,-0.3,0.3);
    pgtree->Draw("((ef_fasercal_y+4.16*ef_ecal_y+38.76*ef_hcal_y)/m_jet_energy)>>h11");
    h11->SetLineColor(kBlack); // Set h11 color to black
    h11->GetXaxis()->SetTitle("normalized py");
    h11->GetYaxis()->SetTitle("Entries");
    TCanvas *c11 = new TCanvas("c11", "normalized py", 800, 600);
    h11->Draw();
    c11->SetLogy();
    c11->Modified();
    c11->Update();
    c11->SaveAs("comp_py_normalized.png");
}