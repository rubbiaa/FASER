// Macro to plot beam spot and extract center from flux tree
void plot_beam_spot() {
    // Open the ROOT file
    TFile *f = TFile::Open("events_light_4x4.root");
    //TFile *f = TFile::Open("events_charm_4x4.root");
    if (!f || f->IsZombie()) {
        //std::cerr << "Error: Could not open file events_light_4x4.root" << std::endl;
        std::cerr << "Error: Could not open file events_charm_4x4.root" << std::endl;
        return;
    }
    
    // Get the flux tree
    TTree *flux = (TTree*)f->Get("flux");
    if (!flux) {
        std::cerr << "Error: Could not find flux tree" << std::endl;
        return;
    }
    
    std::cout << "=== Beam Spot Analysis ===" << std::endl;
    std::cout << "Number of entries: " << flux->GetEntries() << std::endl;
    
    // Create canvas
    TCanvas *c1 = new TCanvas("c1", "Beam Spot", 1200, 500);
    c1->Divide(3, 1);
    
    // Plot 1: 2D beam spot (vtxy vs vtxx)
    c1->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    flux->Draw("vtxy:vtxx>>h2d(100,-2.5,2.5,100,-2.5,2.5)", "", "colz");
    TH2F *h2d = (TH2F*)gDirectory->Get("h2d");
    if (h2d) {
        h2d->SetTitle("Beam Spot Distribution;vtxx (m);vtxy (m)");
        h2d->GetXaxis()->SetTitleSize(0.045);
        h2d->GetYaxis()->SetTitleSize(0.045);
    }
    
    // Plot 2: X projection
    c1->cd(2);
    flux->Draw("vtxx>>hx(100,-2.5,2.5)");
    TH1F *hx = (TH1F*)gDirectory->Get("hx");
    if (hx) {
        hx->SetTitle("Horizontal Beam Distribution;vtxx (m);Entries");
        hx->SetLineColor(kBlue);
        hx->SetLineWidth(2);
    }
    
    // Plot 3: Y projection
    c1->cd(3);
    flux->Draw("vtxy>>hy(100,-2.5,2.5)");
    TH1F *hy = (TH1F*)gDirectory->Get("hy");
    if (hy) {
        hy->SetTitle("Vertical Beam Distribution;vtxy (m);Entries");
        hy->SetLineColor(kRed);
        hy->SetLineWidth(2);
    }
    
    // Calculate beam spot center and RMS with weights
    Double_t vtxx, vtxy, wgt;
    flux->SetBranchAddress("vtxx", &vtxx);
    flux->SetBranchAddress("vtxy", &vtxy);
    flux->SetBranchAddress("wgt", &wgt);
    
    // Unweighted calculations
    Double_t sum_x = 0, sum_y = 0;
    Double_t sum_x2 = 0, sum_y2 = 0;
    Long64_t nEntries = flux->GetEntries();
    
    // Weighted calculations
    Double_t sum_w = 0;
    Double_t sum_wx = 0, sum_wy = 0;
    Double_t sum_wx2 = 0, sum_wy2 = 0;
    
    for (Long64_t i = 0; i < nEntries; i++) {
        flux->GetEntry(i);
        
        // Unweighted
        sum_x += vtxx;
        sum_y += vtxy;
        sum_x2 += vtxx * vtxx;
        sum_y2 += vtxy * vtxy;
        
        // Weighted
        sum_w += wgt;
        sum_wx += wgt * vtxx;
        sum_wy += wgt * vtxy;
        sum_wx2 += wgt * vtxx * vtxx;
        sum_wy2 += wgt * vtxy * vtxy;
    }
    
    // Unweighted results
    Double_t mean_x = sum_x / nEntries;
    Double_t mean_y = sum_y / nEntries;
    Double_t rms_x = TMath::Sqrt(sum_x2 / nEntries - mean_x * mean_x);
    Double_t rms_y = TMath::Sqrt(sum_y2 / nEntries - mean_y * mean_y);
    
    // Weighted results
    Double_t mean_x_weighted = sum_wx / sum_w;
    Double_t mean_y_weighted = sum_wy / sum_w;
    Double_t rms_x_weighted = TMath::Sqrt(sum_wx2 / sum_w - mean_x_weighted * mean_x_weighted);
    Double_t rms_y_weighted = TMath::Sqrt(sum_wy2 / sum_w - mean_y_weighted * mean_y_weighted);
    
    std::cout << "\n=== Beam Spot Center (Unweighted) ===" << std::endl;
    std::cout << "Center X: " << mean_x << " m" << std::endl;
    std::cout << "Center Y: " << mean_y << " m" << std::endl;
    std::cout << "RMS X:    " << rms_x << " m" << std::endl;
    std::cout << "RMS Y:    " << rms_y << " m" << std::endl;
    
    std::cout << "\n=== Beam Spot Center (Weighted) ===" << std::endl;
    std::cout << "Center X: " << mean_x_weighted << " m" << std::endl;
    std::cout << "Center Y: " << mean_y_weighted << " m" << std::endl;
    std::cout << "RMS X:    " << rms_x_weighted << " m" << std::endl;
    std::cout << "RMS Y:    " << rms_y_weighted << " m" << std::endl;
    std::cout << "Sum of weights: " << sum_w << std::endl;
    
    // Add markers at the weighted center on the 2D plot
    c1->cd(1);
    TMarker *marker = new TMarker(mean_x_weighted, mean_y_weighted, 29);
    marker->SetMarkerColor(kWhite);
    marker->SetMarkerSize(3);
    //marker->Draw();
    
    TMarker *marker2 = new TMarker(mean_x_weighted, mean_y_weighted, 5);
    marker2->SetMarkerColor(kBlack);
    marker2->SetMarkerSize(2);
    //marker2->Draw();
    
    // Add text with weighted center coordinates
    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.035);
    text->SetTextColor(kBlack);
    //text->DrawLatex(0.15, 0.88, "Weighted:");
    //text->DrawLatex(0.15, 0.85, Form("Center: (%.3f, %.3f) m", mean_x_weighted, mean_y_weighted));
    text->DrawLatex(0.15, 0.85, Form("Center: (%.3f, %.3f) m", mean_x, mean_y));
    //text->DrawLatex(0.15, 0.83, Form("RMS: (%.3f, %.3f) m", rms_x_weighted, rms_y_weighted));
    
    c1->Update();
    c1->SaveAs("beam_spot.png");
    c1->SaveAs("beam_spot.pdf");
    
    std::cout << "\nPlots saved as beam_spot.png and beam_spot.pdf" << std::endl;
    
    // Clean up
    // f->Close();
}
