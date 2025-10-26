#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TGraph.h>
#include <iostream>
#include <vector>

void plot_muon_reconstruction() {
    // Set ROOT style
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    
    // Open the ROOT file
    TFile* file = TFile::Open("build/scifi_hits.root");
    if (!file || file->IsZombie()) {
        std::cout << "Error: Cannot open file build/scifi_hits.root" << std::endl;
        return;
    }
    
    TTree* tree = (TTree*)file->Get("Hits");
    if (!tree) {
        std::cout << "Error: Cannot find tree 'Hits'" << std::endl;
        file->Close();
        return;
    }
    
    // Declare variables for reading the tree
    std::vector<int>* muon_trackID = nullptr;
    std::vector<int>* muon_pdg = nullptr;
    std::vector<double>* muon_truth_px = nullptr;
    std::vector<double>* muon_truth_py = nullptr;
    std::vector<double>* muon_truth_pz = nullptr;
    std::vector<double>* muon_reco_px = nullptr;
    std::vector<double>* muon_reco_py = nullptr;
    std::vector<double>* muon_reco_pz = nullptr;
    std::vector<int>* muon_nhits = nullptr;
    std::vector<int>* muon_nstations = nullptr;
    std::vector<double>* muon_chi2 = nullptr;
    std::vector<double>* muon_ndf = nullptr;
    std::vector<bool>* muon_fit_success = nullptr;
    
    // Set branch addresses
    tree->SetBranchAddress("muon_trackID", &muon_trackID);
    tree->SetBranchAddress("muon_pdg", &muon_pdg);
    tree->SetBranchAddress("muon_truth_px", &muon_truth_px);
    tree->SetBranchAddress("muon_truth_py", &muon_truth_py);
    tree->SetBranchAddress("muon_truth_pz", &muon_truth_pz);
    tree->SetBranchAddress("muon_reco_px", &muon_reco_px);
    tree->SetBranchAddress("muon_reco_py", &muon_reco_py);
    tree->SetBranchAddress("muon_reco_pz", &muon_reco_pz);
    tree->SetBranchAddress("muon_nhits", &muon_nhits);
    tree->SetBranchAddress("muon_nstations", &muon_nstations);
    tree->SetBranchAddress("muon_chi2", &muon_chi2);
    tree->SetBranchAddress("muon_ndf", &muon_ndf);
    tree->SetBranchAddress("muon_fit_success", &muon_fit_success);
    
    // Create histograms
    TH1F* h_truth_p = new TH1F("h_truth_p", "Truth Momentum;p_{truth} [GeV];Entries", 100, 0, 1000);
    TH1F* h_reco_p = new TH1F("h_reco_p", "Reconstructed Momentum;p_{reco} [GeV];Entries", 100, 0, 1000);
    TH2F* h_truth_vs_reco = new TH2F("h_truth_vs_reco", "Truth vs Reco Momentum;p_{truth} [GeV];p_{reco} [GeV]", 100, 0, 1000, 100, 0, 1000);
    
    TH1F* h_res_px = new TH1F("h_res_px", "p_{x} Resolution;(p_{x,reco} - p_{x,truth})/p_{x,truth};Entries", 100, -1, 1);
    TH1F* h_res_py = new TH1F("h_res_py", "p_{y} Resolution;(p_{y,reco} - p_{y,truth})/p_{y,truth};Entries", 100, -1, 1);
    TH1F* h_res_pz = new TH1F("h_res_pz", "p_{z} Resolution;(p_{z,reco} - p_{z,truth})/p_{z,truth};Entries", 100, -1, 1);
    TH1F* h_res_p = new TH1F("h_res_p", "Total Momentum Resolution;(p_{reco} - p_{truth})/p_{truth};Entries", 100, -1, 1);
    
    TH1F* h_nhits = new TH1F("h_nhits", "Number of Hits per Track;N_{hits};Entries", 50, 0, 50);
    TH1F* h_nstations = new TH1F("h_nstations", "Number of Stations per Track;N_{stations};Entries", 15, 0, 15);
    TH1F* h_chi2_ndf = new TH1F("h_chi2_ndf", "Track Fit Quality;#chi^{2}/NDF;Entries", 100, 0, 10);
    
    // Process events
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Processing " << nEntries << " events..." << std::endl;
    
    int total_muons = 0;
    int successful_fits = 0;
    
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        if (i % 1000 == 0) std::cout << "Processing event " << i << std::endl;
        
        // Process muons in this event
        for (size_t j = 0; j < muon_trackID->size(); j++) {
            total_muons++;
            
            // Calculate truth momentum magnitude
            double truth_px = (*muon_truth_px)[j];
            double truth_py = (*muon_truth_py)[j];
            double truth_pz = (*muon_truth_pz)[j];
            double truth_p = TMath::Sqrt(truth_px*truth_px + truth_py*truth_py + truth_pz*truth_pz);
            
            // Calculate reco momentum magnitude
            double reco_px = (*muon_reco_px)[j];
            double reco_py = (*muon_reco_py)[j];
            double reco_pz = (*muon_reco_pz)[j];
            double reco_p = TMath::Sqrt(reco_px*reco_px + reco_py*reco_py + reco_pz*reco_pz);
            
            // Fill basic histograms
            h_truth_p->Fill(truth_p);
            h_nhits->Fill((*muon_nhits)[j]);
            h_nstations->Fill((*muon_nstations)[j]);
            
            // Check if fit was successful
            if ((*muon_fit_success)[j]) {
                successful_fits++;
                
                h_reco_p->Fill(reco_p);
                h_truth_vs_reco->Fill(truth_p, reco_p);
                
                // Calculate resolutions (only for successful fits and non-zero truth values)
                if (TMath::Abs(truth_px) > 1e-6) {
                    double res_px = (reco_px - truth_px) / truth_px;
                    h_res_px->Fill(res_px);
                }
                if (TMath::Abs(truth_py) > 1e-6) {
                    double res_py = (reco_py - truth_py) / truth_py;
                    h_res_py->Fill(res_py);
                }
                if (TMath::Abs(truth_pz) > 1e-6) {
                    double res_pz = (reco_pz - truth_pz) / truth_pz;
                    h_res_pz->Fill(res_pz);
                }
                if (truth_p > 1e-6) {
                    double res_p = (reco_p - truth_p) / truth_p;
                    h_res_p->Fill(res_p);
                }
                
                // Fill fit quality
                if ((*muon_ndf)[j] > 0) {
                    h_chi2_ndf->Fill((*muon_chi2)[j] / (*muon_ndf)[j]);
                }
            }
        }
    }
    
    std::cout << "Summary:" << std::endl;
    std::cout << "Total muons: " << total_muons << std::endl;
    std::cout << "Successful fits: " << successful_fits << std::endl;
    std::cout << "Fit success rate: " << (total_muons > 0 ? 100.0 * successful_fits / total_muons : 0) << "%" << std::endl;
    
    // Create canvas and plot
    TCanvas* c1 = new TCanvas("c1", "Muon Reconstruction Analysis", 1200, 800);
    c1->Divide(3, 2);
    
    // Plot 1: Truth vs Reco momentum
    c1->cd(1);
    h_truth_p->SetLineColor(kBlue);
    h_truth_p->SetFillColor(kBlue);
    h_truth_p->SetFillStyle(3004);
    h_truth_p->Draw();
    
    h_reco_p->SetLineColor(kRed);
    h_reco_p->SetFillColor(kRed);
    h_reco_p->SetFillStyle(3005);
    h_reco_p->Draw("same");
    
    TLegend* leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg1->AddEntry(h_truth_p, "Truth", "f");
    leg1->AddEntry(h_reco_p, "Reconstructed", "f");
    leg1->Draw();
    
    // Plot 2: Truth vs Reco 2D correlation
    c1->cd(2);
    h_truth_vs_reco->Draw("colz");
    
    // Add diagonal line
    TGraph* diagonal = new TGraph(2);
    diagonal->SetPoint(0, 0, 0);
    diagonal->SetPoint(1, 1000, 1000);
    diagonal->SetLineColor(kRed);
    diagonal->SetLineWidth(2);
    diagonal->Draw("same");
    
    // Plot 3: Total momentum resolution
    c1->cd(3);
    h_res_p->SetLineColor(kBlack);
    h_res_p->Fit("gaus", "Q");
    h_res_p->Draw();
    
    // Plot 4: px, py, pz resolutions
    c1->cd(4);
    h_res_px->SetLineColor(kRed);
    h_res_py->SetLineColor(kGreen);
    h_res_pz->SetLineColor(kBlue);
    
    h_res_px->Draw();
    h_res_py->Draw("same");
    h_res_pz->Draw("same");
    
    TLegend* leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg2->AddEntry(h_res_px, "p_{x} resolution", "l");
    leg2->AddEntry(h_res_py, "p_{y} resolution", "l");
    leg2->AddEntry(h_res_pz, "p_{z} resolution", "l");
    leg2->Draw();
    
    // Plot 5: Number of hits and stations
    c1->cd(5);
    h_nhits->SetLineColor(kBlue);
    h_nhits->Draw();
    
    // Plot 6: Fit quality
    c1->cd(6);
    h_chi2_ndf->SetLineColor(kMagenta);
    h_chi2_ndf->Draw();
    
    c1->Update();
    c1->SaveAs("muon_reconstruction_analysis.png");
    c1->SaveAs("muon_reconstruction_analysis.pdf");
    
    // Create second canvas for detailed resolution plots
    TCanvas* c2 = new TCanvas("c2", "Resolution Details", 800, 600);
    c2->Divide(2, 2);
    
    c2->cd(1);
    h_res_px->Fit("gaus", "Q");
    h_res_px->Draw();
    
    c2->cd(2);
    h_res_py->Fit("gaus", "Q");
    h_res_py->Draw();
    
    c2->cd(3);
    h_res_pz->Fit("gaus", "Q");
    h_res_pz->Draw();
    
    c2->cd(4);
    h_nstations->Draw();
    
    c2->Update();
    c2->SaveAs("muon_resolution_details.png");
    c2->SaveAs("muon_resolution_details.pdf");
    
    std::cout << "Plots saved as muon_reconstruction_analysis.png/pdf and muon_resolution_details.png/pdf" << std::endl;
    
    file->Close();
}