#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include <vector>

void analyze_muon_data() {
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
    
    // Statistics variables
    int total_muons = 0;
    int successful_fits = 0;
    double sum_truth_p = 0, sum_reco_p = 0;
    double sum_res_p = 0, sum_res_p_sq = 0;
    int valid_resolution_entries = 0;
    
    // Process events
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Analyzing " << nEntries << " events..." << std::endl;
    
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        // Process muons in this event
        for (size_t j = 0; j < muon_trackID->size(); j++) {
            total_muons++;
            
            // Calculate truth momentum magnitude
            double truth_px = (*muon_truth_px)[j];
            double truth_py = (*muon_truth_py)[j];
            double truth_pz = (*muon_truth_pz)[j];
            double truth_p = TMath::Sqrt(truth_px*truth_px + truth_py*truth_py + truth_pz*truth_pz);
            
            sum_truth_p += truth_p;
            
            std::cout << "Event " << i << ", Muon " << j << ":" << std::endl;
            std::cout << "  Track ID: " << (*muon_trackID)[j] << ", PDG: " << (*muon_pdg)[j] << std::endl;
            std::cout << "  Truth momentum: (" << truth_px << ", " << truth_py << ", " << truth_pz << ") = " << truth_p << " GeV" << std::endl;
            std::cout << "  Hits: " << (*muon_nhits)[j] << ", Stations: " << (*muon_nstations)[j] << std::endl;
            std::cout << "  Fit success: " << ((*muon_fit_success)[j] ? "Yes" : "No") << std::endl;
            
            // Check if fit was successful
            if ((*muon_fit_success)[j]) {
                successful_fits++;
                
                // Calculate reco momentum magnitude
                double reco_px = (*muon_reco_px)[j];
                double reco_py = (*muon_reco_py)[j];
                double reco_pz = (*muon_reco_pz)[j];
                double reco_p = TMath::Sqrt(reco_px*reco_px + reco_py*reco_py + reco_pz*reco_pz);
                
                sum_reco_p += reco_p;
                
                std::cout << "  Reco momentum: (" << reco_px << ", " << reco_py << ", " << reco_pz << ") = " << reco_p << " GeV" << std::endl;
                std::cout << "  Chi2/NDF: " << (*muon_chi2)[j] << "/" << (*muon_ndf)[j];
                if ((*muon_ndf)[j] > 0) {
                    std::cout << " = " << ((*muon_chi2)[j] / (*muon_ndf)[j]);
                }
                std::cout << std::endl;
                
                // Calculate resolution
                if (truth_p > 1e-6) {
                    double res_p = (reco_p - truth_p) / truth_p;
                    sum_res_p += res_p;
                    sum_res_p_sq += res_p * res_p;
                    valid_resolution_entries++;
                    std::cout << "  Momentum resolution: " << res_p * 100 << "%" << std::endl;
                }
            }
            std::cout << std::endl;
            
            // Only show first few examples to avoid too much output
            if (total_muons >= 10) break;
        }
        if (total_muons >= 10) break;
    }
    
    // Print summary statistics
    std::cout << "\n=== SUMMARY STATISTICS ===" << std::endl;
    std::cout << "Total events processed: " << nEntries << std::endl;
    std::cout << "Total muons: " << total_muons << std::endl;
    std::cout << "Successful fits: " << successful_fits << std::endl;
    std::cout << "Fit success rate: " << (total_muons > 0 ? 100.0 * successful_fits / total_muons : 0) << "%" << std::endl;
    
    if (total_muons > 0) {
        std::cout << "Average truth momentum: " << sum_truth_p / total_muons << " GeV" << std::endl;
    }
    
    if (successful_fits > 0) {
        std::cout << "Average reco momentum: " << sum_reco_p / successful_fits << " GeV" << std::endl;
    }
    
    if (valid_resolution_entries > 0) {
        double mean_res = sum_res_p / valid_resolution_entries;
        double rms_res = TMath::Sqrt(sum_res_p_sq / valid_resolution_entries - mean_res * mean_res);
        std::cout << "Mean momentum resolution: " << mean_res * 100 << "%" << std::endl;
        std::cout << "RMS momentum resolution: " << rms_res * 100 << "%" << std::endl;
    }
    
    file->Close();
}