// Compute neutrino pseudorapidity coverage from FASER GENIE ntuples
//
// Usage:
//   etaCoverage <genieroot> <outputRootFile>
// ------------------------------------------------------------

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>

#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>


// ----------------- Geometry (in mm) -----------------

// Upstream veto
constexpr double Z_FRONT_TARGET_MAX_MM = -1150.0;

// 3DCAL region
constexpr double Z_3DCAL_MIN_MM       = -1150.0;
constexpr double Z_3DCAL_MAX_MM       =  1150.0;

// rear ECAL (from GDML)
constexpr double Z_REARECAL_MIN_MM    = 1205.0;
constexpr double Z_REARECAL_MAX_MM    = 1550.0;

// AHCAL / rear HCAL
constexpr double Z_AHCAL_MIN_MM       = 1551.0;
constexpr double Z_AHCAL_MAX_MM       = 2450.0;

// Muon spectrometer
constexpr double Z_MUSPEC_MIN_MM      = 2451.0;
constexpr double Z_MUSPEC_MAX_MM      = 4350.0;

// ------------------------------------------------------------

void ensurePlotDir()
{
    const char *dirname = "plots";
    struct stat st{};
    if (stat(dirname, &st) == -1) {
        mkdir(dirname, 0755);
        std::cout << "Created directory 'plots/'\n";
    } else {
        std::cout << "Directory 'plots/' already exists\n";
    }
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        std::cout << "Usage: " << argv[0]
                  << " <genieroot> <outputRootFile>\n";
        std::cout << "  <genieroot>      Input FASER GENIE ROOT file(s), tree gFaser\n";
        std::cout << "  <outputRootFile> Output ROOT file with eta histograms\n";
        return 1;
    }

    std::string inputPattern  = argv[1];
    std::string outputRoot    = argv[2];

    ensurePlotDir();

    TChain *tree = new TChain("gFaser");
    tree->Add(inputPattern.c_str());
    Long64_t n_entries = tree->GetEntries();
    if (n_entries <= 0) {
        std::cerr << "ERROR: no entries found in tree gFaser from " << inputPattern << "\n";
        return 1;
    }

    std::cout << "Input GENIE file(s): " << inputPattern << "\n";
    std::cout << "Total entries in gFaser: " << n_entries << "\n";

    // Branches
    Double_t vx = 0, vy = 0, vz = 0; // in meters
    Int_t n = 0;
    std::vector<std::string> *name = nullptr;
    std::vector<int> *pdgc = nullptr;
    std::vector<int> *status = nullptr;
    std::vector<int> *firstMother = nullptr;
    std::vector<int> *lastMother = nullptr;
    std::vector<int> *firstDaughter = nullptr;
    std::vector<int> *lastDaughter = nullptr;
    std::vector<double> *px = nullptr;
    std::vector<double> *py = nullptr;
    std::vector<double> *pz = nullptr;
    std::vector<double> *E = nullptr;
    std::vector<double> *m = nullptr;
    std::vector<double> *M = nullptr;

    TBranch *b_pdgc   = nullptr;
    TBranch *b_status = nullptr;

    tree->SetBranchAddress("vx", &vx);
    tree->SetBranchAddress("vy", &vy);
    tree->SetBranchAddress("vz", &vz);
    tree->SetBranchAddress("n", &n);
    tree->SetBranchAddress("name", &name);
    tree->SetBranchAddress("pdgc", &pdgc, &b_pdgc);
    tree->SetBranchAddress("status", &status, &b_status);
    tree->SetBranchAddress("firstMother", &firstMother);
    tree->SetBranchAddress("lastMother", &lastMother);
    tree->SetBranchAddress("firstDaughter", &firstDaughter);
    tree->SetBranchAddress("lastDaughter", &lastDaughter);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("E", &E);
    tree->SetBranchAddress("m", &m);
    tree->SetBranchAddress("M", &M);

    int    nbins = 80;
    double etamin = 4.0;
    double etamax = 12.0;
    int    nbins_E   = 80;
    double Emin      = 0.0;
    double Emax      = 5000.0; // 5 TeV
    // Histograms for eta
    TH1D *h_eta_all      = new TH1D("h_eta_all",      "Neutrino #eta (all regions);#eta;Events",         nbins, etamin, etamax);
    TH1D *h_eta_3dcal    = new TH1D("h_eta_3dcal",    "Neutrino #eta (3DCAL);#eta;Events",               nbins, etamin, etamax);
    TH1D *h_eta_rearecal = new TH1D("h_eta_rearecal", "Neutrino #eta (rear ECAL);#eta;Events",           nbins, etamin, etamax);
    TH1D *h_eta_ahcal    = new TH1D("h_eta_ahcal",    "Neutrino #eta (AHCAL);#eta;Events",               nbins, etamin, etamax);
    TH1D *h_eta_muspec   = new TH1D("h_eta_muspec",   "Neutrino #eta (Muon Spectrometer);#eta;Events",   nbins, etamin, etamax);
    TH2D *h_E_vs_eta_3dcal = new TH2D("h_E_vs_eta_3dcal", "Neutrino E vs #eta (3DCAL);#eta;E_{#nu} [GeV]", nbins, etamin, etamax, nbins_E,   Emin,   Emax);
    TH2D *h_E_vs_eta_ahcal = new TH2D("h_E_vs_eta_ahcal", "Neutrino E vs #eta (AHCAL);#eta;E_{#nu} [GeV]", nbins, etamin, etamax, nbins_E,   Emin,   Emax);
    
    Long64_t total_after_vtx = 0;
    Long64_t total_with_nu   = 0;

    for (Long64_t ievt = 0; ievt < n_entries; ++ievt)
    {
        tree->GetEntry(ievt);

        if (ievt % 10000 == 0)
            std::cout << "Processing event " << ievt << " / " << n_entries << "...\n";

        // Convert vertex to mm
        double z_mm = vz * 1e3;

        // Upstream veto
        if (z_mm < Z_FRONT_TARGET_MAX_MM)
            continue;

        ++total_after_vtx;

        // Find primary neutrino in the list
        int nu_index = -1;
        for (int i = 0; i < n; ++i) {
            int pdg = pdgc->at(i);
            int apdg = std::abs(pdg);
            if (apdg == 12 || apdg == 14 || apdg == 16) {
                nu_index = i;
                break;
            }
        }

        if (nu_index < 0) {
            continue; // no neutrino found, skip
        }

        ++total_with_nu;

        double pxNu = px->at(nu_index);
        double pyNu = py->at(nu_index);
        double pzNu = pz->at(nu_index);

        double pt = std::sqrt(pxNu*pxNu + pyNu*pyNu);
        double theta = std::atan2(pt, pzNu); // angle wrt z-axis

        double ENu  = E->at(nu_index);

        if (theta <= 0.0 || theta >= M_PI) {
            continue;
        }

        double eta = -std::log(std::tan(theta/2.0));

        h_eta_all->Fill(eta);

        // Region classification by vertex z
        bool in_3dcal    = (z_mm >= Z_3DCAL_MIN_MM    && z_mm < Z_3DCAL_MAX_MM);
        bool in_rearecal = (z_mm >= Z_REARECAL_MIN_MM && z_mm < Z_REARECAL_MAX_MM);
        bool in_ahcal    = (z_mm >= Z_AHCAL_MIN_MM    && z_mm < Z_AHCAL_MAX_MM);
        bool in_muspec   = (z_mm >= Z_MUSPEC_MIN_MM   && z_mm < Z_MUSPEC_MAX_MM);

        if (in_3dcal)
        {
            h_eta_3dcal->Fill(eta);
            h_E_vs_eta_3dcal->Fill(eta, ENu);
        } 
        else if (in_rearecal)
            h_eta_rearecal->Fill(eta);
        else if (in_ahcal)
        {
            h_eta_ahcal->Fill(eta);
            h_E_vs_eta_ahcal->Fill(eta, ENu);
        }
        else if (in_muspec)
            h_eta_muspec->Fill(eta);
    }

    std::cout << "\n==========================================\n";
    std::cout << "Events after front-target veto (z > " << Z_FRONT_TARGET_MAX_MM << " mm): "
              << total_after_vtx << "\n";
    std::cout << "Events with identified incoming neutrino: " << total_with_nu << "\n";
    std::cout << "==========================================\n";

    // Save histograms to ROOT file
    TFile *fout = new TFile(outputRoot.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "ERROR: could not create output ROOT file " << outputRoot << "\n";
        return 1;
    }

    //h_eta_all->Write();
    h_eta_3dcal->Write();
    //h_eta_rearecal->Write();
    h_eta_ahcal->Write();
    //h_eta_muspec->Write();
    
    h_E_vs_eta_3dcal->Write();
    h_E_vs_eta_ahcal->Write();

    fout->Close();
    std::cout << "Eta histograms saved to " << outputRoot << "\n";

    // Quick PNGs
    gStyle->SetOptStat(1110);

    auto save1D = [](TH1 *h, const char *name) {
        TCanvas *c = new TCanvas();
        h->Draw();
        std::string out = std::string("plots/") + name;
        c->SaveAs(out.c_str());
        delete c;
    };

    //save1D(h_eta_all,      "eta_all.png");
    save1D(h_eta_3dcal,    "eta_3dcal.png");
    //save1D(h_eta_rearecal, "eta_rearecal.png");
    save1D(h_eta_ahcal,    "eta_ahcal.png");
    //save1D(h_eta_muspec,   "eta_muspec.png");

    gStyle->SetOptStat(0);

    auto save2D = [](TH2D *h, const char *name) {
        TCanvas *c = new TCanvas();
        h->Draw("COLZ");
        std::string out = std::string("plots/") + name;
        c->SaveAs(out.c_str());
        delete c;
    };

    save2D(h_E_vs_eta_3dcal, "E_vs_eta_3dcal.png");
    save2D(h_E_vs_eta_ahcal, "E_vs_eta_ahcal.png");
    
    std::cout << "Done.\n";

    return 0;
}
