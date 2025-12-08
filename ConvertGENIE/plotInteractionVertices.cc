// Plot neutrino interaction vertices (GENIE gFaser ntuples)
// Saves PNG plots into a "plots/" directory
//
// Usage:
//   plotInteractionVertices <genieroot> <outputRootFile>
//
// Example:
//   plotInteractionVertices FASERCal_tilted_5deg_1abInv_500.gfaser.root vertices.root
//
// This program:
//  - Reads neutrino interaction vertex from gFaser
//  - Converts (vx,vy,vz) from m → mm
//  - Fills histograms for: all, 3DCAL, AHCAL, muon spectrometer
//  - Dumps PNG plots into ./plots/
//  - Writes all histograms into outputRootFile
//
// ------------------------------------------------------------

#include <iostream>
#include <string>
#include <sys/stat.h>  
#include <sys/types.h>

#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TBox.h>
#include <TLine.h>
#include <TLegend.h>
// ----------------- Geometry (in mm) -----------------

constexpr double Z_FRONT_TARGET_MAX_MM = -1150.0;

// 3DCAL region
constexpr double Z_3DCAL_MIN_MM = -1150.0;
constexpr double Z_3DCAL_MAX_MM =  1150.0;

// ECAL region
constexpr double Z_ECAL_MIN_MM  =  1205.0;
constexpr double Z_ECAL_MAX_MM  =  1550.0;

// AHCAL region
constexpr double Z_AHCAL_MIN_MM = 1551.0;
constexpr double Z_AHCAL_MAX_MM = 2450.0;

// Muon spectrometer region
constexpr double Z_MS_MIN_MM    = 2451.0;
constexpr double Z_MS_MAX_MM    = 4350.0;

// ------------------------------------------------------------

void ensurePlotDirectory()
{
    const char *dirname = "plots";
    struct stat st = {0};

    if (stat(dirname, &st) == -1) {
        mkdir(dirname, 0755);
        std::cout << "Created directory 'plots/'" << std::endl;
    } else {
        std::cout << "Directory 'plots/' already exists" << std::endl;
    }
}

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cout << "Usage: " << argv[0]
                  << " <genieroot> <outputRootFile>" << std::endl;
        return 1;
    }

    std::string rootinputString = argv[1];
    std::string outputRootFile  = argv[2];

    ensurePlotDirectory();

    // Load GENIE files
    TChain *tree = new TChain("gFaser");
    tree->Add(rootinputString.c_str());
    Long64_t n_entries = tree->GetEntries();
    if (n_entries == 0) {
        std::cerr << "ERROR: no entries found." << std::endl;
        return 1;
    }

    // GENIE branches
    Double_t vx, vy, vz;
    Int_t n;
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

    TBranch *b_pdgc = nullptr;
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

    // Histograms
    TH1D *h_z_all     = new TH1D("h_z_all",     "Interaction z (all regions);z [mm];Events",          200, -3000., 7000.);
    TH1D *h_z_3dcal   = new TH1D("h_z_3dcal",   "Interaction z (3DCAL);z [mm];Events",                200, -3000., 7000.);
    TH1D *h_z_ecal    = new TH1D("h_z_ecal",    "Interaction z (ECAL);z [mm];Events",                 200, -3000., 7000.);
    TH1D *h_z_ahcal   = new TH1D("h_z_ahcal",   "Interaction z (AHCAL);z [mm];Events",                200, -3000., 7000.);
    TH1D *h_z_muspec  = new TH1D("h_z_muspec",  "Interaction z (MuonSpectrometer);z [mm];Events",     200, -3000., 7000.);

    TH2D *h_xz_all    = new TH2D("h_xz_all",    "x vs z (all);z [mm];x [mm]",                         200, -3000., 7000., 100, -2000., 2000.);
    TH2D *h_xz_3dcal  = new TH2D("h_xz_3dcal",  "x vs z (3DCAL);z [mm];x [mm]",                       200, -3000., 7000., 100, -2000., 2000.);
    TH2D *h_xz_ecal   = new TH2D("h_xz_ecal",   "x vs z (ECAL);z [mm];x [mm]",                        200, -3000., 7000., 100, -2000., 2000.);
    TH2D *h_xz_ahcal  = new TH2D("h_xz_ahcal",  "x vs z (AHCAL);z [mm];x [mm]",                       200, -3000., 7000., 100, -2000., 2000.);
    TH2D *h_xz_muspec = new TH2D("h_xz_muspec", "x vs z (MuonSpec);z [mm];x [mm]",                    200, -3000., 7000., 100, -2000., 2000.);

    // XY plots (transverse view)
    TH2D *h_xy_all      = new TH2D("h_xy_all",      "x vs y (all);x [mm];y [mm]", 100, -500., 1500., 100, -500., 1500.);
    TH2D *h_xy_3dcal    = new TH2D("h_xy_3dcal",    "x vs y (3DCAL);x [mm];y [mm]", 100, -500., 1500., 100, -500., 1500.);
    TH2D *h_xy_ecal     = new TH2D("h_xy_ecal",     "x vs y (rear ECAL);x [mm];y [mm]", 100, -500., 1500., 100, -500., 1500.);
    TH2D *h_xy_ahcal    = new TH2D("h_xy_ahcal",    "x vs y (AHCAL);x [mm];y [mm]", 100, -500., 1500., 100, -500., 1500.);
    TH2D *h_xy_muspec   = new TH2D("h_xy_muspec",   "x vs y (MuonSpec);x [mm];y [mm]", 100, -500., 1500., 100, -500., 1500.);
    // Loop over events
    for (Long64_t ievt = 0; ievt < n_entries; ++ievt)
    {
        tree->GetEntry(ievt);

        double x_mm = vx * 1e3;
        double y_mm = vy * 1e3;
        double z_mm = vz * 1e3;

        if (z_mm < Z_FRONT_TARGET_MAX_MM)
            continue;

        h_z_all->Fill(z_mm);
        h_xz_all->Fill(z_mm, x_mm);
        h_xy_all->Fill(x_mm, y_mm);

        if (z_mm >= Z_3DCAL_MIN_MM && z_mm < Z_3DCAL_MAX_MM) {
            h_z_3dcal->Fill(z_mm);
            h_xz_3dcal->Fill(z_mm, x_mm);
            h_xy_3dcal->Fill(x_mm, y_mm);
        }
        if (z_mm >= Z_ECAL_MIN_MM && z_mm < Z_ECAL_MAX_MM) {
            h_z_ecal->Fill(z_mm);
            h_xz_ecal->Fill(z_mm, x_mm);
            h_xy_ecal->Fill(x_mm, y_mm);
        }
        if (z_mm >= Z_AHCAL_MIN_MM && z_mm < Z_AHCAL_MAX_MM) {
            h_z_ahcal->Fill(z_mm);
            h_xz_ahcal->Fill(z_mm, x_mm);
            h_xy_ahcal->Fill(x_mm, y_mm);
        }
        if (z_mm >= Z_MS_MIN_MM && z_mm < Z_MS_MAX_MM) {
            h_z_muspec->Fill(z_mm);
            h_xz_muspec->Fill(z_mm, x_mm);
            h_xy_muspec->Fill(x_mm, y_mm);
        }
    }

    // Save histograms
    TFile *fout = new TFile(outputRootFile.c_str(), "RECREATE");
    h_z_all->Write();
    h_z_3dcal->Write();
    h_z_ecal->Write();
    h_z_ahcal->Write();
    h_z_muspec->Write();
    h_xz_all->Write();
    h_xz_3dcal->Write();
    h_xz_ecal->Write();
    h_xz_ahcal->Write();
    h_xz_muspec->Write();
    h_xy_all->Write();
    h_xy_3dcal->Write();
    h_xy_ecal->Write();
    h_xy_ahcal->Write();
    h_xy_muspec->Write();
    fout->Close();

    /*
    // Save PNGs to plots/
    gStyle->SetOptStat(1110);

    auto savePlot = [&](TH1 *h, const char *filename) {
        TCanvas *c = new TCanvas();
        h->Draw();
        std::string out = std::string("plots/") + filename;
        c->SaveAs(out.c_str());
        delete c;
    };

    auto savePlot2D = [&](TH2 *h, const char *filename) {
        TCanvas *c = new TCanvas();
        h->Draw("COLZ");
        std::string out = std::string("plots/") + filename;
        c->SaveAs(out.c_str());
        delete c;
    };

    savePlot(h_z_all,    "vertex_z_all.png");
    savePlot(h_z_3dcal,  "vertex_z_3dcal.png");
    savePlot(h_z_ecal,   "vertex_z_ecal.png");
    savePlot(h_z_ahcal,  "vertex_z_ahcal.png");
    savePlot(h_z_muspec, "vertex_z_muspec.png");

    savePlot2D(h_xz_all,    "vertex_xz_all.png");
    savePlot2D(h_xz_3dcal,  "vertex_xz_3dcal.png");
    savePlot2D(h_xz_ecal,   "vertex_xz_ecal.png");
    savePlot2D(h_xz_ahcal,  "vertex_xz_ahcal.png");
    savePlot2D(h_xz_muspec, "vertex_xz_muspec.png");

    std::cout << "Plots saved in ./plots/" << std::endl;
    std::cout << "ROOT file saved as " << outputRootFile << std::endl;
*/
// ----------------- Plots with fiducial overlays -----------------

    gStyle->SetOptStat(1110);

    // 1D z: ALL + colored bands
    {
        TCanvas *c = new TCanvas("c_z_all", "z (all) with fiducial volumes", 900, 600);
        h_z_all->SetLineWidth(2);
        h_z_all->SetLineColor(kBlack);
        h_z_all->Draw();

        double ymax = h_z_all->GetMaximum() * 1.20;

        // TBoxes: [z_min, z_max] × [0, ymax]
        TBox *b_3dcal    = new TBox(Z_3DCAL_MIN_MM,    0., Z_3DCAL_MAX_MM,    ymax);
        TBox *b_ecal = new TBox(Z_ECAL_MIN_MM, 0., Z_ECAL_MAX_MM, ymax);
        TBox *b_ahcal    = new TBox(Z_AHCAL_MIN_MM,    0., Z_AHCAL_MAX_MM,    ymax);
        TBox *b_muspec   = new TBox(Z_MS_MIN_MM,       0., Z_MS_MAX_MM,       ymax);

        b_3dcal->SetFillColor(kYellow);
        b_3dcal->SetFillStyle(3004);
        b_3dcal->Draw("same");

        b_ecal->SetFillColor(kBlue-9);
        b_ecal->SetFillStyle(3005);
        b_ecal->Draw("same");

        b_ahcal->SetFillColor(kRed-9);
        b_ahcal->SetFillStyle(3004);
        b_ahcal->Draw("same");

        b_muspec->SetFillColor(kGreen-9);
        b_muspec->SetFillStyle(3005);
        b_muspec->Draw("same");

        // Draw histogram again on top of boxes
        h_z_all->Draw("same");

        TLegend *leg = new TLegend(0.15, 0.70, 0.45, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h_z_all,    "All interactions",   "l");
        leg->AddEntry(b_3dcal,    "3DCAL",              "f");
        leg->AddEntry(b_ecal, " ECAL",          "f");
        leg->AddEntry(b_ahcal,    "AHCAL (rear HCAL)",  "f");
        leg->AddEntry(b_muspec,   "Muon spectrometer",  "f");
        leg->Draw();

        c->SaveAs("plots/vertex_z_all_with_fiducials.png");
        delete c;
    }

    // 2D x–z ALL: vertical lines at boundaries
    {
        TCanvas *c = new TCanvas("c_xz_all", "x vs z (all) with fiducial lines", 900, 600);
        h_xz_all->Draw("COLZ");

        double yMin = h_xz_all->GetYaxis()->GetXmin();
        double yMax = h_xz_all->GetYaxis()->GetXmax();

        auto drawBoundary = [&](double z, int color, int style) {
            TLine *ln = new TLine(z, yMin, z, yMax);
            ln->SetLineColor(color);
            ln->SetLineWidth(2);
            ln->SetLineStyle(style);
            ln->Draw("same");
        };

        // Boundaries between regions
        drawBoundary(Z_3DCAL_MIN_MM,    kGray+1, 2);
        drawBoundary(Z_3DCAL_MAX_MM,    kGray+1, 2);
        drawBoundary(Z_ECAL_MIN_MM, kBlue+1, 2);
        drawBoundary(Z_ECAL_MAX_MM, kBlue+1, 2);
        drawBoundary(Z_AHCAL_MIN_MM,    kRed+1,  2);
        drawBoundary(Z_AHCAL_MAX_MM,    kRed+1,  2);
        drawBoundary(Z_MS_MIN_MM,       kGreen+2,2);
        drawBoundary(Z_MS_MAX_MM,       kGreen+2,2);

        // Legend (just a text block explaining colors)
        TLegend *leg = new TLegend(0.15, 0.70, 0.45, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry((TObject*)0, "Vertical lines:", "");
        leg->AddEntry((TObject*)0, "grey: 3DCAL",       "");
        leg->AddEntry((TObject*)0, "blue: rear ECAL",   "");
        leg->AddEntry((TObject*)0, "red : AHCAL",       "");
        leg->AddEntry((TObject*)0, "green: muSpec",     "");
        leg->Draw();

        c->SaveAs("plots/vertex_xz_all_with_fiducials.png");
        delete c;
    }

    // Simple region-only plots (no overlays needed)

    auto savePlot1D = [&](TH1 *h, const char *filename) {
        TCanvas *c = new TCanvas();
        h->Draw();
        std::string out = std::string("plots/") + filename;
        c->SaveAs(out.c_str());
        delete c;
    };

    auto savePlot2D = [&](TH2 *h, const char *filename) {
        TCanvas *c = new TCanvas();
        h->Draw("COLZ");
        std::string out = std::string("plots/") + filename;
        c->SaveAs(out.c_str());
        delete c;
    };

    // 1D z per-region
    savePlot1D(h_z_all,       "vertex_z_all.png");          // without bands
    savePlot1D(h_z_3dcal,     "vertex_z_3dcal.png");
    savePlot1D(h_z_ecal,      "vertex_z_ecal.png");
    savePlot1D(h_z_ahcal,     "vertex_z_ahcal.png");
    savePlot1D(h_z_muspec,    "vertex_z_muspec.png");

    // 2D x–z per-region
    savePlot2D(h_xz_all,      "vertex_xz_all.png");         // without lines
    savePlot2D(h_xz_3dcal,    "vertex_xz_3dcal.png");
    savePlot2D(h_xz_ecal,     "vertex_xz_ecal.png");
    savePlot2D(h_xz_ahcal,    "vertex_xz_ahcal.png");
    savePlot2D(h_xz_muspec,   "vertex_xz_muspec.png");

    // 2D x–y per-region
    savePlot2D(h_xy_all,      "vertex_xy_all.png");
    savePlot2D(h_xy_3dcal,    "vertex_xy_3dcal.png");
    savePlot2D(h_xy_ecal,     "vertex_xy_ecal.png");
    savePlot2D(h_xy_ahcal,    "vertex_xy_ahcal.png");
    savePlot2D(h_xy_muspec,   "vertex_xy_muspec.png");

    std::cout << "Plots saved in ./plots/" << std::endl;
    std::cout << "Done." << std::endl;
    return 0;
}
