// Count neutrino interactions by type in FASER GENIE ntuples
//
// Usage:
//   countInteractions [-3dcal] [-ECAL] [-AHCAL] [-muSpec] <genieroot>
//
// If no region is specified, all interactions downstream of the front target
// (z >= Z_FRONT_TARGET_MAX_MM) are counted.
//


#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <cstring>

#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TBranch.h>

// ----------------- Geometry (in mm) – 3DCAL/ECAL/AHCAL/muSpec -----------------

// Upstream target / shield cut 
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

// ----------------- Region selection flags -----------------
bool select_3dcal  = false;
bool select_ecal   = false;
bool select_ahcal  = false;
bool select_muspec = false;

// ----------------- Simple container for counts -----------------
struct InteractionCounts {
  long cc  = 0;  // charged-current
  long nc  = 0;  // neutral-current
};

int main(int argc, char **argv)
{
  if (argc < 2)
  {
    std::cout << "Usage: " << argv[0]
              << " [-3dcal] [-ECAL] [-AHCAL] [-muSpec] <genieroot>" << std::endl;
    std::cout << std::endl;
    std::cout << "  <genieroot>     The input FASER GENIE ROOT file(s), "
                 "can use wildcards" << std::endl;
    std::cout << "Options (fiducial volume):" << std::endl;
    std::cout << "  -3dcal          Count interactions in 3DCAL region" << std::endl;
    std::cout << "  -ECAL           Count interactions in ECAL region" << std::endl;
    std::cout << "  -AHCAL          Count interactions in AHCAL (rear HCAL) region" << std::endl;
    std::cout << "  -muSpec         Count interactions in muon spectrometer region" << std::endl;
    std::cout << "If no region is specified, all regions downstream of the front target "
                 "are counted." << std::endl;
    return 1;
  }

  int iarg = 1;

  // Parse options
  while (iarg < argc && argv[iarg][0] == '-') {
    if (std::strcmp(argv[iarg], "-3dcal") == 0) {
      select_3dcal = true;
      iarg++;
      continue;
    }
    if (std::strcmp(argv[iarg], "-ECAL") == 0) {
      select_ecal = true;
      iarg++;
      continue;
    }
    if (std::strcmp(argv[iarg], "-AHCAL") == 0) {
      select_ahcal = true;
      iarg++;
      continue;
    }
    if (std::strcmp(argv[iarg], "-muSpec") == 0) {
      select_muspec = true;
      iarg++;
      continue;
    }

    std::cerr << "Unknown option: " << argv[iarg] << std::endl;
    return 1;
  }

  if (iarg >= argc) {
    std::cerr << "Error: missing <genieroot> argument." << std::endl;
    return 1;
  }

  std::string rootinputString = argv[iarg++];

  // ----------------------------------------------------------------------
  // Build the TChain
  // ----------------------------------------------------------------------
  TChain *tree = new TChain("gFaser");
  tree->Add(rootinputString.c_str());
  size_t n_entries = tree->GetEntries();

  if (n_entries == 0) {
    std::cerr << "No entries found in tree gFaser from " << rootinputString << std::endl;
    return 1;
  }

  std::cout << "Input GENIE file(s): " << rootinputString << std::endl;
  std::cout << "Total entries in gFaser: " << n_entries << std::endl;

  // ----------------------------------------------------------------------
  // Set branch addresses (as in your original converter)
  // ----------------------------------------------------------------------
  Double_t vx;
  Double_t vy;
  Double_t vz;
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

  TBranch *b_pdgc    = nullptr;
  TBranch *b_status  = nullptr;

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

  // ----------------------------------------------------------------------
  // Counters
  // ----------------------------------------------------------------------
  long total_events              = 0;  // after fiducial selection
  long total_events_all_regions  = 0;  // just downstream of front target

  long unknown_neutrino_flavour  = 0;
  long unknown_ccnc              = 0;

  // Map key = label string (e.g. "numu", "numubar"), value = counts
  std::map<std::string, InteractionCounts> counts;

  auto get_flavour_label = [](int nu_pdg) -> std::string {
    int apdg = std::abs(nu_pdg);
    bool is_anti = (nu_pdg < 0);

    std::string flav;
    if (apdg == 12) flav = "nue";
    else if (apdg == 14) flav = "numu";
    else if (apdg == 16) flav = "nutau";
    else flav = "unknown";

    if (flav != "unknown" && is_anti) flav += "bar";
    return flav;
  };

  // ----------------------------------------------------------------------
  // Event loop
  // ----------------------------------------------------------------------
  for (size_t ievt = 0; ievt < n_entries; ++ievt)
  {
    tree->GetEntry(ievt);

    if (ievt % 10000 == 0) {
      std::cout << "Processing event " << ievt << " / " << n_entries << "..." << std::endl;
    }

    // Convert z to mm (GENIE gives meters)
    double z_mm = vz * 1e3;

    // Skip upstream "front target" events (same logic as before)
    if (z_mm < Z_FRONT_TARGET_MAX_MM)
      continue;

    total_events_all_regions++;

    // Determine which region the vertex is in
    bool in3dcal   = (z_mm >= Z_3DCAL_MIN_MM && z_mm < Z_3DCAL_MAX_MM);
    bool inEcal    = (z_mm >= Z_ECAL_MIN_MM  && z_mm < Z_ECAL_MAX_MM);
    bool inAHCAL   = (z_mm >= Z_AHCAL_MIN_MM && z_mm < Z_AHCAL_MAX_MM);
    bool inMuonSpec= (z_mm >= Z_MS_MIN_MM    && z_mm < Z_MS_MAX_MM);

    bool keep = false;
    if (!select_3dcal && !select_ecal && !select_ahcal && !select_muspec) {
      // No region specified: keep anything downstream of front target
      keep = true;
    } else {
      if (select_3dcal  && in3dcal)    keep = true;
      if (select_ecal   && inEcal)     keep = true;
      if (select_ahcal  && inAHCAL)    keep = true;
      if (select_muspec && inMuonSpec) keep = true;
    }

    if (!keep) continue;

    total_events++;

    // ------------------------------------------------------------------
    // Identify primary neutrino and determine CC vs NC
    // ------------------------------------------------------------------
    int nu_pdg = 0;
    int nu_index = -1;

    for (int i = 0; i < n; ++i) {
      int pdg = pdgc->at(i);
      if (std::abs(pdg) == 12 || std::abs(pdg) == 14 || std::abs(pdg) == 16) {
        nu_pdg = pdg;
        nu_index = i;
        break; // take the first neutrino we find
      }
    }

    if (nu_index < 0) {
      unknown_neutrino_flavour++;
      continue;
    }

    std::string flav = get_flavour_label(nu_pdg);
    if (flav == "unknown") {
      unknown_neutrino_flavour++;
      continue;
    }

    // Look for a final-state charged lepton: e, mu, tau
    bool is_cc = false;

    for (int i = 0; i < n; ++i) {
      int pdg = pdgc->at(i);
      int apdg = std::abs(pdg);

      // Only consider final-state-like status codes: GENIE uses 1 for final state,
      // but in practice you may want to relax this. Start with status==1.
      int st = status->at(i);
      if (st != 1) continue;

      if (apdg == 11 || apdg == 13 || apdg == 15) {
        // We found an outgoing charged lepton -> CC
        is_cc = true;
        break;
      }
    }

    // Update counts
    InteractionCounts &c = counts[flav];
    if (is_cc) {
      c.cc++;
    } else {
      c.nc++;
    }
  }

  // ----------------------------------------------------------------------
  // Print results
  // ----------------------------------------------------------------------
  std::cout << "\n==========================================" << std::endl;
  std::cout << "Fiducial region selection:" << std::endl;
  if (!select_3dcal && !select_ecal && !select_ahcal && !select_muspec) {
    std::cout << "  All regions downstream of front target (z >= "
              << Z_FRONT_TARGET_MAX_MM << " mm)" << std::endl;
  } else {
    if (select_3dcal)
      std::cout << "  3DCAL   : " << Z_3DCAL_MIN_MM << " to "
                << Z_3DCAL_MAX_MM << " mm" << std::endl;
    if (select_ecal)
      std::cout << "  ECAL    : " << Z_ECAL_MIN_MM  << " to "
                << Z_ECAL_MAX_MM  << " mm" << std::endl;
    if (select_ahcal)
      std::cout << "  AHCAL   : " << Z_AHCAL_MIN_MM << " to "
                << Z_AHCAL_MAX_MM << " mm" << std::endl;
    if (select_muspec)
      std::cout << "  muSpec  : " << Z_MS_MIN_MM    << " to "
                << Z_MS_MAX_MM    << " mm" << std::endl;
  }

  std::cout << "==========================================" << std::endl;
  std::cout << "Total events: " << total_events_all_regions << std::endl;
  std::cout << "Total events after fiducial selection : " << total_events << std::endl;
  std::cout << "==========================================" << std::endl;

  for (const auto &kv : counts) {
    const std::string &flav = kv.first;
    const InteractionCounts &c = kv.second;
    long total_f = c.cc + c.nc;
    std::cout << flav << " interactions: " << total_f
              << "  (CC=" << c.cc << ", NC=" << c.nc << ")" << std::endl;
  }

  if (unknown_neutrino_flavour > 0) {
    std::cout << "Unknown neutrino flavour events: " << unknown_neutrino_flavour << std::endl;
  }
  if (unknown_ccnc > 0) {
    std::cout << "Unknown CC/NC classification events: " << unknown_ccnc << std::endl;
  }

  std::cout << "==========================================" << std::endl;
  std::cout << "Done." << std::endl;

  return 0;
}
