// Count neutrino interactions per detector region using TPOEvent
//
// Regions counted (in mm, along z):
//  - 3DCAL:       [-1150,  1150]
//  - rear ECAL:   [1205,   1550]
//  - AHCAL:       [1551,   2450]
//  - MuSpect:     [2451,   4350]
// Upstream "front target" veto: z < -1150 mm
//
// Usage:
//   CountGENIERegions.exe [-charmonly] [-tauCConly] <genieroot> <run>
//
// Notes:
//   - Uses TPOEvent as in the original converter, but does NOT write
//     any POEvent ROOT file. It only prints event counts per region.
//   - If -charmonly is used, only events with charm are counted.
//   - If -tauCConly is used, only events with a tau lepton are counted.

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <cstring>

#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TBranch.h>

#include "TPOEvent.hh"

bool charm_only = false;
bool tauCC_only = false;

// Geometry / regions (mm)
constexpr double Z_FRONT_TARGET_MAX_MM = -1150.0;  // upstream veto

// 3DCAL region (same as your original cuts)
constexpr double Z_3DCAL_MIN_MM       = -1150.0;
constexpr double Z_3DCAL_MAX_MM       =  1150.0;

// rear ECAL (from GDML)
constexpr double Z_REARECAL_MIN_MM    = 1205.0;
constexpr double Z_REARECAL_MAX_MM    = 1550.0;

// AHCAL / rear HCAL (from GDML)
constexpr double Z_AHCAL_MIN_MM       = 1551.0;
constexpr double Z_AHCAL_MAX_MM       = 2450.0;

// Muon spectrometer (from GDML)
constexpr double Z_MUSPEC_MIN_MM      = 2451.0;
constexpr double Z_MUSPEC_MAX_MM      = 4350.0;

// Per-region stats
struct RegionStats {
  long events      = 0;   // events in this region (after cuts)
  long nueCC       = 0;
  long numuCC      = 0;
  long nutauCC     = 0;
  long NC          = 0;   // all NC 
  long ES          = 0;   // elastic scattering on electrons
};

struct StatSnapshot {
  long nueCC   = 0;
  long numuCC  = 0;
  long nutauCC = 0;
  long NC      = 0;
  long ES      = 0;
};

// get the text printed by fTPOEvent.dump_stats()
StatSnapshot get_dump_stats_output(const std::string &s)
{
  StatSnapshot snap;
  std::istringstream iss(s);
  std::string tok;

  // Example printed lines:
  //  nueCC = 1778 numuCC = 24834 nutauCC = 0 NC = 8663
  //  ES = 70
  while (iss >> tok) {
    if (tok == "nueCC") {
      std::string eq; long v;
      if (iss >> eq >> v) snap.nueCC = v;
    } else if (tok == "numuCC") {
      std::string eq; long v;
      if (iss >> eq >> v) snap.numuCC = v;
    } else if (tok == "nutauCC") {
      std::string eq; long v;
      if (iss >> eq >> v) snap.nutauCC = v;
    } else if (tok == "NC") {
      std::string eq; long v;
      if (iss >> eq >> v) snap.NC = v;
    } else if (tok == "ES") {
      std::string eq; long v;
      if (iss >> eq >> v) snap.ES = v;
    }
  }
  return snap;
}

// call fTPOEvent.dump_stats() and capture its stdout into a string
StatSnapshot snapshot_stats(TPOEvent &ev)
{
  std::ostringstream buffer;
  std::streambuf *old_cout_buf = std::cout.rdbuf();
  std::cout.rdbuf(buffer.rdbuf());

  ev.dump_stats();  // prints something like: "nueCC = ... numuCC = ..."

  std::cout.rdbuf(old_cout_buf);

  return get_dump_stats_output(buffer.str());
}

void count_FASERMC_regions(int run_number, TTree *tree,
                           int min_event, int max_event, int mask /*unused*/)
{
  std::cout << "Counting events per region (via TPOEvent) ..." << std::endl;

  TPOEvent fTPOEvent;

  // FASER GENIE branches
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

  TBranch *b_vx;            //!
  TBranch *b_vy;            //!
  TBranch *b_vz;            //!
  TBranch *b_n;             //!
  TBranch *b_name;          //!
  TBranch *b_pdgc;          //!
  TBranch *b_status;        //!
  TBranch *b_firstMother;   //!
  TBranch *b_lastMother;    //!
  TBranch *b_firstDaughter; //!
  TBranch *b_lastDaughter;  //!
  TBranch *b_px;            //!
  TBranch *b_py;            //!
  TBranch *b_pz;            //!
  TBranch *b_E;             //!
  TBranch *b_m;             //!
  TBranch *b_M;             //!

  tree->SetBranchAddress("vx", &vx, &b_vx);
  tree->SetBranchAddress("vy", &vy, &b_vy);
  tree->SetBranchAddress("vz", &vz, &b_vz);
  tree->SetBranchAddress("n", &n, &b_n);
  tree->SetBranchAddress("name", &name, &b_name);
  tree->SetBranchAddress("pdgc", &pdgc, &b_pdgc);
  tree->SetBranchAddress("status", &status, &b_status);
  tree->SetBranchAddress("firstMother", &firstMother, &b_firstMother);
  tree->SetBranchAddress("lastMother", &lastMother, &b_lastMother);
  tree->SetBranchAddress("firstDaughter", &firstDaughter, &b_firstDaughter);
  tree->SetBranchAddress("lastDaughter", &lastDaughter, &b_lastDaughter);
  tree->SetBranchAddress("px", &px, &b_px);
  tree->SetBranchAddress("py", &py, &b_py);
  tree->SetBranchAddress("pz", &pz, &b_pz);
  tree->SetBranchAddress("E", &E, &b_E);
  tree->SetBranchAddress("m", &m, &b_m);
  tree->SetBranchAddress("M", &M, &b_M);

  // Counters
  long total_events_read        = 0;
  long total_events_after_vtx   = 0;
  long total_events_after_cuts  = 0;

  // Per-region stats
  RegionStats reg_3dcal;
  RegionStats reg_rearecal;
  RegionStats reg_ahcal;
  RegionStats reg_muspec;
  RegionStats reg_other;

  int evt_to_dump = 0;

  for (int event = min_event; event < max_event; ++event)
  {
    tree->GetEntry(event);
    ++total_events_read;

    if (event == min_event)
      fTPOEvent.reset_stats();

    if (event % 1000 == 0)
      std::cout << "Processing event " << event << " ..." << std::endl;

    fTPOEvent.clear_event();
    fTPOEvent.run_number = run_number;

    // Convert vertex from meters to mm
    fTPOEvent.setPrimaryVtx(vx * 1e3, vy * 1e3, vz * 1e3);
    fTPOEvent.use_GENIE_vtx = true; // tell FASERG4 to use this vtx

    double z_mm = fTPOEvent.prim_vx.z();

    // Upstream front target veto (same logic as before)
    if (z_mm < Z_FRONT_TARGET_MAX_MM)
      continue;

    ++total_events_after_vtx;

    // Determine region based on z
    bool in_3dcal    = (z_mm >= Z_3DCAL_MIN_MM    && z_mm < Z_3DCAL_MAX_MM);
    bool in_rearecal = (z_mm >= Z_REARECAL_MIN_MM && z_mm < Z_REARECAL_MAX_MM);
    bool in_ahcal    = (z_mm >= Z_AHCAL_MIN_MM    && z_mm < Z_AHCAL_MAX_MM);
    bool in_muspec   = (z_mm >= Z_MUSPEC_MIN_MM   && z_mm < Z_MUSPEC_MAX_MM);

    RegionStats *R = &reg_other;
    if (in_3dcal)    R = &reg_3dcal;
    else if (in_rearecal) R = &reg_rearecal;
    else if (in_ahcal)    R = &reg_ahcal;
    else if (in_muspec)   R = &reg_muspec;

    // Build PO list from GENIE,
    bool found_tau_lepton = false;

    for (int i = 0; i < n; ++i)
    {
      struct PO aPO;
      aPO.m_track_id = i;
      aPO.m_pdg_id   = pdgc->at(i);
      if (i == 1) {
        fTPOEvent.GENIE_vtx_name = name->at(i);
      }
      aPO.m_status = status->at(i);
      if (aPO.m_status == 0)
        aPO.m_status = 4;

      aPO.m_px     = px->at(i);
      aPO.m_py     = py->at(i);
      aPO.m_pz     = pz->at(i);
      aPO.m_energy = E->at(i);
      aPO.geanttrackID = -1;
      aPO.nparent = 0;

      int iMo = firstMother->at(i);
      if (iMo > -1) {
        aPO.nparent = 1;
        aPO.m_trackid_in_particle[0] = iMo;
      }

      // Only keep status 1 or 4, as before
      if (aPO.m_status == 1 || aPO.m_status == 4)
      {
        fTPOEvent.POs.push_back(aPO);
        fTPOEvent.perform_charmhadron_decay(aPO);
      }

      if (!found_tau_lepton && std::abs(aPO.m_pdg_id) == 15)
      {
        found_tau_lepton = true;
#ifdef _INCLUDE_PYTHIA_
        fTPOEvent.perform_taulepton_decay(aPO);
#endif
      }
    }

    // Kinematics and charm/tau identification
    fTPOEvent.kinematics_event();

    if (charm_only && !fTPOEvent.isCharmed())
      continue;
    if (tauCC_only && !found_tau_lepton)
      continue;

    // If we reached here, event is counted
    ++total_events_after_cuts;

    ++R->events;

    // Snapshot BEFORE
    StatSnapshot before = snapshot_stats(fTPOEvent);
    // Update stats
    fTPOEvent.update_stats();
    // Snapshot AFTER
    StatSnapshot after = snapshot_stats(fTPOEvent);

    StatSnapshot delta;
    delta.nueCC   = after.nueCC   - before.nueCC;
    delta.numuCC  = after.numuCC  - before.numuCC;
    delta.nutauCC = after.nutauCC - before.nutauCC;
    delta.NC      = after.NC      - before.NC;
    delta.ES      = after.ES      - before.ES;

    // Accumulate in region
    R->nueCC   += delta.nueCC;
    R->numuCC  += delta.numuCC;
    R->nutauCC += delta.nutauCC;
    R->NC      += delta.NC;
    R->ES      += delta.ES;

    if (evt_to_dump++ < 5) {
      std::cout << "---- Event " << event
                << "  z_vtx = " << z_mm << " mm" << std::endl;
      if (in_3dcal)    std::cout << "  -> in 3DCAL" << std::endl;
      if (in_rearecal) std::cout << "  -> in rear ECAL" << std::endl;
      if (in_ahcal)    std::cout << "  -> in AHCAL" << std::endl;
      if (in_muspec)   std::cout << "  -> in MuSpect" << std::endl;
    }
  }
  // Summary
  std::cout << "\n==========================================" << std::endl;
  std::cout << "Total events read from GENIE:          " << total_events_read << std::endl;
  std::cout << "After front-target veto (z > " << Z_FRONT_TARGET_MAX_MM << " mm): "
            << total_events_after_vtx << std::endl;
  std::cout << "After charm/tau cuts (if any):         " << total_events_after_cuts << std::endl;

  auto print_region = [](const std::string &name, double zmin, double zmax,
                         const RegionStats &R) {
    std::cout << "------------------------------------------" << std::endl;
    std::cout << name << " [" << zmin << ", " << zmax << "] mm:" << std::endl;
    std::cout << "  events  = " << R.events  << std::endl;
    std::cout << "  nueCC   = " << R.nueCC   << std::endl;
    std::cout << "  numuCC  = " << R.numuCC  << std::endl;
    std::cout << "  nutauCC = " << R.nutauCC << std::endl;
    std::cout << "  NC      = " << R.NC      << std::endl;
    std::cout << "  ES      = " << R.ES      << std::endl;
  };

  print_region("3DCAL",     Z_3DCAL_MIN_MM,    Z_3DCAL_MAX_MM,    reg_3dcal);
  print_region("rear ECAL", Z_REARECAL_MIN_MM, Z_REARECAL_MAX_MM, reg_rearecal);
  print_region("AHCAL",     Z_AHCAL_MIN_MM,    Z_AHCAL_MAX_MM,    reg_ahcal);
  print_region("MuSpect",   Z_MUSPEC_MIN_MM,   Z_MUSPEC_MAX_MM,   reg_muspec);
  print_region("OTHER",    -999999.,           999999.,           reg_other);

  std::cout << "==========================================" << std::endl;
  std::cout << "TPOEvent accumulated stats (all regions combined):" << std::endl;
  fTPOEvent.dump_stats();
}

int main(int argc, char **argv)
{
  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0]
              << " [-charmonly] [-tauCConly] <genieroot> <run>" << std::endl;
    std::cout << "  <genieroot>   Input FASER GENIE root file(s)" << std::endl;
    std::cout << "  <run>         Run number (stored in TPOEvent)" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -charmonly    Process only charm events" << std::endl;
    std::cout << "  -tauCConly    Process only tau CC events" << std::endl;
    return 1;
  }

  int iarg = 1;
  charm_only = false;
  tauCC_only = false;

  if (std::strcmp(argv[iarg], "-charmonly") == 0) {
    charm_only = true;
    iarg++;
  }
  if (std::strcmp(argv[iarg], "-tauCConly") == 0) {
    tauCC_only = true;
    iarg++;
  }

  if (iarg + 1 >= argc) {
    std::cerr << "Error: missing <genieroot> or <run> argument." << std::endl;
    return 1;
  }

  std::string rootinputString = argv[iarg++];
  std::string runString       = argv[iarg++];

  int run_number = 0;
  try {
    run_number = std::stoi(runString);
  } catch (const std::exception &e) {
    std::cerr << "Error converting run number: " << e.what() << std::endl;
    return 1;
  }

  int min_event = 0;
  int max_event = -1;
  int event_mask = 0; // unused, 

  TChain *tree = new TChain("gFaser");
  tree->Add(rootinputString.c_str());
  Long64_t n_entries = tree->GetEntries();
  if (max_event < 0 || max_event > n_entries)
    max_event = n_entries;

  std::cout << "Counting FASERMC events from: " << rootinputString << std::endl;
  std::cout << "Run number: " << run_number << std::endl;
  std::cout << "Entry range: [" << min_event << ", " << max_event << ")" << std::endl;

  count_FASERMC_regions(run_number, tree, min_event, max_event, event_mask);

  std::cout << "I'm done." << std::endl;
  return 0;
}
