// Convert FASER GENIE human output to POevents
//
// A. Rubbia, August 2024
//

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>
#include <map>

#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TH1D.h>
#include <TH2D.h>

#include "TPOEvent.hh"

bool charm_only = false;
bool tauCC_only = false;

struct DetectorStats {
  long n_processed = 0;  // events that passed geometry check
  long n_written   = 0;  // events actually written (after charm/tau filters)
};

void load_geometry(std::string geometryFile){
    // Load the GDML geometry
    TGeoManager::Import(geometryFile.c_str());
    // Check if the geometry was loaded successfully
    if (!gGeoManager) {
        std::cerr << "Error: Could not load geometry from file " << geometryFile << std::endl;
        exit(1);
    }
    std::cout << "Geometry loaded from " << geometryFile << std::endl;
}

void convert_FASERMC(int run_number, TTree *tree, int min_event, int max_event,
                     std::string ROOTOutputFile, int mask, 
                     std::string detector)
{
  std::cout << "Converting events ..." << std::endl;
  std::cout << "Detector selection: " << detector << std::endl;
    
  bool want_3dcal = (detector == "3DCAL" || detector == "ALL");
  bool want_ecal = (detector == "ECAL" || detector == "ALL");
  bool want_ahcal = (detector == "AHCAL" || detector == "ALL");
  bool want_magnet = (detector == "MuonSpec" || detector == "ALL");
  
  TFile *m_rootFile = new TFile(ROOTOutputFile.c_str(), "RECREATE", "", 505); // last is the compression level
  if (!m_rootFile || !m_rootFile->IsOpen())
  {
    throw std::runtime_error("Could not create output ROOT file");
  }
  
    // Create output file
  TFile *outfile = new TFile("interaction_plots.root", "RECREATE");
  outfile->cd();

  // Create histograms
  TH1D *h_z_all = new TH1D("h_z_all", "All Interactions;Z [mm];Events", 400, 0, 8000);
  TH1D *h_x_all = new TH1D("h_x_all", "All Interactions;X [mm];Events", 200, -1000, 1000);
  TH1D *h_y_all = new TH1D("h_y_all", "All Interactions;Y [mm];Events", 200, -1000, 1000);
  TH2D *h_xy_all = new TH2D("h_xy_all", "All Interactions;X [mm];Y [mm]", 100, -1000, 1000, 100, -1000, 1000);
  TH2D *h_xz_all = new TH2D("h_xz_all", "All Interactions;Z [mm];X [mm]", 200, 0, 8000, 100, -1000, 1000);
  TH2D *h_yz_all = new TH2D("h_yz_all", "All Interactions;Z [mm];Y [mm]", 200, 0, 8000, 100, -1000, 1000);
  
  m_rootFile->cd();

  TPOEvent fTPOEvent;
  TPOEvent *branch_POEvent = &fTPOEvent;
  TTree *m_POEventTree = new TTree("POEvent", "POEvent");
  m_POEventTree->Branch("event", &branch_POEvent); // this should be named POEvent

  // FASER GENIE ntuple

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

  // List of branches
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

  int evt_to_dump = 0;
  size_t iseq = 0;
  
  // Initialize statistics counters
  fTPOEvent.reset_stats();

  // Per-detector event counters
  std::map<std::string, DetectorStats> det_stats;
  for (const auto &d : {"3DCAL", "ECAL", "AHCAL", "MuonSpec"})
    det_stats[d] = {};
  
  for (size_t event = min_event; event < max_event; event++)
  {

    tree->GetEntry(event);

    if (event % 1000 == 0)
    {
      std::cout << "Processing event " << event << " ..." << std::endl;
    }

    fTPOEvent.clear_event();
    fTPOEvent.run_number = run_number;

    double x = vx * 1e2; // convert from meters to cm
    double y = vy * 1e2;
    double z = vz * 1e2;
    TGeoNode *node = gGeoManager->FindNode(x, y, z);

    if (!node)
    {
      std::cout << "No volume contains this point.\n";
      exit(1);
    }

    std::string path = gGeoManager->GetPath();
   
    std::string volumename = node->GetName();
    bool in_3dcal = path.find("ContainerPlacement") != std::string::npos;
    bool in_ECAL = path.find("rearCal") != std::string::npos;
    bool in_HCAL = path.find("rearHCal") != std::string::npos;
    bool in_Magnet = path.find("MuonSpectrometer") != std::string::npos;
    bool in_elsewhere = path.find("World") != std::string::npos || path.find("DetectorAssemblyPV") != std::string::npos;
  
    if(!in_3dcal && !in_Magnet && !in_HCAL && !in_ECAL) {
      if(in_elsewhere) continue;
      std::cout << "Primary vertex coordinates (mm): (" << x << ", " << y << ", " << z << ")" << std::endl;
      std::cout << "Primary vertex is in volume: " << node->GetName() << std::endl;
      std::cout << "Path to volume: " << path << std::endl;
      continue; // skip events in front of the magnet
    }

    if(!((want_3dcal && in_3dcal) || (want_ecal && in_ECAL) || (want_ahcal && in_HCAL) || (want_magnet && in_Magnet))) {
      continue; // skip events not in the selected detector
    }
  
    // Tag which detector this vertex is in and count it
    std::string current_det = in_3dcal   ? "3DCAL"
                            : in_ECAL    ? "ECAL"
                            : in_HCAL    ? "AHCAL"
                            :              "MuonSpec";
    det_stats[current_det].n_processed++;
  
    fTPOEvent.setPrimaryVtx(vx * 1e3, vy * 1e3, vz * 1e3); // use original coordinates for G4
    fTPOEvent.use_GENIE_vtx = true;     // tell FASERG4 to use this vtx

    fTPOEvent.event_id = iseq;

    bool found_tau_lepton = false;

    for (size_t i = 0; i < n; i++)
    {
      struct PO aPO;
      aPO.m_track_id = i;
      aPO.m_pdg_id = pdgc->at(i);
      if(i==1) {
        fTPOEvent.GENIE_vtx_name = name->at(i);
      }
      aPO.m_status = status->at(i);
      if (aPO.m_status == 0)
        aPO.m_status = 4;
      aPO.m_px = px->at(i);
      aPO.m_py = py->at(i);
      aPO.m_pz = pz->at(i);
      aPO.m_energy = E->at(i);
      aPO.geanttrackID = -1;
      aPO.nparent = 0;
      int iMo = firstMother->at(i);
      if (iMo > -1)
      {
        aPO.nparent = 1;
        aPO.m_trackid_in_particle[0] = iMo;
      }

      if (aPO.m_status == 1 || aPO.m_status == 4)
      {
        fTPOEvent.POs.push_back(aPO);
        // decay charm hadron if necessary
        fTPOEvent.perform_charmhadron_decay(aPO);
      }

      if (!found_tau_lepton && abs(aPO.m_pdg_id) == 15)
      {
        found_tau_lepton = true;
#ifdef _INCLUDE_PYTHIA_
        fTPOEvent.perform_taulepton_decay(aPO);
#endif
      }

    }

    fTPOEvent.kinematics_event();

    // skip events with no charm 
    if(charm_only && !fTPOEvent.isCharmed()) continue;
    if(tauCC_only && !found_tau_lepton) continue;

    if (evt_to_dump++ < 20 && !found_tau_lepton)
    {
      fTPOEvent.dump_event();
    };
    fTPOEvent.update_stats();
    det_stats[current_det].n_written++;
    // Fill all interactions
    double vtx_x = vx * 1e3; // convert from meters to mm
    double vtx_y = vy * 1e3;
    double vtx_z = vz * 1e3;
    h_z_all->Fill(vtx_z);
    h_x_all->Fill(vtx_x);
    h_y_all->Fill(vtx_y);
    h_xy_all->Fill(vtx_x, vtx_y);
    h_xz_all->Fill(vtx_z, vtx_x);
    h_yz_all->Fill(vtx_z, vtx_y);
    m_POEventTree->Fill();
    iseq++;
  }

  m_POEventTree->Write();
  m_rootFile->Close();
  std::cout << "Done saving..." << std::endl;

  fTPOEvent.dump_stats();

  std::cout << "\n--- Per-detector event counts ---\n";
  std::cout << std::left
            << std::setw(12) << "Detector"
            << std::setw(14) << "Processed"
            << std::setw(14) << "Written" << "\n";
  std::cout << std::string(40, '-') << "\n";
  for (const auto &kv : det_stats) {
    std::cout << std::setw(12) << kv.first
              << std::setw(14) << kv.second.n_processed
              << std::setw(14) << kv.second.n_written << "\n";
  }
  std::cout << std::string(40, '-') << "\n";
  std::cout << "Total number of events written " << iseq << "\n";

 // Write all histograms
  outfile->cd(); 
  h_z_all->Write();
  h_x_all->Write();
  h_y_all->Write();
  h_xy_all->Write();
  h_xz_all->Write();
  h_yz_all->Write();

  // save all histograms as C macro
  h_z_all->SaveAs("z_distribution.C");
  h_x_all->SaveAs("x_distribution.C");
  h_y_all->SaveAs("y_distribution.C");
  h_xy_all->SaveAs("xy_distribution.C");
  h_xz_all->SaveAs("xz_distribution.C");
  h_yz_all->SaveAs("yz_distribution.C");
  outfile->Close();
}

int main(int argc, char **argv)
{
  // get the output file name as the first argument
  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0] << " [options] <genieroot> <run> [detector]" << std::endl;
    std::cout << std::endl;
    std::cout << "  <genieroot>                The input FASER GENIE root files" << std::endl;
    std::cout << "  <run>                      The output run number" << std::endl;
    std::cout << "  [detector]                 Detector selection: 3DCAL, ECAL, AHCAL, or ALL (default: ALL)" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -g <gdmlfile>              GDML geometry file (extracts detector boundaries and tilt)" << std::endl;
    std::cout << "  -charmonly                 Process only charm events " << std::endl;
    std::cout << "  -tauCConly                 Process only tau CC events " << std::endl;
    std::cout << std::endl;
    return 1;
  }

  int iarg = 1;
  
  // Parse options
  std::string gdmlFile = "";
  charm_only = false;
  tauCC_only = false;
  
  while (iarg < argc && argv[iarg][0] == '-') {
    if (strcmp(argv[iarg], "-g") == 0) {
      iarg++;
      if (iarg >= argc) {
        std::cerr << "Error: -g requires a GDML filename" << std::endl;
        return 1;
      }
      gdmlFile = argv[iarg++];
    } else if (strcmp(argv[iarg], "-charmonly") == 0) {
      charm_only = true;
      iarg++;
    } else if (strcmp(argv[iarg], "-tauCConly") == 0) {
      tauCC_only = true;
      iarg++;
    } else {
      std::cerr << "Unknown option: " << argv[iarg] << std::endl;
      return 1;
    }
  }
  
  if(argc - iarg < 2) {
    std::cerr << "Error: Missing required arguments" << std::endl;
    return 1;
  }
  
  std::string rootinputString = argv[iarg++];
  std::string runString = argv[iarg++];
  
  // Optional detector selection (default: ALL)
  std::string detector = "ALL";
  if(iarg < argc) {
    detector = argv[iarg++];
  }
  
  load_geometry(gdmlFile);

  int run_number;
  try
  {
    run_number = std::stoi(runString);
  }
  catch (const std::invalid_argument &e)
  {
    std::cerr << "Invalid argument for run: " << e.what() << std::endl;
    return 1;
  }
  catch (const std::out_of_range &e)
  {
    std::cerr << "Out of range for run: " << e.what() << std::endl;
    return 1;
  }

  int min_event = 0;
  int max_event = -1;
  int event_mask = 0;

  std::ostringstream inputDirFiles;
  inputDirFiles << rootinputString;

  TChain *tree = new TChain("gFaser");
  int nfiles = tree->Add(inputDirFiles.str().c_str());
  if (nfiles == 0) {
    std::cerr << "Error: No files found matching pattern: " << inputDirFiles.str() << std::endl;
    std::cerr << "Please check the file path and try again." << std::endl;
    return 1;
  }
  
  size_t n_entries = tree->GetEntries();
  std::cout << "Found " << nfiles << " file(s) with " << n_entries << " total entries" << std::endl;
  
  if (n_entries == 0) {
    std::cerr << "Error: No entries found in input file(s)" << std::endl;
    return 1;
  }
  
  if (max_event == -1)
  {
    max_event = n_entries;
  }

  std::ostringstream ROOTOutputFile;
  ROOTOutputFile << "FASERMC-PO-Run" << run_number << "-" << min_event << "_" << max_event;
  if (detector != "ALL") {
    ROOTOutputFile << "_" << detector;
  }
  if (event_mask > 0)
  {
    const char *mask = TPOEvent::DecodeEventMask(event_mask);
    ROOTOutputFile << "_" << mask;
  }
  ROOTOutputFile << ".root";

  std::cout << "Converting FASERMC from " << inputDirFiles.str() << std::endl;
  std::cout << "The output file is " << ROOTOutputFile.str() << std::endl;

  convert_FASERMC(run_number, tree, min_event, max_event,
                  ROOTOutputFile.str(), event_mask, detector);

  std::cout
      << "I'm done." << std::endl;
  return 0;
}
