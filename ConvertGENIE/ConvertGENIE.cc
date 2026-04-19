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

#include "TPOEvent.hh"

bool charm_only = false;
bool tauCC_only = false;

// Structure to hold detector geometry information
struct DetectorBoundaries {
  double z_3dcal_min;
  double z_3dcal_max;
  double z_ecal_min;
  double z_ecal_max;
  double z_ahcal_min;
  double z_ahcal_max;
  double tilt_angle;
};

// Parse GDML file to extract detector dimensions and tilt angle
bool parseGDML(const std::string& gdmlFile, DetectorBoundaries& boundaries) {
  std::ifstream file(gdmlFile);
  if (!file.is_open()) {
    std::cerr << "Error: Cannot open GDML file " << gdmlFile << std::endl;
    return false;
  }
  
  std::string line;
  std::map<std::string, double> box_sizes;
  std::map<std::string, double> box_positions;
  
  boundaries.tilt_angle = 0.0; // default
  
  // Read file line by line
  while (std::getline(file, line)) {
    // Look for box definitions with z dimension
    if (line.find("<box") != std::string::npos && line.find("lunit=\"mm\"") != std::string::npos) {
      size_t name_start = line.find("name=\"");
      size_t z_start = line.find("z=\"");
      if (name_start != std::string::npos && z_start != std::string::npos) {
        name_start += 6; // skip 'name="'
        size_t name_end = line.find("\"", name_start);
        std::string name = line.substr(name_start, name_end - name_start);
        
        z_start += 3; // skip 'z="'
        size_t z_end = line.find("\"", z_start);
        std::string z_str = line.substr(z_start, z_end - z_start);
        double z_size = std::stod(z_str);
        box_sizes[name] = z_size;
      }
    }
    
    // Look for physical volume positions
    if (line.find("<position") != std::string::npos && line.find("_pos\"") != std::string::npos) {
      size_t name_start = line.find("name=\"");
      size_t z_start = line.find("z=\"");
      if (name_start != std::string::npos && z_start != std::string::npos) {
        name_start += 6; // skip 'name="'
        size_t name_end = line.find("_pos\"", name_start);
        if (name_end != std::string::npos) {
          std::string name = line.substr(name_start, name_end - name_start);
          
          z_start += 3; // skip 'z="'
          size_t z_end = line.find("\"", z_start);
          std::string z_str = line.substr(z_start, z_end - z_start);
          double z_pos = std::stod(z_str);
          box_positions[name] = z_pos;
        }
      }
    }
    
    // Look for detector assembly rotation
    if (line.find("<rotation") != std::string::npos && 
        line.find("DetectorAssemblyPV") != std::string::npos && 
        line.find("_rot\"") != std::string::npos) {
      size_t y_start = line.find("y=\"");
      if (y_start != std::string::npos) {
        y_start += 3; // skip 'y="'
        size_t y_end = line.find("\"", y_start);
        std::string y_str = line.substr(y_start, y_end - y_start);
        boundaries.tilt_angle = std::stod(y_str);
      }
    }
  }
  
  file.close();
  
  // Calculate detector boundaries
  // 3DCAL: ContainerBox centered at origin
  double z_3dcal = -999;
  for (const auto& kv : box_sizes) {
    if (kv.first.find("ContainerBox") != std::string::npos) {
      z_3dcal = kv.second;
      break;
    }
  }
  
  // ECAL: ContainerEcal at specific position
  double z_ecal_size = -999;
  double z_ecal_pos = -999;
  for (const auto& kv : box_sizes) {
    if (kv.first.find("ContainerEcal") != std::string::npos) {
      z_ecal_size = kv.second;
      break;
    }
  }
  for (const auto& kv : box_positions) {
    if (kv.first.find("rearCal") != std::string::npos && 
        kv.first.find("Abs") == std::string::npos && 
        kv.first.find("Scint") == std::string::npos) {
      z_ecal_pos = kv.second;
      break;
    }
  }
  
  // AHCAL: ContainerHcal at specific position
  double z_hcal_size = -999;
  double z_hcal_pos = -999;
  for (const auto& kv : box_sizes) {
    if (kv.first.find("ContainerHcal") != std::string::npos) {
      z_hcal_size = kv.second;
      break;
    }
  }
  for (const auto& kv : box_positions) {
    if (kv.first.find("rearHCal") != std::string::npos) {
      z_hcal_pos = kv.second;
      break;
    }
  }
  
  // Validate and set boundaries
  if (z_3dcal > 0) {
    boundaries.z_3dcal_min = -z_3dcal / 2.0;
    boundaries.z_3dcal_max = z_3dcal / 2.0;
    std::cout << "GDML: 3DCAL: " << boundaries.z_3dcal_min << " to " << boundaries.z_3dcal_max << " mm" << std::endl;
  } else {
    return false;
  }
  
  if (z_ecal_size > 0 && z_ecal_pos != -999) {
    boundaries.z_ecal_min = z_ecal_pos - z_ecal_size / 2.0;
    boundaries.z_ecal_max = z_ecal_pos + z_ecal_size / 2.0;
    std::cout << "GDML: ECAL: " << boundaries.z_ecal_min << " to " << boundaries.z_ecal_max << " mm" << std::endl;
  } else {
    return false;
  }
  
  if (z_hcal_size > 0 && z_hcal_pos != -999) {
    boundaries.z_ahcal_min = z_hcal_pos - z_hcal_size / 2.0;
    boundaries.z_ahcal_max = z_hcal_pos + z_hcal_size / 2.0;
    std::cout << "GDML: AHCAL: " << boundaries.z_ahcal_min << " to " << boundaries.z_ahcal_max << " mm" << std::endl;
  } else {
    return false;
  }
  
  std::cout << "GDML: Tilt angle: " << boundaries.tilt_angle << " degrees" << std::endl;
  
  return true;
}

// Apply tilt correction to vertex coordinates
// For a detector tilted by angle (in degrees) around Y-axis
void apply_tilt_correction(double &x, double &z, double tilt_deg) {
  double tilt_rad = tilt_deg * M_PI / 180.0;
  double cos_tilt = cos(tilt_rad);
  double sin_tilt = sin(tilt_rad);
  
  // Rotate coordinates: tilt around Y-axis
  double x_new = x * cos_tilt + z * sin_tilt;
  double z_new = -x * sin_tilt + z * cos_tilt;
  
  x = x_new;
  z = z_new;
}

void convert_FASERMC(int run_number, TTree *tree, int min_event, int max_event,
                     std::string ROOTOutputFile, int mask, 
                     std::string detector, const DetectorBoundaries& boundaries)
{
  std::cout << "Converting events ..." << std::endl;
  std::cout << "Detector selection: " << detector << std::endl;
  std::cout << "Tilt angle: " << boundaries.tilt_angle << " degrees" << std::endl;
  
  // Define detector boundaries (in mm) from GDML or defaults
  double zmin = -1e9, zmax = 1e9; // default: accept all
  
  if (detector == "3DCAL") {
    zmin = boundaries.z_3dcal_min;
    zmax = boundaries.z_3dcal_max;
  } else if (detector == "ECAL") {
    zmin = boundaries.z_ecal_min;
    zmax = boundaries.z_ecal_max;
  } else if (detector == "AHCAL") {
    zmin = boundaries.z_ahcal_min;
    zmax = boundaries.z_ahcal_max;
  } else if (detector != "ALL") {
    throw std::runtime_error("Unknown detector: " + detector + ". Use 3DCAL, ECAL, AHCAL, or ALL");
  }

  TFile *m_rootFile = new TFile(ROOTOutputFile.c_str(), "RECREATE", "", 505); // last is the compression level
  if (!m_rootFile || !m_rootFile->IsOpen())
  {
    throw std::runtime_error("Could not create output ROOT file");
  }
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

  for (size_t event = min_event; event < max_event; event++)
  {

    tree->GetEntry(event);

    if (event == 0) fTPOEvent.reset_stats();
    if (event % 1000 == 0)
    {
      std::cout << "Processing event " << event << " ..." << std::endl;
    }

    fTPOEvent.clear_event();
    fTPOEvent.run_number = run_number;
    
    // Apply tilt correction if needed
    double vtx_x = vx * 1e3; // convert from meters to mm
    double vtx_y = vy * 1e3;
    double vtx_z = vz * 1e3;
    
    if (boundaries.tilt_angle != 0.0) {
      apply_tilt_correction(vtx_x, vtx_z, boundaries.tilt_angle);
    }
    
    // Filter by detector region
    if (vtx_z < zmin || vtx_z > zmax) continue;
    
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
    m_POEventTree->Fill();
    iseq++;
  }

  m_POEventTree->Write();
  m_rootFile->Close();
  std::cout << "Done saving..." << std::endl;

  fTPOEvent.dump_stats();
  std::cout << "Total number of events written " << iseq << std::endl;

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
    std::cout << "Default detector Z-ranges (mm) when no GDML provided:" << std::endl;
    std::cout << "  3DCAL:  -1205 to  1205" << std::endl;
    std::cout << "  ECAL:    1215 to  1650" << std::endl;
    std::cout << "  AHCAL:   1660 to  2775" << std::endl;
    std::cout << "  Tilt:    5.0 degrees" << std::endl;
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
  
  // Get detector boundaries from GDML or use defaults
  DetectorBoundaries boundaries;
  if (!gdmlFile.empty()) {
    std::cout << "Parsing GDML file: " << gdmlFile << std::endl;
    if (!parseGDML(gdmlFile, boundaries)) {
      std::cerr << "Error: Failed to parse GDML file" << std::endl;
      return 1;
    }
  } else {
    // Use default boundaries
    std::cout << "Using default detector boundaries and tilt angle" << std::endl;
    boundaries.z_3dcal_min = -1205.0;
    boundaries.z_3dcal_max = 1205.0;
    boundaries.z_ecal_min = 1215.0;
    boundaries.z_ecal_max = 1650.0;
    boundaries.z_ahcal_min = 1660.0;
    boundaries.z_ahcal_max = 2775.0;
    boundaries.tilt_angle = -5.0;
  }

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
  tree->Add(inputDirFiles.str().c_str());
  size_t n_entries = tree->GetEntries();
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
                  ROOTOutputFile.str(), event_mask, detector, boundaries);

  std::cout
      << "I'm done." << std::endl;
  return 0;
}
