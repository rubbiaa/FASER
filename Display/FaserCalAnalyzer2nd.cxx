#include "FaserCalDisplay.h" 
#include <TFile.h>
#include <TTree.h>
#include <TDatabasePDG.h>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <sstream>
#include "TChain.h"
#include <string>

#include "TcalEvent.hh"
#include "TPOEvent.hh"
#include "TPORecoEvent.hh"

void LoadAllRecoEvents(display::FaserCalDisplay* display, int runNumber, int maxEvents, std::string mask_str)
{
  
  const std::string& input_folder_path = "input/";
  std::string input_file_path;
  // Variables to store in the tree
  int eventNumber = 0;
  int totalMultiplicity = 0, chargedMultiplicity = 0, neutralMultiplicity = 0;
  int gammaMultiplicity = 0, neutronMultiplicity = 0;
  int CharmType = 0, CCNC = 0, neutrinoType = 0;
  double neutrinoEnergy = 0;
  double prim_vx = 0, prim_vy = 0, prim_vz = 0;
  double TotalEnergy = 0, VisibleEnergy = 0, RearECalEnergy = 0, RearHCalEnergy = 0, RearMuCalEnergy = 0;
  Float_t JetEne = 0, JetPt = 0, Ptmiss = 0, Evis_true = 0;
  Int_t tauDecayModeTruth = 0;

  std::vector<float> v_x, v_y, v_z;
  std::vector<int> v_ntrks;
  int t_closest_vertex = -1;
  int n_tktracks = 0;
  int n_pstracks = 0;
  int n_clusters = 0;
  int n_vertices = 0; 
  Float_t c_E1 = 0, c_E1T = 0, c_chi2_1 = 0, c_a_1 = 0, c_b_1 = 0, c_E2 = 0;

  TDatabasePDG* pdgDB = TDatabasePDG::Instance();
  display->AddCustomNucleusParticles(); 

TGeoManager::Import("../../GeomGDML/geometry.gdml");


  const std::string& folder_path = "input/";
  std::vector<std::string> file_paths;
  int cnt = 0;

  // Get all file paths in the input directory
  for (const auto& entry : std::filesystem::directory_iterator(folder_path))
    {
      file_paths.push_back(entry.path().string());
    }
  
  for (const std::string& file_path : file_paths)
    {
      std::cout << "Processing file: " << file_path << std::endl;
      cnt++;
      // Extract the base name from the file path
      std::string base_name = std::filesystem::path(file_path).stem().string();
      std::cout << "basename " << base_name << std::endl;    
      // Split the base_name using '_' as a delimiter
      std::istringstream ss(base_name);
      std::string token;
      std::vector<std::string> parts;
      while (std::getline(ss, token, '_'))
        {
	  parts.push_back(token);
        }
      if (parts.size() < 3)
        {
	  std::cerr << "Error: Invalid filename format, unable to parse: " << base_name << std::endl;
	  continue; // Skip this file if the format is incorrect
        }
      try {
	mask_str = (parts.size() == 3) ? "NoMask" : parts[3];
        runNumber = std::stoi(parts[1]);
        eventNumber = std::stoi(parts[2]);
	std::cout << cnt << " Event Number: " << eventNumber << ", Run Number: " << runNumber << ", Mask: " << mask_str << std::endl;
	int imask = 0;
	if (mask_str == "nueCC") imask = 1;
	else if (mask_str == "numuCC") imask = 2;
	else if (mask_str == "nutauCC") imask = 3;
	else if (mask_str == "nuNC") imask = 4;
	else if (mask_str == "nuES") imask = 5;
	
        //TPORecoEvent *fTPORecoEvent = nullptr; // &fTPORecoEvent;
	
	display->fTcalEvent = new TcalEvent();
	display->POevent = new TPOEvent();
	display->fTcalEvent->Load_event(input_folder_path, runNumber, eventNumber, imask, display->POevent);
	
        prim_vx = display->POevent->prim_vx.x();
        prim_vy = display->POevent->prim_vx.y();
        prim_vz = display->POevent->prim_vx.z();
        std::cout << prim_vx << " " << prim_vy << " " << prim_vz << std::endl; 
	//display->MuTagInfo();
	// ========== RESET ALL SHORT-LIVED PARTICLE VARIABLES ==========
	// reset SLP info
	display->fdecay_vx = display->fdecay_vy = display->fdecay_vz = 0.0;
	display->fdecayFlightLength = 0.0;
	display->fCharmEnergy = 0.0;
	display->fCharmCharge = 0;
	display->fnumChargedDaughters = -1;
	display->fnumNeutralDaughters = -1;
	display->fcharmname.clear();
	display->fdecayMode.clear();
	display->fCharm = -1;
	display->fCharmParentID = display->fCharmDaughterID = -1;
	
	// Reset tau variables
	display->ftau_vx = display->ftau_vy = display->ftau_vz = 0.0;
	display->ftauDecayFlightLength = 0.0;
	display->fTauEnergy = 0.0;
	display->fTauCharge = 0;
	display->fnumTauChargedDaughters = -1;
	display->fnumTauNeutralDaughters = -1;
	display->ftauname.clear();
	display->ftauDecayMode.clear();
	display->fTau = -1;
	display->fTauParentID = display->fTauDaughterID = -1;
	
	// Reset generic SLP variables
	display->fSLPParentIDs.clear();
	display->fSLPNames.clear();
	display->fSLPTypes.clear();
	totalMultiplicity = chargedMultiplicity = neutralMultiplicity = 0;
	gammaMultiplicity = neutronMultiplicity = 0;
	CharmType = neutrinoType = 0;
	//
	// display->CharmedEvent();
	
	// ========== ANALYZE SHORT-LIVED PARTICLES ==========
	bool hasSLP = display->ShortLivedParticleEvent(); // This analyzes both charm and tau
	
	if (hasSLP) {
	  std::cout << "Event " << eventNumber << " contains short-lived particles:" << std::endl;
	  for (size_t i = 0; i < display->fSLPTypes.size(); i++) {
	    std::string particleType = "";
	    switch(display->fSLPTypes[i]) {
            case 1: particleType = "Charm"; break;
            case 2: particleType = "Tau"; break;
            case 3: particleType = "Other SLP"; break;
            default: particleType = "Unknown"; break;
	    }
	    std::cout << "  - " << particleType << ": " << display->fSLPNames[i] 
		      << " (Parent ID: " << display->fSLPParentIDs[i] << ")" << std::endl;
	  }
	} // end of hasSLP
	
	  // particle counts
	totalMultiplicity = chargedMultiplicity = neutralMultiplicity = 0;
	gammaMultiplicity = neutronMultiplicity = 0;
	for (size_t i=0;i<display->POevent->n_particles(); ++i) {
	  struct PO &aPO = display->POevent->POs[i];
	  TParticlePDG* particle = pdgDB->GetParticle(aPO.m_pdg_id);
	  int charge = particle ? particle->Charge() : 0;
	  if (aPO.m_status == 1) {
	    if (display->IsCharmed(aPO.m_pdg_id)) CharmType = aPO.m_pdg_id;
	    totalMultiplicity++;
	    if (charge != 0) chargedMultiplicity++; else neutralMultiplicity++;
	    if (aPO.m_pdg_id == 22) gammaMultiplicity++;
	    if (aPO.m_pdg_id == 2112) neutronMultiplicity++;
	  }
	}
	// Get neutrino energy
	neutrinoEnergy = display->POevent->in_neutrino.m_energy;
	neutrinoType = display->POevent->in_neutrino.m_pdg_id;
	CCNC = display->POevent->isCC;	  
	
	JetEne = 0; JetPt = 0; Evis_true = 0; Ptmiss = 0;
	
	display->fTcalEvent -> fTPOEvent -> dump_event();
	
	tauDecayModeTruth = display->POevent->tau_decaymode;
	
	display->GetMuTagInfo();
	// ========== OUTPUT SHORT-LIVED PARTICLE SUMMARY ==========
	if (hasSLP) {
	  std::cout << "SLP Summary for Event " << eventNumber << ":" << std::endl;
	  if (display->fCharm != -1) {
	    std::cout << "  Charm: " << display->fcharmname << " (Type: " << display->fCharm 
		      << ", Prongs: " << display->fnumChargedDaughters << ")" << std::endl;
	  }
	  if (display->fTau != -1) {
	    std::cout << "  Tau: " << display->ftauname << " (Type: " << display->fTau 
		      << ", Prongs: " << display->fnumTauChargedDaughters << ")" << std::endl;
	  }
	}
      }
      catch (const std::exception& e) {
        std::cerr << "Error parsing filename or processing event: " << e.what() << std::endl;
        continue; // Skip this file on error
      }
    } // end of file_paths loop
}



int main(int argc, char** argv)
{

  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <RunNumber> <MaxEvents> <Masks>" << std::endl;
    std::cout << "Masks: nueCC, numuCC, nutauCC, nuNC, nuES, -" << std::endl;
    std::cout << "This version includes short-lived particle (charm and tau) analysis" << std::endl;

    return 1;
  }
  
  int runNumber = std::stoi(argv[1]);
  int maxEvents = std::stoi(argv[2]);
  std::string mask_str = argv[3];

  std::cout << "Starting FaserCalAnalyzer with Short-Lived Particle Analysis" << std::endl;
  std::cout << "Run: " << runNumber << ", Max Events: " << maxEvents << ", Mask: " << mask_str << std::endl;

  display::FaserCalDisplay* disp = new display::FaserCalDisplay();
  LoadAllRecoEvents(disp, runNumber, maxEvents, mask_str);
  delete disp;
  return 0;
}
