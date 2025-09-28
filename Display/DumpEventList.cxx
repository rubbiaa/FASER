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
  const std::string reco_folder_path = "input_reco/";
  std::vector<std::string> reco_file_paths;

  const std::string& input_folder_path = "input/";
  std::string input_file_path;

  for (const auto& entry : std::filesystem::directory_iterator(reco_folder_path)) {
    reco_file_paths.push_back(entry.path().string());
  }
  std::cout << " number of recofiles " << reco_file_paths.size() << std::endl;
  //const size_t max_files = 100;
  //if (reco_file_paths.size() > max_files)
  //  reco_file_paths.resize(max_files);
  
  int imask = 0;
  if (mask_str == "nueCC") imask = 1;
  else if (mask_str == "numuCC") imask = 2;
  else if (mask_str == "nutauCC") imask = 3;
  else if (mask_str == "nuNC") imask = 4;
  else if (mask_str == "nuES") imask = 5;
  else if (mask_str == " ") imask = 0;
  else if (mask_str == "-") imask = 0;

  // Variables to store in the tree
  int eventNumber, totalMultiplicity, chargedMultiplicity, neutralMultiplicity, gammaMultiplicity, neutronMultiplicity;
  int CharmType, CCNC, neutrinoType;
  double neutrinoEnergy;
  double prim_vx, prim_vy, prim_vz;
  double TotalEnergy, VisibleEnergy, RearECalEnergy, RearHCalEnergy,RearMuCalEnergy;
  Int_t max_vertices = 100;
  Int_t n_vertices;
  Float_t v_x[max_vertices], v_y[max_vertices], v_z[max_vertices];
  Int_t v_ntrks[max_vertices]; // number of tracks in vertex
  Int_t t_closest_vertex;   // index of vertex closest to true primary vertex
  Int_t n_tktracks;
  Int_t n_pstracks;
  Int_t n_clusters;
  Float_t c_E1;    // energy most energetic cluster
  Float_t c_E1T;     // transverse energy most energetic cluster
  Float_t c_chi2_1;
  Float_t c_a_1;
  Float_t c_b_1;
  Float_t c_E2;    // energy 2nd most energetic cluster

  TDatabasePDG* pdgDB = TDatabasePDG::Instance();
  display->AddCustomNucleusParticles(); 

  for (const std::string& file_path : reco_file_paths) 
    {
      std::string base_name = std::filesystem::path(file_path).stem().string();
      std::cout << file_path << std::endl;

      TChain *event_tree = new TChain("RecoEvent","READ");
      event_tree->Add(file_path.c_str());
      
      Long_t nentries = event_tree->GetEntries();
      std::cout << "Number of entries " << nentries << std::endl;
      
      // TPORecoEvent class and histograms
      TPORecoEvent *fTPORecoEvent = nullptr; // &fTPORecoEvent;
      event_tree -> SetBranchAddress("TPORecoEvent", &fTPORecoEvent);
      
      //std::cout << "debuggg " << file_path << " " << nentries << std::endl;
      // process events
      int ievent = 0;
      //if(maxEvents == -1) maxEvents = nentries;
      int error = 0;

      TPOEvent myPOevent; // local copy; just a temporary PO event to store stats
      myPOevent.reset_stats();
      
      //std::cout << " max event  " << file_path << " " << maxEvents << std::endl;   

      while (error == 0 && ievent<nentries) 
	{
	  event_tree->GetEntry(ievent);
	  
	  if(ievent % 1000 == 0) {
	    std::cout << "Processing event " << ievent << std::endl;
	  }
	  
	  std::cout << ievent << " ËventNumber_reco: " << fTPORecoEvent->GetPOEvent()-> event_id << " ";
	  eventNumber = fTPORecoEvent->GetPOEvent()-> event_id;
	  
	  input_file_path = input_folder_path + "FASERG4-Tcalevent_" +std::to_string(runNumber)+"_"+std::to_string(eventNumber)+".root";
	  //std::cout << input_file_path << std::endl;
	  // Validate input file before processing
	  TFile* input_check = TFile::Open(input_file_path.c_str(), "READ");
	  if (!input_check || input_check->IsZombie()) {
	    std::cout << eventNumber << " no input: ";
	    //std::cerr << "Error: Cannot open or read input file: " << input_file_path << std::endl;
	    delete input_check;
	    ievent++;
	    continue;
	  }
	  input_check->Close();
	  delete input_check;
	  
	  std::cout << eventNumber << " ok " << std::endl;  

	  ievent++;
	}

    }
  std::cout << "Finished processing " << " files." << std::endl;
}





int main(int argc, char** argv)
{

  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <RunNumber> <MaxEvents> <Masks>" << std::endl;
    std::cout << "Masks: nueCC, numuCC, nutauCC, nuNC, nuES, -" << std::endl;
    return 1;
  }
  
  int runNumber = std::stoi(argv[1]);
  int maxEvents = std::stoi(argv[2]);
  std::string mask_str = argv[3];

  display::FaserCalDisplay* disp = new display::FaserCalDisplay();
  LoadAllRecoEvents(disp, runNumber, maxEvents, mask_str);
  delete disp;
  return 0;
}
