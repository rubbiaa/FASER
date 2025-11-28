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

void LoadAllRecoEvents(int runNumber, int maxEvents, std::string mask_str)
{
  const std::string reco_folder_path = "input_reco/";
  std::vector<std::string> reco_file_paths;

  for (const auto& entry : std::filesystem::directory_iterator(reco_folder_path)) {
    reco_file_paths.push_back(entry.path().string());
  }
  
  int imask = 0;
  if (mask_str == "nueCC") imask = 1;
  else if (mask_str == "numuCC") imask = 2;
  else if (mask_str == "nutauCC") imask = 3;
  else if (mask_str == "nuNC") imask = 4;
  else if (mask_str == "nuES") imask = 5;
  else if (mask_str == " ") imask = 0;
  else if (mask_str == "-") imask = 0;

  TFile* outputFile = new TFile("AnalysisData.root", "RECREATE");
  if (!outputFile || outputFile->IsZombie()) 
    {
      std::cerr << "Error: Could not create output ROOT file." << std::endl;
      return;
    }  
  TTree* tree = new TTree("AnalysisTree", "Storing FaserCalAnalyzer data");
 
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

  //TDatabasePDG* pdgDB = TDatabasePDG::Instance();
  //display->AddCustomNucleusParticles(); 


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
      
      // process events
      int ievent = 0;
      //if(maxEvents == -1) maxEvents = nentries;
      int error = 0;

      TPOEvent myPOevent; // local copy; just a temporary PO event to store stats
      myPOevent.reset_stats();
      
      while (error == 0 && ievent<nentries) 
	    {
	      event_tree->GetEntry(ievent);
    	  if(ievent % 1000 == 0) 
        {
	        std::cout << "Processing event " << ievent << std::endl;
	      }
        //fTPORecoEvent -> GetPOEvent()->dump_event();
        std::cout << "Truth Tau decay Mode: " << fTPORecoEvent -> GetPOEvent()->tau_decaymode << std::endl;
            tauDecayModeTruth = fTPORecoEvent -> GetPOEvent()->tau_decaymode;
	  
	      std::cout << ievent << " EventNumberReco: " << fTPORecoEvent->GetPOEvent()-> event_id << std::endl;
	      eventNumber = fTPORecoEvent->GetPOEvent()-> event_id;

    	  ROOT::Math::XYZVector true_vtx(fTPORecoEvent->GetPOEvent()->prim_vx.x(), fTPORecoEvent->GetPOEvent()->prim_vx.y(), fTPORecoEvent->GetPOEvent()->prim_vx.z());
        std::cout << "True Vertex: " << true_vtx.X() << ", " << true_vtx.Y() << ", " << true_vtx.Z() << std::endl;

	      VisibleEnergy = fTPORecoEvent->GetPOFullRecoEvent()->TotalEvis();
	      TotalEnergy = fTPORecoEvent->GetPOFullRecoEvent()->TotalET();
	      RearECalEnergy = fTPORecoEvent->rearCals.rearCalDeposit;
	      RearHCalEnergy = fTPORecoEvent->rearCals.rearHCalDeposit;
	      RearMuCalEnergy = fTPORecoEvent->rearCals.rearMuCalDeposit;
	  
        std::cout << VisibleEnergy << " "
            << TotalEnergy << " "
            << RearECalEnergy << " "
            << RearHCalEnergy << " "
            << RearMuCalEnergy << " "
            << std::endl;

        
 std::cout << " Inside Muon Spectrometer " << std::endl;
    for (int i = 0; i < fTPORecoEvent->fMuTracks.size(); ++i)
      {
        TMuTrack* muTrack = &fTPORecoEvent->fMuTracks[i];
        if (!muTrack) continue;
        std::cout << "Muon track info: "
                  << "  " << i << " " 
                  << muTrack->ftrackID << " "
                  << muTrack->fPDG << " "
                  << muTrack->fcharge << " "
                  << muTrack->fpos.size() << " "
                  << muTrack->fpx << " "
                  << muTrack->fpy << " "
                  << muTrack->fpz << " "
                  << muTrack->fp << " "
                  << muTrack->fchi2 << " "
                  << muTrack->fnDoF << " "
                  << muTrack->fpval << " "
                  << muTrack->fpErr << " "
                  << muTrack->fipErr << " "
                  << std::endl;
        for (size_t j=0; j<muTrack->fpos.size(); ++j)
        {
          const auto position = muTrack->fpos[j];
          std::cout << j << " Muon Spectrometer hit position " 
                    << muTrack->ftrackID << " "
                    << muTrack->layerID[j] << " "
                    << position.x() << " "
                    << position.y() << " "
                    << position.z() << " "
                    << std::endl;
        }
      }


	  ievent++;
	}

    }

  outputFile->cd();
  tree->Write();
  std::cout << "Finished processing " << " files." << std::endl;

// Print summary statistics
  std::cout << "\n========== ANALYSIS SUMMARY ==========" << std::endl;
  std::cout << "Tree entries: " << tree->GetEntries() << std::endl;
  std::cout << "Output file: AnalysisData.root" << std::endl;
  std::cout << "Branches include: Charm, Tau, and generic SLP information" << std::endl;
  outputFile->Close();

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

 //display::FaserCalDisplay* disp = new display::FaserCalDisplay();
  LoadAllRecoEvents(runNumber, maxEvents, mask_str);
  //delete disp;
  return 0;
}
