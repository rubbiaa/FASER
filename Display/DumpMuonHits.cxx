#include "FaserCalDisplay.h" 
#include <TFile.h>
#include <TTree.h>
#include <TDatabasePDG.h>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <cmath>
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
  
  const size_t max_files = 100;
  if (reco_file_paths.size() > max_files)
    reco_file_paths.resize(max_files);
  
  int imask = 0;
  if (mask_str == "nueCC") imask = 1;
  else if (mask_str == "numuCC") imask = 2;
  else if (mask_str == "nutauCC") imask = 3;
  else if (mask_str == "nuNC") imask = 4;
  else if (mask_str == "nuES") imask = 5;
  else if (mask_str == " ") imask = 0;
  else if (mask_str == "-") imask = 0;

  // Create output ROOT file with TTree
  TFile* outfile = new TFile("muon_hits_analysis.root", "RECREATE");
  TTree* muonTree = new TTree("MuonHits", "Muon Hit Information in 3DCAL and MuonSpectrometer");
  
  // Branch variables
  Int_t eventNumber, trackID, primaryID, parentID, PDG;
  Bool_t isLeading;
  Int_t origin; // 1=primary, 2=charm, 3=light meson/other
  Int_t nu_type; // neutrino PDG code (14=numu, -14=numubar, 12=nue, -12=nuebar, 16=nutau, -16=nutaubar)
  Int_t CCNC; // 1=CC (charged current), 0=NC (neutral current)
  Double_t prim_vx, prim_vy, prim_vz;
  Double_t first_x, first_y, first_z;
  Double_t slope_xz, slope_yz;
  Int_t nhits_3DCAL, nhits_MuonSpec, nstations_MuonSpec;
  Bool_t escapes_3DCAL;
  Double_t last_z_3DCAL, first_z_MuonSpec;
  
  // Set up branches
  muonTree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  muonTree->Branch("trackID", &trackID, "trackID/I");
  muonTree->Branch("primaryID", &primaryID, "primaryID/I");
  muonTree->Branch("parentID", &parentID, "parentID/I");
  muonTree->Branch("PDG", &PDG, "PDG/I");
  muonTree->Branch("isLeading", &isLeading, "isLeading/O");
  muonTree->Branch("origin", &origin, "origin/I");
  muonTree->Branch("nu_type", &nu_type, "nu_type/I");
  muonTree->Branch("CCNC", &CCNC, "CCNC/I");
  muonTree->Branch("prim_vx", &prim_vx, "prim_vx/D");
  muonTree->Branch("prim_vy", &prim_vy, "prim_vy/D");
  muonTree->Branch("prim_vz", &prim_vz, "prim_vz/D");
  muonTree->Branch("first_x", &first_x, "first_x/D");
  muonTree->Branch("first_y", &first_y, "first_y/D");
  muonTree->Branch("first_z", &first_z, "first_z/D");
  muonTree->Branch("slope_xz", &slope_xz, "slope_xz/D");
  muonTree->Branch("slope_yz", &slope_yz, "slope_yz/D");
  muonTree->Branch("nhits_3DCAL", &nhits_3DCAL, "nhits_3DCAL/I");
  muonTree->Branch("nhits_MuonSpec", &nhits_MuonSpec, "nhits_MuonSpec/I");
  muonTree->Branch("nstations_MuonSpec", &nstations_MuonSpec, "nstations_MuonSpec/I");
  muonTree->Branch("escapes_3DCAL", &escapes_3DCAL, "escapes_3DCAL/O");
  muonTree->Branch("last_z_3DCAL", &last_z_3DCAL, "last_z_3DCAL/D");
  muonTree->Branch("first_z_MuonSpec", &first_z_MuonSpec, "first_z_MuonSpec/D");

  TDatabasePDG* pdgDB = TDatabasePDG::Instance();
  display->AddCustomNucleusParticles(); 

  TGeoManager::Import("../../GeomGDML/geometry_tilted_5degree.gdml");

  int total_muons = 0;
  int muons_in_3DCAL = 0;
  int muons_escaping_3DCAL = 0;
  int muons_in_MuonSpec = 0;
  int total_events_processed = 0;

  for (const std::string& file_path : reco_file_paths) 
    {
      // Check if we've already processed enough events
      if (maxEvents != -1 && total_events_processed >= maxEvents) {
        break;
      }
      
      std::string base_name = std::filesystem::path(file_path).stem().string();
      std::cout << file_path << std::endl;

      TChain *event_tree = new TChain("RecoEvent","READ");
      event_tree->Add(file_path.c_str());
      
      Long_t nentries = event_tree->GetEntries();
      std::cout << "Number of entries " << nentries << std::endl;
      
      TPORecoEvent *fTPORecoEvent = nullptr;
      event_tree->SetBranchAddress("TPORecoEvent", &fTPORecoEvent);
      
      int error = 0;
      
      for (int ievent = 0; ievent < nentries; ievent++) 
	{
	  // Check if we've reached the maximum number of events
	  if (maxEvents != -1 && total_events_processed >= maxEvents) {
	    break;
	  }
	  
	  event_tree->GetEntry(ievent);
	  
	  if(total_events_processed % 100 == 0) {
	    std::cout << "Processing event " << total_events_processed << std::endl;
	  }
	  
	  eventNumber = fTPORecoEvent->GetPOEvent()->event_id;
	  
	  input_file_path = input_folder_path + "FASERG4-Tcalevent_" + std::to_string(runNumber) + "_" + std::to_string(eventNumber) + ".root";
	  
	  // Validate input file
	  TFile* input_check = TFile::Open(input_file_path.c_str(), "READ");
	  if (!input_check || input_check->IsZombie()) {
	    std::cerr << "Error: Cannot open file: " << input_file_path << std::endl;
	    delete input_check;
	    total_events_processed++;
	    continue;
	  }
	  input_check->Close();
	  delete input_check;
	  
	  display->fTcalEvent = new TcalEvent();
	  display->POevent = new TPOEvent();
	  display->fTcalEvent->Load_event(input_folder_path, runNumber, eventNumber, imask, display->POevent);
	  
	  prim_vx = display->POevent->prim_vx.x();
	  prim_vy = display->POevent->prim_vx.y();
	  prim_vz = display->POevent->prim_vx.z();
	  
	  // Get neutrino information
	  nu_type = display->POevent->in_neutrino.m_pdg_id;
	  CCNC = display->POevent->isCC;

	  // Identify primary and charm muons
	  int PryMuGeantTrkID = -1;
	  int CharmParentTrkID = -1;
	  
	  for (size_t i = 0; i < display->POevent->n_particles(); i++)
	    {
	      struct PO& aPO = display->POevent->POs[i];
	      if (aPO.m_status == 1 && abs(aPO.m_pdg_id) == 13)
		{
		  PryMuGeantTrkID = aPO.geanttrackID;
		}
	      if (aPO.m_status == 1 && display->IsCharmed(aPO.m_pdg_id))
		{
		  CharmParentTrkID = aPO.geanttrackID;
		}
	    }
	  
	  // Map to store muon track information
	  std::map<int, std::vector<ROOT::Math::XYZVector>> muon_3DCAL_hits;
	  
	  // Process hits in 3DCAL
	  for (const auto& track : display->fTcalEvent->getfTracks())
	    {
	      if (abs(track->fPDG) != 13) continue; // Only muons
	      
	      size_t nhits = track->fhitIDs.size();
	      for (size_t i = 0; i < nhits; i++)
		{
		  long hittype = display->fTcalEvent->getChannelTypefromID(track->fhitIDs[i]);
		  if (hittype != 0) continue; // Only scintillator hits (3DCAL)
		  if (track->fEnergyDeposits[i] < 0.5) continue; // Energy threshold
		  
		  ROOT::Math::XYZVector position = display->fTcalEvent->getChannelXYZfromID(track->fhitIDs[i]);
		  muon_3DCAL_hits[track->ftrackID].push_back(position);
		}
	    }
	  
	  // Now process each muon track and fill the tree
	  for (const auto& track : display->fTcalEvent->getfTracks())
	    {
	      if (abs(track->fPDG) != 13) continue; // Only muons
	      
	      trackID = track->ftrackID;
	      primaryID = track->fprimaryID;
	      parentID = track->fparentID;
	      PDG = track->fPDG;
	      
	      // Determine if leading muon
	      isLeading = (PryMuGeantTrkID == primaryID);
	      
	      // Determine origin
	      if (PryMuGeantTrkID == primaryID) {
		origin = 1; // Primary muon
	      } else if (CharmParentTrkID == parentID) {
		origin = 2; // Charm muon
	      } else {
		origin = 3; // Light meson or other
	      }
	      
	      // Get 3DCAL hits for this track
	      auto& hits_3DCAL = muon_3DCAL_hits[trackID];
	      nhits_3DCAL = hits_3DCAL.size();
	      
	      if (nhits_3DCAL > 0) {
		total_muons++;
		muons_in_3DCAL++;
		
		// Calculate first hit position and last z
		first_x = hits_3DCAL[0].x();
		first_y = hits_3DCAL[0].y();
		first_z = hits_3DCAL[0].z();
		
		last_z_3DCAL = hits_3DCAL[0].z();
		for (const auto& hit : hits_3DCAL) {
		  if (hit.z() > last_z_3DCAL) last_z_3DCAL = hit.z();
		}
		
		// Calculate slope (using first and last hits if multiple hits)
		if (nhits_3DCAL >= 2) {
		  double dz = hits_3DCAL.back().z() - hits_3DCAL[0].z();
		  if (std::abs(dz) > 0.1) {
		    slope_xz = (hits_3DCAL.back().x() - hits_3DCAL[0].x()) / dz;
		    slope_yz = (hits_3DCAL.back().y() - hits_3DCAL[0].y()) / dz;
		  } else {
		    slope_xz = 0.0;
		    slope_yz = 0.0;
		  }
		} else {
		  slope_xz = 0.0;
		  slope_yz = 0.0;
		}
	      } else {
		// No hits in 3DCAL
		first_x = first_y = first_z = 0.0;
		slope_xz = slope_yz = 0.0;
		last_z_3DCAL = 0.0;
	      }
	      
	      // Check MuonSpectrometer (MuTag)
	      nhits_MuonSpec = 0;
	      nstations_MuonSpec = 0;
	      first_z_MuonSpec = -999.0;
	      std::set<int> stations_hit;
	      
	      for (const auto& mutrack : display->fTcalEvent->fMuTagTracks)
		{
		  if (mutrack->ftrackID == trackID)
		    {
		      nhits_MuonSpec = mutrack->pos.size();
		      
		      if (nhits_MuonSpec > 0) {
			first_z_MuonSpec = mutrack->pos[0].z();
			
			// Count unique stations (stationID = layerID / 4)
			for (size_t i = 0; i < mutrack->layerID.size(); i++) {
			  int stationID = mutrack->layerID[i] / 4;
			  stations_hit.insert(stationID);
			}
			nstations_MuonSpec = stations_hit.size();
		      }
		      break;
		    }
		}
	      
	      // Determine if muon escapes 3DCAL
	      escapes_3DCAL = (nhits_MuonSpec > 0);
	      
	      if (escapes_3DCAL) {
		muons_escaping_3DCAL++;
	      }
	      
	      if (nhits_MuonSpec > 0) {
		muons_in_MuonSpec++;
	      }
	      
	      // Fill tree only if muon has hits in either 3DCAL or MuonSpec
	      if (nhits_3DCAL > 0 || nhits_MuonSpec > 0) {
		muonTree->Fill();
		
		// Print summary
		std::cout << "Muon: Event=" << eventNumber 
			  << " TrackID=" << trackID
			  << " Origin=" << origin 
			  << " (1=primary,2=charm,3=other)"
			  << " Hits_3DCAL=" << nhits_3DCAL
			  << " Hits_MuonSpec=" << nhits_MuonSpec
			  << " Stations=" << nstations_MuonSpec
			  << " Escapes=" << (escapes_3DCAL ? "YES" : "NO")
			  << std::endl;
	      }
	    }
	  
	  total_events_processed++;
	  delete display->POevent;
	  delete display->fTcalEvent;
	}
      
      delete event_tree;
    }
  
  // Write tree and close file
  outfile->cd();
  muonTree->Write();
  outfile->Close();
  delete outfile;
  
  // Print summary statistics
  std::cout << "\n========================================" << std::endl;
  std::cout << "MUON STATISTICS SUMMARY" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Total events processed: " << total_events_processed << std::endl;
  std::cout << "Total muons with hits: " << total_muons << std::endl;
  std::cout << "Muons with hits in 3DCAL: " << muons_in_3DCAL << std::endl;
  std::cout << "Muons escaping 3DCAL (reaching MuonSpec): " << muons_escaping_3DCAL 
	    << " (" << (muons_in_3DCAL > 0 ? 100.0*muons_escaping_3DCAL/muons_in_3DCAL : 0) << "%)" << std::endl;
  std::cout << "Muons with hits in MuonSpectrometer: " << muons_in_MuonSpec << std::endl;
  std::cout << "========================================\n" << std::endl;
  std::cout << "Output saved to: muon_hits_analysis.root" << std::endl;
  std::cout << "TTree name: MuonHits" << std::endl;
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
