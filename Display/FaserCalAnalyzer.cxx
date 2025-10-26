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
  int t_closest_vertex = -1;
  int n_tktracks = 0;
  int n_pstracks = 0;
  int n_clusters = 0;
  int n_vertices = 0; 
  Float_t c_E1 = 0, c_E1T = 0, c_chi2_1 = 0, c_a_1 = 0, c_b_1 = 0, c_E2 = 0;

  /*
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

  Float_t JetEne, JetPt; // true jet energy and PT
  Float_t Ptmiss, Evis_true; // true Ptmiss and Evis

  Int_t tauDecayModeTruth; // Truth tau decay mode from POEvent
  */
  // Basic event information
  tree->Branch("Run", &runNumber, "Run/I");
  tree->Branch("Event", &eventNumber, "Event/I");
  
  tree->Branch("Vtx_pry",&prim_vx,"Vtx_pry/D");
  tree->Branch("Vty_pry",&prim_vy,"Vty_pry/D");
  tree->Branch("Vtz_pry",&prim_vz,"Vtz_pry/D");

  tree->Branch("TotalMultiplicity", &totalMultiplicity, "TotalMultiplicity/I");
  tree->Branch("ChargedMultiplicity", &chargedMultiplicity, "ChargedMultiplicity/I");
  tree->Branch("NeutralMultiplicity", &neutralMultiplicity, "NeutralMultiplicity/I");
  tree->Branch("GammaMultiplicity", &gammaMultiplicity, "GammaMultiplicity/I");
  tree->Branch("NeutronMultiplicity", &neutronMultiplicity, "NeutronMultiplicity/I");

  tree->Branch("NeutrinoEnergy", &neutrinoEnergy, "NeutrinoEnergy/D");
  tree->Branch("CCNC", &CCNC, "CCNC/I");
  tree->Branch("NeutrinoType", &neutrinoType, "NeutrinoType/I");

  tree->Branch("JetEne",&JetEne,"JetEne/F");
  tree->Branch("JetPt",&JetPt,"JetPt/F");
  tree->Branch("Evis_true",&Evis_true,"Evis_true/F");
  tree->Branch("Ptmiss",&Ptmiss,"Ptmiss/F");  

  //========== CHARM PARTICLE INFORMATION ==========
  tree->Branch("CharmType", &CharmType, "CharmType/I");
  tree->Branch("CharmProngs", &display->fnumChargedDaughters, "CharmProngs/I");
  tree->Branch("CharmDecayMode", &display->fdecayMode);
  tree->Branch("CharmDecay", &display->fCharm);
  tree->Branch("CharmParticle", &display->fcharmname);
  tree->Branch("CharmEnergy", &display->fCharmEnergy);

  tree->Branch("Vtx_dcy",&display->fdecay_vx,"Vtx_dcy/D");
  tree->Branch("Vty_dcy",&display->fdecay_vy,"Vty_dcy/D");
  tree->Branch("Vtz_dcy",&display->fdecay_vz,"Vtz_dcy/D");
  tree->Branch("FlightLength", &display->fdecayFlightLength, "FlightLength/D");

  // ========== TAU LEPTON INFORMATION ==========
  tree->Branch("TauType", &display->fTau, "TauType/I");
  tree->Branch("TauProngs", &display->fnumTauChargedDaughters, "TauProngs/I");
  tree->Branch("TauDecayMode", &display->ftauDecayMode);
  tree->Branch("TauDecayModeTruth", &tauDecayModeTruth);
  tree->Branch("TauEnergy", &display->fTauEnergy, "TauEnergy/D");
  tree->Branch("TauVtx_dcy", &display->ftau_vx, "TauVtx_dcy/D");
  tree->Branch("TauVty_dcy", &display->ftau_vy, "TauVty_dcy/D");
  tree->Branch("TauVtz_dcy", &display->ftau_vz, "TauVtz_dcy/D");
  tree->Branch("TauFlightLength", &display->ftauDecayFlightLength, "TauFlightLength/D");
 
  // ========== GENERIC SHORT-LIVED PARTICLE INFORMATION ==========
  tree->Branch("NumSLPs", &display->fSLPParentIDs); // Number of short-lived particles
  tree->Branch("SLPTypes", &display->fSLPTypes);     // Types: 1=charm, 2=tau, 3=other
  tree->Branch("SLPNames", &display->fSLPNames);     // Particle names

  // ========== VERTEX AND TRACKING INFORMATION ==========
  /*tree->Branch("n_vertices", &n_vertices);
  tree->Branch("v_x", &v_x, "v_x[n_vertices]/F");
  tree->Branch("v_y", &v_y, "v_y[n_vertices]/F");
  tree->Branch("v_z", &v_z, "v_z[n_vertices]/F");
  tree->Branch("v_ntrks", &v_ntrks, "v_ntrks[n_vertices]/I");
  tree->Branch("t_closest_vertex",&t_closest_vertex);
  tree->Branch("n_tktracks",&n_tktracks);
  tree->Branch("n_pstracks",&n_pstracks);
  tree->Branch("n_clusters",&n_clusters);
  tree->Branch("c_E1", &c_E1);
  tree->Branch("c_E1T", &c_E1T);
  tree->Branch("c_chi2_1", &c_chi2_1);
  tree->Branch("c_a_1", &c_a_1);
  tree->Branch("c_b_1", &c_b_1);
  tree->Branch("c_E2", &c_E2);
*/
  tree->Branch("n_vertices", &n_vertices, "n_vertices/I");
  tree->Branch("v_x", "std::vector<float>", &v_x);
    tree->Branch("v_y", "std::vector<float>", &v_y);
    tree->Branch("v_z", "std::vector<float>", &v_z);
    tree->Branch("v_ntrks", "std::vector<int>", &v_ntrks);
  tree->Branch("t_closest_vertex", &t_closest_vertex, "t_closest_vertex/I");
  tree->Branch("n_tktracks", &n_tktracks, "n_tktracks/I");
  tree->Branch("n_pstracks", &n_pstracks, "n_pstracks/I");
  tree->Branch("n_clusters", &n_clusters, "n_clusters/I");
  tree->Branch("c_E1", &c_E1, "c_E1/F");
  tree->Branch("c_E1T", &c_E1T, "c_E1T/F");
  tree->Branch("c_chi2_1", &c_chi2_1, "c_chi2_1/F");
  tree->Branch("c_a_1", &c_a_1, "c_a_1/F");
  tree->Branch("c_b_1", &c_b_1, "c_b_1/F");
  tree->Branch("c_E2", &c_E2, "c_E2/F");

  // ========== ENERGY INFORMATION ==========
  tree->Branch("VisibleEnergy", &VisibleEnergy,"VisibleEnergy/D");
  tree->Branch("TotalEnergy", &TotalEnergy,"TotalEnergy/D");
  tree->Branch("RearECalEnergy", &RearECalEnergy,"RearECalEnergy/D");
  tree->Branch("RearHCalEnergy", &RearHCalEnergy,"RearHCalEnergy/D");
  tree->Branch("RearMuCalEnergy", &RearMuCalEnergy,"RearMuCalEnergy/D");

  // ========== MUON TAGGING INFORMATION ==========
  tree->Branch("MuTag_alltrk", &display->fmuTag_alltrk);
  tree->Branch("MuTag_mom", &display->fmuTag_p);
  tree->Branch("MuTag_px", &display->fmuTag_px);
  tree->Branch("MuTag_py", &display->fmuTag_py);
  tree->Branch("MuTag_pz", &display->fmuTag_pz);
  tree->Branch("MuTag_ene", &display->fmuTag_E);
  tree->Branch("MuTag_dist", &display->fmuTag_dist);
  tree->Branch("MuTag_muon", &display->fmuTag_muon);
  tree->Branch("MuTag_muonSign", &display->fmuTag_muonSign);
  tree->Branch("MuTag_id", &display->fmuTag_id);
  tree->Branch("PryMu_id", &display->fPryMuonID);
  tree->Branch("CharmMuonDau_id", &display->fCharmDaughterID);

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
      
      // process events
      int ievent = 0;
      //if(maxEvents == -1) maxEvents = nentries;
      int error = 0;

      TPOEvent myPOevent; // local copy; just a temporary PO event to store stats
      myPOevent.reset_stats();
      
      while (error == 0 && ievent<nentries) 
	{
	  event_tree->GetEntry(ievent);
	  
	  if(ievent % 1000 == 0) {
	    std::cout << "Processing event " << ievent << std::endl;
	  }
	  // bool dump_event_cout = (ievent < 20);
	  //if(dump_event_cout) {
            fTPORecoEvent -> GetPOEvent()->dump_event();
            std::cout << "Truth Tau decay Mode: " << fTPORecoEvent -> GetPOEvent()->tau_decaymode << std::endl;
            tauDecayModeTruth = fTPORecoEvent -> GetPOEvent()->tau_decaymode;
            //}
	  
	  std::cout << ievent << " EventNumberReco: " << fTPORecoEvent->GetPOEvent()-> event_id << std::endl;
	  eventNumber = fTPORecoEvent->GetPOEvent()-> event_id;
	  // Clear vectors for new event
            v_x.clear();
            v_y.clear();
            v_z.clear();
            v_ntrks.clear();

            // Fill vertex information
            n_vertices = fTPORecoEvent->fTKVertices.size();
            for (int i = 0; i < n_vertices; ++i) {
                v_x.push_back(fTPORecoEvent->fTKVertices[i].position.x());
                v_y.push_back(fTPORecoEvent->fTKVertices[i].position.y());
                v_z.push_back(fTPORecoEvent->fTKVertices[i].position.z());
                v_ntrks.push_back(fTPORecoEvent->fTKVertices[i].trackIDs.size());
            }

	  // find vertex closest to true primary vertex
	  double min_dist = 1e9;
	  int min_index = -1;
	  ROOT::Math::XYZVector true_vtx(fTPORecoEvent->GetPOEvent()->prim_vx.x(), fTPORecoEvent->GetPOEvent()->prim_vx.y(), fTPORecoEvent->GetPOEvent()->prim_vx.z());
    for (int i=0;i<n_vertices;++i) {
        double dist = (fTPORecoEvent->fTKVertices[i].position - true_vtx).R();
        if (dist < min_dist) { min_dist = dist; t_closest_vertex = i; }
      }
	  // store number of tracks
	  n_tktracks = fTPORecoEvent->fTKTracks.size();
	  // store number of PS tracks
	  n_pstracks = fTPORecoEvent->fTKTracks.size();
	  // store number of PS clusters
	  n_clusters = fTPORecoEvent->n_psclustersX();
	  if(n_clusters>0) {
            c_E1 = fTPORecoEvent->PSClustersX[0].rawenergy/1e3;   // convert to GeV
            ROOT::Math::XYZVector dir = fTPORecoEvent->PSClustersX[0].cog-fTPORecoEvent->PSClustersX[0].vtx;
            c_E1T = fTPORecoEvent->PSClustersX[0].rawenergy/1e3*sqrt(dir.Unit().Perp2());
            c_chi2_1 = fTPORecoEvent->PSClustersX[0].longenergyprofile.chi2_per_ndf;
            if (std::isnan(c_chi2_1)) {
	      c_chi2_1 = -1.0;
            }
            c_a_1 = fTPORecoEvent->PSClustersX[0].longenergyprofile.a;
            c_b_1 = fTPORecoEvent->PSClustersX[0].longenergyprofile.b;
	  } else {
            c_E1 = 0;
            c_chi2_1 = -999;
            c_a_1 = c_b_1 = -999;
	  }
	  if(n_clusters>1) {
            c_E2 = fTPORecoEvent->PSClustersX[1].rawenergy/1e3;   // convert to GeV
	  } else {
            c_E2 = 0;
	  }
	  
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
	 
	  
	   input_file_path = input_folder_path + "FASERG4-Tcalevent_" +std::to_string(runNumber)+"_"+std::to_string(eventNumber)+".root";
	   std::cout << input_file_path << std::endl;
	   // Validate input file before processing
	   TFile* input_check = TFile::Open(input_file_path.c_str(), "READ");
	   if (!input_check || input_check->IsZombie()) {
	     std::cerr << "Error: Cannot open or read input file: " << input_file_path << std::endl;
	     delete input_check;
	     ievent++;
	     continue;
	   }
	   input_check->Close();
	   delete input_check;
	   
	   std::cout << ievent << " EventNumberTCaLEvent: " << fTPORecoEvent->GetPOEvent()-> event_id << " " << eventNumber << std::endl;

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
      }

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
    TVector3 jet = TVector3(fTPORecoEvent->GetPOEvent()->jetpx, fTPORecoEvent->GetPOEvent()->jetpy, fTPORecoEvent->GetPOEvent()->jetpz);
	  JetEne = jet.Mag();
	  JetPt = jet.Pt();
	  Evis_true = fTPORecoEvent->GetPOEvent()->Evis;
	  Ptmiss = fTPORecoEvent->GetPOEvent()->ptmiss;
	  tauDecayModeTruth = fTPORecoEvent->GetPOEvent()->tau_decaymode;

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

	  tree->Fill();

	  ievent++;
	  delete display->POevent;
	  delete display->fTcalEvent;
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

  display::FaserCalDisplay* disp = new display::FaserCalDisplay();
  LoadAllRecoEvents(disp, runNumber, maxEvents, mask_str);
  delete disp;
  return 0;
}
