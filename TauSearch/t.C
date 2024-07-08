// FASER kinematical analysis
// - t.C : code to convert FASER ntuple into event summary tuples
// A. Rubbia May 2024
//

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <map>

#include <cmath>
#include <random>

#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>

// #include "Pythia8/Pythia.h"

#include "faserntuplib.hh"

struct EVENT event;

TFile *tuple_file;
TTree *event_tree;

std::ofstream outFile("error.txt");

void create_histos(std::string outputFile) {
  // Open a ROOT file for writing
  tuple_file = new TFile(outputFile.c_str(), "RECREATE");
  event_tree = new TTree("event_tree", "Event data");

  event_tree->Branch("run_number", &event.run_number);
  event_tree->Branch("event_id", &event.event_id);
  event_tree->Branch("isCC", &event.isCC);
  event_tree->Branch("isES", &event.isES);
  event_tree->Branch("istau", &event.istau);
  event_tree->Branch("tau_decaymode",&event.tau_decaymode);
  event_tree->Branch("n_particles", &event.n_particles);
  event_tree->Branch("n_charged", &event.n_charged);
  event_tree->Branch("prim_vx", &event.prim_vx[0]);
  event_tree->Branch("prim_vy", &event.prim_vx[1]);
  event_tree->Branch("prim_vz", &event.prim_vx[2]);

  event_tree->Branch("in_lepton_pdgid", &event.in_neutrino.m_pdg_id);
  event_tree->Branch("vis_spx", &event.vis_spx);
  event_tree->Branch("vis_spy", &event.vis_spy);
  event_tree->Branch("vis_spz", &event.vis_spz);
  event_tree->Branch("jetpx", &event.jetpx);
  event_tree->Branch("jetpy", &event.jetpy);
  event_tree->Branch("jetpz", &event.jetpz);
  event_tree->Branch("tauvis_px", &event.tauvis_px);
  event_tree->Branch("tauvis_py", &event.tauvis_py);
  event_tree->Branch("tauvis_pz", &event.tauvis_pz);
  event_tree->Branch("tautracklength", &event.tautracklength);
  event_tree->Branch("Evis", &event.Evis);
  event_tree->Branch("ptmiss", &event.ptmiss);
  event_tree->Branch("cost", &event.cost);
  event_tree->Branch("cosf", &event.cosf);

  // tau->1pi candidate
  event_tree->Branch("taupi_cand[3]", event.taupi_cand);
  // tau->rho candidate
  event_tree->Branch("taurho_cand[3]", event.taurho_cand);
};

void fill_histos() {
  // Fill the tree
  event_tree->Fill();
}

void close_histos() {
  tuple_file->cd();
  event_tree->Write();
  tuple_file->Close();
}


void readDataCards(const std::string& filename, std::map<std::string, std::string>& config) {
  std::cout << "Reading datacards " << filename << std::endl;
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return;
  }
  
  std::string line;
  while (std::getline(file, line)) {
    // Skip comment lines
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream iss(line);
    std::cout << line << std::endl;
    std::string key, value;
    if (!(iss >> key >> value)) {
      std::cerr << "Error: Invalid line format: " << line << std::endl;
      continue;
    }
    config[key] = value;
  }
  file.close();
}



void t() {

  std::map<std::string, std::string> config;
  readDataCards("datacards.txt", config);

  std::string inputFile = config["input"];
  std::string outputFile = config["output"];

  // Open the ROOT file

 TChain *tree = new TChain("m_NuMCTruth_tree");
 // tree->Add("./sim/FaserMC-MC22_Genie_all_6invab-200006-00000-s0010-NTUP.root");
 tree->Add(inputFile.c_str());
 
 auto it = config.find("input2");
 if (it != config.end()) {
   tree->Add(it->second.c_str());
 }
 auto it3 = config.find("input3");
 if (it3 != config.end()) {
   tree->Add(it3->second.c_str());
 }

// Set up variables to hold the data and link them to the tree branches
Int_t m_runnumber, m_event_id_MC, m_track_id, m_pdg_id, m_num_in_particle, m_num_out_particle;
Double_t m_px, m_py, m_pz, m_energy, m_kinetic_energy, m_mass;
Float_t m_vx_prod, m_vy_prod, m_vz_prod, m_vx_decay, m_vy_decay, m_vz_decay;
std::vector<int> *m_pdg_in_particle = nullptr, *m_pdg_out_particle = nullptr;
std::vector<int> *m_trackid_in_particle = nullptr, *m_trackid_out_particle = nullptr;
Int_t m_status;

// Branch linking
// all energies are in MeV
tree->SetBranchAddress("m_runnumber", &m_runnumber);
tree->SetBranchAddress("m_event_id_MC", &m_event_id_MC);
tree->SetBranchAddress("m_track_id", &m_track_id);
tree->SetBranchAddress("m_pdg_id", &m_pdg_id);
tree->SetBranchAddress("m_px", &m_px);
tree->SetBranchAddress("m_py", &m_py);
tree->SetBranchAddress("m_pz", &m_pz);
tree->SetBranchAddress("m_energy", &m_energy);
tree->SetBranchAddress("m_kinetic_energy", &m_kinetic_energy);
tree->SetBranchAddress("m_mass", &m_mass);
tree->SetBranchAddress("m_vx_prod", &m_vx_prod);
tree->SetBranchAddress("m_vy_prod", &m_vy_prod);
tree->SetBranchAddress("m_vz_prod", &m_vz_prod);
tree->SetBranchAddress("m_vx_decay", &m_vx_decay);
tree->SetBranchAddress("m_vy_decay", &m_vy_decay);
tree->SetBranchAddress("m_vz_decay", &m_vz_decay);
tree->SetBranchAddress("m_pdg_in_particle", &m_pdg_in_particle);
tree->SetBranchAddress("m_pdg_out_particle", &m_pdg_out_particle);
tree->SetBranchAddress("m_trackid_in_particle", &m_trackid_in_particle);
tree->SetBranchAddress("m_trackid_out_particle", &m_trackid_out_particle);
tree->SetBranchAddress("m_status", &m_status);

 Int_t last_event_id_MC = -1;
 Int_t current_event_first_entry = -1;
 bool found_tau_lepton = false;
 bool got_primvtx = false;
 int tau_lepton_track_id = 0;
 Long64_t event_count = 0;
 Long64_t event_max = 0;
 
 create_histos(outputFile);
 
 std::cout << "Number of entries " << tree->GetEntries() << std::endl;

 // dump event
 bool dump = true;

// Now you can loop over the entries in the tree to read them
 for (Long64_t ientry = 0; ientry < tree->GetEntries() && (event_count < event_max || event_max == 0); ientry++) {
    tree->GetEntry(ientry);

     if (m_event_id_MC != last_event_id_MC) {

       if(event_count % 1000 == 0) {
	 std::cout << "Event #" << event_count << std::endl;
       }

       if(m_event_id_MC < last_event_id_MC) {
	 std::cout << "Event counter has gone back to low value... duplicate events??" << std::endl;
	 TFile *currentFile = tree->GetFile();
	 std::cout << "Entry  is from file: " << currentFile->GetName() << std::endl;
	 //	 break;
       }

       current_event_first_entry = ientry;
       
       event_count++;
       last_event_id_MC = m_event_id_MC;

       // process previous event
       if(event_count>1) {
	 //	 smear_event();
	 kinematics_event();
	 if(event.istau && event.isCC & event.n_taudecay==0) {
	   std::cout << "Could not find tau decay product??" << std::endl;
	   //	 exit(1);
	   dump_event();
	   outFile << "Tau error " << event.run_number << " " << event.event_id << std::endl;
	 } else {
	   fill_histos();
	 }

	 dump = event.istau || event_count < 10;
	 
	 if(dump) {
	   dump_event();
	 };
	 
       }
	 
       clear_event();
       event.run_number = m_runnumber;
       event.event_id = m_event_id_MC;
       found_tau_lepton = false;
       got_primvtx = false;
     }

     if(!got_primvtx && m_status == 1) {
       event.prim_vx[0] = m_vx_prod;
       event.prim_vx[1] = m_vy_prod;
       event.prim_vx[2] = m_vz_prod;
       got_primvtx = true;
     }

     struct PO aPO;
     aPO.m_pdg_id = m_pdg_id;
     aPO.m_track_id = m_track_id;
     aPO.m_px = m_px/1e3;
     aPO.m_py = m_py/1e3;
     aPO.m_pz = m_pz/1e3;
     aPO.m_vx_decay = m_vx_prod-event.prim_vx[0];
     aPO.m_vy_decay = m_vy_prod-event.prim_vx[1];
     aPO.m_vz_decay = m_vz_prod-event.prim_vx[2];
     aPO.nparent = m_trackid_in_particle->size();
     for (int i=0; i<aPO.nparent;i++){
       aPO.m_trackid_in_particle[i] = m_trackid_in_particle->at(i);
     };
     aPO.m_status = m_status;

     if(m_track_id < 20000 && m_status != 3) {
       event.POs[event.n_particles++] = aPO;
     }

     // found charged tau lepton - store decay products
     if(!found_tau_lepton && abs(aPO.m_pdg_id) == 15) {
       found_tau_lepton = true;
       tau_lepton_track_id = m_track_id;
     }

     if(found_tau_lepton) {
       for(int i=0;i<aPO.nparent;i++) {
	 if(aPO.m_trackid_in_particle[i] == tau_lepton_track_id) {
	   event.taudecay[event.n_taudecay++] = aPO;
	   double decaylength = sqrt(aPO.m_vx_decay*aPO.m_vx_decay+aPO.m_vy_decay*aPO.m_vy_decay+aPO.m_vz_decay*aPO.m_vz_decay);
	   event.tautracklength = decaylength;
	 }
       }
     }
     
 }
 
 
 close_histos();

 // Close the file
 outFile.close();

 std::cout << "We're done." << std::endl;
 
}

int main() {

  // Set up random number generator
  std::random_device rd;  // Seed generator
  std::mt19937 gen(rd()); // Mersenne Twister engine
  
  t();
}

