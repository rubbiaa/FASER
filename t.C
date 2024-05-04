// FASER kinematical analysis
// A. Rubbia May 2024
//

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <map>

#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>

// #include "Pythia8/Pythia.h"
#include <TDatabasePDG.h>

#define MAXPARENT 10
struct PO {
  int m_pdg_id;
  int m_track_id;
  double m_px;
  double m_py;
  double m_pz;
  int nparent;
  int m_trackid_in_particle[MAXPARENT];
  int m_status;
};

#define MAXPARTICLES 1000
struct EVENT {
  int run_number;
  int event_id;
  bool isCC;
  bool istau;
  int tau_decaymode; // =1 e, =2 mu, =3 1-prong, =4 rho =5 3-prong, =6 other
  size_t n_particles;
  struct PO in_neutrino;
  struct PO out_lepton;
  struct PO POs[MAXPARTICLES];
  size_t n_taudecay;
  struct PO taudecay[MAXPARTICLES];
  double spx, spy, spz;
  double vis_spx, vis_spy, vis_spz;
  double jetpx, jetpy, jetpz;
  double tauvis_px, tauvis_py, tauvis_pz;
  double Evis, ptmiss;
} event;

void clear_event() {
  event.run_number = event.event_id = -1;
  event.n_particles = 0;
  event.n_taudecay = 0;
  event.tau_decaymode = -1;
  event.isCC = false;
  event.istau = false;
  event.spx=event.spy=event.spz=0;
  event.tauvis_px=event.tauvis_py=event.tauvis_pz=0;
};

bool is_lepton(int pdgid) {
  int pdgidabs = abs(pdgid);
  return (pdgidabs >= 11 && pdgidabs <= 16);
}

bool is_neutrino(int pdgid) {
  int pdgidabs = abs(pdgid);
  return (pdgidabs == 12 || pdgidabs == 14 || pdgidabs == 16);
}

void kinematics_event() {
  bool got_out_lepton = false;
  for (size_t i=0; i<event.n_particles; i++) {
    struct PO aPO = event.POs[i];
    if(aPO.m_status == 4 && i==0) {
      event.in_neutrino = aPO;
      event.istau = (abs(aPO.m_pdg_id) == 16);
    }
    if(!got_out_lepton && aPO.m_status == 1 && is_lepton(aPO.m_pdg_id) ) {
      event.out_lepton = aPO;
      got_out_lepton = true;
    }
    if(aPO.m_status == 1) {
      event.spx += aPO.m_px;
      event.spy += aPO.m_py;
      event.spz += aPO.m_pz;
    }
  }
  event.jetpx = event.spx-event.out_lepton.m_px;
  event.jetpy = event.spy-event.out_lepton.m_py;
  event.jetpz = event.spz-event.out_lepton.m_pz;
  event.vis_spx = event.spx;
  event.vis_spy = event.spy;
  event.vis_spz = event.spz;
  event.isCC = !(event.in_neutrino.m_pdg_id == event.out_lepton.m_pdg_id);
  if(event.istau && event.isCC){
    TDatabasePDG *pdgDB = TDatabasePDG::Instance();
    int nc = 0, nn = 0;
    for(int i=0; i<event.n_taudecay; i++) {
      struct PO aPO = event.taudecay[i];
      TParticlePDG *particle = pdgDB->GetParticle(aPO.m_pdg_id);
      if(aPO.m_status == 1 && !is_neutrino(aPO.m_pdg_id)){
	if(abs(aPO.m_pdg_id) == 11) {
	  event.tau_decaymode = 1;
	}
	if(abs(aPO.m_pdg_id) == 13) {
	  event.tau_decaymode = 2;
	}
	if(particle->Charge() == 0){
	  nn++;
	} else {
	  nc++;
	}
	event.tauvis_px += aPO.m_px;
	event.tauvis_py += aPO.m_py;
	event.tauvis_pz += aPO.m_pz;
      }
    }
    if(event.tau_decaymode < 0){
      event.tau_decaymode = 6;
      if(nc==1 && nn == 0) event.tau_decaymode = 3;
      if(nc==1 && nn == 1) event.tau_decaymode = 4;
      if(nc==3) event.tau_decaymode = 5;
    }
    event.vis_spx = event.tauvis_px + event.jetpx;
    event.vis_spy = event.tauvis_py + event.jetpy;
    event.vis_spz = event.tauvis_pz + event.jetpz;
  }
  event.Evis = sqrt(event.vis_spx*event.vis_spx + event.vis_spy*event.vis_spy + event.vis_spz*event.vis_spz);
  event.ptmiss = sqrt(event.vis_spx*event.vis_spx + event.vis_spy*event.vis_spy);
}


void dump_PO(struct PO aPO,  TDatabasePDG *pdgDB) {
  TParticlePDG *particle = pdgDB->GetParticle(aPO.m_pdg_id);
  std::cout << std::setw(10) << aPO.m_track_id << " " << aPO.m_pdg_id << " ";
  if(particle) {
    std::cout << std::setw(10) << particle->GetName() << " ";
  } else {
    std::cout << std::setw(10) << " ?? ";
  }
  std::cout << std::setw(10) << aPO.m_px << " " << " " << aPO.m_py << " " << aPO.m_pz << " " << aPO.m_status << " ";
  for (size_t j=0; j<aPO.nparent; j++) {
    std::cout << std::setw(10) << aPO.m_trackid_in_particle[j] << " ";
  }
  std::cout << std::endl;
}

void dump_event() {
  double spx=0, spy=0, spz=0;
  TDatabasePDG *pdgDB = TDatabasePDG::Instance();
  if(event.isCC) {
    std::cout << "--- Run " << event.run_number << " Event " << event.event_id << " ---------------------------------------- CC --------------------------------------------" << std::endl;
  } else {
    std::cout << "--- Run " << event.run_number << " Event " << event.event_id << " ------------------------------------------- NC ------------------------------------------" << std::endl;
  }
  for (size_t i=0; i<event.n_particles; i++) {
    struct PO aPO = event.POs[i];
    dump_PO(aPO, pdgDB);
  }
  std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
  std::cout << std::setw(10) << "Outgoing lepton:           ";
  dump_PO(event.out_lepton, pdgDB);
  std::cout << std::setw(10) << "Jet :                      " << event.jetpx << " " << event.jetpy << " " << event.jetpz << std::endl;
  std::cout << std::setw(10) << "Sum final state particles: " << event.spx << " " << event.spy << " " << event.spz << std::endl;
  std::cout << std::setw(10) << "Sum final state particles (VIS): " << event.vis_spx << " " << event.vis_spy << " " << event.vis_spz << std::endl;
  std::cout << std::setw(10) << "Ptmiss = " << event.ptmiss << "  Evis = " << event.Evis << std::endl;
  std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
  if(event.n_taudecay>0) {
    std::cout << "Tau decay mode : " << event.tau_decaymode << std::endl;
    for (size_t i=0; i<event.n_taudecay; i++) {
      struct PO aPO = event.taudecay[i];
      dump_PO(aPO, pdgDB);
    }
    std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
  }
}

TFile *tuple_file;
TTree *event_tree;


void create_histos(std::string outputFile) {
  // Open a ROOT file for writing
  tuple_file = new TFile(outputFile.c_str(), "RECREATE");
  event_tree = new TTree("event_tree", "Event data");

  event_tree->Branch("run_number", &event.run_number);
  event_tree->Branch("event_id", &event.event_id);
  event_tree->Branch("isCC", &event.isCC);
  event_tree->Branch("istau", &event.istau);
  event_tree->Branch("tau_decaymode",&event.tau_decaymode);
  event_tree->Branch("n_particles", &event.n_particles);

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
  event_tree->Branch("Evis", &event.Evis);
  event_tree->Branch("ptmiss", &event.ptmiss);
  
};

void fill_histos() {
  // Fill the tree
  event_tree->Fill();
}

void close_histos() {
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
 int tau_lepton_track_id = 0;
 
 Int_t event_count = 0;
 Int_t event_max = 10000000;
 
 create_histos(outputFile);
 
 std::cout << "Number of entries " << tree->GetEntries() << std::endl;

 // dump event
 bool dump = true;

 std::ofstream outFile("error.txt");

// Now you can loop over the entries in the tree to read them
for (Long64_t ientry = 0; ientry < tree->GetEntries() && event_count < event_max; ientry++) {
    tree->GetEntry(ientry);

     if (m_event_id_MC != last_event_id_MC) {

       if(event_count % 1000 == 0) {
	 std::cout << "Event #" << event_count << std::endl;
       }

       if(m_event_id_MC < last_event_id_MC) {
	 std::cout << "Event counter has gone back to low value... duplicate events??" << std::endl;
	 TFile *currentFile = tree->GetFile();
	 std::cout << "Entry  is from file: " << currentFile->GetName() << std::endl;
	 break;
       }

       current_event_first_entry = ientry;
       
       event_count++;
       last_event_id_MC = m_event_id_MC;

       // process previous event
       if(event_count>1) {
	 kinematics_event();
	 fill_histos();

	 dump = event.istau;
	 
	 if(dump) {
	   dump_event();
	 };
	 
	 if(event.istau && event.n_taudecay==0) {
	   std::cout << "Could not find tau decay product??" << std::endl;
	   //	 exit(1);
	   outFile << event.run_number << " " << event.event_id << std::endl;
	 }
       }
	 
       clear_event();
       event.run_number = m_runnumber;
       event.event_id = m_event_id_MC;
       found_tau_lepton = false;	
     }

     struct PO aPO;
     aPO.m_pdg_id = m_pdg_id;
     aPO.m_track_id = m_track_id;
     aPO.m_px = m_px/1e3;
     aPO.m_py = m_py/1e3;
     aPO.m_pz = m_pz/1e3;
     aPO.nparent = m_trackid_in_particle->size();
     for (int i=0; i<aPO.nparent;i++){
       aPO.m_trackid_in_particle[i] = m_trackid_in_particle->at(i);
     };
     aPO.m_status = m_status;

     if(m_track_id < 20000 && m_status != 3) {
       event.POs[event.n_particles++] = aPO;
     }

     if(event.event_id == 724) {
       std::cout << " shit " << std::endl;
       TDatabasePDG *pdgDB = TDatabasePDG::Instance();
       dump_PO(aPO,pdgDB);
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
  t();
}

