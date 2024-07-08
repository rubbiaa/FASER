#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "PrimaryGeneratorMessenger.hh"
#include "ParticleManager.hh"


#include "typedef.h"

#include "TPOEvent.hh"

class G4ParticleGun;
class G4Event;
class PrimaryGeneratorMessenger;

/**
 * @class PrimaryGeneratorAction
 * @brief Primary generator action class, which is called at the beginning of each event, and generates the primary particles.
 *
 * This class is responsible for generating the primary particles.
 * It is feed with the information from the PrimaryGeneratorMessenger class.
 */
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  
  PrimaryGeneratorAction(ParticleManager *t_particleManager); ///< Constructor
  ~PrimaryGeneratorAction() override; ///< Destructor
  
  void GeneratePrimaries(G4Event *) override; ///< Generate the primary particles, called at the beginning of each event by Geant4 kernel
  
  
  void SetRandomFlag(G4bool flag);
  
  void SetROOTInputFileName(G4String name); ///< Set the name of the ROOT input file for the primary particles
  G4String fROOTInputFileName = "def";    ///< Name of the ROOT input file for the primary particles - default value
  
  void SetFileNumber(G4int number); ///< Set the number of the file to be read, used for larger data sets
  G4int fFileNumber = 0;	    ///< Number of the file to be read, used for larger data sets - default value
  
  void SetNEventsPerFile(G4int number); ///< Set the number of events per file to be read, used for larger data sets
  G4int fNEventsPerFile = 0;	    ///< Number of events per file to be read, used for larger data sets - default value

  void Clear(); ///< Clear the data from the previous event

  const TPOEvent* GetTPOEvent() const { return &fTPOEvent; };
  
private:
  ParticleManager *fParticleManager = nullptr; ///< Particle manager, which is used to generate the primary particles
  
  PrimaryGeneratorMessenger *fMessenger = nullptr;  ///< Messenger class for the primary generator action

  std::vector<G4ParticleGun*> fParticleGuns; // Multiple particle guns

  int NStartEvent;
  int valid_event;
  int n_passed_event;
  
  TChain *tree; // input FASER ntuple chain
  size_t tree_ientry;
  
  Int_t last_event_id_MC = -1;
  Int_t current_event_first_entry = -1;
  bool found_tau_lepton = false;
  bool got_primvtx = false;
  int tau_lepton_track_id = 0;

  // Set up variables to hold the data and link them to the tree branches
  Int_t m_runnumber, m_event_id_MC, m_track_id, m_pdg_id, m_num_in_particle, m_num_out_particle;
  Double_t m_px, m_py, m_pz, m_energy, m_kinetic_energy, m_mass;
  Float_t m_vx_prod, m_vy_prod, m_vz_prod, m_vx_decay, m_vy_decay, m_vz_decay;
  std::vector<int> *m_pdg_in_particle = nullptr, *m_pdg_out_particle = nullptr;
  std::vector<int> *m_trackid_in_particle = nullptr, *m_trackid_out_particle = nullptr;
  Int_t m_status;

  TPOEvent fTPOEvent;

};

#endif
