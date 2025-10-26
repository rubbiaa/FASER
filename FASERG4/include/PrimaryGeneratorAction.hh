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
#include <fstream>
#include <mutex>

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
  
  // added by Umut
  void SetWantMuonBackground(G4bool flag) { fWantMuonBackground = flag; }
  void SetWantSingleParticle(G4bool flag) { fWantSingleParticle = flag; }
 
  
  void SetROOTInputFileName(G4String name); ///< Set the name of the ROOT input file for the primary particles
  G4String fROOTInputFileName = "def";    ///< Name of the ROOT input file for the primary particles - default value
  
  void SetFileNumber(G4int number); ///< Set the number of the file to be read, used for larger data sets
  G4int fFileNumber = 0;	    ///< Number of the file to be read, used for larger data sets - default value
  
  void SetfNStartEvent(G4int value) { fNStartEvent = value; }; ///< Set the number of events per file to be read, used for larger data sets
  G4int fNStartEvent = 0;	    ///< Number of events per file to be read, used for larger data sets - default value

  void Clear(); ///< Clear the data from the previous event

  const TPOEvent* GetTPOEvent() const { return &fTPOEvent; };
  
  // added by Umut: set single-particle momentum (GeV) from messenger
  void SetSingleParticleMomentum(double gev);

private:
  ParticleManager *fParticleManager = nullptr; ///< Particle manager, which is used to generate the primary particles
  
  PrimaryGeneratorMessenger *fMessenger = nullptr;  ///< Messenger class for the primary generator action

  std::vector<G4ParticleGun*> fParticleGuns; // Multiple particle guns

  int valid_event;
  int n_passed_event;
  
  TFile *m_ROOTInputFile = nullptr;
  TTree *m_POEventTree = nullptr;
  size_t tree_ientry;
  TPOEvent fTPOEvent;

  // added by Umut: file to dump generated muons (CSV)
  std::ofstream m_muonDumpFile;
  std::mutex m_muonDumpMutex;
  // runtime control flags (can be set via PrimaryGeneratorMessenger)
  bool fWantMuonBackground = false;
  bool fWantSingleParticle = false;
  // single particle mode momentum (in GeV)
  double fSingleParticleMomentum = 100.0; // default 100 GeV

};

#endif
