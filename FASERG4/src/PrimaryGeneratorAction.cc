#include "PrimaryGeneratorAction.hh"

#include "G4Box.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "DetectorConstruction.hh"

#include "TPOEvent.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(ParticleManager* f_particleManager) : G4VUserPrimaryGeneratorAction()
{

  fMessenger = new PrimaryGeneratorMessenger(this);
  
  fParticleManager = f_particleManager;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	m_ROOTInputFile->Close();
	delete fMessenger;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

	// DEBUG : only primary lepton if CC otherwise random pion
	bool want_particleGun = false; //  true;

	const TPOEvent *branch_POEvent = GetTPOEvent();

	if (m_ROOTInputFile == nullptr) {
		// Open FASERMC PO input ntuple files
		std::string inputFile = fROOTInputFileName;
		m_ROOTInputFile = new TFile(inputFile.c_str(), "READ");

		m_ROOTInputFile->GetObject("POEvent", m_POEventTree);
		if (!m_POEventTree)
		{
			std::cerr << "Error: Cannot find the TTree named 'myTree' in the file." << std::endl;
			m_ROOTInputFile->Close();
			exit(1);
		}

		m_POEventTree->SetBranchAddress("event", &branch_POEvent);

		tree_ientry = fNStartEvent;

		valid_event = 0;

		n_passed_event = 0;
	}

	bool found_tau_lepton = false;
  	bool got_primvtx = false;
  	int tau_lepton_track_id = 0;

	fParticleGuns.clear();

	if(tree_ientry >= m_POEventTree->GetEntries()) {
	    G4cout << "Not enough events in input ntuple..." << G4endl;
		G4RunManager::GetRunManager()->AbortRun();
		return;
	}
	  
	m_POEventTree -> GetEntry(tree_ientry++);

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

	const DetectorConstruction* detector = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

	// Generate primary vertex position
	G4double theta = G4UniformRand() * 2 * M_PI;
	G4double x = 50 * cos(theta);
	G4double y = 50 * sin(theta);
	G4double z = 0;
	// uniformly distributed in the first 10 layers
	G4int wanted_layer = floor(G4UniformRand() * 10);
	XYZVector vtxpos;
	// decide where the event is generated
	if (G4UniformRand() < detector->fTotalWMass / detector->fTotalMass) {
		// Generate an event in the target
		fTPOEvent.setVtxTarget(TPOEvent::kVtx_in_W);
		G4double zfront = -detector->getNumberReplicas() * detector->fSandwichLength / 2.0;
		z = zfront + wanted_layer * detector->fSandwichLength + detector->getScintillatorSizeZ() + 
				G4UniformRand() * detector->gettargetWSizeZ();
	} else {
		// Generate an event in the tungsten
		fTPOEvent.setVtxTarget(TPOEvent::kVtx_in_Scint);
		G4double zfront = -detector->getNumberReplicas() * detector->fSandwichLength / 2.0;
		z = zfront + wanted_layer * detector->fSandwichLength + G4UniformRand() * detector->getScintillatorSizeZ();
	}
	vtxpos.SetX(x);
	vtxpos.SetY(y);
	vtxpos.SetZ(z);
	fParticleManager->setVertexInformation(vtxpos);
	fTPOEvent.setPrimaryVtx(x,y,z);

	int ipo_maxhadron = -1;
	if(want_particleGun) {
		// shift run number to big number
		fTPOEvent.run_number += 1000000000;
		// find most energetic pion
		for (G4int i = 0; i < fTPOEvent.n_particles(); ++i)
		{
			struct PO aPO = fTPOEvent.POs[i];
			if (abs(aPO.m_pdg_id) != 211 && aPO.m_pdg_id != 111)
				continue;
			if(ipo_maxhadron == -1 || aPO.m_energy > fTPOEvent.POs[ipo_maxhadron].m_energy) {
				ipo_maxhadron = i;
			}			
		}
	}
	bool got_pion = false;

	fTPOEvent.dump_event();	

	valid_event++;
	for (G4int i = 0; i < fTPOEvent.n_particles(); ++i)
	{
		struct PO aPO = fTPOEvent.POs[i];

		// run in particle gun mode keeping only one relevant track from event
		if(want_particleGun) {
			if(got_pion) continue;
			if(fTPOEvent.isCC) {
				if(!fTPOEvent.is_lepton(aPO.m_pdg_id)) continue;
			} else {
				if(i != ipo_maxhadron) continue;				
				got_pion = true;
			}
		}

		G4ParticleDefinition *particle = particleTable->FindParticle(aPO.m_pdg_id);

		if (particle != nullptr && aPO.m_status == 1)
		{
//			if(aPO.m_pdg_id != 15) continue;  // TODO/FIXME debug to process only taus
			G4ParticleGun *particleGun = new G4ParticleGun(1);

//			ParticlePDGCode.push_back(aPO.m_pdg_id);

			particleGun->SetParticleDefinition(particle);

			particleGun->SetParticlePosition(G4ThreeVector(x * mm, y * mm, z * mm));
			G4ThreeVector StartMomentum(aPO.m_px * GeV, aPO.m_py * GeV, aPO.m_pz * GeV);

			particleGun->SetParticleMomentum(StartMomentum);
			fParticleGuns.push_back(particleGun);
		}
	}

	int NParticlesIF = fParticleGuns.size();

	std::cout<<"Number of particles: "<<NParticlesIF<<std::endl;
	
	for (auto gun : fParticleGuns) {
	  gun->GeneratePrimaryVertex(anEvent);
	}
	// Check the number of primary particles
	std::cout<<"Number of primary particles: "<<anEvent->GetNumberOfPrimaryVertex()<<std::endl;

	for (auto gun : fParticleGuns) {
    	delete gun;
    }
	fParticleGuns.clear();

	std::cout<<"Valid events: "<<valid_event<<std::endl;	
}


void PrimaryGeneratorAction::SetROOTInputFileName(G4String value) { fROOTInputFileName = value; }

void PrimaryGeneratorAction::SetFileNumber(G4int value) { fFileNumber = value; }

