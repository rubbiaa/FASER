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

	const TPOEvent *branch_POEvent = GetTPOEvent();

	if (m_ROOTInputFile == nullptr) {
		// Open FASERMC PO input ntuple files
		std::string inputFile = fROOTInputFileName;
		m_ROOTInputFile = new TFile(inputFile.c_str());

		m_ROOTInputFile->GetObject("POEvent", m_POEventTree);
		if (!m_POEventTree)
		{
			std::cerr << "Error: Cannot find the TTree named 'myTree' in the file." << std::endl;
			m_ROOTInputFile->Close();
			exit(1);
		}

		m_POEventTree->SetBranchAddress("event", &branch_POEvent);

		tree_ientry = 0;

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
		// Generate an event in the tungsten
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

	fTPOEvent.dump_event();	

	std::vector<XYZVector> ParticlePosition;
	std::vector<XYZVector> ParticleMomentum;
	std::vector<double> ParticleTime;
	std::vector<int> ParticlePDGCode;
	std::vector<double> ParticleEnergy;
	std::vector<double> ParticleKinEnergy;
	std::vector<int> ParticleTrackID;
	std::vector<int> DecayModeFlag;

	ParticlePosition.clear();
	ParticleMomentum.clear();
	ParticleTime.clear();
	ParticlePDGCode.clear();
	ParticleEnergy.clear();
	ParticleKinEnergy.clear();
	ParticleTrackID.clear();
	DecayModeFlag.clear();

	valid_event++;
	for (G4int i = 0; i < fTPOEvent.n_particles(); ++i)
	{
		struct PO aPO = fTPOEvent.POs[i];
		G4ParticleDefinition *particle = particleTable->FindParticle(aPO.m_pdg_id);

		if (particle != nullptr && aPO.m_status == 1)
		{
//			if(aPO.m_pdg_id != 311) continue;  // TODO/FIXME remove
			G4ParticleGun *particleGun = new G4ParticleGun(1);

			ParticlePDGCode.push_back(aPO.m_pdg_id);

			// 	G4cout << "Particle found: " << particle->GetParticleName()
			//  	<< ", mass: " << particle->GetPDGMass()/GeV << " GeV"
			//  	<< ", charge: " << particle->GetPDGCharge() << G4endl;

			particleGun->SetParticleDefinition(particle);

			particleGun->SetParticlePosition(G4ThreeVector(x * mm, y * mm, z * mm));
			G4ThreeVector StartMomentum(aPO.m_px * GeV, aPO.m_py * GeV, aPO.m_pz * GeV);

#if 0

			ParticlePosition.push_back(XYZVector(x * mm, y * mm, z * mm));

			ParticleMomentum.push_back(XYZVector(aPO.m_px * GeV, aPO.m_py * GeV, aPO.m_pz * GeV));

			ParticleTime.push_back(0);
			ParticleEnergy.push_back(aPO.m_energy);
			ParticleKinEnergy.push_back(aPO.m_kinetic_energy);
			ParticleTrackID.push_back(aPO.m_track_id);
			DecayModeFlag.push_back(0);
#endif

			particleGun->SetParticleMomentum(StartMomentum);
			fParticleGuns.push_back(particleGun);
		}
	}

	// still need to add tau decay products if any present
	if (fTPOEvent.n_taudecay() > 0)
	{
		G4cout << "Found a tau lepton..." << G4endl;
		fTPOEvent.dump_event();
		std::cout << "Tau decay mode : " << fTPOEvent.tau_decaymode << std::endl;

		for (size_t i = 0; i < fTPOEvent.n_taudecay(); i++)
		{
			struct PO aPO = fTPOEvent.taudecay[i];
			if (aPO.m_status == 1)
			{
				G4ParticleGun *particleGun = new G4ParticleGun(1);

				G4ParticleDefinition *particle = particleTable->FindParticle(aPO.m_pdg_id);
				ParticlePDGCode.push_back(aPO.m_pdg_id);
				G4cout << "Particle found: " << particle->GetParticleName()
					   << ", mass: " << particle->GetPDGMass() / GeV << " GeV"
					   << ", charge: " << particle->GetPDGCharge() << G4endl;

				particleGun->SetParticleDefinition(particle);

				particleGun->SetParticlePosition(G4ThreeVector((x + aPO.m_vx_decay) * mm, (y + aPO.m_vy_decay) * mm, (z + aPO.m_vz_decay) * mm));
				ParticlePosition.push_back(XYZVector((x + aPO.m_vx_decay) * mm, (y + aPO.m_vy_decay) * mm, (z + aPO.m_vz_decay) * mm));

				G4ThreeVector StartMomentum(aPO.m_px * GeV, aPO.m_py * GeV, aPO.m_pz * GeV);
				ParticleMomentum.push_back(XYZVector(aPO.m_px * GeV, aPO.m_py * GeV, aPO.m_pz * GeV));

				ParticleTime.push_back(0);
				ParticleEnergy.push_back(aPO.m_energy);
				ParticleKinEnergy.push_back(aPO.m_kinetic_energy);
				ParticleTrackID.push_back(aPO.m_track_id);
				DecayModeFlag.push_back(1);

				particleGun->SetParticleMomentum(StartMomentum);
				fParticleGuns.push_back(particleGun);
			}
		}
	}

	int NParticlesIF = ParticlePDGCode.size();

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

void PrimaryGeneratorAction::SetNEventsPerFile(G4int value) { fNEventsPerFile = value; }
