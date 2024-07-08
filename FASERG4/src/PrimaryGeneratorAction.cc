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
  
	// Open FASER input ntuple files
	//	std::string inputFile = "/home/rubbiaa/FASER/FASERDATA/sim/*.root";
	//std::string inputFile = "/home/rubbiaa/FASER/FASERDATA/sim/FaserMC-MC22_Genie_all_6invab-200006-00001-s0010-NTUP.root";
	
	tree_ientry = 0;

	valid_event = 0;

	n_passed_event = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fMessenger;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	int G4EventIt = anEvent->GetEventID();
	if (G4EventIt == 0){ 
		NStartEvent = fFileNumber * fNEventsPerFile;
	}

	std::string inputFile = fROOTInputFileName;
	tree = new TChain("m_NuMCTruth_tree");  

	fParticleGuns.clear();

	tree->Add(inputFile.c_str());

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

	if (tree_ientry==0) std::cout << "Number of entries " << tree->GetEntries() << std::endl;

	fTPOEvent.clear_event();
	found_tau_lepton = false;
	got_primvtx = false;

	tree->GetEntry(tree_ientry);
	if(last_event_id_MC != m_event_id_MC) {
	  last_event_id_MC = m_event_id_MC;
	}
	// find the start position of the event
	while (G4EventIt==0 && n_passed_event < NStartEvent && tree_ientry < tree->GetEntries()) {
		//std::cout<<"Event ID: "<<m_event_id_MC<<"passed"<<std::endl;
		tree_ientry++;
		tree->GetEntry(tree_ientry);
		if(last_event_id_MC != m_event_id_MC) {
	  		last_event_id_MC = m_event_id_MC;
			n_passed_event++;
		}
	} 


	while (m_event_id_MC == last_event_id_MC) {
	  if(tree_ientry >= tree->GetEntries()) {
	    std::cout << "Not enough events in input ntuple..." << std::endl;
	    exit(1);
	  }
	  
	  tree->GetEntry(tree_ientry);
	  if(m_event_id_MC == last_event_id_MC) {
	    tree_ientry++;
	    fTPOEvent.run_number = m_runnumber;
	    fTPOEvent.event_id = m_event_id_MC;

	    if(!got_primvtx && m_status == 1) {
	      fTPOEvent.prim_vx[0] = m_vx_prod;
	      fTPOEvent.prim_vx[1] = m_vy_prod;
	      fTPOEvent.prim_vx[2] = m_vz_prod;
	      got_primvtx = true;
	    }
	    
	    struct PO aPO;
	    aPO.m_pdg_id = m_pdg_id;
	    aPO.m_track_id = m_track_id;
	    aPO.m_px = m_px/1e3;
	    aPO.m_py = m_py/1e3;
	    aPO.m_pz = m_pz/1e3;
		aPO.m_energy = m_energy/1e3;
		aPO.m_kinetic_energy = m_kinetic_energy/1e3;	
	    aPO.m_vx_decay = m_vx_prod-fTPOEvent.prim_vx[0];
	    aPO.m_vy_decay = m_vy_prod-fTPOEvent.prim_vx[1];
	    aPO.m_vz_decay = m_vz_prod-fTPOEvent.prim_vx[2];
	    aPO.nparent = m_trackid_in_particle->size();
	    for (int i=0; i<aPO.nparent;i++){
	      aPO.m_trackid_in_particle[i] = m_trackid_in_particle->at(i);
	    };
	    aPO.m_status = m_status;
		aPO.geanttrackID = -1;
	    
	    if(m_track_id < 20000 && m_status != 3) {
	      fTPOEvent.POs[fTPOEvent.n_particles++] = aPO;
	    }
	    
	    // found charged tau lepton - store decay products
	    if(!found_tau_lepton && abs(aPO.m_pdg_id) == 15) {
	      found_tau_lepton = true;
	      tau_lepton_track_id = m_track_id;
	    }
	    
	    if(found_tau_lepton) {
	      for(int i=0;i<aPO.nparent;i++) {
				if(aPO.m_trackid_in_particle[i] == tau_lepton_track_id) {
				  fTPOEvent.taudecay[fTPOEvent.n_taudecay++] = aPO;
				  double decaylength = sqrt(aPO.m_vx_decay*aPO.m_vx_decay+aPO.m_vy_decay*aPO.m_vy_decay+aPO.m_vz_decay*aPO.m_vz_decay);
				  fTPOEvent.tautracklength = decaylength;
				}
	      }
	    }
	    
	  }	    
	}
	
	if(fTPOEvent.istau && fTPOEvent.n_taudecay==0) {
	  std::cout << "Could not find tau decay product??" << std::endl;
	  //	 exit(1);
	}

	fTPOEvent.kinematics_event();
	fTPOEvent.dump_header();
	fTPOEvent.dump_event();	

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

	const DetectorConstruction* detector = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

	// Generate an event in the tungsten
	G4double theta = G4UniformRand() * 2 * M_PI;
	G4double x = 50 * cos(theta);
	G4double y = 50 * sin(theta);
	G4double zlayer = detector->getScintillatorSizeZ()+detector->gettargetWSizeZ();
	G4double zfront = -detector->getNumberReplicas()*zlayer/2.0;
	G4int wanted_layer = 0; // select first layer
	G4double z = zfront + wanted_layer*zlayer + detector->getScintillatorSizeZ() + G4UniformRand()*detector->gettargetWSizeZ();

	XYZVector vtxpos;
	vtxpos.SetX(x);	vtxpos.SetY(y);vtxpos.SetZ(z);
	fParticleManager->setVertexInformation(vtxpos);

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

	// only focus on nueCC events
	bool isnueCC = fTPOEvent.isCC && abs(fTPOEvent.out_lepton.m_pdg_id) == 11;
	bool event_wanted = true;
	if (event_wanted){
	    valid_event++;
		for (G4int i = 0; i < fTPOEvent.n_particles; ++i)
		{
		  	struct PO aPO = fTPOEvent.POs[i];
			G4ParticleDefinition* particle = particleTable->FindParticle(aPO.m_pdg_id);

		  	if(particle != nullptr && aPO.m_status == 1) {
		  	  	G4ParticleGun* particleGun = new G4ParticleGun(1);

			  	ParticlePDGCode.push_back(aPO.m_pdg_id);

		  	 // 	G4cout << "Particle found: " << particle->GetParticleName()
		  	 //  	<< ", mass: " << particle->GetPDGMass()/GeV << " GeV"
		  	 //  	<< ", charge: " << particle->GetPDGCharge() << G4endl;

		  	  	particleGun->SetParticleDefinition(particle);

		  	  	particleGun->SetParticlePosition(G4ThreeVector(x * mm, y * mm, z * mm));
		  		ParticlePosition.push_back(XYZVector(x * mm, y * mm, z * mm));

		  	  	G4ThreeVector StartMomentum(aPO.m_px * GeV, aPO.m_py * GeV, aPO.m_pz * GeV);
		  		ParticleMomentum.push_back(XYZVector(aPO.m_px * GeV, aPO.m_py * GeV, aPO.m_pz * GeV));

		  		ParticleTime.push_back(0);
		  		ParticleEnergy.push_back(aPO.m_energy);
		  		ParticleKinEnergy.push_back(aPO.m_kinetic_energy);
		  		ParticleTrackID.push_back(aPO.m_track_id);
		  		DecayModeFlag.push_back(0);

		  		particleGun->SetParticleMomentum(StartMomentum);
		  	  	fParticleGuns.push_back(particleGun);
		  	}
		}

		// still need to add tau decay products if any present
	  	if(fTPOEvent.n_taudecay>0) {
	    	G4cout << "Found a tau lepton..." << G4endl;
	    	fTPOEvent.dump_event();
	    	std::cout << "Tau decay mode : " << fTPOEvent.tau_decaymode << std::endl;

	    	for (size_t i=0; i<fTPOEvent.n_taudecay; i++) {
	    	  	struct PO aPO = fTPOEvent.taudecay[i];
	    	  	if(aPO.m_status == 1) {
					G4ParticleGun* particleGun = new G4ParticleGun(1);

					G4ParticleDefinition* particle = particleTable->FindParticle(aPO.m_pdg_id);
					ParticlePDGCode.push_back(aPO.m_pdg_id);
					G4cout << "Particle found: " << particle->GetParticleName()
					       << ", mass: " << particle->GetPDGMass()/GeV << " GeV"
					       << ", charge: " << particle->GetPDGCharge() << G4endl;

					particleGun->SetParticleDefinition(particle);

					particleGun->SetParticlePosition(G4ThreeVector((x+aPO.m_vx_decay) * mm, (y+aPO.m_vy_decay) * mm, (z+aPO.m_vz_decay) * mm));
					ParticlePosition.push_back(XYZVector((x+aPO.m_vx_decay) * mm, (y+aPO.m_vy_decay) * mm, (z+aPO.m_vz_decay) * mm));
					
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
	}

	int NParticlesIF = ParticlePDGCode.size();
	#if 0
	fParticleManager->setVertexInformation(0, NParticlesIF, ParticlePosition, ParticleMomentum, ParticleTime, 
											ParticleEnergy, ParticleKinEnergy, ParticleTrackID, DecayModeFlag,
											ParticlePDGCode, 0);
	#endif

	std::cout<<"Number of particles: "<<NParticlesIF<<std::endl;
	
	for (auto gun : fParticleGuns) {
	  gun->GeneratePrimaryVertex(anEvent);
	}
	// Check the number of primary particles
	std::cout<<"Number of primary particles: "<<anEvent->GetNumberOfPrimaryVertex()<<std::endl;

	delete tree;

	for (auto gun : fParticleGuns) {
    	delete gun;
    }
	fParticleGuns.clear();

	std::cout<<"Valid events: "<<valid_event<<std::endl;
	
#if 0
	std::vector<XYZVector> ParticlePosition;
	std::vector<XYZVector> ParticleMomentum;
	std::vector<double> ParticleTime;
	std::vector<int> ParticlePDGCode;

	for (int partIt = 0; partIt < NParticlesIF; partIt++) {
		auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(ParticlePDGCodeIF->at(partIt));
		fParticleGun->SetParticleDefinition(particleDefinition);

		ParticlePDGCode.push_back(ParticlePDGCodeIF->at(partIt));

		G4ThreeVector StartPosition(ParticlePositionXIF->at(partIt) * mm, ParticlePositionYIF->at(partIt) * mm,
					    ParticlePositionZIF->at(partIt) * mm);
		fParticleGun->SetParticlePosition(StartPosition);

		// Used for storage in the ROOT file
		ParticlePosition.push_back(
		    XYZVector(ParticlePositionXIF->at(partIt) * mm, ParticlePositionYIF->at(partIt) * mm, ParticlePositionZIF->at(partIt) * mm));

		G4ThreeVector StartMomentum(ParticleMomentumXIF->at(partIt) * MeV, ParticleMomentumYIF->at(partIt) * MeV,
					    ParticleMomentumZIF->at(partIt) * MeV);

		fParticleGun->SetParticleMomentum(StartMomentum);

		// Used for storage in the ROOT file
		ParticleMomentum.push_back(
		    XYZVector(ParticleMomentumXIF->at(partIt) * MeV, ParticleMomentumYIF->at(partIt) * MeV, ParticleMomentumZIF->at(partIt) * MeV));

		ParticleTime.push_back(0);  /// Currently this is not used for vertex generation, hence it is set to 0, as all the particles are
					    /// generated by default at t = 0.

		fParticleGun->GeneratePrimaryVertex(anEvent);
	}
	fParticleManager->setVertexInformation(ModeIF, NParticlesIF, ParticlePosition, ParticleMomentum, ParticleTime, ParticlePDGCode, EnuIF);
#endif
	
}


void PrimaryGeneratorAction::SetROOTInputFileName(G4String value) { fROOTInputFileName = value; }

void PrimaryGeneratorAction::SetFileNumber(G4int value) { fFileNumber = value; }

void PrimaryGeneratorAction::SetNEventsPerFile(G4int value) { fNEventsPerFile = value; }
