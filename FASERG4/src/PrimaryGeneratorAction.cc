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
#include "TVector3.h"
#include "TRotation.h"

PrimaryGeneratorAction::PrimaryGeneratorAction(ParticleManager* f_particleManager) : G4VUserPrimaryGeneratorAction()
{

  fMessenger = new PrimaryGeneratorMessenger(this);
  
  fParticleManager = f_particleManager;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	if(m_ROOTInputFile != nullptr) m_ROOTInputFile->Close();
	delete fMessenger;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

	// DEBUG : only primary lepton if CC otherwise random pion
	bool want_particleGun = false; //  true;
	bool want_muon_background = false; // true;
	bool want_single_particle = false; // true;
	bool want_zeropt_jet = true; // true;

	const TPOEvent *branch_POEvent = GetTPOEvent();

	if (m_ROOTInputFile == nullptr && !want_single_particle) {
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

	if(m_POEventTree != nullptr && tree_ientry >= m_POEventTree->GetEntries()) {
	    G4cout << "Not enough events in input ntuple..." << G4endl;
		G4RunManager::GetRunManager()->AbortRun();
		return;
	}
	  
	if(m_POEventTree != nullptr) m_POEventTree -> GetEntry(tree_ientry++);

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

	const DetectorConstruction* detector = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

	XYZVector vtxpos;
	if(fTPOEvent.use_GENIE_vtx) {
		double x = fTPOEvent.prim_vx.x();
		double y = fTPOEvent.prim_vx.y();
		double z = fTPOEvent.prim_vx.z();
		if(want_muon_background) {
		  // put muon in front of the detector
		  // uniform in x and y from -20 to 20 cm
		  x = (G4UniformRand() - 0.5) * 400; // in mm
		  y = (G4UniformRand() - 0.5) * 400; // in mm
		  z = -2000; // in mm, in front of the detector
		} else {
		  std::cout << " Using GENIE vtx:  x=" << x << " y=" << y << " z=" << z << " ";
		}
		vtxpos.SetX(x);
		vtxpos.SetY(y);
		vtxpos.SetZ(z);
		fParticleManager->setVertexInformation(vtxpos);
		std::cout << " Vertex target " << fTPOEvent.GENIE_vtx_name << std::endl;
	} else {
		// Generate primary vertex position
		G4double theta = G4UniformRand() * 2 * M_PI;
		G4double x = 50 * cos(theta);
		G4double y = 50 * sin(theta);
		G4double z = 0;
		// uniformly distributed in the first "n" layers
		G4int maxlayer = 5;
		G4int wanted_layer = floor(G4UniformRand() * maxlayer);
		// decide where the event is generated
		if (G4UniformRand() < detector->fTotalWMass / detector->fTotalMass) {
			// Generate an event in the target
			fTPOEvent.setVtxTarget(TPOEvent::kVtx_in_W);
			G4double zfront = -detector->getNumberReplicas() * detector->fSandwichLength / 2.0;
			z = zfront + wanted_layer * detector->fSandwichLength + detector->getScintillatorSizeZ() + 
					G4UniformRand() * detector->gettargetWSizeZ();
		} else {
			// Generate an event in the scintillator
			fTPOEvent.setVtxTarget(TPOEvent::kVtx_in_Scint);
			G4double zfront = -detector->getNumberReplicas() * detector->fSandwichLength / 2.0;
			z = zfront + wanted_layer * detector->fSandwichLength + G4UniformRand() * detector->getScintillatorSizeZ();
		}
		vtxpos.SetX(x);
		vtxpos.SetY(y);
		vtxpos.SetZ(z);
		fParticleManager->setVertexInformation(vtxpos);
		fTPOEvent.setPrimaryVtx(x,y,z);
	}

	// single particle gun mode
	if(want_single_particle) {

		int popt = 3;
		// 0 = electron
		// 1 = photon
		// 2 = muon
		// 3 = pion
		int particle_options[5] = {11, 22, 13, 211}; // e, gamma, mu+, pi+
		int pdgid = particle_options[popt];
		double momentumMagnitude = 100.0; // in GeV
	
		fTPOEvent.clear_event();
		fTPOEvent.POs.clear();
		fTPOEvent.setVtxTarget(TPOEvent::kVtx_in_Scint);
		vtxpos.SetX(0);
		vtxpos.SetY(0);
		vtxpos.SetZ(-800); // in mm, in front of the detector
		fTPOEvent.setPrimaryVtx(vtxpos.x(), vtxpos.y(), vtxpos.z());
		fParticleManager->setVertexInformation(vtxpos);
		fTPOEvent.run_number = 888*100000 + popt*10000 + int(momentumMagnitude);
		fTPOEvent.event_id = valid_event;
		struct PO aPO;
		aPO.m_pdg_id = pdgid;
		G4ParticleDefinition *particle = particleTable->FindParticle(aPO.m_pdg_id);
		double mass = particle->GetPDGMass()/GeV;
		aPO.m_track_id = 1;
		aPO.m_status = 1;
		aPO.m_px = 0;
		aPO.m_py = 0;
		aPO.m_pz = momentumMagnitude;
		aPO.m_energy = sqrt(aPO.m_px*aPO.m_px + aPO.m_py*aPO.m_py + aPO.m_pz*aPO.m_pz + mass*mass); // in GeV
		aPO.m_vx_decay = 0;
		aPO.m_vy_decay = 0;
		aPO.m_vz_decay = 0;
		aPO.nparent = 0;
		aPO.geanttrackID = -1;
		fTPOEvent.POs.push_back(aPO);
	}

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

	if(want_zeropt_jet) {
		vtxpos.SetX(0);
		vtxpos.SetY(0);
		vtxpos.SetZ(-800); // in mm, in front of the detector
		fTPOEvent.setPrimaryVtx(vtxpos.x(), vtxpos.y(), vtxpos.z());
		fParticleManager->setVertexInformation(vtxpos);
		fTPOEvent.run_number = 9990000;

		TVector3 nHat(fTPOEvent.jetpx, fTPOEvent.jetpy, fTPOEvent.jetpz);
		nHat = nHat.Unit(); // Normalize it
		TVector3 zAxis(0, 0, 1);
		// Compute rotation axis (cross product)
	    TVector3 rotationAxis = nHat.Cross(zAxis);

    	double angle = nHat.Angle(zAxis); // Angle between nHat and z-axis

    	TRotation rot;

		if (rotationAxis.Mag() != 0) {
			rotationAxis = rotationAxis.Unit(); // Normalize the axis
			rot.Rotate(angle, rotationAxis);    // Rotate by angle around rotationAxis
		}

		// loop over PO and rotate momenta
		TVector3 totvec = TVector3(0,0,0);
		for (G4int i = 0; i < fTPOEvent.n_particles(); ++i)
		{
			struct PO &aPO = fTPOEvent.POs[i];
			// check if parent of the particle is the neutrino
			if(aPO.nparent > 0) {
				int parent_trackid = aPO.m_trackid_in_particle[0];
				if(parent_trackid == 0) {
					continue;
				}
			} else {
				continue;
			}
			TVector3 pvec(aPO.m_px, aPO.m_py, aPO.m_pz);
			TVector3 pvec_rotated = rot * pvec;
			totvec += pvec_rotated;
			// print rotated momentum
			std::cout << " Original p: " << pvec.X() << " " << pvec.Y() << " " << pvec.Z() << " ";
			std::cout << " Rotated p: " << pvec_rotated.X() << " " << pvec_rotated.Y() << " " << pvec_rotated.Z() << " ";
			std::cout << std::endl;
			aPO.m_px = pvec_rotated.X();
			aPO.m_py = pvec_rotated.Y();
			aPO.m_pz = pvec_rotated.Z();
		}
		std::cout << " Jet total p: " << totvec.X() << " " << totvec.Y() << " " << totvec.Z() << " ";
		std::cout << " Jet magnitude: " << totvec.Mag() << " Jet pt: " << totvec.Perp() << std::endl;
	}

	fTPOEvent.dump_event();	

	valid_event++;
	if (!want_muon_background)
	{
		for (G4int i = 0; i < fTPOEvent.n_particles(); ++i)
		{
			struct PO aPO = fTPOEvent.POs[i];

			// run in particle gun mode keeping only one relevant track from event
			if (want_particleGun)
			{
				if (got_pion)
					continue;
				if (fTPOEvent.isCC)
				{
					if (!fTPOEvent.is_lepton(aPO.m_pdg_id))
						continue;
				}
				else
				{
					if (i != ipo_maxhadron)
						continue;
					got_pion = true;
				}
			}

			G4ParticleDefinition *particle = particleTable->FindParticle(aPO.m_pdg_id);

			if (particle != nullptr && aPO.m_status == 1)
			{

				if(want_zeropt_jet) {
					// check if parent of the particle is the neutrino
					if(aPO.nparent > 0) {
						int parent_trackid = aPO.m_trackid_in_particle[0];
						if(parent_trackid == 0) {
							continue;
						}
					} 
				}

				//			if(aPO.m_pdg_id != 15) continue;  // TODO/FIXME debug to process only taus
				G4ParticleGun *particleGun = new G4ParticleGun(1);

				//			ParticlePDGCode.push_back(aPO.m_pdg_id);

				particleGun->SetParticleDefinition(particle);

				particleGun->SetParticlePosition(G4ThreeVector(vtxpos.x() * mm, vtxpos.y() * mm, vtxpos.z() * mm));
				G4ThreeVector StartMomentum(aPO.m_px * GeV, aPO.m_py * GeV, aPO.m_pz * GeV);

				particleGun->SetParticleMomentum(StartMomentum);
				fParticleGuns.push_back(particleGun);
			}
		}
	}
	else if(want_muon_background) 
	{
		// generate muon background
		fTPOEvent.run_number = 999;
		G4ParticleDefinition *muon = particleTable->FindParticle("mu-");
		if (muon != nullptr)
		{
			G4ParticleGun *particleGun = new G4ParticleGun(1);
			particleGun->SetParticleDefinition(muon);
			particleGun->SetParticlePosition(G4ThreeVector(vtxpos.x() * mm, vtxpos.y() * mm, vtxpos.z() * mm));
			// Define angular spread (in radians)
			double sigmaTheta = 0.1; // 100 mrad
			// Sample θ from Gaussian centered at 0 with std dev 0.1
			double theta = G4RandGauss::shoot(0.0, sigmaTheta);
			// Sample φ uniformly from 0 to 2π
			double phi = G4UniformRand() * 2.0 * CLHEP::pi;
			// Convert (θ, φ) to Cartesian direction vector
			double px = std::sin(theta) * std::cos(phi);
			double py = std::sin(theta) * std::sin(phi);
			double pz = std::cos(theta);

			// Set direction vector with given momentum magnitude
			double momentumMagnitude = 500.0; // in GeV
			G4ThreeVector StartMomentum(px, py, pz);
			StartMomentum = StartMomentum.unit() * (momentumMagnitude * GeV); // Normalize and scale
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

