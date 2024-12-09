#include "TrackerSD.hh"

TrackerSD::TrackerSD(const G4String& name, ParticleManager* i_ParticleManager) : G4VSensitiveDetector(name) { fParticleManager = i_ParticleManager; }

TrackerSD::~TrackerSD() {}

G4bool TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	// energy deposit
	G4double edep = aStep->GetTotalEnergyDeposit();

	// paritcle Type
	G4String particleType = aStep->GetTrack()->GetDefinition()->GetParticleName();

	if (particleType == "opticalphoton") {
		//TODO Here for now I am just killing all the optical photons created. This might not be what we want ot have in the long term
		aStep->GetTrack()->SetTrackStatus(fStopAndKill);
	}
	else {	// here we deal with a particle in geant4.
		XYZVector position = XYZVector(aStep->GetPreStepPoint()->GetPosition().x(), aStep->GetPreStepPoint()->GetPosition().y(),
					     aStep->GetPreStepPoint()->GetPosition().z());

		XYZVector direction = XYZVector(aStep->GetPreStepPoint()->GetMomentum().x(), aStep->GetPreStepPoint()->GetMomentum().y(),
					      aStep->GetPreStepPoint()->GetMomentum().z());
		float time = aStep->GetPreStepPoint()->GetGlobalTime();
		float energyDeposit = aStep->GetTotalEnergyDeposit();

		// if no energy deposited or track is not a primary, then we skip
		G4Track *track = aStep->GetTrack();
		if(energyDeposit < 1e-6 && track->GetParentID()!=0) {
			return true;
		}

//		std::string processName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
//		Geant4Process process = Geant4ProcessMap[processName];

		// std::string VolumeName = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName();

		G4TouchableHandle touchable = aStep->GetPreStepPoint()->GetTouchableHandle();

    	// Get the physical volume
    	G4VPhysicalVolume* physicalVolume = touchable->GetVolume();

		if (physicalVolume == nullptr) { return true; };

		// Get the logical volume
		G4LogicalVolume* logicalVolume = physicalVolume->GetLogicalVolume();

		// Get the volume name
		G4String volumeName = logicalVolume->GetName();

		// Get the material name
		G4String materialName = logicalVolume->GetMaterial()->GetName();

		// Get the copy number
		G4int copyNumber = touchable->GetCopyNumber();

		G4int motherCopyNumber = -1;
        if (touchable->GetHistoryDepth() > 0) {
            motherCopyNumber = touchable->GetCopyNumber(1); // 1 indicates one level up
        }

		// std::cout << process << std::endl;
		int parentID = aStep->GetTrack()->GetParentID();
		int trackID = aStep->GetTrack()->GetTrackID();
		int pdgCode = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();

		fParticleManager->processParticleHit(aStep->GetTrack(), position, direction, time, energyDeposit, 
			parentID, pdgCode, volumeName, copyNumber, motherCopyNumber);
    
		// determine if the particle is leavingthe world volume or is beeing killed/absorbed
		// check with fStopAndKill
		if (!aStep->GetTrack()->GetNextVolume() || aStep->GetTrack()->GetTrackStatus() == fStopAndKill) {
			// Add another photon hit with the post step point
			XYZVector position = XYZVector(aStep->GetPostStepPoint()->GetPosition().x(), aStep->GetPostStepPoint()->GetPosition().y(),
						     aStep->GetPostStepPoint()->GetPosition().z());
			XYZVector direction =
			    XYZVector(aStep->GetPostStepPoint()->GetMomentumDirection().x(), aStep->GetPostStepPoint()->GetMomentumDirection().y(),
				     aStep->GetPostStepPoint()->GetMomentumDirection().z());
			float time = aStep->GetPostStepPoint()->GetGlobalTime();
			float energyDeposit = aStep->GetTotalEnergyDeposit();
//			std::string processName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
//			Geant4Process process = Geant4ProcessMap[processName];
			int parentID = aStep->GetTrack()->GetParentID();
			int trackID = aStep->GetTrack()->GetTrackID();
			int pdgCode = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
			if (!aStep->GetTrack()->GetNextVolume()) {
				fParticleManager->processParticleHit(aStep->GetTrack(), position, direction, time, 
				energyDeposit, parentID, pdgCode,
								     "OutOfWorld",0,0);
			}
			else {
				std::string VolumeName_ = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();

				G4VPhysicalVolume* physicalVolume = aStep->GetPostStepPoint()->GetPhysicalVolume();

				if (physicalVolume == nullptr) { return true; };

				// Get the logical volume
				G4LogicalVolume* logicalVolume = physicalVolume->GetLogicalVolume();

				// Get the volume name
				G4String volumeName = logicalVolume->GetName();

				// Get the material name
				G4String materialName = logicalVolume->GetMaterial()->GetName();

				// Get the copy number
				G4int copyNumber = touchable->GetCopyNumber();

				G4int motherCopyNumber = -1;
				if (touchable->GetHistoryDepth() > 0) {
					motherCopyNumber = touchable->GetCopyNumber(1); // 1 indicates one level up
				}

				fParticleManager->processParticleHit(aStep->GetTrack(), position, direction, time, 
						energyDeposit, parentID, pdgCode,
								     volumeName,copyNumber, motherCopyNumber);
			}

			// fParticleManager->processParticleHit(trackID, position, direction, time, energyDeposit, process, parentID, pdgCode,
			// VolumeName);
		}
	}

	return true;
}
