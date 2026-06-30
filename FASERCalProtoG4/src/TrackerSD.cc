#include "TrackerSD.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

// Initialize static members
std::ofstream TrackerSD::fMuonTrackingFile;
G4bool TrackerSD::fMuonFileOpen = false;

TrackerSD::TrackerSD(const G4String& name, ParticleManager* i_ParticleManager) : G4VSensitiveDetector(name) { fParticleManager = i_ParticleManager; }

TrackerSD::~TrackerSD() {}

void TrackerSD::OpenMuonTrackingFile(const G4String& filename) {
	if (!fMuonFileOpen) {
		fMuonTrackingFile.open(filename.c_str());
		if (fMuonTrackingFile.is_open()) {
			fMuonFileOpen = true;
			// Write header
			fMuonTrackingFile << "# Muon Tracking Data\n";
			fMuonTrackingFile << "# Format: EventID PDG TrackID Edep[MeV] GlobalX[mm] GlobalY[mm] GlobalZ[mm] "
			                  << "LocalX[mm] LocalY[mm] LocalZ[mm] Layer IsActive\n";
			fMuonTrackingFile << "# IsActive: 1=active region, 0=reflective layer\n";
			G4cout << "Opened muon tracking file: " << filename << G4endl;
		} else {
			G4cerr << "Failed to open muon tracking file: " << filename << G4endl;
		}
	}
}

void TrackerSD::CloseMuonTrackingFile() {
	if (fMuonFileOpen) {
		fMuonTrackingFile.close();
		fMuonFileOpen = false;
		G4cout << "Closed muon tracking file" << G4endl;
	}
}

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

		XYZVector local_position = XYZVector(touchable->GetHistory()->GetTopTransform().TransformPoint(
		    aStep->GetPreStepPoint()->GetPosition()));
			// std::cout << " Local position: " << local_position.X() << ", " << local_position.Y() << ", " << local_position.Z() << std::endl;
	
		G4int motherCopyNumber = -1;
		int stepsUp = 1; // go two levels up to get to the top volume since everything is in the detector assembly
        if (touchable->GetHistoryDepth() > 0) {
            motherCopyNumber = touchable->GetCopyNumber(1); // 1 indicates one level up
			if(volumeName == "rearHCalscintillatorLogical") {
				stepsUp = 2;
			}
			local_position = XYZVector(touchable->GetHistory()->GetTransform(stepsUp)
			.TransformPoint(
				aStep->GetPreStepPoint()->GetPosition()));
			// std::cout << " Mother Local position: " << local_position.X() << ", " << local_position.Y() << ", " << local_position.Z() << std::endl;
        }

		// std::cout << process << std::endl;
		int parentID = aStep->GetTrack()->GetParentID();
		int trackID = aStep->GetTrack()->GetTrackID();
		int pdgCode = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();

		// Check if hit is in active voxel region (not in reflective layer)
		// For 3DCAL scintillator, energy deposits in reflective layers should be set to 0
		if(volumeName == "ScintillatorLogical") {
			const DetectorConstruction* detector = static_cast<const DetectorConstruction*>(
				G4RunManager::GetRunManager()->GetUserDetectorConstruction());
			
			bool isActive = detector->isInActiveVoxelRegion(volumeName, local_position, motherCopyNumber);
			// FOR DEBUG ONLY: Dump info about hits in reflective layer and muon steps for visualization
			// Dump all muon steps for 3D visualization
			if(std::abs(pdgCode) == 13) {  // muon or antimuon
				// Get event ID
				G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
				
				// Write to file if open
				if(fMuonFileOpen) {
					fMuonTrackingFile << eventID << " "
					                  << pdgCode << " "
					                  << trackID << " "
					                  << energyDeposit/MeV << " "
					                  << position.X()/mm << " "
					                  << position.Y()/mm << " "
					                  << position.Z()/mm << " "
					                  << local_position.X()/mm << " "
					                  << local_position.Y()/mm << " "
					                  << local_position.Z()/mm << " "
					                  << motherCopyNumber << " "
					                  << (isActive ? 1 : 0) << "\n";
				}
				
				// Also print to console for immediate feedback
				G4cout << "MUON_STEP: EventID=" << eventID
				       << " PDG=" << pdgCode 
				       << " TrackID=" << trackID 
				       << " Edep=" << energyDeposit/MeV << " MeV"
				       << " GlobalPos[mm]=(" << position.X()/mm << "," 
				       << position.Y()/mm << "," 
				       << position.Z()/mm << ")"
				       << " LocalPos[mm]=(" << local_position.X()/mm << "," 
				       << local_position.Y()/mm << "," 
				       << local_position.Z()/mm << ")"
				       << " Layer=" << motherCopyNumber
				       << " IsActive=" << isActive << G4endl;
			}
			
			if(!isActive) {
				// Hit is in the reflective layer - dump info and set energy deposit to 0
				G4cout << "INACTIVE_REGION_HIT: PDG=" << pdgCode 
				       << " TrackID=" << trackID 
				       << " OrigEdep=" << energyDeposit/MeV << " MeV"
				       << " LocalPos[mm]=(" << local_position.X()/mm << "," 
				       << local_position.Y()/mm << "," 
				       << local_position.Z()/mm << ")"
				       << " Layer=" << motherCopyNumber
				       << " -> Edep set to 0" << G4endl;
				energyDeposit = 0;
			}
		}

		fParticleManager->processParticleHit(aStep->GetTrack(), position, local_position, direction,time, energyDeposit, 
			parentID, pdgCode, volumeName, copyNumber, motherCopyNumber);
    
		// determine if the particle is leavingthe world volume or is beeing killed/absorbed
		// check with fStopAndKill
		if (!aStep->GetTrack()->GetNextVolume() || aStep->GetTrack()->GetTrackStatus() == fStopAndKill) {
			// Add another photon hit with the post step point
			XYZVector position = XYZVector(aStep->GetPostStepPoint()->GetPosition().x(), aStep->GetPostStepPoint()->GetPosition().y(),
						     aStep->GetPostStepPoint()->GetPosition().z());
			XYZVector local_position = XYZVector(touchable->GetHistory()->GetTransform(stepsUp)
			.TransformPoint(
				aStep->GetPreStepPoint()->GetPosition()));
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
				fParticleManager->processParticleHit(aStep->GetTrack(), position, local_position, direction, time, 
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
				//////////////////////////////////////////////////
				// Check if hit is in active voxel region (not in reflective layer)
				if(volumeName == "ScintillatorLogical") {
					const DetectorConstruction* detector = static_cast<const DetectorConstruction*>(
						G4RunManager::GetRunManager()->GetUserDetectorConstruction());
					
					bool isActive = detector->isInActiveVoxelRegion(volumeName, local_position, motherCopyNumber);
					
					// Dump all muon exit steps for 3D visualization
					if(std::abs(pdgCode) == 13) {  // muon or antimuon
						// Get event ID
						G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
						
						// Write to file if open
						if(fMuonFileOpen) {
							fMuonTrackingFile << eventID << " "
							                  << pdgCode << " "
							                  << trackID << " "
							                  << energyDeposit/MeV << " "
							                  << position.X()/mm << " "
							                  << position.Y()/mm << " "
							                  << position.Z()/mm << " "
							                  << local_position.X()/mm << " "
							                  << local_position.Y()/mm << " "
							                  << local_position.Z()/mm << " "
							                  << motherCopyNumber << " "
							                  << (isActive ? 1 : 0) << "\n";
						}
						
						// Also print to console for immediate feedback
						G4cout << "MUON_STEP_EXIT: EventID=" << eventID
						       << " PDG=" << pdgCode 
						       << " TrackID=" << trackID 
						       << " Edep=" << energyDeposit/MeV << " MeV"
						       << " GlobalPos[mm]=(" << position.X()/mm << "," 
						       << position.Y()/mm << "," 
						       << position.Z()/mm << ")"
						       << " LocalPos[mm]=(" << local_position.X()/mm << "," 
						       << local_position.Y()/mm << "," 
						       << local_position.Z()/mm << ")"
						       << " Layer=" << motherCopyNumber
						       << " IsActive=" << isActive << G4endl;
					}
					
					if(!isActive) {
						// Hit is in the reflective layer (particle exit/kill) - dump info and set energy deposit to 0
						G4cout << "INACTIVE_REGION_HIT_EXIT: PDG=" << pdgCode 
						       << " TrackID=" << trackID 
						       << " OrigEdep=" << energyDeposit/MeV << " MeV"
						       << " LocalPos[mm]=(" << local_position.X()/mm << "," 
						       << local_position.Y()/mm << "," 
						       << local_position.Z()/mm << ")"
						       << " Layer=" << motherCopyNumber
						       << " -> Edep set to 0" << G4endl;
						energyDeposit = 0;
					}
				}

				// std::cout << " +++++++++++++++++++++++ particle killed : " << VolumeName_ << std::endl;
				fParticleManager->processParticleHit(aStep->GetTrack(), position, local_position, direction, time, 
						energyDeposit, parentID, pdgCode,
								     volumeName,copyNumber, motherCopyNumber);
			}
									
		

			// fParticleManager->processParticleHit(trackID, position, direction, time, energyDeposit, process, parentID, pdgCode,
			// VolumeName);
		}
	}

	return true;
}
