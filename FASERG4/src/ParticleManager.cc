#include "ParticleManager.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Track.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "G4SystemOfUnits.hh"

ParticleManager::ParticleManager(int number) : m_rootFile(nullptr)
{
}

ParticleManager::~ParticleManager() {}

void ParticleManager::processParticleHit(G4Track *track, XYZVector const& position, XYZVector const &local_position, XYZVector const& direction, double const& time,
					 double const& energydeposit, int const& parentID, int const& pdg,
					 std::string const& VolumeName, int CopyNumber, int MotherCopyNumber)
{
	const G4ParticleDefinition *particledef = track->GetParticleDefinition();
	// special treatment for the rear calorimeter hits
	if(VolumeName == "rearCalscintillatorLogical") {
//		std::cout << VolumeName << " copy=" << CopyNumber << " MotherCopy = " << MotherCopyNumber << std::endl;
		auto it = std::find_if(fTcalEvent->rearCalDeposit.begin(), fTcalEvent->rearCalDeposit.end(),
                       [MotherCopyNumber](const TcalEvent::REARCALDEPOSIT& deposit) {
                           return deposit.moduleID == MotherCopyNumber;
                       });
		if (it != fTcalEvent->rearCalDeposit.end()) {
			it->energyDeposit += energydeposit;
		} else {
			struct TcalEvent::REARCALDEPOSIT rearhit = {MotherCopyNumber, energydeposit};
			fTcalEvent->rearCalDeposit.push_back(rearhit);
		}
		return;
	} else if(VolumeName == "rearHCalscintillatorLogical") {
		if(particledef->GetPDGCharge() == 0) return;
		// special treatment for the rear hadronic calorimeter hits
//		std::cout << VolumeName << " copy=" << CopyNumber << " MotherCopy = " << MotherCopyNumber << std::endl;

		const DetectorConstruction *detector = static_cast<const DetectorConstruction *>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
		const XYZVector pos = local_position;
		G4long moduleID = detector->getHCalChannelIDfromXYZ(CopyNumber, pos);

		auto it = std::find_if(fTcalEvent->rearHCalDeposit.begin(), fTcalEvent->rearHCalDeposit.end(),
                       [moduleID](const TcalEvent::REARCALDEPOSIT& deposit) {
                           return deposit.moduleID == moduleID;
                       });
		if (it != fTcalEvent->rearHCalDeposit.end()) {
			it->energyDeposit += energydeposit;
		} else {
			struct TcalEvent::REARCALDEPOSIT rearhit = {moduleID, energydeposit};
			fTcalEvent->rearHCalDeposit.push_back(rearhit);
		}
		#if 0
		G4cout << "HcalScint: pdgid=" << pdg << " edepo=" << energydeposit << G4endl;
		G4double kineticEnergy = track->GetKineticEnergy();
	    G4ThreeVector momentum = track->GetMomentum();
		G4cout << "  Kinetic Energy: " << kineticEnergy / MeV << " MeV" << G4endl;
		G4cout << "  Momentum: " << momentum / MeV << " MeV/c" << G4endl;
		#endif
		return;
	} else if(VolumeName == "muCalscintillatorLogical" || VolumeName == "SciFiLayerLV") {
		// special treatment for the rear muon calorimeter hits
		if(particledef->GetPDGCharge() == 0) return;
		fTcalEvent->rearMuCalDeposit += energydeposit;
		#if 0
		{
		G4cout << "MucalScint: pdgid=" << pdg << " edepo=" << energydeposit << G4endl;
        G4double kineticEnergy = track->GetKineticEnergy();
	    G4ThreeVector momentum = track->GetMomentum();
		G4cout << "  Kinetic Energy: " << kineticEnergy / MeV << " MeV" << G4endl;
        G4cout << "  Momentum: " << momentum / MeV << " MeV/c" << G4endl;
		G4cout << " CopyNumber: " << CopyNumber << " MotherCopyNumber: " << MotherCopyNumber << G4endl;
		}
		#endif
		int trackID = track->GetTrackID();
		int layerID = CopyNumber; // 1-40
		G4ThreeVector momentum = track->GetMomentum();
		auto it = m_MuTagTrackMap.find(trackID);
		if (it != m_MuTagTrackMap.end())
		{
			ROOT::Math::XYZVector pos;
			pos.SetXYZ(position.x(), position.y(), position.z());
			it->second->pos.push_back(pos);
			ROOT::Math::XYZVector mom;
			mom.SetXYZ(momentum.x()/MeV, momentum.y()/MeV, momentum.z()/MeV);
			it->second->mom.push_back(mom);
			it->second->layerID.push_back(layerID);
		}
		else
		{
			// If the track does not exist, create a new one
			MuTagTrack *newT = new MuTagTrack(trackID);
			newT->fPDG = pdg;
			ROOT::Math::XYZVector pos;
			pos.SetXYZ(position.x(), position.y(), position.z());
			newT->pos.push_back(pos);
			ROOT::Math::XYZVector mom;
			mom.SetXYZ(momentum.x()/MeV, momentum.y()/MeV, momentum.z()/MeV);
			newT->mom.push_back(mom);
			newT->layerID.push_back(layerID);
			m_MuTagTrackMap[trackID] = newT;
		}
		return;
	}
	// special treatment for the magnet volume, store all positions
	if(VolumeName == "ShortCylLogical" || VolumeName == "LongCylLogical") {
		if(particledef->GetPDGCharge() == 0) return;

		int trackID = track->GetTrackID();
		auto it = m_magnetTrackMap.find(trackID);
		if (it != m_magnetTrackMap.end())
		{
			ROOT::Math::XYZVector pos;
			pos.SetXYZ(position.x(), position.y(), position.z());
			it->second->pos.push_back(pos);
		}
		else
		{
			// If the track does not exist, create a new one
			MagnetTrack *newT = new MagnetTrack(trackID);
			newT->fPDG = pdg;
			ROOT::Math::XYZVector pos;
			pos.SetXYZ(position.x(), position.y(), position.z());
			newT->pos.push_back(pos);
			m_magnetTrackMap[trackID] = newT;
		}
		return;
	}
	// normal treatment for all other volumes which should essentially be scintillator voxels or tracker hits
	if (VolumeName == "ScintillatorLogical" || VolumeName == "SiTrackerLogical")
	{
		//std::cout << "Position: " << position.x()/mm << " " << position.y()/mm << " " << position.z()/mm << " mm" << std::endl;
		//std::cout << "Local Position: " << local_position.x()/mm << " " << local_position.y()/mm << " " << local_position.z()/mm << " mm" << std::endl;
	
		int trackID = track->GetTrackID();
		if (energydeposit > 0 || parentID == 0)
		{ // or some energy or is primary
			auto it = m_particleMap.find(trackID);
			if (it != m_particleMap.end())
			{
				it->second->addTotalEnergyDeposit(energydeposit);
				it->second->update(local_position, direction, time, energydeposit, VolumeName, CopyNumber, MotherCopyNumber);
			}
			else
			{
				// If the track does not exist, create a new one
				Track *newParticle = new Track(trackID, parentID, pdg, local_position, direction, time,
											   energydeposit, VolumeName, CopyNumber, MotherCopyNumber);
				m_particleMap[trackID] = newParticle;
			}
		}
	}
}

void ParticleManager::setDetectorInformation(std::string material1, XYZVector size1, std::string material2, XYZVector size2, int NRep)
{

	m_material1 = material1;
	m_material2 = material2;
	m_size1 = size1;
	m_size2 = size2;
	m_nrep = NRep;
//	saveDetector = true;
}

void ParticleManager::beginOfEvent()
{
	// Clear the track map
	m_particleMap.clear();
	m_magnetTrackMap.clear();
	m_MuTagTrackMap.clear();
	fPrimaries.clear();

	PrimaryGeneratorAction* prim = dynamic_cast<PrimaryGeneratorAction*>
		(const_cast<G4VUserPrimaryGeneratorAction*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction()));
	int mask_event = const_cast<TPOEvent*>(prim->GetTPOEvent())->event_mask;
	fTcalEvent = new TcalEvent(prim->GetTPOEvent()->run_number, prim->GetTPOEvent()->event_id, mask_event);

	// set the link to the primary event information
//	PrimaryGeneratorAction* prim = dynamic_cast<PrimaryGeneratorAction*>
//		(const_cast<G4VUserPrimaryGeneratorAction*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction()));
	fTcalEvent -> fTPOEvent = const_cast<TPOEvent*>(prim->GetTPOEvent());

	fTcalEvent -> fTPOEvent -> setPrimaryVtx(primary_vertex_position.X(), primary_vertex_position.Y(), primary_vertex_position.Z()) ;

	// setup the geometry
	const DetectorConstruction* detector = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
	fTcalEvent->geom_detector.fScintillatorSizeX = detector->fScintillatorSizeX;
	fTcalEvent->geom_detector.fScintillatorSizeY = detector->fScintillatorSizeY;
	fTcalEvent->geom_detector.fScintillatorVoxelSize = detector->fScintillatorVoxelSize;
	fTcalEvent->geom_detector.fSiTrackerSizeZ = detector->fSiTrackerSizeZ;
	fTcalEvent->geom_detector.fSiTrackerPixelSize = detector->fSiTrackerPixelSize;
	fTcalEvent->geom_detector.fTargetSizeZ = detector->ftargetWSizeZ;
	fTcalEvent->geom_detector.fSandwichLength = detector->fSandwichLength;
	fTcalEvent->geom_detector.fTotalLength = detector->fTotalLength;
	fTcalEvent->geom_detector.NRep = detector->getNumberReplicas();
	fTcalEvent->geom_detector.fTotalMass = detector->fTotalMass;
	fTcalEvent->geom_detector.fTotalWmass = detector->fTotalWMass;
	fTcalEvent->geom_detector.fTotalScintmass = detector->fTotalScintMass;
	fTcalEvent->geom_detector.rearCalSizeX = 121.2;
	fTcalEvent->geom_detector.rearCalSizeY = 121.2;
    fTcalEvent->geom_detector.rearCalLocZ = detector->fTotalLength/2.0; // when no magnet + 3500.0;
	fTcalEvent->geom_detector.rearCalNxy = 5;
	fTcalEvent->geom_detector.rearHCalSizeX = 720.0; // mm
	fTcalEvent->geom_detector.rearHCalSizeY = 720.0; // mm
	fTcalEvent->geom_detector.rearHCalSizeZ = 23; // mm
	fTcalEvent->geom_detector.rearHCalVoxelSize = detector->fRearHCalVoxelSize;
	fTcalEvent->geom_detector.rearHCalLocZ =  fTcalEvent->geom_detector.rearCalLocZ + 66*6.0; 
	fTcalEvent->geom_detector.rearHCalNxy = 18;
	fTcalEvent->geom_detector.rearHCalLength = detector->fRearHCalLength;
	fTcalEvent->geom_detector.rearHCalNlayer = 40;
	fTcalEvent->geom_detector.rearMuSpectLocZ = detector->fRearMuSpectLocZ;
	fTcalEvent->geom_detector.rearMuSpectSizeZ = detector->fRearMuSpectSizeZ;
	fTcalEvent->geom_detector.fFASERCal_LOS_shiftX = detector->fFASERCal_LOS_shiftX * mm;
	fTcalEvent->geom_detector.fFASERCal_LOS_shiftY = detector->fFASERCal_LOS_shiftY * mm;
	fTcalEvent->geom_detector.fAirGap = detector->fAirGap * mm;
	fTcalEvent->geom_detector.fAlPlateThickness = detector->fAlPlateThickness * mm;
	fTcalEvent->geom_detector.fSiTrackerGap = detector->fSiTrackerGap * mm;

	// clear the rear calorimeter
	fTcalEvent->rearCalDeposit.clear();
	// clear the rear hcal scintillator
	fTcalEvent->rearHCalDeposit.clear();
	// clear the rear muCal scintillator
	fTcalEvent->rearMuCalDeposit = 0;
}

void ParticleManager::endOfEvent(G4Event const* event)
{
	// assign geantTrackID to the TPOEvent entry
	G4int number_primary_vtx = event->GetNumberOfPrimaryVertex();
	G4cout << "Number of primary vertices " << number_primary_vtx << G4endl;
	for (int i = 0; i < number_primary_vtx; i++){
		G4PrimaryVertex *pv = event->GetPrimaryVertex(i);
		G4int number_particles = pv->GetNumberOfParticle();
		for (int j = 0; j < number_particles; j++) {
			G4PrimaryParticle *pp = pv->GetPrimary(j);
			fTcalEvent -> AssignGEANTTrackID(pp->GetTrackID(), pp->GetPDGcode(), 
				pp->GetPx()/CLHEP::GeV, pp->GetPy()/CLHEP::GeV, pp->GetPz()/CLHEP::GeV);
		}
	}

	int numberOfTracks = m_particleMap.size();
	std::cout << "Number of particles: " << numberOfTracks << std::endl;

	if(numberOfTracks > 0) {
		// Fill tCalEvent and dump the primaries
		for (const auto& it : m_particleMap) {
			auto track = it.second;
		//	if(track->getParentID() == 0) track->Dump();

			DigitizedTrack* t = fTcalEvent->addTrack(track->getTrackID());
			t->fparentID = track->getParentID();
			t->fPDG = track->getPDG();
			t->fprimaryID = FindPrimaryParticle(track->getTrackID());
			for (const auto& hit : track->getmhitIDMap()) {
				t->fhitIDs.push_back(hit.first);
				t->fEnergyDeposits.push_back(hit.second);
			}
		}

		// Fill magnet tracks
		for (const auto& it : m_magnetTrackMap) {
			MagnetTrack* magnettrk = it.second;
			fTcalEvent->fMagnetTracks.push_back(magnettrk);
		}
	
		// Fill mutag tracks
		for (const auto& it : m_MuTagTrackMap) {
			MuTagTrack* mutagtrk = it.second;
			fTcalEvent->fMuTagTracks.push_back(mutagtrk);
		}

		// dump RearCal
		#if 0
		for (const auto &it : fTcalEvent->rearCalDeposit) {
			G4cout << " Module " << it.moduleID << " deposit: " << it.energyDeposit << G4endl;
		}
		// dump rear muCal
		G4cout << " Rear Mu Scintillator: " << fTcalEvent->rearMuCalDeposit << G4endl;
		#endif
		
		fTcalEvent->fillTree();
	//	fTcalEvent_number++;
	}
	delete fTcalEvent;
	fTcalEvent = nullptr;

#if 0
	// what is this???
	for( int i = 0; i < m_numberParticleVectors; i++){
		m_particleVectorOutput->at(i).clear();
	}
	int i = 0;
	int vectorIt = -1;
	for (auto& entry : m_particleMap) {
		if(i % m_numberParticlesPerVector == 0){
			vectorIt++;
		}
		m_particleVectorOutput->at(vectorIt).push_back(*(entry.second));
		int trackID = entry.second->getTrackID();
		m_particleIDMap->insert(std::pair<int, int>(trackID, i));
		i++;
	}

	// Write the m_rayVector to the ROOT file
	m_particleTree->Fill();
#endif

	for (auto& entry : m_particleMap) {
		delete entry.second;
	}
	
	m_particleMap.clear();
}

void ParticleManager::beginOfRun()
{
	fPrimaries.clear();
}

void ParticleManager::endOfRun()
{
}

/// @brief Record track into the map to link trackID to primaries and create Track
/// @param track 
void ParticleManager::RecordTrack(const G4Track* track) {
	// are we in the map that will allow us to easily find primary track of any track?
	int trackID = track->GetTrackID();
	auto it = fPrimaries.find(trackID);
	if(it == fPrimaries.end()) {
		struct PrimaryInfo pinfo;
		pinfo.parentID = track->GetParentID();
		if(pinfo.parentID > 0) {
			pinfo.primaryID = FindPrimaryParticle(pinfo.parentID);
		} else {
			pinfo.primaryID = trackID;     // it's myself
		}
		fPrimaries[trackID] = pinfo;		
	}

	// handle case of neutral primary particle that don't deposit energy
	// this is to take care of the case where the neutral particle never 
	// triggers a TrackerSD hit
	G4int parentID = track->GetParentID();
	if(parentID == 0) {
		const G4ParticleDefinition *particle = track->GetParticleDefinition();
		if(particle->GetPDGCharge() == 0){
			auto it = m_particleMap.find(trackID);
			if (it == m_particleMap.end()) {
				G4int pdg = particle->GetPDGEncoding();

				// insert a Track but skip final state neutrinos
				if(!fTcalEvent -> fTPOEvent->is_neutrino(pdg)) {
					XYZVector position(0, 0, 0);
					XYZVector direction(0, 0, 0);
					double Time = 0;
					double energydeposit = 0;

					Track *newParticle = new Track(trackID, parentID, pdg,
												   position, direction, Time, energydeposit, "", 0,0);
					m_particleMap[trackID] = newParticle;
				}
			}
		}
	}
}

G4int ParticleManager::FindPrimaryParticle(G4int trackID)
{
	auto it = fPrimaries.find(trackID);
	if(it != fPrimaries.end()) {
		return it->second.primaryID;
	}
	G4cerr << "FindPrimaryParticle failed ..." << G4endl;
	return -1;
}

