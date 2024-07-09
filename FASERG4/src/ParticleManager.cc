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

void ParticleManager::processParticleHit(int trackID, XYZVector const& position, XYZVector const& direction, double const& time,
					 double const& energydeposit, Geant4Process const& process, int const& parentID, int const& pdg,
					 std::string const& VolumeName, int CopyNumber)
{
	if(energydeposit>0) {
		auto it = m_particleMap.find(trackID);
		if (it != m_particleMap.end()) {
			it->second->addTotalEnergyDeposit(energydeposit);
			it->second->update(position, direction, time, energydeposit, process, VolumeName, CopyNumber);
		}
		else {
			// If the photon does not exist, create a new one
			Track* newParticle = new Track(trackID, parentID, pdg, position, direction, time, energydeposit, process, VolumeName, CopyNumber);
			m_particleMap[trackID] = newParticle;
		}
	}
}

#if 0
void ParticleManager::setVertexInformation(int Mode, int NParticles, std::vector<XYZVector> const& VertexPositions,
					   std::vector<XYZVector> const& VertexMomenta, std::vector<double> const& VertexTimes, 
					   std::vector<double> const& VertexEnergies, std::vector<double> const& VertexKinEnergy,
					   std::vector<int> const& VertexTrackID, std::vector<int> const& VertexDecayMode,
					   std::vector<int> const& VertexPDG, double const& vertexNuEnergy)
{
	m_vertexMode = Mode;
	m_vertexNParticles = NParticles;
	m_vertexPositions->clear();
	m_vertexMomenta->clear();
	m_vertexTimes->clear();
	m_vertexEnergy->clear();
	m_vertexKinEnergy->clear();
	m_vertexTrackID->clear();
	m_vertexDecayMode->clear();
	m_vertexPDG->clear();
	m_vertexNuEnergy = vertexNuEnergy;
	for (auto const& pos : VertexPositions) {
		m_vertexPositions->push_back(pos);
	}
	for (auto const& mom : VertexMomenta) {
		m_vertexMomenta->push_back(mom);
	}
	for (auto const& time : VertexTimes) {
		m_vertexTimes->push_back(time);
	}
	for (auto const& energy : VertexEnergies) {
		m_vertexEnergy->push_back(energy);
	}
	for (auto const& kinEnergy : VertexKinEnergy) {
		m_vertexKinEnergy->push_back(kinEnergy);
	}
	for (auto const& trackID : VertexTrackID) {
		m_vertexTrackID->push_back(trackID);
	}
	for (auto const& decayMode : VertexDecayMode) {
		m_vertexDecayMode->push_back(decayMode);
	}
	for (auto const& pdg : VertexPDG) {
		m_vertexPDG->push_back(pdg);
	}
	m_vertexTree->Fill();
}
#endif


void ParticleManager::setDetectorInformation(std::string material1, XYZVector size1, std::string material2, XYZVector size2, int NRep)
{

	m_material1 = material1;
	m_material2 = material2;
	m_size1 = size1;
	m_size2 = size2;
	m_nrep = NRep;
	saveDetector = true;
}

void ParticleManager::beginOfEvent()
{
	// Clear the track map
	m_particleMap.clear();
	fPrimaries.clear();

	PrimaryGeneratorAction* prim = dynamic_cast<PrimaryGeneratorAction*>
		(const_cast<G4VUserPrimaryGeneratorAction*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction()));
	fTcalEvent = new TcalEvent(prim->GetTPOEvent()->run_number, fTcalEvent_number);

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
	fTcalEvent->geom_detector.fSandwichLength = detector->fSandwichLength;
	fTcalEvent->geom_detector.fTotalLength = detector->fTotalLength;
	fTcalEvent->geom_detector.NRep = detector->getNumberReplicas();
	fTcalEvent->geom_detector.fTotalMass = detector->fTotalMass;
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
		fTcalEvent->fillTree();
		fTcalEvent_number++;
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

void ParticleManager::RecordTrack(const G4Track* track) {
	// are we in?
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

