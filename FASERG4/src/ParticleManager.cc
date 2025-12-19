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

// add by Umut
// Field probing for diagnostics
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagneticField.hh"

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
		//std::cout << VolumeName << " copy=" << CopyNumber << " MotherCopy = " << MotherCopyNumber << std::endl;

		const DetectorConstruction *detector = static_cast<const DetectorConstruction *>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
		const XYZVector pos = local_position;
		//////////
    	// Added by Umut: Check that the local position is within expected bounds
		if(m_verbose) {
			double halfX = detector->fRearHCalSizeX/2.0;
			double halfY = detector->fRearHCalSizeY/2.0;
			if (pos.X() < -halfX || pos.X() > halfX || pos.Y() < -halfY || pos.Y() > halfY) {
				G4cout << "[WARN] rearHCAL 'local_position' out of local bounds: "
				<< "pos=(" << pos.X() << "," << pos.Y() << "," << pos.Z() << ") "
				<< " expected X in [" << -halfX << "," << halfX << "]"
				<< " Y in [" << -halfY << "," << halfY << "]" << G4endl;
			}
		}
		/////////
		G4long moduleID = detector->getHCalChannelIDfromXYZ(CopyNumber, pos);
		// DEBUG DUMP for moduleID == 0 (ix=0, iy=0, iz=0)
        if (moduleID == 0 && m_verbose) {
            auto &g = fTcalEvent->geom_detector;
            // Decode indices from moduleID (ix + iy*1000 + iz*1000000)
            long ix = moduleID % 1000;
            long iy = (moduleID / 1000) % 1000;
            long iz = (moduleID / 1000000LL) % 1000;

            // Expected local center using both shift conventions
			double expX_FASER = g.fFASERCal_LOS_shiftX + (ix - g.rearHCalNxy/2.0 + 0.5) * g.rearHCalVoxelSize;
			double expY_FASER = g.fFASERCal_LOS_shiftY + (iy - g.rearHCalNxy/2.0 + 0.5) * g.rearHCalVoxelSize;
			double expZ       = g.rearHCalLocZ + (iz - g.rearHCalNlayer/2.0 + 0.5) * g.rearHCalSizeZ;

			double expX_rear  = g.frearHCal_LOS_shiftX + (ix - g.rearHCalNxy/2.0 + 0.5) * g.rearHCalVoxelSize;
			double expY_rear  = g.frearHCal_LOS_shiftY + (iy - g.rearHCalNxy/2.0 + 0.5) * g.rearHCalVoxelSize;

			// Deltas of the (supposed) local hit vs expected centers
			double dXFASER = pos.X() - expX_FASER;
			double dYFASER = pos.Y() - expY_FASER;
			double dXrear  = pos.X() - expX_rear;
			double dYrear  = pos.Y() - expY_rear;

			// Bounds check for local
			double halfX = g.rearHCalSizeX/2.0;
			double halfY = g.rearHCalSizeY/2.0;
			bool inLocalBounds = (pos.X() >= -halfX && pos.X() <= halfX &&
								pos.Y() >= -halfY && pos.Y() <= halfY);

			// World position from Geant4 tracking step
			const G4ThreeVector worldPosG4(position.x(), position.y(), position.z());

			// Simple Y-tilt rotation (det-frame -> world) using rear HCal shift center
			double theta = g.fTiltAngleY;   // radians
			double c = std::cos(-theta);    // match TcalEvent’s sign convention
			double s = std::sin(-theta);
			double x_det = expX_rear;
			double y_det = expY_rear;
			double z_det = expZ;
			double x_rot =  c * x_det + s * z_det;
			double z_rot = -s * x_det + c * z_det;
			double y_rot =  y_det;

			// World position per display reconstruction (TcalEvent): matrix if available, else Y-tilt helper
			ROOT::Math::XYZVector worldFromTcal = fTcalEvent->getChannelXYZRearHCal(moduleID);

			// Deltas to the Geant4 world
			double dWx_rot = worldPosG4.x()/mm - x_rot;
			double dWy_rot = worldPosG4.y()/mm - y_rot;
			double dWz_rot = worldPosG4.z()/mm - z_rot;

			double dWx_tcal = worldPosG4.x()/mm - worldFromTcal.x();
			double dWy_tcal = worldPosG4.y()/mm - worldFromTcal.y();
			double dWz_tcal = worldPosG4.z()/mm - worldFromTcal.z();

			#if 0
			G4cout
			<< "\n[HCal DEBUG] moduleID=" << moduleID
			<< " ix=" << ix << " iy=" << iy << " iz=" << iz << "\n"
			<< "  Region ranges: ix[0," << (g.rearHCalNxy-1) << "]"
			<< " iy[0," << (g.rearHCalNxy-1) << "]"
			<< " iz[0," << (g.rearHCalNlayer-1) << "]\n"
			<< "  SizeX=" << g.rearHCalSizeX << " mm  SizeY=" << g.rearHCalSizeY << " mm"
			<< "  VoxelSize=" << g.rearHCalVoxelSize << " mm  LayerThickness=" << g.rearHCalSizeZ << " mm\n"
			<< "  CopyNumber(=iz from SD)=" << CopyNumber << "  MotherCopy=" << MotherCopyNumber
			<< "  Volume=" << VolumeName << "\n"
			<< "  LocalPos (from SD)  = (" << pos.X() << ", " << pos.Y() << ", " << pos.Z() << ") mm"
			<< "  inLocalBounds=" << (inLocalBounds ? "yes" : "NO") << "\n"
			<< "  ExpectedLocal(FASERCalLOS)= (" << expX_FASER << ", " << expY_FASER << ", " << expZ << ") mm"
			<< "  d=( " << dXFASER << ", " << dYFASER << " )\n"
			<< "  ExpectedLocal(rearHCalLOS)= (" << expX_rear  << ", " << expY_rear  << ", " << expZ << ") mm"
			<< "  d=( " << dXrear  << ", " << dYrear  << " )\n"
			<< "  WorldPos (G4 step)         = (" << worldPosG4.x()/mm << ", "
													<< worldPosG4.y()/mm << ", "
													<< worldPosG4.z()/mm << ") mm\n"
			<< "  World (Y-tilt on det ctr)  = (" << x_rot << ", " << y_rot << ", " << z_rot << ") mm"
			<< "  Δ(G4 - Y-tilt)=(" << dWx_rot << ", " << dWy_rot << ", " << dWz_rot << ") mm\n"
			<< "  World (TcalEvent recon)    = (" << worldFromTcal.x() << ", "
													<< worldFromTcal.y() << ", "
													<< worldFromTcal.z() << ") mm"
			<< "  Δ(G4 - Tcal)=(" << dWx_tcal << ", " << dWy_tcal << ", " << dWz_tcal << ") mm\n"
			<< "  Note: TcalEvent uses TGeo matrix if available; else a Y-tilt fallback.\n"
			<< "--------------------------------------------------------------\n"
			<< G4endl;
			#endif

        }
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
		#if 0
		// Added by Umut: store hit truth info for rear HCal
		TcalEvent::REARHCALHITTRUTH th{moduleID, track->GetTrackID(), pdg,
                               position.x(), position.y(), position.z(),
                               energydeposit};
		fTcalEvent->rearHCalTruth.push_back(th);
		// Added by Umut for debugging RearHCal truth hits
		//std::cout << " ------> RearHCalTruth: moduleID=" << moduleID 
		//		  << " trackID=" << track->GetTrackID() 
		//		  << " pdg=" << pdg 
		//		  << " x=" << position.x() 
		//		  << " y=" << position.y()
		//		  << " z=" << position.z() 
		//		  << " edep=" << energydeposit << std::endl;
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
		int layerID = CopyNumber; // 0-43
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
		// added by Umut
		// Check probe magnetic field at hit (if requested)
		if (m_verbose) {
			G4ThreeVector gp(position.x(), position.y(), position.z());
			const G4FieldManager* fm = G4TransportationManager::GetTransportationManager()->GetFieldManager();
			G4ThreeVector B(0.,0.,0.);
			if (fm) {
				const G4MagneticField* mf = dynamic_cast<const G4MagneticField*>(fm->GetDetectorField());
				if (mf) {
					G4double point[4] = {gp.x(), gp.y(), gp.z(), 0.};
					G4double Bfield[3] = {0.,0.,0.};
					mf->GetFieldValue(point, Bfield);
					B = G4ThreeVector(Bfield[0], Bfield[1], Bfield[2]);
				}
			}
			G4cout << "SPECT_HIT: SciFi/MuTag pdg=" << pdg << " edep=" << energydeposit << " B[T]=" << B.x()/tesla << "," << B.y()/tesla << "," << B.z()/tesla << " pos[mm]=" << gp.x()/mm <<","<< gp.y()/mm <<","<< gp.z()/mm << G4endl;
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
		// added by Umut
		// Check probe magnetic field at hit (if requested)
		if (m_verbose) {
			G4ThreeVector gp(position.x(), position.y(), position.z());
			const G4FieldManager* fm = G4TransportationManager::GetTransportationManager()->GetFieldManager();
			G4ThreeVector B(0.,0.,0.);
			if (fm) {
				const G4MagneticField* mf = dynamic_cast<const G4MagneticField*>(fm->GetDetectorField());
				if (mf) {
					G4double point[4] = {gp.x(), gp.y(), gp.z(), 0.};
					G4double Bfield[3] = {0.,0.,0.};
					mf->GetFieldValue(point, Bfield);
					B = G4ThreeVector(Bfield[0], Bfield[1], Bfield[2]);
				}
			}
			G4cout << "SPECT_HIT: Magnet pdg=" << pdg << " edep=" << energydeposit << " B[T]=" << B.x()/tesla << "," << B.y()/tesla << "," << B.z()/tesla << " pos[mm]=" << gp.x()/mm <<","<< gp.y()/mm <<","<< gp.z()/mm << " Vol=" << VolumeName << G4endl;
		}	
		return;
	}
	// normal treatment for all other volumes which should essentially be scintillator voxels or tracker hits
	if (VolumeName == "ScintillatorLogical" || VolumeName == "SiTrackerLogical")
	{
		// we save the absolute position in the world coordinate system
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
	fTcalEvent->geom_detector.NRep_SiTracker = detector->fNumberRep_SiTracker;
	fTcalEvent->geom_detector.frearHCal_LOS_shiftX = detector->fRearHCal_LOS_shiftX * mm;
	fTcalEvent->geom_detector.frearHCal_LOS_shiftY = detector->fRearHCal_LOS_shiftY * mm;
	fTcalEvent->geom_detector.fRearMuSpect_LOS_shiftX = detector->fRearMuSpect_LOS_shiftX * mm;
	fTcalEvent->geom_detector.fRearMuSpect_LOS_shiftY = detector->fRearMuSpect_LOS_shiftY * mm;
	//UMUT: copy the GEANT4 tilt angle (radians) so it's saved in the geom branch
    fTcalEvent->geom_detector.fTiltAngleY = detector->GetTiltAngleY();

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

