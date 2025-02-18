#include "DetectorConstruction.hh"
#include "G4ios.hh"
#include "G4GDMLParser.hh"

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

DetectorConstruction::DetectorConstruction(ParticleManager* photonManager)
{
	fMessenger = new DetectorMessenger(this);

	fParticleManager = photonManager;
}

DetectorConstruction::~DetectorConstruction()
{
	delete fStepLimit;
	delete fMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	// Define materials
	DefineMaterials();

	// Define volumes
	return DefineVolumes();
}

void DetectorConstruction::DefineMaterials()
{
	// Material definition

	G4NistManager* nistManager = G4NistManager::Instance();

	// Air defined using NIST Manager
	if (fWorldMaterial == nullptr) {
		fWorldMaterial = nistManager->FindOrBuildMaterial("G4_AIR");
		fAir_MPT = new G4MaterialPropertiesTable();
		fAir_MPT->AddProperty("RINDEX", fPhotonEnergyAir, fRefractiveIndex_Air, nEntriesAir);
		fWorldMaterial->SetMaterialPropertiesTable(fAir_MPT);
	}

	// Build PVT
	// check if the material is defined already, if not define it
	// This check is necessary because this method is invoked in the messeneger call funciton to set the material
	// If this check is not done, it will be invoked twice and the second time it will cause an error
	if (fPolyvinyltoluene == nullptr) {
		// Build PVT (PolyVinylToluene) from C and H elements
		fHydrogen = new G4Element("Hydrogen", "H", 1., 1.01 * g / mole);
		fCarbon = new G4Element("Carbon", "C", 6., 12.01 * g / mole);

		// PVT (PolyVinylToluene, C_9H_10)
		double HAtomsPerVolume = 10;
		double CAtomsPerVolume = 9;

		double HPerCent = HAtomsPerVolume / (HAtomsPerVolume + CAtomsPerVolume);
		double CPerCent = CAtomsPerVolume / (HAtomsPerVolume + CAtomsPerVolume);

		fPolyvinyltoluene = new G4Material("PVT", 1.03 * g / cm3, 2);
		fPolyvinyltoluene->AddElement(fHydrogen, HPerCent);
		fPolyvinyltoluene->AddElement(fCarbon, CPerCent);

		fPolyvinyltoluene_MPT = new G4MaterialPropertiesTable();
		fPolyvinyltoluene_MPT->AddConstProperty("SCINTILLATIONYIELD", fLightYield / MeV);  // The light yield of the scintillator
		fPolyvinyltoluene_MPT->AddConstProperty(
		    "RESOLUTIONSCALE",
		    1.0);  // This determines the width of the distribution of number of photons produced. If the number is lower than 10,
			   // it is sampled from a poission. If it is larger, a gaussian with a widht of 1/sqrt(N)*resScale is used
		fPolyvinyltoluene_MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", fScintillationDecayTime * ns);  // The decay time

		fPolyvinyltoluene_MPT->AddProperty(
		    "RINDEX", fPhotonEnergyPVT, fRefractiveIndex_PVT,
		    nEntriesPVT);  // The refractive index of the scintillator, see in the header file for the valuesTODO
		fPolyvinyltoluene_MPT->AddProperty(
		    "ABSLENGTH", fPhotonEnergyPVT, fAbsorption_PVT,
		    nEntriesPVT);  // The absorption length of the scintillator, see in the header file for the valuesTODO
		fPolyvinyltoluene_MPT->AddProperty(
		    "SCINTILLATIONCOMPONENT1", fPhotonEnergyPVT, fScintillation_PVT,
		    nEntriesPVT);  // The scintillation spectrum of the scintillator, see in the header file for the valuesTODO
		fPolyvinyltoluene->SetMaterialPropertiesTable(fPolyvinyltoluene_MPT);
	
	// 0.898e-2 g/cm^2/MeV
	// M.Hirschberg et al., IEEE Trans. Nuc. Sci. 39 (1992) 511
    // SCSN-38: kB = (0.806 +/- 0.012)E-2 g/cm^2/MeV
    // SCSN-28: kB = (0.877 +/- 0.03)E-2 g/cm^2/MeV
    // GS 2003: kB = (0.844 +/- 0.015)E-2 g/cm^2/MeV
    // NE 102A: kB = (0.882 +/- 0.012)E-2 g/cm^2/MeV
    // NE 102A: kB = (0.888 +/- 0.025)E-2 g/cm^2/MeV
    // KSTI 390: kB = (1.09 +/- 0.015)E-2 g/cm^2/MeV
    // Simple average: (0.898 +/- 0.099)E-2 g/cm^2/MeV
    // Average excluding KSTI 390: (0.859 +/- 0.034)E-2 g/cm^2/MeV
    //
    // Also: M. Bongrand AIP Conf. Proc. 807 (2007) 14: kB = 9E-2 g/cm^2/MeV
    //
    // Wikipedia give a value of 0.126 mm/MeV for polystyrene based
    // scintillator, and this value seems to pop up in various places.  It
    // comes from a measurement of 1 mm scintillating fibers.  This value also
    // appears in Geant4 Examples.
    //
    // Leverington, Anelli, Campana, and Rosellini, arxiv:1106.5649 (2011))
    //
    // As best I can tell, the value was not fit to data, and is simply what
    // was used in their simulation.  I don't find any supporting
    // information in the paper.

		// The Birks constant of the scintillator
		double birks_constant = (0.898e-2 * g / cm / cm / MeV);
		fPolyvinyltoluene->GetIonisation()->SetBirksConstant(8.718e-3 * cm / MeV);	
	}

	// Print materials
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
	// Sizes of the principal geometrical components (solids)
	// The values are given via the messenger, same as the units



	G4int NRep = getNumberReplicas();
	G4double sizetargetWX = gettargetWSizeX();
	G4double sizetargetWY = gettargetWSizeY();
	G4double sizetargetWZ = gettargetWSizeZ();
	G4double sizeScintillatorX = getScintillatorSizeX();
	G4double sizeScintillatorY = getScintillatorSizeY();
	G4double sizeScintillatorZ = getScintillatorSizeZ();

	G4double sizeX = std::max(sizetargetWX, sizeScintillatorX);
	G4double sizeY = std::max(sizetargetWY, sizeScintillatorY);

	G4double WorldSizeX = 2.5* sizeX;
	G4double WorldSizeY = 2.5* sizeY;
	G4double WorldSizeZ = 6*m; // 1.2* NRep*(sizetargetWZ + sizeScintillatorZ);

	G4cout << "Size of the world " << WorldSizeX << " " << WorldSizeY << " " << WorldSizeZ << " mm" << G4endl;

	// get the maximum size of the detector
	G4double maxDetectorSize = std::max(std::max(WorldSizeX, WorldSizeY), WorldSizeZ);

	// Definitions of Solids, Logical Volumes, Physical Volumes

	G4GeometryManager::GetInstance()->SetWorldMaximumExtent(maxDetectorSize);

	auto worldSolid = new G4Box("world",					      // its name
				    WorldSizeX / 2, WorldSizeY / 2, WorldSizeZ );  // its size
	auto worldLV = new G4LogicalVolume(worldSolid,				      // its solid
					   fWorldMaterial,			      // its material
					   "World");				      // its name

	//  Must place the World Physical volume unrotated at (0,0,0).
	//
	auto worldPV = new G4PVPlacement(nullptr,	   // no rotation
					 G4ThreeVector( 0,0,0*cm), 
					 worldLV,	   // its logical volume
					 "World",	   // its name
					 nullptr,	   // its mother  volume
					 false,		   // no boolean operations
					 0,		   // copy number
					 fCheckOverlaps);  // checking overlaps

	// target is composed of W or Copper
	G4Material * G4_W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
	G4Material * G4_Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
	G4Material * G4_graphite = G4NistManager::Instance()->FindOrBuildMaterial("G4_GRAPHITE");
	G4Material * G4_Pb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");

	G4Material * G4_Target = G4_W;

	G4double zLocation = 0*cm;
	CreateFaserCal(zLocation, fPolyvinyltoluene, G4_Target, G4ThreeVector(sizeScintillatorX, sizeScintillatorY, 
	sizeScintillatorZ),G4ThreeVector(sizetargetWX, sizetargetWY, sizetargetWZ),worldLV, NRep);
	fParticleManager->setDetectorInformation(fPolyvinyltoluene->GetName(),
	XYZVector(sizeScintillatorX, sizeScintillatorY, sizeScintillatorZ), G4_Target->GetName(),  
	XYZVector(sizetargetWX, sizetargetWY, sizetargetWZ), NRep);

//	G4double sizeZ = (sizeScintillatorZ + sizetargetWZ + fNumberRep_SiTracker*fSiTrackerSizeZ)*NRep;
	G4double sizeZ = fTotalLength;

//	CreateFrontTarget(-sizeZ/2.0-20.0*cm, worldLV);

	zLocation += sizeZ/2.0;

#if magnet
	CreateMagnetSystem(zLocation, worldLV);
	zLocation += 350.0*cm;    // length of the magnet system
#endif

	CreateRearCal(zLocation, worldLV);

	G4double sizeZmu = zLocation + 66*6*mm;
	CreateRearHCalMuTag(sizeZmu, worldLV);

	// Save the geometry of the detector
	G4GDMLParser parser;
	parser.Write("geometry.gdml", worldPV->GetLogicalVolume());

	return worldPV;
}

// This function makes the block of scintillator sensitive to the particles
void DetectorConstruction::ConstructSDandField()
{
	// Sensitive detectors

	G4String FASERCalSDname = "/FaserCalSD";
	auto aTrackerSD = new TrackerSD(FASERCalSDname, fParticleManager);
	G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
	SetSensitiveDetector("ScintillatorLogical", aTrackerSD, true);
	SetSensitiveDetector("trackerSiLogical", aTrackerSD, false);
	SetSensitiveDetector("rearCalscintillatorLogical", aTrackerSD, false);
	SetSensitiveDetector("rearHCalscintillatorLogical", aTrackerSD, false);
	SetSensitiveDetector("muCalscintillatorLogical", aTrackerSD, false);

#if magnet
	SetSensitiveDetector("ShortCylLogical", aTrackerSD, false);
	SetSensitiveDetector("LongCylLogical", aTrackerSD, false);
#endif

	//SetSensitiveDetector("World", aTrackerSD, true);

	// Create global magnetic field messenger.
	// Uniform magnetic field is then created automatically if
	// the field value is not zero.
	G4ThreeVector fieldValue = G4ThreeVector();
	fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
	fMagFieldMessenger->SetVerboseLevel(1);

	// Register the field messenger for deleting
	G4AutoDelete::Register(fMagFieldMessenger);
}

// These set functions are used to set the values of the detector size and the material of the scintillator via the messenger
void DetectorConstruction::SetMaxStep(G4double maxStep)
{
	if ((fStepLimit) && (maxStep > 0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

void DetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps) { fCheckOverlaps = checkOverlaps; }

void DetectorConstruction::SetLightYield(G4double lightYield)
{
	fLightYield = lightYield;
	fScintillatorMaterial->GetMaterialPropertiesTable()->AddConstProperty("SCINTILLATIONYIELD", fLightYield / MeV);
}

void DetectorConstruction::SetScintillationDecayTime(G4double scintillatorDecayTime)
{
	fScintillationDecayTime = scintillatorDecayTime;
	fScintillatorMaterial->GetMaterialPropertiesTable()->AddConstProperty("SCINTILLATIONTIMECONSTANT", fScintillationDecayTime / ns);
}

void DetectorConstruction::SetScintillatorSizeX(G4double detectorSizeX) { fScintillatorSizeX = detectorSizeX; }
void DetectorConstruction::SetScintillatorSizeY(G4double detectorSizeY) { fScintillatorSizeY = detectorSizeY; }
void DetectorConstruction::SetScintillatorSizeZ(G4double detectorSizeZ) { fScintillatorSizeZ = detectorSizeZ; }

void DetectorConstruction::SettargetWSizeX(G4double detectorSizeX) { ftargetWSizeX = detectorSizeX; }
void DetectorConstruction::SettargetWSizeY(G4double detectorSizeY) { ftargetWSizeY = detectorSizeY; }
void DetectorConstruction::SettargetWSizeZ(G4double detectorSizeZ) { ftargetWSizeZ = detectorSizeZ; }
void DetectorConstruction::SetNumberReplicas(G4int rep) { fNumberReplicas = rep; }



void DetectorConstruction::SetScintillatorMaterial(G4String materialName)
{
	G4Material* pttoMaterial = nullptr;

	if (materialName ==
	    "PVT") {  // PVT is a custom material, hence it is not defined in the G4NistManager and we have to take care of it ourselves
		DetectorConstruction::DefineMaterials();
		G4cout << fPolyvinyltoluene->GetName() << G4endl;
		pttoMaterial = fPolyvinyltoluene;
	}
	else {
		G4NistManager* nistManager = G4NistManager::Instance();
		pttoMaterial = nistManager->FindOrBuildMaterial(materialName);
	}

	if (fScintillatorMaterial != pttoMaterial) {
		if (pttoMaterial) {
			fScintillatorMaterial = pttoMaterial;
			if (fLogicScintillator) fLogicScintillator->SetMaterial(fScintillatorMaterial);
			G4cout << G4endl << "----> The scintillator is made of " << materialName << G4endl;
		}
		else {
			G4cout << G4endl << "-->  WARNING from SetTargetMaterial : " << materialName << " not found" << G4endl;
		}
	}
}

void DetectorConstruction::CreateFaserCal(G4double zLocation, G4Material* material1, G4Material* material2, G4ThreeVector size1, 
		G4ThreeVector size2, G4LogicalVolume* parent, G4int NRep){

	G4Material * G4_Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

	// add precision tracker of composed of Si
	G4Material * G4_Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
	G4cout << "Thickness of a single Silicon tracker layer" << fSiTrackerSizeZ << " mm " << G4endl;
	G4cout << "Number of silicon tracker layer " << fNumberRep_SiTracker << G4endl;

	G4double sizeX = std::max(size1.getX(), size2.getX());
	G4double sizeY = std::max(size1.getY(), size2.getY());
	G4double sizeZ = size1.getZ() + size2.getZ() + 	fNumberRep_SiTracker*fSiTrackerSizeZ;

	fSandwichLength = fAlPlateThickness + size1.getZ() + size2.getZ() + fAlPlateThickness + fAirGap;

	G4cout << "Total size of one sandwich layer " << fSandwichLength*mm << " mm " << G4endl;

	G4cout << "Number of layers " << NRep << G4endl;

	fTotalLength = NRep*fSandwichLength;
	G4cout << "Total length " << fTotalLength << " mm" << G4endl;

	G4double density = material2->GetDensity()/(g/cm3);  // Density in g/cm^3
	fTotalWMass = sizeX*sizeY*size2.getZ()*density*1e-3*NRep*1e-3;
	G4cout << "Total mass target (W, Cu, ...) " << fTotalWMass << " kg" << G4endl;

	fTotalScintMass = sizeX*sizeY*size1.getZ()*1.03e-3*NRep*1e-3;
	G4cout << "Total mass scint " << fTotalScintMass << " kg" << G4endl;

	fTotalMass = fTotalWMass + fTotalScintMass;
	G4cout << "Total mass target+scint " << fTotalMass << " kg" << G4endl;

	G4double radlen1 = material1->GetRadlen()/(mm);
	G4double interlen1 = material1->GetNuclearInterLength()/(mm);
	G4double radlen2 = material2->GetRadlen()/(mm);
	G4double interlen2 = material2->GetNuclearInterLength()/(mm);
	G4cout << "Radiation length scintillator layer " << size1.getZ()/radlen1 << G4endl;
	G4cout << "Interaction length scintillator layer " << size1.getZ()/interlen1 << G4endl;
	G4cout << "Radiation length target layer " << size2.getZ()/radlen2 << G4endl;
	G4cout << "Interaction length target layer " << size2.getZ()/interlen2 << G4endl;

	G4cout << "Total rad length " << NRep*(size1.getZ()/radlen1+size2.getZ()/radlen2) << G4endl;
	G4cout << "Total interac length " << NRep*(size1.getZ()/interlen1+size2.getZ()/interlen2) << G4endl;

	// TODO Setup via macro commands - currently set to large values for simulation speed.
	G4UserLimits* userLimits_Scint = new G4UserLimits();
	userLimits_Scint->SetMaxAllowedStep(2*mm); // necessary; should be tuned to size of voxel
	G4UserLimits* userLimits_targetW = new G4UserLimits();
//	userLimits_targetW->SetMaxAllowedStep(1*mm);   // is this really necessary??

	// In this replica solid we place the actual materials. This replica solid is then replicated NRep times.
    G4Box * replicaSolid = new G4Box("ReplicaSolid" , sizeX/2, sizeY/2, fSandwichLength/2);

    G4Box* scintillatorSolid = new G4Box("ScintillatorSlab", size1.getX() / 2, size1.getY() / 2, size1.getZ() / 2);
	G4Box* targetWSolid = nullptr;
	if(size2.getZ() > 0)
	    targetWSolid = new G4Box("targetWSlab", size2.getX() / 2, size2.getY() / 2, size2.getZ() / 2);
    G4Box* trackerSiSolid = new G4Box("trackerSiSlab", size2.getX() / 2, size2.getY() / 2, fSiTrackerSizeZ / 2);
    G4Box* AlPlateSolid = new G4Box("AlPlateSlab", size2.getX() / 2, size2.getY() / 2, fAlPlateThickness / 2);

	G4LogicalVolume* replicaLogic = new G4LogicalVolume(replicaSolid, fWorldMaterial, "replicaLogical");

	G4LogicalVolume* scintillatorLogic = new G4LogicalVolume(scintillatorSolid, material1, "ScintillatorLogical");
	scintillatorLogic->SetUserLimits(userLimits_Scint);
	G4LogicalVolume* targetWLogic = nullptr;
	if(size2.getZ() > 0) {
		targetWLogic = new G4LogicalVolume(targetWSolid, material2, "targetWLogical");
		targetWLogic->SetUserLimits(userLimits_targetW);
	}
	G4LogicalVolume* trackerSiLogic = new G4LogicalVolume(trackerSiSolid, G4_Si, "trackerSiLogical");
	G4LogicalVolume* AlPlateLogic = new G4LogicalVolume(AlPlateSolid, G4_Al, "AlPlateLogical");


    //We place the scinitllator and targetW into the replica
	double zShift = -fSandwichLength/2.0 + fAlPlateThickness/2.0;
	new G4PVPlacement(0, G4ThreeVector(0,0,zShift), AlPlateLogic, "AlPlate", replicaLogic, false, 0, true);
	zShift += fAlPlateThickness/2.0 + size2.getZ()/2.0;
	if(size2.getZ() > 0)
	    new G4PVPlacement(0, G4ThreeVector(0,0,zShift), targetWLogic, "targetW", replicaLogic, false, 1, true);
	zShift += size2.getZ()/2.0 + size1.getZ()/2.0;
    new G4PVPlacement(0, G4ThreeVector(0,0,zShift), scintillatorLogic, "Scintillator", replicaLogic, false, 1, true);
	zShift += size1.getZ()/2.0 + fAlPlateThickness/2.0;
	new G4PVPlacement(0, G4ThreeVector(0,0,zShift), AlPlateLogic, "AlPlate", replicaLogic, false, 1, true);

	zShift = fSandwichLength/2.0;
	#if 0
    new G4PVPlacement(0, G4ThreeVector(0,0,zShift - fSiTrackerGap - fSiTrackerSizeZ/2), 
					trackerSiLogic, "trackerSi", replicaLogic, false, 1, true);
    new G4PVPlacement(0, G4ThreeVector(0,0,zShift - fSiTrackerSizeZ/2), trackerSiLogic, "trackerSi", replicaLogic, false, 2, true);
	#endif

    // Create container, which hosts Nrep replicas of our layer. This container then is set inside the world volume
    G4Box* containerSolid = new G4Box("ContainerBox", sizeX/2, sizeY / 2, NRep*fSandwichLength / 2);

    G4LogicalVolume* containerLogic = new G4LogicalVolume(containerSolid, fWorldMaterial, "ContainerLogical");
//    new G4PVReplica("replica", replicaLogic, containerLogic, kZAxis, NRep, sizeZ, 0);
//    new G4PVPlacement(0, G4ThreeVector(0,0,NRep*sizeZ/2), containerLogic, "ContainerPlacement", parent,  false, 0, true);
	// make nRep copies
	for(int i = 0; i < NRep; i++) {
		double zshift = i*fSandwichLength-fSandwichLength*NRep/2.0+fSandwichLength/2.0;
		std::cout << "Placing first replica at " << zshift << std::endl;
		new G4PVPlacement(0, G4ThreeVector(0,0,zshift), replicaLogic, "replica", containerLogic, false, i, true);
	}
    new G4PVPlacement(0, G4ThreeVector(
		fFASERCal_LOS_shiftX,fFASERCal_LOS_shiftY,zLocation), 
		containerLogic, "ContainerPlacement", parent,  false, 0, true);
}

static int getchannelIDerrorcount = 0;

G4long DetectorConstruction::getChannelIDfromXYZ(std::string const& VolumeName, int CopyNumber, int MotherCopyNumber, XYZVector const& position) const {

	G4double dx = position.X()-fFASERCal_LOS_shiftX+fScintillatorSizeX/2.0;
	G4double dy = position.Y()-fFASERCal_LOS_shiftY+fScintillatorSizeY/2.0;
	G4double epsilon = 1e-4;   // avoid rounding errors at volume boundary
	G4double dz = position.Z()+fTotalLength/2.0+epsilon;

	// sanity check
	if((dx < 0 || dx > fScintillatorSizeX) || (dy < 0 || dy > fScintillatorSizeY) || (dz < 0 || dz > fTotalLength)) {
		getchannelIDerrorcount++;
		if(getchannelIDerrorcount < 100) {
			G4cerr << "ERROR : getCHannelIDfromXYZ problem dx:" << dx << " dy:" << dy << " dz:" << dz << G4endl;
		}
		return 0;
	}

//	G4long ilayer = floor(dz / fSandwichLength);
//	if( ilayer != CopyNumber) {
//		G4cerr << "ERROR getChannelIDfromXYZ " << VolumeName << " - ilayer " << ilayer << " copyNumber " << CopyNumber << G4endl;
//		return 0;
//	}

	G4long ilayer = MotherCopyNumber;

	if(VolumeName == "ScintillatorLogical") {

		G4long ix = floor(dx / fScintillatorVoxelSize);
		G4long iy = floor(dy / fScintillatorVoxelSize);
		G4long iz = floor((dz-ilayer*fSandwichLength-fAlPlateThickness-ftargetWSizeZ)/ fScintillatorVoxelSize);
		// sanity check
		if (ix < 0 || ix > 999 || iy < 0 || iy > 999 || iz < 0 || iz > 999 || ilayer < 0) {
			if(iz>1000)
				G4cerr << "ERROR : getCHannelIDfromXYZ problem ix:" << ix << " iy:" << iy << " iz:" << iz << G4endl;
			return 0;
		}
		G4long ID = ix + iy*1000 + iz*1000000 + ilayer*1000000000;
		return ID;
	} if(VolumeName == "trackerSiLogical") {
		G4long ix = floor(dx / fSiTrackerPixelSize);
		G4long iy = floor(dy / fSiTrackerPixelSize);
		G4long iz = ilayer;
		// sanity check
		if (ix < 0 || ix > 9999 || iy < 0 || iy > 9999 || iz < 0 || iz > 999) {
			G4cerr << "ERROR : getCHannelIDfromXYZ problem ix:" << ix << " iy:" << iy << " iz:" << iz << G4endl;
			return 0;
		}
		G4long ID = ix + iy*10000 + iz*100000000 + CopyNumber*10000000000LL + 100000000000LL;
		return ID;
	} else {
		G4cerr << "Don't know how to handle volume " << VolumeName << G4endl;
		return 0;
	}
}

void DetectorConstruction::CreateMagnetSystem(G4double zLocation, G4LogicalVolume* parent) {

	G4Material *G4_air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");

	G4double innerRadius = 0.0*cm;
	G4double outerRadius = 10.0*cm;
	G4double halfHeight1 = 100.0*cm/2.0;  
	G4double startAngle = 0.0*deg;
	G4double spanningAngle = 360.0*deg;
	G4Tubs* shortsolidCylinder = new G4Tubs("MagnetShort", innerRadius, outerRadius, halfHeight1, startAngle, spanningAngle);
	G4LogicalVolume* shortCylLogic = new G4LogicalVolume(shortsolidCylinder, G4_air, "ShortCylLogical");

	G4double halfHeight2 = 150.0*cm/2.0;  
	G4Tubs* longsolidCylinder = new G4Tubs("MagnetLong", innerRadius, outerRadius, halfHeight2, startAngle, spanningAngle);
	G4LogicalVolume* longCylLogic = new G4LogicalVolume(longsolidCylinder, G4_air, "LongCylLogical");

		// Create a uniform magnetic field (1 Tesla along the Z-axis)
	G4MagneticField* magField = new G4UniformMagField(G4ThreeVector(0., 0.57*tesla, 0.));

	// Create the field manager and assign the magnetic field
	G4FieldManager* fieldManager = new G4FieldManager();
	fieldManager->SetDetectorField(magField);

	// Define the equation of motion and stepper
	G4Mag_UsualEqRhs* equationOfMotion = new G4Mag_UsualEqRhs(magField);
	G4int nvar = 8;
	G4ClassicalRK4* stepper = new G4ClassicalRK4(equationOfMotion, nvar);
	G4ChordFinder* chordFinder = new G4ChordFinder(magField, 1.0e-2*mm, stepper);
	fieldManager->SetChordFinder(chordFinder);

	shortCylLogic->SetFieldManager(fieldManager, true);
	longCylLogic->SetFieldManager(fieldManager, true);

	G4double z = zLocation;
	z += halfHeight2;
	new G4PVPlacement(0, G4ThreeVector(0,0,z), longCylLogic, "Magnet", parent, false, 0, true);
	z += halfHeight2 + halfHeight1;
	new G4PVPlacement(0, G4ThreeVector(0,0,z), shortCylLogic, "Magnet", parent, false, 1, true);
	z += halfHeight1 + halfHeight1;
	new G4PVPlacement(0, G4ThreeVector(0,0,z), shortCylLogic, "Magnet", parent, false, 2, true);
}


void DetectorConstruction::CreateRearCal(G4double zLocation, G4LogicalVolume* parent) {
	G4double sizeX = 121.2*mm;
	G4double sizeY = 121.2*mm;
	G4double sizeZ_Pb = 2*mm;
	G4double sizeZ_PS = 4*mm;
	G4int nlayer = 66;
	G4Material * G4_Pb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
	G4Material* plasticScintillator = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

	G4Box* scintillatorSolid = new G4Box("PSSlab", sizeX / 2, sizeY / 2, sizeZ_PS / 2);
	G4LogicalVolume* scintillatorLogic = new G4LogicalVolume(scintillatorSolid, plasticScintillator, "rearCalscintillatorLogical");
	G4Box* absorberSolid = new G4Box("PbSlab", sizeX / 2, sizeY / 2, sizeZ_Pb / 2);
	G4LogicalVolume* absorberLogic = new G4LogicalVolume(absorberSolid, G4_Pb, "absorberLogical");

	G4double sizeZ_module = nlayer * (sizeZ_Pb + sizeZ_PS);
	G4Box* rearCal_module = new G4Box("rearCalmodule", sizeX / 2, sizeY / 2, sizeZ_module / 2);

	G4LogicalVolume *rearCal_moduleLogic = new G4LogicalVolume(rearCal_module, fWorldMaterial, "rearCalmoduleLogical");

	for (int i = 0; i < nlayer; i++){
		double z = i*(sizeZ_Pb+sizeZ_PS)-sizeZ_module/2.0+sizeZ_Pb/2.0;
		new G4PVPlacement(0, G4ThreeVector(0,0,z), absorberLogic, "rearCalAbs", rearCal_moduleLogic, false, i, true);
		z += sizeZ_Pb/2.0+sizeZ_PS/2.0;
		new G4PVPlacement(0, G4ThreeVector(0,0,z), scintillatorLogic, "rearCalScint", rearCal_moduleLogic, false, i, true);
	}

	// now place 5x5 matrix
	for (int i = 0; i<25; i++) {
		double x = (i%5)*sizeX - 2.5*sizeX + sizeX/2.0 + fFASERCal_LOS_shiftX;
		double y = (i/5)*sizeY - 2.5*sizeY + sizeY/2.0 + fFASERCal_LOS_shiftY;
		new G4PVPlacement(0, G4ThreeVector(x,y,zLocation+sizeZ_module/2.0), rearCal_moduleLogic, "rearCal", parent, false, i, true);
	}
}

void DetectorConstruction::CreateRearHCalMuTag(G4double zLocation, G4LogicalVolume* parent) {
	// dimensions of the rear HCal
	G4double sizeZ_Pb = 9*cm;
	G4double sizeZ_PS = 1*cm;
	G4int nlayer = 9;
	// dimensions of the neutron absorber
	G4double sizeZ_nabs = 10*cm;  /// polyethylene slab
	// dimensions of the rear MuTag
	G4double sizeX = 60*cm;
	G4double sizeY = 60*cm;
	G4double sizeZ = 4*cm;

	G4Material * G4_Pb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
	G4Material* polyethylene = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYETHYLENE");
	G4Material* plasticScintillator = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

	double sizeZ_HCalmodule = sizeZ_Pb + sizeZ_PS;
	double total_sizeZ = nlayer*sizeZ_HCalmodule;

	double Pb_mass = nlayer*sizeX*sizeY*sizeZ_Pb*(G4_Pb->GetDensity()/(kg/cm3));
	G4cout << "Total mass of RearMuCal absorber " << Pb_mass << " kg " << G4endl;

	G4Box* absorberSolid = new G4Box("PbSlab", sizeX / 2, sizeY / 2, sizeZ_Pb / 2);
	G4LogicalVolume* absorberLogic = new G4LogicalVolume(absorberSolid, G4_Pb, "absorberLogical");

	G4Box* nabsorberSolid = new G4Box("NeutronAbsSlab", sizeX / 2, sizeY / 2, sizeZ_nabs / 2);
	G4LogicalVolume* nabsorberLogic = new G4LogicalVolume(nabsorberSolid, polyethylene, "neutabsorberLogical");

	G4Box* HscintillatorSolid = new G4Box("HCALPSSlab", sizeX / 2, sizeY / 2, sizeZ_PS / 2);
	G4LogicalVolume* HscintillatorLogic = new G4LogicalVolume(HscintillatorSolid, plasticScintillator, "rearHCalscintillatorLogical");

	G4Box* scintillatorSolid = new G4Box("PSSlab", sizeX / 2, sizeY / 2, sizeZ / 2);
	G4LogicalVolume* scintillatorLogic = new G4LogicalVolume(scintillatorSolid, plasticScintillator, "muCalscintillatorLogical");

	double x = fFASERCal_LOS_shiftX;
	double y = fFASERCal_LOS_shiftY;
	for(int i = 0; i < nlayer; i++) {
		double z = zLocation + i*sizeZ_HCalmodule + sizeZ_Pb/2.0;
		new G4PVPlacement(0, G4ThreeVector(x,y,z), absorberLogic, "rearHCalAbs", parent, false, i, true);
		z += sizeZ_Pb/2.0 + sizeZ_PS/2.0;
		new G4PVPlacement(0, G4ThreeVector(x,y,z), HscintillatorLogic, "rearHCalScint", parent, false, i, true);
	}
//	new G4PVPlacement(0, G4ThreeVector(0,0,z), absorberLogic, "rearMuCalAbs", parent, false, 0, true);

	double z = zLocation + total_sizeZ + sizeZ_nabs/2.0;
	new G4PVPlacement(0, G4ThreeVector(x,y,z), nabsorberLogic, "rearMuCalNAbs", parent, false, 0, true);

	z = zLocation + total_sizeZ + sizeZ_nabs + sizeZ/2.0;
	new G4PVPlacement(0, G4ThreeVector(x,y,z), scintillatorLogic, "rearMuCalScint", parent, false, 0, true);
}

void DetectorConstruction::CreateFrontTarget(G4double zLocation, G4LogicalVolume *parent) {
	G4double sizeX = 48*cm;
	G4double sizeY = 48*cm;
	G4double sizeZ = 20*cm;
	G4Material * G4_Pb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
	double Pb_mass = sizeX*sizeY*sizeZ*(G4_Pb->GetDensity()/(kg/cm3));
	G4cout << "Total mass of Front target " << Pb_mass << " g " << G4endl;

	G4Box* absorberSolid = new G4Box("PbSlab", sizeX / 2, sizeY / 2, sizeZ / 2);
	G4LogicalVolume* absorberLogic = new G4LogicalVolume(absorberSolid, G4_Pb, "frontTargetLogical");
	double z = zLocation + sizeZ/2.0;
	new G4PVPlacement(0, G4ThreeVector(0,0,z), absorberLogic, "frontTarget", parent, false, 0, true);
}
