#include "DetectorConstruction.hh"
#include "G4ios.hh"
#include "G4GDMLParser.hh"

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


	G4double WorldSizeX = 1.2* sizeX;
	G4double WorldSizeY = 1.2* sizeY;
	G4double WorldSizeZ = 1.2* NRep*(sizetargetWZ + sizeScintillatorZ);

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
					 G4ThreeVector( 0,0,0),  // at (0,0,0)
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

	G4Material * G4_Target = fWorldMaterial; // G4_graphite; // G4_Cu; // fWorldMaterial;

	CreateFaserNu(fPolyvinyltoluene, G4_Target, G4ThreeVector(sizeScintillatorX, sizeScintillatorY, 
	sizeScintillatorZ),G4ThreeVector(sizetargetWX, sizetargetWY, sizetargetWZ),worldLV, NRep);
	fParticleManager->setDetectorInformation(fPolyvinyltoluene->GetName(),
	XYZVector(sizeScintillatorX, sizeScintillatorY, sizeScintillatorZ), G4_Target->GetName(),  
	XYZVector(sizetargetWX, sizetargetWY, sizetargetWZ), NRep);

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

void DetectorConstruction::CreateFaserNu(G4Material* material1, G4Material* material2, G4ThreeVector size1, 
		G4ThreeVector size2, G4LogicalVolume* parent, G4int NRep){

	// add precision tracker of composed of Si
	G4Material * G4_Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
	G4cout << "Thickness of a single Silicon tracker layer" << fSiTrackerSizeZ << " mm " << G4endl;
	G4cout << "Number of silicon tracker layer " << fNumberRep_SiTracker << G4endl;

	G4double sizeX = std::max(size1.getX(), size2.getX());
	G4double sizeY = std::max(size1.getY(), size2.getY());
	G4double sizeZ = size1.getZ() + size2.getZ() + 	fNumberRep_SiTracker*fSiTrackerSizeZ;

	fSandwichLength = sizeZ;

	G4cout << "Total size of one sandwich layer " << sizeZ*mm << " mm " << G4endl;

	G4cout << "Number of layers " << NRep << G4endl;

	fTotalLength = NRep*sizeZ;
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
	userLimits_Scint->SetMaxAllowedStep(1*mm);
	G4UserLimits* userLimits_targetW = new G4UserLimits();
//	userLimits_targetW->SetMaxAllowedStep(1*mm);   // is this really necessary??

	// In this replica solid we place the two actual materials. This replica solid is then replicated NRep times.
    G4Box * replicaSolid = new G4Box("ReplicaSolid" , sizeX/2, sizeY/2, sizeZ/2);

    G4Box* scintillatorSolid = new G4Box("ScintillatorSlab", size1.getX() / 2, size1.getY() / 2, size1.getZ() / 2);
    G4Box* targetWSolid = new G4Box("targetWSlab", size2.getX() / 2, size2.getY() / 2, size2.getZ() / 2);
    G4Box* trackerSiSolid = new G4Box("trackerSiSlab", size2.getX() / 2, size2.getY() / 2, fSiTrackerSizeZ / 2);

	G4LogicalVolume* replicaLogic = new G4LogicalVolume(replicaSolid, fWorldMaterial, "replicaLogical");

	G4LogicalVolume* scintillatorLogic = new G4LogicalVolume(scintillatorSolid, material1, "ScintillatorLogical");
	scintillatorLogic->SetUserLimits(userLimits_Scint);
	G4LogicalVolume* targetWLogic = new G4LogicalVolume(targetWSolid, material2, "targetWLogical");
	targetWLogic->SetUserLimits(userLimits_targetW);
	G4LogicalVolume* trackerSiLogic = new G4LogicalVolume(trackerSiSolid, G4_Si, "trackerSiLogical");


    //We place the scinitllator and targetW into the replica
    new G4PVPlacement(0, G4ThreeVector(0,0,-sizeZ/2 + size1.getZ()/2), scintillatorLogic, "Scintillator", replicaLogic, false, 1, true);
    new G4PVPlacement(0, G4ThreeVector(0,0,sizeZ/2 - size2.getZ()/2 - fSiTrackerSizeZ), targetWLogic, "targetW", replicaLogic, false, 1, true);
    new G4PVPlacement(0, G4ThreeVector(0,0,sizeZ/2 - fSiTrackerSizeZ - size2.getZ() - fSiTrackerSizeZ/2), 
					trackerSiLogic, "trackerSi", replicaLogic, false, 1, true);
    new G4PVPlacement(0, G4ThreeVector(0,0,sizeZ/2 - fSiTrackerSizeZ/2), trackerSiLogic, "trackerSi", replicaLogic, false, 2, true);

    // Create container, which hosts Nrep replicas of our layer. This container then is set inside the world volume
    G4Box* containerSolid = new G4Box("ContainerBox", sizeX/2, sizeY / 2, NRep*sizeZ / 2);

    G4LogicalVolume* containerLogic = new G4LogicalVolume(containerSolid, fWorldMaterial, "ContainerLogical");
    new G4PVReplica("replica", replicaLogic, containerLogic, kZAxis, NRep, sizeZ, 0);
//    new G4PVPlacement(0, G4ThreeVector(0,0,NRep*sizeZ/2), containerLogic, "ContainerPlacement", parent,  false, 0, true);
    new G4PVPlacement(0, G4ThreeVector(0,0,0), containerLogic, "ContainerPlacement", parent,  false, 0, true);
}

G4long DetectorConstruction::getChannelIDfromXYZ(std::string const& VolumeName, int CopyNumber, int MotherCopyNumber, XYZVector const& position) const {

	G4double dx = position.X()+fScintillatorSizeX/2.0;
	G4double dy = position.Y()+fScintillatorSizeY/2.0;
	G4double epsilon = 1e-4;   // avoid rounding errors at volume boundary
	G4double dz = position.Z()+fTotalLength/2.0+epsilon;

	// sanity check
	if((dx < 0 || dx > fScintillatorSizeX) || (dy < 0 || dy > fScintillatorSizeY) || (dz < 0 || dz > fTotalLength)) {
		G4cerr << "ERROR : getCHannelIDfromXYZ problem dx:" << dx << " dy:" << dy << " dz:" << dz << G4endl;
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
		G4long iz = floor((dz-ilayer*fSandwichLength)/ fScintillatorVoxelSize);
		// sanity check
		if (ix < 0 || ix > 999 || iy < 0 || iy > 999 || iz < 0 || iz > 999 || ilayer < 0) {
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
