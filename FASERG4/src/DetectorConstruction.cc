#include "DetectorConstruction.hh"
#include "G4ios.hh"
#include "G4GDMLParser.hh"
#include <cstdio> // added by UMUT

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolumeStore.hh"


#include "MuonDetMagneticField.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

DetectorConstruction::DetectorConstruction(ParticleManager* photonManager)
{
	fMessenger = new DetectorMessenger(this);

	fParticleManager = photonManager;
	// UMUT: initialize tilt
    fTiltAngleY = 0.*deg;  // <<< tiltY
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

	// FASER calorimeter only flag
	fonlyFaserCal = false;    // true if only the FASER calorimeter is constructed

	G4int NRep = getNumberReplicas();
	G4double sizetargetWX = gettargetWSizeX();
	G4double sizetargetWY = gettargetWSizeY();
	G4double sizetargetWZ = gettargetWSizeZ();
	G4double sizeScintillatorX = getScintillatorSizeX();
	G4double sizeScintillatorY = getScintillatorSizeY();
	G4double sizeScintillatorZ = getScintillatorSizeZ();

	G4double sizeX = fRearMuSpectSizeX;
	G4double sizeY = fRearMuSpectSizeY;

	//G4double WorldSizeX = 2.5* sizeX;
	//G4double WorldSizeY = 2.5* sizeY;
	//G4double WorldSizeZ = 9*m; // 1.2* NRep*(sizetargetWZ + sizeScintillatorZ);
	// UMUT: tilt the assembly by 5° around Y, the corners move a bit in X and Z.
	G4double WorldSizeX = 5* sizeX;
	G4double WorldSizeY = 5* sizeY;
	G4double WorldSizeZ = 12.*m; // 1.2* NRep*(sizetargetWZ + sizeScintillatorZ);
	
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
	
	// -------------------------------------------------
	// UMUT: global detector rotation
    G4RotationMatrix* detRot = nullptr;                      // <<< tiltY
    if (fTiltAngleY != 0.) {                                 // <<< tiltY
        detRot = new G4RotationMatrix();                     // <<< tiltY
        detRot->rotateY(fTiltAngleY);                        // <<< tiltY
    }   	//	
	// DETECTOR ASSEMBLY (the whole detector lives inside this LV)
    //
    auto detSolid = new G4Box("DetectorAssemblySolid",
                              WorldSizeX/2, WorldSizeY/2, WorldSizeZ/2);
    auto detLV = new G4LogicalVolume(detSolid, fWorldMaterial, "DetectorAssemblyLV");
	// ------------------------------------------------
	// target is composed of W or Copper
	G4Material * G4_W = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
	G4Material * G4_Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
	G4Material * G4_graphite = G4NistManager::Instance()->FindOrBuildMaterial("G4_GRAPHITE");
	G4Material * G4_Pb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");

	G4Material * G4_Target = G4_W;

	G4double zLocation = 0*cm;
	//CreateFaserCal(zLocation, fPolyvinyltoluene, G4_Target, G4ThreeVector(sizeScintillatorX, sizeScintillatorY, 
	//sizeScintillatorZ),G4ThreeVector(sizetargetWX, sizetargetWY, sizetargetWZ),worldLV, NRep);
	CreateFaserCal(zLocation, fPolyvinyltoluene, G4_Target, G4ThreeVector(sizeScintillatorX, sizeScintillatorY, 
	sizeScintillatorZ),G4ThreeVector(sizetargetWX, sizetargetWY, sizetargetWZ),detLV, NRep);
	fParticleManager->setDetectorInformation(fPolyvinyltoluene->GetName(),
	XYZVector(sizeScintillatorX, sizeScintillatorY, sizeScintillatorZ), G4_Target->GetName(),  
	XYZVector(sizetargetWX, sizetargetWY, sizetargetWZ), NRep);

//	G4double sizeZ = (sizeScintillatorZ + sizetargetWZ + fNumberRep_SiTracker*fSiTrackerSizeZ)*NRep;
	G4double sizeZ = fTotalLength;

//	CreateFrontTarget(-sizeZ/2.0-20.0*cm, worldLV);

	zLocation += sizeZ/2.0;

#if magnet
	// CreateMagnetSystem(zLocation, worldLV);
	CreateMagnetSystem(zLocation, detLV);
	zLocation += 350.0*cm;    // length of the magnet system
#endif

	if (!fonlyFaserCal)
	{
		//CreateRearCal(zLocation, worldLV);
		CreateRearCal(zLocation, detLV);

		G4double locZHcal = zLocation + fRearCalSizeZ;
		fRearHCal_LOS_shiftX = fFASERCal_LOS_shiftX + (fRearHCalSizeX - fECalSizeX) / 2.0;
		fRearHCal_LOS_shiftY = fFASERCal_LOS_shiftY + (fRearHCalSizeY - fECalSizeY) / 2.0;
		std::cout << "HCal shift x " << fRearHCal_LOS_shiftX / cm << " cm,  y " << fRearHCal_LOS_shiftY / cm << " cm" << std::endl;
		//CreateRearHCal(locZHcal, worldLV);
		CreateRearHCal(locZHcal, detLV);

		G4double locMuSpect = locZHcal + fRearHCalLength;
		fRearMuSpect_LOS_shiftX = fFASERCal_LOS_shiftX + (fRearMuSpectSizeX - fECalSizeX) / 2.0;
		fRearMuSpect_LOS_shiftY = fFASERCal_LOS_shiftY + (fRearMuSpectSizeY - fECalSizeY) / 2.0;		//CreateRearMuSpectrometer(locMuSpect, worldLV);
		//CreateRearMuSpectrometer(locMuSpect, worldLV);
		CreateRearMuSpectrometer(locMuSpect, detLV);
	}

	// UMUT: Place the detector assembly into the world with the global rotation
	new G4PVPlacement(detRot,                      // <<< tiltY
	                  G4ThreeVector(0,0,0), // given in cm
	                  detLV,
	                  "DetectorAssemblyPV",
	                  worldLV,
	                  false,
	                  0,
	                  fCheckOverlaps);

	// Save the geometry of the detector
	G4GDMLParser parser;
	//UMUT: If a previous geometry file exists, remove it so G4GDML doesn't abort
	// (parser.Write throws if the file already exists)
	std::remove("geometry_tilted_5degree.gdml");
	parser.Write("geometry_tilted_5degree.gdml", worldPV->GetLogicalVolume());

	// Print the total mass of the detector
	G4cout << "----------------------------------" << G4endl;
	//G4cout << "Total mass of the detector : " << worldLV->GetMass() / kg << " kg" << G4endl;
	G4cout << "Total mass of the detector : " << detLV->GetMass() / kg << " kg" << G4endl;
	G4cout << "----------------------------------" << G4endl;
	// loop over all logical volumes hanging from worldLV and print their masses
	G4cout << "Masses of logical volumes : " << G4endl;
	// Collect all logical volumes in the geometry
	std::vector<G4LogicalVolume*> fLogicalVolumes;
	const auto& lvStore = *G4LogicalVolumeStore::GetInstance();
	for (const auto& lv : lvStore) {
		fLogicalVolumes.push_back(lv);
	}
	for (const auto& volume : fLogicalVolumes) {
			G4cout << volume->GetName() << " : " << volume->GetMass() / kg << " kg" << G4endl;
	}
	G4cout << "----------------------------------" << G4endl;

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
	if(!fonlyFaserCal){
		SetSensitiveDetector("rearCalscintillatorLogical", aTrackerSD, false);
		SetSensitiveDetector("rearHCalscintillatorLogical", aTrackerSD, false);
		//	SetSensitiveDetector("muCalscintillatorLogical", aTrackerSD, false);
		SetSensitiveDetector("SciFiLayerLV", aTrackerSD, false);
	}

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

// UMUT: tilt setter
void DetectorConstruction::SetTiltAngleY(G4double angle)         // <<< tiltY
{
    fTiltAngleY = angle;                                        // <<< tiltY
}

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
	#if 1
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

	// position is given in the local coordinate system of the volume
	G4double dx = position.X()+fScintillatorSizeX/2.0-fFASERCal_LOS_shiftX ;
	G4double dy = position.Y()+fScintillatorSizeY/2.0-fFASERCal_LOS_shiftY ;
	G4double epsilon = 1e-4;   // avoid rounding errors at volume boundary
	G4double dz = position.Z()+fTotalLength/2.0-epsilon;
	//G4double dz = position.Z() + fTotalLength/2.0;
	// sanity check
	if((dx < 0 || dx > fScintillatorSizeX) || (dy < 0 || dy > fScintillatorSizeY) || (dz < 0 || dz > fTotalLength)) {
		getchannelIDerrorcount++;
		if(getchannelIDerrorcount < 100) {
			G4cerr << "ERROR : getCHannelIDfromXYZ problem dx:" << dx << " dy:" << dy << " dz:" << dz << G4endl;
			G4cerr << "VolumeName: " << VolumeName << " CopyNumber: " << CopyNumber << " MotherCopyNumber: " << MotherCopyNumber << G4endl;
			G4cerr << "Position: " << position.X() << ", " << position.Y() << ", " << position.Z() << G4endl;
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

		G4long ix = floor((dx / fScintillatorVoxelSize) - epsilon);
		G4long iy = floor((dy / fScintillatorVoxelSize) - epsilon);
		G4long iz = floor((dz-ilayer*fSandwichLength-fAlPlateThickness-ftargetWSizeZ)/ fScintillatorVoxelSize);
		//const G4double eps = 1e-9;
		//G4double zLocal = dz - ilayer*fSandwichLength - fAlPlateThickness - ftargetWSizeZ;
		//G4long iz = (G4long) std::floor((zLocal / fScintillatorVoxelSize) - eps);
		//G4long Nz = (G4long) std::floor((fScintillatorSizeZ / fScintillatorVoxelSize) + eps);
		//if (iz < 0) iz = 0;
		//else if (iz >= Nz) iz = Nz - 1;
		
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

G4long DetectorConstruction::getHCalChannelIDfromXYZ(int CopyNumber, XYZVector const& position) const {
	
	// position is given in the local coordinate system of the volume
	G4double dx = position.X()+fRearHCalSizeX/2.0;
	G4double dy = position.Y()+fRearHCalSizeY/2.0;
	G4double epsilon = 1e-4;   // avoid rounding errors at volume boundary

	// sanity check
	if((dx < 0 || dx > fRearHCalSizeX) || (dy < 0 || dy > fRearHCalSizeY)) {
		getchannelIDerrorcount++;
		if(getchannelIDerrorcount < 100) {
			G4cerr << "ERROR : getHCalCHannelIDfromXYZ problem dx:" << dx << " dy:" << dy << G4endl;
		}
		return 0;
	}

	G4long iz = CopyNumber;
	G4long ix = floor(dx / fRearHCalVoxelSize);
	G4long iy = floor(dy / fRearHCalVoxelSize);
	// sanity check
	if (ix < 0 || ix > 999 || iy < 0 || iy > 999 || iz < 0 || iz > 999) {
		if(iz>1000)
			G4cerr << "ERROR : getHCalCHannelIDfromXYZ problem ix:" << ix << " iy:" << iy << " iz:" << iz << G4endl;
		return 0;
	}
	G4long ID = ix + iy*1000L + iz*1000000LL;
	return ID;
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
	// dimensions of the rear calorimeter
	int nMat = 5;
	fECalSizeX = sizeX*nMat;
	fECalSizeY = sizeY*nMat;
	G4double sizeZ_Pb = 2*mm;
	G4double sizeZ_PS = 4*mm;
	G4int nlayer = 66;
	fRearCalSizeZ = nlayer*(sizeZ_Pb+sizeZ_PS);
	G4cout << "Total length of the rear calorimeter " << fRearCalSizeZ << " mm " << G4endl;

	// Pb and plastic scintillator
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

void DetectorConstruction::CreateRearHCal(G4double zLocation, G4LogicalVolume* parent) {
	// dimensions of the rear HCal
	G4double sizeZ_Fe = 2*cm;
	G4double sizeZ_PS = 0.3*cm;
	G4int nlayer = 40;
	// dimensions of the neutron absorber
	G4double sizeZ_nabs = 10*cm;  /// polyethylene slab
	// dimensions of the rear HCal
	G4double sizeX = 72*cm;
	G4double sizeY = 72*cm;
	G4double sizeZ = 4*cm;

	G4Material * G4_Fe = G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");
	G4Material* polyethylene = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYETHYLENE");
	G4Material* plasticScintillator = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

	double sizeZ_HCalmodule = sizeZ_Fe + sizeZ_PS;
	double total_sizeZ = nlayer*sizeZ_HCalmodule;

	fRearHCalLength = total_sizeZ;

	double Fe_mass = nlayer*sizeX*sizeY*sizeZ_Fe*(G4_Fe->GetDensity()/(kg/cm3));
	G4cout << "Total mass of RearMuCal absorber " << Fe_mass << " kg " << G4endl;

	G4Box* absorberSolid = new G4Box("PbSlab", sizeX / 2, sizeY / 2, sizeZ_Fe / 2);
	G4LogicalVolume* absorberLogic = new G4LogicalVolume(absorberSolid, G4_Fe, "absorberLogical");

	G4Box* nabsorberSolid = new G4Box("NeutronAbsSlab", sizeX / 2, sizeY / 2, sizeZ_nabs / 2);
	G4LogicalVolume* nabsorberLogic = new G4LogicalVolume(nabsorberSolid, polyethylene, "neutabsorberLogical");

	G4Box* HscintillatorSolid = new G4Box("HCALPSSlab", sizeX / 2, sizeY / 2, sizeZ_PS / 2);
	G4LogicalVolume* HscintillatorLogic = new G4LogicalVolume(HscintillatorSolid, plasticScintillator, "rearHCalscintillatorLogical");

//	G4Box* scintillatorSolid = new G4Box("PSSlab", sizeX / 2, sizeY / 2, sizeZ / 2);
//	G4LogicalVolume* scintillatorLogic = new G4LogicalVolume(scintillatorSolid, plasticScintillator, "muCalscintillatorLogical");

    // Create container, which hosts Nrep replicas of our layer. This container then is set inside the world volume
    G4Box* HcalSolid = new G4Box("ContainerHcal", sizeX/2, sizeY / 2, total_sizeZ / 2);
    G4LogicalVolume* HCalcontainerLogic = new G4LogicalVolume(HcalSolid, fWorldMaterial, "HCalContainerLogical");

	for(int i = 0; i < nlayer; i++) {
		double z = i*sizeZ_HCalmodule + sizeZ_Fe/2.0 - total_sizeZ/2.0;
		new G4PVPlacement(0, G4ThreeVector(0,0,z), absorberLogic, "rearHCalAbs", HCalcontainerLogic, false, i, true);
		z += sizeZ_Fe/2.0 + sizeZ_PS/2.0;
		new G4PVPlacement(0, G4ThreeVector(0,0,z), HscintillatorLogic, "rearHCalScint", HCalcontainerLogic, false, i, true);
	}
	double x = fRearHCal_LOS_shiftX;
	double y = fRearHCal_LOS_shiftY;
	new G4PVPlacement(0, G4ThreeVector(x,y,zLocation + total_sizeZ/2.0), HCalcontainerLogic, "rearHCal", parent, false, 0, true);

	#if 0
	double z = zLocation + total_sizeZ + sizeZ_nabs/2.0;
	new G4PVPlacement(0, G4ThreeVector(x,y,z), nabsorberLogic, "rearMuCalNAbs", parent, false, 0, true);

	z = zLocation + total_sizeZ + sizeZ_nabs + sizeZ/2.0;
	new G4PVPlacement(0, G4ThreeVector(x,y,z), scintillatorLogic, "rearMuCalScint", parent, false, 0, true);
	#endif
}

void DetectorConstruction::CreateRearMuSpectrometer(G4double zLocation, G4LogicalVolume* parent) {
	auto nist = G4NistManager::Instance();
	G4Material *air = nist->FindOrBuildMaterial("G4_AIR");
	G4Material *steel = nist->FindOrBuildMaterial("G4_Fe");
	G4Material *plastic = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
	G4Material *aluminum = nist->FindOrBuildMaterial("G4_Al");

	// Parameters
	const int nMagnets = 10;
	const int nSciFiGroups = nMagnets + 1;
	const int nSciFiPerGroup = 4;
	const int nSciFiPlanes = nSciFiGroups * nSciFiPerGroup;
	G4double magnetSizeX = fRearMuSpectSizeX;
	G4double magnetSizeY = fRearMuSpectSizeY;
	G4double scifilayerThickness = 2.5 * mm;
	G4double magnetThickness = 150. * mm; // Thickness of each magnet changed from 100 to 150
	G4double gapBeforeMagnet = 10. * mm;
	G4double gapAfterMagnet = 10. * mm;

	fRearMuSpectLocZ = zLocation;

	// Calculate total length needed
	G4double totalSciFiThickness = nSciFiPlanes * scifilayerThickness;
	G4double totalMagnetThickness = nMagnets * magnetThickness;
	G4double totalGapBefore = nMagnets * gapBeforeMagnet;
	G4double totalGapAfter = nMagnets * gapAfterMagnet;
	G4double totalLength = totalSciFiThickness + totalMagnetThickness + totalGapBefore + totalGapAfter;
	fRearMuSpectSizeZ = totalLength;
	G4cout << "Total length of the muon spectrometer " << totalLength / mm << " mm " << G4endl;

	G4double zStart = -totalLength / 2;

	double x = fRearMuSpect_LOS_shiftX;
	double y = fRearMuSpect_LOS_shiftY;
	double z = zLocation + totalLength / 2;
	G4Box* muonSpectrometerBox = new G4Box("MuonSpectrometer", magnetSizeX / 2, magnetSizeY / 2, totalLength / 2);
	G4LogicalVolume* muonSpectrometerLV = new G4LogicalVolume(muonSpectrometerBox, air, "MuonSpectrometerLV");
	new G4PVPlacement(0, G4ThreeVector(x,y,z), muonSpectrometerLV, "MuonSpectrometer", parent, false, 0, true);

	int scifiLayerID = 0;
	G4double zPos = zStart + scifilayerThickness / 2;

	G4Box *scifiBox = new G4Box("SciFiLayer", magnetSizeX / 2, magnetSizeY / 2, scifilayerThickness / 2);
	G4LogicalVolume *scifiLV = new G4LogicalVolume(scifiBox, plastic, "SciFiLayerLV");

			for (int group = 0; group < nSciFiGroups; ++group)
	{
		// Place 4 SciFi layers
		for (int l = 0; l < nSciFiPerGroup; ++l)
		{
			G4String physName = "SciFiLayer_" + std::to_string(scifiLayerID);
			new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zPos), scifiLV, physName, muonSpectrometerLV, false, scifiLayerID);

			auto scifiVis = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); // green
			scifiVis->SetForceSolid(true);
			scifiLV->SetVisAttributes(scifiVis);
			// scifiLV->SetSensitiveDetector(scifiSD);

			zPos += scifilayerThickness;
			scifiLayerID++;
		}

		// After last group, don't place a magnet
		if (group == nMagnets)
			break;

		// Gap before magnet
		zPos += gapBeforeMagnet;

		// Place magnet
		G4Box *solidMagnet = new G4Box("Magnet", magnetSizeX / 2, magnetSizeY / 2, magnetThickness / 2);
		// Create slits in the magnet
		G4double slitWidth = magnetSizeX / 2.0;
		G4double slitHeight = 20. * mm;
		G4double slitPosition = 250. * mm;
		// Create two slits symmetrically positioned above and below the center
		G4Box *slitBox = new G4Box("Slit", slitWidth / 2.0, slitHeight / 2.0, magnetThickness / 2);

		G4SubtractionSolid *magnetMinusSlit1 = new G4SubtractionSolid(
			"MagnetMinusSlit1", solidMagnet, slitBox, nullptr, G4ThreeVector(0, slitPosition, 0));
		G4SubtractionSolid *magnetWithSlits = new G4SubtractionSolid(
			"MagnetWithSlits", magnetMinusSlit1, slitBox, nullptr, G4ThreeVector(0, -slitPosition, 0));

		G4LogicalVolume *magnetLV = new G4LogicalVolume(magnetWithSlits, steel, "MagnetLV");
		new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zPos + magnetThickness / 2 - scifilayerThickness / 2), magnetLV, "Magnet", muonSpectrometerLV, false, group);

		// Magnetic field on the magnet volume
		auto field = new MuonMagneticField();
		field->SetSlitPosition(slitPosition);
		auto fieldManager = new G4FieldManager(field);
		fieldManager->CreateChordFinder(field);
		magnetLV->SetFieldManager(fieldManager, true);
		
		////////////// added umut
		// If particle manager verbosity enabled, probe the field at magnet center and slit positions
		if (fParticleManager && fParticleManager->verbose()) {
			G4ThreeVector centerGlobal(0, 0, zPos + magnetThickness / 2 - scifilayerThickness / 2);
			G4double slitPos = slitPosition;
			G4ThreeVector topSlitGlobal(0, slitPos, zPos + magnetThickness / 2 - scifilayerThickness / 2);
			G4ThreeVector botSlitGlobal(0, -slitPos, zPos + magnetThickness / 2 - scifilayerThickness / 2);
			// Query the field attached to this logical volume via its field manager
			const G4FieldManager* fm = magnetLV->GetFieldManager();
			if (fm) {
				const G4MagneticField* mf = dynamic_cast<const G4MagneticField*>(fm->GetDetectorField());
				if (mf) {
					G4double p[4]; G4double B[3];
					p[0]=centerGlobal.x(); p[1]=centerGlobal.y(); p[2]=centerGlobal.z(); p[3]=0.;
					mf->GetFieldValue(p, B);
					G4cout << "MAGNET_FIELD_TEST: magnet=" << group << " center B[T]=" << B[0]/tesla << "," << B[1]/tesla << "," << B[2]/tesla << " pos[mm]=" << centerGlobal.x()/mm << "," << centerGlobal.y()/mm << "," << centerGlobal.z()/mm << G4endl;
					p[0]=topSlitGlobal.x(); p[1]=topSlitGlobal.y(); p[2]=topSlitGlobal.z();
					mf->GetFieldValue(p, B);
					G4cout << "MAGNET_FIELD_TEST: magnet=" << group << " topslit B[T]=" << B[0]/tesla << "," << B[1]/tesla << "," << B[2]/tesla << " pos[mm]=" << topSlitGlobal.x()/mm << "," << topSlitGlobal.y()/mm << "," << topSlitGlobal.z()/mm << G4endl;
					p[0]=botSlitGlobal.x(); p[1]=botSlitGlobal.y(); p[2]=botSlitGlobal.z();
					mf->GetFieldValue(p, B);
					G4cout << "MAGNET_FIELD_TEST: magnet=" << group << " botslit B[T]=" << B[0]/tesla << "," << B[1]/tesla << "," << B[2]/tesla << " pos[mm]=" << botSlitGlobal.x()/mm << "," << botSlitGlobal.y()/mm << "," << botSlitGlobal.z()/mm << G4endl;
				}
			}
		}
		// Magnet (gray)
		auto magnetVis = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)); // RGB: gray
		magnetVis->SetForceSolid(true);
		magnetLV->SetVisAttributes(magnetVis);

		#if 0
		// Place aluminum strips as before, using the magnet center z
		G4double magnetZ = zPos + magnetThickness / 2 - scifilayerThickness / 2;

		int nTopStrips = 11;
		G4double stripWidth = 25. * mm;
		G4double stripHalfLength = 125. * mm;
		G4double stripThickness = 2. * mm;

		for (int j = 0; j < nTopStrips; ++j)
		{
			G4bool placeFront = (j % 2 == 0);
			G4double xOffset = -250. * mm + j * 50. * mm;
			G4double yOffset = 375. * mm;
			G4Box *aluStripSolid = new G4Box("AluStrip", stripWidth, stripHalfLength, stripThickness);
			G4LogicalVolume *aluStripLV = new G4LogicalVolume(aluStripSolid, aluminum, "AluStripLV");
			auto aluVis = new G4VisAttributes(G4Colour(0.8, 0.6, 0.4));
			aluVis->SetForceSolid(true);
			aluStripLV->SetVisAttributes(aluVis);
			G4double zPlacement = magnetZ + (placeFront ? -magnetThickness / 2 - 2. * mm : magnetThickness / 2 + 2. * mm);
			G4ThreeVector pos(xOffset, yOffset, zPlacement);
			new G4PVPlacement(nullptr, pos, aluStripLV, "AluStrip", muonSpectrometerLV, false, group * 1000 + j);
		}
		for (int j = 0; j < nTopStrips; ++j)
		{
			G4bool placeFront = (j % 2 == 0);
			G4double xOffset = -250. * mm + j * 50. * mm;
			G4double yOffset = -375. * mm;
			G4Box *aluStripSolid = new G4Box("AluStrip", stripWidth, stripHalfLength, stripThickness);
			G4LogicalVolume *aluStripLV = new G4LogicalVolume(aluStripSolid, aluminum, "AluStripLV");
			auto aluVis = new G4VisAttributes(G4Colour(0.8, 0.6, 0.4));
			aluVis->SetForceSolid(true);
			aluStripLV->SetVisAttributes(aluVis);
			G4double zPlacement = magnetZ + (placeFront ? -magnetThickness / 2 - 2. * mm : magnetThickness / 2 + 2. * mm);
			G4ThreeVector pos(xOffset, yOffset, zPlacement);
			new G4PVPlacement(nullptr, pos, aluStripLV, "AluStrip", muonSpectrometerLV, false, group * 1000 + j + 100);
		}
		int nMiddleStrips = 12;
		G4double MiddlestripHalfLength = 250. * mm;
		for (int j = 1; j < nMiddleStrips; ++j)
		{
			G4bool placeFront = (j % 2 == 0);
			G4double xOffset = -300. * mm + j * 50. * mm;
			G4double yOffset = 0;
			G4Box *aluStripSolid = new G4Box("AluStrip", stripWidth, MiddlestripHalfLength, stripThickness);
			G4LogicalVolume *aluStripLV = new G4LogicalVolume(aluStripSolid, aluminum, "AluStripLV");
			auto aluVis = new G4VisAttributes(G4Colour(0.8, 0.6, 0.4));
			aluVis->SetForceSolid(true);
			aluStripLV->SetVisAttributes(aluVis);
			G4double zPlacement = magnetZ + (placeFront ? -magnetThickness / 2 - 2. * mm : magnetThickness / 2 + 2. * mm);
			G4ThreeVector pos(xOffset, yOffset, zPlacement);
			new G4PVPlacement(nullptr, pos, aluStripLV, "AluStrip", muonSpectrometerLV, false, group * 1000 + j + 200);
		}
#endif
		// Advance zPos past the magnet
		zPos += magnetThickness;

		// Gap after magnet
		zPos += gapAfterMagnet;
	}
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
