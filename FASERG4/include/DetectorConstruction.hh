#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "DetectorMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PVReplica.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4UserLimits.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "ParticleManager.hh"
#include "TrackerSD.hh"
#include "globals.hh"
#include "tls.hh"
#include "G4UserLimits.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;

class DetectorMessenger;

/**
 * @class DetectorConstruction
 * @brief Builds the detector geometry
 * @details This class is responsible for constructing the detector geometry
 * and materials. It also sets the user limits for the simulation.
 * Any kind of material is defined here and potential fields are set.
 */
class DetectorConstruction : public G4VUserDetectorConstruction {
    public:
	DetectorConstruction(ParticleManager* photonManager);
	~DetectorConstruction() override;

	G4VPhysicalVolume* Construct() override;
	void ConstructSDandField() override;

	// Set methods for scintillator
	void SetScintillatorMaterial(G4String materialChoice);
	void SetMaxStep(G4double);
	void SetCheckOverlaps(G4bool);
	void SetLightYield(G4double);
	void SetScintillationDecayTime(G4double);
	void SetScintillatorSizeX(G4double);
	void SetScintillatorSizeY(G4double);
	void SetScintillatorSizeZ(G4double);
	void SetVoxelSize(G4double voxel) { fScintillatorVoxelSize = voxel;};

	G4double fScintillatorSizeX = 1000 * mm;  // The size of the detector in mm, default 1000 mm, set by macro
	G4double fScintillatorSizeY = 1000 * mm;  // The size of the detector in mm, default 1000 mm, set by macro
	G4double fScintillatorSizeZ = 1000 * mm;  // The size of the detector in mm, default 1000 mm, set by macro
	G4double fScintillatorVoxelSize = 5 * mm;

	G4double fFASERCal_LOS_shiftX = 0 * cm;
	G4double fFASERCal_LOS_shiftY = 0 * cm;

	double getScintillatorSizeX() const { return fScintillatorSizeX; }
	double getScintillatorSizeY() const { return fScintillatorSizeY; }
	double getScintillatorSizeZ() const { return fScintillatorSizeZ; }

	void SettargetWSizeX(G4double);
	void SettargetWSizeY(G4double);
	void SettargetWSizeZ(G4double);

	G4double ftargetWSizeX = 1000 * mm;  // The size of the detector in mm, default 1000 mm, set by macro
	G4double ftargetWSizeY = 1000 * mm;  // The size of the detector in mm, default 1000 mm, set by macro
	G4double ftargetWSizeZ = 1000 * mm;  // The size of the detector in mm, default 1000 mm, set by macro

	double gettargetWSizeX() const { return ftargetWSizeX; }
	double gettargetWSizeY() const { return ftargetWSizeY; }
	double gettargetWSizeZ() const { return ftargetWSizeZ; }

	G4double fSiTrackerSizeZ = 0.2 * mm;
	G4double fSiTrackerPixelSize = 0.1 * mm;
	G4double fSiTrackerGap = 4 * mm;
	G4int fNumberReplicas = 1;
	G4int fNumberRep_SiTracker = 2;

	G4double fAirGap = 3*cm;
	G4double fAlPlateThickness = 0.5*cm;

	G4double fTotalLength;      // of the full detector
	G4double fSandwichLength;   // of the given sandwich of scint+W+Silicon
	G4double fTotalMass;
	G4double fTotalWMass;
	G4double fTotalScintMass;

	// read ECAL
	G4double fECalSizeX = 0;
	G4double fECalSizeY = 0;
	G4double fRearCalSizeZ;

	// rear HCAL
	G4double fRearHCalSizeX = 720 * mm;
	G4double fRearHCalSizeY = 720 * mm;
	G4double fRearHCalVoxelSize = 40 * mm;
	G4double fRearHCalLength;
	G4double fRearHCal_LOS_shiftX = 0 * cm;
	G4double fRearHCal_LOS_shiftY = 0 * cm;

	// read muon spectrometer
	G4double fRearMuSpectSizeX = 1000 * mm;
	G4double fRearMuSpectSizeY = 1000 * mm;
	G4double fRearMuSpectLocZ; // in mm
	G4double fRearMuSpectSizeZ; // in mm
	G4double fRearMuSpect_LOS_shiftX = 0 * cm;
	G4double fRearMuSpect_LOS_shiftY = 0 * cm;

	void SetNumberReplicas(G4int);
	G4int getNumberReplicas() const {return fNumberReplicas ;}

	// channel ID
	G4long getChannelIDfromXYZ(std::string const& VolumeName, int CopyNumber, int MotherCopyVolume, XYZVector const& position) const;

	G4long getHCalChannelIDfromXYZ(int CopyNumber, XYZVector const& position) const;

    private:
	// methods
	void DefineMaterials();		     ///< Define materials used in the detector
	G4VPhysicalVolume* DefineVolumes();  ///< Define volumes used in the detector

	// static data members
	static G4ThreadLocal G4GlobalMagFieldMessenger* fMagFieldMessenger;  ///< Magnetic field messenger (Currently not used)
	ParticleManager* fParticleManager = nullptr;  ///< Pointer to the particle manager, which takes care of the tracking of all the partilces
	// magnetic field messenger
	// data members

	/// @brief Creates the FASERcal geometry
	/// @param material1 = scintillator geometry
	/// @param material2 = target geometry
	/// @param size1 = thickness scintillator
	/// @param size2 = thickness target
	/// @param parent 
	/// @param Nrep = number of layers
	void CreateFaserCal(G4double zLocation, G4Material* material1, G4Material* material2, G4ThreeVector size1, G4ThreeVector size2, G4LogicalVolume* parent, G4int Nrep);

	/// @brief Create the Magnet system
	void CreateMagnetSystem(G4double zLocation, G4LogicalVolume* parent);

	/// @brief Create the RearCalorimeter
	void CreateRearCal(G4double zLocation, G4LogicalVolume* parent);

	/// @brief Create the rear HCAL and rear muon tagger
	/// @param zLocation 
	/// @param parent 
	void CreateRearHCal(G4double zLocation, G4LogicalVolume* parent);

	/// @brief Create the rear muon spectrometer
	void CreateRearMuSpectrometer(G4double zLocation, G4LogicalVolume* parent);

	/// @brief Create the front Pb target
	/// @param zLocation 
	/// @param parent 
	void CreateFrontTarget(G4double zLocation, G4LogicalVolume *parent);

	G4LogicalVolume* fLogicScintillator = nullptr;	///< Logical volume of the scintillator

	G4Element* fCarbon = nullptr;		  ///< Element Carbon to build PVT
	G4Element* fHydrogen = nullptr;		  ///< Element Hydrogen to build PVT
	G4Material* fPolyvinyltoluene = nullptr;  ///< Material PVT

	G4Material* fScintillatorMaterial = nullptr;  ///< Material of the scintillator (Default PVT, can be set via macro)
	G4Material* fWorldMaterial = nullptr;	      ///< Material of the world (Air)

	// Material property tables are used to define optical properties
	G4MaterialPropertiesTable* fPolyvinyltoluene_MPT = nullptr;
	G4MaterialPropertiesTable* fAir_MPT = nullptr;

	G4UserLimits* fStepLimit = nullptr;  ///< pointer to user step limits

	DetectorMessenger* fMessenger = nullptr;  ///< messenger

	G4bool fCheckOverlaps = true;  ///< option to activate checking of volumes overlaps

	// This section contains the  optical properties of the materials. Long lists of numbers, sorry for that.

	// Air
	const static int nEntriesAir = 2;
	double fPhotonEnergyAir[nEntriesAir] = {0.01 * eV, 100 * eV};
	double fRefractiveIndex_Air[nEntriesAir] = {1.0, 1.0};

	// The scintillator moddeled here is EJ-262 from Eljen Technology
	const static int nEntriesPVT = 301;
	G4double fScintillationDecayTime = 2.1;	 // The decay time of the scintillator in ns
	G4double fLightYield = 8700;		 // The light yield of the scintillator in photons/MeV, can be also set by the macro file
	double fPhotonEnergyPVT[nEntriesPVT] = {
	    2.0664 * eV,  2.06813 * eV, 2.06985 * eV, 2.07158 * eV, 2.07331 * eV, 2.07505 * eV, 2.07679 * eV, 2.07853 * eV, 2.08027 * eV,
	    2.08202 * eV, 2.08377 * eV, 2.08552 * eV, 2.08728 * eV, 2.08903 * eV, 2.0908 * eV,	2.09256 * eV, 2.09433 * eV, 2.0961 * eV,
	    2.09787 * eV, 2.09965 * eV, 2.10143 * eV, 2.10321 * eV, 2.10499 * eV, 2.10678 * eV, 2.10857 * eV, 2.11037 * eV, 2.11217 * eV,
	    2.11397 * eV, 2.11577 * eV, 2.11758 * eV, 2.11939 * eV, 2.1212 * eV,  2.12302 * eV, 2.12484 * eV, 2.12666 * eV, 2.12848 * eV,
	    2.13031 * eV, 2.13214 * eV, 2.13398 * eV, 2.13582 * eV, 2.13766 * eV, 2.1395 * eV,	2.14135 * eV, 2.1432 * eV,  2.14506 * eV,
	    2.14691 * eV, 2.14877 * eV, 2.15064 * eV, 2.1525 * eV,  2.15437 * eV, 2.15625 * eV, 2.15812 * eV, 2.16 * eV,    2.16189 * eV,
	    2.16377 * eV, 2.16566 * eV, 2.16756 * eV, 2.16945 * eV, 2.17135 * eV, 2.17326 * eV, 2.17516 * eV, 2.17707 * eV, 2.17898 * eV,
	    2.1809 * eV,  2.18282 * eV, 2.18474 * eV, 2.18667 * eV, 2.1886 * eV,  2.19053 * eV, 2.19247 * eV, 2.19441 * eV, 2.19635 * eV,
	    2.1983 * eV,  2.20025 * eV, 2.20221 * eV, 2.20416 * eV, 2.20612 * eV, 2.20809 * eV, 2.21006 * eV, 2.21203 * eV, 2.214 * eV,
	    2.21598 * eV, 2.21796 * eV, 2.21995 * eV, 2.22194 * eV, 2.22393 * eV, 2.22593 * eV, 2.22793 * eV, 2.22993 * eV, 2.23194 * eV,
	    2.23395 * eV, 2.23596 * eV, 2.23798 * eV, 2.24 * eV,    2.24203 * eV, 2.24406 * eV, 2.24609 * eV, 2.24813 * eV, 2.25017 * eV,
	    2.25221 * eV, 2.25426 * eV, 2.25631 * eV, 2.25836 * eV, 2.26042 * eV, 2.26249 * eV, 2.26455 * eV, 2.26662 * eV, 2.2687 * eV,
	    2.27077 * eV, 2.27285 * eV, 2.27494 * eV, 2.27703 * eV, 2.27912 * eV, 2.28122 * eV, 2.28332 * eV, 2.28542 * eV, 2.28753 * eV,
	    2.28964 * eV, 2.29176 * eV, 2.29388 * eV, 2.296 * eV,   2.29813 * eV, 2.30026 * eV, 2.3024 * eV,  2.30454 * eV, 2.30668 * eV,
	    2.30883 * eV, 2.31098 * eV, 2.31314 * eV, 2.3153 * eV,  2.31746 * eV, 2.31963 * eV, 2.3218 * eV,  2.32398 * eV, 2.32616 * eV,
	    2.32834 * eV, 2.33053 * eV, 2.33272 * eV, 2.33492 * eV, 2.33712 * eV, 2.33932 * eV, 2.34153 * eV, 2.34375 * eV, 2.34596 * eV,
	    2.34819 * eV, 2.35041 * eV, 2.35264 * eV, 2.35488 * eV, 2.35711 * eV, 2.35936 * eV, 2.3616 * eV,  2.36386 * eV, 2.36611 * eV,
	    2.36837 * eV, 2.37063 * eV, 2.3729 * eV,  2.37518 * eV, 2.37745 * eV, 2.37974 * eV, 2.38202 * eV, 2.38431 * eV, 2.38661 * eV,
	    2.38891 * eV, 2.39121 * eV, 2.39352 * eV, 2.39583 * eV, 2.39815 * eV, 2.40047 * eV, 2.40279 * eV, 2.40513 * eV, 2.40746 * eV,
	    2.4098 * eV,  2.41214 * eV, 2.41449 * eV, 2.41685 * eV, 2.4192 * eV,  2.42157 * eV, 2.42393 * eV, 2.42631 * eV, 2.42868 * eV,
	    2.43106 * eV, 2.43345 * eV, 2.43584 * eV, 2.43823 * eV, 2.44063 * eV, 2.44304 * eV, 2.44545 * eV, 2.44786 * eV, 2.45028 * eV,
	    2.4527 * eV,  2.45513 * eV, 2.45757 * eV, 2.46 * eV,    2.46245 * eV, 2.46489 * eV, 2.46735 * eV, 2.4698 * eV,  2.47227 * eV,
	    2.47473 * eV, 2.47721 * eV, 2.47968 * eV, 2.48217 * eV, 2.48465 * eV, 2.48715 * eV, 2.48964 * eV, 2.49214 * eV, 2.49465 * eV,
	    2.49716 * eV, 2.49968 * eV, 2.5022 * eV,  2.50473 * eV, 2.50726 * eV, 2.5098 * eV,	2.51234 * eV, 2.51489 * eV, 2.51745 * eV,
	    2.52 * eV,	  2.52257 * eV, 2.52514 * eV, 2.52771 * eV, 2.53029 * eV, 2.53287 * eV, 2.53546 * eV, 2.53806 * eV, 2.54066 * eV,
	    2.54327 * eV, 2.54588 * eV, 2.54849 * eV, 2.55112 * eV, 2.55374 * eV, 2.55638 * eV, 2.55901 * eV, 2.56166 * eV, 2.56431 * eV,
	    2.56696 * eV, 2.56962 * eV, 2.57229 * eV, 2.57496 * eV, 2.57763 * eV, 2.58032 * eV, 2.583 * eV,   2.5857 * eV,  2.5884 * eV,
	    2.5911 * eV,  2.59381 * eV, 2.59653 * eV, 2.59925 * eV, 2.60198 * eV, 2.60471 * eV, 2.60745 * eV, 2.61019 * eV, 2.61294 * eV,
	    2.6157 * eV,  2.61846 * eV, 2.62123 * eV, 2.624 * eV,   2.62678 * eV, 2.62957 * eV, 2.63236 * eV, 2.63516 * eV, 2.63796 * eV,
	    2.64077 * eV, 2.64359 * eV, 2.64641 * eV, 2.64924 * eV, 2.65207 * eV, 2.65491 * eV, 2.65775 * eV, 2.66061 * eV, 2.66346 * eV,
	    2.66633 * eV, 2.6692 * eV,	2.67207 * eV, 2.67496 * eV, 2.67784 * eV, 2.68074 * eV, 2.68364 * eV, 2.68655 * eV, 2.68946 * eV,
	    2.69238 * eV, 2.69531 * eV, 2.69824 * eV, 2.70118 * eV, 2.70413 * eV, 2.70708 * eV, 2.71004 * eV, 2.713 * eV,   2.71597 * eV,
	    2.71895 * eV, 2.72194 * eV, 2.72493 * eV, 2.72793 * eV, 2.73093 * eV, 2.73394 * eV, 2.73696 * eV, 2.73998 * eV, 2.74301 * eV,
	    2.74605 * eV, 2.7491 * eV,	2.75215 * eV, 2.7552 * eV};

	double fScintillation_PVT[nEntriesPVT] = {
	    0.000132212, 0.000137236, 0.00014747,  0.000157705, 0.000162895, 0.000166889, 0.000170883, 0.000174877, 0.000178871, 0.000182865,
	    0.000186859, 0.000190853, 0.000194847, 0.000198842, 0.000202836, 0.00020683,  0.000210824, 0.000215131, 0.000224228, 0.000233326,
	    0.000242424, 0.000251521, 0.000260619, 0.000269726, 0.000279173, 0.000288621, 0.000298068, 0.000307516, 0.000316963, 0.000326411,
	    0.000335858, 0.000345306, 0.000353371, 0.000358333, 0.000363296, 0.000368258, 0.00037322,  0.000378183, 0.000383145, 0.000388107,
	    0.00039307,	 0.000398032, 0.000402994, 0.000410071, 0.000423718, 0.000437364, 0.00045101,  0.000464657, 0.000478303, 0.00049195,
	    0.000505596, 0.00051886,  0.000530557, 0.000542254, 0.000553951, 0.000565648, 0.000577345, 0.000589042, 0.000600739, 0.000612436,
	    0.000624132, 0.00063618,  0.00064831,  0.00066044,	0.00067257,  0.0006847,	  0.00069683,  0.000708961, 0.000721091, 0.000733221,
	    0.000746496, 0.000760143, 0.000773789, 0.000787435, 0.000801082, 0.000814728, 0.000828375, 0.000842021, 0.000855668, 0.000869314,
	    0.000888661, 0.000910495, 0.00093233,  0.000954164, 0.000975998, 0.000997833, 0.00101967,  0.0010415,   0.00106334,	 0.00108517,
	    0.00110914,	 0.00113434,  0.00115953,  0.00118473,	0.00120992,  0.00123511,  0.00126031,  0.0012855,   0.00131069,	 0.00135379,
	    0.00139845,	 0.00144311,  0.00148777,  0.00153243,	0.00157709,  0.00162176,  0.00167329,  0.0017306,   0.00178792,	 0.00184523,
	    0.00190255,	 0.00195986,  0.00201718,  0.0020665,	0.00211426,  0.00216202,  0.00220979,  0.00225755,  0.00230531,	 0.00235307,
	    0.00240084,	 0.00246414,  0.00253157,  0.002599,	0.00266642,  0.00273385,  0.00280085,  0.00286346,  0.00292608,	 0.00298869,
	    0.0030513,	 0.00311392,  0.00316808,  0.00320902,	0.00324995,  0.00329089,  0.00333183,  0.00337277,  0.00341627,	 0.00346444,
	    0.0035126,	 0.00356077,  0.00360893,  0.00365709,	0.00368783,  0.00370938,  0.00373093,  0.00375247,  0.00377402,	 0.00379557,
	    0.00381834,	 0.0038442,   0.00387005,  0.00389591,	0.00392177,  0.00394762,  0.00397348,  0.00400896,  0.00404535,	 0.00408174,
	    0.00411813,	 0.00415452,  0.00419091,  0.00424727,	0.0043064,   0.00436554,  0.00442467,  0.00448381,  0.00454294,	 0.00462617,
	    0.0047139,	 0.00480163,  0.00488936,  0.00497902,	0.00507979,  0.00518057,  0.00528134,  0.00538212,  0.0055003,	 0.00562312,
	    0.00574593,	 0.00589558,  0.0060798,   0.00626403,	0.00646027,  0.00666041,  0.00686056,  0.00705035,  0.00723644,	 0.00742253,
	    0.00760998,	 0.00782412,  0.00803826,  0.00825241,	0.00846655,  0.00868359,  0.00890194,  0.00912028,  0.00933862,	 0.00952061,
	    0.00968437,	 0.00984812,  0.0100119,   0.0101834,	0.0103723,   0.0105613,	  0.0107502,   0.0109392,   0.0111215,	 0.0113034,
	    0.0114854,	 0.0115869,   0.0116824,   0.0116827,	0.0116745,   0.0116663,	  0.0116421,   0.0116057,   0.0115694,	 0.0114687,
	    0.0113152,	 0.0111617,   0.0110251,   0.0108911,	0.0107571,   0.0106178,	  0.0104541,   0.0102903,   0.0101266,	 0.00994031,
	    0.00974926,	 0.00955821,  0.00930423,  0.0090352,	0.0087581,   0.00847152,  0.00819106,  0.00791813,  0.0076452,	 0.00737227,
	    0.00708113,	 0.00676531,  0.00645732,  0.00624239,	0.00602746,  0.00577343,  0.00545956,  0.0051797,   0.00495454,	 0.00472937,
	    0.00445705,	 0.00416463,  0.00389601,  0.00367085,	0.00344568,  0.00321462,  0.00298068,  0.00274885,  0.00252141,	 0.00229397,
	    0.0020815,	 0.00190286,  0.00172422,  0.00154557,	0.00140341,  0.00128374,  0.00116407,  0.0010444,   0.00092875,	 0.000824541,
	    0.000720332, 0.000616122, 0.000520421, 0.000431099, 0.000341777, 0.000252455, 0.000203456, 0.000159788, 0.000116119, 7.24502e-05,
	    2.87816e-05, 2.23978e-05, 2.23978e-05, 2.23978e-05, 2.23978e-05, 2.23978e-05, 2.23978e-05, 2.23978e-05, 2.23978e-05, 2.23978e-05,
	    2.23978e-05, 2.23978e-05, 2.23978e-05, 2.23978e-05, 2.23978e-05, 2.23978e-05, 2.23978e-05, 2.23978e-05, 2.23978e-05, 2.23978e-05,
	    2.23978e-05};

	double fAbsorption_PVT[nEntriesPVT] = {
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm, 250 * cm,
	    250 * cm, 250 * cm};
	double fRefractiveIndex_PVT[nEntriesPVT] = {
	    1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62,
	    1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62,
	    1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62,
	    1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62,
	    1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62,
	    1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62,
	    1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62,
	    1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62,
	    1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62,
	    1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62,
	    1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62,
	    1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62,
	    1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62,
	    1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62};	// 1.62
    public:
	double getLightYield() const { return fLightYield; }
	double getScintillationDecayTime() const { return fScintillationDecayTime; }
	std::vector<double> getRefractiveIndex() const { return std::vector<double>(std::begin(fRefractiveIndex_PVT), std::end(fRefractiveIndex_PVT)); }
	std::vector<double> getAbsorptionLength() const { return std::vector<double>(std::begin(fAbsorption_PVT), std::end(fAbsorption_PVT)); }
	std::vector<double> getEmissionSpectrum() const { return std::vector<double>(std::begin(fScintillation_PVT), std::end(fScintillation_PVT)); }
	std::vector<double> getPhotonEnergy() const { return std::vector<double>(std::begin(fPhotonEnergyPVT), std::end(fPhotonEnergyPVT)); }


};

#endif
