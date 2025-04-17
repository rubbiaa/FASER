#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIdirectory.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;


/**
 * @class DetectorMessenger
 * @brief Messenger for the detector construction
 *
 * This class is used to define the commands that can be used to modify the
 * detector geometry.
 * It implements the following commands:
 * - /PLATON/scint/material Sets the material of the plastic scintillator(PVT as only option currently)
 * - /PLATON/scint/stepMax The maximal step length in the scintillator
 * - /PLATON/scint/lightYield The light yield of the scintillator in photons/MeV (default 8700 from EJ-262)
 * - /PLATON/scint/decayTime The decay time of the scintillator in ns (default 2.1 from EJ-262)
 * - /PLATON/scint/sizeX The size of the scintillator in x direction in mm (the cube is always centered around 0,0,0)
 * - /PLATON/scint/sizeY The size of the scintillator in y direction in mm (the cube is always centered around 0,0,0)
 * - /PLATON/scint/sizeZ The size of the scintillator in z direction in mm (the cube is always centered around 0,0,0)
 */


class DetectorConstruction;
class DetectorMessenger : public G4UImessenger {
    public:
	DetectorMessenger(DetectorConstruction*); ///< Constructor
	~DetectorMessenger() override; ///< Destructor

	void SetNewValue(G4UIcommand*, G4String); ///< Set the new value of one fo the properties set by the messenger. This function is called by the geometry initialization and calls the Setters of the DetectorConstruction class.

    private:
	DetectorConstruction* fDetectorConstruction = nullptr; ///< Pointer to the detector construction class

	G4UIdirectory* fDirectory = nullptr; ///< Directory for the commands (for book keeping stuff)
	G4UIdirectory* fDetDirectory = nullptr; ///< Directory for the detector commands

	G4UIcmdWithAString* fScintMatCmd = nullptr; ///< Command to set the target material

	G4UIcmdWithADoubleAndUnit* fStepMaxCmd = nullptr; ///< Command to set the maximal step length in the scintillator

	G4UIcmdWithADouble* fLightYieldCmd = nullptr; ///< Command to set the light yield of the scintillator
	G4UIcmdWithADouble* fScintillationDecayTimeCmd = nullptr; ///< Command to set the decay time of the scintillator

	G4UIcmdWithADoubleAndUnit* fScintillatorSizeXCmd = nullptr; ///< Command to set the size of the scintillator in x direction
	G4UIcmdWithADoubleAndUnit* fScintillatorSizeYCmd = nullptr;	///< Command to set the size of the scintillator in y direction
	G4UIcmdWithADoubleAndUnit* fScintillatorSizeZCmd = nullptr;	///< Command to set the size of the scintillator in z direction

	G4UIcmdWithADoubleAndUnit* fScintillatorGranuCmd = nullptr;	///< Command to set the size of the scintillator in z direction

	G4UIcmdWithABool* fSaveScintillatorInformation = nullptr; ///< Command to set if the detector information should be saved

	G4UIcmdWithADoubleAndUnit* ftargetWSizeXCmd = nullptr; 
	G4UIcmdWithADoubleAndUnit* ftargetWSizeYCmd = nullptr;	
	G4UIcmdWithADoubleAndUnit* ftargetWSizeZCmd = nullptr;	

	G4UIcmdWithADoubleAndUnit* fLOSShiftXCmd = nullptr;	
	G4UIcmdWithADoubleAndUnit* fLOSShiftYCmd = nullptr;	

	G4UIcmdWithAnInteger * fNumberReplicasCmd = nullptr;



};


#endif
