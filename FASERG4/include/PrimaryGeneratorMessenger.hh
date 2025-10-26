#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
// added by Umut
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "PrimaryGeneratorAction.hh"

class PrimaryGeneratorAction;

/**
 * @class PrimaryGeneratorMessenger
 * @brief Messenger for the primary generator action
 *
 * This class is used to define the commands that can be used to control the
 * primary generator action.
 * It contains the following commands:
 * - /generator/rootinputfilename filename : set the name of the ROOT file
 * - /generator/filenumber number : Determines where to read in the ROOT file
 * - /generator/neventsperfile number : Determines how many events to read in and process
 */

class PrimaryGeneratorMessenger: public G4UImessenger
{
  public:

    PrimaryGeneratorMessenger(PrimaryGeneratorAction*); ///< Constructor
   ~PrimaryGeneratorMessenger(); ///< Destructor

    void SetNewValue(G4UIcommand*, G4String); ///< Set a new value for the variables

  private:

    PrimaryGeneratorAction* fAction; ///< Pointer to the primary generator action
    G4UIdirectory* fGunDir;	  ///< Pointer to the UI directory
    G4UIcmdWithAString* fROOTInputFileNameCmd; ///< Input command for the ROOT file name
    G4UIcmdWithAnInteger* fFileNumberCmd; ///< Input command for the file number
    G4UIcmdWithAnInteger* fNStartEvent; ///< Input command for the number of events per file
  // added by Umut
    G4UIcmdWithABool* fWantMuonBackground; ///< Input command for muon background option
    G4UIcmdWithABool* fWantSingleParticle; ///< Input command for single particle option
    G4UIcmdWithADoubleAndUnit* fSingleMomentumCmd; ///< Command to set single particle momentum (with unit)
};

#endif
