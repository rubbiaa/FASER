#ifndef MuonDISMessenger_h
#define MuonDISMessenger_h 1

// Macro/UI control for MuonDISPhysics, following the same G4UImessenger
// pattern as PrimaryGeneratorMessenger / DetectorMessenger elsewhere in this
// app.
//
// IMPORTANT: unlike /generator/... commands (which may run any time before
// /run/beamOn), every command here configures the physics constructor and
// MUST be issued in your macro BEFORE /run/initialize -- ConstructProcess()
// (where the MuonDIS process actually gets installed on mu-/mu+) runs during
// initialization, and macro commands issued afterwards will have no effect
// on it. All commands are therefore restricted to G4State_PreInit.

#include "G4UImessenger.hh"
#include "globals.hh"

class MuonDISPhysics;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;

class MuonDISMessenger : public G4UImessenger {
public:
  explicit MuonDISMessenger(MuonDISPhysics* physics);
  ~MuonDISMessenger() override;

  void SetNewValue(G4UIcommand* command, G4String newValue) override;

private:
  MuonDISPhysics* fPhysics{};
  G4UIdirectory* fDirectory{};
  G4UIcmdWithABool* fEnableCmd{};
  G4UIcmdWithADouble* fCrossSectionBiasCmd{};
  G4UIcmdWithADouble* fQ2MinCmd{};
  G4UIcmdWithAString* fInteractionLogCmd{};
  G4UIcmdWithABool* fDebugCmd{};
};

#endif
