#include "MuonDISMessenger.hh"
#include "MuonDISPhysics.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

MuonDISMessenger::MuonDISMessenger(MuonDISPhysics* physics) : fPhysics(physics) {
  fDirectory = new G4UIdirectory("/physics/muondis/");
  fDirectory->SetGuidance("Muon deep-inelastic-scattering physics process control.");
  fDirectory->SetGuidance(
      "IMPORTANT: every command below must be issued BEFORE /run/initialize -- they "
      "configure the physics constructor, which is locked in during initialization.");

  fEnableCmd = new G4UIcmdWithABool("/physics/muondis/enable", this);
  fEnableCmd->SetGuidance(
      "Enable the Pythia8-driven MuonDIS muon-nuclear process (default: disabled).");
  fEnableCmd->SetParameterName("enable", false);
  fEnableCmd->AvailableForStates(G4State_PreInit);

  fCrossSectionBiasCmd = new G4UIcmdWithADouble("/physics/muondis/crossSectionBias", this);
  fCrossSectionBiasCmd->SetGuidance(
      "Bias factor applied to the muon-nuclear cross section for primary muons only "
      "(default: 150).");
  fCrossSectionBiasCmd->SetParameterName("bias", false);
  fCrossSectionBiasCmd->AvailableForStates(G4State_PreInit);

  fQ2MinCmd = new G4UIcmdWithADouble("/physics/muondis/q2min", this);
  fQ2MinCmd->SetGuidance(
      "Minimum Q^2 in GeV^2 passed to Pythia8 (PhaseSpace:Q2Min) for the DIS process "
      "(default: 1.0). Review/tune for your analysis.");
  fQ2MinCmd->SetParameterName("q2min", false);
  fQ2MinCmd->AvailableForStates(G4State_PreInit);

  fInteractionLogCmd = new G4UIcmdWithAString("/physics/muondis/interactionLog", this);
  fInteractionLogCmd->SetGuidance(
      "CSV path for exact MuonDIS Geant4 interaction points; leave unset to disable.");
  fInteractionLogCmd->SetParameterName("path", false);
  fInteractionLogCmd->AvailableForStates(G4State_PreInit);

  fDebugCmd = new G4UIcmdWithABool("/physics/muondis/debug", this);
  fDebugCmd->SetGuidance("Emit verbose MuonDIS process/Pythia8 installation logs.");
  fDebugCmd->SetParameterName("debug", false);
  fDebugCmd->AvailableForStates(G4State_PreInit);
}

MuonDISMessenger::~MuonDISMessenger() {
  delete fEnableCmd;
  delete fCrossSectionBiasCmd;
  delete fQ2MinCmd;
  delete fInteractionLogCmd;
  delete fDebugCmd;
  delete fDirectory;
}

void MuonDISMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if (command == fEnableCmd) {
    fPhysics->SetEnabled(fEnableCmd->GetNewBoolValue(newValue));
  } else if (command == fCrossSectionBiasCmd) {
    fPhysics->SetCrossSectionBias(fCrossSectionBiasCmd->GetNewDoubleValue(newValue));
  } else if (command == fQ2MinCmd) {
    fPhysics->SetQ2Min(fQ2MinCmd->GetNewDoubleValue(newValue));
  } else if (command == fInteractionLogCmd) {
    fPhysics->SetInteractionLogPath(newValue);
  } else if (command == fDebugCmd) {
    fPhysics->SetDebug(fDebugCmd->GetNewBoolValue(newValue));
  }
}
