#include "DetectorMessenger.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* det) : fDetectorConstruction(det)
{
	G4cout << "DetectorMessenger::DetectorMessenger" << G4endl;
	fDirectory = new G4UIdirectory("/FASER/");
	fDirectory->SetGuidance("UI commands specific to this example.");

	fDetDirectory = new G4UIdirectory("/FASER/scint/");
	fDetDirectory->SetGuidance("Detector construction control");

	fScintMatCmd = new G4UIcmdWithAString("/FASER/scint/material", this);
	fScintMatCmd->SetGuidance("Select Material of the Target.");
	fScintMatCmd->SetParameterName("choice", false);
	fScintMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fLightYieldCmd = new G4UIcmdWithADouble("/FASER/scint/lightYield", this);
	fLightYieldCmd->SetGuidance("Set scintillation light yield");
	fLightYieldCmd->SetParameterName("lightYield", false);
	fLightYieldCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fScintillationDecayTimeCmd = new G4UIcmdWithADouble("/FASER/scint/decayTime", this);
	fScintillationDecayTimeCmd->SetGuidance("Set scintillation decay time");
	fScintillationDecayTimeCmd->SetParameterName("decayTime", false);
	fScintillationDecayTimeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fScintillatorSizeXCmd = new G4UIcmdWithADoubleAndUnit("/FASER/scint/sizeX", this);
	fScintillatorSizeXCmd->SetGuidance("Set scintillator size");
	fScintillatorSizeXCmd->SetParameterName("size", false);
	fScintillatorSizeXCmd->SetUnitCategory("Length");
	fScintillatorSizeXCmd->SetDefaultUnit("mm");
	fScintillatorSizeXCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fScintillatorSizeYCmd = new G4UIcmdWithADoubleAndUnit("/FASER/scint/sizeY", this);
	fScintillatorSizeYCmd->SetGuidance("Set scintillator size");
	fScintillatorSizeYCmd->SetParameterName("size", false);
	fScintillatorSizeYCmd->SetUnitCategory("Length");
	fScintillatorSizeYCmd->SetDefaultUnit("mm");
	fScintillatorSizeYCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fScintillatorSizeZCmd = new G4UIcmdWithADoubleAndUnit("/FASER/scint/sizeZ", this);
	fScintillatorSizeZCmd->SetGuidance("Set scintillator size");
	fScintillatorSizeZCmd->SetParameterName("size", false);
	fScintillatorSizeZCmd->SetUnitCategory("Length");
	fScintillatorSizeZCmd->SetDefaultUnit("mm");
	fScintillatorSizeZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fScintillatorGranuCmd = new G4UIcmdWithADoubleAndUnit("/FASER/scint/voxel", this);
	fScintillatorGranuCmd->SetGuidance("Set scintillator readout voxel size");
	fScintillatorGranuCmd->SetParameterName("size", false);
	fScintillatorGranuCmd->SetUnitCategory("Length");
	fScintillatorGranuCmd->SetDefaultUnit("mm");
	fScintillatorGranuCmd->AvailableForStates(G4State_PreInit, G4State_Idle);


	ftargetWSizeXCmd = new G4UIcmdWithADoubleAndUnit("/FASER/targetW/sizeX", this);
	ftargetWSizeXCmd->SetGuidance("Set targetW size");
	ftargetWSizeXCmd->SetParameterName("size", false);
	ftargetWSizeXCmd->SetUnitCategory("Length");
	ftargetWSizeXCmd->SetDefaultUnit("mm");
	ftargetWSizeXCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	ftargetWSizeYCmd = new G4UIcmdWithADoubleAndUnit("/FASER/targetW/sizeY", this);
	ftargetWSizeYCmd->SetGuidance("Set targetW size");
	ftargetWSizeYCmd->SetParameterName("size", false);
	ftargetWSizeYCmd->SetUnitCategory("Length");
	ftargetWSizeYCmd->SetDefaultUnit("mm");
	ftargetWSizeYCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	ftargetWSizeZCmd = new G4UIcmdWithADoubleAndUnit("/FASER/targetW/sizeZ", this);
	ftargetWSizeZCmd->SetGuidance("Set targetW size");
	ftargetWSizeZCmd->SetParameterName("size", false);
	ftargetWSizeZCmd->SetUnitCategory("Length");
	ftargetWSizeZCmd->SetDefaultUnit("mm");
	ftargetWSizeZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fNumberReplicasCmd = new G4UIcmdWithAnInteger("/FASER/layers", this);
	fNumberReplicasCmd->SetGuidance("Set number of targetW-scintillator layers");
	fNumberReplicasCmd->SetParameterName("size", false);
	fNumberReplicasCmd->AvailableForStates(G4State_PreInit, G4State_Idle);



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
	delete fScintMatCmd;
	delete fStepMaxCmd;
	delete fLightYieldCmd;
	delete fScintillationDecayTimeCmd;
	delete fScintillatorSizeXCmd;
	delete fScintillatorSizeYCmd;
	delete fScintillatorSizeZCmd;
	delete fScintillatorGranuCmd;
	delete ftargetWSizeXCmd;
	delete ftargetWSizeYCmd;
	delete ftargetWSizeZCmd;
	delete fNumberReplicasCmd;
	delete fDirectory;
	delete fDetDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if (command == fScintMatCmd) {
		fDetectorConstruction->SetScintillatorMaterial(newValue);
	}

	if (command == fStepMaxCmd) {
		fDetectorConstruction->SetMaxStep(fStepMaxCmd->GetNewDoubleValue(newValue));
	}

	if (command == fLightYieldCmd) {
		fDetectorConstruction->SetLightYield(fLightYieldCmd->GetNewDoubleValue(newValue));
	}

	if (command == fScintillationDecayTimeCmd) {
		fDetectorConstruction->SetScintillationDecayTime(fScintillationDecayTimeCmd->GetNewDoubleValue(newValue));
	}

	if (command == fScintillatorSizeXCmd) {
		fDetectorConstruction->SetScintillatorSizeX(fScintillatorSizeXCmd->GetNewDoubleValue(newValue));
	}

	if (command == fScintillatorSizeYCmd) {
		fDetectorConstruction->SetScintillatorSizeY(fScintillatorSizeYCmd->GetNewDoubleValue(newValue));
	}

	if (command == fScintillatorSizeZCmd) {
		fDetectorConstruction->SetScintillatorSizeZ(fScintillatorSizeZCmd->GetNewDoubleValue(newValue));
	}

	if (command == fScintillatorGranuCmd) {
		fDetectorConstruction->SetVoxelSize(fScintillatorGranuCmd->GetNewDoubleValue(newValue));
	}

	if (command == ftargetWSizeXCmd) {
		fDetectorConstruction->SettargetWSizeX(ftargetWSizeXCmd->GetNewDoubleValue(newValue));
	}

	if (command == ftargetWSizeYCmd) {
		fDetectorConstruction->SettargetWSizeY(ftargetWSizeYCmd->GetNewDoubleValue(newValue));
	}

	if (command == ftargetWSizeZCmd) {
		fDetectorConstruction->SettargetWSizeZ(ftargetWSizeZCmd->GetNewDoubleValue(newValue));
	}

	if (command == fNumberReplicasCmd) {
		fDetectorConstruction->SetNumberReplicas(fNumberReplicasCmd->GetNewIntValue(newValue));
	}




}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
