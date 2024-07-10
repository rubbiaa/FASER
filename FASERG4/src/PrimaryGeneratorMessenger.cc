#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* Gun) : fAction(Gun)
{
	fGunDir = new G4UIdirectory("/generator/");
	fGunDir->SetGuidance("generator control");

	fROOTInputFileNameCmd = new G4UIcmdWithAString("/generator/rootinputfilename", this);
	fROOTInputFileNameCmd->SetGuidance("Select ROOT input file name.");
	fROOTInputFileNameCmd->SetParameterName("ROOTInputFileName", false);
	fROOTInputFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fFileNumberCmd = new G4UIcmdWithAnInteger("/generator/filenumber", this);
	fFileNumberCmd->SetGuidance("Select file number.");
	fFileNumberCmd->SetParameterName("FileNumber", false);
	fFileNumberCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fNEventsPerFileCmd = new G4UIcmdWithAnInteger("/generator/neventsperfile", this);
	fNEventsPerFileCmd->SetGuidance("Select number of events per file.");
	fNEventsPerFileCmd->SetParameterName("NEventsPerFile", false);
	fNEventsPerFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
	delete fROOTInputFileNameCmd;
	delete fFileNumberCmd;
	delete fNEventsPerFileCmd;
	delete fGunDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if (command == fROOTInputFileNameCmd) {
		fAction->SetROOTInputFileName(newValue);
	}

	if (command == fFileNumberCmd) {
		fAction->SetFileNumber(fFileNumberCmd->GetNewIntValue(newValue));
	}

	if (command == fNEventsPerFileCmd) {
		fAction->SetNEventsPerFile(fNEventsPerFileCmd->GetNewIntValue(newValue));
	}
}