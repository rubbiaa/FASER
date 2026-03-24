#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"
// added by Umut
#include "G4SystemOfUnits.hh"


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

	fNStartEvent = new G4UIcmdWithAnInteger("/generator/startevent", this);
	fNStartEvent->SetGuidance("Select start event.");
	fNStartEvent->SetParameterName("StartEvent", false);
	fNStartEvent->AvailableForStates(G4State_PreInit, G4State_Idle);

	fWantMuonBackground = new G4UIcmdWithABool("/generator/wantMuonBackground", this);
	fWantMuonBackground->SetGuidance("Enable muon background mode (transport only)");
	fWantMuonBackground->SetParameterName("wantMuonBackground", false);
	fWantMuonBackground->AvailableForStates(G4State_PreInit, G4State_Idle);
	
	fWantSingleParticle = new G4UIcmdWithABool("/generator/wantSingleParticle", this);
	fWantSingleParticle->SetGuidance("Enable single particle mode (use only 1 particle from event)");
	fWantSingleParticle->SetParameterName("wantSingleParticle", false);
	fWantSingleParticle->AvailableForStates(G4State_PreInit, G4State_Idle);

	// Umut::adding for single particle momentum command
	fSingleMomentumCmd = new G4UIcmdWithADoubleAndUnit("/generator/singleMomentum", this);
	fSingleMomentumCmd->SetGuidance("Set momentum magnitude for single-particle mode (with unit, e.g. 50 GeV)");
	fSingleMomentumCmd->SetParameterName("SingleMomentum", false);
	fSingleMomentumCmd->SetDefaultUnit("GeV");
	fSingleMomentumCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fSingleParticleNameCmd = new G4UIcmdWithAString("/generator/singleParticleName", this);
	fSingleParticleNameCmd->SetGuidance("Set single-particle species, e.g. neutron, lambda0, k0, k0s, k0l");
	fSingleParticleNameCmd->SetParameterName("SingleParticleName", false);
	fSingleParticleNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fUseNeutralHadronLogSpectrumCmd = new G4UIcmdWithABool("/generator/useNeutralHadronLogSpectrum", this);
	fUseNeutralHadronLogSpectrumCmd->SetGuidance("Enable log-uniform kinetic energy spectrum for neutral hadrons");
	fUseNeutralHadronLogSpectrumCmd->SetParameterName("UseNeutralHadronLogSpectrum", false);
	fUseNeutralHadronLogSpectrumCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fNeutralHadronLogEminCmd = new G4UIcmdWithADoubleAndUnit("/generator/neutralHadronLogEmin", this);
	fNeutralHadronLogEminCmd->SetGuidance("Set minimum kinetic energy for neutral-hadron log spectrum");
	fNeutralHadronLogEminCmd->SetParameterName("NeutralHadronLogEmin", false);
	fNeutralHadronLogEminCmd->SetDefaultUnit("GeV");
	fNeutralHadronLogEminCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fNeutralHadronLogEmaxCmd = new G4UIcmdWithADoubleAndUnit("/generator/neutralHadronLogEmax", this);
	fNeutralHadronLogEmaxCmd->SetGuidance("Set maximum kinetic energy for neutral-hadron log spectrum");
	fNeutralHadronLogEmaxCmd->SetParameterName("NeutralHadronLogEmax", false);
	fNeutralHadronLogEmaxCmd->SetDefaultUnit("GeV");
	fNeutralHadronLogEmaxCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
	delete fROOTInputFileNameCmd;
	delete fFileNumberCmd;
	delete fNStartEvent;
	delete fWantMuonBackground;
	delete fWantSingleParticle;	
	delete fSingleMomentumCmd;
	delete fSingleParticleNameCmd;
	delete fUseNeutralHadronLogSpectrumCmd;
	delete fNeutralHadronLogEminCmd;
	delete fNeutralHadronLogEmaxCmd;
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

	if (command == fNStartEvent) {
		fAction->SetfNStartEvent(fNStartEvent->GetNewIntValue(newValue));
	}

	if (command == fWantMuonBackground) {
		fAction->SetWantMuonBackground(fWantMuonBackground->GetNewBoolValue(newValue));
	}
	
	if (command == fWantSingleParticle) {
		fAction->SetWantSingleParticle(fWantSingleParticle->GetNewBoolValue(newValue));
	}

	if (command == fSingleParticleNameCmd) {
		fAction->SetSingleParticleName(newValue);
	}

	if (command == fUseNeutralHadronLogSpectrumCmd) {
		fAction->SetUseNeutralHadronLogSpectrum(fUseNeutralHadronLogSpectrumCmd->GetNewBoolValue(newValue));
	}

	if (command == fNeutralHadronLogEminCmd) {
		double valueGeV = fNeutralHadronLogEminCmd->GetNewDoubleValue(newValue) / GeV;
		fAction->SetNeutralHadronLogEmin(valueGeV);
	}

	if (command == fNeutralHadronLogEmaxCmd) {
		double valueGeV = fNeutralHadronLogEmaxCmd->GetNewDoubleValue(newValue) / GeV;
		fAction->SetNeutralHadronLogEmax(valueGeV);
	}
	
	// added by Umut
	if (command == fSingleMomentumCmd) {
		// log the raw string passed by UI before conversion to double
		G4cout << "PrimaryGeneratorMessenger: /generator/singleMomentum raw newValue='" << newValue << "'" << G4endl;
		double val = fSingleMomentumCmd->GetNewDoubleValue(newValue);
		double val_raw = val; // in internal units (e.g. MeV)
		double val_GeV = val_raw / GeV;
		G4cout << "PrimaryGeneratorMessenger: parsed raw value = " << val_raw << " (internal units), converted to " << val_GeV << " GeV" << G4endl;
		fAction->SetSingleParticleMomentum(val_GeV);
	}
}
