#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "FTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4SteppingVerbose.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "ParticleManager.hh"
#include "Randomize.hh"
#include "RunAction.hh"
#include "TFile.h"
#include "TH2F.h"


int main(int argc, char** argv)
{
	// get the output file name as the first argument
	if (argc != 3) {
		G4cout << "Usage: " << argv[0] << " <macro|vis>  <output.root>" << G4endl;
		return 1;
	}

	// Detect interactive mode (if no arguments) and define UI session

	G4UIExecutive* ui = nullptr;
	std::string fileName = argv[2];
	if (argv[1] == "vis") {
		ui = new G4UIExecutive(argc - 1, argv);
	}

	// Optionally: choose a different Random engine...
	G4Random::setTheEngine(new CLHEP::MTwistEngine);
	//G4long seed = time(NULL);
	G4long seed = 123456789;
	G4Random::setTheSeed(seed);

	// use G4SteppingVerboseWithUnits
	G4int precision = 4;
	G4SteppingVerbose::UseBestUnit(precision);

	// Create the ParticleManager
	auto particleManager = new ParticleManager(10, fileName);

	// Construct the default run manager
	//
	auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
	//
	auto detectorConstruction = new DetectorConstruction(particleManager);
	runManager->SetUserInitialization(detectorConstruction);

	G4StepLimiterPhysics* stepLimitPhys = new G4StepLimiterPhysics();
	auto physicsList = new FTFP_BERT;
	physicsList->RegisterPhysics(new G4StepLimiterPhysics());
	//	physicsList->RegisterPhysics(new G4OpticalPhysics());
	physicsList->RegisterPhysics(stepLimitPhys);
	runManager->SetUserInitialization(physicsList);
	// Set the ParticleManager in RunAction and EventAction
	//auto runAction = new RunAction(photonManager);
	//auto eventAction = new EventAction(photonManager, runAction);

	// Set the user actions
	auto actionInitialization = new ActionInitialization(particleManager);
	runManager->SetUserInitialization(actionInitialization);

	// Need to call the GeometryHasBeenModified() method as the macros can change the geometry
	runManager->GeometryHasBeenModified();

	// Initialize visualization with the default graphics system
	auto visManager = new G4VisExecutive(argc, argv);
	// Constructors can also take optional arguments:
	// - a graphics system of choice, eg. "OGL"
	// - and a verbosity argument - see /vis/verbose guidance.
	// auto visManager = new G4VisExecutive(argc, argv, "OGL", "Quiet");
	// auto visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();

	// Get the pointer to the User Interface manager
	auto UImanager = G4UImanager::GetUIpointer();

	// Process macro or start UI session
	//
	if (!ui) {
		// batch mode
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UImanager->ApplyCommand(command + fileName);
	}
	else {
		// interactive mode
		UImanager->ApplyCommand("/control/execute init_vis.mac");
		if (ui->IsGUI()) {
			UImanager->ApplyCommand("/control/execute gui.mac");
		}
		ui->SessionStart();
		delete ui;
	}

	// Job termination
	// Free the store: user actions, physics_list and detector_description are
	// owned and deleted by the run manager, so they should not be deleted
	// in the main() program !
	//
	delete visManager;
	delete runManager;
	delete particleManager;
}
