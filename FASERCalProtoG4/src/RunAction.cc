#include "RunAction.hh"
#include "TrackerSD.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(ParticleManager* photonManager): fParticleManager(photonManager)
{
  // set printing event number per each 100 events
  G4RunManager::GetRunManager()->SetPrintProgress(1);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  fParticleManager->beginOfRun();
  
  // Open muon tracking file
  TrackerSD::OpenMuonTrackingFile("muon_tracking.dat");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* )
{
   std::cout << "RunAction::EndOfRunAction" << std::endl;
  fParticleManager->endOfRun();
  
  // Close muon tracking file
  TrackerSD::CloseMuonTrackingFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
