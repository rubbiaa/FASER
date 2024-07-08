#include "EventAction.hh"

EventAction::EventAction(ParticleManager* photonManager) : fParticleManager(photonManager)
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
	std::cout << "EventAction::BeginOfEventAction" << std::endl;
	fParticleManager->beginOfEvent();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
    std::cout << "EventAction::EndOfEventAction" << std::endl;

   fParticleManager->endOfEvent(event);
}
