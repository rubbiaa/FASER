#include "ActionInitialization.hh"


ActionInitialization::ActionInitialization(ParticleManager* pm)
    : G4VUserActionInitialization(), fParticleManager(pm) {}  // Initialize the member variable in the constructor

void ActionInitialization::BuildForMaster() const { SetUserAction(new RunAction(fParticleManager)); }

void ActionInitialization::Build() const
{
	SetUserAction(new PrimaryGeneratorAction(fParticleManager));
	SetUserAction(new RunAction(fParticleManager));
	SetUserAction(new EventAction(fParticleManager));
	SetUserAction(new TrackingAction(fParticleManager));
}
