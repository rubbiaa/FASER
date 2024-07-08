#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "ParticleManager.hh"

#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"



/**
 * @class ActionInitialization
 * @brief Initialize all user actions.
 */
class ActionInitialization : public G4VUserActionInitialization {
    private:
	ParticleManager* fParticleManager;  // The particle manager is used to store all particles and their tracks. It is used in the RunAction to
					    // write the output file and in the EventAction to fill the trees inside the root file.

    public:
	ActionInitialization(ParticleManager* photonManager);  // Constructor taking ParticleManager as an argument
	~ActionInitialization() override = default; // Destructor

	void BuildForMaster() const override; // Build for master thread
	void Build() const override; // Build for worker threads
};


#endif
