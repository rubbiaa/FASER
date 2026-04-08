#ifndef RunAction_h
#define RunAction_h 1

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UserRunAction.hh"

#include "ParticleManager.hh"

#include "globals.hh"

class G4Run;


/**
 * @class RunAction
 * @brief RunAction class, which is called at the beginning and end of each run.
 * Used to trigger the ParticleManager to write out the particle information.
 * Called inside the main simulation code
 */
class RunAction : public G4UserRunAction {
    private:
	ParticleManager* fParticleManager;  ///< Pointer to the ParticleManager
	

    public:
	RunAction(ParticleManager* photonManager);  ///< Constructor, which takes a pointer to the ParticleManager
	~RunAction() override = default;	    ///< Destructor

	void BeginOfRunAction(const G4Run* run) override;  ///< Called at the beginning of each run, calls ParticleManager::beginOfRun
	void EndOfRunAction(const G4Run* run) override;	   ///< Called at the end of each run, calls ParticleManager::endOfRun

};


#endif
