#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "ParticleManager.hh"
#include "globals.hh"


#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "RunAction.hh"


/**
 * @brief Handles the details of a single event
 *
 * This class takes care on the ways how a event is handled.
 * Two main functions can be used and are called by the G4 kernel:
 * BeginOfEventAction() and EndOfEventAction().
 * They do what they say, setting up everything for the event and cleaning up afterowrds.
 * Here they are mostly used to reset the particle manager and fill the trees inside the particle manager to save the simulation data.
 */
class EventAction : public G4UserEventAction {
    private:
	ParticleManager* fParticleManager; ///< Pointer to the particle manager

    public:
	EventAction(ParticleManager* photonManager); ///< Constructor, sets the pointer to the particle manager
	~EventAction() override = default; ///< Destructor

	void BeginOfEventAction(const G4Event*) override; ///< Called at the beginning of an event - calls particle Managers begin of event action
	void EndOfEventAction(const G4Event*) override; ///< Called at the end of an event - calls particle Managers end of event action
};


#endif
