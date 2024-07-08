#ifndef TrackerSD_h
#define TrackerSD_h 1

#include <vector>

#include "G4VSensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"


#include "ParticleManager.hh"
#include "Geant4Process.hh"
#include "typedef.h"

class G4Step;
class G4HCofThisEvent;

/**
 * @class TrackerSD
 * @brief Tracker sensitive detector class
 *
 * The hits are accounted in hits in ProcessHits() function which is called
 * by Geant4 kernel at each step. A hit is created with each step with non zero
 * energy deposit.
 */

class TrackerSD : public G4VSensitiveDetector {
    public:
	/**
	 * @brief Constructor.
	 * @param name Name of the sensitive detector.
	 * @param pParticleManager Particle manager, which is used to store and organize the hits.
	 */
	TrackerSD(const G4String& name, ParticleManager* pParticleManager);
	~TrackerSD();  ///< Destructor.

	/**
	 * @brief Processes the hits.
	 * @details This function is called by Geant4 kernel at each step. It internally calls the ParticleManager class functions to store the hits
	 * of the particles and photons
	 * @param step Step object.
	 * @param history History of the step.
	 */
	G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;

    private:
	ParticleManager* fParticleManager = nullptr;  ///< Particle manager, which is used to store and organize the hits.
};


#endif
