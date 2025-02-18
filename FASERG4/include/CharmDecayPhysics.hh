#ifndef CharmDecayerPhysics_H
#define CharmDecayerPhysics_H

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class G4Decay;

class CharmDecayPhysics: public G4VPhysicsConstructor
{
  public:
    CharmDecayPhysics(G4int verb=1);
    virtual ~CharmDecayPhysics();

  protected:
    // methods
    // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
