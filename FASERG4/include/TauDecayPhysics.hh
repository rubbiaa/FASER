#ifndef TauDecayerPhysics_H
#define TauDecayerPhysics_H

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class G4Decay;

class TauDecayPhysics: public G4VPhysicsConstructor
{
  public:
    TauDecayPhysics(G4int verb=1);
    virtual ~TauDecayPhysics();

  protected:
    // methods
    // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
