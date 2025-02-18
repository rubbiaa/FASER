#ifndef CustomCharmDecay_H
#define CustomCharmDecay_H

#include "G4VExtDecayer.hh"
#include "globals.hh"

class G4Track;
class G4DecayProducts;

class CustomCharmDecay : public G4VExtDecayer
{
  
   public:

      //ctor & dtor
      CustomCharmDecay();
      virtual ~CustomCharmDecay();

      virtual G4DecayProducts* ImportDecayProducts(const G4Track&);
    
   private:

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
