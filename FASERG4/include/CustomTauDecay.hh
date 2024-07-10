#ifndef CustomTauDecay_H
#define CustomTauDecay_H

#include "G4VExtDecayer.hh"
#include "globals.hh"

class G4Track;
class G4DecayProducts;

class CustomTauDecay : public G4VExtDecayer
{
  
   public:

      //ctor & dtor
      CustomTauDecay();
      virtual ~CustomTauDecay();

      virtual G4DecayProducts* ImportDecayProducts(const G4Track&);
    
   private:

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
