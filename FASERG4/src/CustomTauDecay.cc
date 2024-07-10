#include "CustomTauDecay.hh"

#include "G4DynamicParticle.hh"
#include "G4ParticleTable.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

#include "G4RunManager.hh"
#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CustomTauDecay::CustomTauDecay()
   : G4VExtDecayer("CustomTauDecay")
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CustomTauDecay::~CustomTauDecay()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DecayProducts* CustomTauDecay::ImportDecayProducts(const G4Track& track)
{

    G4cout << "CustomTauDecay ImportDecayProducts" << G4endl;

    G4DecayProducts* dproducts = nullptr;   
   
    G4ParticleDefinition* pd = track.GetDefinition();
    int    pdgid   = pd->GetPDGEncoding();

   G4cout << pdgid << " " << track.GetMomentum().x() / CLHEP::GeV <<
                " " << track.GetMomentum().y() / CLHEP::GeV  <<
                " " << track.GetMomentum().z() / CLHEP::GeV  <<
                " " << track.GetDynamicParticle()->GetTotalEnergy() / CLHEP::GeV <<
                " " << pd->GetPDGMass() / CLHEP::GeV << G4endl;

    const PrimaryGeneratorAction* primaryGenAction = 
    dynamic_cast<const PrimaryGeneratorAction*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
    const TPOEvent *TPOevent = primaryGenAction->GetTPOEvent();
    TPOevent->dump_event();

   // create & fill up decay products
    dproducts = new G4DecayProducts(*(track.GetDynamicParticle()));
    
    // get decay products from the POEvent - assume the particle hasn't lost too much energy in dE/dx
    for (size_t i=0; i<TPOevent->n_taudecay(); i++) {
        struct PO aPO = TPOevent->taudecay[i];
        if(TPOevent->is_neutrino(aPO.m_pdg_id))continue;
        G4ParticleDefinition* pddec = 
            G4ParticleTable::GetParticleTable()->FindParticle(aPO.m_pdg_id);
        if ( !pddec ) continue; // maybe we should print out a warning !
        G4ThreeVector momentum = G4ThreeVector( aPO.m_px * CLHEP::GeV,
                                                aPO.m_py * CLHEP::GeV,
                                                aPO.m_pz * CLHEP::GeV ); 
        dproducts->PushProducts( new G4DynamicParticle( pddec, momentum) ); 
    };

#if 0   
   // create G4DynamicParticle out of each fDecayer->event entry (except the 1st one)
   // and push into dproducts
   
   for ( int ip=npart_before_decay; ip<npart_after_decay; ++ip )
   {
      
      // only select final state decay products (direct or via subsequent decays);
      // skip all others
      //
      // NOTE: in general, final state decays products will have 
      //       positive status code between 91 and 99 
      //       (in case such information could be of interest in the future)
      //
      if ( fDecayer->event[ip].status() < 0 ) continue;
            
      // should we also skip neutrinos ???
      // if so, skip products with abs(fDecayer->event[ip].id()) of 12, 14, or 16
            
      G4ParticleDefinition* pddec = 
         G4ParticleTable::GetParticleTable()->FindParticle( fDecayer->event[ip].id() );
      if ( !pddec ) continue; // maybe we should print out a warning !
      G4ThreeVector momentum = G4ThreeVector( fDecayer->event[ip].px() * CLHEP::GeV,
                                              fDecayer->event[ip].py() * CLHEP::GeV,
                                              fDecayer->event[ip].pz() * CLHEP::GeV ); 
      dproducts->PushProducts( new G4DynamicParticle( pddec, momentum) ); 
   }
#endif
   
   return dproducts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
