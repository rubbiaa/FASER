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
    G4bool istau = std::abs(pdgid) == 15;

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
    G4ThreeVector decay_position = track.GetPosition();

    if (istau)
    {
       // get decay products from the POEvent - assume the particle hasn't lost too much energy in dE/dx
       for (size_t i = 0; i < TPOevent->n_taudecay(); i++)
       {
          struct PO aPO = TPOevent->taudecay[i];
          // insert decay position
          aPO.m_vx_decay = decay_position.x();
          aPO.m_vy_decay = decay_position.y();
          aPO.m_vz_decay = decay_position.z();
          if (TPOevent->is_neutrino(aPO.m_pdg_id))
             continue;
          G4ParticleDefinition *pddec =
              G4ParticleTable::GetParticleTable()->FindParticle(aPO.m_pdg_id);
          if (!pddec)
             continue; // maybe we should print out a warning !
          G4ThreeVector momentum = G4ThreeVector(aPO.m_px * CLHEP::GeV,
                                                 aPO.m_py * CLHEP::GeV,
                                                 aPO.m_pz * CLHEP::GeV);
          dproducts->PushProducts(new G4DynamicParticle(pddec, momentum));
       }
    }
    else
    {
       // get decay products from the POEvent - assume the particle hasn't lost too much energy in dE/dx
       for (size_t i = 0; i < TPOevent->charmdecay.size(); i++)
       {
          struct PO aPO = TPOevent->charmdecay[i];
          // insert decay position
          aPO.m_vx_decay = decay_position.x();
          aPO.m_vy_decay = decay_position.y();
          aPO.m_vz_decay = decay_position.z();
          if (TPOevent->is_neutrino(aPO.m_pdg_id))
             continue;
          G4ParticleDefinition *pddec =
              G4ParticleTable::GetParticleTable()->FindParticle(aPO.m_pdg_id);
          if (!pddec)
             continue; // maybe we should print out a warning !
          G4ThreeVector momentum = G4ThreeVector(aPO.m_px * CLHEP::GeV,
                                                 aPO.m_py * CLHEP::GeV,
                                                 aPO.m_pz * CLHEP::GeV);
          dproducts->PushProducts(new G4DynamicParticle(pddec, momentum));
       };
    }

   return dproducts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
