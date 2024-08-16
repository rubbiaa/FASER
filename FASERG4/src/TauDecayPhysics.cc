#include "TauDecayPhysics.hh"
#include "CustomTauDecay.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Decay.hh"
#include "G4DecayTable.hh"

// factory
//
#include "G4PhysicsConstructorFactory.hh"
//
// register it with contructor factory
//
G4_DECLARE_PHYSCONSTR_FACTORY(TauDecayPhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TauDecayPhysics::TauDecayPhysics(G4int)
  : G4VPhysicsConstructor("TauDecayPhysics")
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TauDecayPhysics::~TauDecayPhysics() 
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TauDecayPhysics::ConstructParticle()
{
   // Nothing needs to be done here
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TauDecayPhysics::ConstructProcess()
{
   // Adding external decayer to G4Decay process (per each thread).
   // G4Decay will use the external decayer if G4Decay process is
   // assigned to an unstable particle and that particle does not
   // have its decay table.

   // Loop over all particles instantiated and remove already-assigned
   // decay table for tau's and B+/- so that they will decay through
   // the external decayer (Pythia8).

   // NOTE: The extDecayer will be deleted in G4Decay destructor
   
   CustomTauDecay* extDecayer = new CustomTauDecay();
   G4bool setOnce = true;

   auto particleIterator=GetParticleIterator();
   particleIterator->reset();
   while ((*particleIterator)())
   {    
      G4ParticleDefinition* particle = particleIterator->value();

      // remove native/existing decay table for tau's 
      // so that G4Decay will use the external decayer
      if ( std::abs(particle->GetPDGEncoding()) == 15 )
      {
        if ( particle->GetDecayTable() )
        {
          particle->DumpTable();
          delete particle->GetDecayTable();
          particle->SetDecayTable(nullptr);
          particle->DumpTable();
        }
      }

      if(setOnce)
      // One G4Decay object is shared by all unstable particles (per thread).
      // Thus, we set the external decayer only once.
      {
        G4ProcessManager* pmanager = particle->GetProcessManager();    
        G4ProcessVector* processVector = pmanager->GetProcessList();
        for ( size_t i=0; i<processVector->length(); ++i ) 
        {    
           G4Decay* decay = dynamic_cast<G4Decay*>((*processVector)[i]);
           if ( decay ) 
           {
             decay->SetExtDecayer(extDecayer);
             setOnce = false;
           }
        }
      }              
   }

   return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
