#ifndef MuonDISNuclearWrapperProcess_h
#define MuonDISNuclearWrapperProcess_h 1

// Ported near-verbatim from the Calypso MuonDIS Geant4 extension -- this class
// has no Athena/HepMC3 dependency, only Geant4, so it carries over unchanged.
//
// Wraps a pair of G4MuonNuclearProcess instances (one "biased" -- cross-section
// boosted and using the Pythia8-driven MuonDISMuonNuclearModel -- and one
// "unbiased" -- the original stock Geant4 process) and dispatches per-track to
// whichever is appropriate: biased for primary muons (ParentID()==0) only,
// unbiased for everything else (secondaries, muons produced later in the
// shower). This keeps the MuonDIS enhancement scoped to the primary beam muon,
// as intended.

#include "G4WrapperProcess.hh"

class G4MuonNuclearProcess;
class G4Track;
class G4Step;

class MuonDISNuclearWrapperProcess final : public G4WrapperProcess {
public:
  MuonDISNuclearWrapperProcess(const G4String& name, G4MuonNuclearProcess* biasedProcess,
                               G4MuonNuclearProcess* unbiasedProcess, bool enableDebug);
  ~MuonDISNuclearWrapperProcess() override;

  G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                G4double previousStepSize,
                                                G4ForceCondition* condition) override;
  G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& stepData) override;

  G4bool IsApplicable(const G4ParticleDefinition& particle) override;
  void BuildPhysicsTable(const G4ParticleDefinition& particle) override;
  void PreparePhysicsTable(const G4ParticleDefinition& particle) override;
  G4bool StorePhysicsTable(const G4ParticleDefinition* particle, const G4String& directory,
                           G4bool ascii = false) override;
  G4bool RetrievePhysicsTable(const G4ParticleDefinition* particle, const G4String& directory,
                              G4bool ascii = false) override;
  void StartTracking(G4Track* track) override;
  void EndTracking() override;
  void SetProcessManager(const G4ProcessManager* procMan) override;
  const G4ProcessManager* GetProcessManager() override;
  void ResetNumberOfInteractionLengthLeft() override;
  void SetMasterProcess(G4VProcess* masterP) override;

private:
  G4bool useBiasedProcess(const G4Track& track) const;
  G4MuonNuclearProcess* selectProcess(const G4Track& track) const;
  void demoteSurvivingPrimary(const G4Track& track, const G4Step& stepData,
                              G4VParticleChange& particleChange) const;

  G4MuonNuclearProcess* m_biasedProcess{};
  G4MuonNuclearProcess* m_unbiasedProcess{};
  const G4ProcessManager* m_processManager{};
  bool m_enableDebug{false};
};

#endif
