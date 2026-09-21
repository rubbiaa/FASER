#ifndef MuonDISMuonNuclearModel_h
#define MuonDISMuonNuclearModel_h 1

// Adapted from the Calypso MuonDIS package's MuonDISMuonNuclearModel: same
// G4HadronicInteraction role (replaces the standard muon-nuclear final state
// for the primary muon), but the final state now comes from an on-the-fly
// Pythia8 call (MuonDISPythiaGenerator) instead of a lookup into a
// pre-generated POWHEG/HepMC3 event library.

#include "G4HadronicInteraction.hh"

#include <iosfwd>
#include <utility>

class MuonDISMuonNuclearModel : public G4HadronicInteraction {
public:
  explicit MuonDISMuonNuclearModel(const G4String& modelName = "MuonDISMuonNuclearModel",
                                   bool enableDebug = false);
  ~MuonDISMuonNuclearModel() override;

  G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
                                 G4Nucleus& targetNucleus) override;
  void ModelDescription(std::ostream& outFile) const override;
  const std::pair<G4double, G4double> GetFatalEnergyCheckLevels() const override;

private:
  bool m_enableDebug{false};
};

#endif
