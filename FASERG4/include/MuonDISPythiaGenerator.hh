#ifndef MuonDISPythiaGenerator_h
#define MuonDISPythiaGenerator_h 1

// On-the-fly muon-nucleon deep inelastic scattering event generation using
// Pythia8, replacing the pre-generated POWHEG/HepMC3 event library approach
// from the original Calypso MuonDIS package.
//
// *** PLEASE REVIEW THE PHYSICS SETTINGS BELOW BEFORE TRUSTING THE OUTPUT ***
// This is a first cut, based on Pythia8's generic electroweak process
// library (verified against the online Pythia8 manual: "Electroweak
// Processes" and "Beam Parameters" pages), not yet cross-checked against a
// running Pythia8 install (Pythia8 is not present in the sandbox this was
// written in). Things you should verify once you can build & run this:
//   - WeakBosonExchange:ff2ff(t:gmZ) is the standard t-channel photon/Z
//     exchange process (verified: Pythia8 manual, process code 211). It
//     needs no lepton PDF since the incoming lepton itself is a hard-process
//     leg, not resolved into partons.
//   - Charged-current (W exchange, mu -> nu_mu) is intentionally left OFF
//     (negligible at these Q^2/energies for a charged lepton beam); flip on
//     WeakBosonExchange:ff2ff(t:W) below if your physics case needs it.
//   - PhaseSpace:Q2Min (default 1.0 GeV^2 here, configurable via
//     /physics/muondis/q2min) keeps the process away from the divergent
//     photoproduction region; pick a value appropriate to your analysis.
//   - Beams:frameType = 3 (verified: Pythia8 manual, "Beam Parameters",
//     option 3) lets beam A be given an arbitrary 3-momentum -- this is used
//     to set beam A to the *actual* Geant4 muon energy/direction for every
//     single interaction, so Pythia8's lab-frame final state comes out
//     already aligned with the physical track. No POWHEG-style
//     generator-frame rotation is needed (contrast with the original
//     Calypso MuonDISMuonNuclearModel::rotateToMuonDirection()).
//   - Because both species (isoscalar proton/neutron pick) and kinematics
//     change essentially every call, this reinitializes Pythia8 (init())
//     for every interaction rather than using the variable-energy
//     setKinematics() machinery (which, per the manual, is documented for
//     energy changes at fixed frameType/species, not a full re-aim of beam
//     A's direction). If MuonDIS interaction rates make this a bottleneck,
//     look at Beams:allowVariableEnergy plus the setKinematics() overload
//     matching frameType=3 in your installed Pythia8's Pythia.h.
//
// Nuclear target treatment: per your instructions, the target nucleon (beam
// B) is picked event-by-event as a free proton or neutron with probability
// Z/A and (A-Z)/A respectively (isoscalar mix), at rest, ignoring Fermi
// motion, binding energy, shadowing and the EMC effect. Pythia8 ships no
// nuclear PDFs itself.

#include "G4ThreeVector.hh"

#include <memory>
#include <vector>

namespace Pythia8 {
  class Pythia;
}

class MuonDISPythiaGenerator {
public:
  struct GeneratedParticle {
    int pdgId{0};
    double px{0.0}, py{0.0}, pz{0.0}, e{0.0};  // GeV, in the same lab frame as the muon direction passed in
  };

  struct GeneratedEvent {
    bool valid{false};
    double q2GeV2{0.0};        // -tHat() of the hard 2->2 process; diagnostic only
    int targetNucleonPdg{0};   // 2212 or 2112 -- the isoscalar pick actually used this call
    std::vector<GeneratedParticle> finalState;
  };

  static MuonDISPythiaGenerator& instance();

  MuonDISPythiaGenerator(const MuonDISPythiaGenerator&) = delete;
  MuonDISPythiaGenerator& operator=(const MuonDISPythiaGenerator&) = delete;

  /// Must be called (idempotent) before the first generate(); safe to call again to
  /// change debug/Q2Min before they take effect on the next re-init.
  void configure(bool enableDebug, double q2MinGeV2);

  /// Generate one muon-nucleon DIS event.
  /// @param muonPdgId      13 (mu-) or -13 (mu+)
  /// @param muonEnergyGeV  total lab energy of the incoming muon, in GeV
  /// @param muonDirection  unit vector of the incoming muon direction, Geant4 lab frame
  /// @param targetZ, targetA  the struck nucleus, used only to pick the isoscalar nucleon
  GeneratedEvent generate(int muonPdgId, double muonEnergyGeV,
                          const G4ThreeVector& muonDirection, int targetZ, int targetA);

private:
  MuonDISPythiaGenerator() = default;
  ~MuonDISPythiaGenerator();

  void ensureInitialized(int muonPdgId, int nucleonPdgId, double muonEnergyGeV,
                         const G4ThreeVector& muonDirection);
  int chooseNucleon(int targetZ, int targetA) const;

  std::unique_ptr<Pythia8::Pythia> m_pythia;
  bool m_settingsApplied{false};
  bool m_enableDebug{false};
  double m_q2MinGeV2{1.0};
  int m_lastMuonPdg{0};
  int m_lastNucleonPdg{0};
};

#endif
