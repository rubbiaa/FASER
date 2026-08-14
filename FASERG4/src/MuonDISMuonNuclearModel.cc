#include "MuonDISMuonNuclearModel.hh"

#include "MuonDISInteractionRecorder.hh"
#include "MuonDISPythiaGenerator.hh"

#include "G4DynamicParticle.hh"
#include "G4HadFinalState.hh"
#include "G4IonTable.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <cfloat>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

[[noreturn]] void failModel(const std::string& message) {
  G4Exception("MuonDISMuonNuclearModel", "MuonDISModelFailure", FatalException, message.c_str());
  throw std::runtime_error(message);
}

}  // namespace

MuonDISMuonNuclearModel::MuonDISMuonNuclearModel(const G4String& modelName, bool enableDebug)
    : G4HadronicInteraction(modelName), m_enableDebug(enableDebug) {
  SetMinEnergy(0.0);
  SetMaxEnergy(DBL_MAX);
}

MuonDISMuonNuclearModel::~MuonDISMuonNuclearModel() = default;

G4HadFinalState* MuonDISMuonNuclearModel::ApplyYourself(const G4HadProjectile& aTrack,
                                                        G4Nucleus& targetNucleus) {
  theParticleChange.Clear();

  const double muEnergyGeV = aTrack.GetTotalEnergy() / GeV;
  const G4ThreeVector muDirection = aTrack.Get4Momentum().vect().unit();
  const int muonPdgId = aTrack.GetDefinition()->GetPDGEncoding();

  const int targetZ = targetNucleus.GetZ_asInt();
  const int targetA = targetNucleus.GetA_asInt();

  MuonDISPythiaGenerator::GeneratedEvent generated;
  constexpr int kMaxAttempts = 5;
  for (int attempt = 0; attempt < kMaxAttempts && !generated.valid; ++attempt) {
    generated = MuonDISPythiaGenerator::instance().generate(muonPdgId, muEnergyGeV, muDirection,
                                                            targetZ, targetA);
  }
  if (!generated.valid) {
    failModel("MuonDISMuonNuclearModel: Pythia8 failed to produce a usable DIS final state "
              "after several attempts.");
  }

  theParticleChange.SetStatusChange(stopAndKill);
  theParticleChange.SetEnergyChange(0.0);
  theParticleChange.SetLocalEnergyDeposit(0.0);
  theParticleChange.SetMomentumChange(G4ThreeVector());

  auto* particleTable = G4ParticleTable::GetParticleTable();
  std::vector<int> disFinalPdgs;
  disFinalPdgs.reserve(generated.finalState.size());

  if (m_enableDebug) {
    G4cout << "MuonDISMuonNuclearModel: E_mu=" << muEnergyGeV << " GeV, Q2~" << generated.q2GeV2
           << " GeV^2, target nucleon pdg=" << generated.targetNucleonPdg
           << ", final-state count=" << generated.finalState.size() << G4endl;
  }

  for (const auto& particle : generated.finalState) {
    disFinalPdgs.push_back(particle.pdgId);
    auto* definition = particleTable->FindParticle(particle.pdgId);
    if (!definition) {
      if (m_enableDebug) {
        G4cout << "MuonDISMuonNuclearModel: skipping unknown PDG " << particle.pdgId << G4endl;
      }
      continue;
    }
    const G4ThreeVector momentum(particle.px * GeV, particle.py * GeV, particle.pz * GeV);
    theParticleChange.AddSecondary(new G4DynamicParticle(definition, momentum));
  }

  auto* recoilDefinition = G4IonTable::GetIonTable()->GetIon(targetZ, targetA, 0.0);
  if (recoilDefinition) {
    theParticleChange.AddSecondary(new G4DynamicParticle(recoilDefinition, G4ThreeVector()));
  }

  if (theParticleChange.GetNumberOfSecondaries() <= 0) {
    failModel("MuonDISMuonNuclearModel: no Geant4 secondaries were created.");
  }

  G4ThreeVector interactionPosition;
  double interactionGlobalTime = 0.0;
  int trackId = -1;
  int parentTrackId = -1;
  std::string volumeName;
  int geantEventId = -1;
  bool hasInteractionLocalPosition = false;
  G4ThreeVector interactionLocalPosition;
  const bool hasInteractionPosition =
      MuonDISInteractionRecorder::instance().currentInteraction(geantEventId,
                                                                interactionPosition,
                                                                interactionGlobalTime,
                                                                trackId,
                                                                parentTrackId,
                                                                volumeName,
                                                                hasInteractionLocalPosition,
                                                                interactionLocalPosition);

  MuonDISInteractionRecorder::instance().record({
      geantEventId,
      muEnergyGeV,
      generated.q2GeV2,
      generated.targetNucleonPdg,
      trackId,
      parentTrackId,
      hasInteractionPosition,
      interactionPosition,
      hasInteractionLocalPosition,
      interactionLocalPosition,
      interactionGlobalTime,
      volumeName,
      targetZ,
      targetA,
      static_cast<int>(generated.finalState.size()),
      theParticleChange.GetNumberOfSecondaries(),
      disFinalPdgs,
  });

  return &theParticleChange;
}

void MuonDISMuonNuclearModel::ModelDescription(std::ostream& outFile) const {
  outFile << "MuonDIS muon-nuclear interaction model backed by on-the-fly Pythia8 "
             "lepton-nucleon DIS (WeakBosonExchange:ff2ff(t:gmZ)).";
}

const std::pair<G4double, G4double> MuonDISMuonNuclearModel::GetFatalEnergyCheckLevels() const {
  // Loosened from the G4 default: the recoil nucleus is added at rest (no
  // nuclear-medium/Fermi-motion recoil calculation), so exact energy-momentum
  // conservation against the *nuclear* target is only approximate even though
  // Pythia8 conserves it exactly against the free nucleon it was given.
  return {1e6, 10000 * GeV};
}
