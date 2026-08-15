#include "MuonDISMuonNuclearModel.hh"

#include "MuonDISInteractionRecorder.hh"
#include "MuonDISPythiaGenerator.hh"

#include "G4DynamicParticle.hh"
#include "G4HadFinalState.hh"
#include "G4IonTable.hh"
#include "G4ParticleTable.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "PrimaryGeneratorAction.hh"
#include "TPOEvent.hh"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

[[noreturn]] void failModel(const std::string& message) {
  G4Exception("MuonDISMuonNuclearModel", "MuonDISModelFailure", FatalException, message.c_str());
  throw std::runtime_error(message);
}

// Obtain write access to the current event's TPOEvent truth structure. TPOEvent is owned by
// PrimaryGeneratorAction (one instance per worker thread) and GetTPOEvent() only exposes a
// const accessor everywhere else in this codebase; ParticleManager.cc already works around
// this with the same const_cast, which we reuse here rather than inventing a new accessor.
TPOEvent* currentTPOEvent() {
  auto* runManager = G4RunManager::GetRunManager();
  if (!runManager) {
    return nullptr;
  }
  const auto* primaryGenAction =
      dynamic_cast<const PrimaryGeneratorAction*>(runManager->GetUserPrimaryGeneratorAction());
  if (!primaryGenAction) {
    return nullptr;
  }
  return const_cast<TPOEvent*>(primaryGenAction->GetTPOEvent());
}

// Append the Pythia8 DIS final-state particles to the truth PO list -- the same TPOEvent.POs
// vector that PrimaryGeneratorAction::GeneratePrimaries fills for the incoming primary muon.
// Each new PO is linked back to the muon's own PO entry via nparent/m_trackid_in_particle[0],
// following the parent-bookkeeping convention already used for tau/charm decay products in
// TPOEvent::perform_charmhadron_decay(). geanttrackID is left at -1: unlike the primary muon
// (whose G4 track ID is known up front, from PrimaryGeneratorAction), these particles don't
// have a GEANT4 track ID yet -- they only become real G4 tracks after being pushed via
// G4ParticleChange::AddSecondary() above, and G4 doesn't assign their track IDs until they are
// popped off the tracking stack. Matching them back to these PO entries (PDG + momentum, the
// way TcalEvent::AssignGEANTTrackID already does for primary-vertex particles) would need a
// similar hook at tracking time -- e.g. extending ParticleManager::RecordTrack() /
// TrackingAction::PreUserTrackingAction() -- which is not yet wired up for MuonDIS secondaries.
void recordDISFinalStateTruth(const MuonDISPythiaGenerator::GeneratedEvent& generated,
                              int muonTrackId, bool hasInteractionPosition,
                              const G4ThreeVector& interactionPosition) {
  TPOEvent* tpoEvent = currentTPOEvent();
  if (!tpoEvent) {
    return;
  }

  int parentPOTrackId = muonTrackId;
  bool foundParent = false;
  for (const auto& candidate : tpoEvent->POs) {
    if (candidate.m_track_id == muonTrackId) {
      foundParent = true;
      break;
    }
  }
  if (!foundParent && !tpoEvent->POs.empty()) {
    // Fall back to the first PO (the primary muon, in the muon-background running mode this
    // extension currently targets) if the track ID recorded by MuonDISNuclearWrapperProcess
    // could not be matched -- better to link to *something* sensible than to silently drop
    // the parent reference.
    parentPOTrackId = tpoEvent->POs.front().m_track_id;
  }

  const double vx = hasInteractionPosition ? interactionPosition.x() / mm : 0.0;
  const double vy = hasInteractionPosition ? interactionPosition.y() / mm : 0.0;
  const double vz = hasInteractionPosition ? interactionPosition.z() / mm : 0.0;

  for (const auto& particle : generated.finalState) {
    const double mass2 = particle.e * particle.e - particle.px * particle.px -
                         particle.py * particle.py - particle.pz * particle.pz;
    const double mass = std::sqrt(std::max(mass2, 0.0));

    PO aPO{};
    aPO.m_pdg_id = particle.pdgId;
    aPO.m_track_id = static_cast<int>(tpoEvent->POs.size()) + 1;
    aPO.m_status = 1;
    aPO.m_px = particle.px;
    aPO.m_py = particle.py;
    aPO.m_pz = particle.pz;
    aPO.m_energy = particle.e;
    aPO.m_kinetic_energy = particle.e - mass;
    aPO.m_vx_decay = vx;
    aPO.m_vy_decay = vy;
    aPO.m_vz_decay = vz;
    aPO.nparent = 1;
    aPO.m_trackid_in_particle[0] = parentPOTrackId;
    aPO.geanttrackID = -1;
    tpoEvent->POs.push_back(aPO);
  }

  // Recompute in_neutrino/out_lepton/jet/Evis/ptmiss/DIS-kinematics now that the DIS final
  // state has been appended -- otherwise these stay at whatever kinematics_event() last
  // computed (or their cleared defaults) until something else happens to call it again, and
  // PrimaryGeneratorAction::GeneratePrimaries's dump_event() call for the *next* event would
  // print stale/default values for this one.
  tpoEvent->kinematics_event();
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

  // Fetch the interaction context recorded by MuonDISNuclearWrapperProcess (the muon's own G4
  // track/parent IDs, interaction position/time/volume) up front: it is needed both to link the
  // truth PO entries below back to the muon's PO and for the CSV interaction log further down,
  // so we read it once and reuse it for both rather than calling currentInteraction() twice.
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

  // Record the full Pythia8 final state as generator-level truth, the same way
  // PrimaryGeneratorAction::GeneratePrimaries records the incoming primary muon itself. This
  // covers every particle Pythia8 reported, even any whose PDG code G4ParticleTable could not
  // resolve above (and which therefore did not become a real G4 secondary) -- truth-level POs
  // are meant to capture the full generator output, not just what G4 could instantiate.
  recordDISFinalStateTruth(generated, trackId, hasInteractionPosition, interactionPosition);

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
      static_cast<int>(theParticleChange.GetNumberOfSecondaries()),
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
