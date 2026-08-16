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
#include <string>
#include <vector>

namespace {

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
                              const G4ThreeVector& interactionPosition,
                              const std::string& volumeName,
                              double muEnergyGeV, const G4ThreeVector& muMomentumGeV) {
  // For MuonDIS events only: record the actual Info::x2() Pythia8 used to accept this event, so
  // it can be directly cross-checked against the truth-level TPOEvent::xBj recomputed below by
  // kinematics_event() from in_neutrino/out_lepton -- a separate calculation that should agree
  // closely but is not literally the same quantity.
  TPOEvent* tpoEvent = currentTPOEvent();
  if (!tpoEvent) {
    return;
  }

  int parentPOTrackId = muonTrackId;
  PO* parentPO = nullptr;
  for (auto& candidate : tpoEvent->POs) {
    if (candidate.m_track_id == muonTrackId) {
      parentPO = &candidate;
      break;
    }
  }
  if (!parentPO && !tpoEvent->POs.empty()) {
    // Fall back to the first PO (the primary muon, in the muon-background running mode this
    // extension currently targets) if the track ID recorded by MuonDISNuclearWrapperProcess
    // could not be matched -- better to link to *something* sensible than to silently drop
    // the parent reference.
    parentPO = &tpoEvent->POs.front();
  }
  if (parentPO) {
    parentPOTrackId = parentPO->m_track_id;

    // The incoming muon PO was filled once at event-generation time (PrimaryGeneratorAction),
    // recording the muon's momentum where it was injected in front of the detector. By the
    // time this DIS interaction actually fires, the real G4 track has typically lost some
    // energy (ionization etc. along the flight path) -- TPOEvent::kinematics_event()'s
    // nu/Q2/x/y only make sense if in_neutrino and the generated final state are evaluated
    // at the same instant, so overwrite the muon PO with its actual momentum right at the
    // interaction vertex: the same 4-momentum just handed to Pythia8 as beam A. Without this,
    // in_neutrino stayed at the launch-time value while out_lepton/jet came from a final
    // state balanced against the (lower-energy, possibly redirected) vertex momentum, which
    // could make Q2 inconsistent with nu/W2 (even yielding an unphysical negative W2).
    parentPO->m_px = muMomentumGeV.x();
    parentPO->m_py = muMomentumGeV.y();
    parentPO->m_pz = muMomentumGeV.z();
    parentPO->m_energy = muEnergyGeV;
    parentPO->m_kinetic_energy = muEnergyGeV - parentPO->m_mass();
  }

  const double vx = hasInteractionPosition ? interactionPosition.x() / mm : 0.0;
  const double vy = hasInteractionPosition ? interactionPosition.y() / mm : 0.0;
  const double vz = hasInteractionPosition ? interactionPosition.z() / mm : 0.0;

  // The primary vertex PrimaryGeneratorAction set at the start of the event is just where the
  // muon was injected in front of the detector, not where the physics interaction actually
  // happened. Move it to the real DIS interaction point, the same way a GENIE-derived event's
  // prim_vx records where the neutrino actually interacted rather than where it originated.
  if (hasInteractionPosition) {
    tpoEvent->setPrimaryVtx(vx, vy, vz);
    if (volumeName == "targetW") {
      tpoEvent->setVtxTarget(TPOEvent::kVtx_in_W);
    } else if (volumeName == "Scintillator") {
      tpoEvent->setVtxTarget(TPOEvent::kVtx_in_Scint);
    }
  }

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
  tpoEvent->pythiaXbj = generated.pythiaX2;
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
    // Pythia8 occasionally fails to init()/generate a usable DIS final state for some
    // (E_mu, nucleon) combinations -- observed at E_mu ~ 2.1 TeV, well within the real
    // FASERnu muon flux's range (the extreme lab-frame asymmetry between a multi-TeV muon
    // beam and an at-rest nucleon apparently trips up Pythia8's arbitrary-momentum-beam
    // (frameType=3) initialization for some configurations). This is not a programming
    // error worth aborting an entire production run over -- this used to call
    // G4Exception(..., FatalException, ...), which killed the whole job the first time any
    // muon in the sampled flux hit such an energy. Instead, treat this specific biased step
    // as a no-op: the muon survives with its kinematics unchanged, exactly as
    // MuonDISNuclearWrapperProcess::demoteSurvivingPrimary() already expects for a step
    // where the biased process didn't actually produce a hard interaction -- it converts
    // the surviving primary back to an ordinary (unbiased) muon and tracking continues
    // normally, so the muon simply gets another chance to undergo MuonDIS further along its
    // path.
    static G4int failureCount = 0;
    ++failureCount;
    G4cout << "MuonDISMuonNuclearModel: Pythia8 failed to produce a usable DIS final state for "
              "E_mu=" << muEnergyGeV << " GeV after " << kMaxAttempts << " attempts (failure #"
           << failureCount << "); muon survives this step unchanged." << G4endl;
    theParticleChange.Clear();
    theParticleChange.SetStatusChange(isAlive);
    theParticleChange.SetEnergyChange(aTrack.GetKineticEnergy());
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
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
    // Same reasoning as the Pythia8-failure fallback above: a valid-looking generated final
    // state whose PDG codes G4ParticleTable/G4IonTable could not resolve into any real
    // secondary is a rare generator/table mismatch, not something worth crashing a production
    // run over. Undo the stopAndKill/zero-momentum commitment made above and let the muon
    // survive this step unchanged instead.
    static G4int noSecondaryFailureCount = 0;
    ++noSecondaryFailureCount;
    G4cout << "MuonDISMuonNuclearModel: generated DIS final state produced no resolvable Geant4 "
              "secondaries (failure #" << noSecondaryFailureCount
           << "); muon survives this step unchanged." << G4endl;
    theParticleChange.Clear();
    theParticleChange.SetStatusChange(isAlive);
    theParticleChange.SetEnergyChange(aTrack.GetKineticEnergy());
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }

  // Record the full Pythia8 final state as generator-level truth, the same way
  // PrimaryGeneratorAction::GeneratePrimaries records the incoming primary muon itself. This
  // covers every particle Pythia8 reported, even any whose PDG code G4ParticleTable could not
  // resolve above (and which therefore did not become a real G4 secondary) -- truth-level POs
  // are meant to capture the full generator output, not just what G4 could instantiate.
  const G4ThreeVector muMomentumGeV = aTrack.Get4Momentum().vect() / GeV;
  recordDISFinalStateTruth(generated, trackId, hasInteractionPosition, interactionPosition, volumeName,
                           muEnergyGeV, muMomentumGeV);

  MuonDISInteractionRecorder::instance().record({
      geantEventId,
      muEnergyGeV,
      generated.q2GeV2,
      generated.pythiaX2,
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
