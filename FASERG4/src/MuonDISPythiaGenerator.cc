#include "MuonDISPythiaGenerator.hh"

#include <Pythia8/Pythia.h>

#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "globals.hh"

#include <cmath>
#include <stdexcept>

namespace {
// PDG muon mass, GeV (matches G4MuonMinus/G4MuonPlus to good precision; not
// pulled from Geant4 directly since this header avoids depending on
// G4ParticleDefinition).
constexpr double kMuonMassGeV = 0.1056583755;
}  // namespace

MuonDISPythiaGenerator& MuonDISPythiaGenerator::instance() {
  static MuonDISPythiaGenerator generator;
  return generator;
}

MuonDISPythiaGenerator::~MuonDISPythiaGenerator() = default;

void MuonDISPythiaGenerator::configure(bool enableDebug, double q2MinGeV2) {
  m_enableDebug = enableDebug;
  m_q2MinGeV2 = q2MinGeV2;
}

int MuonDISPythiaGenerator::chooseNucleon(int targetZ, int targetA) const {
  if (targetA <= 0) {
    return 2212;  // fall back to a free proton if target bookkeeping is missing
  }
  const double protonFraction = static_cast<double>(targetZ) / static_cast<double>(targetA);
  return (G4UniformRand() < protonFraction) ? 2212 : 2112;
}

void MuonDISPythiaGenerator::ensureInitialized(int muonPdgId, int nucleonPdgId,
                                               double muonEnergyGeV,
                                               const G4ThreeVector& muonDirection) {
  if (!m_pythia) {
    m_pythia = std::make_unique<Pythia8::Pythia>();
  }

  if (!m_settingsApplied) {
    // Keep Pythia8's own banner/summary output out of the way; MuonDISPhysics
    // already logs what it installs, and per-event listings would be far too
    // verbose for a production run.
    m_pythia->readString("Print:quiet = true");
    m_pythia->readString("Init:showProcesses = false");
    m_pythia->readString("Init:showMultipartonInteractions = false");
    m_pythia->readString("Init:showChangedSettings = false");
    m_pythia->readString("Init:showChangedParticleData = false");
    m_pythia->readString("Next:numberCount = 0");
    m_pythia->readString("Next:numberShowInfo = 0");
    m_pythia->readString("Next:numberShowProcess = 0");
    m_pythia->readString("Next:numberShowEvent = 0");

    // Neutral-current t-channel gamma*/Z exchange: the standard Pythia8 route
    // to lepton-nucleon DIS (see MuonDISPythiaGenerator.hh for the caveats
    // and the charged-current toggle).
    m_pythia->readString("WeakBosonExchange:ff2ff(t:gmZ) = on");

    // Recommended scale choice for this class of t-channel process (Pythia8
    // manual, Electroweak Processes).
    m_pythia->readString("SigmaProcess:factorScale2 = 6");
    m_pythia->readString("SigmaProcess:renormScale2 = 6");

    // Arbitrary 3-momentum beams (Pythia8 manual, Beam Parameters, option 3)
    // so beam A can track the true Geant4 muon direction every interaction.
    m_pythia->readString("Beams:frameType = 3");

    m_settingsApplied = true;
  }

  const double muonMomentumGeV =
      std::sqrt(std::max(muonEnergyGeV * muonEnergyGeV - kMuonMassGeV * kMuonMassGeV, 0.0));
  const G4ThreeVector p = muonDirection.unit() * muonMomentumGeV;

  m_pythia->settings.mode("Beams:idA", muonPdgId);
  m_pythia->settings.mode("Beams:idB", nucleonPdgId);
  m_pythia->settings.parm("Beams:pxA", p.x());
  m_pythia->settings.parm("Beams:pyA", p.y());
  m_pythia->settings.parm("Beams:pzA", p.z());
  // Beam B (the struck nucleon) at rest in this lab frame -- i.e. this frame
  // *is* the target-nucleon rest frame, consistent with a fixed-target setup.
  m_pythia->settings.parm("Beams:pxB", 0.0);
  m_pythia->settings.parm("Beams:pyB", 0.0);
  m_pythia->settings.parm("Beams:pzB", 0.0);
  m_pythia->settings.parm("PhaseSpace:Q2Min", m_q2MinGeV2);

  // Beam species/kinematics are latched at init() time; since the muon
  // direction and energy differ essentially every call here, we simply
  // re-init() every time (see header comment on setKinematics as a possible
  // optimization once this is running).
  const bool ok = m_pythia->init();
  if (!ok) {
    throw std::runtime_error(
        "MuonDISPythiaGenerator: Pythia8 init() failed for the requested muon-nucleon DIS beam "
        "configuration (E_mu=" + std::to_string(muonEnergyGeV) + " GeV, nucleon pdg=" +
        std::to_string(nucleonPdgId) + ").");
  }
  m_lastMuonPdg = muonPdgId;
  m_lastNucleonPdg = nucleonPdgId;
}

MuonDISPythiaGenerator::GeneratedEvent MuonDISPythiaGenerator::generate(
    int muonPdgId, double muonEnergyGeV, const G4ThreeVector& muonDirection, int targetZ,
    int targetA) {
  GeneratedEvent result;

  const int nucleonPdgId = chooseNucleon(targetZ, targetA);

  try {
    ensureInitialized(muonPdgId, nucleonPdgId, muonEnergyGeV, muonDirection);
  } catch (const std::exception& error) {
    G4cout << "MuonDISPythiaGenerator: " << error.what() << G4endl;
    return result;
  }

  if (!m_pythia->next()) {
    if (m_enableDebug) {
      G4cout << "MuonDISPythiaGenerator: Pythia8 next() failed for this interaction, retrying "
                "once." << G4endl;
    }
    if (!m_pythia->next()) {
      return result;  // caller decides how many times to retry / whether to abort
    }
  }

  result.targetNucleonPdg = nucleonPdgId;
  // -tHat() is the standard proxy for the (spacelike) momentum-transfer Q^2
  // of a t-channel 2->2 process; this is a diagnostic field only and does not
  // feed back into the injected kinematics below.
  result.q2GeV2 = -m_pythia->info.tHat();

  for (int i = 0; i < m_pythia->event.size(); ++i) {
    const Pythia8::Particle& particle = m_pythia->event[i];
    if (!particle.isFinal()) {
      continue;
    }
    GeneratedParticle gp;
    gp.pdgId = particle.id();
    gp.px = particle.px();
    gp.py = particle.py();
    gp.pz = particle.pz();
    gp.e = particle.e();
    result.finalState.push_back(gp);
  }

  result.valid = !result.finalState.empty();
  return result;
}
