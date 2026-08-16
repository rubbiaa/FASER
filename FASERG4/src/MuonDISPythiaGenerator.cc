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

// Average isoscalar (p/n) nucleon mass, GeV -- MUST match TPOEvent::kinematics_event()'s own
// kNucleonMassGeV, since observedXBj() below exists specifically to compute the same quantity
// TPOEvent will independently recompute downstream, for a consistent /physics/muondis/xbjmin cut.
constexpr double kNucleonMassGeV = 0.93892;

// /physics/muondis/xbjmin is enforced by re-drawing next() (cheap -- the beam is already
// init()'d) until an event clears the cut; there is no cap on that search itself (see
// generate()). This cap instead guards against a genuinely broken/exhausted beam
// configuration where Pythia8's next() itself keeps failing outright (unrelated to the xbj
// value of any generated event) -- without it, that pathological case would spin forever.
constexpr int kMaxConsecutiveNextFailures = 25;

// Bjorken x as it will actually be RECONSTRUCTED downstream from the observable final state
// (exactly the formula TPOEvent::kinematics_event() applies to in_neutrino/out_lepton), i.e.
// AFTER Pythia8's full parton shower/hadronization -- NOT Pythia8's own Info::x2(), which is
// the momentum fraction at the hard 2->2 vertex, BEFORE any subsequent final-state QED/QCD
// radiation off the outgoing muon leg. That radiation lowers the observed outgoing-muon energy,
// which raises the observed nu and therefore LOWERS the observed x relative to Info::x2() -- so
// enforcing /physics/muondis/xbjmin against Info::x2() alone does not guarantee the truth-level
// TPOEvent::xBj actually shown/analyzed downstream clears the cut (confirmed: every accepted
// event's Info::x2() was comfortably above the configured cutoff, while the corresponding
// TPOEvent::xBj was several times smaller and often below it). Returns -1 if no final-state
// particle with the same PDG code as the incoming muon (event[1], per the frameType=3 beam-A
// convention set in ensureInitialized()) can be found, or if the reconstructed nu is not
// positive (kinematically degenerate/unphysical draw) -- both treated by the caller as "does
// not clear the cut", triggering another next() draw.
double observedXBj(const Pythia8::Event& event) {
  if (event.size() < 2) {
    return -1.0;
  }
  const Pythia8::Particle& muonIn = event[1];
  double bestE = -1.0;
  int bestIndex = -1;
  for (int i = 2; i < event.size(); ++i) {
    const Pythia8::Particle& particle = event[i];
    if (!particle.isFinal() || particle.id() != muonIn.id()) {
      continue;
    }
    if (particle.e() > bestE) {
      bestE = particle.e();
      bestIndex = i;
    }
  }
  if (bestIndex < 0) {
    return -1.0;
  }
  const Pythia8::Particle& muonOut = event[bestIndex];
  const double nu = muonIn.e() - muonOut.e();
  if (nu <= 0.0) {
    return -1.0;
  }
  const double dpx = muonIn.px() - muonOut.px();
  const double dpy = muonIn.py() - muonOut.py();
  const double dpz = muonIn.pz() - muonOut.pz();
  const double q2 = dpx * dpx + dpy * dpy + dpz * dpz - nu * nu;
  return q2 / (2.0 * kNucleonMassGeV * nu);
}
}  // namespace

MuonDISPythiaGenerator& MuonDISPythiaGenerator::instance() {
  static MuonDISPythiaGenerator generator;
  return generator;
}

MuonDISPythiaGenerator::~MuonDISPythiaGenerator() = default;

void MuonDISPythiaGenerator::configure(bool enableDebug, double q2MinGeV2,
                                       const std::string& pdfSetPath, double xBjMin) {
  m_enableDebug = enableDebug;
  m_q2MinGeV2 = q2MinGeV2;
  m_pdfSetPath = pdfSetPath;
  m_xBjMin = xBjMin;
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

    // Optional non-default PDF for the struck nucleon (beam B), e.g. one of the grids bundled
    // for the muDIS charm-asymmetry study (see MuonDISPythiaGenerator.hh and
    // README_MuonDIS.md). Loaded natively via Pythia8's built-in LHAGrid1 reader -- no
    // external LHAPDF6 install needed. Left at Pythia8's own built-in proton PDF if empty
    // (the default, set via /physics/muondis/pdfSet if you want to change it).
    if (!m_pdfSetPath.empty()) {
      m_pythia->readString("PDF:pSet = LHAGrid1:" + m_pdfSetPath);
    }

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

  // Enforce the minimum Bjorken-x cut (/physics/muondis/xbjmin), if configured, against the
  // RECONSTRUCTED (post-shower) x -- see observedXBj() above for why this must NOT be
  // Info::x2() (the pre-shower, hard-vertex value): the two can differ by several-fold once
  // final-state radiation off the outgoing muon is included, and it's the reconstructed value
  // that TPOEvent::xBj (what's actually displayed/analyzed downstream) will match. Re-generating
  // via next() on the already-init()'d beam is cheap (unlike ensureInitialized(), which re-inits
  // Pythia8's full beam setup), so retry the cut here rather than pushing every attempt out to
  // the caller's much coarser retry loop.
  if (m_xBjMin > 0.0) {
    int xBjAttempt = 0;
    int consecutiveNextFailures = 0;
    while (observedXBj(m_pythia->event) < m_xBjMin) {
      ++xBjAttempt;
      if (!m_pythia->next()) {
        ++consecutiveNextFailures;
        if (consecutiveNextFailures >= kMaxConsecutiveNextFailures) {
          if (m_enableDebug) {
            G4cout << "MuonDISPythiaGenerator: giving up after " << consecutiveNextFailures
                   << " consecutive Pythia8 next() failures while hunting for an event with "
                      "xBj >= " << m_xBjMin << " (" << xBjAttempt << " draw(s) total)."
                   << G4endl;
          }
          return result;  // caller treats this like any other failed-to-generate attempt
        }
        continue;  // next() itself failed -- retry, does not count against the xBj search
      }
      consecutiveNextFailures = 0;
      if (m_enableDebug && xBjAttempt % 1000 == 0) {
        G4cout << "MuonDISPythiaGenerator: still hunting for xBj >= " << m_xBjMin << " after "
               << xBjAttempt << " draw(s) (most recent reconstructed x="
               << observedXBj(m_pythia->event) << ")." << G4endl;
      }
    }
  }

  result.targetNucleonPdg = nucleonPdgId;
  // -tHat() is the standard proxy for the (spacelike) momentum-transfer Q^2
  // of a t-channel 2->2 process; this is a diagnostic field only and does not
  // feed back into the injected kinematics below.
  result.q2GeV2 = -m_pythia->info.tHat();
  // Record the same RECONSTRUCTED (post-shower) Bjorken x used to enforce /physics/muondis/
  // xbjmin above (always populated, whether or not the cut is active), so callers get a value
  // that should closely match the independently-recomputed truth-level TPOEvent::xBj -- unlike
  // Info::x2(), which is the pre-shower hard-vertex value and can differ substantially (see
  // observedXBj() above).
  result.pythiaX2 = observedXBj(m_pythia->event);

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
