#include "MuonDISPhysics.hh"

#include "MuonDISInteractionRecorder.hh"
#include "MuonDISMessenger.hh"
#include "MuonDISMuonNuclearModel.hh"
#include "MuonDISNuclearWrapperProcess.hh"
#include "MuonDISPythiaGenerator.hh"

#include "G4MuonMinus.hh"
#include "G4MuonNuclearProcess.hh"
#include "G4MuonPlus.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "globals.hh"

// factory registration, mirroring TauDecayPhysics in this directory
#include "G4PhysicsConstructorFactory.hh"
G4_DECLARE_PHYSCONSTR_FACTORY(MuonDISPhysics);

#include <cmath>
#include <exception>

MuonDISPhysics::MuonDISPhysics(G4int) : G4VPhysicsConstructor("MuonDISPhysics") {
  fMessenger = new MuonDISMessenger(this);
}

MuonDISPhysics::~MuonDISPhysics() {
  delete fMessenger;
}

void MuonDISPhysics::ConstructParticle() {
  // Nothing to do here -- standard mu-/mu+ already exist from the base physics list.
}

void MuonDISPhysics::ConstructProcess() {
  if (!fEnabled) {
    G4cout << "MuonDISPhysics: disabled (enable with /physics/muondis/enable true BEFORE "
              "/run/initialize). Standard Geant4 muon-nuclear process left untouched."
           << G4endl;
    return;
  }

  try {
    MuonDISInteractionRecorder::instance().configure(fInteractionLogPath);
  } catch (const std::exception& error) {
    G4Exception("MuonDISPhysics", "MuonDISInteractionLogFailure", FatalException, error.what());
    return;
  }

  MuonDISPythiaGenerator::instance().configure(fEnableDebug, fQ2Min, fPdfSetPath, fXBjMin);

  auto patchParticle = [&](G4ParticleDefinition* particle) {
    if (!particle) {
      return;
    }

    auto* manager = particle->GetProcessManager();
    if (!manager) {
      G4cout << "MuonDISPhysics: no process manager for " << particle->GetParticleName()
             << G4endl;
      return;
    }

    G4MuonNuclearProcess* muonNuclear = nullptr;
    G4ProcessVector* processList = manager->GetProcessList();
    for (size_t index = 0; index < processList->length(); ++index) {
      auto* process = (*processList)[index];
      muonNuclear = dynamic_cast<G4MuonNuclearProcess*>(process);
      if (muonNuclear) {
        break;
      }
    }

    if (!muonNuclear) {
      G4cout << "MuonDISPhysics: failed to find G4MuonNuclearProcess for "
             << particle->GetParticleName() << G4endl;
      return;
    }

    manager->RemoveProcess(muonNuclear);

    auto* unbiasedProcess = muonNuclear;
    auto* biasedProcess = new G4MuonNuclearProcess();
    auto* model = new MuonDISMuonNuclearModel(
        "MuonDISMuonNuclearModel_" + particle->GetParticleName(), fEnableDebug);
    biasedProcess->RegisterMe(model);

    if (fCrossSectionBias > 0.0 && std::abs(fCrossSectionBias - 1.0) > 1e-12) {
      biasedProcess->BiasCrossSectionByFactor(fCrossSectionBias);
    }

    auto* replacement = new MuonDISNuclearWrapperProcess(
        "MuonDISNuclearWrapperProcess_" + particle->GetParticleName(), biasedProcess,
        unbiasedProcess, fEnableDebug);
    manager->AddDiscreteProcess(replacement);

    G4cout << "MuonDISPhysics: installed Pythia8-driven MuonDIS process for "
           << particle->GetParticleName() << " with CrossSectionBias=" << fCrossSectionBias
           << " (primary muons only; secondaries and post-DIS outgoing muons remain unbiased)."
           << G4endl;
  };

  patchParticle(G4MuonMinus::MuonMinus());
  patchParticle(G4MuonPlus::MuonPlus());
}
