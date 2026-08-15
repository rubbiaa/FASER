#include "MuonDISNuclearWrapperProcess.hh"

#include "MuonDISInteractionRecorder.hh"

#include "G4AffineTransform.hh"
#include "G4DynamicParticle.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4MuonMinus.hh"
#include "G4MuonNuclearProcess.hh"
#include "G4MuonPlus.hh"
#include "G4ParticleChange.hh"
#include "G4ProcessManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TouchableHistory.hh"
#include "globals.hh"

#include <stdexcept>
#include <string>

namespace {

const G4ParticleDefinition* standardMuonDefinition(const G4Track& track) {
  return (track.GetDefinition()->GetPDGCharge() < 0.0)
             ? static_cast<const G4ParticleDefinition*>(G4MuonMinus::MuonMinus())
             : static_cast<const G4ParticleDefinition*>(G4MuonPlus::MuonPlus());
}

class InteractionScope {
public:
  InteractionScope(int geantEventId,
                   const G4ThreeVector& position, double globalTime,
                   int trackId, int parentTrackId, const std::string& volumeName,
                   bool hasLocalPosition, const G4ThreeVector& localPosition)
      : m_recorder(MuonDISInteractionRecorder::instance()) {
    m_recorder.setCurrentInteraction(geantEventId, position, globalTime,
                                     trackId, parentTrackId, volumeName,
                                     hasLocalPosition, localPosition);
  }

  ~InteractionScope() {
    m_recorder.clearCurrentInteraction();
  }

private:
  MuonDISInteractionRecorder& m_recorder;
};

int currentGeantEventId() {
  auto* eventManager = G4EventManager::GetEventManager();
  if (!eventManager) {
    return -1;
  }
  const auto* currentEvent = eventManager->GetConstCurrentEvent();
  return currentEvent ? currentEvent->GetEventID() : -1;
}

bool localInteractionPosition(const G4Step& stepData, G4ThreeVector& localPosition) {
  const auto* touch = dynamic_cast<const G4TouchableHistory*>(
      stepData.GetPreStepPoint()->GetTouchable());
  if (!touch) {
    return false;
  }
  const G4AffineTransform transformation = touch->GetHistory()->GetTopTransform();
  localPosition = transformation.TransformPoint(stepData.GetPostStepPoint()->GetPosition());
  return true;
}

}  // namespace

MuonDISNuclearWrapperProcess::MuonDISNuclearWrapperProcess(const G4String& name,
                                                           G4MuonNuclearProcess* biasedProcess,
                                                           G4MuonNuclearProcess* unbiasedProcess,
                                                           bool enableDebug)
    : G4WrapperProcess(name),
      m_biasedProcess(biasedProcess),
      m_unbiasedProcess(unbiasedProcess),
      m_enableDebug(enableDebug) {
  if (!m_biasedProcess || !m_unbiasedProcess) {
    throw std::runtime_error(
        "MuonDISNuclearWrapperProcess requires both biased and unbiased processes.");
  }

  SetProcessType(m_unbiasedProcess->GetProcessType());
  SetProcessSubType(m_unbiasedProcess->GetProcessSubType());
}

MuonDISNuclearWrapperProcess::~MuonDISNuclearWrapperProcess() {
  // Do NOT delete m_biasedProcess / m_unbiasedProcess here. Every G4VProcess
  // object (including these two, and this wrapper itself) self-registers
  // into Geant4's global/thread-local G4ProcessTable when constructed, and
  // G4ProcessTable::Clear() (invoked from G4TaskRunManager's destructor at
  // shutdown) independently deletes every process it tracks -- including
  // m_biasedProcess and m_unbiasedProcess directly, separately from this
  // wrapper. Explicitly deleting them here as well raced with that and
  // produced a double-free segfault inside G4ProcessTable::~G4ProcessTable()
  // at program exit (observed 2026-08, ryzen01 build). G4ProcessTable owns
  // their cleanup; we just drop our raw pointers.
  m_biasedProcess = nullptr;
  m_unbiasedProcess = nullptr;
}

G4double MuonDISNuclearWrapperProcess::PostStepGetPhysicalInteractionLength(
    const G4Track& track, G4double previousStepSize, G4ForceCondition* condition) {
  return selectProcess(track)->PostStepGetPhysicalInteractionLength(track, previousStepSize,
                                                                    condition);
}

G4VParticleChange* MuonDISNuclearWrapperProcess::PostStepDoIt(const G4Track& track,
                                                              const G4Step& stepData) {
  if (!useBiasedProcess(track)) {
    return m_unbiasedProcess->PostStepDoIt(track, stepData);
  }

  const auto& post = *stepData.GetPostStepPoint();
  const auto* volume = post.GetPhysicalVolume();
  G4ThreeVector localPosition;
  const bool hasLocalPosition = localInteractionPosition(stepData, localPosition);
  InteractionScope interactionScope(currentGeantEventId(),
                                    post.GetPosition(), post.GetGlobalTime(),
                                    track.GetTrackID(), track.GetParentID(),
                                    volume ? volume->GetName() : "",
                                    hasLocalPosition, localPosition);
  auto* particleChange = m_biasedProcess->PostStepDoIt(track, stepData);

  demoteSurvivingPrimary(track, stepData, *particleChange);
  return particleChange;
}

G4bool MuonDISNuclearWrapperProcess::IsApplicable(const G4ParticleDefinition& particle) {
  return m_unbiasedProcess->IsApplicable(particle);
}

void MuonDISNuclearWrapperProcess::BuildPhysicsTable(const G4ParticleDefinition& particle) {
  m_unbiasedProcess->BuildPhysicsTable(particle);
  m_biasedProcess->BuildPhysicsTable(particle);
}

void MuonDISNuclearWrapperProcess::PreparePhysicsTable(const G4ParticleDefinition& particle) {
  m_unbiasedProcess->PreparePhysicsTable(particle);
  m_biasedProcess->PreparePhysicsTable(particle);
}

G4bool MuonDISNuclearWrapperProcess::StorePhysicsTable(const G4ParticleDefinition* particle,
                                                       const G4String& directory, G4bool ascii) {
  return m_unbiasedProcess->StorePhysicsTable(particle, directory, ascii) &&
         m_biasedProcess->StorePhysicsTable(particle, directory, ascii);
}

G4bool MuonDISNuclearWrapperProcess::RetrievePhysicsTable(const G4ParticleDefinition* particle,
                                                          const G4String& directory,
                                                          G4bool ascii) {
  return m_unbiasedProcess->RetrievePhysicsTable(particle, directory, ascii) &&
         m_biasedProcess->RetrievePhysicsTable(particle, directory, ascii);
}

void MuonDISNuclearWrapperProcess::StartTracking(G4Track* track) {
  selectProcess(*track)->StartTracking(track);
}

void MuonDISNuclearWrapperProcess::EndTracking() {
  m_unbiasedProcess->EndTracking();
  m_biasedProcess->EndTracking();
}

void MuonDISNuclearWrapperProcess::SetProcessManager(const G4ProcessManager* procMan) {
  m_processManager = procMan;
  m_unbiasedProcess->SetProcessManager(procMan);
  m_biasedProcess->SetProcessManager(procMan);
}

const G4ProcessManager* MuonDISNuclearWrapperProcess::GetProcessManager() {
  return m_processManager;
}

void MuonDISNuclearWrapperProcess::ResetNumberOfInteractionLengthLeft() {
  m_unbiasedProcess->ResetNumberOfInteractionLengthLeft();
  m_biasedProcess->ResetNumberOfInteractionLengthLeft();
}

void MuonDISNuclearWrapperProcess::SetMasterProcess(G4VProcess* masterP) {
  m_unbiasedProcess->SetMasterProcess(masterP);
  m_biasedProcess->SetMasterProcess(masterP);
}

G4bool MuonDISNuclearWrapperProcess::useBiasedProcess(const G4Track& track) const {
  return track.GetParentID() == 0;
}

G4MuonNuclearProcess* MuonDISNuclearWrapperProcess::selectProcess(const G4Track& track) const {
  return useBiasedProcess(track) ? m_biasedProcess : m_unbiasedProcess;
}

void MuonDISNuclearWrapperProcess::demoteSurvivingPrimary(const G4Track& track,
                                                          const G4Step& stepData,
                                                          G4VParticleChange& particleChange) const {
  if (particleChange.GetTrackStatus() == fStopAndKill) {
    return;
  }

  const auto* stdMuon = standardMuonDefinition(track);
  const auto& post = *stepData.GetPostStepPoint();
  auto* replacement =
      new G4Track(new G4DynamicParticle(stdMuon, post.GetMomentumDirection(),
                                        post.GetKineticEnergy()),
                  post.GetGlobalTime(), post.GetPosition());
  replacement->SetTouchableHandle(track.GetTouchableHandle());
  replacement->SetWeight(track.GetWeight());

  particleChange.AddSecondary(replacement);
  particleChange.ProposeTrackStatus(fStopAndKill);

  if (m_enableDebug) {
    G4cout << "MuonDISNuclearWrapperProcess: converted surviving primary "
           << stdMuon->GetParticleName()
           << " to an unbiased secondary after the biased DIS step." << G4endl;
  }
}
