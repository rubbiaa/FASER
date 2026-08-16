#include "MuonDISInteractionRecorder.hh"

#include "G4SystemOfUnits.hh"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace {

struct InteractionContext {
  bool active{false};
  int geantEventId{-1};
  G4ThreeVector position{};
  double globalTime{0.0};
  int trackId{-1};
  int parentTrackId{-1};
  std::string volumeName{};
  bool hasLocalPosition{false};
  G4ThreeVector localPosition{};
};

thread_local InteractionContext currentInteractionContext;

std::string joinPdgs(const std::vector<int>& pdgs) {
  std::ostringstream output;
  for (size_t index = 0; index < pdgs.size(); ++index) {
    if (index != 0) {
      output << ';';
    }
    output << pdgs[index];
  }
  return output.str();
}

std::string csvQuote(const std::string& value) {
  std::string output = "\"";
  for (const char character : value) {
    if (character == '"') {
      output += "\"\"";
    } else {
      output += character;
    }
  }
  output += '"';
  return output;
}

}  // namespace

MuonDISInteractionRecorder& MuonDISInteractionRecorder::instance() {
  static MuonDISInteractionRecorder recorder;
  return recorder;
}

void MuonDISInteractionRecorder::configure(const std::string& outputPath) {
  std::scoped_lock lock{m_mutex};
  m_outputPath = outputPath;
  m_nextRecordIndex = 0;
  if (m_outputPath.empty()) {
    return;
  }

  std::ofstream output(m_outputPath, std::ios::out | std::ios::trunc);
  if (!output) {
    throw std::runtime_error("MuonDISInteractionRecorder: failed to open " + m_outputPath);
  }
  output
      << "interaction_index,geant_event_id,geant_muon_energy_GeV,q2_GeV2,pythia_x2,"
      << "target_nucleon_pdg,"
      << "track_id,parent_track_id,has_interaction_position,"
      << "interaction_x_mm,interaction_y_mm,interaction_z_mm,interaction_t_ns,"
      << "has_interaction_local_position,"
      << "interaction_local_x_mm,interaction_local_y_mm,interaction_local_z_mm,"
      << "volume_name,target_Z,target_A,dis_final_count,geant_secondary_count,dis_final_pdgs\n";
}

void MuonDISInteractionRecorder::setCurrentInteraction(int geantEventId,
                                                       const G4ThreeVector& position,
                                                       double globalTime,
                                                       int trackId,
                                                       int parentTrackId,
                                                       const std::string& volumeName,
                                                       bool hasLocalPosition,
                                                       const G4ThreeVector& localPosition) {
  currentInteractionContext.active = true;
  currentInteractionContext.geantEventId = geantEventId;
  currentInteractionContext.position = position;
  currentInteractionContext.globalTime = globalTime;
  currentInteractionContext.trackId = trackId;
  currentInteractionContext.parentTrackId = parentTrackId;
  currentInteractionContext.volumeName = volumeName;
  currentInteractionContext.hasLocalPosition = hasLocalPosition;
  currentInteractionContext.localPosition = localPosition;
}

void MuonDISInteractionRecorder::clearCurrentInteraction() {
  currentInteractionContext = {};
}

bool MuonDISInteractionRecorder::currentInteraction(int& geantEventId,
                                                    G4ThreeVector& position,
                                                    double& globalTime,
                                                    int& trackId,
                                                    int& parentTrackId,
                                                    std::string& volumeName,
                                                    bool& hasLocalPosition,
                                                    G4ThreeVector& localPosition) const {
  if (!currentInteractionContext.active) {
    return false;
  }
  geantEventId = currentInteractionContext.geantEventId;
  position = currentInteractionContext.position;
  globalTime = currentInteractionContext.globalTime;
  trackId = currentInteractionContext.trackId;
  parentTrackId = currentInteractionContext.parentTrackId;
  volumeName = currentInteractionContext.volumeName;
  hasLocalPosition = currentInteractionContext.hasLocalPosition;
  localPosition = currentInteractionContext.localPosition;
  return true;
}

void MuonDISInteractionRecorder::record(const Record& record) {
  std::scoped_lock lock{m_mutex};
  if (m_outputPath.empty()) {
    return;
  }

  std::ofstream output(m_outputPath, std::ios::out | std::ios::app);
  if (!output) {
    throw std::runtime_error("MuonDISInteractionRecorder: failed to append to " + m_outputPath);
  }

  const unsigned long long recordIndex = m_nextRecordIndex++;
  output << std::setprecision(17)
         << recordIndex << ','
         << record.geantEventId << ','
         << record.geantMuonEnergyGeV << ','
         << record.q2GeV2 << ','
         << record.pythiaX2 << ','
         << record.targetNucleonPdg << ','
         << record.trackId << ','
         << record.parentTrackId << ','
         << (record.hasInteractionPosition ? 1 : 0) << ','
         << record.interactionPosition.x() / mm << ','
         << record.interactionPosition.y() / mm << ','
         << record.interactionPosition.z() / mm << ','
         << record.interactionTimeNs / ns << ','
         << (record.hasInteractionLocalPosition ? 1 : 0) << ','
         << record.interactionLocalPosition.x() / mm << ','
         << record.interactionLocalPosition.y() / mm << ','
         << record.interactionLocalPosition.z() / mm << ','
         << csvQuote(record.volumeName) << ','
         << record.targetZ << ','
         << record.targetA << ','
         << record.disFinalCount << ','
         << record.geantSecondaryCount << ','
         << joinPdgs(record.disFinalPdgs) << '\n';
}
