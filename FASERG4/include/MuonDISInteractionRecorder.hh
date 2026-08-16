#ifndef MuonDISInteractionRecorder_h
#define MuonDISInteractionRecorder_h 1

// Ported from the Calypso MuonDIS Geant4 extension (Simulation/G4Extensions/MuonDIS)
// with no functional changes other than swapping the POWHEG-library fields for the
// on-the-fly Pythia8 diagnostic fields (Q2, chosen target nucleon).
//
// Thread-safety note: this is a process-wide singleton with a std::mutex-guarded
// append; the "current interaction" context (set by MuonDISNuclearWrapperProcess
// right before calling the biased process) is thread_local. This app currently
// always runs with /run/numberOfThreads 1, so this has not been exercised
// multi-threaded -- revisit if that changes.

#include "G4ThreeVector.hh"

#include <mutex>
#include <string>
#include <vector>

class MuonDISInteractionRecorder {
public:
  struct Record {
    int geantEventId{-1};
    double geantMuonEnergyGeV{0.0};
    double q2GeV2{0.0};
    double pythiaX2{0.0};  // Info::x2() from MuonDISPythiaGenerator -- Bjorken x actually used
                           // to accept this event; cross-check against TPOEvent::xBj.
    int targetNucleonPdg{0};
    int trackId{-1};
    int parentTrackId{-1};
    bool hasInteractionPosition{false};
    G4ThreeVector interactionPosition{};
    bool hasInteractionLocalPosition{false};
    G4ThreeVector interactionLocalPosition{};
    double interactionTimeNs{0.0};
    std::string volumeName{};
    int targetZ{0};
    int targetA{0};
    int disFinalCount{0};
    int geantSecondaryCount{0};
    std::vector<int> disFinalPdgs{};
  };

  static MuonDISInteractionRecorder& instance();

  void configure(const std::string& outputPath);
  void setCurrentInteraction(int geantEventId,
                             const G4ThreeVector& position, double globalTime,
                             int trackId, int parentTrackId,
                             const std::string& volumeName,
                             bool hasLocalPosition,
                             const G4ThreeVector& localPosition);
  void clearCurrentInteraction();
  bool currentInteraction(int& geantEventId,
                          G4ThreeVector& position, double& globalTime,
                          int& trackId, int& parentTrackId,
                          std::string& volumeName,
                          bool& hasLocalPosition,
                          G4ThreeVector& localPosition) const;
  void record(const Record& record);
  const std::string& outputPath() const { return m_outputPath; }
  bool enabled() const { return !m_outputPath.empty(); }

private:
  MuonDISInteractionRecorder() = default;

  mutable std::mutex m_mutex{};
  std::string m_outputPath{};
  unsigned long long m_nextRecordIndex{0};
};

#endif
