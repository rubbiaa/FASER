#ifndef MuonDISPhysics_h
#define MuonDISPhysics_h 1

// Registers the MuonDIS muon-nuclear process replacement, following exactly
// the same pattern as TauDecayPhysics in this same directory: a
// G4VPhysicsConstructor added directly to the physics list in faserps.cc
// (physicsList->RegisterPhysics(new MuonDISPhysics())), rather than the
// Athena AthAlgTool + Gaudi::Property wrapper the original Calypso
// MuonDISPhysicsTool used (there is no Gaudi/Athena in this app).
//
// Disabled by default -- opt in from your run macro with, BEFORE
// /run/initialize:
//   /physics/muondis/enable true
//   /physics/muondis/crossSectionBias 150
//   /physics/muondis/q2min 1.0
//   /physics/muondis/interactionLog muondis_interactions.csv
//   /physics/muondis/pdfSet input/NNPDF40_nnlo_as_01180_charmasy_0000.dat
//   /physics/muondis/xbjmin 0.01
// See MuonDISMessenger.hh for the full command list and MuonDISPythiaGenerator.hh
// for the physics settings these commands feed into.

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class MuonDISMessenger;

class MuonDISPhysics : public G4VPhysicsConstructor {
public:
  explicit MuonDISPhysics(G4int verb = 1);
  ~MuonDISPhysics() override;

  void ConstructParticle() override;
  void ConstructProcess() override;

  // Called from MuonDISMessenger; must happen before /run/initialize to have
  // any effect (see class comment).
  void SetEnabled(G4bool value) { fEnabled = value; }
  void SetCrossSectionBias(G4double value) { fCrossSectionBias = value; }
  void SetInteractionLogPath(const G4String& path) { fInteractionLogPath = path; }
  void SetDebug(G4bool value) { fEnableDebug = value; }
  void SetQ2Min(G4double value) { fQ2Min = value; }
  void SetPdfSet(const G4String& path) { fPdfSetPath = path; }
  void SetXBjMin(G4double value) { fXBjMin = value; }

private:
  G4bool fEnabled{false};
  G4double fCrossSectionBias{150.0};
  G4String fInteractionLogPath{""};
  G4bool fEnableDebug{false};
  G4double fQ2Min{1.0};
  G4String fPdfSetPath{""};
  G4double fXBjMin{0.0};
  MuonDISMessenger* fMessenger{};
};

#endif
