#ifndef MDTSD_h
#define MDTSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include "globals.hh"

class G4Step;
class G4TouchableHistory;

/**
 * @brief Data structure for MDT (Monitored Drift Tube) hit information
 * 
 * MDT measures:
 * - X position along the tube (from wire position)
 * - Drift radius (perpendicular distance from wire in Y-Z plane)
 * - Drift time (what is actually measured)
 * Note: Y and Z are NOT measured individually - only radial distance!
 */
struct MDTHit {
  G4ThreeVector tubeCenter;
  G4double hitX;              // local X along tube (= global X - tubeCenter.x())
  G4double trueDriftRadius;
  G4double driftAngle;
  G4ThreeVector truePosition;     // raw Geant4 pre-step position
  G4ThreeVector closestPosition;  // closest approach to MDT wire
  G4ThreeVector trueMomentum;     // raw Geant4 pre-step momentum

  G4int trackID;
  G4int pdgID;
  G4int stationID;
  G4int planeID;
  G4int tubeID;

  // Optional:
  G4double edep;
  G4double driftTime;
};

/**
 * @brief Sensitive Detector for ATLAS-style MDT (Monitored Drift Tubes)
 * 
 * Simulates drift tube response with:
 * - 30 mm outer diameter aluminum tubes
 * - 400 µm wall thickness
 * - ~80 µm spatial resolution
 * - Tubes oriented horizontally along X-axis
 * - Staggered rows in Y-direction (honeycomb pattern)
 */
class MDTSD : public G4VSensitiveDetector {
public:
  MDTSD(const G4String& name);
  virtual ~MDTSD();
  
  virtual void   Initialize(G4HCofThisEvent* hce) override;
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
  
  const std::vector<MDTHit>& GetHits() const { return fHits; }
  void ClearHits() { fHits.clear(); }
    
  static void SetVerbose(int v) { fVerbose = v; }
  static int  GetVerbose()      { return fVerbose; }

  
private:
  std::vector<MDTHit> fHits;
  static int fVerbose;
  
};

#endif // MDTSD_h
