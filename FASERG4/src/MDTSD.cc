#include "MDTSD.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHistory.hh"
#include "G4TouchableHandle.hh"
#include "G4HCofThisEvent.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include <cmath>
#include <limits>

MDTSD::MDTSD(const G4String& name)
  : G4VSensitiveDetector(name)
{}

MDTSD::~MDTSD()
{}

int MDTSD::fVerbose = 1;

void MDTSD::Initialize(G4HCofThisEvent*)
{
  fHits.clear();
}

G4bool MDTSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  if (!step) return false;

  G4double edep = step->GetTotalEnergyDeposit();

  // Ignore tiny deposits
  if (edep < 0.1 * keV) return false;

  G4StepPoint* pre  = step->GetPreStepPoint();
  G4StepPoint* post = step->GetPostStepPoint();
  if (!pre || !post) return false;

  G4Track* track = step->GetTrack();
  if (!track) return false;

  G4TouchableHandle touchable = pre->GetTouchableHandle();
  if (!touchable) return false;

  // ------------------------------------------------------------
  // Get global tube/wire center
  //
  // The sensitive volume is the gas volume inside the MDT tube.
  // Its local origin is the wire/tube center.
  // ------------------------------------------------------------
  G4ThreeVector tubeCenter =
    touchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0., 0., 0.));

  // ------------------------------------------------------------
  // Copy-number decoding
  //
  // Outer tube copy number:
  // copyNumber = (stationID) * 10000 + planeID * 1000 + tubeID
  //
  // Since the gas volume is placed inside the tube, use parent
  // copy number: GetCopyNumber(1)
  // ------------------------------------------------------------
  G4int copyNumber = touchable->GetCopyNumber(1);

  G4int tubeID    = copyNumber % 1000;
  G4int planeID   = (copyNumber / 1000) % 10;
  G4int stationID = copyNumber / 10000;

  G4int trackID = track->GetTrackID();
  G4int pdgID   = track->GetDefinition()->GetPDGEncoding();

  // Optional: keep only muons
  // if (std::abs(pdgID) != 13) return false;

  // ------------------------------------------------------------
  // True step positions
  // ------------------------------------------------------------
  G4ThreeVector p0 = pre->GetPosition();
  G4ThreeVector p1 = post->GetPosition();

  G4ThreeVector truePos = p0;
  G4ThreeVector trueMom = pre->GetMomentum();

  // ------------------------------------------------------------
  // MDT wire is along global X (assumption: tubes oriented along X).
  // If the geometry changes tube orientation this must be updated.
  // Therefore the drift radius is calculated in the Y-Z plane.
  //
  // We compute the minimum distance between the Geant4 step segment
  // and the wire center in the Y-Z plane.
  // ------------------------------------------------------------
  G4double y0 = p0.y() - tubeCenter.y();
  G4double z0 = p0.z() - tubeCenter.z();

  G4double y1 = p1.y() - tubeCenter.y();
  G4double z1 = p1.z() - tubeCenter.z();

  G4double dy = y1 - y0;
  G4double dz = z1 - z0;

  G4double denom = dy * dy + dz * dz;

  G4double t = 0.0;

  if (denom > 1e-12 * mm * mm) {
    t = -(y0 * dy + z0 * dz) / denom;

    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
  }

  // Closest point on the step segment to the wire center in Y-Z plane
  // This is the point of closest approach to the wire
  // close to what we see in the real detector, where the drift radius is measured.
  G4double yClosest = y0 + t * dy;
  G4double zClosest = z0 + t * dz;

  G4double trueDriftRadius =
    std::sqrt(yClosest * yClosest + zClosest * zClosest);

  G4double driftAngle = std::atan2(zClosest, yClosest);

  G4ThreeVector closestPos(
    p0.x() + t * (p1.x() - p0.x()),
    tubeCenter.y() + yClosest,
    tubeCenter.z() + zClosest
  );

  // hitX: position along the tube in local (tube-centered) X coordinates
  G4double hitX = closestPos.x() - tubeCenter.x();

  // Simple reference drift time
  const G4double driftVelocity = 0.025 * mm / ns;
  G4double driftTime = trueDriftRadius / driftVelocity;

  // ------------------------------------------------------------
  // Merge Geant4 steps into one MDT hit per:
  //
  //   trackID + stationID + planeID + tubeID
  //
  // Keep the minimum drift radius.
  // ------------------------------------------------------------
  for (auto& oldHit : fHits) {

    if (oldHit.trackID   == trackID   &&
        oldHit.stationID == stationID &&
        oldHit.planeID   == planeID   &&
        oldHit.tubeID    == tubeID) 
        {
          // Keep the closest approach to the wire
          if (trueDriftRadius < oldHit.trueDriftRadius) 
          {
            oldHit.trueDriftRadius = trueDriftRadius;
            oldHit.driftAngle      = driftAngle;
            // MDT measurement-relevant position
            oldHit.closestPosition = closestPos;
            // Raw Geant4 truth step position
            oldHit.truePosition    = truePos;

            oldHit.trueMomentum    = trueMom;
            oldHit.hitX            = hitX;
            oldHit.driftTime       = driftTime;
          }

        // Optional: if MDTHit has edep member, accumulate it:
        oldHit.edep += edep;
        return true;
      }
  }

  // ------------------------------------------------------------
  // First hit in this tube for this track
  // ------------------------------------------------------------
  MDTHit hit;

  hit.tubeCenter      = tubeCenter;
  hit.hitX            = hitX;
  hit.trueDriftRadius = trueDriftRadius;
  hit.driftAngle      = driftAngle;
  hit.truePosition    = truePos;
  hit.trueMomentum    = trueMom;
  hit.closestPosition = closestPos;

  hit.trackID   = trackID;
  hit.pdgID     = pdgID;
  hit.stationID = stationID;
  hit.planeID   = planeID;
  hit.tubeID    = tubeID;

  // Optional, if your MDTHit struct has these:
  hit.edep = edep;
  hit.driftTime = driftTime;

  fHits.push_back(hit);
  
  if (fVerbose >= 1) {
    G4cout
      << "MDT Hit:"
      << " station=" << stationID
      << " plane="   << planeID
      << " tube="    << tubeID
      << " trackID=" << trackID
      << " pdgID="   << pdgID
      << " r="       << trueDriftRadius/mm << " mm"
      << G4endl;
  }
  if (fVerbose >= 2) {
    G4cout
      << "  rawPos=("
      << p0.x()/mm << ", "
      << p0.y()/mm << ", "
      << p0.z()/mm << ") mm"
      << " closestPos=("
      << closestPos.x()/mm << ", "
      << closestPos.y()/mm << ", "
      << closestPos.z()/mm << ") mm"
      << G4endl;
  }

  if(fVerbose >= 3)  
    G4cout << "MDT Hit: "
           << "station=" << stationID
           << " plane=" << planeID
           << " tube=" << tubeID
           << " trackID=" << trackID
           << " pdgID=" << pdgID
           << " hitX=" << hitX / mm << " mm"
           << " trueDriftRadius=" << trueDriftRadius / mm << " mm"
           << " driftTime=" << driftTime / ns << " ns"
           << " driftAngle=" << driftAngle << " rad"
           << " tubeCenter=("
           << tubeCenter.x() / mm << ", "
           << tubeCenter.y() / mm << ", "
           << tubeCenter.z() / mm << ") mm"
           << " edep=" << edep / keV << " keV"
          << G4endl;

  return true;
}
