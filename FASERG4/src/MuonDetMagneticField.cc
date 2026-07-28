#include "MuonDetMagneticField.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>

#include "G4TransportationManager.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

void MuonMagneticField::GetFieldValue(const G4double point[4], G4double* Bfield) const {
    // point[0] = x, point[1] = y, point[2] = z

    // Assume magnet box extends: x ∈ [-500, 500] mm, y ∈ [-500, 500] mm, z ∈ [–25, +25] mm (centered)
    G4double x = point[0];
    G4double y = point[1];
    G4double z = point[2];

    // return logical volume position for debugging
    // Use G4Navigator to find the volume at the given global point.
    // This is the correct approach for magnetic field calculation,
    // as it doesn't rely on the current track's G4TouchableHistory.
    G4String volumeName = "Unknown";
    G4ThreeVector globalPoint(x, y, z);

    #if 0
    // Get the navigator for the current world (assuming standard geometry)
    G4Navigator *navigator = G4TransportationManager::GetTransportationManager()
                                 ->GetNavigatorForTracking();

    if (navigator)
    {
      // Locate the volume at the global point
      G4VPhysicalVolume *physVol = navigator->LocateGlobalPointAndSetup(globalPoint);

      if (physVol)
      {
        volumeName = physVol->GetLogicalVolume()->GetName();
      }
    }

    // For debugging: print or log volumeName if needed
    // This debug output is now correct for a magnetic field calculation context.
    G4cout << "MagneticField calucation in volume: " << volumeName
           << " x = " << x / CLHEP::mm
           << " y = " << y / CLHEP::mm
           << " z = " << z / CLHEP::mm
           << G4endl;
#endif
    Bfield[0] = Bfield[1] = Bfield[2] = 0.0;

    // Assume ±1.5 Tesla in steel, depending on y (top/bottom vs center).
    // Rotation around Y leaves y unchanged, so the global y coordinate
    // already equals the tilted assembly's local y — no correction needed here.
    G4double Blocal = 0.0;
    if (std::abs(y) >= slitposition * mm && std::abs(y) <= 2 * slitposition * mm) {
      // Top or Bottom: +1.5 T
      Blocal = +1.5 * tesla;
    } else if (std::abs(y) < slitposition * mm) {
      // Middle: –1.5 T
      Blocal = -1.5 * tesla;
    }

    // The field points along the Fe slab's local +x axis, which is tilted
    // by tiltAngleY (rad) around Y with respect to the global frame (same
    // rotation as new G4RotationMatrix()->rotateY(fTiltAngleY) used to place
    // the detector assembly). Rotate the local field vector into global
    // coordinates so it stays aligned with the iron even when tilted.
    Bfield[0] = Blocal * std::cos(tiltAngleY);
    Bfield[2] = Blocal * std::sin(tiltAngleY);
}
/*
// previously used code
void MagneticField::GetFieldValue(const G4double point[4], G4double* Bfield) const {
    // point[0] = x, point[1] = y, point[2] = z

    // Assume magnet box extends: x ∈ [-500, 500] mm, y ∈ [-500, 500] mm, z ∈ [–25, +25] mm (centered)
    G4double x = point[0];
    G4double y = point[1];
    G4double z = point[2];

    // Default: no field
    Bfield[0] = Bfield[1] = Bfield[2] = 0.0;

    // Assume ±1.5 Tesla in steel, depending on y (top/bottom vs center)
    if (std::abs(y) >= 250. * mm && std::abs(y) <= 500. * mm) {
      // Top or Bottom: +1.5 T
      Bfield[0] = +1.5 * tesla;
    } else if (std::abs(y) < 250. * mm) {
      // Middle: –1.5 T
      Bfield[0] = -1.5 * tesla;
    }

    
}
*/