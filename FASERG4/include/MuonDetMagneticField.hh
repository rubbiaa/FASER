#ifndef MAGNETICFIELD_HH
#define MAGNETICFIELD_HH

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"

class MuonMagneticField : public G4MagneticField {
public:
    MuonMagneticField() = default;
    ~MuonMagneticField() = default;

    G4double slitposition = 0.0; // position of the slit along y in mm
    G4double tiltAngleY = 0.0;   // detector assembly tilt around Y, in radians (same as fTiltAngleY)

    void SetSlitPosition(G4double pos) { slitposition = pos; }
    void SetTiltAngleY(G4double angle) { tiltAngleY = angle; }
    virtual void GetFieldValue(const G4double point[4], G4double* Bfield) const override;
};

#endif

