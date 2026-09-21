#include "MuonDetMagneticField.hh"
#include "MuonSpectrometerField.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>

void MuonMagneticField::GetFieldValue(const G4double point[4], G4double* Bfield) const {
    // point[0] = x, point[1] = y, point[2] = z
    const G4double y = point[1];

    Bfield[0] = Bfield[1] = Bfield[2] = 0.0;

    // Translate to the assembly's local y (relative to magnet center):
    // the field boundaries (slit position, +/-1.5T regions) are defined
    // in the assembly's LOCAL frame, not the global frame, and the
    // detector can be shifted in Y in global coordinates.
    const double y_local_cm = (y - centreY) / cm;

    // See CoreUtils/MuonSpectrometerField.hh for what this shared
    // function computes and why: this is the SAME field model GenFit's
    // reconstruction uses (CoreUtils/GenMagneticField.hh), so simulation
    // and reconstruction can no longer silently disagree about it.
    FASER::MuonSpectrometerFieldParams params;
    params.slitPositionCm = slitposition / cm;
    params.tiltAngleRad = tiltAngleY;
    params.fieldMagnitudeKG = 15.0; // 1.5 T, matching GenFit's field model
    // rampHalfWidthCm left at its default (0): G4 reproduces the exact
    // hard step at the slit boundary. GenFit's own field additionally
    // ramps this over a small window purely for its Runge-Kutta
    // stepper's numerical stability - see rampHalfWidthCm's doc comment
    // in CoreUtils/MuonSpectrometerField.hh.

    double Bx_kG = 0.0, By_kG = 0.0, Bz_kG = 0.0;
    FASER::ComputeMuonSpectrometerField(y_local_cm, params, Bx_kG, By_kG, Bz_kG);

    // kGauss (GenFit's native unit) -> G4 internal field units: 1 T = 10 kG.
    Bfield[0] = (Bx_kG / 10.0) * tesla;
    Bfield[1] = (By_kG / 10.0) * tesla;
    Bfield[2] = (Bz_kG / 10.0) * tesla;
}
