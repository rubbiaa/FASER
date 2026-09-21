#ifndef _MUONSPECTROMETERFIELD_HH_
#define _MUONSPECTROMETERFIELD_HH_ 1

// Single, framework-independent implementation of the FASER muon
// spectrometer's dipole field model, shared by:
//
//   - the Geant4 simulation:     FASERG4/src/MuonDetMagneticField.cc
//                                (class MuonMagneticField)
//   - the GenFit reconstruction: CoreUtils/GenMagneticField.hh
//                                (class GenMagneticField, both its plain
//                                "Magnet" and "MDTMagnet" branches)
//
// Both are thin adapters: they convert their own position/unit
// conventions (G4: mm and Tesla; GenFit: cm and kGauss) to this
// function's (local y in cm, tilt in radians, field in kGauss) and
// convert the result back. The field MODEL - where the +/-1.5T
// boundaries sit and how the field vector rotates with the assembly's
// tilt - is written ONCE, here, so simulation and reconstruction cannot
// silently disagree about it the way they did before this file existed
// (see CoreUtils/GenMagneticField.hh's own history: its plain "Magnet"
// branch was missing the final tilt-rotation that both G4's field and
// GenFit's own "MDTMagnet" branch already applied).
//
// Tests/MuonSpectrometerFieldTest.cc is the gtest that cross-checks the
// G4 and GenFit adapters actually agree - run it whenever this file or
// either adapter changes.
namespace FASER {

struct MuonSpectrometerFieldParams {
    // Half-gap (cm) between the assembly's local y=0 plane and the
    // +/-fieldMagnitudeKG boundary. The modelled field envelope spans
    // |y_local| in [0, 2*slitPositionCm]; beyond that, the field is zero.
    double slitPositionCm = 25.0;

    // Assembly tilt around the global Y axis, in radians. Same sign
    // convention as DetectorConstruction's fTiltAngleY: the local field
    // direction (Blocal, 0, 0) is rotated into the global frame as
    // (Blocal*cos(tiltAngleRad), 0, Blocal*sin(tiltAngleRad)), matching
    // the G4RotationMatrix()->rotateY(fTiltAngleY) used to place the
    // detector assembly itself.
    double tiltAngleRad = 0.0;

    // |B| (kGauss) in each of the two field regions. 15 kG == 1.5 T,
    // the value both FASERG4 and GenFit have always used.
    double fieldMagnitudeKG = 15.0;

    // Half-width (cm) of a linear ramp centred on the slit boundary,
    // replacing the physically-real hard step there with a smooth
    // blend. 0 (the default - and what FASERG4's adapter uses)
    // reproduces the exact hard step G4 simulates. GenFit's adapters
    // pass a small nonzero value (currently 0.5, on the "MDTMagnet"
    // branch only) purely to keep its Runge-Kutta stepper numerically
    // stable while integrating across the discontinuity - this is a
    // reconstruction-side numerical accommodation, not a physical
    // difference from the field G4 actually simulates, which is why it
    // lives here as an explicit, visible parameter instead of a second
    // hand-copied implementation.
    double rampHalfWidthCm = 0.0;
};

// Evaluates the muon-spectrometer dipole field model at a given
// transverse position.
//
// y_local_cm: position along the magnet assembly's LOCAL y axis. A
// rotation around the global Y axis never changes y, so "local" and
// "global" y coincide once the caller has subtracted the assembly's own
// global-frame Y shift; that subtraction is each adapter's job, not
// this function's - see MuonMagneticField::GetFieldValue's "y -
// centreY" and GenMagneticField::get()'s "position.Y() -
// rearMuSpec_LOS_shiftY" for the two existing examples.
//
// Only y_local_cm enters this decision: the modelled assembly is
// symmetric in x, and its z-extent is enforced upstream by each caller
// (the G4 volume boundary that GetFieldValue is only ever invoked
// within; GenFit's magnet_z_ranges_cm_/mdt_magnet_z_ranges_cm_ gating)
// before this function is ever reached.
//
// Writes Bx_kG, By_kG, Bz_kG (kGauss) and returns true if y_local_cm
// falls inside the modelled |y_local| <= 2*slitPositionCm envelope;
// returns false and zeroes the field otherwise.
bool ComputeMuonSpectrometerField(double y_local_cm,
                                   const MuonSpectrometerFieldParams& params,
                                   double& Bx_kG, double& By_kG, double& Bz_kG);

} // namespace FASER

#endif
