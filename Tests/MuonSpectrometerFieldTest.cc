// Cross-checks that the Geant4 simulation's muon-spectrometer field
// (MuonMagneticField, FASERG4/src/MuonDetMagneticField.cc) and the
// GenFit reconstruction's muon-spectrometer field (GenMagneticField,
// CoreUtils/GenMagneticField.hh) agree - both are thin adapters around
// the one shared field model in CoreUtils/MuonSpectrometerField.hh, and
// this is what actually proves the two adapters wire it up identically
// rather than just trusting that by inspection. Run this after touching
// either adapter or the shared field model itself.
#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "MuonSpectrometerField.hh"
#include "MuonDetMagneticField.hh" // G4 adapter
#include "GenMagneticField.hh"     // GenFit adapter
#include "G4SystemOfUnits.hh"

namespace {

constexpr double kSlitPositionCm = 25.0;
constexpr double kTiltDeg = 10.0; // nonzero, to exercise the tilt rotation on both sides
constexpr double kTiltRad = kTiltDeg * M_PI / 180.0;
constexpr double kShiftYCm = 3.0; // nonzero global-frame Y shift of the assembly

// Queries a G4 MuonMagneticField (mm in, Tesla out) at a position given
// in cm and returns the result in kGauss, GenFit's native unit, so the
// two sides can be compared directly.
void QueryG4FieldKG(const MuonMagneticField& field, double x_cm, double y_cm, double z_cm,
                     double& Bx_kG, double& By_kG, double& Bz_kG) {
    const G4double point[4] = {x_cm * 10.0, y_cm * 10.0, z_cm * 10.0, 0.0}; // cm -> mm
    G4double Bfield[3];
    field.GetFieldValue(point, Bfield);
    Bx_kG = (Bfield[0] / tesla) * 10.0; // Tesla -> kGauss (1 T = 10 kG)
    By_kG = (Bfield[1] / tesla) * 10.0;
    Bz_kG = (Bfield[2] / tesla) * 10.0;
}

// y_local offsets (cm) sampled around both boundaries (the slit at
// slitPositionCm, and the outer envelope edge at 2*slitPositionCm) plus
// a point clearly beyond the modelled envelope, mirrored to +/-.
const std::vector<double> kYLocalOffsetsCm = {
    0.0, 12.5, 24.9, 25.0, 25.1, 37.5, 49.9, 50.0, 50.1, 60.0,
    -12.5, -24.9, -25.0, -25.1, -37.5, -49.9, -50.0, -50.1, -60.0,
};

} // namespace

// The shared field model itself, independent of either adapter: sanity
// check the sign/magnitude at a few unambiguous points with zero tilt.
TEST(MuonSpectrometerField, SanityAtZeroTilt) {
    FASER::MuonSpectrometerFieldParams params;
    params.slitPositionCm = kSlitPositionCm;

    double Bx_kG, By_kG, Bz_kG;

    ASSERT_TRUE(FASER::ComputeMuonSpectrometerField(0.0, params, Bx_kG, By_kG, Bz_kG));
    EXPECT_NEAR(Bx_kG, -15.0, 1e-9); // middle region: -1.5 T
    EXPECT_NEAR(Bz_kG, 0.0, 1e-9);   // zero tilt: field stays along local x

    ASSERT_TRUE(FASER::ComputeMuonSpectrometerField(37.5, params, Bx_kG, By_kG, Bz_kG));
    EXPECT_NEAR(Bx_kG, +15.0, 1e-9); // top/bottom region: +1.5 T

    EXPECT_FALSE(FASER::ComputeMuonSpectrometerField(60.0, params, Bx_kG, By_kG, Bz_kG));
    EXPECT_NEAR(Bx_kG, 0.0, 1e-9); // beyond the modelled envelope: zero field
}

// G4's field and GenFit's plain "Magnet" branch (event-station gating,
// no ramp on either side) must agree EXACTLY, including at the slit
// boundary itself, since neither smooths the step there.
TEST(MuonSpectrometerFieldConsistency, MagnetBranchMatchesG4Exactly) {
    MuonMagneticField g4field;
    g4field.SetSlitPosition(kSlitPositionCm * 10.0); // cm -> mm
    g4field.SetTiltAngleY(kTiltRad);
    g4field.SetCentreY(kShiftYCm * 10.0); // cm -> mm

    GenMagneticField genField;
    genField.SetSlitPosition(kSlitPositionCm);
    genField.SetEventStationZsCm({-1000.0, 1000.0}); // wide "on" envelope; z=0 stays far from either dead-band
    // rearMuSpec_tilt_deg's sign convention is the negative of
    // fTiltAngleY/kTiltRad's (see GenMagneticField.hh and
    // TPORecoEvent.cc) - pass -kTiltDeg so the class's own internal
    // sign flip nets back to +kTiltRad, matching the G4 field above.
    genField.SetRearMuSpectShift(0.0, kShiftYCm, 0.0, -kTiltDeg);
    genfit::AbsBField* genFieldBase = &genField; // get() is private on GenMagneticField itself

    const double x_cm = 7.0, z_cm = 0.0; // arbitrary and nonzero: the field model doesn't depend on x/z
    for (double yLocal_cm : kYLocalOffsetsCm) {
        const double y_cm = yLocal_cm + kShiftYCm;

        double g4Bx_kG, g4By_kG, g4Bz_kG;
        QueryG4FieldKG(g4field, x_cm, y_cm, z_cm, g4Bx_kG, g4By_kG, g4Bz_kG);

        const TVector3 genB_kG = genFieldBase->get(TVector3(x_cm, y_cm, z_cm));

        SCOPED_TRACE(::testing::Message() << "y_local_cm=" << yLocal_cm);
        EXPECT_NEAR(g4Bx_kG, genB_kG.X(), 1e-9);
        EXPECT_NEAR(g4By_kG, genB_kG.Y(), 1e-9);
        EXPECT_NEAR(g4Bz_kG, genB_kG.Z(), 1e-9);
    }
}

// G4's field and GenFit's "MDTMagnet" branch must agree everywhere
// EXCEPT inside the +/-0.5 cm ramp GenFit's adapter deliberately adds
// around the slit boundary for its Runge-Kutta stepper's numerical
// stability - see MuonSpectrometerFieldParams::rampHalfWidthCm's doc
// comment. Inside that window we only check GenFit's value stays a
// physically sane blend between the two field extremes, since it is
// *meant* to differ from G4's exact hard step there.
TEST(MuonSpectrometerFieldConsistency, MDTMagnetBranchMatchesG4OutsideRampZone) {
    constexpr double kRampHalfWidthCm = 0.5; // must match GenMagneticField.hh's MDT branch

    MuonMagneticField g4field;
    g4field.SetSlitPosition(kSlitPositionCm * 10.0);
    g4field.SetTiltAngleY(kTiltRad);
    g4field.SetCentreY(kShiftYCm * 10.0);

    GenMagneticField genField;
    genField.SetSlitPosition(kSlitPositionCm);
    genField.SetMDTMagnetZRangesCm({{-1000.0, 1000.0}});
    genField.SetRearMuSpectShift(0.0, kShiftYCm, 0.0, -kTiltDeg);
    genfit::AbsBField* genFieldBase = &genField;

    const double x_cm = -4.0, z_cm = 0.0;
    for (double yLocal_cm : kYLocalOffsetsCm) {
        const double y_cm = yLocal_cm + kShiftYCm;

        double g4Bx_kG, g4By_kG, g4Bz_kG;
        QueryG4FieldKG(g4field, x_cm, y_cm, z_cm, g4Bx_kG, g4By_kG, g4Bz_kG);

        const TVector3 genB_kG = genFieldBase->get(TVector3(x_cm, y_cm, z_cm));

        SCOPED_TRACE(::testing::Message() << "y_local_cm=" << yLocal_cm);
        if (std::abs(std::abs(yLocal_cm) - kSlitPositionCm) <= kRampHalfWidthCm) {
            // Inside GenFit's deliberate smoothing window: just check
            // it's a sane blend (magnitude never exceeds the field
            // strength), not equal to G4's hard step.
            const double magnitude_kG = std::hypot(genB_kG.X(), genB_kG.Z());
            EXPECT_LE(magnitude_kG, 15.0 + 1e-9);
            continue;
        }
        EXPECT_NEAR(g4Bx_kG, genB_kG.X(), 1e-9);
        EXPECT_NEAR(g4By_kG, genB_kG.Y(), 1e-9);
        EXPECT_NEAR(g4Bz_kG, genB_kG.Z(), 1e-9);
    }
}
