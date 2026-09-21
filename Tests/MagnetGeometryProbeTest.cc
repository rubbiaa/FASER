// Unit tests for FASER::ProbeMagnetSlit (CoreUtils/MagnetGeometryProbe.hh),
// independent of any G4/GDML machinery: builds small ROOT TGeo shapes
// directly, by hand, to check the probe finds the right answer on a
// shape shaped just like the real MDT magnet (a block with two
// symmetric slits), correctly flags shapes that don't match that
// pattern, and - the actual point of this whole exercise - correctly
// DETECTS a real mismatch between the field model's assumed envelope
// and the block's true physical edge, rather than just reproducing the
// nominal numbers by coincidence.
#include <gtest/gtest.h>

#include "MagnetGeometryProbe.hh"

#include <TGeoManager.h>
#include <TGeoBBox.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TGeoMatrix.h>

namespace {

// A fresh TGeoManager is required before building any TGeo shapes (they
// register themselves with the current one). One process-wide instance,
// shared by every TEST below, is enough - none of these tests navigate
// or close the geometry, they only call TGeoShape::Contains() directly
// on hand-built shapes.
void EnsureGeoManager() {
    if (!gGeoManager) {
        new TGeoManager("MagnetGeometryProbeTestGeom", "for MagnetGeometryProbeTest");
    }
}

// Builds a solid block of half-extents (halfX, halfY, halfZ) cm with two
// slits of half-extents (slitHalfX, slitHalfHeight, halfZ) cm cut in at
// y = +-slitPositionCm - exactly the shape
// FASERG4::DetectorConstruction::CreateMuSpectWithMDT builds (just named
// uniquely per test so repeated TGeo shape names never collide).
TGeoShape* BuildBlockWithSlits(const char* nameTag, double halfX, double halfY, double halfZ,
                                double slitHalfX, double slitHalfHeight, double slitPositionCm) {
    EnsureGeoManager();
    std::string block = std::string("block_") + nameTag;
    std::string slit  = std::string("slit_") + nameTag;
    std::string s1    = std::string("minusSlit1_") + nameTag;
    std::string s2    = std::string("withSlits_") + nameTag;

    TGeoBBox* blockShape = new TGeoBBox(block.c_str(), halfX, halfY, halfZ);
    TGeoBBox* slitShape  = new TGeoBBox(slit.c_str(), slitHalfX, slitHalfHeight, halfZ);

    auto* trUp   = new TGeoTranslation(0, slitPositionCm, 0);
    auto* trDown = new TGeoTranslation(0, -slitPositionCm, 0);

    auto* sub1 = new TGeoSubtraction(blockShape, slitShape, nullptr, trUp);
    auto* minusSlit1 = new TGeoCompositeShape(s1.c_str(), sub1);

    auto* sub2 = new TGeoSubtraction(minusSlit1, slitShape, nullptr, trDown);
    auto* withSlits = new TGeoCompositeShape(s2.c_str(), sub2);

    return withSlits;
}

} // namespace

// Nominal case: exactly the real MDT magnet's numbers (in cm: 50/50 half
// extent, 25/1 slit half-extent, slit centred at 25cm) - envelope
// (2*slitPositionCm=50cm) matches the block's real half-extent (50cm) by
// construction, same as production today.
TEST(MagnetGeometryProbe, NominalMDTMagnetShapeMatchesExpectedNumbers) {
    TGeoShape* shape = BuildBlockWithSlits("nominal", 50.0, 50.0, 20.0, 25.0, 1.0, 25.0);
    FASER::MagnetSlitProbeResult result = FASER::ProbeMagnetSlit(shape);

    ASSERT_TRUE(result.ok);
    EXPECT_NEAR(result.slitPositionCm, 25.0, 0.1);
    EXPECT_NEAR(result.slitHalfWidthCm, 1.0, 0.15);
    EXPECT_NEAR(result.blockHalfExtentYCm, 50.0, 1e-9);
    EXPECT_NEAR(result.asymmetryCm, 0.0, 0.1);

    // The actual point of the whole exercise: this must hold, or the
    // field's zero-cutoff envelope doesn't line up with the real iron
    // edge.
    EXPECT_NEAR(2.0 * result.slitPositionCm, result.blockHalfExtentYCm, 0.2);
}

// The real regression this probe exists to catch: someone widens the
// block (say, to give the MDT tubes more Y coverage) without touching
// the slit position. The probe must report the TRUE, now-mismatched
// numbers - not silently agree with 2*slitPositionCm just because that
// happened to be true before.
TEST(MagnetGeometryProbe, DetectsEnvelopeGeometryMismatch) {
    // Same slit as the real magnet (still centred at 25cm), but the
    // block is now 70cm half-extent instead of 50cm: the field's
    // envelope (2*25=50cm) would incorrectly go to zero 20cm before the
    // iron actually ends.
    TGeoShape* shape = BuildBlockWithSlits("widened_block", 50.0, 70.0, 20.0, 25.0, 1.0, 25.0);
    FASER::MagnetSlitProbeResult result = FASER::ProbeMagnetSlit(shape);

    ASSERT_TRUE(result.ok);
    EXPECT_NEAR(result.slitPositionCm, 25.0, 0.1);
    EXPECT_NEAR(result.blockHalfExtentYCm, 70.0, 1e-9);

    // This is the failure SetupMDTMagneticField's own runtime check is
    // built to catch and warn about: envelope != real half-extent.
    EXPECT_GT(std::abs(2.0 * result.slitPositionCm - result.blockHalfExtentYCm), 10.0);
}

// Asymmetric slits (deliberately mis-built, to exercise the asymmetry
// diagnostic): +y slit at 20cm, -y slit at 30cm.
TEST(MagnetGeometryProbe, DetectsAsymmetricSlits) {
    EnsureGeoManager();
    TGeoBBox* blockShape = new TGeoBBox("block_asym", 50.0, 50.0, 20.0);
    TGeoBBox* slitShape  = new TGeoBBox("slit_asym", 25.0, 1.0, 20.0);
    auto* trUp   = new TGeoTranslation(0, 20.0, 0);
    auto* trDown = new TGeoTranslation(0, -30.0, 0);
    auto* sub1 = new TGeoSubtraction(blockShape, slitShape, nullptr, trUp);
    auto* minusSlit1 = new TGeoCompositeShape("minusSlit1_asym", sub1);
    auto* sub2 = new TGeoSubtraction(minusSlit1, slitShape, nullptr, trDown);
    TGeoShape* shape = new TGeoCompositeShape("withSlits_asym", sub2);

    FASER::MagnetSlitProbeResult result = FASER::ProbeMagnetSlit(shape);

    ASSERT_TRUE(result.ok);
    EXPECT_NEAR(result.asymmetryCm, 10.0, 0.2); // |20 - 30|
}

// A plain solid block with no slit at all must be reported as NOT
// matching this probe's expected pattern, rather than returning some
// meaningless number - callers must be able to trust ok=false as "don't
// use this value".
TEST(MagnetGeometryProbe, PlainBlockWithNoSlitIsNotOk) {
    EnsureGeoManager();
    TGeoShape* shape = new TGeoBBox("plain_block", 50.0, 50.0, 20.0);
    FASER::MagnetSlitProbeResult result = FASER::ProbeMagnetSlit(shape);
    EXPECT_FALSE(result.ok);
}

// A null shape must be handled safely, not crash.
TEST(MagnetGeometryProbe, NullShapeIsNotOk) {
    FASER::MagnetSlitProbeResult result = FASER::ProbeMagnetSlit(nullptr);
    EXPECT_FALSE(result.ok);
}
