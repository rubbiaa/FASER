#ifndef _MAGNETGEOMETRYPROBE_HH_
#define _MAGNETGEOMETRYPROBE_HH_ 1

class TGeoShape;

// Probes a magnet solid's ACTUAL shape - as loaded into ROOT, normally
// from the very same GDML file FASERG4 exports (see
// DetectorConstruction::Construct()'s parser.Write("FASERCAL_V10.gdml",
// ...) and Batch/BatchReco.cc's TGeoManager::Import(...) of it) - to find
// where it physically transitions from solid iron to the horizontal slit
// gap that splits the muon spectrometer's field into two
// oppositely-magnetized halves, and how far the block itself extends in
// y.
//
// Why this exists: FASER::MuonSpectrometerFieldParams::slitPositionCm
// (CoreUtils/MuonSpectrometerField.hh) decides where GenFit's
// reconstruction thinks the field flips sign, and 2*slitPositionCm
// decides where it thinks the field stops entirely. That value used to
// be a hand-typed literal (GenMagneticField.cc's SetupMDTMagneticField
// called SetSlitPosition(25.0)) with zero structural connection to the
// number FASERG4::DetectorConstruction actually built the iron block and
// its slit cutouts with - so a future change to the block's geometry
// could silently desynchronize the field model from the geometry it's
// supposed to describe, with nothing to catch it. This probe closes that
// gap by deriving the value from the real, currently-loaded shape every
// time, instead of trusting the two sides to be kept in sync by hand.
namespace FASER {

struct MagnetSlitProbeResult {
    // false if `shape` didn't look like a plain solid block with a
    // single horizontal slit cut symmetrically into each side (e.g. a
    // null shape, or one with no detectable gap) - callers MUST check
    // this and fall back to a known-good default rather than trusting a
    // probe that didn't find what it expected.
    bool ok = false;

    // |y_local| (cm) where solid material ends and the slit gap begins:
    // the value to feed to MuonSpectrometerFieldParams::slitPositionCm.
    // Averaged between the +y and -y slits.
    double slitPositionCm = 0.0;

    // Half-height (cm) of the detected gap itself, for diagnostics only
    // (the field model treats the slit as a zero-width flip, not a
    // modelled gap - see MuonSpectrometerField.cc).
    double slitHalfWidthCm = 0.0;

    // Half-extent (cm) of the block's own bounding box along y: the
    // block's true physical outer edge. The field model's own envelope
    // cutoff (2*slitPositionCm) is SUPPOSED to coincide with this -
    // callers should compare the two and warn (not silently ignore) if
    // they don't.
    double blockHalfExtentYCm = 0.0;

    // |(+y slit position) - (-y slit position)|, cm: how far the two
    // slits deviate from being perfectly symmetric about y=0. Small
    // asymmetry can be a real (if tiny) construction detail; a large one
    // means this probe's assumptions (a single symmetric slit pair) may
    // not actually hold for this shape.
    double asymmetryCm = 0.0;
};

// Samples `shape` along its own local y axis (the same TGeoShape
// returned by TGeoVolume::GetShape(), evaluated in ITS OWN local frame -
// no placement matrix applied, exactly like the existing z-half-extent
// lookup in GenMagneticField.cc's SetupMDTMagneticField) in `stepCm`
// increments outward from its center, on both the +y and -y sides,
// looking for the solid -> gap -> solid pattern a slit cuts into an
// otherwise solid block.
MagnetSlitProbeResult ProbeMagnetSlit(TGeoShape* shape, double stepCm = 0.1);

} // namespace FASER

#endif
