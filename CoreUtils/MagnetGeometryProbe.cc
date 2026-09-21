#include "MagnetGeometryProbe.hh"

#include <TGeoShape.h>

#include <algorithm>
#include <cmath>

namespace FASER {

namespace {

// Walks from the shape's center towards +y (direction=+1) or -y
// (direction=-1), starting from a point already known to be solid,
// looking for a single solid -> gap -> solid transition. On success,
// fills [gapStartCm, gapEndCm] with the (unsigned, center-relative)
// distances where the gap begins and ends and returns true.
bool FindGapAlongDirection(TGeoShape* shape, double xCenter, double yCenter,
                           double zCenter, double yHalfExtent, double stepCm,
                           int direction, double& gapStartCm, double& gapEndCm)
{
    double probe[3] = {xCenter, yCenter, zCenter};

    double ySolidEdge = yCenter;
    double y = yCenter;
    bool foundGapStart = false;
    while (std::abs(y - yCenter) < yHalfExtent) {
        y += direction * stepCm;
        probe[1] = y;
        if (!shape->Contains(probe)) { foundGapStart = true; break; }
        ySolidEdge = y;
    }
    if (!foundGapStart) return false; // stayed solid all the way to the block's own edge

    double yGapEnd = y;
    bool foundGapEnd = false;
    while (std::abs(y - yCenter) < yHalfExtent) {
        y += direction * stepCm;
        probe[1] = y;
        if (shape->Contains(probe)) { yGapEnd = y; foundGapEnd = true; break; }
    }
    if (!foundGapEnd) return false; // gap never closes before the block's own edge

    gapStartCm = std::abs(ySolidEdge - yCenter);
    gapEndCm   = std::abs(yGapEnd - yCenter);
    if (gapEndCm < gapStartCm) std::swap(gapStartCm, gapEndCm);
    return true;
}

} // namespace

MagnetSlitProbeResult ProbeMagnetSlit(TGeoShape* shape, double stepCm)
{
    MagnetSlitProbeResult result;
    if (!shape || stepCm <= 0.0) return result;

    // TGeoShape::GetAxisRange's axis indices are 1=X, 2=Y, 3=Z (exactly
    // as GenMagneticField.cc's SetupMDTMagneticField already uses axis 3
    // for the z half-extent) - this is the shape's own bounding box, in
    // its own local frame, before any placement matrix.
    double xlo, xhi, ylo, yhi, zlo, zhi;
    shape->GetAxisRange(1, xlo, xhi);
    shape->GetAxisRange(2, ylo, yhi);
    shape->GetAxisRange(3, zlo, zhi);

    const double xCenter = 0.5 * (xlo + xhi);
    const double yCenter = 0.5 * (ylo + yhi);
    const double zCenter = 0.5 * (zlo + zhi);
    const double yHalfExtent = 0.5 * (yhi - ylo);
    if (yHalfExtent <= 0.0) return result;

    // The block's own center must be solid, or this isn't the "solid
    // block with a slit cut into each side" shape this probe expects.
    double probeCenter[3] = {xCenter, yCenter, zCenter};
    if (!shape->Contains(probeCenter)) return result;

    double posStart, posEnd, negStart, negEnd;
    const bool foundPos = FindGapAlongDirection(shape, xCenter, yCenter, zCenter,
                                                 yHalfExtent, stepCm, +1, posStart, posEnd);
    const bool foundNeg = FindGapAlongDirection(shape, xCenter, yCenter, zCenter,
                                                 yHalfExtent, stepCm, -1, negStart, negEnd);
    if (!foundPos || !foundNeg) return result;

    const double posCenter = 0.5 * (posStart + posEnd);
    const double negCenter = 0.5 * (negStart + negEnd);

    result.ok = true;
    result.slitPositionCm = 0.5 * (posCenter + negCenter);
    result.slitHalfWidthCm = 0.25 * ((posEnd - posStart) + (negEnd - negStart));
    result.blockHalfExtentYCm = yHalfExtent;
    result.asymmetryCm = std::abs(posCenter - negCenter);
    return result;
}

} // namespace FASER
