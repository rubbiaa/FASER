#include "MuonSpectrometerField.hh"
#include <cmath>

namespace FASER {

bool ComputeMuonSpectrometerField(double y_local_cm,
                                   const MuonSpectrometerFieldParams& params,
                                   double& Bx_kG, double& By_kG, double& Bz_kG)
{
    Bx_kG = By_kG = Bz_kG = 0.0;

    const double y_abs = std::abs(y_local_cm);
    const double lo = params.slitPositionCm - params.rampHalfWidthCm;
    const double hi = params.slitPositionCm + params.rampHalfWidthCm;

    double Blocal_kG;
    if (y_abs < lo) {
        // Middle region.
        Blocal_kG = -params.fieldMagnitudeKG;
    } else if (y_abs < hi) {
        // Linear ramp across the slit boundary. When rampHalfWidthCm==0
        // (FASERG4's default), lo==hi==slitPositionCm: the branch above
        // already covers y_abs<lo, and the branch below covers
        // y_abs>=hi, so this condition can never be true and the exact
        // hard step G4 simulates is reproduced bit-for-bit.
        const double t = (y_abs - lo) / (hi - lo);
        Blocal_kG = -params.fieldMagnitudeKG + 2.0 * params.fieldMagnitudeKG * t;
    } else if (y_abs <= 2.0 * params.slitPositionCm) {
        // Top/bottom region.
        Blocal_kG = params.fieldMagnitudeKG;
    } else {
        // Outside the modelled envelope entirely.
        return false;
    }

    // The field points along the assembly's local +x axis; rotate that
    // local vector (Blocal, 0, 0) into the global frame by the SAME
    // rotation used to place the assembly, so B stays aligned with the
    // iron even when the assembly is tilted.
    Bx_kG = Blocal_kG * std::cos(params.tiltAngleRad);
    Bz_kG = Blocal_kG * std::sin(params.tiltAngleRad);
    return true;
}

} // namespace FASER
