#ifndef _GENMAGNETICFIELD_HH_
#define _GENMAGNETICFIELD_HH_ 1

#include <TVector3.h>

#include <AbsBField.h>

class GenMagneticField : public genfit::AbsBField {
public:
    GenMagneticField() = default;
    virtual ~GenMagneticField() = default;

    double slitposition = 25; // position of the slit along y in cm
    void SetSlitPosition(double pos) { slitposition = pos; }

    double rearMuSpectLocZ = 0.0; // in cm
    double rearMuSpectSizeZ = 0.0; // in cm
    void SetRearMuSpectGeometry(double locZ, double sizeZ) { rearMuSpectLocZ = locZ; rearMuSpectSizeZ = sizeZ; }

    double rearMuSpec_LOS_shiftX = 0.0; // in cm
    double rearMuSpec_LOS_shiftY = 0.0; // in cm
    void SetRearMuSpectShift(double shiftX, double shiftY) { rearMuSpec_LOS_shiftX = shiftX; rearMuSpec_LOS_shiftY = shiftY; }

    TVector3 get(const TVector3 &position) const override
    {
        // Implement the magnetic field calculation here
        // std::cout << "Calculating magnetic field at position: "
        //          << "x=" << position.X() << ", y=" << position.Y() << ", z=" << position.Z() << std::endl;

        // Check if the position is within the muon spectrometer region
        if (position.Z() >= rearMuSpectLocZ && position.Z() <= (rearMuSpectLocZ + rearMuSpectSizeZ))
        {
            // Adjust position for LOS shifts
            TVector3 localPos = position;
            localPos.SetX(position.X() - rearMuSpec_LOS_shiftX);
            localPos.SetY(position.Y() - rearMuSpec_LOS_shiftY);
            if (std::abs(localPos.Y()) >= slitposition && std::abs(localPos.Y()) <= 2 * slitposition)
            {
                // Top or Bottom: +1.5 T
                return TVector3(+15.0, 0, 0); // 15 kGauss
            }
            else if (std::abs(localPos.Y()) < slitposition)
            {
                // Middle: -1.5 T
                return TVector3(-15.0, 0, 0); // -15 kGauss
            }
        }
        std::cout << "Outside magnet region, no field." << std::endl;
        return TVector3(0, 1e-3, 0); // No field outside the magnet
    }
};
#endif
