#ifndef _GENMAGNETICFIELD_HH_
#define _GENMAGNETICFIELD_HH_ 1

#include <TVector3.h>

#include <AbsBField.h>

class GenMagneticField : public genfit::AbsBField {
public:
    GenMagneticField() = default;
    virtual ~GenMagneticField() = default;

    TVector3 get(const TVector3& position) const override {
        // Implement the magnetic field calculation here
        std::cout << "Calculating magnetic field at position: "
                  << "x=" << position.X() << ", y=" << position.Y() << ", z=" << position.Z() << std::endl;

        double slitposition = 25; // in cm
        if(std::abs(position.Y()) >= slitposition && std::abs(position.Y()) <= 2 * slitposition) {
            // Top or Bottom: +1.5 T
            return TVector3(+15.0, 0, 0); // 15 kGauss
        } else if (std::abs(position.Y()) < slitposition) {
            // Middle: -1.5 T
            return TVector3(-15.0, 0, 0); // -15 kGauss
        }
    return TVector3(0, 0, 0); // No field outside the magnet
    }

};
#endif
