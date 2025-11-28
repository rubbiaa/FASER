#ifndef _GENMAGNETICFIELD_HH_
#define _GENMAGNETICFIELD_HH_ 1

#include <iostream>
#include <TVector3.h>

#include <AbsBField.h>

//////////(o^o)///////////
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
//////////(o^o)///////////

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

    //////////(o^o)///////////
    // Geometry-driven field control: field off at tracking stations, on in-between
    // Units for all z values stored here are cm to match GenFit conventions
    void BuildGeometryMaps(int verbose = 0);
    bool IsInMagnet(double z_cm) const; // true if z is within a known magnet z-range
    bool IsNearStation(double z_cm, double eps_cm = 0.1) const; // within eps of a SciFi layer
    bool HasAnyMagnet() const { return !magnet_z_ranges_cm_.empty(); }
    // Tunables
    double station_dead_half_thickness_cm_ = 0.1; // 1 mm dead band around layer

    // Event-driven station z positions (cm). If provided, prefer these to gate the field
    void SetEventStationZsCm(const std::vector<double>& zs_cm) {
        event_station_z_cm_ = zs_cm;
        std::sort(event_station_z_cm_.begin(), event_station_z_cm_.end());
    }
    const std::vector<double>& GetEventStationZsCm() const { return event_station_z_cm_; }
    //////////(o^o)///////////

private:
    std::vector<std::pair<double,double>> magnet_z_ranges_cm_;
    std::vector<double> scifi_layer_z_cm_;
    //////////(o^o)///////////
    std::vector<double> event_station_z_cm_;
    //////////(o^o)///////////
    TVector3 get(const TVector3 &position) const override
    {
        // Mapping magnetic field here
        //std::cout << "Magnetic field at position: "
        //           << "x=" << position.X() << ", y=" << position.Y() << ", z=" << position.Z() << std::endl;

        // ////////    ///////////
        // Determine if magnetic field should be active at this z (cm)
        const double z_cm = position.Z();
        bool fieldOnHere = false;
        // Prefer per-event station z envelope when available
        if (!event_station_z_cm_.empty()) {
            const double zmin = event_station_z_cm_.front();
            const double zmax = event_station_z_cm_.back();
            if (z_cm > zmin && z_cm < zmax) fieldOnHere = true;
            // kill field near any provided station planes
            for (double zl : event_station_z_cm_) { if (std::abs(zl - z_cm) <= station_dead_half_thickness_cm_) { fieldOnHere = false; break; } }
        } else if (!magnet_z_ranges_cm_.empty()) {
            // Use explicit magnet volume z-ranges from geometry
            for (const auto& rng : magnet_z_ranges_cm_) {
                if (z_cm > rng.first && z_cm < rng.second) { fieldOnHere = true; break; }
            }
            // Enforce no field in a dead band around any station plane
            if (IsNearStation(z_cm, station_dead_half_thickness_cm_)) fieldOnHere = false;
        } else if (!scifi_layer_z_cm_.empty()) {
            // Fallback: define field regions in-between SciFi layers
            const double zmin = *std::min_element(scifi_layer_z_cm_.begin(), scifi_layer_z_cm_.end());
            const double zmax = *std::max_element(scifi_layer_z_cm_.begin(), scifi_layer_z_cm_.end());
            if (z_cm >= zmin && z_cm <= zmax) fieldOnHere = true;
            // Enforce no field near layers
            if (IsNearStation(z_cm, station_dead_half_thickness_cm_)) fieldOnHere = false;
        }

        if (fieldOnHere) {
            //////////(o^o)///////////
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
        // No field at stations or outside magnet regions
        return TVector3(0, 1e-3, 0);
    }
};
#endif