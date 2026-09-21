#ifndef _GENMAGNETICFIELD_HH_
#define _GENMAGNETICFIELD_HH_ 1

#include <iostream>
#include <TVector3.h>

#include <AbsBField.h>

#include <TGeoManager.h>
#include <TGeoNode.h>
#include <string>

#include "MuonSpectrometerField.hh"

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
    double rearMuSpec_LOS_shiftZ = 0.0; // in cm
    double rearMuSpec_tilt_deg = 0.0;   // rotation around y-axis in degrees
    void SetRearMuSpectShift(double shiftX, double shiftY, double shiftZ = 0.0, double tilt_deg = 0.0) { 
        rearMuSpec_LOS_shiftX = shiftX; 
        rearMuSpec_LOS_shiftY = shiftY; 
        rearMuSpec_LOS_shiftZ = shiftZ;
        rearMuSpec_tilt_deg = tilt_deg;
    }

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

    // Precomputed MDT magnet z-ranges in GLOBAL cm (zmin,zmax per magnet).
    // Set once per event from the cached geometry so that get() never has to
    // call gGeoManager->FindNode() while GenFit's RK stepper is mid-propagation:
    // GenFit's own TGeoMaterialInterface concurrently drives gGeoManager's
    // (stateful, single) navigator for material stepping, and an independent
    // FindNode() call from inside the field functor can desynchronize that
    // navigator's current-node bookkeeping, which manifests as the Kalman fit
    // diverging/non-converging mid-track.
    void SetMDTMagnetZRangesCm(const std::vector<std::pair<double,double>>& ranges_cm) {
        mdt_magnet_z_ranges_cm_ = ranges_cm;
    }

private:
    std::vector<std::pair<double,double>> magnet_z_ranges_cm_;
    std::vector<double> scifi_layer_z_cm_;
    //////////(o^o)///////////
    std::vector<double> event_station_z_cm_;
    std::vector<std::pair<double,double>> mdt_magnet_z_ranges_cm_;
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
        /////////////////////////////////////////
        // Both branches below describe the SAME physical field model -
        // FASER::ComputeMuonSpectrometerField (CoreUtils/
        // MuonSpectrometerField.hh), also used by FASERG4's
        // MuonMagneticField::GetFieldValue - just at two different
        // magnet placements ("Magnet" here vs "MDTMagnet" just below),
        // exactly like FASERG4's DetectorConstruction places the same
        // MuonMagneticField class at both. rearMuSpec_tilt_deg is set to
        // -fTiltAngleY (deg) (see TPORecoEvent.cc), so recover the
        // detector-frame tilt in radians once, matching FASERG4's own
        // convention.
        const double tiltAngleRad = -rearMuSpec_tilt_deg * M_PI / 180.0;

        // Dedicated MDT magnet field.
        // Uses precomputed global z-ranges (set once per event via
        // SetMDTMagnetZRangesCm) instead of a live gGeoManager->FindNode()
        // lookup: see the comment on SetMDTMagnetZRangesCm for why a live
        // lookup here is unsafe during GenFit propagation.
        for (const auto& rng : mdt_magnet_z_ranges_cm_) {
            if (z_cm > rng.first && z_cm < rng.second) {
                // Rotation around Y leaves Y unchanged, so global Y
                // already equals the tilted assembly's local Y -- no
                // correction needed beyond the shift.
                const double y_local_cm = position.Y() - rearMuSpec_LOS_shiftY;
                FASER::MuonSpectrometerFieldParams params;
                params.slitPositionCm = slitposition;
                params.tiltAngleRad = tiltAngleRad;
                // Smooth the central/outer field boundary with a 1 cm linear
                // ramp (+/-0.5 cm) to avoid a step-function discontinuity
                // that can destabilise the Runge-Kutta integrator for
                // tracks near |y_local|=slitposition. Numerical
                // accommodation only - see rampHalfWidthCm's doc comment.
                params.rampHalfWidthCm = 0.5;
                double Bx_kG, By_kG, Bz_kG;
                if (FASER::ComputeMuonSpectrometerField(y_local_cm, params, Bx_kG, By_kG, Bz_kG))
                    return TVector3(Bx_kG, By_kG, Bz_kG);
                return TVector3(0.0, 0.0, 0.0);
            }
        }

        if (fieldOnHere) {
            // Rotation around Y leaves Y unchanged, so global Y already
            // equals the tilted assembly's local Y -- no correction
            // needed beyond the shift (see the MDT branch above, and
            // FASERG4::MuonMagneticField::GetFieldValue, for the same
            // simplification).
            const double y_local_cm = position.Y() - rearMuSpec_LOS_shiftY;
            FASER::MuonSpectrometerFieldParams params;
            params.slitPositionCm = slitposition;
            params.tiltAngleRad = tiltAngleRad;
            // rampHalfWidthCm left at its default (0): matches FASERG4's
            // exact hard step at the slit boundary for this magnet type.
            double Bx_kG, By_kG, Bz_kG;
            if (FASER::ComputeMuonSpectrometerField(y_local_cm, params, Bx_kG, By_kG, Bz_kG))
                return TVector3(Bx_kG, By_kG, Bz_kG);
        }
        // No field at stations or outside magnet regions
        return TVector3(0, 1e-3, 0);
    }
};
#endif