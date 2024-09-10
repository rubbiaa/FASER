#ifndef _TTKTRACK_
#define _TTKTRACK_ 1

#include <TObject.h>
#include <vector>
#include <TH2D.h>
#include <Math/Vector3D.h>
#include <TVector3.h>

#include <TDatabasePDG.h>

#include "TcalEvent.hh"
#include "TPOEvent.hh"
#include "TPSCluster.hh"

/// @brief Holds a reconstructed track from the precise pixel tracker
class TTKTrack : public TObject {
public:

   struct TRACKHIT {
        long ID;
        ROOT::Math::XYZVector point;
        float eDeposit;
    };

    std::vector<TRACKHIT> tkhit;
    TVector3 centroid;
    TVector3 direction;
    double SSR;

    /// @brief Function to sort hits by the Z coordinate of the hits
    void SortHitsByZ() {
        std::sort(tkhit.begin(), tkhit.end(), [](const TRACKHIT& a, const TRACKHIT& b) {
            return a.point.z() < b.point.z();
        });
    }

    TTKTrack() {};

    /// @brief The direction defined by two first hits of the track
    /// @return the direction
    TVector3 direction2Hits();

    /// @brief Fit line through hits
    /// @param centroid - returns the centroid of the track
    /// @return the direction of the track
    TVector3 fitLineThroughHits(TVector3& centroid);

    ClassDef(TTKTrack, 1)
};


#endif
