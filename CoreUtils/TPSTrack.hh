#ifndef _TPSTRACK_
#define _TPSTRACK_ 1

#include <TObject.h>
#include <vector>
#include <algorithm>  // For std::find_if
#include <Math/Vector3D.h>
#include <TVector3.h>

#include "TcalEvent.hh"

/// @brief Holds a reconstructed track from the plastic scintillator detector
class TPSTrack : public TObject {
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

    TPSTrack() {};

    /// @brief Check if ID belongs to track
    struct TRACKHIT *VoxelBelongstoTrack(long ID)
    {
        auto it = std::find_if(tkhit.begin(), tkhit.end(),
                               [ID](const auto &hit)
                               { return hit.ID == ID; });
        // If an element is found, return it; otherwise return nullptr
        return (it != tkhit.end()) ? &(*it) : nullptr;
    }

    /// @brief Function to sort hits by the Z coordinate of the hits
    inline void SortHitsByZ() {
        std::sort(tkhit.begin(), tkhit.end(), [](const TRACKHIT& a, const TRACKHIT& b) {
            return a.point.z() < b.point.z();
        });
    }

    /// @brief The direction defined by two first hits of the track
    /// @return the direction
    TVector3 direction2Hits(TVector3& centroid);

    /// @brief Fit line through hits
    /// @param centroid - returns the centroid of the track
    /// @return the direction of the track
    TVector3 fitLineThroughHits(TVector3& centroid);

    void Dump() const {
        std::cout << "PSTrack: " << tkhit.size() << " hits ";
        std::cout << "direction " << direction.x() << " " << direction.y() << " " << direction.z() << std::endl;
        for (const auto &it : tkhit) {
            std::cout << "   " << it.ID << " eDeposit: " << it.eDeposit << std::endl;
        }        
    }
    ClassDef(TPSTrack, 1)
};


#endif
