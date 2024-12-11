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

// genfit
#include <RKTrackRep.h>
#include <Track.h>

/// @brief Holds a reconstructed track from the precise pixel tracker
class TTKTrack : public TObject {
public:

   struct TRACKHIT {
        long ID;
        int type; // 0 - scintillator, 1 - tracker
        ROOT::Math::XYZVector point;
        float eDeposit;
    };

    std::vector<TRACKHIT> tkhit;
    TVector3 centroid;
    TVector3 direction;
    double SSR;
    int trackID;                         // the unique track ID (<0 if not yet assigned)
    int vertexID;                        //the vertex ID the track belongs to (<0 if not associated to a vertex)

    genfit::Track *fitTrack;                //! the genfit track

    /// @brief Function to sort hits by the Z coordinate of the hits
    void SortHitsByZ() {
        std::sort(tkhit.begin(), tkhit.end(), [](const TRACKHIT& a, const TRACKHIT& b) {
            return a.point.z() < b.point.z();
        });
    }

    TTKTrack() : TObject(), fitTrack(0), vertexID(-1), trackID(-1) { 
        //std::cout << "TTKTrack::TTKTrack - constructor . " << this << std::endl; 
        };
    TTKTrack(const TTKTrack &t);
    virtual ~TTKTrack() { 
//        std::cout << "TTKTrack::~TTKTrack - destructor . " << this << std::endl;
        if(fitTrack!=nullptr) delete fitTrack;
    };

    /// @brief Compute the distance between a point and a line
    double pointLineDistance(const ROOT::Math::XYZVector& point, const TVector3& direction, const TVector3& centroid);

    /// @brief Merge two tracks
    /// @param track2 - the track to merge
    void MergeTracks(TTKTrack &track2);

    /// @brief The direction defined by two first hits of the track
    /// @return the direction
    TVector3 direction2Hits();

    /// @brief Fit line through hits
    /// @param centroid - returns the centroid of the track
    /// @return the direction of the track
    TVector3 fitLineThroughHits(TVector3& centroid);

    /// @brief Use GenFit to fit the track
    void GenFitTrackFit();

    /// @brief Extrapolate the track to a given Z coordinate
    TVector3 extrapolateTracktoZ(double z, int &failed);

    /// @brief Get the minimum and maximum module of the track
    void GetMinMaxModule(int &minLayer, int &maxLayer);

    /// @brief Update the position with the fitted values
    void UpdateFittedPosition();

    /// @brief Dump the track information
    void Dump(int verbose = 1) const;

    ClassDef(TTKTrack, 1)
};


#endif
