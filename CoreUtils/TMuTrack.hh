#ifndef _TMUTRACK_
#define _TMUTRACK_ 1

#include <TObject.h>
#include <Math/Vector3D.h>
#include <TVector3.h>

// genfit
#include <RKTrackRep.h>
#include <Track.h>

class TMuTrack : public TObject {
public:
    int ftrackID;       // the unique track ID
    int fPDG;         // the PDG code of the particle
    float fcharge;    // the charge of the particle
    std::vector<ROOT::Math::XYZVector> fpos;
    std::vector<int> layerID; // the layer index in the muon spectrometer (0-43)

    genfit::Track *fitTrack;                //! the genfit track
    double fpx, fpy, fpz, fp; // fitted momentum at the first point
    double fpErr;            // error on fitted momentum at the first point
    double fchi2;              // chi2 of the fit
    int fnDoF;                // nDoF of the fit
    double fpval;              // p-value of the fit
    double fipErr;            // error on fitted inverse momentum at the first point

    TMuTrack() : ftrackID(-1), fcharge(0), fitTrack(nullptr) {}
    virtual ~TMuTrack() {
        if (fitTrack) {
            delete fitTrack;
            fitTrack = nullptr;
        }
    }

    void Dump(int verbose = 1) const
    {
        std::cout << "TMuTrack: TrackID: " << ftrackID << " PDG: " << fPDG << " Charge: " << fcharge << " npoints: " << fpos.size() << std::endl;
        if (verbose > 1)
        {
            for (size_t i = 0; i < fpos.size(); i++)
            {
                std::cout << "  Point: " << i << " LayerID: " << layerID[i] << " Pos: (" << fpos[i].x() << "," << fpos[i].y() << "," << fpos[i].z() << ")" << std::endl;
            }
        }
    };

    /// @brief Use GenFit to fit the track
    void GenFitTrackFit(int verbose, double detectorResolutionPSmm);

    ClassDef(TMuTrack,2) // A reconstructed track in the muon spectrometer
};

#endif
