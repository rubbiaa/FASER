#ifndef _TMUTRACK_
#define _TMUTRACK_ 1

#include <TObject.h>
#include <Math/Vector3D.h>
#include <TVector3.h>

// genfit
#include <RKTrackRep.h>
#include <Track.h>


struct MDTMeas {
    // x_mm is not measured
    double x_mm = 0.0;
    double y_mm = 0.0;
    double z_mm = 0.0;

    double wireX_mm = 0.0;
    double wireY_mm = 0.0;
    double wireZ_mm = 0.0;
    double r_meas_mm = 0.0;
    double r_true_mm = 0.0;

    int stationID = 0;
    int planeID = 0;
    int tubeID = 0;
    int side = 0;

    // localX_mm is arbitrary; localY/localZ define the MDT 1D measurement.
    double localX_mm = 0.0;
    double localY_mm = 0.0;
    double localZ_mm = 0.0;
    // GLOBAL directions of MDT-local axes.
    // u = measured direction = MDT-local Y.
    // v = unmeasured wire direction = MDT-local X.
    // w = downstream direction = MDT-local Z.
    double uX = 0.0, uY = 1.0, uZ = 0.0;
    double vX = 1.0, vY = 0.0, vZ = 0.0;
    double wX = 0.0, wY = 0.0, wZ = 1.0;
};


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
    double fQOverP;            // fitted inverse momentum at the first point
    double fQOverPErr;         // error on fitted inverse momentum at the first point
    
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
    bool GenFitMDTFit(const std::vector<MDTMeas>& meas, int pdg, double seedMomentumGeV = 10.0, int verbose = 0);
    
    
    // ////////    ///////////
    /// Fast circle fit (Taubin-style) in the bending plane (y-z) under Bx field
    /// Fills fpx,fpy,fpz,fp as a lightweight alternative to GenFit. Units:
    /// - Inputs: positions in mm (fpos)
    /// - Field queried in kGauss via GenFit FieldManager
    /// - Outputs: momentum components in GeV/c
    void CircleFitTaubin(int verbose, double detectorResolutionPSmm);
    // ////////    ///////////

    ClassDef(TMuTrack,2) // A reconstructed track in the muon spectrometer
};

#endif
