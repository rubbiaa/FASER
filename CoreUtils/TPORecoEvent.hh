#ifndef _TPORECOEVENT_
#define _TPORECOEVENT_ 1

#include <TObject.h>
#include <vector>
#include <TH2D.h>
#include <Math/Vector3D.h>
#include <TVector3.h>
#include <TMatrixDSym.h>
#include <Math/SMatrix.h>

#include <TDatabasePDG.h>

#include "TcalEvent.hh"
#include "TPOEvent.hh"
#include "TPSCluster.hh"
#include "TTKTrack.hh"
#include "TPSTrack.hh"

/// @brief TPORec holds a reconstructed particle object
class TPORec : public TObject {
public:

   struct CALENERGIES {
        double em;             // in GeV
        double had;            // in GeV
        double Ecompensated;   // in GeV
        ROOT::Math::XYZVector cog;
        ROOT::Math::XYZVector Eflow;
    };

    int POID;                          // the primary track in POEvent (the index)
    std::vector<int> fGEANTTrackIDs;   //! all the geant track id that belong to this PORec
    std::vector<DigitizedTrack*> DTs;  //! all the DigitizedTracks that belong to this POREC
    std::vector<struct CALENERGIES> fEnergiesCogs; //! the energies and COG of each Digitized track

    struct CALENERGIES fTotal;         // the cumulative energies for the primary

    struct TRACKHIT {
        long ID;
        ROOT::Math::XYZVector point;
        float eDeposit;
    };

    struct TRACK {
        std::vector<TRACKHIT> tkhit;
        TVector3 centroid;
        TVector3 direction;
        double SSR;
    };

    std::vector<TRACK> fTracks;       // all the reconstructed tracks that belong to this PORec

    // constructors & destructors
    TPORec() = default;
    TPORec(int id) : POID(id) {};
    virtual ~TPORec() = default;

    // computed quantities
    double TotalEvis() { return fTotal.Ecompensated; };   // Total visible energy
    double TotalET() { return sqrt(fTotal.Eflow.Perp2()); };// Total transverse energy

    ClassDef(TPORec,1)
};

class TPORecoEvent : public TObject {
private:

    /// @brief Private function to compute energies and COG belonging from a Digitized Track hits.
    struct TPORec::CALENERGIES computeEnergiesAndCOG(DigitizedTrack *dt);   //! (no ROOT I/O output)

    /// @brief The vector that holds all the PORec (Reconstructed POs) in the event
    std::vector<class TPORec*> fPORecs;                                     //! (no ROOT I/O output)

    TPORec *fPOFullEvent = nullptr;                   // the kinematics of the full event (TRUTH)
    TPORec *fPOFullRecoEvent = nullptr;                   // the kinematics of the full event (RECO)

    TcalEvent* fTcalEvent;                            //! Reference to the TCAL event
    TPOEvent* fTPOEvent;                              // Reference to the TPOEvent

// static    TVector3 fitLineThroughPoints(const struct TPORec::TRACK &track, TVector3& centroid);
public:

   /// @brief A hit in a two dimension plastic scintillator view
    struct PSHIT2D {
        float electromagneticity;       // =0 if hadronic, =1 if electromagnetic
        int  ntracks;                   // number of overlapping tracks
        float Edeposited;               // deposited energy
    };

    std::map<long, PSHIT2D> PShitmapX;      //! the X-Z view hit map
    std::map<long, PSHIT2D> PShitmapY;      //! the Y-Z view hit map
    std::map<int, std::map<long, PSHIT2D>> PShitmapsZ; //! the X-Y views per layer

    /// @brief Return the x,y,z position of a 2D hit (one coordinate x, or y should always be ignored)
    void pshit2d_position(long ID, double &fix, double &fiy, double &fiz);

    TH2D* xviewPS = nullptr;                          //! 2Dview scintillator X-Z
    TH2D* yviewPS = nullptr;                          //! 2Dview scintillator Y-Z
    std::vector<TH2D*> zviewPS;                                //! 2Dview scintillator X-Z view for 50 planes
    TH2D* xviewPS_em = nullptr;                          //! 2Dview scintillator X-Z
    TH2D* yviewPS_em = nullptr;                          //! 2Dview scintillator Y-Z
    TH2D* xviewPS_had = nullptr;                          //! 2Dview scintillator X-Z
    TH2D* yviewPS_had = nullptr;                          //! 2Dview scintillator Y-Z0
    TH2D* xviewPS_eldepo = nullptr;                          //! 2Dview scintillator X-Z
    TH2D* yviewPS_eldepo = nullptr;                          //! 2Dview scintillator Y-Z

    std::vector<TPSCluster> PSClustersX;                // 2Dview clusters XZ
    std::vector<TPSCluster> PSClustersY;                // 2Dview clusters YZ
    size_t n_psclustersX() { return PSClustersX.size(); };    // number of reconstructed cluster in XZ view
    size_t n_psclustersY() { return PSClustersY.size(); };    // number of reconstructed cluster in YZ view

    struct PSVOXEL3D {
        long ID;
        float RawEnergy;            // MeV
        bool ghost;
        bool member_of_TKtrack;       // is this voxel part of a TKtrack
    };

    struct Voxel {
        float value;
        Voxel() : value(0) {}
    };

    std::map<long, struct PSVOXEL3D> PSvoxelmap;

    /// @brief All reconstruced TKTracks in event
    std::vector<TTKTrack> fTKTracks;

    struct TTKVertex {
        int vertexID;
        ROOT::Math::XYZVector position;
        TMatrixDSym covariance;
        int ntracks;
        int ndof;
        double chi2;
        std::vector<int> trackIDs;  // the track IDs associated to this vertex

        // Constructor
        TTKVertex()
           : vertexID(-1), covariance(3), ntracks(0), ndof(0), chi2(0.0) {} // Initialize covariance with size 3
    };

    std::vector<TTKVertex> fTKVertices;

    /// @brief All reconstruced TPSTracks in event
    std::vector<TPSTrack> fPSTracks;

    struct FASERCAL {
        int ModuleID;
        double EDeposit;   // total energy in FaserCal module in GeV
    };
    std::vector<struct FASERCAL> faserCals;

    /// @brief Structure to hold rear calorimeter and mutag deposited energies
    struct REARCALS {
        double rearCalDeposit;   // total energy in RearCal in GeV
        double rearHCalDeposit;  // total energy in hCal in GeV
        double rearMuCalDeposit; // total energy in muTag in MeV
        std::vector<struct TcalEvent::REARCALDEPOSIT> rearCalModule; // individual module deposits
        std::vector<struct TcalEvent::REARCALDEPOSIT> rearHCalModule; // individual module deposits
    };
    struct REARCALS rearCals;

    // @brief A copy of the geometry originally stored in TCalEvent
    struct TcalEvent::GEOM_DETECTOR geom_detector;

    TPORecoEvent();
    TPORecoEvent(TcalEvent* c, TPOEvent* p);
    virtual ~TPORecoEvent();

    int verbose = 0;                            //! controls amount of debug information
    bool multiThread = true;                   //! controls if multi-threading is used

    struct RECOCONFIG {
        // resolution fudge factor used in TTKTrack Genfit fitting for PS voxel hits
        double psvoxel_fudge_factor;

        // energy compensation factors for the PS voxel hits
        double alpha;
        double beta;

        double findpattern_max_hit_layers;
        double findpattern_dist_min_cut;
        double findpattern_parallel_cut;
        double findpattern_mindZ_fudge;
    //    double findpattern_cut_SSR_merge;
        double findpattern_parallel_cut_merge;

        double extendtracks_closest_voxel_cut;
        double extendtracks_dist2_perp_voxel_cut;

        // genfit track minimum pVal 
        double genfit_min_pVal;
        double genfit_min_pMom;

        size_t findvtx_cut_max_trk;
        double findvtx_chi2ndf_cut;
        double findvtx_trk_dist_cut;
        double findvtx_merge_dist_cut;

        double clusters_threshold_2dhit;
        double clusters_eps;
        double clusters_minPts;
        double clusters_threshold_cluster;

        double PS3D_nvox_max_after_iteration;
        double PS3D_total_score_min_break;
        double PS3D_ehit_threshold;
        double PS3D_evox_threshold;
        int    PS3D_nvox_per_layer_max;

        int PSFilter_max_number_track_seeds;
        double PSFilter_closest_voxel_cut;
        double PSFilter_parallel_cut;
        double PSFilter_mindZcut;

    } recoConfig;

    /// @brief Reconstruct the FASERG4 simulated event to the PORec
    void ReconstructTruth();

    /// @brief Full track reconstruction based on pixel and scintillator hits
    void TrackReconstruct();

    /// @brief Pattern recognition for the tracker
    void FindPatternTracks();

    /// @brief Extend the tracks using voxels in the plastic scintillator
    void ExtendTracks();

    /// @brief Find the tracks vertices using GenFit2
    void FindTrackVertices();

    /// @brief Reconstruct all the tracks associated to the PORec (call this after Reconstruct)
    void TrackReconstructTruth();

    /// @brief Reconstruct the 2D plastic scintillator views XZ and YZ
    void Reconstruct2DViewsPS();

    /// @brief Reconstruct all 2D clusters for the xz and the yz views (view=0 for XZ, and view=1 for YZ)
    void ReconstructClusters(int view);

    /// @brief Reconstruct 3D voxels from 2D views in plastic scintillator
    void Reconstruct3DPS(int maxIter = 150);
    void Reconstruct3DPS_2(int maxIter = 150);
    void reconstruct3DPS_module(int maxIter, int imodule, std::vector<std::vector<std::vector<Voxel>>> &V,
        std::vector<std::vector<float>> &XZ, std::vector<std::vector<float>> &YZ, std::vector<std::vector<std::vector<float>>> &XY,
        std::vector<int>& nvox_per_layer, int nvox_per_layer_max);
    void Reconstruct3DPS_Eflow();

    /// @brief Recontruct particle tracks from 3D PS voxels
    void PSVoxelParticleFilter();

    /// @brief Reconstruct FASERCAL and rear calorimeters and rear mu tag
    void ReconstructRearCals();

    /// @brief Dump PORecs to the screen
    void Dump();

    /// @brief Dump reconstructed TTKTracks to the screen
    void DumpReconstructedTracks();

    /// @brief Returns the vector of Reconstructed POs
    std::vector<class TPORec*> GetPORecs() { return fPORecs;};

    /// @brief Returns the kinematic quantities of the full TRUTH event
    /// @return TPORec of the full event
    TPORec *GetPOFullEvent() { return fPOFullEvent; };

    /// @brief Returns the kinematic quantities of the full event (RECONSTRUCTED)
    /// @return TPORec of the full event
    TPORec *GetPOFullRecoEvent() { return fPOFullRecoEvent; };

    /// @brief Returns the truth MC information
    /// @return TPOEvent pointer of the MC truth event
    TPOEvent *GetPOEvent() { return fTPOEvent; };

    // Reconstructed event summary variables

    /// @brief Number of charged particles at the primary
    int primary_n_charged;

    /// @brief Number of scintillator hits left by tau
    int nhits_tau;

    /// @brief Number of hits found in first tracker layers
    int nhits_tracker_first;

    /// @brief Fill the 2D (x-z) and (y-z) views of the Scintillator detector
    void Fill2DViewsPS();
    TH2D* Get2DViewXPS() { return xviewPS; };
    TH2D* Get2DViewYPS() { return yviewPS; };

    ClassDef(TPORecoEvent,3)
};

#endif
