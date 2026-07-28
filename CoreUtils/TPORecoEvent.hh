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
#include "TMuTrack.hh"
#include "GenMagneticField.hh"
#include "TPORecoEvent.hh"

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
    struct CALENERGIES fTotal_fasercal;        // total energy in FASER calorimeter (compensated)
    struct CALENERGIES fTotal_ecal;           // total calorimetric energy
    struct CALENERGIES fTotal_hcal;           // total hadronic calorimetric energy

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
    double TotalEvis() { return sqrt(fTotal.Eflow.Mag2()); };   // Total visible energy
    double TotalET() { return sqrt(fTotal.Eflow.Perp2()); };// Total transverse energy

    ClassDef(TPORec,1)
};

class TPORecoEvent : public TObject {
private:

    /// @brief Private function to compute energies and COG belonging from a Digitized Track hits.
    struct TPORec::CALENERGIES computeEnergiesAndCOG(DigitizedTrack *dt);   //! (no ROOT I/O output)

    //////////////////////////////////////////////////////////
    // FASERCAL DETECTOR RESPONSE FUNCTIONS
    // Convert energy deposits to photoelectrons (PE)
    // 
    // Fiber geometry: 11×48×48 voxels per layer
    // - Each voxel contributes light to 3 fibers (X, Y, Z directions)
    // - Fibers run line-by-line through voxels collecting light
    // - All fibers read out at positive (+) end with SiPMs
    // - Light attenuates with distance along fiber to readout
    //////////////////////////////////////////////////////////
    
    /// @brief Apply FASERCal detector response to convert energy deposits to PE
    /// @details Processes all tracks and voxel hits, computing PE for X/Y/Z fibers.
    ///          Includes optical crosstalk between neighboring voxels.
    /// @param tracks Input digitized tracks with energy deposits
    /// @return Vector of voxel PE responses (one per voxel with signal)
    std::vector<TcalEvent::FASERCALVOXELRESPONSE> applyFaserCalDetectorResponse(
        const std::vector<DigitizedTrack*>& tracks);
    
    /// @brief Compute PE for a single voxel from energy deposit
    /// @details Calculates light propagation from voxel to SiPM readout at +X, +Y, +Z ends.
    ///          Uses dual-component attenuation model (bulk + Rayleigh scattering).
    /// @param channelID Voxel channel ID (encoded ix, iy, iz, ilayer)
    /// @param energyDepositMeV Energy deposited in MeV
    /// @param fiberPE Output array for X,Y,Z fiber PE [out]
    /// @param applyStatistics If true (default), apply Poisson fluctuations (per
    ///        recoConfig.faserCal_applyPoissonStatistics) before returning. Pass
    ///        false to get the raw mean PE, e.g. when the caller needs to split
    ///        the mean across several voxels (crosstalk) before fluctuating each
    ///        target independently.
    /// @return Total PE (sum of all three fibers)
    double computeFaserCalDirectPE(long channelID, double energyDepositMeV,
                                   std::array<double, 3>& fiberPE,
                                   bool applyStatistics = true);
    
    /// @brief Apply optical leakage and accumulate PE with crosstalk
    /// @details Light can leak through voxel faces to neighbors (1% per face).
    ///          Voxels without direct hits can still receive PE from neighbors.
    /// @param voxelPEMap Map to accumulate total PE per voxel [in/out]
    /// @param voxelPEFibersMap Map to accumulate PE per fiber direction [in/out]
    /// @param channelID Central voxel ID where energy was deposited
    /// @param energyDepositMeV Energy deposited in central voxel (MeV)
    void accumulateFaserCalPEWithCrosstalk(
        std::map<long, double>& voxelPEMap,
        std::map<long, std::array<double, 3>>& voxelPEFibersMap,
        long channelID, double energyDepositMeV);
    
    /// @brief Helper function to decode FASERCal voxel ID
    /// @param id Encoded channel ID
    /// @param ix X voxel index [out]
    /// @param iy Y voxel index [out]
    /// @param iz Z voxel index within layer [out]
    /// @param ilayer Layer/module index [out]
    /// @return true if valid scintillator ID (hittype==0)
    bool decodeFaserCalScintID(long id, int& ix, int& iy, int& iz, int& ilayer);
    
    /// @brief Helper function to encode FASERCal voxel ID
    /// @return Encoded channel ID
    long encodeFaserCalScintID(int ix, int iy, int iz, int ilayer);
    
    /// @brief Helper function to check if voxel indices are valid
    /// @return true if indices within detector bounds
    bool isValidFaserCalVoxelIndex(int ix, int iy, int iz, int ilayer);

    /// @brief The vector that holds all the PORec (Reconstructed POs) in the event
    std::vector<class TPORec*> fPORecs;                                     // (keep ROOT I/O output)

    TPORec *fPOFullEvent = nullptr;                   // the kinematics of the full event (TRUTH)
    TPORec *fPOFullRecoEvent = nullptr;                   // the kinematics of the full event (RECO)

    TcalEvent* fTcalEvent;                            //! Reference to the TCAL event
    TPOEvent* fTPOEvent;                              // Reference to the TPOEvent

    GenMagneticField* fMagField = nullptr;         //! The magnetic field of the FASER detector

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

    TH2D* xviewPS = nullptr;                          // 2Dview scintillator X-Z
    TH2D* yviewPS = nullptr;                          // 2Dview scintillator Y-Z
    std::vector<TH2D*> zviewPS;                                // 2Dview scintillator X-Z view for all planes
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

    //////////////////////////////
    // added by Umut
    /// Added to access hits information in clusters
    const std::vector<TPSCluster>& GetPSClusters(int view) const {
        return (view == 0) ? PSClustersX : PSClustersY;
    }
    std::vector<TPSCluster>& GetPSClusters(int view) {
        return (view == 0) ? PSClustersX : PSClustersY;
    }
  
    std::vector<TPSCluster> PSClusters3D;
    const std::vector<TPSCluster>* GetPSClusters3D() const { return &PSClusters3D; }
    void Reconstruct3DClusters();
    ///////////////////////////////////

    struct PSVOXEL3D {
        long ID;
        float RawEnergy;            // MeV
        bool ghost;
        bool member_of_TKtrack;       // is this voxel part of a TKtrack
        std::vector<int> pdgs;       // all the PDG codes of tracks contributing to this voxel
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
        std::vector<struct TcalEvent::REARCALDEPOSIT> rearHCalModule; //! (NO ROOT I/O) individual module deposits
    };
    struct REARCALS rearCals;
    // Umut: to understand whats happening at rear hadron calorimeter
    void DumpRearHCalTruth(int maxPrint = 200, bool uniquePerModuleAndTrack = true);
    
    // muon tracks
    std::vector<TMuTrack> fMuTracks;

    // FASERCal detector response (PE values per voxel)
    std::vector<TcalEvent::FASERCALVOXELRESPONSE> faserCalVoxelResponse;

    /// @brief FASERCal fiber channel PE distributions (what each SiPM sees)
    /// Each fiber collects light from all voxels along its path
    struct FIBERCHANNEL {
        int channel_id;      // Unique channel identifier
        int coord1, coord2;  // Position indices (Y,Z for X-fiber; X,Z for Y-fiber; X,Y for Z-fiber)
        int layer;           // Layer index (relevant for Z-fibers)
        double totalPE;      // Total PE collected by this fiber (sum from all voxels)
        int nVoxelsHit;      // Number of voxels along this fiber with energy deposits
    };
    std::vector<FIBERCHANNEL> faserCalFiberChannelsX;  // X-fibers: indexed by (Y, Z) - SAVED TO ROOT
    std::vector<FIBERCHANNEL> faserCalFiberChannelsY;  // Y-fibers: indexed by (X, Z) - SAVED TO ROOT  
    std::vector<FIBERCHANNEL> faserCalFiberChannelsZ;  // Z-fibers: indexed by (X, Y, layer) - SAVED TO ROOT
    
    // @brief A copy of the geometry originally stored in TCalEvent
    struct TcalEvent::GEOM_DETECTOR geom_detector;

    TPORecoEvent();
    TPORecoEvent(TcalEvent* c, TPOEvent* p);
    virtual ~TPORecoEvent();

    friend std::ostream& operator<<(std::ostream& os, const TPORecoEvent& evt) {
    os << "TPORecoEvent: nMuTracks = " << evt.fMuTracks.size();
    // loop over muon tracks and dump them
    for (const auto& muTrack : evt.fMuTracks) {
        os << "\n" << muTrack.ftrackID << ": PDG=" << muTrack.fPDG << ", Charge=" << muTrack.fcharge
           << ", nPoints=" << muTrack.fpos.size() << ", Momentum=" << muTrack.fp << " GeV" << ", Chi2=" << muTrack.fchi2;
    }
    return os;
}

    int verbose = 0;                            //! controls amount of debug information
    bool multiThread = true;                   //! controls if multi-threading is used

    struct RECOCONFIG {
        // resolution fudge factor used in TTKTrack Genfit fitting for PS voxel hits
        double psvoxel_fudge_factor;

        // energy compensation factors for the PS voxel hits
        double alpha;
        double beta;
        double gamma; // for ECAL energy in total energy flow
        double delta; // for HCAL energy in total energy flow
        // e/pi compensation at voxel level
        double PS_EoverEhad_electron;
        double PS_EoverEhad_muon;
        double PS_EoverEhad_photon;
        double PS_EoverEhad_hadron;

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

        //////////////////////////////////////////////////////////
        // FASERCAL DETECTOR RESPONSE CONFIGURATION
        // Parameters for converting energy deposits to photoelectrons
        // 
        // FIBER GEOMETRY (11 x 48 x 48 voxels per layer):
        // - X-fibers: 48×48 fibers running along X (11 voxels/fiber), readout at +X
        // - Y-fibers: 11×48 fibers running along Y (48 voxels/fiber), readout at +Y  
        // - Z-fibers: 11×48 fibers running along Z (48 voxels/fiber per layer), readout at +Z
        // 
        // Each voxel contributes light to 3 fibers (one in each direction).
        // Light propagates along fiber to SiPM at positive (+) end with attenuation.
        //////////////////////////////////////////////////////////
        
        // Scintillation and light collection
        double faserCal_scintPhotonYieldPerMeV = 8000.0;  // photons/MeV
        double faserCal_globalFiberCapture = 0.10;        // 10% combined light collection efficiency
        double faserCal_fiberTrappingEfficiency = 0.05;   // 5% trapping in fiber
        double faserCal_sensorPDE = 0.25;                 // 25% SiPM photon detection efficiency
        
        // Dual-component attenuation model for each fiber direction (mm)
        // Model: Attenuation = f_short * exp(-d/L_short) + (1-f_short) * exp(-d/L_long)
        // Short component: fast attenuation (bulk absorption, defects)
        // Long component: slow attenuation (Rayleigh scattering, high-quality transmission)
        
        // Short attenuation lengths (from SuperFGD) https://arxiv.org/html/2603.14921v1
        double faserCal_fiberAttenuationLengthShortX = 350.0; // in mm
        double faserCal_fiberAttenuationLengthShortY = 350.0;
        double faserCal_fiberAttenuationLengthShortZ = 430.0;
        
        // Long attenuation lengths (from SuperFGD) https://arxiv.org/html/2603.14921v1
        double faserCal_fiberAttenuationLengthLongX = 4000.0; // in mm
        double faserCal_fiberAttenuationLengthLongY = 4000.0;
        double faserCal_fiberAttenuationLengthLongZ = 5000.0;
        
        // Fraction of light in short component (typically 0.2-0.4) https://arxiv.org/html/2603.14921v1
        double faserCal_fiberShortComponentFractionX = 0.29;
        double faserCal_fiberShortComponentFractionY = 0.29;
        double faserCal_fiberShortComponentFractionZ = 0.33;
        
        // Optical crosstalk between neighboring voxels, measured from data as
        // XT_face = sum(light leakage hits in that neighbor) / sum(light yield
        // in the central/seed cube), separately for the 4 transverse-plane
        // neighbors (x+1, x-1, y+1, y-1). VALUES ARE IN PERCENT:
        //   XT y+1 = 0.2660%, XT y-1 = 0.2784%, XT x+1 = 0.2818%, XT x-1 = 0.2398%
        // Despite Mylar foils, data shows nonzero crosstalk in X too (contrary
        // to the original "no leakage in X" prototype assumption), at a level
        // similar to Y. NOTE: the data ratio is likely an UNDERESTIMATE - the
        // 1 p.e. detection threshold suppresses low-PE leakage counts, so the
        // true crosstalk may be somewhat higher.
        //
        // These are RATIOS to the central signal (PE_neighbor / PE_central), not
        // fractions of a shared total (the 4 ratios sum to ~1.07%). They are
        // converted to model fractions (frac_face, used directly in addNeighbor
        // below) via: frac_face = ratio_face * (1 - S), S = R/(1+R), where R is
        // the sum of the 4 ratios (as fractions, i.e. divided by 100). Solving
        // gives centralFrac ~= 0.989 (~98.9% of light stays in the seed cube,
        // ~1.1% leaks out total), and per-face fractions ~0.24%-0.28%.
        // faserCal_opticalLeakageSideX/Y below are set to the per-axis AVERAGE
        // of the +/- fractions (the code applies one value to both +/-
        // neighbors on a given axis).
        //
        // faserCal_opticalLeakageZ is for the ±iz neighbors (in-layer,
        // along-fiber depth) - a different axis, NOT covered by the
        // transverse-plane crosstalk study above (no direct data). Rather than
        // the old ad hoc 1% guess, it is set to the same order of magnitude as
        // the measured X/Y values (average of SideX, SideY) as a placeholder
        // until Z crosstalk is measured directly.
        double faserCal_opticalLeakageSideX = 0.00258;  // leakage to ±X neighbors (avg of x+1=0.279%, x-1=0.237%)
        double faserCal_opticalLeakageSideY = 0.00269;  // leakage to ±Y neighbors (avg of y+1=0.263%, y-1=0.276%)
        double faserCal_opticalLeakageZ = 0.00264;       // leakage to ±Z (in-layer fiber-depth) neighbors - placeholder = avg(SideX, SideY), not yet measured directly
        
        // Statistical fluctuations
        bool faserCal_applyPoissonStatistics = true;  // Apply Poisson fluctuations to PE

        // Electronic noise / resolution at the SiPM channel level.
        // From data/MC comparison of the 2D view PE spectra: the width of the
        // Gaussian peak in data exceeds MC by these amounts (in PE), attributed
        // to SiPM+electronics noise not modeled by Poisson photon statistics alone.
        // Applied once per aggregated fiber channel (not per voxel), added in
        // quadrature-equivalent as an independent zero-mean Gaussian smearing.
        bool faserCal_applyElectronicNoise = true;
        double faserCal_electronicNoiseSigmaX = 5.11;  // PE, from XZ view data/MC comparison
        double faserCal_electronicNoiseSigmaY = 6.24;  // PE, from YZ view data/MC comparison
        double faserCal_electronicNoiseSigmaZ = 0.0;   // PE, not yet measured

        // Enable/disable detector response (set false to use energy deposits directly)
        bool faserCal_applyDetectorResponse = true;

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
        std::vector<int>& nvox_per_layer, int nvox_per_layer_max, int nzlayer);
    void Reconstruct3DPS_Eflow();

    /// @brief Recontruct particle tracks from 3D PS voxels
    void PSVoxelParticleFilter();

    /// @brief Reconstruct FASERCAL and rear calorimeters and rear mu tag
    void ReconstructRearCals();
    
    /// @brief Apply FASERCal detector response at reconstruction level
    /// @details Converts energy deposits to photoelectrons using detector response model.
    /// This replaces any PE data from Geant4, allowing parameter tuning without re-simulation.
    /// Results are stored in faserCalVoxelResponse member variable
    void ApplyFaserCalDetectorResponse();

    /// @brief Compute fiber channel PE distributions from voxel responses
    /// @details Aggregates PE from all voxels along each fiber to show what each SiPM sees.
    ///          Must be called after ApplyFaserCalDetectorResponse().
    ///          Results stored in faserCalFiberChannelsX/Y/Z member variables.
    void ComputeFaserCalFiberChannels();

    /// @brief Dump FASERCal fiber channel PE distributions
    /// @param maxChannels Maximum number of channels to print per direction (0=all)
    /// @param sortByPE If true, sort channels by PE (highest first)
    void DumpFaserCalFiberChannels(int maxChannels = 50, bool sortByPE = true);

    /// @brief Create histograms of fiber channel PE distributions
    /// @param prefix Histogram name prefix (default: "h_fiber")
    /// @return Map of histogram pointers: "X_PE", "Y_PE", "Z_PE", "X_2D", "Y_2D", "Z_2D"
    std::map<std::string, TH1*> PlotFaserCalFiberChannels(const std::string& prefix = "h_fiber");

    /// @brief Create 2D map of fiber PE for a given direction
    /// @param direction 0=X-fibers, 1=Y-fibers, 2=Z-fibers
    /// @param name Histogram name
    /// @param title Histogram title
    /// @return TH2D showing PE distribution in detector coordinates
    TH2D* CreateFiberPEMap(int direction, const std::string& name, const std::string& title);

    /// @brief Reconstruct the muon spectrometer's tracks
    void ReconstructMuonSpectrometer(); // added by Umut
    void ReconstructMuonSpectrometer_obs(); 

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

    void ReconstructMDT();
    std::vector<double> GetMDTMagnetCentersZ() const;
    std::vector<ROOT::Math::XYZVector> GetMDTMagnetCentersGlobal() const;



    ClassDef(TPORecoEvent,4)
};

#endif
