#ifndef _TPORECOEVENT_
#define _TPORECOEVENT_ 1

#include <TObject.h>
#include <vector>
#include <TH2D.h>

#include <TDatabasePDG.h>

#include "TcalEvent.hh"
#include "TPOEvent.hh"

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

    TPORec *fPOFullEvent = nullptr;                   // the kinematics of the full event

    TcalEvent* fTcalEvent;                            //! Reference to the TCAL event
    TPOEvent* fTPOEvent;                              //! Reference to the TPOEvent

// static    TVector3 fitLineThroughPoints(const struct TPORec::TRACK &track, TVector3& centroid);
public:

   /// @brief A hit in a two dimension plastic scintillator view
    struct PSHIT2D {
        float electromagneticity;       // =0 if hadronic, =1 if electromagnetic
        int  ntracks;                   // number of overlapping tracks
        float Edeposited;               // deposited energy
    };

    std::map<long, PSHIT2D> PShitmapX;      // the X-Z view hit map
    std::map<long, PSHIT2D> PShitmapY;      // the Y-Z view hit map

    TH2D* xviewPS = nullptr;                          //! 2Dview scintillator X-Z
    TH2D* yviewPS = nullptr;                          //! 2Dview scintillator Y-Z
    TH2D* xviewPS_em = nullptr;                          //! 2Dview scintillator X-Z
    TH2D* yviewPS_em = nullptr;                          //! 2Dview scintillator Y-Z
    TH2D* xviewPS_had = nullptr;                          //! 2Dview scintillator X-Z
    TH2D* yviewPS_had = nullptr;                          //! 2Dview scintillator Y-Z0
    TH2D* xviewPS_eldepo = nullptr;                          //! 2Dview scintillator X-Z
    TH2D* yviewPS_eldepo = nullptr;                          //! 2Dview scintillator Y-Z

    TPORecoEvent() : fTcalEvent(0), fTPOEvent(0) {};
    TPORecoEvent(TcalEvent* c, TPOEvent* p);
    virtual ~TPORecoEvent();

    /// @brief Reconstruct the FASERG4 simulated event to the PORec
    void Reconstruct();

    /// @brief Reconstruct all the tracks associated to the PORec (call this after Reconstruct)
    void TrackReconstruct();

    /// @brief Dump PORecs to the screen
    void Dump();

    /// @brief Returns the vector of Reconstructed POs
    std::vector<class TPORec*> GetPORecs() { return fPORecs;}; 

    /// @brief Returns the kinematic quantities of the full event
    /// @return TPORec of the full event
    TPORec *GetPOFullEvent() { return fPOFullEvent; };

    /// @brief Returns the truth MC information
    /// @return TPOEvent pointer of the MC truth event
    TPOEvent *GetPOEvent() { return fTPOEvent; };

    // Reconstructed event summary variables

    /// @brief Number of charged particles at the primary
    int primary_n_charged;

    /// @brief Number of scintillator hits left by tau
    int nhits_tau;

    /// @brief Number of hits left found in first tracker layers
    int nhits_tracker_first;

    /// @brief Fill the 2D (x-z) and (y-z) views of the Scintillator detector
    void Fill2DViewsPS();
    TH2D* Get2DViewXPS() { return xviewPS; };
    TH2D* Get2DViewYPS() { return yviewPS; };

    ClassDef(TPORecoEvent,1)
};

#endif