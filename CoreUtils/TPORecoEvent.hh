#ifndef _TPORECOEVENT_
#define _TPORECOEVENT_ 1

#include <TObject.h>
#include <vector>

#include <TDatabasePDG.h>

#include "TcalEvent.hh"
#include "TPOEvent.hh"

class TPORec : public TObject {
public:

   struct CALENERGIES {
        double em;             // in GeV
        double had;            // in GeV
        double Ecompensated;   // in GeV
        ROOT::Math::XYZVector cog;
        ROOT::Math::XYZVector Eflow;
    };

    int POID;                          // the primary track in POEvent
    std::vector<int> fGEANTTrackIDs;   //! all the geant track id that belong to this PORec
    std::vector<DigitizedTrack*> DTs;  //! all the DigitizedTracks that belong to this POREC
    std::vector<struct CALENERGIES> fEnergiesCogs; //! the energies and COG of each Digitized track
    
    struct CALENERGIES fTotal;         // the cumulative energies for the primary

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

public:

    TPORecoEvent() : fTcalEvent(0), fTPOEvent(0) {};
    TPORecoEvent(TcalEvent* c, TPOEvent* p);
    virtual ~TPORecoEvent();

    /// @brief Reconstruct the FASERG4 simulated event to the PORec
    void Reconstruct();

    /// @brief Dump PORecs to the screen
    void Dump();

    /// @brief Returns the vector of Reconstructed POs
    std::vector<class TPORec*> GetPORecs() { return fPORecs;}; 

    /// @brief Returns the kinematic quantities of the full event
    /// @return TPORec of the full event
    TPORec *GetPOFullEvent() { return fPOFullEvent; };

    // Reconstructed event summary variables

    /// @brief Number of charged particles at the primary
    int primary_n_charged;

    /// @brief Number of scintillator hits left by tau
    int nhits_tau;

    /// @brief Number of hits left found in first tracker layers
    int nhits_tracker_first;

    ClassDef(TPORecoEvent,1)
};

#endif