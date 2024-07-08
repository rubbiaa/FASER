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
    };

    int POID;                          // the primary track in POEvent
    std::vector<int> fGEANTTrackIDs;   // all the geant track id that belong to this PORec
    std::vector<struct CALENERGIES> fEnergiesCogs;
    
    struct CALENERGIES fTotal;         // the cumulative energies for the primary

    TPORec() = default;
    virtual ~TPORec() = default;

    TPORec(int id) : POID(id) {};

    ClassDef(TPORec,1)
};

class TPORecoEvent : public TObject {

public:

    std::vector<class TPORec*> fPORecs;

    TcalEvent* fTcalEvent;
    TPOEvent* fTPOEvent;

    TPORecoEvent(TcalEvent* c, TPOEvent* p);
    virtual ~TPORecoEvent();

    void Reconstruct();
    void Dump();

    struct TPORec::CALENERGIES computeEnergiesAndCOG(DigitizedTrack *dt);

    ClassDef(TPORecoEvent,1)
};

#endif