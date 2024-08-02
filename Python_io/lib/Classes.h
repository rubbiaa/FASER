// DigitizedTrack.h
#include <vector>
#include "TObject.h"

class DigitizedTrack : public TObject {
public:
    DigitizedTrack() { fhitIDs.clear(); fEnergyDeposits.clear(); }
    virtual ~DigitizedTrack() {}

    int ftrackID;            // the GEANT4 track id of this track
    int fparentID;           // the GEANT4 parent of this track
    int fprimaryID;          // the GEANT4 track id of the primary that generated this track
    int fPDG;
    std::vector<long> fhitIDs;
    std::vector<double> fEnergyDeposits;

    ClassDef(DigitizedTrack, 1)
};

class TcalEvent : public TObject {

public:
    TcalEvent() {}
    virtual ~TcalEvent() {}

    struct GEOM_DETECTOR {
        Double_t fScintillatorSizeX;  // in mm
        Double_t fScintillatorSizeY;  // in mm
        Double_t fScintillatorVoxelSize; // in mm
        Double_t fSiTrackerSizeZ; // in mm
        Double_t fSiTrackerPixelSize; // in mm
        Double_t fSandwichLength; // in mm
        Double_t fTotalLength; // in mm
        Int_t NRep;
        Double_t fTotalMass;  // in kg
        Double_t fTotalWmass; // in kg
        Double_t fTotalScintmass; // in kg

        ClassDef(GEOM_DETECTOR, 1)
    };

    GEOM_DETECTOR geom_detector;

    ClassDef(TcalEvent, 1)
};
