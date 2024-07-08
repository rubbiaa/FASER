#ifndef CALEVENT_H
#define CALEVENT_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <Math/Vector3D.h>

#include "TObject.h"
#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"

#include "TPOEvent.hh"

/**
 * @class calEvent
 * @brief Object to store the information of FASER calo event
 * @details 
 */

class DigitizedTrack : public TObject {

    public:
    DigitizedTrack() {};
//    DigitizedTrack(int trackID):ftrackID(trackID) {};

    virtual ~DigitizedTrack() {};

    int ftrackID;            // the GEANT4 track id of this track
    int fparentID;           // the GEANT4 parent of this track
    int fprimaryID;          // the GEANT4 track id of the primary that generated this track
    int fPDG;
    std::vector<long> fhitIDs;
    std::vector<double> fEnergyDeposits;

    ClassDef(DigitizedTrack, 1)
};

class TcalEvent : public TObject {
    TFile* m_rootFile;
    TTree* m_calEventTree;

public:
    TcalEvent();
    TcalEvent(long event_number);
    virtual ~TcalEvent();

    int Load_event(std::string base_path, int ievent, TPOEvent *POevent);

    TPOEvent* fTPOEvent;
    void AssignGEANTTrackID(int G4TrackID, int PDGcode, double px, double py, double pz);

    std::vector<DigitizedTrack*> fTracks;
    std::vector<DigitizedTrack*> getfTracks() { return fTracks;};

    DigitizedTrack* addTrack(int trackID);
    void fillTree();

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
    };
 
    struct GEOM_DETECTOR geom_detector;
    long getChannelTypefromID(long ID) const;
    ROOT::Math::XYZVector getChannelXYZfromID(long ID) const;

    ClassDef(TcalEvent, 1)
};

#endif

