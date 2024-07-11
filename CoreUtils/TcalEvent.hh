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
    DigitizedTrack() { fhitIDs.clear(); fEnergyDeposits.clear(); };
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

private:
    TFile* m_rootFile;
    TTree* m_calEventTree;
    /// @brief The digitized tracks associated to this event
    std::vector<DigitizedTrack*> fTracks;
    
    int verbose = 1;

public:
    TcalEvent();
    TcalEvent(int run_number, long event_number, int event_mask);
    virtual ~TcalEvent();

    /// @brief Sets verbosity level
    /// @param verb = 0 quiet, =1 dump info to cout
    void SetVerbose(int verb) { verbose = verb; };

    /// @brief Load FASERG4 event
    /// @param base_path 
    /// @param run_number 
    /// @param ievent 
    /// @param event_mask 
    /// @param POevent 
    /// @return 
    int Load_event(std::string base_path, int run_number, int ievent, int event_mask, TPOEvent *POevent);

    TPOEvent* fTPOEvent;
    void AssignGEANTTrackID(int G4TrackID, int PDGcode, double px, double py, double pz);

    std::vector<DigitizedTrack*> getfTracks() { return fTracks;};

    DigitizedTrack* addTrack(int trackID);
    void fillTree();

    /// @brief Structure contains all basic geometry parameters (see FASERG4 DetectorConstruction class)
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
    };

    /// @brief The summary of the detector geometry
    struct GEOM_DETECTOR geom_detector;

    /// @brief Returns type ID of a hit
    /// @param ID The hit ID (see FASERG4 DetectorConstruction class)
    /// @return =0 if scintillator hit, =1 if silicon tracker hit
    long getChannelTypefromID(long ID) const;

    /// @brief Returns (x,y,z) absolute position for a given hit ID
    /// @param ID The hit ID
    /// @return The (x,y,z) absolute position of the hit
    ROOT::Math::XYZVector getChannelXYZfromID(long ID) const;

    ClassDef(TcalEvent, 1)
};

#endif

