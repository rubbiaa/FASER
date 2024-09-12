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
    TPOEvent* fTPOEvent;

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

    void AssignGEANTTrackID(int G4TrackID, int PDGcode, double px, double py, double pz);

    std::vector<DigitizedTrack*> getfTracks() { return fTracks;};

    DigitizedTrack* addTrack(int trackID);
    
    void fillTree();

    struct REARCALDEPOSIT {
        Int_t moduleID;
        Double_t energyDeposit;
    };
    std::vector<struct REARCALDEPOSIT> rearCalDeposit;    // energy deposited in rear calorimeter

    double rearMuCalDeposit;     // energy in the rear MuCal scintillator
        
    /// @brief Structure contains all basic geometry parameters (see FASERG4 DetectorConstruction class)
    struct GEOM_DETECTOR {
        Double_t fScintillatorSizeX;  // in mm
        Double_t fScintillatorSizeY;  // in mm
        Double_t fScintillatorVoxelSize; // in mm
        Double_t fSiTrackerSizeZ; // in mm
        Double_t fSiTrackerPixelSize; // in mm
        Double_t fTargetSizeZ; // in mm
        Int_t NRep_SiTracker;     // number of Si layers
        Double_t fSandwichLength; // in mm
        Double_t fTotalLength; // in mm
        Int_t NRep;
        Double_t fTotalMass;  // in kg
        Double_t fTotalWmass; // in kg
        Double_t fTotalScintmass; // in kg
        Double_t rearCalSizeX; // in mm
        Double_t rearCalSizeY; // in mm
        Double_t rearCalLocZ;  // in mm
        Int_t rearCalNxy = 5;   // number of modules in x, and y
    };

    /// @brief The summary of the detector geometry
    struct GEOM_DETECTOR geom_detector;

    /// @brief Returns type ID of a hit
    /// @param ID The hit ID (see FASERG4 DetectorConstruction class)
    /// @return =0 if scintillator hit, =1 if silicon tracker hit
    long getChannelTypefromID(long ID) const;

    /// @brief Returns the layer number of a hit
    /// @param ID The hit ID (see FASERG4 DetectorConstruction class)
    /// @return The hit layer
    long getChannelLayerfromID(long ID) const;

    /// @brief Returns the precise tracker "layer" (or copy of volume)
    /// @param ID The hit ID (see FASERG4 DetectorConstruction class)
    /// @return The precise tracker layer 
    long getChannelCopyfromID(long ID) const;

    /// @brief Returns z coordinate of layer
    /// @param layer the layer index
    /// @param iz the z position index within the layer
    /// @return z coordinate
    inline double getZofLayer(long ilayer, long iz) const;

    /// @brief Returns (x,y,z) absolute position for a given hit ID (the center of the hit/voxel)
    /// @param ID The hit ID
    /// @return The (x,y,z) absolute position of the hit
    ROOT::Math::XYZVector getChannelXYZfromID(long ID) const;

    ROOT::Math::XYZVector getChannelXYZRearCal(int moduleID) const;

    ClassDef(TcalEvent, 3)
};

#endif

