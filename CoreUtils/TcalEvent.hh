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

class MagnetTrack : public TObject {
    public:

    MagnetTrack() { };
    MagnetTrack(int trackID) { ftrackID = trackID; pos.clear(); };

    int ftrackID;           // the GEANT4 track id
    int fPDG;               // the PDG id of this track
    std::vector<ROOT::Math::XYZVector> pos;

    ClassDef(MagnetTrack,1)
};

class MuTagTrack : public TObject {
    public:

    MuTagTrack() { };
    MuTagTrack(int trackID) { ftrackID = trackID; pos.clear(); };
    ~MuTagTrack() {};

    int ftrackID;           // the GEANT4 track id
    int fPDG;               // the PDG id of this track
    std::vector<ROOT::Math::XYZVector> pos;
    std::vector<ROOT::Math::XYZVector> mom;
    
    ClassDef(MuTagTrack,1)
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
    
    /// @brief The truth tracks in the magnet to provide nice event display of MC
    std::vector<MagnetTrack*> fMagnetTracks;

    /// @Brief The truch tracks in the muTag
    std::vector<MuTagTrack*> fMuTagTracks;

    void fillTree();

    struct REARCALDEPOSIT {
        Int_t moduleID;
        Double_t energyDeposit;
    };
    std::vector<struct REARCALDEPOSIT> rearCalDeposit;    // energy deposited in rear calorimeter
    std::vector<struct REARCALDEPOSIT> rearHCalDeposit;     // energy in the rear HCal scintillator

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
        Double_t rearHCalSizeX; // in mm
        Double_t rearHCalSizeY; // in mm
        Double_t rearHCalSizeZ; // in mm
        Double_t rearHCalLocZ;  // in mm
        Int_t rearHCalNxy = 9;   // number of modules in z   
        Double_t fFASERCal_LOS_shiftX; // in mm
        Double_t fFASERCal_LOS_shiftY; // in mm 
        Double_t fAirGap; // in mm
        Double_t fAlPlateThickness; // in mm
        Double_t fSiTrackerGap; // in mm
    };

    /// @brief The summary of the detector geometry
    struct GEOM_DETECTOR geom_detector;

    /// @brief Returns type ID of a hit
    /// @param ID The hit ID (see FASERG4 DetectorConstruction class)
    /// @return =0 if scintillator hit, =1 if silicon tracker hit
    static inline long getChannelTypefromID(long ID) {
        return ID / 100000000000LL;
    }   

    /// @brief Returns the module number of a hit
    /// @param ID The hit ID (see FASERG4 DetectorConstruction class)
    /// @return The module of the hit
    static inline long getChannelModulefromID(long ID) {
        long hittype = ID / 100000000000LL;
        if (hittype == 0)
        { // hit in scintillator
            long ilayer = (ID / 1000000000);
            return ilayer;
        }
        else if (hittype == 1)
        {
            long ilayer = (ID / 100000000) % 100;
            return ilayer;
        }
        return 0;
    }

    /// @brief Returns the "layer" within module
    static inline int getChannelnzfromID(long ID) {
        return (ID / 1000000) % 1000;
    }

    /// @brief Returns the precise tracker "layer" (or copy of volume)
    /// @param ID The hit ID (see FASERG4 DetectorConstruction class)
    /// @return The precise tracker layer 
    static inline long getChannelCopyfromID(long ID) {
        long icopy = (ID / 10000000000LL) % 10;
        return icopy;
    }

    /// @brief Returns z coordinate of layer in the scintillator cubes
    /// @param layer the layer index
    /// @param iz the z position index within the layer
    /// @return z coordinate
    inline double getZofLayer(long ilayer, long iz) const {
        double z = ilayer * geom_detector.fSandwichLength + iz * geom_detector.fScintillatorVoxelSize
        - (geom_detector.NRep * geom_detector.fSandwichLength) / 2.0
        + geom_detector.fScintillatorVoxelSize/2.0;
        z += geom_detector.fAlPlateThickness + geom_detector.fTargetSizeZ;
        return z;
    }

    /// @brief Returns (x,y,z) absolute position for a given hit ID (the center of the hit/voxel)
    /// @param ID The hit ID
    /// @return The (x,y,z) absolute position of the hit
    ROOT::Math::XYZVector getChannelXYZfromID(long ID) const;

    ROOT::Math::XYZVector getChannelXYZRearCal(int moduleID) const;

    ROOT::Math::XYZVector getChannelXYZRearHCal(int moduleID) const;

    ClassDef(TcalEvent, 4)
};

#endif

