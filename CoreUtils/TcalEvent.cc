#include "TcalEvent.hh"

#include <TChain.h>

ClassImp(TcalEvent)

TcalEvent::TcalEvent()
{
    fTPOEvent = nullptr;
    m_rootFile = nullptr;
}

/// @brief Create a TcalEvent with a given event number for OUTPUT
/// @param event_number 
TcalEvent::TcalEvent(long event_number) : TcalEvent()
{
    std::ostringstream fileNameStream;
    fileNameStream << "output/tcalevent_" << event_number << ".root";
    std::string m_rootOutputFileName = fileNameStream.str();

    m_rootFile = new TFile(m_rootOutputFileName.c_str(), "RECREATE", "", 505); // last is the compression level
    if (!m_rootFile || !m_rootFile->IsOpen())
    {
        throw std::runtime_error("Could not create ROOT file");
    }
    m_rootFile->cd();

    m_calEventTree = new TTree("calEvent", "calEvent");
    m_calEventTree->Branch("tracks", &fTracks);
    m_calEventTree->Branch("event", &fTPOEvent);
    m_calEventTree->Branch("geom", &geom_detector);

    //    fTracks = new std::vector<DigitizedTrack*>;
}

TcalEvent::~TcalEvent()
{
    if(m_rootFile != nullptr){
        m_rootFile->cd();
        m_calEventTree->Write();
        m_rootFile->Close();
    }

    // Delete all tracks
    for (auto track : fTracks)
    {
        delete track;
    }
    fTracks.clear();
    //    delete fTracks;
}

/// @brief Load event from G4 output for INPUT
/// @param base_path 
/// @param ievent 
/// @return error = 0, ok, error = 1, file not found
int TcalEvent::Load_event(std::string base_path, int ievent, TPOEvent *POevent) {
    std::string extension = ".root";

    // Create the filename based on ievent
    std::ostringstream filename;
    filename << base_path << ievent << extension;

    // Print the filename to verify
    std::cout << "Loading file: " << filename.str() << " ..... ";

    TChain *event_tree = new TChain("calEvent");
    int nFiles = event_tree->Add(filename.str().c_str());

    Long_t nentries = event_tree->GetEntries();
    std::cout << "Number of entries " << nentries << std::endl;

    if(nentries == 0){
        return 1;
    }

    // Set the branch address
    std::vector<DigitizedTrack*> *t = &fTracks;
    event_tree->SetBranchAddress("tracks", &t);

//    const TPOEvent *POevent = new TPOEvent();
    event_tree -> SetBranchAddress("event", &POevent);
//    fTPOEvent = const_cast<TPOEvent*>(POevent);
    fTPOEvent = POevent;

    struct TcalEvent::GEOM_DETECTOR *g_d = &geom_detector;
    event_tree -> SetBranchAddress("geom", &g_d);
 
    // Read the first entry
    event_tree->GetEntry(0);

    // Use the loaded data (example)
    std::cout << "Loaded event data for event 0" << std::endl;
    std::cout << " digitized tracks " << t->size() << std::endl;

    delete event_tree;

    return 0;
}

void TcalEvent::AssignGEANTTrackID(int G4TrackID, int PDGcode, double px, double py, double pz)
{
    //    std::cout << "Assigning GEANT4 track ID " << G4TrackID << std::endl;

    for (int i = 0; i < fTPOEvent->n_particles; i++)
    {
        const struct PO *aPO = &fTPOEvent->POs[i];
        if (aPO->geanttrackID > -1)
            continue;
        if (aPO->m_pdg_id != PDGcode)
            continue;
        double dpx = fabs(px - aPO->m_px);
        double dpy = fabs(py - aPO->m_py);
        double dpz = fabs(pz - aPO->m_pz);
        if (dpx < 1e-6 && dpy < 1e-6 && dpz < 1e-6)
        {
            //           std::cout << ".... Assign to PO number " << i << std::endl;
            fTPOEvent->setGEANT4TrackID(i, G4TrackID);
            break;
        }
    }
}

/// @brief Return the type of the hit 
/// @param ID The hit ID
/// @return =0 for scintillator, = 1 for silicon tracker hit
long TcalEvent::getChannelTypefromID(long ID) const {
    return ID / 100000000000LL;
}

ROOT::Math::XYZVector TcalEvent::getChannelXYZfromID(long ID) const
{
    long hittype = ID / 100000000000LL;
    if(hittype == 0) {        // hit in scintillator
        long ix = ID % 1000;
        long iy = (ID / 1000) % 1000;
        long iz = (ID / 1000000) % 1000;
        long ilayer = (ID / 1000000000);

        double x = ix * geom_detector.fScintillatorVoxelSize - geom_detector.fScintillatorSizeX / 2.0;
        double y = iy * geom_detector.fScintillatorVoxelSize - geom_detector.fScintillatorSizeY / 2.0;
        double z = ilayer * geom_detector.fSandwichLength + iz * geom_detector.fScintillatorVoxelSize 
            - (geom_detector.NRep * geom_detector.fSandwichLength) / 2.0;
        return ROOT::Math::XYZVector(x, y, z);
    } else if (hittype == 1) {

        double fSiTrackerSizeZ = 0.2;
	    double fSiTrackerPixelSize = 0.1;  // TODO: read from geometry

        long ix = ID % 10000;
        long iy = (ID / 10000) % 10000;
        long ilayer = (ID / 100000000) % 1000;
        double x = ix * fSiTrackerPixelSize - geom_detector.fScintillatorSizeX / 2.0;
        double y = iy * fSiTrackerPixelSize - geom_detector.fScintillatorSizeY / 2.0;
        double z = ilayer * geom_detector.fSandwichLength + geom_detector.fSandwichLength - fSiTrackerSizeZ
            - (geom_detector.NRep * geom_detector.fSandwichLength) / 2.0;
        return ROOT::Math::XYZVector(x, y, z);
    } else {
        std::cerr << " TcalEvent::getChannelXYZfromID - hit of unknown type" << hittype << std::endl;
        return ROOT::Math::XYZVector(0,0,0);
    }
}

void TcalEvent::fillTree()
{
    m_rootFile->cd();
    m_calEventTree->Fill();
}

ClassImp(DigitizedTrack)

DigitizedTrack *TcalEvent::addTrack(int trackID)
{
    DigitizedTrack *digitizedTrack = new DigitizedTrack();
    digitizedTrack->ftrackID = trackID;
    fTracks.push_back(digitizedTrack);
    return digitizedTrack;
}



