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
TcalEvent::TcalEvent(int run_number, long event_number, int event_mask) : TcalEvent()
{
    std::ostringstream fileNameStream;
    fileNameStream << "output/FASERG4-Tcalevent_" << run_number << "_" << event_number ;
    if(event_mask>0) {
        const char *mask = TPOEvent::DecodeEventMask(event_mask);
        fileNameStream << "_" << mask;
    }
    fileNameStream << ".root";
    std::string m_rootOutputFileName = fileNameStream.str();

    m_rootFile = new TFile(m_rootOutputFileName.c_str(), "RECREATE", "", 505); // last is the compression level
    if (!m_rootFile || !m_rootFile->IsOpen())
    {
        throw std::runtime_error("Could not create ROOT file");
    }
    m_rootFile->cd();

    m_calEventTree = new TTree("calEvent", "calEvent");
    m_calEventTree->Branch("tracks", &fTracks);
    m_calEventTree->Branch("magnetracks", &fMagnetTracks);
    m_calEventTree->Branch("event", &fTPOEvent);    // this should be labelled POEvent !
    m_calEventTree->Branch("geom", &geom_detector);
    m_calEventTree->Branch("rearcal", &rearCalDeposit);
    m_calEventTree->Branch("rearmucal", &rearMuCalDeposit);

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

    // Delete all magnet tracks
    for (auto track : fMagnetTracks)
    {
        delete track;
    }
    fMagnetTracks.clear();

    //    delete fTracks;
}

/// @brief Load event from G4 output for INPUT
/// @param base_path
/// @param ievent
/// @return error = 0, ok, error = 1, file not found
int TcalEvent::Load_event(std::string base_path, int run_number, int ievent,
                            int event_mask, TPOEvent *POevent) {
    std::string extension = ".root";

    // Create the filename based on ievent
    std::ostringstream filename;
    filename << base_path << "FASERG4-Tcalevent_" << run_number << "_" << ievent;
    if(event_mask>0) {
        const char *mask = TPOEvent::DecodeEventMask(event_mask);
        filename << "_" << mask;
    }
    filename << extension;

    if(verbose > 0) {
       // Print the filename to verify
       std::cout << "Loading file: " << filename.str() << " ..... ";
    }

    TFile *m_rootFile = new TFile(filename.str().c_str(), "READ"); 
    if (!m_rootFile || !m_rootFile->IsOpen())
    {
        return 1;
    }
    m_rootFile->cd();
    TTree *event_tree;
    m_rootFile->GetObject("calEvent",event_tree);

    Long_t nentries = event_tree->GetEntries();
    if(verbose > 0) std::cout << "Number of entries " << nentries << std::endl;

    // Set the branch address
    std::vector<DigitizedTrack*> *t = &fTracks;
    event_tree->SetBranchAddress("tracks", &t);

    std::vector<MagnetTrack*> *mt = &fMagnetTracks;
    event_tree->SetBranchAddress("magnetracks", &mt);

//    const TPOEvent *POevent = new TPOEvent();
    event_tree -> SetBranchAddress("event", &POevent);
//    fTPOEvent = const_cast<TPOEvent*>(POevent);
    fTPOEvent = POevent;

    struct TcalEvent::GEOM_DETECTOR *g_d = &geom_detector;
    event_tree -> SetBranchAddress("geom", &g_d);

    std::vector<struct REARCALDEPOSIT> *g_r = &rearCalDeposit;
    event_tree ->SetBranchAddress("rearcal", &g_r);

    event_tree -> SetBranchAddress("rearmucal", &rearMuCalDeposit);

    // Read the first entry
    event_tree->GetEntry(0);

    if(verbose > 0) {
        std::cout << "Loaded event data for event 0" << std::endl;
        std::cout << " digitized tracks " << t->size() << std::endl;
    }

    delete event_tree;
    m_rootFile -> Close();
    delete m_rootFile;

    return 0;
}

void TcalEvent::AssignGEANTTrackID(int G4TrackID, int PDGcode, double px, double py, double pz)
{
    //    std::cout << "Assigning GEANT4 track ID " << G4TrackID << std::endl;

    for (int i = 0; i < fTPOEvent->n_particles(); i++)
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

inline double TcalEvent::getZofLayer(long ilayer, long iz) const {
    double z = ilayer * geom_detector.fSandwichLength + iz * geom_detector.fScintillatorVoxelSize
            - (geom_detector.NRep * geom_detector.fSandwichLength) / 2.0
            + geom_detector.fScintillatorVoxelSize/2.0;
    return z;
}

ROOT::Math::XYZVector TcalEvent::getChannelXYZfromID(long ID) const
{
    long hittype = ID / 100000000000LL;
    if(hittype == 0) {        // hit in scintillator
        long ix = ID % 1000;
        long iy = (ID / 1000) % 1000;
        long iz = (ID / 1000000) % 1000;
        long ilayer = (ID / 1000000000);

        double x = ix * geom_detector.fScintillatorVoxelSize - geom_detector.fScintillatorSizeX / 2.0
            + geom_detector.fScintillatorVoxelSize/2.0;
        double y = iy * geom_detector.fScintillatorVoxelSize - geom_detector.fScintillatorSizeY / 2.0
            + geom_detector.fScintillatorVoxelSize/2.0;
//        double z = ilayer * geom_detector.fSandwichLength + iz * geom_detector.fScintillatorVoxelSize
//            - (geom_detector.NRep * geom_detector.fSandwichLength) / 2.0
//            + geom_detector.fScintillatorVoxelSize/2.0;
        double z = getZofLayer(ilayer, iz);
        return ROOT::Math::XYZVector(x, y, z);
    } else if (hittype == 1) {

        long ix = ID % 10000;
        long iy = (ID / 10000) % 10000;
        long ilayer = (ID / 100000000) % 100;
        long icopy = (ID / 10000000000LL) % 10;
        double x = ix * geom_detector.fSiTrackerPixelSize - geom_detector.fScintillatorSizeX / 2.0;
        double y = iy * geom_detector.fSiTrackerPixelSize - geom_detector.fScintillatorSizeY / 2.0;
        double z = ilayer * geom_detector.fSandwichLength + geom_detector.fSandwichLength
            - (geom_detector.NRep * geom_detector.fSandwichLength) / 2.0;
        if(icopy == 1) z -= geom_detector.fTargetSizeZ;
        return ROOT::Math::XYZVector(x, y, z);
    } else {
        std::cerr << " TcalEvent::getChannelXYZfromID - hit of unknown type" << hittype << std::endl;
        return ROOT::Math::XYZVector(0,0,0);
    }
}

ROOT::Math::XYZVector TcalEvent::getChannelXYZRearCal(int moduleID) const {
    double x = (moduleID%geom_detector.rearCalNxy - geom_detector.rearCalNxy/2.0 + 0.5)*geom_detector.rearCalSizeX;
    double y = (moduleID/geom_detector.rearCalNxy - geom_detector.rearCalNxy/2.0 + 0.5)*geom_detector.rearCalSizeY;
    double z = geom_detector.rearCalLocZ;
    return ROOT::Math::XYZVector(x, y, z);
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


