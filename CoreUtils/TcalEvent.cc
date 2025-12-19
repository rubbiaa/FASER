#include <sstream>
#include "TcalEvent.hh"

#include <TChain.h>
#include <TGeoManager.h>
#include <TGeoBBox.h>

ClassImp(TcalEvent)

TcalEvent::TcalEvent() : TObject(), fTracks(), fMagnetTracks(), fMuTagTracks() , fTPOEvent(nullptr), m_rootFile(nullptr)
{
    geom_detector = {0};
    rearCalDeposit = {};
    rearHCalDeposit = {};
    rearMuCalDeposit = {};
    // Umut: to understand whats happening at rear hadron calorimeter
    //rearHCalTruth = {};

    if(gGeoManager == nullptr) {
//        std::cerr << "Warning: gGeoManager is null! Cannot initialize TcalEvent TGeom nodes." << std::endl;
        return;
    }

    // Find the TGeoNode corresponding to the rear HCal
    TGeoNode *node = gGeoManager->GetTopNode();
    // loop over all daughters
    int ndaughters = node->GetNdaughters();
    for(int i=0; i<ndaughters; i++) {
        TGeoNode *dnode = node->GetDaughter(i);
        std::string name = dnode->GetName();
        std::cout << "Daughter " << i << " name: " << name << std::endl;
//        if(name.find("rearHCal")!=std::string::npos) {
        if(name.find("DetectorAssembly")!=std::string::npos) {
            fdetectorAssemblyTGeomNode = dnode;
            TGeoMatrix *matrix = dnode->GetMatrix();
            //matrix->Print();
            break;
        }
    }
    // check if found detector assembly node
    if(fdetectorAssemblyTGeomNode == nullptr) {
        std::cerr << "FATAL error: Could not find DetectorAssembly TGeoNode in Geometry!" << std::endl;
        throw std::runtime_error("Could not find DetectorAssembly TGeoNode");
    }

    // now iterate over all nodes to find rearHCal node
    ndaughters = fdetectorAssemblyTGeomNode->GetNdaughters();
    for(int i=0; i<ndaughters; i++) {
        TGeoNode *dnode = fdetectorAssemblyTGeomNode->GetDaughter(i);
        std::string name = dnode->GetName();
        std::cout << "  Daughter " << i << " name: " << name << std::endl;
        if(name.find("rearHCal")!=std::string::npos) {
            frearHCalTGeomNode = dnode;
            TGeoMatrix *matrix = dnode->GetMatrix();
            //matrix->Print();
            break;
        }
    }
    // check if found rearHCal node
    if(frearHCalTGeomNode == nullptr) {
        std::cerr << "FATAL error: Could not find rearHCal TGeoNode in Geometry!" << std::endl;
        throw std::runtime_error("Could not find rearHCal TGeoNode");
    }

    // Find the TGeoNode corresponding to the rear Muon Spectrometer
    TGeoNode *nodem = gGeoManager->GetTopNode();
    // loop over all daughters
    ndaughters = node->GetNdaughters();
    for(int i=0; i<ndaughters; i++) {
        TGeoNode *dnode = node->GetDaughter(i);
        std::string name = dnode->GetName();
        if(name.find("MuonSpectrometer")!=std::string::npos) {
            frearMuSpectTGeomNode = dnode;
            TGeoMatrix *matrix = dnode->GetMatrix();
            // matrix->Print();
            break;
        }
    }
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
    m_calEventTree->Branch("mutagtracks", &fMuTagTracks);
    m_calEventTree->Branch("event", &fTPOEvent);    // this should be labelled POEvent !
    m_calEventTree->Branch("geom", &geom_detector);
    m_calEventTree->Branch("rearcal", &rearCalDeposit);
    m_calEventTree->Branch("rearhcal", &rearHCalDeposit);
    m_calEventTree->Branch("rearmucal", &rearMuCalDeposit);
    // Umut: to understand whats happening at rear hadron calorimeter
    //m_calEventTree->Branch("rearhcaltruth", &rearHCalTruth);

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

    // Delete all magnet tracks
    for (auto track : fMuTagTracks)
    {
        delete track;
    }
    fMuTagTracks.clear();
    //    delete fTracks;
}

/// @brief Load event from G4 output for INPUT
/// @param base_path
/// @param ievent
/// @return error = 0, ok, error = 1, file not found, error=2 no entries found
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
    if (!m_rootFile || m_rootFile-> IsZombie() || !m_rootFile->IsOpen())
    {
      std::cerr << "Error : failed to open file " << filename.str() << std::endl;
        return 1;
    }
    m_rootFile->cd();
    TTree *event_tree = nullptr;
    m_rootFile->GetObject("calEvent",event_tree);
    if(!event_tree) {
      std::cerr << "Error : TTree 'calEvent' not found in " << filename.str() << std::endl;
      return 1;
    }

    Long_t nentries = event_tree->GetEntries();
    if(verbose > 0) std::cout << "Number of entries " << nentries << std::endl;
    if(nentries < 1) {
        std::cerr << "Number of entries " << nentries << std::endl;
        m_rootFile -> Close();
        delete m_rootFile;
        return 2;
    }

    // Set the branch address
    std::vector<DigitizedTrack*> *t = &fTracks;
    event_tree->SetBranchAddress("tracks", &t);

    std::vector<MagnetTrack*> *mt = &fMagnetTracks;
    event_tree->SetBranchAddress("magnetracks", &mt);

    std::vector<MuTagTrack*> *mut = &fMuTagTracks;
    event_tree->SetBranchAddress("mutagtracks", &mut);

//    const TPOEvent *POevent = new TPOEvent();
    event_tree -> SetBranchAddress("event", &POevent);
//    fTPOEvent = const_cast<TPOEvent*>(POevent);
    fTPOEvent = POevent;

    struct TcalEvent::GEOM_DETECTOR *g_d = &geom_detector;
    event_tree -> SetBranchAddress("geom", &g_d);

    std::vector<struct REARCALDEPOSIT> *g_r = &rearCalDeposit;
    event_tree ->SetBranchAddress("rearcal", &g_r);

    std::vector<struct REARCALDEPOSIT> *g_h = &rearHCalDeposit;
    event_tree -> SetBranchAddress("rearhcal", &g_h);

    event_tree -> SetBranchAddress("rearmucal", &rearMuCalDeposit);

    // Umut: to understand whats happening at rear hadron calorimeter
    //std::vector<struct REARHCALHITTRUTH>* g_htruth = &rearHCalTruth;
    //event_tree->SetBranchAddress("rearhcaltruth", &g_htruth);

    // Read the first entry
    event_tree->GetEntry(0);

    if(verbose > 0) {
        std::cout << "Loaded event data for event 0" << std::endl;
        std::cout << " digitized tracks " << t->size() << std::endl;
    }

    delete event_tree;
    m_rootFile -> Close();
    delete m_rootFile;

    // check geometry
        const char *volumeName = "ContainerLogical";
    TGeoVolume *targetVolume = gGeoManager->GetVolume(volumeName);

    if (targetVolume)
    {
        //std::cout << "Found volume: " << targetVolume->GetName() << std::endl;
        // You can now draw it, modify it, or use it for positioning
        // targetVolume->Draw();
        // check its dimensions
        Double_t dims[3] = {0.0, 0.0, 0.0};
        TGeoShape *shape = targetVolume->GetShape();
        if (shape) {
            // If it's a box, use the typed accessors (ROOT shapes are in cm)
            TGeoBBox *bbox = dynamic_cast<TGeoBBox*>(shape);
            if (bbox) {
                // GetDx/GetDy/GetDz return half-lengths in cm -> convert to mm and full length
                dims[0] = 2.0 * bbox->GetDX() * 10.0; // X in mm
                dims[1] = 2.0 * bbox->GetDY() * 10.0; // Y in mm
                dims[2] = 2.0 * bbox->GetDZ() * 10.0; // Z in mm
            } 
        }
        // compare to what is stored in geom_detector
        if (fabs(dims[0] - geom_detector.fScintillatorSizeX) > 1e-3 ||
            fabs(dims[1] - geom_detector.fScintillatorSizeY) > 1e-3 ||
            fabs(dims[2] - geom_detector.fTotalLength) > 1e-3)
        {
            std::cerr << "Dimensions of " << volumeName << ": "
                      << "X: " << dims[0] << " mm, "
                      << "Y: " << dims[1] << " mm, "
                      << "Z: " << dims[2] << " mm" << std::endl;
            std::cerr << "FATAL error: Geometry mismatch for " << volumeName << std::endl;
            std::cerr << " TGeoManager has X: " << dims[0] << " Y: " << dims[1] << " Z: " << dims[2] << std::endl;
            std::cerr << " geom_detector has X: " << geom_detector.fScintillatorSizeX
                      << " Y: " << geom_detector.fScintillatorSizeY
                      << " Z: " << geom_detector.fTotalLength << std::endl;
            exit(1);
        }
    }
    else
    {
        std::cerr << "Volume '" << volumeName << "' not found." << std::endl;
        exit(1);
    }

    //std::cout << "[TcalEvent] TcalEvent initialized and event loaded." << std::endl;
    //std::cout << "Check LoS shifts: FASERCal ("
    //          << geom_detector.fFASERCal_LOS_shiftX << ", "
    //          << geom_detector.fFASERCal_LOS_shiftY << "), RearHCal ("
    //          << geom_detector.frearHCal_LOS_shiftX << ", "
    //          << geom_detector.frearHCal_LOS_shiftY << "), RearMuSpect ("
    //          << geom_detector.fRearMuSpect_LOS_shiftX << ", "
    //          << geom_detector.fRearMuSpect_LOS_shiftY << "), Tilt Angle Y ("
    //          << geom_detector.fTiltAngleY << " rad)"
    //          << std::endl;

    //DumpFaserCalNodeCandidates();   

    // return success
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

// Added by UMUT
// to rotate detector coordinates to world coordinates
ROOT::Math::XYZVector TcalEvent::DetToWorld(const ROOT::Math::XYZVector& detPos) const
{
    // Detector frame -> world frame
    // Same rotation as in Geant4:
    // new G4RotationMatrix()->rotateY(fTiltAngleY);

    //constexpr double kTiltAngleDeg = 5.0; 
    //const double theta = kTiltAngleDeg * TMath::DegToRad();
    double theta = -1*geom_detector.fTiltAngleY; // given in radians
//
    const double c = std::cos(theta);
    const double s = std::sin(theta);

    const double x = detPos.X();
    const double y = detPos.Y();
    const double z = detPos.Z();

    // Rotation around Y (right–handed Geant4 / ROOT convention)
    const double x_rot =  c * x + s * z;
    const double z_rot = -s * x + c * z;

    return ROOT::Math::XYZVector(x_rot, y, z_rot);
}

// UMUT: use TGEO instead doing by hand: localPos in mm, theta in radians, optional translation (mm)
ROOT::Math::XYZVector ApplyYRotationWithTGeo(const ROOT::Math::XYZVector &localPos_mm, double theta_rad, 
                                             double tx_mm = 0.0, double ty_mm = 0.0, double tz_mm = 0.0)
{
    double theta_deg = -1 * theta_rad * 180.0 / M_PI; // to match the tilt in geometry i need to invert the angle

    TGeoRotation rot;
    rot.RotateY(theta_deg);
    TGeoCombiTrans combi(tx_mm/10.0, ty_mm/10.0, tz_mm/10.0, &rot);
    double local_cm[3] = { localPos_mm.x() / 10.0, localPos_mm.y() / 10.0, localPos_mm.z() / 10.0 };
    double global_cm[3] = {0,0,0};
    combi.LocalToMaster(local_cm, global_cm);
    return ROOT::Math::XYZVector(global_cm[0]*10.0, global_cm[1]*10.0, global_cm[2]*10.0);
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
        // i needed to comment the following lines to match the positions correctly
        x += geom_detector.fFASERCal_LOS_shiftX;
        y += geom_detector.fFASERCal_LOS_shiftY;
        double z = getZofLayer(ilayer, iz);
         ROOT::Math::XYZVector localPos(x, y, z);
        return ApplyYRotationWithTGeo(localPos, geom_detector.fTiltAngleY /*rad*/, 0.0, 0.0, 0.0);
    } else if (hittype == 1) {

        long ix = ID % 10000;
        long iy = (ID / 10000) % 10000;
        long ilayer = (ID / 100000000) % 100;
        long icopy = (ID / 10000000000LL) % 10;
        double x = ix * geom_detector.fSiTrackerPixelSize - geom_detector.fScintillatorSizeX / 2.0;
        double y = iy * geom_detector.fSiTrackerPixelSize - geom_detector.fScintillatorSizeY / 2.0;
        // i needed to comment the following lines to match the positions correctly
        x += geom_detector.fFASERCal_LOS_shiftX;
        y += geom_detector.fFASERCal_LOS_shiftY;
        double z = ilayer * geom_detector.fSandwichLength + geom_detector.fSandwichLength
            - (geom_detector.NRep * geom_detector.fSandwichLength) / 2.0;
        if(icopy == 1) z -= geom_detector.fAirGap;
        
        ROOT::Math::XYZVector localPos(x, y, z);
        return ApplyYRotationWithTGeo(localPos, geom_detector.fTiltAngleY /*rad*/, 0.0, 0.0, 0.0);
    } else {
        std::cerr << " TcalEvent::getChannelXYZfromID - hit of unknown type" << hittype << std::endl;
        return ROOT::Math::XYZVector(0,0,0);
    }
}

ROOT::Math::XYZVector TcalEvent::getChannelXYZRearCal(int moduleID) const { // the following works well
    double x = (moduleID%geom_detector.rearCalNxy - geom_detector.rearCalNxy/2.0 + 0.5)*geom_detector.rearCalSizeX;
    double y = (moduleID/geom_detector.rearCalNxy - geom_detector.rearCalNxy/2.0 + 0.5)*geom_detector.rearCalSizeY;
    x += geom_detector.fFASERCal_LOS_shiftX;
    y += geom_detector.fFASERCal_LOS_shiftY;
    double z = geom_detector.rearCalLocZ;
    return ApplyYRotationWithTGeo(ROOT::Math::XYZVector(x,y,z), geom_detector.fTiltAngleY /*rad*/, 0.0, 0.0, 0.0);
}

ROOT::Math::XYZVector TcalEvent::getChannelXYZRearHCal(int moduleID) const
{
    long ix = moduleID % 1000;
    long iy = (moduleID / 1000) % 1000;
    long iz = (moduleID / 1000000LL) % 1000;
    //    double x = geom_detector.frearHCal_LOS_shiftX + (ix - geom_detector.rearHCalNxy/2.0 + 0.5)*geom_detector.rearHCalVoxelSize;
    //    double y = geom_detector.frearHCal_LOS_shiftY + (iy - geom_detector.rearHCalNxy/2.0 + 0.5)*geom_detector.rearHCalVoxelSize;
    //    double z = geom_detector.rearHCalLocZ + iz * geom_detector.rearHCalSizeZ + geom_detector.rearHCalSizeZ/2.0;
    //    std::cout << "HCal moduleID " << moduleID << " ix " << ix << " iy " << iy << " iz " << iz << " x " << x << " y " << y << " z " << z << std::endl;
    std::cout << "[HCAL DECODE] moduleID=" << moduleID
          << " ix_dec=" << ix
          << " iy_dec=" << iy
          << " iz_dec=" << iz << std::endl;
    // Compute local coordinates 
    double x = (ix - geom_detector.rearHCalNxy / 2.0 + 0.5) * geom_detector.rearHCalVoxelSize;
    double y = (iy - geom_detector.rearHCalNxy / 2.0 + 0.5) * geom_detector.rearHCalVoxelSize;
    double z = (iz - geom_detector.rearHCalNlayer / 2.0 + 0.5) * geom_detector.rearHCalSizeZ;
    //std::cout << "HCal moduleID " << moduleID << " ix " << ix << " iy " << iy << " iz " << iz << " x " << x << " y " << y << " z " << z << " " << geom_detector.rearHCalLocZ << " " << geom_detector.rearHCalNlayer << " " << geom_detector.rearHCalVoxelSize << std::endl;
    
    ROOT::Math::XYZVector localPos(x, y, z);

    // Get the transformation matrix of the HCAL mother node
    TGeoMatrix *matrix = frearHCalTGeomNode->GetMatrix();
    //matrix->Print();

    // Transform the point using the matrix
    double local[3] = {localPos.x() / 10.0, localPos.y() / 10.0, localPos.z() / 10.0}; // Convert mm to cm for TGeo
    double global[3] = {0, 0, 0};
    matrix->LocalToMaster(local, global);
    ROOT::Math::XYZVector globalPos(global[0] * 10.0, global[1] * 10.0, global[2] * 10.0);
    std::cout << "Local (" << local[0] << ", " << local[1] << ", " << local[2] << ") -> Global ("
              << global[0] << ", " << global[1] << ", " << global[2] << ")" << std::endl;

    // now go from HCAL local to the detector assembly local
    TGeoMatrix *parentMatrix = fdetectorAssemblyTGeomNode->GetMatrix();
    //parentMatrix->Print();
    double hcal_local[3] = {globalPos.x() / 10.0, globalPos.y() / 10.0, globalPos.z() / 10.0}; // Convert mm to cm for TGeo
    double hcal_global[3] = {0, 0, 0};
    parentMatrix->LocalToMaster(hcal_local, hcal_global);

    ROOT::Math::XYZVector finalGlobalPos(hcal_global[0] * 10.0, hcal_global[1] * 10.0, hcal_global[2] * 10.0);
    std::cout << "HCal Local (" << hcal_local[0] << ", " << hcal_local[1] << ", " << hcal_local[2] << ") -> World ("
              << hcal_global[0] << ", " << hcal_global[1] << ", " << hcal_global[2] << ")" << std::endl;    

/*    // now go from detector assembly local to world
    TGeoMatrix *worldMatrix = fdetectorAssemblyTGeomNode->GetMatrix();
   // worldMatrix->Inverse();
    worldMatrix->Print(); **/
    double assembly_local[3] = {hcal_global[0], hcal_global[1], hcal_global[2]}; // Convert mm to cm for TGeo
    double assembly_global[3] = {assembly_local[0], assembly_local[1], assembly_local[2]};
    //worldMatrix->LocalToMaster(assembly_local, assembly_global);

    ROOT::Math::XYZVector worldGlobalPos(assembly_global[0] * 10.0, assembly_global[1] * 10.0, assembly_global[2] * 10.0);
    std::cout << "Detector Assembly Local (" << assembly_local[0] << ", " << assembly_local[1] << ", " << assembly_local[2] << ") -> World ("
              << assembly_global[0] << ", " << assembly_global[1] << ", " << assembly_global[2] << ")" << std::endl;

    return worldGlobalPos;
}

void TcalEvent::fillTree()
{
    m_rootFile->cd();
    m_calEventTree->Fill();
}

void TcalEvent::DumpFaserCalNodeCandidates() const {
    if (gGeoManager == nullptr) return;
    TGeoVolume* top = gGeoManager->GetTopVolume();
    TGeoIterator it(top);
    TGeoNode* n;
    std::cout << "[TcalEvent] Listing candidate nodes (name : volume) and matrix info:\n";
    while ((n = (TGeoNode*)it.Next())) {
        if (!n || !n->GetVolume()) continue;
        std::string vname = n->GetVolume()->GetName();
        if (vname.find("Container") != std::string::npos || vname.find("rear") != std::string::npos || vname.find("HCal") != std::string::npos) {
            std::cout << " node: " << n->GetName() << "  volume: " << vname << "  ndaughters=" << n->GetNdaughters() << "\n";
            if (n->GetMatrix()) {
                n->GetMatrix()->Print();
            }
        }
    }
}


ClassImp(DigitizedTrack)

DigitizedTrack *TcalEvent::addTrack(int trackID)
{
    DigitizedTrack *digitizedTrack = new DigitizedTrack();
    digitizedTrack->ftrackID = trackID;
    fTracks.push_back(digitizedTrack);
    return digitizedTrack;
}


