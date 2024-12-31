#include <TCanvas.h>
#include <TGeoManager.h>
#include "TGeoVolume.h"
#include "TGeoSphere.h"
#include "TGeoBBox.h"
#include "TGeoMedium.h"
#include <TView3D.h>
#include <TChain.h>
#include <TText.h>
#include <TGTextEntry.h>
#include <TStyle.h>
#include <TControlBar.h>
#include <TButton.h>
#include <TGTab.h>
#include <TGListTree.h>
#include <TGClient.h>

#include "MyMainFrame.h"
#include "TPORecoEvent.hh"

MyMainFrame::MyMainFrame(int run_number, int ieve, int mask, bool pre, const TGWindow *p, UInt_t w, UInt_t h) : fTcalEvent(0) {

    process_reco_event = pre;

// create window
    fMain = new TGMainFrame(p, w, h);

    TGTab *tab = new TGTab(fMain, w, h);
    TGCompositeFrame *tab1 = tab->AddTab("Event");
    TGCompositeFrame *tab2 = tab->AddTab("2DPSView");
    TGCompositeFrame *tab2b = tab->AddTab("2DPSViewZ");
    TGCompositeFrame *tab3 = tab->AddTab("2DPSView_emhad");
    TGCompositeFrame *tab4 = tab->AddTab("eldepo");
    fMain->AddFrame(tab, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    // Create a horizontal frame to contain the canvas and tree
    TGHorizontalFrame *hFrameMain = new TGHorizontalFrame(tab1, 1600, 600);
    tab1->AddFrame(hFrameMain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    // Create a TGListTree on the left
    TGCanvas *listCanvas = new TGCanvas(hFrameMain, 300, 600);
    listTree = new TGListTree(listCanvas->GetViewPort(), kHorizontalFrame);
    listCanvas->SetContainer(listTree);
    hFrameMain->AddFrame(listCanvas, new TGLayoutHints(kLHintsLeft | kLHintsExpandY));

    // Create an embedded canvas
    fCanvas = new TRootEmbeddedCanvas("EmbeddedCanvas", tab1, 1300, 600);
    hFrameMain->AddFrame(fCanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    // Create a horizontal frame to contain the toggle buttons
    TGHorizontalFrame *hFrame = new TGHorizontalFrame(tab1);
    fButton = new TGTextButton(hFrame, "Toggle prim_em");
    fButton->Connect("Clicked()", "MyMainFrame", this, "toggle_prim_em()");
    hFrame->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));
    fButton = new TGTextButton(hFrame, "Toggle prim_had");
    fButton->Connect("Clicked()", "MyMainFrame", this, "toggle_prim_had()");
    hFrame->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));
    fButton = new TGTextButton(hFrame, "Toggle sec_em");
    fButton->Connect("Clicked()", "MyMainFrame", this, "toggle_sec_em()");
    hFrame->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));
    fButton = new TGTextButton(hFrame, "Toggle sec_had");
    fButton->Connect("Clicked()", "MyMainFrame", this, "toggle_sec_had()");
    hFrame->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));
    fButton = new TGTextButton(hFrame, "Toggle reco_track");
    fButton->Connect("Clicked()", "MyMainFrame", this, "toggle_reco_track()");
    hFrame->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));
    fButton = new TGTextButton(hFrame, "Toggle reco_pstrack");
    fButton->Connect("Clicked()", "MyMainFrame", this, "toggle_recon_ps_tracks()");
    hFrame->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));
    fButton = new TGTextButton(hFrame, "Toggle reco_voxel");
    fButton->Connect("Clicked()", "MyMainFrame", this, "toggle_reco_voxels()");
    hFrame->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));
    fButton = new TGTextButton(hFrame, "Toggle color fakes");
    fButton->Connect("Clicked()", "MyMainFrame", this, "toggle_color_fakes()");
    hFrame->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));
    fButton = new TGTextButton(hFrame, "ONLY RECO");
    fButton->Connect("Clicked()", "MyMainFrame", this, "only_reco()");
    hFrame->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));

    f_fullreco_CheckBox = new TGCheckButton(hFrame, "Full Reco");
    f_fullreco_CheckBox->Connect("Toggled(Bool_t)", "MyMainFrame", this, "on_fullreco_toggle(Bool_t)");
    hFrame->AddFrame(f_fullreco_CheckBox, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));

    // Create a horizontal frame to contain the zoom and sideview buttons
    TGHorizontalFrame *hFrame2 = new TGHorizontalFrame(tab1);
    // Create a button
    fButton = new TGTextButton(hFrame2, "Side View");
    fButton->Connect("Clicked()", "MyMainFrame", this, "SideView()");
    hFrame2->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));
    fButton = new TGTextButton(hFrame2, "Zoom vtx");
    fButton->Connect("Clicked()", "MyMainFrame", this, "HandleButton()");
    hFrame2->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));
    fButton = new TGTextButton(hFrame2, "Zoom in");
    fButton->Connect("Clicked()", "MyMainFrame", this, "ZoomIn()");
    hFrame2->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));
    fButton = new TGTextButton(hFrame2, "Zoom out");
    fButton->Connect("Clicked()", "MyMainFrame", this, "ZoomOut()");
    hFrame2->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));
    fButton = new TGTextButton(hFrame2, "Move up");
    fButton->Connect("Clicked()", "MyMainFrame", this, "MoveUp()");
    hFrame2->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));
    fButton = new TGTextButton(hFrame2, "Move down");
    fButton->Connect("Clicked()", "MyMainFrame", this, "MoveDown()");
    hFrame2->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));
    fButton = new TGTextButton(hFrame2, "Move left");
    fButton->Connect("Clicked()", "MyMainFrame", this, "MoveLeft()");
    hFrame2->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));
    fButton = new TGTextButton(hFrame2, "Move right");
    fButton->Connect("Clicked()", "MyMainFrame", this, "MoveRight()");
    hFrame2->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));

    // Add the horizontal frame to the main frame
    tab1->AddFrame(hFrame, new TGLayoutHints(kLHintsCenterX | kLHintsBottom, 5, 5, 3, 4));

   // Add the horizontal frame to the main frame
    tab1->AddFrame(hFrame2, new TGLayoutHints(kLHintsCenterX | kLHintsBottom, 5, 5, 3, 4));

    fCanvas_2DPSview = new TRootEmbeddedCanvas("EmbeddedCanvas2", tab2, 1200, 600);;
    tab2->AddFrame(fCanvas_2DPSview, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    fCanvas_2DPSviewZ = new TRootEmbeddedCanvas("EmbeddedCanvas2b", tab2b, 1200, 600);;
    tab2b->AddFrame(fCanvas_2DPSviewZ, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    fCanvas_2DPSview_emhad = new TRootEmbeddedCanvas("EmbeddedCanvas3", tab3, 1200, 600);;
    tab3->AddFrame(fCanvas_2DPSview_emhad, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    fCanvas_eldepo = new TRootEmbeddedCanvas("EmbeddedCanvas4", tab4, 1200, 600);;
    tab4->AddFrame(fCanvas_eldepo, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    // Create a horizontal frame to contain the next and goto event buttons
    TGHorizontalFrame *hFrame3 = new TGHorizontalFrame(fMain);
    fButton = new TGTextButton(hFrame3, "Next Event");
    fButton->Connect("Clicked()", "MyMainFrame", this, "next_event()");
    hFrame3->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsBottom, 5, 5, 3, 4));
    fButton = new TGTextButton(hFrame3, "Goto Event:");
    fButton->Connect("Clicked()", "MyMainFrame", this, "goto_event()");
    textNextEventEntry = new TGTextEntry(hFrame3, new TGTextBuffer(50));
    textNextEventEntry->Connect("ReturnPressed()", "MyMainFrame", this, "goto_event()");
    hFrame3->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5, 5, 3, 4));
    hFrame3->AddFrame(textNextEventEntry, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));
    // Add the horizontal frame to the main frame
    fMain->AddFrame(hFrame3, new TGLayoutHints(kLHintsCenterX | kLHintsBottom, 5, 5, 3, 4));

    fMain->SetWindowName("The FASERCal event display");
    fMain->MapSubwindows();
    fMain->Resize(fMain->GetDefaultSize());
    fMain->MapWindow();

    range_min[0] = 12.5; range_min[1] = 12.5; range_min[2] = -200.0;
    range_max[0] = 25.0; range_max[1] = 25.0; range_max[2] = 200.0;

    // load event
    ievent = ieve;
    opened_reco_event = false;

    if(!process_reco_event) 
        Load_event(run_number, ievent, mask);
    else
        Load_Recoevent(run_number, ievent);

    TCanvas *canvas = fCanvas->GetCanvas();
    canvas->cd();
    // Draw the geometry
    Draw_event();
    Update_ListTree();

    toggle_primary_em=
    toggle_primary_had=
    toggle_secondary_em=
    toggle_secondary_had=false;

    toggle_reconstructed_tracks = true;
    toggle_reconstructed_ps_tracks = true;

    toggle_color_fake_voxel = false;
    toggle_reco_voxel = true;
}

// Destructor
MyMainFrame::~MyMainFrame() {
    fMain->Cleanup();
    delete fMain;
}

void MyMainFrame::Load_event(int run_number, int ievent, int mask) {

    std::string base_path = "input/";

    // Create an instance of TcalEvent and TPOEvent
    fTcalEvent = new TcalEvent();
    POevent = new TPOEvent();

    event_mask = mask;
    int error = fTcalEvent -> Load_event(base_path, run_number, ievent, mask, POevent);
    if(error !=0) return;
    std::cout << "Transverse size " << fTcalEvent->geom_detector.fScintillatorSizeX << " mm " << std::endl;
    std::cout << "Total size of one sandwich layer " << fTcalEvent->geom_detector.fSandwichLength << " mm " << std::endl;
	std::cout << "Number of layers " << fTcalEvent->geom_detector.NRep << std::endl;
    std::cout << "Voxel size " << fTcalEvent->geom_detector.fScintillatorVoxelSize << " mm " << std::endl;

    std::cout << " copied digitized tracks " << fTcalEvent->getfTracks().size() << std::endl;

    fTcalEvent -> fTPOEvent -> dump_event();

    fPORecoEvent = new TPORecoEvent(fTcalEvent, fTcalEvent->fTPOEvent);
    fPORecoEvent->verbose = 3;
    std::cout << "Start reconstruction of PORecs..." << std::endl;
    fPORecoEvent -> ReconstructTruth();
    std::cout << "Start reconstruction of clusters..." << std::endl;
    fPORecoEvent -> Reconstruct2DViewsPS();
    if(f_fullreco_CheckBox->IsOn() || true) {
        fPORecoEvent->Reconstruct3DPS_2();
    }
    fPORecoEvent -> ReconstructRearCals();
    fPORecoEvent -> ReconstructClusters(0);    // this is very slow
    fPORecoEvent -> Reconstruct3DPS_Eflow();
    std::cout << "Start reconstruction of tracks..." << std::endl;
    if(f_fullreco_CheckBox->IsOn() || true) {
        fPORecoEvent->TrackReconstruct();
    }
    fPORecoEvent->PSVoxelParticleFilter();
    fPORecoEvent -> Dump();

    // fill 2D maps
    fPORecoEvent -> Fill2DViewsPS();
}

void MyMainFrame::Load_Recoevent(int run_number, int ievent) {
   if(!opened_reco_event){
       std::ostringstream inputfilename;
       inputfilename << "../Batch/Batch-TPORecevent_" << run_number << "_*_*.root";

       reco_event_tree = new TChain("RecoEvent", "READ");
       reco_event_tree->Add(inputfilename.str().c_str());

       Long_t nentries = reco_event_tree->GetEntries();
       std::cout << "Number of entries " << nentries << std::endl;

       // TPORecoEvent class and histograms
       fPORecoEvent = nullptr; // &fTPORecoEvent;
       reco_event_tree->SetBranchAddress("TPORecoEvent", &fPORecoEvent);
       opened_reco_event = true;

       // dummy tcalEvent
       fTcalEvent = new TcalEvent();
   }
   reco_event_tree->GetEntry(ievent);
   POevent = fPORecoEvent -> GetPOEvent();
   fTcalEvent -> geom_detector = fPORecoEvent -> geom_detector;
   fTcalEvent -> rearCalDeposit.clear();
   for (const auto &it : fPORecoEvent -> rearCals.rearCalModule) {
       fTcalEvent -> rearCalDeposit.push_back(it);
   }
   fTcalEvent -> fTPOEvent = POevent;
}

void MyMainFrame::Update_ListTree() {
    // loop over all children and delete them
    TGListTreeItem *child = listTree->GetFirstItem();
    while (child) {
        listTree->DeleteItem(child);
        child = listTree->GetFirstItem();
    }
    trackMap.clear();
    vertexMap.clear();

    // Populate the tree with some items
    TGListTreeItem *root = listTree->AddItem(nullptr, "FASERCal Event");
    TGListTreeItem *item1 = listTree->AddItem(root, "Tracks");
    TGListTreeItem *item2 = listTree->AddItem(root, "Vertices");
    listTree->OpenItem(root); // Expand the root node
    
    // Add the tracks from TPORecoEvent
    for (const auto &track : fPORecoEvent->fTKTracks) {
        std::ostringstream trackname;
        trackname << "Track " << track.trackID;
        // print fit status, chi2 and pval
        if(track.fitTrack != nullptr) {
            trackname << "(p=" << track.fitTrack->getFittedState().getMom().Mag();
            trackname << ",chi2=" << track.fitTrack->getFitStatus()->getChi2() << ", pval=" << track.fitTrack->getFitStatus()->getPVal() << ")";
        }
        TGListTreeItem *trackitem = listTree->AddItem(item1, trackname.str().c_str());
        trackMap[trackname.str()] = track.trackID;
        for (const auto &hit : track.tkhit) {
            std::ostringstream hitname;
            hitname << "Hit " << hit.ID;
            listTree->AddItem(trackitem, hitname.str().c_str());
        }
    }

    // Add the vertices from TPORecoEvent
    for (const auto &vertex : fPORecoEvent->fTKVertices) {
        std::ostringstream vtxname;
        vtxname << "Vertex " << vertex.vertexID;
        // add position and chi2 to the name
        vtxname << " (" << vertex.position.X() << ", " << vertex.position.Y() << ", " << vertex.position.Z() << ") ";
        vtxname << "chi2=" << vertex.chi2;
        TGListTreeItem *vtxitem = listTree->AddItem(item2, vtxname.str().c_str());
        vertexMap[vtxname.str()] = vertex.vertexID;
        for (const auto &trackid : vertex.trackIDs) {
            std::ostringstream trackname;
            trackname << "Track " << trackid;
            listTree->AddItem(vtxitem, trackname.str().c_str());
        }
    }

    listTree->Connect("Clicked(TGListTreeItem*,Int_t)", "MyMainFrame", this, "OnTrackSelected(TGListTreeItem*,Int_t)");
}

void MyMainFrame::OnTrackSelected(TGListTreeItem *entry, Int_t btn) {
    if (entry) {
        std::string name = entry->GetText();
        if (trackMap.find(name) != trackMap.end()) {
            selectedTrackID = trackMap[name];
            std::cout << "Selected track " << selectedTrackID << std::endl;
            // dump track
            for (const auto &track : fPORecoEvent->fTKTracks) {
                if (track.trackID == selectedTrackID) {
                    track.Dump();
                    break;
                }
            }
        } else 
            selectedTrackID = -1;
        if (vertexMap.find(name) != vertexMap.end()) {
            selectedVertexID = vertexMap[name];
            std::cout << "Selected vertex " << selectedVertexID << std::endl;
        } else
            selectedVertexID = -1;

        // Redraw the event
        toggle_reconstructed_tracks = false;
        toggle_reco_track();
    }
}

void MyMainFrame::Draw_event() {

    TCanvas *canvas = fCanvas->GetCanvas();
    canvas->cd();
    canvas->Clear();

    // Draw the geometry
    gGeoManager->GetTopVolume()->Draw("gl");
    SideView();
    double wdx, wdy, wdz;
    if (TGeoBBox *box = dynamic_cast<TGeoBBox*>(gGeoManager->GetTopVolume()->GetShape())) {
        wdx = box->GetDX();
        wdy = box->GetDY();
        wdz = box->GetDZ()/2.0;
        printf("Top volume dimensions: DX = %f, DY = %f, DZ = %f\n", wdx, wdy, wdz);
    } else {
        exit(1);
    }
    bigbox = new TGeoBBox("bigbox", wdx, wdy, wdz);

    TGeoMedium *air = gGeoManager->GetMedium("AIR");
    primary_em = new TGeoVolume("primary_em", bigbox, air);
    primary_had = new TGeoVolume("primary_had", bigbox, air);
    secondary_em = new TGeoVolume("secondary_em", bigbox, air);
    secondary_had = new TGeoVolume("secondary_had", bigbox, air);
    si_tracker = new TGeoVolume("si_tracker", bigbox, air);
    si_vertices = new TGeoVolume("si_vertices", bigbox, air);
    rearcal = new TGeoVolume("rearcal", bigbox, air);
    rearHcal = new TGeoVolume("rearHcal", bigbox, air);
    rearMucal = new TGeoVolume("rearMucal", bigbox, air);

    TGeoMaterial *matAluminum = new TGeoMaterial("Aluminum", 26.98, 13, 2.7);
    TGeoMedium *aluminum = new TGeoMedium("Aluminum", 2, matAluminum);
    double voxelsize = fTcalEvent->geom_detector.fScintillatorVoxelSize/10.0;
    TGeoShape *box = new TGeoBBox("box", voxelsize/2.0,voxelsize/2.0,voxelsize/2.0);

    TGeoShape *trackerhitbox = new TGeoBBox("box", 0.1/2.0,0.1/2.0,0.1/2.0);

    for (const auto &track : fTcalEvent->getfTracks())
    {
        //        std::cout << track->ftrackID << std::endl;
        size_t nhits = track->fhitIDs.size();
        //        std::cout << nhits << std::endl;
        //        if(track->fparentID > 0) continue;
        for (size_t i = 0; i < nhits; i++)
        {

            long hittype = fTcalEvent->getChannelTypefromID(track->fhitIDs[i]);

            // apply energy cut on scintillator voxel
            if (hittype == 0 && track->fEnergyDeposits[i] < 0.5)
                continue;
            // apply cut on pixel hit
            if(hittype == 1 && track->fEnergyDeposits[i] < 1e-3)continue;

            ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(track->fhitIDs[i]);
            // Create a translation matrix for the hit position

            if (hittype == 0)
            {
                //                position += ROOT::Math::XYZVector(2.5,2.5,2.5);    // in mm
                TGeoTranslation *trans = new TGeoTranslation(position.X() / 10.0,
                                                             position.Y() / 10.0, position.Z() / 10.0);
                TGeoVolume *hitVolume = new TGeoVolume("HitVolume", box, air);
                hitVolume->SetLineColor(kRed);
                if (fabs(track->fPDG) == 11)
                {
                    hitVolume->SetLineColor(kBlue); // electromagnetic is blue
                }
                else if (fabs(track->fPDG) == 13)
                {
                    hitVolume->SetLineColor(kGreen); // muons
                }
                else if (fabs(track->fPDG) == 15)
                {
                    hitVolume->SetLineColor(kCyan); // taus
                }

                // Add the hit volume to the top volume with the translation
                if (track->fparentID == 0)
                {
                    if (fabs(track->fPDG) == 11 || fabs(track->fPDG) == 13 || fabs(track->fPDG) == 15)
                    {
                        primary_em->AddNode(hitVolume, i, trans);
                    }
                    else
                    {
                        primary_had->AddNode(hitVolume, i, trans);
                    }
                }
                else
                {
                    if (fabs(track->fPDG) == 11 || fabs(track->fPDG) == 13)
                    {
                        secondary_em->AddNode(hitVolume, i, trans);
                    }
                    else
                    {
                        secondary_had->AddNode(hitVolume, i, trans);
                    }
                }
            }
            else if (hittype == 1)
            {
                //                position += ROOT::Math::XYZVector(0.05,0.05,-0.2);
                TGeoTranslation *trans = new TGeoTranslation(position.X() / 10.0,
                                                             position.Y() / 10.0, position.Z() / 10.0);
                TGeoVolume *hitVolume = new TGeoVolume("TrackerHitVolume", trackerhitbox, air);
                hitVolume->SetLineColor(kRed);
                si_tracker->AddNode(hitVolume, i, trans);
            }
            else
            {
                std::cout << " Unknown type of hit " << std::endl;
            }
        }
    }

    // draw magnet MC tracks
    polylineMagnetTracks.clear();
    for (const auto &it : fTcalEvent->fMagnetTracks)
    {
        int nhits = it->pos.size();
        if (nhits == 0)
            continue;
        double *x = (double *)malloc(nhits * sizeof(double));
        double *y = (double *)malloc(nhits * sizeof(double));
        double *z = (double *)malloc(nhits * sizeof(double));
        int idx = 0;
        for (const auto &hit : it->pos)
        {
            x[idx] = hit.x() / 10.0;
            y[idx] = hit.y() / 10.0;
            z[idx] = hit.z() / 10.0;
            idx++;
        }
        TPolyLine3D *trackpoly = new TPolyLine3D(nhits, x, y, z);
        if (std::abs(it->fPDG) == 13)
        {
            trackpoly->SetLineColor(kGreen);
        }
        else
            trackpoly->SetLineColor(kBlack);
        trackpoly->SetLineWidth(1);
        trackpoly->Draw("same");
        polylineTracks.push_back(trackpoly);
    }

    // draw rear Ecal hits
    for (const auto &it : fTcalEvent->rearCalDeposit)
    {
        ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZRearCal(it.moduleID);
        double zBox = it.energyDeposit / 1e2; // 1cm is 1 GeV
        TGeoShape *box = new TGeoBBox("rearcalbox", fTcalEvent->geom_detector.rearCalSizeX / 20.0,
                                      fTcalEvent->geom_detector.rearCalSizeY / 20.0, zBox / 20.0);
        TGeoVolume *hitVolume = new TGeoVolume("RearCalVolume", box, air);
        hitVolume->SetLineColor(kBlue);
        TGeoTranslation *trans = new TGeoTranslation(position.X() / 10.0,
                                                     position.Y() / 10.0, (position.Z() + zBox / 2.0) / 10.0);
        rearcal->AddNode(hitVolume, it.moduleID, trans);
    }

    // draw rear Hcal hits
    for (const auto &it : fTcalEvent->rearHCalDeposit)
    {
        ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZRearHCal(it.moduleID);
        double zBox = it.energyDeposit*10.0; // 1mm is 100 MeV
        TGeoShape *box = new TGeoBBox("rearhcalbox", fTcalEvent->geom_detector.rearHCalSizeX / 20.0,
                                      zBox / 20.0, 
                                      fTcalEvent->geom_detector.rearHCalSizeZ / 20.0);
        TGeoVolume *hitVolume = new TGeoVolume("RearHCalVolume", box, air);
        hitVolume->SetLineColor(kRed);
        TGeoTranslation *trans = new TGeoTranslation(position.X() / 10.0,
                                                     (position.Y() + zBox / 2.0) / 10.0, position.Z() / 10.0);
        rearHcal->AddNode(hitVolume, it.moduleID, trans);
    }

    // draw rear Muon Cal hits
    if(true) {
        ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZRearHCal(10);
        double zBox = fTcalEvent->rearMuCalDeposit*10.0; // 1mm is 100 MeV
        TGeoShape *box = new TGeoBBox("rearmucalbox", fTcalEvent->geom_detector.rearHCalSizeX / 20.0,
                                      zBox / 20.0, 
                                      fTcalEvent->geom_detector.rearHCalSizeZ / 20.0);
        TGeoVolume *hitVolume = new TGeoVolume("RearMuCalVolume", box, air);
        hitVolume->SetLineColor(kMagenta);
        TGeoTranslation *trans = new TGeoTranslation(position.X() / 10.0,
                                                (position.Y()+zBox/2.0) / 10.0, position.Z() / 10.0);
        rearMucal->AddNode(hitVolume, 0, trans);
    }

    fCanvas->GetCanvas()->cd();
    delete runText;
    runText = new TText(0.05, 0.95, Form("Seq. Event: %d - Run: %d Event: %d", 
        ievent, POevent->run_number, POevent->event_id));
    runText->SetNDC();
    runText->SetTextSize(0.03);
    runText->Draw();

    delete eventypeText;
    std::ostringstream eventtype;
    int pdgin = POevent->in_neutrino.m_pdg_id;
    switch(pdgin) {
        case -12:
            eventtype << "antinu_e";
            break;
        case 12:
            eventtype << "nu_e";
            break;
        case -14:
            eventtype << "antinu_mu";
            break;
        case  14:
            eventtype << "nu_mu";
            break;
        case -16:
            eventtype << "antinu_tau";
            break;
        case  16:
            eventtype << "nu_tau";
            break;
        default:
            eventtype << " ?? ";
    }
    if(POevent->isES()){
        eventtype << " ES ";
    } else if(POevent->isCC) {
        eventtype << " CC ";
    } else {
        eventtype << " NC ";
    }
    if(POevent->istau) {
        switch(POevent->tau_decaymode) {
            case 1:
                eventtype << "  tau->e  ";
                break;
            case 2:
                eventtype << "  tau->mu  ";
                break;
            case 3:
                eventtype << "  tau->1-prong  ";
                break;
            case 4:
                eventtype << "  tau->rho  ";
                break;
            case 5:
                eventtype << "  tau->3-prong  ";
                break;
            default:
                eventtype << "  tau->??  ";
        }
    }
    eventtype << "    Vtx: x=" << POevent->prim_vx.x() << " y=" << POevent->prim_vx.y() << " z=" << POevent->prim_vx.z() << "  ";
    eventtype << " Target " << POevent->GENIE_vtx_name;
    eventtype << " " << Form("       Etrue:%6.2f GeV", POevent->in_neutrino.m_energy);
    if(POevent->isCharmed()) {
        eventtype << "  Charmed";
    }
    eventypeText = new TText(0.05, 0.9, eventtype.str().c_str());
    eventypeText->SetNDC();
    eventypeText->SetTextSize(0.03);
    eventypeText->Draw();

    delete energyText;
    std::ostringstream energies;
    if(fPORecoEvent->GetPOFullRecoEvent()!=nullptr) {
        energies << Form("Evis:%6.2f GeV", fPORecoEvent->GetPOFullRecoEvent()->TotalEvis());
        energies << Form("  ET:%6.2f GeV", fPORecoEvent->GetPOFullRecoEvent()->TotalET());
    }
    energies << Form("  RearCal:%6.2f GeV", fPORecoEvent->rearCals.rearCalDeposit);
    energies << Form("  RearHCal:%6.2f GeV", fPORecoEvent->rearCals.rearHCalDeposit);
    energies << Form("  RearMuCal:%6.2f MeV", fPORecoEvent->rearCals.rearMuCalDeposit);
    energyText = new TText(0.05, 0.85, energies.str().c_str());
    energyText->SetNDC();
    energyText->SetTextSize(0.03);
    energyText->Draw();

    Draw_event_reco_tracks();
    
    // Draw reconstructed voxels
    Draw_event_reco_voxel(bigbox, air, box);
    gGeoManager->GetTopVolume()->AddNode(ps_reco_voxel,1);

    // Draw reconstructed PSTracks
    ps_tracks = new TGeoVolume("ps_tracks", bigbox, air);
    int ivox = 1;
    Int_t colors[] = {kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink};
    for (const auto& it : fPORecoEvent->fPSTracks) {
        Int_t kcolor = colors[(ivox/6)%6]+ivox%6-3;
        size_t nhits = it.tkhit.size();
        for ( size_t i = 0; i < nhits; i++) {
            long ID = it.tkhit[i].ID;
            ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(ID);
            TGeoTranslation *trans = new TGeoTranslation(position.X() / 10.0, 
                position.Y() / 10.0, position.Z() / 10.0);
            TGeoVolume *hitVolume = new TGeoVolume("HitVolume", box, air);
            hitVolume->SetLineColor(kcolor); 
            ps_tracks->AddNode(hitVolume, ivox++, trans);
        }
    }

    if(toggle_secondary_em)
        gGeoManager->GetTopVolume()->AddNode(secondary_em,1);
    if(toggle_secondary_em)
        gGeoManager->GetTopVolume()->AddNode(secondary_had,1);
    if(toggle_primary_em)
        gGeoManager->GetTopVolume()->AddNode(primary_em,1);
    if(toggle_primary_had)
        gGeoManager->GetTopVolume()->AddNode(primary_had,1);

    gGeoManager->GetTopVolume()->AddNode(si_tracker,1);
    gGeoManager->GetTopVolume()->AddNode(si_vertices,1);
    gGeoManager->GetTopVolume()->AddNode(ps_tracks,1);
    gGeoManager->GetTopVolume()->AddNode(rearcal,1);
    gGeoManager->GetTopVolume()->AddNode(rearHcal,1);
    gGeoManager->GetTopVolume()->AddNode(rearMucal,1);

    TCanvas *c1 = fCanvas_2DPSview->GetCanvas();
    c1->Clear();
    c1->cd();
    c1->Divide(2, 1);
    c1->cd(1);
    gPad->SetLogz();
    gStyle->SetOptStat(0);  // Disable the statistics box
    if(fPORecoEvent -> Get2DViewXPS() != nullptr)
        fPORecoEvent -> Get2DViewXPS() -> Draw("COLZ");
    c1->cd(2);
    gPad->SetLogz();
    gStyle->SetOptStat(0);  // Disable the statistics box
    if(fPORecoEvent -> Get2DViewYPS() != nullptr)
        fPORecoEvent -> Get2DViewYPS() -> Draw("COLZ");
    c1->Modified();
    c1->Update();

    TCanvas *c2 = fCanvas_2DPSview_emhad->GetCanvas();
    c2->Clear();
    c2->Divide(2, 2);
    c2->cd(1);
    gPad->SetLogz();
    gStyle->SetOptStat(0);  // Disable the statistics box
    if(fPORecoEvent -> xviewPS_em != nullptr)
        fPORecoEvent -> xviewPS_em -> Draw("COLZ");
    c2->cd(2);
    gPad->SetLogz();
    gStyle->SetOptStat(0);  // Disable the statistics box
    if(fPORecoEvent -> yviewPS_em != nullptr)
        fPORecoEvent -> yviewPS_em -> Draw("COLZ");
    c2->cd(3);
    gPad->SetLogz();
    gStyle->SetOptStat(0);  // Disable the statistics box
    if(fPORecoEvent -> xviewPS_had != nullptr)
        fPORecoEvent -> xviewPS_had -> Draw("COLZ");
    c2->cd(4);
    gPad->SetLogz();
    gStyle->SetOptStat(0);  // Disable the statistics box
    if(fPORecoEvent -> yviewPS_had != nullptr)
        fPORecoEvent -> yviewPS_had -> Draw("COLZ");
    c2->Modified();
    c2->Update();

    TCanvas *c3 = fCanvas_2DPSviewZ->GetCanvas();
    c3->Clear();
    c3->Divide(5, 4);
    for (int i=0; i<20; i++) {
        c3->cd(i+1);
        gPad->SetLogz();
        gStyle->SetOptStat(0);  // Disable the statistics box
        if(fPORecoEvent -> zviewPS.size() > 0)
            fPORecoEvent -> zviewPS[i] -> Draw("COLZ");
    }
    c3->Modified();
    c3->Update();

    TCanvas *c4 = fCanvas_eldepo->GetCanvas();
    c4->Clear();
    c4->Divide(2, 1);
    c4->cd(1);
    gPad->SetLogz();
    gStyle->SetOptStat(0);  // Disable the statistics box
    if(fPORecoEvent -> xviewPS_eldepo != nullptr)
    {
        fPORecoEvent -> xviewPS_eldepo -> GetXaxis() -> SetTitle("Electromagneticity");
        fPORecoEvent -> xviewPS_eldepo -> GetYaxis() -> SetTitle("Deposited energy (MeV)");
        fPORecoEvent -> xviewPS_eldepo -> Draw("COLZ");
    }
    c4->cd(2);
    gPad->SetLogz();
    gStyle->SetOptStat(0);  // Disable the statistics box
    if(fPORecoEvent -> yviewPS_eldepo != nullptr)
    {
        fPORecoEvent -> yviewPS_eldepo -> GetXaxis() -> SetTitle("Electromagneticity");
        fPORecoEvent -> yviewPS_eldepo -> GetYaxis() -> SetTitle("Deposited energy (MeV)");
        fPORecoEvent -> yviewPS_eldepo -> Draw("COLZ");
    }
    c4->Modified();
    c4->Update();
}

void MyMainFrame::Draw_event_reco_tracks() {
    // draw tracks
    #if 0 // done in ROOT when canvas is cleared
    for(auto &it : polylineTracks) {
        delete it;
    }
    #endif
    polylineTracks.clear();
    polylineVertices.clear();

    if(toggle_reconstructed_tracks) {
#if 0
        for (auto &it : fPORecoEvent->GetPORecs())
        {
            // plot all track of each PORec
            //            struct TPORec::TRACK itrk = it->fTracks[0];
            for (auto &itrk : it->fTracks)
#endif
            for (auto &itrk : fPORecoEvent -> fTKTracks)
            {
                int nhits = itrk.tkhit.size();
                if (nhits == 0)
                    continue;
                double *x = (double *)malloc(nhits * sizeof(double));
                double *y = (double *)malloc(nhits * sizeof(double));
                double *z = (double *)malloc(nhits * sizeof(double));
                int idx = 0;
                for (auto &itrk : itrk.tkhit)
                {
                    x[idx] = itrk.point.x() / 10.0;
                    y[idx] = itrk.point.y() / 10.0;
                    z[idx] = itrk.point.z() / 10.0;
                    idx++;
                }
                TPolyLine3D *trackpoly = new TPolyLine3D(nhits, x, y, z);
                trackpoly->SetLineColor(kBlack);
                trackpoly->SetLineWidth(2);
                if(itrk.trackID == selectedTrackID) {
                    trackpoly->SetLineColor(kRed);
                    trackpoly->SetLineWidth(4);
                }
                trackpoly->Draw("same");
                polylineTracks.push_back(trackpoly);
            }
#if 0
        }
#endif
    }

    if(toggle_reconstructed_tracks) {
        // draw all vertices
        for (auto &it : fPORecoEvent->fTKVertices)
        {
            double x = it.position.x() / 10.0;
            double y = it.position.y() / 10.0;
            double z = it.position.z() / 10.0;
            TPolyMarker3D *polyMarker = new TPolyMarker3D(1);
            polyMarker->SetPoint(0, x, y, z);
            // Set marker properties
            polyMarker->SetMarkerColor(kRed);
            polyMarker->SetMarkerStyle(20); // Full circle
            polyMarker->SetMarkerSize(1.5);
            if(it.vertexID == selectedVertexID) {
                polyMarker->SetMarkerColor(kBlue);
                polyMarker->SetMarkerSize(3.0);
            }
            polyMarker->Draw("same");
            polylineVertices.push_back(polyMarker);
        }
    }
}

void MyMainFrame::Draw_event_reco_voxel(TGeoShape *bigbox, TGeoMedium *air, TGeoShape *box) {
    // Draw reconstructed voxels
    ps_reco_voxel = new TGeoVolume("ps_reco_voxel", bigbox, air);

    int i = 1;
    for (auto& it : fPORecoEvent->PSvoxelmap) {
        long ID = it.first;
        double ehit = it.second.RawEnergy;
        if(ehit < 0.5) continue;
        ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(ID);
//        position += ROOT::Math::XYZVector(2.5,2.5,2.5);    // in mm
        TGeoTranslation *trans = new TGeoTranslation(position.X() / 10.0, 
                position.Y() / 10.0, position.Z() / 10.0);
        TGeoVolume *hitVolume = new TGeoVolume("HitVolume", box, air);
        int color = kBlack;
        if(ehit > 10.0) color = kGreen;
        if(ehit > 100.0) color = kYellow;
        if(it.second.ghost) {
            hitVolume->SetLineColor(toggle_color_fake_voxel ? kMagenta : color); 
        } else {
            hitVolume->SetLineColor(color); 
        }
        ps_reco_voxel->AddNode(hitVolume, i++, trans);
    }
}

void MyMainFrame::Next_Event(int ievent) {

    TCanvas *canvas = fCanvas->GetCanvas();

    // remove the previous event
    TGeoNode *nodeToRemove1 = gGeoManager->GetTopVolume()->FindNode("primary_em_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove1);
    TGeoNode *nodeToRemove2 = gGeoManager->GetTopVolume()->FindNode("primary_had_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove2);
    TGeoNode *nodeToRemove3 = gGeoManager->GetTopVolume()->FindNode("secondary_em_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove3);
    TGeoNode *nodeToRemove4 = gGeoManager->GetTopVolume()->FindNode("secondary_had_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove4);
 
    TGeoNode *nodeToRemove5 = gGeoManager->GetTopVolume()->FindNode("si_tracker_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove5);
    TGeoNode *nodeToRemove11 = gGeoManager->GetTopVolume()->FindNode("si_vertices_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove11);

    TGeoNode *nodeToRemove6 = gGeoManager->GetTopVolume()->FindNode("ps_reco_voxel_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove6);
    TGeoNode *nodeToRemove8 = gGeoManager->GetTopVolume()->FindNode("ps_tracks_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove8);

    TGeoNode *nodeToRemove7 = gGeoManager->GetTopVolume()->FindNode("rearcal_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove7);
    TGeoNode *nodeToRemove9 = gGeoManager->GetTopVolume()->FindNode("rearHcal_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove9);
    TGeoNode *nodeToRemove10 = gGeoManager->GetTopVolume()->FindNode("rearMucal_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove10);

    int run_number = fTcalEvent->fTPOEvent->run_number;
    if(!process_reco_event) {
        delete POevent;
        delete fTcalEvent;
    }
    delete fPORecoEvent;
    fPORecoEvent = nullptr;

    if(!process_reco_event) 
        Load_event(run_number, ievent, event_mask);
    else
        Load_Recoevent(run_number, ievent);

    toggle_primary_em=
    toggle_primary_had=
    toggle_secondary_em=
    toggle_secondary_had=false;

    toggle_reconstructed_tracks = true;
    toggle_reconstructed_ps_tracks = true;
    toggle_reco_voxel = true;

    // remove selected track and vertex
    selectedTrackID = -1;
    selectedVertexID = -1;

    Draw_event();
    Update_ListTree();
    canvas->Modified();
    canvas->Update();
}


// Function to handle button click
void MyMainFrame::HandleButton() {
    // get the first hit
    double zvtx = fTcalEvent->fTPOEvent->prim_vx.z();
    ZoomToPosition(0,0,zvtx/10.0);
}

// Function to handle goto event button click
void MyMainFrame::goto_event() {
    const char* text = textNextEventEntry->GetText();
    printf("Text Entry Content: %s\n", text);
    try {
        ievent = std::stoi(text);
        Next_Event(ievent); 
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid input: not an integer" << std::endl;
    } catch (const std::out_of_range& e) {
        std::cerr << "Invalid input: out of range" << std::endl;
    }
}

// Function to handle next eventbutton click
void MyMainFrame::next_event() {
    Next_Event(++ievent);  
}

void MyMainFrame::toggle_prim_em() {
    TCanvas *canvas = fCanvas->GetCanvas();
    toggle_primary_em = !toggle_primary_em;
    if(toggle_primary_em) {
        gGeoManager->GetTopVolume()->AddNode(primary_em,1);
    } else {
        TGeoNode *nodeToRemove = gGeoManager->GetTopVolume()->FindNode("primary_em_1");
        gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove);
    }
    canvas->Modified();
    canvas->Update();
}

void MyMainFrame::toggle_prim_had() {
    TCanvas *canvas = fCanvas->GetCanvas();
    toggle_primary_had = !toggle_primary_had;
    if(toggle_primary_had) {
        gGeoManager->GetTopVolume()->AddNode(primary_had,1);
    } else {
        TGeoNode *nodeToRemove = gGeoManager->GetTopVolume()->FindNode("primary_had_1");
        gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove);
    }
    canvas->Modified();
    canvas->Update();
}
void MyMainFrame::toggle_sec_em() {
    TCanvas *canvas = fCanvas->GetCanvas();
    toggle_secondary_em = !toggle_secondary_em;
    if(toggle_secondary_em) {
        gGeoManager->GetTopVolume()->AddNode(secondary_em,1);
    } else {
        TGeoNode *nodeToRemove = gGeoManager->GetTopVolume()->FindNode("secondary_em_1");
        gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove);
    }
    canvas->Modified();
    canvas->Update();
}
void MyMainFrame::toggle_sec_had() {
    TCanvas *canvas = fCanvas->GetCanvas();
    toggle_secondary_had = !toggle_secondary_had;
    if(toggle_secondary_had) {
        gGeoManager->GetTopVolume()->AddNode(secondary_had,1);
    } else {
        TGeoNode *nodeToRemove = gGeoManager->GetTopVolume()->FindNode("secondary_had_1");
        gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove);
    }
    canvas->Modified();
    canvas->Update();
}
void MyMainFrame::toggle_reco_track() {
    TCanvas *canvas = fCanvas->GetCanvas();
    toggle_reconstructed_tracks = !toggle_reconstructed_tracks;
/*    if(toggle_reconstructed_tracks) {
        gGeoManager->GetTopVolume()->AddNode(si_tracker,1);
    } else {
        TGeoNode *nodeToRemove = gGeoManager->GetTopVolume()->FindNode("si_tracker_1");
        gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove);
    }
    */
    canvas->cd();
    for(auto &it : polylineTracks) {
        delete it;
    }
    for(auto &it : polylineVertices) {
        delete it;
    }
    polylineTracks.clear();
    polylineVertices.clear();
    Draw_event_reco_tracks();
    canvas->Modified();
    canvas->Update();
}
void MyMainFrame::toggle_color_fakes() {
    TCanvas *canvas = fCanvas->GetCanvas();
    toggle_color_fake_voxel = !toggle_color_fake_voxel;
    double voxelsize = fTcalEvent->geom_detector.fScintillatorVoxelSize/10.0;
    TGeoMedium *air = gGeoManager->GetMedium("AIR");
    TGeoShape *box = new TGeoBBox("box", voxelsize/2.0,voxelsize/2.0,voxelsize/2.0);
    TGeoNode *nodeToRemove = gGeoManager->GetTopVolume()->FindNode("ps_reco_voxel_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove);
    Draw_event_reco_voxel(bigbox, air, box);
    gGeoManager->GetTopVolume()->AddNode(ps_reco_voxel,1);
    canvas->Modified();
    canvas->Update();
}
void MyMainFrame::toggle_recon_ps_tracks() {
    TCanvas *canvas = fCanvas->GetCanvas();
    toggle_reconstructed_ps_tracks = !toggle_reconstructed_ps_tracks;
    if(toggle_reconstructed_ps_tracks) {
        gGeoManager->GetTopVolume()->AddNode(ps_tracks,1);
    } else {
        TGeoNode *nodeToRemove = gGeoManager->GetTopVolume()->FindNode("ps_tracks_1");
        gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove);
    }
    canvas->Modified();
    canvas->Update();
}
void MyMainFrame::toggle_reco_voxels() {
    TCanvas *canvas = fCanvas->GetCanvas();
    toggle_reco_voxel = !toggle_reco_voxel;
    if(toggle_reco_voxel) {
        std::cout << "Adding ps reco voxel..." << std::endl;
        gGeoManager->GetTopVolume()->AddNode(ps_reco_voxel,1);
    } else {
        std::cout << "Removing ps reco voxel..." << std::endl;
        TGeoNode *nodeToRemove = gGeoManager->GetTopVolume()->FindNode("ps_reco_voxel_1");
        gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove);
    }
    canvas->Modified();
    canvas->Update();
}

void MyMainFrame::only_reco() {
    TCanvas *canvas = fCanvas->GetCanvas();
    // remove the previous event
    TGeoNode *nodeToRemove1 = gGeoManager->GetTopVolume()->FindNode("primary_em_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove1);
    TGeoNode *nodeToRemove2 = gGeoManager->GetTopVolume()->FindNode("primary_had_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove2);
    TGeoNode *nodeToRemove3 = gGeoManager->GetTopVolume()->FindNode("secondary_em_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove3);
    TGeoNode *nodeToRemove4 = gGeoManager->GetTopVolume()->FindNode("secondary_had_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove4);
 
 /*   TGeoNode *nodeToRemove5 = gGeoManager->GetTopVolume()->FindNode("si_tracker_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove5);*/

    toggle_primary_em=
    toggle_primary_had=
    toggle_secondary_em=
    toggle_secondary_had=false; 
    toggle_reconstructed_tracks = false;
    toggle_reco_voxel = true;
    canvas->Modified();
    canvas->Update();
}

void MyMainFrame::SideView() {
    TCanvas *canvas = fCanvas->GetCanvas();
    TView *view = (TView *)canvas->GetView();
    if(view == nullptr) return;
    view->SetPsi(90);
    view->SetRange(range_min, range_max);
    canvas->Modified();
    canvas->Update();
}

void MyMainFrame::ZoomToPosition(Double_t x, Double_t y, Double_t z) {
    TCanvas *canvas = fCanvas->GetCanvas();    
    TView *view = (TView *)canvas->GetView();
    view->SetPsi(0);
    view->SetRange(0,0,z-20,0.1,0.1,z+60);
    canvas->Modified();
    canvas->Update();
}

void MyMainFrame::ZoomIn() {
    TCanvas *canvas = fCanvas->GetCanvas();    
    TView *view = (TView *)canvas->GetView();
    view->ZoomIn();
    canvas->Modified();
    canvas->Update();
}
void MyMainFrame::ZoomOut() {
    TCanvas *canvas = fCanvas->GetCanvas();    
    TView *view = (TView *)canvas->GetView();
    view->ZoomOut();
    canvas->Modified();
    canvas->Update();
}
void MyMainFrame::MoveUp() {
    TCanvas *canvas = fCanvas->GetCanvas();    
    TView *view = (TView *)canvas->GetView();
    view->MoveWindow('u');
    canvas->Modified();
    canvas->Update();
}
void MyMainFrame::MoveDown() {
    TCanvas *canvas = fCanvas->GetCanvas();    
    TView *view = (TView *)canvas->GetView();
    view->MoveWindow('i');
    canvas->Modified();
    canvas->Update();
}
void MyMainFrame::MoveLeft() {
    TCanvas *canvas = fCanvas->GetCanvas();    
    TView *view = (TView *)canvas->GetView();
    view->MoveWindow('l');
    canvas->Modified();
    canvas->Update();
}
void MyMainFrame::MoveRight() {
    TCanvas *canvas = fCanvas->GetCanvas();    
    TView *view = (TView *)canvas->GetView();
    view->MoveWindow('h');
    canvas->Modified();
    canvas->Update();
}
void MyMainFrame::on_fullreco_toggle(Bool_t state)
{
    if (state)
    {
        std::cout << "Checkbox is checked." << std::endl;
    }
    else
    {
        std::cout << "Checkbox is unchecked." << std::endl;
    }
}

ClassImp(MyMainFrame)
