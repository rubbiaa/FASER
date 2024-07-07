#include <TCanvas.h>
#include <TGeoManager.h>
#include "TGeoVolume.h"
#include <TPolyMarker3D.h>
#include "TGeoSphere.h"
#include "TGeoBBox.h"
#include "TGeoMedium.h"
#include <TView3D.h>
#include <TChain.h>
#include <TText.h>

#include "MyMainFrame.h"
#include "TPORecoEvent.hh"

MyMainFrame::MyMainFrame(int ieve, const TGWindow *p, UInt_t w, UInt_t h) {

    // load event
    ievent = ieve;
    Load_event(ievent);

    fMain = new TGMainFrame(p, w, h);

    // Create an embedded canvas
    fCanvas = new TRootEmbeddedCanvas("EmbeddedCanvas", fMain, 800, 600);
    fMain->AddFrame(fCanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    // Create a button
    fButton = new TGTextButton(fMain, "&Zoom vtx");
    fButton->Connect("Clicked()", "MyMainFrame", this, "HandleButton()");
    fMain->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsBottom, 5, 5, 3, 4));
    fButton = new TGTextButton(fMain, "Toggle prim_em");
    fButton->Connect("Clicked()", "MyMainFrame", this, "toggle_prim_em()");
    fMain->AddFrame(fButton, new TGLayoutHints(kLHintsBottom, 5, 5, 3, 4));
    fButton = new TGTextButton(fMain, "Toggle prim_had");
    fButton->Connect("Clicked()", "MyMainFrame", this, "toggle_prim_had()");
    fMain->AddFrame(fButton, new TGLayoutHints(kLHintsBottom, 5, 5, 3, 4));
    fButton = new TGTextButton(fMain, "Toggle sec_em");
    fButton->Connect("Clicked()", "MyMainFrame", this, "toggle_sec_em()");
    fMain->AddFrame(fButton, new TGLayoutHints(kLHintsBottom, 5, 5, 3, 4));
    fButton = new TGTextButton(fMain, "Toggle sec_had");
    fButton->Connect("Clicked()", "MyMainFrame", this, "toggle_sec_had()");
    fMain->AddFrame(fButton, new TGLayoutHints(kLHintsBottom, 5, 5, 3, 4));
    fButton = new TGTextButton(fMain, "Next");
    fButton->Connect("Clicked()", "MyMainFrame", this, "next_event()");
    fMain->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsBottom, 5, 5, 3, 4));

    fMain->SetWindowName("The FASERkine event display");
    fMain->MapSubwindows();
    fMain->Resize(fMain->GetDefaultSize());
    fMain->MapWindow();

    TCanvas *canvas = fCanvas->GetCanvas();
    canvas->cd();
    // Draw the geometry
    gGeoManager->GetTopVolume()->Draw("ogl");

    Draw_event();

    toggle_primary_em=
    toggle_primary_had=
    toggle_secondary_em=
    toggle_secondary_had=true;
}

// Destructor
MyMainFrame::~MyMainFrame() {
    fMain->Cleanup();
    delete fMain;
}

void MyMainFrame::Load_event(int ievent) {

    std::string base_path = "input/tcalevent_";

    // Create an instance of TcalEvent and TPOEvent
    fTcalEvent = new TcalEvent();
    POevent = new TPOEvent();

    fTcalEvent -> Load_event(base_path, ievent, POevent);
    std::cout << "Transverse size " << fTcalEvent->geom_detector.fScintillatorSizeX << " mm " << std::endl;
    std::cout << "Total size of one sandwich layer " << fTcalEvent->geom_detector.fTotalLength << " mm " << std::endl;
	std::cout << "Number of layers " << fTcalEvent->geom_detector.NRep << std::endl;
    std::cout << "Voxel size " << fTcalEvent->geom_detector.fScintillatorVoxelSize << " mm " << std::endl;

    std::cout << " copied digitized tracks " << fTcalEvent->fTracks.size() << std::endl;

    fTcalEvent -> fTPOEvent -> dump_event();

    fPORecoEvent = new TPORecoEvent(fTcalEvent, fTcalEvent->fTPOEvent);
    fPORecoEvent -> Reconstruct();
    fPORecoEvent -> Dump();
}

void MyMainFrame::Draw_event() {

    double wdx, wdy, wdz;
    if (TGeoBBox *box = dynamic_cast<TGeoBBox*>(gGeoManager->GetTopVolume()->GetShape())) {
        wdx = box->GetDX();
        wdy = box->GetDY();
        wdz = box->GetDZ()/2.0;
        printf("Top volume dimensions: DX = %f, DY = %f, DZ = %f\n", wdx, wdy, wdz);
    } else {
        exit(1);
    }
    TGeoShape *bigbox = new TGeoBBox("bigbox", wdx, wdy, wdz);

    TGeoMedium *air = gGeoManager->GetMedium("AIR");
    primary_em = new TGeoVolume("primary_em", bigbox, air);
    primary_had = new TGeoVolume("primary_had", bigbox, air);
    secondary_em = new TGeoVolume("secondary_em", bigbox, air);
    secondary_had = new TGeoVolume("secondary_had", bigbox, air);

    TGeoShape *hitShape = new TGeoSphere("HitShape", 0, 0.5);
 
    TGeoMaterial *matAluminum = new TGeoMaterial("Aluminum", 26.98, 13, 2.7);
    TGeoMedium *aluminum = new TGeoMedium("Aluminum", 2, matAluminum);
    TGeoShape *box = new TGeoBBox("box", 0.5/2.0,0.5/2.0,0.5/2.0);

    for (const auto& track : fTcalEvent->fTracks) {
//        std::cout << track->ftrackID << std::endl;
        size_t nhits = track->fhitIDs.size();
//        std::cout << nhits << std::endl;
//        if(track->fparentID > 0) continue;
        for ( size_t i = 0; i < nhits; i++) {
            if(track->fEnergyDeposits[i] < 0.5)continue;
            XYZVector position = fTcalEvent->getChannelXYZfromID(track->fhitIDs[i]);
//            marker->SetNextPoint(position.X()/10.0-12.5, position.Y()/10.0-12.5, position.Z()/10.0-60.0);

            TGeoVolume *hitVolume = new TGeoVolume("HitVolume", box, air);
            hitVolume->SetLineColor(kRed); 
            if(fabs(track->fPDG) == 11){
                hitVolume->SetLineColor(kBlue); // electromagnetic is blue
            } else if(fabs(track->fPDG) == 13){
                hitVolume->SetLineColor(kGreen); // muons
            }

        // Create a translation matrix for the hit position
            TGeoTranslation *trans = new TGeoTranslation(position.X() / 10.0, position.Y() / 10.0, position.Z() / 10.0);

        // Add the hit volume to the top volume with the translation
            if(track->fparentID == 0) {
                if(fabs(track->fPDG) == 11) {
                    primary_em->AddNode(hitVolume, i, trans);
                } else {
                    primary_had->AddNode(hitVolume, i, trans);
                }
            } else {
                if(fabs(track->fPDG) == 11) {
                    secondary_em->AddNode(hitVolume, i, trans);
                } else {
                    secondary_had->AddNode(hitVolume, i, trans);
                }
            }

        }
    }

    gGeoManager->GetTopVolume()->AddNode(secondary_em,1);
    gGeoManager->GetTopVolume()->AddNode(secondary_had,1);
    gGeoManager->GetTopVolume()->AddNode(primary_em,1);
    gGeoManager->GetTopVolume()->AddNode(primary_had,1);
//    gGeoManager->GetTopVolume()->Print();

    delete runText;
    runText = new TText(0.05, 0.95, Form("Run: %d Event: %d", 
        POevent->run_number, POevent->event_id));
    runText->SetNDC();
    runText->SetTextSize(0.03);
    runText->Draw();

    delete eventypeText;
    std::ostringstream eventtype;
    int pdgin = POevent->in_neutrino.m_pdg_id;
    switch(pdgin) {
        case -12:
        case 12:
            eventtype << "nu_e";
            break;
        case -14:
        case  14:
            eventtype << "nu_mu";
            break;
        case -16:
        case  16:
            eventtype << "nu_tau";
            break;
        default:
            eventtype << " ?? ";
    }
    if(POevent->isCC) {
        eventtype << " CC ";
    } else {
        eventtype << " NC ";
    }
    eventypeText = new TText(0.05, 0.9, eventtype.str().c_str());
    eventypeText->SetNDC();
    eventypeText->SetTextSize(0.03);
    eventypeText->Draw();

    delete energyText;
    energyText = new TText(0.05, 0.85, Form("%f GeV", 
        POevent->in_neutrino.m_energy));
    energyText->SetNDC();
    energyText->SetTextSize(0.03);
    energyText->Draw();
}


// Function to handle button click
void MyMainFrame::HandleButton() {
    // get the first hit
    XYZVector pos = fTcalEvent->getChannelXYZfromID(0);
    ZoomToPosition(0,0,pos.Z()/10.0);
}

// Function to handle button click
void MyMainFrame::next_event() {
    TCanvas *canvas = fCanvas->GetCanvas();
    TGeoNode *nodeToRemove1 = gGeoManager->GetTopVolume()->FindNode("primary_em_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove1);
    TGeoNode *nodeToRemove2 = gGeoManager->GetTopVolume()->FindNode("primary_had_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove2);
    TGeoNode *nodeToRemove3 = gGeoManager->GetTopVolume()->FindNode("secondary_em_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove3);
    TGeoNode *nodeToRemove4 = gGeoManager->GetTopVolume()->FindNode("secondary_had_1");
    gGeoManager->GetTopVolume()->RemoveNode(nodeToRemove4);

    delete POevent;
    delete fTcalEvent;
    delete fPORecoEvent;

    Load_event(++ievent);
    Draw_event();
    canvas->Modified();
    canvas->Update();
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

void MyMainFrame::ZoomToPosition(Double_t x, Double_t y, Double_t z) {
    TCanvas *canvas = fCanvas->GetCanvas();
    
    // Define the range for the view manually
    Double_t viewRange[6] = {-10, -10, -10, 10, 10, 10}; // xmin, ymin, zmin, xmax, ymax, zmax
    
    TView *view = (TView *)canvas->GetView();
    view->SetRange(0,0,z-10,0.1,0.1,z+30);
    canvas->Modified();
    canvas->Update();
}

ClassImp(MyMainFrame)
