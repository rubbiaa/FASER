#include <TCanvas.h>
#include <TGeoManager.h>
#include "TGeoVolume.h"
#include <TPolyMarker3D.h>
#include "TGeoSphere.h"
#include "TGeoBBox.h"
#include "TGeoMedium.h"
#include <TView3D.h>


#include "MyMainFrame.h"

MyMainFrame::MyMainFrame(TcalEvent *ft, const TGWindow *p, UInt_t w, UInt_t h) : fTcalEvent(ft) {
    fMain = new TGMainFrame(p, w, h);

    // Create an embedded canvas
    fCanvas = new TRootEmbeddedCanvas("EmbeddedCanvas", fMain, 800, 600);
    fMain->AddFrame(fCanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    // Create a button
    fButton = new TGTextButton(fMain, "&Click me");
    fButton->Connect("Clicked()", "MyMainFrame", this, "HandleButton()");
    fMain->AddFrame(fButton, new TGLayoutHints(kLHintsCenterX | kLHintsBottom, 5, 5, 3, 4));

    fMain->SetWindowName("The FASERkine event display");
    fMain->MapSubwindows();
    fMain->Resize(fMain->GetDefaultSize());
    fMain->MapWindow();

    TCanvas *canvas = fCanvas->GetCanvas();
    canvas->cd();
    // Draw the geometry
    gGeoManager->GetTopVolume()->Draw("ogl");

    TGeoShape *hitShape = new TGeoSphere("HitShape", 0, 0.5);
    TGeoMedium *air = gGeoManager->GetMedium("AIR");

    TGeoMaterial *matAluminum = new TGeoMaterial("Aluminum", 26.98, 13, 2.7);
    TGeoMedium *aluminum = new TGeoMedium("Aluminum", 2, matAluminum);
    TGeoShape *box = new TGeoBBox("box", 0.5/2.0,0.5/2.0,0.5/2.0);

    for (const auto& track : fTcalEvent->fTracks) {
        TPolyMarker3D* marker = new TPolyMarker3D();
//        std::cout << track->ftrackID << std::endl;
        size_t nhits = track->fhitIDs.size();
//        std::cout << nhits << std::endl;
//        if(track->fparentID > 0) continue;
        for ( size_t i = 0; i < nhits; i++) {
            if(track->fEnergyDeposits[i] < 0.5)continue;
            XYZVector position = fTcalEvent->getChannelXYZfromID(track->fhitIDs[i]);
            marker->SetNextPoint(position.X()/10.0-12.5, position.Y()/10.0-12.5, position.Z()/10.0-60.0);

            TGeoVolume *hitVolume = new TGeoVolume("HitVolume", box, air);
            hitVolume->SetLineColor(kBlue); // Set color to blue

        // Create a translation matrix for the hit position
            TGeoTranslation *trans = new TGeoTranslation(position.X() / 10.0-12.5, position.Y() / 10.0-12.5, position.Z() / 10.0);

        // Add the hit volume to the top volume with the translation
            gGeoManager->GetTopVolume()->AddNode(hitVolume, i, trans);

        }
        marker->SetMarkerColor(kRed);
        if(fabs(track->fPDG) == 11) marker->SetMarkerColor(kBlue);
        marker->SetMarkerSize(0.5);
        marker->SetMarkerStyle(4);  // Full circle
//        marker->Draw("same");
    }

}

// Destructor
MyMainFrame::~MyMainFrame() {
    fMain->Cleanup();
    delete fMain;
}

// Function to handle button click
void MyMainFrame::HandleButton() {
    printf("Button clicked!\n");
    ZoomToPosition(0,0,-30);
}

void MyMainFrame::ZoomToPosition(Double_t x, Double_t y, Double_t z) {
    TCanvas *canvas = fCanvas->GetCanvas();
    
    // Define the range for the view manually
    Double_t viewRange[6] = {-10, -10, -10, 10, 10, 10}; // xmin, ymin, zmin, xmax, ymax, zmax
    
    TView *view = (TView *)canvas->GetView();
    view->SetRange(0,0,-70,0.1,0.1,-40);
    canvas->Update();
}

ClassImp(MyMainFrame)
