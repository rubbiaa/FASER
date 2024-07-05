#include <TApplication.h>
#include <TCanvas.h>
#include <TGeoManager.h>
#include <TChain.h>
#include <TSystem.h>
#include <TPolyMarker3D.h>
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TGeoSphere.h"
#include "TGeoBBox.h"

// #include "CreateGeometry.h"
// #include "EventData.h"

#include "TcalEvent.hh"

TcalEvent* fTcalEvent;

void load_geometry() {
    // Load the GDML geometry
    TGeoManager::Import("/home/rubbiaa/faserps/FASERPS/geometry.gdml");

    // Draw the geometry
//    gGeoManager->GetTopVolume()->Draw("ogl");
}

void load_event() {
    TChain *event_tree = new TChain("calEvent");
    event_tree->Add("/home/rubbiaa/faserps/FASERPS/output/tcalevent_2.root");
    std::cout << "Number of entries " << event_tree->GetEntries() << std::endl;

    // Create an instance of TcalEvent
    fTcalEvent = new TcalEvent();

    // Set the branch address
    std::vector<DigitizedTrack*> *t;
    t = new std::vector<DigitizedTrack*>;
    event_tree->SetBranchAddress("tracks", &t);
 
    // Read the first entry
    event_tree->GetEntry(0);

    // Use the loaded data (example)
    std::cout << "Loaded event data for event 0" << std::endl;
    std::cout << " digitized tracks " << t->size() << std::endl;

    for (const auto& track : *t) {
        fTcalEvent->fTracks.push_back(track);
    }
    std::cout << " copied digitized tracks " << fTcalEvent->fTracks.size() << std::endl;
}

int main(int argc, char** argv) {

    TApplication app("app", &argc, argv);

    load_geometry();

    // load event
    load_event();

    // Draw the event data
 //   eventData.Draw();

    // Draw the geometry
    gGeoManager->GetTopVolume()->Draw("ogl");

    TGeoShape *hitShape = new TGeoSphere("HitShape", 0, 0.5);
    TGeoMedium *air = gGeoManager->GetMedium("AIR");

    for (const auto& track : fTcalEvent->fTracks) {
        TPolyMarker3D* marker = new TPolyMarker3D();
//        std::cout << track->ftrackID << std::endl;
        size_t nhits = track->fhitIDs.size();
//        std::cout << nhits << std::endl;
//        if(track->fparentID > 0) continue;
        for ( size_t i = 0; i < nhits; i++) {
            if(track->fEnergyDeposits[i] < 0.5)continue;
            XYZVector position = fTcalEvent->getChannelXYZfromID(track->fhitIDs[i]);
            marker->SetNextPoint(position.X()/10.0-12.5, position.Y()/10.0-12.5, position.Z()/10.0);

//        TGeoVolume *hitVolume = new TGeoVolume("HitVolume", hitShape, air);
//        hitVolume->SetLineColor(kBlue); // Set color to blue

        // Create a translation matrix for the hit position
//        TGeoTranslation *trans = new TGeoTranslation(position.X() / 10.0, position.Y() / 10.0, position.Z() / 10.0);

        // Add the hit volume to the top volume with the translation
//        gGeoManager->GetTopVolume()->AddNode(hitVolume, i, trans);

        }
        marker->SetMarkerColor(kRed);
        if(fabs(track->fPDG) == 11) marker->SetMarkerColor(kBlue);
        marker->SetMarkerSize(0.5);
        marker->SetMarkerStyle(4);  // Full circle
        marker->Draw("same");
    }
    // Run the application
    app.Run();

    return 0;
}
