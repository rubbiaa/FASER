#include <TApplication.h>
#include <TCanvas.h>
#include <TGeoManager.h>
#include <TChain.h>
#include <TSystem.h>
#include "TGeoVolume.h"
#include "TGeoMatrix.h"

#include "TcalEvent.hh"
#include "TPOEvent.hh"

#include "MyMainFrame.h"

TcalEvent* fTcalEvent;

void load_geometry() {
    // Load the GDML geometry
    TGeoManager::Import("/home/rubbiaa/faserps/FASERPS/geometry.gdml");

    // Draw the geometry
//    gGeoManager->GetTopVolume()->Draw("ogl");
}

void load_event() {
    TChain *event_tree = new TChain("calEvent");
    event_tree->Add("/home/rubbiaa/faserps/FASERPS/output/tcalevent_0.root");
    std::cout << "Number of entries " << event_tree->GetEntries() << std::endl;

    // Create an instance of TcalEvent
    fTcalEvent = new TcalEvent();

    // Set the branch address
    std::vector<DigitizedTrack*> *t;
    t = new std::vector<DigitizedTrack*>;
    event_tree->SetBranchAddress("tracks", &t);

    TPOEvent *POevent;
    POevent = new TPOEvent();
    event_tree -> SetBranchAddress("event", &POevent);
    fTcalEvent->fTPOEvent = POevent;
 
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

    fTcalEvent -> fTPOEvent -> dump_event();

    new MyMainFrame(fTcalEvent, gClient->GetRoot(), 800, 600);

    // Run the application
    app.Run();

    return 0;
}
