#include <TApplication.h>
#include <TCanvas.h>
#include <TGeoManager.h>

#include <TSystem.h>
#include "TGeoVolume.h"
#include "TGeoMatrix.h"

#include "TcalEvent.hh"
#include "TPOEvent.hh"

#include "MyMainFrame.h"

#include <sstream>
#include <iostream>

void load_geometry() {
    // Load the GDML geometry
    TGeoManager::Import("../GeomGDML/geometry.gdml");

    // Draw the geometry
//    gGeoManager->GetTopVolume()->Draw("ogl");
}


int main(int argc, char** argv) {

    TApplication app("app", &argc, argv);

    load_geometry();

    new MyMainFrame(0, gClient->GetRoot(), 800, 600);

    // Run the application
    app.Run();

    return 0;
}
