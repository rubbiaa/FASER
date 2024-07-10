// Interactive event display of FASERCAL events
//
// A. Rubbia, July 2024
//

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

    // get the run number as the first argument
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " <run>" << std::endl;
		return 1;
	}

    std::string runString = argv[1];
    int run_number;

    try {
        run_number = std::stoi(runString);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument for run: " << e.what() << std::endl;
    } catch (const std::out_of_range& e) {
        std::cerr << "Out of range for run: " << e.what() << std::endl;
    }

    TApplication app("app", &argc, argv);

    load_geometry();

    new MyMainFrame(run_number, 0, gClient->GetRoot(), 1200, 600);

    // Run the application
    app.Run();

    return 0;
}
