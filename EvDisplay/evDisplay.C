// Interactive event display of FASERCAL events
//
// A. Rubbia, July 2024
//

#include <sstream>
#include <iostream>

#include <TApplication.h>
#include <TCanvas.h>
#include <TGeoManager.h>

#include <TSystem.h>
#include "TGeoVolume.h"
#include "TGeoMatrix.h"

#include "TcalEvent.hh"
#include "TPOEvent.hh"

#include "MyMainFrame.h"

// #include <TEveManager.h>

void load_geometry() {
    // Load the GDML geometry
    TGeoManager::Import("../GeomGDML/geometry.gdml");

#if 0
    TEveManager::Create();
    TEveGeoTopNode *FASERCAL = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
    gEve->AddGlobalElement(FASERCAL);
    gEve->FullRedraw3D(kTrue)
#endif
}


int main(int argc, char** argv) {

    // get the run number as the first argument
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " <run> [mask]" << std::endl;
        std::cout << "   <run>                     Run number" << std::endl;
        std::cout << "   mask                      To process only specific events (def=none): ";
        std::cout << "  nueCC, numuCC, nutauCC, or nuNC" << std::endl;
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

    int event_mask = 0;
    if(argc>1) {
        int mask = TPOEvent::EncodeEventMask(argv[2]);
        if(mask>0) {
            event_mask = mask;
        } else {
            std::cerr << "Unknown mask " << argv[2] << std::endl;
            exit(1);            
        }
    }

    TApplication app("app", &argc, argv);

    load_geometry();

    new MyMainFrame(run_number, 0, event_mask, gClient->GetRoot(), 1200, 600);

    // Run the application
    app.Run();

    return 0;
}
