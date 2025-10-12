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

#include <TEveManager.h>
#include <TEveGeoNode.h>

// Function to recursively list volumes
void ListVolumes(TGeoVolume* volume, int depth = 0) {
    // Print the volume name with indentation for hierarchy
    for (int i = 0; i < depth; ++i) std::cout << "  "; // Indent based on depth
    std::cout << volume->GetName();
    std::cout << "  Volume type: " << volume->GetShape()->GetTitle() << std::endl;
    std::cout << "  Volume material: " << volume->GetMaterial()->GetName() << std::endl;
    std::cout << "  Number of daughters: " << volume->GetNdaughters() << std::endl;
    std::cout << "  Volume dimensions: ";
    volume->GetShape()->Print(); // Print shape details
    std::cout << std::endl;

    // Loop through the list of daughter volumes
    int nDaughters = volume->GetNdaughters();
    for (int i = 0; i < nDaughters; ++i) {
        TGeoVolume* daughter = volume->GetNode(i)->GetVolume();
        ListVolumes(daughter, depth + 1);  // Recursive call with increased depth
    }
}

void load_geometry(std::string geometryFile) {
    // Load the GDML geometry
    TGeoManager::Import(geometryFile.c_str());

    TGeoVolume* topVolume = gGeoManager->GetTopVolume();
    topVolume->SetTransparency(50);  // Value from 0 (opaque) to 100 (fully transparent)

    // List all volumes in the geometry
//    std::cout << "Listing all volumes in the geometry:" << std::endl;
//    ListVolumes(topVolume);

    TGeoVolume *volume1 = gGeoManager->FindVolumeFast("ShortCylLogical");
//    volume1->SetLineColor(kRed);      // Set color (optional)
//    volume1->SetVisibility(kTRUE);    // Ensure it is visible
//    volume1->SetDrawOption("wireframe");
#if 0
    TEveManager::Create();

//    GeoTopVolume *topVolume = dynamic_cast<GeoTopVolume*>(gGeoManager->GetTopVolume());
    TEveGeoNode *FASERCAL = new TEveGeoNode(gGeoManager->GetTopNode());
    gEve->AddGlobalElement(FASERCAL);
    gEve->FullRedraw3D(true);
#endif
}


int main(int argc, char** argv) {

    // get the run number as the first argument
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " [-g <geometryfile>] [-r] <run> [mask]" << std::endl;
        std::cout << "   -g <geometryfile>         Load a specific geometry file" << std::endl;
        std::cout << "   -r                        Open reconstructed files" << std::endl;
        std::cout << "   <run>                     Run number" << std::endl;
        std::cout << "   mask                      To process only specific events (def=none): ";
        std::cout << "  nueCC, numuCC, nutauCC, nuNC or nuES" << std::endl;
		return 1;
	}

    bool pre = false;
    int idx = 1;

    if(strcmp(argv[idx], "-r")==0) {
        pre = true;
        idx++;
    }

    // -g <geometryfile>  to load a specific geometry file
    std::string geometryFile = "../GeomGDML/geometry.gdml";
    if (strcmp(argv[idx], "-g") == 0) {
        idx++;
        if (argc < idx + 1) {
            std::cerr << "Error: -g option requires a geometry file argument." << std::endl;
            return 1;
        }
        geometryFile = argv[idx];
        idx++;
    } 

    // check number of arguments on command line
    if (argc < idx + 1) {
        std::cerr << "Error: Missing run number argument." << std::endl;
        return 1;
    }

    std::string runString = argv[idx];
    int run_number;

    try {
        run_number = std::stoi(runString);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument for run: " << e.what() << std::endl;
    } catch (const std::out_of_range& e) {
        std::cerr << "Out of range for run: " << e.what() << std::endl;
    }

    int event_mask = 0;
    if(argc>idx+1) {
        int mask = TPOEvent::EncodeEventMask(argv[idx]);
        if(mask>0) {
            event_mask = mask;
        } else {
            std::cerr << "Unknown mask " << argv[idx] << std::endl;
            exit(1);            
        }
    }

    TApplication app("app", &argc, argv);

    load_geometry(geometryFile);
#if 0
    TGeoVolume *specificVolume = gGeoManager->FindVolumeFast("rearCalmoduleLogical");
    if (specificVolume) {
    specificVolume->SetLineColor(kRed);      
    specificVolume->SetTransparency(50);       
    }
#endif

    new MyMainFrame(run_number, 0, event_mask, pre, gClient->GetRoot(), 1700, 800);

    // Run the application
    app.Run();

    return 0;
}
