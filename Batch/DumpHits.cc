// Dump hits of FASERCAL events
//
// A. Rubbia, July 2024
//

#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <TGeoManager.h>

#include "TcalEvent.hh"
#include "TPORecoEvent.hh"
#include "TTauSearch.hh"

void load_geometry()
{
    // Load the GDML geometry
    TGeoManager::Import("../GeomGDML/geometry.gdml");
}

int main(int argc, char **argv)
{

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <run> [maxevent] [mask]" << std::endl;
        std::cout << "   <run>                     Run number" << std::endl;
        std::cout << "   maxevent                  Maximum number of events to process (def=-1)" << std::endl;
        std::cout << "   mask                      To process only specific events (def=none): ";
        std::cout << "  nueCC, numuCC, nutauCC, or nuNC" << std::endl;
        return 1;
    }

    // get the run number as the first argument
    std::string runString = argv[1];
    int run_number;

    try
    {
        run_number = std::stoi(runString);
    }
    catch (const std::invalid_argument &e)
    {
        std::cerr << "Invalid argument for run: " << e.what() << std::endl;
        exit(1);
    }
    catch (const std::out_of_range &e)
    {
        std::cerr << "Out of range for run: " << e.what() << std::endl;
        exit(1);
    }

    int max_event = -1;
    if (argc > 2)
    {
        try
        {
            max_event = std::stoi(argv[2]);
        }
        catch (const std::invalid_argument &e)
        {
            std::cerr << "Invalid argument for maxevent: " << e.what() << std::endl;
            exit(1);
        }
        catch (const std::out_of_range &e)
        {
            std::cerr << "Out of range for maxevent: " << e.what() << std::endl;
            exit(1);
        }
    }

    int event_mask = 0;
    if (argc > 3)
    {
        int mask = TPOEvent::EncodeEventMask(argv[3]);
        if (mask > 0)
        {
            event_mask = mask;
        }
        else
        {
            std::cerr << "Unknown mask " << argv[3] << std::endl;
            exit(1);
        }
    }

    load_geometry();

    std::string base_path = "input/";

    std::ostringstream filename;
    filename << "Batch-TPORecevent_" << run_number;
    if (event_mask > 0)
    {
        const char *mask = TPOEvent::DecodeEventMask(event_mask);
        filename << "_" << mask;
    }
    filename << ".root";

    // process events
    int ievent = 0;
    if (max_event == -1)
        max_event = 99999999;
    int error = 0;

    while (error == 0 && ievent < max_event)
    {

        if (ievent % 1000 == 0)
        {
            std::cout << "Processing event " << ievent << std::endl;
        }
        bool dump_event_cout = (ievent < 20);

        // Create an instance of TcalEvent and TPOEvent
        TcalEvent *fTcalEvent = new TcalEvent();
        TPOEvent *POevent = new TPOEvent();

        if (!dump_event_cout)
            fTcalEvent->SetVerbose(0);
        error = fTcalEvent->Load_event(base_path, run_number, ievent++, event_mask, POevent);
        if (error != 0)
            break;

        if (dump_event_cout)
            fTcalEvent->fTPOEvent->dump_event();

        TPORecoEvent *fPORecoEvent = new TPORecoEvent(fTcalEvent, fTcalEvent->fTPOEvent);
        fPORecoEvent->Reconstruct();
        if (dump_event_cout)
            fPORecoEvent->Dump();

        TPORecoEvent *branch_TPORecoEvent = fPORecoEvent;

        for (const auto &track : fTcalEvent->getfTracks())
        {
            //        std::cout << track->ftrackID << std::endl;
            size_t nhits = track->fhitIDs.size();
            //        std::cout << nhits << std::endl;
            //        if(track->fparentID > 0) continue;
            for (size_t i = 0; i < nhits; i++)
            {

                long hittype = fTcalEvent->getChannelTypefromID(track->fhitIDs[i]);

                // apply energy cut and keep only scintillator voxel
                if (hittype != 0 || track->fEnergyDeposits[i] < 0.5)
                    continue;

                ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(track->fhitIDs[i]);
                std::cout << std::setw(10) << track->ftrackID << " ";
                std::cout << std::setw(10) << track->fPDG << " ";
                std::cout << std::setw(10) << track->fparentID << " ";
                std::cout << std::setw(10) << track->fparentID << " ";
                std::cout << std::setw(10) << position.x() << " " << std::setw(10) << position.y() << " " << std::setw(10) <<  position.z() << " ";
                std::cout << std::setw(10) << track->fEnergyDeposits[i];
                std::cout << std::endl;
            }
        }

        delete fPORecoEvent;
        delete POevent;
        delete fTcalEvent;
    }

    std::cout << "I'm done." << std::endl;

    return 0;
}
