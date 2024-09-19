// Create file Masks of FASERCAL events
//
// A. Rubbia, September 2024
//

#include <string>

#include "TFile.h"
#include "TChain.h"

#include "TPORecoEvent.hh"
#include "TFileMask.hh"

int main(int argc, char** argv) {

	if (argc < 2) {
	std::cout << "Usage: " << argv[0] << " <run> " << std::endl;
        std::cout << "   <run>                     Run number" << std::endl;
		return 1;
	}

   // get the run number as the first argument
    std::string runString = argv[1];
    int run_number;

    try {
        run_number = std::stoi(runString);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument for run: " << e.what() << std::endl;
        exit(1);
    } catch (const std::out_of_range& e) {
        std::cerr << "Out of range for run: " << e.what() << std::endl;
        exit(1);
    }

    std::ostringstream inputfilename;
    inputfilename << "input/Batch-TPORecevent_" << run_number << "_*_*.root";

    TChain *event_tree = new TChain("RecoEvent","READ");
    event_tree->Add(inputfilename.str().c_str());

    Long_t nentries = event_tree->GetEntries();
    std::cout << "Number of entries " << nentries << std::endl;

    // TPORecoEvent class and histograms
    TPORecoEvent *fTPORecoEvent = nullptr; // &fTPORecoEvent;
    event_tree -> SetBranchAddress("TPORecoEvent", &fTPORecoEvent);

    // the file masks
    TFileMask charm_fileMask(run_number, "charm");
    TFileMask tau_fileMask(run_number, "tau");

    // process events
    int ievent = 0;
    int max_event = -1;
    if(max_event == -1) max_event = nentries;
    int error = 0;

    TPOEvent myPOevent; // local copy; just a temporary PO event to store stats
    myPOevent.reset_stats();

    while (error == 0 && ievent<max_event) {

        event_tree->GetEntry(ievent);

        if(ievent % 1000 == 0) {
            std::cout << "Processing event " << ievent << std::endl;
        }
        bool dump_event_cout = (ievent < 20);
        if(dump_event_cout) {
            fTPORecoEvent -> GetPOEvent()->dump_event();
        }
        myPOevent.clone(fTPORecoEvent -> GetPOEvent());
        myPOevent.update_stats();

        int event = fTPORecoEvent->GetPOEvent()->event_id;
        if(fTPORecoEvent -> GetPOEvent() -> isCharmed()) {
            charm_fileMask.addEvent(event);
        }
        if(fTPORecoEvent -> GetPOEvent() -> istau) {
            tau_fileMask.addEvent(event);
        }

        ievent++;
    }

    charm_fileMask.Dump();
    tau_fileMask.Dump();

    myPOevent.dump_stats();

    std::cout << "I'm done." << std::endl;

    return 0;
 
}
