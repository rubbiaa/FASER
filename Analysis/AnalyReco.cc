// Analysis of reconstructed events of FASERCAL events
//
// A. Rubbia, September 2024
//

#include <string>

#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include <TGeoManager.h>

#include "TcalEvent.hh"
#include "TPORecoEvent.hh"
#include "TTauSearch.hh"
#include "TParticleGun.hh"

struct EVENT {
    // Truth variables
    Int_t t_reaction; // reaction 1=nue, 2=numu, 3=nutau, +10 if NC, 20=ES (truth)
    Float_t t_Eneutrino; // true incoming neutrino energy
    // reconstructed quantities
    Int_t n_clusters;
    Float_t c_E1;    // energy most energetic cluster
    Float_t c_chi2_1;
    Float_t c_a_1;
    Float_t c_b_1;
    Float_t c_E2;    // energy 2nd most energetic cluster
    Float_t rear_Cal; // energy deposit in rear Cal
    Float_t rear_MuCal; // energy deposit in rear muCal
};

int main(int argc, char** argv) {

    struct EVENT event;

	if (argc < 2) {
	std::cout << "Usage: " << argv[0] << " <run> [maxevent] [mask]" << std::endl;
        std::cout << "   <run>                     Run number" << std::endl;
        std::cout << "   maxevent                  Maximum number of events to process (def=-1)" << std::endl;
        std::cout << "   mask                      To process only specific events (def=none): ";
        std::cout << "  nueCC, numuCC, nutauCC, nuNC or nuES" << std::endl;
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

    int max_event = -1;
    if(argc>2) {
        try {
            max_event = std::stoi(argv[2]);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument for maxevent: " << e.what() << std::endl;
            exit(1);
        } catch (const std::out_of_range& e) {
            std::cerr << "Out of range for maxevent: " << e.what() << std::endl;
            exit(1);
        }
    }

    int event_mask = 0;
    if(argc>3) {
        int mask = TPOEvent::EncodeEventMask(argv[3]);
        if(mask>0) {
            event_mask = mask;
        } else {
            std::cerr << "Unknown mask " << argv[3] << std::endl;
            exit(1);            
        }
    }

    std::string base_path = "input/";

    std::ostringstream filename;
    filename << "Analysis_" << run_number;
    if(event_mask>0) {
        const char *mask = TPOEvent::DecodeEventMask(event_mask);
        filename << "_" << mask;
    }
    filename << ".root";

    // Print the filename to verify
    std::cout << "Writing TPORecEvent and histograms into file: " << filename.str() << " ..... " << std::endl;

    TFile *m_rootFile = new TFile(filename.str().c_str(), "RECREATE", "", 0); // last is the compression level
    if (!m_rootFile || !m_rootFile->IsOpen())
    {
        throw std::runtime_error("Could not create ROOT file");
    }
    m_rootFile->cd();
    TTree *t = new TTree("Event", "Event");
    t->Branch("t_reaction",&event.t_reaction);
    t->Branch("t_Eneutrino",&event.t_Eneutrino);
    t->Branch("n_clusters",&event.n_clusters);
    t->Branch("c_E1", &event.c_E1);
    t->Branch("c_chi2_1", &event.c_chi2_1);
    t->Branch("c_a_1", &event.c_a_1);
    t->Branch("c_b_1", &event.c_b_1);

    t->Branch("c_E2", &event.c_E2);
    t->Branch("rear_Cal", &event.rear_Cal);
    t->Branch("rear_MuCal", &event.rear_MuCal);

    // charm
    TH1D h_charm_Enuall = TH1D("h_charm_Enuall", "Neutrino energy", 50, 0, 4000.0);
    TH1D h_charm_Enucharmed = TH1D("h_charm_Enucharmed", "Neutrino energy", 50, 0, 4000.0);

    std::ostringstream inputfilename;
    inputfilename << "input/Batch-TPORecevent_" << run_number << "_*_*.root";

#if 0
    TFile *m_inrootFile = new TFile(inputfilename.str().c_str(), "READ"); 
    if (!m_inrootFile || !m_inrootFile->IsOpen())
    {
        return 1;
    }
    m_inrootFile->cd();    
    TTree *event_tree;
    m_inrootFile->GetObject("RecoEvent",event_tree);
#endif
    TChain *event_tree = new TChain("RecoEvent","READ");
    event_tree->Add(inputfilename.str().c_str());

    Long_t nentries = event_tree->GetEntries();
    std::cout << "Number of entries " << nentries << std::endl;

    // TPORecoEvent class and histograms
    TPORecoEvent *fTPORecoEvent = nullptr; // &fTPORecoEvent;
    event_tree -> SetBranchAddress("TPORecoEvent", &fTPORecoEvent);

    // process events
    int ievent = 0;
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
 
        // fill event ntuple
        event.t_Eneutrino = fTPORecoEvent -> GetPOEvent()->in_neutrino.m_energy;
        if(fTPORecoEvent -> GetPOEvent()->isES()) {
            event.t_reaction = 20;
        } else {
            event.t_reaction = fTPORecoEvent -> GetPOEvent()->isCC ? 0 : 10;
            int m_pdg = abs(fTPORecoEvent -> GetPOEvent()->in_neutrino.m_pdg_id);
            switch(m_pdg){
                case 12:
                    event.t_reaction += 1;
                    break;
                case 14:
                    event.t_reaction += 2;
                    break;
                case 16:
                    event.t_reaction += 3;
                    break;
                default:
                break;
            }
        }
        event.n_clusters = fTPORecoEvent->n_psclustersX();
        if(event.n_clusters>0) {
            event.c_E1 = fTPORecoEvent->PSClustersX[0].rawenergy/1e3;   // convert to GeV
            event.c_chi2_1 = fTPORecoEvent->PSClustersX[0].longenergyprofile.chi2_per_ndf;
            event.c_a_1 = fTPORecoEvent->PSClustersX[0].longenergyprofile.a;
            event.c_b_1 = fTPORecoEvent->PSClustersX[0].longenergyprofile.b;
        } else {
            event.c_E1 = 0;
            event.c_chi2_1 = -999;
            event.c_a_1 = event.c_b_1 = -999;
        }
        if(event.n_clusters>1) {
            event.c_E2 = fTPORecoEvent->PSClustersX[1].rawenergy/1e3;   // convert to GeV
        } else {
            event.c_E2 = 0;
        }

        // store energies in rear calorimeters
        event.rear_Cal = fTPORecoEvent->rearCals.rearCalDeposit;
        event.rear_MuCal = fTPORecoEvent->rearCals.rearMuCalDeposit;
    
        // charm 
        double enu = fTPORecoEvent -> GetPOEvent() -> POs[0].m_energy;
        h_charm_Enuall.Fill(enu);
        if(fTPORecoEvent -> GetPOEvent() -> isCharmed()) {
            fTPORecoEvent -> GetPOEvent() -> dump_event();
            h_charm_Enucharmed.Fill(enu);
        }

        t->Fill();
        ievent++;
    }

    m_rootFile->Write();
    m_rootFile->Close();

    myPOevent.dump_stats();

    std::cout << "I'm done." << std::endl;

    return 0;
}
