// Batch reconstruction of FASERCAL events
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

void load_geometry() {
    // Load the GDML geometry
    TGeoManager::Import("../GeomGDML/geometry.gdml");
}

int main(int argc, char** argv) {

	if (argc < 2) {
	std::cout << "Usage: " << argv[0] << " <run> [mask]" << std::endl;
        std::cout << "   <run>                     Run number" << std::endl;
        std::cout << "   mask                      To process only specific events (def=none): ";
        std::cout << "  nueCC, numuCC, nutauCC, or nuNC" << std::endl;
    	return 1;
	}

   // get the run number as the first argument
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
    if(argc>2) {
        int mask = TPOEvent::EncodeEventMask(argv[2]);
        if(mask>0) {
            event_mask = mask;
        } else {
            std::cerr << "Unknown mask " << argv[2] << std::endl;
            exit(1);            
        }
    }

    load_geometry();

    std::string base_path = "input/";

    std::ostringstream filename;
    filename << "Batch-TPORecevent_" << run_number;
    if(event_mask>0) {
        const char *mask = TPOEvent::DecodeEventMask(event_mask);
        filename << "_" << mask;
    }
    filename << ".root";

    // Print the filename to verify
    std::cout << "Writing TPORecEvent into file: " << filename.str() << " ..... ";

    TFile *m_rootFile = new TFile(filename.str().c_str(), "RECREATE", "", 0); // last is the compression level
    if (!m_rootFile || !m_rootFile->IsOpen())
    {
        throw std::runtime_error("Could not create ROOT file");
    }
    m_rootFile->cd();

    // calorimetry histograms
    TH1D e_c_energy = TH1D("e_c_energy", "electrons: compensated energy fraction", 200, -1., 1.);
    TH2D e_c_energy2 = TH2D("e_c_energy2", "electrons: compensated energy fraction", 100, 0., 1000., 200, -1., 1.);
    TH1D pi_c_energy = TH1D("pi_c_energy", "pions: compensated energy fraction", 100, -1., 1.);
    TH2D pi_c_energy2 = TH2D("pi_c_energy2", "pions : compensated energy fraction vs E", 100, 0.,200.,100,-1.,1.);
    TH1D p_c_energy = TH1D("p_c_energy", "protons: compensated energy fraction", 100, -1., 1.);
    TH2D p_c_energy2 = TH2D("p_c_energy2", "protons : compensated energy fraction vs E", 100, 0.,200.,100,-1.,1.);

    // muons histograms
    TH1D mu_dedx = TH1D("mu_dedx", "muons: dE in scintillator voxel (MeV)", 200, 0., 5.);

    // TPORecoEvent class and histograms
    TPORecoEvent *branch_TPORecoEvent = nullptr; // &fTPORecoEvent;
    TTree *m_POEventTree = nullptr;
    TH1D h_fullevent_Evis = TH1D("h_fullevent_Evis", "Full event energy", 200, 0., 2000.);
    TH1D h_fullevent_ET = TH1D("h_fullevent_ET", "Full event transverse energy", 200, 0., 40.);

    // process events
    int ievent = 0;
    int error = 0;

    while (error == 0 && ievent<99999999) {

        if(ievent % 1000 == 0) {
            std::cout << "Processing event " << ievent << std::endl;
        }
        bool dump_event_cout = (ievent < 20);

        // Create an instance of TcalEvent and TPOEvent
        TcalEvent *fTcalEvent = new TcalEvent();
        TPOEvent *POevent = new TPOEvent();

        if(!dump_event_cout) fTcalEvent->SetVerbose(0);
        error = fTcalEvent -> Load_event(base_path, run_number, ievent++, event_mask, POevent);
        if(error != 0) break;
    
        if(dump_event_cout) {
            std::cout << "Transverse size " << fTcalEvent->geom_detector.fScintillatorSizeX << " mm " << std::endl;
            std::cout << "Total size of one sandwich layer " << fTcalEvent->geom_detector.fTotalLength << " mm " << std::endl;
            std::cout << "Number of layers " << fTcalEvent->geom_detector.NRep << std::endl;
            std::cout << "Voxel size " << fTcalEvent->geom_detector.fScintillatorVoxelSize << " mm " << std::endl;

            std::cout << " copied digitized tracks " << fTcalEvent->getfTracks().size() << std::endl;
        }

        if(dump_event_cout) fTcalEvent -> fTPOEvent -> dump_event();

        TPORecoEvent* fPORecoEvent = new TPORecoEvent(fTcalEvent, fTcalEvent->fTPOEvent);
        fPORecoEvent -> Reconstruct();
        if(dump_event_cout) fPORecoEvent -> Dump();

        TPORecoEvent *branch_TPORecoEvent = fPORecoEvent;
        if(m_POEventTree == nullptr) {
            m_POEventTree = new TTree("RecoEvent", "RecoEvent");
            m_POEventTree->Branch("TPORecoEvent", &branch_TPORecoEvent);
        }

        for (auto it : fPORecoEvent->GetPORecs())
        {
            int POID = it->POID;
            struct PO *aPO = &fTcalEvent->fTPOEvent->POs[POID];
            int PDG = aPO->m_pdg_id;
            double POEne = aPO->m_energy;
            double RecoEne = it->fTotal.Ecompensated;
            double f = (RecoEne - POEne) / POEne;
            int ntracks;
            size_t nhits;
            int hittype;
            switch (abs(PDG)) {
            case 11:                    // electrons
                e_c_energy.Fill(f);
                e_c_energy2.Fill(POEne, f);
                break;
            case 13:                    // muons
                ntracks = it->fGEANTTrackIDs.size();
                for (size_t i = 0; i < ntracks; i++) {
                    DigitizedTrack* dt = it->DTs[i];
                    nhits = dt->fEnergyDeposits.size();
                    for (size_t j = 0; j < nhits; j++) {
                        hittype = fTcalEvent -> getChannelTypefromID(dt->fhitIDs[i]);
                        if(hittype!=0)continue;
                        mu_dedx.Fill(dt->fEnergyDeposits[j]);
                    }
                }
                break;
            case 111:                           // pi0
            case 211:                           // pi+-
                pi_c_energy.Fill(f);
                pi_c_energy2.Fill(POEne, f);
                break;
            case 2112:                          // protons
                p_c_energy.Fill(f);
                p_c_energy2.Fill(POEne, f);
                break;
            }
        };

        // full event histograms
        h_fullevent_Evis.Fill(fPORecoEvent->GetPOFullEvent()->TotalEvis());
        h_fullevent_ET.Fill(fPORecoEvent->GetPOFullEvent()->TotalET());

        m_POEventTree -> Fill();

        delete fPORecoEvent;
        delete POevent;
        delete fTcalEvent;       
    }

    m_rootFile->Write();
    m_rootFile->Close();

    std::cout << "I'm done." << std::endl;

    return 0;
}
