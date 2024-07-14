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
#include "TTauSearch.hh"

void load_geometry() {
    // Load the GDML geometry
    TGeoManager::Import("../GeomGDML/geometry.gdml");
}

int main(int argc, char** argv) {

	if (argc < 2) {
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
    std::cout << "Writing TPORecEvent and histograms into file: " << filename.str() << " ..... ";

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

    // dE/dx
    TH1D pi_dedx = TH1D("pi_dedx", "pions: dE in scintillator voxel (MeV)", 200, 0., 5.);
    TH1D p_dedx = TH1D("p_dedx", "protons: dE in scintillator voxel (MeV)", 200, 0., 5.);
    TH1D mu_dedx = TH1D("mu_dedx", "muons: dE in scintillator voxel (MeV)", 200, 0., 5.);
    TH1D tau_dedx = TH1D("tau_dedx", "taus: dE in scintillator voxel (MeV)", 200, 0., 5.);

    // hit level
    TH1D hit_scint_nhits = TH1D("hit_scint_nhits", "Scintillator nhits", 100, 0.,5e5);
    TH1D hit_scint_edepo = TH1D("hit_scint_edepo", "Scintillator edepo MeV", 100, 0.,10.);
    TH1D hit_tracker_nhits = TH1D("hit_tracker_nhits", "Tracker nhits", 100, 0.,100000.);
    TH1D hit_tracker_edepo = TH1D("hit_tracker_edepo", "Tracker edepo MeV", 100, 0.,10.);

    // TPORecoEvent class and histograms
    TPORecoEvent *branch_TPORecoEvent = nullptr; // &fTPORecoEvent;
    TTree *m_POEventTree = nullptr;
    TH1D h_fullevent_Evis = TH1D("h_fullevent_Evis", "Full event energy", 200, 0., 2000.);
    TH1D h_fullevent_ET = TH1D("h_fullevent_ET", "Full event transverse energy", 200, 0., 40.);

    // TauSearches
    TTauSearch fTTauSearch_e;
    TTree *m_tausearch_e_Tree = new TTree("Tausearch_e", "Tausearch_e");
    fTTauSearch_e.Create_Sel_Tree(m_tausearch_e_Tree);

    // process events
    int ievent = 0;
    if(max_event == -1) max_event = 99999999;
    int error = 0;

    while (error == 0 && ievent<max_event) {

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
            m_rootFile->cd();
            m_POEventTree = new TTree("RecoEvent", "RecoEvent");
            m_POEventTree->Branch("TPORecoEvent", &branch_TPORecoEvent);
        }

        // fill hit level histograms (time consuming!)
        std::map<long, double> hitmap_scint;
        std::map<long, double> hitmap_tracker;
        for (auto it : fPORecoEvent->GetPORecs())
        {
            int POID = it->POID;
            struct PO *aPO = &fTcalEvent->fTPOEvent->POs[POID];
            size_t ntracks = it->fGEANTTrackIDs.size();
            for (size_t i = 0; i < ntracks; i++)
            {
                DigitizedTrack *dt = it->DTs[i];
                size_t nhits = dt->fEnergyDeposits.size();
                for (size_t j = 0; j < nhits; j++)
                {
                    long chanID = dt->fhitIDs[j];
                    int hittype = fTcalEvent->getChannelTypefromID(chanID);
                    double energydeposit = dt->fEnergyDeposits[j];
                    if (hittype == 0){
                        auto it = hitmap_scint.find(chanID);
                        if (it != hitmap_scint.end())
                        {
                            it->second += energydeposit;
                        }
                        else
                        {
                            hitmap_scint[chanID] = energydeposit;
                        }
                    } else if (hittype == 1){
                        auto it = hitmap_tracker.find(chanID);
                        if (it != hitmap_tracker.end())
                        {
                            it->second += energydeposit;
                        }
                        else
                        {
                            hitmap_tracker[chanID] = energydeposit;
                        }
                    }
                }
            }
        }
        size_t nhits_scint = hitmap_scint.size();
        size_t nhits_tracker = hitmap_tracker.size();
        hit_scint_nhits.Fill(nhits_scint);
        hit_tracker_nhits.Fill(nhits_tracker);
        for (auto it : hitmap_scint) {
            hit_scint_edepo.Fill(it.second);
        }
        for (auto it : hitmap_tracker) {
            hit_tracker_edepo.Fill(it.second);
        }

        // fill reconstructed track histograms and some additional ntuple variables
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
            DigitizedTrack* dt = it->DTs[0];
            switch (abs(PDG)) {
            case 11:                    // electrons
                e_c_energy.Fill(f);
                e_c_energy2.Fill(POEne, f);
                break;
            case 13:                    // muons
                ntracks = it->fGEANTTrackIDs.size();
                for (size_t i = 0; i < ntracks; i++) {
                    DigitizedTrack* dt = it->DTs[i];
                    if(abs(dt->fPDG) != 13)continue;
                    nhits = dt->fEnergyDeposits.size();
                    for (size_t j = 0; j < nhits; j++) {
                        hittype = fTcalEvent -> getChannelTypefromID(dt->fhitIDs[j]);
                        if(hittype!=0)continue;
                        mu_dedx.Fill(dt->fEnergyDeposits[j]);
                    }
                }
                break;
            case 15:                            // taus
                // only primary track
                nhits = dt->fEnergyDeposits.size();
                for (size_t j = 0; j < nhits; j++) {
                    hittype = fTcalEvent -> getChannelTypefromID(dt->fhitIDs[j]);
                    if(hittype!=0)continue;
                    tau_dedx.Fill(dt->fEnergyDeposits[j]);
                }
                break;
            case 111:                           // pi0
            case 211:                           // pi+-
                if(POEne>2.0) {
                    pi_c_energy.Fill(f);
                    pi_c_energy2.Fill(POEne, f);
                }
                if(PDG == 111)continue;
                nhits = dt->fEnergyDeposits.size();
                for (size_t j = 0; j < nhits; j++) {
                    hittype = fTcalEvent -> getChannelTypefromID(dt->fhitIDs[j]);
                    if(hittype!=0)continue;
                    pi_dedx.Fill(dt->fEnergyDeposits[j]);
                }
                break;
            case 2112:                          // protons
                p_c_energy.Fill(f);
                p_c_energy2.Fill(POEne, f);
                nhits = dt->fEnergyDeposits.size();
                for (size_t j = 0; j < nhits; j++) {
                    hittype = fTcalEvent -> getChannelTypefromID(dt->fhitIDs[j]);
                    if(hittype!=0)continue;
                    p_dedx.Fill(dt->fEnergyDeposits[j]);
                }
                break;
            }
        };

        // full event histograms
        h_fullevent_Evis.Fill(fPORecoEvent->GetPOFullEvent()->TotalEvis());
        h_fullevent_ET.Fill(fPORecoEvent->GetPOFullEvent()->TotalET());

        // Tau Searches !!
        bool found_electron = false;
        bool found_tau_e = false;
        ROOT::Math::XYZVector spx = fPORecoEvent->GetPOFullEvent()->fTotal.Ecompensated*
                fPORecoEvent->GetPOFullEvent()->fTotal.Eflow.Unit();
        for (auto it : fPORecoEvent->GetPORecs())
        {
            int POID = it->POID;
            struct PO *aPO = &fTcalEvent->fTPOEvent->POs[POID];
            int PDG = aPO->m_pdg_id;

            double RecoLeptonEne = it->fTotal.Ecompensated;
            ROOT::Math::XYZVector lepton = RecoLeptonEne*it->fTotal.Eflow.Unit();
            ROOT::Math::XYZVector jet = spx-lepton;
            fTTauSearch_e.SetLepton(lepton);
            fTTauSearch_e.SetJet(jet);
            fTTauSearch_e.Kinematics(fPORecoEvent->primary_n_charged);

            // nueCC background
            if(abs(PDG) == 11 && !found_electron) {
                found_electron = true;
                fTTauSearch_e.Fill_Sel_Tree(m_tausearch_e_Tree);
            }

            // nutau signal -> e
            if(abs(PDG) == 15 && !found_tau_e) {
                if(fPORecoEvent->GetPOEvent()->tau_decaymode==1) {
                    found_tau_e = true;
                    fTTauSearch_e.Fill_Sel_Tree(m_tausearch_e_Tree);
                }
            }
        }

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
