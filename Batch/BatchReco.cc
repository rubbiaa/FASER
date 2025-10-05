// Batch reconstruction of FASERCAL events
//
// A. Rubbia, July 2024
//

#include <string>
#include <csignal>
#include <atomic>
#include <iostream>
#include <fstream>
#include <sys/resource.h>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <TGeoManager.h>

#include "TcalEvent.hh"
#include "TPORecoEvent.hh"
#include "TTauSearch.hh"
#include "TParticleGun.hh"
#include "TMuonSpect.hh"

#include <chrono>

// Global atomic flag to indicate whether the program should continue running
std::atomic<bool> keepRunning(true);
static int n_sigint = 0;

// Signal handler function
void handleSignal(int signal) {
    if (signal == SIGINT || signal == SIGTERM) {
        if(++n_sigint<3){
            std::cout << "\nCaught Ctrl+C (SIGINT) or SIGTERM. Exiting cleanly...\n";
            keepRunning = false; // Update the flag to stop the program loop
        } else {
            std::cout << "\nCaught Ctrl+C (SIGINT) or SIGTERM. Exiting now...\n";
            exit(1);
        }
    }
}

void load_geometry() {
    // Load the GDML geometry
    TGeoManager::Import("../GeomGDML/geometry.gdml");
}

int main(int argc, char** argv) {

	if (argc < 2) {
	std::cout << "Usage: " << argv[0] << " [-mt] <run> [minevent] [maxevent] [mask]" << std::endl;
        std::cout << "   <run>                     Run number" << std::endl;
        std::cout << "   minevent                  Minimum number of events to process (def=-1)" << std::endl;
        std::cout << "   maxevent                  Maximum number of events to process (def=-1)" << std::endl;
        std::cout << "   mask                      To process only specific events (def=none): ";
        std::cout << "  nueCC, numuCC, nutauCC, nuNC or nuES" << std::endl;
		return 1;
	}

    int argv_index = 1;
    bool multiThread_option = false;
    if(argc>argv_index) {
        if(std::string(argv[argv_index]) == "-mt") {
            multiThread_option = true;
            argv_index++;
        }
    }

   // get the run number as the first argument
    std::string runString = argv[argv_index++];
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

    int min_event = 0;
    if(argc>argv_index) {
        try {
            min_event = std::stoi(argv[argv_index++]);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument for minevent: " << e.what() << std::endl;
            exit(1);
        } catch (const std::out_of_range& e) {
            std::cerr << "Out of range for minevent: " << e.what() << std::endl;
            exit(1);
        }
    }

    int max_event = -1;
    if(argc>argv_index) {
        try {
            max_event = std::stoi(argv[argv_index++]);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument for maxevent: " << e.what() << std::endl;
            exit(1);
        } catch (const std::out_of_range& e) {
            std::cerr << "Out of range for maxevent: " << e.what() << std::endl;
            exit(1);
        }
    }
    if(max_event == -1) max_event = 99999999;

    int event_mask = 0;
    if(argc>argv_index) {
        int mask = TPOEvent::EncodeEventMask(argv[argv_index]);
        if(mask>0) {
            event_mask = mask;
        } else {
            std::cerr << "Unknown mask " << argv[argv_index] << std::endl;
            exit(1);            
        }
    }

    load_geometry();

    std::string base_path = "input/";

    std::ostringstream filename;
    filename << "Batch-TPORecevent_" << run_number << "_" << min_event << "_" << max_event;
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
    TH1D pi_c_energy = TH1D("pi_c_energy", "charged pions: compensated energy fraction", 100, -1., 1.);
    TH2D pi_c_energy2 = TH2D("pi_c_energy2", "charged pions : compensated energy fraction vs E", 100, 0.,200.,100,-1.,1.);
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

    // charm
    TH1D h_charm_Enuall = TH1D("h_charm_Enuall", "Neutrino energy", 50, 0, 4000.0);
    TH1D h_charm_Enucharmed = TH1D("h_charm_Enucharmed", "Neutrino energy", 50, 0, 4000.0);
    
    // TauSearches
    TTauSearch fTTauSearch_e;
    TTree *m_tausearch_e_Tree = new TTree("Tausearch_e", "Tausearch_e");
    fTTauSearch_e.Create_Sel_Tree(m_tausearch_e_Tree);

    // TParticleGun
    TParticleGun fTParticleGun;
    TTree *m_particlegun_Tree = new TTree("ParticleGun","ParticleGun");
    fTParticleGun.Create_Sel_Tree(m_particlegun_Tree);

    // TMuonSpectrometer
    TMuonSpectrometer fTMuonSpectrometer;
    TTree *m_muonspect_Tree = new TTree("MuonSpectrometer","MuonSpectrometer");
    fTMuonSpectrometer.Create_Sel_Tree(m_muonspect_Tree);

    // process events
    int ievent = min_event;
    int error = 0;

    // Register the signal handler for SIGINT (Ctrl+C)
    std::signal(SIGINT, handleSignal);
    std::signal(SIGTERM, handleSignal);

    // Total elapsed time to process job
    long long total_time = 0;
    size_t n_events = 0;

    while (keepRunning && error == 0 && ievent<max_event) {

        auto start = std::chrono::high_resolution_clock::now();

        if(ievent % 1000 == 0) {
            std::cout << "Processing event " << ievent << std::endl;
        }
        bool dump_event_cout = (ievent < 20);

        // Create an instance of TcalEvent and TPOEvent
        TcalEvent *fTcalEvent = new TcalEvent();
        TPOEvent *POevent = new TPOEvent();

        if(!dump_event_cout) fTcalEvent->SetVerbose(0);
        error = fTcalEvent -> Load_event(base_path, run_number, ievent++, event_mask, POevent);
        if(error==2){ error = 0; continue; }
        if(error != 0) break;

        // skip empty events (can happen for example for particle guns G4 where no lepton or pion has been found (i.e. NC ES))
        if(POevent -> POs.size() == 0) continue;
    
        if(dump_event_cout) {
            std::cout << "Transverse size " << fTcalEvent->geom_detector.fScintillatorSizeX << " mm " << std::endl;
            std::cout << "Total size of one sandwich layer " << fTcalEvent->geom_detector.fTotalLength << " mm " << std::endl;
            std::cout << "Number of layers " << fTcalEvent->geom_detector.NRep << std::endl;
            std::cout << "Voxel size " << fTcalEvent->geom_detector.fScintillatorVoxelSize << " mm " << std::endl;

            std::cout << " copied digitized tracks " << fTcalEvent->getfTracks().size() << std::endl;
        }

        if(dump_event_cout) { POevent -> dump_event(); }
        else POevent->dump_header();

        // charm 
        double enu = POevent -> POs[0].m_energy;
        h_charm_Enuall.Fill(enu);
        if(POevent -> isCharmed()) {
            POevent -> dump_event();
            h_charm_Enucharmed.Fill(enu);
        }

        // tau decay length
        if(POevent -> istau) {
           POevent -> dump_event();
           struct PO aPO = POevent-> POs[4];    // first decay product
           ROOT::Math::XYZVector decayvtx = ROOT::Math::XYZVector(aPO.m_vx_decay, aPO.m_vy_decay, aPO.m_vz_decay);   
        }

#if 0
        //// 
        delete POevent;
        delete fTcalEvent;   
        continue;
        /////
#endif
        TPORecoEvent* fPORecoEvent = new TPORecoEvent(fTcalEvent, fTcalEvent->fTPOEvent);
        fPORecoEvent -> verbose = 0;   // 0: no output, 1: some output, 2: more output, 3: full output;
        if(!dump_event_cout) fPORecoEvent -> verbose = 0;
        fPORecoEvent -> multiThread = multiThread_option;
        fPORecoEvent -> ReconstructTruth();
        fPORecoEvent -> Reconstruct2DViewsPS();
        fPORecoEvent -> ReconstructClusters(0);
        fPORecoEvent -> Reconstruct3DPS_2();
        fPORecoEvent -> ReconstructRearCals();
        fPORecoEvent -> ReconstructMuonSpectrometer();
        fPORecoEvent -> Reconstruct3DPS_Eflow();
        fPORecoEvent -> TrackReconstruct();

        if(dump_event_cout) { fPORecoEvent -> Dump(); }

        fPORecoEvent -> Fill2DViewsPS();
    
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
            if(POID < 0) continue;
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
            if(POID < 0) continue;
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
            case 211:                           // pi+-
                if(POEne>2.0) {
                    pi_c_energy.Fill(f);
                    pi_c_energy2.Fill(POEne, f);
                }
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
            if(POID < 0) continue;
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
                if(!POevent->isES())                     // skip ES events
                    fTTauSearch_e.Fill_Sel_Tree(m_tausearch_e_Tree);
            }

            // nutau signal -> e
            if(abs(PDG) == 15 && !found_tau_e) {
                if(fPORecoEvent->GetPOEvent()->tau_decaymode==1) {
                    found_tau_e = true;
                    if(!POevent->isES())                     // skip ES events
                        fTTauSearch_e.Fill_Sel_Tree(m_tausearch_e_Tree);
                }
            }
        }

        // particle gun studies
        bool particle_gun = true;
        double emax_cluster = 0;
        if (particle_gun)
        {
            if ((fPORecoEvent->GetPORecs()).size() == 0)
                continue;
            struct TPORec *aPORec = (fPORecoEvent->GetPORecs())[0];
            int POID = aPORec->POID;
            if (POID >= 0)
            {
                struct PO *aPO = &fTcalEvent->fTPOEvent->POs[POID];
                fTParticleGun.features.m_pdg_id = aPO->m_pdg_id;
                fTParticleGun.features.m_energy = aPO->m_energy;

                // fill features of most energetic reconstructed cluster
                if (fPORecoEvent->PSClustersX.size() > 0)
                {
                    TPSCluster *c = &fPORecoEvent->PSClustersX[0]; // most energetic one
                    fTParticleGun.features.ep_chi2_per_ndf = c->longenergyprofile.chi2_per_ndf;
                    fTParticleGun.features.ep_E0 = c->longenergyprofile.E0;
                    fTParticleGun.features.ep_a = c->longenergyprofile.a;
                    fTParticleGun.features.ep_b = c->longenergyprofile.b;
                    fTParticleGun.features.ep_tmax = c->longenergyprofile.tmax;
                    fTParticleGun.features.ep_c = c->longenergyprofile.c;
                }

                fTParticleGun.Fill_Sel_Tree(m_particlegun_Tree);
            }
        }

        // TMuonSpectrometer
        fTMuonSpectrometer.features.ntracks = fPORecoEvent->fMuTracks.size();
        int ntracks = fPORecoEvent->fMuTracks.size();
        if(ntracks>MAXMUTRACKS) ntracks = MAXMUTRACKS;
        for(int i=0; i<ntracks; i++) {
            TMuTrack *aMuTrack = &fPORecoEvent->fMuTracks[i];
            fTMuonSpectrometer.features.charge[i] = aMuTrack->fcharge;
            fTMuonSpectrometer.features.npoints[i] = aMuTrack->fpos.size();
            fTMuonSpectrometer.features.px[i] = aMuTrack->fpx;
            fTMuonSpectrometer.features.py[i] = aMuTrack->fpy;
            fTMuonSpectrometer.features.pz[i] = aMuTrack->fpz;
            fTMuonSpectrometer.features.p[i] = aMuTrack->fp;
            fTMuonSpectrometer.features.chi2[i] = aMuTrack->fchi2;
            fTMuonSpectrometer.features.nDoF[i] = aMuTrack->fnDoF;
            fTMuonSpectrometer.features.pval[i] = aMuTrack->fpval;
        }
        fTMuonSpectrometer.Fill_Sel_Tree(m_muonspect_Tree);

        m_POEventTree -> Fill();


        struct rusage usage;
        double mem_usage = -1;
        if (getrusage(RUSAGE_SELF, &usage) == 0) {
            mem_usage = usage.ru_maxrss;
        } else {
            perror("getrusage failed");
        }

        delete fPORecoEvent;
        delete POevent;
        delete fTcalEvent;       

        // Record the end time
        auto end = std::chrono::high_resolution_clock::now();

        // Calculate the elapsed time in milliseconds
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        std::cerr << "+++++ Run " << run_number << " Event processed " << n_events << " - Elapsed time: " << elapsed.count() << " ms  ";
        std::cerr << "Memory usage: " << usage.ru_maxrss << " KB" << std::endl;

        total_time += elapsed.count();

        // chekc if environment variable is set
        const char *heartbeat_file;
        char *env = std::getenv("HEARTBEAT_FILE");
        if(env) {
            heartbeat_file = env;
        } else {
            heartbeat_file = "heartbeat.txt";
        }
        // write down the heartbeat file
        std::ofstream heartbeat(heartbeat_file);
        heartbeat << "Run " << run_number << " Event " << ievent << " processed - elasped time " << elapsed.count() << " ms - memory usage " << mem_usage << " KB" << std::endl;
        heartbeat.close();

        n_events++;

    }

    m_rootFile->Write();
    m_rootFile->Close();

    std::cerr << "Total elapsed time: " << total_time << " ms\n";
    double average_time = static_cast<double>(total_time) / n_events;

    std::cerr << "Average time per event: " << average_time << " ms\n";
    std::cout << "I'm done." << std::endl;

    return 0;
}
