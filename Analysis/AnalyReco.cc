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

#define EVENT_MAX_VERTICES 100
struct EVENT {
    // index
    Int_t event_number;
    // Truth variables
    Int_t t_reaction; // reaction 1=nue, 2=numu, 3=nutau, +10 if NC, 20=ES (truth)
    Float_t t_Eneutrino; // true incoming neutrino energy
    Float_t t_outlepton; // true outgoing lepton energy
    Float_t t_jetE; // true jet energy
    Float_t t_jetPt; // true jet transverse momentum
    Float_t t_Evis; // true visible energy
    Float_t t_ptmiss; // true missing transverse momentum
    Float_t t_primvx, t_primvy, t_primvz; // true primary vertex
    Int_t t_ischarmed; // true if charmed event
    Int_t t_tau_decaymode; // tau decay mode
    Float_t t_tautracklength; // tau track length
    Float_t t_tauKinAngle; // tau kink angle
    // reconstructed quantities
    Int_t n_vertices;
    Float_t v_x[EVENT_MAX_VERTICES], v_y[EVENT_MAX_VERTICES], v_z[EVENT_MAX_VERTICES];
    Int_t v_ntrks[EVENT_MAX_VERTICES]; // number of tracks in vertex
    Int_t t_closest_vertex;   // index of vertex closest to true primary vertex
    Int_t n_tktracks;
    Int_t n_pstracks;
    Int_t n_clusters;
    Float_t c_E1;    // energy most energetic cluster
    Float_t c_E1T;     // transverse energy most energetic cluster
    Float_t c_chi2_1;
    Float_t c_a_1;
    Float_t c_b_1;
    Float_t c_E2;    // energy 2nd most energetic cluster
    Float_t eflow_E;
    Float_t eflow_Et;
    Float_t faserCal; // total energy deposit in FaserCal
    Float_t faserCal_frac[15]; // fraction of energy in each FaserCal module
    Float_t rear_Cal; // energy deposit in rear Cal
    Float_t rear_HCal; // total energy deposit in rear hCal
    Float_t rear_HCal_frac[9]; // fraction of energy in each rear HCal module
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
    t->Branch("event_number",&event.event_number);
    t->Branch("t_reaction",&event.t_reaction);
    t->Branch("t_Eneutrino",&event.t_Eneutrino);
    t->Branch("t_outlepton",&event.t_outlepton);
    t->Branch("t_jetE",&event.t_jetE);
    t->Branch("t_jetPt",&event.t_jetPt);
    t->Branch("t_Evis",&event.t_Evis);
    t->Branch("t_ptmiss",&event.t_ptmiss);
    t->Branch("t_primvx",&event.t_primvx);
    t->Branch("t_primvy",&event.t_primvy);
    t->Branch("t_primvz",&event.t_primvz);
    t->Branch("t_ischarmed",&event.t_ischarmed);
    t->Branch("t_tau_decaymode",&event.t_tau_decaymode);
    t->Branch("t_tautracklength",&event.t_tautracklength);
    t->Branch("t_tauKinAngle",&event.t_tauKinAngle);
    t->Branch("n_vertices",&event.n_vertices);
    t->Branch("v_x", &event.v_x, "v_x[n_vertices]/F");
    t->Branch("v_y", &event.v_y, "v_y[n_vertices]/F");
    t->Branch("v_z", &event.v_z, "v_z[n_vertices]/F");
    t->Branch("v_ntrks", &event.v_ntrks, "v_ntrks[n_vertices]/I");
    t->Branch("t_closest_vertex",&event.t_closest_vertex);
    t->Branch("n_tktracks",&event.n_tktracks);
    t->Branch("n_pstracks",&event.n_pstracks);
    t->Branch("n_clusters",&event.n_clusters);
    t->Branch("c_E1", &event.c_E1);
    t->Branch("c_E1T", &event.c_E1T);
    t->Branch("c_chi2_1", &event.c_chi2_1);
    t->Branch("c_a_1", &event.c_a_1);
    t->Branch("c_b_1", &event.c_b_1);
    t->Branch("c_E2", &event.c_E2);
    t->Branch("eflow_E", &event.eflow_E);
    t->Branch("eflow_Et", &event.eflow_Et);
    t->Branch("faserCal", &event.faserCal);
    t->Branch("faserCal_frac", &event.faserCal_frac, "faserCal_frac[15]/F");
    t->Branch("rear_Cal", &event.rear_Cal);
    t->Branch("rear_HCal", &event.rear_HCal);
    t->Branch("rear_HCal_frac", &event.rear_HCal_frac, "rear_HCal_frac[9]/F");
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
        event.event_number = fTPORecoEvent -> GetPOEvent()-> event_id;
        // t_reaction = 1=nue, 2=numu, 3=nutau, +10 if NC, 20=ES (truth)
        event.t_Eneutrino = fTPORecoEvent -> GetPOEvent()->in_neutrino.m_energy;
        event.t_outlepton = fTPORecoEvent -> GetPOEvent()->out_lepton.m_energy;
        TVector3 jet = TVector3(fTPORecoEvent -> GetPOEvent()->jetpx, fTPORecoEvent -> GetPOEvent()->jetpy, fTPORecoEvent -> GetPOEvent()->jetpz);
        event.t_jetE = jet.Mag();
        event.t_jetPt = jet.Pt();
        event.t_Evis = fTPORecoEvent -> GetPOEvent()->Evis;
        event.t_ptmiss = fTPORecoEvent -> GetPOEvent()->ptmiss;
        if(fTPORecoEvent -> GetPOEvent()->isES()) {
            event.t_reaction = 20;
        } else {
            event.t_reaction = fTPORecoEvent -> GetPOEvent()->isCC ? 0 : 10;
            int m_pdg = fTPORecoEvent -> GetPOEvent()->in_neutrino.m_pdg_id;
            switch(abs(m_pdg)){
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
            if(m_pdg<0) event.t_reaction = -event.t_reaction;
        }
        event.t_ischarmed = fTPORecoEvent -> GetPOEvent()->isCharmed();
        event.t_tau_decaymode = fTPORecoEvent -> GetPOEvent()->tau_decaymode;
        event.t_tautracklength = fTPORecoEvent -> GetPOEvent()->tauDecaylength();
        event.t_tauKinAngle = fTPORecoEvent -> GetPOEvent()->tauKinkAngle();
        event.t_primvx = fTPORecoEvent -> GetPOEvent()->prim_vx.x();
        event.t_primvy = fTPORecoEvent -> GetPOEvent()->prim_vx.y();
        event.t_primvz = fTPORecoEvent -> GetPOEvent()->prim_vx.z();

        // store vertices
        if(fTPORecoEvent->fTKVertices.size()>EVENT_MAX_VERTICES) {
            std::cerr << "Number of vertices " << fTPORecoEvent->fTKVertices.size() << " exceeds maximum " << EVENT_MAX_VERTICES << std::endl;
            exit(1);
        }
        event.n_vertices = fTPORecoEvent->fTKVertices.size();
        for(int i=0; i<std::min(EVENT_MAX_VERTICES, event.n_vertices); i++) {
            event.v_x[i] = fTPORecoEvent->fTKVertices[i].position.x();
            event.v_y[i] = fTPORecoEvent->fTKVertices[i].position.y();
            event.v_z[i] = fTPORecoEvent->fTKVertices[i].position.z();
            event.v_ntrks[i] = fTPORecoEvent->fTKVertices[i].trackIDs.size();
        }

        // find vertex closest to true primary vertex
        double min_dist = 1e9;
        int min_index = -1;
        ROOT::Math::XYZVector true_vtx(fTPORecoEvent->GetPOEvent()->prim_vx.x(), fTPORecoEvent->GetPOEvent()->prim_vx.y(), fTPORecoEvent->GetPOEvent()->prim_vx.z());
        for(int i=0; i<event.n_vertices; i++) {
            double dist = (fTPORecoEvent->fTKVertices[i].position - true_vtx).R();
            if(dist<min_dist) {
                min_dist = dist;
                min_index = i;
            }
        }
        event.t_closest_vertex = min_index;

        // store number of tracks
        event.n_tktracks = fTPORecoEvent->fTKTracks.size();
        // store number of PS tracks
        event.n_pstracks = fTPORecoEvent->fTKTracks.size();

        // store number of PS clusters
        event.n_clusters = fTPORecoEvent->n_psclustersX();
        if(event.n_clusters>0) {
            event.c_E1 = fTPORecoEvent->PSClustersX[0].rawenergy/1e3;   // convert to GeV
            ROOT::Math::XYZVector dir = fTPORecoEvent->PSClustersX[0].cog-fTPORecoEvent->PSClustersX[0].vtx;
            event.c_E1T = fTPORecoEvent->PSClustersX[0].rawenergy/1e3*sqrt(dir.Unit().Perp2());
            event.c_chi2_1 = fTPORecoEvent->PSClustersX[0].longenergyprofile.chi2_per_ndf;
            if (std::isnan(event.c_chi2_1)) {
                event.c_chi2_1 = -1.0;
            }
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

        // store eflow
        event.eflow_E = fTPORecoEvent->GetPOFullRecoEvent()->TotalEvis();
        event.eflow_Et = fTPORecoEvent->GetPOFullRecoEvent()->TotalET();

        // store energies in FaserCal
        double totalFaserCal = 1e-9; // to avoid division by zero
        int istart = 0;
        for(int i=0; i<fTPORecoEvent->faserCals.size(); i++) {
            double edep = fTPORecoEvent->faserCals[i].EDeposit;
            totalFaserCal += edep;
            if(edep>10.0 && istart == 0) istart = i;
            event.faserCal_frac[i] = 0;
        }
        event.faserCal = totalFaserCal;
        for(int i=0; i<fTPORecoEvent->faserCals.size(); i++) {
            int moduleID = fTPORecoEvent->faserCals[i].ModuleID;
            if(moduleID>=istart)
                event.faserCal_frac[moduleID-istart] = fTPORecoEvent->faserCals[i].EDeposit/totalFaserCal;
        }
        // store energies in rear calorimeters
        event.rear_Cal = fTPORecoEvent->rearCals.rearCalDeposit;
        event.rear_HCal = fTPORecoEvent->rearCals.rearHCalDeposit;
        for(int i=0; i<fTPORecoEvent->rearCals.rearHCalModule.size(); i++) {
            int moduleID = fTPORecoEvent->rearCals.rearHCalModule[i].moduleID;
            event.rear_HCal_frac[moduleID] = fTPORecoEvent->rearCals.rearHCalModule[i].energyDeposit/(event.rear_HCal+1e-9);
        }
        event.rear_MuCal = fTPORecoEvent->rearCals.rearMuCalDeposit;
    
        // charm 
        double enu = fTPORecoEvent -> GetPOEvent() -> POs[0].m_energy;
        h_charm_Enuall.Fill(enu);
        if(fTPORecoEvent -> GetPOEvent() -> isCharmed()) {
//            fTPORecoEvent -> GetPOEvent() -> dump_event();
            h_charm_Enucharmed.Fill(enu);
        }

        t->Fill();
        ievent++;
    }

    m_rootFile->Write();
    m_rootFile->Close();

    myPOevent.dump_stats();

    std::cout << "Number of events processed: " << ievent << std::endl;

    std::cout << "I'm done." << std::endl;

    return 0;
}
