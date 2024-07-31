// Convert FASER official MC events into POevents
//
// A. Rubbia, July 2024
//

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <map>
#include <regex>

#include <TTree.h>
#include <TChain.h>
#include <TFile.h>

#include "TPOEvent.hh"

std::map<int, long> index_events;

int load_map_XML(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile) {
        std::cerr << "Error opening file for reading: " << filename << std::endl;
        return 1;
    }

    std::cout << "Loading map..." << std::endl;

    std::string line;
    std::string key;
    int event;
    long index;
    std::regex keyRegex("<key>([^<]+)</key>");
    std::regex valueRegex("<value>([^<]+)</value>");

    while (std::getline(inFile, line)) {
        std::smatch match;
        if (std::regex_search(line, match, keyRegex)) {
            key = match[1];
        } else if (std::regex_search(line, match, valueRegex)) {
            event = std::stol(key);
            index = std::stol(match[1]);
            index_events[event] = index;
        }
    }

    inFile.close();
    return 0;
}

int load_map_TXT(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile) {
        std::cerr << "Error opening file for reading: " << filename << std::endl;
        return 1;
    }
    std::cout << "Loading map..." << std::endl;

    int event;
    long index;
    while (inFile >> event >> index) {
        index_events[event] = index;
    }
    inFile.close();
    return 0;
}


void save_map_TXT(const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    for (const auto& pair : index_events) {
        outFile << pair.first << " " << pair.second << std::endl;
    }

    outFile.close();
}

void create_map(int run_number, std::string inputDirFiles) {
    std::cout << "Loading FASERMC files " << inputDirFiles << " ..." << std::endl;

    TChain *tree = new TChain("m_NuMCTruth_tree");
    tree->Add(inputDirFiles.c_str());  

    int m_event_id_MC;
    tree->SetBranchAddress("m_event_id_MC", &m_event_id_MC);

    long nentries = tree->GetEntries();
    std::cout << "Number of entries ... " << nentries << std::endl;

    long tree_ientry = 0;
    int last_event_id_MC = -1;

    while (tree_ientry < nentries) {
        tree->GetEntry(tree_ientry);
        if(last_event_id_MC != m_event_id_MC) {
    	    last_event_id_MC = m_event_id_MC;
            index_events[m_event_id_MC] = tree_ientry;
            if(m_event_id_MC % 1000 == 0){
                std::cout << "Event " << m_event_id_MC << " ... index " << tree_ientry << std::endl;
            }
	    }
        tree_ientry++;
    }
    delete tree;
}

void convert_FASERMC(int run_number, std::string inputDirFiles, int min_event, int max_event,
    std::string ROOTOutputFile, int mask) {

    std::cout << "Converting events ..." << std::endl;

    TFile *m_rootFile = new TFile(ROOTOutputFile.c_str(), "RECREATE", "", 505); // last is the compression level
    if (!m_rootFile || !m_rootFile->IsOpen())
    {
        throw std::runtime_error("Could not create output ROOT file");
    }
    m_rootFile->cd();

    TPOEvent fTPOEvent;
    TPOEvent *branch_POEvent = &fTPOEvent;
    TTree *m_POEventTree = new TTree("POEvent", "POEvent");
    m_POEventTree->Branch("event", &branch_POEvent);    // this should be named POEvent

    TChain *tree = new TChain("m_NuMCTruth_tree");
    tree->Add(inputDirFiles.c_str());  

    // FASER MC ntuple
    // Set up variables to hold the data and link them to the tree branches
    Int_t m_runnumber, m_event_id_MC, m_track_id, m_pdg_id, m_num_in_particle, m_num_out_particle;
    Double_t m_px, m_py, m_pz, m_energy, m_kinetic_energy, m_mass;
    Float_t m_vx_prod, m_vy_prod, m_vz_prod, m_vx_decay, m_vy_decay, m_vz_decay;
    std::vector<int> *m_pdg_in_particle = nullptr, *m_pdg_out_particle = nullptr;
    std::vector<int> *m_trackid_in_particle = nullptr, *m_trackid_out_particle = nullptr;
    Int_t m_status;

    // all energies are in MeV
	tree->SetBranchAddress("m_runnumber", &m_runnumber);
	tree->SetBranchAddress("m_event_id_MC", &m_event_id_MC);
	tree->SetBranchAddress("m_track_id", &m_track_id);
	tree->SetBranchAddress("m_pdg_id", &m_pdg_id);
	tree->SetBranchAddress("m_px", &m_px);
	tree->SetBranchAddress("m_py", &m_py);
	tree->SetBranchAddress("m_pz", &m_pz);
	tree->SetBranchAddress("m_energy", &m_energy);
	tree->SetBranchAddress("m_kinetic_energy", &m_kinetic_energy);
	tree->SetBranchAddress("m_mass", &m_mass);
	tree->SetBranchAddress("m_vx_prod", &m_vx_prod);
	tree->SetBranchAddress("m_vy_prod", &m_vy_prod);
	tree->SetBranchAddress("m_vz_prod", &m_vz_prod);
	tree->SetBranchAddress("m_vx_decay", &m_vx_decay);
	tree->SetBranchAddress("m_vy_decay", &m_vy_decay);
	tree->SetBranchAddress("m_vz_decay", &m_vz_decay);
	tree->SetBranchAddress("m_pdg_in_particle", &m_pdg_in_particle);
	tree->SetBranchAddress("m_pdg_out_particle", &m_pdg_out_particle);
	tree->SetBranchAddress("m_trackid_in_particle", &m_trackid_in_particle);
	tree->SetBranchAddress("m_trackid_out_particle", &m_trackid_out_particle);
	tree->SetBranchAddress("m_status", &m_status);

    int evt_to_dump = 0;

    size_t n_entries = tree->GetEntries();

    for (auto it : index_events) {
        int event = it.first;
        if( event < min_event || event > max_event) continue;

        long tree_ientry = it.second;

        if(event % 1000 == 0) {
            std::cout << "Processing event " << event << " ..." << std::endl;
        }

        fTPOEvent.clear_event();
    
        bool must_end = false;
        bool got_primvtx = false;
        bool found_tau_lepton = false;
        int tau_lepton_track_id = 0;
        while (!must_end){
            tree->GetEntry(tree_ientry);
            if(event != m_event_id_MC || tree_ientry >= n_entries) {
                if(fTPOEvent.n_particles() == 0){
                    std::cerr << " Something went wrong .. maybe wrong index? -- empty event" << std::endl;
                }
                must_end = true;
                continue;
		    }
            fTPOEvent.run_number = m_runnumber;
	        fTPOEvent.event_id = m_event_id_MC;

    	    if(!got_primvtx && m_status == 1) {
    	        fTPOEvent.setPrimaryVtx(m_vx_prod, m_vy_prod, m_vz_prod);  // in mm
	            got_primvtx = true;
	        }
	    
    	    struct PO aPO;
	        aPO.m_pdg_id = m_pdg_id;
	        aPO.m_track_id = m_track_id;
	        aPO.m_px = m_px/1e3;
	        aPO.m_py = m_py/1e3;
	        aPO.m_pz = m_pz/1e3;
		    aPO.m_energy = m_energy/1e3;
		    aPO.m_kinetic_energy = m_kinetic_energy/1e3;	
	        aPO.m_vx_decay = m_vx_prod-fTPOEvent.prim_vx.X();
	        aPO.m_vy_decay = m_vy_prod-fTPOEvent.prim_vx.Y();
	        aPO.m_vz_decay = m_vz_prod-fTPOEvent.prim_vx.Z();
	        aPO.nparent = m_trackid_in_particle->size();
	        for (int i=0; i<aPO.nparent;i++){
	            aPO.m_trackid_in_particle[i] = m_trackid_in_particle->at(i);
    	    };
    	    aPO.m_status = m_status;
    		aPO.geanttrackID = -1;
	    
    	    if(m_track_id < 20000 && m_status != 3) {
    	        fTPOEvent.POs.push_back(aPO);
	        }
	    
	        // found charged tau lepton - store decay products
	        if(!found_tau_lepton && abs(aPO.m_pdg_id) == 15) {
	            found_tau_lepton = true;
	            tau_lepton_track_id = m_track_id;   // obsolete
#ifdef _INCLUDE_PYTHIA_
                fTPOEvent.perform_taulepton_decay(aPO);
#endif
	        }

#if 0	    
	        if(found_tau_lepton) {
	            for(int i=0;i<aPO.nparent;i++) {
				    if(aPO.m_trackid_in_particle[i] == tau_lepton_track_id) {
				        fTPOEvent.taudecay.push_back(aPO);
				        double decaylength = sqrt(aPO.m_vx_decay*aPO.m_vx_decay+aPO.m_vy_decay*aPO.m_vy_decay+aPO.m_vz_decay*aPO.m_vz_decay);
				        fTPOEvent.tautracklength = decaylength;
				    }
                }
	        }
#endif

            tree_ientry++;

        } // while
	
    	fTPOEvent.kinematics_event();

    	if(fTPOEvent.istau && fTPOEvent.isCC && fTPOEvent.n_taudecay()==0) {
    	    std::cerr << "Convert_FASERMC: Could not find tau decay product??" << std::endl;
            fTPOEvent.dump_event();
    	}

        // now check for mask
        bool masked = false;
        if(mask>0) {
            bool wanted = false;
            wanted |= (mask == TPOEvent::kMask_nueCC && fTPOEvent.isCC && 
                abs(fTPOEvent.in_neutrino.m_pdg_id) == 12);
            wanted |= (mask == TPOEvent::kMask_numuCC && fTPOEvent.isCC && 
                abs(fTPOEvent.in_neutrino.m_pdg_id) == 14);
            wanted |= (mask == TPOEvent::kMask_nutauCC && fTPOEvent.isCC && 
                abs(fTPOEvent.in_neutrino.m_pdg_id) == 16);
            wanted |= (mask == TPOEvent::kMask_NC && !fTPOEvent.isCC);
            wanted |= (mask == TPOEvent::kMask_ES && fTPOEvent.isES());
            masked = !wanted;
            fTPOEvent.SetEventMask(mask);
        }
        if(!masked) {
            if(evt_to_dump++ < 20) {
        	    fTPOEvent.dump_header();
        	    fTPOEvent.dump_event();	
            };
            m_POEventTree -> Fill();
        }
    } // for

    std::cout << " Processed " << evt_to_dump << " events." << std::endl;
    m_POEventTree->Write();
    m_rootFile->Close();
    std::cout << "Done saving..." << std::endl;
}

bool confirmAction() {
    char response;
    while (true) {
        std::cout << "Do you want to proceed to regenerate the map (very slow)? (y/n): ";
        std::cin >> response;

        // Convert response to lowercase to make the check case-insensitive
        response = std::tolower(response);

        if (response == 'y') {
            return true;
        } else if (response == 'n') {
            return false;
        } else {
            std::cout << "Invalid input. Please enter 'y' for yes or 'n' for no." << std::endl;
        }
    }
}

int main(int argc, char** argv) {
    // get the output file name as the first argument
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " <run> [minevent] [maxevent] [mask]" << std::endl;
        std::cout << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "   minevent                  Minimum event to process (def=0)" << std::endl;
        std::cout << "   maxevent                  Maximum event to process (def=1000):";
        std::cout << "  =-1 for all events" << std::endl;
        std::cout << "   mask                      To process only specific events (def=none): ";
        std::cout << "  nueCC, numuCC, nutauCC, nuNC or nuES" << std::endl;
		return 1;
	}

    std::string runString = argv[1];

    int run_number;
    int min_event = 0;
    int max_event = 1000;

    try {
        run_number = std::stoi(runString);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument for run: " << e.what() << std::endl;
    } catch (const std::out_of_range& e) {
        std::cerr << "Out of range for run: " << e.what() << std::endl;
    }

    if(argc>2){
        try {
            min_event = std::stoi(argv[2]);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument for minevent: " << e.what() << std::endl;
        } catch (const std::out_of_range& e) {
            std::cerr << "Out of range for minevent: " << e.what() << std::endl;
        }
    }

    if(argc>3){
        try {
            max_event = std::stoi(argv[3]);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument for maxevent: " << e.what() << std::endl;
        } catch (const std::out_of_range& e) {
            std::cerr << "Out of range for maxevent: " << e.what() << std::endl;
        }
    }
    if(max_event == -1) {
        max_event = 99999999;
    }

    int event_mask = 0;
    if(argc>4) {
        int mask = TPOEvent::EncodeEventMask(argv[4]);
        if(mask>0) {
            event_mask = mask;
        } else {
            std::cerr << "Unknown mask " << argv[4] << std::endl;
            exit(1);            
        }
    }

    std::ostringstream inputDirFiles;
    inputDirFiles << "../FASERDATA/sim" << run_number << "/*.root";
    std::ostringstream ROOTOutputFile;
    ROOTOutputFile << "FASERMC-PO-Run" << run_number << "-" << min_event << "_" << max_event;
    if(event_mask>0) {
        const char *mask = TPOEvent::DecodeEventMask(event_mask);
        ROOTOutputFile << "_" << mask;
    }
    ROOTOutputFile << ".root";
    std::ostringstream MapFile;
    MapFile << "map_run" << run_number << ".txt";

    std::cout << "Converting FASERMC from " << inputDirFiles.str() << std::endl;
    std::cout << "The event indices map is " << MapFile.str() << std::endl;
    std::cout << "The output file is " << ROOTOutputFile.str() << std::endl;

    if(load_map_TXT(MapFile.str()) != 0) {
        if (confirmAction()) {
            create_map(run_number, inputDirFiles.str());
            save_map_TXT(MapFile.str());
        }
    }

    convert_FASERMC(run_number, inputDirFiles.str(), min_event, max_event, ROOTOutputFile.str(), event_mask);
    
    std::cout << "I'm done." << std::endl;
    return 0;
}
