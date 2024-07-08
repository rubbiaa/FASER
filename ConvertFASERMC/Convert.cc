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

#include "TPOEvent.hh"

std::map<int, long> index_events;

int load_map_XML(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile) {
        std::cerr << "Error opening file for reading: " << filename << std::endl;
        return 1;
    }

    std::string line;
    std::string key;
    int value;
    std::regex keyRegex("<key>([^<]+)</key>");
    std::regex valueRegex("<value>([^<]+)</value>");

    while (std::getline(inFile, line)) {
        std::smatch match;
        if (std::regex_search(line, match, keyRegex)) {
            key = match[1];
        } else if (std::regex_search(line, match, valueRegex)) {
            value = std::stoi(match[1]);
            index_events[key] = value;
        }
    }

    inFile.close();
    return 0;
}


void save_map_XML(const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    outFile << "<map>" << std::endl;
    for (const auto& pair : index_events) {
        outFile << "  <entry>" << std::endl;
        outFile << "    <key>" << pair.first << "</key>" << std::endl;
        outFile << "    <value>" << pair.second << "</value>" << std::endl;
        outFile << "  </entry>" << std::endl;
    }
    outFile << "</map>" << std::endl;

    outFile.close();
}

void create_map(int run_number, std::string inputDirFiles) {
    std::cout << "Loading map of files " << inputDirFiles << " ..." << std::endl;

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
            std::cout << "Event " << m_event_id_MC << " ... index " << tree_ientry << std::endl;
	    }
        tree_ientry++;
    }
    delete tree;
}

void convert_FASERMC(int run_number, std::string inputDirFiles, int min_event, int max_event) {

    std::cout << "Converting events ..." << std::endl;

    TFile *m_rootFile = new TFile(m_rootOutputFileName.c_str(), "RECREATE", "", 505); // last is the compression level
    if (!m_rootFile || !m_rootFile->IsOpen())
    {
        throw std::runtime_error("Could not create output ROOT file");
    }
    m_rootFile->cd();

    TPOEvent *branch_POEvent;
    TTree *m_POEventTree = new TTree("POEvent", "POEvent");
    m_calEventTree->Branch("event", &branch_POEvent);

    TChain *tree = new TChain("m_NuMCTruth_tree");
    tree->Add(inputDirFiles.c_str());  

    // FASER MC ntuple
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

    for (auto it : index_events) {
        int event = it.first;
        if( event < min_event || event > max_event) continue;

        int tree_ientry = it.second;

        TPOEvent fTPOEvent();
    
        bool must_end = false;
        bool got_primvtx = false;
        bool found_tau_lepton = false;
        int tau_lepton_track_id = 0;
        while (!must_end){
            tree->GetEntry(tree_ientry);
            if(event != m_event_id_MC) {
                if(TPOEvent.n_particles == 0){
                    std::cerr << " Something went wrong .. maybe wrong index? -- empty event" << std::endl;
                    must_end = true;
                }
                continue;
		    }
            fTPOEvent.run_number = m_runnumber;
	        fTPOEvent.event_id = m_event_id_MC;

    	    if(!got_primvtx && m_status == 1) {
    	        fTPOEvent.prim_vx[0] = m_vx_prod;
	            fTPOEvent.prim_vx[1] = m_vy_prod;
	            fTPOEvent.prim_vx[2] = m_vz_prod;
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
	        aPO.m_vx_decay = m_vx_prod-fTPOEvent.prim_vx[0];
	        aPO.m_vy_decay = m_vy_prod-fTPOEvent.prim_vx[1];
	        aPO.m_vz_decay = m_vz_prod-fTPOEvent.prim_vx[2];
	        aPO.nparent = m_trackid_in_particle->size();
	        for (int i=0; i<aPO.nparent;i++){
	            aPO.m_trackid_in_particle[i] = m_trackid_in_particle->at(i);
    	    };
    	    aPO.m_status = m_status;
    		aPO.geanttrackID = -1;
	    
    	    if(m_track_id < 20000 && m_status != 3) {
    	        fTPOEvent.POs[fTPOEvent.n_particles++] = aPO;
	        }
	    
	        // found charged tau lepton - store decay products
	        if(!found_tau_lepton && abs(aPO.m_pdg_id) == 15) {
	            found_tau_lepton = true;
	            tau_lepton_track_id = m_track_id;
	        }
	    
	        if(found_tau_lepton) {
	            for(int i=0;i<aPO.nparent;i++) {
				    if(aPO.m_trackid_in_particle[i] == tau_lepton_track_id) {
				        fTPOEvent.taudecay[fTPOEvent.n_taudecay++] = aPO;
				        double decaylength = sqrt(aPO.m_vx_decay*aPO.m_vx_decay+aPO.m_vy_decay*aPO.m_vy_decay+aPO.m_vz_decay*aPO.m_vz_decay);
				        fTPOEvent.tautracklength = decaylength;
				    }
                }
	        }
        } // while
	
    	if(fTPOEvent.istau && fTPOEvent.n_taudecay==0) {
    	    std::cout << "Could not find tau decay product??" << std::endl;
    	}

    	fTPOEvent.kinematics_event();
    	fTPOEvent.dump_header();
    	fTPOEvent.dump_event();	

        // should write to output
        branch_POEvent = &fTPOEvent;
        branch_POEvent -> Fill();
    
        delete fTPOEvent;
    }
    m_calEventTree->Write();
    m_rootFile->Close();
}

bool confirmAction() {
    char response;
    while (true) {
        std::cout << "Do you want to proceed? (y/n): ";
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
	if (argc != 1) {
		std::cout << "Usage: " << argv[0] << " <run> " << std::endl;
		return 1;
	}

//    std::string inputDir = argv[1];
//	std::string ROOToutputFile = argv[2];

    std::string runString("200026");

    int run_number;
    try {
        run_number = std::stoi(runString);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument: " << e.what() << std::endl;
    } catch (const std::out_of_range& e) {
        std::cerr << "Out of range: " << e.what() << std::endl;
    }

    std::ostringstream inputDirFiles;
    inputDirFiles << "../FASERDATA/sim" << run_number << "/*.root";
    std::ostringstream ROOTOutputFile;
    ROOTOutputFile << "run" << run_number << ".root";
    std::ostringstream MapFile;
    MapFile << "map_run" << run_number << ".xml";

    std::cout << "Converting FASERMC from " << inputDirFiles.str() << std::endl;
    std::cout << "The event indices map is " << MapOutputFile.str() << std::endl;

    if(load_map_XML(MapFile.str()) != 0) {
        if (confirmAction()) {
            create_map(run_number, inputDirFiles.str());
            save_map_XML(MapFile.str());
        }
    }

}
