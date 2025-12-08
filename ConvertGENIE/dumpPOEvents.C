#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TClass.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../CoreUtils/TPOEvent.hh"
#include "../CoreUtils/TPOEvent.cc" // include implementation so methods resolve

static TTree* FindPOEventTree(TFile* f, std::string& branchName) {
    if(!f) return nullptr;
    TIter keyIter(f->GetListOfKeys());
    while(auto *key = (TKey*)keyIter()) {
        if(strcmp(key->GetClassName(), "TTree")!=0) continue;
        TTree* t = (TTree*)key->ReadObj();
        std::vector<std::string> candidates = {"event", "poevent", "POEvent", "TPOEvent"};
        for(const auto& b : candidates) {
            if(t->GetBranch(b.c_str())) {
                TBranch* br = t->GetBranch(b.c_str());
                if(br && br->GetClassName() && std::string(br->GetClassName()).find("TPOEvent")!=std::string::npos) {
                    branchName = b; return t; }
            }
        }
        TObjArray* bl = t->GetListOfBranches();
        for(int i=0;i<bl->GetEntries();++i) {
            TBranch* br = (TBranch*)bl->At(i);
            if(br && br->GetClassName() && std::string(br->GetClassName()).find("TPOEvent")!=std::string::npos) {
                branchName = br->GetName(); return t; }
        }
    }
    return nullptr;
}

// Dump full event including particles
static void DumpFullEvent(const TPOEvent* ev, std::ostream& out, bool showPOs=true) {
    if(!ev) return;
    ev->dump_header(out);
    if(showPOs) {
        out << "Particles (trackID PDG name px py pz E status geant parents)" << std::endl;
        TDatabasePDG* pdgDB = TDatabasePDG::Instance();
        for(size_t i=0;i<ev->POs.size(); ++i) {
            ev->dump_PO(ev->POs[i], pdgDB, out);
        }
    }
    out << "Reaction: " << ev->reaction_desc() << "  CC=" << (ev->isCC?"yes":"no") << std::endl;
    out << "In neutrino PDG=" << ev->in_neutrino.m_pdg_id << " Out lepton PDG=" << ev->out_lepton.m_pdg_id << std::endl;
    out << "Summed momentum (spx spy spz) = (" << ev->spx << "," << ev->spy << "," << ev->spz << ")" << std::endl;
    out << "Visible momentum (vis_spx vis_spy vis_spz) = (" << ev->vis_spx << "," << ev->vis_spy << "," << ev->vis_spz << ")" << std::endl;
    out << "jet (px py pz) = (" << ev->jetpx << "," << ev->jetpy << "," << ev->jetpz << ")" << std::endl;
    out << "Evis=" << ev->Evis << " MeV  ptmiss=" << ev->ptmiss << " MeV" << std::endl;
    if(ev->istau) {
        out << "Tau decay mode=" << ev->tau_decaymode << " vis momentum (" << ev->tauvis_px << "," << ev->tauvis_py << "," << ev->tauvis_pz << ")" << std::endl;
    }
    out << "--------------------------------------------------------------------------------" << std::endl;
}

// Helper: parse comma-separated filter list into set
static std::set<std::string> ParseFilter(const std::string& s) {
    std::set<std::string> out;
    std::string token; std::istringstream iss(s);
    while(std::getline(iss, token, ',')) {
        // trim spaces
        size_t b = token.find_first_not_of(" \t\n");
        size_t e = token.find_last_not_of(" \t\n");
        if(b==std::string::npos) continue;
        std::string t = token.substr(b, e-b+1);
        if(!t.empty()) out.insert(t);
    }
    return out;
}

// Macro entry point with reaction filter.
// filterReactions: comma-separated list of reaction_desc strings (e.g. "numuCC,nuNC"). Empty = no filter.
// Example: root -l -b -q 'dumpPOEvents.C("file.root", 10, true, "numuCC,nueCC")'
void dumpPOEvents(const char* filename="FASERMC-PO-Run500-0_99261.root", Long64_t maxEvents=5, bool full=true, const char* filterReactions="") {
    TFile* f = TFile::Open(filename, "READ");
    if(!f || f->IsZombie()) { std::cerr << "Cannot open file " << filename << std::endl; return; }
    std::string branchName; TTree* tree = FindPOEventTree(f, branchName);
    if(!tree) { std::cerr << "No TPOEvent tree in file." << std::endl; f->Close(); return; }
    std::set<std::string> filters = ParseFilter(filterReactions ? filterReactions : "");
    bool useFilter = !filters.empty();
    std::cout << "Dumping up to " << maxEvents << " events from tree " << tree->GetName() << " branch " << branchName;
    if(useFilter) {
        std::cout << " with reaction filter: ";
        for(auto &r : filters) std::cout << r << " ";
    }
    std::cout << std::endl;

    TPOEvent* event=nullptr; tree->SetBranchAddress(branchName.c_str(), &event);
    Long64_t nEntries = tree->GetEntries();
    if(maxEvents>nEntries) maxEvents = nEntries;
    Long64_t dumped = 0;
    Long64_t scanned = 0;
    for(; scanned<nEntries && dumped<maxEvents; ++scanned) {
        tree->GetEntry(scanned);
        if(!event) continue;
        event->kinematics_event();
        const char* react = event->reaction_desc();
        if(useFilter && filters.find(react) == filters.end()) continue; // skip non-matching
        std::cout << "Event index=" << scanned << " run=" << event->run_number << " id=" << event->event_id << std::endl;
        DumpFullEvent(event, std::cout, full);
        ++dumped;
    }
    std::cout << "Summary: scanned=" << scanned << " dumped=" << dumped;
    if(useFilter) std::cout << " (filter applied)";
    std::cout << std::endl;
    f->Close();
}
