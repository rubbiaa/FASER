// Simple CLI combinator for ROOT TTrees using TChain
// Usage examples:
//   CombineFluxes -o output.root input1.root input2.root ...
//   CombineFluxes -t gFaser -o output.root input*.root

#include <TFile.h>
#include <TChain.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <vector>

static void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [-t <treeName>] -o <output.root> <input1.root> [input2.root ...]\n"
              << "  -t <treeName>  : Name of the TTree to merge (default: gFaser)\n"
              << "  -o <output>    : Output ROOT file name (required)\n"
              << "  inputs         : One or more input ROOT files (wildcards allowed)\n";
}

int main(int argc, char** argv) {
    std::string treeName = "gFaser";
    std::string outName;
    std::vector<std::string> inputs;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "-t" && i + 1 < argc) {
            treeName = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            outName = argv[++i];
        } else {
            inputs.push_back(arg);
        }
    }

    if (outName.empty() || inputs.empty()) {
        std::cerr << "Error: missing required arguments.\n";
        print_usage(argv[0]);
        return 1;
    }

    std::cout << "Merging tree '" << treeName << "' from " << inputs.size() << " file(s)\n";

    TChain chain(treeName.c_str());
    for (const auto& f : inputs) {
        int added = chain.Add(f.c_str());
        if (added == 0) {
            std::cerr << "Warning: no files matched or could not add '" << f << "'\n";
        } else {
            std::cout << "  + added: " << f << "\n";
        }
    }

    Long64_t nEntries = chain.GetEntries();
    std::cout << "Total entries in combined chain: " << nEntries << std::endl;

    if (nEntries <= 0) {
        std::cerr << "Error: no entries to write. Aborting.\n";
        return 2;
    }

    TFile outputFile(outName.c_str(), "RECREATE");
    if (outputFile.IsZombie()) {
        std::cerr << "Error: could not create output file '" << outName << "'\n";
        return 3;
    }

    TTree* outputTree = chain.CloneTree();
    if (!outputTree) {
        std::cerr << "Error: failed to clone tree.\n";
        return 4;
    }

    outputTree->Write();
    outputFile.Close();

    std::cout << "Files combined and saved to '" << outName << "'\n";
    return 0;
}



