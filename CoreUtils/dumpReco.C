
#include "TPORecoEvent.hh"

void dumpReco()
{

    gSystem->AddIncludePath("-I../CoreUtils");
    gSystem->AddIncludePath("-I../GenFit-install/include");

    gSystem->Load("libCoreUtilsDict.so");

    TChain *event_tree = new TChain("RecoEvent", "READ");
    event_tree->Add("../Batch/Batch-TPORecevent_301_0_100.root");

    TPORecoEvent *reco = nullptr;
    event_tree->SetBranchAddress("TPORecoEvent", &reco);

    Long_t nentries = event_tree->GetEntries();
    std::cout << "Number of entries " << nentries << std::endl;
    Long_t ientry = 0;
    for (ientry = 0; ientry < nentries; ientry++)
    {
        event_tree->GetEntry(ientry++);

        std::cout << *reco << std::endl;
    }
}