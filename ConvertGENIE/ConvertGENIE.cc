// Convert FASER GENIE human output to POevents
//
// A. Rubbia, August 2024
//

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>

#include <TTree.h>
#include <TChain.h>
#include <TFile.h>

#include "TPOEvent.hh"

void convert_FASERMC(int run_number, TTree *tree, int min_event, int max_event,
                     std::string ROOTOutputFile, int mask)
{
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
  m_POEventTree->Branch("event", &branch_POEvent); // this should be named POEvent

  // FASER GENIE ntuple

  Double_t vx;
  Double_t vy;
  Double_t vz;
  Int_t n;
  std::vector<std::string> *name = nullptr;
  std::vector<int> *pdgc = nullptr;
  std::vector<int> *status = nullptr;
  std::vector<int> *firstMother = nullptr;
  std::vector<int> *lastMother = nullptr;
  std::vector<int> *firstDaughter = nullptr;
  std::vector<int> *lastDaughter = nullptr;
  std::vector<double> *px = nullptr;
  std::vector<double> *py = nullptr;
  std::vector<double> *pz = nullptr;
  std::vector<double> *E = nullptr;
  std::vector<double> *m = nullptr;
  std::vector<double> *M = nullptr;

  // List of branches
  TBranch *b_vx;            //!
  TBranch *b_vy;            //!
  TBranch *b_vz;            //!
  TBranch *b_n;             //!
  TBranch *b_name;          //!
  TBranch *b_pdgc;          //!
  TBranch *b_status;        //!
  TBranch *b_firstMother;   //!
  TBranch *b_lastMother;    //!
  TBranch *b_firstDaughter; //!
  TBranch *b_lastDaughter;  //!
  TBranch *b_px;            //!
  TBranch *b_py;            //!
  TBranch *b_pz;            //!
  TBranch *b_E;             //!
  TBranch *b_m;             //!
  TBranch *b_M;             //!

  tree->SetBranchAddress("vx", &vx);
  tree->SetBranchAddress("vy", &vy);
  tree->SetBranchAddress("vz", &vz);
  tree->SetBranchAddress("n", &n);
  tree->SetBranchAddress("name", &name);
  tree->SetBranchAddress("pdgc", &pdgc, &b_pdgc);
  tree->SetBranchAddress("status", &status, &b_status);
  tree->SetBranchAddress("firstMother", &firstMother);
  tree->SetBranchAddress("lastMother", &lastMother);
  tree->SetBranchAddress("firstDaughter", &firstDaughter);
  tree->SetBranchAddress("lastDaughter", &lastDaughter);
  tree->SetBranchAddress("px", &px);
  tree->SetBranchAddress("py", &py);
  tree->SetBranchAddress("pz", &pz);
  tree->SetBranchAddress("E", &E);
  tree->SetBranchAddress("m", &m);
  tree->SetBranchAddress("M", &M);

  int evt_to_dump = 0;
  size_t iseq = 0;

  for (size_t event = min_event; event < max_event; event++)
  {

    tree->GetEntry(event);

    if (event == 0) fTPOEvent.reset_stats();
    if (event % 1000 == 0)
    {
      std::cout << "Processing event " << event << " ..." << std::endl;
    }

    fTPOEvent.clear_event();
    fTPOEvent.run_number = run_number;
    fTPOEvent.setPrimaryVtx(vx * 1e3, vy * 1e3, vz * 1e3); // convert from meters to mm
    fTPOEvent.use_GENIE_vtx = true;     // tell FASERG4 to use this vtx

    // SKIP VERTICES IN REAR CAL
    if(fTPOEvent.prim_vx.z() > 1533.0) continue;

    fTPOEvent.event_id = iseq++;

    bool found_tau_lepton = false;

    for (size_t i = 0; i < n; i++)
    {
      struct PO aPO;
      aPO.m_track_id = i;
      aPO.m_pdg_id = pdgc->at(i);
      if(i==1) {
        fTPOEvent.GENIE_vtx_name = name->at(i);
      }
      aPO.m_status = status->at(i);
      if (aPO.m_status == 0)
        aPO.m_status = 4;
      aPO.m_px = px->at(i);
      aPO.m_py = py->at(i);
      aPO.m_pz = pz->at(i);
      aPO.m_energy = E->at(i);
      aPO.geanttrackID = -1;
      aPO.nparent = 0;
      int iMo = firstMother->at(i);
      if (iMo > -1)
      {
        aPO.nparent = 1;
        aPO.m_trackid_in_particle[0] = iMo;
      }

      if (aPO.m_status == 1 || aPO.m_status == 4)
      {
        fTPOEvent.POs.push_back(aPO);
      }

      if (!found_tau_lepton && abs(aPO.m_pdg_id) == 15)
      {
        found_tau_lepton = true;
#ifdef _INCLUDE_PYTHIA_
        fTPOEvent.perform_taulepton_decay(aPO);
#endif
      }
    }

    fTPOEvent.kinematics_event();
    if (evt_to_dump++ < 20)
    {
      fTPOEvent.dump_event();
    };
    fTPOEvent.update_stats();
    m_POEventTree->Fill();
  }

  m_POEventTree->Write();
  m_rootFile->Close();
  std::cout << "Done saving..." << std::endl;

  fTPOEvent.dump_stats();
  std::cout << "Total number of events written " << iseq << std::endl;

}

int main(int argc, char **argv)
{
  // get the output file name as the first argument
  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0] << " <genieroot> <run>" << std::endl;
    std::cout << std::endl;
    std::cout << "  <genieroot>                The input FASER GENIE root files" << std::endl;
    std::cout << "  <run>                      The output run number" << std::endl;
#if 0
    std::cout << "Options:" << std::endl;
    std::cout << "   mask                      To process only specific events (def=none): ";
    std::cout << "  nueCC, numuCC, nutauCC, nuNC or nuES" << std::endl;
#endif
    return 1;
  }

  std::string rootinputString = argv[1];
  std::string runString = argv[2];

  int run_number;
  try
  {
    run_number = std::stoi(runString);
  }
  catch (const std::invalid_argument &e)
  {
    std::cerr << "Invalid argument for run: " << e.what() << std::endl;
  }
  catch (const std::out_of_range &e)
  {
    std::cerr << "Out of range for run: " << e.what() << std::endl;
  }

  int min_event = 0;
  int max_event = -1;
  int event_mask = 0;

  std::ostringstream inputDirFiles;
  inputDirFiles << rootinputString;

  TChain *tree = new TChain("gFaser");
  tree->Add(inputDirFiles.str().c_str());
  size_t n_entries = tree->GetEntries();
  if (max_event == -1)
  {
    max_event = n_entries;
  }

  std::ostringstream ROOTOutputFile;
  ROOTOutputFile << "FASERMC-PO-Run" << run_number << "-" << min_event << "_" << max_event;
  if (event_mask > 0)
  {
    const char *mask = TPOEvent::DecodeEventMask(event_mask);
    ROOTOutputFile << "_" << mask;
  }
  ROOTOutputFile << ".root";

  std::cout << "Converting FASERMC from " << inputDirFiles.str() << std::endl;
  std::cout << "The output file is " << ROOTOutputFile.str() << std::endl;

  convert_FASERMC(run_number, tree, min_event, max_event,
                  ROOTOutputFile.str(), event_mask);

  std::cout
      << "I'm done." << std::endl;
  return 0;
}
