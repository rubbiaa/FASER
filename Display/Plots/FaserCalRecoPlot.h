#ifndef FaserCalRecoPlot_h
#define FaserCalRecoPlot_h

#include <TROOT.h>
#include <TTree.h>
#include <TBranch.h>
#include <TNamed.h>
#include <TObject.h>
#include <TFile.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <TTreePerfStats.h>
#include "TVector3.h"

#include "string"
#include "vector"


using namespace std;

class FaserCalRecoPlot : public TObject {
 public :
  // Declaration of leaf types
  static const Int_t kMaxVtx = 100;
  Int_t           Run;
  Int_t           Event;
  Double_t        Vtx_pry;
  Double_t        Vty_pry;
  Double_t        Vtz_pry;
  Int_t           TotalMultiplicity;
  Int_t           ChargedMultiplicity;
  Int_t           NeutralMultiplicity;
  Int_t           GammaMultiplicity;
  Int_t           NeutronMultiplicity;
  Double_t        NeutrinoEnergy;
  Int_t           CCNC;
  Int_t           NeutrinoType;
  Int_t           CharmType;
  Int_t           CharmProngs;
  string          *CharmDecayMode;
  string          *CharmParticle;
  Int_t           CharmDecay;
  Double_t        CharmEnergy;
  Double_t        Vtx_dcy;
  Double_t        Vty_dcy;
  Double_t        Vtz_dcy;
  Double_t        FlightLength;
  Int_t           n_vertices;
  Float_t         v_x[kMaxVtx];   //[n_vertices]
  Float_t         v_y[kMaxVtx];   //[n_vertices]
  Float_t         v_z[kMaxVtx];   //[n_vertices]
  Int_t           v_ntrks[kMaxVtx];   //[n_vertices]
  Int_t           t_closest_vertex;
  Int_t           n_tktracks;
  Int_t           n_pstracks;
  Int_t           n_clusters;
  Float_t         c_E1;
  Float_t         c_E1T;
  Float_t         c_chi2_1;
  Float_t         c_a_1;
  Float_t         c_b_1;
  Float_t         c_E2;
  Double_t        VisibleEnergy;
  Double_t        TotalEnergy;
  Double_t        RearECalEnergy;
  Double_t        RearHCalEnergy;
  Double_t        RearMuCalEnergy;
  vector<int>     *MuTag_alltrk;
  vector<double>  *MuTag_mom;
  vector<double>  *MuTag_px;
  vector<double>  *MuTag_py;
  vector<double>  *MuTag_pz;
  vector<double>  *MuTag_ene;
  vector<double>  *MuTag_dist;
  Int_t           MuTag_muon;
  vector<int>     *MuTag_muonSign;
  vector<double>  *MuTag_id;
  
  // List of branches
  TBranch        *b_Run;   //!
  TBranch        *b_Event;   //!
  TBranch        *b_Vtx_pry;   //!
  TBranch        *b_Vty_pry;   //!
  TBranch        *b_Vtz_pry;   //!
  TBranch        *b_TotalMultiplicity;   //!
  TBranch        *b_ChargedMultiplicity;   //!
  TBranch        *b_NeutralMultiplicity;   //!
  TBranch        *b_GammaMultiplicity;   //!
  TBranch        *b_NeutronMultiplicity;   //!
  TBranch        *b_NeutrinoEnergy;   //!
  TBranch        *b_CCNC;   //!
  TBranch        *b_NeutrinoType;   //!
  TBranch        *b_CharmType;   //!
  TBranch        *b_CharmProngs;   //!
  TBranch        *b_CharmDecayMode;   //!
  TBranch        *b_CharmParticle;   //!
  TBranch        *b_CharmDecay;   //!
  TBranch        *b_CharmEnergy;   //!
  TBranch        *b_Vtx_dcy;   //!
  TBranch        *b_Vty_dcy;   //!
  TBranch        *b_Vtz_dcy;   //!
  TBranch        *b_FlightLength;   //!
  TBranch        *b_n_vertices;   //!
  TBranch        *b_v_x;   //!
  TBranch        *b_v_y;   //!
  TBranch        *b_v_z;   //!
  TBranch        *b_v_ntrks;   //!
  TBranch        *b_t_closest_vertex;   //!
  TBranch        *b_n_tktracks;   //!
  TBranch        *b_n_pstracks;   //!
  TBranch        *b_n_clusters;   //!
  TBranch        *b_c_E1;   //!
  TBranch        *b_c_E1T;   //!
  TBranch        *b_c_chi2_1;   //!
  TBranch        *b_c_a_1;   //!
  TBranch        *b_c_b_1;   //!
  TBranch        *b_c_E2;   //!
  TBranch        *b_VisibleEnergy;   //!
  TBranch        *b_TotalEnergy;   //!
  TBranch        *b_RearECalEnergy;   //!
  TBranch        *b_RearHCalEnergy;   //!
  TBranch        *b_RearMuCalEnergy;   //!
  TBranch        *b_MuTag_alltrk;   //!
  TBranch        *b_MuTag_mom;   //!
  TBranch        *b_MuTag_px;   //!
  TBranch        *b_MuTag_py;   //!
  TBranch        *b_MuTag_pz;   //!
  TBranch        *b_MuTag_ene;   //!
  TBranch        *b_MuTag_dist;   //!
  TBranch        *b_MuTag_muon;   //!
  TBranch        *b_MuTag_muonSign;   //!
  TBranch        *b_MuTag_id;   //!

  FaserCalRecoPlot();
  // ~FaserCalRecoPlot(){};
~FaserCalRecoPlot() {}
  static FaserCalRecoPlot* giveThis();
  static FaserCalRecoPlot* giveThis(TTree* aTree, const std::string& option);
  static void releaseThis();
  void setupRead(const std::string& option = "");
  bool SetBranchAddresses();
  FaserCalRecoPlot(TTree* aTree, const std::string& option);
  TTree* GetInputTree();
  void GetEvent(Long64_t entry);

  // Basic Scalars
  Int_t       RunNumber()                        { return m_run;         };
  Int_t       EventNumber()                      { return m_event;         };
  TVector3    *GetPryVtx()      { return m_pryVtx; };
  TVector3    *GetDcyVtx()      { return m_dcyVtx; };
  TVector3    *GetRecoVtx(Int_t index)      { return m_recoVtx[index]; };
  Int_t     Multiplicity()           { return m_multip; };
  Int_t     Charged()    { return m_nch; };
  Int_t     Neutral()    { return m_nneutral; };
  Int_t     Gamma()    { return m_ngamma; };
  Int_t     Neutron()    { return m_nneutron; };
  Double_t  EneNu()         { return m_nuEne; };
  Int_t     CCNCType()              { return m_intType; };
  Int_t     NuType()           { return m_nuType; };
  Int_t     CharmParticleType()         { return m_charmType; };
  Int_t     CharmProngCount()        { return m_nprong; };
  Int_t      CharmDcy()       { return m_dcy; };
  std::string CharmDecayM()       { return m_dcymode; };
  std::string CharmParticleName()   { return m_charmP; };
  Double_t  CharmEne()            { return m_charmEne; };
  Double_t  FL()           { return m_charmFL; };
  Double_t  VisibleE()               { return m_eVis; };
  Double_t  TotalE()                 { return m_eTot; };
  Double_t  RearECalE()              { return m_rEcalEne; };
  Double_t  RearHCalE()              { return m_rHcalEne; };
  Double_t  RearMuCalE()             { return m_rMucalEne; };
  
  // Cluster Info
  Int_t     ClusterCount()           { return m_recoclstr; };
  Float_t   ClusterE1()              { return m_e1; };
  Float_t   ClusterE1T()             { return m_cE1T; };
  Float_t   ClusterChi2()            { return m_cchi2; };
  Float_t   ClusterParamA()          { return m_ca1; };
  Float_t   ClusterParamB()          { return m_cb1; };
  Float_t   ClusterE2()              { return m_cE2; };
  
  // Vertex Info
  Int_t     VertexCount()            { return m_nvtx; };
  Int_t     RecoTrackCount()         { return m_recotrks; };
  Int_t     PSRecoTrackCount()       { return m_recopstrks; };
  Double_t  ClosestRecoVtxIndex()    { return m_vtx_close; };
  
  // MuTag Data
  Int_t     MuonTagStatus() const         { return m_mutagmuon; };
  
  const std::vector<int>& MuTagTrks()     const { return m_MuTagTrks; };
  const std::vector<double>& MuTagMom()   const { return m_MuTagMom; };
  const std::vector<double>& MuTagEne()   const { return m_MuTagEne; };
  const std::vector<double>& MuTagDist()  const { return m_MuTagDist; };
  const std::vector<int>& MuTagMuSign()     const { return m_MuTagMuSign; };
  const std::vector<double>& MuTagID()    const { return m_MuTagID; };
  
 private:
  TTree* m_treeIn;
  int   m_run, m_event, m_dcy; 
  int m_multip, m_nch, m_nneutral, m_ngamma, m_nneutron;
  int m_intType, m_nuType, m_charmType, m_nprong, m_mutagmuon;
  int m_nvtx, m_ntrks[kMaxVtx], m_recotrks, m_recopstrks, m_recoclstr;
  std::string m_dcymode, m_charmP;
  double m_charmFL ;
  TVector3 *m_pryVtx, *m_dcyVtx, *m_recoVtx[kMaxVtx];
  vector<double> m_MuTagP, m_MuTagMom, m_MuTagEne, m_MuTagDist, m_MuTagID;
  vector<int> m_MuTagTrks, m_MuTagMuSign;
  double m_vtx_close; 
  double m_nuEne, m_charmEne, m_eVis, m_eTot, m_rEcalEne, m_rHcalEne, m_rMucalEne;
  double m_e1, m_cE1T, m_cchi2, m_ca1, m_cb1, m_cE2 ;
  static FaserCalRecoPlot* m_instance;
  ClassDef(FaserCalRecoPlot,0)
    };

// initialisation of the FaserCalRecoPlot pointer
FaserCalRecoPlot* FaserCalRecoPlot::m_instance = 0;
FaserCalRecoPlot::FaserCalRecoPlot(){}

// to get a unique instance 
inline FaserCalRecoPlot* FaserCalRecoPlot::giveThis()
{
  if (0 == m_instance){
    cout << "FaserCalRecoPlot::giveThis error not constructed properly " << endl;
  }
  return m_instance;
}

// Open a TTree for reading (writing option not now)
inline FaserCalRecoPlot* FaserCalRecoPlot::giveThis(TTree* aTree, const std::string& option)
{
  if (0 == m_instance){
    m_instance = new FaserCalRecoPlot(aTree, option);
  } else{
    cout << "FaserCalRecoPlot::giveThis Warning " << aTree->GetTitle() << endl;
  }
  return m_instance;
}

// Delete unique instance
inline void FaserCalRecoPlot::releaseThis() {
  if ( m_instance != 0 ) {
    delete m_instance;
    m_instance = 0;
  }
}

// constructor for one TTree with read/write option
FaserCalRecoPlot::FaserCalRecoPlot(TTree* aTree, const std::string& option){
  if(option=="read"){
    Long64_t nevent = 9999999;
    m_treeIn = aTree;
    setupRead();
  }
}

// setup the input tree for reading
void FaserCalRecoPlot::setupRead(const std::string& option){
  if(!m_treeIn){
    cout << "setupRead error: m_treeIn undefined " << endl;
    exit(1);
  }
  if(SetBranchAddresses()){}else{
    cerr << "TBranch error.."<< endl;
  }
}

// Setting up Branch content of input TTree
bool FaserCalRecoPlot::SetBranchAddresses(){
  CharmDecayMode = 0;
  CharmParticle = 0;
  MuTag_mom = 0;
  MuTag_px = 0;
  MuTag_py = 0;
  MuTag_pz = 0;
  MuTag_ene = 0;
  MuTag_dist = 0;
  MuTag_id = 0;
  MuTag_muonSign = 0;
  MuTag_alltrk = 0;
  // Set branch addresses and branch pointers
  m_treeIn->SetMakeClass(1);
  m_treeIn->SetBranchAddress("Run", &Run, &b_Run);
  m_treeIn->SetBranchAddress("Event", &Event, &b_Event);
  m_treeIn->SetBranchAddress("Vtx_pry", &Vtx_pry, &b_Vtx_pry);
  m_treeIn->SetBranchAddress("Vty_pry", &Vty_pry, &b_Vty_pry);
  m_treeIn->SetBranchAddress("Vtz_pry", &Vtz_pry, &b_Vtz_pry);
  m_treeIn->SetBranchAddress("TotalMultiplicity", &TotalMultiplicity, &b_TotalMultiplicity);
  m_treeIn->SetBranchAddress("ChargedMultiplicity", &ChargedMultiplicity, &b_ChargedMultiplicity);
  m_treeIn->SetBranchAddress("NeutralMultiplicity", &NeutralMultiplicity, &b_NeutralMultiplicity);
  m_treeIn->SetBranchAddress("GammaMultiplicity", &GammaMultiplicity, &b_GammaMultiplicity);
  m_treeIn->SetBranchAddress("NeutronMultiplicity", &NeutronMultiplicity, &b_NeutronMultiplicity);
  m_treeIn->SetBranchAddress("NeutrinoEnergy", &NeutrinoEnergy, &b_NeutrinoEnergy);
  m_treeIn->SetBranchAddress("CCNC", &CCNC, &b_CCNC);
  m_treeIn->SetBranchAddress("NeutrinoType", &NeutrinoType, &b_NeutrinoType);
  m_treeIn->SetBranchAddress("CharmType", &CharmType, &b_CharmType);
  m_treeIn->SetBranchAddress("CharmProngs", &CharmProngs, &b_CharmProngs);
  m_treeIn->SetBranchAddress("CharmDecayMode", &CharmDecayMode, &b_CharmDecayMode);
  m_treeIn->SetBranchAddress("CharmParticle", &CharmParticle, &b_CharmParticle);
  m_treeIn->SetBranchAddress("CharmDecay", &CharmDecay, &b_CharmDecay);
  m_treeIn->SetBranchAddress("CharmEnergy", &CharmEnergy, &b_CharmEnergy);
  m_treeIn->SetBranchAddress("Vtx_dcy", &Vtx_dcy, &b_Vtx_dcy);
  m_treeIn->SetBranchAddress("Vty_dcy", &Vty_dcy, &b_Vty_dcy);
  m_treeIn->SetBranchAddress("Vtz_dcy", &Vtz_dcy, &b_Vtz_dcy);
  m_treeIn->SetBranchAddress("FlightLength", &FlightLength, &b_FlightLength);
  m_treeIn->SetBranchAddress("n_vertices", &n_vertices, &b_n_vertices);
  m_treeIn->SetBranchAddress("v_x", v_x, &b_v_x);
  m_treeIn->SetBranchAddress("v_y", v_y, &b_v_y);
  m_treeIn->SetBranchAddress("v_z", v_z, &b_v_z);
  m_treeIn->SetBranchAddress("v_ntrks", v_ntrks, &b_v_ntrks);
  m_treeIn->SetBranchAddress("t_closest_vertex", &t_closest_vertex, &b_t_closest_vertex);
  m_treeIn->SetBranchAddress("n_tktracks", &n_tktracks, &b_n_tktracks);
  m_treeIn->SetBranchAddress("n_pstracks", &n_pstracks, &b_n_pstracks);
  m_treeIn->SetBranchAddress("n_clusters", &n_clusters, &b_n_clusters);
  m_treeIn->SetBranchAddress("c_E1", &c_E1, &b_c_E1);
  m_treeIn->SetBranchAddress("c_E1T", &c_E1T, &b_c_E1T);
  m_treeIn->SetBranchAddress("c_chi2_1", &c_chi2_1, &b_c_chi2_1);
  m_treeIn->SetBranchAddress("c_a_1", &c_a_1, &b_c_a_1);
  m_treeIn->SetBranchAddress("c_b_1", &c_b_1, &b_c_b_1);
  m_treeIn->SetBranchAddress("c_E2", &c_E2, &b_c_E2);
  m_treeIn->SetBranchAddress("VisibleEnergy", &VisibleEnergy, &b_VisibleEnergy);
  m_treeIn->SetBranchAddress("TotalEnergy", &TotalEnergy, &b_TotalEnergy);
  m_treeIn->SetBranchAddress("RearECalEnergy", &RearECalEnergy, &b_RearECalEnergy);
  m_treeIn->SetBranchAddress("RearHCalEnergy", &RearHCalEnergy, &b_RearHCalEnergy);
  m_treeIn->SetBranchAddress("RearMuCalEnergy", &RearMuCalEnergy, &b_RearMuCalEnergy);
  m_treeIn->SetBranchAddress("MuTag_alltrk", &MuTag_alltrk, &b_MuTag_alltrk);
  m_treeIn->SetBranchAddress("MuTag_mom", &MuTag_mom, &b_MuTag_mom);
  m_treeIn->SetBranchAddress("MuTag_px", &MuTag_px, &b_MuTag_px);
  m_treeIn->SetBranchAddress("MuTag_py", &MuTag_py, &b_MuTag_py);
  m_treeIn->SetBranchAddress("MuTag_pz", &MuTag_pz, &b_MuTag_pz);
  m_treeIn->SetBranchAddress("MuTag_ene", &MuTag_ene, &b_MuTag_ene);
  m_treeIn->SetBranchAddress("MuTag_dist", &MuTag_dist, &b_MuTag_dist);
  m_treeIn->SetBranchAddress("MuTag_muon", &MuTag_muon, &b_MuTag_muon);
  m_treeIn->SetBranchAddress("MuTag_muonSign", &MuTag_muonSign, &b_MuTag_muonSign);
  m_treeIn->SetBranchAddress("MuTag_id", &MuTag_id, &b_MuTag_id);
  m_treeIn->AddBranchToCache("*",kTRUE);
  return true;
}

// accessors for the input TTree
TTree* FaserCalRecoPlot::GetInputTree() {return m_treeIn;};
// get all branch contents of input TTree for entry i
void FaserCalRecoPlot::GetEvent(Long64_t entry)
{
  if (!m_treeIn) {
    cout << "FaserCalRecoPlot::getEntry error" << endl;
    exit(1);
  }
  //
  m_treeIn->GetEntry(entry);
  m_run       = Run; 
  m_event     = Event; 

  m_pryVtx = new TVector3();
  m_pryVtx->SetXYZ(Vtx_pry,Vty_pry,Vtz_pry);

  m_multip    = TotalMultiplicity;
  m_nch       = ChargedMultiplicity;
  m_nneutral  = NeutralMultiplicity;
  m_ngamma    = GammaMultiplicity;
  m_nneutron  = NeutronMultiplicity;
  m_nuEne     = NeutrinoEnergy;
  m_intType   =  CCNC;
  m_nuType    = NeutrinoType;
  m_charmType = CharmType;
  m_nprong    = CharmProngs;
  m_dcymode   = *CharmDecayMode;
  m_charmP    = *CharmParticle;
  m_charmEne  = CharmEnergy;
  m_dcy = CharmDecay;

  m_dcyVtx = new TVector3();
  m_dcyVtx->SetXYZ(Vtx_dcy,Vty_dcy,Vtz_dcy);

  m_charmFL   = FlightLength;
  m_nvtx      = n_vertices;

  for (int i = 0; i< n_vertices; i++)
    {  
      m_recoVtx[i] = new TVector3();
      m_recoVtx[i]->SetXYZ( v_x[i], v_y[i], v_z[i]);
      m_ntrks[i]    = v_ntrks[i];
    }
  m_vtx_close  = t_closest_vertex;

  m_recotrks   = n_tktracks;
  m_recopstrks = n_pstracks;
  m_recoclstr  = n_clusters;

  m_e1         = c_E1;
  m_cE1T       = c_E1T;
  m_cchi2      = c_chi2_1;
  m_ca1        = c_a_1;
  m_cb1        = c_b_1;
  m_cE2        = c_E2;

  m_eVis       = VisibleEnergy;
  m_eTot       = TotalEnergy;
  m_rEcalEne   = RearECalEnergy;
  m_rHcalEne   = RearHCalEnergy;
  m_rMucalEne   = RearMuCalEnergy;
  m_mutagmuon = MuTag_muon;
  
  // Copy from TTree vectors to internal class vectors
  m_MuTagMom    = *MuTag_mom;
  m_MuTagEne    = *MuTag_ene;
  m_MuTagDist   = *MuTag_dist;
  m_MuTagID     = *MuTag_id;
  m_MuTagTrks   = *MuTag_alltrk;
  m_MuTagMuSign = *MuTag_muonSign;


}

#endif
