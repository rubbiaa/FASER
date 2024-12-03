#ifndef _FASERNTUPLIB_
#define _FASERNTUPLIB_ 1

#include <TDatabasePDG.h>

#define MAXPARENT 10
struct PO {
  int m_pdg_id;
  int m_track_id;
  double m_px;
  double m_py;
  double m_pz;
  double m_vx_decay, m_vy_decay, m_vz_decay;
  int nparent;
  int m_trackid_in_particle[MAXPARENT];
  int m_status;
};

#define MAXPARTICLES 1000
struct EVENT {
  int run_number;
  int event_id;
  Float_t prim_vx[3];
  bool isCC;
  bool isES; // true if elastic scattering off electrons
  bool istau;
  int tau_decaymode; // =1 e, =2 mu, =3 1-prong, =4 rho =5 3-prong, =6 other
  Int_t n_particles;
  Int_t n_charged;
  struct PO in_neutrino;
  struct PO out_lepton;
  struct PO POs[MAXPARTICLES];

  // tau decay 
  Int_t n_taudecay;
  struct PO taudecay[MAXPARTICLES];
  Float_t tautracklength;

  // selection
  Float_t taupi_cand[3];
  Float_t taurho_cand[3];

  // kinematical variables
  Float_t spx, spy, spz;
  Float_t vis_spx, vis_spy, vis_spz;
  Float_t jetpx, jetpy, jetpz;
  Float_t tauvis_px, tauvis_py, tauvis_pz;
  Float_t Evis, ptmiss;
  Float_t cost, cosf;

};

extern struct EVENT event;

void clear_event();
bool is_lepton(int pdgid);
bool is_neutrino(int pdgid);
void smear_event();
void kinematics_event();
void dump_PO(struct PO aPO,  TDatabasePDG *pdgDB);
void dump_event();

#endif
