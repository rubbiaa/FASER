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
  double prim_vx[3];
  bool isCC;
  bool istau;
  int tau_decaymode; // =1 e, =2 mu, =3 1-prong, =4 rho =5 3-prong, =6 other
  size_t n_particles;
  size_t n_charged;
  struct PO in_neutrino;
  struct PO out_lepton;
  struct PO POs[MAXPARTICLES];
  size_t n_taudecay;
  struct PO taudecay[MAXPARTICLES];
  double tautracklength;
  double spx, spy, spz;
  double vis_spx, vis_spy, vis_spz;
  double jetpx, jetpy, jetpz;
  double tauvis_px, tauvis_py, tauvis_pz;
  double Evis, ptmiss;
  double cost, cosf;
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
