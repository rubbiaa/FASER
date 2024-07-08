#ifndef _TPOEVENT_
#define _TPOEVENT_ 1

#include <TDatabasePDG.h>

#define MAXPARENT 10
struct PO {
  int m_pdg_id;
  int m_track_id;
  double m_px;
  double m_py;
  double m_pz;
  double m_energy;
  double m_kinetic_energy;
  double m_vx_decay, m_vy_decay, m_vz_decay;
  int nparent;
  int m_trackid_in_particle[MAXPARENT];
  int m_status;
  int geanttrackID;

  ClassDef(PO,1)

};

#define MAXPARTICLES 1000
class TPOEvent : public TObject {

public:
  int run_number;
  int event_id;
  double prim_vx[3];
  bool isCC;
  bool istau;
  int tau_decaymode; // =1 e, =2 mu, =3 1-prong, =4 rho =5 3-prong, =6 other
  size_t n_particles;
  struct PO in_neutrino;
  struct PO out_lepton;
  struct PO POs[MAXPARTICLES];       // TODO: move to vector
  size_t n_taudecay;
  struct PO taudecay[MAXPARTICLES];       // TODO: move to vector
  double tautracklength;
  double spx, spy, spz;
  double vis_spx, vis_spy, vis_spz;
  double jetpx, jetpy, jetpz;
  double tauvis_px, tauvis_py, tauvis_pz;
  double Evis, ptmiss;

  void clear_event();
  bool is_lepton(int pdgid);
  bool is_neutrino(int pdgid);
  void kinematics_event();
  void dump_PO(struct PO aPO,  TDatabasePDG *pdgDB) const;
  void dump_header() const;
  void dump_event() const;

  void setGEANT4TrackID(int iPO, int trackID) { POs[iPO].geanttrackID = trackID; };
  int  findFromGEANT4TrackID(int trackID);
  void setPrimaryVtx(double x, double y, double z) { prim_vx[0]=x; prim_vx[1]=y; prim_vx[2]=z;};

  ClassDef(TPOEvent, 1)
};

#endif