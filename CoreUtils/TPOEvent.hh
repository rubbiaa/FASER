/* 
  FASER ntuple interface library
  A. Rubbia May 2024 
*/

#ifndef _TPOEVENT_
#define _TPOEVENT_ 1

#include <TDatabasePDG.h>
#include <Math/PositionVector3D.h>
#include <Math/Vector3D.h>
#include <Math/Point3D.h>

#define MAXPARENT 10

/// @brief Particle Object structure to hold generator level (truth) information.
struct PO { // : public TObject {
public:
  int m_pdg_id;                  // PDG code 
  int m_track_id;                // track ID within MC
  double m_px;                   // x momentum
  double m_py;                   // y momentum
  double m_pz;                   // z momentum
  double m_energy;               // total energy
  double m_kinetic_energy;       // kinetic energy
  double m_vx_decay, m_vy_decay, m_vz_decay;    // this is mislabeled - it's actually vertex of particle
  int nparent;                   // number of parents
  int m_trackid_in_particle[MAXPARENT]; // list of parents
  int m_status;                  // MC status

  int geanttrackID;              // GEANT4 track id of this primary

#if 1
  /// @brief Compute mass of the PO
  /// @return mass
  double m_mass() const { 
    double mass2 = m_energy*m_energy-m_px*m_px-m_py*m_py-m_pz*m_pz;
    return sqrt(std::max(mass2,0.0));
  };   //!
  double m_charge() const {
    TDatabasePDG *pdgDB = TDatabasePDG::Instance();
    TParticlePDG *particle = pdgDB->GetParticle(m_pdg_id);
    if(particle != nullptr) return particle->Charge();
    return 0;
  };   //!
#endif

//  ClassDef(PO,1)

};

/// @brief Particle Object Event to hold a full generator level (truth) event.
/// Many entries of compute by calling kinematics_event().
class TPOEvent : public TObject
{
private:
  struct stats
  {
    size_t nueCC;
    size_t numuCC;
    size_t nutauCC;
    size_t NC;
    size_t ES;
  } stats;         //!

public:

  // target constants (not saved in ROOT I/O)
  static const int kVtx_in_W = 1;          //!
  static const int kVtx_in_Scint = 2;      //!

  static const int kMask_nueCC = 1;
  static const int kMask_numuCC = 2;
  static const int kMask_nutauCC = 3;
  static const int kMask_NC = 4;
  static const int kMask_ES = 5;

  int run_number;                   // run number
  int event_id;                     // event number
  ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> prim_vx;  // primary vertex in mm
  bool use_GENIE_vtx = false;       // tell FASERG4 to use prim_vx (if true)
  std::string GENIE_vtx_name = "";  // the name of the hit nucleus
  int vtx_target;                   // in which target did the interaction occur
  int event_mask = 0;               // if events are masked (see kMask_... constants)
  bool isCC;                        // event is a charged current
  bool isES() const { return jetpx==0 && jetpy == 0 && jetpz == 0;}; // event is elastic scattering on electron
  bool isCharmed() const;
  bool istau;                       // incoming neutrino is a nutau
  int tau_decaymode;                // =1 e, =2 mu, =3 1-prong, =4 rho =5 3-prong, =6 other
  std::vector<struct PO> POs;       // vector of PO of the event
  std::vector<struct PO> taudecay;  // vector of the tau decay products
  double tauDecaylength();            // tau track length (in mm)
  double tauKinkAngle();            // Kink angle of tau decay products

  // Kinematics of the event
  struct PO in_neutrino;            // PO of incoming neutrino
  struct PO out_lepton;             // PO of outgoing lepton
  double spx, spy, spz;             // Sum of momenta of particles (except final state neutrinos)
  double vis_spx, vis_spy, vis_spz; // same as spx,spy,spz except replaces tau by its decay products excluding neutrinos
  double jetpx, jetpy, jetpz;       // spx,spy,spz minus outgoing lepton
  double tauvis_px, tauvis_py, tauvis_pz; // Sum of tau decay products except neutrinos
  double Evis, ptmiss;              // Visible energy and miss transverse momentum

  /// @brief Main constructor
  TPOEvent() { clear_event(); };

  /// @brief Reset PO structure to empty event
  void clear_event();

  /// @brief Check if PDGid is a lepton (e,mu,tau or neutrinos)
  /// @param pdgid 
  /// @return true if lepton
  static bool is_lepton(int pdgid) {
      int pdgidabs = abs(pdgid);
      return (pdgidabs >= 11 && pdgidabs <= 16);
  }

  /// @brief Check if PDGid is a neutrino (any flavor)
  /// @param pdgid 
  /// @return true if neutrino
  static bool is_neutrino(int pdgid) {
      int pdgidabs = abs(pdgid);
      return (pdgidabs == 12 || pdgidabs == 14 || pdgidabs == 16);
  }

  /// @brief Return number of POs in the event
  /// @return Number of POs
  size_t n_particles() const { return POs.size(); };

  /// @brief Return number of electrically charged particles in the event
  /// @return Number of electrically charged POs.
  size_t n_charged() const;

  /// @brief Perform the decay of the charged tau lepton and store decay products
  /// @param tauPO - the PO of the tau to be decayed
#ifdef _INCLUDE_PYTHIA_
  void perform_taulepton_decay(struct PO tauPO);
#endif

  /// @brief Return number of tau decay products
  /// @return Number of tau decay products
  size_t n_taudecay() const { return taudecay.size(); };

  /// @brief Set primary vertex of the event
  /// @param x in mm
  /// @param y in mm
  /// @param z in mm
  void setPrimaryVtx(double x, double y, double z) { prim_vx.SetCoordinates(x,y,z);};

  void setVtxTarget(int target) { vtx_target = target; };

  /// @brief Compute event kinematics including Evis, ptmiss, etc. Important to call.
  void kinematics_event();

  /// @brief Nicely formatted header of the event on std::cout
  void dump_header(std::ostream& out = std::cout) const;

  /// @brief Nicely formatted dump of the full event on std::cout
  void dump_event(std::ostream& out = std::cout) const;

  /// @brief Nicely formatted dump of a single PO on std::cout
  /// @param aPO 
  /// @param pdgDB 
  void dump_PO(struct PO aPO,  TDatabasePDG *pdgDB, std::ostream& out = std::cout) const;

  /// @brief Return brief description (nueCC, numuCC, nutauCC, NC, etc..)
  /// @return The text describing the reaction
  const char* reaction_desc() const;

  /// @brief Set the GEANT4 track id corresponding to the PO
  /// @param iPO Index of PO in POs vector
  /// @param trackID GEANT4 tracj id
  void setGEANT4TrackID(int iPO, int trackID) { POs[iPO].geanttrackID = trackID; };

  /// @brief Return PO index in POs of the track with the given GEANT4 trackID
  /// @param trackID 
  /// @return 
  int  findFromGEANT4TrackID(int trackID);

  int GetEventMask() const { return event_mask; }
  void SetEventMask(int mask) { event_mask = mask; };
  static int EncodeEventMask(std::string maskname)
  {
    if (maskname == "nueCC")
    {
      return kMask_nueCC;
    }
    else if (maskname == "numuCC")
    {
      return kMask_numuCC;
    }
    else if (maskname == "nutauCC")
    {
      return kMask_nutauCC;
    }
    else if (maskname == "nuNC")
    {
      return kMask_NC;
    }
    else if (maskname == "nuES")
    {
      return kMask_ES;
    }
    return -1;
  };

  static const char *DecodeEventMask(int mask)
  {
    switch (mask)
    {
    case kMask_nueCC:
      return "nueCC";
    case kMask_numuCC:
      return "numuCC";
    case kMask_nutauCC:
      return "nutauCC";
    case kMask_NC:
      return "nuNC";
    case kMask_ES:
      return "nuES";
    }
    return "unkmask";
  }

  void clone(TPOEvent *src);

  /// @brief Event statistics
  void reset_stats();
  void update_stats();
  void dump_stats();

  ClassDef(TPOEvent, 4)
};

#endif
