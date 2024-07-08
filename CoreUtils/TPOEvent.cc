// FASER ntuple interface library
// A. Rubbia May 2024
//

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <map>

#include "TPOEvent.hh"


ClassImp(TPOEvent)
ClassImp(PO)

void TPOEvent::clear_event() {
  run_number = event_id = -1;
  prim_vx[0] = prim_vx[1] = prim_vx[2] = 0;
  n_particles = 0;
  n_taudecay = 0;
  tau_decaymode = -1;
  isCC = false;
  istau = false;
  spx=spy=spz=0;
  tauvis_px=tauvis_py=tauvis_pz=0;
};

bool TPOEvent::is_lepton(int pdgid) {
  int pdgidabs = abs(pdgid);
  return (pdgidabs >= 11 && pdgidabs <= 16);
}

bool TPOEvent::is_neutrino(int pdgid) {
  int pdgidabs = abs(pdgid);
  return (pdgidabs == 12 || pdgidabs == 14 || pdgidabs == 16);
}

void TPOEvent::kinematics_event() {
  bool got_out_lepton = false;
  for (size_t i=0; i<n_particles; i++) {
    struct PO aPO = POs[i];
    if(aPO.m_status == 4 && i==0) {
      in_neutrino = aPO;
      istau = (abs(aPO.m_pdg_id) == 16);
    }
    if(!got_out_lepton && aPO.m_status == 1 && is_lepton(aPO.m_pdg_id) ) {
      out_lepton = aPO;
      got_out_lepton = true;
    }
    if(aPO.m_status == 1 && !is_neutrino(aPO.m_pdg_id)) {
      spx += aPO.m_px;
      spy += aPO.m_py;
      spz += aPO.m_pz;
    }
  }
  jetpx = spx-out_lepton.m_px;
  jetpy = spy-out_lepton.m_py;
  jetpz = spz-out_lepton.m_pz;
  vis_spx = spx;
  vis_spy = spy;
  vis_spz = spz;
  isCC = !(in_neutrino.m_pdg_id == out_lepton.m_pdg_id);
  if(istau && isCC){
    TDatabasePDG *pdgDB = TDatabasePDG::Instance();
    int nc = 0, nn = 0;
    for(int i=0; i<n_taudecay; i++) {
      struct PO aPO = taudecay[i];
      TParticlePDG *particle = pdgDB->GetParticle(aPO.m_pdg_id);
      if(aPO.m_status == 1 && !is_neutrino(aPO.m_pdg_id)){
	if(abs(aPO.m_pdg_id) == 11) {
	  tau_decaymode = 1;
	}
	if(abs(aPO.m_pdg_id) == 13) {
	  tau_decaymode = 2;
	}
	if(particle->Charge() == 0){
	  nn++;
	} else {
	  nc++;
	}
	tauvis_px += aPO.m_px;
	tauvis_py += aPO.m_py;
	tauvis_pz += aPO.m_pz;
      }
    }
    if(tau_decaymode < 0){
      tau_decaymode = 6;
      if(nc==1 && nn == 0) tau_decaymode = 3;
      if(nc==1 && nn == 1) tau_decaymode = 4;
      if(nc==3) tau_decaymode = 5;
    }
    vis_spx = tauvis_px + jetpx;
    vis_spy = tauvis_py + jetpy;
    vis_spz = tauvis_pz + jetpz;
  }
  Evis = sqrt(vis_spx*vis_spx + vis_spy*vis_spy + vis_spz*vis_spz);
  ptmiss = sqrt(vis_spx*vis_spx + vis_spy*vis_spy);
}


void TPOEvent::dump_PO(struct PO aPO,  TDatabasePDG *pdgDB) const {
  TParticlePDG *particle = pdgDB->GetParticle(aPO.m_pdg_id);
  std::cout << std::setw(10) << aPO.m_track_id << " " << std::setw(12) << aPO.m_pdg_id << " ";
  if(particle) {
    std::cout << std::setw(10) << particle->GetName() << " ";
  } else {
    std::cout << "        ?? ";
  }
//  double decaylength = sqrt(aPO.m_vx_decay*aPO.m_vx_decay+aPO.m_vy_decay*aPO.m_vy_decay+aPO.m_vz_decay*aPO.m_vz_decay);
//  std::cout << std::setw(10) << decaylength << " ";
//  std::cout << std::setw(10) << aPO.m_vx_decay << " " << std::setw(10) << aPO.m_vy_decay << " " << std::setw(10) << aPO.m_vz_decay << " ";
  std::cout << std::setw(12) << aPO.m_px << " " << " " << std::setw(12) << aPO.m_py << " " << std::setw(12) << aPO.m_pz 
          << " " << std::setw(12) << aPO.m_energy << " " << std::setw(4) << aPO.m_status << " " << std::setw(4) << aPO.geanttrackID << " ";
  for (size_t j=0; j<aPO.nparent; j++) {
    std::cout << std::setw(10) << aPO.m_trackid_in_particle[j] << " ";
  }
  std::cout << std::endl;
}

void TPOEvent::dump_header() const {
  if(isCC) {
    std::cout << "--- Run " << run_number << " Event " << event_id << " ---------------------------------------- CC --------------------------------------------" << std::endl;
  } else {
    std::cout << "--- Run " << run_number << " Event " << event_id << " ------------------------------------------- NC ------------------------------------------" << std::endl;
  }
}
 
void TPOEvent::dump_event() const {
  double spx=0, spy=0, spz=0;
  TDatabasePDG *pdgDB = TDatabasePDG::Instance();
  dump_header();
  std::cout << " Primary vtx = " << prim_vx[0] << " " << prim_vx[1] << " " << prim_vx[2] << std::endl;
  std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "¨    trackID, pdg_ID, name, px, py, pz, E, status, geant4ID, parents" << std::endl;
  for (size_t i=0; i<n_particles; i++) {
    struct PO aPO = POs[i];
    dump_PO(aPO, pdgDB);
  }
  std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
  std::cout << std::setw(10) << "Outgoing lepton: " << std::endl;
  dump_PO(out_lepton, pdgDB);
  std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
  std::cout << std::setw(10) << "Jet :                      " << jetpx << " " << jetpy << " " << jetpz << std::endl;
  std::cout << std::setw(10) << "Sum final state particles: " << spx << " " << spy << " " << spz << std::endl;
  std::cout << std::setw(10) << "Sum final state particles (VIS): " << vis_spx << " " << vis_spy << " " << vis_spz << std::endl;
  std::cout << std::setw(10) << "Ptmiss = " << ptmiss << "  Evis = " << Evis << std::endl;
  std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
  if(n_taudecay>0) {
    std::cout << "Tau decay mode : " << tau_decaymode << std::endl;
    for (size_t i=0; i<n_taudecay; i++) {
      struct PO aPO = taudecay[i];
      dump_PO(aPO, pdgDB);
    }
    std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
  }
}

int TPOEvent::findFromGEANT4TrackID(int trackID) {
  for (size_t i=0; i<n_particles; i++) {
    if(POs[i].geanttrackID == trackID){
      return i;
    }
  }
  return -1;
}