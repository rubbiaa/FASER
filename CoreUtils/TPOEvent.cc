/* 
  FASER ntuple interface library
  A. Rubbia May 2024 
*/

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <map>

#include "TPOEvent.hh"

// #define _INCLUDE_PYTHIA_ 1
#ifdef _INCLUDE_PYTHIA_
#include "TPythia8.h"
#include <Pythia8/Event.h>
static  TPythia8 *fPythia8 = nullptr;
#endif

ClassImp(TPOEvent)
// ClassImp(PO)

void TPOEvent::clear_event() {
  run_number = event_id = -1;
  setPrimaryVtx(0,0,0);
  vtx_target = 0;
  POs.clear();
  taudecay.clear();
  tau_decaymode = -1;
  isCC = false;
  istau = false;
  spx=spy=spz=0;
  tauvis_px=tauvis_py=tauvis_pz=0;
};

#ifdef _INCLUDE_PYTHIA_
void TPOEvent::perform_taulepton_decay(struct PO tauPO) {
  taudecay.clear();
  istau = false;
  for (size_t i=0; i<n_particles(); i++) {
    struct PO aPO = POs[i];
    if(aPO.m_status == 4 && i==0) {
      istau = (abs(aPO.m_pdg_id) == 16);
    }
  }
  if(!istau) {
    std::cerr << "perform_taulepton_decay: event is not a nutauCC" << std::endl;
    return;
  }
  if(abs(tauPO.m_pdg_id) != 15) {
    std::cerr << "perform_taulepton_decay: PO id is not a charged tau lepton" << std::endl;
    return;
  }
  if(fPythia8 == nullptr) {
    fPythia8 = new TPythia8();

#if 1
    // this is the trick to make Pythia8 work as "decayer"
    //
    fPythia8->ReadString("ProcessLevel:all = off");

    fPythia8->ReadString("ProcessLevel:resonanceDecays=on");

    // shut off Pythia8 (default) verbosty
    //
    fPythia8->ReadString("Init:showAllSettings=false");
    fPythia8->ReadString("Init:showChangedSettings=false");
    fPythia8->ReadString("Init:showAllParticleData=false");
    fPythia8->ReadString("Init:showChangedParticleData=false");
    //
    // specify how many Py8 events to print out, at either level
    // in this particular case print out a maximum of 10 events
    //
    fPythia8->ReadString("Next:numberShowProcess = 0");
    fPythia8->ReadString("Next:numberShowEvent = 10");
#endif

    fPythia8->Pythia8()->init();

    // shut off decays of pi0's as we want Geant4 to handle them
    // if other immediate decay products should be handled by Geant4,
    // their respective decay modes should be shut off as well
    //
    fPythia8->ReadString("111:onMode = off");

  }
#if 1

    fPythia8->Pythia8()->event.reset();

    // Create a Pythia8::Vec4 object to hold the particle's momentum and energy
    double px = tauPO.m_px;
    double py = tauPO.m_py;
    double pz = tauPO.m_pz;
    double e = tauPO.m_energy;
    double mass2 = e*e-px*px-py*py-pz*pz;
    double mass = sqrt(std::max(mass2,0.0));

    // Append the particle to the event
    fPythia8->Pythia8()->event.append(
        tauPO.m_pdg_id,   // PDG ID
        1,                // Status (assuming final state particle)
 //       0, 0,             // Mother indices
        0, 0,             // Color indices
        px, py, pz, e,    // Momentum components and energy
        mass,             // Mass
        0.0,              // scale (default)
        0.0               // polarization (set unpolarized) - TODO: set proper value
    );

#endif

  // TODO: // specify polarization!!

  int npart_before_decay = fPythia8->Pythia8()->event.size();
   
  fPythia8->Pythia8()->next();
   
  int npart_after_decay = fPythia8->Pythia8()->event.size();

 for ( int ip=npart_before_decay; ip<npart_after_decay; ++ip )
   {
            // only select final state decay products (direct or via subsequent decays=
      if ( fPythia8->Pythia8()->event[ip].status() < 0 ) continue;
      struct PO decayProdPO;
      decayProdPO.m_pdg_id = fPythia8->Pythia8()->event[ip].id();
      decayProdPO.m_px = fPythia8->Pythia8()->event[ip].px();
      decayProdPO.m_py = fPythia8->Pythia8()->event[ip].py();
      decayProdPO.m_pz = fPythia8->Pythia8()->event[ip].pz();
      decayProdPO.m_energy = fPythia8->Pythia8()->event[ip].e();
      taudecay.push_back(decayProdPO);
   }

	double decaylength = 0; // sqrt(aPO.m_vx_decay*aPO.m_vx_decay+aPO.m_vy_decay*aPO.m_vy_decay+aPO.m_vz_decay*aPO.m_vz_decay);
	tautracklength = decaylength;

}
#endif

size_t TPOEvent::n_charged() const {
  TDatabasePDG *pdgDB = TDatabasePDG::Instance();
  size_t nc=0;
  for (size_t i=0; i<n_particles(); i++) {
    struct PO aPO = POs[i];
    if(aPO.m_status != 1) continue;
    TParticlePDG *particle = pdgDB->GetParticle(aPO.m_pdg_id);
    if(particle != nullptr && particle->Charge() == 0) nc++;
  }
  return nc;
}

void TPOEvent::kinematics_event() {
  bool got_out_lepton = false;
  spx=spy=spz=0;
  tauvis_px=tauvis_py=tauvis_pz=0;
  for (size_t i=0; i<n_particles(); i++) {
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

  // additiona processing in case of tau decay
  if(istau && isCC){
    TDatabasePDG *pdgDB = TDatabasePDG::Instance();
    int nc = 0, nn = 0;
    for(int i=0; i<n_taudecay(); i++) {
      struct PO aPO = taudecay[i];
      TParticlePDG *particle = pdgDB->GetParticle(aPO.m_pdg_id);
      if(aPO.m_status == 1 && !is_neutrino(aPO.m_pdg_id)){
      	if(abs(aPO.m_pdg_id) == 11) {
	        tau_decaymode = 1;
	      } else if(abs(aPO.m_pdg_id) == 13) {
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

const char* TPOEvent::reaction_desc() const {
  if(isCC) {
    switch(in_neutrino.m_pdg_id) {
      case 12:
        return "nueCC";
      case 14:
        return "numuCC";
      case 16:
        return "nutauCC";
      case -12:
        return "antinueCC";
      case -14:
        return "antinumuCC";
      case -16:
        return "antinutauCC";
    }
  } else {
    return "nuNC";
  }
  return " ?? ";
}

void TPOEvent::dump_header() const {
  const char *reaction = reaction_desc();
  std::cout << "--- Run " << run_number << " Event " << event_id << " ---------------------------------------- " << reaction << " --------------------------------------------" << std::endl;
}
 
void TPOEvent::dump_event() const {
  TDatabasePDG *pdgDB = TDatabasePDG::Instance();
  dump_header();
  std::cout << " Primary vtx = " << prim_vx.x() << " " << prim_vx.y() << " " << prim_vx.z() << " mm ";
  if(vtx_target>0) {
    switch(vtx_target) {
      case kVtx_in_W: 
        std::cout << " - in W target";
        break;
      case kVtx_in_Scint:
        std::cout << " - in Scint target";
        break;
    }
  } 
  std::cout << std::endl;
  std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "¨    trackID, pdg_ID, name, px, py, pz, E, status, geant4ID, parents" << std::endl;
  for (size_t i=0; i<n_particles(); i++) {
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
  if(n_taudecay()>0) {
    std::cout << "Tau decay mode : " << tau_decaymode << std::endl;
    std::cout << "¨    trackID, pdg_ID, name, px, py, pz, E, status, geant4ID, parents" << std::endl;
    for (size_t i=0; i<n_taudecay(); i++) {
      struct PO aPO = taudecay[i];
      dump_PO(aPO, pdgDB);
    }
    std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
  }
}

int TPOEvent::findFromGEANT4TrackID(int trackID) {
  for (size_t i=0; i<n_particles(); i++) {
    if(POs[i].geanttrackID == trackID){
      return i;
    }
  }
  return -1;
}




