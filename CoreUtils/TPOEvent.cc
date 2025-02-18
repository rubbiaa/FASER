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

static void initialize_pythia() {
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

    // same for K0S and K0L
    fPythia8->ReadString("310:onMode = off");
    fPythia8->ReadString("130:onMode = off");

  }
}
#endif

ClassImp(TPOEvent)
// ClassImp(PO)

void TPOEvent::clear_event() {
  run_number = event_id = -1;
  setPrimaryVtx(0,0,0);
  vtx_target = 0;
  POs.clear();
  taudecay.clear();
  charmdecay.clear();
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
  initialize_pythia();
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

}
#endif

bool TPOEvent::isCharmed() const {
  bool charmed = false;
  for (size_t i=0; i<n_particles(); i++) {
    if(charmed) break;
    struct PO aPO = POs[i];
    if(abs(aPO.m_pdg_id)/100 == 4) charmed = true;
    if(abs(aPO.m_pdg_id)/100 == 104) charmed = true;
    if(abs(aPO.m_pdg_id)/1000 == 4) charmed = true;
    if(abs(aPO.m_pdg_id)/1000 == 14) charmed = true;
  }
  return charmed;
}

#ifdef _INCLUDE_PYTHIA_
void TPOEvent::perform_charmhadron_decay(struct PO PO) {
//  charmdecay.clear();
  int pdg = abs(PO.m_pdg_id);
  int nq1 = (pdg/1000)%10;
  int nq2 = (pdg/100)%10;
  bool ischarm = (nq1 == 0 && nq2 == 4) || (nq1 == 4);
  if(!ischarm) {
    return;
  }

  initialize_pythia();
  fPythia8->Pythia8()->event.reset();

  // Create a Pythia8::Vec4 object to hold the particle's momentum and energy
  double px = PO.m_px;
  double py = PO.m_py;
  double pz = PO.m_pz;
  double e = PO.m_energy;
  double mass2 = e * e - px * px - py * py - pz * pz;
  double mass = sqrt(std::max(mass2, 0.0));

  // Append the particle to the event
  fPythia8->Pythia8()->event.append(
      PO.m_pdg_id,   // PDG ID
      1,             // Status (assuming final state particle)
                     //       0, 0,             // Mother indices
      0, 0,          // Color indices
      px, py, pz, e, // Momentum components and energy
      mass,          // Mass
      0.0,           // scale (default)
      0.0            // polarization (set unpolarized) - TODO: set proper value
  );

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
      decayProdPO.nparent = 1;
      decayProdPO.m_trackid_in_particle[0] = PO.m_track_id;
      decayProdPO.geanttrackID = -1;
      decayProdPO.m_status = 1;
      charmdecay.push_back(decayProdPO);
   }

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

double TPOEvent::tauDecaylength() {
  if(!istau || !isCC) return -1;
  double xdecay = taudecay[0].m_vx_decay;
  double ydecay = taudecay[0].m_vy_decay;
  double zdecay = taudecay[0].m_vz_decay;
  double xprim = prim_vx.x();
  double yprim = prim_vx.y();
  double zprim = prim_vx.z();
  return sqrt((xdecay-xprim)*(xdecay-xprim) + (ydecay-yprim)*(ydecay-yprim) + (zdecay-zprim)*(zdecay-zprim));
}

double TPOEvent::tauKinkAngle() {
  if(!istau) return -1;
  double tauvis_p = sqrt(tauvis_px*tauvis_px + tauvis_py*tauvis_py + tauvis_pz*tauvis_pz);
  double tau_p = sqrt(out_lepton.m_px*out_lepton.m_px + out_lepton.m_py*out_lepton.m_py + out_lepton.m_pz*out_lepton.m_pz);
  double cosangle = (tauvis_px*out_lepton.m_px + tauvis_py*out_lepton.m_py + tauvis_pz*out_lepton.m_pz)/(tauvis_p*tau_p);
  return acos(cosangle);
}

void TPOEvent::dump_PO(struct PO aPO,  TDatabasePDG *pdgDB, std::ostream& out) const {
  TParticlePDG *particle = pdgDB->GetParticle(aPO.m_pdg_id);
  out << std::setw(10) << aPO.m_track_id << " " << std::setw(12) << aPO.m_pdg_id << " ";
  if(particle) {
    out << std::setw(10) << particle->GetName() << " ";
  } else {
    out << "        ?? ";
  }
//  double decaylength = sqrt(aPO.m_vx_decay*aPO.m_vx_decay+aPO.m_vy_decay*aPO.m_vy_decay+aPO.m_vz_decay*aPO.m_vz_decay);
//  std::cout << std::setw(10) << decaylength << " ";
//  std::cout << std::setw(10) << aPO.m_vx_decay << " " << std::setw(10) << aPO.m_vy_decay << " " << std::setw(10) << aPO.m_vz_decay << " ";
  out << std::setw(12) << aPO.m_px << " " << " " << std::setw(12) << aPO.m_py << " " << std::setw(12) << aPO.m_pz 
          << " " << std::setw(12) << aPO.m_energy << " " << std::setw(4) << aPO.m_status << " " << std::setw(4) << aPO.geanttrackID << " ";
  for (size_t j=0; j<aPO.nparent; j++) {
    out << std::setw(10) << aPO.m_trackid_in_particle[j] << " ";
  }
  out << std::endl;
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

void TPOEvent::dump_header(std::ostream& out) const {
  const char *reaction = reaction_desc();
  out << "--- Run " << run_number << " Event " << event_id << " ---------------------------------------- " << reaction << " --------------------------------------------" << std::endl;
}
 
void TPOEvent::dump_event(std::ostream& out) const {
  TDatabasePDG *pdgDB = TDatabasePDG::Instance();
  dump_header(out);
  out << " Primary vtx = " << prim_vx.x() << " " << prim_vx.y() << " " << prim_vx.z() << " mm ";
  if(use_GENIE_vtx) {
    out << "(this is the true vertex to be used in G4)";
  }
  if(vtx_target>0) {
    switch(vtx_target) {
      case kVtx_in_W: 
        out << " - in W target";
        break;
      case kVtx_in_Scint:
        out << " - in Scint target";
        break;
    }
  }
  if(GENIE_vtx_name != "") {
    out << " Genie has chosen this nucleus target " << GENIE_vtx_name;
  }
  out << std::endl; 
  out << "--------------------------------------------------------------------------------------------" << std::endl;
  out << "¨    trackID, pdg_ID, name, px, py, pz, E, status, geant4ID, parents" << std::endl;
  for (size_t i=0; i<n_particles(); i++) {
    struct PO aPO = POs[i];
    dump_PO(aPO, pdgDB, out);
  }
  out << "--------------------------------------------------------------------------------------------" << std::endl;
  out << std::setw(10) << "Outgoing lepton: " << std::endl;
  dump_PO(out_lepton, pdgDB, out);
  out << "--------------------------------------------------------------------------------------------" << std::endl;
  out << std::setw(10) << "Jet :                      " << jetpx << " " << jetpy << " " << jetpz << std::endl;
  out << std::setw(10) << "Sum final state particles: " << spx << " " << spy << " " << spz << std::endl;
  out << std::setw(10) << "Sum final state particles (VIS): " << vis_spx << " " << vis_spy << " " << vis_spz << std::endl;
  out << std::setw(10) << "Ptmiss = " << ptmiss << "  Evis = " << Evis << std::endl;
  out << "--------------------------------------------------------------------------------------------" << std::endl;
  if(n_taudecay()>0) {
    out << "Tau decay mode : " << tau_decaymode << std::endl;
    out << "¨    trackID, pdg_ID, name, px, py, pz, E, status, geant4ID, parents" << std::endl;
    for (size_t i=0; i<n_taudecay(); i++) {
      struct PO aPO = taudecay[i];
      dump_PO(aPO, pdgDB, out);
    }
    out << "--------------------------------------------------------------------------------------------" << std::endl;
  }
  if(charmdecay.size()>0) {
    out << "Charm decays products : " << std::endl;
    out << "¨    trackID, pdg_ID, name, px, py, pz, E, status, geant4ID, parents" << std::endl;
    for (size_t i=0; i<charmdecay.size(); i++) {
      struct PO aPO = charmdecay[i];
      dump_PO(aPO, pdgDB, out);
    }
    out << "--------------------------------------------------------------------------------------------" << std::endl;
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

void TPOEvent::clone(TPOEvent *src) {
  run_number = src->run_number;
  event_id = src->event_id;
  isCC = src->isCC;
  istau = src->istau;
  in_neutrino = src->in_neutrino;
  out_lepton = src->out_lepton;
  jetpx = src->jetpx;
  jetpy = src->jetpy;
  jetpz = src->jetpz;
}

void TPOEvent::reset_stats() {
  stats.ES = stats.nueCC = stats.numuCC = stats.nutauCC = stats.NC = 0;
}

void TPOEvent::update_stats()
{
  if (isES())
  {
    stats.ES++;
  }
  else
  {
    if (isCC)
    {
      int in_lepton_pdgid = in_neutrino.m_pdg_id;
      if (abs(in_lepton_pdgid) == 12)
      {
        stats.nueCC++;
      }
      if (abs(in_lepton_pdgid) == 14)
      {
        stats.numuCC++;
      }
      if (abs(in_lepton_pdgid) == 16)
      {
        stats.nutauCC++;
      }
    } else {
      stats.NC++;
    }
  }
}

void TPOEvent::dump_stats()
{
  std::cout << " nueCC = " << stats.nueCC;
  std::cout << " numuCC = " << stats.numuCC;
  std::cout << " nutauCC = " << stats.nutauCC;
  std::cout << " NC = " << stats.NC << std::endl;
  std::cout << " ES = " << stats.ES << std::endl;
}
