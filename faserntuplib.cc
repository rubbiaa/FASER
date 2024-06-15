// FASER ntuple interface library
// A. Rubbia May 2024
//

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <map>

#include <cmath>
#include <random>

#include "faserntuplib.hh"

void clear_event() {
  event.run_number = event.event_id = -1;
  event.prim_vx[0] = event.prim_vx[1] = event.prim_vx[2] = 0;
  event.n_particles = event.n_charged = 0;
  event.n_taudecay = 0;
  event.tau_decaymode = -1;
  event.isCC = false;
  event.isES = false;
  event.istau = false;
  event.spx=event.spy=event.spz=0;
  event.tauvis_px=event.tauvis_py=event.tauvis_pz=0;
  event.Evis = event.ptmiss = -1;
  event.cost = event.cosf = -999;
  event.taupi_cand[0] = event.taupi_cand[1] = event.taupi_cand[2] = 0;
  event.taurho_cand[0] = event.taurho_cand[1] = event.taurho_cand[2] = 0;
};

bool is_lepton(int pdgid) {
  int pdgidabs = abs(pdgid);
  return (pdgidabs >= 11 && pdgidabs <= 16);
}

double compute_qtlep(double p[3], double jet[3]) {
  // projection of momentum on hadronic jet
  double pproj = (p[0]*jet[0] + p[1]*jet[1] + p[2]*jet[2])/sqrt(jet[0]*jet[0]+jet[1]*jet[1]+jet[2]*jet[2]);
  double p2 = p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
  return sqrt(p2-pproj*pproj);
}

bool is_neutrino(int pdgid) {
  int pdgidabs = abs(pdgid);
  return (pdgidabs == 12 || pdgidabs == 14 || pdgidabs == 16);
}

void smear_PO(struct PO *aPO, TDatabasePDG *pdgDB) {
  // skip neutrinos
  if(is_neutrino(aPO->m_pdg_id)) {
    return;
  }

  TParticlePDG *particle = pdgDB->GetParticle(aPO->m_pdg_id);
  //  dump_PO(*aPO, pdgDB);
  // particle not found? probably heavy nucleus recoil.. ignore
  if(particle == nullptr) {
    aPO->m_px = 0;
    aPO->m_py = 0;
    aPO->m_pz = 0;
    return;
  }
  double m = particle->Mass();
  double px = aPO->m_px;
  double py = aPO->m_py;
  double pz = aPO->m_pz;
  double p2 = px * px + py * py + pz * pz;
  double p = sqrt(p2);
  double e = sqrt(p2+m*m);

  // default smearing factors
  double eres_stoch = 0.5;
  double eres_const = 0.1;
  double angle_res = 1.0; // in degrees

  // Define smearing factors based on particle ID
  switch (abs(aPO->m_pdg_id)) {
  case 11: // Electron
  case 22: // Gamma
  case 111: // pi0
    eres_stoch = 0.10;
    eres_const = 0.05;
    break;
  case 13: // Muon
    eres_stoch = 0;
    eres_const = 0.2;
    break;
  default:
    // Keep default values
    break;
  }
  
  double eneres = sqrt(pow(eres_stoch/sqrt(e),2)+pow(eres_const,2));
  std::normal_distribution<> d(0, eneres); // Normal distribution with mean 0 and std deviation eneres
  
  // Apply Gaussian smearing to energy
  // Set up random number generator
  std::random_device rd;  // Seed generator
  std::mt19937 gen(rd()); // Mersenne Twister engine
  double smeared_e = e + e * d(gen);

  // Angular smearing
  double degree_to_radian = M_PI / 180.0;
  double angle_stddev = angle_res * degree_to_radian; // 1 degree in radians
  std::normal_distribution<> angle_d(0, angle_stddev); // Normal distribution for angles
  
  // Generate random angles for smearing
  double delta_theta = angle_d(gen);
  double delta_phi = angle_d(gen);
  
  // Current direction in spherical coordinates
  
  double theta = acos(pz / p);
  double phi = atan2(py, px);
  
  // Smear angles
  theta += delta_theta;
  phi += delta_phi;

  double smeared_p = smeared_e*smeared_e - m*m;
  if(smeared_p>0) {
    smeared_p = sqrt(smeared_p);
  } else {
    smeared_p = 0.1;
  }
  
  // Convert back to Cartesian coordinates
  aPO->m_px = smeared_p * sin(theta) * cos(phi);
  aPO->m_py = smeared_p * sin(theta) * sin(phi);
  aPO->m_pz = smeared_p * cos(theta);

  //  std::cout << "After smearing" << std::endl;
  //  dump_PO(*aPO, pdgDB);
}

void smear_event() {
  TDatabasePDG *pdgDB = TDatabasePDG::Instance();
  for (size_t i=0; i<event.n_particles; i++) {
    struct PO *aPO = &event.POs[i];
    if(aPO->m_status == 1) {
      smear_PO(aPO, pdgDB);
    };
  }
}

void kinematics_event() {
  TDatabasePDG *pdgDB = TDatabasePDG::Instance();
  bool got_out_lepton = false;
  for (size_t i=0; i<event.n_particles; i++) {
    struct PO aPO = event.POs[i];
    if(aPO.m_status == 4 && i==0) {
      event.in_neutrino = aPO;
      event.istau = (abs(aPO.m_pdg_id) == 16);
    }
    if(!got_out_lepton && aPO.m_status == 1 && is_lepton(aPO.m_pdg_id) ) {
      event.out_lepton = aPO;
      got_out_lepton = true;
    }
    // should add all the final state particles, except neutrinos and taus
    if(aPO.m_status == 1 && !is_neutrino(aPO.m_pdg_id) && abs(aPO.m_pdg_id) != 15) {
      TParticlePDG *particle = pdgDB->GetParticle(aPO.m_pdg_id);
      if(particle != nullptr && particle->Charge()!=0) {
	event.n_charged++;
      }
      event.spx += aPO.m_px;
      event.spy += aPO.m_py;
      event.spz += aPO.m_pz;
    }
  }
  event.isCC = !(event.in_neutrino.m_pdg_id == event.out_lepton.m_pdg_id);

  if(event.isCC && !event.istau) {
    event.jetpx = event.spx-event.out_lepton.m_px;
    event.jetpy = event.spy-event.out_lepton.m_py;
    event.jetpz = event.spz-event.out_lepton.m_pz;
  } else {
    event.jetpx = event.spx;
    event.jetpy = event.spy;
    event.jetpz = event.spz;
  }
  if(event.jetpx == 0 && event.jetpy == 0 && event.jetpz == 0){
    event.isES = true;
    std::cout << "jet momentum = 0" << std::endl;
    dump_event();
  }
  event.vis_spx = event.spx;
  event.vis_spy = event.spy;
  event.vis_spz = event.spz;

  // check decay mode of the tau
  if(event.istau && event.isCC){
    int nc = 0, nn = 0;
    for(int i=0; i<event.n_taudecay; i++) {
      struct PO aPO = event.taudecay[i];
      TParticlePDG *particle = pdgDB->GetParticle(aPO.m_pdg_id);
      if(aPO.m_status == 1 && !is_neutrino(aPO.m_pdg_id)){
	if(abs(aPO.m_pdg_id) == 11) {
	  event.tau_decaymode = 1;
	}
	if(abs(aPO.m_pdg_id) == 13) {
	  event.tau_decaymode = 2;
	}
	if(particle->Charge() == 0){
	  nn++;
	} else {
	  nc++;
	}
	event.tauvis_px += aPO.m_px;
	event.tauvis_py += aPO.m_py;
	event.tauvis_pz += aPO.m_pz;
      }
    }
    if(event.tau_decaymode < 0){
      event.tau_decaymode = 6;
      if(nc==1 && nn == 0) event.tau_decaymode = 3;
      if(nc==1 && nn == 1) event.tau_decaymode = 4;
      if(nc==3) event.tau_decaymode = 5;
    }
    event.vis_spx = event.tauvis_px + event.jetpx;
    event.vis_spy = event.tauvis_py + event.jetpy;
    event.vis_spz = event.tauvis_pz + event.jetpz;
    // tau -> pi decay mode - save candidate
    if(event.tau_decaymode == 3) {
      event.taupi_cand[0] = event.tauvis_px;
      event.taupi_cand[1] = event.tauvis_py;
      event.taupi_cand[2] = event.tauvis_pz;
    }
    // tau -> rho decay mode - save candidate
    if(event.tau_decaymode == 4) {
      event.taurho_cand[0] = event.tauvis_px;
      event.taurho_cand[1] = event.tauvis_py;
      event.taurho_cand[2] = event.tauvis_pz;
    }
  }
  event.Evis = sqrt(event.vis_spx*event.vis_spx + event.vis_spy*event.vis_spy + event.vis_spz*event.vis_spz);
  event.ptmiss = sqrt(event.vis_spx*event.vis_spx + event.vis_spy*event.vis_spy);

  // extra kinematics
  if(abs(event.in_neutrino.m_pdg_id) == 16) {
    event.cost = (event.tauvis_px*event.jetpx + event.tauvis_py*event.jetpy)/(sqrt(event.tauvis_px*event.tauvis_px+event.tauvis_py*event.tauvis_py)*sqrt(event.jetpx*event.jetpx+event.jetpy*event.jetpy));
    double ptmissx = -(event.tauvis_px+event.jetpx);
    double ptmissy = -(event.tauvis_py+event.jetpy);
    event.cosf = (event.tauvis_px*ptmissx + event.tauvis_py*ptmissy)/(sqrt(event.tauvis_px*event.tauvis_px+event.tauvis_py*event.tauvis_py)*sqrt(ptmissx*ptmissx+ptmissy*ptmissy));
  } else {
    double lep_px = event.vis_spx-event.jetpx;
    double lep_py = event.vis_spy-event.jetpy;
    event.cost = (lep_px*event.jetpx + lep_py*event.jetpy)/(sqrt(lep_px*lep_px+lep_py*lep_py)*sqrt(event.jetpx*event.jetpx+event.jetpy*event.jetpy));
    double ptmissx = -event.vis_spx;
    double ptmissy = -event.vis_spy;
    event.cosf = (lep_px*ptmissx + lep_py*ptmissy)/(sqrt(lep_px*lep_px+lep_py*lep_py)*sqrt(ptmissx*ptmissx+ptmissy*ptmissy));
  }

  // for NC events, search for the tau->pi candidate
  if(!event.isCC && !event.isES) {
    double jet[3];
    jet[0] = event.jetpx;
    jet[1] = event.jetpy;
    jet[2] = event.jetpz;
    double qtlepmax = 0;
    for (size_t i=0; i<event.n_particles; i++) {
      struct PO aPO = event.POs[i];
      // should consider all charged final state particles, except neutrinos and taus
      if(aPO.m_status == 1 && !is_neutrino(aPO.m_pdg_id) && abs(aPO.m_pdg_id) != 15) {
	TParticlePDG *particle = pdgDB->GetParticle(aPO.m_pdg_id);
	if(particle != nullptr && particle->Charge()!=0) {
	  double p[3];
	  p[0] = aPO.m_px;
	  p[1] = aPO.m_py;
	  p[2] = aPO.m_pz;
	  double qtlep = compute_qtlep(p, jet);
	  if(qtlep > qtlepmax) {
	    event.taupi_cand[0] = aPO.m_px;
	    event.taupi_cand[1] = aPO.m_py;
	    event.taupi_cand[2] = aPO.m_pz;
	    qtlepmax = qtlep;
	  }
	}	    
      }
    }
  }

  // for NC events, search for the tau->rho candidate
  if(!event.isCC && !event.isES) {
    double jet[3];
    jet[0] = event.jetpx;
    jet[1] = event.jetpy;
    jet[2] = event.jetpz;
    double qtlepmax = 0;
    for (size_t i=0; i<event.n_particles; i++) {
      struct PO aPO = event.POs[i];
      // should consider all charged final state particles, except neutrinos and taus
      if(aPO.m_status == 1 && !is_neutrino(aPO.m_pdg_id) && abs(aPO.m_pdg_id) != 15) {
	TParticlePDG *particle = pdgDB->GetParticle(aPO.m_pdg_id);
	if(particle != nullptr && particle->Charge()!=0) {
	  double p1[3];
	  p1[0] = aPO.m_px;
	  p1[1] = aPO.m_py;
	  p1[2] = aPO.m_pz;
	  for (size_t j=0; j<event.n_particles; j++) {
	    struct PO aPO2 = event.POs[j];
	    // should consider only pi0
	    if(aPO2.m_status == 1 && aPO2.m_pdg_id == 111) {
	      double p2[3];
	      p2[0] = aPO2.m_px;
	      p2[1] = aPO2.m_py;
	      p2[2] = aPO2.m_pz;
	      double ptot[3];
	      ptot[0] = p1[0]+p2[0];
	      ptot[1] = p1[1]+p2[1];
	      ptot[2] = p1[2]+p2[2];
	      double qtlep = compute_qtlep(ptot, jet);
	      if(qtlep > qtlepmax) {
		event.taurho_cand[0] = ptot[0];
		event.taurho_cand[1] = ptot[1];
		event.taurho_cand[2] = ptot[2];
		qtlepmax = qtlep;
	      }
	    }
	  } 
	}
      }
    }
  }
}


void dump_PO(struct PO aPO,  TDatabasePDG *pdgDB) {
  TParticlePDG *particle = pdgDB->GetParticle(aPO.m_pdg_id);
  std::cout << std::setw(10) << aPO.m_track_id << " " << aPO.m_pdg_id << " ";
  if(particle) {
    std::cout << std::setw(10) << particle->GetName() << " ";
  } else {
    std::cout << std::setw(10) << " ?? ";
  }
  double decaylength = sqrt(aPO.m_vx_decay*aPO.m_vx_decay+aPO.m_vy_decay*aPO.m_vy_decay+aPO.m_vz_decay*aPO.m_vz_decay);
  std::cout << std::setw(10) << "gct=" << decaylength << " ";
  std::cout << "vx=" << aPO.m_vx_decay << " vy = " << aPO.m_vy_decay << " vz = " << aPO.m_vz_decay << " ";
  std::cout << std::setw(10) << aPO.m_px << " " << " " << aPO.m_py << " " << aPO.m_pz << " " << aPO.m_status << " ";
  for (size_t j=0; j<aPO.nparent; j++) {
    std::cout << std::setw(10) << aPO.m_trackid_in_particle[j] << " ";
  }
  std::cout << std::endl;
}

void dump_event() {
  double spx=0, spy=0, spz=0;
  TDatabasePDG *pdgDB = TDatabasePDG::Instance();
  if(event.isCC) {
    std::cout << "--- Run " << event.run_number << " Event " << event.event_id << " ---------------------------------------- CC --------------------------------------------" << std::endl;
  } else {
    std::cout << "--- Run " << event.run_number << " Event " << event.event_id << " ------------------------------------------- NC ------------------------------------------" << std::endl;
  }
  std::cout << " Primary vtx = " << event.prim_vx[0] << " " << event.prim_vx[1] << " " << event.prim_vx[2] << std::endl;
  for (size_t i=0; i<event.n_particles; i++) {
    struct PO aPO = event.POs[i];
    dump_PO(aPO, pdgDB);
  }
  std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
  std::cout << std::setw(10) << "Outgoing lepton:           ";
  dump_PO(event.out_lepton, pdgDB);
  std::cout << std::setw(10) << "Jet :                      " << event.jetpx << " " << event.jetpy << " " << event.jetpz << std::endl;
  std::cout << std::setw(10) << "Sum final state particles: " << event.spx << " " << event.spy << " " << event.spz << std::endl;
  std::cout << std::setw(10) << "Sum final state particles (VIS): " << event.vis_spx << " " << event.vis_spy << " " << event.vis_spz << std::endl;
  std::cout << std::setw(10) << "Ptmiss = " << event.ptmiss << "  Evis = " << event.Evis << std::endl;
  std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
  if(event.n_taudecay>0) {
    std::cout << "Tau decay mode : " << event.tau_decaymode << std::endl;
    for (size_t i=0; i<event.n_taudecay; i++) {
      struct PO aPO = event.taudecay[i];
      dump_PO(aPO, pdgDB);
    }
    std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
  }
  std::cout << std::setw(10) << "   Tau pi cand " << event.taupi_cand[0] << " " << event.taupi_cand[1] << " " << event.taupi_cand[2] << std::endl;
    std::cout << std::setw(10) << "   Tau rho cand " << event.taurho_cand[0] << " " << event.taurho_cand[1] << " " << event.taurho_cand[2] << std::endl;
  std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
}
