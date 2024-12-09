#include "TParticleGun.hh"

#include <TTree.h>

ClassImp(TParticleGun);

void TParticleGun::Create_Sel_Tree(TTree *t) {
  t->Branch("m_pdg_id",&features.m_pdg_id);
  t->Branch("m_energy",&features.m_energy);
  t->Branch("ep_chi2_per_ndf", &features.ep_chi2_per_ndf);
  t->Branch("ep_E0", &features.ep_E0);
  t->Branch("ep_a", &features.ep_a);
  t->Branch("ep_b", &features.ep_b);
  t->Branch("ep_tmax", &features.ep_tmax);
  t->Branch("ep_c", &features.ep_c);
};
