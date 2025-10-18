#include "TParticleGun.hh"

#include <TTree.h>

ClassImp(TParticleGun);

void TParticleGun::Create_Sel_Tree(TTree *t) {
  t->Branch("m_pdg_id",&features.m_pdg_id);
  t->Branch("m_energy",&features.m_energy);
  t->Branch("m_jet_energy",&features.m_jet_energy);
  t->Branch("ef_evis", &features.ef_evis);
  t->Branch("ef_et", &features.ef_et);
  t->Branch("ef_fasercal_x", &features.ef_fasercal_x);
  t->Branch("ef_fasercal_y", &features.ef_fasercal_y);
  t->Branch("ef_fasercal_z", &features.ef_fasercal_z);
  t->Branch("ef_ecal_x", &features.ef_ecal_x);
  t->Branch("ef_ecal_y", &features.ef_ecal_y);
  t->Branch("ef_ecal_z", &features.ef_ecal_z);
  t->Branch("ef_hcal_x", &features.ef_hcal_x);
  t->Branch("ef_hcal_y", &features.ef_hcal_y);
  t->Branch("ef_hcal_z", &features.ef_hcal_z);
  t->Branch("ep_chi2_per_ndf", &features.ep_chi2_per_ndf);
  t->Branch("ep_E0", &features.ep_E0);
  t->Branch("ep_a", &features.ep_a);
  t->Branch("ep_b", &features.ep_b);
  t->Branch("ep_tmax", &features.ep_tmax);
  t->Branch("ep_c", &features.ep_c);
};

int TParticleGun::ProcessEvent(TcalEvent *fTcalEvent, TPORecoEvent *fPORecoEvent)
{
  // fill features from recoEvent
  if ((fPORecoEvent->GetPORecs()).size() == 0)
    return 0;
  struct TPORec *aPORec = (fPORecoEvent->GetPORecs())[0];
  int POID = aPORec->POID;
  if (POID >= 0)
  {
    struct PO *aPO = &fTcalEvent->fTPOEvent->POs[POID];
    features.m_pdg_id = aPO->m_pdg_id;
    features.m_energy = aPO->m_energy;

    // fill features of most energetic reconstructed cluster
    if (fPORecoEvent->PSClustersX.size() > 0)
    {
      TPSCluster *c = &fPORecoEvent->PSClustersX[0]; // most energetic one
      features.ep_chi2_per_ndf = c->longenergyprofile.chi2_per_ndf;
      features.ep_E0 = c->longenergyprofile.E0;
      features.ep_a = c->longenergyprofile.a;
      features.ep_b = c->longenergyprofile.b;
      features.ep_tmax = c->longenergyprofile.tmax;
      features.ep_c = c->longenergyprofile.c;
    }
  }

  // compute true jet kinematics
  TVector3 totvec = TVector3(0, 0, 0);
  TPOEvent &fTPOEvent = *(fTcalEvent->fTPOEvent);
  for (int i = 0; i < fTPOEvent.n_particles(); ++i)
  {
    struct PO &aPO = fTPOEvent.POs[i];
    // check if parent of the particle is the neutrino
    if (aPO.nparent > 0)
    {
      int parent_trackid = aPO.m_trackid_in_particle[0];
      if (parent_trackid == 0)
      {
        continue;
      }
    }
    else
    {
      continue;
    }
    TVector3 pvec(aPO.m_px, aPO.m_py, aPO.m_pz);
    totvec += pvec;
  }
  std::cout << " Jet total p: " << totvec.X() << " " << totvec.Y() << " " << totvec.Z() << " ";
  std::cout << " Jet magnitude: " << totvec.Mag() << " Jet pt: " << totvec.Perp() << std::endl;
  features.m_jet_energy = totvec.Mag();

  // full event features
  features.ef_evis = fPORecoEvent->GetPOFullRecoEvent()->TotalEvis();
  features.ef_et = fPORecoEvent->GetPOFullRecoEvent()->TotalET();
  features.ef_fasercal_x = fPORecoEvent->GetPOFullRecoEvent()->fTotal_fasercal.Eflow.X();
  features.ef_fasercal_y = fPORecoEvent->GetPOFullRecoEvent()->fTotal_fasercal.Eflow.Y();
  features.ef_fasercal_z = fPORecoEvent->GetPOFullRecoEvent()->fTotal_fasercal.Eflow.Z();

  features.ef_ecal_x = fPORecoEvent->GetPOFullRecoEvent()->fTotal_ecal.Eflow.X();
  features.ef_ecal_y = fPORecoEvent->GetPOFullRecoEvent()->fTotal_ecal.Eflow.Y();
  features.ef_ecal_z = fPORecoEvent->GetPOFullRecoEvent()->fTotal_ecal.Eflow.Z();
  features.ef_hcal_x = fPORecoEvent->GetPOFullRecoEvent()->fTotal_hcal.Eflow.X();
  features.ef_hcal_y = fPORecoEvent->GetPOFullRecoEvent()->fTotal_hcal.Eflow.Y();
  features.ef_hcal_z = fPORecoEvent->GetPOFullRecoEvent()->fTotal_hcal.Eflow.Z();
  return 1;
}