#ifndef _TPARTICLEGUN_H
#define _TPARTICLEGUN_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TObject.h"

#include <TTree.h>

#include "TPORecoEvent.hh"
#include "TcalEvent.hh"

class TParticleGun : public TObject {

public:
struct FEATURES {
    Int_t m_pdg_id;        // truth particle id
    Float_t m_energy;      // truth energy
    Float_t m_jet_energy;  // truth jet energy
    // reco eflow
    Float_t ef_evis;    // visible energy
    Float_t ef_et;      // visible transverse energy
    Float_t ef_fasercal_x; // energy in FASER calorimeter (compensated)
    Float_t ef_fasercal_y;
    Float_t ef_fasercal_z;
    Float_t ef_ecal_x;    // calorimetric energy
    Float_t ef_ecal_y;    // calorimetric energy
    Float_t ef_ecal_z;    // calorimetric energy
    Float_t ef_hcal_x;    // hadronic calorimetric energy
    Float_t ef_hcal_y;
    Float_t ef_hcal_z;
    // shower profile parameters
    Float_t ep_chi2_per_ndf; // longitudinal energy profile fit: fitted chi2/ndf
    Float_t ep_E0;         
    Float_t ep_a;
    Float_t ep_b;
    Float_t ep_tmax;
    Float_t ep_c;
};

private:

public:
    ClassDef(TParticleGun, 1)

    struct FEATURES features;

    void Create_Sel_Tree(TTree *t);
    int ProcessEvent(TcalEvent *fTcalEvent, TPORecoEvent *fPORecoEvent);
    void Fill_Sel_Tree(TTree *t) { t->Fill();};

};


#endif
