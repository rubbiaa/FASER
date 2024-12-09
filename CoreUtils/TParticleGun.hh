#ifndef _TPARTICLEGUN_H
#define _TPARTICLEGUN_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TObject.h"

#include <TTree.h>

class TParticleGun : public TObject {

public:
struct FEATURES {
    Int_t m_pdg_id;        // truth particle id
    Float_t m_energy;      // truth energy
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
    void Fill_Sel_Tree(TTree *t) { t->Fill();};

};


#endif
