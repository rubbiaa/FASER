#ifndef _TTAUSEARCH_H
#define _TTAUSEARCH_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <Math/Vector3D.h>

#include "TObject.h"
#include "TVector3.h"

#include <TTree.h>

#include "TcalEvent.hh"
#include "TPORecoEvent.hh"

class TTauSearch : public TObject {

public:
struct KINEMATICS {
    Float_t Evis;
    Float_t ptmiss;
    Float_t plep;
    Float_t ptlep;
    Float_t qtlep;
    Float_t cost;
    Float_t cosf;
    Int_t primary_n_charged;
};

private:
    ROOT::Math::XYZVector lepton;
    ROOT::Math::XYZVector jet;

public:
    ClassDef(TTauSearch, 1)

    struct KINEMATICS kine;

    void Create_Sel_Tree(TTree *t);
    void Fill_Sel_Tree(TTree *t) { t->Fill();};

    void SetLepton(ROOT::Math::XYZVector newlep) { lepton = newlep;};
    void SetJet(ROOT::Math::XYZVector newjet) { jet = newjet;};

    /// @brief Compute the tau search kinematical variables
    void Kinematics(int primary_n_charged);

    int ProcessEvent(TPORecoEvent *fPORecoEvent);
};


#endif
