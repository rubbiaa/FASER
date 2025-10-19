#include "TTauSearch.hh"

#include <TTree.h>

ClassImp(TTauSearch);

void TTauSearch::Create_Sel_Tree(TTree *t) {
//  t->Branch("infid",&sumevent.infid);
  t->Branch("evis",&kine.Evis);
  t->Branch("plep",&kine.plep);
  t->Branch("ptlep",&kine.ptlep);
  t->Branch("qtlep",&kine.qtlep);
  t->Branch("ptmiss",&kine.ptmiss);
  t->Branch("cost",&kine.cost);
  t->Branch("cosf",&kine.cosf);
//  t->Branch("ptfitchi2",&sumevent.ptFitChi2);
  t->Branch("primary_n_charged", &kine.primary_n_charged);
};

void TTauSearch::Kinematics(int primary_n_charged) {

    ROOT::Math::XYZVector spx = lepton+jet;

    ROOT::Math::XYZVector lepperp = lepton;
    lepperp.SetZ(0);

    kine.Evis = sqrt(spx.Mag2());
    kine.plep = sqrt(lepton.Mag2());
    kine.ptlep = sqrt(lepperp.Mag2());
    kine.ptmiss = sqrt(spx.Perp2());

    // projection of momentum on hadronic jet
    double pproj = (lepton.Dot(jet))/sqrt(jet.Mag2());
    double p2 = lepton.Mag2();
    double qtlep2 = p2-pproj*pproj;
    if(qtlep2 > 0) {
        kine.qtlep = sqrt(qtlep2);
    } else if(qtlep2 > 1e-3) {
        kine.qtlep = 0;
    } else {
 //       std::cerr << "Error in compute qtlep " << qtlep2 << std::endl;
        kine.qtlep = 0;
    }

    ROOT::Math::XYZVector jetperp = jet;
    jetperp.SetZ(0);
    ROOT::Math::XYZVector spxperp = spx;
    spxperp.SetZ(0);

    kine.cost = (lepperp.Dot(jetperp))/(sqrt(lepperp.Mag2()*jetperp.Mag2()));
    kine.cosf = -(lepperp.Dot(spxperp))/(sqrt(lepperp.Mag2()*spxperp.Mag2()));

    kine.primary_n_charged = primary_n_charged;
}

int TTauSearch::ProcessEvent(TPORecoEvent *fPORecoEvent)
{
    // fill kinematics from recoEvent
    if ((fPORecoEvent->GetPORecs()).size() == 0)
        return 0;

        fPORecoEvent->Dump();

        TPOEvent *POevent = fPORecoEvent->GetPOEvent();

        bool found_electron = false;
        bool found_tau_e = false;
        ROOT::Math::XYZVector spx = fPORecoEvent->GetPOFullEvent()->fTotal.Ecompensated*
                fPORecoEvent->GetPOFullEvent()->fTotal.Eflow.Unit();
        for (auto it : fPORecoEvent->GetPORecs())
        {
            int POID = it->POID;
            if(POID < 0) continue;
            struct PO *aPO = &POevent->POs[POID];
            int PDG = aPO->m_pdg_id;

            double RecoLeptonEne = it->fTotal.Ecompensated;
            ROOT::Math::XYZVector lepton = RecoLeptonEne*it->fTotal.Eflow.Unit();
            ROOT::Math::XYZVector jet = spx-lepton;
            SetLepton(lepton);
            SetJet(jet);
            Kinematics(fPORecoEvent->primary_n_charged);

            // nueCC background
            if(abs(PDG) == 11 && !found_electron) {
                found_electron = true;
                if(!POevent->isES())                     // skip ES events
                    return 1;
            }

            // nutau signal -> e
            if(abs(PDG) == 15 && !found_tau_e) {
                if(fPORecoEvent->GetPOEvent()->tau_decaymode==1) {
                    found_tau_e = true;
                    if(!POevent->isES())                     // skip ES events
                    return 1;
                }
            }
        }

    return 0;
}
