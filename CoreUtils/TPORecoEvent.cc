#include "TPORecoEvent.hh"

ClassImp(TPORec);
ClassImp(TPORecoEvent);

TPORecoEvent::TPORecoEvent(TcalEvent* c, TPOEvent* p) : fTcalEvent(c), fTPOEvent(p) { 
};

TPORecoEvent::~TPORecoEvent() {
    for(auto it : fPORecs) {
        delete it;
    }
    delete fPOFullEvent;
    fPORecs.clear();
}


void TPORecoEvent::Reconstruct() {
    fPORecs.clear();

    std::cout << "Starting reconstruction... " << fTcalEvent->fTracks.size() << " tracks to process" << std::endl;
    // loop over all digitized tracks and create primaries
    for (auto it : fTcalEvent->fTracks) {
        if(it -> fparentID ==0) {
            int POID = fTPOEvent->findFromGEANT4TrackID(it->ftrackID);

            // skip final state neutrinos
            if(fTPOEvent->is_neutrino(fTPOEvent->POs[POID].m_pdg_id)) continue;

            TPORec* aPORec = new TPORec(POID);
            aPORec->POID = POID;
            aPORec->fGEANTTrackIDs.push_back(it->ftrackID);
            aPORec->DTs.push_back(it);
            struct TPORec::CALENERGIES calene = computeEnergiesAndCOG(it);
            aPORec->fEnergiesCogs.push_back(calene);
            fPORecs.push_back(aPORec);
        }
    }

    // now loop over secondaries and add to the primary
    for (auto it : fTcalEvent->fTracks) {
        if(it -> fparentID !=0){
            // find primary
            int primaryID = it->fprimaryID;
            if(primaryID>-1){
                for (auto itRecs : fPORecs) {
                    if(itRecs->fGEANTTrackIDs[0] == primaryID) {
                        itRecs->fGEANTTrackIDs.push_back(it->ftrackID);
                        itRecs->DTs.push_back(it);
                        struct TPORec::CALENERGIES calene = computeEnergiesAndCOG(it);
                        itRecs->fEnergiesCogs.push_back(calene);
                        break;
                    }
                }
            } else {
                // this is a hanging secondary - it means that some intermediate particle did 
                // deposit energy anywhere
                std::cout << " This shouldn't happen...." << std::endl;
                exit(1);
            }
        }
    }

    // now sum all quantities belowing to a given primary applying compensation
    double alpha = 1.0/(1.0-0.341)*0.98;
    double beta = 3.0;
    for(auto it : fPORecs) {
        int ntracks = it->fGEANTTrackIDs.size();
        it->fTotal.em = 0;
        it->fTotal.had = 0;
        it->fTotal.cog.SetCoordinates(0,0,0);
        for (int i = 0; i < ntracks; i++) {
            it->fTotal.em += it->fEnergiesCogs[i].em;
            it->fTotal.had += it->fEnergiesCogs[i].had;
            it->fTotal.cog += it->fEnergiesCogs[i].cog*(it->fEnergiesCogs[i].em+it->fEnergiesCogs[i].had);
        }
        double Eraw = it->fTotal.em+it->fTotal.had;
        if(Eraw > 0) {
            it->fTotal.cog /= Eraw;
        }
        // apply compensation
        it->fTotal.Ecompensated = it->fTotal.em*alpha+it->fTotal.had*beta;

        // if the energy is less than 2 GeV then we should go for the integration of dE/dx
        double Threshold_for_dEdx = 2.0;
        if(it->fTotal.Ecompensated>0 && it->fTotal.Ecompensated < Threshold_for_dEdx) {
            double EKin = 0;
            for (int i = 0; i < ntracks; i++) {
                EKin += it->fEnergiesCogs[i].em + it->fEnergiesCogs[i].had;               
            }
            // get the mass from the PO information - assumes perfect particle ID!
            double px = fTPOEvent->POs[it->POID].m_px;
            double py = fTPOEvent->POs[it->POID].m_py;
            double pz = fTPOEvent->POs[it->POID].m_pz;
            double E = fTPOEvent->POs[it->POID].m_energy;
            double mass = sqrt(E*E-px*px-py*py-pz*pz);
            double Ereco = EKin + mass;
            it->fTotal.Ecompensated = Ereco;
        }

        // now compute the energy compensated eflow relative to primary vertex
        ROOT::Math::XYZVector primary(fTPOEvent->prim_vx.X(), fTPOEvent->prim_vx.Y(), fTPOEvent->prim_vx.Z());
        ROOT::Math::XYZVector direction = it->fTotal.cog - primary;
        it->fTotal.Eflow = it->fTotal.Ecompensated *direction.Unit();

    }

    // Now compute full event quantities
    fPOFullEvent = new TPORec(-1);
    fPOFullEvent->fTotal.Ecompensated = 0;
    fPOFullEvent->fTotal.cog.SetCoordinates(0,0,0);
    for(auto it : fPORecs) {
        fPOFullEvent->fTotal.Ecompensated += it->fTotal.Ecompensated;
        fPOFullEvent->fTotal.cog += it->fTotal.cog*it->fTotal.Ecompensated;
    }
    fPOFullEvent->fTotal.cog /= fPOFullEvent->fTotal.Ecompensated;
    // now compute the energy compensated eflow relative to primary vertex
    ROOT::Math::XYZVector primary(fTPOEvent->prim_vx.X(), fTPOEvent->prim_vx.Y(), fTPOEvent->prim_vx.Z());
    ROOT::Math::XYZVector direction = fPOFullEvent->fTotal.cog - primary;
    fPOFullEvent->fTotal.Eflow = fPOFullEvent->fTotal.Ecompensated * direction.Unit();
}

struct TPORec::CALENERGIES TPORecoEvent::computeEnergiesAndCOG(DigitizedTrack *dt) {
    struct TPORec::CALENERGIES result;

    result.em = result.had = 0;
    result.Ecompensated = 0;
    double cogx = 0;
    double cogy = 0;
    double cogz = 0;

    bool isEM = (abs(dt->fPDG) == 11);

    size_t nhits = dt->fhitIDs.size();
    for ( size_t i = 0; i < nhits; i++) {
        ROOT::Math::XYZVector position = fTcalEvent -> getChannelXYZfromID(dt->fhitIDs[i]);
        double ehit = dt->fEnergyDeposits[i]/1e3;
        if(isEM) {
            result.em += ehit;
        } else {
            result.had += ehit;
        }
        cogx += position.X()*ehit;
        cogy += position.Y()*ehit;
        cogz += position.Z()*ehit;
    }
    double etot = result.em + result.had;
    if(etot>0) {
        result.cog.SetX(cogx/etot);
        result.cog.SetY(cogy/etot);
        result.cog.SetZ(cogz/etot);

        // compute energy flow relative to primary vertex
        ROOT::Math::XYZVector primary(fTPOEvent->prim_vx.X(), fTPOEvent->prim_vx.Y(), fTPOEvent->prim_vx.Z());
        ROOT::Math::XYZVector direction = result.cog - primary;
        result.Eflow = etot*direction.Unit();

    } else {
        result.cog.SetCoordinates(0,0,0);   
        result.Eflow.SetCoordinates(0,0,0);
    }

    return result;
}

void TPORecoEvent::Dump() {
    TDatabasePDG *pdgDB = TDatabasePDG::Instance();
    fTPOEvent->dump_header();
    for(auto it : fPORecs) {
        fTPOEvent->dump_PO(fTPOEvent->POs[it->POID], pdgDB);
        int ntracks = it->fGEANTTrackIDs.size();
        std::cout << "RECO>>" << std::setw(10) << ntracks << " tracks ";
        std::cout << "em: " << std::setw(10) << it->fTotal.em << " had: " << std::setw(10) << it->fTotal.had << " ";
        std::cout << "comp: " << std::setw(10) << it->fTotal.Ecompensated << " ";
        std::cout << std::endl;
        std::cout << "EFLOW>> ";
        std::cout << " Ex: " << it->fTotal.Eflow.X();
        std::cout << " Ey: " << it->fTotal.Eflow.Y();
        std::cout << " Ez: " << it->fTotal.Eflow.Z();
        std::cout << "      (COG: " << it->fTotal.cog.X() << " " << it->fTotal.cog.Y() << " " << it->fTotal.cog.Z() << ")";
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "FULL EFLOW>> ";
    std::cout << " Ex: " << fPOFullEvent->fTotal.Eflow.X();
    std::cout << " Ey: " << fPOFullEvent->fTotal.Eflow.Y();
    std::cout << " Ez: " << fPOFullEvent->fTotal.Eflow.Z();
//    std::cout << " cogx: " << fPOFullEvent->fTotal.cog.X() << " ";
//    std::cout << " cogy: " << fPOFullEvent->fTotal.cog.Y() << " ";
//    std::cout << " cogz: " << fPOFullEvent->fTotal.cog.Z() << " ";
    std::cout << std::endl;
    std::cout << "FULL EVENT>> Ene: " << fPOFullEvent->fTotal.Ecompensated << " ";
    double ETmiss = sqrt(fPOFullEvent->fTotal.Eflow.Perp2());
    std::cout << " ETmiss: " << ETmiss;
    std::cout << std::endl;
    fTPOEvent->dump_header();
}