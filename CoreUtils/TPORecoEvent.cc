#include "TPORecoEvent.hh"
#include "DBScan.hh"

#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <TDecompSVD.h>
#include <TVector3.h>
#include <TCanvas.h>

#include <random>

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

    std::cout << "Starting reconstruction... " << fTcalEvent->getfTracks().size() << " G4 tracks to process" << std::endl;
    // loop over all digitized tracks and create primaries
    for (auto it : fTcalEvent->getfTracks()) {
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
    for (auto it : fTcalEvent->getfTracks()) {
        if(it -> fparentID !=0){
            // find primary
            int primaryID = it->fprimaryID;
            if(primaryID>-1){
                bool foundPORec = false;
                for (auto itRecs : fPORecs) {
                    if(itRecs->fGEANTTrackIDs[0] == primaryID) {
                        itRecs->fGEANTTrackIDs.push_back(it->ftrackID);
                        itRecs->DTs.push_back(it);
                        struct TPORec::CALENERGIES calene = computeEnergiesAndCOG(it);
                        itRecs->fEnergiesCogs.push_back(calene);
                        foundPORec = true;
                        break;
                    }
                }
                if (!foundPORec)
                { // hanging secondary - for example, if primary did not leave any
                  // signal in sensitive detectors (e.g. tau decaying within W)
                    // first insert a primary
                    int POID = fTPOEvent->findFromGEANT4TrackID(primaryID);
                    TPORec *aPORec = new TPORec(POID);
                    aPORec->POID = POID;
                    aPORec->fGEANTTrackIDs.push_back(primaryID);
                    DigitizedTrack *itprim = new DigitizedTrack();
                    itprim -> ftrackID = primaryID;
                    itprim -> fparentID = 0;
                    itprim -> fprimaryID = primaryID;
                    itprim -> fPDG = 0; // FIXME
                    aPORec->DTs.push_back(itprim);
                    struct TPORec::CALENERGIES calene = computeEnergiesAndCOG(itprim);
                    aPORec->fEnergiesCogs.push_back(calene);
                    // now add secondary
                    aPORec->fGEANTTrackIDs.push_back(it->ftrackID);
                    aPORec->DTs.push_back(it);
                    struct TPORec::CALENERGIES calene2 = computeEnergiesAndCOG(it);
                    aPORec->fEnergiesCogs.push_back(calene2);
                    fPORecs.push_back(aPORec);
                }
            }
            else
            {
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
        if(std::isnan(Eraw)) {
            std::cerr << "TPORecoEvent::Reconstruct - raw energy is nan!" << std::endl;
        }
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
            double mass = sqrt(std::max(E*E-px*px-py*py-pz*pz,0.0));

            double Ereco = EKin + mass;
            if(std::isnan(Ereco)) {
                std::cerr << "TPORecoEvent::Reconstruct - Ereco is nan!" << std::endl;
            }
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
    if(std::isnan(fPOFullEvent->fTotal.Ecompensated)) {
        std::cerr << "TPORecoEvent::Reconstruct - compensated total energy is nan!" << std::endl;
    }
    fPOFullEvent->fTotal.cog /= fPOFullEvent->fTotal.Ecompensated;
    // now compute the energy compensated eflow relative to primary vertex
    ROOT::Math::XYZVector primary(fTPOEvent->prim_vx.X(), fTPOEvent->prim_vx.Y(), fTPOEvent->prim_vx.Z());
    ROOT::Math::XYZVector direction = fPOFullEvent->fTotal.cog - primary;
    fPOFullEvent->fTotal.Eflow = fPOFullEvent->fTotal.Ecompensated * direction.Unit();

    // some additional event summary variables
    primary_n_charged = fTPOEvent->n_charged();
    nhits_tau = 0;
    nhits_tracker_first = 0;
    for (auto it : fPORecs)
    {
        int POID = it->POID;
        struct PO *aPO = &fTcalEvent->fTPOEvent->POs[POID];
        int PDG = aPO->m_pdg_id;
        // consider only primary track
        DigitizedTrack *dt = it->DTs[0];
        size_t nhits = dt->fEnergyDeposits.size();
        int hittype;
        long minlayer = 999999999;
        // first search for the minimal layer
        for (size_t i = 0; i < nhits; i++)
        {
            // loop over hits in the tracker
            hittype = fTcalEvent->getChannelTypefromID(dt->fhitIDs[i]);
            if (hittype != 1) continue;
            long layer = fTcalEvent->getChannelLayerfromID(dt->fhitIDs[i]);
            if(layer < minlayer) minlayer = layer;
        }
        for (size_t i = 0; i < nhits; i++)
        {
            // loop over hits in the tracker
            hittype = fTcalEvent->getChannelTypefromID(dt->fhitIDs[i]);
            if (hittype != 1) continue;
            long layer = fTcalEvent->getChannelLayerfromID(dt->fhitIDs[i]);
            if(abs(layer - minlayer) <=1) {
                nhits_tracker_first++;
            }
        }
        switch (abs(PDG))
        {
        case 15: // taus
            // compute number of hits in scintillator
            for (size_t i = 0; i < nhits; i++)
            {
                hittype = fTcalEvent->getChannelTypefromID(dt->fhitIDs[i]);
                if (hittype != 0)
                    continue;
                nhits_tau++;
            }
        }
    }
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
        double ehit = dt->fEnergyDeposits[i]/1e3;   // *CLHEP::MeV
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
    if(std::isnan(etot)) {
        std::cerr << "TPORecoEvent::computeEnergiesAndCOG - compensated energy is nan!" << std::endl;
    }
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

static TVector3 fitLineThroughPoints(const struct TPORec::TRACK &track, TVector3& centroid) {
    std::vector<TPORec::TRACKHIT> tkhit = track.tkhit;
    int N = tkhit.size();
    // Calculate the centroid of the points
    centroid.SetXYZ(0, 0, 0);
    for (const auto& hit : tkhit) {
        centroid += TVector3(hit.point.x(), hit.point.y(), hit.point.z());
    }
    centroid *= (1.0 / N);

    // Compute covariance matrix
    TMatrixDSym covariance(3);
    for (const auto& hit : tkhit) {
        TVector3 centered = TVector3(hit.point.x(), hit.point.y(), hit.point.z()) - centroid;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                covariance(i, j) += centered[i] * centered[j];
            }
        }
    }

 //   std::cout << " Cov: "; covariance.Print();

    // Perform Singular Value Decomposition (SVD) to find the best-fit line direction
    TDecompSVD svd(covariance);
    TMatrixD eigenVectors = svd.GetU();
//    std::cout << "U: "; eigenVectors.Print();
//    std::cout << "sigma: "; svd.GetSig().Print();

    // The direction of the best-fit line is given by the eigenvector with the largest eigenvalue
    double fac = 1.0;
    if(eigenVectors(2, 0)<0) fac = -1.0;
    TVector3 direction(fac*eigenVectors(0, 0), fac*eigenVectors(1, 0), fac*eigenVectors(2, 0));

    return direction;
}

// Calculate the perpendicular distance from a point to a line
static double pointLineDistance(const ROOT::Math::XYZVector& point, const TVector3& direction, const TVector3& centroid) {
    TVector3 pointVec(point.x(), point.y(), point.z());
    TVector3 pointToCentroid = pointVec - centroid;
    TVector3 crossProduct = pointToCentroid.Cross(direction);
    return crossProduct.Mag() / direction.Mag();
}

static double calculateSSR(const struct TPORec::TRACK &track, const TVector3& direction, const TVector3& centroid) {
    std::vector<TPORec::TRACKHIT> tkhit = track.tkhit;
    double ssr = 0.0;
    for (const auto& hit : tkhit) {
        double distance = pointLineDistance(hit.point, direction, centroid);
        ssr += distance * distance;
    }
    return ssr;
}

void TPORecoEvent::TrackReconstruct() {

    for (auto it : fPORecs)
    {
        it->fTracks.clear();
        // consider only primary track
        //        DigitizedTrack *dt = it->DTs[0];
        int idx = 0;
        for (auto &dt : it->DTs)
        {
            struct TPORec::TRACK track;
            size_t nhits = dt->fEnergyDeposits.size();
            for (size_t i = 0; i < nhits; i++)
            {
                long ID = dt->fhitIDs[i];
                long hittype = fTcalEvent->getChannelTypefromID(ID);
                if (hittype != 1)
                    continue; // only tracker hits
                ROOT::Math::XYZVector point = fTcalEvent->getChannelXYZfromID(ID);
                float ehit = dt->fEnergyDeposits[i];
                struct TPORec::TRACKHIT tkhit = {ID, point, ehit};
                track.tkhit.push_back(tkhit);
            }

            // now fit track of primary
            if (nhits > 1 && idx ==0)
            {
                track.direction = fitLineThroughPoints(track, track.centroid);
                track.SSR = calculateSSR(track, track.direction, track.centroid);
            }
            else
            {
                track.direction = {0, 0, 0};
                track.SSR = -1;
            }
            it->fTracks.push_back(track);
            idx++;
        }
    }
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
        std::cout << "TRACK>> ";
        std::cout << "ntrack = " << it->fTracks.size() << " ";
        if(it->fTracks.size() > 0){
        std::cout << " tracks: ";
        //for (auto itrk : it->fTracks)
        struct TPORec::TRACK itrk = it->fTracks[0]; 
        {
            std::cout << "nhit=" << itrk.tkhit.size() << " ";
            std::cout << "dir: " << itrk.direction.x() << " " << itrk.direction.y() << " " << itrk.direction.z() << " ";
            std::cout << "SSR: " << itrk.SSR << " ";
        }
        }
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
    std::cout << "FULL EVENT>> Ene: " << fPOFullEvent->TotalEvis() << " ";;
    std::cout << " ET: " << fPOFullEvent->TotalET();
    std::cout << std::endl;
    fTPOEvent->dump_header();
}

void TPORecoEvent::Reconstruct2DViewsPS() {

    PShitmapX.clear();
    PShitmapY.clear();
    for (auto &it : PShitmapsZ) it.second.clear();
    PShitmapsZ.clear();

    int nx = fTcalEvent->geom_detector.fScintillatorSizeX / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int ny = fTcalEvent->geom_detector.fScintillatorSizeY / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int nzlayer = fTcalEvent->geom_detector.fSandwichLength / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int nztot = fTcalEvent->geom_detector.NRep * nzlayer;

    for (auto it : fPORecs)
    {
        int ntracks = it->fGEANTTrackIDs.size();
        for (size_t i = 0; i < ntracks; i++)
        {
            // only primaries
          ////  if(i>0) continue;
            /////////////////
            DigitizedTrack *dt = it->DTs[i];
            size_t nhits = dt->fEnergyDeposits.size();
            for (size_t j = 0; j < nhits; j++)
            {
                long ID = dt->fhitIDs[j];
                long hittype = ID / 100000000000LL;
                if (hittype != 0)
                    continue;
                long ix = ID % 1000;
                long iy = (ID / 1000) % 1000;
                long iz = (ID / 1000000) % 1000;
                long ilayer = (ID / 1000000000);
                float ehit = dt->fEnergyDeposits[j];
                int fPDG = dt->fPDG;
                float electromagnetic = 0;
                if (fabs(dt->fPDG) == 11)
                    electromagnetic = 1.0;
                // XZ view
                long IDX = ix + iz * 1000000 + ilayer * 1000000000;
                auto hitX = PShitmapX.find(IDX);
                if (hitX != PShitmapX.end())
                {
                    hitX->second.Edeposited += ehit;
                    hitX->second.electromagneticity += electromagnetic;
                    hitX->second.ntracks++;
                }
                else
                {
                    PSHIT2D hitx = {electromagnetic, 1, ehit};
                    PShitmapX[IDX] = hitx;
                }
                // YZ view
                long IDY = iy * 1000 + iz * 1000000 + ilayer * 1000000000;
                auto hitY = PShitmapY.find(IDY);
                if (hitY != PShitmapY.end())
                {
                    hitY->second.Edeposited += ehit;
                    hitY->second.electromagneticity += electromagnetic;
                    hitY->second.ntracks++;
                }
                else
                {
                    PSHIT2D hity = {electromagnetic, 1, ehit};
                    PShitmapY[IDY] = hity;
                }
                // XY view per layer
                long IDZ = ix + iy * 1000 + ilayer * 1000000000;
                auto layerhitmap = PShitmapsZ.find(ilayer);
                if (layerhitmap != PShitmapsZ.end()) {
                    auto hitZ = layerhitmap->second.find(IDZ);
                    if (hitZ != layerhitmap->second.end()) {
                        hitZ->second.Edeposited += ehit;
                        hitZ->second.electromagneticity += electromagnetic;
                        hitZ->second.ntracks++;
                    } else {
                        PSHIT2D hitz = {electromagnetic, 1, ehit};
                        layerhitmap->second[IDZ] = hitz;
                    }
                } else {
                    std::map<long, PSHIT2D> xymap;
                    PSHIT2D hitz = {electromagnetic, 1, ehit};
                    xymap[IDZ] = hitz;
                    PShitmapsZ[ilayer] = xymap;
                }
            }
        }
    }

    // renormalize the electromagneticity
    for (auto& it : PShitmapX)
    {
        it.second.electromagneticity /= float(it.second.ntracks);
    }
    for (auto& it : PShitmapY)
    {
        it.second.electromagneticity /= float(it.second.ntracks);
    }
    for (auto& itm : PShitmapsZ)
    {
        for (auto& it : itm.second) {
            it.second.electromagneticity /= float(it.second.ntracks);
        }
    }
}

void TPORecoEvent::Fill2DViewsPS() {
    int nx = fTcalEvent->geom_detector.fScintillatorSizeX / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int ny = fTcalEvent->geom_detector.fScintillatorSizeY / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int nzlayer = fTcalEvent->geom_detector.fSandwichLength / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int nztot = fTcalEvent->geom_detector.NRep * nzlayer;
    // fill histrograms
    xviewPS = (TH2D*)gDirectory->Get("xviewPS");
    if(xviewPS != nullptr) {
        xviewPS->Reset(); 
        yviewPS = (TH2D*)gDirectory->Get("yviewPS");
        yviewPS->Reset(); 
        for (int i=0; i<50; i++) {
            std::string hname = "zviewPS_" + std::to_string(i);
            zviewPS[i] = (TH2D*)gDirectory->Get(hname.c_str());
            zviewPS[i]->Reset(); 
        }
        xviewPS_em = (TH2D*)gDirectory->Get("xviewPS_em");
        xviewPS_em->Reset(); 
        yviewPS_em = (TH2D*)gDirectory->Get("yviewPS_em");
        yviewPS_em->Reset(); 
        xviewPS_had = (TH2D*)gDirectory->Get("xviewPS_had");
        xviewPS_had->Reset(); 
        yviewPS_had = (TH2D*)gDirectory->Get("yviewPS_had");
        yviewPS_had->Reset(); 
        xviewPS_eldepo = (TH2D*)gDirectory->Get("xviewPS_eldepo");
        xviewPS_eldepo->Reset(); 
        yviewPS_eldepo = (TH2D*)gDirectory->Get("yviewPS_eldepo");
        yviewPS_eldepo->Reset(); 
    } else {
    xviewPS = new TH2D("xviewPS", "Scintillator x-view", nztot, 0, nztot, nx, 0, nx);
    yviewPS = new TH2D("yviewPS", "Scintillator y-view", nztot, 0, nztot, ny, 0, ny);
    for (int i=0; i<50; i++) {
        std::string hname = "zviewPS_" + std::to_string(i);
        zviewPS[i] = new TH2D(hname.c_str(), "Scintillator y-view", nx, 0, nx, ny, 0, ny);
    }
    xviewPS_em = new TH2D("xviewPS_em", "Scintillator x-view - EM", nztot, 0, nztot, nx, 0, nx);
    yviewPS_em = new TH2D("yviewPS_em", "Scintillator y-view - EM", nztot, 0, nztot, ny, 0, ny);
    xviewPS_had = new TH2D("xviewPS_had", "Scintillator x-view - HAD", nztot, 0, nztot, nx, 0, nx);
    yviewPS_had = new TH2D("yviewPS_had", "Scintillator y-view - HAD", nztot, 0, nztot, ny, 0, ny);
    xviewPS_eldepo = new TH2D("xviewPS_eldepo", "Scintillator x-view", 11, 0.,1.1,100.,0.,25.);
    yviewPS_eldepo = new TH2D("yviewPS_eldepo", "Scintillator y-view", 11, 0.,1.1,100.,0.,25.);
    }

    double electromagneticity_threshold = 0.8;
    for (auto it : PShitmapX)
    {
        long ID = it.first;
        double ehit = it.second.Edeposited;
        double elec = it.second.electromagneticity;
        long ix = ID % 1000;
        long iz = (ID / 1000000) % 1000;
        long ilayer = (ID / 1000000000);
        double fix = ix + 0.5;
        double fiz = ilayer * nzlayer + iz + 0.5;
        xviewPS->Fill(fiz, fix, ehit);
        if(elec>electromagneticity_threshold) {
            xviewPS_em->Fill(fiz, fix, ehit);
        } else {
            xviewPS_had->Fill(fiz, fix, ehit);
        }
        xviewPS_eldepo -> Fill(elec, ehit);
    }

    for (auto it : PShitmapY)
    {
        long ID = it.first;
        double ehit = it.second.Edeposited;
        double elec = it.second.electromagneticity;
        long iy = (ID / 1000) % 1000;
        long iz = (ID / 1000000) % 1000;
        long ilayer = (ID / 1000000000);
        double fiy = iy + 0.5;
        double fiz = ilayer * nzlayer + iz + 0.5;
        yviewPS->Fill(fiz, fiy, ehit);
        if(elec>electromagneticity_threshold) {
            yviewPS_em->Fill(fiz, fiy, ehit);
        } else {
            yviewPS_had->Fill(fiz, fiy, ehit);
        }
        yviewPS_eldepo -> Fill(elec, ehit);
    }

    for (auto itm : PShitmapsZ)
    {
        for (auto it : itm.second) {
            long ID = it.first;
            double ehit = it.second.Edeposited;
            long ix = ID % 1000;
            long iy = (ID / 1000) % 1000;
            long ilayer = (ID / 1000000000);
            double fix = ix + 0.5;
            double fiy = iy + 0.5;
            zviewPS[ilayer]->Fill(fix, fiy, ehit);
        }
    }

#if 0
    // Visualize the results
    TCanvas *c1 = new TCanvas("xyviews", "xyviews", 800, 600);
    c1->Divide(2,2);
    for (int i=0; i<4; i++) {
        c1->cd(i+1);
        zviewPS[i+3]->Draw();
    }
    c1->Modified();
    c1->Update();
#endif

}

void TPORecoEvent::pshit2d_position(long ID, double &fix, double &fiy, double &fiz) {
    int nzlayer = fTcalEvent->geom_detector.fSandwichLength / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    long ix = ID % 1000;
    long iy = (ID / 1000) % 1000;
    long iz = (ID / 1000000) % 1000;
    long ilayer = (ID / 1000000000);
    fix = ix + 0.5;
    fiy = iy + 0.5;
    fiz = ilayer * nzlayer + iz + 0.5;
}

/// @brief Reconstruct all 2D clusters for the xz and the yz views
/// @param view = 0 for XZ, and .ne.0 for YZ
void TPORecoEvent::ReconstructClusters(int view, bool verbose) {

    double threshold_2dhit = 2.0; // MeV
    double eps = 5; // mm
    int minPts = 10; // minimum number of points

    DBScan dbscan;

    std::map<int, class TPSCluster> *PSClusters = (view==0) ? &PSClustersX : &PSClustersY;
    PSClusters->clear();

    std::vector<DBScan::Point> points;

    // XZ or YZ view view
    for (auto it : (view==0) ? PShitmapX : PShitmapY)
    {
        long ID = it.first;
        double ehit = it.second.Edeposited;
        if(ehit < threshold_2dhit) continue;
        double fix, fiy, fiz;
        pshit2d_position(ID, fix, fiy, fiz);
        DBScan::Point p = {ID, ehit, (view==0) ? fix : fiy, fiz};
        points.push_back(p);
    }
    dbscan.scan(points, eps, minPts);
    for (const auto& point : points) {
        if(point.clusterID == 0)continue;

        TPSCluster::PSCLUSTERHIT hit = {point.ID, (float)point.ehit}; 
        auto c = PSClusters->find(point.clusterID);
        if (c != PSClusters->end())
        {
            c->second.hits.push_back(hit);
            c->second.rawenergy += point.ehit;
        }
        else
        {
            TPSCluster newc(view, fTcalEvent);
            newc.clusterID = point.clusterID;
            newc.rawenergy = point.ehit;
            newc.hits.push_back(hit);
            (*PSClusters)[point.clusterID] = newc;
        }
    }

    // now compute the energy compensated eflow relative to primary vertex
    ROOT::Math::XYZVector primary(fTPOEvent->prim_vx.X(), fTPOEvent->prim_vx.Y(), fTPOEvent->prim_vx.Z());

    for (auto& c : *PSClusters) {
        std::cout << "Cluster ID:" << c.second.clusterID << " nhits=" << c.second.hits.size();
        std::cout << " rawEnergy(MeV): " << c.second.rawenergy; 
        std::cout << std::endl;        
        if(c.second.rawenergy < 10*1e3) continue;
        c.second.ComputeCOG();
        c.second.setVtx(primary.x(), primary.y(), primary.z());
        c.second.ComputeLongProfile(verbose);
    }
}

void TPORecoEvent::Reconstruct3DPS(int maxIter) {

    std::random_device rd;  // Seed for the random number generator
    std::mt19937 gen(rd());  // Mersenne Twister random number generator
 
    double ehit_threshold = 0.5; // MeV
    int nvox_per_layer_max = 150;

    int nx = fTcalEvent->geom_detector.fScintillatorSizeX / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int ny = fTcalEvent->geom_detector.fScintillatorSizeY / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int nzlayer = fTcalEvent->geom_detector.fSandwichLength / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int nrep = fTcalEvent->geom_detector.NRep;
    int nztot =  nrep * nzlayer;

    std::uniform_int_distribution<> rnd_layer(0, nzlayer-1);
    std::uniform_int_distribution<> rnd_nx(0, nx-1);
    std::uniform_int_distribution<> rnd_ny(0, ny-1);

    // the maximum number layer (from 0 to nRep) that is reconstructed
    int maxLayer = 25;

    // Step 0: organize hits for easy access
    std::vector<std::vector<float>> XZ(nx, std::vector<float>(nztot, 0.0));
    std::vector<std::vector<float>> YZ(ny, std::vector<float>(nztot, 0.0));
    std::vector<std::vector<std::vector<float>>> XY(nrep, 
            std::vector<std::vector<float>>(nx, std::vector<float>(ny, 0.0)));

#if 0
// debug
    nztot = 1;
    nzlayer = 1;
    nx = 5;
    ny = 5;

    XZ[2][0] = 1.0;
    XZ[1][0] = 1.0;
//    XZ[4][0] = 1.0;

    YZ[0][0] = 1.0;
//    YZ[1][0] = 1.0;
    YZ[3][0] = 1.0;
//    YZ[4][0] = 1.0;
#endif

    #if 1
    for (auto it : PShitmapX)
    {
        long ID = it.first;
        double ehit = it.second.Edeposited;
        if(ehit < ehit_threshold) continue;
        long ix = ID % 1000;
        long iz = (ID / 1000000) % 1000;
        long ilayer = (ID / 1000000000);
        if(ilayer > maxLayer) continue;
        long izz = ilayer * nzlayer + iz;
        if(ix < nx && izz < nztot) {
            XZ[ix][izz] = ehit;
        }
    }
    for (auto it : PShitmapY)
    {
        long ID = it.first;
        double ehit = it.second.Edeposited;
        if(ehit < ehit_threshold) continue;
        long iy = (ID / 1000) % 1000;
        long iz = (ID / 1000000) % 1000;
        long ilayer = (ID / 1000000000);
        if(ilayer > maxLayer) continue;
        long izz = ilayer * nzlayer + iz;
        if(iy < ny && izz < nztot) {
            YZ[iy][izz] = ehit;
        }
    }
    for (auto itm : PShitmapsZ) {
        for (auto it: itm.second) {
            long ID = it.first;
            double ehit = it.second.Edeposited;
            if (ehit < ehit_threshold)
                continue;
            long ix = ID % 1000;
            long iy = (ID / 1000) % 1000;
            long ilayer = (ID / 1000000000);
            if (ilayer > maxLayer)
                continue;
            if(iy < ny) {
                XY[ilayer][ix][iy] = ehit;
            }
        }
    }
#endif

    struct Voxel {
        float value;
        Voxel() : value(0) {}
    };

    std::vector<std::vector<std::vector<Voxel>>> V(
        nx, std::vector<std::vector<Voxel>>(
                   ny, std::vector<Voxel>(
                               nztot, Voxel())));

    std::vector<std::vector<std::vector<Voxel>>> V_min(
        nx, std::vector<std::vector<Voxel>>(
                   ny, std::vector<Voxel>(
                               nztot, Voxel())));

    double min_score = 1e9;
    double prev_score = 1e9;

    // Step 1: Initial assignment based on projections
    for (int z = 0; z < nztot; ++z)
    {
        long ilayer = z / nrep;
        for (int x = 0; x < nx; ++x)
        {
            if (XZ[x][z] > 0)
            {
                for (int y = 0; y < ny; ++y)
                {
                    if (YZ[y][z] > 0)
                    {
                        if(XY[ilayer][x][y] > 0) {
                            V[x][y][z].value = std::min(XZ[x][z], std::min(YZ[y][z], XY[ilayer][x][y]));
                        }
                    }
                }
            }
        }
    }

     // decide which layers should be used for reconstructing 3D voxels
    int nvox_per_layer[nztot];
    for (int z = 0; z < nztot; ++z)
    {
        int nvox_layer = 0;
        for (int x = 0; x < nx; ++x)
            for (int y = 0; y < ny; ++y) {
                if(V[x][y][z].value>0) {
                    nvox_layer++;
                }
            }
        nvox_per_layer[z] = nvox_layer;
    }

    // Step 2: Iterative Refinement to match projections
    for (int iter = 0; iter < maxIter; ++iter) {

    for (int z = 0; z < nztot; ++z)
    {
        int nvox_layer = 0;
        for (int x = 0; x < nx; ++x)
            for (int y = 0; y < ny; ++y) {
                if(V[x][y][z].value>0) {
                    nvox_layer++;
                }
            }
        nvox_per_layer[z] = nvox_layer;
    }

        if (verbose > 4)
        {
            // print plane
            std::cout << "ITERATION " << iter << " ---------------------------------" << std::endl;
            for (int x = 0; x < nx; ++x)
            {
                double sumXZ = 0;
                for (int y = 0; y < ny; ++y)
                {
                    std::cout << V[x][y][0].value << " ";
                    sumXZ += V[x][y][0].value;
                }
                std::cout << "   s=" << sumXZ;
                std::cout << "   a=" << XZ[x][0];
                std::cout << std::endl;
            }
            for (int y = 0; y < nx; ++y)
            {
                double sumYZ = 0;
                for (int x = 0; x < nx; ++x)
                {
                    sumYZ += V[x][y][0].value;
                }
                std::cout << sumYZ << " ";
            }
            std::cout << ": s" << std::endl;
            for (int y = 0; y < nx; ++y)
            {
                std::cout << YZ[y][0] << " ";
            }
            std::cout << ": a" << std::endl;
        }

        int adjusted = 0;
        double score = 0;

//        verbose = 4; // FIXME: DEBUG

        // adjust XY
        for (int ilayer = 0; ilayer < nrep; ilayer++)
        {
            // count number of voxel in this layer
            int sum_nvox = 0;
            for (int z = 0; z < nzlayer; z++){
                int izz = ilayer * nzlayer + z;
                sum_nvox += nvox_per_layer[izz];
            }
            if(sum_nvox > nvox_per_layer_max) continue;

            for (int x = 0; x < nx; ++x)
            {
                for (int y = 0; y < ny; ++y)
                {
                    double sumXY = 0;
                    for (int z = 0; z < nzlayer; z++)
                    {
                        int izz = ilayer * nzlayer + z;
                        sumXY += V[x][y][izz].value;
                    }
                    float difference = XY[ilayer][x][y] - sumXY;
                    score += abs(difference);
                    if(verbose>3 && abs(difference)>1e-3) 
                        std::cout << " ilayer: " << ilayer << " difference=" << difference << std::endl;
                    int niter = 0;
                    while(abs(difference)>1e-3 && niter++ < 2*nzlayer) {
                        int z = rnd_layer(gen);
                        int izz = ilayer * nzlayer + z;

                        // check difference of two 2D projections
 //                       double diff = fabs(XZ[x][izz] - YZ[y][izz]);
 //                       if(diff < 1e-3) continue;  // this is likely not a fake

                        float adjust = std::min(difference, XY[ilayer][x][y] - V[x][y][izz].value);
                        if(V[x][y][izz].value + adjust < 0) adjust = -V[x][y][izz].value;
                        if(verbose>3 && abs(adjust)>1e-3) 
                            std::cout << "    z: " << z << "  adjust " << adjust << std::endl;
                        V[x][y][izz].value += adjust;
                        difference -= adjust;
                        if(adjust!=0) adjusted++;
                    }
                }
            }
        }

        // adjust XZ and YZ
        for (int z = 0; z < nztot; ++z) {
            if(nvox_per_layer[z] > nvox_per_layer_max) continue;
            
            // Adjust XZ projection
            for (int x = 0; x < nx; ++x) {
                double sumXZ = 0;
                for (int y = 0; y < ny; ++y) {
                    sumXZ += V[x][y][z].value;
                }
                float difference = XZ[x][z] - sumXZ;
                score += abs(difference);
                if(verbose>3 && abs(difference)>1e-3) std::cout << " x: " << x << " difference=" << difference << std::endl;
                int niter = 0;
                while(abs(difference)>1e-3 && niter++ < ny) {
                    int y = rnd_ny(gen);
                    float adjust = std::min(difference, YZ[y][z] - V[x][y][z].value);
                    if(V[x][y][z].value + adjust < 0) adjust = -V[x][y][z].value;
                    if(verbose>3 && abs(adjust)>1e-3) std::cout << "    y: " << y << "  adjust " << adjust << std::endl;
                    V[x][y][z].value += adjust;
                    difference -= adjust;
                    if(adjust!=0) adjusted++;
                }
            }

            if (verbose > 4)
            {

                // print plane
                std::cout << "AFTER XZ " << iter << " ---------------------------------" << std::endl;
                for (int x = 0; x < nx; ++x)
                {
                    double sumXZ = 0;
                    for (int y = 0; y < ny; ++y)
                    {
                        std::cout << V[x][y][0].value << " ";
                        sumXZ += V[x][y][0].value;
                    }
                    std::cout << "   s=" << sumXZ;
                    std::cout << "   a=" << XZ[x][0];
                    std::cout << std::endl;
                }
                for (int y = 0; y < nx; ++y)
                {
                    double sumYZ = 0;
                    for (int x = 0; x < nx; ++x)
                    {
                        sumYZ += V[x][y][0].value;
                    }
                    std::cout << sumYZ << " ";
                }
                std::cout << ": s" << std::endl;
                for (int y = 0; y < nx; ++y)
                {
                    std::cout << YZ[y][0] << " ";
                }
                std::cout << ": a" << std::endl;
            }

           // Adjust YZ projection
            for (int y = 0; y < ny; ++y) {
                double sumYZ = 0;
                for (int x = 0; x < nx; ++x) {
                    sumYZ += V[x][y][z].value;
                }
                float difference = YZ[y][z] - sumYZ;
                score += abs(difference);
                if(verbose>3 && abs(difference)>1e-3)std::cout << " y: " << y << " difference=" << difference << std::endl;
                int niter = 0;
                while(abs(difference)>1e-3 && niter++ < nx) {
                    int x = rnd_nx(gen);
                    float adjust = std::min(difference, XZ[x][z] - V[x][y][z].value);
                    if(V[x][y][z].value + adjust < 0) adjust = -V[x][y][z].value;
                    if(verbose>3 && abs(adjust)>1e-3) std::cout << "    x: " << x << "  adjust " << adjust << std::endl;
                    V[x][y][z].value += adjust;
                    difference -= adjust;
                    if(adjust!=0) adjusted++;
                }
            }
        }

        if(verbose>0 && adjusted>0) std::cout << "Iteration " << iter << ": " << adjusted << 
        " voxels adjusted. Score = " << score << std::endl;

        // save solution of with minimum score
        if(score < min_score){
        for (int z = 0; z < nztot; ++z)
            for (int x = 0; x < nx; ++x)
                for (int y = 0; y < ny; ++y)
                {
                    V_min[x][y][z].value = V[x][y][z].value;
                }
            min_score = score;
        }

        // if score degraded, switch back to best answer
        if(score > prev_score) {
        for (int z = 0; z < nztot; ++z)
            for (int x = 0; x < nx; ++x)
                for (int y = 0; y < ny; ++y)
                {
                    V[x][y][z].value = V_min[x][y][z].value;
                }
        }
        prev_score = score;

        if(score == 0) break;
    }

    std::cout << " We keep best score solution:" << min_score << std::endl;

    PSvoxelmap.clear();
    for (int z = 0; z < nztot; ++z) {
        // 
            if(nvox_per_layer[z] > nvox_per_layer_max) continue;
        //
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                float ehit = V_min[x][y][z].value;
                if(ehit>0) {
                    long iz = z % nzlayer;
                    long ilayer = z / nzlayer;
    	        	long ID = x + y*1000 + iz*1000000 + ilayer*1000000000;
                    struct PSVOXEL3D v = {ehit, true};
                    PSvoxelmap[ID] = v;
                }
            }
        }
    }

    // now flag real hits
    for (const auto& track : fTcalEvent->getfTracks()) {
        size_t nhits = track->fhitIDs.size();
        for ( size_t i = 0; i < nhits; i++) {
            auto v = PSvoxelmap.find(track->fhitIDs[i]);
            if (v != PSvoxelmap.end()) {
                v->second.ghost = false;
            }
        }
    }

    // some stats on hits
    size_t ntot = 0;
    size_t fakes = 0;
    size_t nreal = 0;
    double ave_e_ghost = 0;
    double ave_e_real = 0;

    TH1D h_diff_fake = TH1D("h_diff_fake", "difference XZ vs YZ fakes", 100, 0., 30.);
    TH1D h_diff_real = TH1D("h_diff_real", "difference XZ vs YZ real", 100, 0., 30.);

    for (const auto& v : PSvoxelmap) {
        ntot++;
        long ID = v.first;
        long ix = ID % 1000;
        long iy = (ID / 1000) % 1000;
        long iz = (ID / 1000000) % 1000;
        long ilayer = (ID / 1000000000);
        if(ilayer > maxLayer) continue;
        long izz = ilayer * nzlayer + iz;
        double eXZ = XZ[ix][izz];
        double eYZ = YZ[iy][izz];
        double diff = fabs(eXZ - eYZ);

        if(v.second.ghost) {
            fakes++;
            ave_e_ghost += v.second.RawEnergy;
            h_diff_fake.Fill(diff);
        } else {
            ave_e_real += v.second.RawEnergy;
            nreal++;
            h_diff_real.Fill(diff);
        }
    }
    std::cout << " STATS: " << ntot << " hits " << fakes << " ghosts." << std::endl;
    std::cout << " Avg energy: " << ave_e_real/float(nreal) << " ghosts: " << ave_e_ghost/float(fakes);
    std::cout << std::endl;

#if 0
    // Visualize the results
    TCanvas *c1 = new TCanvas("3dreco", "3dreco", 800, 600);
    c1->Divide(1,2);
    c1->cd(1);
    h_diff_real.Draw();
    c1->cd(2);
    h_diff_fake.Draw();
    c1->Modified();
    c1->Update();
    c1->SaveAs("3dreco.png");
#endif
}
