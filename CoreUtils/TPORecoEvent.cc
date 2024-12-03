#include "TPORecoEvent.hh"
#include "DBScan.hh"
#include "TTKTrack.hh"
#include "TPSTrack.hh"

#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <TDecompSVD.h>
#include <TVector3.h>
#include <TCanvas.h>

// Genfit
#include <GFRaveVertexFactory.h>
#include <GFRaveVertex.h>

#include <random>
#include <thread>

ClassImp(TPORec);
ClassImp(TPORecoEvent);

TPORecoEvent::TPORecoEvent(TcalEvent* c, TPOEvent* p) : fTcalEvent(c), fTPOEvent(p), fPOFullEvent(0), fPOFullRecoEvent(0) {
    // copy geometry
    geom_detector = fTcalEvent->geom_detector;

    // empty histogram pointers
	for(int i =0; i < 50; i++){
		zviewPS.push_back(nullptr);
	}
};

TPORecoEvent::~TPORecoEvent() {
    for(auto it : fPORecs) {
        delete it;
    }
    delete fPOFullEvent;
    delete fPOFullRecoEvent;
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
            long layer = fTcalEvent->getChannelModulefromID(dt->fhitIDs[i]);
            if(layer < minlayer) minlayer = layer;
        }
        for (size_t i = 0; i < nhits; i++)
        {
            // loop over hits in the tracker
            hittype = fTcalEvent->getChannelTypefromID(dt->fhitIDs[i]);
            if (hittype != 1) continue;
            long layer = fTcalEvent->getChannelModulefromID(dt->fhitIDs[i]);
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
    int nrep = fTcalEvent->geom_detector.NRep;

    // fill map of hits
    std::map <int, std::vector<struct TTKTrack::TRACKHIT>> hitMap;

    for (const auto &track : fTcalEvent->getfTracks())
    {
        size_t nhits = track->fhitIDs.size();
        for (size_t i = 0; i < nhits; i++)
        {
            long ID = track->fhitIDs[i];
            float ehit = track->fEnergyDeposits[i];
            long hittype = fTcalEvent->getChannelTypefromID(ID);
            if (hittype != 1)
                continue;
            int ilayer = fTcalEvent->getChannelModulefromID(ID);
            ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(track->fhitIDs[i]);
            struct TTKTrack::TRACKHIT hit = {ID, position, ehit};
            auto hitlayer = hitMap.find(ilayer);
            if(hitlayer != hitMap.end()) {
                hitlayer->second.push_back(hit);
            } else {
                std::vector<struct TTKTrack::TRACKHIT> hits;
                hits.push_back(hit);
                hitMap[ilayer] = hits;
            }
        }
    }

    fTKTracks.clear();
    // now match doublets in each layer
    for(auto &it : hitMap) {
        int ilayer = it.first;
        for (const auto &hit1 : it.second) {
            int icopy1 = fTcalEvent->getChannelCopyfromID(hit1.ID);
            if(icopy1 != 1) continue;
            double distmin = 1e9;
            struct TTKTrack::TRACKHIT hitmin;
            for (const auto &hit2 : it.second) {
                int icopy2 = fTcalEvent->getChannelCopyfromID(hit2.ID);
                if(icopy2 != 2) continue;
                ROOT::Math::XYZVector distV = (hit2.point - hit1.point);
                // check if hits are close in XY plane
                double dist = distV.X()*distV.X() + distV.Y()*distV.Y();
//                double dist = (hit2.point - hit1.point).R();
                if(dist < distmin) {
                    distmin = dist;
                    hitmin = hit2;
                }
            }

            // FIXME: define cut
            if(distmin>5.0) continue;

            // now create a TKTrack with the doublet if hits are close
            TTKTrack trk;
            trk.tkhit.push_back(hit1);
            trk.tkhit.push_back(hitmin);
            trk.direction = trk.direction2Hits();
            ROOT::Math::XYZVector centroid = (hit1.point + hitmin.point) * 0.5;
            trk.centroid.SetXYZ(centroid.X(), centroid.Y(), centroid.Z());
            fTKTracks.push_back(trk);
            // remove hit min from the list of hits
            it.second.erase(std::remove_if(it.second.begin(), it.second.end(),
                                           [&hitmin](const TTKTrack::TRACKHIT &hit)
                                           {
                                               // Use the ID to match the hit to be removed
                                               return hit.ID == hitmin.ID;
                                           }),
                            it.second.end());
        }
    }

    // now sort doublets by z coordinate
    std::sort(fTKTracks.begin(), fTKTracks.end(), [](const TTKTrack &a, const TTKTrack &b)
              { return a.centroid.Z() < b.centroid.Z(); });

    // now loop over tracks and merge segments
    for (auto &trk : fTKTracks) {
        trk.SortHitsByZ();
    }

    double parallel_cut = 0.01; // FIXME: adjust value
    double mindZcut = fTcalEvent->geom_detector.fTargetSizeZ*1.1;

    for (size_t i = 0; i < fTKTracks.size(); i++) {
        TTKTrack &track1 = fTKTracks[i];
        // for subsequent iterations, consider only tracks with more than 2 hits
        size_t besttrk = -1;
        double mindZ = 1e9;
        for( size_t j = 0; j < fTKTracks.size(); j++) {
            if (i==j) continue;
            TTKTrack &track2 = fTKTracks[j];
            // check if segments are parallel
            TVector3 normDir1 = track1.direction.Unit();
            TVector3 normDir2 = track2.direction.Unit();
            double dotProduct = normDir1.Dot(normDir2);
            if (std::abs(std::abs(dotProduct) - 1.0) > parallel_cut) continue;
            // ensure that segments belong to different planes
            struct TTKTrack::TRACKHIT hit1 = track1.tkhit.back();
            struct TTKTrack::TRACKHIT hit2 = track2.tkhit.front();
            double dz = std::abs(hit1.point.z()-hit2.point.z());
            if(dz < mindZcut) continue;
/*            // ensure that segments don't start in the same plane
            struct TTKTrack::TRACKHIT hit1b = track1.tkhit.back();
            struct TTKTrack::TRACKHIT hit2b = track2.tkhit.back();
            dz = std::abs(hit1b.point.z()-hit2b.point.z());
            if(dz < mindZcut) continue;
            */
            // ensure that segments are parallel to main line joining hits
            ROOT::Math::XYZVector normDir = (hit2.point - hit1.point).Unit();
            dotProduct = std::min(std::abs(normDir.Dot(normDir1)), std::abs(normDir.Dot(normDir2)));
            if (std::abs(dotProduct - 1.0) > parallel_cut) continue;
            if(dz < mindZ) {
                mindZ = dz;
                besttrk = j;
            }
        }
        // if match found, then merge the tracks
        if(besttrk != -1) {
            // add hits from track2 to track1
            for (const auto& h : fTKTracks[besttrk].tkhit) {
                track1.tkhit.push_back(h);
            }
            track1.direction = track1.fitLineThroughHits(track1.centroid);
            fTKTracks.erase(fTKTracks.begin() + besttrk);
            if(besttrk < i) i--;
        }
    }

    // get rid of tracks with less than 3 hits
    fTKTracks.erase(std::remove_if(fTKTracks.begin(), fTKTracks.end(),
                                   [](const TTKTrack &track)
                                   {
                                       return track.tkhit.size() < 3;
                                   }),
                    fTKTracks.end());

    // always order hits in the track by increasing "z"
    for (auto &trk : fTKTracks) {
        trk.SortHitsByZ();
    }

    // merge broken tracks
    double cut_SSR_merge = 10; // FIXME: adjust value
    for (size_t i = 0; i < fTKTracks.size(); i++) {
        TTKTrack &track1 = fTKTracks[i];
        for( size_t j = 0; j < fTKTracks.size(); j++) {
            if (i==j) continue;
            TTKTrack &track2 = fTKTracks[j];
            TTKTrack *tempTrack = new TTKTrack(track1);
            tempTrack->MergeTracks(track2);
            // fit line through hits
            tempTrack->direction = tempTrack->fitLineThroughHits(tempTrack->centroid);
            if(verbose > 5) {
                std::cout << "Merging tracks: " << i << " " << j << std::endl;
                track1.Dump();
                track2.Dump();
                tempTrack->Dump();
            }
            if(tempTrack->SSR < cut_SSR_merge) {
                fTKTracks[i] = *tempTrack;
                fTKTracks.erase(fTKTracks.begin() + j);
                if(j < i) i--;
                break;
            }
            delete tempTrack;
        }
    }

    // now fit tracks with GENFIT2
    for (auto &trk : fTKTracks) {
        trk.GenFitTrackFit();
        if(verbose > 5) {
            trk.fitTrack->Print();
        }
    }

    if(verbose > 2) {
        DumpReconstructedTracks();
    }
}

void TPORecoEvent::DumpReconstructedTracks() {
    for (auto &trk : fTKTracks) {
        trk.Dump();
    }
}

void TPORecoEvent::FitTrackVertices() {
    genfit::GFRaveVertexFactory vertexFactory(2);
    vertexFactory.setMethod("kalman-smoothing:1");

    std::vector<genfit::Track*> tracks;
    std::vector<genfit::GFRaveVertex*> vertices;

    for (auto &trk : fTKTracks) {
        tracks.push_back(trk.fitTrack);
    }

    // vertexing
    vertexFactory.findVertices(&vertices, tracks);

    // print reconstructed vertices
    for (auto &vertex : vertices) {
        vertex->Print();
    }
}

void TPORecoEvent::TrackReconstructTruth() {
// truth tracks
    for (auto it : fPORecs)
    {
        it->fTracks.clear();
        struct PO aPO = fTPOEvent->POs[it->POID];
        if(aPO.m_charge() == 0) continue;
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
    std::cout << std::endl;
    std::cout << "FULL EVENT>> Ene: " << fPOFullEvent->TotalEvis() << " ";
    std::cout << " ET: " << fPOFullEvent->TotalET();
    std::cout << std::endl;
    if(fPOFullRecoEvent!=nullptr) {
        std::cout << "FULL RECO EFLOW>> ";
        std::cout << " Ex: " << fPOFullRecoEvent->fTotal.Eflow.X();
        std::cout << " Ey: " << fPOFullRecoEvent->fTotal.Eflow.Y();
        std::cout << " Ez: " << fPOFullRecoEvent->fTotal.Eflow.Z();
        std::cout << std::endl;
        std::cout << "FULL RECO EVENT>> Ene: " << fPOFullRecoEvent->TotalEvis() << " ";
        std::cout << " ET: " << fPOFullRecoEvent->TotalET();
        std::cout << std::endl;
    }
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
            zviewPS.at(i) = (TH2D*)gDirectory->Get(hname.c_str());
            zviewPS.at(i)->Reset();
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
    xviewPS = new TH2D("xviewPS", "Scintillator xz-view", nztot, 0, nztot, nx, 0, nx);
    yviewPS = new TH2D("yviewPS", "Scintillator yz-view", nztot, 0, nztot, ny, 0, ny);
    for (int i=0; i<50; i++) {
        std::string hname = "zviewPS_" + std::to_string(i);
        zviewPS.at(i) = new TH2D(hname.c_str(), hname.c_str(), nx, 0, nx, ny, 0, ny);
    }
    xviewPS_em = new TH2D("xviewPS_em", "Scintillator xz-view - EM", nztot, 0, nztot, nx, 0, nx);
    yviewPS_em = new TH2D("yviewPS_em", "Scintillator yz-view - EM", nztot, 0, nztot, ny, 0, ny);
    xviewPS_had = new TH2D("xviewPS_had", "Scintillator xz-view - HAD", nztot, 0, nztot, nx, 0, nx);
    yviewPS_had = new TH2D("yviewPS_had", "Scintillator yz-view - HAD", nztot, 0, nztot, ny, 0, ny);
    xviewPS_eldepo = new TH2D("xviewPS_eldepo", "Scintillator xz-view", 11, 0.,1.1,100.,0.,25.);
    yviewPS_eldepo = new TH2D("yviewPS_eldepo", "Scintillator yz-view", 11, 0.,1.1,100.,0.,25.);
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
            zviewPS.at(ilayer)->Fill(fix, fiy, ehit);
        }
    }

#if 0
    std::cout << zviewPS[0] << std::endl;
    TCanvas *c1 = new TCanvas("","", 1000, 1000);
    zviewPS[2]->Draw("COLZ");
    c1->SaveAs("tempfile.png");
#endif

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
void TPORecoEvent::ReconstructClusters(int view) {

    double threshold_2dhit = 2.0; // MeV
    double eps = 5; // mm
    int minPts = 10; // minimum number of points

    DBScan dbscan;

    std::map<int, class TPSCluster> PSClustersMap;

    std::vector<DBScan::Point> points;

    // XZ or YZ view
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
        auto c = PSClustersMap.find(point.clusterID);
        if (c != PSClustersMap.end())
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
            PSClustersMap[point.clusterID] = newc;
        }
    }

    // now compute the energy compensated eflow relative to primary vertex
    ROOT::Math::XYZVector primary(fTPOEvent->prim_vx.X(), fTPOEvent->prim_vx.Y(), fTPOEvent->prim_vx.Z());

    // store in vector
    std::vector<TPSCluster> *PSClusters = (view==0) ? &PSClustersX : &PSClustersY;
    PSClusters->clear();
    for (auto& c : PSClustersMap) {
        if (c.second.rawenergy < 10 * 1e3)
            continue;
        PSClusters->push_back(c.second);
    }

    for (auto &c : *PSClusters)
    {
        c.ComputeCOG();
        c.setVtx(primary.x(), primary.y(), primary.z());
        c.ComputeLongProfile(verbose);
    }

    // Sort PSClusters by rawEnergy in descending order (highest energy first)
    std::sort(PSClusters->begin(), PSClusters->end(), [](const TPSCluster& a, const TPSCluster& b) {
    return a.rawenergy > b.rawenergy;
    });

    if(verbose > 0) {
        for (const auto &c : *PSClusters)
        {
            std::cout << "Cluster ID:" << c.clusterID << " nhits=" << c.hits.size();
            std::cout << " rawEnergy(MeV): " << c.rawenergy;
            std::cout << std::endl;
        }
    }
}

void TPORecoEvent::Reconstruct3DPS(int maxIter) {

    std::random_device rd;  // Seed for the random number generator
    std::mt19937 gen(rd());  // Mersenne Twister random number generator

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine rng(seed);

    double ehit_threshold = 0.5; // MeV
    double evox_threshold = 0.5; // MeV
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
            if(ix < nx && iy < ny) {
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

#if 1
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
#endif

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

    struct VOXEL
    {
        int x, y, z, layer;
        float adjust;
        float score;
    };
    std::vector<struct VOXEL> voxels;

    // Step 2: Iterative Refinement to match projections
    for (int iter = 0; iter < maxIter; ++iter) {

        for (int z = 0; z < nztot; ++z)
        {
            int nvox_layer = 0;
            for (int x = 0; x < nx; ++x)
                for (int y = 0; y < ny; ++y)
                {
                    if (V[x][y][z].value > 0)
                    {
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
        double total_score = 0;

//        verbose = 4; // FIXME: DEBUG

        // adjust XY
        for (int ilayer = 0; ilayer < nrep; ilayer++)
        {
            // count number of voxel in this module
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
                    total_score += fabs(difference);
                    if(verbose>3 && fabs(difference)>1e-3)
                        std::cout << " ilayer: " << ilayer << " difference=" << difference << std::endl;

                    if(fabs(difference) < 1e-3)continue;

                    voxels.clear();
                    for (int z = 0; z < nzlayer; z++) {
                        int izz = ilayer * nzlayer + z;
                        float adjust = XY[ilayer][x][y] - V[x][y][izz].value;
                        if(V[x][y][izz].value + adjust < 0) 
                            adjust = -V[x][y][izz].value;

                        double eXZ = XZ[x][izz];
                        double eYZ = YZ[y][izz];
                        double diff = std::max(fabs(eXZ - eYZ),0.1);

                        float score = fabs(adjust) + diff;
                        struct VOXEL vox = {x,y,izz,ilayer,adjust,score};
                        voxels.push_back(vox);
                    }

                    // Sorting by score in descending order
//                    std::sort(voxels.begin(), voxels.end(), [](const VOXEL &a, const VOXEL &b)
//                              { return a.score > b.score; });
//                    std::shuffle(voxels.begin(), voxels.end(), rng);

                    // Compute cumulative distribution (CDF)
                    std::vector<float> cdf(voxels.size());
                    cdf[0] = voxels[0].score; // Start with the first score
                    for (size_t i = 1; i < voxels.size(); ++i)
                    {
                        cdf[i] = cdf[i - 1] + voxels[i].score;
                    }
                    // Get the total sum of scores
                    float integ_score = cdf.back();
                    std::uniform_real_distribution<float> dist(0.0, integ_score);

                    int niter = 0;
                    while(fabs(difference)>1e-3 && niter < voxels.size()) {
                        float random_value = dist(rng);
                    // Use binary search (std::upper_bound) to find the voxel corresponding to the random value
                        auto it = std::upper_bound(cdf.begin(), cdf.end(), random_value);
                        int index = std::distance(cdf.begin(), it);
                        if(index >= voxels.size()) continue;
                        auto vox = voxels[index];
                        int x = vox.x;
                        int y = vox.y;
                        int izz = vox.z;
                        vox.adjust = XY[ilayer][x][y] - V[x][y][izz].value;  // recompute at each iteration
                        float adjust = std::min(difference, vox.adjust);
#if 0
                        if(fabs(difference) > fabs(vox.adjust)) {
                            adjust = vox.adjust;                            
                        } else {
                            adjust = difference;
                        }
#endif
                        if(verbose>3 && fabs(adjust)>1e-3)
                            std::cout << "    z: " << izz << "  adjust " << adjust << std::endl;
                        V[x][y][izz].value += adjust;
                        difference -= adjust;
                        if(adjust!=0) adjusted++;
                        niter++;
                    }
#if 0
                    for (const VOXEL& vox : voxels) {
                        if(abs(difference)<1e-3) break;
                        int x = vox.x;
                        int y = vox.y;
                        int izz = vox.z;
                        float adjust = std::min(difference, vox.adjust);
                        if(verbose>3 && abs(adjust)>1e-3)
                            std::cout << "    z: " << izz << "  adjust " << adjust << std::endl;
                        V[x][y][izz].value += adjust;
                        difference -= adjust;
                        if(adjust!=0) adjusted++;
                    }
#endif

#if 0
                    int niter = 0;
                    while(abs(difference)>1e-3 && niter++ < 2*nzlayer) {
                        int z = rnd_layer(gen);

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
#endif
                }
            }
        }

        // adjust XZ and YZ
        for (int z = 0; z < nztot; ++z) {
            if(nvox_per_layer[z] > nvox_per_layer_max) continue;

            long ilayer = z / nrep;

            // Adjust XZ projection
            for (int x = 0; x < nx; ++x) {
                double sumXZ = 0;
                for (int y = 0; y < ny; ++y) {
                    sumXZ += V[x][y][z].value;
                }
                float difference = XZ[x][z] - sumXZ;
                total_score += fabs(difference);
                if(verbose>3 && fabs(difference)>1e-3) std::cout << " x: " << x << " difference=" << difference << std::endl;

#if 0
                if (fabs(difference) < 1e-3)
                    continue;

                voxels.clear();
                for (int y = 0; y < ny; y++)
                {
                    float adjust = YZ[y][z] - V[x][y][z].value;
                    if (V[x][y][z].value + adjust < 0)
                        adjust = -V[x][y][z].value;
                    float score = fabs(adjust);
                    if(score>0){
                        struct VOXEL vox = {x, y, z, 0, adjust, score};
                        voxels.push_back(vox);
                    }
                }
#endif
#if 0
                // Sorting by score in descending order
                std::sort(voxels.begin(), voxels.end(), [](const VOXEL &a, const VOXEL &b)
                          { return a.score > b.score; });
#endif
//                std::shuffle(voxels.begin(), voxels.end(), rng);
#if 0
                bool stop = false;
                while (!stop && !voxels.empty())
                {
                    // Compute cumulative distribution (CDF)
                    std::vector<float> cdf(voxels.size());
                    cdf[0] = voxels[0].score; // Start with the first score
                    for (size_t i = 1; i < voxels.size(); ++i)
                    {
                        cdf[i] = cdf[i - 1] + voxels[i].score;
                    }
                    // Get the total sum of scores
                    float integ_score = cdf.back();
                    std::uniform_real_distribution<float> dist(0.0, integ_score);

                    float random_value = dist(rng);
                    // Use binary search (std::upper_bound) to find the voxel corresponding to the random value
                    auto it = std::upper_bound(cdf.begin(), cdf.end(), random_value);
                    int index = std::distance(cdf.begin(), it);
                    auto vox = voxels[index];
                    int x = vox.x;
                    int y = vox.y;
                    int z = vox.z;
                    float adjust = std::min(difference, vox.adjust);
                    V[x][y][z].value += adjust;
                    difference -= adjust;
                    if (adjust != 0)
                        adjusted++;
                    voxels.erase(voxels.begin()+index);
//                    stop = true;
                }
#endif
#if 1
                int niter = 0;
                while(fabs(difference)>1e-3 && niter++ < ny) {
                    int y = rnd_ny(gen);
                    float adjust = std::min(difference, YZ[y][z] - V[x][y][z].value);
                    if(V[x][y][z].value + adjust < 0) adjust = -V[x][y][z].value;
                    if(verbose>3 && abs(adjust)>1e-3) std::cout << "    y: " << y << "  adjust " << adjust << std::endl;
                    V[x][y][z].value += adjust;
                    difference -= adjust;
                    if(adjust!=0) adjusted++;
                }
#endif
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
                total_score += fabs(difference);
                if(verbose>3 && fabs(difference)>1e-3)std::cout << " y: " << y << " difference=" << difference << std::endl;
                int niter = 0;
                while(fabs(difference)>1e-3 && niter++ < nx) {
                    int x = rnd_nx(gen);
                    float adjust = std::min(difference, XZ[x][z] - V[x][y][z].value);
                    if(V[x][y][z].value + adjust < 0) adjust = -V[x][y][z].value;
                    if(verbose>3 && fabs(adjust)>1e-3) std::cout << "    x: " << x << "  adjust " << adjust << std::endl;
                    V[x][y][z].value += adjust;
                    difference -= adjust;
                    if(adjust!=0) adjusted++;
                }
            }
        }

        if(verbose>1 && adjusted>0) std::cout << "Iteration " << iter << ": " << adjusted <<
        " voxels adjusted. Score = " << total_score << std::endl;

#if 0
        // do some cleanup removing hits with no XY hit
        for (int z = 0; z < nztot; ++z) {
            long ilayer = z / nrep;
            for (int x = 0; x < nx; ++x) {
                for (int y = 0; y < ny; ++y) {
                    if(V[x][y][z].value > 1e-3) {
                        if(XY[ilayer][x][y] < 0.5) {
//                            V[x][y][z].value = 0;
                        }
                    }
                    double neigh = 0;
                    if(x>0) neigh += V[x-1][y][z].value;
                    if(x<nx-1) neigh += V[x+1][y][z].value;
                    if(y>0) neigh += V[x][y-1][z].value;
                    if(y<ny-1) neigh += V[x][y+1][z].value;
                    if(neigh<0.5) {
                            V[x][y][z].value = 0;
                    }
                }

            }
        }
#endif

        // save solution of with minimum score
        if(total_score < min_score){
        for (int z = 0; z < nztot; ++z)
            for (int x = 0; x < nx; ++x)
                for (int y = 0; y < ny; ++y)
                {
                    V_min[x][y][z].value = V[x][y][z].value;
                }
            min_score = total_score;
        }

        // if score degraded, switch back to best answer
        if(total_score > prev_score) {
        for (int z = 0; z < nztot; ++z)
            for (int x = 0; x < nx; ++x)
                for (int y = 0; y < ny; ++y)
                {
                    V[x][y][z].value = V_min[x][y][z].value;
                }
        }
        prev_score = total_score;

        if(total_score == 0) break;
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

                // cut on minimum energy
                if(ehit<evox_threshold)
                    continue;

#if 0
                // remove totally isolated hits
                double neigh = 0;
                if (x > 0)
                    neigh += V_min[x - 1][y][z].value;
                if (x < nx - 1)
                    neigh += V_min[x + 1][y][z].value;
                if (y > 0)
                    neigh += V_min[x][y - 1][z].value;
                if (y < ny - 1)
                    neigh += V_min[x][y + 1][z].value;
                if (neigh < 1.0)
                    continue;
#endif

                long ilayer = z / nzlayer;
#if 0
                double eXZ = XZ[x][z];
                double eYZ = YZ[y][z];
                double eXY = XY[ilayer][x][y];
                double diff = fabs(eXZ - eYZ);
                if(diff > 10.0) continue;
#endif
                long iz = z % nzlayer;
                long ID = x + y * 1000 + iz * 1000000 + ilayer * 1000000000;
                struct PSVOXEL3D v = {ID, ehit, true};
                PSvoxelmap[ID] = v;
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
    TH1D h_XY_fake = TH1D("h_XY_fake", "XY hit fakes", 100, 0., 30.);
    TH1D h_XY_real = TH1D("h_XY_real", "XY hit real", 100, 0., 30.);

    for (const auto& v : PSvoxelmap) {
        ntot++;
        long ID = v.first;
        long ix = ID % 1000;
        long iy = (ID / 1000) % 1000;
        long iz = (ID / 1000000) % 1000;
        long ilayer = fTcalEvent->getChannelModulefromID(ID);
        if(ilayer > maxLayer) continue;
        long izz = ilayer * nzlayer + iz;
        double eXZ = XZ[ix][izz];
        double eYZ = YZ[iy][izz];
        double eXY = XY[ilayer][ix][iy];
        double diff = fabs(eXZ - eYZ);
//        diff = v.second.RawEnergy;

        if(v.second.ghost) {
            fakes++;
            ave_e_ghost += v.second.RawEnergy;
            h_diff_fake.Fill(diff);
            h_XY_fake.Fill(eXY);
        } else {
            ave_e_real += v.second.RawEnergy;
            nreal++;
            h_diff_real.Fill(diff);
            h_XY_real.Fill(eXY);
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
//    h_XY_real.Draw();
    h_diff_real.Draw();
    c1->cd(2);
    h_diff_fake.Draw();
//    h_XY_fake.Draw();
    c1->Modified();
    c1->Update();
    c1->SaveAs("3dreco.png");
#endif
}

void TPORecoEvent::reconstruct3DPS_module(int maxIter, int imodule, std::vector<std::vector<std::vector<Voxel>>> &V,
    std::vector<std::vector<float>> &XZ, std::vector<std::vector<float>> &YZ, std::vector<std::vector<std::vector<float>>> &XY,
    std::vector<int>& nvox_per_layer, int nvox_per_layer_max) {

    int nvox_max_after_iteration = 25;  // after this iteration limit the number of voxels in module

    int nx = fTcalEvent->geom_detector.fScintillatorSizeX / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int ny = fTcalEvent->geom_detector.fScintillatorSizeY / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int nzlayer = fTcalEvent->geom_detector.fSandwichLength / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int nrep = fTcalEvent->geom_detector.NRep;
    int nztot =  nrep * nzlayer;

    std::random_device rd;  // Seed for the random number generator
    std::mt19937 gen(rd());  // Mersenne Twister random number generator
    std::uniform_int_distribution<> rnd_nx(0, nx-1);
    std::uniform_int_distribution<> rnd_ny(0, ny-1);

    double total_score_min_break = 10;

    double total_score = 0;
    for (int iter = 0; iter < maxIter; ++iter)
    {
        int adjusted = 0;
        total_score = 0;
        // decide which layers should be used for reconstructing 3D voxels
        for (int z = imodule * nzlayer; z < (imodule + 1) * nzlayer; ++z)
        {
            int nvox_layer = 0;
            for (int x = 0; x < nx; ++x)
                for (int y = 0; y < ny; ++y)
                {
                    if (V[x][y][z].value > 0)
                    {
                        nvox_layer++;
                    }
                }
            nvox_per_layer[z] = nvox_layer;
        }

        // count number of voxel in this module
        int sum_nvox = 0;
        for (int z = 0; z < nzlayer; z++)
        {
            int izz = imodule * nzlayer + z;
            sum_nvox += nvox_per_layer[izz];
        }
        if (iter > nvox_max_after_iteration && sum_nvox > nvox_per_layer_max)
            continue;

        double module_score = 0;
        for (int iz = 0; iz < nzlayer; ++iz)
        {

            int z = imodule * nzlayer + iz;

#if 0
                // THIS MAKES RESULTS WORSE... more fake hits
                // if first iteration, copy result from previous layer as starting point
                if(iter ==0 && z > 0) {
                for (int x = 0; x < nx; ++x)
                    for (int y = 0; y < ny; ++y)
                    {
                        V[x][y][z].value = V[x][y][z-1].value;
                    }
                }
#endif

            for (int iterl = 0; iterl < nx * ny; ++iterl)
            {

                int x = rnd_nx(gen);
                int y = rnd_ny(gen);

                double sumXZ = 0;
                for (int iy = 0; iy < ny; ++iy)
                {
                    sumXZ += V[x][iy][z].value;
                }
                double sumYZ = 0;
                for (int ix = 0; ix < nx; ++ix)
                {
                    sumYZ += V[ix][y][z].value;
                }
                long ilayer = imodule; // z / nzlayer;
                double sumXY = 0;
                for (int iz = 0; iz < nzlayer; iz++)
                {
                    int izz = ilayer * nzlayer + iz;
                    sumXY += V[x][y][izz].value;
                }

                double tolerance = 1.0;
                double learning_rate = 0.01;
                int max_iters = 1000;
                double prev_value = 0;

                double adjust = 0;

                for (int itermin = 0; itermin < max_iters; itermin++)
                {
                    double a1 = sumXY - XY[ilayer][x][y];
                    double a2 = sumXZ - XZ[x][z];
                    double a3 = sumYZ - YZ[y][z];
#define SQR(a) ((a) * (a))
                    //                        double chi2 = pow(a1 + adjust, 2) + pow(a2 + adjust, 2) + pow(a3 + adjust, 2);
                    double chi2 = SQR(a1 + adjust) + SQR(a2 + adjust) + SQR(a3 + adjust);

                    //                   std::cout << itermin << "- adjust: " << adjust << " chi2: " << chi2 << std::endl;

                    if (std::abs(chi2 - prev_value) < tolerance)
                        break;
                    double adjust2 = adjust + tolerance;
                    //                        double chi2_2 = pow(a1 + adjust2, 2) + pow(a2 + adjust2, 2) + pow(a3 + adjust2, 2);
                    double chi2_2 = SQR(a1 + adjust2) + SQR(a2 + adjust2) + SQR(a3 + adjust2);
                    double gradient = (chi2_2 - chi2) / tolerance;
                    adjust -= learning_rate * gradient;
                    prev_value = chi2;
                }

                if (std::abs(adjust) == 0)
                    continue;

                V[x][y][z].value += adjust;
                total_score += std::abs(adjust);
                module_score += std::abs(adjust);
                if (V[x][y][z].value < 0)
                    V[x][y][z].value = 0;
                adjusted++;
            }
        }
        if (verbose > 4)
        {
            std::cout << " module " << imodule << " has " << sum_nvox << " voxels " << " score " << module_score << std::endl;
            std::cout << "Module " << imodule << " - Iteration " << iter << ": " << adjusted << " voxels adjusted. Score = " << total_score << std::endl;
        }
        if (adjusted == 0 || total_score < total_score_min_break)
            break;
    }
    if (verbose > 1)
        std::cout << "Module " << imodule << " - Score = " << total_score << std::endl;
}

void TPORecoEvent::Reconstruct3DPS_2(int maxIter) {

    maxIter = 150; // 300;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine rng(seed);

    double ehit_threshold = 0.5; // MeV
    double evox_threshold = 0.5; // MeV

    int nvox_per_layer_max = 3000;      // maximum number of voxels in module

    int nx = fTcalEvent->geom_detector.fScintillatorSizeX / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int ny = fTcalEvent->geom_detector.fScintillatorSizeY / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int nzlayer = fTcalEvent->geom_detector.fSandwichLength / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int nrep = fTcalEvent->geom_detector.NRep;
    int nztot =  nrep * nzlayer;

    std::uniform_int_distribution<> rnd_layer(0, nzlayer-1);

    // the maximum number layer (from 0 to nRep) that is reconstructed
    int maxLayer = 25;

    // Step 0: organize hits for easy access
    std::vector<std::vector<float>> XZ(nx, std::vector<float>(nztot, 0.0));
    std::vector<std::vector<float>> YZ(ny, std::vector<float>(nztot, 0.0));
    std::vector<std::vector<std::vector<float>>> XY(nrep,
            std::vector<std::vector<float>>(nx, std::vector<float>(ny, 0.0)));

    if(verbose>0)
       std::cout << "Start reconstruction of 3D voxels..." << std::endl;

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
            if(ix < nx && iy < ny) {
                XY[ilayer][ix][iy] = ehit;
            }
        }
    }

    std::vector<std::vector<std::vector<Voxel>>> V(
        nx, std::vector<std::vector<Voxel>>(
                   ny, std::vector<Voxel>(
                               nztot, Voxel())));

    std::vector<std::vector<std::vector<Voxel>>> V_min(
        nx, std::vector<std::vector<Voxel>>(
                   ny, std::vector<Voxel>(
                               nztot, Voxel())));

    std::vector<int> nvox_per_layer(nztot);

    std::vector<std::thread> threads;
    for (int imodule = 0; imodule < nrep; imodule++)
    {
 //       reconstruct3DPS_module(maxIter, imodule, V, XZ, YZ, XY, nvox_per_layer, nvox_per_layer_max);
        threads.emplace_back(&TPORecoEvent::reconstruct3DPS_module, this, maxIter, imodule, 
            std::ref(V), std::ref(XZ), std::ref(YZ), std::ref(XY),
            std::ref(nvox_per_layer), nvox_per_layer_max);
    }

    for (auto& th : threads)
    {
        th.join();
    }

    PSvoxelmap.clear();
    for (int z = 0; z < nztot; ++z) {
        //
        if(nvox_per_layer[z] > nvox_per_layer_max) continue;
        //
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                float ehit = V[x][y][z].value;

                // cut on minimum energy
                if(ehit<evox_threshold)
                    continue;

                // remove totally isolated hits
                double neigh = 0;
                if (x > 0)
                    neigh += V[x - 1][y][z].value;
                if (x < nx - 1)
                    neigh += V[x + 1][y][z].value;
                if (y > 0)
                    neigh += V[x][y - 1][z].value;
                if (y < ny - 1)
                    neigh += V[x][y + 1][z].value;
                if (z > 0) 
                    neigh += V[x][y][z-1].value;
                if (z < nztot - 1)
                    neigh += V[x][y][z+1].value;

                if (neigh < 1.0)
                    continue;

                long ilayer = z / nzlayer;
                long iz = z % nzlayer;
                long ID = x + y * 1000 + iz * 1000000 + ilayer * 1000000000;
                struct PSVOXEL3D v = {ID, ehit, true};
                PSvoxelmap[ID] = v;
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
    TH1D h_XY_fake = TH1D("h_XY_fake", "XY hit fakes", 100, 0., 30.);
    TH1D h_XY_real = TH1D("h_XY_real", "XY hit real", 100, 0., 30.);

    int ntotl[nrep];
    int nfake[nrep];
    for(int imodule = 0; imodule < nrep; imodule++){
        ntotl[imodule] = 0;
        nfake[imodule] = 0;
    }
    for (const auto& v : PSvoxelmap) {
        ntot++;
        long ID = v.first;
        long ix = ID % 1000;
        long iy = (ID / 1000) % 1000;
        long iz = (ID / 1000000) % 1000;
        long ilayer = fTcalEvent->getChannelModulefromID(ID);
        ntotl[ilayer]++;
        if(ilayer > maxLayer) continue;
        long izz = ilayer * nzlayer + iz;
        double eXZ = XZ[ix][izz];
        double eYZ = YZ[iy][izz];
        double eXY = XY[ilayer][ix][iy];
        double diff = fabs(eXZ - eYZ);
//        diff = v.second.RawEnergy;

        if(v.second.ghost) {
            fakes++;
            nfake[ilayer]++;
            ave_e_ghost += v.second.RawEnergy;
            h_diff_fake.Fill(diff);
            h_XY_fake.Fill(eXY);
        } else {
            ave_e_real += v.second.RawEnergy;
            nreal++;
            h_diff_real.Fill(diff);
            h_XY_real.Fill(eXY);
        }
    }
    for(int imodule = 0; imodule < nrep; imodule++){
        double frac = ntotl[imodule] > 0 ? nfake[imodule]*100.0/float(ntotl[imodule]) : -1;
        std::cout << " Module " << imodule << " " << ntotl[imodule] << " hits " << nfake[imodule] << " fakes ";
        std::cout << frac << "%" << std::endl;
    }
    std::cout << " STATS: " << ntot << " hits " << fakes << " ghosts." << std::endl;
    std::cout << " Avg energy: " << ave_e_real/float(nreal) << " ghosts: " << ave_e_ghost/float(fakes);
    std::cout << std::endl;
}


void TPORecoEvent::Reconstruct3DPS_Eflow() {
    struct TPORec::CALENERGIES result;
    double cogx = 0;
    double cogy = 0;
    double cogz = 0;
    double etot = 0;
    for (const auto& v : PSvoxelmap) {
        long ID = v.first;
        ROOT::Math::XYZVector position = fTcalEvent -> getChannelXYZfromID(ID);
        double ehit = v.second.RawEnergy;
        cogx += position.X()*ehit;
        cogy += position.Y()*ehit;
        cogz += position.Z()*ehit;
        etot += ehit;
    }
    if(etot>0) {
        result.cog.SetX(cogx/etot);
        result.cog.SetY(cogy/etot);
        result.cog.SetZ(cogz/etot);

        // compute energy flow relative to primary vertex
        ROOT::Math::XYZVector primary(fTPOEvent->prim_vx.X(), fTPOEvent->prim_vx.Y(), fTPOEvent->prim_vx.Z());
        ROOT::Math::XYZVector direction = result.cog - primary;
        result.Eflow = (etot/1e3)*direction.Unit(); // convert to GeV

    } else {
        result.cog.SetCoordinates(0,0,0);
        result.Eflow.SetCoordinates(0,0,0);
    }
    fPOFullRecoEvent = new TPORec(-1);
    fPOFullRecoEvent->fTotal.cog = result.cog;
    fPOFullRecoEvent->fTotal.Eflow = result.Eflow;
    fPOFullRecoEvent->fTotal.Ecompensated = etot / 1e3; // in GeV
}

void TPORecoEvent::PSVoxelParticleFilter() {

    int max_number_track_seeds = 100;

    if(verbose>0) std::cout << "Start PS voxel particle filter..." << std::endl;

    int nmodules = fTcalEvent->geom_detector.NRep;
    int nzlayer = fTcalEvent->geom_detector.fSandwichLength / fTcalEvent->geom_detector.fScintillatorVoxelSize;

    double closest_voxel_cut = fTcalEvent->geom_detector.fScintillatorVoxelSize*2.0;
    closest_voxel_cut = closest_voxel_cut*closest_voxel_cut;

    std::map<int, std::vector<struct PSVOXEL3D>> PSvoxelmapModule;
    // organize voxels by module number
    for (const auto &it : PSvoxelmap) {
        long ID = it.first;
        long module = fTcalEvent->getChannelModulefromID(ID);
        auto lay = PSvoxelmapModule.find(module);
        if(lay != PSvoxelmapModule.end()) {
            lay->second.push_back(it.second);
        } else {
            std::vector <struct PSVOXEL3D> voxels;
            voxels.push_back(it.second);
            PSvoxelmapModule[module] = voxels;
        }
    }

    std::vector<TPSTrack> track_seeds;
    track_seeds.clear();

    // loop over modules
    for (int i = 0; i < nmodules; i++) {
        if(verbose>1) std::cout << " ------------------- Module " << i << " - track seeds: " << track_seeds.size() << std::endl;
        //
        if(track_seeds.size() > max_number_track_seeds)
            continue;
        //
        auto module = PSvoxelmapModule.find(i);
        if(module == PSvoxelmapModule.end()) continue;
        // loop over layers within module in reverse order
        for (int layer = nzlayer-1; layer > 0; layer--) {
            for (const auto &hit : module->second) {
                long iz = (hit.ID / 1000000) % 1000;
                ROOT::Math::XYZVector position = fTcalEvent -> getChannelXYZfromID(hit.ID);
                if(iz != layer) continue;

                TPSTrack track_seed;
                TPSTrack *track;

                // check if this voxel already belong to a track seed
                bool found_seed = false;
                for (auto &ts : track_seeds) {
                    struct TPSTrack::TRACKHIT *voxel = ts.VoxelBelongstoTrack(hit.ID);
                    if(voxel != nullptr) {
                        track = &ts;
                        found_seed = true;
                        break;
                    }
                }

                // check if it's touching a track
                if (!found_seed)
                {
                    for (auto &ts : track_seeds)
                    {
                        if (ts.VoxelTouchesTrack(hit.ID))
                        {
                            track = &ts;
                            found_seed = true;
                            break;
                        }
                    }
                }

                if(!found_seed) {
                    // create PS track seed
                    track = &track_seed;
                    struct TPSTrack::TRACKHIT vox = {hit.ID, position, hit.RawEnergy};
                    track->tkhit.push_back(vox);
                }

                int prev_layer = layer - 1;
                // find closest voxel in previous layer
                double transv_dist = 1e9;
                struct PSVOXEL3D *closest_vox = nullptr;
                for (auto &prev_hit : module->second)
                {
                    long prev_iz = (prev_hit.ID / 1000000) % 1000;
                    if (prev_iz != prev_layer)
                        continue;
                    ROOT::Math::XYZVector prev_position = fTcalEvent->getChannelXYZfromID(prev_hit.ID);
                    ROOT::Math::XYZVector diff = position - prev_position;
                    double d2 = diff.x() * diff.x() + diff.y() * diff.y();
                    if (d2 < transv_dist)
                    {
                        transv_dist = d2;
                        closest_vox = &prev_hit;
                    }
                }
//                std::cout << track->direction.x() << " " << track->direction.x()*track->direction.x() << std::endl;
//                std::cout << track->direction.y() << " " << track->direction.y()*track->direction.y() << std::endl;
                double rt = track->direction.x()*track->direction.x() + track->direction.y()*track->direction.y();
//                std::cout << rt << std::endl;
                double tanangle = sqrt(rt) / (std::max(track->direction.z(),1e-3));
//                std::cout << tanangle << std::endl;
                double cut = closest_voxel_cut*std::max(1.0,tanangle);
//                std::cout << cut << std::endl;
#if 0
                if(rt>0) {
                    std::cout << " track angle found - ";
                    track->Dump();
                }
#endif
                if (transv_dist < cut)
                {
                    ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(closest_vox->ID);
                    struct TPSTrack::TRACKHIT vox2 = {closest_vox->ID, position, closest_vox->RawEnergy};
                    track->tkhit.push_back(vox2);
                    int nhits = track->tkhit.size();
                    track->direction.SetXYZ(0,0,0);
                    track->SortHitsByZ();
                    if (nhits == 2)
                        track->direction = track->direction2Hits(track->centroid);
                    if (nhits > 2)
                        track->direction = track->fitLineThroughHits(track->centroid);
                }
                if(!found_seed)
                    track_seeds.push_back(*track);
            }
        }
        // we should now have all track seeds from the given layer
        for (auto &tk : track_seeds)
        {
            int nhits = tk.tkhit.size();
            tk.direction.SetXYZ(0,0,0);
            tk.SortHitsByZ();
            if(nhits==2) tk.direction = tk.direction2Hits(tk.centroid);
            if(nhits>2) tk.direction = tk.fitLineThroughHits(tk.centroid);
        }

        // now loop over tracks and merge segments
        double parallel_cut = 0.01; // FIXME: adjust value
        double mindZcut = fTcalEvent->geom_detector.fSandwichLength;

        for (size_t i = 0; i < track_seeds.size(); i++)
        {
            TPSTrack &track1 = track_seeds[i];
            size_t besttrk = -1;
            double mindZ = 1e9;
            for (size_t j = 0; j < track_seeds.size(); j++)
            {
                if (i == j)
                    continue;
                TPSTrack &track2 = track_seeds[j];
                // check if segments are parallel
                TVector3 normDir1 = track1.direction.Unit();
                TVector3 normDir2 = track2.direction.Unit();
                double dotProduct = normDir1.Dot(normDir2);
                if (std::abs(std::abs(dotProduct) - 1.0) > parallel_cut)
                    continue;
                // ensure that segments belong to different planes
                struct TPSTrack::TRACKHIT hit1 = track1.tkhit.back();
                struct TPSTrack::TRACKHIT hit2 = track2.tkhit.front();
                double dz = std::abs(hit1.point.z() - hit2.point.z());
                if (dz < mindZcut)
                    continue;
                // ensure that segments are parallel to main line joining hits
                ROOT::Math::XYZVector normDir = (hit2.point - hit1.point).Unit();
                dotProduct = std::min(std::abs(normDir.Dot(normDir1)), std::abs(normDir.Dot(normDir2)));
                if (std::abs(dotProduct - 1.0) > parallel_cut)
                    continue;
                if (dz < mindZ)
                {
                    mindZ = dz;
                    besttrk = j;
                }
            }
            // if match found, then merge the tracks
            if (besttrk != -1)
            {
                // add hits from track2 to track1
                for (const auto &h : track_seeds[besttrk].tkhit)
                {
                    track1.tkhit.push_back(h);
                }
                track1.direction = track1.fitLineThroughHits(track1.centroid);
                track_seeds.erase(track_seeds.begin() + besttrk);
                if (besttrk < i)
                    i--;
            }
        }

    }

    // always order hits in the track by increasing "z"
    for (auto &trk : track_seeds)
    {
        trk.SortHitsByZ();
    }

    // copy into permanent storage
    fPSTracks.clear();
    for (auto &trk : track_seeds) {
        if(trk.tkhit.size()<2) continue;               // ignore isolated hits
        fPSTracks.push_back(trk);
        if (verbose > 2)
            trk.Dump();
    }
}

void TPORecoEvent::ReconstructRearCals() {
    double eDeposit = 0;
    rearCals.rearCalModule.clear();
    for (const auto &it : fTcalEvent->rearCalDeposit) {
        eDeposit += it.energyDeposit;
        struct TcalEvent::REARCALDEPOSIT im = {it.moduleID, it.energyDeposit};
        rearCals.rearCalModule.push_back(im);
    }
    rearCals.rearCalDeposit = eDeposit/1e3;  // convert to GeV
    rearCals.rearMuCalDeposit = fTcalEvent->rearMuCalDeposit; // in MeV
    if(verbose>0){
        std::cout << "Rear CalDeposit " << rearCals.rearCalDeposit << std::endl;
        std::cout << "Rear MuCalDeposit " << rearCals.rearMuCalDeposit << std::endl;
    }
}
