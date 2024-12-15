#include "TPORecoEvent.hh"
#include "DBScan.hh"
#include "TTKTrack.hh"
#include "TPSTrack.hh"

#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <TDecompSVD.h>
#include <TVector3.h>
#include <TCanvas.h>

// GENFIT
#include <GFRaveVertexFactory.h>
#include <GFRaveVertex.h>
#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <MaterialEffects.h>
#include <TGeoMaterialInterface.h>

#include <random>
#include <thread>

#include <set>

#include <signal.h>

ClassImp(TPORec);
ClassImp(TPORecoEvent);

static int TPORecoEvent_initGenfit = 0;
static int TPORecoEvent_configPrinted = 0;

// Function to print a table of RECOCONFIG values
static void printRecoConfig(const struct TPORecoEvent::RECOCONFIG& config) {
    std::cout << std::left << std::setw(40) << "Variable" 
              << std::setw(20) << "Value" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    // Print each variable and its value
    std::cout << std::setw(40) << "psvoxel_fudge_factor" << config.psvoxel_fudge_factor << std::endl;
    std::cout << std::setw(40) << "alpha" << config.alpha << std::endl;
    std::cout << std::setw(40) << "beta" << config.beta << std::endl;

    std::cout << std::setw(40) << "findpattern_max_hit_layers" << config.findpattern_max_hit_layers << std::endl;
    std::cout << std::setw(40) << "findpattern_dist_min_cut" << config.findpattern_dist_min_cut << std::endl;
    std::cout << std::setw(40) << "findpattern_parallel_cut" << config.findpattern_parallel_cut << std::endl;
    std::cout << std::setw(40) << "findpattern_mindZ_fudge" << config.findpattern_mindZ_fudge << std::endl;
    std::cout << std::setw(40) << "findpattern_parallel_cut_merge" << config.findpattern_parallel_cut_merge << std::endl;

    std::cout << std::setw(40) << "extendtracks_closest_voxel_cut" << config.extendtracks_closest_voxel_cut << std::endl;
    std::cout << std::setw(40) << "extendtracks_dist2_perp_voxel_cut" << config.extendtracks_dist2_perp_voxel_cut << std::endl;

    std::cout << std::setw(40) << "genfit_min_pVal" << config.genfit_min_pVal << std::endl;
    std::cout << std::setw(40) << "genfit_min_pMom" << config.genfit_min_pMom << std::endl;

    std::cout << std::setw(40) << "findvtx_cut_max_trk" << config.findvtx_cut_max_trk << std::endl;
    std::cout << std::setw(40) << "findvtx_chi2ndf_cut" << config.findvtx_chi2ndf_cut << std::endl;
    std::cout << std::setw(40) << "findvtx_trk_dist_cut" << config.findvtx_trk_dist_cut << std::endl;
    std::cout << std::setw(40) << "findvtx_merge_dist_cut" << config.findvtx_merge_dist_cut << std::endl;

    std::cout << std::setw(40) << "clusters_threshold_2dhit" << config.clusters_threshold_2dhit << std::endl;
    std::cout << std::setw(40) << "clusters_eps" << config.clusters_eps << std::endl;
    std::cout << std::setw(40) << "clusters_minPts" << config.clusters_minPts << std::endl;
    std::cout << std::setw(40) << "clusters_threshold_cluster" << config.clusters_threshold_cluster << std::endl;

    std::cout << std::setw(40) << "PS3D_nvox_max_after_iteration" << config.PS3D_nvox_max_after_iteration << std::endl;
    std::cout << std::setw(40) << "PS3D_total_score_min_break" << config.PS3D_total_score_min_break << std::endl;
    std::cout << std::setw(40) << "PS3D_ehit_threshold" << config.PS3D_ehit_threshold << std::endl;
    std::cout << std::setw(40) << "PS3D_evox_threshold" << config.PS3D_evox_threshold << std::endl;
    std::cout << std::setw(40) << "PS3D_nvox_per_layer_max" << config.PS3D_nvox_per_layer_max << std::endl;

    std::cout << std::setw(40) << "PSFilter_max_number_track_seeds" << config.PSFilter_max_number_track_seeds << std::endl;
    std::cout << std::setw(40) << "PSFilter_closest_voxel_cut" << config.PSFilter_closest_voxel_cut << std::endl;
    std::cout << std::setw(40) << "PSFilter_parallel_cut" << config.PSFilter_parallel_cut << std::endl;
    std::cout << std::setw(40) << "PSFilter_mindZcut" << config.PSFilter_mindZcut << std::endl;
    std::cout << std::string(60, '-') << std::endl;
}

TPORecoEvent:: TPORecoEvent() : TObject(), fTcalEvent(0), fTPOEvent(0), fPOFullEvent(0), fPOFullRecoEvent(0) {
    // initialize genfit
    if(!TPORecoEvent_initGenfit) {
        TPORecoEvent_initGenfit = 1;
        std::cout << "Initializing Genfit" << std::endl;
        genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
        genfit::FieldManager::getInstance()->init(new genfit::ConstField(0. ,1e-4, 0.)); 
    }
}

TPORecoEvent::TPORecoEvent(TcalEvent* c, TPOEvent* p) : TPORecoEvent() {
    fTcalEvent = c;
    fTPOEvent = p;
    // copy geometry
    geom_detector = fTcalEvent->geom_detector;

    // empty histogram pointers
	for(int i =0; i < 50; i++){
		zviewPS.push_back(nullptr);
	}

        // set default values for the recoConfig
    recoConfig.psvoxel_fudge_factor = 1.25/sqrt(12);  // fudge factor for voxel size

    recoConfig.alpha = 1.0/(1.0-0.341)*0.98;
    recoConfig.beta = 3.0;

    recoConfig.findpattern_max_hit_layers = 1000;
    recoConfig.findpattern_dist_min_cut = 15.0;
    recoConfig.findpattern_parallel_cut = 0.01;
    recoConfig.findpattern_mindZ_fudge = 1.1;
//    recoConfig.findpattern_cut_SSR_merge = 10.0;
    recoConfig.findpattern_parallel_cut_merge = 0.01;

    recoConfig.extendtracks_closest_voxel_cut = fTcalEvent->geom_detector.fScintillatorVoxelSize*2.0;
    recoConfig.extendtracks_dist2_perp_voxel_cut = fTcalEvent->geom_detector.fScintillatorVoxelSize*fTcalEvent->geom_detector.fScintillatorVoxelSize; 

    recoConfig.genfit_min_pVal = 0.01;
    recoConfig.genfit_min_pMom = 1e-3;

    recoConfig.findvtx_cut_max_trk = 100;
    recoConfig.findvtx_chi2ndf_cut = 10000;  // FIXME: tune value
    recoConfig.findvtx_trk_dist_cut = 1000.0; // in millimeters - FIXME: tune value
    recoConfig.findvtx_merge_dist_cut = 10;  // in centimeters - FIXME: tune value

    recoConfig.clusters_threshold_2dhit = 2.0; // MeV
    recoConfig.clusters_eps = 5; // in mm
    recoConfig.clusters_minPts = 10; // minimum number of points to form a cluster
    recoConfig.clusters_threshold_cluster = 10*1e3; // MeV

    recoConfig.PS3D_nvox_max_after_iteration = 25; // after this iteration limit the number of voxels in module
    recoConfig.PS3D_total_score_min_break = 10.0;
    recoConfig.PS3D_ehit_threshold = 0.5; // MeV
    recoConfig.PS3D_evox_threshold = 0.5; // MeV
    recoConfig.PS3D_nvox_per_layer_max = 3000; // maximum number of voxels per module

    recoConfig.PSFilter_max_number_track_seeds = 1000;
    recoConfig.PSFilter_closest_voxel_cut = fTcalEvent->geom_detector.fScintillatorVoxelSize*2.0;

    recoConfig.PSFilter_parallel_cut = 0.01;
    recoConfig.PSFilter_mindZcut = fTcalEvent->geom_detector.fSandwichLength;

    if(!TPORecoEvent_configPrinted) {
        printRecoConfig(recoConfig);
        TPORecoEvent_configPrinted = 1;
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

    if(!TPORecoEvent_configPrinted) {
        std::cerr << "TPORecoEvent::Reconstruct - configuration not initialized!" << std::endl;
        exit(1);
    }

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
        it->fTotal.Ecompensated = it->fTotal.em*recoConfig.alpha+it->fTotal.had*recoConfig.beta;

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

    if(!TPORecoEvent_configPrinted) {
        std::cerr << "TPORecoEvent::TrackReconstruct - configuration not initialized!" << std::endl;
        exit(1);
    }

    if(verbose > 0) {
        std::cout << "TPORecoEvent::TrackReconstruct - start" << std::endl;
    }
    FindPatternTracks();
    if(verbose > 2) {
        DumpReconstructedTracks();
    }
    FindTrackVertices();
    ExtendTracks();

    std::vector<TTKTrack> tempTracks;
    for (auto &trk : fTKTracks) {
        if(trk.tkhit.size() < 3) continue;
        // reject tracks with bad chi2
        if(trk.fitTrack->getFitStatus()->getChi2() > 0 && trk.fitTrack->getFitStatus()->getPVal() < recoConfig.genfit_min_pVal) continue;
        tempTracks.push_back(trk);
    }
    fTKTracks.clear();
    for (auto &trk : tempTracks) {
        fTKTracks.push_back(trk);
    }
}

void TPORecoEvent::FindPatternTracks() {
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
            if(ehit < 1e-3) {
 //               std::cout << "TrackReconstruct: Energy deposit is too low: " << ehit << std::endl;
                continue;
            }
            int ilayer = fTcalEvent->getChannelModulefromID(ID);
            ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(track->fhitIDs[i]);
            struct TTKTrack::TRACKHIT hit = {ID, (int)hittype, position, ehit};
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

    // dump the number of hits in each layer
    if(verbose>3) {
        for(auto it : hitMap) {
            std::cout << "Layer: " << it.first << " Hits: " << it.second.size() << std::endl;
        }
    }

    std::vector<TTKTrack*> tempTracks;

    // now match doublets in each layer
    for(auto &it : hitMap) {
        int ilayer = it.first;

        // skip layers with too many hits
        if(hitMap[ilayer].size() > recoConfig.findpattern_max_hit_layers) continue;

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

            if(verbose > 4)
                std::cout << "Layer: " << ilayer << " Hit1: " << hit1.ID << " Hit2: " << hitmin.ID << " Dist: " << sqrt(distmin) << std::endl;
            if(distmin>recoConfig.findpattern_dist_min_cut) continue;

            // now create a TKTrack with the doublet if hits are close
            TTKTrack *trk = new TTKTrack();
            trk->tkhit.push_back(hit1);
            trk->tkhit.push_back(hitmin);
            trk->direction = trk->direction2Hits();
            ROOT::Math::XYZVector centroid = (hit1.point + hitmin.point) * 0.5;
            trk->centroid.SetXYZ(centroid.X(), centroid.Y(), centroid.Z());
            tempTracks.push_back(trk);
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
    std::sort(tempTracks.begin(), tempTracks.end(), [](const TTKTrack *a, const TTKTrack *b)
              { return a->centroid.Z() < b->centroid.Z(); });

    // now loop over tracks and merge segments
    for (auto trk : tempTracks) {
        trk->SortHitsByZ();
    }

    // dump temp tracks
    if(verbose > 4) {
        for (auto trk : tempTracks) {
            trk->Dump();
        }
    }

    double mindZcut = fTcalEvent->geom_detector.fTargetSizeZ*recoConfig.findpattern_mindZ_fudge;

    for (size_t i = 0; i < tempTracks.size(); i++) {
        TTKTrack *track1 = tempTracks[i];
        // for subsequent iterations, consider only tracks with more than 2 hits
        size_t besttrk = -1;
        double mindZ = 1e9;
        for( size_t j = 0; j < tempTracks.size(); j++) {
            if (i==j) continue;
            TTKTrack *track2 = tempTracks[j];
            // check if segments are parallel
            TVector3 normDir1 = track1->direction.Unit();
            TVector3 normDir2 = track2->direction.Unit();
            double dotProduct = normDir1.Dot(normDir2);
            if (std::abs(std::abs(dotProduct) - 1.0) > recoConfig.findpattern_parallel_cut) continue;
            // ensure that segments belong to different planes
            struct TTKTrack::TRACKHIT hit1 = track1->tkhit.back();
            struct TTKTrack::TRACKHIT hit2 = track2->tkhit.front();
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
            if (std::abs(dotProduct - 1.0) > recoConfig.findpattern_parallel_cut) continue;
            if(dz < mindZ) {
                mindZ = dz;
                besttrk = j;
            }
        }
        // if match found, then merge the tracks
        if(besttrk != -1) {
            // add hits from track2 to track1
            for (const auto& h : tempTracks[besttrk]->tkhit) {
                track1->tkhit.push_back(h);
            }
            track1->direction = track1->fitLineThroughHits(track1->centroid);
            // erase best track
            delete tempTracks[besttrk];
            tempTracks.erase(tempTracks.begin() + besttrk);
            if(besttrk < i) i--;
        }
    }

    // get rid of tracks with less than 3 hits
    #if 0
    fTKTracks.erase(std::remove_if(fTKTracks.begin(), fTKTracks.end(),
                                   [](const TTKTrack &track)
                                   {
                                       return track.tkhit.size() < 3;
                                   }),
                    fTKTracks.end());
    #endif

    // always order hits in the track by increasing "z"
    for (auto trk : tempTracks) {
        trk->SortHitsByZ();
    }

    if(verbose > 0)
        std::cout << "Total number of tracks before merging: " << tempTracks.size() << std::endl;
    // dump temp tracks
    if(verbose > 4) {
        for (auto trk : tempTracks) {
            trk->Dump();
        }
    }

    #if 0
    // merge broken tracks
    double cut_SSR_merge = 10; // FIXME: adjust value
    double cut_parallel_merge = 0.01; // FIXME: adjust value
    for (size_t i = 0; i < tempTracks.size(); i++) {
        TTKTrack *track1 = tempTracks[i];
        for( size_t j = 0; j < tempTracks.size(); j++) {
            if (i==j) continue;
            TTKTrack *track2 = tempTracks[j];

            // check if segments are parallel
            TVector3 normDir1 = track1->direction.Unit();
            TVector3 normDir2 = track2->direction.Unit();
            double dotProduct = normDir1.Dot(normDir2);
            if (std::abs(std::abs(dotProduct) - 1.0) > cut_parallel_merge) continue;

            TTKTrack *tempTrack = new TTKTrack(*track1);
            tempTrack->MergeTracks(*track2);
            // fit line through hits
            tempTrack->direction = tempTrack->fitLineThroughHits(tempTrack->centroid);
            if(verbose > 0) {
                std::cout << "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM" << std::endl;
                std::cout << "Trying to merge tracks: " << i << " " << j << "SSR: " << tempTrack->SSR << std::endl;
                track1->Dump();
                track2->Dump();
                tempTrack->Dump();
            }
            if(tempTrack->SSR < cut_SSR_merge) {
                tempTracks[i] = tempTrack;
                tempTracks.erase(tempTracks.begin() + j);
                delete track2;
                if(j < i) i--;
                std::cout << "---- MERGED";
            } else {
                delete tempTrack;
            }
            std::cout << std::endl;
        }
    }
    #endif

    // order tracks by increasing "z"
    std::sort(tempTracks.begin(), tempTracks.end(), [](const TTKTrack *a, const TTKTrack *b)
              { return a->centroid.Z() < b->centroid.Z(); });

    if(verbose > 0)
        std::cout << "Total number of tracks before fitting: " << tempTracks.size() << std::endl;
    // always order hits in the track by increasing "z"
    for (auto trk : tempTracks) {
        trk->SortHitsByZ();
    }

    // now fit tracks with GENFIT2
    for (auto trk : tempTracks) {
        trk->GenFitTrackFit(fTcalEvent->geom_detector.fScintillatorVoxelSize*recoConfig.psvoxel_fudge_factor);
        // print fit result
        if(verbose > 0) {
            std::cout << "Fitted track - ";
            TVector3 momf = trk->fitTrack->getFittedState().getMom();
            std::cout << "Momentum: " << momf.Mag() << std::endl;
        }
        if(verbose > 3) {
            trk->fitTrack->Print();
        }
    }

    if(verbose > 0)
        std::cout << "Total number of tracks after fitting: " << tempTracks.size() << std::endl;

    // remove tracks with bad fitted pval (pval < 0.01) or low momentum
    for (size_t i = 0; i < tempTracks.size(); i++) {
        TTKTrack *track = tempTracks[i];
        if(!track->fitTrack->getFitStatus()->isFitConverged() ||
           (track->fitTrack->getFitStatus()->getChi2() > 0 && track->fitTrack->getFitStatus()->getPVal() < recoConfig.genfit_min_pVal) ||
           track->fitTrack->getFittedState().getMom().Mag() < recoConfig.genfit_min_pMom) {
            if(verbose > 0) {
                std::cout << "Removing track with bad fit: " << i << std::endl;
                track->Dump(verbose);
            }
            delete track;
            tempTracks.erase(tempTracks.begin() + i);
            i--;
        }
    }

    if(verbose > 0)
        std::cout << "Total number of tracks after quality cuts: " << tempTracks.size() << std::endl;
    // dump temp tracks
    if(verbose > 3) {
        for (auto trk : tempTracks) {
            trk->Dump();
        }
    }


    // now merge broken tracks
    double cut_merge_mindZ = fTcalEvent->geom_detector.fTargetSizeZ*recoConfig.findpattern_mindZ_fudge;
    for (size_t i = 0; i < tempTracks.size(); i++) {
        for( size_t j = 0; j < tempTracks.size(); j++) {
            if (i==j) continue;
            TTKTrack *track1 = tempTracks[i];
            TTKTrack *track2 = tempTracks[j];

            struct TTKTrack::TRACKHIT hit1 = track1->tkhit.front();
            struct TTKTrack::TRACKHIT hit2 = track2->tkhit.front();
//            if(hit1.ID == 110325340785L && hit2.ID == 110526731650L) raise(SIGTRAP);
            double dz = std::abs(hit1.point.z()-hit2.point.z());
            if(dz < cut_merge_mindZ) continue;

            // check if segments are parallel
            TVector3 normDir1 = track1->direction.Unit();
            TVector3 normDir2 = track2->direction.Unit();
            double dotProduct = normDir1.Dot(normDir2);
            if (std::abs(std::abs(dotProduct) - 1.0) > recoConfig.findpattern_parallel_cut_merge) continue;

            TTKTrack *tempTrack = new TTKTrack(*track1);
            tempTrack->MergeTracks(*track2);
            tempTrack->GenFitTrackFit(fTcalEvent->geom_detector.fScintillatorVoxelSize*recoConfig.psvoxel_fudge_factor);
            if(verbose > 0) {
                std::cout << "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM"<<std::endl;
                std::cout << "Trying to merge tracks: " << i << " " << j << std::endl;
                track1->Dump();
                track2->Dump();
                tempTrack->Dump();
            }
            if(tempTrack->fitTrack->getFitStatus()->isFitConverged() &&
               (tempTrack->fitTrack->getFitStatus()->getChi2() > 0 && tempTrack->fitTrack->getFitStatus()->getPVal() > recoConfig.genfit_min_pVal) &&
               tempTrack->fitTrack->getFittedState().getMom().Mag() > recoConfig.genfit_min_pMom) {
                if(verbose > 0)
                    std::cout << "---- MERGED" << std::endl;
                tempTrack->direction = tempTrack->fitLineThroughHits(tempTrack->centroid);
                tempTracks[i] = tempTrack;
                tempTracks.erase(tempTracks.begin() + j);
//                delete track1;  // CRASHES
                delete track2;
                if(j < i) i--;
            } else {
                delete tempTrack;
                if(verbose > 0)
                    std::cout << "FAILED" << std::endl;
            }
        }
    }
                

    // now copy temporary tracks to fTKTracks
    fTKTracks.clear();
    int trackID = 0;
    for (auto trk : tempTracks) {
        trk->trackID = trackID++;
        fTKTracks.push_back(*trk);
        delete trk;
    }
}

void TPORecoEvent::ExtendTracks() {

    if(verbose>0) std::cout << "Start extendingtracks..." << std::endl;

    int nmodules = fTcalEvent->geom_detector.NRep;
    int nzlayer = fTcalEvent->geom_detector.fSandwichLength / fTcalEvent->geom_detector.fScintillatorVoxelSize;

    double closest_voxel_cut = recoConfig.extendtracks_closest_voxel_cut*recoConfig.extendtracks_closest_voxel_cut;

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

    // loop over all tracks
    for (auto &trk : fTKTracks) {
        int min, max;
        trk.GetMinMaxModule(min, max);
        int module = min;
        int nmodules = fTcalEvent->geom_detector.NRep;
        bool stop_extend = false;
        while(module < nmodules && !stop_extend) {
            // loop over all layers in module
            for (int iz = 0; iz < nzlayer; iz++) {
                // get the closest voxel to the track
                struct PSVOXEL3D closestVoxel;
                int failed(0);
                double z = fTcalEvent->getZofLayer(module, iz);
                if(verbose > 3) {
                    std::cout << "Module: " << module << " Layer: " << iz << " Z: " << z << " ";
                    std::cout << "Track: fittrak: " << trk.fitTrack << std::endl;
                }
                TVector3 pos = trk.extrapolateTracktoZ(z, failed);
                if(failed) continue;
                if(verbose>3)
                    std::cout << "Successfully extrapolated to Z: " << z << " Position: " << pos.X() << " " << pos.Y() << " " << pos.Z() << std::endl;
                // cut on the fiducial volume
                if(std::abs(pos.X()) > fTcalEvent->geom_detector.fScintillatorSizeX/2.0) stop_extend = true;
                if(std::abs(pos.Y()) > fTcalEvent->geom_detector.fScintillatorSizeY/2.0) stop_extend = true;
                if(stop_extend) continue;
                for (const auto &voxel : PSvoxelmapModule[module]) {
                    if(iz != fTcalEvent->getChannelnzfromID(voxel.ID)) 
                        continue;
                    // get position of voxel
                    ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(voxel.ID);
                    double dist2 = (pos.X()-position.X())*(pos.X()-position.X()) + (pos.Y()-position.Y())*(pos.Y()-position.Y());
                    if(dist2 > recoConfig.extendtracks_dist2_perp_voxel_cut) 
                        continue;
                    if(verbose > 1)
                        std::cout << "Trying to add closest voxel: " << voxel.ID << " Dist: " << sqrt(dist2) << std::endl;
                    float ehit = -1; // FIXME: get energy from voxel
                    struct TTKTrack::TRACKHIT hit = {voxel.ID, (int)TcalEvent::getChannelTypefromID(voxel.ID), position, ehit};

                    // check if track already has a voxel in this layer
                    if(verbose>1)
                        std::cout << "Checking if another at the same layer is already in track" << std::endl;
                    bool found = false;
                    double dist2existing = 1e9;
                    long IDexisting = -1;
                    for (const auto &hit : trk.tkhit) {
                        // check if this hit is in the same module
                        if(fTcalEvent->getChannelModulefromID(hit.ID) != module) continue;
                        if(fTcalEvent->getChannelnzfromID(hit.ID) == iz) {
                            // get position of voxel
                            ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(hit.ID);
                            dist2existing = (pos.X()-position.X())*(pos.X()-position.X()) + (pos.Y()-position.Y())*(pos.Y()-position.Y());
                            IDexisting = hit.ID;
                            found = true;
                            break;
                        }
                    }
                    if(found) {
                        if(verbose > 1)
                            std::cout << "Found existing hit in the same layer: " << IDexisting << " Dist: " << sqrt(dist2existing) << std::endl;
                        if(dist2 < dist2existing) {
                            if(verbose > 1)
                                std::cout << "Replacing hit in the same layer: " << IDexisting << " Dist: " << sqrt(dist2existing) << 
                                    " with new one " << hit.ID << " Dist: " << sqrt(dist2) << std::endl;
                            trk.tkhit.erase(std::remove_if(trk.tkhit.begin(), trk.tkhit.end(),
                                                           [IDexisting](const TTKTrack::TRACKHIT &hit2)
                                                           {
                                                               return IDexisting == hit2.ID;
                                                           }),
                                           trk.tkhit.end());
                            found = false; // force adding the new hit
                        }
                    }
                    if(found) continue;
                    if(verbose > 1)
                        std::cout << "Adding hit " << hit.ID << " to track" << std::endl;
                    trk.tkhit.push_back(hit);
                    // now mark voxel as used by track
                    // find voxel in PSvoxelmap
                    auto it = PSvoxelmap.find(voxel.ID);
                    if(it != PSvoxelmap.end()) {
                        it->second.member_of_TKtrack = true;
                    }
                }
            }
        module++;
        }
        // sort hits by z
        trk.SortHitsByZ();
        if(verbose > 1){
            std::cout << "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB" << std::endl;
            std::cout << "Track " << trk.trackID << " - Before refit" << std::endl;
            trk.Dump();
        }
        delete trk.fitTrack;
        trk.GenFitTrackFit(fTcalEvent->geom_detector.fScintillatorVoxelSize*recoConfig.psvoxel_fudge_factor);
        if(verbose > 1){
            std::cout << "Track " << trk.trackID << "After refit" << std::endl;
            trk.Dump();
            // trk.fitTrack->Print();
        }
        trk.UpdateFittedPosition();
    }

    // dump temp tracks
    if(verbose > 1) {
        std::cout << "Extended tracks: " << fTKTracks.size() << std::endl;
        for (auto &trk : fTKTracks) {
            trk.Dump();
        }
    }

}

void TPORecoEvent::DumpReconstructedTracks() {
    // dump tracks
    std::cout << "Total number of tracks: " << fTKTracks.size() << std::endl;
    for (auto &trk : fTKTracks) {
        trk.Dump(verbose);
    }
    // dump vertices
    std::cout << "Total number of vertices: " << fTKVertices.size() << std::endl;
    for (auto &vertex : fTKVertices) {
        std::cout << "Vertex: " << vertex.position.X() << " " << vertex.position.Y() << " " << vertex.position.Z() << " ";
//        std::cout << "Covariance: " << vertex.covariance(0, 0) << " " << vertex.covariance(1, 1) << " " << vertex.covariance(2, 2) << std::endl;
        std::cout << "Errors on position: " << sqrt(vertex.covariance(0, 0)) << " " << sqrt(vertex.covariance(1, 1)) << " " << sqrt(vertex.covariance(2, 2)) << " ";
        std::cout << "Ndf: " << vertex.ndof << ", Chi2: " << vertex.chi2 << ", Ntracks: " << vertex.ntracks << std::endl;
    }
}

void TPORecoEvent::FindTrackVertices() {
    genfit::GFRaveVertexFactory vertexFactory(5);
    vertexFactory.setMethod("kalman-smoothing:1");

    fTKVertices.clear();
    // reset all tracks vertexID
    for (auto &trk : fTKTracks) {
        trk.vertexID = -1;
    }

    std::vector<TTKTrack*> tempTracks;
    for (auto &trk : fTKTracks) {
        if(trk.tkhit.size() < 3) continue;
        // reject tracks with bad chi2
        if(trk.fitTrack->getFitStatus()->getChi2() > 0 && trk.fitTrack->getFitStatus()->getPVal() < recoConfig.genfit_min_pVal) continue;
        tempTracks.push_back(&trk);
    }
    std::cout << "+V+V+ Total number of tracks: " << tempTracks.size() << " retained out of " << fTKTracks.size() << std::endl;

    // order tracks by increasing "z" of the first hit of the track
    std::sort(tempTracks.begin(), tempTracks.end(), [](const TTKTrack *a, const TTKTrack *b)
                { return a->tkhit.front().point.z() < b->tkhit.front().point.z(); });

    // if there are less than 2 tracks, then no vertexing
    if(tempTracks.size() < 3) return;

    size_t ntracks = tempTracks.size();
#if 1
    if(ntracks > recoConfig.findvtx_cut_max_trk) {
        ntracks = recoConfig.findvtx_cut_max_trk;
        std::cerr << "+V+V+ Number of tracks considered for vertexing cut to: " << ntracks << std::endl;
    }
#endif

    struct TUPLET {
        std::set<int> tracks;
        TVector3 position;          // in cm
        double chi2;
        int ndf;
        double chi2ndf;            // chi2/ndf of the vertex
        // FIXME: add covariance matrix
        TMatrixDSym covariance;

        TUPLET()
           : covariance(3), chi2(0.0) {} // Initialize covariance with size 3
    };
    std::vector <TUPLET> tuplets;

    // loop over all possible combinations of 3 tracks
    for (int i = 0; i < ntracks-2; i++) {
        for (int j = i+1; j < ntracks-1; j++) {
//            if(tempTracks[i]->vertexID != -1) continue;
//            if(tempTracks[j]->vertexID != -1) continue;
            // first try to fit a vertex with two tracks
            // prepare tracks for vertexing    
            std::vector<genfit::Track*> tracks;
            std::vector<genfit::GFRaveVertex*> vertices;
            tempTracks[i]->fitTrack->checkConsistency();
            tempTracks[j]->fitTrack->checkConsistency();
            tracks.push_back(tempTracks[i]->fitTrack);
            tracks.push_back(tempTracks[j]->fitTrack);
            // vertexing
            if(verbose > 0) {
                std::cout << "****************************************************************************" << std::endl;
                std::cout << "Calling findVertices with two tracks: " << tempTracks[i]->trackID << " " << tempTracks[j]->trackID << std::endl;  
            }
            vertexFactory.findVertices(&vertices, tracks);
            if(verbose > 0) {
                std::cout << "number of vertices: " << vertices.size() << std::endl;
                // print reconstructed vertices
                for (auto &vertex : vertices) {
                    std::cout << "GFRaveVertex of two tracks " << tempTracks[i]->trackID << " " << tempTracks[j]->trackID << std::endl;
                    std::cout << "Position: "; vertex->getPos().Print();
//                    std::cout << "Covariance: "; vertex->getCov().Print();
                    std::cout << "Ndf: " << vertex->getNdf() << ", Chi2: " << vertex->getChi2() << ", Id: " << vertex->getId() << "\n";
                }
            }
            // apply quality cuts
            bool good_vertex = true;
            if(vertices.size() == 0) 
                good_vertex = false;
            else {
                // negative chi2 means that vertexing failed due to numerical reasons
                if(vertices[0]->getChi2() < 0) good_vertex = false;
                if(vertices[0]->getChi2()/vertices[0]->getNdf() > recoConfig.findvtx_chi2ndf_cut) good_vertex = false;
                // if position is outside of the detector then break
                if(std::abs(vertices[0]->getPos().X())*10.0 > fTcalEvent->geom_detector.fScintillatorSizeX/2.0) good_vertex = false;
                if(std::abs(vertices[0]->getPos().Y())*10.0 > fTcalEvent->geom_detector.fScintillatorSizeY/2.0) good_vertex = false;
                // check z position of the vertex witihn the detector
                if(vertices[0]->getPos().Z()*10.0 > fTcalEvent->geom_detector.rearCalLocZ) good_vertex = false;
                if(vertices[0]->getPos().Z()*10.0 < -fTcalEvent->geom_detector.fTotalLength/2.0) good_vertex = false;
            }
            if(!good_vertex) {
                // delete vertices
                for (auto &vertex : vertices) {
                    delete vertex;
                }
                continue;
            }
    
            for (int k = j+1; k < ntracks; k++) {
//                if(tempTracks[i]->vertexID != -1) continue;
//                if(tempTracks[j]->vertexID != -1) continue;
//                if(tempTracks[k]->vertexID != -1) continue;
                // prepare tracks for vertexing    
                std::vector<genfit::Track*> tracks;
                std::vector<genfit::GFRaveVertex*> vertices;
                tracks.push_back(tempTracks[i]->fitTrack);
                tracks.push_back(tempTracks[j]->fitTrack);
                tracks.push_back(tempTracks[k]->fitTrack);
                tempTracks[i]->fitTrack->checkConsistency();
                tempTracks[j]->fitTrack->checkConsistency();
                tempTracks[k]->fitTrack->checkConsistency();
                // vertexing
                if(verbose > 4) {
                    std::cout << "+++++++++++++++++++++++++++++++ 33333333333 +++++++++++++++++++++++++++++++" << std::endl;
                    std::cout << "Calling findVertices with three tracks: " << tempTracks[i]->trackID << " " << 
                        tempTracks[j]->trackID << " " << tempTracks[k]->trackID << std::endl;
                }
                vertexFactory.findVertices(&vertices, tracks);
                if(verbose > 4) {
                    std::cout << "number of vertices: " << vertices.size() << std::endl;
                    // print reconstructed vertices
                    for (auto &vertex : vertices) {
                        std::cout << "GFRaveVertex of three tracks " << tempTracks[i]->trackID << " " << tempTracks[j]->trackID << " " << tempTracks[k]->trackID << std::endl;
                        std::cout << "Position: "; vertex->getPos().Print();
//                        std::cout << "Covariance: "; vertex->getCov().Print();
                        std::cout << "Ndf: " << vertex->getNdf() << ", Chi2: " << vertex->getChi2() << ", Id: " << vertex->getId() << "\n";
                    }
                }
                good_vertex = true;
                if(vertices.size() == 0) 
                    good_vertex = false;
                else {
                // negative chi2 means that vertexing failed due to numerical reasons
                    if(vertices[0]->getChi2() < 0) good_vertex = false;
                    if(vertices[0]->getChi2()/vertices[0]->getNdf() > recoConfig.findvtx_chi2ndf_cut) good_vertex = false;
                    // if position is outside of the detector then break
                    if(std::abs(vertices[0]->getPos().X())*10.0 > fTcalEvent->geom_detector.fScintillatorSizeX/2.0) good_vertex = false;
                    if(std::abs(vertices[0]->getPos().Y())*10.0 > fTcalEvent->geom_detector.fScintillatorSizeY/2.0) good_vertex = false;
                    // check z position of the vertex witihn the detector
                    if(vertices[0]->getPos().Z()*10.0 > fTcalEvent->geom_detector.rearCalLocZ) good_vertex = false;
                    if(vertices[0]->getPos().Z()*10.0 < -fTcalEvent->geom_detector.fTotalLength/2.0) good_vertex = false;

                    // check that vertex is close to the three tracks
                    ROOT::Math::XYZVector vtxpos(vertices[0]->getPos().X()*10.0, vertices[0]->getPos().Y()*10.0, vertices[0]->getPos().Z()*10.0);
                    ROOT::Math::XYZVector dist1 = vtxpos - tempTracks[i]->tkhit[0].point;
                    ROOT::Math::XYZVector dist2 = vtxpos - tempTracks[j]->tkhit[0].point;
                    ROOT::Math::XYZVector dist3 = vtxpos - tempTracks[k]->tkhit[0].point;
                    if(verbose > 3)
                        std::cout << "Distances: " << dist1.R() << " " << dist2.R() << " " << dist3.R() << std::endl;
                    if(dist1.R() > recoConfig.findvtx_trk_dist_cut || 
                        dist2.R() > recoConfig.findvtx_trk_dist_cut || 
                        dist3.R() > recoConfig.findvtx_trk_dist_cut) good_vertex = false;
                }
                if(good_vertex) {
                    // add vertex to the list of triplets
                    TUPLET aTuplet;
                    aTuplet.tracks = {i, j, k};
                    aTuplet.position = vertices[0]->getPos();
                    aTuplet.chi2 = vertices[0]->getChi2();
                    aTuplet.ndf = vertices[0]->getNdf();
                    aTuplet.chi2ndf = vertices[0]->getChi2()/vertices[0]->getNdf();
                    aTuplet.covariance = vertices[0]->getCov();
                    tuplets.push_back(aTuplet);
                    // mark tracks as used
//                    tempTracks[i]->vertexID = vertices[0]->getId();
//                    tempTracks[j]->vertexID = vertices[0]->getId();
//                    tempTracks[k]->vertexID = vertices[0]->getId();
                }
                // delete vertices
                for (auto &vertex : vertices) {
                    delete vertex;
                }
            }
        }
    }

    // dump all elements of set in tuplets
    if(verbose > 0){
        for (auto &tuplet : tuplets)
        {
            std::cout << "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" << std::endl;
            std::cout << "Tuplet: ";
            for (auto &track : tuplet.tracks)
            {
                std::cout << track << " ";
            }
            std::cout << "Chi2/ndf: " << tuplet.chi2 << std::endl;
            std::cout << "Position: " << tuplet.position.X() << " " << tuplet.position.Y() << " " << tuplet.position.Z() << std::endl;
            std::cout << std::endl;
            std::cout << "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" << std::endl;
        }
    }

    // now merge vertices that have at least two tracks in common
    for(int iter = 0; iter < 2; iter++) {

        for(int i = 0; i < tuplets.size(); i++) {
            for (int j = 0; j < tuplets.size(); j++) {
                if(i == j) continue;
                // check if two tuplets have at least two tracks in common
                // compute intersection of two sets
                std::set<int> intersection;
                std::set_intersection(tuplets[i].tracks.begin(), tuplets[i].tracks.end(),
                                    tuplets[j].tracks.begin(), tuplets[j].tracks.end(),
                                    std::inserter(intersection, intersection.begin()));
                // compute distance between vertices
                TVector3 dist = tuplets[i].position - tuplets[j].position;
                double distance = dist.Mag();
                // FIXME: replace distance cut by cut on distance/uncertainty
                // uncertainty on diffrence is given by sum of covariance matrices
                if(intersection.size() > 1 && distance < recoConfig.findvtx_merge_dist_cut) {
                    if(verbose > 0)
                        std::cout << "Merge i: " << i << " j: " << j << " - Distance between vertices: " << distance << std::endl;
                    // merge two tuplets
                    std::set<int> merged;
                    std::set_union(tuplets[i].tracks.begin(), tuplets[i].tracks.end(),
                                tuplets[j].tracks.begin(), tuplets[j].tracks.end(),
                                std::inserter(merged, merged.begin()));
                    // replace old tuplet with new one
                    TUPLET aTuplet;
                    aTuplet.tracks = merged;
                    aTuplet.position = (tuplets[i].position + tuplets[j].position) * 0.5;
                    aTuplet.chi2 = -1;
                    aTuplet.ndf = 0;
                    aTuplet.chi2ndf = 0.0;
                    aTuplet.covariance = tuplets[i].covariance + tuplets[j].covariance;
                    tuplets[i] = aTuplet;
                    tuplets.erase(tuplets.begin() + j);
                    j--;
                }
            }
        }

        // now fit vertices with tracks from tuplets
        for (int i = 0; i < tuplets.size(); i++) {
            auto &tuplet = tuplets[i];
            // prepare tracks for vertexing    
            std::vector<genfit::Track*> tracks;
            std::vector<genfit::GFRaveVertex*> vertices;
            for(auto &track : tuplet.tracks) 
                tracks.push_back(tempTracks[track]->fitTrack);
            // vertexing

            if(verbose > 0) {
                std::cout << "+++++++++++++++++++++++++++++++ MMMMMMMMMMMM +++++++++++++++++++++++++++++++" << std::endl;
                std::cout << "Calling findVertices with " << tracks.size() << " tracks: ";
                for(auto &track : tuplet.tracks) 
                    std::cout << track << " ";
                std::cout << std::endl; 
            }

            vertexFactory.findVertices(&vertices, tracks);

            if(verbose > 0) {
                std::cout << "number of vertices: " << vertices.size() << std::endl;
                // print reconstructed vertices
                for (auto &vertex : vertices) {
                    std::cout << "FINAL GFRaveVertex of several tracks\n";
                    std::cout << "Position: "; vertex->getPos().Print();
                    std::cout << "Covariance: "; vertex->getCov().Print();
                    std::cout << "Ndf: " << vertex->getNdf() << ", Chi2: " << vertex->getChi2() << ", Id: " << vertex->getId() << "\n";
                    std::cout << "Number of tracks: " << vertex->getNTracks() << "\n";
                }
            }

            bool good_vertex = true;
            if(vertices.size() == 0) 
                good_vertex = false;
            else {
                // negative chi2 means that vertexing failed due to numerical reasons
                if(vertices[0]->getChi2() < 0) good_vertex = false;
                if(vertices[0]->getChi2()/vertices[0]->getNdf() > recoConfig.findvtx_chi2ndf_cut) good_vertex = false;
                // if position is outside of the detector then break
                if(std::abs(vertices[0]->getPos().X())*10.0 > fTcalEvent->geom_detector.fScintillatorSizeX/2.0) good_vertex = false;
                if(std::abs(vertices[0]->getPos().Y())*10.0 > fTcalEvent->geom_detector.fScintillatorSizeY/2.0
                ) good_vertex = false;
                // check z position of the vertex witihn the detector
                if(vertices[0]->getPos().Z()*10.0 > fTcalEvent->geom_detector.rearCalLocZ) good_vertex = false;
                if(vertices[0]->getPos().Z()*10.0 < -fTcalEvent->geom_detector.fTotalLength/2.0) good_vertex = false;
            }

            // if chi2 is bad then skip
            if(!good_vertex) {
                // erase tuplet that failed
                if(verbose > 0)
                    std::cout << "Erasing tuplet that failed to refit: " << i << std::endl;
                tuplets.erase(tuplets.begin() + i);
                i--;
                continue;
            }

            // update position of the vertex
            tuplet.position = vertices[0]->getPos();
            // update covariance matrix
            tuplet.covariance = vertices[0]->getCov();
            // update chi2 and ndf
            tuplet.chi2 = vertices[0]->getChi2();
            tuplet.ndf = vertices[0]->getNdf();
            tuplet.chi2ndf = vertices[0]->getChi2()/vertices[0]->getNdf();
            // delete vertices
            for (auto &vertex : vertices) {
                delete vertex;
            }
        }
    }

    if(verbose > 0){
        for (auto &tuplet : tuplets)
        {
            std::cout << "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" << std::endl;
            std::cout << "Tuplet: ";
            for (auto &track : tuplet.tracks)
            {
                std::cout << track << " ";
            }
            std::cout << "Chi2/ndf: " << tuplet.chi2 << std::endl;
            std::cout << "Position: " << tuplet.position.X() << " " << tuplet.position.Y() << " " << tuplet.position.Z() << std::endl;
            std::cout << std::endl;
            std::cout << "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" << std::endl;
        }
    }

    // now fill the list of vertices
    int vertexID = 0;
    for (auto &tuplet : tuplets) {
        struct TTKVertex aVertex;
        aVertex.vertexID = vertexID++;
        aVertex.position = tuplet.position;
        aVertex.position *= 10.0; // convert to mm
        aVertex.covariance = tuplet.covariance;
        aVertex.covariance *= 100.0; // convert to mm
        aVertex.chi2 = tuplet.chi2;
        aVertex.ndof = tuplet.ndf;
        // add tracks to the vertex
        for(auto &trackID : tuplet.tracks) 
            aVertex.trackIDs.push_back(tempTracks[trackID]->trackID);
        fTKVertices.push_back(aVertex);
    }

    // now loop over tracks and assign vertexID
    for (auto &trk : fTKTracks) {
        for (auto &tk : tempTracks) {
            if(trk.trackID == tk->trackID) {
                trk.vertexID = tk->vertexID;
                break;
            }
        }
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

    DBScan dbscan;

    std::map<int, class TPSCluster> PSClustersMap;

    std::vector<DBScan::Point> points;

    // XZ or YZ view
    for (auto it : (view==0) ? PShitmapX : PShitmapY)
    {
        long ID = it.first;
        double ehit = it.second.Edeposited;
        if(ehit < recoConfig.clusters_threshold_2dhit) continue;
        double fix, fiy, fiz;
        pshit2d_position(ID, fix, fiy, fiz);
        DBScan::Point p = {ID, ehit, (view==0) ? fix : fiy, fiz};
        points.push_back(p);
    }
    dbscan.scan(points, recoConfig.clusters_eps, recoConfig.clusters_minPts);
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
        if (c.second.rawenergy < recoConfig.clusters_threshold_cluster)
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

///// OBSOLETE - use Reconstruct3DPS_2

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
                struct PSVOXEL3D v = {ID, ehit, true, false};
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

    int nx = fTcalEvent->geom_detector.fScintillatorSizeX / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int ny = fTcalEvent->geom_detector.fScintillatorSizeY / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int nzlayer = fTcalEvent->geom_detector.fSandwichLength / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int nrep = fTcalEvent->geom_detector.NRep;
    int nztot =  nrep * nzlayer;

    std::random_device rd;  // Seed for the random number generator
    std::mt19937 gen(rd());  // Mersenne Twister random number generator
    std::uniform_int_distribution<> rnd_nx(0, nx-1);
    std::uniform_int_distribution<> rnd_ny(0, ny-1);

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
        if (iter > recoConfig.PS3D_nvox_max_after_iteration && sum_nvox > nvox_per_layer_max)
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
        if (adjusted == 0 || total_score < recoConfig.PS3D_total_score_min_break)
            break;
    }
    if (verbose > 1)
        std::cout << "Module " << imodule << " - Score = " << total_score << std::endl;
}

void TPORecoEvent::Reconstruct3DPS_2(int maxIter) {

    maxIter = 150; // 300;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine rng(seed);

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
        if(ehit < recoConfig.PS3D_ehit_threshold) continue;
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
        if(ehit < recoConfig.PS3D_ehit_threshold) continue;
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
            if (ehit < recoConfig.PS3D_ehit_threshold)
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

    if(multiThread) {
        std::vector<std::thread> threads;
        for (int imodule = 0; imodule < nrep; imodule++)
        {
    //       reconstruct3DPS_module(maxIter, imodule, V, XZ, YZ, XY, nvox_per_layer, nvox_per_layer_max);
            threads.emplace_back(&TPORecoEvent::reconstruct3DPS_module, this, maxIter, imodule, 
                std::ref(V), std::ref(XZ), std::ref(YZ), std::ref(XY),
                std::ref(nvox_per_layer), recoConfig.PS3D_nvox_per_layer_max);
        }

        for (auto& th : threads)
        {
            th.join();
        }
    } else {
        for (int imodule = 0; imodule < nrep; imodule++)
        {
            reconstruct3DPS_module(maxIter, imodule, V, XZ, YZ, XY, nvox_per_layer, recoConfig.PS3D_nvox_per_layer_max);
        }
    }

    PSvoxelmap.clear();
    for (int z = 0; z < nztot; ++z) {
        //
        if(nvox_per_layer[z] > recoConfig.PS3D_nvox_per_layer_max) continue;
        //
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                float ehit = V[x][y][z].value;

                // cut on minimum energy
                if(ehit<recoConfig.PS3D_evox_threshold)
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
    if(verbose > 0) {
        for(int imodule = 0; imodule < nrep; imodule++){
            double frac = ntotl[imodule] > 0 ? nfake[imodule]*100.0/float(ntotl[imodule]) : -1;
            std::cout << " Module " << imodule << " " << ntotl[imodule] << " hits " << nfake[imodule] << " fakes ";
            std::cout << frac << "%" << std::endl;
        }
        std::cout << " STATS: " << ntot << " hits " << fakes << " ghosts." << std::endl;
        std::cout << " Avg energy: " << ave_e_real/float(nreal) << " ghosts: " << ave_e_ghost/float(fakes);
        std::cout << std::endl;
    }
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

    if(verbose>0) std::cout << "Start PS voxel particle filter..." << std::endl;

    int nmodules = fTcalEvent->geom_detector.NRep;
    int nzlayer = fTcalEvent->geom_detector.fSandwichLength / fTcalEvent->geom_detector.fScintillatorVoxelSize;

    double closest_voxel_cut2 = recoConfig.PSFilter_closest_voxel_cut*recoConfig.PSFilter_closest_voxel_cut;

    std::map<int, std::vector<struct PSVOXEL3D>> PSvoxelmapModule;
    // organize voxels by module number and skip voxels that belong to a TK track
    for (const auto &it : PSvoxelmap) {
        long ID = it.first;

        // check if this voxel already belong to a TK track
        if(it.second.member_of_TKtrack) {
            if(verbose>3)
                std::cout << "Voxel " << ID << " already belongs to a track." << std::endl;
            continue;
        }

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
        if(track_seeds.size() > recoConfig.PSFilter_max_number_track_seeds)
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
                double transv_dist2 = 1e9;
                struct PSVOXEL3D *closest_vox = nullptr;
                for (auto &prev_hit : module->second)
                {
                    long prev_iz = (prev_hit.ID / 1000000) % 1000;
                    if (prev_iz != prev_layer)
                        continue;
                    ROOT::Math::XYZVector prev_position = fTcalEvent->getChannelXYZfromID(prev_hit.ID);
                    ROOT::Math::XYZVector diff = position - prev_position;
                    double d2 = diff.x() * diff.x() + diff.y() * diff.y();
                    if (d2 < transv_dist2)
                    {
                        transv_dist2 = d2;
                        closest_vox = &prev_hit;
                    }
                }
//                std::cout << track->direction.x() << " " << track->direction.x()*track->direction.x() << std::endl;
//                std::cout << track->direction.y() << " " << track->direction.y()*track->direction.y() << std::endl;
                double rt = track->direction.x()*track->direction.x() + track->direction.y()*track->direction.y();
//                std::cout << rt << std::endl;
                double tanangle = sqrt(rt) / (std::max(track->direction.z(),1e-3));
//                std::cout << tanangle << std::endl;
                double cut2 = closest_voxel_cut2*std::max(1.0,tanangle*tanangle);
//                std::cout << cut << std::endl;
#if 0
                if(rt>0) {
                    std::cout << " track angle found - ";
                    track->Dump();
                }
#endif
                if (transv_dist2 < cut2)
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
                if (std::abs(std::abs(dotProduct) - 1.0) > recoConfig.PSFilter_parallel_cut)
                    continue;
                // ensure that segments belong to different planes
                struct TPSTrack::TRACKHIT hit1 = track1.tkhit.back();
                struct TPSTrack::TRACKHIT hit2 = track2.tkhit.front();
                double dz = std::abs(hit1.point.z() - hit2.point.z());
                if (dz < recoConfig.PSFilter_mindZcut)
                    continue;
                // ensure that segments are parallel to main line joining hits
                ROOT::Math::XYZVector normDir = (hit2.point - hit1.point).Unit();
                dotProduct = std::min(std::abs(normDir.Dot(normDir1)), std::abs(normDir.Dot(normDir2)));
                if (std::abs(dotProduct - 1.0) > recoConfig.PSFilter_parallel_cut)
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
        if (verbose > 3)
            trk.Dump();
    }
}

void TPORecoEvent::ReconstructRearCals() {

    // number of modules; loop over all modules and collect deposited energy
    int nmodules = fTcalEvent->geom_detector.NRep;
    faserCals.clear();
    for(int im = 0; im < nmodules; im++) {
        double eDeposit = 0;
        for (const auto it : fPORecs) {
            int ntracks = it->fGEANTTrackIDs.size();
            for (size_t i = 0; i< ntracks; i++) {
                DigitizedTrack *dt = it->DTs[i];
                size_t nhits = dt->fhitIDs.size();
                for (size_t j = 0; j < nhits; j++) {
                    long ID = dt->fhitIDs[j];
                    long hittype = TcalEvent::getChannelTypefromID(ID);
                    if(hittype != 0) continue; // only scintillator hits
                    long module = fTcalEvent->getChannelModulefromID(ID);
                    if(module == im) {
                        eDeposit += dt->fEnergyDeposits[j];
                    }
                }
            }
        }
        struct FASERCAL faserCal;
        faserCal.ModuleID = im;
        faserCal.EDeposit = eDeposit;
        faserCals.push_back(faserCal);
    }

    double eDeposit = 0;
    double eDepositHCal = 0;
    rearCals.rearCalModule.clear();
    rearCals.rearHCalModule.clear();
    for (const auto &it : fTcalEvent->rearCalDeposit) {
        eDeposit += it.energyDeposit;
        struct TcalEvent::REARCALDEPOSIT im = {it.moduleID, it.energyDeposit};
        rearCals.rearCalModule.push_back(im);
    }
    for (const auto &it : fTcalEvent->rearHCalDeposit) {
        eDepositHCal += it.energyDeposit;
        struct TcalEvent::REARCALDEPOSIT im = {it.moduleID, it.energyDeposit};
        rearCals.rearHCalModule.push_back(im);
    }
    rearCals.rearCalDeposit = eDeposit/1e3;  // convert to GeV
    rearCals.rearHCalDeposit = eDepositHCal/1e3;  // convert to GeV
    rearCals.rearMuCalDeposit = fTcalEvent->rearMuCalDeposit; // in MeV
    if(verbose>0){
        for (const auto &it : faserCals) {
            std::cout << "Module " << it.ModuleID << " - EDeposit " << it.EDeposit << std::endl;
        }
        std::cout << "Rear CalDeposit " << rearCals.rearCalDeposit << std::endl;
        std::cout << "Rear HCalDeposit " << rearCals.rearHCalDeposit << " - ";
        for(auto &it : rearCals.rearHCalModule) {
            std::cout << it.energyDeposit << " ";
        }
        std::cout << std::endl;   
        std::cout << "Rear MuCalDeposit " << rearCals.rearMuCalDeposit << std::endl;
    }
}
