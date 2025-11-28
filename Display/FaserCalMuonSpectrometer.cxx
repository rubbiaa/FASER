#include "FaserCalDisplay.h" 
#include <TFile.h>
#include <TTree.h>
#include <TDatabasePDG.h>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <sstream>
#include "TChain.h"
#include <string>
#include <random>
#include <cmath>
#include <set>

#include "TcalEvent.hh"
#include "TPOEvent.hh"
#include "TPORecoEvent.hh"
#include "TMuTrack.hh"
#include "GenMagneticField.hh"

// GENFIT includes
#include <ConstField.h>
#include <FieldManager.h>
#include <MaterialEffects.h>
#include <TGeoMaterialInterface.h>

struct MuonRecoData {
    int trackID;
    int pdg;
    double truth_px, truth_py, truth_pz, truth_p;
    // Backward-compat single-set (will mirror GenFit results)
    double reco_px, reco_py, reco_pz, reco_p;
    int nhits, nstations;
    double chi2, ndf;
    bool fit_success;
    // GenFit results
    double genfit_px, genfit_py, genfit_pz, genfit_p;
    double genfit_chi2, genfit_ndf, genfit_pval;
    bool genfit_fit_success;
    // Taubin results
    double taubin_px, taubin_py, taubin_pz, taubin_p;
    double taubin_p_err;
    double taubin_chi2, taubin_ndf, taubin_pval;
    int taubin_charge;
    bool taubin_fit_success;
};

std::vector<MuonRecoData> ReconstructMuonSpectrometer(display::FaserCalDisplay* display, int eventNumber) 
  {
    std::vector<MuonRecoData> recoResults;
    // Initialize GenFit if not already done
    static bool genfit_initialized = false;
    if (!genfit_initialized) {
        std::cout << "Initializing GenFit for muon track fitting..." << std::endl;
        genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
        // Initialize magnetic field 
        GenMagneticField* magField = new GenMagneticField();
        //////////(o^o)///////////
        // Build geometry-driven field maps so B=0 at stations and |B|=1.5T in between
        int venv = 0; if (const char* s = std::getenv("MS_VERBOSE")) { try { venv = std::stoi(s); } catch(...) { venv = 1; } }
        magField->BuildGeometryMaps(venv);
        //////////(o^o)///////////
        genfit::FieldManager::getInstance()->init(magField);
        
        genfit_initialized = true;
        std::cout << "GenFit initialization complete." << std::endl;
    }

    std::cout << " ------------ MUTAG truth tracks" << std::endl;
    for (const auto &it : display->fTcalEvent->fMuTagTracks) 
    {
      std::cout << it->ftrackID << " " << it->fPDG << " Mom:" << it->mom[0].x() << " " << it->mom[0].y() << " " << it->mom[0].z() << std::endl;
    }

    std::vector<TMuTrack> fMuTracks;

   // Build temp tracks: smear each raw hit by 100 μm (x,y) BEFORE station averaging
    std::vector<TMuTrack> tempMuTracks;
    {
      std::random_device rd_hits;
      std::mt19937 gen_hits(rd_hits());
      std::normal_distribution<> smear_xy_mm(0.0, 0.1); // 0.1 mm = 100 μm
      for (const auto &it : display->fTcalEvent->fMuTagTracks) 
      {
        if(it->fPDG != 13 && it->fPDG != -13) continue;
        TMuTrack mutrack;
        mutrack.ftrackID = it->ftrackID;
        mutrack.fPDG = it->fPDG;
        mutrack.layerID = it->layerID;
        mutrack.fpos.reserve(it->pos.size());
        for (size_t ih = 0; ih < it->pos.size(); ++ih) {
          ROOT::Math::XYZVector p = it->pos[ih];
          p.SetX(p.x() + smear_xy_mm(gen_hits));
          p.SetY(p.y() + smear_xy_mm(gen_hits));
          // z is not smeared
          mutrack.fpos.push_back(p);
        }
        tempMuTracks.push_back(std::move(mutrack));
      }
    }
    // now loop over tempMuTracks and create MuTrack with average hits in the same station
    fMuTracks.clear();
    for(auto &it : tempMuTracks) {
        TMuTrack mutrack;
        mutrack.ftrackID = it.ftrackID;
        mutrack.fPDG = it.fPDG;
        // loop over all hits and average position over hits in the same station
        size_t nhits = it.fpos.size();
        int currentStation = -1;
        int nstation = 0;
        int hitsInCurrentStation = 0;
        for(size_t i = 0; i < nhits; i++) {
            int stationID = it.layerID[i] / 4;
            ROOT::Math::XYZVector pos = it.fpos[i];
            if(stationID == currentStation) {
                // average position with existing position in current station
                hitsInCurrentStation++;
                std::cout << "  Averaging hit " << hitsInCurrentStation << " in station " << stationID 
                         << " for track " << mutrack.ftrackID << std::endl;
                std::cout << "    Before averaging: (" << mutrack.fpos[nstation-1].x() << ", " 
                         << mutrack.fpos[nstation-1].y() << ", " << mutrack.fpos[nstation-1].z() << ")" << std::endl;
                std::cout << "    New hit position: (" << pos.x() << ", " << pos.y() << ", " << pos.z() << ")" << std::endl;
                
                mutrack.fpos[nstation-1].SetX( (mutrack.fpos[nstation-1].x()*float(hitsInCurrentStation-1) + pos.x())/float(hitsInCurrentStation) );
                mutrack.fpos[nstation-1].SetY( (mutrack.fpos[nstation-1].y()*float(hitsInCurrentStation-1) + pos.y())/float(hitsInCurrentStation) );
                mutrack.fpos[nstation-1].SetZ( (mutrack.fpos[nstation-1].z()*float(hitsInCurrentStation-1) + pos.z())/float(hitsInCurrentStation) );
                
                std::cout << "    After averaging:  (" << mutrack.fpos[nstation-1].x() << ", " 
                         << mutrack.fpos[nstation-1].y() << ", " << mutrack.fpos[nstation-1].z() << ")" << std::endl;
            } else {
                // new station - add the first hit position
                std::cout << "  New station " << stationID << " for track " << mutrack.ftrackID 
                         << " at position (" << pos.x() << ", " << pos.y() << ", " << pos.z() << ")" << std::endl;
                mutrack.fpos.push_back(pos);
                mutrack.layerID.push_back(stationID * 4);
                currentStation = stationID;
                nstation++;
                hitsInCurrentStation = 1;
            }
        }
        fMuTracks.push_back(mutrack);        
        // Summary of averaging for this track
        std::cout << "Track " << mutrack.ftrackID << " averaging complete: " 
                 << mutrack.fpos.size() << " averaged positions from " << nhits << " original hits" << std::endl;
        for (size_t j = 0; j < mutrack.fpos.size(); j++) {
            int stationID = mutrack.layerID[j] / 4;
            std::cout << "  Station " << stationID << ": averaged position (" 
                     << mutrack.fpos[j].x() << ", " << mutrack.fpos[j].y() << ", " << mutrack.fpos[j].z() << ")" << std::endl;
        }
    }
  // No smearing here: already applied per-hit before averaging
  // Field unit sanity check at first hit (Tesla) --> GenFit expects kGauss; convert to Tesla for print.
  if (!fMuTracks.empty() && !fMuTracks.front().fpos.empty()) {
    auto &trk0 = fMuTracks.front();
    TVector3 pos_cm(trk0.fpos.front().x()/10.0, trk0.fpos.front().y()/10.0, trk0.fpos.front().z()/10.0);
    genfit::AbsBField* field = genfit::FieldManager::getInstance()->getField();
    //////////(o^o)///////////
    // Provide per-event station Zs to the magnetic field for correct gating when geometry is misaligned
    if (auto* gfield = dynamic_cast<GenMagneticField*>(field)) {
      std::vector<double> zs_cm;
      // Use the first muon track's averaged station positions as representatives
      for (const auto& p : trk0.fpos) zs_cm.push_back(p.z()/10.0);
      std::sort(zs_cm.begin(), zs_cm.end());
      gfield->SetEventStationZsCm(zs_cm);
    }
    //////////(o^o)///////////
    TVector3 BkG = field->get(pos_cm);
    std::cout << "[Field] First-hit field (Tesla): Bx=" << BkG.X()/10.0
          << " By=" << BkG.Y()/10.0
          << " Bz=" << BkG.Z()/10.0 << std::endl;
  }
    // No smearing here: already applied per-hit before averaging
    // dump fMuTracks
    // Verbosity for GenFit prints: set MS_VERBOSE=1 (or higher) in env to enable
    int verbose = 0;
    if (const char* v = std::getenv("MS_VERBOSE")) {
      try { verbose = std::stoi(v); } catch (...) { verbose = 1; }
    }
    //////////(o^o)///////////
    // just for diagnostic: compute sum B_perp · dl along the averaged hit polyline and a quick p estimate
    int diag_bdl = 0; if (const char* d = std::getenv("MS_DIAG_BDL")) { try { diag_bdl = std::stoi(d); } catch(...) { diag_bdl = 1; } }
    auto computeSigmaBdl = [&](const TMuTrack& trk) {
      struct BdlResult { double bdl_kG_cm{0.0}; double bdl_T_m{0.0}; double theta_rad{0.0}; double p_est_GeV{0.0}; } res;
      if (trk.fpos.size() < 2) return res;
      genfit::AbsBField* field = genfit::FieldManager::getInstance()->getField();
      // Deflection angle from first and last segment
      auto vdir = [&](size_t i0, size_t i1) {
        ROOT::Math::XYZVector d = trk.fpos[i1] - trk.fpos[i0];
        double n = d.R();
        if (n == 0) return ROOT::Math::XYZVector(0,0,1);
        return (1.0/n) * d; // unit
      };
      ROOT::Math::XYZVector u_in  = vdir(0, 1);
      ROOT::Math::XYZVector u_out = vdir(trk.fpos.size()-2, trk.fpos.size()-1);
      double dot = u_in.Dot(u_out);
      dot = std::max(-1.0, std::min(1.0, dot));
      res.theta_rad = std::acos(dot);
      // Integrate along polyline
      for (size_t i = 0; i + 1 < trk.fpos.size(); ++i) {
        ROOT::Math::XYZVector p0_mm = trk.fpos[i];
        ROOT::Math::XYZVector p1_mm = trk.fpos[i+1];
        ROOT::Math::XYZVector seg = p1_mm - p0_mm;
        double seg_len_mm = seg.R();
        if (seg_len_mm <= 0) continue;
        int nSteps = std::max(1, (int)std::ceil(seg_len_mm / 2.0)); // about 2 mm per step
        ROOT::Math::XYZVector u_seg = (1.0/seg_len_mm) * seg; // unit in mm basis
        for (int s = 0; s < nSteps; ++s) {
          double t0 = (double)s / nSteps;
          double t1 = (double)(s+1) / nSteps;
          // Midpoint for field query (convert to cm)
          ROOT::Math::XYZVector pmid_mm = p0_mm + (t0 + t1)*0.5 * seg;
          TVector3 pos_cm(pmid_mm.x()/10.0, pmid_mm.y()/10.0, pmid_mm.z()/10.0);
          TVector3 B_kG = field->get(pos_cm);
          // Convert to XYZVector for cross product with direction
          ROOT::Math::XYZVector B(B_kG.X(), B_kG.Y(), B_kG.Z());
          // v_hat ~ segment direction in mm; use same unitless direction
          ROOT::Math::XYZVector vhat = u_seg;
          ROOT::Math::XYZVector BcrossV = ROOT::Math::XYZVector(
            B.y()*vhat.z() - B.z()*vhat.y(),
            B.z()*vhat.x() - B.x()*vhat.z(),
            B.x()*vhat.y() - B.y()*vhat.x());
          double Bperp_kG = BcrossV.R();
          double dl_cm = (t1 - t0) * seg_len_mm / 10.0; // step length in cm
          res.bdl_kG_cm += Bperp_kG * dl_cm;
        }
      }
      res.bdl_T_m = res.bdl_kG_cm * 1e-3; // 1 kG*cm = 1e-3 T*m
      if (res.theta_rad > 1e-6) res.p_est_GeV = 0.3 * (res.bdl_T_m) / res.theta_rad;
      return res;
    };
    std::cout << " ------------ Spectrometer reconstructed muon tracks" << std::endl;
    for (const auto &it : fMuTracks) 
    {
      std::cout << it.ftrackID << " " << it.fPDG << " Smearred Positions: " << it.fpos.size() << " hits" << std::endl;
      for (size_t i = 0; i < it.fpos.size(); i++) {
        std::cout << "   Hit " << i << " LayerID: " << it.layerID[i] << " Pos: (" << it.fpos[i].x() << "," << it.fpos[i].y() << "," << it.fpos[i].z() << ")" << std::endl;
      }
    }

  for(auto &it : fMuTracks) {
        // Count unique stations for this track
        std::set<int> stationSet;
        for (size_t i = 0; i < it.layerID.size(); i++) {
            int stationID = it.layerID[i] / 4;
            stationSet.insert(stationID);
        }
        int nstations = stationSet.size();
        //
        // Fitting if we have at least 3 stations
        bool attemptFit = (nstations >= 3);
        
        // Containers for both fits
        double g_px=0, g_py=0, g_pz=0, g_p=0, g_chi2=-1, g_ndf=0, g_pval=0; bool g_ok=false;
        double t_px=0, t_py=0, t_pz=0, t_p=0, t_perr=0, t_chi2=-1, t_ndf=0, t_pval=0; int t_q=0; bool t_ok=false;

        if (attemptFit) {
          std::cout << "Fitting track " << it.ftrackID << " with " << nstations << " stations..." << std::endl;
          // Run Taubin first
          it.CircleFitTaubin(verbose, 0.1);
          t_px = it.fpx; t_py = it.fpy; t_pz = it.fpz; t_p = it.fp; t_perr = it.fpErr;
          t_chi2 = it.fchi2; t_ndf = it.fnDoF; t_pval = it.fpval; t_q = (int)it.fcharge;
          t_ok = (t_ndf > 0 ? (t_chi2 / t_ndf) < 20.0 : true);
          // Run GenFit next
          it.GenFitTrackFit(verbose, 0.1);
          g_px = it.fpx; g_py = it.fpy; g_pz = it.fpz; g_p = it.fp;
          g_chi2 = it.fchi2; g_ndf = it.fnDoF; g_pval = it.fpval;
          g_ok = (it.fpval > 0.01);
          if (diag_bdl > 0) {
            auto res = computeSigmaBdl(it);
            std::cout << "[Diag] Σ B_perp·dl = " << res.bdl_T_m << " T·m (" << res.bdl_kG_cm
                    << " kG·cm), deflection θ = " << res.theta_rad
                    << " rad, p_est ≈ " << res.p_est_GeV << " GeV/c" << std::endl;
          }
        } else {
            std::cout << "Skipping track " << it.ftrackID << " - only " << nstations << " stations (minimum 3 required)" << std::endl;
        }
        // Collect reconstruction data
        MuonRecoData recoData;
        recoData.trackID = it.ftrackID;
        recoData.pdg = it.fPDG;
        
        // Get truth momentum from original tracks
        for (const auto &truthTrack : display->fTcalEvent->fMuTagTracks) {
            if (truthTrack->ftrackID == it.ftrackID) {
        // Convert MeV/c to GeV/c for truth
        recoData.truth_px = truthTrack->mom[0].x() / 1000.0;
        recoData.truth_py = truthTrack->mom[0].y() / 1000.0;
        recoData.truth_pz = truthTrack->mom[0].z() / 1000.0;
                recoData.truth_p = sqrt(recoData.truth_px*recoData.truth_px + 
                                       recoData.truth_py*recoData.truth_py + 
                                       recoData.truth_pz*recoData.truth_pz);
                break;
            }
        }
        
        // Store hit information
        recoData.nhits = it.fpos.size();
        recoData.nstations = nstations;
        
    // Store fit results (only valid if fit was attempted). Mirror GenFit into legacy reco_* fields
    if (attemptFit) {
      recoData.genfit_px = g_px;
      recoData.genfit_py = g_py;
      recoData.genfit_pz = g_pz;
      recoData.genfit_p  = g_p;
      recoData.genfit_chi2 = g_chi2;
      recoData.genfit_ndf  = g_ndf;
      recoData.genfit_pval = g_pval;
      recoData.genfit_fit_success = g_ok;

      recoData.taubin_px = t_px;
      recoData.taubin_py = t_py;
      recoData.taubin_pz = t_pz;
      recoData.taubin_p  = t_p;
      recoData.taubin_p_err = t_perr;
      recoData.taubin_chi2 = t_chi2;
      recoData.taubin_ndf  = t_ndf;
      recoData.taubin_pval = t_pval;
      recoData.taubin_charge = t_q;
      recoData.taubin_fit_success = t_ok;

      // Legacy fields mirror GenFit
      recoData.reco_px = g_px;
      recoData.reco_py = g_py;
      recoData.reco_pz = g_pz;
      recoData.reco_p  = g_p;
      recoData.chi2 = g_chi2;
      recoData.ndf  = g_ndf;
      recoData.fit_success = g_ok;
    } else {
            // Set default values for tracks that weren't fitted
            recoData.reco_px = 0.0;
            recoData.reco_py = 0.0;
            recoData.reco_pz = 0.0;
            recoData.reco_p = 0.0;
            recoData.chi2 = -1.0;
            recoData.ndf = 0;
            recoData.fit_success = false;
      recoData.genfit_px = recoData.genfit_py = recoData.genfit_pz = recoData.genfit_p = 0.0;
      recoData.genfit_chi2 = -1.0; recoData.genfit_ndf = 0; recoData.genfit_pval = 0.0; recoData.genfit_fit_success = false;
      recoData.taubin_px = recoData.taubin_py = recoData.taubin_pz = recoData.taubin_p = 0.0;
      recoData.taubin_p_err = 0.0; recoData.taubin_chi2 = -1.0; recoData.taubin_ndf = 0; recoData.taubin_pval = 0.0; recoData.taubin_fit_success = false; recoData.taubin_charge = 0;
        }
        
        recoResults.push_back(recoData);
    }
  std::cout << " ------------ MUTAG reconstructed muon tracks (Event " << eventNumber << ")" << std::endl;
for (const auto &it : recoResults) 
{
  std::cout << "Event " << eventNumber << " - Muon track info: "
    << "trackID: " << it.trackID << " "
    << "pdg: " << it.pdg << " "
    << "reco_px: " << it.reco_px << " GeV/c "
    << "reco_py: " << it.reco_py << " GeV/c "
    << "reco_pz: " << it.reco_pz << " GeV/c "
    << "reco_p: " << it.reco_p << " GeV/c "
    << "truth_p: " << it.truth_p << " GeV/c "
          << "chi2: " << it.chi2 << " "
          << "ndf: " << it.ndf << " "
          << "nhits: " << it.nhits << " "
          << "nstations: " << it.nstations << " "
          << "fit_success: " << it.fit_success << " "
          << std::endl;
}
    return recoResults;
}


void LoadAllEvents(display::FaserCalDisplay* display, int runNumber, int maxEvents, std::string mask_str)
{
  
  //const std::string& input_folder_path = "input/";
  std::string input_file_path;
  // Variables to store in the tree
  int eventNumber = 0;

  TDatabasePDG* pdgDB = TDatabasePDG::Instance();
  display->AddCustomNucleusParticles(); 

  const std::string& input_folder_path = "input/";
  std::vector<std::string> file_paths;
  int cnt = 0;

  std::vector<double> fx, fy, fz;
  std::vector<double> fpx, fpy, fpz;
  std::vector<int> fpdg, ftrackID, flayerID, fstationID;
  int feventID;

  // Muon reconstruction variables
  std::vector<int> fmuon_trackID, fmuon_pdg;
  std::vector<double> fmuon_truth_px, fmuon_truth_py, fmuon_truth_pz, fmuon_truth_p;
  std::vector<double> fmuon_reco_px, fmuon_reco_py, fmuon_reco_pz, fmuon_reco_p;
  std::vector<int> fmuon_nhits, fmuon_nstations;
  std::vector<double> fmuon_chi2, fmuon_ndf;
  std::vector<bool> fmuon_fit_success;

  TFile *RootFile = new TFile("scifi_hits_all.root", "RECREATE");
  TTree *fTree = new TTree("Hits", "SciFi Hits");

  fTree->Branch("eventID", &feventID);
  fTree->Branch("trackID", &ftrackID);
  fTree->Branch("pdg", &fpdg);
  fTree->Branch("stationID", &fstationID);
  fTree->Branch("layerID", &flayerID);
  fTree->Branch("x", &fx);
  fTree->Branch("y", &fy);
  fTree->Branch("z", &fz);
  fTree->Branch("px", &fpx);
  fTree->Branch("py", &fpy);
  fTree->Branch("pz", &fpz);

  // Muon reconstruction branches
  fTree->Branch("muon_trackID", &fmuon_trackID);
  fTree->Branch("muon_pdg", &fmuon_pdg);
  fTree->Branch("muon_truth_px", &fmuon_truth_px);
  fTree->Branch("muon_truth_py", &fmuon_truth_py);
  fTree->Branch("muon_truth_pz", &fmuon_truth_pz);
  fTree->Branch("muon_truth_p", &fmuon_truth_p);
  fTree->Branch("muon_reco_px", &fmuon_reco_px);
  fTree->Branch("muon_reco_py", &fmuon_reco_py);
  fTree->Branch("muon_reco_pz", &fmuon_reco_pz);
  fTree->Branch("muon_reco_p", &fmuon_reco_p);
  fTree->Branch("muon_nhits", &fmuon_nhits);
  fTree->Branch("muon_nstations", &fmuon_nstations);
  fTree->Branch("muon_chi2", &fmuon_chi2);
  fTree->Branch("muon_ndf", &fmuon_ndf);
  fTree->Branch("muon_fit_success", &fmuon_fit_success);

  // New explicit branches for GenFit
  std::vector<double> fgenfit_px, fgenfit_py, fgenfit_pz, fgenfit_p;
  std::vector<double> fgenfit_chi2, fgenfit_ndf, fgenfit_pval;
  std::vector<bool>   fgenfit_fit_success;
  fTree->Branch("genfit_px", &fgenfit_px);
  fTree->Branch("genfit_py", &fgenfit_py);
  fTree->Branch("genfit_pz", &fgenfit_pz);
  fTree->Branch("genfit_p",  &fgenfit_p);
  fTree->Branch("genfit_chi2", &fgenfit_chi2);
  fTree->Branch("genfit_ndf",  &fgenfit_ndf);
  fTree->Branch("genfit_pval", &fgenfit_pval);
  fTree->Branch("genfit_fit_success", &fgenfit_fit_success);

  // New explicit branches for Taubin
  std::vector<double> ftaubin_px, ftaubin_py, ftaubin_pz, ftaubin_p, ftaubin_p_err;
  std::vector<double> ftaubin_chi2, ftaubin_ndf, ftaubin_pval;
  std::vector<int>    ftaubin_charge;
  std::vector<bool>   ftaubin_fit_success;
  fTree->Branch("taubin_px", &ftaubin_px);
  fTree->Branch("taubin_py", &ftaubin_py);
  fTree->Branch("taubin_pz", &ftaubin_pz);
  fTree->Branch("taubin_p",  &ftaubin_p);
  fTree->Branch("taubin_p_err", &ftaubin_p_err);
  fTree->Branch("taubin_chi2", &ftaubin_chi2);
  fTree->Branch("taubin_ndf",  &ftaubin_ndf);
  fTree->Branch("taubin_pval", &ftaubin_pval);
  fTree->Branch("taubin_charge", &ftaubin_charge);
  fTree->Branch("taubin_fit_success", &ftaubin_fit_success);


  TGeoManager::Import("../../GeomGDML/geometry.gdml");

  // Get all file paths in the input directory
  for (const auto& entry : std::filesystem::directory_iterator(input_folder_path))
    {
      file_paths.push_back(entry.path().string());
    }
  
  int eventsProcessed = 0;
  bool done = false;
  for (const std::string& file_path : file_paths)
    {
      if (done) break;
      std::cout << "Processing file: " << file_path << std::endl;
      cnt++;
      // Extract the base name from the file path
      std::string base_name = std::filesystem::path(file_path).stem().string();
      std::cout << "basename " << base_name << std::endl;    
      // Split the base_name using '_' as a delimiter
      std::istringstream ss(base_name);
      std::string token;
      std::vector<std::string> parts;
      while (std::getline(ss, token, '_'))
        {
	        parts.push_back(token);
        }
        if (parts.size() < 3)
          {
	          std::cerr << "Error: Invalid filename format, unable to parse: " << base_name << std::endl;
	          continue; // Skip this file if the format is incorrect
          }
        try 
        {
	        mask_str = (parts.size() == 3) ? "NoMask" : parts[3];
          int fileRun = std::stoi(parts[1]);
          int fileEvent = std::stoi(parts[2]);
          // Only process files that match requested runNumber unless runNumber <= 0 meaning "all runs"
          if (runNumber > 0 && fileRun != runNumber) {
            continue;
          }
          eventNumber = fileEvent;
	        std::cout << cnt << " Event Number: " << eventNumber << ", Run Number: " << runNumber << ", Mask: " << mask_str << std::endl;
	        int imask = 0;
	        if (mask_str == "nueCC") imask = 1;
	        else if (mask_str == "numuCC") imask = 2;
	        else if (mask_str == "nutauCC") imask = 3;
	        else if (mask_str == "nuNC") imask = 4;
	        else if (mask_str == "nuES") imask = 5;
          //TPORecoEvent *fTPORecoEvent = nullptr; // &fTPORecoEvent;
          
          display->fTcalEvent = new TcalEvent();
          display->POevent = new TPOEvent();
          // Use the actual file's run number when loading, even if the CLI runNumber was 0 (ALL)
          display->fTcalEvent->Load_event(input_folder_path, fileRun, eventNumber, imask, display->POevent);
          
          //display->fTcalEvent -> fTPOEvent -> dump_event();

          // Clear vectors for this event to avoid accumulating data across events
          ftrackID.clear();
          fpdg.clear();
          flayerID.clear();
          fstationID.clear();
          fpx.clear();
          fpy.clear();
          fpz.clear();
          fx.clear();
          fy.clear();
          fz.clear();

	        
          for (const auto &it : display->fTcalEvent->fMuTagTracks)
	         {

	            for(size_t i = 0; i < it->pos.size(); ++i)
	              {
		              const auto position = it->pos[i];
		              const auto momentum = it->mom[i];
              		std::cout << "MuSpectrometer hit position " 
                  << eventNumber << " "
                  << it->ftrackID << " "
                  << it->fPDG << " "
                  << it->layerID[i] << " "
                  << momentum.x() << " "
                  << momentum.y() << " "
                  << momentum.z() << " "
                  << position.x() << " "
                  << position.y() << " "
                  << position.z() << " "
                  << std::endl;
                  int stationID = it->layerID[i] / 4 + 1; // from 1 to 10
                  feventID = eventNumber;
                  ftrackID.push_back(it->ftrackID);
                  fpdg.push_back(it->fPDG);
                  flayerID.push_back(it->layerID[i]);
                  fstationID.push_back(stationID);
                  fpx.push_back(momentum.x());
                  fpy.push_back(momentum.y());
                  fpz.push_back(momentum.z());
                  fx.push_back(position.x());
                  fy.push_back(position.y());
                  fz.push_back(position.z());
                }
            }
            ///////////// Reconstruct the muon spectrometer
            std::vector<MuonRecoData> recoResults = ReconstructMuonSpectrometer(display, eventNumber);
            
            // Fill muon reconstruction data
            fmuon_trackID.clear();
            fmuon_pdg.clear();
            fmuon_truth_px.clear();
            fmuon_truth_py.clear();
            fmuon_truth_pz.clear();
            fmuon_truth_p.clear();
            fmuon_reco_px.clear();
            fmuon_reco_py.clear();
            fmuon_reco_pz.clear();
            fmuon_reco_p.clear();
            fmuon_nhits.clear();
            fmuon_nstations.clear();
            fmuon_chi2.clear();
            fmuon_ndf.clear();
            fmuon_fit_success.clear();
            fgenfit_px.clear(); fgenfit_py.clear(); fgenfit_pz.clear(); fgenfit_p.clear();
            fgenfit_chi2.clear(); fgenfit_ndf.clear(); fgenfit_pval.clear(); fgenfit_fit_success.clear();
            ftaubin_px.clear(); ftaubin_py.clear(); ftaubin_pz.clear(); ftaubin_p.clear(); ftaubin_p_err.clear();
            ftaubin_chi2.clear(); ftaubin_ndf.clear(); ftaubin_pval.clear(); ftaubin_charge.clear(); ftaubin_fit_success.clear();
            
      for (const auto& recoData : recoResults) {
                fmuon_trackID.push_back(recoData.trackID);
                fmuon_pdg.push_back(recoData.pdg);
                fmuon_truth_px.push_back(recoData.truth_px);
                fmuon_truth_py.push_back(recoData.truth_py);
                fmuon_truth_pz.push_back(recoData.truth_pz);
    fmuon_truth_p.push_back(recoData.truth_p);
        // Legacy reco mirrors GenFit
        fmuon_reco_px.push_back(recoData.reco_px);
        fmuon_reco_py.push_back(recoData.reco_py);
        fmuon_reco_pz.push_back(recoData.reco_pz);
  fmuon_reco_p.push_back(recoData.reco_p);
                fmuon_nhits.push_back(recoData.nhits);
                fmuon_nstations.push_back(recoData.nstations);
        fmuon_chi2.push_back(recoData.chi2);
        fmuon_ndf.push_back(recoData.ndf);
        fmuon_fit_success.push_back(recoData.fit_success);

        // GenFit results
        fgenfit_px.push_back(recoData.genfit_px);
        fgenfit_py.push_back(recoData.genfit_py);
        fgenfit_pz.push_back(recoData.genfit_pz);
        fgenfit_p.push_back(recoData.genfit_p);
        fgenfit_chi2.push_back(recoData.genfit_chi2);
        fgenfit_ndf.push_back(recoData.genfit_ndf);
        fgenfit_pval.push_back(recoData.genfit_pval);
        fgenfit_fit_success.push_back(recoData.genfit_fit_success);

        // Taubin results
        ftaubin_px.push_back(recoData.taubin_px);
        ftaubin_py.push_back(recoData.taubin_py);
        ftaubin_pz.push_back(recoData.taubin_pz);
        ftaubin_p.push_back(recoData.taubin_p);
        ftaubin_p_err.push_back(recoData.taubin_p_err);
        ftaubin_chi2.push_back(recoData.taubin_chi2);
        ftaubin_ndf.push_back(recoData.taubin_ndf);
        ftaubin_pval.push_back(recoData.taubin_pval);
        ftaubin_charge.push_back(recoData.taubin_charge);
        ftaubin_fit_success.push_back(recoData.taubin_fit_success);
            }
            /////////////
          // Fill the tree once per event (after processing all muon tracks)
          if (!ftrackID.empty()) {
            fTree->Fill();
            eventsProcessed++;
            if (maxEvents > 0 && eventsProcessed >= maxEvents) {
              done = true;
            }
          }
        } 
      catch (const std::exception& e) {
        std::cerr << "Error parsing filename or processing event: " << e.what() << std::endl;
        continue; // Skip this file on error
      }
    } // end of file_paths loop

  // Write and close the ROOT file after processing all events
  RootFile->cd();
  fTree->Write("", TObject::kOverwrite);
  RootFile->Close();
  delete RootFile;
}



int main(int argc, char** argv)
{

  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <RunNumber> <MaxEvents> <Masks>" << std::endl;
    std::cout << "Masks: nueCC, numuCC, nutauCC, nuNC, nuES, -" << std::endl;
    std::cout << "This version includes short-lived particle (charm and tau) analysis" << std::endl;

    return 1;
  }
  
  int runNumber = std::stoi(argv[1]);
  int maxEvents = std::stoi(argv[2]);
  std::string mask_str = argv[3];

  std::cout << "Starting FaserCalMuonSpectrometer for Muon tracks in spectrometer" << std::endl;
  std::cout << "Run: " << (runNumber>0?std::to_string(runNumber):std::string("ALL"))
            << ", Max Events: " << (maxEvents>0?std::to_string(maxEvents):std::string("ALL"))
            << ", Mask: " << mask_str << std::endl;

  display::FaserCalDisplay* disp = new display::FaserCalDisplay();
  LoadAllEvents(disp, runNumber, maxEvents, mask_str);
  delete disp;
  return 0;
}
