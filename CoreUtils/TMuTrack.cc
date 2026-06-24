#include "TMuTrack.hh"
#include "GenMagneticField.hh"

// genfit
#include <RKTrackRep.h>
#include <Track.h>
#include <TrackPoint.h>
#include <PlanarMeasurement.h>
#include <KalmanFitterRefTrack.h>
#include <FitStatus.h>
// ////////    ///////////
#include <cstdlib> // getenv for MS_FORCE_SEED_10GEV toggle
// ////////    ///////////
// ////////    ///////////
#include <FieldManager.h>

#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TGeoMedium.h>
#include <TGeoMaterial.h>
// ////////    ///////////

// in GenFit Momentum in GeV/c, length in cm and magnetic field in kGauss!

ClassImp(TMuTrack)

void TMuTrack::GenFitTrackFit(int verbose, double detectorResolutionPSmm) {

    const int pdg = 13;    

// trackrep
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

    // start values for the fit, e.g. from pattern recognition
    // position of first hit of the track
    TVector3 pos(fpos[0].x()/10.0, fpos[0].y()/10.0, fpos[0].z()/10.0);
    double pmom = 10.0; //10.0*1e3;  // in GeV;
    //
    //////////(o^o)///////////
    // CoreUtils/TMuTrack.cc to force seed momentum to 10 GeV when needed:
    //if (!getenv("MS_FORCE_SEED_10GEV")) {
        // keep existing curvature-based override of pmom
        // (only adjust pmom when env is NOT set)
    //}
    // Export MS_FORCE_SEED_10GEV=1 in the environment to guarantee pmom=10 GeV.
    TVector3 mom(0, 0, pmom);
    // Seed direction from first two hits when available; default is +z
    TVector3 direction(0, 0, 1);
    if (fpos.size() > 1) {
        TVector3 hit1(fpos[0].x(), fpos[0].y(), fpos[0].z());
        TVector3 hit2(fpos[1].x(), fpos[1].y(), fpos[1].z());
        TVector3 diff = hit2 - hit1;
        if (diff.Mag() > 0) direction = diff.Unit();
    }
    mom = direction * pmom;
    if (verbose>0) {
        std::cout << "Seed direction: (" << direction.X() << ", " << direction.Y() << ", " << direction.Z() << ")" << std::endl;
    }
    if (verbose>0) {
        std::cout << "Seed momentum p0 (GeV/c): " << pmom << std::endl;
    }
    //////////(o^o)///////////

    // create track
    fitTrack = new genfit::Track(rep, pos, mom);

    const int detId(0); // detector ID
    int planeId(0); // detector plane ID
    int hitId(0); // hit ID

 // Scifi resolution: convert from mm to cm (GenFit uses cm)
    // detectorResolutionPSmm is in mm (default 0.1 mm = 100 μm)
    double detectorResolution(detectorResolutionPSmm/10.0); // Convert mm to cm: 0.1 mm → 0.01 cm
    if (verbose > 0) {
        std::cout << "Detector resolution: " << detectorResolutionPSmm << " mm = " 
                  << detectorResolutionPSmm*1000 << " μm (GenFit internal: " 
                  << detectorResolution << " cm)" << std::endl;
    }   
    TMatrixDSym hitCov(2);
    hitCov.UnitMatrix();
    hitCov *= detectorResolution*detectorResolution;

    int nhits = fpos.size();
    // loop over hits
    for (size_t i = 0; i < nhits; i++) {
        TVectorD hitCoords(2);
        hitCoords[0] = fpos[i].x()/10.0;
        hitCoords[1] = fpos[i].y()/10.0;
        genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, 
            detId, ++hitId, nullptr);
        measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,fpos[i].z()/10.0),
                     TVector3(1,0,0), TVector3(0,1,0))), ++planeId);
        fitTrack->insertPoint(new genfit::TrackPoint(measurement, fitTrack));
    }

    fitTrack->checkConsistency();

     // init fitter
    genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();
    fitter->setMaxIterations(10);  // Allow up to 10 iterations for convergence
    fitter->setRelChi2Change(0.001);  // Stop when chi2 changes less than 0.1%

    // do the fit
    try {
        //fitter->processTrack(fitTrack);
        fitter->processTrackWithRep(fitTrack, rep);  // CRITICAL: Use processTrackWithRep, not processTrack
    }    
    catch(genfit::Exception& e){
        std::cerr << e.what();
        std::cerr << "Exception when track fitting with GENFIT" << std::endl;
        delete fitter;
        return;  // Exit on exception  
    }
    // Check if fit converged
    if (!fitTrack->getFitStatus()->isFitConverged()) {
        if (verbose > 0) {
            std::cerr << "WARNING: Fit did not converge!" << std::endl;
        }
        fchi2 = -1;
        fnDoF = 0;
        fpval = 0;
        fp = 0;
        fpx = fpy = fpz = 0;
        fpErr = 0;
        fcharge = 0;
        delete fitter;
        return;
    }

    if(verbose > 3) {
        std::cout << "After fit: " << std::endl;
        //fitTrack->getFittedState().Print();
        fitTrack->getFittedState(0).Print();  // Explicitly get state at first point
        fitTrack->Print();
    }

    // compute momentum at the first point
    double chi2 = fitTrack->getFitStatus()->getChi2();
    double pval = fitTrack->getFitStatus()->getPVal();
    fchi2 = chi2;
    fnDoF = fitTrack->getFitStatus()->getNdf();
    fpval = pval;
    if (verbose>0) {
        std::cout << "Track fit results: chi2 = " << chi2 << " nDoF = " << fnDoF << " pval = " << pval << std::endl;
    }
    // print momentum 
    //genfit::MeasuredStateOnPlane state = fitTrack->getFittedState();
    genfit::MeasuredStateOnPlane state = fitTrack->getFittedState(0); // Explicitly get state at first point
    TVector3 p = state.getMom();
    if (verbose>0) {
        std::cout << "fitted momentum (GeV/c): " << p.Mag() << " px: " << p.X() << " py: " << p.Y() << " pz: " << p.Z();
        std::cout << " pT: " << sqrt(p.X()*p.X()+p.Y()*p.Y()) << std::endl;
    }
    fpx = p.X();
    fpy = p.Y();
    fpz = p.Z();
    fp = p.Mag();
    fcharge = state.getCharge();
    // compute error on momentum 
    fipErr = sqrt(state.getCov()(0,0));
    fpErr = sqrt(state.getMomVar()); 
    if (verbose>0) {
        std::cout << "fitted momentum error (GeV/c): " << fpErr << std::endl;
    }

    // clean up (don't delete rep because it is owned by the track)
    delete fitter;
}

// ////////    ///////////
// Helper: compute circle from three 2D points (returns success flag)
static bool circleFrom3Points(double y1, double z1,
                              double y2, double z2,
                              double y3, double z3,
                              double &cy, double &cz, double &R) {
    // Using perpendicular bisector intersection
    double a = y1 - y2;
    double b = z1 - z2;
    double c = y1 - y3;
    double d = z1 - z3;
    double e = ((y1*y1 - y2*y2) + (z1*z1 - z2*z2)) / 2.0;
    double f = ((y1*y1 - y3*y3) + (z1*z1 - z3*z3)) / 2.0;
    double det = a*d - b*c;
    if (std::fabs(det) < 1e-9) return false;
    cy = (d*e - b*f) / det;
    cz = (-c*e + a*f) / det;
    R = std::sqrt((cy - y1)*(cy - y1) + (cz - z1)*(cz - z1));
    return std::isfinite(R) && R > 0;
}

// Added analytic alternative to GenFit
void TMuTrack::CircleFitTaubin(int verbose, double detectorResolutionPSmm) {
    // Defaults
    fpx = fpy = fpz = fp = 0.0;
    fchi2 = -1.0;
    fnDoF = 0;
    fpval = 0.0;
    fipErr = 0.0;
    fpErr = 0.0;
    // Initialize charge from PDG sign (mu-: 13 -> -1, mu+: -13 -> +1); may refine after fit
    fcharge = (fPDG == 13 ? -1.f : +1.f);

    const size_t N = fpos.size();
    if (N < 3) return;

    // Build y-z points in meters (m) for better conditioning
    std::vector<double> y_m; y_m.reserve(N);
    std::vector<double> z_m; z_m.reserve(N);
    for (size_t i = 0; i < N; ++i) {
        y_m.push_back(fpos[i].y() / 1000.0);
        z_m.push_back(fpos[i].z() / 1000.0);
    }

    // Compute centroid
    double meanY = 0.0, meanZ = 0.0;
    for (size_t i = 0; i < N; ++i) { meanY += y_m[i]; meanZ += z_m[i]; }
    meanY /= (double)N; meanZ /= (double)N;

    // Centered coordinates and moments
    double Suu=0, Svv=0, Suv=0, Suuu=0, Svvv=0, Suuv=0, Suvv=0;
    for (size_t i = 0; i < N; ++i) {
        double u = y_m[i] - meanY;
        double v = z_m[i] - meanZ;
        double uu = u*u, vv = v*v;
        Suu  += uu;
        Svv  += vv;
        Suv  += u*v;
        Suuu += uu*u;
        Svvv += vv*v;
        Suuv += uu*v;
        Suvv += u*vv;
    }

    // Solve linear system for circle center in (u,v): A * [uc; vc] = b
    double A11 = 2*Suu, A12 = 2*Suv;
    double A21 = 2*Suv, A22 = 2*Svv;
    double b1  = Suuu + Suvv;
    double b2  = Svvv + Suuv;
    double det = A11*A22 - A12*A21;

    // Fallback: if nearly singular, try 3-point circle
    double uc=0, vc=0;
    if (std::fabs(det) > 1e-16) {
        uc = ( b1*A22 - b2*A12) / det;
        vc = (-b1*A21 + b2*A11) / det;
    } else {
        // ////////    ///////////
        // Fallback to 3-point circle fit
        size_t i0 = 0, i1 = N/2, i2 = N-1;
        double cy=0, cz=0, Rm=0;
        bool ok = circleFrom3Points(y_m[i0], z_m[i0], y_m[i1], z_m[i1], y_m[i2], z_m[i2], cy, cz, Rm);
        if (!ok) return;
        uc = cy - meanY; vc = cz - meanZ;
    }

    // Circle center in original coordinates
    double yc = meanY + uc;
    double zc = meanZ + vc;
    // Radius estimate
    double R_m = std::sqrt(uc*uc + vc*vc + (Suu + Svv) / (double)N);

    // Residuals and chi2 using sigma from detector resolution
    double sigma_m = std::max(1e-6, detectorResolutionPSmm / 1000.0); // mm -> m, guard
    double chi2 = 0.0; int nd = 0; double sigma_acc = 0.0;
    for (size_t i = 0; i < N; ++i) {
        double dy = y_m[i] - yc;
        double dz = z_m[i] - zc;
        double r  = std::sqrt(dy*dy + dz*dz);
        double res = (r - R_m);
        chi2 += (res*res) / (sigma_m*sigma_m);
        sigma_acc += res*res;
        ++nd;
    }
    fnDoF = std::max(0, nd - 3);
    fchi2 = chi2;
    // RMS of radial residuals (meters)
    double sigma_R_m = std::sqrt(sigma_acc / std::max(1, fnDoF));

    // Local field at middle point (Tesla, assume dominant Bx)
    genfit::AbsBField* field = genfit::FieldManager::getInstance()->getField();
    double Bx_T = 0.0, Bmag_T = 0.0;
    if (field) {
        size_t imid = N/2;
        TVector3 pos_cm(fpos[imid].x()/10.0, fpos[imid].y()/10.0, fpos[imid].z()/10.0);
        TVector3 BkG = field->get(pos_cm);
        Bx_T = BkG.X() / 10.0;
        Bmag_T = std::sqrt(BkG.X()*BkG.X()+BkG.Y()*BkG.Y()+BkG.Z()*BkG.Z())/10.0;
    }
    double Beff_T = (Bmag_T > 1e-4 ? std::fabs(Bx_T) : 1.5); // if field off, assume 1.5 T

    // Momentum and its simple error propagation from radius
    double p_mag = 0.3 * Beff_T * R_m;           // GeV/c
    double p_err = 0.3 * Beff_T * sigma_R_m;     // GeV/c
    fp = p_mag; fpErr = p_err;

    // Charge from sagitta sign and Bx sign
    int charge_sign = (fPDG == 13 ? -1 : +1);
    if (N >= 3 && std::fabs(Bx_T) > 1e-3) {
        double y1 = y_m.front();
        double y2 = y_m[N/2];
        double y3 = y_m.back();
        double sagitta_y = y2 - 0.5*(y1 + y3);
        charge_sign = (sagitta_y * Bx_T > 0) ? -1 : +1;
    }
    fcharge = (float)charge_sign;

    // Tangent direction at middle (approx with neighbor points)
    TVector3 dir(0, 0, 1);
    size_t i_mid = N/2;
    size_t i_prev = (i_mid > 0 ? i_mid - 1 : i_mid);
    size_t i_next = (i_mid + 1 < N ? i_mid + 1 : i_mid);
    TVector3 dmm(0, (fpos[i_next].y() - fpos[i_prev].y()), (fpos[i_next].z() - fpos[i_prev].z()));
    if (dmm.Mag() > 0) dir = dmm.Unit();

    // 3D momentum components; assume negligible px bending
    fpx = 0.0;
    fpy = dir.Y() * p_mag;
    fpz = dir.Z() * p_mag;

    // Simple p-value heuristic by reduced chi2
    fpval = (fnDoF > 0 && (fchi2 / fnDoF) < 20.0) ? 0.5 : 0.0;

    if (verbose>0) {
        std::cout << "[Taubin] yc= " << yc << " m, zc= " << zc << " m, R= " << R_m
                  << " m, |Bx|~" << Beff_T << " T, p ≈ " << p_mag << " ± " << p_err
                  << " GeV/c, chi2/ndf = " << (fnDoF>0? fchi2/fnDoF : fchi2)
                  << ", q= " << fcharge << std::endl;
    }
}
// ////////    ///////////
bool TMuTrack::GenFitMDTFit(const std::vector<MDTMeas>& meas,
                            int pdg, 
                            double seedMomentumGeV,
                            int verbose)
{
    if (meas.size() < 5) {
        std::cerr << "[GenFitMDTFit] not enough MDT measurements: " << meas.size() << std::endl;
        return false;
    }
    // Sort measurements along z.
    // This avoids relying on the input order.
    std::vector<MDTMeas> sortedMeas = meas;
    std::sort(sortedMeas.begin(), sortedMeas.end(),
              [](const MDTMeas& a, const MDTMeas& b) {
                  return a.z_mm < b.z_mm;
              });

    int fitPdg = pdg;
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(fitPdg);

    // Seed position: first measurement, in cm.
    TVector3 p0(sortedMeas.front().x_mm / 10.0,
                sortedMeas.front().y_mm / 10.0,
                sortedMeas.front().z_mm / 10.0);

    // Build a local straight seed from the first few hits, not from the full bent track.
    const size_t nSeed = std::min<size_t>(8, sortedMeas.size());

    double sumZ = 0.0, sumX = 0.0, sumY = 0.0;
    double sumZZ = 0.0, sumZX = 0.0, sumZY = 0.0;
    for (size_t i = 0; i < nSeed; ++i) {
        const double z = sortedMeas[i].z_mm / 10.0; // cm
        const double x = sortedMeas[i].x_mm / 10.0; // cm
        const double y = sortedMeas[i].y_mm / 10.0; // cm

        sumZ  += z;
        sumX  += x;
        sumY  += y;
        sumZZ += z * z;
        sumZX += z * x;
        sumZY += z * y;
    }

    const double denom = nSeed * sumZZ - sumZ * sumZ;

    double tx = 0.0; // dx/dz
    double ty = 0.0; // dy/dz

    if (std::fabs(denom) > 1e-12) {
        tx = (nSeed * sumZX - sumZ * sumX) / denom;
        ty = (nSeed * sumZY - sumZ * sumY) / denom;
    } else {
        std::cerr << "[GenFitMDTFit] bad seed line denominator" << std::endl;
        delete rep;
        return false;
    }

    // Since the particle goes mainly along +z:
    TVector3 dir(tx, ty, 1.0);
    if (dir.Mag() <= 0.0) {
        std::cerr << "[GenFitMDTFit] bad seed direction" << std::endl;
        delete rep;
        return false;
    }

    dir = dir.Unit();

    TVector3 momSeed = dir * seedMomentumGeV;

    if (verbose > 0) {
        std::cout << "[GenFitMDTFit] gGeoManager = " << gGeoManager << std::endl;
        if (gGeoManager && gGeoManager->GetTopVolume()) {
            std::cout << "[GenFitMDTFit] Geometry top volume = "
                    << gGeoManager->GetTopVolume()->GetName()
                    << std::endl;
        }
        auto* field = genfit::FieldManager::getInstance()->getField();
        std::cout << "[GenFitMDTFit] Field pointer = " << field << std::endl;
        for (const auto& m : meas) {
            TVector3 pos_cm(m.x_mm/10.0, m.y_mm/10.0, m.z_mm/10.0);
            TVector3 BkG = field ? field->get(pos_cm) : TVector3(0,0,0);
            TGeoNode* node = nullptr;
            const char* matName = "NULL";
            if (gGeoManager) {
                node = gGeoManager->FindNode(pos_cm.X(), pos_cm.Y(), pos_cm.Z());
                if (node && node->GetVolume() && node->GetVolume()->GetMedium()) {
                    matName = node->GetVolume()->GetMedium()->GetMaterial()->GetName();
                }
            }
            std::cout << Form(
                "[GenFitMDTFit] z=%+.1f mm  B=(%+.2f,%+.2f,%+.2f) kG  material=%s",
                m.z_mm, BkG.X(), BkG.Y(), BkG.Z(), matName
            ) << std::endl;
        }
    }

    std::cout << "[GenFitMDTFit] seed p = " << momSeed.Mag() << " GeV\n";
    std::cout << "[GenFitMDTFit] seed pos cm = "
          << p0.X() << " "
          << p0.Y() << " "
          << p0.Z() << "\n";
    std::cout << "[GenFitMDTFit] seed mom GeV = "
          << momSeed.X() << " "
          << momSeed.Y() << " "
          << momSeed.Z() << "\n";

    // create track
    genfit::Track* trk = new genfit::Track(rep, p0, momSeed);
    // Measurement covariance.
    //
    // For Bx field, bending is in y(z).
    // y is the precision coordinate.
    // x is along the tube/wire direction and should be weakly constrained.    const double sigmaDriftMm = 0.080;
    const double sigmaDriftMm = 0.080;       // MDT precision, 80 um
    const double sigmaY_cm = 0.1; //sigmaDriftMm / 10.0;
    const double sigmaX_cm = 10.0;          // very weak x constraint

    TMatrixDSym cov(2);
    cov.Zero();
    cov(0,0) = sigmaX_cm * sigmaX_cm; // x is along the tube/wire direction, weakly constrained
    cov(1,1) = sigmaY_cm * sigmaY_cm; // y is the precision coordinate, strongly constrained

    for (size_t i = 0; i < sortedMeas.size(); ++i) {
        TVectorD hitCoords(2);

        // Local plane coordinates:
        // coordinate 0 -> x, weak
        // coordinate 1 -> y, precision/bending coordinate
        hitCoords[0] = sortedMeas[i].x_mm / 10.0;
        hitCoords[1] = sortedMeas[i].y_mm / 10.0;

        auto* measurement = new genfit::PlanarMeasurement(
            hitCoords, cov, 0, static_cast<int>(i), nullptr);

        genfit::SharedPlanePtr plane(new genfit::DetPlane(
            TVector3(0.0, 0.0, sortedMeas[i].z_mm / 10.0),
            TVector3(1.0, 0.0, 0.0),   // local u = global x
            TVector3(0.0, 1.0, 0.0)    // local v = global y
        ));

        measurement->setPlane(plane, static_cast<int>(i));
        trk->insertPoint(new genfit::TrackPoint(measurement, trk));
    }
    trk->checkConsistency();

    auto* field = genfit::FieldManager::getInstance()->getField();

    if (verbose > 0 && field) {
        std::cout << "\n[GenFitMDTFit] Field scan between first and last MDT hits\n";

        const double x_cm = p0.X();
        const double y_cm = p0.Y();
        const double zMin_cm = sortedMeas.front().z_mm / 10.0;
        const double zMax_cm = sortedMeas.back().z_mm  / 10.0;

        for (double z_cm = zMin_cm; z_cm <= zMax_cm; z_cm += 5.0) {
            TVector3 pos_cm(x_cm, y_cm, z_cm);
            TVector3 BkG = field->get(pos_cm);

            if (std::fabs(BkG.X()) > 1e-3 ||
                std::fabs(BkG.Y()) > 1e-3 ||
                std::fabs(BkG.Z()) > 1e-3) {
                std::cout << Form(
                    "[GenFitMDTFit] FIELD z=%+.1f cm  B=(%+.3f,%+.3f,%+.3f) kG",
                    z_cm, BkG.X(), BkG.Y(), BkG.Z()
                ) << std::endl;
            }
        }
        std::cout << "[GenFitMDTFit] End field scan\n\n";
    }

    // init fitter
    // Use KalmanFitterRefTrack for reference track fitting
    genfit::KalmanFitterRefTrack fitter;
    fitter.setMaxIterations(50);
    fitter.setRelChi2Change(1e-4);

    try {
        //fitter.processTrack(trk);
        fitter.processTrackWithRep(trk, rep);
    }
    catch (genfit::Exception& e) {
        std::cerr << "[GenFitMDTFit] GenFit exception: "
                  << e.what() << std::endl;
        delete trk;
        return false;
    }
    // Check fit status
    genfit::FitStatus* status = trk->getFitStatus(rep);

    if (!status || !status->isFitConverged()) {
        std::cerr << "[GenFitMDTFit] fit did not converge" << std::endl;
        delete trk;
        return false;
    }
    // Extract fitted state at the first measurement
    genfit::MeasuredStateOnPlane state = trk->getFittedState();
    TVector3 pfit = state.getMom();
    fpx = pfit.X();
    fpy = pfit.Y();
    fpz = pfit.Z();
    fp  = pfit.Mag();
    fcharge = state.getCharge();
    fchi2   = status->getChi2();
    fnDoF   = status->getNdf();
    fpval   = status->getPVal();
    fitTrack = trk;

    fpos.clear();
    layerID.clear();
    for (const auto& m : sortedMeas) {
        fpos.emplace_back(m.x_mm, m.y_mm, m.z_mm);
        layerID.push_back(
            m.stationID * 10000
          + m.planeID   * 1000
          + m.tubeID
        );
    }
    if (verbose) {
        std::cout << "[GenFitMDTFit] fitted momentum:"
                  << " px=" << fpx
                  << " py=" << fpy
                  << " pz=" << fpz
                  << " p="  << fp
                  << " q="  << fcharge
                  << std::endl;

        std::cout << "[GenFitMDTFit] fitted slopes:"
                  << " dx/dz=" << fpx / fpz
                  << " dy/dz=" << fpy / fpz
                  << std::endl;
    }

    return true;
}





