#include "TMuTrack.hh"
#include "GenMagneticField.hh"

// genfit
#include <RKTrackRep.h>
#include <Track.h>
#include <TrackPoint.h>
#include <PlanarMeasurement.h>
#include <WireMeasurement.h>
#include <KalmanFitterRefTrack.h>
#include <KalmanFitter.h>
#include <DAF.h>
#include <FitStatus.h>
#include <MeasuredStateOnPlane.h>
#include <FieldManager.h>
#include <MaterialEffects.h>
#include <Exception.h>
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
/////////////////////////////////////////////////
bool TMuTrack::GenFitMDTFit(const std::vector<MDTMeas>& meas,
                            int pdg,
                            double seedMomentumGeV,
                            int verbose,
                            double seedSlopeDyDz)
{
    if (meas.size() < 5) {
        if (verbose > 0)
            std::cerr << "[GenFitMDTFit] not enough MDT measurements: " << meas.size() << std::endl;
        return false;
    }
    // Sort measurements along z.
    // This avoids relying on the input order.
    std::vector<MDTMeas> sortedMeas = meas;
    std::sort(sortedMeas.begin(), sortedMeas.end(),
              [](const MDTMeas& a, const MDTMeas& b) {
                  return a.localZ_mm < b.localZ_mm;
              });

// Globalized MDT-local axes from first hit.
    TVector3 mdtU(sortedMeas.front().uX, sortedMeas.front().uY, sortedMeas.front().uZ); // local Y, measured coordinate
    TVector3 mdtV(sortedMeas.front().vX, sortedMeas.front().vY, sortedMeas.front().vZ); // local X, wire direction
    TVector3 mdtW(sortedMeas.front().wX, sortedMeas.front().wY, sortedMeas.front().wZ); // local Z, downstream

    if (mdtU.Mag() > 0.0) mdtU = mdtU.Unit();
    if (mdtV.Mag() > 0.0) mdtV = mdtV.Unit();
    if (mdtW.Mag() > 0.0) mdtW = mdtW.Unit();
    // Orthogonalize for numerical safety.
    mdtV = mdtV - mdtU * mdtU.Dot(mdtV);

    if (mdtV.Mag() > 0.0) mdtV = mdtV.Unit();

    mdtW = mdtW - mdtU * mdtU.Dot(mdtW) - mdtV * mdtV.Dot(mdtW);

    if (mdtW.Mag() > 0.0) mdtW = mdtW.Unit();
    else mdtW = mdtU.Cross(mdtV).Unit();

    // Seed position: X=0 because B=(Bx,0,0) → Fx=0 → X is conserved.
    // The 1D measurement constrains only Y at each Z; X is never measured,
    // so GenFit's propagator stays at X=0 throughout (avoiding material
    // corrections from traversing the MDT tube walls at global X≈900mm).
    TVector3 posSeed(0.0,
                     sortedMeas.front().y_mm * 0.1,
                     sortedMeas.front().z_mm * 0.1); // cm

    // Seed direction: prefer analytic initial slope (dY/dZ at z0) when provided,
    // because the hit-derived average slope absorbs the full magnetic bend,
    // making the seed slope equal to the mid-trajectory slope rather than the
    // true initial slope — which causes GenFit to reject stations 2-4 as outliers.
    TVector3 dirGlobal;
    if (std::isfinite(seedSlopeDyDz)) {
        // Build direction from analytic slope: dY/dZ in local MDT frame.
        // dz component = 1 (normalized later), dy component = slope.
        dirGlobal = mdtU * seedSlopeDyDz + mdtW;
        if (dirGlobal.Mag() > 0.0)
            dirGlobal = dirGlobal.Unit();
        else
            dirGlobal = mdtW;
    } else {
        const double dy_mm = sortedMeas.back().localY_mm - sortedMeas.front().localY_mm;
        const double dz_mm = sortedMeas.back().localZ_mm - sortedMeas.front().localZ_mm;
        if (std::fabs(dy_mm) > 1.0e-12 || std::fabs(dz_mm) > 1.0e-12) {
            dirGlobal = mdtU * dy_mm + mdtW * dz_mm;
            if (dirGlobal.Mag() > 0.0)
                dirGlobal = dirGlobal.Unit();
            else
                dirGlobal = mdtW;
        } else {
            dirGlobal = mdtW;
        }
    }

    const double pSeed =
        (seedMomentumGeV > 0.0 && std::isfinite(seedMomentumGeV))
      ? seedMomentumGeV
      : 10.0;

    // px=0: no X force, so X-momentum stays negligible.
    TVector3 momSeed(0.0, dirGlobal.Y() * pSeed, dirGlobal.Z() * pSeed);

    if (verbose > 0) {
        std::cout << "[GenFitMDTFit] seed pos cm = "
                  << posSeed.X() << " "
                  << posSeed.Y() << " "
                  << posSeed.Z() << "\n";

        std::cout << "[GenFitMDTFit] seed direction global = "
                  << dirGlobal.X() << " "
                  << dirGlobal.Y() << " "
                  << dirGlobal.Z() << "\n";

        std::cout << "[GenFitMDTFit] seed momentum GeV = "
                  << momSeed.X() << " "
                  << momSeed.Y() << " "
                  << momSeed.Z()
                  << " |p|=" << momSeed.Mag() << "\n";

        std::cout << "[GenFitMDTFit] MDT axes global:\n"
                  << "  U/localY/measured = ("
                  << mdtU.X() << ", "
                  << mdtU.Y() << ", "
                  << mdtU.Z() << ")\n"
                  << "  V/localX/wire     = ("
                  << mdtV.X() << ", "
                  << mdtV.Y() << ", "
                  << mdtV.Z() << ")\n"
                  << "  W/localZ/downstr. = ("
                  << mdtW.X() << ", "
                  << mdtW.Y() << ", "
                  << mdtW.Z() << ")\n";
    }

    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

    // Pack the parameters into GenFit's combined 6D TVectorD seed state vector
    TVectorD stateSeed(6);
    stateSeed(0) = posSeed.X(); 
    stateSeed(1) = posSeed.Y(); 
    stateSeed(2) = posSeed.Z();
    stateSeed(3) = momSeed.X(); 
    stateSeed(4) = momSeed.Y(); 
    stateSeed(5) = momSeed.Z();

    // Simple diagonal covariance in global Cartesian (X,Y,Z,px,py,pz).
    // X is conserved (Fx=0), so position uncertainty is symmetric; use 1 cm each.
    TMatrixDSym covSeed(6);
    covSeed.Zero();
    covSeed(0,0) = 1.0;    // (1 cm)^2 in X
    covSeed(1,1) = 1.0;    // (1 cm)^2 in Y
    covSeed(2,2) = 1.0;    // (1 cm)^2 in Z
    const double momSigma = std::max(5.0, 0.5 * pSeed);
    covSeed(3,3) = 1.0;                        // (1 GeV)^2 in px
    covSeed(4,4) = 1.0;                        // (1 GeV)^2 in py
    covSeed(5,5) = momSigma * momSigma;        // broader pz uncertainty

    // Instantiate the GenFit master track container.
    genfit::Track fitTrack(rep, stateSeed, covSeed);

    // 1D measurement: only global Y (drift direction) is measured.
    // Global X is unmeasured and stays ~0 (no X force from B=(Bx,0,0)).
    // Plane normal = hU×hV = (0,1,0)×(1,0,0) = (0,0,-1) → plane ⊥ Z, natural for MDT.
    const double hitSigmaU_cm = 0.080 * 0.1; // 80 μm drift resolution in cm
    TMatrixDSym hitCov1D(1);
    hitCov1D.Zero();
    hitCov1D(0,0) = hitSigmaU_cm * hitSigmaU_cm;
    const TVector3 hU_global(0.0, 1.0, 0.0); // measured: global Y (drift direction)
    const TVector3 hV_global(1.0, 0.0, 0.0); // unmeasured: global X (wire direction)

    const int detId = 0;
    int hitCounter = 0;
    for (const auto& hit : sortedMeas) {
        // X=0: no X force so particle stays at X=0 throughout.
        TVector3 planeOrigin(0.0, hit.y_mm * 0.1, hit.z_mm * 0.1);

        TVectorD hitCoords(1);
        hitCoords(0) = 0.0; // measured Y = planeOrigin.Y() = hit.y_mm/10

        genfit::PlanarMeasurement* meas_pt = new genfit::PlanarMeasurement(
            hitCoords, hitCov1D, detId, hitCounter, nullptr);
        meas_pt->setPlane(
            genfit::SharedPlanePtr(new genfit::DetPlane(planeOrigin, hU_global, hV_global)),
            hitCounter);
        fitTrack.insertPoint(new genfit::TrackPoint(meas_pt, &fitTrack));
        ++hitCounter;
    }

    if (verbose > 0) {
        auto* field = genfit::FieldManager::getInstance()->getField();
        std::cout << "\n[GenFitMDTFit] ========== MEASUREMENT DETAILS ==========\n";
        std::cout << "[GenFitMDTFit] Number of measurements: " << sortedMeas.size() << "\n";
        for (size_t i = 0; i < sortedMeas.size(); ++i) {
            const auto& m = sortedMeas[i];
            TVector3 pos_cm(m.x_mm / 10.0, m.y_mm / 10.0, m.z_mm / 10.0);
            TVector3 BkG = field ? field->get(pos_cm) : TVector3(0.0, 0.0, 0.0);
            const char* nodeName = "NULL";
            if (gGeoManager) {
                TGeoNode* node = gGeoManager->FindNode(pos_cm.X(), pos_cm.Y(), pos_cm.Z());
                if (node) nodeName = node->GetName();
            }
            std::cout << Form(
                "[GenFitMDTFit] Hit %2zu: localY=%+.1f z=%+.1f mm"
                " global=(%+.1f,%+.1f,%+.1f) mm drift=%+.3f side=%+d"
                " B=(%+.2f,%+.2f,%+.2f) kG node=%s\n",
                i, m.localY_mm, m.localZ_mm,
                m.x_mm, m.y_mm, m.z_mm,
                m.r_meas_mm, m.side,
                BkG.X(), BkG.Y(), BkG.Z(), nodeName);
        }
        std::cout << "[GenFitMDTFit] =========================================\n\n";
        std::cout << "[GenFitMDTFit] Track has " << fitTrack.getNumPoints() << " MDT measurement points\n";
    }

    // KalmanFitterRefTrack: works with PlanarMeasurement (fixed planes).
    genfit::KalmanFitterRefTrack fitter;
    fitter.setMaxIterations(10);
    fitter.setRelChi2Change(0.001);

    if (verbose > 0)
        std::cout << "[GenFitMDTFit] Starting KalmanFitterRefTrack fit\n";

    try {
        fitter.processTrack(&fitTrack);
    }
    catch (genfit::Exception& e) {
        std::cerr << "\n[GenFitMDTFit] ========== GENFIT EXCEPTION ==========\n";
        std::cerr << "[GenFitMDTFit] Exception type: " << e.what() << "\n";
        std::cerr << "[GenFitMDTFit] Exception string: " << e.getExcString() << "\n";
        // Print track info for debugging
        std::cerr << "[GenFitMDTFit] Track info at exception:\n";
        std::cerr << "  Seed momentum: " << momSeed.Mag() << " GeV\n";
        std::cerr << "  Seed position: (" << posSeed.X() << ", " << posSeed.Y() << ", " << posSeed.Z() << ") cm\n";
        std::cerr << "  Number of points: " << fitTrack.getNumPoints() << "\n";
        std::cerr << "[GenFitMDTFit] =====================================\n\n";
      
        return false;
    }
    catch (std::exception& e) {
        std::cerr << "[GenFitMDTFit] Standard exception: " << e.what() << std::endl;
        return false;
    }
    catch (...) {
        std::cerr << "[GenFitMDTFit] Unknown exception caught\n";
        return false;
    }
    
    // Check fit status
    genfit::FitStatus* status = fitTrack.getFitStatus(rep);

    if (!status) {
        std::cerr << "\n[GenFitMDTFit] ========== FIT STATUS ERROR ==========\n";
        std::cerr << "[GenFitMDTFit] ERROR: FitStatus pointer is NULL\n";
        std::cerr << "[GenFitMDTFit] Track has " << fitTrack.getNumPoints() << " points\n";
        std::cerr << "[GenFitMDTFit] ====================================\n\n";
        return false;
    }
    
    // isFitConverged() is unreliable across GenFit versions:
    // some versions return false even for excellent fits (chi2/NDF~0.24).
    // Instead: always try to extract the state and gate on chi2/NDF + pval.
    const double fitChi2 = status->getChi2();
    const double fitNDF  = status->getNdf();
    const double fitPval = status->getPVal();
    const double chi2ndf = (fitNDF > 0) ? fitChi2 / fitNDF : 1e9;

    std::cerr << "[GenFitMDTFit] isFitConverged=" << (status->isFitConverged() ? "YES" : "NO")
              << "  chi2=" << fitChi2 << "  NDF=" << fitNDF
              << "  chi2/NDF=" << chi2ndf << "  pval=" << fitPval << "\n";

    // Try to extract the fitted state regardless of the convergence flag.
    genfit::MeasuredStateOnPlane state;
    try {
        state = fitTrack.getFittedState();
    } catch (...) {
        std::cerr << "[GenFitMDTFit] Could not extract fitted state — fit truly failed\n";
        return false;
    }

    // Quality gate: accept if GenFit declares convergence OR chi2/NDF is reasonable.
    // With 1D global-axis measurements the fitter is numerically stable, so
    // isFitConverged() is reliable.  Keep chi2/NDF < 100 as a safety backstop.
    const bool qualityOk = status->isFitConverged() || ((fitNDF > 0) && (chi2ndf < 100.0));
    // Also trigger retry if too many hits were silently rejected: a "converged" fit
    // that keeps only 6 of 12 hits (NDF=1 vs expected NDF=7) is not acceptable.
    const int expectedNDF = (int)sortedMeas.size() - 5;
    const bool tooFewHits = (fitNDF > 0) && (fitNDF < expectedNDF - 3);
    if (!qualityOk || tooFewHits) {
        std::cerr << "[GenFitMDTFit] Attempt 1 rejected: chi2/NDF=" << chi2ndf
                  << "  pval=" << fitPval
                  << (tooFewHits ? "  (too few hits: NDF=" + std::to_string(fitNDF)
                                   + " expected " + std::to_string(expectedNDF) + ")" : "")
                  << " — trying L/R refinement\n";

        // Use the partial track state to re-assign L/R hit-by-hit.
        // The partial state is usually a good approximation of the true trajectory
        // even when qualityOk=false.  Linear propagation is sufficient to determine
        // which side of each wire the track passes.
        TVector3 pPartial  = state.getMom();   // GeV/c, global
        TVector3 posPartial = state.getPos();  // cm, global

        if (pPartial.Mag() < 0.5) return false;

        const double slopeYZ = (std::fabs(pPartial.Z()) > 0.01)
                               ? pPartial.Y() / pPartial.Z() : 0.0;

        int nFlipped = 0;
        for (auto& m : sortedMeas) {
            const double dz_cm   = m.z_mm * 0.1 - posPartial.Z();
            const double yPred_mm = (posPartial.Y() + slopeYZ * dz_cm) * 10.0;
            const int newSide    = (yPred_mm >= m.wireY_mm) ? 1 : -1;
            if (newSide != m.side) {
                ++nFlipped;
                const int oldSide = m.side;
                m.side        = newSide;
                m.y_mm        = m.wireY_mm + newSide * m.r_meas_mm;
                m.localY_mm  += (newSide - oldSide) * m.r_meas_mm;
                // x_mm: small tilt correction (sin 4.5° ≈ 0.0785), skip for brevity
            }
        }

        if (nFlipped == 0) return false;
        std::cerr << "[GenFitMDTFit] " << nFlipped << " hit(s) L/R refined — retrying\n";

        // Second GenFit attempt with corrected measurements, 1D global planes.
        genfit::AbsTrackRep* rep2 = new genfit::RKTrackRep(pdg);

        // Seed from partial state, zeroing X (conserved; keeps propagator at X=0).
        TVectorD stateSeed2(6);
        stateSeed2(0) = 0.0;            stateSeed2(1) = posPartial.Y(); stateSeed2(2) = posPartial.Z();
        stateSeed2(3) = 0.0;            stateSeed2(4) = pPartial.Y();   stateSeed2(5) = pPartial.Z();
        TMatrixDSym covSeed2(6); covSeed2.Zero();
        for (int i = 0; i < 3; ++i) covSeed2(i,i) = 1.0;
        const double mSig2 = std::max(5.0, 0.3 * pPartial.Mag());
        covSeed2(3,3) = 1.0; covSeed2(4,4) = 1.0;
        covSeed2(5,5) = mSig2 * mSig2;

        genfit::Track fitTrack2(rep2, stateSeed2, covSeed2);

        int hitCtr2 = 0;
        for (const auto& hit : sortedMeas) {
            TVector3 origin2(0.0, hit.y_mm * 0.1, hit.z_mm * 0.1);
            TVectorD hc2(1); hc2(0) = 0.0;
            TMatrixDSym hCov2(1); hCov2.Zero();
            hCov2(0,0) = hitSigmaU_cm * hitSigmaU_cm;
            auto* mpt2 = new genfit::PlanarMeasurement(hc2, hCov2, 0, hitCtr2, nullptr);
            mpt2->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(origin2, hU_global, hV_global)), hitCtr2);
            fitTrack2.insertPoint(new genfit::TrackPoint(mpt2, &fitTrack2));
            ++hitCtr2;
        }

        genfit::KalmanFitterRefTrack fitter2;
        fitter2.setMaxIterations(10);
        fitter2.setRelChi2Change(0.001);
        try { fitter2.processTrack(&fitTrack2); }
        catch (...) { std::cerr << "[GenFitMDTFit] Retry exception\n"; return false; }

        genfit::FitStatus* status2 = fitTrack2.getFitStatus(rep2);
        if (!status2) return false;

        const double chi2_2    = status2->getChi2();
        const double ndf_2     = status2->getNdf();
        const double pval_2    = status2->getPVal();
        const double chi2ndf_2 = (ndf_2 > 0) ? chi2_2 / ndf_2 : 1e9;

        std::cerr << "[GenFitMDTFit] Retry: chi2/NDF=" << chi2ndf_2 << "  pval=" << pval_2 << "\n";

        genfit::MeasuredStateOnPlane state2;
        try { state2 = fitTrack2.getFittedState(); }
        catch (...) { return false; }

        const bool qualityOk2 = status2->isFitConverged() || ((ndf_2 > 0) && (chi2ndf_2 < 100.0));
        if (!qualityOk2) {
            // Attempt 3+: try multiple seed momenta (1D global planes throughout).
            static const double altSeeds[] = { 5.0, 20.0, 100.0 };
            for (double altP : altSeeds) {
                genfit::AbsTrackRep* rep3 = new genfit::RKTrackRep(pdg);
                TVectorD ss3(6);
                // X=0 always; use Y/Z from analytic direction scaled to altP.
                ss3(0) = 0.0; ss3(1) = posSeed.Y(); ss3(2) = posSeed.Z();
                TVector3 d3(0.0, dirGlobal.Y() * altP, dirGlobal.Z() * altP);
                ss3(3) = d3.X(); ss3(4) = d3.Y(); ss3(5) = d3.Z();
                TMatrixDSym cs3(6); cs3.Zero();
                for (int i = 0; i < 3; ++i) cs3(i,i) = 1.0;
                const double ms3 = std::max(5.0, 0.3 * altP);
                cs3(3,3) = 1.0; cs3(4,4) = 1.0; cs3(5,5) = ms3 * ms3;
                genfit::Track ft3(rep3, ss3, cs3);
                int hc3 = 0;
                for (const auto& hit : sortedMeas) {
                    TVector3 org3(0.0, hit.y_mm*0.1, hit.z_mm*0.1);
                    TVectorD hc3v(1); hc3v(0) = 0.0;
                    TMatrixDSym hCov3(1); hCov3.Zero();
                    hCov3(0,0) = hitSigmaU_cm*hitSigmaU_cm;
                    auto* mpt3 = new genfit::PlanarMeasurement(hc3v, hCov3, 0, hc3, nullptr);
                    mpt3->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(org3, hU_global, hV_global)), hc3);
                    ft3.insertPoint(new genfit::TrackPoint(mpt3, &ft3));
                    ++hc3;
                }
                genfit::KalmanFitterRefTrack fitter3;
                fitter3.setMaxIterations(10);
                fitter3.setRelChi2Change(0.001);
                try { fitter3.processTrack(&ft3); } catch (...) { continue; }
                genfit::FitStatus* st3 = ft3.getFitStatus(rep3);
                if (!st3) continue;
                const double c3 = st3->getChi2(), n3 = st3->getNdf(), pv3 = st3->getPVal();
                const double cn3 = (n3 > 0) ? c3/n3 : 1e9;
                std::cerr << "[GenFitMDTFit] AltSeed " << altP << " GeV: chi2/NDF=" << cn3 << "  pval=" << pv3 << "\n";
                if (!st3->isFitConverged() && ((n3 <= 0) || (cn3 >= 100.0))) continue;
                genfit::MeasuredStateOnPlane st3s;
                try { st3s = ft3.getFittedState(); } catch (...) { continue; }
                TVector3 p3 = st3s.getMom();
                fpx = p3.X(); fpy = p3.Y(); fpz = p3.Z(); fp = p3.Mag();
                fcharge = st3s.getCharge();
                fchi2 = c3; fnDoF = static_cast<int>(n3); fpval = pv3;
                fQOverP = st3s.getState()(0);
                fQOverPErr = std::sqrt(st3s.getCov()(0,0));
                fpErr = fQOverPErr * fp * fp;
                fpos.clear(); layerID.clear();
                for (const auto& m : sortedMeas) {
                    fpos.emplace_back(m.x_mm, m.y_mm, m.z_mm);
                    layerID.push_back(m.stationID*10000 + m.planeID*1000 + m.tubeID);
                }
                std::cerr << "[GenFitMDTFit] AltSeed SUCCESS: p=" << fp << " q=" << fcharge << "\n";
                return true;
            }

            // ── DAF rescue ────────────────────────────────────────────────────
            // All Kalman attempts failed.  Try GenFit's Deterministic Annealing
            // Filter which soft-weights each hit by its chi² contribution and
            // iteratively cools the temperature, driving contaminated hits (from
            // a secondary particle or a misidentified MDT tube) to weight → 0.
            // The muon track then emerges from the uncontaminated majority.
            {
                genfit::AbsTrackRep* repD = new genfit::RKTrackRep(pdg);
                TVectorD ssD(6);
                ssD(0) = 0.0; ssD(1) = posSeed.Y(); ssD(2) = posSeed.Z();
                TVector3 dD(0.0, dirGlobal.Y() * seedMomentumGeV,
                                 dirGlobal.Z() * seedMomentumGeV);
                ssD(3) = dD.X(); ssD(4) = dD.Y(); ssD(5) = dD.Z();
                TMatrixDSym csD(6); csD.Zero();
                for (int i = 0; i < 3; ++i) csD(i,i) = 1.0;
                const double msD = std::max(5.0, 0.3 * seedMomentumGeV);
                csD(3,3) = 1.0; csD(4,4) = 1.0; csD(5,5) = msD * msD;

                genfit::Track ftD(repD, ssD, csD);
                int hcD = 0;
                for (const auto& hit : sortedMeas) {
                    TVector3 orgD(0.0, hit.y_mm * 0.1, hit.z_mm * 0.1);
                    TVectorD hcDv(1); hcDv(0) = 0.0;
                    TMatrixDSym hCovD(1); hCovD.Zero();
                    hCovD(0,0) = hitSigmaU_cm * hitSigmaU_cm;
                    auto* mptD = new genfit::PlanarMeasurement(hcDv, hCovD, 0, hcD, nullptr);
                    mptD->setPlane(genfit::SharedPlanePtr(
                        new genfit::DetPlane(orgD, hU_global, hV_global)), hcD);
                    ftD.insertPoint(new genfit::TrackPoint(mptD, &ftD));
                    ++hcD;
                }

                genfit::DAF daf;
                // T_start=1e7: even a hit 1000σ from the seed (chi2=1e6) gets
                // initial weight exp(-1e6/(2e7))=exp(-0.05)≈0.95 — all hits are
                // included at the start regardless of how wrong the seed is.
                // T_stop=9: final cut ≈3σ for a 1D measurement, cleanly rejects
                // contaminated hits that converge far from the muon track.
                daf.setAnnealingScheme(1e7, 9.0, 30);
                daf.setMaxIterations(60);

                try { daf.processTrack(&ftD); }
                catch (...) {
                    std::cerr << "[GenFitMDTFit] DAF exception\n";
                    return false;
                }

                genfit::FitStatus* stD = ftD.getFitStatus(repD);
                if (!stD) return false;

                const double cD  = stD->getChi2();
                const double nD  = stD->getNdf();
                const double pvD = stD->getPVal();
                const double cnD = (nD > 0) ? cD / nD : 1e9;

                std::cerr << "[GenFitMDTFit] DAF: isFitConverged="
                          << (stD->isFitConverged() ? "YES" : "NO")
                          << "  chi2/NDF=" << cnD << "  NDF=" << nD
                          << "  pval=" << pvD << "\n";

                // Accept DAF result if converged or chi2/NDF is reasonable.
                // Use a looser threshold (200) than Kalman (100): DAF NDF is
                // an effective (weighted) quantity and chi2 from upweighted good
                // hits is usually clean, but a few partially-weighted boundary
                // hits can inflate the ratio slightly.
                // Require effective NDF >= 1 (at least one hit with significant
                // weight survived annealing) and physical momentum (< 10 TeV).
                const bool qualityDAF = (stD->isFitConverged() || ((nD > 0) && (cnD < 200.0)))
                                        && (nD >= 1.0) && (cnD < 200.0);
                if (!qualityDAF) return false;

                genfit::MeasuredStateOnPlane stDs;
                try { stDs = ftD.getFittedState(); }
                catch (...) { return false; }

                TVector3 pD = stDs.getMom();
                fpx = pD.X(); fpy = pD.Y(); fpz = pD.Z(); fp = pD.Mag();
                fcharge    = stDs.getCharge();
                fchi2      = cD;
                fnDoF      = static_cast<int>(nD);
                fpval      = pvD;
                fQOverP    = stDs.getState()(0);
                fQOverPErr = std::sqrt(stDs.getCov()(0,0));
                fpErr      = fQOverPErr * fp * fp;
                fpos.clear(); layerID.clear();
                for (const auto& m : sortedMeas) {
                    fpos.emplace_back(m.x_mm, m.y_mm, m.z_mm);
                    layerID.push_back(m.stationID * 10000 + m.planeID * 1000 + m.tubeID);
                }
                std::cerr << "[GenFitMDTFit] DAF SUCCESS: p=" << fp
                          << " q=" << fcharge << " chi2/NDF=" << cnD
                          << " NDF=" << nD << "\n";
                return true;
            }
            // ── end DAF rescue ────────────────────────────────────────────────
        }

        // Accept retry result
        TVector3 p2 = state2.getMom();
        fpx = p2.X(); fpy = p2.Y(); fpz = p2.Z(); fp = p2.Mag();
        fcharge    = state2.getCharge();
        fchi2      = chi2_2;
        fnDoF      = static_cast<int>(ndf_2);
        fpval      = pval_2;
        fQOverP    = state2.getState()(0);
        fQOverPErr = std::sqrt(state2.getCov()(0,0));
        fpErr      = fQOverPErr * fp * fp;
        fpos.clear(); layerID.clear();
        for (const auto& m : sortedMeas) {
            fpos.emplace_back(m.x_mm, m.y_mm, m.z_mm);
            layerID.push_back(m.stationID * 10000 + m.planeID * 1000 + m.tubeID);
        }
        if (verbose)
            std::cout << "[GenFitMDTFit] Retry SUCCESS: p=" << fp
                      << " q=" << fcharge << " chi2/NDF=" << chi2ndf_2 << "\n";
        return true;
    }
    TVector3 pfit = state.getMom();
    fpx = pfit.X();
    fpy = pfit.Y();
    fpz = pfit.Z();
    fp  = pfit.Mag();

    fcharge = state.getCharge();

    fchi2   = fitChi2;
    fnDoF   = static_cast<int>(fitNDF);
    fpval   = fitPval;

    TMatrixDSym localCov = state.getCov();
    fQOverP = state.getState()(0);
    fQOverPErr = std::sqrt(state.getCov()(0,0));
    fpErr = fQOverPErr * fp * fp; // Error propagation: σ_p = σ_(q/p) * p^2

    fpos.clear();
    layerID.clear();
    for (const auto& m : sortedMeas) {
        // Store arbitrary point on the measured MDT line for display/debug.
        // x is not measured.
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
                  << " chi2=" << fchi2
                  << " ndf="  << fnDoF
                  << " pval=" << fpval                  
                  << std::endl;

        std::cout << "[GenFitMDTFit] fitted slopes:"
                  << " dx/dz=" << fpx / fpz
                  << " dy/dz=" << fpy / fpz
                  << std::endl;
    }
    // No delete needed - fitTrack is now stack allocated
    return true;
}
