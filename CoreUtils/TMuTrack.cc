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

// ============================================================================
/// GenFitMDTFit: Clean, optimized MDT fitting with correct charge handling
/// ============================================================================
///
/// Design goals:
///   1. Single-hypothesis fitting (fits PDG code independently)
///   2. Clear charge sign convention (documented, consistent)
///   3. Parametrized quality gates (easy to tune for different configurations)
///   4. Separated concerns (L/R assignment → fitting → result extraction)
///   5. Proper error handling and logging
///
/// Physics:
///   - Magnetic field: Bx = -1.5 T (constant in MDT region)
///   - Charge convention: PDG 13 = µ⁻ (charge -1), PDG -13 = µ⁺ (charge +1)
///   - Bending: F = q(v × B), so µ⁻ curves +Y, µ⁺ curves -Y
///   - qSignAlt: +1 for µ⁻, -1 for µ⁺ (used in L/R prediction formula)
///
/// Stages:
///   Stage 1: Simple Kalman fit with Phase-3.5 L/R assignment
///   Stage 2: L/R refinement using partial state (if Stage 1 poor)
///   Stage 3: Curvature-corrected L/R at multiple momentum guesses
///   Stage 4: DAF rescue path (for contaminated hits)
///
/// ============================================================================
// ============================================================================
// Quality Gates & Parameters (tune here for different configurations)
// ============================================================================

struct FitQualityGates {
    // Momentum validity
    double pxSaneMaxRatio = 0.5;           // max |px| / |pz| (drift dominated)
    double pSaneMinMom = 1.0;              // floor momentum [GeV]
    double pSaneMaxMultiplier = 50.0;      // max p / seed_p ratio
    double pSaneMinMultiplier = 0.05;      // min p / seed_p ratio
    double maxRelMomErr = 1.0;             // reject if sigma_p/p >= this (degenerate fit,
                                            // 2026-08-16: covariance-based catch for the
                                            // "sloppy direction" failure mode -- a near-straight
                                            // high-momentum track can report chi2/NDF~0.2 while
                                            // being essentially unconstrained in curvature; chi2
                                            // alone cannot see this, but the fit's own covariance
                                            // (fpErr) does. Empirically on 5210-event validation,
                                            // all 36 known q*p>2000 GeV/c outliers had fpErr/p in
                                            // [1.04, 11.5], while well-behaved tracks cluster below
                                            // ~0.3-0.5 with a clean cliff after ~0.7.

    // Chi2/NDF thresholds
    double chi2NdfStage1Gate = 100.0;      // acceptance for Stage 1 Kalman
    double chi2NdfStage2Gate = 100.0;      // acceptance for Stage 2 refinement
    double chi2NdfDAFGate = 200.0;         // looser for DAF (weighted hits)

    // Hit count validation
    double expectedNdfMargin = 3;          // allow NDF < (N-5) by this amount

    // L/R refinement parameters
    double minPartialMom = 0.5;            // min momentum to use partial state

    // Curvature correction seeds
    std::vector<double> altMomentaGeV = {5.0, 20.0, 100.0};
    double B_field_magnitude = 1.5;        // Tesla

    // DAF annealing
    double daf_T_start = 1e7;
    double daf_T_stop = 9.0;
    int daf_n_anneal = 30;
    int daf_max_iter = 60;
};

// ============================================================================
// Helper: Sanity checks for momentum
// ============================================================================

static bool pxSane(const TVector3& mom, double maxRatio = 0.5) {
    return std::fabs(mom.X()) < maxRatio * std::fabs(mom.Z()) + 1.0;
}

static bool pSane(const TVector3& mom, double seedMom, double minMult = 0.05, double maxMult = 50.0) {
    double seedRef = std::max(1.0, seedMom);
    double pmag = mom.Mag();
    return pmag > std::max(1.0, minMult * seedRef) && pmag < maxMult * seedRef;
}

// Relative momentum uncertainty from the fit's own covariance:
// sigma_p/p = sigma_(q/p) * p  (since |q/p| = 1/p). Unlike pSane(), which only
// checks the central value against the seed, this checks whether the fit
// actually constrains momentum at all -- catching degenerate/"sloppy
// direction" solutions that chi2/NDF alone cannot see.
static bool pRelErrSane(double qOverPErr, double pmag, double maxRelErr = 1.0) {
    if (pmag <= 0.0) return false;
    double relErr = qOverPErr * pmag;
    return relErr < maxRelErr;
}

// ============================================================================
// Main Fit Function
// ============================================================================

bool TMuTrack::GenFitMDTFit(
    const std::vector<MDTMeas>& meas,
    int pdg,
    double seedMomentumGeV,
    int verbose,
    double seedSlopeDyDz)
{
    // Charge sign convention
    const int expectedCharge = (pdg == 13) ? -1 : +1;
    const bool isMuPlus = (pdg == -13);
    const double qSignAlt = isMuPlus ? -1.0 : +1.0;  // for L/R formula

    if (verbose > 0) {
        std::cout << "[GenFitMDTFit] Starting fit for PDG=" << pdg
                  << " (charge=" << expectedCharge << "), seed p=" << seedMomentumGeV
                  << " GeV\n";
    }

    // Validate input
    if (meas.size() < 5) {
        if (verbose > 0) std::cerr << "[GenFitMDTFit] Too few measurements: " << meas.size() << "\n";
        return false;
    }

    // Initialize quality gates
    FitQualityGates gates;

    // ========================================================================
    // STAGE 1: Simple Kalman Fit with Phase-3.5 L/R Assignment
    // ========================================================================

    // Sort measurements by Z
    std::vector<MDTMeas> sortedMeas = meas;
    std::sort(sortedMeas.begin(), sortedMeas.end(),
              [](const MDTMeas& a, const MDTMeas& b) { return a.localZ_mm < b.localZ_mm; });

    // Extract MDT coordinate system from first hit
    TVector3 mdtU(sortedMeas.front().uX, sortedMeas.front().uY, sortedMeas.front().uZ);
    TVector3 mdtV(sortedMeas.front().vX, sortedMeas.front().vY, sortedMeas.front().vZ);
    TVector3 mdtW(sortedMeas.front().wX, sortedMeas.front().wY, sortedMeas.front().wZ);

    // Normalize and orthogonalize
    if (mdtU.Mag() > 0) mdtU = mdtU.Unit();
    if (mdtV.Mag() > 0) mdtV = mdtV.Unit();
    mdtV = mdtV - mdtU * mdtU.Dot(mdtV);
    if (mdtV.Mag() > 0) mdtV = mdtV.Unit();
    mdtW = mdtW - mdtU * mdtU.Dot(mdtW) - mdtV * mdtV.Dot(mdtW);
    if (mdtW.Mag() > 0) mdtW = mdtW.Unit();
    else mdtW = mdtU.Cross(mdtV).Unit();

    // Seed position
    TVector3 posSeed(sortedMeas.front().x_mm * 0.1,
                     sortedMeas.front().y_mm * 0.1,
                     sortedMeas.front().z_mm * 0.1);

    // Seed direction: prefer analytic slope if provided
    TVector3 dirGlobal;
    if (std::isfinite(seedSlopeDyDz)) {
        dirGlobal = mdtU * seedSlopeDyDz + mdtW;
        if (dirGlobal.Mag() > 0) dirGlobal = dirGlobal.Unit();
        else dirGlobal = mdtW;
    } else {
        // Hit-derived average slope
        const double dy_mm = sortedMeas.back().localY_mm - sortedMeas.front().localY_mm;
        const double dz_mm = sortedMeas.back().localZ_mm - sortedMeas.front().localZ_mm;
        if (std::fabs(dy_mm) > 1.0e-12 || std::fabs(dz_mm) > 1.0e-12) {
            dirGlobal = mdtU * dy_mm + mdtW * dz_mm;
            if (dirGlobal.Mag() > 0) dirGlobal = dirGlobal.Unit();
            else dirGlobal = mdtW;
        } else {
            dirGlobal = mdtW;
        }
    }

    // Seed momentum
    const double pSeed = (seedMomentumGeV > 0.0 && std::isfinite(seedMomentumGeV))
                       ? seedMomentumGeV : 10.0;
    TVector3 momSeed = dirGlobal * pSeed;

    if (verbose > 2) {
        std::cout << "[GenFitMDTFit] Seed: pos=(" << posSeed.X() << "," << posSeed.Y()
                  << "," << posSeed.Z() << ") cm, dir=(" << dirGlobal.X() << "," << dirGlobal.Y()
                  << "," << dirGlobal.Z() << "), p=" << momSeed.Mag() << " GeV\n";
    }

    // Create GenFit track
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

    TVectorD stateSeed(6);
    stateSeed(0) = posSeed.X();
    stateSeed(1) = posSeed.Y();
    stateSeed(2) = posSeed.Z();
    stateSeed(3) = momSeed.X();
    stateSeed(4) = momSeed.Y();
    stateSeed(5) = momSeed.Z();

    TMatrixDSym covSeed(6);
    covSeed.Zero();
    covSeed(0, 0) = covSeed(1, 1) = covSeed(2, 2) = 1.0;  // 1 cm position uncertainty
    covSeed(3, 3) = covSeed(4, 4) = 1.0;                  // 1 GeV momentum uncertainty
    const double momSigma = std::max(5.0, 0.5 * pSeed);
    covSeed(5, 5) = momSigma * momSigma;

    genfit::Track fitTrack(rep, stateSeed, covSeed);

    // Add measurements (1D: only Y is measured, X and Z define the plane)
    const double hitSigmaU_cm = 0.080 * 0.1;  // 80 μm drift resolution
    TMatrixDSym hitCov1D(1);
    hitCov1D(0, 0) = hitSigmaU_cm * hitSigmaU_cm;

    for (size_t i = 0; i < sortedMeas.size(); ++i) {
        const auto& hit = sortedMeas[i];
        TVector3 planeOrigin(hit.x_mm * 0.1, hit.y_mm * 0.1, hit.z_mm * 0.1);
        TVectorD hitCoords(1);
        hitCoords(0) = 0.0;  // measured coordinate on plane

        auto* meas_pt = new genfit::PlanarMeasurement(hitCoords, hitCov1D, 0, i, nullptr);
        meas_pt->setPlane(genfit::SharedPlanePtr(
            new genfit::DetPlane(planeOrigin, mdtU, mdtV)), i);
        fitTrack.insertPoint(new genfit::TrackPoint(meas_pt, &fitTrack));
    }

    // Fit with Kalman
    genfit::KalmanFitterRefTrack fitter;
    fitter.setMaxIterations(10);
    fitter.setRelChi2Change(0.001);

    try {
        fitter.processTrack(&fitTrack);
    } catch (genfit::Exception& e) {
        if (verbose > 0) std::cerr << "[GenFitMDTFit] Kalman exception: " << e.what() << "\n";
        return false;
    }

    // Extract results
    genfit::FitStatus* status = fitTrack.getFitStatus(rep);
    if (!status) {
        if (verbose > 0) std::cerr << "[GenFitMDTFit] No fit status\n";
        return false;
    }

    double chi2_1 = status->getChi2();
    double ndf_1 = status->getNdf();
    double pval_1 = status->getPVal();

    genfit::MeasuredStateOnPlane state_1;
    try {
        state_1 = fitTrack.getFittedState();
    } catch (...) {
        if (verbose > 0) std::cerr << "[GenFitMDTFit] Could not extract fitted state\n";
        return false;
    }

    double chi2ndf_1 = (ndf_1 > 0) ? chi2_1 / ndf_1 : 1e9;
    // AND, never OR, between "GenFit says converged" and "chi2/NDF is sane".
    // isFitConverged() reflects numerical iteration stability, not physical
    // sanity: a degenerate near-straight/near-infinite-momentum solution can
    // report converged=true while chi2/NDF is in the thousands. An OR here
    // lets that garbage through untouched -- this is precisely the bug that
    // put a p=3.3e11 GeV, chi2/NDF=7505 track into production data (event
    // 398, track 1, 2026-08-15 investigation). Require both, always.
    bool qualityOk_1 = status->isFitConverged()
                    && (chi2ndf_1 < gates.chi2NdfStage1Gate)
                    && pxSane(state_1.getMom())
                    && pSane(state_1.getMom(), pSeed)
                    && pRelErrSane(std::sqrt(state_1.getCov()(0,0)), state_1.getMom().Mag(), gates.maxRelMomErr);

    if (verbose > 1) {
        std::cout << "[GenFitMDTFit] Stage 1 Kalman: chi2/ndf=" << chi2ndf_1
                  << " pval=" << pval_1 << " quality=" << (qualityOk_1 ? "OK" : "POOR") << "\n";
    }

    // If Stage 1 successful, use results
    if (qualityOk_1) {
        TVector3 p = state_1.getMom();
        fpx = p.X();
        fpy = p.Y();
        fpz = p.Z();
        fp = p.Mag();
        fcharge = state_1.getCharge();
        fchi2 = chi2_1;
        fnDoF = static_cast<int>(ndf_1);
        fpval = pval_1;
        fQOverP = state_1.getState()(0);
        fQOverPErr = std::sqrt(state_1.getCov()(0,0));
        fpErr = fQOverPErr * fp * fp; // Error propagation: sigma_p = sigma_(q/p) * p^2

        fpos.clear();
        layerID.clear();
        for (const auto& m : sortedMeas) {
            fpos.emplace_back(m.x_mm, m.y_mm, m.z_mm);
            layerID.push_back(m.stationID * 10000 + m.planeID * 1000 + m.tubeID);
        }

        if (verbose > 0) {
            std::cout << "[GenFitMDTFit] Stage 1 SUCCESS: p=" << fp << " GeV, q=" << fcharge
                      << " chi2/ndf=" << chi2ndf_1 << "\n";
        }
        return true;
    }

    // ========================================================================
    // STAGE 2: L/R Refinement (if Stage 1 poor)
    // ========================================================================

    if (verbose > 1) std::cout << "[GenFitMDTFit] Attempting Stage 2: L/R refinement\n";

    TVector3 pPartial = state_1.getMom();
    TVector3 posPartial = state_1.getPos();

    if (pPartial.Mag() >= gates.minPartialMom) {
        const double slopeYZ = (std::fabs(pPartial.Z()) > 0.01)
                             ? pPartial.Y() / pPartial.Z() : 0.0;

        int nFlipped = 0;
        for (auto& m : sortedMeas) {
            const double dz_cm = m.z_mm * 0.1 - posPartial.Z();
            const double yPred_mm = (posPartial.Y() + slopeYZ * dz_cm) * 10.0;
            const int newSide = (yPred_mm >= m.wireY_mm) ? 1 : -1;
            if (newSide != m.side) {
                nFlipped++;
                const int oldSide = m.side;
                m.side = newSide;
                m.y_mm = m.wireY_mm + newSide * m.r_meas_mm;
                m.localY_mm += (newSide - oldSide) * m.r_meas_mm;
            }
        }

        if (nFlipped > 0) {
            if (verbose > 1) std::cout << "[GenFitMDTFit] " << nFlipped << " hits flipped, retrying fit\n";

            // Re-fit with corrected L/R
            genfit::Track fitTrack2(new genfit::RKTrackRep(pdg), stateSeed, covSeed);
            for (size_t i = 0; i < sortedMeas.size(); ++i) {
                const auto& hit = sortedMeas[i];
                TVector3 planeOrigin(hit.x_mm * 0.1, hit.y_mm * 0.1, hit.z_mm * 0.1);
                TVectorD hitCoords(1);
                hitCoords(0) = 0.0;

                auto* meas_pt = new genfit::PlanarMeasurement(hitCoords, hitCov1D, 0, i, nullptr);
                meas_pt->setPlane(genfit::SharedPlanePtr(
                    new genfit::DetPlane(planeOrigin, mdtU, mdtV)), i);
                fitTrack2.insertPoint(new genfit::TrackPoint(meas_pt, &fitTrack2));
            }

            try {
                fitter.processTrack(&fitTrack2);
            } catch (...) {
                if (verbose > 1) std::cout << "[GenFitMDTFit] Stage 2 fit exception\n";
                // Fall through to Stage 3
                sortedMeas = meas;
                goto stage3;
            }

            genfit::FitStatus* status2 = fitTrack2.getFitStatus();
            if (status2) {
                double chi2_2 = status2->getChi2();
                double ndf_2 = status2->getNdf();
                double chi2ndf_2 = (ndf_2 > 0) ? chi2_2 / ndf_2 : 1e9;

                try {
                    auto state2 = fitTrack2.getFittedState();
                    // AND, not OR -- see Stage 1 comment above.
                    bool qualityOk_2 = status2->isFitConverged()
                                    && (chi2ndf_2 < gates.chi2NdfStage2Gate)
                                    && pxSane(state2.getMom())
                                    && pSane(state2.getMom(), pSeed)
                                    && pRelErrSane(std::sqrt(state2.getCov()(0,0)), state2.getMom().Mag(), gates.maxRelMomErr);

                    if (qualityOk_2) {
                        TVector3 p = state2.getMom();
                        fpx = p.X();
                        fpy = p.Y();
                        fpz = p.Z();
                        fp = p.Mag();
                        fcharge = state2.getCharge();
                        fchi2 = chi2_2;
                        fnDoF = static_cast<int>(ndf_2);
                        fpval = status2->getPVal();
                        fQOverP = state2.getState()(0);
                        fQOverPErr = std::sqrt(state2.getCov()(0,0));
                        fpErr = fQOverPErr * fp * fp; // Error propagation: sigma_p = sigma_(q/p) * p^2

                        fpos.clear();
                        layerID.clear();
                        for (const auto& m : sortedMeas) {
                            fpos.emplace_back(m.x_mm, m.y_mm, m.z_mm);
                            layerID.push_back(m.stationID * 10000 + m.planeID * 1000 + m.tubeID);
                        }

                        if (verbose > 0) {
                            std::cout << "[GenFitMDTFit] Stage 2 SUCCESS: p=" << fp << " GeV, q=" << fcharge
                                      << " chi2/ndf=" << chi2ndf_2 << "\n";
                        }
                        return true;
                    }
                } catch (...) {}
            }
        }
    }

    // ========================================================================
    // STAGE 3: Curvature-Corrected L/R at Multiple Seeds
    // ========================================================================

stage3:
    if (verbose > 1) std::cout << "[GenFitMDTFit] Attempting Stage 3: curvature-corrected L/R\n";

    // Restore Phase-3.5 L/R as baseline for alternative seeds
    sortedMeas = meas;
    std::sort(sortedMeas.begin(), sortedMeas.end(),
              [](const MDTMeas& a, const MDTMeas& b) { return a.localZ_mm < b.localZ_mm; });

    for (double altP : gates.altMomentaGeV) {
        // Straight-line fit through wire centers
        double zMeanW = 0;
        for (const auto& mm : sortedMeas) zMeanW += mm.wireZ_mm;
        zMeanW /= static_cast<double>(sortedMeas.size());

        double sNw = 0, sZw = 0, sZ2w = 0, sYw = 0, sZYw = 0;
        for (const auto& mm : sortedMeas) {
            double z = mm.wireZ_mm - zMeanW;
            sNw += 1;
            sZw += z;
            sZ2w += z * z;
            sYw += mm.wireY_mm;
            sZYw += z * mm.wireY_mm;
        }
        double detW = sNw * sZ2w - sZw * sZw;
        double slpW = (std::fabs(detW) > 1e-9) ? (sNw * sZYw - sZw * sYw) / detW : 0.0;
        double y0W = (std::fabs(detW) > 1e-9) ? (sYw - slpW * sZw) / sNw : 0.0;

        // Helical prediction with curvature
        for (auto& mm : sortedMeas) {
            double zRel = mm.wireZ_mm - zMeanW;
            double K = gates.B_field_magnitude * zRel * zRel * 5.0e-4;
            double yPred = y0W + slpW * zRel + (qSignAlt / altP) * K;
            mm.side = (yPred >= mm.wireY_mm) ? 1 : -1;
            mm.y_mm = mm.wireY_mm + mm.side * mm.r_meas_mm;
        }

        // Fit with this L/R assignment
        genfit::Track fitTrack3(new genfit::RKTrackRep(pdg), stateSeed, covSeed);
        for (size_t i = 0; i < sortedMeas.size(); ++i) {
            const auto& hit = sortedMeas[i];
            TVector3 planeOrigin(hit.x_mm * 0.1, hit.y_mm * 0.1, hit.z_mm * 0.1);
            TVectorD hitCoords(1);
            hitCoords(0) = 0.0;

            auto* meas_pt = new genfit::PlanarMeasurement(hitCoords, hitCov1D, 0, i, nullptr);
            meas_pt->setPlane(genfit::SharedPlanePtr(
                new genfit::DetPlane(planeOrigin, mdtU, mdtV)), i);
            fitTrack3.insertPoint(new genfit::TrackPoint(meas_pt, &fitTrack3));
        }

        try {
            fitter.processTrack(&fitTrack3);
        } catch (...) {
            continue;
        }

        genfit::FitStatus* status3 = fitTrack3.getFitStatus();
        if (!status3) continue;

        double chi2_3 = status3->getChi2();
        double ndf_3 = status3->getNdf();
        double chi2ndf_3 = (ndf_3 > 0) ? chi2_3 / ndf_3 : 1e9;

        try {
            auto state3 = fitTrack3.getFittedState();
            // AND, not OR -- see Stage 1 comment above.
            bool qualityOk_3 = status3->isFitConverged()
                            && (chi2ndf_3 < gates.chi2NdfStage1Gate)
                            && pxSane(state3.getMom())
                            && pSane(state3.getMom(), pSeed)
                            && pRelErrSane(std::sqrt(state3.getCov()(0,0)), state3.getMom().Mag(), gates.maxRelMomErr);

            if (qualityOk_3) {
                TVector3 p = state3.getMom();
                fpx = p.X();
                fpy = p.Y();
                fpz = p.Z();
                fp = p.Mag();
                fcharge = state3.getCharge();
                fchi2 = chi2_3;
                fnDoF = static_cast<int>(ndf_3);
                fpval = status3->getPVal();
                fQOverP = state3.getState()(0);
                fQOverPErr = std::sqrt(state3.getCov()(0,0));
                fpErr = fQOverPErr * fp * fp; // Error propagation: sigma_p = sigma_(q/p) * p^2

                fpos.clear();
                layerID.clear();
                for (const auto& m : sortedMeas) {
                    fpos.emplace_back(m.x_mm, m.y_mm, m.z_mm);
                    layerID.push_back(m.stationID * 10000 + m.planeID * 1000 + m.tubeID);
                }

                if (verbose > 0) {
                    std::cout << "[GenFitMDTFit] Stage 3 (altP=" << altP << " GeV) SUCCESS: p="
                              << fp << " GeV, q=" << fcharge << " chi2/ndf=" << chi2ndf_3 << "\n";
                }
                return true;
            }
        } catch (...) {}
    }

    // ========================================================================
    // STAGE 4: DAF Rescue Path
    // ========================================================================

    if (verbose > 1) std::cout << "[GenFitMDTFit] Attempting Stage 4: DAF rescue\n";

    sortedMeas = meas;

    genfit::Track fitTrackDAF(new genfit::RKTrackRep(pdg), stateSeed, covSeed);
    for (size_t i = 0; i < sortedMeas.size(); ++i) {
        const auto& hit = sortedMeas[i];
        TVector3 planeOrigin(hit.x_mm * 0.1, hit.y_mm * 0.1, hit.z_mm * 0.1);
        TVectorD hitCoords(1);
        hitCoords(0) = 0.0;

        auto* meas_pt = new genfit::PlanarMeasurement(hitCoords, hitCov1D, 0, i, nullptr);
        meas_pt->setPlane(genfit::SharedPlanePtr(
            new genfit::DetPlane(planeOrigin, mdtU, mdtV)), i);
        fitTrackDAF.insertPoint(new genfit::TrackPoint(meas_pt, &fitTrackDAF));
    }

    genfit::DAF daf;
    daf.setAnnealingScheme(gates.daf_T_start, gates.daf_T_stop, gates.daf_n_anneal);
    daf.setMaxIterations(gates.daf_max_iter);

    try {
        daf.processTrack(&fitTrackDAF);
    } catch (...) {
        if (verbose > 0) std::cerr << "[GenFitMDTFit] DAF exception\n";
        return false;
    }

    genfit::FitStatus* statusDAF = fitTrackDAF.getFitStatus();
    if (!statusDAF) {
        if (verbose > 0) std::cerr << "[GenFitMDTFit] DAF: no fit status\n";
        return false;
    }

    double chi2_DAF = statusDAF->getChi2();
    double ndf_DAF = statusDAF->getNdf();
    double pval_DAF = statusDAF->getPVal();
    double chi2ndf_DAF = (ndf_DAF > 0) ? chi2_DAF / ndf_DAF : 1e9;

    // AND, not OR -- see Stage 1 comment above.
    bool qualityOk_DAF = statusDAF->isFitConverged()
                      && (chi2ndf_DAF < gates.chi2NdfDAFGate)
                      && ndf_DAF >= 1.0;

    if (verbose > 1) {
        std::cout << "[GenFitMDTFit] DAF: chi2/ndf=" << chi2ndf_DAF << " quality="
                  << (qualityOk_DAF ? "OK" : "POOR") << "\n";
    }

    if (!qualityOk_DAF) {
        if (verbose > 0) std::cerr << "[GenFitMDTFit] DAF quality failed\n";
        return false;
    }

    try {
        auto stateDAF = fitTrackDAF.getFittedState();
        TVector3 p = stateDAF.getMom();

        if (!pxSane(p) || !pSane(p, pSeed) ||
            !pRelErrSane(std::sqrt(stateDAF.getCov()(0,0)), p.Mag(), gates.maxRelMomErr)) {
            if (verbose > 0) std::cerr << "[GenFitMDTFit] DAF momentum sanity failed\n";
            return false;
        }

        fpx = p.X();
        fpy = p.Y();
        fpz = p.Z();
        fp = p.Mag();
        fcharge = stateDAF.getCharge();
        fchi2 = chi2_DAF;
        fnDoF = static_cast<int>(ndf_DAF);
        fpval = pval_DAF;
        fQOverP = stateDAF.getState()(0);
        fQOverPErr = std::sqrt(stateDAF.getCov()(0,0));
        fpErr = fQOverPErr * fp * fp; // Error propagation: sigma_p = sigma_(q/p) * p^2

        fpos.clear();
        layerID.clear();
        for (const auto& m : sortedMeas) {
            fpos.emplace_back(m.x_mm, m.y_mm, m.z_mm);
            layerID.push_back(m.stationID * 10000 + m.planeID * 1000 + m.tubeID);
        }

        if (verbose > 0) {
            std::cout << "[GenFitMDTFit] DAF SUCCESS: p=" << fp << " GeV, q=" << fcharge
                      << " chi2/ndf=" << chi2ndf_DAF << "\n";
        }
        return true;
    } catch (...) {
        if (verbose > 0) std::cerr << "[GenFitMDTFit] DAF state extraction failed\n";
        return false;
    }
}
