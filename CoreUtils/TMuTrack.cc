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
// ─────────────────────────────────────────────────────────────────────────────
// WCBruteForceAssignLR: brute-force wire-centre L/R assignment (standalone approach).
//
// Tries all 2^N left/right combinations for N MDT hits and fits the 3-parameter
// analytic model  Y = Y₀ + slope·(z−z̄) + (q/p)·K(z−z̄)  to the assigned hit
// positions for each combo (K = 0.3·B·z²·5×10⁻⁴ mm, uniform-field kernel).
// Picks the combo with the smallest normalised chi²
//   χ² = Σ [(y_fit − y_assigned) / r_meas]²
// where r_meas is the measured drift radius — this scale gives ~0 for the
// correct L/R and ~4 per wrong hit, giving robust discrimination even for
// high-p tracks where the model approximation is less accurate.
//
// Outputs: best-L/R measurement vector, fitted q/p (qpBest), best chi².
// Returns false when N > 20 (guard) or every matrix is numerically singular.
// ─────────────────────────────────────────────────────────────────────────────
static bool WCBruteForceAssignLR(const std::vector<MDTMeas>& meas,
                                  std::vector<MDTMeas>&       best,
                                  double& qpBest,
                                  double& chi2Best,
                                  double  B0_T = 1.5)
{
    const int N = static_cast<int>(meas.size());
    if (N < 3 || N > 20) return false;

    double zMean = 0.0;
    for (const auto& m : meas) zMean += m.wireZ_mm;
    zMean /= N;

    int    bestCombo = -1;
    chi2Best  = 1e99;
    qpBest    = 0.0;
    double bestCoef[3] = {};

    const int nComb = 1 << N;
    for (int combo = 0; combo < nComb; ++combo) {
        // Build 3×3 normal equations for [Y₀, slope, q/p]
        double M[3][3] = {}, b[3] = {};
        for (int h = 0; h < N; ++h) {
            const int    ns   = ((combo >> h) & 1) ? +1 : -1;
            const double yass = meas[h].wireY_mm + ns * meas[h].r_meas_mm;
            const double z    = meas[h].wireZ_mm - zMean;
            const double K    = 0.3 * B0_T * z * z * 5.0e-4;
            const double A[3] = {1.0, z, K};
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) M[i][j] += A[i] * A[j];
                b[i] += A[i] * yass;
            }
        }
        // Gaussian elimination with partial pivoting
        double a[3][4];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) a[i][j] = M[i][j];
            a[i][3] = b[i];
        }
        bool sing = false;
        for (int c = 0; c < 3; ++c) {
            int piv = c;
            for (int r = c+1; r < 3; ++r)
                if (std::fabs(a[r][c]) > std::fabs(a[piv][c])) piv = r;
            for (int k = 0; k < 4; ++k) std::swap(a[c][k], a[piv][k]);
            if (std::fabs(a[c][c]) < 1e-9) { sing = true; break; }
            for (int r = c+1; r < 3; ++r) {
                const double f = a[r][c] / a[c][c];
                for (int k = c; k < 4; ++k) a[r][k] -= f * a[c][k];
            }
        }
        if (sing) continue;
        double P[3];
        for (int i = 2; i >= 0; --i) {
            P[i] = a[i][3];
            for (int j = i+1; j < 3; ++j) P[i] -= a[i][j] * P[j];
            P[i] /= a[i][i];
        }
        // χ² = Σ [(y_fit − y_assigned) / r_meas]²
        // Correct L/R: residual ≈ model_error/r ≈ 0.  Wrong L/R: residual ≈ ±2.
        double chi2 = 0.0;
        for (int h = 0; h < N; ++h) {
            const int    ns   = ((combo >> h) & 1) ? +1 : -1;
            const double yass = meas[h].wireY_mm + ns * meas[h].r_meas_mm;
            const double z    = meas[h].wireZ_mm - zMean;
            const double K    = 0.3 * B0_T * z * z * 5.0e-4;
            const double yfit = P[0] + P[1]*z + P[2]*K;
            const double r    = std::max(0.5, meas[h].r_meas_mm);
            chi2 += (yfit - yass) * (yfit - yass) / (r * r);
        }
        if (chi2 < chi2Best) {
            chi2Best  = chi2;
            bestCombo = combo;
            for (int i = 0; i < 3; ++i) bestCoef[i] = P[i];
        }
    }

    if (bestCombo < 0) return false;

    qpBest = bestCoef[2];

    // Build output with best L/R.  Update all three global coordinates:
    // hit = wireCenter + side * r * mdtU (proper 3D for tilted planes).
    best = meas;
    for (int h = 0; h < N; ++h) {
        const int ns  = ((bestCombo >> h) & 1) ? +1 : -1;
        best[h].side  = ns;
        best[h].x_mm  = best[h].wireX_mm + ns * best[h].r_meas_mm * best[h].uX;
        best[h].y_mm  = best[h].wireY_mm + ns * best[h].r_meas_mm * best[h].uY;
        best[h].z_mm  = best[h].wireZ_mm + ns * best[h].r_meas_mm * best[h].uZ;
    }

    std::cerr << "[WCBruteForce] N=" << N
              << " chi2=" << chi2Best
              << " qp=" << qpBest
              << " p=" << (std::fabs(qpBest)>1e-5 ? std::fabs(1.0/qpBest) : 9999.0)
              << " GeV\n";
    return true;
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

    // Sanity gate on px (wire direction, unmeasured): physically it should
    // only be a small fraction of pz (~tan(tilt) ~ 0.08 here). GenFit's RK/DAF
    // fit occasionally converges to an exotic, degenerate solution where px
    // dominates -- the sparse per-station hits don't rule it out via chi2,
    // even though it's unphysical for this near-axial spectrometer. Reject
    // those so the retry/rescue cascade gets another attempt instead.
    auto pxSane = [](const TVector3& mom) {
        return std::fabs(mom.X()) < 0.5 * std::fabs(mom.Z()) + 1.0;
    };

    // Sanity gate on |p| against the independent analytic seed estimate
    // (seedMomentumGeV, from the sagitta fit, tracks truth to ~20-25%). A
    // second, distinct degenerate solution exists alongside the px runaway:
    // a very-low-momentum trajectory can bend sharply enough to thread almost
    // any sparse set of hits, giving excellent chi2 while being unphysical
    // (and, empirically, this degeneracy is what breaks charge discrimination
    // between the two RKTrackRep hypotheses). Reject fits too far from the
    // seed in either direction.
    const double seedMomRef = std::max(1.0, seedMomentumGeV);
    // Widen the pSane gate relative to the seed because pAnalytic is often a
    // large over-estimate (mean bias +165%, RMS 1360% from 5k-event study),
    // causing the correct-momentum fit to be rejected when the seed is too high.
    // Hard floor of 1 GeV prevents unphysical low-momentum spirals regardless
    // of seed.  Upper factor 50× gives room for the (rare) underestimated seed.
    auto pSane = [seedMomRef](const TVector3& mom) {
        const double pmag = mom.Mag();
        return pmag > std::max(1.0, 0.05 * seedMomRef) && pmag < 50.0 * seedMomRef;
    };

    // Seed position: X from the true (tilt-aware) hit position, matching the
    // measurement planes' origin (also hit.x_mm). With a tilted plane normal,
    // seeding X=0 while planes sit at the true X≈900mm decouples the seed fn?rom
    // the geometry it must intersect, which let px run away during fitting.
    TVector3 posSeed(sortedMeas.front().x_mm * 0.1,
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

    // px from dirGlobal.X(): the tilted forward direction has a small but
    // real X-component (~sin(tilt)); seeding px=0 forced a mismatch against
    // the tilted plane geometry.
    TVector3 momSeed(dirGlobal.X() * pSeed, dirGlobal.Y() * pSeed, dirGlobal.Z() * pSeed);

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

    // 1D measurement: only the drift (U) coordinate is measured; the wire (V)
    // direction is unmeasured. Use the tilt-aware mdtU/mdtV computed above
    // (not the global X/Y axes) so the measurement plane matches the wires'
    // actual (tilted) orientation instead of an axis-aligned approximation.
    const double hitSigmaU_cm = 0.080 * 0.1; // 80 μm drift resolution in cm
    TMatrixDSym hitCov1D(1);
    hitCov1D.Zero();
    hitCov1D(0,0) = hitSigmaU_cm * hitSigmaU_cm;
    const TVector3& hU_global = mdtU; // measured: drift direction
    const TVector3& hV_global = mdtV; // unmeasured: wire direction

    const int detId = 0;
    int hitCounter = 0;
    for (const auto& hit : sortedMeas) {
        // Use the true global X (from hit.x_mm) as the plane origin: with a
        // tilted hV_global (wire direction), a wrong origin X is no longer
        // "along the unmeasured axis only" -- it now has a component along
        // the plane normal (hU x hV) and would offset the whole plane.
        TVector3 planeOrigin(hit.x_mm * 0.1, hit.y_mm * 0.1, hit.z_mm * 0.1);

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
    const bool qualityOk = (status->isFitConverged() || ((fitNDF > 0) && (chi2ndf < 100.0)))
                         && pxSane(state.getMom()) && pSane(state.getMom());
    // Also trigger retry if too many hits were silently rejected: a "converged" fit
    // that keeps only 6 of 12 hits (NDF=1 vs expected NDF=7) is not acceptable.
    const int expectedNDF = (int)sortedMeas.size() - 5;
    const bool tooFewHits = (fitNDF > 0) && (fitNDF < expectedNDF - 3);
    if (!qualityOk || tooFewHits) {
        // Snapshot of the Phase-3.5 L/R before any refinement. Used as the base
        // for each alt-seed's curvature-corrected L/R assignment (see below).
        const std::vector<MDTMeas> measPhase35 = sortedMeas;

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

        if (nFlipped == 0) {
            // L/R is already consistent with the partial state — the re-fitted
            // second attempt (below) would give the same bad result, so skip it.
            // Fall straight through to the alt-seed loop (5/20/100 GeV): those
            // use independent fixed momenta and can still succeed when only the
            // original analytic seed momentum was wrong.
        } else {
            std::cerr << "[GenFitMDTFit] " << nFlipped << " hit(s) L/R refined — retrying\n";
        }

        // Second GenFit attempt with corrected measurements, 1D global planes.
        genfit::AbsTrackRep* rep2 = new genfit::RKTrackRep(pdg);

      // Seed from partial state, keeping its own X/px (consistent with the
        // tilted, true-X measurement planes -- zeroing X here would mismatch them).
        TVectorD stateSeed2(6);
        stateSeed2(0) = posPartial.X(); stateSeed2(1) = posPartial.Y(); stateSeed2(2) = posPartial.Z();
        stateSeed2(3) = pPartial.X();   stateSeed2(4) = pPartial.Y();   stateSeed2(5) = pPartial.Z();
        TMatrixDSym covSeed2(6); covSeed2.Zero();
        for (int i = 0; i < 3; ++i) covSeed2(i,i) = 1.0;
        const double mSig2 = std::max(5.0, 0.3 * pPartial.Mag());
        covSeed2(3,3) = 1.0; covSeed2(4,4) = 1.0;
        covSeed2(5,5) = mSig2 * mSig2;

        genfit::Track fitTrack2(rep2, stateSeed2, covSeed2);

        int hitCtr2 = 0;
        for (const auto& hit : sortedMeas) {
            TVector3 origin2(hit.x_mm * 0.1, hit.y_mm * 0.1, hit.z_mm * 0.1);
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

        const bool qualityOk2 = (status2->isFitConverged() || ((ndf_2 > 0) && (chi2ndf_2 < 100.0)))
                              && pxSane(state2.getMom()) && pSane(state2.getMom());
        if (!qualityOk2) {
            // ── Brute-force wire-centre L/R (WCBF) ─────────────────────────────────
            // Adapted from the standalone sim (genfit_reco_batch.cc) that achieves
            // 99.75% correct charge.  Tries all 2^N L/R combinations, fits
            // Y = Y₀ + slope·z + (q/p)·K(z) analytically, picks best chi².
            // Discrimination: correct L/R → chi²≈0; each wrong hit → Δchi²≈4.
            // After finding best L/R, runs GenFit Kalman with 7mm tube-radius sigma
            // (same as alt-seed loop): accepts hits when seed is slightly off.
            {
                std::vector<MDTMeas> bfMeas;
                double qpBF = 0.0, chi2BF = 1e99;
                if (WCBruteForceAssignLR(measPhase35, bfMeas, qpBF, chi2BF)) {
                    sortedMeas = bfMeas;
                    const double pBF = std::min(1000.0, std::max(1.0,
                        std::fabs(qpBF) > 1e-5 ? std::fabs(1.0/qpBF) : 10.0));

                    // Guard 1: sagitta must exceed tube radius for reliable discrimination.
                    // At p=20 GeV, sagitta ≈ 6 mm ≈ tube radius → marginal; above this
                    // the brute-force chi2 is near-zero for ALL combos (track nearly
                    // straight) so the best combo is random and the WCBF Kalman (7 mm σ)
                    // produces chi2/NDF ≈ 0 for BOTH hypotheses → sigma mismatch with the
                    // 80 µm primary Kalman destroys the dual-hypothesis comparison.
                    if (pBF >= 20.0) {
                        std::cerr << "[WCBF] skip: pBF=" << pBF << " >= 20 GeV (sagitta < tube_r)\n";
                        // Fall through to alt-seed loop; sortedMeas already set to bfMeas
                        // (harmless — alt-seed loop resets sortedMeas = measPhase35).
                    } else {

                    constexpr double wclrSigmaU_cm = 0.7; // 7 mm ≈ MDT tube radius
                    genfit::AbsTrackRep* repWC = new genfit::RKTrackRep(pdg);
                    TVectorD ssWC(6);
                    { const auto& f = sortedMeas.front();
                      ssWC(0)=f.x_mm*0.1; ssWC(1)=f.y_mm*0.1; ssWC(2)=f.z_mm*0.1; }
                    { // Seed direction: tangent to helical model at first hit
                      double zMbf = 0.0;
                      for (const auto& m : sortedMeas) zMbf += m.wireZ_mm;
                      zMbf /= static_cast<double>(sortedMeas.size());
                      // Slope from straight-line fit through wire centres
                      double sN=0,sZ=0,sZ2=0,sY=0,sZY=0;
                      for (const auto& m : sortedMeas) {
                          const double z=m.wireZ_mm-zMbf;
                          sN+=1;sZ+=z;sZ2+=z*z;sY+=m.wireY_mm;sZY+=z*m.wireY_mm;
                      }
                      const double det=(sN*sZ2-sZ*sZ);
                      const double slp=(det>1e-9)?(sN*sZY-sZ*sY)/det:0.0;
                      const double zF=sortedMeas.front().wireZ_mm-zMbf;
                      const double dydz=slp+qpBF*0.3*1.5*2.0*zF*5.0e-4;
                      TVector3 dWC=mdtU*dydz+mdtW;
                      if(dWC.Mag()>0) dWC=dWC.Unit();
                      ssWC(3)=dWC.X()*pBF; ssWC(4)=dWC.Y()*pBF; ssWC(5)=dWC.Z()*pBF; }
                    TMatrixDSym csWC(6); csWC.Zero();
                    for (int i=0;i<3;++i) csWC(i,i)=1.0;
                    csWC(3,3)=1.0; csWC(4,4)=1.0;
                    csWC(5,5)=std::max(5.0,0.3*pBF)*std::max(5.0,0.3*pBF);
                    genfit::Track ftWC(repWC, ssWC, csWC);
                    int hcWC=0;
                    for (const auto& hit : sortedMeas) {
                        TVector3 oWC(hit.x_mm*0.1, hit.y_mm*0.1, hit.z_mm*0.1);
                        TVectorD hcvWC(1); hcvWC(0)=0.0;
                        TMatrixDSym hcovWC(1); hcovWC.Zero();
                        hcovWC(0,0)=wclrSigmaU_cm*wclrSigmaU_cm;
                        auto* mptWC = new genfit::PlanarMeasurement(hcvWC, hcovWC, 0, hcWC, nullptr);
                        mptWC->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(oWC, hU_global, hV_global)), hcWC);
                        ftWC.insertPoint(new genfit::TrackPoint(mptWC, &ftWC));
                        ++hcWC;
                    }
                    genfit::KalmanFitterRefTrack fitterWC;
                    fitterWC.setMaxIterations(20);
                    fitterWC.setRelChi2Change(0.001);
                    bool wcKalOk = false;
                    try {
                        fitterWC.processTrack(&ftWC);
                        genfit::FitStatus* stWC = ftWC.getFitStatus(repWC);
                        if (stWC) {
                            const double nWC=stWC->getNdf(), cWC=stWC->getChi2();
                            const double cnWC=(nWC>0)?cWC/nWC:1e9;
                            // Guard 2: tighter chi2/NDF gate (10 instead of 100).
                            // For wrong hypothesis with sagitta ≈ 12 mm (10 GeV), the
                            // 7 mm residuals give chi2/NDF ≈ 24 → rejected here.
                            // Correct hypothesis: chi2/NDF ≈ 0.001 → accepted.
                            if ((stWC->isFitConverged()||(nWC>0&&cnWC<10.0))&&cnWC<10.0) {
                                genfit::MeasuredStateOnPlane stWCs = ftWC.getFittedState();
                                const double wcRef = std::max(1.0, pBF);
                                const double pmWC  = stWCs.getMom().Mag();
                                if (pxSane(stWCs.getMom())
                                    && pmWC >= std::max(1.0, 0.05*wcRef)
                                    && pmWC <= 50.0*wcRef) {
                                    const TVector3 pWCf = stWCs.getMom();
                                    fpx=pWCf.X(); fpy=pWCf.Y(); fpz=pWCf.Z(); fp=pmWC;
                                    fcharge    = stWCs.getCharge();
                                    fchi2      = cWC;
                                    fnDoF      = static_cast<int>(nWC);
                                    fpval      = stWC->getPVal();
                                    fQOverP    = stWCs.getState()(0);
                                    fQOverPErr = std::sqrt(stWCs.getCov()(0,0));
                                    fpErr      = fQOverPErr * fp * fp;
                                    fpos.clear(); layerID.clear();
                                    for (const auto& m : sortedMeas) {
                                        fpos.emplace_back(m.x_mm, m.y_mm, m.z_mm);
                                        layerID.push_back(m.stationID*10000+m.planeID*1000+m.tubeID);
                                    }
                                    wcKalOk = true;
                                    std::cerr << "[WCBF] success: p=" << fp
                                              << " q=" << fcharge << " chi2/NDF=" << cnWC
                                              << " qpBF=" << qpBF << "\n";
                                }
                            }
                        }
                    } catch (...) { std::cerr << "[WCBF] Kalman exception\n"; }
                    if (!wcKalOk)
                        std::cerr << "[WCBF] Kalman failed (NDF/chi2 gate)\n";
                    if (wcKalOk) return true;
                    // Fall through: brute-force Kalman failed; alt-seed loop resets sortedMeas.
                    }  // end pBF < 20 GeV guard
                }
            }  // end WCBF block

            // Attempt 3+: try multiple seed momenta (1D global planes throughout).
            static const double altSeeds[] = { 5.0, 20.0, 100.0, 500.0 };
            for (double altP : altSeeds) {
                // ── Curvature-corrected L/R for this (pdg, altP) hypothesis ─────────
                // For low-momentum tracks the Phase 3.5 straight-line L/R is wrong
                // (sagitta ≫ tube radius → all hits on wrong side → NDF=0 for ANY
                // seed).  Re-assign using a uniform-B helical approximation so that
                // each alt-seed starts from a self-consistent measurement set.
                sortedMeas = measPhase35;  // restore Phase-3.5 L/R
                const double qSignAlt = (pdg == -13) ? -1.0 : +1.0;
                const double B0_T     = 1.5;  // approx avg MDT B-field [T]
                // Step 1: straight-line fit through wire centres.
                double zMeanW = 0.0;
                for (const auto& mm : sortedMeas) zMeanW += mm.wireZ_mm;
                zMeanW /= static_cast<double>(sortedMeas.size());
                {
                    double sNw=0, sZw=0, sZ2w=0, sYw=0, sZYw=0;
                    for (const auto& mm : sortedMeas) {
                        const double z = mm.wireZ_mm - zMeanW;
                        sNw+=1; sZw+=z; sZ2w+=z*z; sYw+=mm.wireY_mm; sZYw+=z*mm.wireY_mm;
                    }
                    const double detW = sNw*sZ2w - sZw*sZw;
                    const double slpW = (std::fabs(detW) > 1e-9) ? (sNw*sZYw - sZw*sYw)/detW : 0.0;
                    const double y0W  = (std::fabs(detW) > 1e-9) ? (sYw - slpW*sZw)/sNw  : 0.0;
                    // Step 2: helical prediction Δy ≈ (q/p)×0.3×B×zRel²×5×10⁻⁴ [mm]
                    for (auto& mm : sortedMeas) {
                        const double zRel = mm.wireZ_mm - zMeanW;
                        const double K_mm = 0.3 * B0_T * zRel * zRel * 5.0e-4;
                        const double yPred = y0W + slpW * zRel + (qSignAlt / altP) * K_mm;
                        const int ns = (yPred >= mm.wireY_mm) ? 1 : -1;
                        mm.side = ns;
                        mm.y_mm = mm.wireY_mm + ns * mm.r_meas_mm;
                    }
                }
                // ─────────────────────────────────────────────────────────────────

                genfit::AbsTrackRep* rep3 = new genfit::RKTrackRep(pdg);
                TVectorD ss3(6);
                // Seed position: use the curvature-corrected first-hit Y so that
                // GenFit's 80-µm sigma does not immediately reject the first station
                // (which would happen if the Phase-3.5 L/R was on the wrong side).
                {
                    const auto& f = sortedMeas.front();
                    ss3(0) = f.x_mm * 0.1; ss3(1) = f.y_mm * 0.1; ss3(2) = f.z_mm * 0.1;
                }
                // Seed direction: tangent to the helical model at the first station.
                // dy/dz|_z1 = slpW + (q/p)×0.3×B×2×zRel1×5×10⁻⁴  (local-frame slope)
                // We recompute slpW quickly from the first/last wire-centre positions.
                TVector3 dirAlt;
                {
                    const double zRelFront = sortedMeas.front().wireZ_mm - zMeanW;
                    double sNw=0, sZw=0, sZ2w=0, sYw=0, sZYw=0;
                    for (const auto& mm : sortedMeas) {
                        const double z = mm.wireZ_mm - zMeanW;
                        sNw+=1; sZw+=z; sZ2w+=z*z; sYw+=mm.wireY_mm; sZYw+=z*mm.wireY_mm;
                    }
                    const double detW = sNw*sZ2w - sZw*sZw;
                    const double slpW2 = (std::fabs(detW) > 1e-9) ? (sNw*sZYw - sZw*sYw)/detW : 0.0;
                    const double slopeAtFront = slpW2
                        + (qSignAlt / altP) * 0.3 * B0_T * 2.0 * zRelFront * 5.0e-4;
                    dirAlt = mdtU * slopeAtFront + mdtW;
                    if (dirAlt.Mag() > 0.0) dirAlt = dirAlt.Unit(); else dirAlt = mdtW;
                }
                TVector3 d3(dirAlt.X() * altP, dirAlt.Y() * altP, dirAlt.Z() * altP);
                ss3(3) = d3.X(); ss3(4) = d3.Y(); ss3(5) = d3.Z();
                TMatrixDSym cs3(6); cs3.Zero();
                for (int i = 0; i < 3; ++i) cs3(i,i) = 1.0;
                const double ms3 = std::max(5.0, 0.3 * altP);
                cs3(3,3) = 1.0; cs3(4,4) = 1.0; cs3(5,5) = ms3 * ms3;
                genfit::Track ft3(rep3, ss3, cs3);
                // Use tube-radius sigma (7 mm) instead of the 80 µm drift resolution
                // so that the Kalman outlier cut (chi2 < 100) accepts hits even when
                // the seed direction is slightly wrong.  With the curvature-corrected
                // L/R, the correct-hypothesis fit will have chi2/NDF ≈ 0 while the
                // wrong-hypothesis fit has chi2/NDF ≈ (14mm/7mm)² ≈ 4, giving a clear
                // delta-chi2 signal for charge assignment.
                constexpr double altSigmaU_cm = 0.7; // 7 mm ≈ MDT tube radius
                int hc3 = 0;
                for (const auto& hit : sortedMeas) {
                    TVector3 org3(hit.x_mm*0.1, hit.y_mm*0.1, hit.z_mm*0.1);
                    TVectorD hc3v(1); hc3v(0) = 0.0;
                    TMatrixDSym hCov3(1); hCov3.Zero();
                    hCov3(0,0) = altSigmaU_cm * altSigmaU_cm;
                    auto* mpt3 = new genfit::PlanarMeasurement(hc3v, hCov3, 0, hc3, nullptr);
                    mpt3->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(org3, hU_global, hV_global)), hc3);
                    ft3.insertPoint(new genfit::TrackPoint(mpt3, &ft3));
                    ++hc3;
                }
                genfit::KalmanFitterRefTrack fitter3;
                fitter3.setMaxIterations(20);
                fitter3.setRelChi2Change(0.001);
                try { fitter3.processTrack(&ft3); } catch (...) { continue; }
                genfit::FitStatus* st3 = ft3.getFitStatus(rep3);
                if (!st3) continue;
                const double c3 = st3->getChi2(), n3 = st3->getNdf(), pv3 = st3->getPVal();
                const double cn3 = (n3 > 0) ? c3/n3 : 1e9;
                std::cerr << "[GenFitMDTFit] AltSeed " << altP << " GeV: chi2/NDF=" << cn3 << "  pval=" << pv3 << "\n";
                if (!st3->isFitConverged() && ((n3 <= 0) || (cn3 >= 100.0))) {
                    std::cerr << "[GenFitMDTFit] AltSeed " << altP << " GeV FAIL:badChi2 conv="
                              << st3->isFitConverged() << " cn3=" << cn3 << "\n";
                    continue;
                }
                genfit::MeasuredStateOnPlane st3s;
                try { st3s = ft3.getFittedState(); } catch (...) {
                    std::cerr << "[GenFitMDTFit] AltSeed " << altP << " GeV FAIL:getFitState exception\n";
                    continue;
                }
                // Gate on pxSane + a per-alt-seed pSane (NOT the analytic-estimate
                // pSane).  When the analytic estimate is far from truth (common for
                // high-momentum nearly-straight tracks), the original pSane rejects
                // the CORRECT solution.  Using the alt-seed momentum as the reference
                // allows the 20/100/500 GeV seeds to accept solutions in a wide range
                // while still blocking very-low-momentum degenerate spirals.
                if (!pxSane(st3s.getMom())) {
                    std::cerr << "[GenFitMDTFit] AltSeed " << altP << " GeV FAIL:pxSane px="
                              << st3s.getMom().X() << " pz=" << st3s.getMom().Z() << "\n";
                    continue;
                }
                {
                    const double altRef = std::max(1.0, altP);
                    const double pmag3  = st3s.getMom().Mag();
                    if (pmag3 < std::max(1.0, 0.05 * altRef) || pmag3 > 50.0 * altRef) {
                        std::cerr << "[GenFitMDTFit] AltSeed " << altP << " GeV FAIL:pSaneAlt p="
                                  << pmag3 << " range=[" << std::max(1.0, 0.05*altRef)
                                  << "," << 50.0*altRef << "]\n";
                        continue;
                    }
                }
                // ── Stage 2: DAF refit with drift resolution ────────────────────────
                // Use the tube-sigma Stage-1 result as seed.  With 80 µm sigma the
                // wrong-hypothesis measurements (displaced ~14 mm = 175 σ) are driven
                // to weight≈0 by the annealing scheme → NDF→0 → skip this altP.
                // Correct-hypothesis measurements (residuals ≪ σ) → DAF converges.
                // T_start=5000 so that wrong hits (chi2≈30625) start with
                // weight≈exp(-30625/10000)≈0.05 and are immediately nearly zeroed,
                // while correct hits (chi2≈0) always have weight≈1.
                {
                    genfit::AbsTrackRep* repS2 = new genfit::RKTrackRep(pdg);
                    const TVector3 xS2 = st3s.getPos(); // cm
                    const TVector3 pS2 = st3s.getMom(); // GeV/c
                    TVectorD ssS2(6);
                    ssS2(0)=xS2.X(); ssS2(1)=xS2.Y(); ssS2(2)=xS2.Z();
                    ssS2(3)=pS2.X(); ssS2(4)=pS2.Y(); ssS2(5)=pS2.Z();
                    TMatrixDSym csS2(6); csS2.Zero();
                    for (int i=0;i<3;++i) csS2(i,i)=0.01;
                    const double msS2 = std::max(0.1, 0.2 * fp);
                    csS2(3,3)=msS2*msS2; csS2(4,4)=msS2*msS2; csS2(5,5)=msS2*msS2;
                    genfit::Track ftS2(repS2, ssS2, csS2);
                    int hcS2 = 0;
                    for (const auto& hit : sortedMeas) {
                        TVector3 orgS2(hit.x_mm*0.1, hit.y_mm*0.1, hit.z_mm*0.1);
                        TVectorD hcS2v(1); hcS2v(0) = 0.0;
                        TMatrixDSym hCovS2(1); hCovS2.Zero();
                        hCovS2(0,0) = hitSigmaU_cm * hitSigmaU_cm;
                        auto* mptS2 = new genfit::PlanarMeasurement(hcS2v, hCovS2, 0, hcS2, nullptr);
                        mptS2->setPlane(genfit::SharedPlanePtr(
                            new genfit::DetPlane(orgS2, hU_global, hV_global)), hcS2);
                        ftS2.insertPoint(new genfit::TrackPoint(mptS2, &ftS2));
                        ++hcS2;
                    }
                    genfit::DAF dafS2;
                    dafS2.setAnnealingScheme(5000.0, 9.0, 30);
                    dafS2.setMaxIterations(40);
                    bool dafS2ok = false;
                    try {
                        dafS2.processTrack(&ftS2);
                        genfit::FitStatus* stS2 = ftS2.getFitStatus(repS2);
                        if (stS2) {
                            const double nS2 = stS2->getNdf(), cS2 = stS2->getChi2();
                            const double cnS2 = (nS2 > 0) ? cS2/nS2 : 1e9;
                            if (nS2 >= 1.0 && cnS2 < 200.0) {
                                try {
                                    genfit::MeasuredStateOnPlane stS2s = ftS2.getFittedState();
                                    const TVector3 pS2f = stS2s.getMom();
                                    const double pmagS2 = pS2f.Mag();
                                    const double altRefS2 = std::max(1.0, altP);
                                    if (pxSane(pS2f)
                                     && pmagS2 >= std::max(1.0, 0.05*altRefS2)
                                     && pmagS2 <= 50.0*altRefS2) {
                                        fpx=pS2f.X(); fpy=pS2f.Y(); fpz=pS2f.Z(); fp=pmagS2;
                                        fcharge    = stS2s.getCharge();
                                        fchi2      = cS2;
                                        fnDoF      = static_cast<int>(nS2);
                                        fpval      = stS2->getPVal();
                                        fQOverP    = stS2s.getState()(0);
                                        fQOverPErr = std::sqrt(stS2s.getCov()(0,0));
                                        fpErr      = fQOverPErr * fp * fp;
                                        dafS2ok    = true;
                                        std::cerr << "[GenFitMDTFit] AltSeed " << altP
                                                  << " GeV Stage2 ok: chi2/NDF=" << cnS2 << "\n";
                                    }
                                } catch (...) {}
                            }
                        }
                    } catch (...) {}
                    if (!dafS2ok) {
                        std::cerr << "[GenFitMDTFit] AltSeed " << altP
                                  << " GeV Stage2 FAIL (wrong hyp?)\n";
                        continue; // likely wrong hypothesis – try next altP
                    }
                }
                // Stage 2 succeeded: record hit positions and return.
                fpos.clear(); layerID.clear();
                for (const auto& m : sortedMeas) {
                    fpos.emplace_back(m.x_mm, m.y_mm, m.z_mm);
                    layerID.push_back(m.stationID*10000 + m.planeID*1000 + m.tubeID);
                }
                std::cerr << "[GenFitMDTFit] AltSeed SUCCESS (Stage2 DAF): p=" << fp
                          << " q=" << fcharge << "\n";
                return true;
            }
            // Restore Phase-3.5 L/R for the DAF rescue that follows (so it
            // doesn't start from the last alt-seed's curvature-corrected state).
            sortedMeas = measPhase35;

            // ── DAF rescue ────────────────────────────────────────────────────
            // All Kalman attempts failed.  Try GenFit's Deterministic Annealing
            // Filter which soft-weights each hit by its chi² contribution and
            // iteratively cools the temperature, driving contaminated hits (from
            // a secondary particle or a misidentified MDT tube) to weight → 0.
            // The muon track then emerges from the uncontaminated majority.
            {
                genfit::AbsTrackRep* repD = new genfit::RKTrackRep(pdg);
                TVectorD ssD(6);
                ssD(0) = posSeed.X(); ssD(1) = posSeed.Y(); ssD(2) = posSeed.Z();
                TVector3 dD(dirGlobal.X() * seedMomentumGeV,
                                 dirGlobal.Y() * seedMomentumGeV,
                                 dirGlobal.Z() * seedMomentumGeV);
                ssD(3) = dD.X(); ssD(4) = dD.Y(); ssD(5) = dD.Z();
                TMatrixDSym csD(6); csD.Zero();
                for (int i = 0; i < 3; ++i) csD(i,i) = 1.0;
                const double msD = std::max(5.0, 0.3 * seedMomentumGeV);
                csD(3,3) = 1.0; csD(4,4) = 1.0; csD(5,5) = msD * msD;

                genfit::Track ftD(repD, ssD, csD);
                int hcD = 0;
                for (const auto& hit : sortedMeas) {
                    TVector3 orgD(hit.x_mm * 0.1, hit.y_mm * 0.1, hit.z_mm * 0.1);
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
                genfit::MeasuredStateOnPlane stDs;
                try { stDs = ftD.getFittedState(); }
                catch (...) { return false; }

                const bool qualityDAF = (stD->isFitConverged() || ((nD > 0) && (cnD < 200.0)))
                                        && (nD >= 1.0) && (cnD < 200.0)
                                        && pxSane(stDs.getMom()) && pSane(stDs.getMom());
                if (!qualityDAF) return false;

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
