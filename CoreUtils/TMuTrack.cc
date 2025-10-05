#include "TMuTrack.hh"
#include "GenMagneticField.hh"

// genfit
#include <RKTrackRep.h>
#include <Track.h>
#include <TrackPoint.h>
#include <PlanarMeasurement.h>
#include <KalmanFitterRefTrack.h>
#include <FitStatus.h>

ClassImp(TMuTrack)

void TMuTrack::GenFitTrackFit(int verbose, double detectorResolutionPSmm) {

    const int pdg = 13;    

// trackrep
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

    // start values for the fit, e.g. from pattern recognition
    // position of first hit of the track
    TVector3 pos(fpos[0].x()/10.0, fpos[0].y()/10.0, fpos[0].z()/10.0);
    double pmom = 10.0*1e3;  // in MeV??
    TVector3 mom(0, 0, pmom);

    // create track
    fitTrack = new genfit::Track(rep, pos, mom);

    const int detId(0); // detector ID
    int planeId(0); // detector plane ID
    int hitId(0); // hit ID

    // Scifi resolution
    double detectorResolution(0.1); // resolution of planar detectors in cm 
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

    // do the fit
    try {
        fitter->processTrack(fitTrack);
    }    
    catch(genfit::Exception& e){
        std::cerr << e.what();
        std::cerr << "Exception when track fitting with GENFIT" << std::endl;
      }

    if(verbose > 3) {
        std::cout << "After fit: " << std::endl;
        fitTrack->getFittedState().Print();
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
    genfit::MeasuredStateOnPlane state = fitTrack->getFittedState();
    TVector3 p = state.getMom();
    if (verbose>0) {
        std::cout << "fitted momentum (GeV/c): " << p.Mag() << " px: " << p.X() << " py: " << p.Y() << " pz: " << p.Z();
        std::cout << " pT: " << sqrt(p.X()*p.X()+p.Y()*p.Y()) << std::endl;
    }
    fpx = p.X();
    fpy = p.Y();
    fpz = p.Z();
    fp = p.Mag();

    // clean up (don't delete rep because it is owned by the track)
    delete fitter;
}
