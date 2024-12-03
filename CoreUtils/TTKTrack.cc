#include <iostream>
#include <cmath>

#include "TTKTrack.hh"

#include <TVector3.h>
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TMinuit.h"
#include "Fit/FitResult.h"
#include <TMatrixDSym.h>
#include <TDecompSVD.h>

// genfit
#include <RKTrackRep.h>
#include <Track.h>
#include <TrackPoint.h>
#include <PlanarMeasurement.h>
#include <KalmanFitterRefTrack.h>

ClassImp(TTKTrack);

TTKTrack::TTKTrack(const TTKTrack &t) : fitTrack(0) {
    tkhit = t.tkhit;
    centroid = t.centroid;
    direction = t.direction;
    SSR = t.SSR;
}

double TTKTrack::pointLineDistance(const ROOT::Math::XYZVector& point, const TVector3& direction, const TVector3& centroid) {
    TVector3 pointVec(point.x(), point.y(), point.z());
    TVector3 pointToCentroid = pointVec - centroid;
    TVector3 crossProduct = pointToCentroid.Cross(direction);
    return crossProduct.Mag() / direction.Mag();
}

void TTKTrack::MergeTracks(TTKTrack &track2) {
    // add hits from track2 to track1
    for (const auto& h : track2.tkhit) {
        tkhit.push_back(h);
    }
    direction = fitLineThroughHits(centroid);
}

TVector3 TTKTrack::direction2Hits() {
    TVector3 direction;
    direction.SetXYZ(0,0,0);
    if(tkhit.size()>1) {
        struct TRACKHIT hit1 = tkhit[0];
        struct TRACKHIT hit2 = tkhit[1];
        direction.SetXYZ(hit2.point.x() - hit1.point.x(), 
                         hit2.point.y() - hit1.point.y(), 
                         hit2.point.z() - hit1.point.z());
        direction = direction.Unit();
    }
    SSR = -1;
    return direction;
}

TVector3 TTKTrack::fitLineThroughHits(TVector3& centroid) {
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

    // compute SSR
    SSR = 0.0;
    for (const auto& hit : tkhit) {
        double distance = pointLineDistance(hit.point, direction, centroid);
        SSR += distance * distance;
    }

    return direction;
}

void TTKTrack::GenFitTrackFit() {

    const int pdg = 13;       
// trackrep
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

    // start values for the fit, e.g. from pattern recognition
    // position of first hit of the track
    TVector3 pos(tkhit[0].point.x(), tkhit[0].point.y(), tkhit[0].point.z());
    double pmom = 100.0*1e3;  // in MeV
    TVector3 mom(pmom*direction.x(), pmom*direction.y(), pmom*direction.z());

    // create track
    fitTrack = new genfit::Track(rep, pos, mom);

    const int detId(0); // detector ID
    int planeId(0); // detector plane ID
    int hitId(0); // hit ID

    double detectorResolution(0.1); // resolution of planar detectors
    TMatrixDSym hitCov(2);
    hitCov.UnitMatrix();
    hitCov *= detectorResolution*detectorResolution;

    for (const auto &hit : tkhit) {
        TVectorD hitCoords(2);
        hitCoords[0] = hit.point.x();
        hitCoords[1] = hit.point.y();
        genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, nullptr);
        measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,hit.point.z()),
                     TVector3(1,0,0), TVector3(0,1,0))), ++planeId);
        fitTrack->insertPoint(new genfit::TrackPoint(measurement, fitTrack));
    }

    fitTrack->checkConsistency();

     // init fitter
    genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();

    // do the fit
    fitter->processTrack(fitTrack);

    // print fit result
//    fitTrack.getFittedState().Print();
//    fitTrack.Print();
}

void TTKTrack::Dump() const {
    std::cout << "TKTrack: " << tkhit.size() << " hits; SSR = " << SSR << std::endl;
    std::cout << "direction " << direction.x() << " " << direction.y() << " " << direction.z() << std::endl;
    for (const auto &it : tkhit) {
        std::cout << "   " << it.ID << " x:" << Form("%5.2f",it.point.x()) << 
        " y:" << Form("%5.2f",it.point.y()) << " z:" << Form("%5.2f",it.point.z())
        << " eDeposit: " << it.eDeposit << std::endl;
    }        
    if(fitTrack) {
        fitTrack->Print();
    }
}
