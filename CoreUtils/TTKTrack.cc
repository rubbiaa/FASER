#include <iostream>
#include <cmath>

#include "TTKTrack.hh"
#include "TcalEvent.hh"

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
#include <FitStatus.h>

ClassImp(TTKTrack);

int TTKTrackRefenceID = 0;

TTKTrack::TTKTrack() : TObject(), fitTrack(0), vertexID(-1), trackID(-1) { 
    TTKTrackRefenceID++;
        //std::cout << "TTKTrack::TTKTrack - constructor . " << this << std::endl; 
};

TTKTrack::TTKTrack(const TTKTrack &t) : TTKTrack() {
    tkhit = t.tkhit;
    centroid = t.centroid;
    direction = t.direction;
    SSR = t.SSR;
    trackID = t.trackID;
    vertexID = t.vertexID;
    if(t.fitTrack) {
        fitTrack = new genfit::Track(*t.fitTrack); // FIXME: this creates leaks; should use clone???
//        std::cout << "TTKTrack::TTKTrack - copy constructor fitrack" << fitTrack << std::endl;
    }
}

TTKTrack::~TTKTrack() { 
//        std::cout << "TTKTrack::~TTKTrack - destructor . " << this << std::endl;
    if(fitTrack!=nullptr) {
        fitTrack->deleteFitterInfo();
        // delete track represemntations
        for (int i=0; i<fitTrack->getNumReps(); ++i) {
            fitTrack->deleteTrackRep(i);
        }
        delete fitTrack;
    }
    TTKTrackRefenceID--;
//    std::cout << "TTKTrack::~TTKTrack - destructor . " << this << " " << TTKTrackRefenceID << std::endl;
};

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

void TTKTrack::GenFitTrackFit(double detectorResolutionPSmm) {

    const int pdg = 13;       
// trackrep
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

    // start values for the fit, e.g. from pattern recognition
    // position of first hit of the track
    TVector3 pos(tkhit[0].point.x()/10.0, tkhit[0].point.y()/10.0, tkhit[0].point.z()/10.0);
    double pmom = 10.0*1e0;  // in MeV??
    TVector3 mom(pmom*direction.x(), pmom*direction.y(), pmom*direction.z());

    // create track
    fitTrack = new genfit::Track(rep, pos, mom);

    const int detId(0); // detector ID
    int planeId(0); // detector plane ID
    int hitId(0); // hit ID

    // tracker resolution
    double detectorResolution(0.01); // 0.1); // resolution of planar detectors in cm 
    TMatrixDSym hitCov(2);
    hitCov.UnitMatrix();
    hitCov *= detectorResolution*detectorResolution;
    // scintillator resolution
    double detectorResolutionPS = detectorResolutionPSmm/10.0; // resolution of planar detectors in cm
    TMatrixDSym hitCovPS(2);
    hitCovPS.UnitMatrix();
    hitCovPS *= detectorResolutionPS*detectorResolutionPS;

    for (const auto &hit : tkhit) {
        TVectorD hitCoords(2);
        hitCoords[0] = hit.point.x()/10.0;
        hitCoords[1] = hit.point.y()/10.0;
        genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, 
            (hit.type == 0) ? hitCovPS  : hitCov, 
            detId, ++hitId, nullptr);
        measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,hit.point.z()/10.0),
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

    // fitTrack->getFittedState().Print();
    // fitTrack->Print();
    delete fitter;
}

TVector3 TTKTrack::extrapolateTracktoZ(double z, int &failed) {
    failed = 0;
    if(!fitTrack) {
        std::cerr << "TTKTrack::extrapolateTracktoZ - no fitTrack" << std::endl;
        failed = 1;
        return TVector3(0,0,0);
    }
    // check that track is fitted
    if(fitTrack->getFitStatus()->isFitted() == false) {
        std::cerr << "TTKTrack::extrapolateTracktoZ - track is not fitted" << std::endl;
        failed = 1;
        return TVector3(0,0,0);
    }

    int planeId(0);

    genfit::StateOnPlane state = fitTrack->getFittedState();
    genfit::AbsTrackRep *rep = fitTrack->getCardinalRep();
    if(!rep) {
        std::cerr << "TTKTrack::extrapolateTracktoZ - no rep" << std::endl;
        failed = 1;
        return TVector3(0,0,0);
    }
    try {
        double length = rep -> extrapolateToPlane(
        state, genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,z/10.0),
                        TVector3(1,0,0), TVector3(0,1,0))));
    } catch(genfit::Exception& e){
        std::cerr << e.what();
        std::cerr << "Exception when track extrapolation to z " << z << std::endl;
        failed = 1;
        return TVector3(0,0,0);
    }
/*    if(length < 0) {
        std::cerr << "TTKTrack::extrapolateTracktoZ - extrapolation failed" << std::endl;
        return TVector3(0,0,0);
    }
    */
    TVector3 pos = state.getPos();
    return TVector3(pos.x()*10.0, pos.y()*10.0, pos.z()*10.0);
}

void TTKTrack::GetMinMaxModule(int &minLayer, int &maxLayer) {
    minLayer = 1000;
    maxLayer = -1;
    for (const auto &hit : tkhit) {
        long ID = hit.ID;
        long ilayer = TcalEvent::getChannelModulefromID(ID);
        if(ilayer < minLayer) minLayer = ilayer;
        if(ilayer > maxLayer) maxLayer = ilayer;
    }
}

void TTKTrack::UpdateFittedPosition() {
    if(!fitTrack) {
        std::cerr << "TTKTrack::UpdateFittedPosition - no fitTrack" << std::endl;
        return;
    }
    if(fitTrack->getFitStatus()->isFitted() == false) {
        std::cerr << "TTKTrack::UpdateFittedPosition - track is not fitted" << std::endl;
        return;
    }
    // try to refit the complete track
    int failed(0);
    for (auto &hit : tkhit) {
        TVector3 pos = extrapolateTracktoZ(hit.point.z(), failed);
        if(!failed) {
            hit.point.SetXYZ(pos.x(), pos.y(), pos.z());
        } else {
            std::cerr << "TTKTrack::UpdateFittedPosition - extrapolation failed" << std::endl;
            break;
        }
    }
}

void TTKTrack::Dump(int verbose) const {
    std::cout << "TKTrack " << trackID << ": " << tkhit.size() << " hits; SSR = " << SSR << std::endl;
    std::cout << "direction " << direction.x() << " " << direction.y() << " " << direction.z() << std::endl;
    for (const auto &it : tkhit) {
        std::cout << "   " << it.ID << " x:" << Form("%5.2f",it.point.x()) << 
        " y:" << Form("%5.2f",it.point.y()) << " z:" << Form("%5.2f",it.point.z())
        << " eDeposit: " << it.eDeposit << std::endl;
        if(it.eDeposit==0) {
            std::cout << "ERROR: zero energy deposit" << std::endl;
        }
    }        
    if(fitTrack) {
//        std::cout << "track " << this << "GenFit track: " << fitTrack;
//        std::cout << "cardinal rep : " << fitTrack->getCardinalRep() << std::endl;
        if(verbose > 4) fitTrack->Print();
        genfit::FitStatus* getFitStatus = fitTrack->getFitStatus();
        if(getFitStatus) {
            getFitStatus->Print();
        }
    }
}
