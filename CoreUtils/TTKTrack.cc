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

ClassImp(TTKTrack);

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

    return direction;
}
