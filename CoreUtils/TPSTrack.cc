#include <iostream>
#include <cmath>

#include "TPSTrack.hh"

#include <TVector3.h>
#include "TMath.h"
#include "TMinuit.h"
#include "Fit/FitResult.h"
#include <TMatrixDSym.h>
#include <TDecompSVD.h>

ClassImp(TPSTrack);

TVector3 TPSTrack::direction2Hits(TVector3& centroid) {
    TVector3 direction;
    direction.SetXYZ(0,0,0);
    if(tkhit.size()>1) {
        struct TRACKHIT hit1 = tkhit[0];
        struct TRACKHIT hit2 = tkhit[1];
        direction.SetXYZ(hit2.point.x() - hit1.point.x(), 
                         hit2.point.y() - hit1.point.y(), 
                         hit2.point.z() - hit1.point.z());
        direction = direction.Unit();
        ROOT::Math::XYZVector centroid = (hit1.point + hit2.point) * 0.5;
        centroid.SetXYZ(centroid.X(), centroid.Y(), centroid.Z());
    }
    return direction;
}

TVector3 TPSTrack::fitLineThroughHits(TVector3& centroid) {
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

bool TPSTrack::VoxelTouchesTrack(long ID) {
    bool touches = false;
    long ix = ID % 1000;
    long iy = (ID / 1000) % 1000;
    long iz = (ID / 1000000) % 1000;
    long ilayer = (ID / 1000000000);
    for (const auto& hit : tkhit) {
        long ID2 = hit.ID;
        long ix2 = ID2 % 1000;
        long iy2 = (ID2 / 1000) % 1000;
        long iz2 = (ID2 / 1000000) % 1000;
        long ilayer2 = (ID2 / 1000000000);
        if(ilayer!=ilayer2 || iz != iz2) continue;
        int dx = abs(ix-ix2);
        int dy = abs(iy-iy2);
        if( dx < 2 && dy < 2 ) {
            touches = true;
            break;
        }
    }
    return touches;
}
