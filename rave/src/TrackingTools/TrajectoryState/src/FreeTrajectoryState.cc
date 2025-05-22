#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToCartesian.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCartesianToCurvilinear.h"

#include <cmath>

// implementation of non-trivial methods of FreeTrajectoryState

// Warning: these methods violate constness

// convert curvilinear errors to cartesian
void FreeTrajectoryState::createCartesianError() const{
  
  JacobianCurvilinearToCartesian curv2Cart(theGlobalParameters);
  const AlgebraicMatrix65& jac = curv2Cart.jacobian();

  ((FreeTrajectoryState*)this)->theCartesianError = 
    ROOT::Math::Similarity(jac, theCurvilinearError.matrix());

  ((FreeTrajectoryState*)this)->theCartesianErrorValid = true;
}

// convert cartesian errors to curvilinear
void FreeTrajectoryState::createCurvilinearError() const{
  
  JacobianCartesianToCurvilinear cart2Curv(theGlobalParameters);
  const AlgebraicMatrix56& jac = cart2Curv.jacobian();
  
  ((FreeTrajectoryState*)this)->theCurvilinearError = 
    ROOT::Math::Similarity(jac, theCartesianError.matrix());
  ((FreeTrajectoryState*)this)->theCurvilinearErrorValid = true;
} 

// check if trajectory can reach given radius

bool FreeTrajectoryState::canReach(double radius) const {
  GlobalPoint x = position();
  GlobalVector p = momentum().unit();
  double rho = transverseCurvature()*p.perp();
  double rx = rho*x.x();
  double ry = rho*x.y();
  double rr = rho*radius;
  double ax = p.x()*rx + p.y()*ry;
  double ay = p.x()*ry - p.y()*rx + p.perp2();
  double cospsi = (.5*(rx*rx + ry*ry - rr*rr) + ay)/sqrt(ax*ax + ay*ay);
  return fabs(cospsi) <= 1.;
}







