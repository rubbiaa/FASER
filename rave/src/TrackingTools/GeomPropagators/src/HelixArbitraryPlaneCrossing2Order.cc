#include "TrackingTools/GeomPropagators/interface/HelixArbitraryPlaneCrossing2Order.h"

#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include <cmath>
#include <cfloat>

HelixArbitraryPlaneCrossing2Order::HelixArbitraryPlaneCrossing2Order(const PositionType& point,
								     const DirectionType& direction,
								     const float curvature,
								     const PropagationDirection propDir) :
  theX0(point.x()),
  theY0(point.y()),
  theZ0(point.z()),
  theRho(curvature),
  thePropDir(propDir)
{
  //
  // Components of direction vector (with correct normalisation)
  //
  double px = direction.x();
  double py = direction.y();
  double pz = direction.z();
  double pt = px*px+py*py;
  double p = sqrt(pt+pz*pz);
  pt = sqrt(pt);
  theCosPhi0 = px/pt;
  theSinPhi0 = py/pt;
  theCosTheta = pz/p;
  theSinTheta = pt/p;
}

//
// Propagation status and path length to intersection
//
std::pair<bool,double>
HelixArbitraryPlaneCrossing2Order::pathLength(const Plane& plane) {
  //
  // get local z-vector in global co-ordinates and
  // distance to starting point
  //
  GlobalVector normalToPlane = plane.normalVector();
  double nPx = normalToPlane.x();
  double nPy = normalToPlane.y();
  double nPz = normalToPlane.z();
  double cP = plane.localZ(GlobalPoint(theX0,theY0,theZ0));
  //
  // coefficients of 2nd order equation to obtain intersection point
  // in approximation (without curvature-related factors).
  //
  double ceq1 = theRho*(nPx*theSinPhi0-nPy*theCosPhi0);
  double ceq2 = nPx*theCosPhi0 + nPy*theSinPhi0 + nPz*theCosTheta/theSinTheta;
  double ceq3 = cP;
  //
  // Check for degeneration to linear equation (zero
  //   curvature, forward plane or direction perp. to plane)
  //
  double dS1,dS2;
  if ( fabs(ceq1)>FLT_MIN ) {
    double deq1 = ceq2*ceq2;
    double deq2 = ceq1*ceq3;
    if ( fabs(deq1)<FLT_MIN || fabs(deq2/deq1)>1.e-6 ) {
      //
      // Standard solution for quadratic equations
      //
      double deq = deq1+2*deq2;
      if ( deq<0. )  return std::pair<bool,double>(false,0);
      double ceq = -0.5*(ceq2+(ceq2>0?1:-1)*sqrt(deq));
      dS1 = -2*ceq/ceq1/theSinTheta;
      dS2 = ceq3/ceq/theSinTheta;
    }
    else {
      //
      // Solution by expansion of sqrt(1+deq)
      //
      double ceq = ceq2/ceq1/theSinTheta;
      double deq = deq2/deq1;
      deq *= (1-deq/2);
      dS1 = -ceq*deq;
      dS2 = ceq*(2+deq);
    }
  }
  else {
    //
    // Special case: linear equation
    //
    dS1 = dS2 = -ceq3/ceq2*(1/theSinTheta);
  }
  //
  // Choose and return solution
  //
  return solutionByDirection(dS1,dS2);
}
//
// Position after a step of path length s (2nd order)
//
HelixPlaneCrossing::PositionType
HelixArbitraryPlaneCrossing2Order::position (double s) const {
  // use double precision result
  PositionTypeDouble pos = positionInDouble(s);
  return PositionType(pos.x(),pos.y(),pos.z());
}
//
// Position after a step of path length s (2nd order) (in double precision)
//
HelixArbitraryPlaneCrossing2Order::PositionTypeDouble
HelixArbitraryPlaneCrossing2Order::positionInDouble (double s) const {
  // based on path length in the transverse plane
  double st = s*theSinTheta;
  return PositionTypeDouble(theX0+(theCosPhi0-st*0.5*theRho*theSinPhi0)*st,
			    theY0+(theSinPhi0+st*0.5*theRho*theCosPhi0)*st,
			    theZ0+st*theCosTheta/theSinTheta);
}
//
// Direction after a step of path length 2 (2nd order) (in double precision)
//
HelixPlaneCrossing::DirectionType
HelixArbitraryPlaneCrossing2Order::direction (double s) const {
  // use double precision result
  DirectionTypeDouble dir = directionInDouble(s);
  return DirectionType(dir.x(),dir.y(),dir.z());
}
//
// Direction after a step of path length 2 (2nd order)
//
HelixArbitraryPlaneCrossing2Order::DirectionTypeDouble
HelixArbitraryPlaneCrossing2Order::directionInDouble (double s) const {
  // based on delta phi
  double dph = s*theRho*theSinTheta;
  return DirectionTypeDouble(theCosPhi0-(theSinPhi0+0.5*theCosPhi0*dph)*dph,
			     theSinPhi0+(theCosPhi0-0.5*theSinPhi0*dph)*dph,
			     theCosTheta/theSinTheta);
}
//
// Choice of solution according to propagation direction
//
std::pair<bool,double>
HelixArbitraryPlaneCrossing2Order::solutionByDirection(const double dS1,
						       const double dS2) const {
  bool valid = false;
  double path = 0;
  if ( thePropDir == anyDirection ) {
    valid = true;
    path = smallestPathLength(dS1,dS2);
  }
  else {
    // use same logic for both directions (invert if necessary)
    int propSign = thePropDir==alongMomentum ? 1 : -1;
    double s1(propSign*dS1);
    double s2(propSign*dS2);
    // sort
    if ( s1 > s2 ) {
      double aux = s1;
      s1 = s2;
      s2 = aux;
    }
    // choose solution (if any with positive sign)
    if ( s1<0 && s2>=0 ) {
      // First solution in backward direction: choose second one.
      valid = true;
      path = propSign*s2;
    }
    else if ( s1>=0 ) {
      // First solution in forward direction: choose it (s2 is further away!).
      valid = true;
      path = propSign*s1;
    }
  }
  return std::pair<bool,double>(valid,path);
}
