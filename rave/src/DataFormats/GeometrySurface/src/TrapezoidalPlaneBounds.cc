////#include \"Utilities/Configuration/interface/Architecture.h"

#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"
#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include <cmath>

TrapezoidalPlaneBounds::TrapezoidalPlaneBounds( float be, float te, 
						float a, float t) : 
  hbotedge(be), htopedge(te), hapothem(a), hthickness(t) {

  // pre-compute offset of triangle vertex and tg of (half) opening
  // angle of the trapezoid for faster inside() implementation.

  offset = a * (te+be) / (te-be);  // check sign if te < be !!! 
  tan_a = te / (offset + a);
}

bool TrapezoidalPlaneBounds::inside( const Local2DPoint& p) const {
  return fabs(p.y()) <= hapothem && 
    fabs(p.x()/(p.y()+offset)) <= tan_a;
}

bool TrapezoidalPlaneBounds::inside( const Local3DPoint& p) const {
  return fabs(p.y()) <= hapothem &&
    fabs(p.x()/(p.y()+offset)) <= tan_a &&
    fabs(p.z()) <= hthickness;
}

bool TrapezoidalPlaneBounds::inside( const Local3DPoint& p,
				     const LocalError& err, float scale) const {
  TrapezoidalPlaneBounds tmp( hbotedge + sqrt(err.xx())*scale,
			      htopedge + sqrt(err.xx())*scale,
			      hapothem + sqrt(err.yy())*scale,
			      hthickness);
  return tmp.inside(p);
}
  

const std::vector<float> TrapezoidalPlaneBounds::parameters() const { 
  std::vector<float> vec(4);
  // Same order as geant3 for constructor compatibility
  vec[0] = hbotedge;
  vec[1] = htopedge;
  vec[3] = hapothem;
  vec[2] = hthickness;
  return vec;
}
