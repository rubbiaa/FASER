#ifndef Geom_TangentPlane_H
#define Geom_TangentPlane_H

#include "DataFormats/GeometrySurface/interface/Plane.h"
//#include "Utilities/GenUtil/interface/ReferenceCountingPointer.h"

/** Plane tangent to a more general surface (e.g. cylinder).
 *  To be constructed by the "parent" surface.
 */

class TangentPlane : public Plane {
public:
  TangentPlane (const PositionType& pos, 
		const RotationType& rot, 
		const Surface* parent) :
    Surface(pos,rot), Plane(pos,rot),
    theParent(parent) {}

  /// access to original surface
  const Surface& parentSurface() {return *theParent;}

private:
  ConstReferenceCountingPointer<Surface> theParent;

};
#endif
