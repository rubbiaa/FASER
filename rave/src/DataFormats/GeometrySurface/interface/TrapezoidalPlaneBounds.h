#ifndef Geom_TrapezoidalPlaneBounds_H
#define Geom_TrapezoidalPlaneBounds_H


#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometrySurface/interface/Bounds.h"
#include <vector>

/**
 *  Trapezoidal detector active volume.
 *  Local Coordinate system coincides with center of the box
 *  with y axis being the symmetry axis along the height
 *  and pointing in the direction of top_edge.
 */

class TrapezoidalPlaneBounds : public Bounds {
public:

  /** constructed from:
   *    half bottom edge (smaller side width)
   *    half top edge    (larger side width)
   *    half apothem (distance from top to bottom sides, 
   *                  measured perpendicularly to them)
   *    half thickness
   */
  TrapezoidalPlaneBounds( float be, float te, float a, float t);
  

  /** apothem (full, not half)*/
  virtual float length() const    { return 2 * hapothem;}

  /** largest width (full, not half)*/
  virtual float width()  const    { return 2 * std::max( hbotedge, htopedge);}

  /** thickness (full, not half)*/
  virtual float thickness() const { return 2 * hthickness;}

  /** Width at half length. Useful for e.g. pitch definition.
   */
  virtual float widthAtHalfLength() const {return hbotedge+htopedge;}

  virtual bool inside( const Local2DPoint& p) const;

  virtual bool inside( const Local3DPoint& p) const;

  virtual bool inside( const Local3DPoint& p, const LocalError& err, float scale) const;

  virtual bool inside( const Local2DPoint& p, const LocalError& err, float scale) const {
    return Bounds::inside(p,err,scale);
  }

  /** returns the 4 parameters needed for construction, in the order
   * ( half bottom edge, half top edge, half thickness, half apothem).
   * Beware! This order is different from the one in the constructor!
   */
  virtual const std::vector<float> parameters() const;

  virtual Bounds* clone() const { 
    return new TrapezoidalPlaneBounds(*this);
  }

private:
  // persistent part
  float hbotedge;
  float htopedge;
  float hapothem;
  float hthickness;

  // transient part 
  float offset;
  float tan_a;
};

#endif // Geom_TrapezoidalPlaneBounds_H
