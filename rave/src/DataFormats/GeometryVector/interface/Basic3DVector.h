#ifndef Basic3DVector_H
#define Basic3DVector_H

#include "DataFormats/GeometryVector/interface/Basic2DVector.h"
#include "DataFormats/GeometryVector/interface/Theta.h"
#include "DataFormats/GeometryVector/interface/Phi.h"
#include "DataFormats/GeometryVector/interface/PreciseFloatType.h"
#include "DataFormats/GeometryVector/interface/CoordinateSets.h"
#include <iosfwd>
#include <iostream>
#include <cmath>

template < typename T> 
class Basic3DVector {
public:

  typedef T                                   ScalarType;
  typedef Geom::Cylindrical2Cartesian<T>      Cylindrical;
  typedef Geom::Spherical2Cartesian<T>        Spherical;
  typedef Spherical                           Polar; // synonym
    
  /** default constructor uses default constructor of T to initialize the 
   *  components. For built-in floating-point types this means initialization 
   * to zero
   */
  Basic3DVector() : theX( T()), theY( T()), theZ( T()) {}

  /// Copy constructor from same type. Should not be needed but for gcc bug 12685
  Basic3DVector( const Basic3DVector & p) : 
    theX(p.x()), theY(p.y()), theZ(p.z()) {}

  /// constructor from 2D vector (X and Y from 2D vector, z set to zero)
  Basic3DVector( const Basic2DVector<T> & p) : 
    theX(p.x()), theY(p.y()), theZ(0) {}

  /** Explicit constructor from other (possibly unrelated) vector classes 
   *  The only constraint on the argument type is that it has methods
   *  x(), y() and z(), and that these methods return a type convertible to T.
   *  Examples of use are
   *   <BR> construction from a Basic3DVector with different precision
   *   <BR> construction from a Hep3Vector
   *   <BR> construction from a coordinate system converter 
   */
  template <class OtherPoint> 
  explicit Basic3DVector( const OtherPoint& p) : 
    theX(p.x()), theY(p.y()), theZ(p.z()) {}

  /// construct from cartesian coordinates
  Basic3DVector( const T& x, const T& y, const T& z) : 
    theX(x), theY(y), theZ(z) {}

  /** Deprecated construct from polar coordinates, use 
   *  <BR> Basic3DVector<T>( Basic3DVector<T>::Polar( theta, phi, r))
   *  instead. 
   */
  template <typename U>
  Basic3DVector( const Geom::Theta<U>& theta, 
		 const Geom::Phi<U>& phi, const T& r) {
    Polar p( theta.value(), phi.value(), r);
    theX = p.x(); theY = p.y(); theZ = p.z();
  }

  /// Cartesian x coordinate
  T x() const { return theX;}

  /// Cartesian y coordinate
  T y() const { return theY;}

  /// Cartesian z coordinate
  T z() const { return theZ;}

  /// The vector magnitude squared. Equivalent to vec.dot(vec)
  T mag2() const { return theX*theX + theY*theY + theZ*theZ;}

  /// The vector magnitude. Equivalent to sqrt(vec.mag2())
  T mag() const  { return sqrt( mag2());}

  /// Squared magnitude of transverse component 
  T perp2() const { return x()*x() + y()*y();}

  /// Magnitude of transverse component 
  T perp() const { return sqrt( perp2());}

  /// Another name for perp()
  T transverse() const { return perp();}

  /** Azimuthal angle. The value is returned in radians, in the range (-pi,pi].
   *  Same precision as the system atan2(x,y) function.
   *  The return type is Geom::Phi<T>, see it's documentation.
   */ 
  Geom::Phi<T> phi() const {return Geom::Phi<T>(atan2(theY,theX));}

  /** Polar angle. The value is returned in radians, in the range [0,pi]
   *  Same precision as the system atan2(x,y) function.
   *  The return type is Geom::Phi<T>, see it's documentation.
   */ 
  Geom::Theta<T> theta() const {return Geom::Theta<T>(atan2(perp(),z()));}

  /** Pseudorapidity. 
   *  Does not check for zero transverse component; in this case the behavior 
   *  is as for divide-by zero, i.e. system-dependent.
   */
  // T eta() const { return -log( tan( theta()/2.));} 
  T eta() const { T x(z()/perp()); return log(x+sqrt(x*x+1));} // faster 

  /** Unit vector parallel to this.
   *  If mag() is zero, a zero vector is returned.
   */
  Basic3DVector unit() const {
    T my_mag = mag();
    return my_mag == 0 ? *this : *this / my_mag;
  }

  /** Operator += with a Basic3DVector of possibly different precision.
   */
  template <class U> 
  Basic3DVector& operator+= ( const Basic3DVector<U>& p) {
    theX += p.x();
    theY += p.y();
    theZ += p.z();
    return *this;
  } 

  /** Operator -= with a Basic3DVector of possibly different precision.
   */
  template <class U> 
  Basic3DVector& operator-= ( const Basic3DVector<U>& p) {
    theX -= p.x();
    theY -= p.y();
    theZ -= p.z();
    return *this;
  } 

  /// Unary minus, returns a vector with components (-x(),-y(),-z())
  Basic3DVector operator-() const { return Basic3DVector(-x(),-y(),-z());}

  /// Scaling by a scalar value (multiplication)
  Basic3DVector& operator*= ( const T& t) {
    theX *= t;
    theY *= t;
    theZ *= t;
    return *this;
  } 

  /// Scaling by a scalar value (division)
  Basic3DVector& operator/= ( const T& t) {
    theX /= t;
    theY /= t;
    theZ /= t;
    return *this;
  } 

  /// Scalar product, or "dot" product, with a vector of same type.
  T dot( const Basic3DVector& v) const { 
    return x()*v.x() + y()*v.y() + z()*v.z();
  }

  /** Scalar (or dot) product with a vector of different precision.
   *  The product is computed without loss of precision. The type
   *  of the returned scalar is the more precise of the scalar types 
   *  of the two vectors.
   */
  template <class U> 
  typename PreciseFloatType<T,U>::Type dot( const Basic3DVector<U>& v) const { 
    return x()*v.x() + y()*v.y() + z()*v.z();
  }

  /// Vector product, or "cross" product, with a vector of same type.
  Basic3DVector cross( const Basic3DVector& v) const {
    return Basic3DVector( y()*v.z() - v.y()*z(), 
			  z()*v.x() - v.z()*x(), 
			  x()*v.y() - v.x()*y());
  }

  /** Vector (or cross) product with a vector of different precision.
   *  The product is computed without loss of precision. The type
   *  of the returned vector is the more precise of the types 
   *  of the two vectors.   
   */
  template <class U> 
  Basic3DVector<typename PreciseFloatType<T,U>::Type> 
  cross( const Basic3DVector<U>& v) const {
    return Basic3DVector<typename PreciseFloatType<T,U>::Type>( y()*v.z() - v.y()*z(), 
								z()*v.x() - v.z()*x(), 
								x()*v.y() - v.x()*y());
  }

private:
  T theX;
  T theY;
  T theZ;
};

/// simple text output to standard streams
template <class T>
inline std::ostream & operator<<( std::ostream& s, const Basic3DVector<T>& v) {
  return s << " (" << v.x() << "," << v.y() << "," << v.z() << ") ";
} 

/// vector sum and subtraction of vectors of possibly different precision
template <class T, class U>
inline Basic3DVector<typename PreciseFloatType<T,U>::Type>
operator+( const Basic3DVector<T>& a, const Basic3DVector<U>& b) {
  typedef Basic3DVector<typename PreciseFloatType<T,U>::Type> RT;
  return RT(a.x()+b.x(), a.y()+b.y(), a.z()+b.z());
}

template <class T, class U>
inline Basic3DVector<typename PreciseFloatType<T,U>::Type>
operator-( const Basic3DVector<T>& a, const Basic3DVector<U>& b) {
  typedef Basic3DVector<typename PreciseFloatType<T,U>::Type> RT;
  return RT(a.x()-b.x(), a.y()-b.y(), a.z()-b.z());
}

/// scalar product of vectors of same precision
template <class T>
inline T operator*( const Basic3DVector<T>& v1, const Basic3DVector<T>& v2) {
  return v1.dot(v2);
}

/// scalar product of vectors of different precision
template <class T, class U>
inline typename PreciseFloatType<T,U>::Type operator*( const Basic3DVector<T>& v1, 
						       const Basic3DVector<U>& v2) {
  return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
}

/** Multiplication by scalar, does not change the precision of the vector.
 *  The return type is the same as the type of the vector argument.
 */
template <class T, class Scalar>
inline Basic3DVector<T> operator*( const Basic3DVector<T>& v, const Scalar& s) {
  T t = static_cast<T>(s);
  return Basic3DVector<T>(v.x()*t, v.y()*t, v.z()*t);
}

/// Same as operator*( Vector, Scalar)
template <class T, class Scalar>
inline Basic3DVector<T> operator*( const Scalar& s, const Basic3DVector<T>& v) {
  T t = static_cast<T>(s);
  return Basic3DVector<T>(v.x()*t, v.y()*t, v.z()*t);
}

/** Division by scalar, does not change the precision of the vector.
 *  The return type is the same as the type of the vector argument.
 */
template <class T, class Scalar>
inline Basic3DVector<T> operator/( const Basic3DVector<T>& v, const Scalar& s) {
  T t = static_cast<T>(s);
  return Basic3DVector<T>(v.x()/t, v.y()/t, v.z()/t);
}

#endif // Basic3DVector_H


