#ifndef _TRACKER_FREETRAJECTORYSTATE_H_
#define _TRACKER_FREETRAJECTORYSTATE_H_

// base trajectory state class
// track parameters and error covariance matrix

#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "DataFormats/TrajectoryState/interface/TrackCharge.h"
#include "TrackingTools/TrajectoryParametrization/interface/TrajectoryStateExceptions.h"

#include <iosfwd>

/** A 6-dimensional state vector of a helix given at some point in 
 *  space along the helix, and the associated error matrix.
 *  The error can be obtained in two different parametrizations:
 *  CurvilinearTrajectoryError and CartesianTrajectoryError
 *  (see descriptions of corresponding classes).
 *  The FreeTrajectoryState can be instantiated with either of these
 *  error parametrisations; it converts from one to the other internally.
 */

class FreeTrajectoryState {
public:
// construst
//default constructor - needed for Persistency

  FreeTrajectoryState():
    theCartesianErrorValid(false),
    theCurvilinearErrorValid(false), theTrackId(-1) {};

  FreeTrajectoryState(const GlobalTrajectoryParameters& aGlobalParameters) :
    theGlobalParameters(aGlobalParameters),
    theCartesianErrorValid(false),
    theCurvilinearErrorValid(false),
    theTrackId(-1)
  {
  }

  FreeTrajectoryState(const GlobalPoint& aX,
                      const GlobalVector& aP,
                      TrackCharge aCharge, 
                      const MagneticField* fieldProvider) :
    theGlobalParameters(aX, aP, aCharge, fieldProvider),
    theCartesianErrorValid(false),
    theCurvilinearErrorValid(false), theTrackId(-1)
  {
  }

  FreeTrajectoryState(const GlobalTrajectoryParameters& aGlobalParameters,
                      const CartesianTrajectoryError& aCartesianError) :
    theGlobalParameters(aGlobalParameters),
    theCartesianError(aCartesianError),
    theCartesianErrorValid(true),
    theCurvilinearErrorValid(false),
    theTrackId(-1)
  {
  }
  FreeTrajectoryState(const GlobalTrajectoryParameters& aGlobalParameters,
                      const CurvilinearTrajectoryError& aCurvilinearError) :
    theGlobalParameters(aGlobalParameters),
    theCurvilinearError(aCurvilinearError),
    theCartesianErrorValid(false),
    theCurvilinearErrorValid(true),
    theTrackId(-1)
  {
  }
  FreeTrajectoryState(const GlobalTrajectoryParameters& aGlobalParameters,
                      const CartesianTrajectoryError& aCartesianError,
                      const CurvilinearTrajectoryError& aCurvilinearError) :
    theGlobalParameters(aGlobalParameters),
    theCartesianError(aCartesianError),
    theCurvilinearError(aCurvilinearError),
    theCartesianErrorValid(true),
    theCurvilinearErrorValid(true),
    theTrackId(-1)
  {
  }

  int trackId() const {
    return theTrackId;
  }

  void setTrackId( int id ) {
    theTrackId=id;
  }
// access
// propagate access to parameters
  GlobalPoint position() const {
    return theGlobalParameters.position();
  }
  GlobalVector momentum() const {
    return theGlobalParameters.momentum();
  }
  TrackCharge charge() const {
    return theGlobalParameters.charge();
  }
  double signedInverseMomentum() const {
    return theGlobalParameters.signedInverseMomentum();
  }
  double transverseCurvature() const {
    return theGlobalParameters.transverseCurvature();
  }
// direct access
  bool hasCartesianError() const {return theCartesianErrorValid;}
  bool hasCurvilinearError() const {return theCurvilinearErrorValid;}
  bool hasError() const {
    return theCurvilinearErrorValid || theCartesianErrorValid;
  }
  const GlobalTrajectoryParameters& parameters() const {
    return theGlobalParameters;
  }
  const CartesianTrajectoryError& cartesianError() const {
    if (!hasError()) throw TrajectoryStateException(
      "FreeTrajectoryState: attempt to access errors when none available");
    if (!theCartesianErrorValid)
      createCartesianError();
    return theCartesianError;
  }
  const CurvilinearTrajectoryError& curvilinearError() const {
    if (!hasError()) throw TrajectoryStateException(
      "FreeTrajectoryState: attempt to access errors when none available");
    if (!theCurvilinearErrorValid)
      createCurvilinearError();
    return theCurvilinearError;
  }
  void rescaleError(double factor) {
    bool zeroField = parameters().magneticFieldInInverseGeV(GlobalPoint(0,0,0)).mag2()==0;
    if (zeroField) {
      if (theCartesianErrorValid){
	if (!theCurvilinearErrorValid) createCurvilinearError();
	theCurvilinearError.zeroFieldScaling(factor*factor);
	createCartesianError();
      }else
	if (theCurvilinearErrorValid) theCurvilinearError.zeroFieldScaling(factor*factor);
    } else{
      if (theCartesianErrorValid){
	theCartesianError *= (factor*factor);
      }
      if (theCurvilinearErrorValid){
	theCurvilinearError *= (factor*factor);
      }
    }
  }

  void setCartesianError(const CartesianTrajectoryError &err) {
        theCartesianError = err; theCartesianErrorValid = true;
  }
  void setCartesianError(const AlgebraicSymMatrix66 &err) {
        theCartesianError = CartesianTrajectoryError(err); theCartesianErrorValid = true;
  }
  void setCurvilinearError(const CurvilinearTrajectoryError &err) {
        theCurvilinearError = err; theCurvilinearErrorValid = true;
  }
  void setCurvilinearError(const AlgebraicSymMatrix55 &err) {
        theCurvilinearError = CurvilinearTrajectoryError(err); theCurvilinearErrorValid = true;
  }
// properties
  bool canReach(double radius) const;
private:
// convert curvilinear errors to cartesian
  void createCartesianError() const;
// convert cartesian errors to curvilinear
  void createCurvilinearError() const;

private:
  GlobalTrajectoryParameters  theGlobalParameters;
  CartesianTrajectoryError    theCartesianError;
  CurvilinearTrajectoryError  theCurvilinearError;
  bool                        theCartesianErrorValid;
  bool                        theCurvilinearErrorValid;
  int                         theTrackId; // which track does this belong to
};

std::ostream& operator<<(std::ostream& os, const FreeTrajectoryState& fts);

#endif
