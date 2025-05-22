#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "DataFormats/TrackReco/interface/Track.h" 
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h" 
#include "DataFormats/GeometrySurface/interface/Surface.h" 
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
// #include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"

using namespace SurfaceSideDefinition;

/*
PTrajectoryStateOnDet* 
TrajectoryStateTransform::persistentState( const TrajectoryStateOnSurface& ts,
					   unsigned int detid) const
{
  AlgebraicSymMatrix55 m = ts.localError().matrix();
  
  int dim = 5; /// should check if corresponds to m

  float localErrors[15];
  int k = 0;
  for (int i=0; i<dim; i++) {
    for (int j=0; j<=i; j++) {
      localErrors[k++] = m(i,j);
    }
  }
  int surfaceSide = static_cast<int>(ts.surfaceSide());

  return new PTrajectoryStateOnDet( ts.localParameters(),
				    localErrors, detid,
				    surfaceSide);
}

TrajectoryStateOnSurface 
TrajectoryStateTransform::transientState( const PTrajectoryStateOnDet& ts,
					  const Surface* surface,
					  const MagneticField* field) const
{
  int dim = 5;
  AlgebraicSymMatrix55 m;
  const std::vector<float> &errs = ts.errorMatrix();  
  int k = 0;
  for (int i=0; i<dim; i++) {
    for (int j=0; j<=i; j++) {
      m(i,j) = errs[k++];       // NOTE: here we do a cast float => double.     
    }
  }
  

  return TrajectoryStateOnSurface( ts.parameters(),
				   LocalTrajectoryError(m), 
				   *surface, field,
				   static_cast<SurfaceSide>(ts.surfaceSide()));

}
*/

FreeTrajectoryState TrajectoryStateTransform::initialFreeState( const reco::Track& tk,
							      const MagneticField* field) const
{
  Basic3DVector<float> pos( tk.vertex());
  GlobalPoint gpos( pos);
  Basic3DVector<float> mom( tk.momentum());
  GlobalVector gmom( mom);
  GlobalTrajectoryParameters par( gpos, gmom, tk.charge(), field);
  CurvilinearTrajectoryError err( tk.covariance());
  return FreeTrajectoryState( par, err);
}

/*
FreeTrajectoryState TrajectoryStateTransform::innerFreeState( const reco::Track& tk,
							      const MagneticField* field) const
{
  Basic3DVector<float> pos( tk.innerPosition());
  GlobalPoint gpos( pos);
  Basic3DVector<float> mom( tk.innerMomentum());
  GlobalVector gmom( mom);
  GlobalTrajectoryParameters par( gpos, gmom, tk.charge(), field);
  CurvilinearTrajectoryError err( tk.extra()->innerStateCovariance());
  return FreeTrajectoryState( par, err);
}


FreeTrajectoryState TrajectoryStateTransform::outerFreeState( const reco::Track& tk,
							      const MagneticField* field) const
{
  Basic3DVector<float> pos( tk.outerPosition());
  GlobalPoint gpos( pos);
  Basic3DVector<float> mom( tk.outerMomentum());
  GlobalVector gmom( mom);
  GlobalTrajectoryParameters par( gpos, gmom, tk.charge(), field);
  CurvilinearTrajectoryError err( tk.extra()->outerStateCovariance());
  return FreeTrajectoryState( par, err);
}
*/

/*
#include <iostream>
using namespace std;

TrajectoryStateOnSurface TrajectoryStateTransform::innerStateOnSurface( const reco::Track& tk, 
									const TrackingGeometry& geom,
									const MagneticField* field) const
{
  cout << "[TrajectoryStateTransform] aiyye need implementation in " << __LINE__ << endl;
  return TrajectoryStateOnSurface();
  const Surface& surface = geom.idToDet( DetId( tk.extra()->innerDetId()))->surface();
  return TrajectoryStateOnSurface( innerFreeState( tk, field), surface);
}

TrajectoryStateOnSurface TrajectoryStateTransform::outerStateOnSurface( const reco::Track& tk, 
									const TrackingGeometry& geom,
									const MagneticField* field) const
{
  cout << "[TrajectoryStateTransform] aiyye need implementation in " << __LINE__ << endl;
  return TrajectoryStateOnSurface();
  const Surface& surface = geom.idToDet( DetId( tk.extra()->outerDetId()))->surface();
  return TrajectoryStateOnSurface( outerFreeState( tk, field), surface);
}
*/
