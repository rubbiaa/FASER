#include "RecoVertex/VertexTools/interface/PerigeeRefittedTrackState.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "RaveTools/Converters/interface/PropagatorSingleton.h"

AlgebraicVector3 PerigeeRefittedTrackState::momentumVector() const
{
 if (!momentumVectorAvailable) {
    momentumAtVertex[0] = theState.perigeeParameters().vector()[0];
    momentumAtVertex[1] = theState.perigeeParameters().theta();
    momentumAtVertex[2] = theState.perigeeParameters().phi();
    momentumVectorAvailable = true;
  }
  return momentumAtVertex;
}

std::vector< PerigeeRefittedTrackState::RefCountedRefittedTrackState > 
PerigeeRefittedTrackState::components() const
{
  std::vector<RefCountedRefittedTrackState> result; result.reserve(1);
  result.push_back(RefCountedRefittedTrackState( 
  				const_cast<PerigeeRefittedTrackState*>(this)));
  return result;
}

PerigeeRefittedTrackState::RefCountedRefittedTrackState
PerigeeRefittedTrackState::stateWithNewWeight (const double newWeight) const
{
  return RefCountedRefittedTrackState(
  		new PerigeeRefittedTrackState(theState, newWeight) );
}

TrajectoryStateOnSurface
PerigeeRefittedTrackState::trajectoryStateOnSurface(const Surface & surface) const
{
  // AnalyticalPropagator thePropagator(&(theState.theState().parameters().magneticField()), anyDirection);
  TrajectoryStateOnSurface tsos = PropagatorSingleton::Instance()->propagator()->propagate(freeTrajectoryState(), surface);
  return TrajectoryStateOnSurface (tsos.globalParameters(),
  	tsos.curvilinearError(), surface ,weight()) ;
} 

TrajectoryStateOnSurface
PerigeeRefittedTrackState::trajectoryStateOnSurface(const Surface & surface,
				const Propagator & propagator) const
{
  std::auto_ptr<Propagator> thePropagator( propagator.clone());
  thePropagator->setPropagationDirection(anyDirection);

  TrajectoryStateOnSurface tsos = thePropagator->propagate(freeTrajectoryState(), surface);
  return TrajectoryStateOnSurface (tsos.globalParameters(),
  	tsos.curvilinearError(), surface ,weight()) ;
}

reco::TransientTrack PerigeeRefittedTrackState::transientTrack() const
{
  TransientTrackFromFTSFactory factory;
  return factory.build(freeTrajectoryState());
}
