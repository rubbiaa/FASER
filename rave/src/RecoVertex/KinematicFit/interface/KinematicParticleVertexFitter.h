#ifndef KinematicParticleVertexFitter_H
#define KinematicParticleVertexFitter_H

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "RecoVertex/VertexTools/interface/LinearizationPointFinder.h"
#include "RecoVertex/VertexTools/interface/VertexTrackFactory.h"
#include "RecoVertex/KinematicFit/interface/InputSort.h"
#include "RecoVertex/VertexPrimitives/interface/VertexFitter.h"

/**
 * Class creating a kinematic particle out
 * of set of daughter particles. Daughter
 * particles are supposed to come from a 
 * common vertex. Arbitrary VerexFitter
 * can be used to fit the common vertex 
 * and refit the daughter particles with the
 * knowledge of vertex. The Kinematic Vertex is also
 * created and the resulting KinematicParticle points
 * on it.
 * 
 * Kirill Prokofiev, December 2002
 */

class KinematicParticleVertexFitter
{

public:
 
/**
 * Constructor with LMSLinearizationPointFinder used as default.
 *
 */ 
 KinematicParticleVertexFitter();
 
 
/*
 * Constructor with the LinearizationPointFinder
 * Linearization point finder should have an 
 * ability to find point out of set of FreeTrajectoryStates
 * LMSLinearizationPointFinder is used as default.
 */
//  KinematicParticleVertexFitter(const LinearizationPointFinder& finder);
 
 ~KinematicParticleVertexFitter();
 
/**
 * Fit method taking set of particles, fitting them to the
 * common vertex and creating tree out of them.
 * Input particles can belong to kinmaticTrees.
 * In such a case it should be TOP particle of
 * corresponding tree.
 */  
 RefCountedKinematicTree  fit(vector<RefCountedKinematicParticle> particles,
      bool usebeamspot=false, const reco::BeamSpot & s = reco::BeamSpot()  ) const;
 
private:
		    
//widely used common tools 
  //SequentialKinematicVertexFitter * fitter; use the default KalmanFilter instead!
  VertexFitter<6> * fitter;
  LinearizationPointFinder * pointFinder; 
  VertexTrackFactory<6> * vFactory;
};
#endif
