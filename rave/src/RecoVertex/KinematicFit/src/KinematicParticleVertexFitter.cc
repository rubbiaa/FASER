#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
// #include "Vertex/LinearizationPointFinders/interface/LMSLinearizationPointFinder.h"
#include "RecoVertex/KinematicFit/interface/FinalTreeBuilder.h"
#include "RecoVertex/VertexTools/interface/SequentialVertexSmoother.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanSmoothedVertexChi2Estimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanTrackToTrackCovCalculator.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "RecoVertex/LinearizationPointFinders/interface/DefaultLinearizationPointFinder.h"
#include "RecoVertex/VertexTools/interface/SequentialVertexFitter.h"
#include "DataFormats/CLHEP/interface/Migration.h"

using namespace std;

KinematicParticleVertexFitter::KinematicParticleVertexFitter()
{ 

  pointFinder =  new DefaultLinearizationPointFinder();
  vFactory = new VertexTrackFactory<6>();

  KalmanVertexTrackUpdator<6> vtu;
  KalmanSmoothedVertexChi2Estimator<6> vse;
  KalmanTrackToTrackCovCalculator<6> covCalc;
  SequentialVertexSmoother<6> smoother(vtu, vse, covCalc);
  edm::ParameterSet pSet;
  pSet.addParameter<double>("maxDistance", 0.01);
  pSet.addParameter<int>("maxNbrOfIterations", 10); //10
  fitter 
    = new SequentialVertexFitter<6>(pSet, *pointFinder, KalmanVertexUpdator<6>(),
				 smoother, ParticleKinematicLinearizedTrackStateFactory());
}

KinematicParticleVertexFitter::~KinematicParticleVertexFitter()
{
 delete vFactory;
 delete pointFinder;
 delete fitter;
}
 
RefCountedKinematicTree KinematicParticleVertexFitter::fit(
    vector<RefCountedKinematicParticle> particles, bool usebeamspot,
    const reco::BeamSpot & spot ) const
{
 typedef ReferenceCountingPointer<VertexTrack<6> > RefCountedVertexTrack;
//sorting the input 
 if(particles.size()<2) throw VertexException("KinematicParticleVertexFitter::input states are less than 2"); 
 InputSort iSort;
 pair<vector<RefCountedKinematicParticle>, vector<FreeTrajectoryState> > input = iSort.sort(particles);
 vector<RefCountedKinematicParticle> & newPart = input.first;
 vector<FreeTrajectoryState> & freeStates = input.second;

 GlobalPoint linPoint = pointFinder->getLinearizationPoint(freeStates);
  
// cout<<"Linearization point found"<<endl; 
 
//making initial veretx seed with lin point as position and a fake error
 AlgebraicSymMatrix33 we;
 we(0,0)=we(1,1)=we(2,2) = 10000.;
 GlobalError error(we);
 VertexState state(linPoint, error);
 
//vector of Vertex Tracks to fit
 vector<RefCountedVertexTrack> ttf; 
 for(vector<RefCountedKinematicParticle>::const_iterator i = newPart.begin();i != newPart.end();i++)
 {ttf.push_back(vFactory->vertexTrack((*i)->particleLinearizedTrackState(linPoint),state,1.));}

// //debugging code to check neutrals: 
//  for(vector<RefCountedVertexTrack>::const_iterator i = ttf.begin(); i!=ttf.end(); i++)
//  {
// //   cout<<"predicted state momentum error"<<(*i)->linearizedTrack()->predictedStateMomentumError()<<endl;
// //  cout<<"Momentum jacobian"<<(*i)->linearizedTrack()->momentumJacobian() <<endl;
//  //  cout<<"predicted state momentum "<<(*i)->linearizedTrack()->predictedStateMomentum()<<endl;
// //   cout<<"constant term"<<(*i)->linearizedTrack()->constantTerm()<<endl;
// 
//  }
//
 /* cout << "[KinematicParticleVertexFitter] now really use beamspot: " 
      << usebeamspot << endl; */
 CachingVertex<6> vtx;
 if ( usebeamspot )
 {
   // cout << "[KinematicParticleVertexFitter] spot=" << spot << endl;
   // vtx = fitter->vertex(ttf ); 
   vtx = fitter->vertex(ttf, spot ); 
 } else {
   vtx = fitter->vertex(ttf); 
 }

 FinalTreeBuilder tBuilder;
 return tBuilder.buildTree(vtx, newPart); 
}
