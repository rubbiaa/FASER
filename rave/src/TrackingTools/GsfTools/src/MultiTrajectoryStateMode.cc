#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateMode.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "TrackingTools/GsfTools/interface/MultiGaussianStateTransform.h"
#include "TrackingTools/GsfTools/interface/MultiGaussianState1D.h"
#include "TrackingTools/GsfTools/interface/GaussianSumUtilities1D.h"

bool
MultiTrajectoryStateMode::momentumFromModeCartesian (const TrajectoryStateOnSurface tsos,
						     GlobalVector& momentum) const
{
  //
  // clear result vector and check validity of the TSOS
  //
  momentum = GlobalVector(0.,0.,0.);
  if ( !tsos.isValid() ) {
    edm::LogInfo("MultiTrajectoryStateMode") << "Cannot calculate mode from invalid TSOS";
    return false;
  }
  //  
  // 1D mode computation for px, py and pz
  // 
  std::vector<TrajectoryStateOnSurface> components(tsos.components());
  unsigned int numb = components.size();
  // vectors of components in x, y and z
  std::vector<SingleGaussianState1D> pxStates; pxStates.reserve(numb);
  std::vector<SingleGaussianState1D> pyStates; pyStates.reserve(numb);
  std::vector<SingleGaussianState1D> pzStates; pzStates.reserve(numb);
  // iteration over components
  for ( std::vector<TrajectoryStateOnSurface>::const_iterator ic=components.begin();
	ic!=components.end(); ++ic ) {
    // extraction of parameters and variances
    GlobalVector mom(ic->globalMomentum());
    AlgebraicSymMatrix66 cov(ic->cartesianError().matrix());
    pxStates.push_back(SingleGaussianState1D(mom.x(),cov(3,3),ic->weight()));
    pyStates.push_back(SingleGaussianState1D(mom.y(),cov(4,4),ic->weight()));
    pzStates.push_back(SingleGaussianState1D(mom.z(),cov(5,5),ic->weight()));
  }
  //
  // transformation in 1D multi-states and creation of utility classes
  //
  MultiGaussianState1D pxState(pxStates);
  MultiGaussianState1D pyState(pyStates);
  MultiGaussianState1D pzState(pzStates);
  GaussianSumUtilities1D pxUtils(pxState);
  GaussianSumUtilities1D pyUtils(pyState);
  GaussianSumUtilities1D pzUtils(pzState);
  //
  // cartesian momentum vector from modes
  //
  momentum = GlobalVector(pxUtils.mode().mean(),pyUtils.mode().mean(),pzUtils.mode().mean());
  return true;
}

bool
MultiTrajectoryStateMode::momentumFromModeLocal (const TrajectoryStateOnSurface tsos,
						 GlobalVector& momentum) const
{
  //
  // clear result vector and check validity of the TSOS
  //
  momentum = GlobalVector(0.,0.,0.);
  if ( !tsos.isValid() ) {
    edm::LogInfo("MultiTrajectoryStateMode") << "Cannot calculate mode from invalid TSOS";
    return false;
  }
  //  
  // mode computation for local co-ordinates q/p, dx/dz, dy/dz
  //
  double qpMode(0);
  double dxdzMode(0);
  double dydzMode(0);
  //
  // first 3 elements of local parameters = q/p, dx/dz, dy/dz
  //
  for ( unsigned int iv=0; iv<3; ++iv ) {
    // extraction of multi-state using helper class
    MultiGaussianState1D state1D = MultiGaussianStateTransform::multiState1D(tsos,iv);
    GaussianSumUtilities1D utils(state1D);
    // mode (in case of failure: mean)
    double result = utils.mode().mean();
    if ( !utils.modeIsValid() )  result = utils.mean();
    if ( iv==0 )  qpMode = result;
    else if ( iv==1 )  dxdzMode = result;
    else  dydzMode = result;
  }
  // local momentum vector from dx/dz, dy/dz and q/p + sign of local pz
  LocalVector localP(dxdzMode,dydzMode,1.);
  localP *= tsos.localParameters().pzSign()/fabs(qpMode)
    /sqrt(dxdzMode*dxdzMode+dydzMode*dydzMode+1.);
  // conversion to global coordinates
  momentum = tsos.surface().toGlobal(localP);
  return true;
}

bool
MultiTrajectoryStateMode::momentumFromModePPhiEta (const TrajectoryStateOnSurface tsos,
						   GlobalVector& momentum) const
{
  //
  // clear result vector and check validity of the TSOS
  //
  momentum = GlobalVector(0.,0.,0.);
  if ( !tsos.isValid() ) {
    edm::LogInfo("MultiTrajectoryStateMode") << "Cannot calculate mode from invalid TSOS";
    return false;
  }
  //  
  // 1D mode computation for p, phi, eta
  // 
  std::vector<TrajectoryStateOnSurface> components(tsos.components());
  unsigned int numb = components.size();
  // vectors of components in p, phi and eta
  std::vector<SingleGaussianState1D> pStates; pStates.reserve(numb);
  std::vector<SingleGaussianState1D> phiStates; phiStates.reserve(numb);
  std::vector<SingleGaussianState1D> etaStates; etaStates.reserve(numb);
  // covariances in cartesian and p-phi-eta and jacobian
  AlgebraicMatrix33 jacobian;
  AlgebraicSymMatrix33 covCart;
  AlgebraicSymMatrix33 covPPhiEta;
  // iteration over components
  for ( std::vector<TrajectoryStateOnSurface>::const_iterator ic=components.begin();
	ic!=components.end(); ++ic ) {
    // parameters
    GlobalVector mom(ic->globalMomentum());
    double px = mom.x();
    double py = mom.y();
    double pz = mom.z();
    double p = mom.mag();
    double pt2 = mom.perp2();
    double phi = mom.phi();
    double eta = mom.eta();
    // jacobian
    jacobian(0,0) = px/p;
    jacobian(0,1) = py/p;
    jacobian(0,2) = pz/p;
    jacobian(1,0) = py/pt2;
    jacobian(1,1) = -px/pt2;
    jacobian(1,2) = 0;
    jacobian(2,0) = px*pz/(pt2*p);
    jacobian(2,1) = py*pz/(pt2*p);
    jacobian(2,2) = -1./p;
    // extraction of the momentum part from the 6x6 cartesian error matrix
    // and conversion to p-phi-eta
    covCart = ic->cartesianError().matrix().Sub<AlgebraicSymMatrix33>(3,3);
    covPPhiEta = ROOT::Math::Similarity(jacobian,covCart);
    pStates.push_back(SingleGaussianState1D(p,covPPhiEta(0,0),ic->weight()));
    phiStates.push_back(SingleGaussianState1D(phi,covPPhiEta(1,1),ic->weight()));
    etaStates.push_back(SingleGaussianState1D(eta,covPPhiEta(2,2),ic->weight()));
  }
  //
  // transformation in 1D multi-states and creation of utility classes
  //
  MultiGaussianState1D pState(pStates);
  MultiGaussianState1D phiState(phiStates);
  MultiGaussianState1D etaState(etaStates);
  GaussianSumUtilities1D pUtils(pState);
  GaussianSumUtilities1D phiUtils(phiState);
  GaussianSumUtilities1D etaUtils(etaState);
  //
  // parameters from mode (in case of failure: mean)
  //
  double p = pUtils.modeIsValid() ? pUtils.mode().mean() : pUtils.mean();
  double phi = phiUtils.modeIsValid() ? phiUtils.mode().mean() : phiUtils.mean();
  double eta = etaUtils.modeIsValid() ? etaUtils.mode().mean() : etaUtils.mean();
  double theta = 2*atan(exp(-eta));
  double tanth2 = exp(-eta);
  double pt = p*2*tanth2/(1+tanth2*tanth2);  // p*sin(theta)
  double pz = p*(1-tanth2*tanth2)/(1+tanth2*tanth2);  // p*cos(theta)
  // conversion to a cartesian momentum vector
  momentum = GlobalVector(pt*cos(phi),pt*sin(phi),pz);
  return true;
}

