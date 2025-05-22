#include "RecoVertex/TrimmedKalmanVertexFinder/interface/ConfigurableTrimmedVertexFinder.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

using namespace reco;

namespace {
 struct CompareTwoTracks {
    int operator() ( const reco::TransientTrack & a,
                     const reco::TransientTrack & b ) {
            return ( a.impactPointState().globalMomentum().perp() >                                                 b.impactPointState().globalMomentum().perp() ) ;
    };
  };
}

ConfigurableTrimmedVertexFinder::ConfigurableTrimmedVertexFinder(
  const VertexFitter<5> * vf,
  const VertexUpdator<5> * vu,
  const VertexTrackCompatibilityEstimator<5> * ve)
  : theClusterFinder(vf, vu, ve), theVtxFitProbCut(0.01),
    theTrackCompatibilityToPV(0.05), theTrackCompatibilityToSV(0.01),
    theMaxNbOfVertices(0)
{
  // default pt cut is 1.5 GeV
  theFilter.setPtCut(1.5);
}

void ConfigurableTrimmedVertexFinder::setParameters ( const edm::ParameterSet & s )
{
  theFilter.setPtCut(s.getParameter<double>("ptCut"));
  theTrackCompatibilityToPV = s.getParameter<double>("trackCompatibilityToPVcut");
  theTrackCompatibilityToSV = s.getParameter<double>("trackCompatibilityToSVcut");
  theVtxFitProbCut = s.getParameter<double>("vtxFitProbCut");
  theMaxNbOfVertices =  s.getParameter<int>("maxNbOfVertices");
}


vector<TransientVertex> ConfigurableTrimmedVertexFinder::vertices(
  const vector<TransientTrack> & ptracks) const
{
  vector<TransientTrack> tracks = ptracks;
  sort ( tracks.begin(), tracks.end(), CompareTwoTracks() );
  vector<TransientTrack> remaining;

  return vertices(tracks, remaining, reco::BeamSpot(), false );

}

vector<TransientVertex> ConfigurableTrimmedVertexFinder::vertices(
  const vector<TransientTrack> & ptracks, const reco::BeamSpot & spot ) const
{
  vector<TransientTrack> remaining;
  vector<TransientTrack> tracks = ptracks;
  sort ( tracks.begin(), tracks.end(), CompareTwoTracks() );
  return vertices ( tracks, remaining, spot, true );
}

vector<TransientVertex> ConfigurableTrimmedVertexFinder::vertices(
  const vector<TransientTrack> & ptracks, vector<TransientTrack> & unused,
  const reco::BeamSpot & spot, bool use_spot ) const
{
  vector<TransientTrack> tracks = ptracks;
  sort ( tracks.begin(), tracks.end(), CompareTwoTracks() );
  resetEvent(tracks);
  analyseInputTracks(tracks);

  vector<TransientTrack> filtered;
  for (vector<TransientTrack>::const_iterator it = tracks.begin();
       it != tracks.end(); it++) {
    if (theFilter(*it)) {
      filtered.push_back(*it);
    }
    else {
      unused.push_back(*it);
    }
  }

  vector<TransientVertex> all = vertexCandidates(filtered, unused,
      spot, use_spot );

  analyseVertexCandidates(all);

  vector<TransientVertex> sel = clean(all);

  analyseFoundVertices(sel);

  return sel;

}


vector<TransientVertex> ConfigurableTrimmedVertexFinder::vertexCandidates(
  const vector<TransientTrack> & tracks, vector<TransientTrack> & unused,
  const reco::BeamSpot & spot, bool use_spot ) const
{

  vector<TransientVertex> cand;

  vector<TransientTrack> remain = tracks;

  while (true) {

    float tkCompCut = (cand.size() == 0 ?
		       theTrackCompatibilityToPV
		       : theTrackCompatibilityToSV);

    //    cout << "PVR:compat cut " << tkCompCut << endl;
    theClusterFinder.setTrackCompatibilityCut(tkCompCut);
    //    cout << "PVCF:compat cut after setting "
    //	 << theClusterFinder.trackCompatibilityCut() << endl;

    vector<TransientVertex> newVertices;
    if ( cand.size() == 0 && use_spot )
    {
      newVertices = theClusterFinder.vertices(remain, spot );
    } else {
      newVertices = theClusterFinder.vertices(remain);
    }
    if (newVertices.empty()) break;

    analyseClusterFinder(newVertices, remain);

    for (vector<TransientVertex>::const_iterator iv = newVertices.begin();
         iv != newVertices.end(); iv++) {
      if ( iv->originalTracks().size() > 1 ) {
        cand.push_back(*iv);
      }
      else {
        // candidate has too few tracks - get them back into the vector
        for ( vector< TransientTrack >::const_iterator trk
		= iv->originalTracks().begin();
              trk != iv->originalTracks().end(); ++trk ) {
          unused.push_back ( *trk );
        }
      }
    }

    // when max number of vertices reached, stop
    if (theMaxNbOfVertices != 0) {
      if (cand.size() >= (unsigned int) theMaxNbOfVertices) break;
    }
  }

  for (vector<TransientTrack>::const_iterator it = remain.begin();
       it != remain.end(); it++) {
    unused.push_back(*it);
  }

  return cand;
}


vector<TransientVertex>
ConfigurableTrimmedVertexFinder::clean(const vector<TransientVertex> & candidates) const
{
  vector<TransientVertex> sel;
  for (vector<TransientVertex>::const_iterator i = candidates.begin();
       i != candidates.end(); i++) {

    if (ChiSquaredProbability((*i).totalChiSquared(), (*i).degreesOfFreedom())
	> theVtxFitProbCut) { sel.push_back(*i); }
  }

  return sel;
}
