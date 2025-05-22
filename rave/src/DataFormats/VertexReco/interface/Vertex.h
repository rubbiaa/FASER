#ifndef VertexReco_Vertex_h
#define VertexReco_Vertex_h
/** \class reco::Vertex
 *  
 * A reconstructed Vertex providing position, error, chi2, ndof 
 * and reconstrudted tracks.
 * The vertex can be valid, fake, or invalid.
 * A valid vertex is one which has been obtained from a vertex fit of tracks, 
 * and all data is meaningful
 * A fake vertex is a vertex which was not made out of a proper fit with
 * tracks, but still has a position and error (chi2 and ndof are null). 
 * For a primary vertex, it could simply be the beam line.
 * A fake vertex is considered valid.
 * An invalid vertex has no meaningful data.
 *
 * \author Luca Lista, INFN
 *
 * \version $Id: Vertex.h,v 1.33 2008/11/06 17:06:23 elmer Exp $
 *
 */
#include <Rtypes.h>
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/RefToBase.h" 

namespace reco {

  class Track;

  class Vertex {
  public:
    /// The iteratator for the vector<TrackRef>
    typedef std::vector<TrackBaseRef >::const_iterator trackRef_iterator;
    /// point in the space
    typedef math::XYZPoint Point;
    /// error matrix dimension
    enum { dimension = 3 };
    /// covariance error matrix (3x3)
    typedef math::Error<dimension>::type Error;
    /// covariance error matrix (3x3)
    typedef math::Error<dimension>::type CovarianceMatrix;
    /// matix size
    enum { size = dimension * ( dimension + 1 ) / 2 };
    /// index type
    typedef unsigned int index;
    /// default constructor - The vertex will not be valid. Position, error,
    /// chi2, ndof will have random entries, and the vectors of tracks will be empty
    /// Use the isValid method to check that your vertex is valid. 
    Vertex() { validity_ = false;}
    /// Constructor for a fake vertex.
    Vertex( const Point &, const Error &);
    /// constructor for a valid vertex, with all data
    Vertex( const Point &, const Error &, double chi2, double ndof, size_t size );
    /// Tells whether the vertex is valid.
    bool isValid() const {return validity_;}
    /// Tells whether a Vertex is fake, i.e. not a vertex made out of a proper
    /// fit with tracks.
    /// For a primary vertex, it could simply be the beam line.
    bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    /// add a reference to a Track
    void add( const TrackBaseRef & r, float w=1.0 );
    /// add the original a Track(reference) and the smoothed Track
    void add( const TrackBaseRef & r, const Track & refTrack, float w=1.0 );
    void removeTracks();
    ///returns the weight with which a Track has contributed to the vertex-fit.
    float trackWeight ( const TrackBaseRef & r ) const;
    ///returns the weight with which a Track has contributed to the vertex-fit.
    float trackWeight ( const TrackRef & r ) const;
    /// first iterator over tracks
    trackRef_iterator tracks_begin() const;
    /// last iterator over tracks
    trackRef_iterator tracks_end() const;
    /// number of tracks
    size_t tracksSize() const;
    /// chi-squares
    double chi2() const { return chi2_; }
    /** Number of degrees of freedom
     *  Meant to be Double32_t for soft-assignment fitters: 
     *  tracks may contribute to the vertex with fractional weights.
     *  The ndof is then = to the sum of the track weights.
     *  see e.g. CMS NOTE-2006/032, CMS NOTE-2004/002
     */
    double ndof() const { return ndof_; }
    /// chi-squared divided by n.d.o.f.
    double normalizedChi2() const { return chi2_ / ndof_; }
    /// position 
    const Point & position() const { return position_; }
    /// x coordinate 
    double x() const { return position_.X(); }
    /// y coordinate 
    double y() const { return position_.Y(); }
    /// y coordinate 
    double z() const { return position_.Z(); }
    /// error on x
    double xError() const { return sqrt( covariance(0, 0) ); }
    /// error on y
    double yError() const { return sqrt( covariance(1, 1) ); }
    /// error on z
    double zError() const { return sqrt( covariance(2, 2) ); }
    /// (i, j)-th element of error matrix, i, j = 0, ... 2
    // Note that:
    //   double error( int i, int j ) const 
    // is OBSOLETE, use covariance(i, j)
    double covariance( int i, int j ) const { 
      return covariance_[ idx( i, j ) ]; 
    }
    /// return SMatrix
    CovarianceMatrix covariance() const { Error m; fill( m ); return m; }
    /// return SMatrix
    Error error() const { Error m; fill( m ); return m; }
    /// fill SMatrix
    void fill( CovarianceMatrix & v ) const;

    /// Checks whether refitted tracks are stored.
    bool hasRefittedTracks() const {return !refittedTracks_.empty();}
    
    /// Returns the original track which corresponds to a particular refitted Track
    /// Throws an exception if now refitted tracks are stored ot the track is not found in the list
    TrackBaseRef originalTrack(const Track & refTrack) const;

    /// Returns the refitted track which corresponds to a particular original Track
    /// Throws an exception if now refitted tracks are stored ot the track is not found in the list
    Track refittedTrack(const TrackBaseRef & track) const;

    /// Returns the refitted track which corresponds to a particular original Track
    /// Throws an exception if now refitted tracks are stored ot the track is not found in the list
    Track refittedTrack(const TrackRef & track) const;

    /// Returns the container of refitted tracks
    const std::vector<Track> & refittedTracks() const { return refittedTracks_;}

  private:
    class TrackEqual {
      public:
	TrackEqual( const Track & t) : track_( t ) { }
	bool operator()( const Track & t ) const { return t.pt()==track_.pt();}
      private:
	const Track & track_;
    };
    /// chi-sqared
    Double32_t chi2_;
    /// number of degrees of freedom
    Double32_t ndof_;
    /// position
    Point position_;
    /// covariance matrix (3x3) as vector
    Double32_t covariance_[ size ];
    /// reference to tracks
    std::vector<TrackBaseRef > tracks_;
    /// The vector of refitted tracks
    std::vector<Track> refittedTracks_;
    std::vector<float> weights_;
    /// tells wether the vertex is really valid.
    bool validity_;


    /// position index
    index idx( index i, index j ) const {
      int a = ( i <= j ? i : j ), b = ( i <= j ? j : i );
      return b * ( b + 1 ) / 2 + a;
    }
  };
  
}

#endif
