#ifndef MBA_LINK__
#define MBA_LINK__

#include <vector>
#include <mba/MBA.h>
#include <LSsystem.h>
#include <mba/UCButils.h>
#include <mba/PointAccessUtils.h>
#include <GoTools/geometry/SplineSurface.h>
#include <GoTools/geometry/SplineInterpolator.h>
#include <GoTools/geometry/GeometryTools.h>
#include "util_mba_link.hpp"

using std::vector;

class MBALink {

public:

  //! Default constructor
  MBALink () = default;

  /*! Constructor with dimension 1
   *
   */
  MBALink ( boost::shared_ptr< vector<double> > u,
            boost::shared_ptr< vector<double> > v,
            boost::shared_ptr< vector<double> > z,
            boost::shared_ptr<Go::SplineSurface> surface_low
          );

  /*! Constructor with dimension 3
   *
   */
  MBALink ( boost::shared_ptr< vector<double> > u,
            boost::shared_ptr< vector<double> > v,
            boost::shared_ptr< vector<double> > x,
            boost::shared_ptr< vector<double> > y,
            boost::shared_ptr< vector<double> > z,
            boost::shared_ptr<Go::SplineSurface> surface_low
          );

  //! Desctructor
  ~MBALink () {};

  //! Compute parameter and create UCB and Go surfaces
  void
  MBAalg ( int m0 = 1,
           int n0 = 1,
           int h = 7,
           int smoothing_steps = 0
         );

  //! Compute the parameters with the least square procedure and create UCB and Go surfaces
  void
  LSalg ( int m0 = 1,
          int n0 = 1,
          int h = 7,
          int smoothing_steps = 0,
          int n = 100,
          double smoothing_fact = 1.          
        );

  //! Return u points
  boost::shared_ptr< vector<double> >
  get_u () const {
    return u_;
  };

  //! Return v points
  boost::shared_ptr< vector<double> >
  get_v () const {
    return u_;
  };

  //! Return x points
  boost::shared_ptr< vector<double> >
  get_x () const {
    return x_;
  };

  //! Return y points
  boost::shared_ptr< vector<double> >
  get_y () const {
    return y_;
  };

  //! Return z points
  boost::shared_ptr< vector<double> >
  get_z () const {
    return z_;
  };

  //! Return the dimension
  int
  get_dim () const {
    return dim_;
  };

  //! Return the fittet surface
  boost::shared_ptr<Go::SplineSurface>
  get_surface () const {
    return go_surface_;
  };

//   /*! Convert from UCBspl::SplineSurface to Go::SplineSurface
//    *
//    */
//   void Convert ();

private:

  //! Compute the MBA approximation and return a shared pointer to a SplineSurface
  void
  UCB2Go ( vector<UCBspl::SplineSurface>& surf );

  //! Compute residual between Hi-Fi points and Lo-Fi prediction
  void
  ComputeResidual ();

  //! Vector of u parameters
  boost::shared_ptr< vector<double> > u_;
  //! Vector of v parameters
  boost::shared_ptr< vector<double> > v_;
  //! Vector of x points
  boost::shared_ptr< vector<double> > x_;
  //! Vector of y points
  boost::shared_ptr< vector<double> > y_;
  //! Vector of z points
  boost::shared_ptr< vector<double> > z_;
  //! Dimension of the points (1 or 3)
  int dim_;
  //! Vector of x residuals points
  boost::shared_ptr< vector<double> > x_res_;
  //! Vector of y residuals points
  boost::shared_ptr< vector<double> > y_res_;
  //! Vector of z residuals points
  boost::shared_ptr< vector<double> > z_res_;
  //! MBA object
  vector<MBA> mba_;
  //! LS object
  vector< boost::shared_ptr<LSsystem> > lssys_;
  //! UCBspl::SplineSurface
  vector<UCBspl::SplineSurface> ucb_surface_;
  //! Go::SplineSurface
  boost::shared_ptr<Go::SplineSurface> go_surface_;
  //! Low fidelity Go::SplineSurface
  boost::shared_ptr<Go::SplineSurface> surface_low_;

};
// } // End of namespace
#endif
