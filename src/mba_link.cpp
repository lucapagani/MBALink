#include "mba_link.hpp"
#include <boost/make_shared.hpp>

// Constructor
MBALink::MBALink ( boost::shared_ptr< vector<double> > u,
                   boost::shared_ptr< vector<double> > v,
                   boost::shared_ptr< vector<double> > z,
                   boost::shared_ptr<Go::SplineSurface> surface_low
                 ) : u_ ( u ), v_ ( v ), z_ ( z ), dim_ ( 1 ), surface_low_ ( surface_low )
{
  ucb_surface_.reserve ( dim_ );
  mba_.resize ( dim_ );
  lssys_.resize ( dim_ );

  z_res_.reset ( new vector<double> );
  
  MBALink::ComputeResidual ();

  mba_[0].init ( u_, v_, z_res_ );
}

// Constructor
MBALink::MBALink ( boost::shared_ptr< vector<double> > u,
                   boost::shared_ptr< vector<double> > v,
                   boost::shared_ptr< vector<double> > x,
                   boost::shared_ptr< vector<double> > y,
                   boost::shared_ptr< vector<double> > z,
                   boost::shared_ptr<Go::SplineSurface> surface_low
                 ) : u_ ( u ), v_ ( v ), x_ ( x ), y_ ( y ), z_ ( z ), dim_ ( 3 ), surface_low_ ( surface_low )
{
  ucb_surface_.reserve ( dim_ );
  mba_.resize ( dim_ );
  lssys_.resize ( dim_ );

  x_res_.reset ( new vector<double> );
  y_res_.reset ( new vector<double> );
  z_res_.reset ( new vector<double> );

  MBALink::ComputeResidual ();

  mba_[0].init ( u_, v_, x_res_ );
  mba_[1].init ( u_, v_, y_res_ );
  mba_[2].init ( u_, v_, z_res_ );
}

//! Compute parameter and create UCB and Go surfaces
void
MBALink::MBAalg ( int m_0,
                  int n_0,
                  int h,
                  int smoothing_steps
                )
{
  // Estimate the parameters
//   std::for_each ( mba_.begin (), mba_.end (), [&] ( MBA mba ) {
//     mba.MBAalg ( m_0, n_0, h, smoothing_steps );
//   } );

  double umin = std::min ( *std::min_element ( u_->begin (), u_->end() ), surface_low_->startparam_u() );
  double vmin = std::min ( *std::min_element ( v_->begin (), v_->end() ), surface_low_->startparam_v() );
  double umax = std::max ( *std::max_element ( u_->begin (), u_->end() ), surface_low_->endparam_u() );
  double vmax = std::max ( *std::max_element ( v_->begin (), v_->end() ), surface_low_->endparam_v() );

  surface_low_->setParameterDomain ( umin, umax, vmin, vmax );    
  
  for ( auto it = mba_.begin(); it != mba_.end(); ++it ) {
    it->setDomain ( umin, vmin, umax, vmax );
    it->MBAalg ( m_0, n_0, h, smoothing_steps );
  }

  // Fill the ucb_surface_ vector
  for ( auto it = mba_.begin (); it != mba_.end (); ++it )
    ucb_surface_.push_back ( it->getSplineSurface () );

  MBALink::UCB2Go ( ucb_surface_ );
}

//! Compute the parameters with the least square procedure and create UCB and Go surfaces
void
MBALink::LSalg ( int m_0,
                 int n_0,
                 int h,
                 int smoothing_steps,
                 int n,
                 double smoothing_fact
               )
{
  MBAalg ( m_0, n_0, h, smoothing_steps );

  vector< boost::shared_ptr<GenMatrixType> > phi;
  phi.reserve ( dim_ );

  for ( auto it = mba_.begin(); it != mba_.end(); ++it )
    phi.push_back ( it->PHI () );

  int no_x = phi[0]->noX ();
  int no_y = phi[0]->noY ();

  if ( dim_ == 1 ) {
    lssys_[0].reset ( new LSsystem ( u_, v_, z_res_, no_x, no_y, phi[0], smoothing_fact ) );
    lssys_[0]->buildEqSystem();
    lssys_[0]->relaxCG ( n );
  } else {
    lssys_[0].reset ( new LSsystem ( u_, v_, x_res_, no_x, no_y, phi[0], smoothing_fact ) );
    lssys_[1].reset ( new LSsystem ( u_, v_, y_res_, no_x, no_y, phi[1], smoothing_fact ) );
    lssys_[2].reset ( new LSsystem ( u_, v_, z_res_, no_x, no_y, phi[2], smoothing_fact ) );

    lssys_[0]->buildEqSystem();
    lssys_[1]->buildEqSystem();
    lssys_[2]->buildEqSystem();

    lssys_[0]->relaxCG ( n );
    lssys_[1]->relaxCG ( n );
    lssys_[2]->relaxCG ( n );
  }

  MBALink::UCB2Go ( ucb_surface_ );

}

//! Compute the MBA approximation and return a shared pointer to a SplineSurface
void
MBALink::UCB2Go ( vector<UCBspl::SplineSurface>& surf )
{
  // Obtain the matrix of coefficients
  vector<boost::shared_ptr< GenMatrixType > > phi;
  vector<int> noX, noY;

  phi.reserve ( dim_ );

  boost::shared_ptr< GenMatrixType > temp ( new GenMatrixType );

  for ( size_t i = 0; i < surf.size (); ++i ) {
//     phi[i] = surf[i].getCoefficients ();
//     phi[i].reset ( new GenMatrixType );
    temp->init ( *surf[i].getCoefficients () );
    phi.push_back ( temp );
    noX.push_back ( phi[i]->noX () );
    noY.push_back ( phi[i]->noY () );
  }

  // Create an empty spline interpolator
  Go::SplineInterpolator interpol;
  interpol.setFreeConditions ();

  // Create vector of parameters and predicted points on a regular grid noX x noY
  vector<double> u ( *std::max_element ( noX.begin (), noX.end () ), 0. ), v ( *std::max_element ( noY.begin (), noY.end () ), 0. ), point;

  double umin = surf[0].umin (), umax = surf[0].umax ();
  double dx = ( umax - umin ) / ( u.size () - 1 );

  double vmin = surf[0].vmin (), vmax = surf[0].vmax ();
  double dy = ( vmax - vmin ) / ( v.size () - 1 );

//   double umin = surface_low_->startparam_u(), umax = surface_low_->endparam_u();
//   double dx = ( umax - umin ) / ( u.size () - 1 );
//
//   double vmin = surface_low_->startparam_v(), vmax = surface_low_->endparam_v();
//   double dy = ( vmax - vmin ) / ( v.size () - 1 );

  /*  double umin = std::max ( surface_low_->startparam_u(), surf[0].umin () ), umax = std::min ( surface_low_->endparam_u(), surf[0].umax () );
    double dx = ( umax - umin ) / ( u.size () - 1 );

    double vmin = std::max ( surface_low_->startparam_v(), surf[0].vmin () ), vmax = std::min ( surface_low_->endparam_v(), surf[0].vmax () );
    double dy = ( vmax - vmin ) / ( v.size () - 1 );  */

  for ( size_t i = 0; i < u.size(); ++i )
    u[i] = umin + dx * i;

  for ( size_t i = 0; i < v.size(); ++i )
    v[i] = vmin + dy * i;

  point.reserve ( u.size () * v.size () * dim_ );

  for ( size_t i = 0; i < v.size(); ++i )
    for ( size_t j = 0; j < u.size(); ++j )
      for ( size_t k = 0; k < dim_; ++k )
        point.push_back ( surf[k].f ( u[j], v[i] ) );

  // Approximate the MBA surface with a SplineSurface
  Go::SplineSurface go_surface;
  go_surface.interpolate ( interpol, interpol, u.size (), v.size (), dim_, &u[0], &v[0], &point[0] );

//   vector<double> num_tol ( 4, 0. );
//   num_tol[0] = fabs ( go_surface.startparam_u() - surface_low_->startparam_u() );
//   num_tol[1] = fabs ( go_surface.endparam_u() - surface_low_->endparam_u() );
//   num_tol[2] = fabs ( go_surface.startparam_v() - surface_low_->startparam_v() );
//   num_tol[3] = fabs ( go_surface.endparam_v() - surface_low_->endparam_v() );

//   shared_ptr<Go::SplineSurface> go_surface_ptr = Go::GeometryTools::surfaceSum ( go_surface, 1, *surface_low_, 1, *std::max_element ( num_tol.begin(), num_tol.end () ) );
  shared_ptr<Go::SplineSurface> go_surface_ptr = Go::GeometryTools::surfaceSum ( go_surface, 1, *surface_low_, 1 );
  go_surface_ = make_shared_ptr ( go_surface_ptr );
}

//! Compute residual between Hi-Fi points and Lo-Fi prediction
void
MBALink::ComputeResidual ()
{

  Go::Point point;
  z_res_->resize ( u_->size () );

  if ( dim_ == 1 ) {

    for ( size_t i = 0; i < u_->size (); ++i ) {
      surface_low_->point ( point, ( *u_ ) [i], ( *v_ ) [i] );
      ( *z_res_ ) [i] = ( *z_ ) [i] - point[0];
    }

  } else {

    x_res_->resize ( u_->size () );
    y_res_->resize ( u_->size () );
    for ( size_t i = 0; i < u_->size (); ++i ) {
      surface_low_->point ( point, ( *u_ ) [i], ( *v_ ) [i] );
      ( *x_res_ ) [i] = ( *x_ ) [i] - point[0];
      ( *y_res_ ) [i] = ( *y_ ) [i] - point[1];
      ( *z_res_ ) [i] = ( *z_ ) [i] - point[2];
    }

  }
}

