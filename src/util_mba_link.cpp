#include "util_mba_link.hpp"

using std::vector;

//! Compute the MBA approximation and return a SplineSurface
Go::SplineSurface
UCB2Go ( UCBspl::SplineSurface surf_ucb )
{
  // Obtain the matrix of coefficients
  boost::shared_ptr< GenMatrixType > phi = surf_ucb.getCoefficients ();
  int noX = phi->noX (), noY = phi->noY ();

  // Create an empty spline interpolator
  Go::SplineInterpolator interpol;
  interpol.setFreeConditions ();

  // Create vector of parameters and predicted points on a regular grid noX x noY
  vector<double> u ( noX, 0. ), v ( noY, 0. ), point;

  double umin = surf_ucb.umin (), umax = surf_ucb.umax ();
  double dx = ( umax - umin ) / ( noX - 1 );

  double vmin = surf_ucb.vmin (), vmax = surf_ucb.vmax ();
  double dy = ( vmax - vmin ) / ( noY - 1 );

  for ( size_t i = 0; i < noX; ++i )
    u[i] = umin + dx * i;

  for ( size_t i = 0; i < noY; ++i )
    v[i] = vmin + dy * i;

  point.reserve ( noX * noY );

  for ( size_t i = 0; i < noY; ++i )
    for ( size_t j = 0; j < noX; ++j )
      point.push_back ( surf_ucb.f ( u[j], v[i] ) );

  // Approximate the MBA surface with a SplineSurface
  Go::SplineSurface surf_go;
  surf_go.interpolate ( interpol, interpol, u.size (), v.size (), 1, &u[0], &v[0], &point[0] );
  
  return surf_go;
};

//! Compute the MBA approximation and return a SplineSurface
Go::SplineSurface
UCB2Go ( vector<UCBspl::SplineSurface> surf )
{
  // Obtain the matrix of coefficients
  vector<boost::shared_ptr< GenMatrixType > > phi;
  vector<int> noX, noY;

  int dim = 3;
  
  phi.reserve ( dim );
  
  for ( size_t i = 0; i < surf.size (); ++i ) {
    phi.push_back ( surf[i].getCoefficients () );
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

  for ( size_t i = 0; i < u.size(); ++i )
    u[i] = umin + dx * i;

  for ( size_t i = 0; i < v.size(); ++i )
    v[i] = vmin + dy * i;

  point.reserve ( u.size () * v.size () * dim );

  for ( size_t i = 0; i < v.size(); ++i )
    for ( size_t j = 0; j < u.size(); ++j )
      for ( size_t k = 0; k < dim; ++k )
        point.push_back ( surf[k].f ( u[j], v[i] ) );

  // Approximate the MBA surface with a SplineSurface
  Go::SplineSurface surf_go;      
  surf_go.interpolate ( interpol, interpol, u.size (), v.size (), dim, &u[0], &v[0], &point[0] ); 
  
  return surf_go;
}


