#ifndef UTIL_MBA_LINK__
#define UTIL_MBA_LINK__

#include <vector>
#include <mba/MBA.h>
#include <mba/UCButils.h>
#include <mba/PointAccessUtils.h>
#include <GoTools/geometry/SplineSurface.h>
#include <GoTools/geometry/SplineInterpolator.h>

using std::vector;

//! Compute the MBA approximation and return a SplineSurface
Go::SplineSurface
UCB2Go ( UCBspl::SplineSurface surf_ucb );

//! Compute the MBA approximation and return a SplineSurface
Go::SplineSurface
UCB2Go ( vector<UCBspl::SplineSurface> surf_ucb );

//! Convert std::shared_ptr to boost::shared_ptr
template<typename T>
boost::shared_ptr<T> make_shared_ptr ( std::shared_ptr<T>& ptr )
{
  return boost::shared_ptr<T> ( ptr.get(), [ptr] ( T* ) mutable {ptr.reset();} );
}

#endif
