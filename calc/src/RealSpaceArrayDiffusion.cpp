#include "RealSpaceArrayDiffusion.h"
#include <iostream>
#include <stdlib.h>

using namespace std;


/*
 * Constructor.
 */

RealSpaceArrayDiffusion::
RealSpaceArrayDiffusion(const RealSpaceGridHandler &realSGH):
  RealSpaceArrayWithVector<double>(realSGH.getSize()),
  _x(realSGH.getSize()), _realSGH(realSGH)
{
  for(int i=0; i<_x.size(); i++){
    _x[i] = _realSGH.getAt(i);
    _a[i] = 0.0;
  }

  _csi = WeightedCubicSplineInterpolator(&_a[0], &_x[0],
					 _realSGH.getSize());
}


/*
 * Destructor.
 */

RealSpaceArrayDiffusion::~RealSpaceArrayDiffusion()
{
}


/*
 * Interpolation for PoissonSolver (implementation for 
 * RealSpaceInterpolator).
 */

double RealSpaceArrayDiffusion::interpolate(double x) const
{
  return _csi.interpolate(x);
}


/*
 * Update the interpolator.
 */

void RealSpaceArrayDiffusion::updateInterpolator()
{
  _csi.update();
}
