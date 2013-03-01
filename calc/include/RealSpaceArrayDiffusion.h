#pragma once
#include <Interpolator/CubicSpline.h>
#include <RealSpaceArray/RealSpaceArrayWithVector.h>
#include <Transistor2D/RealSpaceInterpolator.h>
#include "RealSpaceGridHandler.h"

using namespace std;


/*
 * Class for real-space arrays that contain quantities
 * like e-field or concentration. 
 */

class RealSpaceArrayDiffusion:
  public RealSpaceArrayWithVector<double>,
  public RealSpaceInterpolator
{
public:

  RealSpaceArrayDiffusion(const RealSpaceGridHandler &realSGH);
  ~RealSpaceArrayDiffusion();
  double interpolate(double x) const;
  void updateInterpolator();
  double getX(int i) const;


protected:

  const RealSpaceGridHandler &_realSGH;
  WeightedCubicSplineInterpolator _csi;
  vector<double> _x;

  RealSpaceArrayDiffusion(const RealSpaceArrayDiffusion &a);
  RealSpaceArrayDiffusion &operator=(const RealSpaceArrayDiffusion &a);

};
