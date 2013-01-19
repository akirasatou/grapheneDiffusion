#pragma once
#include "RealSpaceArrayDiffusion.h"


/*
 * Class for 2D charge density that is passed to PoissonSolver.
 * It Bundles concentrations for electrons, holes, and doping.
 */

class ChargeDensity2D: public RealSpaceInterpolator
{
public:

  ChargeDensity2D(RealSpaceArrayDiffusion &SigmaElectron, 
		  RealSpaceArrayDiffusion &SigmaHole,
		  RealSpaceArrayDiffusion &Sigma0);

  double interpolate(double x) const;
  void updateInterpolator();


private:

  RealSpaceArrayDiffusion &_SigmaElectron, &_SigmaHole, &_Sigma0;

};
