#pragma once

#include "RealSpaceGridHandler.h"
#include "RealSpaceArrayDiffusion.h"
#include "FermiDistrGraphene.h"
#include <Transistor2D/PoissonSolver2D.h>

class PoissonDiffusionMediator
{

public:

  PoissonDiffusionMediator(const RealSpaceGridHandler &realSGH,
			   PoissonSolver2D &poisson,
			   const RealSpaceArrayDiffusion &SigmaDope,
			   const FermiDistrGraphene &fermiDistr);
  void setInitialSolutions(const RealSpaceArrayDiffusion &mue0,
			   const RealSpaceArrayDiffusion &muh0,
			   const RealSpaceArrayDiffusion &Ex0,
			   const RealSpaceArrayDiffusion &dEx_dx0);
  void getCurrentSolutions(RealSpaceArrayDiffusion &mue,
			   RealSpaceArrayDiffusion &muh,
			   RealSpaceArrayDiffusion &Ex) const;
  void updateSolutions();

  void setTime(double t, double dt);
  void setInitialSolutionsInNonlinearIteration();
  void updateSolutionsInNonlinearIteration(const RealSpaceArrayDiffusion &mue, const RealSpaceArrayDiffusion &muh);


private:
    
  RealSpaceArrayDiffusion _mue_n, _mue_n1_l, _mue_n1_l1;
  RealSpaceArrayDiffusion _muh_n, _muh_n1_l, _muh_n1_l1;
  RealSpaceArrayDiffusion _Ex_n, _Ex_n1_l, _Ex_n1_l1;
  RealSpaceArrayDiffusion _dEx_dx_n, _dEx_dx_n1_l, _dEx_dx_n1_l1;
  RealSpaceArrayDiffusion _SigmaDope;
  PoissonSolver2D &_poisson;
  double _t, _dt;
  const RealSpaceGridHandler &_realSGH;

  const FermiDistrGraphene _fermiDistr;

};
