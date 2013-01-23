#pragma once

#include "RealSpaceGridHandler.h"
#include "RealSpaceArrayDiffusion.h"


class DiffusionSolutionHolder
{

public:

  DiffusionSolutionHolder(const RealSpaceGridHandler &realSGH);

  void setInitialSolution(const RealSpaceArrayDiffusion &mue0,
			  const RealSpaceArrayDiffusion &muh0,
			  const RealSpaceArrayDiffusion &Ex0);


private:
  
  RealSpaceArrayDiffusion _mue_n, _muh_n, _Ex_n;
  RealSpaceArrayDiffusion _mue_n1_l, _muh_n1_l, _Ex_n1_l;

};
