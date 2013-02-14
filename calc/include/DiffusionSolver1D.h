#pragma once
#include "DiffusionSolver1DDescriptor.h"
#include "RealSpaceGridHandler.h"
#include "DiffusionABCalculator.h"
#include "PoissonDiffusionMediator.h"
#include "ResidualAndJacobianDiffusion.h"
#include <Transistor2D/PoissonSolver2D.h>
#include <libmesh.h>
#include <mesh.h>
#include <equation_systems.h>
//#include <mesh_refinement.h>
#include <string>
#include <vector>


/*
 * Class solving the Diffusion equation using libMesh.
 */

class DiffusionSolver1D
{

public:

  DiffusionSolver1D(const DiffusionSolver1DDescriptor &difDsc,
		    const DiffusionABCalculator &ab,
		    PoissonDiffusionMediator &pdm,
		    double xl, double xr);
  ~DiffusionSolver1D();
  RealSpaceGridHandler getRealSGH() const;
  void solveStep(double t, double dt);


private:

  static const std::string _className;
  std::string _sysName;
  Mesh _mesh;
  EquationSystems _es;
  //MeshRefinement _mRef;
  const DiffusionSolver1DDescriptor &_difDsc;
  const DiffusionABCalculator &_ab;
  PoissonDiffusionMediator &_pdm;
  ResidualAndJacobianDiffusion *_rj;
  double _Xl, _Xr;
  int _nSteps;

  enum BoundaryID { BoundaryIDNone, BoundaryIDLeft, BoundaryIDRight };

  void _setBoundaryID();
  void _setX();

};
