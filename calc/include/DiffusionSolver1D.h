#pragma once
#include "DiffusionSolver1DDescriptor.h"
#include "RealSpaceGridHandler.h"
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
		    PoissonSolver2D &poisson, double xl, double xr);
  ~DiffusionSolver1D();
  RealSpaceGridHandler getRealSGH();
  void solveStep();


private:

  static const std::string _className;
  std::string _sysName;
  Mesh _mesh;
  EquationSystems _es;
  //MeshRefinement _mRef;
  const DiffusionSolver1DDescriptor &_difDsc;
  PoissonSolver2D &_poisson;
  double _Xl, _Xr;
  int _nrSteps;

  enum BoundaryID { BoundaryIDNone, BoundaryIDLeft, BoundaryIDRight };

  void _setBoundaryID();
  void _setX();

};
