#pragma once
#include "nonlinear_implicit_system.h"
#include "mesh_base.h"
#include "sparse_matrix.h"
#include "DiffusionABCalculator.h"
#include "PoissonDiffusionMediator.h"
#include <Transistor2D/PoissonSolver2D.h>


/*
 * Class that implements 'ComputerResidualandJacobian' class to
 * compute the residual 'R' and Jacobian 'J'.
 */

class ResidualAndJacobianDiffusion:
  public NonlinearImplicitSystem::ComputeResidual,
  public NonlinearImplicitSystem::ComputeJacobian
{
public:

  ResidualAndJacobianDiffusion(const DiffusionABCalculator &ab,
			       PoissonDiffusionMediator &pdm,
			       const MeshBase &meshBase,
			       const DofMap &dofMap);
  void residual(const NumericVector<Number> &U,
                NumericVector<Number> &R,
		NonlinearImplicitSystem &sys);
  void jacobian(const NumericVector<Number> &U,
		SparseMatrix<Number> &J, 
		NonlinearImplicitSystem &sys);
  void setTime(double t, double dt);

  virtual ~ResidualAndJacobianDiffusion();


private:

  const DiffusionABCalculator &_ab;
  PoissonDiffusionMediator &_pdm;
  const MeshBase *_meshBaseRef;
  const DofMap *_dofMapRef;
  double _t, _dt;

  void _setNextSolutionsInNonlinearIteration(const NumericVector<Number> &U,
					     NonlinearImplicitSystem &sys);

};
