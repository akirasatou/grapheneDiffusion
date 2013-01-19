#pragma once
#include "nonlinear_implicit_system.h"
#include "mesh_base.h"
#include "sparse_matrix.h"
#include "DiffusionSolver1DDescriptor.h"
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

  ResidualAndJacobianDiffusion(const DiffusionSolver1DDescriptor &difDsc, PoissonSolver2D &poisson);
  void residual(const NumericVector<Number> &X,
                NumericVector<Number> &R,
		NonlinearImplicitSystem &sys);
  void jacobian(const NumericVector<Number> &X,
		SparseMatrix<Number> &J, 
		NonlinearImplicitSystem &sys);

  virtual ~ResidualAndJacobianDiffusion();


private:

  const DiffusionSolver1DDescriptor &_difDsc;
  PoissonSolver2D &_poisson;

};
