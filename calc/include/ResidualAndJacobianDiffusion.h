#pragma once
#include "nonlinear_implicit_system.h"
#include "elem.h"
#include "mesh_base.h"
#include "sparse_matrix.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include "DiffusionABCalculator.h"
#include "PoissonDiffusionMediator.h"
#include "RealSpaceArrayDiffusion.h"
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
			       const DofMap &dofMap,
			       const RealSpaceGridHandler &realSGH);
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
  RealSpaceArrayDiffusion _mue, _dmue_dx, _d2mue_dx2;
  RealSpaceArrayDiffusion _muh, _dmuh_dx, _d2muh_dx2;
  RealSpaceArrayDiffusion _Ex, _dEx_dx;
  RealSpaceGridHandler _realSGH;
  double _t, _dt;

  void _setNextSolutionsInNonlinearIteration(const NumericVector<Number> &U,
					     NonlinearImplicitSystem &sys);

  void _addK(DenseSubMatrix<Number> &K, const Elem *elem, int iElem,
	     const std::vector<unsigned int> &dofInd, int s) const;
  void _add_dKdU_U(DenseSubMatrix<Number> &dKdU, const Elem *elem, 
		   int iElem, const NumericVector<Number> &U,
		   const std::vector<unsigned int> &dofInd,
		   int s) const;
  void _addF(DenseSubVector<Number> &F, const Elem *elem, int iElem,
	     const std::vector<unsigned int> &dofInd, int s) const;

};
