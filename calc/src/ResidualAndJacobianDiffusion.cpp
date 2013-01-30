#include "ResidualAndJacobianDiffusion.h"

using namespace std;


/*
 * Constructor.
 */

ResidualAndJacobianDiffusion::
ResidualAndJacobianDiffusion(const DiffusionABCalculator &ab,
			     PoissonDiffusionMediator &pdm,
			     const MeshBase &meshBase,
			     const DofMap &dofMap):
  _ab(ab), _pdm(pdm), _meshBaseRef(&meshBase),
  _dofMapRef(&dofMap), _doOnceFlag(false)
{
}


/*
 * Destructor.
 */

ResidualAndJacobianDiffusion::~ResidualAndJacobianDiffusion()
{
}


/*
 * Compute the residual 'R' and Jacobian 'J' for the current 
 * solution vector 'X'.
 */

void ResidualAndJacobianDiffusion::
residual(const NumericVector<Number> &U, NumericVector<Number> &R,
	 NonlinearImplicitSystem &sys)
{
  _setNextSolutionsInNonlinearIteration(U, sys);
}

void ResidualAndJacobianDiffusion::
jacobian(const NumericVector<Number> &U, SparseMatrix<Number> &J, 
	 NonlinearImplicitSystem &sys)
{
  _setNextSolutionsInNonlinearIteration(U, sys);
}


/*
 * Set current time and increment to the next time step.
 */

void ResidualAndJacobianDiffusion::setTime(double t, double dt)
{
  _t = t;
  _dt = dt;
}


/*
 * Called both in residual() and jacobian(), but the main portion
 * is run in only one of them.
 */

void ResidualAndJacobianDiffusion::
_setNextSolutionsInNonlinearIteration(const NumericVector<Number> &U,
				      NonlinearImplicitSystem &sys)
{
  
  // If called already, reset the flag and exit.

  if( _doOnceFlag ){
    _doOnceFlag = false;
    return;
  }


  // Main portion.

  _doOnceFlag = true;


  // Save current solutions in nonlinear iteration and 

  _pdm.shiftSolutionsInNonlinearIteration();


  // Calculate $\mu_{r,n+1}^{(l+1)}$ from $U$.

  // Calculate $E_{x,n+1}^{(l+1)}$ from $\mu_{r,n+1}^{(l+1)}$.

  _pdm.updateSolutionsInNonlinearIteration();

}
