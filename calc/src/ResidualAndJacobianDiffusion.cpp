#include "ResidualAndJacobianDiffusion.h"

using namespace std;


/*
 * Constructor.
 */

ResidualAndJacobianDiffusion::
ResidualAndJacobianDiffusion(const DiffusionSolver1DDescriptor &difDsc, PoissonSolver2D &poisson):
  _difDsc(difDsc), _poisson(poisson)
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
residual(const NumericVector<Number> &X, NumericVector<Number> &R,
	 NonlinearImplicitSystem &sys)
{
}

void ResidualAndJacobianDiffusion::
jacobian(const NumericVector<Number> &X, SparseMatrix<Number> &J, 
	 NonlinearImplicitSystem &sys)
{
}
