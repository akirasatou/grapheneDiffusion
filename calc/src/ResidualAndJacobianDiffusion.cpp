#include "ResidualAndJacobianDiffusion.h"

using namespace std;


/*
 * Constructor.
 */

ResidualAndJacobianDiffusion::
ResidualAndJacobianDiffusion(const DiffusionABCalculator &ab,
			     DiffusionSolutionHolder &dsh,
			     const MeshBase &meshBase,
			     const DofMap &dofMap,
			     double dt):
  _ab(ab), _dsh(dsh), _meshBaseRef(&meshBase),
  _dofMapRef(&dofMap), _dt(dt)
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
}

void ResidualAndJacobianDiffusion::
jacobian(const NumericVector<Number> &U, SparseMatrix<Number> &J, 
	 NonlinearImplicitSystem &sys)
{
}
