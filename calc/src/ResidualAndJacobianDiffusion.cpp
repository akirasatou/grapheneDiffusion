#include "ResidualAndJacobianDiffusion.h"
#include "fe.h"
#include "elem.h"
#include "side.h"
#include "quadrature_gauss.h"
#include "boundary_info.h"
#include "dense_vector.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include <algorithm>
#include <numeric>

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
  _dofMapRef(&dofMap)
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
  R.zero();

  // Numeric ids.

  const unsigned int id_e = sys.variable_number("mu_e");
  const unsigned int id_h = sys.variable_number("mu_h");


  // The definition of matrices and vectors for $\mu_{e}$ and 
  // $\mu_{h}$:
  // U^{e} = [U_{ee}^{e}, U_{hh}^{e}]^{t}
  // F^{e} = [F_{ee}^{e}, F_{hh}^{e}]^{t}
  // K^{e} = |K_{ee}^{e} K_{eh}^{e}|
  //         |K_{he}^{e} K_{hh}^{e}| where K_{eh}^{e}=K_{he}^{e}=0.

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe, Ue, Re;

  DenseSubMatrix<Number> Kee(Ke), Khh(Ke); // Keh(Ke), Khe(Ke), 
  DenseSubVector<Number> Fee(Fe), Fhh(Fhh);


  // Loop for each element.

  MeshBase::const_element_iterator el = _meshBaseRef->active_elements_begin();
  const MeshBase::const_element_iterator end_el = _meshBaseRef->active_elements_end();

  for(int iElem=0; el!=end_el; ++el, ++iElem){
    std::vector<unsigned int> dofInd, dofInd_e, dofInd_h;

    _dofMapRef->dof_indices(*el, dofInd);
    _dofMapRef->dof_indices(*el, dofInd_e, id_e);
    _dofMapRef->dof_indices(*el, dofInd_h, id_h);

    const unsigned int nDofs = dofInd.size();
    const unsigned int nDofs_e = dofInd_e.size();
    const unsigned int nDofs_h = dofInd_h.size();
    

    // Calculate the matrix Ke.

    Ke.resize(nDofs, nDofs);
    Ke.zero();
    Kee.reposition(id_e*nDofs_e, id_e*nDofs_e, nDofs_e, nDofs_e);
    //Keh.reposition(id_e*nDofs_e, id_h*nDofs_h, nDofs_e, nDofs_h);
    //Khe.reposition(id_h*nDofs_h, id_e*nDofs_e, nDofs_h, nDofs_e);
    Khh.reposition(id_h*nDofs_h, id_h*nDofs_h, nDofs_h, nDofs_h);
    //Keh.zero();
    //Khe.zero();

    // assembleMatrixInElement(_t, Kee, *el, dofInd_e, -1);
    // assembleMatrixInElement(_t, Khh, *el, dofInd_h, +1);


    // Setup Ue.

    Ue.resize(nDofs);
    Ue.zero();
    for(unsigned int i=0; i<nDofs; i++) Ue(i) = U(dofInd[i]);


    // Re = Ke*Ue.

    Re.resize(nDofs);
    Re.zero();

    Ke.vector_mult(Re, Ue);


    // Calculate the force vector Fe.

    Fe.resize(nDofs);
    Fe.zero();
    Fee.reposition(id_e*nDofs_e, nDofs_e);
    Fhh.reposition(id_e*nDofs_h, nDofs_h);

    //_assmNSS.assembleVectorInElement(_t, Fee, *el, dofInd_e, -1);
    //_assmNSS.assembleVectorInElement(_t, Fhh, *el, dofInd_h, +1);

    for(int i=0; i<nDofs; i++) Re(i) -= Fe(i);


    // Re = Ke*Ue-Fe.

    _dofMapRef->constrain_element_vector(Re, dofInd);
    R.add_vector(Re, dofInd);
  }
}

void ResidualAndJacobianDiffusion::
jacobian(const NumericVector<Number> &U, SparseMatrix<Number> &J, 
	 NonlinearImplicitSystem &sys)
{

  // Update the current solution in nonlinear iteration using the
  // current solution vector $U$ obtained in FEM. 
  // Called in jacobian() but not in residual() because it seems that
  // jacobian() is called first in the nonlinear iteration in libMesh.
  // Besides, residual() is called twice in one iteration.

  _setNextSolutionsInNonlinearIteration(U, sys);


  // Numeric ids.

  const unsigned int id_e = sys.variable_number("mu_e");
  const unsigned int id_h = sys.variable_number("mu_h");


  // The definition of matrices and vectors for $\mu_{e}$ and 
  // $\mu_{h}$:
  // U^{e} = [U_{ee}^{e}, U_{hh}^{e}]^{t}
  // K^{e} = |K_{ee}^{e} K_{eh}^{e}|
  //         |K_{he}^{e} K_{hh}^{e}| where K_{eh}^{e}=K_{he}^{e}=0.

  DenseMatrix<Number> Je;
  DenseVector<Number> Ue;

  DenseSubMatrix<Number> Jee(Je), Jhh(Je); // Jeh(Je), Jhe(Je),
  DenseVector<Number> Uee(Ue), Uhh(Ue);


  // Loop for each element.
 
  MeshBase::const_element_iterator el = _meshBaseRef->active_elements_begin();
  const MeshBase::const_element_iterator end_el = _meshBaseRef->active_elements_end();

  for(int iElem=0; el!=end_el; ++el, ++iElem){
    std::vector<unsigned int> dofInd, dofInd_e, dofInd_h;

    _dofMapRef->dof_indices(*el, dofInd);
    _dofMapRef->dof_indices(*el, dofInd_e, id_e);
    _dofMapRef->dof_indices(*el, dofInd_h, id_h);

    const unsigned int nDofs = dofInd.size();
    const unsigned int nDofs_e = dofInd_e.size();
    const unsigned int nDofs_h = dofInd_h.size();


    // Je = Ke.

    Je.resize(nDofs, nDofs);
    Je.zero();
    Jee.reposition(id_e*nDofs_e, id_e*nDofs_e, nDofs_e, nDofs_e);
    //Jeh.reposition(id_e*nDofs_e, id_h*nDofs_h, nDofs_e, nDofs_h);
    //Jhe.reposition(id_h*nDofs_h, id_e*nDofs_e, nDofs_h, nDofs_e);
    Jhh.reposition(id_h*nDofs_h, id_h*nDofs_h, nDofs_h, nDofs_h);
    //Jeh.zero();
    //Jhe.zero();

    // assembleMatrixInElement(_t, Jee, *el, dofInd_e, -1);
    // assembleMatrixInElement(_t, Jhh, *el, dofInd_h, +1);


    // Setup Ue.

    Ue.resize(nDofs);
    Ue.zero();
    for(unsigned int i=0; i<nDofs; i++) Ue(i) = U(dofInd[i]);


    // Je = Ke-dKe/dUe*Ue.

    // add_dKe_dUe(_t, Jee, *el, *dofInd_e, -1);
    // add_dKe_dUe(_t, Jhh, *el, *dofInd_h, +1);
  }
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
 * Update the current solution in nonlinear iteration using the
 * current solution vector $U$ obtained in FEM:
 * $\mu_{r,n+1}^{(l+1)}$, $E_{x,n+1}^{(l+1)}$, and 
 * $dE_{x,n+1}^{(l+1)}$.
 */

void ResidualAndJacobianDiffusion::
_setNextSolutionsInNonlinearIteration(const NumericVector<Number> &U,
				      NonlinearImplicitSystem &sys)
{

  // Save current solutions in nonlinear iteration and 

  _pdm.shiftSolutionsInNonlinearIteration();


  // Calculate $\mu_{r,n+1}^{(l+1)}$ from $U$.

  // Calculate $E_{x,n+1}^{(l+1)}$ from $\mu_{r,n+1}^{(l+1)}$.

  _pdm.updateSolutionsInNonlinearIteration();

}
