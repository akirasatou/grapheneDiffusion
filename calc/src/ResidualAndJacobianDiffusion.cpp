#include <PhysicalConstants.h>
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
			     const DofMap &dofMap,
			     const RealSpaceGridHandler &realSGH):
  _ab(ab), _pdm(pdm), _meshBaseRef(&meshBase), _dofMapRef(&dofMap),
  _mue(realSGH), _dmue_dx(realSGH), _d2mue_dx2(realSGH),
  _muh(realSGH), _dmuh_dx(realSGH), _d2muh_dx2(realSGH),
  _Ex(realSGH), _dEx_dx(realSGH), _realSGH(realSGH)
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

  DenseSubMatrix<Number> Kee(Ke), Khh(Ke);
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
    Khh.reposition(id_h*nDofs_h, id_h*nDofs_h, nDofs_h, nDofs_h);

    _addK(Kee, *el, iElem, -1);
    _addK(Khh, *el, iElem, +1);


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

    _addF(Fee, *el, iElem, -1);
    _addF(Fhh, *el, iElem, +1);

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

  _setNextSolutionsNI(U, sys);


  // Numeric ids.

  const unsigned int id_e = sys.variable_number("mu_e");
  const unsigned int id_h = sys.variable_number("mu_h");


  // The definition of matrices and vectors for $\mu_{e}$ and 
  // $\mu_{h}$:
  // U^{e} = [U_{ee}^{e}, U_{hh}^{e}]^{t}
  // K^{e} = |K_{ee}^{e} K_{eh}^{e}|
  //         |K_{he}^{e} K_{hh}^{e}| where K_{eh}^{e}=K_{he}^{e}=0.

  DenseMatrix<Number> Je;
  DenseSubMatrix<Number> Jee(Je), Jhh(Je); // Jeh(Je), Jhe(Je),


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
    Jhh.reposition(id_h*nDofs_h, id_h*nDofs_h, nDofs_h, nDofs_h);

    _addK(Jee, *el, iElem, -1);
    _addK(Jhh, *el, iElem, +1);


    // Je = Ke+dKe/dUe*Ue.

    _add_dKdU_U(Jee, *el, iElem, U, dofInd_e, -1);
    _add_dKdU_U(Jhh, *el, iElem, U, dofInd_h, +1);
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
_setNextSolutionsNI(const NumericVector<Number> &U,
		    NonlinearImplicitSystem &sys)
{

  // Calculate $\mu_{r,n+1}^{(l+1)}$ from $U$.

  const unsigned int dim = _meshBaseRef->mesh_dimension();
  FEType feType = _dofMapRef->variable_type(0);
    
  AutoPtr<FEBase> fe(FEBase::build(dim, feType));
  QGauss qrule(dim, FIFTH);
  fe->attach_quadrature_rule(&qrule);
  const vector<vector<Real> > &phi = fe->get_phi();
  const vector<vector<Real> > &dphidx = fe->get_dphidx();
  const vector<vector<Real> > &d2phidx2 = fe->get_d2phidx2();

  const unsigned int id_e = sys.variable_number("mu_e");
  const unsigned int id_h = sys.variable_number("mu_h");

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

    vector<double> tmp = _realSGH.getPointsInElem(iElem);
    vector<Point> qps(tmp.size()), v;

    for(int i=0; i<qps.size(); i++) qps[i] = Point(tmp[i]);

    FE<1, LAGRANGE>::inverse_map(*el, qps, v);
    fe->reinit(*el, &v);

    for(unsigned int j=0; j<qps.size(); j++){
      int i = _realSGH.getPointInvID(iElem, j);

      _mue.setAt(i, 0.0);
      _dmue_dx.setAt(i, 0.0);
      _d2mue_dx2.setAt(i, 0.0);

      for(unsigned int k=0; k<dofInd_e.size(); k++){
	_mue.addAt(i, phi[k][j]*U(dofInd_e[k]));
	_dmue_dx.addAt(i, dphidx[k][j]*U(dofInd_e[k]));
	_d2mue_dx2.addAt(i, d2phidx2[k][j]*U(dofInd_e[k]));
      }

      _muh.setAt(i, 0.0);
      _dmuh_dx.setAt(i, 0.0);
      _d2muh_dx2.setAt(i, 0.0);

      for(unsigned int k=0; k<dofInd_h.size(); k++){
	_muh.addAt(i, phi[k][j]*U(dofInd_h[k]));
	_dmuh_dx.addAt(i, dphidx[k][j]*U(dofInd_h[k]));
	_d2muh_dx2.addAt(i, d2phidx2[k][j]*U(dofInd_h[k]));
      }
    }

  }


  // Calculate $E_{x,n+1}^{(l+1)}$ from $\mu_{r,n+1}^{(l+1)}$.

  _pdm.updateSolutionsNI(_mue, _dmue_dx, _d2mue_dx2,
			 _muh, _dmuh_dx, _d2muh_dx2);

}


/*
 * Add the stiffness matrix $K$ for the element.
 */

void ResidualAndJacobianDiffusion::
_addK(DenseSubMatrix<Number> &K, const Elem *elem, int iElem,
      int s) const
{
  const unsigned int dim = _meshBaseRef->mesh_dimension();
  FEType feType = _dofMapRef->variable_type(0);
    
  AutoPtr<FEBase> fe(FEBase::build(dim, feType));
  QGauss qrule(dim, FIFTH);
  fe->attach_quadrature_rule(&qrule);
    
  AutoPtr<FEBase> feFace(FEBase::build(dim, feType));
  QGauss qface(dim-1, FIFTH);
  feFace->attach_quadrature_rule(&qface);

  fe->reinit(elem);

  const vector<Real> &JxW = fe->get_JxW();
  const vector<Point> &r = fe->get_xyz();
  const vector<vector<Real> > &phi = fe->get_phi();
  const vector<vector<Real> > &dphidx = fe->get_dphidx();
  const vector<vector<Real> > &d2phidx2 = fe->get_d2phidx2();

  for(unsigned int qp=0; qp<qrule.n_points(); qp++){
    int ir = _realSGH.getPointInvID(iElem, qp);
    double C = 0.0;
    double mu_n1_l1, dmu_dx_n1_l1;
    double Ex_n1_l1 = _pdm.get_Ex_n1_l1(ir);

    if( s == -1 ){
      mu_n1_l1 = _pdm.get_mue_n1_l1(ir);
      dmu_dx_n1_l1 = _pdm.get_dmue_dx_n1_l1(ir);
    }
    else {
      mu_n1_l1 = _pdm.get_muh_n1_l1(ir);
      dmu_dx_n1_l1 = _pdm.get_dmuh_dx_n1_l1(ir);
    }

    double A = _ab.calcA(mu_n1_l1);
    double B = _ab.calcB(mu_n1_l1);

    for(unsigned int i=0; i<phi.size(); i++){
      for(unsigned int j=0; j<phi.size(); j++){
	double tmp = 0.0;

	tmp += phi[j][qp];
	tmp += 0.5*s*_dt*A*(e*Ex_n1_l1-s*dmu_dx_n1_l1)*dphidx[j][qp];
	tmp += -0.5*_dt*B*d2phidx2[j][qp];

	K(i, j) += JxW[qp]*phi[i][qp]*tmp;
      }
    }
  }

}


/*
 * Add the $\sum_{k}\frac{\p K_{ik}^{e}}{\p U_{j}^{e}}U_{k}^{e}$
 * for the element.
 */

void ResidualAndJacobianDiffusion::
_add_dKdU_U(DenseSubMatrix<Number> &K, const Elem *elem, 
	    int iElem, const NumericVector<Number> &U,
	    const std::vector<unsigned int> &dofInd, int s) const
{
  const unsigned int dim = _meshBaseRef->mesh_dimension();
  FEType feType = _dofMapRef->variable_type(0);
    
  AutoPtr<FEBase> fe(FEBase::build(dim, feType));
  QGauss qrule(dim, FIFTH);
  fe->attach_quadrature_rule(&qrule);
    
  AutoPtr<FEBase> feFace(FEBase::build(dim, feType));
  QGauss qface(dim-1, FIFTH);
  feFace->attach_quadrature_rule(&qface);

  fe->reinit(elem);

  const vector<Real> &JxW = fe->get_JxW();
  const vector<Point> &r = fe->get_xyz();
  const vector<vector<Real> > &phi = fe->get_phi();
  const vector<vector<Real> > &dphidx = fe->get_dphidx();
  const vector<vector<Real> > &d2phidx2 = fe->get_d2phidx2();

  for(unsigned int qp=0; qp<qrule.n_points(); qp++){
    int ir = _realSGH.getPointInvID(iElem, qp);
    double mu_n1_l1, dmu_dx_n1_l1;
    double Ex_n1_l1 = _pdm.get_Ex_n1_l1(ir);

    if( s == -1 ){
      mu_n1_l1 = _pdm.get_mue_n1_l1(ir);
      dmu_dx_n1_l1 = _pdm.get_dmue_dx_n1_l1(ir);
    }
    else {
      mu_n1_l1 = _pdm.get_muh_n1_l1(ir);
      dmu_dx_n1_l1 = _pdm.get_dmuh_dx_n1_l1(ir);
    }

    double A = _ab.calcA(mu_n1_l1);
    double dA_dmu = _ab.calc_dA_dmu(mu_n1_l1);
    double dB_dmu = _ab.calc_dB_dmu(mu_n1_l1);

    for(unsigned int i=0; i<phi.size(); i++){
      for(unsigned int j=0; j<phi.size(); j++){
	double a = 0.5*s*_dt*phi[j][qp]*dA_dmu*(e*Ex_n1_l1-s*dmu_dx_n1_l1);

	for(unsigned int k=0; k<phi.size(); k++){
	  double tmp = 0.0;
	  double Uk = U(dofInd[k]); // ?

	  tmp += a*dphidx[k][qp];
	  tmp -= 0.5*_dt*A*dphidx[j][qp]*dphidx[k][qp];
	  tmp -= 0.5*_dt*dB_dmu*phi[j][qp]*d2phidx2[k][qp];

	  K(i, j) += JxW[qp]*phi[i][qp]*tmp*Uk;
	}
      }
    }
  }


}


/*
 * Add the load vector $F$ for the element.
 */

void ResidualAndJacobianDiffusion::
_addF(DenseSubVector<Number> &F, const Elem *elem, int iElem,
      int s) const
{
  const unsigned int dim = _meshBaseRef->mesh_dimension();
  FEType feType = _dofMapRef->variable_type(0);
    
  AutoPtr<FEBase> fe(FEBase::build(dim, feType));
  QGauss qrule(dim, FIFTH);
  fe->attach_quadrature_rule(&qrule);
    
  AutoPtr<FEBase> feFace(FEBase::build(dim, feType));
  QGauss qface(dim-1, FIFTH);
  feFace->attach_quadrature_rule(&qface);

  fe->reinit(elem);

  const vector<Real> &JxW = fe->get_JxW();
  const vector<Point> &r = fe->get_xyz();
  const vector<vector<Real> > &phi = fe->get_phi();

  for(unsigned int qp=0; qp<qrule.n_points(); qp++){
    int ir = _realSGH.getPointInvID(iElem, qp);
    double C = 0.0;
    double mu_n, dmu_dx_n, d2mu_dx2_n, mu_n1_l1;
    double Ex_n = _pdm.get_Ex_n(ir), dEx_dx_n = _pdm.get_dEx_dx_n(ir);
    double dEx_dx_n1_l1 = _pdm.get_dEx_dx_n1_l1(ir);

    if( s == -1 ){
      mu_n = _pdm.get_mue_n(ir);
      dmu_dx_n = _pdm.get_dmue_dx_n(ir);
      d2mu_dx2_n = _pdm.get_d2mue_dx2_n(ir);
      mu_n1_l1 = _pdm.get_mue_n1_l1(ir);
    }
    else {
      mu_n = _pdm.get_muh_n(ir);
      dmu_dx_n = _pdm.get_dmuh_dx_n(ir);
      d2mu_dx2_n = _pdm.get_d2muh_dx2_n(ir);
      mu_n1_l1 = _pdm.get_muh_n1_l1(ir);
    }

    C += mu_n-0.5*s*_dt*_ab.calcA(mu_n)*(e*Ex_n-s*dmu_dx_n)*dmu_dx_n;
    C += 0.5*_dt*_ab.calcB(mu_n)*(d2mu_dx2_n-s*e*Ex_n);
    C += -0.5*s*_dt*_ab.calcB(mu_n1_l1)*e*dEx_dx_n1_l1;

    for(unsigned int i=0; i<phi.size(); i++){
      F(i) += JxW[qp]*phi[i][qp]*C;
    }
  }
}
