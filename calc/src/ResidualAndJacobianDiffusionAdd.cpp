#include <PhysicalConstants.h>
#include <PhysicalUnits.h>
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
#include <assert.h>

using namespace std;


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
  fe->reinit(elem);

  const vector<Real> &JxW = fe->get_JxW();
  const vector<vector<Real> > &phi = fe->get_phi();
  const vector<vector<Real> > &dphidx = fe->get_dphidx();
  const vector<vector<Real> > &d2phidx2 = fe->get_d2phidx2();

  for(unsigned int qp=0; qp<qrule.n_points(); qp++){
    int ir = _realSGH.getPointInvID(iElem, qp);
    double mu_n1_l1, dmu_dx_n1_l1;
    double eEx_n1_l1 = 0.0;//e*_pdm.get_Ex_n1_l1(ir)/_eExNorm;

    if( s == -1 ){
      mu_n1_l1 = _pdm.get_mue_n1_l1(ir)/_muNorm;
      dmu_dx_n1_l1 = _pdm.get_dmue_dx_n1_l1(ir)/_dmudxNorm;
    }
    else {
      mu_n1_l1 = _pdm.get_muh_n1_l1(ir)/_muNorm;
      dmu_dx_n1_l1 = _pdm.get_dmuh_dx_n1_l1(ir)/_dmudxNorm;
    }

    double A = _ab.calcA(_muNorm*mu_n1_l1)/_ANorm;
    double B = _ab.calcB(_muNorm*mu_n1_l1)/_BNorm;

    for(unsigned int i=0; i<phi.size(); i++){
      for(unsigned int j=0; j<phi.size(); j++){
	double tmp = 0.0;

	tmp += phi[j][qp];
	tmp += 0.5*s*(_dt/_tNorm)*A*(eEx_n1_l1-s*dmu_dx_n1_l1)*dphidx[j][qp];
	tmp += -0.5*(_dt/_tNorm)*B*d2phidx2[j][qp];

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
  fe->reinit(elem);

  const vector<Real> &JxW = fe->get_JxW();
  const vector<vector<Real> > &phi = fe->get_phi();
  const vector<vector<Real> > &dphidx = fe->get_dphidx();
  const vector<vector<Real> > &d2phidx2 = fe->get_d2phidx2();

  for(unsigned int qp=0; qp<qrule.n_points(); qp++){
    int ir = _realSGH.getPointInvID(iElem, qp);
    double mu_n1_l1, dmu_dx_n1_l1;
    double eEx_n1_l1 = 0.0;//e*_pdm.get_Ex_n1_l1(ir)/(_eExNorm/_muNorm);

    if( s == -1 ){
      mu_n1_l1 = _pdm.get_mue_n1_l1(ir)/_muNorm;
      dmu_dx_n1_l1 = _pdm.get_dmue_dx_n1_l1(ir)/_dmudxNorm;
    }
    else {
      mu_n1_l1 = _pdm.get_muh_n1_l1(ir)/_muNorm;
      dmu_dx_n1_l1 = _pdm.get_dmuh_dx_n1_l1(ir)/_dmudxNorm;
    }

    double A = _ab.calcA(_muNorm*mu_n1_l1)/_ANorm;
    double dA_dmu = _ab.calc_dA_dmu(_muNorm*mu_n1_l1)/(_ANorm/_muNorm);
    double dB_dmu = _ab.calc_dB_dmu(_muNorm*mu_n1_l1)/(_BNorm/_muNorm);

    for(unsigned int i=0; i<phi.size(); i++){
      for(unsigned int j=0; j<phi.size(); j++){
	double a = 0.5*s*(_dt/_tNorm)*phi[j][qp]*dA_dmu*(eEx_n1_l1-s*dmu_dx_n1_l1);

	for(unsigned int k=0; k<phi.size(); k++){
	  double tmp = 0.0;
	  double Uk = U(dofInd[k]); // ?

	  tmp += a*dphidx[k][qp];
	  tmp -= 0.5*(_dt/_tNorm)*A*dphidx[j][qp]*dphidx[k][qp];
	  tmp -= 0.5*(_dt/_tNorm)*dB_dmu*phi[j][qp]*d2phidx2[k][qp];

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
  fe->reinit(elem);

  const vector<Real> &JxW = fe->get_JxW();
  const vector<vector<Real> > &phi = fe->get_phi();

  for(unsigned int qp=0; qp<qrule.n_points(); qp++){
    int ir = _realSGH.getPointInvID(iElem, qp);

    double mu_n, dmu_dx_n, d2mu_dx2_n, mu_n1_l1;
    double eEx_n = 0.0;//e*_pdm.get_Ex_n(ir)/_eExNorm;
    double edEx_dx_n = 0.0;//e*_pdm.get_dEx_dx_n(ir)/(_eExNorm/_muNorm);
    double edEx_dx_n1_l1 = 0.0;//e*_pdm.get_dEx_dx_n1_l1(ir)/(_eExNorm/_muNorm);

    if( s == -1 ){
      mu_n = _pdm.get_mue_n(ir)/_muNorm;
      dmu_dx_n = _pdm.get_dmue_dx_n(ir)/_dmudxNorm;
      d2mu_dx2_n = _pdm.get_d2mue_dx2_n(ir)/_d2mudx2Norm;
      mu_n1_l1 = _pdm.get_mue_n1_l1(ir)/_muNorm;
    }
    else {
      mu_n = _pdm.get_muh_n(ir)/_muNorm;
      dmu_dx_n = _pdm.get_dmuh_dx_n(ir)/_dmudxNorm;
      d2mu_dx2_n = _pdm.get_d2muh_dx2_n(ir)/_d2mudx2Norm;
      mu_n1_l1 = _pdm.get_muh_n1_l1(ir)/_muNorm;
    }

    double An = _ab.calcA(_muNorm*mu_n)/_ANorm;
    double Bn = _ab.calcB(_muNorm*mu_n)/_BNorm;
    double Bn1 = _ab.calcB(_muNorm*mu_n1_l1)/_BNorm;
    double C = 0.0;

    C += mu_n;
    C += -0.5*s*(_dt/_tNorm)*An*(eEx_n-s*dmu_dx_n)*dmu_dx_n;
    C += 0.5*(_dt/_tNorm)*Bn*d2mu_dx2_n;
    C += -0.5*s*(_dt/_tNorm)*Bn*edEx_dx_n;
    C += -0.5*s*(_dt/_tNorm)*Bn1*edEx_dx_n1_l1;

    for(unsigned int i=0; i<phi.size(); i++){
      F(i) += JxW[qp]*phi[i][qp]*C;
    }
  }
}
