#include "DiffusionSolver1DFD.h"
#include "RealSpaceArrayDiffusion.h"
#include <PhysicalUnits.h>
#include <algorithm>

using namespace std;


/*
 * Solve the diffusion equation for the next time step.
 */

double DiffusionSolver1DFD::solveStep(double t)
{

  // Calculate the adaptive time step to reduce the numerical error.

  double dt = _calcNextTimeStep(t);


  // Set the solutions of the previous time step as the start of
  // the nonlinear iteration: $\mu_{r, n+1}^{(l=0)} <- \mu_{r,n}$.

  _pdm.setInitialSolutionsNI();


  // Solve the diffusion equation by nonlinear Newton iteration.

  RealSpaceArrayDiffusion mue(getRealSGH()), muh(getRealSGH());
  const int size = mue.getSize(), lMax = 100;
  const double tol = 1e-10;
  int l;

  //cerr << "nonlinear iteration" << endl;

  for(l=1; l<=lMax; l++){
    double err = 0.0;

    // For electrons.

    for(int i=0; i<size; i++) _U[i] = _pdm.get_mue_n1_l1(i)/_muNorm;
    _calcRJ(_U, dt, -1);
    _J.solveGauss(_dU, _R);
    for(int i=0; i<size; i++) mue.setAt(i, _muNorm*(_U[i]-_dU[i]));
    for(int i=0; i<size; i++) err = max(err, fabs(_dU[i]));


    // For holes.

    for(int i=0; i<size; i++) _U[i] = _pdm.get_muh_n1_l1(i)/_muNorm;
    _calcRJ(_U, dt, +1);
    _J.solveGauss(_dU, _R);
    for(int i=0; i<size; i++) muh.setAt(i, _muNorm*(_U[i]-_dU[i]));
    for(int i=0; i<size; i++) err = max(err, fabs(_dU[i]));


    // Update the solutions and calculate the new field.

    _pdm.updateSolutionsNI(mue, muh);


    // Convergence check.
    /*
    cerr << "l=" << l << endl;
    cerr << "err=" << err << endl;
    */
    if( err < tol ) break;
  }

  if( l > lMax ){
    cerr << _className << "::solveStep: nonlinear iteration did not ";
    cerr << "converge: t=" << t << endl;
    exit(1);
  }


  // Update the solutions to the next time step.

  _pdm.updateSolutions();

  return dt;
}


/*
 * Calculate the residual $R$ and Jacobian $J$.
 */

void DiffusionSolver1DFD::_calcRJ(const Vector &U, double dt, int sr)
{
  const int size = U.size();

  // R = KU-F.

  for(int i=0; i<size; i++){
    _R[i] = -_calcF(i, dt, sr);

    for(int j=i-1; j<=i+1; j++){
      _R[i] += _calcK(i, j, U, dt, sr)*U[_get_jd(j, size)];
    }
  }


  // J = K+(\partial K/\partial U)U.

  _J.setZeros();

  for(int i=0; i<size; i++){
    for(int j=i-1; j<=i+1; j++){
      _J.setAt(i, _get_jd(j, size), _calcK(i, j, U, dt, sr));
      _J.addAt(i, _get_jd(j, size), _calc_dKdU_U(i, j, U, dt, sr));
    }
  }
}


/*
 * Circular index.
 */

int DiffusionSolver1DFD::_get_jd(int j, int size)
{
  if(j < 0) return size+j;
  if(j >= size) return j-size;
  return j;
}


/*
 * Calculate F.
 * Be cautious about the normalization.
 */

double DiffusionSolver1DFD::_calcF(int i, double dt, int sr) const
{
  double mu_n, dmu_dx_n, d2mu_dx2_n, mu_n1_l1;
  double eEx_n = e*_pdm.get_Ex_n(i)/_eExNorm;
  double edEx_dx_n = e*_pdm.get_dEx_dx_n(i)/(_eExNorm/_xNorm);
  double edEx_dx_n1_l1 = e*_pdm.get_dEx_dx_n1_l1(i)/(_eExNorm/_xNorm);

  if( sr == -1 ){
    mu_n = _pdm.get_mue_n(i)/_muNorm;
    dmu_dx_n = _pdm.get_dmue_dx_n(i)/_dmudxNorm;
    d2mu_dx2_n = _pdm.get_d2mue_dx2_n(i)/_d2mudx2Norm;
    mu_n1_l1 = _pdm.get_mue_n1_l1(i)/_muNorm;
  }
  else {
    mu_n = _pdm.get_muh_n(i)/_muNorm;
    dmu_dx_n = _pdm.get_dmuh_dx_n(i)/_dmudxNorm;
    d2mu_dx2_n = _pdm.get_d2muh_dx2_n(i)/_d2mudx2Norm;
    mu_n1_l1 = _pdm.get_muh_n1_l1(i)/_muNorm;
  }

  double An = _ab.calcA(_muNorm*mu_n)/_ANorm;
  double Bn = _ab.calcB(_muNorm*mu_n)/_BNorm;
  double Bn1 = _ab.calcB(_muNorm*mu_n1_l1)/_BNorm;
  double C = 0.0;

  C += mu_n;
  C += -0.5*sr*(dt/_tNorm)*An*(eEx_n-sr*dmu_dx_n)*dmu_dx_n;
  C += 0.5*(dt/_tNorm)*Bn*d2mu_dx2_n;
  C += -0.5*sr*(dt/_tNorm)*Bn*edEx_dx_n;
  C += -0.5*sr*(dt/_tNorm)*Bn1*edEx_dx_n1_l1;

  return (_dx/_xNorm)*(_dx/_xNorm)*C;
}


/*
 * Calculate K.
 * Be cautious about the normalization.
 */

double DiffusionSolver1DFD::
_calcK(int i, int j, const Vector &U, double dt, int sr) const
{
  if( abs(i-j) > 1 ) return 0.0;

  const int size = U.size();
  double mu, mu_ip1, mu_im1;
  double eEx = e*_pdm.get_Ex_n1_l1(i)/_eExNorm;

  if( sr == -1 ){
    mu = _pdm.get_mue_n1_l1(i)/_muNorm;
    mu_ip1 = _pdm.get_mue_n1_l1(_get_jd(i+1, size))/_muNorm;
    mu_im1 = _pdm.get_mue_n1_l1(_get_jd(i-1, size))/_muNorm;
  }
  else {
    mu = _pdm.get_muh_n1_l1(i)/_muNorm;
    mu_ip1 = _pdm.get_muh_n1_l1(_get_jd(i+1, size))/_muNorm;
    mu_im1 = _pdm.get_muh_n1_l1(_get_jd(i-1, size))/_muNorm;
  }
  
  double A = _ab.calcA(_muNorm*mu)/_ANorm;
  double B = _ab.calcB(_muNorm*mu)/_BNorm;

  if( i == j ){
    return (_dx/_xNorm)*(_dx/_xNorm)+(dt/_tNorm)*B;
  }

  double r;
  int sd = j-i;

  r = sd*0.5*sr*0.5*(dt/_tNorm)*(_dx/_xNorm)*A;
  r *= (_dx/_xNorm)*eEx-0.5*sr*(mu_ip1-mu_im1);
  r += -0.5*(dt/_tNorm)*B;

  return r;
}


/*
 * Calculate \frac{\partial K}{\partial U} U.
 * Be cautious about the normalization.
 */

double DiffusionSolver1DFD::
_calc_dKdU_U(int i, int j, const Vector &U, double dt, int sr) const
{
  if( abs(i-j) > 1 ) return 0.0;

  const int size = U.size();
  double mu, mu_ip1, mu_im1;
  double eEx = e*_pdm.get_Ex_n1_l1(i)/_eExNorm;

  if( sr == -1 ){
    mu = _pdm.get_mue_n1_l1(i)/_muNorm;
    mu_ip1 = _pdm.get_mue_n1_l1(_get_jd(i+1, size))/_muNorm;
    mu_im1 = _pdm.get_mue_n1_l1(_get_jd(i-1, size))/_muNorm;
  }
  else {
    mu = _pdm.get_muh_n1_l1(i)/_muNorm;
    mu_ip1 = _pdm.get_muh_n1_l1(_get_jd(i+1, size))/_muNorm;
    mu_im1 = _pdm.get_muh_n1_l1(_get_jd(i-1, size))/_muNorm;
  }
  
  double A = _ab.calcA(_muNorm*mu)/_ANorm;
  double dA_dmu = _ab.calc_dA_dmu(_muNorm*mu)/(_ANorm/_muNorm);
  double dB_dmu = _ab.calc_dB_dmu(_muNorm*mu)/(_BNorm/_muNorm);

  double dKdU[3];
  double dKdU_U = 0.0;

  if( j == i ){
    dKdU[1] = (dt/_tNorm)*dB_dmu;
    
    for(int s=-1; s<=1; s+=2){
      dKdU[s+1] = s*0.5*sr*0.5*(dt/_tNorm)*(_dx/_xNorm)*dA_dmu;
      dKdU[s+1] *= (_dx/_xNorm)*eEx-0.5*sr*(mu_ip1-mu_im1);
      dKdU[s+1] += -0.5*(dt/_tNorm)*dB_dmu;
    }
  }
  else {
    dKdU[1] = 0.0;

    int sd = j-i;

    for(int s=-1; s<=1; s+=2){
      dKdU[s+1] = -s*sd*0.25*0.5*(dt/_tNorm)*(_dx/_xNorm)*A;
    }
  }

  for(int k=i-1; k<=i+1; k++){
    dKdU_U += dKdU[k-(i-1)]*U[_get_jd(k, size)];
  }

  return dKdU_U;
}


/*
 * Calculate the time step to reduce the numerical error.
 */

double DiffusionSolver1DFD::_calcNextTimeStep(double t) const
{
  // Maximum allowed change in $\mu$.

  double dmuMax = meV2J(0.1);


  // Find the maximum $F$ in the diffusion equation.

  double Fmax = 0.0;

  for(int i=0; i<_difDsc.getNx(); i++){
    double Fe = fabs(_calcFInDiffusionEq(i, -1));
    double Fh = fabs(_calcFInDiffusionEq(i, +1));

    Fmax = max(Fmax, max(Fe, Fh));
  }


  // Set the upper bound on the time step.

  return  min(dmuMax/Fmax, fs2s(100));
}


/*
 * Calculate $F$ in the diffusion equation.
 */

double DiffusionSolver1DFD::_calcFInDiffusionEq(int i, int sr) const
{
  double mu_n, dmu_dx_n, d2mu_dx2_n;
  double eEx_n = e*_pdm.get_Ex_n(i);
  double edEx_dx_n = e*_pdm.get_dEx_dx_n(i);

  if( sr == -1 ){
    mu_n = _pdm.get_mue_n(i);
    dmu_dx_n = _pdm.get_dmue_dx_n(i);
    d2mu_dx2_n = _pdm.get_d2mue_dx2_n(i);
  }
  else {
    mu_n = _pdm.get_muh_n(i);
    dmu_dx_n = _pdm.get_dmuh_dx_n(i);
    d2mu_dx2_n = _pdm.get_d2muh_dx2_n(i);
  }

  double An = _ab.calcA(mu_n);
  double Bn = _ab.calcB(mu_n);
  double r = 0.0;

  r += -sr*_ab.calcA(mu_n)*(eEx_n-sr*dmu_dx_n)*dmu_dx_n;
  r += _ab.calcB(mu_n)*d2mu_dx2_n-sr*Bn*edEx_dx_n;

  return r;
}
