#include "PoissonDiffusionMediator.h"
#include "ChargeDensity2D.h"
#include <stdio.h>
#include <PhysicalUnits.h>

using namespace std;


/*
 * Constructor.
 */

PoissonDiffusionMediator::
PoissonDiffusionMediator(const RealSpaceGridHandler &realSGH,
			 PoissonSolver2D &poisson,
			 const RealSpaceArrayDiffusion &SigmaDope,
			 const FermiDistrGraphene &fermiDistr):
  _mue_n(realSGH), _mue_n1_l(realSGH), _mue_n1_l1(realSGH),
  _dmue_dx_n(realSGH), _dmue_dx_n1_l(realSGH),
  _dmue_dx_n1_l1(realSGH),
  _d2mue_dx2_n(realSGH), _d2mue_dx2_n1_l(realSGH),
  _d2mue_dx2_n1_l1(realSGH),
  _muh_n(realSGH), _muh_n1_l(realSGH), _muh_n1_l1(realSGH),
  _dmuh_dx_n(realSGH), _dmuh_dx_n1_l(realSGH),
  _dmuh_dx_n1_l1(realSGH),
  _d2muh_dx2_n(realSGH), _d2muh_dx2_n1_l(realSGH),
  _d2muh_dx2_n1_l1(realSGH),
  _Ex_n(realSGH), _Ex_n1_l(realSGH), _Ex_n1_l1(realSGH),
  _dEx_dx_n(realSGH), _dEx_dx_n1_l(realSGH), _dEx_dx_n1_l1(realSGH),
  _SigmaDope(realSGH), _poisson(poisson),
  _nrStepsNI(0),
  _realSGH(realSGH), _fermiDistr(fermiDistr)
{
  for(int i=0; i<_SigmaDope.getSize(); i++){
    _SigmaDope.setAt(i, SigmaDope.getAt(i));
  }
}


/*
 * Set the initial solutions $\mu_{r,0}$, $E_{x,0}$,
 * and $dE_{x,0}/dx$.
 */

void PoissonDiffusionMediator::
setInitialSolutions(const RealSpaceArrayDiffusion &mue0,
		    const RealSpaceArrayDiffusion &muh0,
		    const RealSpaceArrayDiffusion &Ex0,
		    const RealSpaceArrayDiffusion &dEx_dx0)
{
  const int n = _mue_n.getSize();

  for(int i=0; i<n; i++){
    _mue_n.setAt(i, mue0.getAt(i));
    _muh_n.setAt(i, muh0.getAt(i));
    _Ex_n.setAt(i, Ex0.getAt(i));
    _dEx_dx_n.setAt(i, dEx_dx0.getAt(i));
  }

  // Need to calculate $d\mu_{r}/dx$ and $d^{2}\mu_{r}/dx^{2}$
  // from $\mu_{r}$.

  for(int i=0; i<n; i++){
    double mue_ip1, mue_im1, muh_ip1, muh_im1;
    double h_i, h_im1;

    if( i == 0 ){
      mue_im1 = _mue_n.getAt(n-1);
      muh_im1 = _muh_n.getAt(n-1);
      h_im1 = _realSGH.getAt(n-1)-_realSGH.getAt(n-2);
    }
    else {
      mue_im1 = _mue_n.getAt(i-1);
      muh_im1 = _muh_n.getAt(i-1);
      h_im1 = _realSGH.getAt(i)-_realSGH.getAt(i-1);
    }

    if( i == n-1 ){
      mue_ip1 = _mue_n.getAt(0);
      muh_ip1 = _muh_n.getAt(0);
      h_i = _realSGH.getAt(1)-_realSGH.getAt(0);
    }
    else {
      mue_ip1 = _mue_n.getAt(i+1);
      muh_ip1 = _muh_n.getAt(i+1);
      h_i = _realSGH.getAt(i+1)-_realSGH.getAt(i);
    }

    _dmue_dx_n.setAt(i, (mue_ip1-mue_im1)/(h_i+h_im1));
    _dmuh_dx_n.setAt(i, (muh_ip1-muh_im1)/(h_i+h_im1));

    double d2mue_dx2, d2muh_dx2;

    d2mue_dx2 = (mue_ip1-_mue_n.getAt(i))/h_i;
    d2mue_dx2 -= (_mue_n.getAt(i)-mue_im1)/h_im1;
    d2mue_dx2 *= 0.5*(h_i+h_im1);
    _d2mue_dx2_n.setAt(i, d2mue_dx2);

    d2muh_dx2 = (muh_ip1-_muh_n.getAt(i))/h_i;
    d2muh_dx2 -= (_muh_n.getAt(i)-muh_im1)/h_im1;
    d2muh_dx2 *= 0.5*(h_i+h_im1);
    _d2muh_dx2_n.setAt(i, d2muh_dx2);
  }

}


/*
 * Get the current solutions $\mu_{r,n}$, $E_{x,n}$.
 */

void PoissonDiffusionMediator::
getCurrentSolutions(RealSpaceArrayDiffusion &mue,
		    RealSpaceArrayDiffusion &muh,
		    RealSpaceArrayDiffusion &Ex) const
{
  const int n = _mue_n.getSize();

  for(int i=0; i<n; i++){
    mue.setAt(i, _mue_n.getAt(i));
    muh.setAt(i, _muh_n.getAt(i));
    Ex.setAt(i, _Ex_n.getAt(i));
  }
}


/*
 * Set the current time and the increment to the next time step.
 */

void PoissonDiffusionMediator::setTime(double t, double dt)
{
  _t = t;
  _dt = dt;
}


/*
 * Update the solutions to the next time step.
 */

void PoissonDiffusionMediator::updateSolutions()
{
  const int n = _mue_n.getSize();

  for(int i=0; i<n; i++){
    _mue_n.setAt(i, _mue_n1_l1.getAt(i));
    _dmue_dx_n.setAt(i, _dmue_dx_n1_l1.getAt(i));
    _d2mue_dx2_n.setAt(i, _d2mue_dx2_n1_l1.getAt(i));

    _muh_n.setAt(i, _muh_n1_l1.getAt(i));
    _dmuh_dx_n.setAt(i, _dmuh_dx_n1_l1.getAt(i));
    _d2muh_dx2_n.setAt(i, _d2muh_dx2_n1_l1.getAt(i));

    _Ex_n.setAt(i, _Ex_n1_l1.getAt(i));
    _dEx_dx_n.setAt(i, _dEx_dx_n1_l1.getAt(i));
  }
}


/*
 * Set the solutions of the previous time step as the start of
 * the nonlinear iteration: $\mu_{r, n+1}^{(l=0)} <- \mu_{r,n}$.
 */

void PoissonDiffusionMediator::
setInitialSolutionsNI()
{
  const int n = _mue_n.getSize();

  for(int i=0; i<n; i++){
    _mue_n1_l.setAt(i, _mue_n.getAt(i));
    _dmue_dx_n1_l.setAt(i, _dmue_dx_n.getAt(i));
    _d2mue_dx2_n1_l.setAt(i, _d2mue_dx2_n.getAt(i));

    _muh_n1_l.setAt(i, _muh_n.getAt(i));
    _dmuh_dx_n1_l.setAt(i, _dmuh_dx_n.getAt(i));
    _d2muh_dx2_n1_l.setAt(i, _d2muh_dx2_n.getAt(i));

    _Ex_n1_l.setAt(i, _Ex_n.getAt(i));
    _dEx_dx_n1_l.setAt(i, _dEx_dx_n.getAt(i));

    _mue_n1_l1.setAt(i, _mue_n.getAt(i));
    _dmue_dx_n1_l1.setAt(i, _dmue_dx_n.getAt(i));
    _d2mue_dx2_n1_l1.setAt(i, _d2mue_dx2_n.getAt(i));

    _muh_n1_l1.setAt(i, _muh_n.getAt(i));
    _dmuh_dx_n1_l1.setAt(i, _dmuh_dx_n.getAt(i));
    _d2muh_dx2_n1_l1.setAt(i, _d2muh_dx2_n.getAt(i));

    _Ex_n1_l1.setAt(i, _Ex_n.getAt(i));
    _dEx_dx_n1_l1.setAt(i, _dEx_dx_n.getAt(i));
  }

  _nrStepsNI = 0;
}


/*
 * Update the solutions in nonlinear iteration.
 */

void PoissonDiffusionMediator::
updateSolutionsNI(const RealSpaceArrayDiffusion &mue,
		  const RealSpaceArrayDiffusion &dmue_dx,
		  const RealSpaceArrayDiffusion &d2mue_dx2,
		  const RealSpaceArrayDiffusion &muh,
		  const RealSpaceArrayDiffusion &dmuh_dx,
		  const RealSpaceArrayDiffusion &d2muh_dx2)
{

  // Shift the old solutions.

  const int n = _mue_n.getSize();

  for(int i=0; i<n; i++){
    _mue_n1_l.setAt(i, _mue_n1_l1.getAt(i));
    _dmue_dx_n1_l.setAt(i, _dmue_dx_n1_l1.getAt(i));
    _d2mue_dx2_n1_l.setAt(i, _d2mue_dx2_n1_l1.getAt(i));

    _muh_n1_l.setAt(i, _muh_n1_l1.getAt(i));
    _dmuh_dx_n1_l.setAt(i, _dmuh_dx_n1_l1.getAt(i));
    _d2muh_dx2_n1_l.setAt(i, _d2muh_dx2_n1_l1.getAt(i));

    _Ex_n1_l.setAt(i, _Ex_n1_l1.getAt(i));
    _dEx_dx_n1_l.setAt(i, _dEx_dx_n1_l1.getAt(i));
  }


  // Set the current $mu_{r}$ and their derivatives.

  for(int i=0; i<n; i++){
    _mue_n1_l1.setAt(i, mue.getAt(i));
    _dmue_dx_n1_l1.setAt(i, dmue_dx.getAt(i));
    _d2mue_dx2_n1_l1.setAt(i, d2mue_dx2.getAt(i));

    _muh_n1_l1.setAt(i, mue.getAt(i));
    _dmuh_dx_n1_l1.setAt(i, dmuh_dx.getAt(i));
    _d2muh_dx2_n1_l1.setAt(i, d2muh_dx2.getAt(i));
  }


  // Calculate the charge density from $\mu_{r,n+1}^{(l+1)}$.

  RealSpaceArrayDiffusion se(_realSGH), sh(_realSGH);
  ChargeDensity2D rho2D(se, sh, _SigmaDope);

  for(int i=0; i<_realSGH.getSize(); i++){
    se.setAt(i, _fermiDistr.calcConcentrationFermiExact(_mue_n1_l1.getAt(i)));
    sh.setAt(i, _fermiDistr.calcConcentrationFermiExact(_muh_n1_l1.getAt(i)));
  }

  rho2D.updateInterpolator();


  // Calculate $E_{x,n+1}^{(l+1)}$ and $dE_{x,n+1}^{(l+1)}/dx$.

  _poisson.solveAndCalcField2DEG(_t, _Ex_n1_l1, _dEx_dx_n1_l1, rho2D);


  FILE *f_mue, *f_muh, *f_Ex;
  char filename[100];

  sprintf(filename, "dat/mue-t=%05.1ffs-NI=%d.dat", s2fs(_t), _nrStepsNI);
  f_mue = fopen(filename, "w");
  sprintf(filename, "dat/muh-t=%05.1ffs-NI=%d.dat", s2fs(_t), _nrStepsNI);
  f_muh = fopen(filename, "w");
  sprintf(filename, "dat/Ex-t=%05.1ffs-NI=%d.dat", s2fs(_t), _nrStepsNI);
  f_Ex = fopen(filename, "w");

  for(int i=0; i<n; i++){
    fprintf(f_mue, "%g %g\n", _realSGH.getAt(i), J2meV(_mue_n1_l1.getAt(i)));
    fprintf(f_muh, "%g %g\n", _realSGH.getAt(i), J2meV(_muh_n1_l1.getAt(i)));
    fprintf(f_Ex, "%g %g\n", _realSGH.getAt(i), _Ex_n1_l1.getAt(i));
  }

  fclose(f_mue);
  fclose(f_muh);
  fclose(f_Ex);

  _nrStepsNI++;

  if( _nrStepsNI == 10 ) exit(1);
}
