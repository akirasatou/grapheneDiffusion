#include "PoissonDiffusionMediator.h"

using namespace std;


/*
 * Constructor.
 */

PoissonDiffusionMediator::
PoissonDiffusionMediator(const RealSpaceGridHandler &realSGH,
			 PoissonSolver2D &poisson):
  _mue_n(realSGH), _mue_n1_l(realSGH), _mue_n1_l1(realSGH),
  _muh_n(realSGH), _muh_n1_l(realSGH), _muh_n1_l1(realSGH),
  _Ex_n(realSGH), _Ex_n1_l(realSGH), _Ex_n1_l1(realSGH),
  _dEx_dx_n(realSGH), _dEx_dx_n1_l(realSGH), _dEx_dx_n1_l1(realSGH),
  _poisson(poisson)
{
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
    _muh_n.setAt(i, _muh_n1_l1.getAt(i));
    _Ex_n.setAt(i, _Ex_n1_l1.getAt(i));
    _dEx_dx_n.setAt(i, _dEx_dx_n1_l1.getAt(i));
  }
}


/*
 * Set the solutions of the previous time step as the start of
 * the nonlinear iteration: $\mu_{r, n+1}^{(l=0)} <- \mu_{r,n}$.
 */

void PoissonDiffusionMediator::
setInitialSolutionsInNonlinearIteration()
{
  const int n = _mue_n.getSize();

  for(int i=0; i<n; i++){
    _mue_n1_l.setAt(i, _mue_n.getAt(i));
    _muh_n1_l.setAt(i, _muh_n.getAt(i));
    _Ex_n1_l.setAt(i, _Ex_n.getAt(i));
  }
}


/*
 * Update the solutions in nonlinear iteration.
 */

void PoissonDiffusionMediator::
updateSolutionsInNonlinearIteration(const RealSpaceArrayDiffusion &mue,
				    const RealSpaceArrayDiffusion &muh)
{

  // Shift the old solutions.

  const int n = _mue_n.getSize();

  for(int i=0; i<n; i++){
    _mue_n1_l.setAt(i, _mue_n1_l1.getAt(i));
    _muh_n1_l.setAt(i, _muh_n1_l1.getAt(i));
    _Ex_n1_l.setAt(i, _Ex_n1_l1.getAt(i));
    _dEx_dx_n1_l.setAt(i, _dEx_dx_n1_l1.getAt(i));

    _mue_n1_l1.setAt(i, mue.getAt(i));
    _muh_n1_l1.setAt(i, muh.getAt(i));
  }


  // Calculate the charge density from $\mu_{r,n+1}^{(l+1)}$.

  // Calculate $E_{x,n+1}^{(l+1)}$ and $dE_{x,n+1}^{(l+1)}/dx$.

  // rho2D
  //_poisson.solveAndCalcField2DEG(_t, _Ex_n1_l1, _dEx_dx_n1_l1,
  //rho2D);


}
