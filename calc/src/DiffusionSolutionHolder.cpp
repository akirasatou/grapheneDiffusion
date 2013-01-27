#include "DiffusionSolutionHolder.h"

using namespace std;


/*
 * Constructor.
 */

DiffusionSolutionHolder::
DiffusionSolutionHolder(const RealSpaceGridHandler &realSGH):
  _mue_n(realSGH), _muh_n(realSGH), _Ex_n(realSGH),
  _dEx_dx_n(realSGH),
  _mue_n1_l(realSGH), _muh_n1_l(realSGH), _Ex_n1_l(realSGH),
  _dEx_dx_n1_l(realSGH)
{
}


/*
 * Set the initial solutions $\mu_{r,0}$, $E_{x,0}$,
 * and $dE_{x,0}/dx$.
 */

void DiffusionSolutionHolder::
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

void DiffusionSolutionHolder::
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
