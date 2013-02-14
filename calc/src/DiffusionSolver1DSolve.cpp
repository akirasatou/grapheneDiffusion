#include "DiffusionSolver1D.h"
#include <nonlinear_implicit_system.h>
#include <transient_system.h>

using namespace std;


/*
 * Solve the diffusion equation for the next time step.
 */

void DiffusionSolver1D::solveStep(double t, double dt)
{
  TransientNonlinearImplicitSystem &sys = _es.get_system<TransientNonlinearImplicitSystem>(_sysName);


  // Update the time.

  sys.time = t;
  _pdm.setTime(t, dt);
  _rj->setTime(t, dt);


  // Save the solution of the previous time step.

  *sys.old_local_solution = *sys.current_local_solution;


  // Set the solutions of the previous time step as the start of
  // the nonlinear iteration: $\mu_{r, n+1}^{(l=0)} <- \mu_{r,n}$.

  _pdm.setInitialSolutionsInNonlinearIteration();


  // Solve the diffusion equation by nonlinear Newton iteration.

  // sys.solve();


  // Update the solutions to the next time step.

  _pdm.updateSolutions();
  _nSteps++;
}
