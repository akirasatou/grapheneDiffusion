#include "DiffusionSolver1D.h"
#include <nonlinear_implicit_system.h>
#include <transient_system.h>

using namespace std;


/*
 * Solve the diffusion equation for the next time step.
 */

void DiffusionSolver1D::solveStep()
{
  TransientNonlinearImplicitSystem &sys = _es.get_system<TransientNonlinearImplicitSystem>(_sysName);


  // Update the time.

  _nrSteps++;
  sys.time = _nrSteps*_difDsc.get_dt();


  // Save the solution of the previous time step.

  *sys.old_local_solution = *sys.current_local_solution;


  // Solve the diffusion equation by nonlinear Newton iteration.

  sys.solve();

}
