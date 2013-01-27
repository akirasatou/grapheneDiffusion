#include "GrapheneTransportSolver1D.h"
#include "ChargeDensity2D.h"
#include <PhysicalConstants.h>
#include <PhysicalUnits.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include <Transistor2D/ChargeDensityCalculator.h>

using namespace std;


/*
 * Find the SCF steady-state distributions of the concentrations as
 * the initial distributions, using the Newton method provided by 
 * libmesh for solving the nonlinear Poisson equation.
 */

void GrapheneTransportSolver1D::_setInitialSteadyStateSCF()
{
  char filehead[100];
  const int nrNodesMax = 7000;

  // Ferm-level shift due to the doping.

  RealSpaceArrayDiffusion Ef0(_realSGH);


  for(int i=0; i<Ef0.getSize(); i++){
    Ef0.setAt(i, _fermiDistr.calcEfEHExact(_SigmaDope.getAt(i)));
  }


  Ef0.updateInterpolator();


  // Initialize the nonlinear solver in PoissonSolver2D.

  ChargeDensityCalculatorGraphene cdcg(_difDsc.getT(),
				       _SigmaDope, Ef0);
  double tol = 1e-3;

  _poisson.initNSS(cdcg, tol);


  // Solver the nonlinear Poisson.

  RealSpaceArrayDiffusion pot(_realSGH), se(_realSGH), sh(_realSGH);
  ChargeDensity2D rho2D(se, sh, _SigmaDope);

  _poisson.solveNSSAndCalcPotential2DEG(pot);

  for(int i=0; i<_realSGH.getSize(); i++){
    double E = Ef0.getAt(i)+e*pot.getAt(i);
    
    se.setAt(i, _fermiDistr.calcConcentrationFermiExact(E));
    sh.setAt(i, _fermiDistr.calcConcentrationFermiExact(-E));

    _SigmaElectron.setAt(i, se.getAt(i));
    _SigmaHole.setAt(i, sh.getAt(i));
    _muElectron.setAt(i, E);
    _muHole.setAt(i, E);
  }


  // Solve the linear Poisson once to have the initial solution vector
  // for the nonlinear Poisson. The time is slightly shifted from 0
  // to deal with abrupt gate-voltage off: $V_{g0}\theta(-t)$.

  RealSpaceArrayDiffusion Ex(_realSGH), dEx_dx(_realSGH);

  for(int i=0; i<Ex.getSize(); i++){
    Ex.setAt(i, 0.0);
    dEx_dx.setAt(i, 0.0);
  }

  rho2D.updateInterpolator();
  _poisson.solveAndCalcField2DEG(_difDsc.get_dt()/1000, Ex, dEx_dx,
				 rho2D);


  // Set the initial solutions to the holder.

  _dsh.setInitialSolutions(_muElectron, _muHole, Ex, dEx_dx);


  // Output the field and concentration.

  outputField2DEG(_SSDir.c_str(), "Ex-SS", Ex);
  outputConcentration2DEG(_SSDir.c_str(), "conc-SS", se, sh);
  outputFermiLevel2DEG(_SSDir.c_str(), "mu-SS");

  MPI_Barrier(MPI_COMM_WORLD);
}
