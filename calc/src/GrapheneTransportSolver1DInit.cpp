#include "GrapheneTransportSolver1D.h"
#include "ChargeDensity2D.h"
#include <math.h>
#include <PhysicalConstants.h>
#include <PhysicalUnits.h>
#include <iostream>
#include <mpi.h>
#include <assert.h>

using namespace std;


/*
 * Initialization of the field reference points
 * and refinement of the mesh.
 */

void GrapheneTransportSolver1D::_initPoissonSolver()
{
  char filehead[100];

  // Initialize the field calculation function in PoissonSolver2D.

  _poisson.initSamplePoints2DEG(_realSGH);


  // Mesh refinement before the time-step iteration.
  // Ficticious electron distribution is set to generate 
  // a finer mesh around the 2DEG.
  // !Do not set too high concentration because it induces large
  // numerical errors when calculating with the real distribution.!

  const int nr = 2;
  const double dsAmp = 0.1;
  const double SigmaArt = 1e16;

  RealSpaceArrayDiffusion se(_realSGH), sh(_realSGH), sDope(_realSGH);

  for(int n=0; n<nr && _poisson.getNrNodes()<16000; n++){
  
  // Refinement using the first perturbed concentration.

    for(int i=0; i<_realSGH.getSize(); i++){
      double a = (_realSGH.getAt(i)-_realSGH.getXl())/(_realSGH.getXr()-_realSGH.getXl());
      
      se.setAt(i, SigmaArt*(1.0+dsAmp*sin(2*M_PI*a)));
      sh.setAt(i, 0.0);
      sDope.setAt(i, 0.0);
    }
    
    ChargeDensity2D rho2D(se, sh, sDope);
    double t0 = 0.0-fs2s(1);
    
    rho2D.updateInterpolator();
    _poisson.solve(t0, rho2D);
    _refineMesh();
    
    sprintf(filehead, "phiInit-n=%d", 2*n);
    outputPotential(_SSDir.c_str(), filehead, true);

    
    // Refinement using the second perturbed concentration.
    
    if(_poisson.getNrNodes() >= 16000) break;
    
    for(int i=0; i<_realSGH.getSize(); i++){
      double a = (_realSGH.getAt(i)-_realSGH.getXl())/(_realSGH.getXr()-_realSGH.getXl());
      
      se.setAt(i, SigmaArt*(1.0+dsAmp*cos(2*M_PI*a)));
      sh.setAt(i, 0.0);
      sDope.setAt(i, 0.0);
    }
    
    rho2D.updateInterpolator();
    _poisson.solve(t0, rho2D);
    _refineMesh();
    
    sprintf(filehead, "phiInit-n=%d", 2*n+1);
    outputPotential(_SSDir.c_str(), filehead, true);
  }
}
