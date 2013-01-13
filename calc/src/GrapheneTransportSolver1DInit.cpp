#include "GrapheneTransportSolver1D.h"
//#include "ChargeDensity2D.h"
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
  /*
  _poisson.initSamplePoints2DEG(_realSGH.getXiPtrPoisson(),
				_realSGH.getSizePoisson());


  if(_bltDsc.getFieldModel() == SCF){

    // Mesh refinement before the time-step iteration.
    // Ficticious electron distribution is set to generate 
    // a finer mesh around the 2DEG.
    // !Do not set too high concentration because it induces large
    // numerical errors when calculating with the real distribution.!
    
    //const int nr = 10;
    const int nr = 2;
    const double dsAmp = 0.5;
    const double Sigma0min = 1e15;
    Concentration se(_realSGH), sh(_realSGH);

    for(int n=0; n<nr && _poisson.getNrNodes()<16000; n++){

      // Refinement using the first perturbed concentration.

      for(int i=-1; i<_realSGH.getSizeBlt()+1; i++){
	double a = (_realSGH.getXiBlt(i)-_realSGH.getXlBlt())/(_realSGH.getXrBlt()-_realSGH.getXlBlt());
	double ds = 0.0;
	
	if(_bltDsc.getTypeBCLeft() == BCPeriodic){
	  ds = dsAmp*sin(2*M_PI*a);
	}
	else {
	  ds = dsAmp*sin(M_PI*(0.9*a+0.05));
	}
	
	if(fabs(_Sigma0xi.getAtBlt(i)) < Sigma0min){
	  ds *= Sigma0min;
	}
	else {
	  ds *= _Sigma0xi.getAtBlt(i);
	}
	
	se.setAtBlt(i, _Sigma0xi.getAtBlt(i)+ds);
	sh.setAtBlt(i, 0.0);
      }

      ChargeDensity2D rho2D(se, sh, _Sigma0xi);

      rho2D.updateInterpolator();
      _poisson.solve(rho2D);
      _refineMesh();

      sprintf(filehead, "phiInit-n=%d", 2*n);
      outputPotential(".", filehead, true);


      // Refinement using the second perturbed concentration.
    
      if(_poisson.getNrNodes() >= 16000) break;

      for(int i=-1; i<_realSGH.getSizeBlt()+1; i++){
	double a = (_realSGH.getXiBlt(i)-_realSGH.getXlBlt())/(_realSGH.getXrBlt()-_realSGH.getXlBlt());
	double ds = 0.0;
	
	if(_bltDsc.getTypeBCLeft() == BCPeriodic){
	  ds = dsAmp*cos(2*M_PI*a);
	}
	else {
	  ds = dsAmp*cos(M_PI*(0.9*a+0.05));
	}
	
	if(fabs(_Sigma0xi.getAtBlt(i)) < Sigma0min){
	  ds *= Sigma0min;
	}
	else {
	  ds *= _Sigma0xi.getAtBlt(i);
	}
	
	se.setAtBlt(i, _Sigma0xi.getAtBlt(i)+ds);
	sh.setAtBlt(i, 0.0);
      }
      
      rho2D.updateInterpolator();
      _poisson.solve(rho2D);
      _refineMesh();

      sprintf(filehead, "phiInit-n=%d", 2*n+1);
      outputPotential(".", filehead, true);
    }
  }
  */
}
