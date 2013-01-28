#include "GrapheneTransportSolver1D.h"
//#include "ChargeDensity2D.h"
#include <math.h>
#include <assert.h>
#include <PhysicalConstants.h>
#include <PhysicalUnits.h>
#include <iostream>

using namespace std;


/*
 * Constructor.
 */

GrapheneTransportSolver1D::
GrapheneTransportSolver1D(const PoissonSolver2DDescriptor &poiDsc,
			  const DiffusionSolver1DDescriptor &difDsc,
			  const OutputDirectoryManager &odm):
  _poiDsc(poiDsc), _difDsc(difDsc), _poisson(poiDsc),
  _diffusion(_difDsc, _ab, _dsh, _poiDsc.getXlPoisson(), _poiDsc.getXrPoisson()),
  _realSGH(_diffusion.getRealSGH()),
  _ab(_difDsc.getT(), _difDsc.get_alpha()), _dsh(_realSGH, _poisson),
  _SigmaElectron(_realSGH), _SigmaHole(_realSGH),
  _SigmaDope(_realSGH), _Ex(_realSGH),
  _muElectron(_realSGH), _muHole(_realSGH),
  _fermiDistr(_difDsc.getT()), 
  _nSteps(0), _odm(odm)
{

  // Check the validity of descriptors.

  if( _poiDsc.getTypeBCLeft() != PoissonSolver2DBC::BCPeriodic ){
    cerr << "GrapheneTransportSolver1D::GrapheneTransportSolver1D: ";
    cerr << "BCs must be periodic for Poisson." << endl;
    exit(1);
  }

  if( fabs(_poiDsc.getL()-_difDsc.getLc()) > nm2m(0.01) ){
    cerr << "GrapheneTransportSolver1D::GrapheneTransportSolver1D: ";
    cerr << "channel lengths for Poisson and Diffusion are ";
    cerr << "different." << endl;
    exit(1);
  }


  // Create directories for output.

  if( _difDsc.toOutputSS() ){
    _odm.push(string("SS")); _SSDir = _odm.pop();
  }
  if( _difDsc.toOutputConcentration() ){
    _odm.push(string("conc")); _concDir = _odm.pop();
  }
  if( _difDsc.toOutputPotential2D() ){
    _odm.push(string("pot2D")); _pot2DDir = _odm.pop();
  }
  if( _difDsc.toOutputField() ){
    _odm.push(string("field")); _fieldDir = _odm.pop();
  }
  if( _difDsc.toOutputVelocity() ){
    _odm.push(string("velocity")); _velDir = _odm.pop();
  }


  // Set uniform doping concentration.

  for(int i=0; i<_SigmaDope.getSize(); i++){
    _SigmaDope.setAt(i, _difDsc.getSigma0());
  }


  // Initialize PoissonSolver.

  _initPoissonSolver();


  // Initial SCF distribution.

  _setInitialSteadyStateSCF();

}

/*
double GrapheneTransportSolver1D::
_calcEfHeavilyDoped(double Sigma0hd, bool toAccountHole)
{
  if( toAccountHole ){
    return _fermiDistr.calcEfEHWithGrid(Sigma0hd);
  }
  return _fermiDistr.calcEfWithGrid(Sigma0hd);
}
*/

/*
 * Solve the diffusion equation for the next time step.
 */

void GrapheneTransportSolver1D::solveStep()
{
  // Solve.

  double t = getTime();

  _diffusion.solveStep(t, _difDsc.get_dt());

  cerr << t << endl;
  for(int i=0; i<_poiDsc.getGates().size(); i++){
    cerr << _poiDsc.getGates()[i]->getVoltage(t) << " ";
  }
  cerr << endl;


  // Get the current solutions from the holder.

  _dsh.getCurrentSolutions(_muElectron, _muHole, _Ex);


  // Output.

  int n_output_step = (int)round(_difDsc.get_tOutputStep()/_difDsc.get_dt());
  int n_outputBin_step = (int)round(_difDsc.get_tOutputBinStep()/_difDsc.get_dt());
  char filehead[300];

  /*
  if(_nSteps%n_outputBin_step == 0){
    if( _difDsc.toOutputConcentration() ){
      sprintf(filehead, "conc-t=%04.0ffs", s2fs(t));
      outputConcentration2DEGBin(_concDir.c_str(), filehead);
    }
  }
  */
  
  if(_nSteps%n_output_step == 0){
    if( _difDsc.toOutputFermiLevel() ){
      sprintf(filehead, "mu-t=%04.0ffs", s2fs(t));
      outputFermiLevel2DEG(_concDir.c_str(), filehead);
    }
    if( _difDsc.toOutputConcentration() ){
      sprintf(filehead, "conc-t=%04.0ffs", s2fs(t));
      outputConcentration2DEG(_concDir.c_str(), filehead);
    }
    if( _difDsc.toOutputPotential2D() ){
      sprintf(filehead, "phi-t=%04.0ffs", s2fs(t));
      outputPotential(_pot2DDir.c_str(), filehead);
    }
    if( _difDsc.toOutputField() ){
      sprintf(filehead, "field-t=%04.0ffs", s2fs(t));
      outputField2DEG(_fieldDir.c_str(), filehead);
    }
    if( _difDsc.toOutputVelocity() ){
      sprintf(filehead, "vel-t=%04.0ffs", s2fs(t));
      outputVelocity(_velDir.c_str(), filehead);
    }
  }

  //MPI_Barrier(_gridParam.getWorld());


  // Increment the number of steps.

  _nSteps++;
}


/*
 * Return the current time.
 */

double GrapheneTransportSolver1D::getTime() const
{
  return _nSteps*_difDsc.get_dt();
}


/*
 * Refine the mesh in PoissonSolver2D and reset reference points
 * of the field.
 */

void GrapheneTransportSolver1D::_refineMesh()
{
  _poisson.refineMesh();
  _poisson.initSamplePoints2DEG(_realSGH);
}
