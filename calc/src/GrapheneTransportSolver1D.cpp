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
  _diffusion(_difDsc, _poisson, _poiDsc.getXlPoisson(), _poiDsc.getXrPoisson()),
  _realSGH(_diffusion.getRealSGH()),
/*
  _gridParam(_bltDsc.getNx(), _bltDsc.getNeps(), _bltDsc.getNtheta(),
	     _bltDsc.get_dt(), MPI::COMM_WORLD),
  _Efhd(_calcEfHeavilyDoped(_bltDsc.getSigma0hd(), _bltDsc.toAccountHole())),
  _Ex(_realSGH), _SigmaElectron(_realSGH), _SigmaHole(_realSGH),
  _Sigma0xi(_realSGH),
*/
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


  /*
  _initPoissonSolver();
  */

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


  // Initial SCF distribution.

  //_setInitialSteadyStateSCF();

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
  /*
  double t = getTime();
  double dt = _gridParam.get_dt();


  // Iteration of updates.

  for(int num=0; num<3; num++){
    _bltElectron.calcConcentration(_SigmaElectron);
    _SigmaElectron.updateInterpolator();

    _bltHole.calcConcentration(_SigmaHole);
    _SigmaHole.updateInterpolator();

    _solveAndCalcField2DEG(_Ex, _SigmaElectron, _SigmaHole,
			   _Sigma0xi);

    //for(int i=0; i<_Ex.getSize(); i++) _Ex.setAt(i, 0.0);

    _bltElectron.update(num, _Ex);
    _bltHole.update(num, _Ex);
  }
  */


  // Output.
  /*
  int n_output_step = (int)round(_bltDsc.get_tOutputStep()/dt);
  int n_outputBin_step = (int)round(_bltDsc.get_tOutputBinStep()/dt);
  char filehead[300];

  if(_nSteps%n_outputBin_step == 0){
    if( _bltDsc.toOutputConcentration() ){
      sprintf(filehead, "conc-t=%04.0ffs", s2fs(t));
      outputConcentration2DEGBin(_concDir.c_str(), filehead);
    }
  }
  
  if(_nSteps%n_output_step == 0){
    if( _bltDsc.toOutputConcentration() ){
      sprintf(filehead, "conc-t=%04.0ffs", s2fs(t));
      outputConcentration2DEG(_concDir.c_str(), filehead);
    }
    if( _bltDsc.toOutputPotential2D() ){
      sprintf(filehead, "phi-t=%04.0ffs", s2fs(t));
      outputPotential(_pot2DDir.c_str(), filehead);
    }
    if( _bltDsc.toOutputField() ){
      sprintf(filehead, "field-t=%04.0ffs", s2fs(t));
      outputField2DEG(_fieldDir.c_str(), filehead);
    }
    if( _bltDsc.toOutputVelocity() ){
      sprintf(filehead, "vel-t=%04.0ffs", s2fs(t));
      outputVelocity(_velDir.c_str(), filehead);
    }
  }

  MPI_Barrier(_gridParam.getWorld());
  */

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
/*
void GrapheneTransportSolver1D::_refineMesh()
{
  _poisson.refineMesh();
  _poisson.initSamplePoints2DEG(_realSGH.getXiPtrPoisson(),
				_realSGH.getSizePoisson());
}
*/

/*
 * Solve the Poisson equation and calculate the field in 2DEG.
 */
/*
void GrapheneTransportSolver1D::
_solveAndCalcField2DEG(Field &Ex, const Concentration &SigmaElectron,
		       const Concentration &SigmaHole,
		       const Concentration &Sigma0)
{
  // Initialization. Do not delete this because the field
  // at the fictitious highly-doped regions is set here.

  for(int i=0; i<Ex.getSizeBltExtended(); i++){
    Ex.setAtBltExtended(i, 0.0);
  }


  // Gradual channel approximation.

  if(_bltDsc.getFieldModel() == GCA){

    // With bottom gates is not supported.
    assert(_poiDsc.getBottomGates().size() == 0);

    double epsAvg = (_poiDsc.getEpsilonUpper2DEG()+_poiDsc.getEpsilonLower2DEG())/2.0;
    double a = e*_poiDsc.getWg()/epsAvg;
    
    for(int i=0; i<Ex.getSize(); i++){
      double tmp = SigmaElectron.getDerivativeAt(i);

      tmp -= SigmaHole.getDerivativeAt(i);
      Ex.setAt(i, a*tmp);
    }
  }

  // Self-consistent field calculation.

  else if(_bltDsc.getFieldModel() == SCF){
    ChargeDensity2D rho2D(_SigmaElectron, _SigmaHole, _Sigma0xi);

    //_poisson.solveAndCalcField2DEG(_Ex, rho2D, true);
    _poisson.solveAndCalcField2DEG(_Ex, rho2D);

    //for(int i=0; i<_Ex.getSize(); i++) _Ex.setAt(i, 0.0);
  }
}
 */

/*
 * Solve the Poisson equation and calculate the potential in 2DEG 
 * at the reference points and put it into pot.
 */
/*
void GrapheneTransportSolver1D::
_solveAndCalcPotential2DEG(Potential &pot,
			   const Concentration &SigmaElectron,
			   const Concentration &SigmaHole,
			   const Concentration &Sigma0)
{
  // Initialization. Do not delete this because the field
  // at the fictitious highly-doped regions is set here.

  for(int i=0; i<pot.getSizeBltExtended(); i++){
    pot.setAtBltExtended(i, 0.0);
  }


  // Gradual channel approximation.

  if(_bltDsc.getFieldModel() == GCA){
    double epsAvg = (_poiDsc.getEpsilonUpper2DEG()+_poiDsc.getEpsilonLower2DEG())/2.0;
    double a = e*_poiDsc.getWg()/epsAvg;
    
    for(int i=0; i<pot.getSize(); i++){
      double tmp = SigmaElectron.getAt(i);

      tmp -= SigmaHole.getAt(i);
      tmp -= Sigma0.getAt(i);
      pot.setAt(i, -a*tmp);
    }
  }

  // Self-consistent field calculation.

  else if(_bltDsc.getFieldModel() == SCF){
    ChargeDensity2D rho2D(_SigmaElectron, _SigmaHole, _Sigma0xi);

    _poisson.solveAndCalcPotential2DEG(pot, rho2D, true);
  }

}
*/
