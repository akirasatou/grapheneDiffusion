#pragma once
#include "DiffusionSolver1DDescriptor.h"
/*
#include "DiffusionSolver1D.h"
#include "Field.h"
#include "Concentration.h"
#include "Potential.h"
#include "RealSpaceGridHandler.h"
#include "PhaseSpaceGridHandler.h"
#include "FermiDistrGraphene.h"
*/
#include <PhysicalConstants.h>
#include <PhysicalUnits.h>
#include <Transistor2D/PoissonSolver2D.h>
#include <OutputDirectoryManager/OutputDirectoryManager.h>
//#include "GridParameter.h"

using namespace std;


/*
 * Class solving the Diffusion equation with self-consistent 
 * Poisson equation. Kind of facade" DiffusionnSolver1D and
 * PoissonSolver2D, keeping the algorithm of how they are used in
 * cooperation and providing output methods.
 */

class GrapheneTransportSolver1D
{
  
public:

  GrapheneTransportSolver1D(const PoissonSolver2DDescriptor &poiDsc,
			    const DiffusionSolver1DDescriptor &difDsc,
			    const OutputDirectoryManager &odm);
  
  void solveStep();
  double getTime() const;

  //double calcTotalConcentrationDiff() const;
  //double getSigma0() const;
  //inline int getNx() const { return _realSGH.getSizeBlt(); }

  void outputPotential(const char *dname, const char *fhead,
		       bool toOutputMesh=false, 
		       const std::string &format="") const;
  void outputConcentration2DEG(const char *dname, const char *fhead) const;

  /*
  void outputConcentration2DEG(const char *dname, const char *fhead,
			       const Concentration &se,
			       const Concentration &sh) const;
  */
  void outputConcentration2DEGBin(const char *dname, const char *fhead) const;

  void outputPotential2DEG(const char *dname, const char *fhead) const;
  /*
  void outputPotential2DEG(const char *dname, const char *fhead,
			   const Potential &pot) const;
  */
  void outputField2DEG(const char *dname, const char *fhead) const;
  /*
  void outputField2DEG(const char *dname, const char *fhead,
		       const Field &field) const;
  */
  void outputVelocity(const char *dname, const char *fhead);


private:

  const PoissonSolver2DDescriptor &_poiDsc;
  const DiffusionSolver1DDescriptor &_difDsc;
  PoissonSolver2D _poisson;
  /*
  GridParameter _gridParam;
  RealSpaceGridHandler _realSGH;
  PhaseSpaceGridHandler _phaseSGH;
  FermiDistrGraphene _fermiDistr;
  DiffusionSolver1D _bltElectron, _bltHole;
  Field _Ex;
  Concentration _SigmaElectron, _SigmaHole, _Sigma0xi;
  */
  int _nSteps;

  OutputDirectoryManager _odm;
  std::string _SSDir, _concDir, _pot2DDir, _fieldDir, _velDir;

  void _setDopingProfile();
  void _initPoissonSolver();

  /*
  void _refineMesh();
  void _setInitialSteadyStateLocalFermi();
  void _setInitialSteadyStateSCFFermi(double rVg=1.0,
				      double rVgPrev=-1.0);
  void _setInitialSteadyStateSCFFermiNewton();

  void _updateConcentrations(Concentration &se, Concentration &sh, 
			     Concentration &sePrev,
			     Concentration &shPrev,
			     const Potential &pot, double Ef0, 
			     double maxChangeFraction) const;
  void _solveAndCalcField2DEG(Field &Ex,
			      const Concentration &SigmaElectron,
			      const Concentration &SigmaHole,
			      const Concentration &Sigma0);
  void _solveAndCalcPotential2DEG(Potential &pot,
				  const Concentration &SigmaElectron,
				  const Concentration &SigmaHole,
				  const Concentration &Sigma0);
  void _calcInitialGuess(Concentration &se, Concentration &sh,
			 double Ef0) const;
  void _refineMeshSS(int nrNodesMax, bool toOutputMesh);
  */
};
