#pragma once

#include "DiffusionSolver1DFD.h"
#include "RealSpaceGridHandler.h"
#include "RealSpaceArrayDiffusion.h"
#include "Concentration.h"
#include "FermiDistrGraphene.h"
#include "DiffusionABCalculator.h"
#include "PoissonDiffusionMediator.h"
/*
#include "Field.h"
#include "Concentration.h"
#include "Potential.h"
*/
#include <PhysicalConstants.h>
#include <PhysicalUnits.h>
#include <Transistor2D/PoissonSolver2D.h>
#include <OutputDirectoryManager/OutputDirectoryManager.h>
//#include "GridParameter.h"


/*
 * Class solving the Diffusion equation with self-consistent 
 * Poisson equation. Kind of facade" DiffusionSolver1D and
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

  void outputConcentration2DEG(const char *dname, const char *fhead,
			       const RealSpaceArrayDiffusion &se,
			       const RealSpaceArrayDiffusion &sh) const;
  void outputConcentration2DEGBin(const char *dname, const char *fhead) const;
  void outputFermiLevel2DEG(const char *dname, const char *fhead) const;
  void outputPotential2DEG(const char *dname, const char *fhead) const;
  void outputPotential2DEG(const char *dname, const char *fhead,
			   const RealSpaceArrayDiffusion &pot) const;
  void outputField2DEG(const char *dname, const char *fhead) const;
  void outputField2DEG(const char *dname, const char *fhead,
		       const RealSpaceArrayDiffusion &field) const;
  void outputVelocity(const char *dname, const char *fhead);


private:

  const PoissonSolver2DDescriptor &_poiDsc;
  const DiffusionSolver1DDescriptor &_difDsc;
  PoissonSolver2D _poisson;
  DiffusionSolver1DFD _diffusion;
  RealSpaceGridHandler _realSGH;
  DiffusionABCalculator _ab;

  Concentration _SigmaElectron, _SigmaHole, _SigmaDope;
  RealSpaceArrayDiffusion _Ex, _muElectron, _muHole;

  FermiDistrGraphene _fermiDistr;

  PoissonDiffusionMediator _pdm;

  /*
  GridParameter _gridParam;
  */
  int _nSteps;
  double _t;

  OutputDirectoryManager _odm;
  std::string _SSDir, _muDir, _concDir, _pot2DDir, _fieldDir, _velDir;

  void _initPoissonSolver();
  void _refineMesh();
  void _setInitialSteadyStateSCF();

};
