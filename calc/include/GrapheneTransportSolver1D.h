#pragma once

#include "DiffusionSolver1D.h"
#include "RealSpaceGridHandler.h"
#include "RealSpaceArrayDiffusion.h"
#include "FermiDistrGraphene.h"
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

  void outputConcentration2DEG(const char *dname, const char *fhead,
			       const RealSpaceArrayDiffusion &se,
			       const RealSpaceArrayDiffusion &sh) const;
  void outputConcentration2DEGBin(const char *dname, const char *fhead) const;

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
  DiffusionSolver1D _diffusion;
  RealSpaceGridHandler _realSGH;

  RealSpaceArrayDiffusion _SigmaElectron, _SigmaHole, _SigmaDope;
  RealSpaceArrayDiffusion _Ex, _muElectron, _muHole;

  FermiDistrGraphene _fermiDistr;

  /*
  GridParameter _gridParam;
  */
  int _nSteps;

  OutputDirectoryManager _odm;
  std::string _SSDir, _concDir, _pot2DDir, _fieldDir, _velDir;

  void _initPoissonSolver();
  void _refineMesh();
  void _setInitialSteadyStateSCF();

};
