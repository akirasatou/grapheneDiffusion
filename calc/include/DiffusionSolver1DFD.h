#pragma once
#include "DiffusionSolver1DDescriptor.h"
#include "RealSpaceGridHandler.h"
#include "DiffusionABCalculator.h"
#include "PoissonDiffusionMediator.h"
#include <Transistor2D/PoissonSolver2D.h>
#include <Vector/Vector.h>
//#include <Matrix/SquareMatrix.h>
#include "CyclicTridiagonalMatrix.h"
#include <string>
#include <vector>


/*
 * Class solving the Diffusion equation using libMesh.
 */

class DiffusionSolver1DFD
{

public:

  DiffusionSolver1DFD(const DiffusionSolver1DDescriptor &difDsc,
		      const DiffusionABCalculator &ab,
		      PoissonDiffusionMediator &pdm,
		      double xl, double xr);
  ~DiffusionSolver1DFD();
  RealSpaceGridHandler getRealSGH() const;
  void solveStep(double t);


private:

  static const std::string _className;
  std::string _sysName;
  const DiffusionSolver1DDescriptor &_difDsc;
  const DiffusionABCalculator &_ab;
  PoissonDiffusionMediator &_pdm;
  Vector _R, _U, _dU;
  CyclicTridiagonalMatrix _J;

  double _dx, _Xl, _Xr;

  double _tNorm, _xNorm, _muNorm, _dmudxNorm, _d2mudx2Norm;
  double _eExNorm, _ANorm, _BNorm;

  double _calcNextTimeStep(double t) const;
  double _calcFInDiffusionEq(int i, int sr) const;
  void _calcRJ(const Vector &U, double dt, int sr);
  double _calcF(int i, double dt, int sr) const;
  double _calcK(int i, int j, const Vector &U, double dt, 
		int sr) const;
  double _calc_dKdU_U(int i, int j, const Vector &U, double dt, 
		      int sr) const;
  static int _get_jd(int j, int size);
  
};
