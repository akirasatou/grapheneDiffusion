#include "DiffusionSolver1DFD.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <PhysicalUnits.h>

using namespace std;

#define isEqual(a, b) (fabs((a)/(b)-1.0) < 1e-5)


/*
 * Class name for error output.
 */

const string DiffusionSolver1DFD::_className("DiffusionSolver1DFD");


/*
 * Constructor.
 */

DiffusionSolver1DFD::
DiffusionSolver1DFD(const DiffusionSolver1DDescriptor &difDsc,
		    const DiffusionABCalculator &ab,
		    PoissonDiffusionMediator &pdm,
		    double xl, double xr):
  _sysName(difDsc.getFileHeadStr()), _difDsc(difDsc), _ab(ab),
  _pdm(pdm), _R(_difDsc.getNx()), _U(_difDsc.getNx()),
  _dU(_difDsc.getNx()), _J(_difDsc.getNx(), 0.0),
  _dx(_difDsc.getLc()/_difDsc.getNx()), _Xl(xl), _Xr(xr),
  _tNorm(fs2s(1)), _xNorm(micro2m(1)), _muNorm(meV2J(1)),
  _dmudxNorm(_muNorm/_xNorm), _d2mudx2Norm(_muNorm/(_xNorm*_xNorm)),
  _ANorm(_xNorm*_xNorm/(_muNorm*_tNorm)),
  _BNorm(_xNorm*_xNorm/_tNorm), _eExNorm(_muNorm/_xNorm)
{
}


/*
 * Destructor.
 */

DiffusionSolver1DFD::~DiffusionSolver1DFD()
{

}


/*
 * Return RealSpaceGridHandler object.
 */

RealSpaceGridHandler DiffusionSolver1DFD::getRealSGH() const
{
  vector<double> x(_difDsc.getNx());

  for(int i=0; i<x.size(); i++){
    x[i] = i*_difDsc.getLc()/_difDsc.getNx()+_Xl;
  }

  return RealSpaceGridHandler(x);
}
