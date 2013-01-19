#include "ChargeDensity2D.h"
#include <PhysicalConstants.h>

using namespace std;


ChargeDensity2D::
ChargeDensity2D(RealSpaceArrayDiffusion &SigmaElectron, 
		RealSpaceArrayDiffusion &SigmaHole,
		RealSpaceArrayDiffusion &Sigma0):
  _SigmaElectron(SigmaElectron), _SigmaHole(SigmaHole),
  _Sigma0(Sigma0)
{
}

double ChargeDensity2D::interpolate(double x) const
{
  double rho = 0.0;

  rho += _SigmaElectron.interpolate(x);
  rho -= _SigmaHole.interpolate(x);
  rho -= _Sigma0.interpolate(x);

  return -e*rho;
}

void ChargeDensity2D::updateInterpolator()
{
  _SigmaElectron.updateInterpolator();
  _SigmaHole.updateInterpolator();
  _Sigma0.updateInterpolator();
}
