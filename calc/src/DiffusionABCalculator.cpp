#include "DiffusionABCalculator.h"
#include <math.h>
#include <PhysicalConstants.h>
#include <Graphene/Constants.h>

using namespace std;


/*
 * Constructor.
 */

DiffusionABCalculator::DiffusionABCalculator(double T, double alpha):
  _T(T), _alpha(alpha), _ookBT(1.0/(kB*T))
{
}


/*
 * Calculate the quantities A, B, $dA/d\mu$, and $dB/d\mu$.
 */

double DiffusionABCalculator::calcA(double mu) const
{
  double c = 0.25*_alpha*vF*vF*hbar*_ookBT*_ookBT;
  double x = mu*_ookBT;

  return c/(log(1+exp(x))*(1+cosh(x)));
}

double DiffusionABCalculator::calcB(double mu) const
{
  double c = 0.5*_alpha*vF*vF*hbar*_ookBT;
  double x = mu*_ookBT;

  return c/(log(1+exp(x))*(1+exp(-x)));
}

double DiffusionABCalculator::calc_dA_dmu(double mu) const
{
  double c = 0.25*_alpha*vF*vF*hbar*_ookBT*_ookBT*_ookBT;
  double x = mu*_ookBT;
  double expx = exp(x), expmx = 1/expx;
  double oolog = 1.0/log(1+expx);
  double oocosh = 1.0/(1+cosh(x));

  return -c*oolog*oocosh*(oolog/(1+expmx)+sinh(x)*oocosh);
}

double DiffusionABCalculator::calc_dB_dmu(double mu) const
{
  double c = 0.5*_alpha*vF*vF*hbar*_ookBT*_ookBT;
  double x = mu*_ookBT;
  double expx = exp(x), expmx = 1/expx;
  double oolog = 1.0/log(1+expx);

  return -c*oolog/((1+expmx)*(1+expmx))*(-oolog+expmx);
}
