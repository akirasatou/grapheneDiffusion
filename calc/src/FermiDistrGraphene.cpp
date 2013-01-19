#include "FermiDistrGraphene.h"
#include <math.h>
#include <Integrator/Integrator1D.h>
#include <PhysicalConstants.h>
#include <PhysicalUnits.h>
#include <Graphene/Constants.h>
#include <DistrFunction.h>

using namespace std;

/*
 * Constructor.
 */

FermiDistrGraphene::FermiDistrGraphene(double T):
  _T(T), _kBT(kB*T), _ookBT(1.0/_kBT)
{
}


/*
 * Fermi distribution function.
 */

double FermiDistrGraphene::getFermiDistr(double eps, double Ef) const
{
  return 1.0/(1+exp((eps-Ef)*_ookBT));
}


/*
 * Calculate the concentration from the Fermi level or vice versa.
 * The interal of the Fermi function is evaluated (almost) exactly.
 */

double FermiDistrGraphene::calcEfExact(double Sigma) const
{
  double Efb, Efe, Efm;
  double Sigma_m;

  for(Efb=-meV2J(20), Efe=0.0; ; Efb=Efe, Efe+=meV2J(20)){
    if( calcConcentrationFermiExact(Efe) >= Sigma ) break;
  }

  int maxCnt = 10000;
  int cnt = 0;
  
  do {
    Efm = (Efb+Efe)/2.0;
    Sigma_m = calcConcentrationFermiExact(Efm);
    
    if(Sigma < Sigma_m) Efe = Efm;
    else Efb = Efm;

    if(cnt++ > maxCnt) break;
  }while(fabs(1.0-Sigma_m/Sigma) > 1e-9);

  return Efm;
}

double FermiDistrGraphene::
calcEfEHExact(double Sigma) const
{
  double Efb, Efe, Efm;
  double Sigma_m;

  Efb = -eV2J(2.0);
  Efe = eV2J(2.0);

  /*
  for(Efb=-meV2J(20), Efe=0.0; ; Efb=Efe, Efe+=meV2J(20)){
    if( calcConcentrationFermiWithGrid(Efe) >= Sigma ) break;
  }
  */

  int maxCnt = 10000;
  int cnt = 0;

  do {
    Efm = (Efb+Efe)/2.0;

    // Electron.
    Sigma_m = calcConcentrationFermiExact(Efm);

    // Hole.
    Sigma_m -= calcConcentrationFermiExact(-Efm);

    if(Sigma < Sigma_m) Efe = Efm;
    else Efb = Efm;

    if(cnt++ > maxCnt) break;
  }while(fabs(1.0-Sigma_m/Sigma) > 1e-9);

  return Efm;
}

double FermiDistrGraphene::
calcConcentrationFermiExact(double Ef) const
{
  double s = FermiIntegral(1, Ef, _T);

  return 2.0/M_PI*pow(_kBT/(vF*hbar), 2)*s;
}
