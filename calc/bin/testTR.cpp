#include <stdio.h>
#include <math.h>
#include <PhysicalConstants.h>
#include <Graphene/Constants.h>
#include <PhysicalUnits.h>

using namespace std;

int main()
{
  double R = 2*sqrt(2)*pow(e, 3./2.)/(pi*pi*pow(hbar, 3./2.)*sqrt(vF));
  double E = 1e8;

  R *= pow(E, 3./2.);

  printf("%g\n", R*1e-4);

  return 0;
}

