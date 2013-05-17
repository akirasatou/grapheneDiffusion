#include <stdio.h>
#include <math.h>
#include <PhysicalConstants.h>
#include <PhysicalUnits.h>
#include <Graphene/Constants.h>

using namespace std;

double mueAppr(double L, double Lg, double d, double Vg, double eps)
{
  return vF*hbar*sqrt(pi*eps*eps0*Vg*Lg/(e*d*L));
}

int main(int argc, char **argv)
{
  FILE *fp;
  char filename[100];

  double LgTab[4] = {0.5, 0.75, 1, 1.25};
  double L = micro2m(3), d = micro2m(1);
  double eps = 4;

  for(int i=0; i<4; i++){
    double Lg = micro2m(LgTab[i]);

    sprintf(filename, "../dat/VgLgDep/L=3000/mue-Lg=%g-appr.dat",
	    m2nm(Lg));
    fp = fopen(filename, "w");

    for(double Vg=0.0; Vg<=50.01; Vg+=0.1){
      fprintf(fp, "%g %g\n", Vg, J2meV(mueAppr(L, Lg, d, Vg, eps)));
    }

    fclose(fp);
  }

  return 0;
}
