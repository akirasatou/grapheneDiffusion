#include <stdio.h>
#include <math.h>
#include <PhysicalConstants.h>
#include <PhysicalUnits.h>

using namespace std;

double sigma_inter(double E, double mu, double T)
{
  double f = 1.0/(1.0+exp((E-mu)/(kB*T)));

  return e*e/(4*pi*eps0*4*hbar)*(1-2*f);
}

double calc_gw(double f, double d, double epsilon, double mu, 
	       double T)
{
  double w = 2*pi*f;
  double wd_c = w*d/c;
  double kwd = sqrt(epsilon*wd_c*wd_c-pi*pi/4.0);

  return -2*pi/kwd*(wd_c)*(sigma_inter(hbar*w/2, mu, T)/c)/d;
}

int main()
{
  const double LgTab[4] = {nm2m(1250), nm2m(1000), nm2m(750), nm2m(500)};
  const double f = THz2Hz(37.5), d = nm2m(1000), epsilon = 4, T = 300;

  for(int i=0; i<4; i++){
    double Lg = LgTab[i];
    FILE *fin, *fout;
    char filename[1000];

    sprintf(filename, "../dat/VgLgDep/L=3000/mue-Lg=%g.dat",
	    m2nm(Lg));
    fin = fopen(filename, "r");

    sprintf(filename, "../dat/VgLgDep/L=3000/gw-Lg=%g.dat", m2nm(Lg));
    fout = fopen(filename, "w");

    while( true ){
      double Vg, mu;

      if( fscanf(fin, "%lf %lf", &Vg, &mu) == EOF ) break;

      mu = meV2J(mu);

      fprintf(fout, "%g %g\n", Vg, calc_gw(f, d, epsilon, mu, T)*1e-2);
    }

    fclose(fin);
    fclose(fout);
  }

  return 0;
}
