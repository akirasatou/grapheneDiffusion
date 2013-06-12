#include <stdio.h>
#include <math.h>
#include <PhysicalConstants.h>
#include <PhysicalUnits.h>
#include <Complex/Complex.h>

using namespace std;

double sigma_inter(double E, double mu, double T)
{
  double f = 1.0/(1.0+exp((E-mu)/(kB*T)));

  return e*e/(4*pi*eps0*4*hbar)*(1-2*f);
}

Complex calc_kw(double f, double d, double epsilon, double mu, 
		double T)
{
  double w = 2*pi*f;
  double wd_c = w*d/c;
  double delta = 4.0*(sigma_inter(hbar*w/2, mu, T)/c)*wd_c;
  Complex a = 0.5*pi;
  Complex b = -I*2*pi*delta;

  Complex tmp = pi*sqrt(a*a+b);
  Complex kw_d = sqrt(epsilon*wd_c*wd_c-0.25*(0.5*pi*pi+b+tmp));
  Complex kappa_d = sqrt(epsilon*wd_c*wd_c-kw_d*kw_d);

  return kw_d/d;
}

int main()
{
  const double LgTab[4] = {nm2m(1250), nm2m(1000), nm2m(750), nm2m(500)};
  //const double f = THz2Hz(37.5), d = nm2m(1000), epsilon = 4, T = 300;
  const double d = nm2m(1000), epsilon = 4, T = 300;

  for(int i=0; i<4; i++){
    double Lg = LgTab[i];

    for(int j=70; j<=90; j+=5){
      FILE *fin, *fre, *fim;
      char filename[1000];
      double lambda = micro2m(j/10.0);
      double f = c/lambda;
      
      sprintf(filename, "../dat/VgLgDep/L=3000/mue-Lg=%g.dat",
	      m2nm(Lg));
      fin = fopen(filename, "r");

      sprintf(filename, "../dat/VgLgDep/L=3000/kw-Lg=%g-lambda=%g.dat",
	      m2nm(Lg), m2micro(lambda));
      fre = fopen(filename, "w");

      sprintf(filename, "../dat/VgLgDep/L=3000/gw-Lg=%g-lambda=%g.dat",
	      m2nm(Lg), m2micro(lambda));
      fim = fopen(filename, "w");

      while( true ){
	double Vg, mu;
	
	if( fscanf(fin, "%lf %lf", &Vg, &mu) == EOF ) break;
	
	mu = meV2J(mu);
	
	Complex kw = calc_kw(f, d, epsilon, mu, T);

	fprintf(fre, "%g %g\n", Vg, Re(kw)*1e-2);
	fprintf(fim, "%g %g\n", Vg, -Im(kw)*1e-2);
      }

      fclose(fin);
      fclose(fre);
      fclose(fim);
    }
  }

  return 0;
}
