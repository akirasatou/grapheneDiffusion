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
  const double epsilon = 4;
  const double lambdaTab[3] = {micro2m(12), micro2m(15), micro2m(30)};
  const double WgBeginTab[3] = {nm2m(1400), nm2m(1800), nm2m(3500)};
  const double WgEndTab[3] = {nm2m(1600), nm2m(2000), nm2m(4000)};
  const double TTab[2] = {77, 300};

  // Wg dep.

  for(int i=0; i<3; i++){
    double lambda = lambdaTab[i];
    double f = c/lambda;

    for(int j=0; j<2; j++){
      double T = TTab[j];

      for(int k=0; k<3; k++){
	double mu = meV2J(40)+meV2J(10)*k;

	FILE *fre, *fim;
	char filename[1000];

	sprintf(filename, "../dat/VgWgLgDep-grounded/kw-lambda=%g-mu=%g-T=%g.dat",
		m2micro(lambda), J2meV(mu), T);
	fre = fopen(filename, "w");
      
	sprintf(filename, "../dat/VgWgLgDep-grounded/gw-lambda=%g-mu=%g-T=%g.dat",
		m2micro(lambda), J2meV(mu), T);
	fim = fopen(filename, "w");
      
	double Wgb = WgBeginTab[i], Wge = WgEndTab[i];
	double WgStep = (Wge-Wgb)/300;
      

	for(double Wg=Wgb; Wg<=Wge+nm2m(0.1); Wg+=WgStep){

	
	  Complex kw = calc_kw(f, Wg, epsilon, mu, T);
	
	  fprintf(fre, "%g %g\n", m2micro(Wg), Re(kw)*1e-2);
	  fprintf(fim, "%g %g\n", m2micro(Wg), -Im(kw)*1e-2);
	}
      
	fclose(fre);
	fclose(fim);
      }
    }
  }


  // mu dep.

  const double WgTab[3] = {nm2m(1510), nm2m(1890), nm2m(3780)};

  for(int i=0; i<3; i++){
    double lambda = lambdaTab[i];
    double f = c/lambda;
    double Wg = WgTab[i];

    for(int j=0; j<2; j++){
      double T = TTab[j];

      FILE *fre, *fim;
      char filename[1000];

      sprintf(filename, "../dat/VgWgLgDep-grounded/kw-lambda=%g-Wg=%g-T=%g.dat",
	      m2micro(lambda), m2nm(Wg), T);
      fre = fopen(filename, "w");
      
      sprintf(filename, "../dat/VgWgLgDep-grounded/gw-lambda=%g-Wg=%g-T=%g.dat",
	      m2micro(lambda), m2nm(Wg), T);
      fim = fopen(filename, "w");

      for(double mu=meV2J(20); mu<=meV2J(80); mu+=meV2J(0.1)){
	Complex kw = calc_kw(f, Wg, epsilon, mu, T);
	
	fprintf(fre, "%g %g\n", J2meV(mu), Re(kw)*1e-2);
	fprintf(fim, "%g %g\n", J2meV(mu), -Im(kw)*1e-2);
      }
      
      fclose(fre);
      fclose(fim);
    }
  }

  return 0;
}
