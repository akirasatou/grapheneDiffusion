#include <stdio.h>
#include <math.h>
#include <PhysicalConstants.h>
#include <PhysicalUnits.h>
#include <Graphene/Constants.h>
#include <Complex/Complex.h>
#include "Conductivity.h"

using namespace std;


Complex calc_kw(double f, double d, Complex epsilon, double mu, 
		double T, double v_tau, int nGL)
{
  double w = 2*pi*f;
  double wd_c = w*d/c;
  Complex delta = 4.0*(nGL*sigma_qe_taup(hbar*w, mu, T, v_tau)/c)*wd_c;
  //Complex delta = 4.0*(nGL*sigma_inter_qe(hbar*w, mu, T)/c)*wd_c;
  Complex a = 0.5*pi;
  Complex b = -I*2*pi*delta;

  Complex tmp = pi*sqrt(a*a+b);
  Complex kw_d = sqrt(epsilon*wd_c*wd_c-0.25*(0.5*pi*pi+b+tmp));
  Complex kappa_d = sqrt(epsilon*wd_c*wd_c-kw_d*kw_d);

  return kw_d/d;
}

int main()
{
  const int nKappa = 2, nLambda = 3, nT = 3, nnGL = 2;

  const double lambdaTab[nLambda] = {micro2m(12), micro2m(15), micro2m(30)};
  const double WgBeginTab[nLambda] = {nm2m(1500), nm2m(1880), nm2m(3750)};
  const double WgEndTab[nLambda] = {nm2m(1550), nm2m(1940), nm2m(3850)};
  const double muTab[nLambda] = {meV2J(70), meV2J(60), meV2J(40)};
  const double TTab[nT] = {20, 77, 300};
  const double kappaTab[nnGL][nKappa] = {{2e-3, 4e-3}, {9e-3, 2.5e-2}};
  const int nGLTab[nnGL] = {1, 4};
  const double n = 2.0;
  const double v_tau = 0.1*vF;


  // Conductivity
  
  for(int iLambda=0; iLambda<nLambda; iLambda++){
    double mu = muTab[iLambda];

    for(int iT=0; iT<nT; iT++){
      double T = TTab[iT];


      FILE *fre, *fim;
      char filename[1000];
	  
      sprintf(filename, "../dat/gain/sigma-re-mu=%g-T=%g.dat",
	      J2meV(mu), T);
      fre = fopen(filename, "w");

      sprintf(filename, "../dat/gain/sigma-im-mu=%g-T=%g.dat",
	      J2meV(mu), T);
      fim = fopen(filename, "w");

      for(double f=THz2Hz(10); f<=THz2Hz(25); f+=THz2Hz(0.5)){
	double w = 2*pi*f;
	Complex sigma = sigma_intra_qe_taup(hbar*w, mu, T, v_tau);
	double sigma0 = e*e/((4*pi*eps0)*4*hbar);

	fprintf(fre, "%g %g\n", Hz2THz(f), Re(sigma)/sigma0);
	fprintf(fim, "%g %g\n", Hz2THz(f), Im(sigma)/sigma0);
      }

      fclose(fre);
      fclose(fim);
    }
  }


  // Wg dep.

  for(int iLambda=0; iLambda<nLambda; iLambda++){
    double lambda = lambdaTab[iLambda];
    double mu = muTab[iLambda];
    double Wgb = WgBeginTab[iLambda], Wge = WgEndTab[iLambda];
    double WgStep = (Wge-Wgb)/300;

    double f = c/lambda;

    for(int iT=0; iT<nT; iT++){
      double T = TTab[iT];

      for(int inGL=0; inGL<nnGL; inGL++){
	int nGL = nGLTab[inGL];
	
	for(int iKappa=0; iKappa<nKappa; iKappa++){
	  double kappa = kappaTab[inGL][iKappa];
	  Complex epsilon = Complex(n*n-kappa*kappa, 2*n*kappa);

	  FILE *fre, *fim;
	  char filename[1000];
	  
	  sprintf(filename, "../dat/gain/kw-WgDep-lambda=%g-nGL=%d-kappa=%g-T=%g.dat",
		  m2micro(lambda), nGL, kappa, T);
	  fre = fopen(filename, "w");

	  sprintf(filename, "../dat/gain/gw-WgDep-lambda=%g-nGL=%d-kappa=%g-T=%g.dat",
		  m2micro(lambda), nGL, kappa, T);
	  fim = fopen(filename, "w");
	  
	  for(double Wg=Wgb; Wg<=Wge+nm2m(0.1); Wg+=WgStep){
	    Complex kw = calc_kw(f, Wg, epsilon, mu, T, v_tau, nGL);
	    
	    fprintf(fre, "%g %g\n", m2micro(Wg), Re(kw)*1e-2);
	    fprintf(fim, "%g %g\n", m2micro(Wg), -Im(kw)*1e-2);
	  }
	 
	  fclose(fre);
	  fclose(fim);
	}
      }
    }
  }


  // mu dep.

  const double WgTab[nLambda] = {nm2m(1510), nm2m(1890), nm2m(3780)};

  for(int iLambda=0; iLambda<nLambda; iLambda++){
    double lambda = lambdaTab[iLambda];
    double Wg = WgTab[iLambda];
    double f = c/lambda;

    for(int iT=0; iT<nT; iT++){
      double T = TTab[iT];

      for(int inGL=0; inGL<nnGL; inGL++){
	int nGL = nGLTab[inGL];
	
	for(int iKappa=0; iKappa<nKappa; iKappa++){
	  double kappa = kappaTab[inGL][iKappa];
	  Complex epsilon = Complex(n*n-kappa*kappa, 2*n*kappa);

	  FILE *fre, *fim;
	  char filename[1000];
	  
	  sprintf(filename, "../dat/gain/kw-muDep-lambda=%g-nGL=%d-kappa=%g-T=%g.dat",
		  m2micro(lambda), nGL, kappa, T);
	  fre = fopen(filename, "w");

	  sprintf(filename, "../dat/gain/gw-muDep-lambda=%g-nGL=%d-kappa=%g-T=%g.dat",
		  m2micro(lambda), nGL, kappa, T);
	  fim = fopen(filename, "w");
	  
	  for(double mu=meV2J(20); mu<=meV2J(80); mu+=meV2J(0.1)){
	    Complex kw = calc_kw(f, Wg, epsilon, mu, T, v_tau, nGL);
	    
	    fprintf(fre, "%g %g\n", J2meV(mu), Re(kw)*1e-2);
	    fprintf(fim, "%g %g\n", J2meV(mu), -Im(kw)*1e-2);
	  }
	  
	  fclose(fre);
	  fclose(fim);
	}
      }
    }
  }

  return 0;
}
