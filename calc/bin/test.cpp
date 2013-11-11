#include "Conductivity.h"
#include <PhysicalConstants.h>
#include <PhysicalUnits.h>
#include <Graphene/Constants.h>
#include <OutputDirectoryManager/OutputDirectoryManager.h>
#include <stdio.h>

using namespace std;

int main(){
  double fBegin = THz2Hz(1);
  double fEnd = THz2Hz(5);
  //double epsfTab[3] = {meV2J(40), meV2J(60), meV2J(70)};
  double epsfTab[3] = {meV2J(10), meV2J(20), meV2J(30)};
  double tauTab[3] = {10e-12, 1e-12, 0.1e-12};
  double vtauTab[3] = {0.001*vF, 0.01*vF, 0.1*vF};
  FILE *fp;
  char filename[100];

  OutputDirectoryManager odm("../dat/conductivity");

  for(int iepsf=0; iepsf<3; iepsf++){
    double epsf = epsfTab[iepsf];

    sprintf(filename, "%s/inter-epsf=%g.dat",
	    odm.getDirectory().c_str(), J2meV(epsf));
    fp = fopen(filename, "w");

    for(double f=fBegin; f<fEnd; f+=(fEnd-fBegin)/100){
      double sigma = sigma_inter_qe(h*f, epsf, 300);

      sigma /= e*e/(4*hbar);

      fprintf(fp, "%g %g\n", Hz2THz(f), sigma);
    }

    fclose(fp);

    for(int itau=0; itau<3; itau++){
      double tau = tauTab[itau];

      sprintf(filename, "%s/intra-epsf=%g-tau=%g.dat",
	      odm.getDirectory().c_str(), J2meV(epsf), s2ps(tau));
      fp = fopen(filename, "w");

      for(double f=fBegin; f<fEnd; f+=(fEnd-fBegin)/100){
	double sigma = sigma_intra_qe(h*f, epsf, 300, tau);

	sigma /= e*e/(4*hbar);

	fprintf(fp, "%g %g\n", Hz2THz(f), sigma);
      }

      fclose(fp);
    }

    for(int ivtau=0; ivtau<3; ivtau++){
      double vtau = vtauTab[ivtau];

      sprintf(filename, "%s/intra-epsf=%g-vtau=%g.dat",
	      odm.getDirectory().c_str(), J2meV(epsf), vtau);
      fp = fopen(filename, "w");

      for(double f=fBegin; f<fEnd; f+=(fEnd-fBegin)/100){
	double sigma = sigma_intra_qe_taup(h*f, epsf, 300, vtau);

	sigma /= e*e/(4*hbar);

	fprintf(fp, "%g %g\n", Hz2THz(f), sigma);
      }

      fclose(fp);
    }

    for(int ivtau=0; ivtau<3; ivtau++){
      double vtau = vtauTab[ivtau];

      sprintf(filename, "%s/total-epsf=%g-vtau=%g.dat",
	      odm.getDirectory().c_str(), J2meV(epsf), vtau);
      fp = fopen(filename, "w");

      for(double f=fBegin; f<fEnd; f+=(fEnd-fBegin)/100){
	double sigma = sigma_qe_taup(h*f, epsf, 300, vtau);

	sigma /= e*e/(4*hbar);

	fprintf(fp, "%g %g\n", Hz2THz(f), sigma);
      }

      fclose(fp);
    }
  }
}
