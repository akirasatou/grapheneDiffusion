#include "Conductivity.h"
#include <math.h>
#include <PhysicalConstants.h>
#include <Graphene/Constants.h>
#include <Integrator/Integrator1D.h>

using namespace std;


Complex sigma_qe(double hw, double epsf, double Te, double tau_m){
  return sigma_intra_qe(hw, epsf, Te, tau_m)+sigma_inter_qe(hw, epsf, Te);
}

Complex sigma_qe_taup(double hw, double epsf, double Te, double v_tau){
  return sigma_intra_qe_taup(hw, epsf, Te, v_tau)+sigma_inter_qe(hw, epsf, Te);
}

double sigma_inter_qe(double hw, double epsf, double Te){
  double f = 1/(exp((hw/2-epsf)/(kB*Te))+1);

  return e*e/((4*pi*eps0)*4*hbar)*(1-2*f);
}

double _integrand(double x, double x0){
  double t = exp(x-x0);

  return x*t/((1+t)*(1+t));
}

Complex sigma_intra_qe(double hw, double epsf, double Te, double tau_m){
  Complex r;
  double x0 = epsf/(kB*Te);

  r = integrate(_integrand, 0.0, x0+20, x0);
  //r = log(1+exp(x0));
  //r *= 2*e*e*kB*Te/(pi*hbar*hbar)*tau_m/(1+tau_m*tau_m*hw*hw/hbar/hbar);
  r *= 2*e*e*kB*Te/((4*pi*eps0)*pi*hbar*hbar)/(1-I*tau_m*hw/hbar);

  return r;
}

double _integrand2(double x, double x0, double a){
  double t = exp(x-x0);
  double x2 = x*x;

  return x2/(x2+a*a)*t/((1+t)*(1+t));
}

double _integrand3(double x, double x0, double a){
  double t = exp(x-x0);
  double x2 = x*x;

  return x*a/(x2+a*a)*t/((1+t)*(1+t));
}

Complex sigma_intra_qe_taup(double hw, double epsf, double Te, double v_tau){
  Complex r;
  double x0 = epsf/(kB*Te);
  double b = (hbar*vF)/(kB*Te*v_tau);
  //double a2 = pow(b*hw/hbar, 2.0);
  double a = b*hw/hbar;

  r = integrate(_integrand2, 0.0, x0+20, x0, a);
  r += I*integrate(_integrand3, 0.0, x0+20, x0, a);
  r *= b*2*e*e*kB*Te/((4*pi*eps0)*pi*hbar*hbar);

  return r;
}

/*
double sigma_intra_im(double hw, double epsf, double Te, double tau_m){
  double r;
  double x0 = epsf/(kB*Te);
  double wt = tau_m*hw/hbar;

  r = integrate(_integrand, 0.0, x0+10, x0);
  r *= 2*e*e*kB*Te/(pi*hbar*hbar)*tau_m*wt/(1+wt*wt);

  return r;
}
*/
