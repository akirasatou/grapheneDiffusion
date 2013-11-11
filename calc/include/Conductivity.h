#pragma once

#include <Complex/Complex.h>

Complex sigma_qe(double hw, double epsf, double Te, double tau_m);
Complex sigma_qe_taup(double hw, double epsf, double Te, double v_tau);
double sigma_inter_qe(double hw, double epsf, double Te);
Complex sigma_intra_qe(double hw, double epsf, double Te, double tau_m);
Complex sigma_intra_qe_taup(double hw, double epsf, double Te, double v_tau);

//double sigma_intra_im(double hw, double epsf, double Te, double tau_m);
