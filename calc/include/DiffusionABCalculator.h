#pragma once


/*
 * Class that calculates the quantities A, B, $dA/d\mu$, and 
 * $dB/d\mu$.
 */

class DiffusionABCalculator
{

public:
  
  DiffusionABCalculator(double T, double alpha);

  double calcA(double mu) const;
  double calcB(double mu) const;
  double calc_dA_dmu(double mu) const;
  double calc_dB_dmu(double mu) const;


private:

  double _T, _alpha, _ookBT;

};
