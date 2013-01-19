#pragma once


/*
 * Class for Fermi distribution in graphene.
 */

class FermiDistrGraphene
{

public:

  FermiDistrGraphene(double T);

  double getFermiDistr(double eps, double Ef) const;

  double calcConcentrationFermiExact(double Ef) const;
  double calcEfEHExact(double Sigma) const;
  double calcEfExact(double Sigma) const;


private:

  double _T, _kBT, _ookBT;

};
