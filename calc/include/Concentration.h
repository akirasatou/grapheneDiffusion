#pragma once
#include "RealSpaceArrayDiffusion.h"


using namespace std;


/*
 * Class for carrier concentration.
 */

class Concentration: public RealSpaceArrayDiffusion
{
public:

  Concentration(const RealSpaceGridHandler &realSGH);
  ~Concentration();

  double calcTotalConcentration() const;


private:

  double *_xTab, *_integrandX;

  Concentration(const Concentration &s);
  Concentration & operator=(const Concentration &s);
};
