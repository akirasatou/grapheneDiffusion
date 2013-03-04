#include "Concentration.h"
#include <Integrator/Integrator1D.h>
#include <PhysicalConstants.h>

using namespace std;


/*
 * Constructor.
 */

Concentration::Concentration(const RealSpaceGridHandler &realSGH):
  RealSpaceArrayDiffusion(realSGH)
{
  _xTab = new double[_realSGH.getSize()+1];
  _integrandX = new double[_realSGH.getSize()+1];
}


/*
 * Destructor.
 */

Concentration::~Concentration()
{
  delete [] _xTab;
  delete [] _integrandX;
}


/*
 * Calculate the total concentration (in m$^{-1}$).
 */

double Concentration::calcTotalConcentration() const
{
  int N = _realSGH.getSize();


  for(int i=0; i<N; i++){
    _xTab[i] = _realSGH.getAt(i);
    _integrandX[i] = getAt(i);
  }
  _xTab[N] = _xTab[N-1]+(_xTab[N-1]-_xTab[N-2]);
  _integrandX[N] = getAt(0);

  return integrate(_integrandX, _xTab, N+1);
}
