#include <iostream>
#include <PhysicalConstants.h>
#include <PhysicalUnits.h>

using namespace std;

double Sigma(double Vg, double Wg)
{
  return 4*eps0*Vg/(e*Wg);
}

int main()
{
  double Wg = micro2m(0.3);

  cout << 1 << " " << Sigma(10, Wg)/1e16 << endl;

  return 0;
}
