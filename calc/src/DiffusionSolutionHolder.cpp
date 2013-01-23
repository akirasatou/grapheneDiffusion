#include "DiffusionSolutionHolder.h"

using namespace std;


/*
 * Constructor.
 */

DiffusionSolutionHolder::
DiffusionSolutionHolder(const RealSpaceGridHandler &realSGH):
  _mue_n(realSGH), _muh_n(realSGH), _Ex_n(realSGH),
  _mue_n1_l(realSGH), _muh_n1_l(realSGH), _Ex_n1_l(realSGH)
{
}
