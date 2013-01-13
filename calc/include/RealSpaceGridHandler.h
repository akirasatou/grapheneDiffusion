#pragma once
#include <vector>
#include <utility>


/*
 * Class handling the real-space points for DiffusionSolver1D.
 */

class RealSpaceGridHandler
{

public:

  RealSpaceGridHandler(const std::vector<std::vector<double> > &x);

  int getSize() const;
  double getX(int i) const;
  std::pair<int, int> getPointID(int i) const;
  int getPointInvID(int iElem, int qp) const;


private:

  std::vector<double> _x;
  std::vector<std::pair<int, int> > _ID;
  std::vector<std::vector<int> > _invID;

};
