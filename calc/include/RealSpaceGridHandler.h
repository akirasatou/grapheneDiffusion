#pragma once
#include <vector>
#include <utility>
#include <RealSpaceArray/RealSpaceArray.h>


/*
 * Class handling the real-space points for DiffusionSolver1D.
 */

class RealSpaceGridHandler: public RealSpaceArray<double>
{

public:

  RealSpaceGridHandler(const std::vector<std::vector<double> > &x);

  void setAt(int i, double val);
  void addAt(int i, double val);
  double getAt(int i) const;
  int getSize() const;

  std::pair<int, int> getPointID(int i) const;
  int getPointInvID(int iElem, int qp) const;


private:

  std::vector<double> _x;
  std::vector<std::pair<int, int> > _ID;
  std::vector<std::vector<int> > _invID;

};
