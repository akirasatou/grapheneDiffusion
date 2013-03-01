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

  RealSpaceGridHandler(const std::vector<double> &x);
  RealSpaceGridHandler(const std::vector<std::vector<double> > &x);

  void setAt(int i, double val);
  void addAt(int i, double val);
  double getAt(int i) const;
  int getSize() const;

  double getXl() const;
  double getXr() const;
  std::pair<int, int> getPointID(int i) const;
  std::vector<double> getPointsInElem(int iElem) const;
  int getPointInvID(int iElem, int qp) const;


private:

  std::vector<double> _x; // x-grid in ascending order.
  std::vector<std::pair<int, int> > _ID; // ID (iElem, qp) for x.
  std::vector<std::vector<int> > _invID; // Index of x for ID (iElem, qp).

};
