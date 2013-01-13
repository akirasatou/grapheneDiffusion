#include "RealSpaceGridHandler.h"
#include <stdlib.h>
#include <iostream>
#include <algorithm>

using namespace std;


/*
 * Constructor.
 */

RealSpaceGridHandler::
RealSpaceGridHandler(const std::vector<std::vector<double> > &x)
{
  vector<pair<double, pair<int, int> > > tmp;

  _invID.resize(x.size());

  for(int i=0; i<x.size(); i++){

    _invID[i].resize(x[i].size());

    for(int j=0; j<x[i].size(); j++){
      pair<double, pair<int, int> > p;

      p.first = x[i][j];
      p.second.first = i;
      p.second.second = j;

      tmp.push_back(p);
    }
  }

  sort(tmp.begin(), tmp.end());

  _x.resize(tmp.size());
  _ID.resize(tmp.size());

  for(int i=0; i<tmp.size(); i++){
    _x[i] = tmp[i].first;
    _ID[i] = tmp[i].second;
    _invID[_ID[i].first][_ID[i].second] = i;
  }
}


/*
 * Accessors.
 */

int RealSpaceGridHandler::getSize() const
{
  return _x.size();
}

double RealSpaceGridHandler::getX(int i) const
{
  if( i<0 || i>=_x.size() ){
    cerr << "RealSpaceGridHandler::getX: index out of range: ";
    cerr << i << endl;
    exit(1);
  }

  return _x[i];
}

pair<int, int> RealSpaceGridHandler::getPointID(int i) const
{
  if( i<0 || i>=_x.size() ){
    cerr << "RealSpaceGridHandler::getPointID: index out of range: ";
    cerr << i << endl;
    exit(1);
  }

  return _ID[i];
}
