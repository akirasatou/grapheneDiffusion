#pragma once

#include <Vector/Vector.h>
#include <string>

class CyclicTridiagonalMatrix
{

public:

  CyclicTridiagonalMatrix(int n, double v=0.0);
  ~CyclicTridiagonalMatrix();

  int size() const;
  void setZeros();
  void setAt(int i, int j, double v);
  void addAt(int i, int j, double v);
  double getAt(int i, int j) const;
  void solveGauss(Vector &x, const Vector &b);


private:

  static const std::string _className;  
  int _size;
  double *_dl, *_d, *_du, _alpha, _beta;

  CyclicTridiagonalMatrix();
  CyclicTridiagonalMatrix(const CyclicTridiagonalMatrix &);
  const CyclicTridiagonalMatrix &operator=(const CyclicTridiagonalMatrix &);

};
