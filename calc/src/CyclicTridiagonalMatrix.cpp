#include "CyclicTridiagonalMatrix.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

/*
 * Some tricks to use LAPACK routines.
 */


#ifdef USE_LIBSCALAPACK_UCASE
# define LAPACK_DGTSV DGTSV
#elif USE_LIBSCALAPACK_LCASE
# define LAPACK_DGTSV dgtsv
#elif USE_LIBSCALAPACK_LCASE_WITH_UNDERSCORE
# define LAPACK_DGTSV dgtsv_
#endif

extern "C" {
  void LAPACK_DGTSV(int *n, int *nrhs, double *dl, double *d,
		    double *du, double *b, int *ldb, int *info);
};


/*
 * Static members.
 */

const string CyclicTridiagonalMatrix::_className("CyclicTridiagonalMatrix");


/*
 * Constructor.
 */

CyclicTridiagonalMatrix::CyclicTridiagonalMatrix(int n, double v):
  _size(n)
{
  _dl = new double[n-1];
  _d = new double[n];
  _du = new double[n-1];

  for(int i=0; i<n-1; i++){
    _dl[i] = _d[i] = _du[i] = v;
  }
  _d[n-1] = _alpha = _beta = v;
}


/*
 * Destructor.
 */

CyclicTridiagonalMatrix::~CyclicTridiagonalMatrix()
{
  delete [] _dl;
  delete [] _d;
  delete [] _du;
}


/*
 * Accessors.
 */

int CyclicTridiagonalMatrix::size() const
{
  return _size;
}

void CyclicTridiagonalMatrix::setZeros()
{
  for(int i=0; i<_size-1; i++){
    _dl[i] = _d[i] = _du[i] = 0.0;
  }
  _d[_size-1] = _alpha = _beta = 0.0;
}

void CyclicTridiagonalMatrix::setAt(int i, int j, double v)
{
  if( i<0 || this->_size<=i || j<0 || this->_size<=j ){
    cerr << _className << "::setAt: index out of range" << endl;
    exit(1);
  }

  if( j == i ){ _d[i] = v; return; }
  if( j == i+1 ){ _du[i] = v; return; }
  if( j == i-1 ){ _dl[j] = v; return; }

  if( i==0 && j==_size-1 ){ _beta = v; return; }
  if( i==_size-1 && j==0 ){ _alpha = v; return; }

  cerr << i << " " << j << endl;
  cerr << _className << "::setAt: index out of range" << endl;
  exit(1);
}

void CyclicTridiagonalMatrix::addAt(int i, int j, double v)
{
  if( i<0 || this->_size<=i || j<0 || this->_size<=j ){
    cerr << _className << "::addAt: index out of range" << endl;
    exit(1);
  }

  if( j == i ){ _d[i] += v; return; }
  if( j == i+1 ){ _du[i] += v; return; }
  if( j == i-1 ){ _dl[j] += v; return; }

  if( i==0 && j==_size-1 ){ _beta += v; return; }
  if( i==_size-1 && j==0 ){ _alpha += v; return; }

  cerr << _className << "::addAt: index out of range" << endl;
  exit(1);
}

double CyclicTridiagonalMatrix::getAt(int i, int j) const
{
  if( i<0 || this->_size<=i || j<0 || this->_size<=j ){
    cerr << _className << "::getAt: index out of range" << endl;
    exit(1);
  }

  if( j == i ) return _d[i];
  if( j == i+1 ) return _du[i];
  if( j == i-1 ) return _dl[j];

  if( i==0 && j==_size-1 ) return _beta;
  if( i==_size-1 && j==0 ) return _alpha;

  cerr << _className << "::getAt: index out of range" << endl;
  exit(1);
}


/*
 * Solve the system of linear equations [a]{x} = {b}
 * using Gaussian elimination method with pivoting followed by
 * the application of the Sherman-Morrison formula 
 * (see "Numerical Recipe"). The solution is stored in x.
 * !Caution! this method modifies the content of *this.
 */

void CyclicTridiagonalMatrix::solveGauss(Vector &x, const Vector &b)
{
  double *bu = new double[2*_size];
  int nrhs = 2, ldb = _size, info;
  double gamma = -_d[0];
  double v1, vn;


  // Setup RHS vectors $b$ and $u$ and auxiliary vector $v$.

  for(int i=0; i<_size; i++){
    bu[i+0] = b[i];
    bu[i+_size] = 0.0;
  }

  bu[0+_size] = gamma;
  bu[_size-1+_size] = _alpha;

  v1 = 1.0;
  vn = _beta/gamma;


  // Replace the first and last diagonal elements.

  _d[0] -= gamma;
  _d[_size-1] -= _alpha*_beta/gamma;


  // Calculate $y = A^{-1}b$ and $z = A^{-1}u$.

  LAPACK_DGTSV(&_size, &nrhs, _dl, _d, _du, bu, &ldb, &info);


  // Calculate $x$ from $y$, $z$, and $v$.

  double vy = v1*bu[0+0]+vn*bu[_size-1+0];
  double vz = v1*bu[0+_size]+vn*bu[_size-1+_size];
  double a = vy/(1.0+vz);

  for(int i=0; i<_size; i++){
    x[i] = bu[i+0]-a*bu[i+_size];
  }
}
