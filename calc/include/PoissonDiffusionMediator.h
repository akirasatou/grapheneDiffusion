#pragma once

#include "RealSpaceGridHandler.h"
#include "RealSpaceArrayDiffusion.h"
#include "FermiDistrGraphene.h"
#include <Transistor2D/PoissonSolver2D.h>

class PoissonDiffusionMediator
{

public:

  PoissonDiffusionMediator(const RealSpaceGridHandler &realSGH,
			   PoissonSolver2D &poisson,
			   const RealSpaceArrayDiffusion &SigmaDope,
			   const FermiDistrGraphene &fermiDistr);
  void setInitialSolutions(const RealSpaceArrayDiffusion &mue0,
			   const RealSpaceArrayDiffusion &muh0,
			   const RealSpaceArrayDiffusion &Ex0,
			   const RealSpaceArrayDiffusion &dEx_dx0);
  void getCurrentSolutions(RealSpaceArrayDiffusion &mue,
			   RealSpaceArrayDiffusion &muh,
			   RealSpaceArrayDiffusion &Ex) const;
  void updateSolutions();

  void setTime(double t, double dt);
  void setInitialSolutionsNI();
  void updateSolutionsNI(const RealSpaceArrayDiffusion &mue,
			 const RealSpaceArrayDiffusion &muh);
  void updateSolutionsNI(const RealSpaceArrayDiffusion &mue,
			 const RealSpaceArrayDiffusion &dmue_dx,
			 const RealSpaceArrayDiffusion &d2mue_dx2,
			 const RealSpaceArrayDiffusion &muh,
			 const RealSpaceArrayDiffusion &dmuh_dx,
			 const RealSpaceArrayDiffusion &d2muh_dx2);

  inline double get_mue_n(int i) const { return _mue_n.getAt(i); }
  inline double get_mue_n1_l(int i) const { return _mue_n1_l.getAt(i); }
  inline double get_mue_n1_l1(int i) const { return _mue_n1_l1.getAt(i); }
  inline double get_dmue_dx_n(int i) const { return _dmue_dx_n.getAt(i); }
  inline double get_dmue_dx_n1_l(int i) const { return _dmue_dx_n1_l.getAt(i); }
  inline double get_dmue_dx_n1_l1(int i) const { return _dmue_dx_n1_l1.getAt(i); }
  inline double get_d2mue_dx2_n(int i) const { return _d2mue_dx2_n.getAt(i); }
  inline double get_d2mue_dx2_n1_l(int i) const { return _d2mue_dx2_n1_l.getAt(i); }
  inline double get_d2mue_dx2_n1_l1(int i) const { return _d2mue_dx2_n1_l1.getAt(i); }

  inline double get_muh_n(int i) const { return _muh_n.getAt(i); }
  inline double get_muh_n1_l(int i) const { return _muh_n1_l.getAt(i); }
  inline double get_muh_n1_l1(int i) const { return _muh_n1_l1.getAt(i); }
  inline double get_dmuh_dx_n(int i) const { return _dmuh_dx_n.getAt(i); }
  inline double get_dmuh_dx_n1_l(int i) const { return _dmuh_dx_n1_l.getAt(i); }
  inline double get_dmuh_dx_n1_l1(int i) const { return _dmuh_dx_n1_l1.getAt(i); }
  inline double get_d2muh_dx2_n(int i) const { return _d2muh_dx2_n.getAt(i); }
  inline double get_d2muh_dx2_n1_l(int i) const { return _d2muh_dx2_n1_l.getAt(i); }
  inline double get_d2muh_dx2_n1_l1(int i) const { return _d2muh_dx2_n1_l1.getAt(i); }
  inline double get_Ex_n(int i) const { return _Ex_n.getAt(i); }
  inline double get_Ex_n1_l(int i) const { return _Ex_n1_l.getAt(i); }
  inline double get_Ex_n1_l1(int i) const { return _Ex_n1_l1.getAt(i); }
  inline double get_dEx_dx_n(int i) const { return _dEx_dx_n.getAt(i); }
  inline double get_dEx_dx_n1_l(int i) const { return _dEx_dx_n1_l.getAt(i); }
  inline double get_dEx_dx_n1_l1(int i) const { return _dEx_dx_n1_l1.getAt(i); }


private:
    
  RealSpaceArrayDiffusion _mue_n, _mue_n1_l, _mue_n1_l1;
  RealSpaceArrayDiffusion _dmue_dx_n, _dmue_dx_n1_l, _dmue_dx_n1_l1;
  RealSpaceArrayDiffusion _d2mue_dx2_n, _d2mue_dx2_n1_l, _d2mue_dx2_n1_l1;
  RealSpaceArrayDiffusion _muh_n, _muh_n1_l, _muh_n1_l1;
  RealSpaceArrayDiffusion _dmuh_dx_n, _dmuh_dx_n1_l, _dmuh_dx_n1_l1;
  RealSpaceArrayDiffusion _d2muh_dx2_n, _d2muh_dx2_n1_l, _d2muh_dx2_n1_l1;
  RealSpaceArrayDiffusion _Ex_n, _Ex_n1_l, _Ex_n1_l1;
  RealSpaceArrayDiffusion _dEx_dx_n, _dEx_dx_n1_l, _dEx_dx_n1_l1;
  RealSpaceArrayDiffusion _SigmaDope;
  PoissonSolver2D &_poisson;
  double _t, _dt;
  int _nrStepsNI;
  const RealSpaceGridHandler &_realSGH;
  const FermiDistrGraphene _fermiDistr;

  static double calc_dmudx(const RealSpaceArrayDiffusion &mu, int i);
  static double calc_d2mudx2(const RealSpaceArrayDiffusion &mu, int i);

};
