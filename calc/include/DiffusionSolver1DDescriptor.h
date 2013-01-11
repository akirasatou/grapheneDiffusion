#pragma once
#include <vector>
#include <string>
#include <sstream>
#include <SystemDescriptor/SystemDescriptor.h>

using namespace std;


/*
 * Class reading and holding the parameters.
 * Initializing this class calls gmsh, outputs a '.msh' file, 
 * and passes it to libmesh.
 */

class DiffusionSolver1DDescriptor
{

 public:

  DiffusionSolver1DDescriptor(const string &fileHead);

  inline const string getFileHeadStr() const { return _fileHeadStr; }

  inline int getNx() const { return _Nx; }
  inline double get_dt() const { return _dt; }

  inline double getLc() const { return _Lc; }
  inline double getSigma0() const { return _Sigma0; }
  inline double getT() const { return _T; }
  inline double get_tMax() const { return _tMax; }

  inline double get_tOutputStep() const { return _tOutputStep; }
  inline double get_tOutputBinStep() const { return _tOutputBinStep; }
  inline bool toOutputSS() const { return _toOutputSS; }
  inline bool toOutputConcentration() const { return _toOutputConcentration; }
  inline bool toOutputPotential2D() const { return _toOutputPotential2D; }
  inline bool toOutputField() const { return _toOutputField; }
  inline bool toOutputVelocity() const { return _toOutputVelocity; }

  ~DiffusionSolver1DDescriptor();


 private:

  static const string className;

  SystemDescriptor _sysDsc;
  string _fileHeadStr, _className;
  void _readDsc();

  int _Nx;
  double _dt;
  void _registerMeshSection();
  void _setMeshSection();

  double _Lc, _Sigma0, _T, _alpha, _tMax;
  void _registerSimulationSection();
  void _setSimulationSection();

  double _tOutputStep, _tOutputBinStep;
  bool _toOutputSS, _toOutputConcentration, _toOutputPotential2D;
  bool _toOutputField, _toOutputVelocity;
  void _registerOutputSection();
  void _setOutputSection();

};
