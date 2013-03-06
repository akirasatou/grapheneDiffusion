#pragma once
#include <vector>
#include <string>
#include <sstream>
#include <SystemDescriptor/SystemDescriptor.h>


/*
 * Class reading and holding the parameters.
 * Initializing this class calls gmsh, outputs a '.msh' file, 
 * and passes it to libmesh.
 */

class DiffusionSolver1DDescriptor
{

 public:

  DiffusionSolver1DDescriptor(const std::string &fileHead);

  inline const std::string getFileHeadStr() const { return _fileHeadStr; }

  inline int getNx() const { return _Nx; }
  inline double get_dtMax() const { return _dtMax; }
  inline double get_dmuMax() const { return _dmuMax; }

  inline double getLc() const { return _Lc; }
  inline double getSigma0() const { return _Sigma0; }
  inline double getT() const { return _T; }
  inline double get_alpha() const { return _alpha; }
  inline double get_tMax() const { return _tMax; }

  inline int get_nOutputStep() const { return _nOutputStep; }
  inline int get_nOutputBinStep() const { return _nOutputBinStep; }
  inline bool toOutputSS() const { return _toOutputSS; }
  inline bool toOutputFermiLevel() const { return _toOutputFermiLevel; }
  inline bool toOutputConcentration() const { return _toOutputConcentration; }
  inline bool toOutputPotential2D() const { return _toOutputPotential2D; }
  inline bool toOutputField() const { return _toOutputField; }
  inline bool toOutputVelocity() const { return _toOutputVelocity; }

  ~DiffusionSolver1DDescriptor();


 private:

  static const std::string className;

  SystemDescriptor _sysDsc;
  std::string _fileHeadStr, _className;
  void _readDsc();

  int _Nx;
  double _dtMax, _dmuMax;
  void _registerMeshSection();
  void _setMeshSection();

  double _Lc, _Sigma0, _T, _alpha, _tMax;
  void _registerSimulationSection();
  void _setSimulationSection();

  int _nOutputStep, _nOutputBinStep;
  bool _toOutputSS, _toOutputFermiLevel, _toOutputConcentration;
  bool _toOutputPotential2D, _toOutputField, _toOutputVelocity;
  void _registerOutputSection();
  void _setOutputSection();

};
