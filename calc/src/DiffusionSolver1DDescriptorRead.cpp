#include "DiffusionSolver1DDescriptor.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <PhysicalConstants.h>
#include <PhysicalUnits.h>

using namespace std;


/*
 * Read a descriptor file for DiffusionSolver.
 *
 * Mesh parameters (compulsory):
 *
 * Nx --- The (approximate) number of points in the $x$-mesh: Nx>=50.
 * dtMax --- Maximum time step bounding the adaptive one (in fs).
 * dmuMax --- Maximum allowed change in $\mu$ by one time step 
 *            (in meV).
 *
 *
 * Simulation parameters (compulsory):
 *
 * Lc --- Length of the 2DEG (in nm). Can be shorter than or
 *        equal to the length of the Poisson box.
 * Sigma0 --- "Donor" concentration in the real 2DEG (in 10$^{12} 
 *            cm$^{-2}$). + for donor and - for acceptor.
 * T --- lattice temperature (in K).
 * alpha --- The dimensionless parameter for the momentum relaxation
 *           time ($\tau = \hbar\alpha/\vF p$):
 *           $\alpha \sim 10-10^{3}$ for $\tau = 0.1-10$ ps. 
 * tMax --- Maximum time to stop the time step (in fs), i.e.,
 *          the simulation time is 0<=t<=tMax.
 *
 *
 * Output parameters (compulsory):
 *
 * nOutputStep --- Number of time steps by which the ascii output are 
 *                 generated.
 * nOutputBinStep --- Number of time steps by which the binary output
 *                    are generated:
 *                    Default value: nOutputStep.
 * toOutput* --- whether to output concentration, field, and so on.
 *               Default value: false.
 * 
 */

void DiffusionSolver1DDescriptor::_readDsc()
{
  // Register sections.

  _registerMeshSection();
  _registerSimulationSection();
  _registerOutputSection();


  // Read the descriptor.

  _sysDsc.readSections();


  // Check the inter-parameter and inter-section validity
  // and set parameters of DiffusionSolver.

  _setMeshSection();
  _setSimulationSection();
  _setOutputSection();
}


/*
 * Register sections.
 */

void DiffusionSolver1DDescriptor::_registerMeshSection()
{
  Section &secMesh = _sysDsc.registerSection("Mesh");

  secMesh.registerParameterIntNoDefault("Nx");
  secMesh.getParameter("Nx").setPVC(new PVCIntGEQ(50));

  secMesh.registerParameterDoubleNoDefault("dtMax");
  secMesh.getParameter("dtMax").setPVC(new PVCDoubleGEQ(0.0));
  secMesh.getParameter("dtMax").setPPP(new PPPDoubleMultiply(fs2s(1)));

  secMesh.registerParameterDoubleNoDefault("dmuMax");
  secMesh.getParameter("dmuMax").setPVC(new PVCDoubleGEQ(0.0));
  secMesh.getParameter("dmuMax").setPPP(new PPPDoubleMultiply(meV2J(1)));
}

void DiffusionSolver1DDescriptor::_registerSimulationSection()
{
  Section &secSim = _sysDsc.registerSection("Simulation");

  secSim.registerParameterDoubleNoDefault("Sigma0");
  secSim.getParameter("Sigma0").setPPP(new PPPDoubleMultiply(1e16));

  secSim.registerParameterDoubleNoDefault("Lc");
  secSim.getParameter("Lc").setPVC(new PVCDoubleGEQ(0.0));
  secSim.getParameter("Lc").setPPP(new PPPDoubleMultiply(nm2m(1)));

  secSim.registerParameterDoubleNoDefault("T");
  secSim.getParameter("T").setPVC(new PVCDoubleGEQ(0.0));

  secSim.registerParameterDoubleNoDefault("alpha");
  secSim.getParameter("alpha").setPVC(new PVCDoubleGEQ(0.0));

  secSim.registerParameterDoubleNoDefault("tMax");
  secSim.getParameter("tMax").setPVC(new PVCDoubleGEQ(0.0));
  secSim.getParameter("tMax").setPPP(new PPPDoubleMultiply(fs2s(1)));
}

void DiffusionSolver1DDescriptor::_registerOutputSection()
{
  Section &secOutput = _sysDsc.registerSection("Output");

  secOutput.registerParameterIntNoDefault("nOutputStep");
  secOutput.getParameter("nOutputStep").setPVC(new PVCIntGEQ(1));

  secOutput.registerParameterIntWithDefault("nOutputBinStep", 0);
  secOutput.getParameter("nOutputBinStep").setPVC(new PVCIntGEQ(1));

  const char *boolTab[2] = {"true", "false"};

  secOutput.registerParameterStringWithDefault("toOutputSS", "false");
  secOutput.getParameter("toOutputSS").setPVC(new PVCStringInTab(boolTab, 2));
  secOutput.registerParameterStringWithDefault("toOutputFermiLevel", "false");
  secOutput.getParameter("toOutputFermiLevel").setPVC(new PVCStringInTab(boolTab, 2));
  secOutput.registerParameterStringWithDefault("toOutputConcentration", "false");
  secOutput.getParameter("toOutputConcentration").setPVC(new PVCStringInTab(boolTab, 2));
  secOutput.registerParameterStringWithDefault("toOutputPotential2D", "false");
  secOutput.getParameter("toOutputPotential2D").setPVC(new PVCStringInTab(boolTab, 2));
  secOutput.registerParameterStringWithDefault("toOutputField", "false");
  secOutput.getParameter("toOutputField").setPVC(new PVCStringInTab(boolTab, 2));
  secOutput.registerParameterStringWithDefault("toOutputVelocity", "false");
  secOutput.getParameter("toOutputVelocity").setPVC(new PVCStringInTab(boolTab, 2));
}


/*
 * Check the inter-parameter and inter-section validity
 * and set parameters of DiffusionSolver.
 */

void DiffusionSolver1DDescriptor::_setMeshSection()
{
  if( !_sysDsc.isRead("Mesh") ){
    cerr << className << ": Unread section: 'Mesh'." << endl;
    exit(1);
  }

  Section &secMesh = _sysDsc.getSection("Mesh");
  
  _Nx = secMesh.getParameter("Nx").getInt();
  _dtMax = secMesh.getParameter("dtMax").getDouble();
}

void DiffusionSolver1DDescriptor::_setSimulationSection()
{
  if( !_sysDsc.isRead("Simulation") ){
    cerr << className << ": Unread section: 'Simulation'." << endl;
    exit(1);
  }

  Section &secSim = _sysDsc.getSection("Simulation");
  
  _Sigma0 = secSim.getParameter("Sigma0").getDouble();
  _Lc = secSim.getParameter("Lc").getDouble();
  _T = secSim.getParameter("T").getDouble();
  _alpha = secSim.getParameter("alpha").getDouble();
  _tMax = secSim.getParameter("tMax").getDouble();
}

void DiffusionSolver1DDescriptor::_setOutputSection()
{
  if( !_sysDsc.isRead("Output") ){
    cerr << className << ": Unread section: 'Output'." << endl;
    exit(1);
  }

  Section &secOutput = _sysDsc.getSection("Output");
  
  _nOutputStep = secOutput.getParameter("nOutputStep").getInt();
  _nOutputBinStep = secOutput.getParameter("nOutputBinStep").getInt();

  if( _nOutputBinStep == 0 ){
    _nOutputBinStep = _nOutputStep;
  }

  _toOutputSS = (secOutput.getParameter("toOutputSS").getString()=="true");
  _toOutputFermiLevel = (secOutput.getParameter("toOutputFermiLevel").getString()=="true");
  _toOutputConcentration = (secOutput.getParameter("toOutputConcentration").getString()=="true");
  _toOutputPotential2D = (secOutput.getParameter("toOutputPotential2D").getString()=="true");
  _toOutputField = (secOutput.getParameter("toOutputField").getString()=="true");
  _toOutputVelocity = (secOutput.getParameter("toOutputVelocity").getString()=="true");
}
