#include "DiffusionSolver1D.h"
#include <mesh_generation.h>
#include <edge_edge3.h>
#include <nonlinear_implicit_system.h>
#include <nonlinear_solver.h>
#include <transient_system.h>
#include <periodic_boundaries.h>
#include <dof_map.h>
#include <fe.h>
#include <elem.h>
#include <side.h>
#include <quadrature_gauss.h>
#include <boundary_info.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

#define isEqual(a, b) (fabs((a)/(b)-1.0) < 1e-5)


/*
 * Class name for error output.
 */

const string DiffusionSolver1D::_className("DiffusionSolver1D");


/*
 * Constructor.
 */

DiffusionSolver1D::
DiffusionSolver1D(const DiffusionSolver1DDescriptor &difDsc,
		  const DiffusionABCalculator &ab,
		  DiffusionSolutionHolder &dsh, double xl, double xr):
  _sysName(difDsc.getFileHeadStr()), _mesh(1), _es(_mesh),
  //_mRef(_mesh),
  _difDsc(difDsc), _ab(ab), _dsh(dsh),
  _Xl(xl), _Xr(xr)
{

  // Generate a uniform mesh. Each element has 3 interior points.

  MeshTools::Generation::build_line(_mesh, _difDsc.getNx()/3,
				    _Xl, _Xr, EDGE3);
  
  
  // Set up the equation system.

  TransientNonlinearImplicitSystem &sys = _es.add_system<TransientNonlinearImplicitSystem>(_sysName);

  sys.add_variable("mu_e", SECOND, LAGRANGE);
  sys.add_variable("mu_h", SECOND, LAGRANGE);


  // Do not use automatic calls for assembly function.

  sys.assemble_before_solve = false;


  // Set residual_and_jacobian_object.

  _rj = new ResidualAndJacobianDiffusion(_ab, _dsh, _es.get_mesh(),
					 sys.get_dof_map());
  sys.nonlinear_solver->residual_object = _rj;
  sys.nonlinear_solver->jacobian_object = _rj;

  
  // Initialize the equation systems.

  _es.init();


  // Setup for the periodic boundary.

  _setBoundaryID();

  DofMap &dofMap = sys.get_dof_map();
  PeriodicBoundary pbc(RealVectorValue(_difDsc.getLc()));

  pbc.myboundary = BoundaryIDLeft;
  pbc.pairedboundary = BoundaryIDRight;

  dofMap.add_periodic_boundary(pbc);


  // Setup some parameters in TransientNonlinearSystem.

  sys.time = 0.0;


  // Setup some parameters in EquationSystems.

  _es.parameters.set<Real>("dt") = _difDsc.get_dt();
  _es.parameters.set<Real>("nonlinear solver relative step tolerance") = 1e-5;
  _es.parameters.set<unsigned int>("nonlinear solver maximum iterations") = 100;

}


/*
 * Destructor.
 */

DiffusionSolver1D::~DiffusionSolver1D()
{
  if( _rj != (ResidualAndJacobianDiffusion *)NULL ){
    delete _rj;
  }
}


/*
 * Set IDs to boundary points.
 */

void DiffusionSolver1D::_setBoundaryID()
{
  TransientNonlinearImplicitSystem &sys = _es.get_system<TransientNonlinearImplicitSystem>(_sysName);
  const DofMap &dofMap = sys.get_dof_map();
  const unsigned int dim = _mesh.mesh_dimension();
  FEType feType = dofMap.variable_type(0);
  AutoPtr<FEBase> feFace(FEBase::build(dim, feType));
  QGauss qface(dim-1, FIFTH);
  feFace->attach_quadrature_rule(&qface);
  const std::vector<Point> &rf = feFace->get_xyz();

  MeshBase::const_element_iterator el = _mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_elements_end();
  
  
  for(; el!=end_el; ++el){
    const Elem *elem = *el;
    std::vector<unsigned int> dofInd;

    dofMap.dof_indices(elem, dofInd);

    for(unsigned int side=0; side<elem->n_sides(); side++){

      if(elem->neighbor(side) != NULL) continue;

      int id = _mesh.boundary_info->boundary_id(elem, side);

      feFace->reinit(elem, side);

      int boundaryID = BoundaryIDNone;

      for(unsigned int qp=0; qp<qface.n_points(); qp++){
	Real xf = rf[qp](0);

	if( isEqual(xf, _Xl) ){
	  boundaryID = BoundaryIDLeft;
	}
	else if( isEqual(xf, _Xr) ){
	  boundaryID = BoundaryIDRight;
	}
      }

      if(boundaryID != BoundaryIDNone){
	_mesh.boundary_info->add_side(elem, side, boundaryID);
      }
    }
  }
}


/*
 * Return RealSpaceGridHandler object.
 */

RealSpaceGridHandler DiffusionSolver1D::getRealSGH()
{
  TransientNonlinearImplicitSystem &sys = _es.get_system<TransientNonlinearImplicitSystem>(_sysName);
  const DofMap &dofMap = sys.get_dof_map();
  const unsigned int dim = _mesh.mesh_dimension();
  FEType feType = dofMap.variable_type(0);
  AutoPtr<FEBase> fe(FEBase::build(dim, feType));
  QGauss qrule(dim, FIFTH);
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Point> &r = fe->get_xyz();

  MeshBase::const_element_iterator el = _mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_elements_end();
  
  vector<vector<double> > x;

  x.resize(_mesh.n_active_elem());

  for(int iElem=0; el!=end_el; ++el, ++iElem){
    const Elem *elem = *el;
    std::vector<unsigned int> dofInd;

    dofMap.dof_indices(elem, dofInd);
    fe->reinit(elem);

    x[iElem].resize(qrule.n_points());

    for(unsigned int qp=0; qp<qrule.n_points(); qp++){
      x[iElem][qp] = r[qp](0);
    }
  }

  return RealSpaceGridHandler(x);
}
