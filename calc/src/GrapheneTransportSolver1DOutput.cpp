#include "GrapheneTransportSolver1D.h"
#include <Graphene/Constants.h>
#include <stdio.h>

using namespace std;


/*
 * Output the potential distribution in .pos format (and the mesh
 * in .msh format if toOutputMesh=true).
 */

void GrapheneTransportSolver1D::
outputPotential(const char *dname, const char *fhead,
		bool toOutputMesh, const string &format) const
{
  _poisson.outputPotential2D(dname, fhead, toOutputMesh, format);
}


/*
 * Output the electron concentration in the 2DEG.
 * Note that this function generates files with extension ".dat".
 */

void GrapheneTransportSolver1D::
outputConcentration2DEG(const char *dname, const char *fhead) const
{
  //if( !_gridParam.isOutputProc() ) return;

  FILE *fp;
  char filename[300];

  sprintf(filename, "%s/%s-e.dat", dname, fhead);
  fp = fopen(filename, "w");

  for(int i=0; i<_realSGH.getSize(); i++){
    fprintf(fp, "%g %g\n", _realSGH.getAt(i),
	    _SigmaElectron.getAt(i));
  }

  fclose(fp);

  sprintf(filename, "%s/%s-h.dat", dname, fhead);
  fp = fopen(filename, "w");
    
  for(int i=0; i<_realSGH.getSize(); i++){
    fprintf(fp, "%g %g\n", _realSGH.getAt(i), _SigmaHole.getAt(i));
  }
  
  fclose(fp);
}

void GrapheneTransportSolver1D::
outputConcentration2DEG(const char *dname, const char *fhead,
			const RealSpaceArrayDiffusion &se,
			const RealSpaceArrayDiffusion &sh) const
{
  //if( !_gridParam.isOutputProc() ) return;

  FILE *fp;
  char filename[300];

  sprintf(filename, "%s/%s-e.dat", dname, fhead);
  fp = fopen(filename, "w");

  for(int i=0; i<_realSGH.getSize(); i++){
    fprintf(fp, "%g %g\n", _realSGH.getAt(i), se.getAt(i));
  }

  fclose(fp);

  sprintf(filename, "%s/%s-h.dat", dname, fhead);
  fp = fopen(filename, "w");
  
  for(int i=0; i<_realSGH.getSize(); i++){
    fprintf(fp, "%g %g\n", _realSGH.getAt(i), sh.getAt(i));
  }
  
  fclose(fp);
}


/*
 * Output the electron concentration in the 2DEG in binary format.
 * Note that this function generates files with extension ".bin".
 */

void GrapheneTransportSolver1D::
outputConcentration2DEGBin(const char *dname, const char *fhead) const
{
  //if( !_gridParam.isOutputProc() ) return;

  FILE *fp;
  char filename[300];

  sprintf(filename, "%s/%s-e.bin", dname, fhead);
  fp = fopen(filename, "wb");

  vector<double> buf(_realSGH.getSize());

  for(int i=0; i<_realSGH.getSize(); i++){
    buf[i] = _SigmaElectron.getAt(i);
  }

  fwrite(&buf[0], sizeof(double), buf.size(), fp);
  fclose(fp);

  sprintf(filename, "%s/%s-h.bin", dname, fhead);
  fp = fopen(filename, "wb");

  for(int i=0; i<_realSGH.getSize(); i++){
    buf[i] = _SigmaHole.getAt(i);
  }

  fwrite(&buf[0], sizeof(double), buf.size(), fp);
  fclose(fp);
}


/*
 * Output the potential in the 2DEG.
 * Note that this function generates files with extension ".dat".
 */
/*
void GrapheneTransportSolver1D::
outputPotential2DEG(const char *dname, const char *fhead) const
{
  FILE *fp;
  char filename[300];
  Potential pot(_realSGH);

  // Gradual channel approximation.

  if(_poiDsc.getFieldModel() == FieldModelCGA){

    // With bottom gates is not supported.
    assert(_poiDsc.getBottomGates().size() == 0);
    
    for(int i=0; i<pot.getSize(); i++){
      pot.setAt(i, _pot1DCGA[i]);
    }
  }

  // SCF potential.

  else if(_poiDsc.getFieldModel() == FieldModelSCF){
    _poisson.calcPotential2DEG(pot);
  }

  
  sprintf(filename, "%s/%s.dat", dname, fhead);
  fp = fopen(filename, "w");

  for(int i=0; i<_Ex.getSize(); i++){
    fprintf(fp, "%g %g\n", _realSGH.getXiPoisson(i),
	    pot.getAt(i));
  }

  fclose(fp);
}
*/

void GrapheneTransportSolver1D::
outputPotential2DEG(const char *dname, const char *fhead,
		    const RealSpaceArrayDiffusion &pot) const
{
  //if( !_gridParam.isOutputProc() ) return;

  FILE *fp;
  char filename[300];

  sprintf(filename, "%s/%s.dat", dname, fhead);
  fp = fopen(filename, "w");

  for(int i=0; i<pot.getSize(); i++){
    fprintf(fp, "%g %g\n", _realSGH.getAt(i), pot.getAt(i));
  }

  fclose(fp);
}


/*
 * Output the field in the 2DEG.
 * Note that this function generates files with extension ".dat".
 */

void GrapheneTransportSolver1D::
outputField2DEG(const char *dname, const char *fhead) const
{
  //if( !_gridParam.isOutputProc() ) return;

  FILE *fp;
  char filename[300];

  sprintf(filename, "%s/%s.dat", dname, fhead);
  fp = fopen(filename, "w");

  for(int i=0; i<_Ex.getSize(); i++){
    fprintf(fp, "%g %g\n", _realSGH.getAt(i), _Ex.getAt(i));
  }

  fclose(fp);
}

void GrapheneTransportSolver1D::
outputField2DEG(const char *dname, const char *fhead,
		const RealSpaceArrayDiffusion &field) const
{
  //if( !_gridParam.isOutputProc() ) return;

  FILE *fp;
  char filename[300];

  sprintf(filename, "%s/%s.dat", dname, fhead);
  fp = fopen(filename, "w");

  for(int i=0; i<field.getSize(); i++){
    fprintf(fp, "%g %g\n", _realSGH.getAt(i), field.getAt(i));
  }

  fclose(fp);
}


/*
 * Output the velocity field $v_{x}(x)$ in the 2DEG.
 * Note that this function generates files with extension ".dat".
 */

void GrapheneTransportSolver1D::
outputVelocity(const char *dname, const char *fhead)
{
  /*
  FILE *fp;
  char filename[300];
  RealSpaceArrayWENO vel(_realSGH);

  // Calculation is parallel but output is serial.

  _bltElectron.calcVelocity(vel);

  if( _gridParam.isOutputProc() ){
    sprintf(filename, "%s/%s-e.dat", dname, fhead);
    fp = fopen(filename, "w");
    
    for(int i=0; i<vel.getSizeBlt(); i++){
      fprintf(fp, "%g %e\n", _realSGH.getXiBlt(i), vel.getAtBlt(i));
    }
    
    fclose(fp);
  }

  if( _bltDsc.toAccountHole() ){
    _bltHole.calcVelocity(vel);

    if( _gridParam.isOutputProc() ){
      sprintf(filename, "%s/%s-h.dat", dname, fhead);
      fp = fopen(filename, "w");
      
      for(int i=0; i<vel.getSizeBlt(); i++){
	fprintf(fp, "%g %e\n", _realSGH.getXiBlt(i), vel.getAtBlt(i));
      }
      
      fclose(fp);
    }
  }
  */
}
