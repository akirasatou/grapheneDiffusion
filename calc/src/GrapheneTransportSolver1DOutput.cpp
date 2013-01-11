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
  /*
  if(_bltDsc.getFieldModel() == SCF){
    _poisson.outputPotential2D(dname, fhead, toOutputMesh, format);
  }
  else if(_bltDsc.getFieldModel() == GCA){
    char filename[1000];

    sprintf(filename, "%s/%s.dat", dname, fhead);

    FILE *fp = fopen(filename, "w");

    fprintf(fp, "FieldModelGCA does not have 2D potential.\n");
    fclose(fp);
  }
  */
}


/*
 * Output the electron concentration in the 2DEG.
 * Note that this function generates files with extension ".dat".
 */

void GrapheneTransportSolver1D::
outputConcentration2DEG(const char *dname, const char *fhead) const
{
  /*
  //if( !_gridParam.isOutputProc() ) return;

  FILE *fp;
  char filename[300];

  sprintf(filename, "%s/%s-e.dat", dname, fhead);
  fp = fopen(filename, "w");

  for(int i=0; i<_realSGH.getSizeBlt(); i++){
    fprintf(fp, "%g %.15e\n", _realSGH.getXiBlt(i),
	    _SigmaElectron.getAtBlt(i));
  }

  fclose(fp);

  if( _bltDsc.toAccountHole() ){
    sprintf(filename, "%s/%s-h.dat", dname, fhead);
    fp = fopen(filename, "w");
    
    for(int i=0; i<_realSGH.getSizeBlt(); i++){
      fprintf(fp, "%g %g\n", _realSGH.getXiBlt(i),
	      _SigmaHole.getAtBlt(i));
    }
    
    fclose(fp);
  }
  */
}

/*
void GrapheneTransportSolver1D::
outputConcentration2DEG(const char *dname, const char *fhead,
			const Concentration &se,
			const Concentration &sh) const
{
  if( !_gridParam.isOutputProc() ) return;

  FILE *fp;
  char filename[300];

  sprintf(filename, "%s/%s-e.dat", dname, fhead);
  fp = fopen(filename, "w");

  for(int i=0; i<_realSGH.getSizeBlt(); i++){
    fprintf(fp, "%g %g\n", _realSGH.getXiBlt(i),
	    se.getAtBlt(i));
  }

  fclose(fp);

  if( _bltDsc.toAccountHole() ){
    sprintf(filename, "%s/%s-h.dat", dname, fhead);
    fp = fopen(filename, "w");
    
    for(int i=0; i<_realSGH.getSizeBlt(); i++){
      fprintf(fp, "%g %g\n", _realSGH.getXiBlt(i),
	      sh.getAtBlt(i));
    }
    
    fclose(fp);
  }
}
*/

/*
 * Output the electron concentration in the 2DEG in binary format.
 * Note that this function generates files with extension ".bin".
 */

void GrapheneTransportSolver1D::
outputConcentration2DEGBin(const char *dname, const char *fhead) const
{
  /*
  //if( !_gridParam.isOutputProc() ) return;

  FILE *fp;
  char filename[300];

  sprintf(filename, "%s/%s-e.bin", dname, fhead);
  fp = fopen(filename, "wb");

  vector<double> buf(_realSGH.getSizeBlt());

  for(int i=0; i<_realSGH.getSizeBlt(); i++){
    buf[i] = _SigmaElectron.getAtBlt(i);
  }

  fwrite(&buf[0], sizeof(double), buf.size(), fp);
  fclose(fp);

  if( _bltDsc.toAccountHole() ){
    sprintf(filename, "%s/%s-h.bin", dname, fhead);
    fp = fopen(filename, "wb");

    for(int i=0; i<_realSGH.getSizeBlt(); i++){
      buf[i] = _SigmaHole.getAtBlt(i);
    }

    fwrite(&buf[0], sizeof(double), buf.size(), fp);
    fclose(fp);
  }
  */
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
 /*
void GrapheneTransportSolver1D::
outputPotential2DEG(const char *dname, const char *fhead,
		    const Potential &pot) const
{
  if( !_gridParam.isOutputProc() ) return;

  FILE *fp;
  char filename[300];

  sprintf(filename, "%s/%s.dat", dname, fhead);
  fp = fopen(filename, "w");

  for(int i=0; i<pot.getSize(); i++){
    fprintf(fp, "%g %g\n", _realSGH.getXiPoisson(i),
	    pot.getAt(i));
  }

  fclose(fp);
}
 */


/*
 * Output the field in the 2DEG.
 * Note that this function generates files with extension ".dat".
 */

void GrapheneTransportSolver1D::
outputField2DEG(const char *dname, const char *fhead) const
{
  /*
  if( !_gridParam.isOutputProc() ) return;

  FILE *fp;
  char filename[300];

  sprintf(filename, "%s/%s.dat", dname, fhead);
  fp = fopen(filename, "w");

  for(int i=0; i<_Ex.getSize(); i++){
    fprintf(fp, "%g %g\n", _realSGH.getXiPoisson(i),
	    _Ex.getAt(i));
  }

  fclose(fp);
  */
}

/*
void GrapheneTransportSolver1D::
outputField2DEG(const char *dname, const char *fhead,
		const Field &field) const
{
  if( !_gridParam.isOutputProc() ) return;

  FILE *fp;
  char filename[300];

  sprintf(filename, "%s/%s.dat", dname, fhead);
  fp = fopen(filename, "w");

  for(int i=0; i<field.getSize(); i++){
    fprintf(fp, "%g %g\n", _realSGH.getXiPoisson(i),
	    field.getAt(i));
  }

  fclose(fp);
}
*/


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
