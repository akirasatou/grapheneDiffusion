#include "DiffusionSolver1DDescriptor.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <PhysicalConstants.h>
#include <PhysicalUnits.h>
#include <iostream>
#include <fstream>

using namespace std;


/*
 * Class name for error output.
 */

const string DiffusionSolver1DDescriptor::
className("DiffusionSolver1DDescriptor");


/*
 * Constructor.
 */

DiffusionSolver1DDescriptor::
DiffusionSolver1DDescriptor(const string &fileHead):
  _fileHeadStr(fileHead), _sysDsc(fileHead+".dif")
{
  // Read the descriptor.

  _readDsc();
}


/*
 * Destructor.
 */

DiffusionSolver1DDescriptor::~DiffusionSolver1DDescriptor()
{
}
