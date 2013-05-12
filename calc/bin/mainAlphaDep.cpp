#include "DiffusionSolver1DDescriptor.h"
#include "GrapheneTransportSolver1D.h"
#include <PhysicalUnits.h>
#include <Transistor2D/PoissonSolver2DDescriptor.h>
#include <OutputDirectoryManager/OutputDirectoryManager.h>
#include <OutputDirectoryManager/ValueStringGenerator.h>
#include "libmesh.h"
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fparser.hh>

using namespace std;


void divideMPICommWorld(MPI_Comm &world, int &nrWorlds, int &worldID);


/*
 * Main.
 */

int main(int argc, char *argv[])
{
  // Extract the descriptor name (without extension).
  
  if(argc < 3){
    cerr << "main: Provide the descriptor names (without extension)";
    cerr << endl;
    exit(1);
  }

  string dscNamePoi = argv[1];
  string dscNameDif = argv[2];


  // Initialize MPI.

  MPI_Init(&argc, &argv);


  // Initialize system-descriptor.

  PoissonSolver2DDescriptor poiDsc(dscNamePoi);
  DiffusionSolver1DDescriptor difDsc(dscNameDif);
  

  // Initialize LibMesh.

  MPI_Comm world;
  int nrWorlds, worldID;

  divideMPICommWorld(world, nrWorlds, worldID);

  //LibMeshInit init(argc, argv, MPI_COMM_WORLD);
  LibMeshInit init(argc, argv, world);


  // Create parent directories for output.

  OutputDirectoryManager odm("../dat/alphaDep");

  odm.push(ValueStringGenerator("alpha", "default").getValueString(difDsc.get_alpha()));
  odm.push(ValueStringGenerator("Vg", "default").getValueString(poiDsc.getGates()[0]->getVoltage(-fs2s(100))));


  // Initialize the solver.

  GrapheneTransportSolver1D solver(poiDsc, difDsc, odm);


  // Time step.

  while( solver.getTime() <= difDsc.get_tMax() ){
    solver.solveStep();
  }

  return 0;
}

void divideMPICommWorld(MPI_Comm &world, int &nrWorlds, int &worldID)
{
  const int nrProcsPerWorld = 1;
  int size, rank;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Stay calm in one world.

  if(size < nrProcsPerWorld){
    world = MPI_COMM_WORLD;
    nrWorlds = 1;
    worldID = 0;

    return;
  }


  // # procs must be multiple of 'nrProcsPerWorld'.

  else if (size%nrProcsPerWorld != 0){
    if(rank == 0){
      fprintf(stderr, "main: # procs must be smaller than 16 ");
      fprintf(stderr, "or multiple of %d.\n", nrProcsPerWorld);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    exit(1);
  }


  // Divide the world into pieces!

  nrWorlds = size/nrProcsPerWorld;
  worldID = rank/nrProcsPerWorld;

  MPI_Group group, subGroup;
  vector<int> ranksInSubGroup(nrProcsPerWorld);
  for(int i=0; i<nrProcsPerWorld; i++){
    ranksInSubGroup[i] = nrProcsPerWorld*worldID+i;
  }

  MPI_Comm_group(MPI_COMM_WORLD, &group);
  MPI_Group_incl(group, nrProcsPerWorld, &ranksInSubGroup[0],
                 &subGroup);
  MPI_Comm_create(MPI_COMM_WORLD, subGroup, &world);
}
