// Include the autotools-provided configuration macros
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mpi.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "madam.h"

int main( int argc, char **argv ) {

  MPI_Comm comm = MPI_COMM_WORLD;
  int flag;
  int ntasks=0, rank=0;
  int err;

  MPI_Initialized( &flag );
  if ( flag ) {
    printf( "ERROR: MPI was already initialized\n" );
    return -1;
  }
  
  if ( MPI_Init( &argc, &argv ) ) {
    printf( "ERROR: Failed to initialize MPI\n" );
    return -1;
  }

  if ( MPI_Comm_size( comm, &ntasks ) ) {
    printf( "ERROR: Failed get MPI communicator size\n" );
    return -1;
  }

  if ( MPI_Comm_rank( comm, &rank ) ) {
    printf( "ERROR: Failed to get MPI rank\n" );
    return -1;
  }

  if ( MPI_Finalize() ) {
    printf( "ERROR: Failed finalize MPI\n" );
    return -1;
  }

  // FIXME: add the actual test here

  return 0;
}

