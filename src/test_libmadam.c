// Include the autotools-provided configuration macros
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mpi.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "madam.h"

#define NDET 4
#define NPERIOD 4

int main( int argc, char **argv ) {

  MPI_Comm comm = MPI_COMM_WORLD;
  int flag;
  int ntasks=0, rank=0;
  int err;
  int nside=8;
  int npix=12*nside*nside;
  int fsample=32.5;

  char *parstring = "base_first=1.0;fsample=32.5;nside_map=8;nside_cross=8;nside_submap=8;write_map=T;write_binmap=T;write_matrix=T;write_wcov=T;write_hits=T;kfilter=T;path_output=./maps/";
  char *detstring = "LFI27M;LFI27S;LFI28M;LFI28S";
  long ndet=NDET;
  double detweights[NDET] = {1, 1, 1, 1};
  long nsamp = 1000;
  long nnz = 1;
  double timestamps[nsamp];
  long pixels[nsamp];
  double pixweights[nsamp*nnz];
  double signal[ndet*nsamp];
  long nperiod=NPERIOD;
  long periods[NPERIOD] = {0, nsamp/4, nsamp/2, 3*nsamp/4};
  long npsd[NDET] = {1, 1, 1, 1};
  double psdstarts[NDET] = {0, 0, 0, 0};
  long npsdtot=ndet;
  long npsdbin=10;
  double psdfreqs[npsdbin];
  long npsdval=npsdbin*npsdtot;
  double psdvals[npsdval];

  long i, j, k, l, kk;
  

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

  k = 0;
  kk = 0;
  for (i=0; i<nsamp; ++i) {
    timestamps[i] = i + rank*nsamp;
    for (j=0; j<ndet; ++j) {
      pixels[k] = (k+rank*nsamp*ndet) % npix;
      signal[k] = pixels[k];
      for (l=0; l<nnz; ++l) {
	pixweights[kk] = 1.0;
	++kk;
      }
      ++k;
    }
  }

  for ( i=0; i<npsdbin; ++i ) psdfreqs[i] = i * fsample / npsdbin;
  
  for ( i=0; i<npsdval; ++i) psdvals[i] = 1;

  destripe( (long)comm, parstring, ndet, detstring, detweights, nsamp, nnz, timestamps, pixels, pixweights, signal, nperiod, periods, npsd, npsdtot, psdstarts, npsdbin, psdfreqs, npsdval, psdvals );

  return 0;
}

