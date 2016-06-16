#!/usr/bin/env python

import ctypes
import os
import sys

_libdir = os.path.dirname( __file__ )

path = os.path.join( _libdir, '.libs', 'libmadam.so' )
try:
    #_madam = ctypes.CDLL( path )
except:
    try:
        #_madam = ctypes.CDLL( path.replace('.so','.dylib') )
    except:
        print('Failed to load libmadam at ', path)
        raise

def dict2parstring( d ):

    s = ''
    for key, value in d.items():
        s += '{} = {};'.format( key, value )

    return s

def dets2detstring( dets ):

    s = ''
    for d in dets:
        s += '{};'.format( d )

    return s

if __name__ == '__main__':

    """
    When madam.py is ran directly, it performs a quick check of the code by making a hit map and a binned map
    """

    from mpi4py import MPI
    import numpy as np
    import healpy as hp
    
    comm = MPI.COMM_WORLD
    itask = comm.Get_rank()
    ntask = comm.Get_size()

    fcomm = comm.py2f()

    if itask == 0: print('Running with ', ntask, ' MPI tasks')

    print('Calling Madam')

    nside = 8
    npix = 12*nside**2
    fsample = 32.5

    pars = {}
    pars[ 'base_first' ] = 1.0
    pars[ 'fsample' ] = fsample
    pars[ 'nside_map' ] = nside
    pars[ 'nside_cross' ] = nside
    pars[ 'nside_submap' ] = nside
    pars[ 'write_map' ] = True
    pars[ 'write_binmap' ] = True
    pars[ 'write_matrix' ] = True
    pars[ 'write_wcov' ] = True
    pars[ 'write_hits' ] = True
    pars[ 'kfilter' ] = True
    pars[ 'path_output' ] = './maps/'

    parstring = dict2parstring( pars )

    dets = ['LFI27M', 'LFI27S', 'LFI28M', 'LFI28S']
    detstring = dets2detstring( dets )

    ndet = len(dets)

    weights = np.ones( ndet, dtype='double' )

    nsamp = 1000 # number of time ordered data samples

    nnz = 1 # number or non zero pointing weights, typically 3 for IQU

    timestamps = np.zeros( nsamp, dtype='double' )
    timestamps[:] = np.arange( nsamp ) + itask*nsamp

    pixels = np.zeros( ndet * nsamp, dtype=np.int64 )
    pixels[:] = np.arange( len(pixels) ) % npix

    pixweights = np.zeros( ndet * nsamp * nnz, dtype='double' )
    pixweights[:] = 1

    signal = np.zeros( ndet * nsamp, dtype='double' )
    signal[:] = pixels
    signal[:] += np.random.randn( nsamp * ndet )

    nperiod = 4 # number of pointing periods

    periods = np.zeros( nperiod, dtype=np.int64 )
    periods[1] = int( nsamp * .25 )
    periods[2] = int( nsamp * .50 )
    periods[3] = int( nsamp * .75 )

    npsd = np.ones( ndet, dtype=np.int64 )
    npsdtot = np.sum( npsd )
    psdstarts = np.zeros( npsdtot )
    npsdbin = 10
    psdfreqs = np.arange(npsdbin) * fsample / npsdbin
    npsdval = npsdbin * npsdtot
    psdvals = np.ones( npsdval )

    for i in range(2): # Ensure we can successfully call Madam twice with different inputs
        _madam.destripe(
            fcomm,
            ctypes.c_char_p( parstring ),
            ctypes.c_long( ndet ),
            ctypes.c_char_p( detstring ),
            weights.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
            ctypes.c_long(nsamp),
            ctypes.c_long(nnz),
            timestamps.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
            pixels.ctypes.data_as( ctypes.POINTER( ctypes.c_long ) ),
            pixweights.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
            signal.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
            ctypes.c_long(nperiod),
            periods.ctypes.data_as( ctypes.POINTER( ctypes.c_long ) ),
            npsd.ctypes.data_as( ctypes.POINTER( ctypes.c_long ) ),
            ctypes.c_long(npsdtot),
            psdstarts.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
            ctypes.c_long(npsdbin),
            psdfreqs.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
            ctypes.c_long(npsdval),
            psdvals.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
        )
        nside *= 2
        pars[ 'nside_map' ] = nside
        parstring = dict2parstring( pars )

    hmap = np.zeros( npix, dtype=np.int64 )
    bmap = np.zeros( npix, dtype=np.float64 )

    for p, s in zip( pixels, signal ):
        hmap[ p ] += 1
        bmap[ p ] += s

    hmap_tot = np.zeros( npix, dtype=np.int64 )
    bmap_tot = np.zeros( npix, dtype=np.float64 )    

    comm.Reduce( hmap, hmap_tot, op=MPI.SUM, root=0 )
    comm.Reduce( bmap, bmap_tot, op=MPI.SUM, root=0 )

    if itask == 0:
        good = hmap_tot != 0
        bmap_tot[ good ] /= hmap_tot[ good ]

        hmap = hmap_tot.astype( np.int32 )
        bmap = bmap_tot.astype( np.float32 )

        hp.write_map( 'hits.fits', hmap, nest=True )
        hp.write_map( 'binned.fits', bmap, nest=True )

        madam_hmap = hp.read_map( 'maps/madam_hmap.fits', nest=True )
        madam_bmap= hp.read_map( 'maps/madam_bmap.fits', nest=True )

        good = madam_bmap != hp.UNSEEN

        hitdiff = np.std( (madam_hmap - hmap)[good] )
        bindiff = np.std( (madam_bmap - bmap)[good] )

        if hitdiff != 0:
            print('Hit map check FAILED: hit map difference RMS ', hitdiff)
        else:
            print('Hit map check PASSED')

        if bindiff != 0:
            print('Binned map check FAILED: Binned map difference RMS ', bindiff)
        else:
            print('Binned map check PASSED')

    print('Done')

