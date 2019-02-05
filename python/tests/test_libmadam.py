# Copyright (c) 2019 by the parties listed in the AUTHORS
# file.  All rights reserved.  Use of this source code is governed
# by a BSD-style license that can be found in the LICENSE file.

from __future__ import division
from __future__ import print_function

from mpi4py import MPI

import ctypes as ct
import os
import sys
import shutil
from unittest import TestCase

import numpy as np
import numpy.testing as npt

import libmadam_wrapper as madam


class MadamTest(TestCase):

    np.random.seed(123456)
    comm = MPI.COMM_WORLD
    itask = comm.Get_rank()
    ntask = comm.Get_size()

    try:
        import healpy as hp
    except:
        hp = None
        if itask == 0:
            print("Warning: unable to import healpy. Output maps are not checked.")

    if itask == 0:
        print("Running with ", ntask, " MPI tasks")

    print("Calling Madam")

    nside = 8
    npix = 12 * nside ** 2
    fsample = 32.5
    nsamp = 1000  # number of time ordered data samples
    nnz = 1  # number or non zero pointing weights, typically 3 for IQU

    pars = {}
    pars["base_first"] = 1.0
    pars["fsample"] = fsample
    pars["nside_map"] = nside
    pars["nside_cross"] = nside // 2
    pars["nside_submap"] = nside // 4
    pars["write_map"] = True
    pars["write_binmap"] = True
    pars["write_matrix"] = True
    pars["write_wcov"] = True
    pars["write_hits"] = True
    pars["write_leakmatrix"] = True
    pars["kfilter"] = True
    pars["diagfilter"] = 0
    pars["file_root"] = "madam_pytest"
    pars["path_output"] = "./pymaps/"
    pars["iter_max"] = 100
    pars["nsubchunk"] = 2
    pars["allreduce"] = True
    pars["info"] = 0

    # pars[ 'detset' ] = ['LFI27 : LFI27M, LFI27S',
    #                    'LFI28 : LFI28M, LFI28S']
    pars["survey"] = [
        "hm1 : {} - {}".format(0, nsamp / 2),
        #'hm2 : {} - {}'.format(nsamp/2, nsamp),
        #'odd : {} - {}, {} - {}'.format(0, nsamp/4, nsamp/2, 3*nsamp/4),
        #'even : {} - {}, {} - {}'.format(nsamp/4, nsamp/2, 3*nsamp/4, nsamp)
    ]
    pars["bin_subsets"] = True

    if itask == 0:
        shutil.rmtree("pymaps", ignore_errors=True)
        os.mkdir("pymaps")
        for fn in ["pymaps/madam_pytest_hmap.fits", "pymaps/madam_pytest_bmap.fits"]:
            if os.path.isfile(fn):
                print("Removing old {}".format(fn))
                os.remove(fn)

    dets = ["LFI27M", "LFI27S", "LFI28M", "LFI28S"]

    ndet = len(dets)

    weights = np.ones(ndet, dtype=np.float64)

    timestamps = np.zeros(nsamp, dtype=madam.TIMESTAMP_TYPE)
    timestamps[:] = np.arange(nsamp) + itask * nsamp

    pixels = np.zeros(ndet * nsamp, dtype=madam.PIXEL_TYPE)
    pixels[:] = np.arange(len(pixels)) % npix

    pixweights = np.zeros(ndet * nsamp * nnz, dtype=madam.WEIGHT_TYPE)
    pixweights[:] = 1

    signal = np.zeros(ndet * nsamp, dtype=madam.SIGNAL_TYPE)
    signal[:] = pixels
    signal[:] += np.random.randn(nsamp * ndet)

    nperiod = 4  # number of pointing periods

    periods = np.zeros(nperiod, dtype=np.int64)
    periods[1] = int(nsamp * 0.25)
    periods[2] = int(nsamp * 0.50)
    periods[3] = int(nsamp * 0.75)

    npsd = np.ones(ndet, dtype=np.int64)
    npsdtot = np.sum(npsd)
    psdstarts = np.zeros(npsdtot)
    npsdbin = 10
    psdfreqs = np.arange(npsdbin) * fsample / npsdbin
    npsdval = npsdbin * npsdtot
    psdvals = np.ones(npsdval)

    # Reference maps for checking

    hmap = np.zeros(npix, dtype=np.int64)
    bmap = np.zeros(npix, dtype=np.float64)

    for p, s in zip(pixels, signal):
        hmap[p] += 1
        bmap[p] += s

    hmap_tot = np.zeros(npix, dtype=np.int64)
    bmap_tot = np.zeros(npix, dtype=np.float64)

    comm.Reduce(hmap, hmap_tot, op=MPI.SUM, root=0)
    comm.Reduce(bmap, bmap_tot, op=MPI.SUM, root=0)

    # Ensure we can successfully call Madam twice with different inputs
    for i in range(2):
        madam.destripe(
            comm,
            pars,
            dets,
            weights,
            timestamps,
            pixels,
            pixweights,
            signal,
            periods,
            npsd,
            psdstarts,
            psdfreqs,
            psdvals,
        )
        nside *= 2
        pars["nside_map"] = nside

    if itask == 0 and hp is not None:
        good = hmap_tot != 0
        bmap_tot[good] /= hmap_tot[good]

        hmap = hmap_tot.astype(np.int32)
        bmap = bmap_tot.astype(np.float32)

        # hp.write_map("hits.fits", hmap, nest=True, overwrite=True)
        # hp.write_map("binned.fits", bmap, nest=True, overwrite=True)

        madam_hmap = hp.read_map("pymaps/madam_pytest_hmap.fits", nest=True)
        madam_bmap = hp.read_map("pymaps/madam_pytest_bmap.fits", nest=True)

        npt.assert_allclose(madam_hmap, hmap)
        npt.assert_allclose(madam_bmap, bmap)

    if itask == 0:
        shutil.rmtree("pymaps")

    print("Done")
