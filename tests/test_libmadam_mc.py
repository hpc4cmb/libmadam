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

import libmadam_wrapper as madam


class MadamMCTest(TestCase):

    np.random.seed(9876543)
    comm = MPI.COMM_WORLD
    itask = comm.Get_rank()
    ntask = comm.Get_size()

    try:
        import healpy as hp
    except:
        hp = None
        if itask == 0:
            print("Warning: unable to import healpy. Output maps are not checked.")

    fcomm = comm.py2f()

    if itask == 0:
        print("Running with ", ntask, " MPI tasks")

    print("Calling Madam")

    nside = 8
    npix = 12 * nside ** 2
    fsample = 32.5

    pars = {}
    pars["base_first"] = 1.0
    pars["nsubchunk"] = 1
    pars["fsample"] = fsample
    pars["nside_map"] = nside
    pars["nside_cross"] = nside
    pars["nside_submap"] = nside
    pars["write_map"] = True
    pars["write_binmap"] = True
    pars["write_matrix"] = True
    pars["write_wcov"] = True
    pars["write_hits"] = True
    pars["write_leakmatrix"] = True
    pars["kfilter"] = True
    pars["file_root"] = "madam_pytest"
    pars["path_output"] = "./mc_maps/"
    pars["mcmode"] = False

    if itask == 0:
        shutil.rmtree("mc_maps", ignore_errors=True)
        os.mkdir("mc_maps")
        for fn in ["mc_maps/madam_pytest_hmap.fits", "mc_maps/madam_pytest_bmap.fits"]:
            if os.path.isfile(fn):
                os.remove(fn)

    parstring = madam.dict2parstring(pars)

    dets = ["LFI27M", "LFI27S", "LFI28M", "LFI28S"]
    detstring = madam.dets2detstring(dets)

    ndet = len(dets)

    weights = np.ones(ndet, dtype="double")

    nsamp = 1000  # number of time ordered data samples

    nnz = 1  # number or non zero pointing weights, typically 3 for IQU

    timestamps = np.zeros(nsamp, dtype=madam.TIMESTAMP_TYPE)
    timestamps[:] = np.arange(nsamp) + itask * nsamp

    pixels = np.zeros(ndet * nsamp, dtype=madam.PIXEL_TYPE)
    pixels[:] = np.arange(len(pixels)) % npix

    pixweights = np.zeros(ndet * nsamp * nnz, dtype=madam.WEIGHT_TYPE)
    pixweights[:] = 1

    signal = np.zeros(ndet * nsamp, dtype=madam.SIGNAL_TYPE)
    signal[:] = pixels
    signal[:] += np.random.randn(nsamp * ndet)

    signal2 = np.zeros(ndet * nsamp, dtype=madam.SIGNAL_TYPE)
    signal2[:] = pixels % 12
    signal2[:] += np.random.randn(nsamp * ndet)

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

    signal_copy = signal.copy()

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

    # Basic Madam call with signal # 1

    madam.destripe(
        fcomm,
        parstring,
        ndet,
        detstring,
        weights,
        nsamp,
        nnz,
        timestamps,
        pixels,
        pixweights,
        signal,
        nperiod,
        periods,
        npsd,
        npsdtot,
        psdstarts,
        npsdbin,
        psdfreqs,
        npsdval,
        psdvals,
    )

    pars["mcmode"] = True
    parstring = madam.dict2parstring(pars)
    signal[:] = signal2

    # Basic Madam call with signal # 2, this time, cache the parameters, filters e.t.c.

    madam.destripe(
        fcomm,
        parstring,
        ndet,
        detstring,
        weights,
        nsamp,
        nnz,
        timestamps,
        pixels,
        pixweights,
        signal,
        nperiod,
        periods,
        npsd,
        npsdtot,
        psdstarts,
        npsdbin,
        psdfreqs,
        npsdval,
        psdvals,
    )

    signal[:] = signal_copy

    # MC Madam call again with signal # 1, This should produce the same result as the first call.

    outpath = "./mc_maps/2/"
    if itask == 0:
        if not os.path.isdir(outpath):
            os.mkdir(outpath)
            for fn in [
                "mc_maps/2/madam_pytest_hmap.fits",
                "mc_maps/2/madam_pytest_bmap.fits",
            ]:
                if os.path.isfile(fn):
                    os.remove(fn)
    outpath = outpath.encode("ascii")

    madam.destripe_with_cache(
        fcomm, ndet, nsamp, nnz, timestamps, pixels, pixweights, signal, outpath
    )

    # Test clearing the caches twice

    madam.clear_caches()

    madam.clear_caches()

    if itask == 0 and hp is not None:
        good = hmap_tot != 0
        bmap_tot[good] /= hmap_tot[good]

        hmap = hmap_tot.astype(np.int32)
        bmap = bmap_tot.astype(np.float32)

        try:
            hp.write_map("hits.fits", hmap, nest=True)
        except:
            hp.write_map("hits.fits", hmap, nest=True, overwrite=True)
        try:
            hp.write_map("binned.fits", bmap, nest=True)
        except:
            hp.write_map("binned.fits", bmap, nest=True, overwrite=True)

        madam_hmap = hp.read_map("mc_maps/madam_pytest_hmap.fits", nest=True)
        madam_bmap = hp.read_map("mc_maps/madam_pytest_bmap.fits", nest=True)

        good = hmap != 0

        hitdiff = np.std((madam_hmap - hmap)[good])
        bindiff = np.std((madam_bmap - bmap)[good])

        if hitdiff != 0:
            print("Hit map check FAILED: hit map difference RMS ", hitdiff)
            sys.exit(-1)
        else:
            print("Hit map check PASSED")

        if bindiff != 0:
            print("Binned map check FAILED: Binned map difference RMS ", bindiff)
            sys.exit(-1)
        else:
            print("Binned map check PASSED")

    if itask == 0:
        shutil.rmtree("mc_maps")

    print("Done")
