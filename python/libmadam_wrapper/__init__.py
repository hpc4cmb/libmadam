# Copyright (c) 2019 by the parties listed in the AUTHORS
# file.  All rights reserved.  Use of this source code is governed
# by a BSD-style license that can be found in the LICENSE file.

from __future__ import division
from __future__ import print_function

import ctypes as ct
import ctypes.util as ctu
import os
import sys
import glob
import shutil

import numpy as np
import numpy.ctypeslib as npc

TIMESTAMP_TYPE = np.float64
SIGNAL_TYPE = np.float64
PIXEL_TYPE = np.int32
WEIGHT_TYPE = np.float32
PSD_TYPE = np.float64

REPEATED_KEYS = ["detset", "detset_nopol", "survey"]


try:
    _madam = ct.CDLL("libmadam.so")
except OSError:
    path = ctu.find_library("madam")
    if path is not None:
        _madam = ct.CDLL(path)

available = _madam is not None


_madam.destripe.restype = None
_madam.destripe.argtypes = [
    ct.c_int,  # fcomm
    ct.c_char_p,  # parstring
    ct.c_long,  # ndet
    ct.c_char_p,  # detstring
    npc.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),
    ct.c_long,  # nsamp
    ct.c_long,  # nnz
    npc.ndpointer(dtype=TIMESTAMP_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    npc.ndpointer(dtype=PIXEL_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    npc.ndpointer(dtype=WEIGHT_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    npc.ndpointer(dtype=SIGNAL_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    ct.c_long,  # nperiod
    npc.ndpointer(dtype=np.int64, ndim=1, flags="C_CONTIGUOUS"),
    npc.ndpointer(dtype=np.int64, ndim=1, flags="C_CONTIGUOUS"),
    ct.c_long,  # npsdtot
    npc.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),
    ct.c_long,  # npsdbin
    npc.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),
    ct.c_long,  # npsdval
    npc.ndpointer(dtype=PSD_TYPE, ndim=1, flags="C_CONTIGUOUS"),
]


_madam.destripe_with_cache.restype = None
_madam.destripe_with_cache.argtypes = [
    ct.c_int,  # fcomm
    ct.c_long,  # ndet
    ct.c_long,  # nsamp
    ct.c_long,  # nnz
    npc.ndpointer(dtype=TIMESTAMP_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    npc.ndpointer(dtype=PIXEL_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    npc.ndpointer(dtype=WEIGHT_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    npc.ndpointer(dtype=SIGNAL_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    ct.c_char_p,  # outpath
]


def dict2parstring(d):
    """
    Translate parameter dictionary into a parameter string
    """
    s = ""
    for key, value in d.items():
        if key in REPEATED_KEYS:
            for separate_value in value:
                s += "{} = {};".format(key, separate_value)
        else:
            s += "{} = {};".format(key, value)
    return s.encode("ascii")


def dets2detstring(dets):
    """
    Translate a container of detector names into a detector name string
    """
    s = ""
    for d in dets:
        s += "{};".format(d)
    return s.encode("ascii")


def destripe(
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
):
    """
    Solve and subtract baseline offsets from signal
    """
    if not available:
        raise RuntimeError("No libmadam available, cannot destripe")
    fcomm = comm.py2f()
    parstring = dict2parstring(pars)
    ndet = len(dets)
    detstring = dets2detstring(dets)
    nsamp = timestamps.size
    nperiod = periods.size
    nnz = pixweights.size // nsamp // ndet
    npsdbin = psdfreqs.size
    npsdval = psdvals.size
    npsdtot = npsdval // npsdbin
    _madam.destripe(
        fcomm,
        parstring,
        ndet,
        detstring,
        weights,
        nsamp,
        nnz,
        np.ascontiguousarray(timestamps.reshape(-1), dtype=TIMESTAMP_TYPE),
        np.ascontiguousarray(pixels.reshape(-1), dtype=PIXEL_TYPE),
        np.ascontiguousarray(pixweights.reshape(-1), dtype=WEIGHT_TYPE),
        np.ascontiguousarray(signal.reshape(-1), dtype=SIGNAL_TYPE),
        nperiod,
        periods,
        npsd,
        npsdtot,
        psdstarts,
        npsdbin,
        psdfreqs,
        npsdval,
        np.ascontiguousarray(psdvals.reshape(-1), dtype=PSD_TYPE),
    )

    return


def destripe_with_cache(comm, timestamps, pixels, pixweights, signal, outpath):
    """
    Solve and subtract baseline offsets from signal in Monte Carlo mode
    """
    if not available:
        raise RuntimeError("No libmadam available, cannot destripe")
    fcomm = comm.py2f()
    nsamp = timestamps.size
    ndet = signal.size // nsamp
    nnz = pixweights.size // nsamp // ndet
    _madam.destripe_with_cache(
        fcomm,
        ndet,
        nsamp,
        nnz,
        np.ascontiguousarray(timestamps.reshape(-1), dtype=TIMESTAMP_TYPE),
        np.ascontiguousarray(pixels.reshape(-1), dtype=PIXEL_TYPE),
        np.ascontiguousarray(pixweights.reshape(-1), dtype=WEIGHT_TYPE),
        np.ascontiguousarray(signal.reshape(-1), dtype=SIGNAL_TYPE),
        outpath,
    )
    return


def clear_caches():
    """
    Clear Monte Carlo caches after a destripe() call with mcmode=True
    """
    _madam.clear_caches()
