#!/bin/python
# -*- coding: utf-8 -*-

import ctypes
import os
import platform
from numpy.ctypeslib import ndpointer
from numpy import zeros, fromstring

plat_info = dict(plat=platform.system())
if plat_info['plat'] == 'Windows':
    plat_info['lib'] = './autocorrelation_shared.dll'
    plat_info['com'] = '(mingw-w64 under cygwin) x86_64-w64-mingw32-gcc.exe -std=c99 -O3 -fPIC -shared autocorrelation.c -o autocorrelation_shared.dll -Wall -lmpfr -lgmp -fopenmp'
else:
    plat_info['lib'] = './autocorrelation_shared.so'
    plat_info['com'] = 'gcc -O3 -fPIC -shared autocorrelation.c -o autocorrelation_shared.so -Wall -lmpfr -lgmp -fopenmp'


if not os.path.isfile(plat_info['lib']):
    raise IOError("{lib} is missing. To compile on {plat}:\n{com}\n".format(**plat_info))

lib = ctypes.cdll[plat_info['lib']]

def aCorrUpTo(x, k, n=None):
    """
    First k lags of x autocorrelation with optional nth bit bitmask.

    If n is None, it calculates it on the whole bytes. 
    Otherwise it calculates it on the nth bit of each byte.
    """
    assert n in [None]+range(0,8) and k>0, \
           "Invalid n or k. Condition is: n in [None]+range(0,8) and k>0"
    n = int(n) # If n = 3.0 it passes the assertion so we make it 3
    fct = lib.aCorrUpTo if n is None else lib.aCorrUpToBit
    fct.argtype = (ndpointer(dtype=ctypes.c_uint8, shape=(len(x),)),
                   ctypes.c_uint64, 
                   ndpointer(dtype=ctypes.c_double, shape=(k,)),
                   ctypes.c_int64)

    r = zeros(k).tostring()
    fct(x, ctypes.c_uint64(len(x)), r, k) if n is None else fct(x, ctypes.c_uint64(len(x)), r, k, n)

    return fromstring(r)



