#!/bin/python
# -*- coding: utf-8 -*-

import ctypes
import os
import subprocess
from numpy.ctypeslib import ndpointer  

if not os.path.isfile("./autocorrelation.so"):
    print("autocorrelation.so is missing. To compile on unix:\ngcc -O3 -fPIC -shared autocorrelation.c -o autocorrelation.so -Wall -lmpfr -lgmp -fopenmp -DUNIX\n")

lib = ctypes.cdll["./autocorrelation.so"]

def aCorrUpTo(x, k, n=None):
    """
    First k lags of x autocorrelation with optional nth bit bitmask.

    If n is None, it calculates it on the whole bytes. 
    Otherwise it calculates it on the nth bit of each byte.
    """
    assert (n is None or (0 <= n <= 7)) and k>0, \
           "Invalid n or k. Condition is: (n is None or (0 <= n <= 7)) and k>0"
    fct = lib.aCorrUpTo if n is None else lib.aCorrUpToBit
    fct.restype = ndpointer(dtype=ctypes.c_double, shape=(k,))
    return fct(x, len(x), k) if n is None else fct(x, len(x), k, n)

