#!/bin/python
# -*- coding: utf-8 -*-

import os, sys
from autocorrelation import aCorrUpTo
from numpy import empty, savez, load, fromfile, uint8
from pylab import *

# Tips for the proper approach:
# os.environ["PATH"] += os.pathsep + path
# os.path.dirname(os.path.abspath(__file__))+'123.dll'

def NOP(x):
    return x

def loadfile(fname, format='numpy'):
    if format == 'string':
        with open(fname,'rb') as f:
            x = f.read()
        return x
    if format == 'numpy':
        return fromfile(fname, dtype=uint8)

def get_corrs_from_file(fname, k):
    """
    Get byte-corrs and every bit-corrs from filename for the kth first lags.
    
    Also write results to disk and use them if they already exist.
    """
    basename = os.path.splitext(os.path.basename(fname))[0]
    par = dict(fname=fname, 
            npzname='{}_autocorrelations_{}.npz'.format(basename, k),
            k=k)
    if not os.path.isfile(par['npzname']):
        x = loadfile(par['fname'])
        corrs = empty((9,k)) # Initializing array
        for i,n in enumerate([None]+range(8)[::-1]):   #Bytes corr, then bit corrs from MSB to LSB.
            corrs[i] = aCorrUpTo(x, k, n)
        savez(par['npzname'], corrs)

    else:
        corrs = load(par['npzname'])['arr_0']

    return corrs


def plot_corrs(corrs, xmin=None, xmax=None, ymin=None, ymax=None, save=False, op=NOP, log=False, frame=False, raster=False, dpi=100, loc=0):
    if log:
        oldop = op
        op = lambda x: abs(oldop(x))
    labels = ('Whole Byte',
              '7th bit (MSB)',
              '6th bit',
              '5th bit',
              '4th bit',
              '3rd bit',
              '2nd bit',
              '1st bit',
              '0th bit (LSB)')
    f = plt.figure()
    ax = f.add_subplot(111)
    for c,l in zip(corrs, labels):
        ax.plot(op(c), label=l, rasterized = raster)
    ax.legend(frameon=frame, loc = loc)
    
    if log:
        ax.set_yscale('log')
    
    ax.set_xlim(ax.get_xlim()[0] if xmin is None else xmin, ax.get_xlim()[1] if xmax is None else xmax)
    ax.set_ylim(ax.get_ylim()[0] if ymin is None else ymin, ax.get_ylim()[1] if ymax is None else ymax)
    
        
    ax.set_xlabel('Lag (byte)')
    ax.set_ylabel('Autocorrelation {}'.format('(log)' if log else ''))    

    draw_if_interactive()

    if save is not False:
        fname = save+'-Correlations'
        if log:
            fname += '_log'
        if xmax is not None:
            fname += '_xmax{}'.format(xmax)
        fname += '_0.pdf'
        f.savefig(fname, bbox_inches='tight', dpi = dpi)

    return f,ax
    
def plot_spectrum(Y, SampleRate, log=False, save=False, raster=False, xmin=None, xmax=None, ymin=None, ymax=None, op=lambda x:20*log(abs(x)), ycomment=False):
    Sp = fftshift(fft(append(Y[-1:0:-1],Y))) # Power Spectral Density (Wiener–Khintchine). Hypothèse autocorelation symmétrique.
    X = fftfreq(len(Sp), 1/SampleRate)
    X = append(X[len(X)/2 + 1 if len(Sp)%2 else 0:], X[:len(X)/2 + 1 if len(Sp)%2 else 0])/1e6   # MHz

    f = plt.figure()
    ax = f.add_subplot(111)
    ax.plot(X, op(Sp), rasterized=raster) # Default, op : abs(real(Sp))
    
    if log:
        ax.set_yscale('log')
    
    ax.set_xlim(ax.get_xlim()[0] if xmin is None else xmin, ax.get_xlim()[1] if xmax is None else xmax)
    ax.set_ylim(ax.get_ylim()[0] if ymin is None else ymin, ax.get_ylim()[1] if ymax is None else ymax)
    
    ax.set_xlabel('Frequency (MHz)')
    ax.set_ylabel('Spectrum {} {}'.format('(log)' if log else '', ycomment if ycomment else ''))

    draw_if_interactive()

    if save:
        fname = save+'-Spectrum'
        if log:
            fname += '_log'
        if xmax is not None:
            fname += '_xmax{}'.format(xmax)
        fname += '_0.pdf'
        f.savefig(fname, bbox_inches='tight', dpi = 100)

    return f,ax
    
    
    

            

    
