#!/usr/bin/env python

import scipy
from scipy import fftpack
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import radialProfile as rp

if len(sys.argv) != 4: 
    print 'Usage: ./grid2pk.py <gridfilename> <boxsize> <outbasename>'
    print '       for example the command'
    print '         ./grid2pk.py grid.bin 100 test'
    print '       will read in N^3 32-bit floats from grid.bin'
    print '       (with N determined from the filesize) and output the'
    print '       power spectrum into the file outbasename.npz'
    sys.exit()

# Command line arguments
datafile = sys.argv[1]
boxsize  = sys.argv[2]
outfile  = sys.argv[3]

boxsize  = np.float(boxsize) 

# Determine dimension of grid data assuming it's a cube of float32s
fsize = os.stat(datafile).st_size
ifile = open(sys.argv[1])

dsize = fsize/4
n3    = int(round(dsize**(1./3.)))
print 'N = ',n3

dk  = 2*np.pi/boxsize
d3k = boxsize**(-3)

# Read in grid data
datareal = np.fromfile(ifile,count=dsize,dtype=np.float32)
datareal = np.reshape(datareal,[n3,n3,n3])

# FFT
datacplx  = fftpack.fftn(datareal)
datacplx  = fftpack.fftshift(datacplx)
datacplx /= dsize

# 3D power spectrum
p3d  = np.abs(datacplx)**2 
p3d /= d3k

# 1D spherically averaged power spectrum
p1d = rp.sphericalAverage(p3d)
k1d = np.arange(len(p1d))+1.5 # Add 1 here for the Nyquist plane
k1d = k1d * dk

# Save
np.savez(outfile+'.npz',k=k1d,pofk=p1d)
