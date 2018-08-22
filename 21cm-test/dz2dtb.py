import numpy as np
import sys
import os.path
import matplotlib.pyplot as plt
from evolution import *
    
def getSize(filename):
    fileobject=open(filename)
    fileobject.seek(0,2) # move the cursor to the end of the file
    size = fileobject.tell()
    fileobject.close()
    return size

def int_cube_root(n):
    c = int(n**(1/3.))
    if c**3 == n:
        return c
    elif (c+1)**3 == n:
        return c+1

if len(sys.argv) != 5:
    sys.exit("usage: python dz2dtb.py <densityfile> <zreionfile> <outputfile> <redshift>")

dfile = sys.argv[1]
zfile = sys.argv[2]
ofile = sys.argv[3]
z     = float(sys.argv[4])

nfloat = int(getSize(dfile)/4.)
n      = int_cube_root(nfloat)

dlin = np.fromfile(dfile,dtype=np.float32,count=n**3)
zrei = np.fromfile(zfile,dtype=np.float32,count=n**3)

xHI = 0 * zrei
xHI[z>zrei] = 1

dtb = dtbofz(z) * xHI * (1+growth(z)*dlin)

dtbslice = np.reshape(dtb,(n,n,n))[:,:,0]
plt.imshow(dtbslice)
plt.savefig('dtbslice.pdf')

if os.path.isfile(ofile):
    sys.exit("refusing to overwrite "+ofile+", exiting...")
else:
    dtb.astype(np.float32).tofile(ofile)

