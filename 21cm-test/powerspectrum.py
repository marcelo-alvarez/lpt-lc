import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from evolution import *

boxsize  = '2e3'
redshift = '8'
z        =  8
code     = 'L2Gpch_N0512'

dtbfile    =   'dtb_z'  + redshift + '_' + code
deltafile  = 'delta_z'  + '0'      + '_' + code
zreionfile =  'zreion'  +            '_' + code

dtbpkfile   = dtbfile   + '_pk'
deltapkfile = deltafile + '_pk'

# make 21-cm cube at z=8 in mK
#os.system('rm -f' + ' ' + dtbfile)
#os.system('python dz2dtb.py' + ' ' + deltafile + ' ' + zreionfile + ' ' + dtbfile + ' ' + redshift)

# get power spectrum of linear density contrast and 21-cm cubes
#os.system('python grid2pk.py' + ' ' + deltafile + ' ' + boxsize + ' ' + deltapkfile)
#os.system('python grid2pk.py' + ' ' +   dtbfile + ' ' + boxsize + ' ' +   dtbpkfile)

data_lin = np.load(deltapkfile + '.npz')
data_dtb = np.load(dtbpkfile   + '.npz')

dtb      = np.fromfile(dtbfile,dtype=np.float32)
lin      = np.fromfile(deltafile,dtype=np.float32)

print dtb.var()**0.5,lin.var()**0.5

k = data_lin['k']

pklin = data_lin['pofk']
pkdtb = data_dtb['pofk']

D  = growth(z) / growth(0.0001)
T0 = dtbofz(z)

pkdtb_lin = pklin * D**2 * T0**2 # this is what dtb power spectrum would be if neutral hydrogen traced matter

bias = np.sqrt(pkdtb / pkdtb_lin)

plt.semilogx(k,bias)
plt.show()
