.. cal-sky documentation master file, created by
   sphinx-quickstart on Wed Apr 11 18:59:35 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2

**Astropnysical Models**
------------------

Cosmic infrared background
^^^^^^^^^

Kinetic Sunyaev-Zel'dovich and Optical Depth
^^^^^^^^^

Patchy Reionization
''''''''''

Late time
'''''''''' 

21-cm intensity maps
^^^^^^^^^

Patchy Reionization
''''''''''''''

The 21-cm signal from patchy reionization is obtained, in the $T_s >> T_{cmb}$ limit, by using Lagrangian perturbation 
theory to displace mass cells to their Eulerian positions:
$${\mathbf x}_c = {\mathbf q}_c + {\mathbf s}^{(1)}_c({\mathbf q}_c)$$. 
Each mass element originate from a comoving volume of $a_{\rm latt}^3$. 
For a sufficiently small frequency bin labeled with index $i$ of size $\delta{\nu}_i$ 
in a pixel $p$ subtending a solid angle $\Delta\Omega_p$, the mean differential 
brightness temperature in the corresponding frequency voxel is given by summing 
all the values of the comoving voxels that overlap the frequency voxel, weighted by the 
ratio of comoving voxel volume to frequency voxel, 
$$\delta{T}_{b,ip}=\sum_c f_c \frac{V_c}{V_{ip}}$$, where the frequency voxel volume is 
$$V_{ip}=\Delta\Omega_p\chi_i^2\frac{d\chi}{d\nu}_i\delta{\nu}_i$$. 
The mean over a finite frequency bin, $\Delta{\nu}_j$, corresponds to a sum over infinitesimal ones, 
$$\delta{T}_{b,jp}=\frac{1}{\Delta{\nu}_j}\sum_i \Delta{\nu}_i\delta{T}_{b,ip}$$. Combining all 
these expressions together, we obtain the weighting scheme to average over comoving voxel 
when binning into finite size frequency voxels $jp$: 
$$\delta{T}_{b,jp}=\frac{V_c}{\Delta\Omega_p\Delta{\nu}_j}\sum_c \delta{T}_{b,c}\left(\frac{d\nu}{d\chi}$$


HI in galaxies
''''''''''''''''

Inhomogeneous heating
'''''''''''''''''''''''

Lyman-alpha pumping
'''''''''''''''''''''

Lensing
^^^^^^^






