.. cal-sky documentation master file, created by
   sphinx-quickstart on Wed Apr 11 18:59:35 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2

**Quickstart**
-------------

*Download*
^^^^^^^^^^

The easiest way to obtain the code is by cloning from github:
    ``git clone git@github.com:marcelo-alvarez/cal-sky.git``

It's also possible to get download the zip version here_

*Install*
^^^^^^^^^

#. create a new copy of the file Make.mach.generic and rename to your machine:
    ``cp Make.mach.generic Make.mach.\<your-machine\>``
    
#. edit the file Make.mach.\<your-machine\> to reflect your libraries and C/C++ compiler flags (you need to have mpi, gsl, fftw-2.1.5, and cfitsio installed)
   
#. run configure (which just concatenates two files to make the Makefile):
    ``./configure \<your-machine\>``
    
#. compile:
    ``make``

*Test*
^^^^^^^

#. go to the test directory and run the test:
   ::
     ``cd test``
     ``./run_test.sh``




 .. _here: https://github.com/marcelo-alvarez/cal-sky/archive/master.zip
