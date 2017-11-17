#!/bin/bash
source ~/.bashrc

if [ `uname` == "Darwin" ] ; then TMPDIR=/tmp;  fi
nproc=4
testdata=delta_500Mpc_n256
if [ ! -f $testdata ] ; then curl http://cita.utoronto.ca/~malvarez/data/$testdata > $testdata; fi

# CMB Lensing
echo '------- CMB lensing test --------'
mapnum=8
testname='cmb'
#mpirun -np $nproc ../bin/lin2map -P param.lin2map -v -D $testdata -N 256 -B 500 -p 500 -x -250 -y -250 -z -250 -m $mapnum -o $testname

echo '------- 21cm test --------'
mapnum=8
testname='dtb'
mpirun -np $nproc ../bin/lin2map -P param.lin2map -v -D $testdata -N 256 -B 500 -p 500 -x -250 -y -250 -z -250 -m $mapnum -o $testname

echo '------- 21cm test --------'
mapnum=8
testname='dtb_static'
mpirun -np $nproc ../bin/lin2map -P param.lin2map -v -D $testdata -N 256 -B 500 -p 500 -x -250 -y -250 -z -250 -m $mapnum -o $testname -e 0 

