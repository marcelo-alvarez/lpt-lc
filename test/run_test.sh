#!/bin/bash

if [ `uname` == "Darwin" ] ; then TMPDIR=/tmp;  fi

nproc=2

boxsize=2000
boxsizeot=-`echo $boxsize | awk '{print $1/2}'`
nres=512
testdata=delta_z0_L2Gpch_N0512	
testdataz=zreion_L2Gpch_N0512

if [ ! -f $testdata ] ; then curl http://cita.utoronto.ca/~malvarez/data/$testdata > $testdata; fi

if [ ! -f $testdataz ] ; then curl http://cita.utoronto.ca/~malvarez/data/$testdataz > $testdataz; fi

echo '------- 21cm test --------'
mapnum=8
nchunk=1
testname='21cm'
mpirun -n $nproc ../bin/lin2map -P param.lin2map -v -D $testdata -C $nchunk -N $nres -B $boxsize -p $boxsize -x $boxsizeot -y $boxsizeot -z $boxsizeot -m $mapnum -o $testname -R $testdataz

echo '------- cmb test with evolution high-z --------'
mapnum=6
testname='cmb'
mpirun -n $nproc ../bin/lin2map -P param.lin2map -v -D $testdata -N $nres -B $boxsize -p $boxsize -x $boxsizeot -y $boxsizeot -z $boxsizeot -m $mapnum -o $testname -R $testdataz




