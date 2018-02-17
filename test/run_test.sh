#!/bin/bash

if [ $(uname) == 'Darwin' ] ; then export TMPDIR='/tmp'; fi
echo $TMPDIR

source ~/.bashrc

boxsize=2000
boxsizeot=-`echo $boxsize | awk '{print $1/2}'`
nres=32
testdata=delta_z0_L2Gpch_N0512	
testdataz=zreion_L2Gpch_N0512

if [ ! -f $testdata ] ; then wget http://cita.utoronto.ca/~malvarez/data/$testdata ; fi

if [ ! -f $testdataz ] ; then wget http://cita.utoronto.ca/~malvarez/data/$testdataz ; fi
echo '------- 21cm test (no fits serial) --------'
mapnum=8
nchunk=1
binary_only=1
testname='21cm'
nproc=1
mpirun -n $nproc ../bin/lin2map -P param.lin2map -v -D $testdata -C $nchunk -N $nres -B $boxsize -p $boxsize -x $boxsizeot -y $boxsizeot -z $boxsizeot -m $mapnum -o $testname -R $testdataz -b $binary_only

if [ ! -f $testdataz ] ; then wget http://cita.utoronto.ca/~malvarez/data/$testdataz ; fi
echo '------- 21cm test (no fits parallel) --------'
mapnum=8
nchunk=1
binary_only=1
testname='21cm'
nproc=1
#mpirun -n $nproc ../bin/lin2map -P param.lin2map -v -D $testdata -C $nchunk -N $nres -B $boxsize -p $boxsize -x $boxsizeot -y $boxsizeot -z $boxsizeot -m $mapnum -o $testname -R $testdataz -b $binary_only

echo '------- 21cm test (with fits serial) --------'
mapnum=8
nchunk=1
binary_only=0
testname='21cm'
nproc=1
#mpirun -n $nproc ../bin/lin2map -P param.lin2map -v -D $testdata -C $nchunk -N $nres -B $boxsize -p $boxsize -x $boxsizeot -y $boxsizeot -z $boxsizeot -m $mapnum -o $testname -R $testdataz -b $binary_only

echo '------- 21cm test (with fits parallel) --------'
mapnum=8
nchunk=1
binary_only=0
testname='21cm'
nproc=2
#mpirun -n $nproc ../bin/lin2map -P param.lin2map -v -D $testdata -C $nchunk -N $nres -B $boxsize -p $boxsize -x $boxsizeot -y $boxsizeot -z $boxsizeot -m $mapnum -o $testname -R $testdataz -b $binary_only

#echo '------- cmb test with evolution high-z --------'
#mapnum=6
#testname='cmb'
#mpirun -n $nproc ../bin/lin2map -P param.lin2map -v -D $testdata -N $nres -B $boxsize -p $boxsize -x $boxsizeot -y $boxsizeot -z $boxsizeot -m $mapnum -o $testname -R $testdataz




