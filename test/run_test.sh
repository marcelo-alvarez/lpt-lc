#!/bin/bash

if [ $(uname) == 'Darwin' ] ; then export TMPDIR='/tmp'; fi
echo $TMPDIR

source ~/.bashrc

boxsize=2000
boxsizeot=-`echo $boxsize | awk '{print $1/2}'`
nres=512
testdata=delta_z0_L2Gpch_N0512	
testdataz=zreion_L2Gpch_N0512

if [ ! -f $testdata ]  ; then wget http://cita.utoronto.ca/~malvarez/data/$testdata  ; fi
if [ ! -f $testdataz ] ; then wget http://cita.utoronto.ca/~malvarez/data/$testdataz ; fi

echo '------- 21cm test high-z --------'
mapnum=8
nchunk=1
binary_only=0
testname='21cm_high-z'
nproc=4
#mpirun -n $nproc ../bin/lin2map -P param.high-z -v -D $testdata -C $nchunk -N $nres -B $boxsize -p $boxsize -x $boxsizeot -y $boxsizeot -z $boxsizeot -m $mapnum -o $testname -R $testdataz -b $binary_only
cp weights.dat weights_high-z.dat

echo '------- 21cm test low-z --------'
mapnum=8
nchunk=1
binary_only=0
testname='21cm_low-z'
nproc=4
mpirun -n $nproc ../bin/lin2map -P param.low-z -v -D $testdata -C $nchunk -N $nres -B $boxsize -p $boxsize -x $boxsizeot -y $boxsizeot -z $boxsizeot -m $mapnum -o $testname -b $binary_only
cp weights.dat weights_low-z.dat
