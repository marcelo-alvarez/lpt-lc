#!/bin/bash
source ~/.bashrc

if [ `uname` == "Darwin" ] ; then TMPDIR=/tmp;  fi
nproc=16

#boxsize=500
#nres=256
#testdata=delta_500Mpc_n256

boxsize=2000
boxsizeot=-`echo $boxsize | awk '{print $1/2}'`
echo $boxsizeot
nres=512
testdata=delta_z0_L2Gpch_N0512	
testdataz=zreion_L2Gpch_N0512

if [ ! -f $testdata ] ; then curl http://cita.utoronto.ca/~malvarez/data/$testdata > $testdata; fi

if [ ! -f $testdataz ] ; then curl http://cita.utoronto.ca/~malvarez/data/$testdataz > $testdataz; fi

# CMB Lensing
echo '------- CMB lensing test --------'
mapnum=8
testname='cmb'
#mpirun -np $nproc ../bin/lin2map -P param.lin2map -v -D $testdata -N 256 -B 500 -p 500 -x -250 -y -250 -z -250 -m $mapnum -o $testname

echo '------- 21cm test --------'
mapnum=8
testname='dtb'
#mpirun -np $nproc ../bin/lin2map -P param.lin2map -v -D $testdata -N 256 -B 500 -p 500 -x -250 -y -250 -z -250 -m $mapnum -o $testname

echo '------- 21cm test no evolution --------'
mapnum=8
testname='dtb_static'
#mpirun -np $nproc ../bin/lin2map -P param.lin2map -v -D $testdata -N 256 -B 500 -p 500 -x -250 -y -250 -z -250 -m $mapnum -o $testname -e 0 

echo '------- 21cm test no evolution high-z --------'
mapnum=8
nchunk=5
testname='dtb_static'
mpirun -n $nproc ../bin/lin2map -P param.lin2map -v -D $testdata -C $nchunk -N $nres -B $boxsize -p $boxsize -x $boxsizeot -y $boxsizeot -z $boxsizeot -m $mapnum -o $testname -e 0 

echo '------- 21cm test with evolution high-z --------'
mapnum=8
nchunk=1
testname='dtb_patchy'
#mpirun -n $nproc ../bin/lin2map -P param.lin2map -v -D $testdata -C $nchunk -N $nres -B $boxsize -p $boxsize -x $boxsizeot -y $boxsizeot -z $boxsizeot -m $mapnum -o $testname -R $testdataz



