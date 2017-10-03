#!/bin/bash

testdata=delta_500Mpc_n256
if [ ! -f $testdata ] ; then curl http://cita.utoronto.ca/~malvarez/data/$testdata > $testdata; fi 
mpirun -np 4 ../bin/lin2map -P param.lin2map -v -D $testdata -N 256 -B 500 -p 500 -x 250 -y 250 -z 250 -m 1

