#!/bin/sh
nproc=5
bin=./calc_checksum
bindir=`pwd`
ntrunc=39
#for ntrunc in 39 79 119 159 239 319 639 1279 2559 5119 10239; do
  nlat=`expr \( ${ntrunc} + 1 \) \* 3 / 2`
  mmax=`printf %0.5d ${ntrunc}`
  jmax=`printf %0.5d ${nlat}`
  if [ ! -d ../run/T${ntrunc} ]; then
    mkdir ../run/T${ntrunc}
  fi
  cd ../run/T${ntrunc}
  echo ${ntrunc} > ./ntrunc.txt
  mpirun -np ${nproc} ${bindir}/${bin}
  cd ..
#done
