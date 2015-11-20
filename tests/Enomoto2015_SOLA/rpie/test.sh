#!/bin/sh
nproc=5
ntrunc=39
bin=rpie2
#for ntrunc in 39 79 119 159 239 319 639 1279 2559 5119; do
  echo ${ntrunc}
  mpirun -np ${nproc} ./${bin} ${ntrunc}
  mv ${bin}.txt ${bin}_T${ntrunc}.txt
  rm -f ${bin}.txt
#done
