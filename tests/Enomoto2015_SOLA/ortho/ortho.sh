#!/bin/sh
#ntrunc=1279
ntrunc=39
bin=ortho
DATADIR=.
#for ntrunc in 39 79 119 159 239 319 639 1279 2559 5119 10239; do
  echo ${ntrunc}
  ./${bin} ${ntrunc} > ${DATADIR}/ortho_T${ntrunc}.txt
#done
