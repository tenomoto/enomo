#!/bin/sh
prefix="checksum_T"
PNGDIR=.
ntrunc=39
#for ntrunc in 39 79 119 159 239 319 639 1279 2559 5119 10239; do
  echo ${ntrunc}
  ncl ntrunc=${ntrunc} checksum.ncl
  mv ${prefix}${ntrunc}.000001.png ${PNGDIR}/${prefix}${ntrunc}_alf.png
  mv ${prefix}${ntrunc}.000002.png ${PNGDIR}/${prefix}${ntrunc}_alfx.png
  mv ${prefix}${ntrunc}.000003.png ${PNGDIR}/${prefix}${ntrunc}_alff.png
#done
