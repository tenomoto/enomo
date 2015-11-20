#!/bin/sh
ntrunc=39
PNGDIR=.
#for ntrunc in 39 79 119 159 239 319 639 1279 2559 5119 10239; do
  ncl ntrunc=${ntrunc} plot.ncl > /dev/null
#  ncl ntrunc=${ntrunc} plot2.ncl > /dev/null
#done
#mv *.png ${PNGDIR}
