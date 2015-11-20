#!/bin/sh
ofile=rpie2_meanmax.txt
if [ -f ${ofile} ]; then
  rm -f ${ofile}
fi
ntrunc=39
#for ntrunc in 39 79 119 159 239 319 639 1279 2559 5119 10239; do
  awk -f meanmax.awk rpie2_T${ntrunc}.txt >> ${ofile}
#done
