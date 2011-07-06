#!/bin/sh
# Author: M.P. Belhorn
# Date: 2010.10.17
# Description: Submits all SGMC for an experiment.

for exp in 43 53 67 69 71; do
  for stream in {0..8}; do
    bsub -q a ./basf.SGMC.u5S.OR.sh ${exp} ${stream}
  done
done

