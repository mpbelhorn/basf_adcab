#!/bin/bash
# File: sub.GNMC.u5S.OR.index.sh
# Author: M.P. Belhorn
# Date: 2011.05.09

if [ $# -ne 3 ]
  then
    echo "Wrong number of arguments."
    echo " Usage: $0 exp type stream"
    echo " exp: 43 53 67 69 71"
    echo " type: bsbs charm nonbsbs uds"
    echo " stream: 0 1 2 3 4 5"
    exit 1
  else
    EXPNO=$1
    EVENTTYPE=$2
    STREAMNO=$3
fi

shopt -s nullglob
bsub -q index ./basf.GNMC.u5S.OR.index.e${EXPNO}.sh ${EVENTTYPE} ${STREAMNO}
