#!/bin/sh
# File: sub-DATATYPE-uNS-RUNTYPE-eEXPNOrRUNNO.sh
# Author: M.P. Belhorn
# Date: 2010.10.17
# Description: This script submits a sample (MC or DATA) for processing. Check
#   output filename variables for description of data.

if [ $# -ne 3 ]
  then
    echo "Not enough arguments..."
    echo " Usage: $0 EXPNO EVENTTYPE STREAMNO"
    exit 1
  else
    bsub -q a ./basf.GNMC.u5S.OR.sh $1 $2 $3
fi

