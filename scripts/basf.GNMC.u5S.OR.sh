#!/bin/sh
# File: basf.GNMC.u5S.OR.sh
# Author: M.P. Belhorn
# Date: 2011.02.16
# Description: This script submits a GNMC sample for processing. It is
#   to be called by a sub.GNMC.sh script.

if [ $# -ne 3 ]
  then
    echo "Not enough arguments..."
    echo " Usage: $0 EXPNO EVENTTYPE STREAMNO"
    exit 1
  else
    EXPNO=$1
    EVENTTYPE=$2
    STREAMNO=$3
fi
  

AMOD=Adcab
ADIR=/bwf/g64home/mbelhorn/analysis/adcab

# The following variables set the output filenames of the run.
# DATATYPE is a four letter identifier for the three data types:
#   (G)e(N)eric(MC) (S)i(G)nal(MC) (DATA).
# REGIME is a single number N to give the Y(NS) energy regime.
#   N = 1 2 3 4 5
# RUNTYPE is a two letter identifier for the processed MDST files.
#   (SC)an (O)n(R)esonance (CO)ntinuum
# EXPNO is the experiment number.
#   Inlude leading 0 if necessary.
# RUNNO is the starting run number.
#   0 for MC
DATATYPE=GNMC
REGIME=5
RUNTYPE=OR
RUNNO=ALL

# The following variables should be fixed.
OUTPUTNAME=${AMOD}.${DATATYPE}.u${REGIME}S.${RUNTYPE}.e${EXPNO}r${RUNNO}.${EVENTTYPE}.${STREAMNO}
HBKFILE=${OUTPUTNAME}.hbk

basf << EOF >& ${ADIR}/output/logs/${OUTPUTNAME}.log
path create main
path add_module main fix_mdst
path create analysis
path add_module analysis ${AMOD}
path add_condition main <=:0:KILL
path add_condition main >:0:analysis
initialize
histogram define ${ADIR}/output/hbks/${HBKFILE}
process_url http://bweb3/montecarlo.php?ex=${EXPNO}&rs=1&re=9999&ty=evtgen-${EVENTTYPE}&dt=5S_onresonance&bl=caseB&dv=zfserv&st=${STREAMNO}
output close
terminate
EOF
