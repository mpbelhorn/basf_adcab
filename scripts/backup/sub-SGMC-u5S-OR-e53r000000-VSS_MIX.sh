#!/bin/sh
# File: sub-DATATYPE-uNS-RUNTYPE-eEXPNOrRUNNO.sh
# Author: M.P. Belhorn
# Date: 2010.10.17
# Description: This script submits a sample (MC or DATA) for processing. Check
#   output filename variables for description of data.

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
DATATYPE=SGMC
REGIME=5
RUNTYPE=OR
EXPNO=53
RUNNO=000000
#     ###### <- Match six digits in RUNNO for consistancy.

# The following variables are only used for SIGMC 
SIGMCPATH=/bwf/g64home/mbelhorn/analysis/MC/signalmc-u5s/mdst
SIGMCMDST=Bsdilepton.VSS_MIX-e53-n10000
SIGMCMODEL=VSS_MIX

# The following variables should be fixed.
OUTPUTNAME=${AMOD}.${DATATYPE}-u${REGIME}S-${RUNTYPE}-e${EXPNO}r${RUNNO}
HBKFILE=${OUTPUTNAME}.hbk

basf << EOF >& ${ADIR}/output/${OUTPUTNAME}.log
path create main
path add_module main fix_mdst
path create analysis
path add_module analysis ${AMOD}
path add_condition main <=:0:KILL
path add_condition main >:0:analysis
initialize
histogram define ${ADIR}/output/${HBKFILE}
process_event ${SIGMCPATH}/${SIGMCMDST}-f0.mdst
process_event ${SIGMCPATH}/${SIGMCMDST}-f1.mdst
process_event ${SIGMCPATH}/${SIGMCMDST}-f2.mdst
process_event ${SIGMCPATH}/${SIGMCMDST}-f3.mdst
process_event ${SIGMCPATH}/${SIGMCMDST}-f4.mdst
output close
terminate
EOF
