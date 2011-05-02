#!/bin/sh
# File: sub-DATATYPE-uNS-RUNTYPE-eEXPNOrRUNNO.sh
# Author: M.P. Belhorn
# Date: 2010.10.17
# Description: This script submits a sample (MC or DATA) for processing. Check
#   output filename variables for description of data.

AMOD=Adcab
ANALYSISDIR=/bwf/g64home/mbelhorn/analysis/adcab
OUTPUTDIR=${ANALYSISDIR}/output
DATATYPE=SGMC
SIGMCPATH=/bwf/g68home/louvot/belle/fs/mcprod/mdstfile
SIGMCMDST='evtgen-u5stobb-0-evt_RL-DECAY_BDlnu.dec-20101108_1839_new_e000067.081110_1846'
OUTPUTNAME=${AMOD}.${DATATYPE}-diagnostic

basf << EOF >& ${OUTPUTDIR}/logs/test.log
path create main
path add_module main fix_mdst
path create analysis
path add_module analysis Adcab
path add_condition main <=:0:KILL
path add_condition main >:0:analysis
initialize
histogram define ${OUTPUTDIR}/hbks/${OUTPUTNAME}.hbk
process_event ${SIGMCPATH}/${SIGMCMDST}.mdst
output close
terminate
EOF
