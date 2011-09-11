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
SIGMCMDST0=evtgen-u5stobb-0-evt_RL-DECAY_BDlnu.dec-20101108_1839_new_e000067.081110_1846
SIGMCMDST1=evtgen-u5stobb-1-evt_RL-DECAY_BDlnu.dec-20101108_1845_new_e000067.081110_1846
OUTPUTNAME=${AMOD}.${DATATYPE}.u5S.OR.e67rALL.bb-louvot

basf<< EOF >& ${OUTPUTDIR}/logs/${OUTPUTNAME}.log
path create main
module register fix_mdst
path add_module main fix_mdst

path create analysis
module register Adcab
path add_module analysis Adcab
module put_parameter Adcab JPsi_Veto_OS_Only\0

path add_condition main <=:0:KILL
path add_condition main >:0:analysis

initialize
histogram define ${OUTPUTDIR}/hbks/${OUTPUTNAME}.hbk
process_event ${SIGMCPATH}/${SIGMCMDST0}.mdst
process_event ${SIGMCPATH}/${SIGMCMDST1}.mdst
terminate
EOF
