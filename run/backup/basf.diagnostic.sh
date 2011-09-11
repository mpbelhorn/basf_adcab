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
DATAPATH=/bwf/g64home/mbelhorn/analysis/MC/signalmc-u5s-mix/mdst
DATAMDST=fs19299.s0.e53.mBdBd.n0019123.mdst
OUTPUTNAME=diagnostic

basf<< EOF >& ${OUTPUTDIR}/logs/${OUTPUTNAME}.log
path create main
module register fix_mdst
path add_module main fix_mdst

path create analysis
module register Adcab
path add_module analysis Adcab
module put_parameter Adcab JPsi_Veto_OS_Only\0
module put_parameter Adcab Verbose_Log\2

path add_condition main <=:0:KILL
path add_condition main >:0:analysis

initialize
histogram define ${OUTPUTDIR}/hbks/${OUTPUTNAME}.hbk
process_event ${DATAPATH}/${DATAMDST}
terminate
EOF
