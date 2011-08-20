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
DATAPATH=/bwf/g89home/g0mc/skim/skim5S/index/dilepskim/5S_onresonance/e000069/evtgen-nonbsbs/s00
DATAMDST=dilepskim-e000069r000012r000122-s00-evtgen-nonbsbs-5S_onresonance-b20090127_0910.index
OUTPUTNAME=diagnostic

basf<< EOF >& ${OUTPUTDIR}/logs/${OUTPUTNAME}.log
path create main
module register fix_mdst
path add_module main fix_mdst

path create analysis
module register Adcab
path add_module analysis Adcab
module put_parameter Adcab JPsi_Veto_OS_Only\0
module put_parameter Adcab Verbose_Log\0

path add_condition main <=:0:KILL
path add_condition main >:0:analysis

initialize
histogram define ${OUTPUTDIR}/hbks/${OUTPUTNAME}.hbk
process_event ${DATAPATH}/${DATAMDST}
terminate
EOF
