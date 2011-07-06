#!/bin/sh
# File: sub.GNMC.u5S.OR.index.sh
# Author: M.P. Belhorn
# Date: 2011.05.09

if [ $# -ne 3 ]
  then
    echo "Wrong number of arguments."
    echo " Usage: $0 exp stream"
    echo " exp: 43 53 67 69 71"
    echo " stream: 0-8"
    exit 1
  else
    EXPNO=$1
    STREAMNO=$2
    MDSTPATH=/bwf/g64home/mbelhorn/analysis/MC/signalmc-u5s-mix/mdst
fi

ANALYSISDIR=/bwf/g64home/mbelhorn/analysis/adcab
OUTPUTDIR=${ANALYSISDIR}/output
OUTPUTNAME=Adcab.SGMC.OR.s0${STREAMNO}.e${EXPNO}
echo "Submitting job for ${OUTPUTNAME}"

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
process_event ${MDSTPATH}/Y5S_to_BsBs_to_semileptonic-e${EXPNO}-f${STREAMNO}.mdst
process_event ${MDSTPATH}/Y5S_to_BdBdX_to_semileptonic-e${EXPNO}-f${STREAMNO}.mdst
process_event ${MDSTPATH}/Y5S_to_BdBuX_to_semileptonic-e${EXPNO}-f${STREAMNO}.mdst
process_event ${MDSTPATH}/Y5S_to_BuBuX_to_semileptonic-e${EXPNO}-f${STREAMNO}.mdst
terminate
EOF
