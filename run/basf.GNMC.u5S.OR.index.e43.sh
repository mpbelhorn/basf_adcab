#!/bin/sh
# File: sub.GNMC.u5S.OR.index.sh
# Author: M.P. Belhorn
# Date: 2011.05.09

if [ $# -ne 2 ]
  then
    echo "Wrong number of arguments."
    echo " Usage: $0 type stream"
    echo " type: bsbs charm nonbsbs uds"
    echo " stream: 0 1 2 3 4 5"
    exit 1
  else
    EXPNO=43
    EVENTTYPE=$1
    STREAMNO=$2
    MDSTPATH=/bwf/g89home/g0mc/skim/skim5S/index/dilepskim/5S_onresonance/e0000${EXPNO}/evtgen-${EVENTTYPE}/s0${STREAMNO}
fi

ANALYSISDIR=/bwf/g64home/mbelhorn/analysis/adcab
OUTPUTDIR=${ANALYSISDIR}/output
OUTPUTNAME=Adcab.GNMC.OR.s0${STREAMNO}.e${EXPNO}.${EVENTTYPE}
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
process_event ${MDSTPATH}/dilepskim-e000043r001013r001026-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000043r001028r001034-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
terminate
EOF
