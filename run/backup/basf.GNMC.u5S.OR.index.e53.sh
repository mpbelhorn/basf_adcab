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
    EXPNO=53
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
process_event ${MDSTPATH}/dilepskim-e000053r000001r000022-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000023r000033-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000034r000041-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000045r000061-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000062r000068-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000069r000081-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000082r000097-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000098r000110-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000113r000127-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000128r000141-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000142r000161-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000162r000169-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000170r000187-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000190r000198-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000199r000206-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000208r000212-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000213r000220-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000221r000235-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000053r000241r000272-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
terminate
EOF
