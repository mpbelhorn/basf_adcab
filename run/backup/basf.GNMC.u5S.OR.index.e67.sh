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
    EXPNO=67
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
process_event ${MDSTPATH}/dilepskim-e000067r000001r000022-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000098r000115-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000117r000148-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000149r000167-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000168r000190-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000191r000224-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000225r000251-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000252r000272-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000273r000299-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000300r000321-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000322r000340-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000341r000371-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000372r000391-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000394r000408-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000409r000419-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000420r000431-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000432r000441-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000456r000483-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000484r000488-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000493r000507-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000508r000571-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000572r000591-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000592r000602-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000603r000663-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000664r000678-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000679r000686-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000067r000688r000696-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
terminate
EOF
