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
    EXPNO=71
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
process_event ${MDSTPATH}/dilepskim-e000071r000027r000049-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r000050r000068-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r000069r000094-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r000095r000113-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r000114r000137-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r000138r000171-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r000172r000204-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r000205r000221-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r002001r002028-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r002029r002041-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r002042r002050-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r002051r002065-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r002076r002083-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r002086r002099-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r002100r002113-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r002114r002126-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r002127r002143-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r002144r002161-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r002162r002194-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r002195r002207-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r002211r002232-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000071r002233r002244-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
terminate
EOF
