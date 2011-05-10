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
    EXPNO=69
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
process_event ${MDSTPATH}/dilepskim-e000069r000012r000122-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000124r000189-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000191r000234-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000236r000274-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000275r000305-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000308r000333-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000334r000356-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000358r000385-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000387r000405-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000407r000424-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000425r000440-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000442r000466-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000467r000494-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000496r000507-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000509r000531-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000532r000541-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000542r000554-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000556r000599-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000600r000621-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000623r000636-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000637r000649-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000652r000683-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000684r000697-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000698r000710-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000711r000733-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000734r000756-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000757r000790-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000791r000819-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000892r000909-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000913r000933-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000934r000964-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000965r000980-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000981r000988-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r000989r001004-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r001005r001020-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r001022r001039-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r001040r001083-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r001085r001109-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r001116r001122-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r001125r001130-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r001131r001156-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r001157r001191-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r001193r001211-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r001212r001221-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r001222r001283-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
process_event ${MDSTPATH}/dilepskim-e000069r001289r001309-s0${STREAMNO}-evtgen-${EVENTTYPE}-5S_onresonance-b20090127_0910.index
terminate
EOF
