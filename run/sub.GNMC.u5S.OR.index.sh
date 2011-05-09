#!/bin/sh
# File: sub.GNMC.u5S.OR.index.sh
# Author: M.P. Belhorn
# Date: 2011.05.09

if [ $# -ne 4 ]
  then
    echo "Wrong number of arguments."
    echo " Usage: $0 exp type stream filename"
    echo " exp: 43 53 67 69 71"
    echo " type: bsbs charm nonbsbs uds"
    echo " stream: 0 1 2 3 4 5"
    exit 1
  else
    EXPNO=$1
    EVENTTYPE=$2
    STREAMNO=$3
    MDSTFILE=$4
    MDSTPATH=/bwf/g89home/g0mc/skim/skim5S/index/dilepskim/5S_onresonance/e0000${EXPNO}/evtgen-${EVENTTYPE}/s0${STREAMNO}
fi

if [ ! -f ${MDSTPATH}/${MDSTFILE}.index ]
  then
    echo " File does not exist. Aborting job sumbmission."
    exit 1
  else
    ANALYSISDIR=/bwf/g64home/mbelhorn/analysis/adcab
    OUTPUTDIR=${ANALYSISDIR}/output
    RUNSTART=`echo ${MDSTFILE} | cut -c 19-24`
    RUNEND=`echo ${MDSTFILE} | cut -c 26-31`
    OUTPUTNAME=Adcab.GNMC.OR.e${EXPNO}.rs${RUNSTART}.re${RUNEND}.${EVENTTYPE}.s0${STREAMNO}
    echo "Submitting job for ${OUTPUTNAME}"
fi

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
process_event ${MDSTPATH}/${MDSTFILE}.mdst
terminate
EOF
