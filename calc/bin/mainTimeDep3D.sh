#!/bin/bash

# Definition of an internal function.

function make3DOutput()
{

DIRNAME=$1
FILENAME=$2
SUFFIX=$3

OUTPUT=${DIRNAME}/${FILENAME}.dat

if [ "${SUFFIX}" != "" ]; then
    OUTPUT=${DIRNAME}/${FILENAME}-${SUFFIX}.dat
fi

if ( test -e ${OUTPUT} ); then
    rm -f ${OUTPUT}
fi

filecount=0
timeskip=5
spaceskip=5
timemax=500

for file in `ls ${DIRNAME}` ; do

    if [ ${filecount} -eq 0 ]; then

	if [ "${SUFFIX}" == "" ]; then
	    time=`echo ${file} | cut -d "-" -f 3 | sed -e "s/t=//g" | sed -e "s/fs.dat$//g"`
	else
	    suf=`echo ${file} | cut -d "-" -f 4 | sed -e "s/.dat$//g"`

	    if [ "${suf}" != "${SUFFIX}" ]; then
		continue
	    fi

	    time=`echo ${file} | cut -d "-" -f 3 | sed -e "s/t=//g" | sed -e "s/fs$//g"`
	fi
	

	if [ `echo "${time} <= ${timemax}" | bc` -eq 0 ]; then
	    break
	fi

	eval "awk 'NR%${spaceskip} == 0 {print ${time} \" \" \$1 \" \" \$2}' ${DIRNAME}/${file}" >> ${OUTPUT}
	echo "" >> ${OUTPUT}
    fi

    filecount=`expr ${filecount} + 1`
    if [ ${filecount} -eq ${timeskip} ]; then filecount=0; fi
done
}


# Main body.

MAINDIR=../dat/timeDep/L=3000/Lg=1000/Vg=50
FIELDDIR=${MAINDIR}/field
CONCDIR=${MAINDIR}/conc
MUDIR=${MAINDIR}/mu


make3DOutput ${FIELDDIR} "field" ""
make3DOutput ${MUDIR} "mu" "e"
make3DOutput ${MUDIR} "mu" "h"
make3DOutput ${CONCDIR} "conc" "e"
make3DOutput ${CONCDIR} "conc" "h"
