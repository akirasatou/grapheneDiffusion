#!/bin/bash

MAINDIR=../dat/timeDep/L=3000/Lg=1000/Vg=50
MUDIR=${MAINDIR}/mu

OUTPUT=${MUDIR}/mu-e.dat

if ( test -e ${OUTPUT} ); then
    rm -f ${OUTPUT}
fi

timemax=100000

for file in `ls ${MUDIR}` ; do
    suf=`echo ${file} | cut -d "-" -f 4 | sed -e "s/.dat$//g"`

    if [ "${suf}" != "e" ]; then
	continue
    fi

    time=`echo ${file} | cut -d "-" -f 3 | sed -e "s/t=//g" | sed -e "s/fs$//g"`

    if [ `echo "${time} <= ${timemax}" | bc` -eq 0 ]; then
	break
    fi

    eval "awk '\$1==\"3.75e-06\" {print ${time} \" \" \$2}' ${MUDIR}/${file}" >> ${OUTPUT}
done
