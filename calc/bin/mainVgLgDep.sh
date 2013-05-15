#!/bin/bash

# Definition of an internal function

outputMuVgDep()
{

Lg=$1
MAINDIR=../dat/VgLgDep/L=3000/Lg=${Lg}/Vg=50
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
}


# Main body.

for Lg in "500 750 1000 1250"; do
    outputMuVgDep ${Lg}
done
