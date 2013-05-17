#!/bin/bash

# Definition of an internal function

outputMuVgDep()
{

Lg=$1
MAINDIR=../dat/VgLgDep/L=3000
OUTPUT=${MAINDIR}/mue-Lg=${Lg}.dat

MAINDIR=${MAINDIR}/Lg=${Lg}

if ( test -e ${OUTPUT} ); then
    rm -f ${OUTPUT}
fi

for Vg in 5 10 15 20 25 30 35 40 45 50; do
    MUDIR=${MAINDIR}/Vg=${Vg}/mu
    file=`ls ${MUDIR} | tail -2 | head -1`
    eval "awk 'NR==1 {print ${Vg} \" \" \$2}' ${MUDIR}/${file}" >> ${OUTPUT}
done
}


# Main body.

for Lg in 500 750 1000 1250; do
    outputMuVgDep ${Lg}
done
