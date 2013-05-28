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

for Vg in $(seq 5 5 150); do
    MUDIR=${MAINDIR}/Vg=${Vg}/mu

    time0=300
    echo ${Vg} `./findMinMuAt ${MUDIR} ${time0} e` >> ${OUTPUT}
done
}


# Main body.

for Lg in 500 750 1000 1250; do
    outputMuVgDep ${Lg}
done
