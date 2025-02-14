#!/bin/bash

# Definition of an internal function

outputMuVgDep()
{

Wg=$1
L=`expr ${Wg} \* 3`
MAINDIR=../dat/VgWgDep/L=${L}
OUTPUT=${MAINDIR}/mue-Wg=${Wg}.dat

MAINDIR=${MAINDIR}/Wg=${Wg}

if ( test -e ${OUTPUT} ); then
    rm -f ${OUTPUT}
fi

for Vg in $(seq 10 10 150); do
    MUDIR=${MAINDIR}/Vg=${Vg}/mu

    time0=300
    echo ${Vg} `./findMinMuAt ${MUDIR} ${time0} e` >> ${OUTPUT}
done
}


# Main body.

for Wg in 1500 2000; do
    outputMuVgDep ${Wg}
done
