#!/bin/bash

# Definition of an internal function

outputMuVgWgLgDep()
{

Wg=$1
Lg=$2
MAINDIR=../dat/VgWgLgDep-grounded/
OUTPUT=${MAINDIR}/mue-Wg=${Wg}-Lg=${Lg}.dat

MAINDIR=${MAINDIR}/Wg=${Wg}/Lg=${Lg}

if ( test -e ${OUTPUT} ); then
    rm -f ${OUTPUT}
fi

for Vg in $(seq 10 10 400); do
    MUDIR=${MAINDIR}/Vg=${Vg}/mu

    time0=300
    echo ${Vg} `./findMinMuAt ${MUDIR} ${time0} e` >> ${OUTPUT}
done
}


# Main body.

WgTab="1200 1400 1600 1800"

for Wg in $(seq 1420 20 1580); do
    Lg=`expr ${Wg} \* 3 / 10`

    outputMuVgWgLgDep ${Wg} ${Lg}
done
