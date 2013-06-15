#!/bin/bash

# Definition of an internal function

outputMuVgDep()
{

Lg=$1
MAINDIR=../dat/VgDep-grounded/
OUTPUT=${MAINDIR}/mue-Lg=${Lg}.dat

MAINDIR=${MAINDIR}/Lg=${Lg}

if ( test -e ${OUTPUT} ); then
    rm -f ${OUTPUT}
fi

for Vg in $(seq 50 50 400); do
    MUDIR=${MAINDIR}/Vg=${Vg}/mu

    time0=300
    echo ${Vg} `./findMinMuAt ${MUDIR} ${time0} e` >> ${OUTPUT}
done
}


# Main body.


for Lg in 300 450 600; do
    outputMuVgDep ${Lg}
done
