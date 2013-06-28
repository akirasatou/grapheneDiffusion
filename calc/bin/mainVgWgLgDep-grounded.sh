#!/bin/bash

# Definition of an internal function

outputMuVgWgLgDep()
{

Wg=$1
Lg=$2
L=$3
MAINDIR=../dat/VgWgLgDep-grounded/
OUTPUT=${MAINDIR}/Wg=${Wg}/mue-L=${L}-Lg=${Lg}.dat

MAINDIR=${MAINDIR}/Wg=${Wg}/L=${L}/Lg=${Lg}

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

Wg="1500"
#LTab="4000 6000 8000"
#LLgTab="4 5 8 10 16 20"
LTab="8000"
LLgTab="40"

for L in ${LTab}; do
    for LLg in ${LLgTab}; do
	
	Lg=`expr ${L} / ${LLg}`
	
	outputMuVgWgLgDep ${Wg} ${Lg} ${L}
    done
done
