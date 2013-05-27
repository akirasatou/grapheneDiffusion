#!/bin/bash

dirname=../dat/timeDep/L=3000/Lg=1000/Vg=50/SS

eps=4
eps0=8.85e-12
d=1e-6
Vg=50
v=1e6
hbar=1.05e-34
e=1.6e-19


# Normalize potential

potFileIn=${dirname}/pot-SS.pos
potFileOut=${dirname}/pot-SS-normalized.pos

#pot0=`eval "python -c \"import math; print math.sqrt(2*math.pi*(${hbar}*${v})**2*(${eps}*${eps0}/${d})*${Vg}/(${e}**3))\""`
pot0=${Vg}

eval "awk 'NR>=8 && NF==1 && \$1!=\"\$EndView\" {print \$1/${pot0}; next} {print}' ${potFileIn}" > ${potFileOut}


# Normalize concentration

concFileIn=${dirname}/conc-SS-e.dat
concFileOut=${dirname}/conc-SS-e-normalized.dat

conc0=`eval "python -c \"import math; print 2*${eps}*${eps0}*${Vg}/(${d}*${e})\""`

eval "awk '{print \$1 \" \" \$2/${conc0}}' ${concFileIn}" > ${concFileOut}
