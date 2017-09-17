#!/bin/sh

max=10
for i in `seq 1 $max`
do
    #source step1_createConfigsPion_directional.csh 1000 $i 10
    #source submitCondorPion_directional.csh 1000 $i 10
    source step1_createConfigsPion_directional.csh 1000 $i 10
    #source submitCondorPion_directional.csh 1000 $i 5
done
