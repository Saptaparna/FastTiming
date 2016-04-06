#!/bin/sh

max=10
for i in `seq 1 $max`
do
    source createConfigs.csh 1000 $i 35
    source submitCondor.csh 1000 $i 35
    source createConfigs.csh 1000 $i 50
    source submitCondor.csh 1000 $i 50
    source createConfigs.csh 1000 $i 65
    source submitCondor.csh 1000 $i 65
    source createConfigs.csh 1000 $i 80
    source submitCondor.csh 1000 $i 80
done
