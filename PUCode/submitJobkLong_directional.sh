#!/bin/sh

max=30
for i in `seq 4 $max`
do
    source analyzer_createConfigskLong_directional.csh $[2*i-1] $[2*i] $i
    cmsRun ConfFile_cfg_kLong_directional_ArabellaSample_PU_Part${i}.py 
done

