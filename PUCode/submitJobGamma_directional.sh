#!/bin/sh

max=30
for i in `seq 7 $max`
do
    source analyzer_createConfigsGamma_directional.csh $[2*i-1] $[2*i] $i
    cmsRun ConfFile_cfg_Gamma_directional_ArabellaSample_PU_Part${i}.py 
done

