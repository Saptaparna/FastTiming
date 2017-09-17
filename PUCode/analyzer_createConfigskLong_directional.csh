#!/bin/csh                                                                                                                                                                         
cp ConfFile_cfg_kLong_directional_ArabellaSample_PU_Template.py ConfFile_cfg_kLong_directional_ArabellaSample_PU_Part${3}.py
sed -i "s/NPART1/$1/g" ConfFile_cfg_kLong_directional_ArabellaSample_PU_Part${3}.py
sed -i "s/NPART2/$2/g" ConfFile_cfg_kLong_directional_ArabellaSample_PU_Part${3}.py
sed -i "s/NPART/$3/g" ConfFile_cfg_kLong_directional_ArabellaSample_PU_Part${3}.py
