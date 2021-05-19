#!/bin/csh                                                                                                                                                                         
cp step2_SingleGammaGun_BatchSubmitTemplate_cfg_directional.py step2_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg_directional.py
sed -i "s/NEVENT/$1/g" step2_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg_directional.py
sed -i "s/NPART/$2/g" step2_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg_directional.py
sed -i "s/Energy/$3/g" step2_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg_directional.py

#cp step2_SingleGammaGun_BatchSubmitTemplate_cfg_directional_eta.py step2_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg_directional_eta.py
#sed -i "s/NEVENT/$1/g" step2_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg_directional_eta.py
#sed -i "s/NPART/$2/g" step2_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg_directional_eta.py
#sed -i "s/Energy/$3/g" step2_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg_directional_eta.py
