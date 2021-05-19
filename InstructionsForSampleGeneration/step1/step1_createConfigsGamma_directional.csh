#!/bin/csh                                                                                                                                                                         
cp step1_SingleGammaGun_BatchSubmitTemplate_cfg.py step1_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg.py
sed -i "s/NEVENT/$1/g" step1_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg.py
sed -i "s/NPART/$2/g" step1_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg.py
sed -i "s/Energy/$3/g" step1_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg.py


#cp step1_SingleGammaGun_BatchSubmitTemplate_cfg_eta.py step1_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg_eta.py
#sed -i "s/NEVENT/$1/g" step1_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg_eta.py
#sed -i "s/NPART/$2/g" step1_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg_eta.py
#sed -i "s/E/$3/g" step1_SingleGammaGun_BatchSubmit_NEVENT${1}_NPART${2}_En${3}_cfg_eta.py
