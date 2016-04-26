#!/bin/csh                                                                                                                                                                         
cp SinglePhotonGun_BatchSubmitTemplate_cfg.py SinglePhotonGun_BatchSubmit_NEVENT${1}_NPART${2}_ENERGY${3}_cfg.py
sed -i "s/NEVENT/$1/g" SinglePhotonGun_BatchSubmit_NEVENT${1}_NPART${2}_ENERGY${3}_cfg.py
sed -i "s/NPART/$2/g" SinglePhotonGun_BatchSubmit_NEVENT${1}_NPART${2}_ENERGY${3}_cfg.py
sed -i "s/ENERGY/$3/g" SinglePhotonGun_BatchSubmit_NEVENT${1}_NPART${2}_ENERGY${3}_cfg.py
