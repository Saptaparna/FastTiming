#!/bin/csh                                                                                                                                                                         
cp SinglePhotonGun35_BatchSubmitTemplate_cfg.py SinglePhotonGun35_BatchSubmit_NEVENT${1}_NPART${2}_ENERGY${3}_cfg.py
sed -i "s/NEVENT/$1/g" SinglePhotonGun35_BatchSubmit_NEVENT${1}_NPART${2}_ENERGY${3}_cfg.py
sed -i "s/NPART/$2/g" SinglePhotonGun35_BatchSubmit_NEVENT${1}_NPART${2}_ENERGY${3}_cfg.py
sed -i "s/ENERGY/$3/g" SinglePhotonGun35_BatchSubmit_NEVENT${1}_NPART${2}_ENERGY${3}_cfg.py
