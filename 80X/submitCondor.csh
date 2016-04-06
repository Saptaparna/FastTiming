#!/bin/tcsh
echo $PWD

echo ${1}
echo ${2}
echo ${3}

cat>Job_${1}_${2}_${3}.csh<<EOF
#!/bin/tcsh
tar -xzf CMSSW_8_1_X_2016-03-27-2300.tar.gz
cd CMSSW_8_1_X_2016-03-27-2300
scram b ProjectRename
source /cvmfs/cms.cern.ch/cmsset_default.csh
cmsenv
cd -
cmsRun  SinglePhotonGun35_BatchSubmit_NEVENT${1}_NPART${2}_ENERGY${3}_cfg.py >& SinglePhotonGun35_BatchSubmit_NEVENT${1}_NPART${2}_ENERGY${3}.txt
EOF

chmod 775 $PWD/Job_${1}_${2}_${3}.csh

cat>condor_${1}_${2}_${3}.jdl<<EOF
universe = vanilla
Executable = Job_${1}_${2}_${3}.csh
Requirements = OpSys == "LINUX" && (Arch != "DUMMY" )
request_disk = 10000000
request_memory = 8100
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
transfer_input_files = /uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_8X/CMSSW_8_1_X_2016-03-27-2300/src/RecoLocalCalo/HGCalRecProducers/test/CMSSW_8_1_X_2016-03-27-2300.tar.gz, /uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_8X/CMSSW_8_1_X_2016-03-27-2300/src/RecoLocalCalo/HGCalRecProducers/test/Job_${1}_${2}_${3}.csh, /uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_8X/CMSSW_8_1_X_2016-03-27-2300/src/RecoLocalCalo/HGCalRecProducers/test/SinglePhotonGun35_BatchSubmit_NEVENT${1}_NPART${2}_ENERGY${3}_cfg.py
notification = Never
Output = CondorJobs/STDOUT_${1}_${2}_${3}.stdout
Error = CondorJobs/STDERR_${1}_${2}_${3}.stderr
Log = CondorJobs/LOG_${1}_${2}_${3}.log
Queue 1
EOF

condor_submit condor_${1}_${2}_${3}.jdl
